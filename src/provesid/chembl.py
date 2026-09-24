"""
ChEMBL Database Interface

This module provides access to the ChEMBL SQLite database for querying chemical compounds,
structures, and properties. ChEMBL is a manually curated database of bioactive molecules
with drug-like properties maintained by EMBL-EBI.

Database tables accessed:

- molecule_dictionary: Primary compound information (ChEMBL ID, names, max_phase, drug
                       classifications, approval status, administration routes)
- molecule_hierarchy: Parent-salt-metabolite relationships for compounds and pro-drugs
- compound_structures: Chemical structures (SMILES, InChI, InChIKey)
- compound_properties: Physicochemical properties (MW, ALogP, HBA, HBD, PSA, etc.)
- molecule_synonyms: Alternative names and synonyms
- chembl_id_lookup: ChEMBL ID to internal ID mappings
- pesticide_classification: Pesticide mechanism classifications (FRAC, HRAC, IRAC)
- pesticide_class_mapping: Links compounds to pesticide classifications

Those eight tables are the whole of what PROVESID reads.  A full ChEMBL release
is ~30 GB across 74 tables, the rest of which is bioactivity data this package
never opens, so [`CheMBL.compact`][provesid.chembl.CheMBL.compact] builds an
extract holding only the eight -- about 2.4 GiB, with identical results.  See
``CheMBL.compact`` for details.

A first download builds that extract on the way in and keeps only it
(``CheMBL(source="sqlite")``, the default), so a machine that has never had
ChEMBL installs 2.4 GiB rather than 27.7 GiB.  ``CheMBL(source="mysql")``
builds the same extract from ChEMBL's 2.1 GB MySQL dump instead, read as a
stream, so the 27.7 GiB release is never written and ~4.5 GiB of free disk is
enough.  ``CheMBL(source="full")`` keeps the whole release, for anyone who
wants the other 66 tables.

ChEMBL documents every table of release N at
https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_<N>/schema_documentation.html,
where N is ``CheMBL().release``.
"""

import os
import re
import glob
import gzip
import zlib
import random
import sqlite3
import hashlib
import tarfile
import logging
import datetime
from typing import Optional, Dict, List, Any, Sequence, Tuple
from urllib.parse import urlsplit
import requests
from tqdm import tqdm
from .datasets import DownloadError, download_file
from .mysqldump import CreateTable, DumpFormatError, read_statements, sqlite_affinity
from .sqlite_client import SQLiteClient
from .utils import user_dataset_path


class ChEMBLError(Exception):
    """
    Custom exception for ChEMBL database errors.

    Raised for a database that is missing and may not be fetched, or that
    cannot be built or verified.

    Examples:
        >>> CheMBL.compact(CheMBL().db_path)     # the installed file is an extract
        Traceback (most recent call last):
        ...
        provesid.chembl.ChEMBLError: ... is already a PROVESID extract; nothing to compact. ...
    """
    pass


class _CountingReader:
    """
    A binary file wrapper that reports how many bytes have been read.

    ``tarfile``'s streaming mode only ever calls ``read``, so this is all it
    takes to drive a progress bar over the compressed bytes of an archive
    that is being decompressed on the fly.
    """

    def __init__(self, handle, on_read):
        self._handle = handle
        self._on_read = on_read

    def read(self, size: int = -1) -> bytes:
        data = self._handle.read(size)
        self._on_read(len(data))
        return data


class CheMBL(SQLiteClient):
    """
    Interface to the ChEMBL SQLite database for chemical compound queries.

    The ChEMBL database contains manually curated bioactive compounds with drug-like
    properties. This class provides methods to search compounds by various identifiers
    and retrieve structural and property information.

    Parameters
    ----------
    db_name : str, optional
        Name of the SQLite database file.  When omitted (the default), the newest
        ``chembl_*.db`` already present in the data directory is reused, and if
        there is none the release advertised by the EBI ``latest/`` directory is
        downloaded and named after the archive (e.g. ``chembl_37.db``).
    auto_download : bool, optional
        If True, automatically download database if not found (default: True)
    db_url : str, optional
        Custom URL for database download.  By default the URL is resolved from
        the EBI ``latest/`` directory listing at download time.
    source : {'sqlite', 'mysql', 'full'}, optional
        How a database that is not yet on disk is acquired (default:
        ``'sqlite'``).  ``'sqlite'`` downloads the release archive, builds the
        PROVESID extract from it and keeps only that -- 2.4 GiB installed.
        ``'mysql'`` builds the same extract from the 2.1 GB MySQL dump without
        ever writing the full release.  ``'full'`` keeps the whole 27.7 GiB
        release.  Ignored when a database is already present: no route
        touches what is on disk.

    Attributes
    ----------
    path : str
        Path to the data directory
    db_path : str
        Full path to the SQLite database file
    db_url : str or None
        Archive URL, once resolved.  ``None`` while an existing local database
        makes a download unnecessary.
    release : int or None
        ChEMBL release number parsed from the database filename.
    source : str
        The acquisition route this instance was constructed with.
    is_compact : bool
        True when the open database is a PROVESID extract built by
        [`compact`][provesid.chembl.CheMBL.compact] rather than a full ChEMBL release.
    provenance : dict or None
        What the extract was built from -- release, source file, build time,
        PROVESID version and per-table row counts.  None for a full release.
    conn : sqlite3.Connection
        This thread's SQLite database connection
    cursor : sqlite3.Cursor
        This thread's cursor for queries

    Examples
    --------
    >>> with CheMBL() as chembl:
    ...     compound = chembl.search_by_chembl_id('CHEMBL25')  # Aspirin
    ...     props = chembl.get_properties(compound['molregno'])
    >>> compound['pref_name'], compound['molregno']
    ('ASPIRIN', 1280)
    >>> print(f"MW: {props['mw_freebase']}")
    MW: 180.16

    Notes
    -----
    The database is large and grows with every release: ChEMBL 37 is ~5.8 GB
    compressed and ~30 GB once extracted. Initial setup downloads and extracts it
    from the EMBL-EBI FTP server.

    Most of that is never read.  [`compact`][provesid.chembl.CheMBL.compact]
    builds a ~2.4 GiB extract holding only the eight tables this class queries,
    answering identically, and a first download builds it on the way in::

        CheMBL()                             # downloads 5.8 GB, installs 2.4 GiB
        CheMBL(source="mysql")               # downloads 2.1 GB, installs 2.4 GiB
        CheMBL(source="full")                # ...or keeps the whole 27.7 GiB
        CheMBL.compact(remove_source=True)   # shrink a release already on disk

    Only the download route is affected: a release already on disk is opened as
    it is, and shrinking it is [`compact`][provesid.chembl.CheMBL.compact]'s
    job rather than something a constructor should do to 27.7 GB unasked.  A
    later ``CheMBL()`` opens the extract in preference to a full release of the
    same number.

    ``latest/`` is a moving directory: it holds only the current release, so the
    archive name changes with every ChEMBL release.  The version is therefore
    resolved from the directory listing rather than pinned, with
    `DEFAULT_DB_URL` as the fallback when the listing cannot be read.

    Connection handling comes from
    [`SQLiteClient`][provesid.sqlite_client.SQLiteClient]: use the class as a context
    manager, or call [`close`][provesid.sqlite_client.SQLiteClient.close] when
    finished, and query it from as many threads as you like --- each gets its
    own connection.
    """

    LATEST_DIR_URL = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/"

    FALLBACK_RELEASE = 37
    """Release used when the ``latest/`` listing cannot be read."""

    DEFAULT_DB_URL = f"{LATEST_DIR_URL}chembl_{FALLBACK_RELEASE}_sqlite.tar.gz"

    # ── Acquisition ───────────────────────────────────────────────────────────

    SOURCES: Tuple[str, ...] = ("sqlite", "mysql", "full")
    """How a missing database is acquired, and what is kept afterwards.

    ``sqlite`` downloads the 5.8 GB release archive, extracts the 27.7 GiB
    database, builds the extract from it and deletes the release; ``full``
    stops after the extraction and keeps all 74 tables.  Both transfer the
    same bytes -- the choice is what stays on disk, not what is fetched.

    ``mysql`` downloads ChEMBL's 2.1 GB MySQL dump instead and reads the
    extract's eight tables straight out of it
    ([`build_from_mysql_dump`][provesid.chembl.CheMBL.build_from_mysql_dump]),
    so the 27.7 GiB release is never written at all.  It ends with the same
    extract as ``sqlite``, for less than half the transfer and a seventh of the
    free disk.
    """

    _ARCHIVE_FORMATS: Dict[str, str] = {
        "sqlite": "sqlite", "full": "sqlite", "mysql": "mysql",
    }
    """Which of ChEMBL's archives each route downloads: the ``NN`` in
    ``chembl_NN_<format>.tar.gz``.
    """

    _COMPACT_HINT_BYTES = 1024 ** 3
    """A full release smaller than this is not worth suggesting
    [`compact`][provesid.chembl.CheMBL.compact] for -- and in the tests, the
    miniature databases are far below it.
    """

    # ── Compaction ────────────────────────────────────────────────────────────

    COMPACT_SUFFIX = "_provesid"
    """Filename marker distinguishing a PROVESID extract from a full release:
    ``chembl_37.db`` is the 30 GB original, ``chembl_37_provesid.db`` the
    extract [`compact`][provesid.chembl.CheMBL.compact] builds from it.
    """

    COMPACT_SCHEMA_VERSION = 1
    """Bumped whenever
    [`COMPACT_TABLES`][provesid.chembl.CheMBL.COMPACT_TABLES] or the kept
    columns change, so an extract built by an older PROVESID can be recognised
    as incomplete instead of failing later with ``no such table``.
    """

    COMPACT_TABLES: Dict[str, Dict[str, Any]] = {
        "molecule_dictionary": {"columns": None, "where": None},
        "compound_structures": {
            "columns": ("molregno", "canonical_smiles", "standard_inchi",
                        "standard_inchi_key"),
            "where": None,
        },
        "compound_properties": {"columns": None, "where": None},
        "molecule_synonyms": {
            "columns": ("molregno", "syn_type", "molsyn_id", "synonyms"),
            "where": None,
        },
        "molecule_hierarchy": {
            "columns": ("molregno", "parent_molregno", "active_molregno"),
            "where": None,
        },
        "chembl_id_lookup": {
            "columns": ("chembl_id", "entity_type", "entity_id", "status",
                        "last_active"),
            "where": "entity_type = 'COMPOUND'",
        },
        "pesticide_classification": {"columns": None, "where": None},
        "pesticide_class_mapping": {"columns": None, "where": None},
    }
    """The only tables this package reads, with the columns kept for each.

    ``columns=None`` keeps every column.  ``where`` restricts the rows: the
    ChEMBL id lookup carries an entry for every entity type (assays, targets,
    documents), and
    [`chembl_id_to_molregno`][provesid.chembl.CheMBL.chembl_id_to_molregno]
    only ever asks about compounds.  ``molfile`` is deliberately absent from
    ``compound_structures``: it is a quarter of the whole database, nothing in
    PROVESID consumes it, and a MOL block is reconstructible from the SMILES
    with RDKit.
    """

    COMPACT_INDEXES: Tuple[Tuple[str, str], ...] = (
        ("ix_lookup_chembl_id", "chembl_id_lookup(chembl_id)"),
        ("ix_md_molregno", "molecule_dictionary(molregno)"),
        ("ix_md_chembl_id", "molecule_dictionary(chembl_id)"),
        ("ix_md_pref_lower", "molecule_dictionary(lower(pref_name))"),
        ("ix_cs_molregno", "compound_structures(molregno)"),
        ("ix_cs_inchikey", "compound_structures(standard_inchi_key)"),
        ("ix_cs_smiles", "compound_structures(canonical_smiles)"),
        ("ix_cs_inchi", "compound_structures(standard_inchi)"),
        ("ix_cp_molregno", "compound_properties(molregno)"),
        ("ix_ms_molregno", "molecule_synonyms(molregno)"),
        ("ix_ms_syn_lower", "molecule_synonyms(lower(synonyms))"),
        ("ix_mh_molregno", "molecule_hierarchy(molregno)"),
        ("ix_pc_class", "pesticide_classification(pest_class_id)"),
        ("ix_pcm_molregno", "pesticide_class_mapping(molregno)"),
        ("ix_pcm_class", "pesticide_class_mapping(pest_class_id)"),
    )
    """Indexes built on the extract, one per lookup the package performs.

    ChEMBL's own indexes are not copied: most of them serve range queries on
    ``compound_properties`` (``alogp``, ``psa``, ``rtb``, …) that this package
    never issues.  The two ``lower(...)`` expression indexes are new, and are
    what makes [`search_by_name`][provesid.chembl.CheMBL.search_by_name] an
    index lookup rather than a scan of all 2.9 M rows: 743 ms per exact lookup
    on a full release, 10 µs here. They cost a few tens of MB.
    """

    PROVENANCE_TABLE = "provesid_provenance"
    """Table recording where an extract came from; absent in a full release."""

    _VERIFY_SAMPLE = 500
    """How many compounds `_verify_compact` re-reads from both databases
    before an extract is trusted enough to delete its source.
    """

    def __init__(
        self,
        db_name: Optional[str] = None,
        auto_download: bool = True,
        db_url: Optional[str] = None,
        data_dir: Optional[str] = None,
        db_path: Optional[str] = None,
        redownload: bool = False,
        source: str = "sqlite",
    ):
        """
        Initialize ChEMBL database interface.

        Parameters
        ----------
        db_name : str, optional
            Database filename.  When omitted, the newest ``chembl_*.db`` in the
            data directory is reused, or the name is derived from the archive
            that is about to be downloaded (e.g. ``chembl_37.db``).
        auto_download : bool, optional
            Auto-download if database missing (default: True)
        db_url : str, optional
            Custom download URL.  Default: resolved from the EBI ``latest/``
            directory listing at download time.
        data_dir : str, optional
            Directory to store the database when ``db_path`` is not provided.
        db_path : str, optional
            Full path to a database file. Overrides ``db_name``/``data_dir``.
        redownload : bool, optional
            If True, force re-download when ``auto_download`` is enabled.
        source : {'sqlite', 'mysql', 'full'}, optional
            How a download is made, and what it leaves on disk (default:
            ``'sqlite'``).  With ``'sqlite'`` the release is compacted into the
            PROVESID extract as the last step of the download and the 27.7 GiB
            original is deleted, so installing ChEMBL costs 2.4 GiB but needs
            33.4 GiB free on the way.  ``'mysql'`` downloads the 2.1 GB MySQL
            dump and builds the same extract from it directly
            ([`build_from_mysql_dump`][provesid.chembl.CheMBL.build_from_mysql_dump]),
            needing ~4.5 GiB free.  With ``'full'`` the whole release is kept.
             Only consulted when a download actually happens.

        Raises
        ------
        FileNotFoundError
            If database not found and auto_download is False
        ValueError
            If ``source`` is not one of [`SOURCES`][provesid.chembl.CheMBL.SOURCES].
        ChEMBLError
            If database connection or validation fails
        """
        self.logger = logging.getLogger(__name__)
        self.db_url = db_url
        self.source = self._validate_source(source)

        if db_path is not None:
            self.db_path = os.path.abspath(os.path.expanduser(db_path))
            self.path = os.path.dirname(self.db_path)
        else:
            self.path = data_dir or user_dataset_path()
            self.db_path = (
                os.path.join(self.path, db_name) if db_name is not None else None
            )

        # Ensure data directory exists
        os.makedirs(self.path, exist_ok=True)

        # With no explicit filename, any release already on disk is good enough —
        # re-downloading tens of GB just because the release number moved is not.
        if self.db_path is None and not redownload:
            self.db_path = self._find_local_database()
            if self.db_path is not None:
                self.logger.info("Using existing ChEMBL database: %s", self.db_path)

        needs_download = (
            redownload or self.db_path is None or not os.path.exists(self.db_path)
        )

        if needs_download:
            if auto_download:
                self.db_url = self.db_url or self.resolve_latest_db_url(
                    archive=self._ARCHIVE_FORMATS[self.source]
                )
                if self.db_path is None:
                    self.db_path = os.path.join(
                        self.path, self._db_name_from_url(self.db_url)
                    )
                if redownload and os.path.exists(self.db_path):
                    self.logger.info(
                        "Forced ChEMBL redownload requested for: %s", self.db_path
                    )
                else:
                    self.logger.info(f"Database not found at {self.db_path}, downloading...")
                self.download_database(url=self.db_url, force=redownload)
                superseded = [
                    p
                    for p in glob.glob(os.path.join(self.path, "chembl_*.db"))
                    if os.path.abspath(p) != os.path.abspath(self.db_path)
                ]
                if superseded:
                    self.logger.info(
                        "Older ChEMBL database(s) left in place and no longer used "
                        "(delete to reclaim tens of GB each): %s",
                        ", ".join(sorted(superseded)),
                    )
            else:
                if self.db_path is None:
                    self.db_path = os.path.join(
                        self.path,
                        self._db_name_from_url(self.db_url or self.DEFAULT_DB_URL),
                    )
                raise FileNotFoundError(
                    f"ChEMBL database not found at {self.db_path}. "
                    f"Set auto_download=True or manually download from "
                    f"{self.db_url or self.LATEST_DIR_URL}"
                )

        release_match = re.search(r"chembl_(\d+)", os.path.basename(self.db_path))
        self.release = int(release_match.group(1)) if release_match else None

        # An extract carries its own provenance; a full release has none.
        self.provenance = self.read_provenance(self.db_path)
        self.is_compact = self.provenance is not None

        # Connect to the database.  One connection per thread, released by
        # close() or by leaving a ``with`` block --- see
        # ``SQLiteClient``.
        try:
            self._open_database(self.db_path)
            self.logger.info(f"Connected to ChEMBL database at {self.db_path}")
        except sqlite3.Error as e:
            raise ChEMBLError(f"Failed to connect to database: {str(e)}")

        if self.is_compact:
            self._warn_if_extract_is_stale()
        else:
            self._suggest_compacting_a_full_release()

    @classmethod
    def _validate_source(cls, source: str) -> str:
        """
        Check an acquisition route name, and say what is wrong with it.

        Parameters
        ----------
        source : str
            The ``source`` argument as given.

        Returns
        -------
        str
            The same value, once it is known to be a route.

        Raises
        ------
        ValueError
            If the name is not in [`SOURCES`][provesid.chembl.CheMBL.SOURCES].
        """
        if source in cls.SOURCES:
            return source
        options = ", ".join(repr(name) for name in cls.SOURCES)
        raise ValueError(
            f"CheMBL(source={source!r}) is not a download route. "
            f"Use one of {options}."
        )

    def _suggest_compacting_a_full_release(self) -> None:
        """
        Mention [`compact`][provesid.chembl.CheMBL.compact] when a large full
        release has been opened.

        The constructor will not shrink a database the user already has --
        deleting 27 GB is asked for, not assumed -- so the only thing left to
        do about it is to say that the option exists, once, where the user is
        already looking.  Silent below `_COMPACT_HINT_BYTES`, since
        there is then nothing worth reclaiming.
        """
        try:
            size = os.path.getsize(self.db_path)
        except OSError:  # pragma: no cover - the file was just opened
            return
        if size < self._COMPACT_HINT_BYTES:
            return
        self.logger.info(
            "This is a full ChEMBL release (%.2f GB). PROVESID reads eight of "
            "its 74 tables; CheMBL.compact(remove_source=True) rebuilds it as "
            "~2.4 GiB with identical results.",
            size / 1e9,
        )

    def _warn_if_extract_is_stale(self) -> None:
        """
        Say so when the open extract predates the tables this version needs.

        An extract is a subset of ChEMBL, so a PROVESID that later reads a ninth
        table would meet ``no such table`` with no hint that the database is
        merely old.  Comparing the stored
        [`COMPACT_SCHEMA_VERSION`][provesid.chembl.CheMBL.COMPACT_SCHEMA_VERSION]
        turns that into an instruction.

        Notes
        -----
        Warns rather than raises: an older extract still answers every query the
        package made when it was built, and refusing to open it would be worse
        than saying what to do about it.
        """
        stored = self.provenance.get("schema_version", "0")
        built_for = int(stored) if stored.isdigit() else 0
        if built_for < self.COMPACT_SCHEMA_VERSION:
            self.logger.warning(
                "ChEMBL extract %s was built for PROVESID extract schema v%s; "
                "this version expects v%d. Rebuild it with "
                "CheMBL.compact(force=True) if a query fails with 'no such "
                "table'.",
                os.path.basename(self.db_path), stored, self.COMPACT_SCHEMA_VERSION,
            )
        else:
            self.logger.info(
                "Using the PROVESID ChEMBL extract (release %s, built %s).",
                self.provenance.get("release", "?"),
                self.provenance.get("built_at", "?"),
            )

    @classmethod
    def resolve_latest_db_url(cls, timeout: float = 30, archive: str = "sqlite") -> str:
        """
        Resolve the download URL of the current ChEMBL archive.

        ``latest/`` only ever holds the newest release, so its archive name
        (``chembl_NN_sqlite.tar.gz``) changes with every ChEMBL release.  This
        reads the directory listing and returns the highest ``NN`` it advertises.

        Parameters
        ----------
        timeout : float, optional
            Timeout in seconds for the listing request (default: 30)
        archive : {'sqlite', 'mysql'}, optional
            Which archive to resolve (default: ``'sqlite'``): the SQLite
            release, or the MySQL dump ``source="mysql"`` reads.

        Returns
        -------
        str
            URL of the newest archive, or the same archive of
            [`FALLBACK_RELEASE`][provesid.chembl.CheMBL.FALLBACK_RELEASE] when
            the listing cannot be read or names no archive.  Never raises — a
            stale pin is better than a hard failure.

        Examples
        --------
        >>> CheMBL.resolve_latest_db_url()                      # doctest: +SKIP
        'https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_37_sqlite.tar.gz'
        >>> CheMBL.resolve_latest_db_url(archive="mysql")       # doctest: +SKIP
        'https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_37_mysql.tar.gz'
        """
        logger = logging.getLogger(__name__)
        pattern = re.compile(rf"chembl_(\d+)_{re.escape(archive)}\.tar\.gz")
        fallback = f"{cls.LATEST_DIR_URL}chembl_{cls.FALLBACK_RELEASE}_{archive}.tar.gz"
        try:
            response = requests.get(cls.LATEST_DIR_URL, timeout=timeout)
            response.raise_for_status()
            releases = {int(m.group(1)) for m in pattern.finditer(response.text)}
        except requests.RequestException as exc:
            logger.warning(
                "Could not read the ChEMBL release listing at %s (%s); "
                "falling back to the pinned release chembl_%d",
                cls.LATEST_DIR_URL, exc, cls.FALLBACK_RELEASE,
            )
            return fallback

        if not releases:
            logger.warning(
                "No chembl_NN_%s.tar.gz found in the listing at %s; "
                "falling back to the pinned release chembl_%d",
                archive, cls.LATEST_DIR_URL, cls.FALLBACK_RELEASE,
            )
            return fallback

        release = max(releases)
        logger.info("Resolved current ChEMBL release: chembl_%d", release)
        return f"{cls.LATEST_DIR_URL}chembl_{release}_{archive}.tar.gz"

    @staticmethod
    def _db_name_from_url(url: str) -> str:
        """
        Derive the local database filename from an archive URL.

        ``.../chembl_37_sqlite.tar.gz`` yields ``chembl_37.db``, matching the
        ``.db`` member inside the archive.  URLs that carry no release number
        fall back to ``chembl.db``.
        """
        archive = os.path.basename(urlsplit(url).path)
        match = re.search(r"chembl_(\d+)", archive)
        return f"chembl_{match.group(1)}.db" if match else "chembl.db"

    def _find_local_database(self) -> Optional[str]:
        """
        Return the best ``chembl_*.db`` already present in the data directory.

        Used when no filename was requested, so an already-downloaded release is
        reused instead of re-downloading tens of GB for a release number that
        merely moved.

        "Best" is the highest release number, and at equal release numbers the
        PROVESID extract in preference to the full database.  Once
        [`compact`][provesid.chembl.CheMBL.compact] has run, a plain
        ``CheMBL()`` should open the 2.6 GB extract rather than the 30 GB
        original that may still sit beside it -- they answer identically, and
        one of them costs a tenth of the page cache.

        Returns
        -------
        str or None
            Path to the preferred database, or None if the data directory holds
            no ChEMBL database.

        Examples
        --------
        With ``chembl_37.db`` and ``chembl_37_provesid.db`` side by side, the
        extract is chosen; with ``chembl_36_provesid.db`` and ``chembl_37.db``,
        release 37 wins, because a newer release beats a smaller file.
        """
        found: List[Tuple[int, int, str]] = []
        for path in glob.glob(os.path.join(self.path, "chembl_*.db")):
            match = re.search(r"chembl_(\d+)", os.path.basename(path))
            release = int(match.group(1)) if match else -1
            found.append((release, 1 if self.is_compact_database(path) else 0, path))
        return max(found)[2] if found else None

    # ── Compaction ────────────────────────────────────────────────────────────

    @classmethod
    def is_compact_database(cls, db_path: str) -> bool:
        """
        Report whether a SQLite file is a PROVESID ChEMBL extract.

        An extract carries a
        [`PROVENANCE_TABLE`][provesid.chembl.CheMBL.PROVENANCE_TABLE]; a full
        ChEMBL release does not.  The check opens the file read-only and never
        raises for a missing or unreadable path — a file that cannot be read
        is, for this purpose, not an extract.

        Parameters
        ----------
        db_path : str
            Path to the SQLite file to inspect.

        Returns
        -------
        bool
            True if the file is an extract built by
            [`compact`][provesid.chembl.CheMBL.compact].

        Examples
        --------
        >>> CheMBL.is_compact_database("/no/such/file.db")
        False
        """
        if not os.path.exists(db_path):
            return False
        conn = None
        try:
            conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
            row = conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name = ?",
                (cls.PROVENANCE_TABLE,),
            ).fetchone()
            return row is not None
        except sqlite3.Error:
            return False
        finally:
            # `with sqlite3.connect(...)` manages a transaction, not the
            # connection, so closing has to be explicit -- and this runs once
            # per candidate file in _find_local_database.
            if conn is not None:
                conn.close()

    @classmethod
    def read_provenance(cls, db_path: str) -> Optional[Dict[str, str]]:
        """
        Read the provenance record written by
        [`compact`][provesid.chembl.CheMBL.compact].

        Parameters
        ----------
        db_path : str
            Path to a ChEMBL extract.

        Returns
        -------
        dict or None
            Every key/value pair from the provenance table, or None when the
            file is a full ChEMBL release rather than an extract.

        Examples
        --------
        >>> prov = CheMBL.read_provenance("chembl_37_provesid.db")  # doctest: +SKIP
        >>> prov["release"], prov["schema_version"]                 # doctest: +SKIP
        ('37', '1')
        """
        if not cls.is_compact_database(db_path):
            return None
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
        try:
            rows = conn.execute(
                f"SELECT key, value FROM {cls.PROVENANCE_TABLE}"
            ).fetchall()
        finally:
            conn.close()
        return {key: value for key, value in rows}

    @classmethod
    def compact_path_for(cls, source_path: str) -> str:
        """
        Derive the extract's filename from a full database's path.

        ``/data/chembl_37.db`` becomes ``/data/chembl_37_provesid.db``, so the
        extract sits beside the release it came from and still carries the
        release number that [`CheMBL`][provesid.chembl.CheMBL] parses out of
        the filename.

        Parameters
        ----------
        source_path : str
            Path to a full ChEMBL database.

        Returns
        -------
        str
            Path the extract should be written to.

        Examples
        --------
        >>> CheMBL.compact_path_for("/data/chembl_37.db")
        '/data/chembl_37_provesid.db'
        """
        directory, filename = os.path.split(os.path.abspath(source_path))
        stem, extension = os.path.splitext(filename)
        return os.path.join(directory, f"{stem}{cls.COMPACT_SUFFIX}{extension}")

    @classmethod
    def compact(
        cls,
        source_path: Optional[str] = None,
        dest_path: Optional[str] = None,
        *,
        data_dir: Optional[str] = None,
        keep_inchi: bool = True,
        remove_source: bool = False,
        force: bool = False,
    ) -> str:
        """
        Build a small ChEMBL extract holding only the tables PROVESID reads.

        A full ChEMBL release is about 30 GB across 74 tables.  This package
        opens eight of them
        ([`COMPACT_TABLES`][provesid.chembl.CheMBL.COMPACT_TABLES]) and reads
        no bioactivity data at all, so almost the whole file is dead weight.
         The extract copies those eight tables, drops the ``molfile`` column --
        a quarter of the entire database on its own, and consumed nowhere --
        keeps only the ``COMPOUND`` rows of ``chembl_id_lookup``, and rebuilds
        just the indexes the package's own queries need.

        Measured on ChEMBL 36: **29.74 GB to 2.60 GB in 31 seconds**, with every
        public method of this class returning the same compounds.

        The extract is written to a temporary file and verified against its
        source before it replaces anything, and ``remove_source`` deletes the
        original only after that verification passes.

        Parameters
        ----------
        source_path : str, optional
            Full ChEMBL database to read.  When omitted, the newest
            ``chembl_*.db`` in ``data_dir`` that is not already an extract.
        dest_path : str, optional
            Where to write the extract.  Defaults to
            [`compact_path_for`][provesid.chembl.CheMBL.compact_path_for] of
            the source, i.e. ``chembl_37.db`` produces
            ``chembl_37_provesid.db`` beside it.
        data_dir : str, optional
            Directory searched for the source and used for the default
            destination.  Defaults to the shared PROVESID dataset directory.
        keep_inchi : bool, optional
            Keep the ``standard_inchi`` column and its index (default: True).
            Dropping it saves a further ~1.0 GB, but
            [`search_by_inchi`][provesid.chembl.CheMBL.search_by_inchi] then
            has no column to match and ``Search`` loses the InChI it currently
            reads straight from ChEMBL.  Leave this alone unless disk is
            genuinely short.
        remove_source : bool, optional
            Delete the full database once the extract verifies (default:
            False).  This is the step that reclaims the ~27 GB; it is off by
            default because deleting 30 GB should be asked for, not assumed.
        force : bool, optional
            Overwrite an existing extract (default: False).

        Returns
        -------
        str
            Path to the extract.

        Raises
        ------
        FileNotFoundError
            If no full ChEMBL database can be found to compact.
        FileExistsError
            If the destination exists and ``force`` is False.
        ChEMBLError
            If the source is itself an extract, if it is missing a table the
            extract needs, or if verification fails.  A failed build leaves the
            source untouched and removes the partial extract.

        Examples
        --------
        Shrink whatever release is already on disk, then reclaim the space:

        >>> path = CheMBL.compact(remove_source=True)            # doctest: +SKIP
        >>> path                                                  # doctest: +SKIP
        '/home/me/.local/share/provesid/chembl_36_provesid.db'
        >>> CheMBL.compact_path_for("/data/chembl_37.db")
        '/data/chembl_37_provesid.db'

        A later ``CheMBL()`` picks up the extract automatically, because
        `_find_local_database` prefers it over a full release of the same
        number.

        Notes
        -----
        The source is opened read-only, so compacting is safe while other
        processes are reading the same file.
        """
        logger = logging.getLogger(__name__)

        source_path = cls._resolve_compact_source(source_path, data_dir)
        if cls.is_compact_database(source_path):
            raise ChEMBLError(
                f"{source_path} is already a PROVESID extract; nothing to compact. "
                "Pass the full chembl_NN.db, or delete the extract and rebuild it."
            )

        dest_path = os.path.abspath(
            os.path.expanduser(dest_path or cls.compact_path_for(source_path))
        )
        if os.path.exists(dest_path) and not force:
            raise FileExistsError(
                f"ChEMBL extract already exists at {dest_path}. "
                "Pass force=True to rebuild it."
            )

        tables = cls._compact_table_spec(keep_inchi)
        cls._require_source_tables(source_path, tables)

        temp_path = f"{dest_path}.tmp"
        for leftover in (temp_path,):
            if os.path.exists(leftover):
                os.remove(leftover)

        source_bytes = os.path.getsize(source_path)
        logger.info(
            "Compacting %s (%.2f GB) into %s",
            source_path, source_bytes / 1e9, dest_path,
        )

        try:
            row_counts = cls._build_compact(source_path, temp_path, tables, keep_inchi)
            cls._verify_compact(source_path, temp_path, tables, row_counts)
        except Exception:
            if os.path.exists(temp_path):
                os.remove(temp_path)
            raise

        os.replace(temp_path, dest_path)
        dest_bytes = os.path.getsize(dest_path)
        logger.info(
            "ChEMBL extract written: %.2f GB, %.1f%% smaller than the source.",
            dest_bytes / 1e9, 100 * (1 - dest_bytes / source_bytes),
        )

        if remove_source:
            os.remove(source_path)
            logger.info(
                "Removed the full ChEMBL database %s, reclaiming %.2f GB.",
                source_path, source_bytes / 1e9,
            )
        else:
            logger.info(
                "The full database is still at %s (%.2f GB). Delete it, or call "
                "compact(remove_source=True), to reclaim that space.",
                source_path, source_bytes / 1e9,
            )

        return dest_path

    @classmethod
    def _resolve_compact_source(
        cls, source_path: Optional[str], data_dir: Optional[str]
    ) -> str:
        """
        Find the full ChEMBL database
        [`compact`][provesid.chembl.CheMBL.compact] should read.

        Parameters
        ----------
        source_path : str or None
            An explicit path, used as-is when given.
        data_dir : str or None
            Directory to search when no path was given.

        Returns
        -------
        str
            Absolute path to a full ChEMBL database.

        Raises
        ------
        FileNotFoundError
            If the explicit path does not exist, or the directory holds no full
            release to compact.
        """
        if source_path is not None:
            source_path = os.path.abspath(os.path.expanduser(source_path))
            if not os.path.exists(source_path):
                raise FileNotFoundError(f"No ChEMBL database at {source_path}")
            return source_path

        directory = data_dir or user_dataset_path()
        candidates = [
            path
            for path in glob.glob(os.path.join(directory, "chembl_*.db"))
            if not cls.is_compact_database(path)
        ]
        if not candidates:
            raise FileNotFoundError(
                f"No full ChEMBL database found in {directory}. "
                "Download one with CheMBL() first, or pass source_path."
            )
        # Highest release number wins, as elsewhere; an unparsable name sorts last.
        def release_of(path: str) -> int:
            match = re.search(r"chembl_(\d+)", os.path.basename(path))
            return int(match.group(1)) if match else -1

        return os.path.abspath(max(candidates, key=release_of))

    @classmethod
    def _compact_table_spec(cls, keep_inchi: bool) -> Dict[str, Dict[str, Any]]:
        """
        Return [`COMPACT_TABLES`][provesid.chembl.CheMBL.COMPACT_TABLES] with
        ``keep_inchi`` applied.

        Parameters
        ----------
        keep_inchi : bool
            When False, ``standard_inchi`` is removed from the columns kept for
            ``compound_structures``.

        Returns
        -------
        dict
            A copy of the table specification; the class attribute is untouched.
        """
        spec = {name: dict(entry) for name, entry in cls.COMPACT_TABLES.items()}
        if not keep_inchi:
            columns = spec["compound_structures"]["columns"]
            spec["compound_structures"]["columns"] = tuple(
                column for column in columns if column != "standard_inchi"
            )
        return spec

    @staticmethod
    def _require_source_tables(
        source_path: str, tables: Dict[str, Dict[str, Any]]
    ) -> None:
        """
        Check the source holds every table and column the extract needs.

        Failing here, before anything is written, gives the caller the missing
        name rather than a partial extract and an opaque SQL error.

        Parameters
        ----------
        source_path : str
            Full ChEMBL database.
        tables : dict
            Table specification from `_compact_table_spec`.

        Raises
        ------
        ChEMBLError
            If a table or a requested column is absent.
        """
        conn = sqlite3.connect(f"file:{source_path}?mode=ro", uri=True)
        try:
            present = {
                row[0]
                for row in conn.execute(
                    "SELECT name FROM sqlite_master WHERE type='table'"
                )
            }
            missing = [name for name in tables if name not in present]
            if missing:
                raise ChEMBLError(
                    f"{source_path} is not a full ChEMBL database: it is missing "
                    f"{', '.join(sorted(missing))}."
                )
            for name, entry in tables.items():
                if entry["columns"] is None:
                    continue
                available = {
                    row[1] for row in conn.execute(f"PRAGMA table_info({name})")
                }
                absent = [c for c in entry["columns"] if c not in available]
                if absent:
                    raise ChEMBLError(
                        f"Table {name} in {source_path} is missing column(s) "
                        f"{', '.join(absent)}."
                    )
        finally:
            conn.close()

    @classmethod
    def _build_compact(
        cls,
        source_path: str,
        dest_path: str,
        tables: Dict[str, Dict[str, Any]],
        keep_inchi: bool,
    ) -> Dict[str, int]:
        """
        Copy the wanted tables into a new database and index them.

        The source is attached read-only and each table is produced by a single
        ``CREATE TABLE ... AS SELECT``, which is both the simplest formulation
        and the fastest: SQLite streams the rows without a Python round-trip.

        Parameters
        ----------
        source_path : str
            Full ChEMBL database.
        dest_path : str
            File to create.  Must not already exist.
        tables : dict
            Table specification from `_compact_table_spec`.
        keep_inchi : bool
            Recorded in the provenance table.

        Returns
        -------
        dict
            Row count per copied table.
        """
        logger = logging.getLogger(__name__)
        conn = sqlite3.connect(dest_path)
        try:
            # This database is rebuilt from scratch on any failure, so durability
            # during the build buys nothing and costs a great deal of time.
            conn.execute("PRAGMA journal_mode = OFF")
            conn.execute("PRAGMA synchronous = OFF")
            conn.execute(
                "ATTACH DATABASE ? AS source", (f"file:{source_path}?mode=ro",)
            )

            row_counts: Dict[str, int] = {}
            for name, entry in tables.items():
                columns = "*" if entry["columns"] is None else ", ".join(entry["columns"])
                where = f" WHERE {entry['where']}" if entry["where"] else ""
                conn.execute(
                    f"CREATE TABLE {name} AS SELECT {columns} FROM source.{name}{where}"
                )
                count = conn.execute(f"SELECT COUNT(*) FROM {name}").fetchone()[0]
                row_counts[name] = count
                logger.info("  copied %s: %s rows", name, f"{count:,}")

            cls._write_provenance(conn, source_path, tables, row_counts, keep_inchi)
            conn.commit()
            conn.execute("DETACH DATABASE source")
            cls._index_and_seal(conn, tables)
        finally:
            conn.close()
        return row_counts

    @classmethod
    def _index_and_seal(
        cls, conn: sqlite3.Connection, tables: Dict[str, Dict[str, Any]]
    ) -> None:
        """
        Finish an extract whose tables are filled: indexes, statistics, VACUUM.

        Shared by both ways of building an extract -- from a SQLite release
        ([`compact`][provesid.chembl.CheMBL.compact]) and from a MySQL dump
        ([`build_from_mysql_dump`][provesid.chembl.CheMBL.build_from_mysql_dump])
        -- so that the two cannot drift apart in anything a query would notice.

        Parameters
        ----------
        conn : sqlite3.Connection
            Open connection to the extract, with every table loaded and no
            other database attached.
        tables : dict
            Table specification actually applied; decides whether the
            ``standard_inchi`` index has a column to be built on.
        """
        logger = logging.getLogger(__name__)
        kept_columns = tables["compound_structures"]["columns"]
        for index_name, target in cls.COMPACT_INDEXES:
            if "standard_inchi)" in target and "standard_inchi" not in kept_columns:
                continue
            conn.execute(f"CREATE INDEX IF NOT EXISTS {index_name} ON {target}")
        logger.info("  built %d indexes", len(cls.COMPACT_INDEXES))
        conn.commit()
        conn.execute("ANALYZE main")
        conn.commit()
        conn.execute("VACUUM")

    @classmethod
    def _write_provenance(
        cls,
        conn: sqlite3.Connection,
        source_path: str,
        tables: Dict[str, Dict[str, Any]],
        row_counts: Dict[str, int],
        keep_inchi: bool,
        source_format: str = "sqlite",
    ) -> None:
        """
        Record what this extract was built from, and by what.

        An extract is a *subset*, so it goes stale differently from a copy: a
        later PROVESID that needs a ninth table has to be able to say so rather
        than fail with ``no such table``.
         [`COMPACT_SCHEMA_VERSION`][provesid.chembl.CheMBL.COMPACT_SCHEMA_VERSION]
        is what makes that possible, and the rest of the record makes a
        database on disk able to answer where it came from.

        Parameters
        ----------
        conn : sqlite3.Connection
            Open connection to the extract being built.
        source_path : str
            Full ChEMBL database the extract was read from.
        tables : dict
            Table specification actually applied.
        row_counts : dict
            Rows copied per table.
        keep_inchi : bool
            Whether ``standard_inchi`` was kept.
        source_format : {'sqlite', 'mysql'}, optional
            What ``source_path`` is: a ChEMBL SQLite release (default), or the
            MySQL dump archive
            [`build_from_mysql_dump`][provesid.chembl.CheMBL.build_from_mysql_dump]
            read.
        """
        from . import __version__

        release_match = re.search(r"chembl_(\d+)", os.path.basename(source_path))
        record = {
            "schema_version": str(cls.COMPACT_SCHEMA_VERSION),
            "provesid_version": __version__,
            "built_at": datetime.datetime.now(datetime.timezone.utc).isoformat(
                timespec="seconds"
            ),
            "release": release_match.group(1) if release_match else "",
            "source_database": os.path.basename(source_path),
            "source_format": source_format,
            "source_bytes": str(os.path.getsize(source_path)),
            "keep_inchi": "1" if keep_inchi else "0",
            "tables": ",".join(sorted(tables)),
        }
        for name, count in row_counts.items():
            record[f"rows.{name}"] = str(count)

        conn.execute(
            f"CREATE TABLE {cls.PROVENANCE_TABLE} "
            "(key TEXT PRIMARY KEY, value TEXT NOT NULL)"
        )
        conn.executemany(
            f"INSERT INTO {cls.PROVENANCE_TABLE} (key, value) VALUES (?, ?)",
            sorted(record.items()),
        )

    @classmethod
    def _verification_sample(cls, source: sqlite3.Connection) -> List[int]:
        """
        Pick molregnos spread across the whole table, cheaply.

        Sampling matters here: ``remove_source`` deletes 30 GB on the strength
        of this check, so the rows compared should not all come from one end of
        the table.  ``ORDER BY random()`` would be unbiased but scans 2.9 M rows;
        drawing random values from the observed molregno range and keeping the
        ones that exist costs a few hundred index seeks and covers the table
        evenly.  The seed is fixed so a failure can be reproduced.

        Parameters
        ----------
        source : sqlite3.Connection
            Open connection to the full database.

        Returns
        -------
        list of int
            Up to `_VERIFY_SAMPLE` molregnos that exist in
            ``compound_structures``.  Fewer when the table is smaller than that,
            in which case every row is returned.
        """
        low, high, total = source.execute(
            "SELECT MIN(molregno), MAX(molregno), COUNT(*) FROM compound_structures"
        ).fetchone()
        if not total:
            return []
        if total <= cls._VERIFY_SAMPLE:
            return [
                row[0]
                for row in source.execute("SELECT molregno FROM compound_structures")
            ]

        rng = random.Random(0)
        found: List[int] = []
        seen = set()
        # Generous attempt budget: molregnos are dense, so most draws hit, but
        # the loop must terminate even on a sparse table.
        for _ in range(cls._VERIFY_SAMPLE * 20):
            if len(found) >= cls._VERIFY_SAMPLE:
                break
            candidate = rng.randint(low, high)
            if candidate in seen:
                continue
            seen.add(candidate)
            row = source.execute(
                "SELECT molregno FROM compound_structures WHERE molregno = ?",
                (candidate,),
            ).fetchone()
            if row is not None:
                found.append(row[0])
        return found

    @classmethod
    def _verify_compact(
        cls,
        source_path: str,
        dest_path: str,
        tables: Dict[str, Dict[str, Any]],
        row_counts: Dict[str, int],
    ) -> None:
        """
        Prove the extract matches its source before anything is deleted.

        Three checks, cheapest first: SQLite's own structural check, a row count
        per table against the source under the same filter, and a sample of
        whole rows re-read from both databases and compared column by column.
        ``remove_source`` destroys 30 GB, so it is worth being sure.

        Parameters
        ----------
        source_path : str
            Full ChEMBL database.
        dest_path : str
            The extract just built.
        tables : dict
            Table specification actually applied.
        row_counts : dict
            Rows the build reported copying.

        Raises
        ------
        ChEMBLError
            On the first check that fails, naming what disagreed.
        """
        logger = logging.getLogger(__name__)
        source = sqlite3.connect(f"file:{source_path}?mode=ro", uri=True)
        dest = sqlite3.connect(f"file:{dest_path}?mode=ro", uri=True)
        source.row_factory = sqlite3.Row
        dest.row_factory = sqlite3.Row
        try:
            structural = dest.execute("PRAGMA quick_check").fetchone()[0]
            if structural != "ok":
                raise ChEMBLError(f"Extract failed SQLite's quick_check: {structural}")

            for name, entry in tables.items():
                where = f" WHERE {entry['where']}" if entry["where"] else ""
                expected = source.execute(
                    f"SELECT COUNT(*) FROM {name}{where}"
                ).fetchone()[0]
                actual = dest.execute(f"SELECT COUNT(*) FROM {name}").fetchone()[0]
                if actual != expected or actual != row_counts[name]:
                    raise ChEMBLError(
                        f"Row count mismatch in {name}: source has {expected:,}, "
                        f"extract has {actual:,}, build reported "
                        f"{row_counts[name]:,}."
                    )

            columns = tables["compound_structures"]["columns"]
            selection = ", ".join(columns)
            query = f"SELECT {selection} FROM compound_structures WHERE molregno = ?"
            compared = 0
            for molregno in cls._verification_sample(source):
                if (tuple(source.execute(query, (molregno,)).fetchone())
                        != tuple(dest.execute(query, (molregno,)).fetchone())):
                    raise ChEMBLError(
                        f"Extract disagrees with its source for molregno {molregno}."
                    )
                compared += 1
            logger.info(
                "  verified: quick_check ok, %d table counts match, %d compounds "
                "compared row for row.", len(tables), compared,
            )
        finally:
            source.close()
            dest.close()

    # ── Building the extract from the MySQL dump (source="mysql") ─────────────

    @classmethod
    def build_from_mysql_dump(
        cls,
        dump_path: str,
        dest_path: Optional[str] = None,
        *,
        keep_inchi: bool = True,
        remove_source: bool = False,
        force: bool = False,
    ) -> str:
        """
        Build the PROVESID extract straight from ChEMBL's MySQL dump.

        The same extract as [`compact`][provesid.chembl.CheMBL.compact] -- the
        eight tables in
        [`COMPACT_TABLES`][provesid.chembl.CheMBL.COMPACT_TABLES], the same
        columns, the same indexes -- read out of ``chembl_NN_mysql.tar.gz`` as
        it is decompressed, rather than out of a 27.7 GiB SQLite release that
        first has to be written to disk.  The other 66 tables stream past
        unparsed.

        ================  ================  =====================
        route             download          free disk at peak
        ================  ================  =====================
        ``sqlite``        5.8 GB            33.4 GiB
        ``mysql``         2.1 GB            ~4.5 GiB
        ================  ================  =====================

        Column types are taken from the dump's ``CREATE TABLE`` statements and
        mapped to the affinity SQLite would give them
        ([`provesid.mysqldump.sqlite_affinity`][provesid.mysqldump.sqlite_affinity]),
        so every value is stored the way
        [`compact`][provesid.chembl.CheMBL.compact] stores it: a ``max_phase``
        of ``4.0`` becomes the integer 4 in both.

        Parameters
        ----------
        dump_path : str
            The ``chembl_NN_mysql.tar.gz`` archive, as downloaded.  It is read
            as a stream and never extracted.
        dest_path : str, optional
            Where to write the extract.  Defaults to
            ``chembl_NN_provesid.db`` beside the archive, the name
            [`compact`][provesid.chembl.CheMBL.compact] would give the same release.
        keep_inchi : bool, optional
            Keep the ``standard_inchi`` column and its index (default: True).
            See [`compact`][provesid.chembl.CheMBL.compact].
        remove_source : bool, optional
            Delete the archive as soon as it has been read (default: False).
            That is before the indexes are built and the file is vacuumed,
            which is what keeps the peak at roughly the archive plus the
            extract, rather than the archive plus two copies of it.
        force : bool, optional
            Overwrite an existing extract (default: False).

        Returns
        -------
        str
            Path to the extract.

        Raises
        ------
        FileNotFoundError
            If ``dump_path`` does not exist.
        FileExistsError
            If the destination exists and ``force`` is False.
        ChEMBLError
            If the archive holds no ``.dmp`` or ``.sql`` file, if it cannot be
            decompressed, if a table or column the extract needs is missing, if
            one of the eight tables has no rows, or if a line of a needed table
            cannot be parsed
            ([`provesid.mysqldump.DumpFormatError`][provesid.mysqldump.DumpFormatError]).
             A failed build removes the partial extract.

        Examples
        --------
        With a dump downloaded by hand:

        >>> path = CheMBL.build_from_mysql_dump("chembl_37_mysql.tar.gz",
        ...                                     remove_source=True)   # doctest: +SKIP
        >>> chembl = CheMBL(db_path=path)                            # doctest: +SKIP

        ``CheMBL(source="mysql")`` downloads the dump and calls this.

        Notes
        -----
        There is no source database to check the result against, so the
        verification [`compact`][provesid.chembl.CheMBL.compact] performs is
        replaced by SQLite's ``quick_check``, the row counts, and the
        requirement that all eight tables were found and none is empty.
         Equivalence with the ``sqlite`` route is a property of the parser and
        was measured once, on real data, with
        [`extract_digest`][provesid.chembl.CheMBL.extract_digest].
        """
        logger = logging.getLogger(__name__)

        dump_path = os.path.abspath(os.path.expanduser(dump_path))
        if not os.path.exists(dump_path):
            raise FileNotFoundError(f"No ChEMBL MySQL dump at {dump_path}")
        if dest_path is None:
            dest_path = cls.compact_path_for(
                os.path.join(os.path.dirname(dump_path),
                             cls._db_name_from_url(dump_path))
            )
        dest_path = os.path.abspath(os.path.expanduser(dest_path))
        if os.path.exists(dest_path) and not force:
            raise FileExistsError(
                f"ChEMBL extract already exists at {dest_path}. "
                "Pass force=True to rebuild it."
            )

        tables = cls._compact_table_spec(keep_inchi)
        temp_path = f"{dest_path}.tmp"
        if os.path.exists(temp_path):
            os.remove(temp_path)

        dump_bytes = os.path.getsize(dump_path)
        logger.info(
            "Building the ChEMBL extract %s from the MySQL dump %s (%.2f GB)",
            dest_path, dump_path, dump_bytes / 1e9,
        )

        try:
            conn = sqlite3.connect(temp_path)
            try:
                conn.execute("PRAGMA journal_mode = OFF")
                conn.execute("PRAGMA synchronous = OFF")
                row_counts = cls._load_mysql_dump(conn, dump_path, tables)
                cls._write_provenance(conn, dump_path, tables, row_counts,
                                      keep_inchi, source_format="mysql")
                conn.commit()
                if remove_source:
                    os.remove(dump_path)
                    logger.info("  removed the dump %s (%.2f GB)",
                                dump_path, dump_bytes / 1e9)
                cls._index_and_seal(conn, tables)
            finally:
                conn.close()
            cls._verify_dump_extract(temp_path, tables, row_counts)
        except Exception:
            if os.path.exists(temp_path):
                os.remove(temp_path)
            raise

        os.replace(temp_path, dest_path)
        logger.info("ChEMBL extract written: %.2f GB.",
                    os.path.getsize(dest_path) / 1e9)
        return dest_path

    @classmethod
    def _load_mysql_dump(
        cls,
        conn: sqlite3.Connection,
        dump_path: str,
        tables: Dict[str, Dict[str, Any]],
    ) -> Dict[str, int]:
        """
        Stream the wanted tables out of a dump archive into ``conn``.

        The archive is opened in ``tarfile``'s streaming mode (``r|gz``), so it
        is decompressed once, front to back, and nothing is written but the
        rows kept.  Each wanted ``CREATE TABLE`` creates its table
        (`_create_table_from_dump`); each ``INSERT`` is projected and
        filtered by the statement that call returns.

        Parameters
        ----------
        conn : sqlite3.Connection
            The extract being built.
        dump_path : str
            ``chembl_NN_mysql.tar.gz``.
        tables : dict
            Table specification from `_compact_table_spec`.

        Returns
        -------
        dict
            Rows kept per table -- after the ``where`` filter, so
            ``chembl_id_lookup`` counts compounds only.

        Raises
        ------
        ChEMBLError
            If the archive is unreadable or holds no dump, or a table is
            missing, lacks a column, or is empty.
        """
        wanted = set(tables)
        dump_columns: Dict[str, Tuple[str, ...]] = {}
        insert_sql: Dict[str, str] = {}
        row_counts = {name: 0 for name in tables}
        dumps_read = 0

        try:
            with open(dump_path, "rb") as raw, tqdm(
                total=os.path.getsize(dump_path), unit="B", unit_scale=True,
                desc="Reading ChEMBL MySQL dump",
            ) as bar, tarfile.open(fileobj=_CountingReader(raw, bar.update),
                                   mode="r|gz") as tar:
                for member in tar:
                    if not (member.isfile()
                            and member.name.endswith((".dmp", ".sql"))):
                        continue
                    dumps_read += 1
                    for statement in read_statements(tar.extractfile(member), wanted):
                        name = statement.table
                        if isinstance(statement, CreateTable):
                            dump_columns[name] = tuple(c.name for c in statement.columns)
                            insert_sql[name] = cls._create_table_from_dump(
                                conn, statement, tables[name], dump_path
                            )
                            continue
                        rows = statement.rows
                        if (statement.columns is not None
                                and statement.columns != dump_columns[name]):
                            order = [statement.columns.index(c)
                                     for c in dump_columns[name]]
                            rows = [tuple(row[i] for i in order) for row in rows]
                        row_counts[name] += conn.executemany(
                            insert_sql[name], rows
                        ).rowcount
        except (DumpFormatError, tarfile.TarError, EOFError, zlib.error,
                gzip.BadGzipFile, UnicodeDecodeError) as exc:
            raise ChEMBLError(f"Could not read the MySQL dump {dump_path}: {exc}") from exc

        if not dumps_read:
            raise ChEMBLError(
                f"{dump_path} holds no .dmp or .sql file; is it ChEMBL's "
                "chembl_NN_mysql.tar.gz?"
            )
        missing = sorted(name for name in tables if name not in dump_columns)
        if missing:
            raise ChEMBLError(
                f"The MySQL dump {dump_path} has no CREATE TABLE for "
                f"{', '.join(missing)}."
            )
        # A table that exists but received nothing most likely means its INSERT
        # lines were written in a shape the reader did not recognise, and so
        # were skipped as another table's would be.  Every one of the eight is
        # populated in every ChEMBL release.
        empty = sorted(name for name, count in row_counts.items() if count == 0)
        if empty:
            raise ChEMBLError(
                f"The MySQL dump {dump_path} yielded no rows for "
                f"{', '.join(empty)}."
            )
        for name, count in row_counts.items():
            logging.getLogger(__name__).info("  loaded %s: %s rows", name, f"{count:,}")
        return row_counts

    @staticmethod
    def _create_table_from_dump(
        conn: sqlite3.Connection,
        statement: CreateTable,
        spec: Dict[str, Any],
        dump_path: str,
    ) -> str:
        """
        Create one extract table from its ``CREATE TABLE`` in the dump.

        Parameters
        ----------
        conn : sqlite3.Connection
            The extract being built.
        statement : CreateTable
            The table's definition as the dump gives it.
        spec : dict
            This table's entry in the table specification: which columns to
            keep (None for all) and which rows (``where``).
        dump_path : str
            Used in the error message only.

        Returns
        -------
        str
            The ``INSERT`` statement that takes one dump row, in the dump's
            column order, and stores the kept columns of it if it passes the
            filter.  Projection and filter are left to SQLite, which reuses the
            same ``where`` text [`compact`][provesid.chembl.CheMBL.compact]
            applies to a release.

        Raises
        ------
        ChEMBLError
            If a column the extract keeps is not in the dump.
        """
        name = statement.table
        available = tuple(column.name for column in statement.columns)
        kept = spec["columns"] or available
        absent = [column for column in kept if column not in available]
        if absent:
            raise ChEMBLError(
                f"Table {name} in {dump_path} is missing column(s) "
                f"{', '.join(absent)}."
            )
        affinity = {column.name: sqlite_affinity(column.type)
                    for column in statement.columns}
        definitions = ", ".join(f'"{column}" {affinity[column]}'.rstrip()
                                for column in kept)
        conn.execute(f"CREATE TABLE {name} ({definitions})")

        if kept == available and not spec["where"]:
            return f"INSERT INTO {name} VALUES ({', '.join('?' * len(kept))})"
        selected = ", ".join(f'"{column}"' for column in kept)
        bound = ", ".join(f'? AS "{column}"' for column in available)
        where = f" WHERE {spec['where']}" if spec["where"] else ""
        return (f"INSERT INTO {name} ({selected}) "
                f"SELECT {selected} FROM (SELECT {bound}){where}")

    @staticmethod
    def _verify_dump_extract(
        dest_path: str,
        tables: Dict[str, Dict[str, Any]],
        row_counts: Dict[str, int],
    ) -> None:
        """
        Check an extract built from a dump before it is moved into place.

        Parameters
        ----------
        dest_path : str
            The extract just built.
        tables : dict
            Table specification actually applied.
        row_counts : dict
            Rows the load reported keeping.

        Raises
        ------
        ChEMBLError
            If SQLite's ``quick_check`` fails, or a table holds a different
            number of rows than were loaded into it.
        """
        dest = sqlite3.connect(f"file:{dest_path}?mode=ro", uri=True)
        try:
            structural = dest.execute("PRAGMA quick_check").fetchone()[0]
            if structural != "ok":
                raise ChEMBLError(f"Extract failed SQLite's quick_check: {structural}")
            for name in tables:
                actual = dest.execute(f"SELECT COUNT(*) FROM {name}").fetchone()[0]
                if actual != row_counts[name]:
                    raise ChEMBLError(
                        f"Row count mismatch in {name}: loaded "
                        f"{row_counts[name]:,}, extract has {actual:,}."
                    )
        finally:
            dest.close()

    @classmethod
    def extract_digest(cls, db_path: str) -> Dict[str, Tuple[int, str]]:
        """
        Fingerprint each extract table: its row count and a content hash.

        Two extracts with equal digests hold the same rows, with the same
        values of the same SQLite types, whatever order the rows are stored
        in -- which is how the ``mysql`` route was checked against the
        ``sqlite`` route, and how anyone can repeat that check.  Columns are
        compared by name, so column order does not matter either.

        Parameters
        ----------
        db_path : str
            A database holding the tables of
            [`COMPACT_TABLES`][provesid.chembl.CheMBL.COMPACT_TABLES] -- an
            extract, or a full release.

        Returns
        -------
        dict
            ``{table: (row_count, hex_digest)}`` for each table in
            [`COMPACT_TABLES`][provesid.chembl.CheMBL.COMPACT_TABLES] present
            in the file.

        Examples
        --------
        >>> a = CheMBL.extract_digest("chembl_36_provesid.db")        # doctest: +SKIP
        >>> b = CheMBL.extract_digest("chembl_36_from_mysql.db")      # doctest: +SKIP
        >>> {t for t in a if a[t] != b[t]}                            # doctest: +SKIP
        set()

        Notes
        -----
        Reads every row of every table: about a minute on a real extract.
        The hash of each row is summed modulo 2**128, which is what makes the
        result independent of row order.  ``typeof`` is part of what is hashed,
        so an integer 180 and a real 180.0 do not collide.
        """
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
        try:
            present = {
                row[0] for row in conn.execute(
                    "SELECT name FROM sqlite_master WHERE type='table'"
                )
            }
            digests: Dict[str, Tuple[int, str]] = {}
            for name, spec in cls.COMPACT_TABLES.items():
                if name not in present:
                    continue
                available = [
                    row[1] for row in conn.execute(f"PRAGMA table_info({name})")
                ]
                columns = sorted(column for column in spec["columns"] or available
                                 if column in available)
                where = f" WHERE {spec['where']}" if spec["where"] else ""
                selection = ", ".join(f'typeof("{c}"), "{c}"' for c in columns)
                total, count = 0, 0
                for row in conn.execute(f"SELECT {selection} FROM {name}{where}"):
                    total += int.from_bytes(
                        hashlib.md5(repr(row).encode("utf-8")).digest(), "big"
                    )
                    count += 1
                digests[name] = (count, f"{total % 2 ** 128:032x}")
            return digests
        finally:
            conn.close()

    def download_database(self, url: Optional[str] = None, force: bool = False):
        """
        Download and extract ChEMBL SQLite database from EMBL-EBI FTP.

        Downloads the compressed tar.gz archive (~5.8 GB for release 37), extracts
        the SQLite database (~27.7 GiB), and validates its integrity by querying
        the molecule_dictionary table.

        With ``source="sqlite"`` (the default) the release is then compacted
        into the PROVESID extract ([`compact`][provesid.chembl.CheMBL.compact])
        and the full database is deleted, leaving ~2.4 GiB on disk and moving
        ``db_path`` onto the extract.  The two together never need more room
        than the extraction already did: the archive is deleted before the
        extract is built.

        With ``source="mysql"`` none of that happens: the 2.1 GB MySQL dump is
        downloaded instead and the same extract is read straight out of it
        ([`build_from_mysql_dump`][provesid.chembl.CheMBL.build_from_mysql_dump]),
        so there is no release to extract, validate or delete.

        The download is resumable. An interrupted transfer leaves a ``.part``
        file beside the archive and the next call continues from it, which
        matters more here than anywhere else in the package: this is the
        largest single download PROVESID makes, and it used to start again from
        zero.

        Parameters
        ----------
        url : str, optional
            Download URL.  Defaults to this instance's ``db_url``, or — when that
            was never set — the release resolved from the EBI ``latest/`` listing.
        force : bool, optional
            If True, re-download even if database exists (default: False)

        Raises
        ------
        ChEMBLError
            If download, extraction, or validation fails

        Notes
        -----
        A compaction that fails after a successful download does *not* raise:
        the full release is in place by then and answers every query, so the
        failure costs disk rather than function.  It is logged as a warning
        naming [`compact`][provesid.chembl.CheMBL.compact], and ``db_path``
        stays on the full release.

        Examples
        --------
        >>> chembl = CheMBL(auto_download=False)      # doctest: +SKIP
        >>> chembl.download_database(force=True)      # doctest: +SKIP

        Both lines are marked ``+SKIP`` deliberately: under
        ``pytest --doctest-modules`` this example would otherwise *execute*, and
        fetch ~5.8 GB into whatever ``db_path`` happens to be -- including over a
        compacted extract.
        """
        url = url or self.db_url or self.resolve_latest_db_url(
            archive=self._ARCHIVE_FORMATS[self.source]
        )
        self.db_url = url

        # Check if already exists
        if os.path.exists(self.db_path) and not force:
            self.logger.info(f"Database already exists at {self.db_path}")
            return

        if self.source == "mysql":
            self._install_from_mysql_dump(url)
            self.logger.info("ChEMBL database download and setup complete")
            return

        # The archive is kept beside the database rather than in a temporary
        # directory, so that an interrupted download's .part file is found and
        # resumed by the next call.  5.8 GB is far too much to fetch twice.
        archive_path = self.db_path + ".tar.gz"

        # The database is extracted under a temporary name and only moved onto
        # db_path once it has answered a query, so a failed extraction cannot
        # destroy the release that is already there -- the same contract the
        # other four downloads now keep through download_file's `verify`.
        staged_path = self.db_path + ".incoming"

        try:
            download_file(
                url,
                archive_path,
                description=f"ChEMBL archive ({os.path.basename(url)})",
                log=self.logger,
            )

            self.logger.info("Download complete. Extracting database...")
            self._extract_database(archive_path, staged_path)

            # Validate database integrity
            self.logger.info("Validating database integrity...")
            try:
                test_conn = sqlite3.connect(staged_path)
                test_cursor = test_conn.cursor()
                test_cursor.execute("SELECT COUNT(*) FROM molecule_dictionary")
                count = test_cursor.fetchone()[0]
                test_conn.close()
                self.logger.info(f"Database validated successfully. Contains {count:,} compounds.")
            except sqlite3.Error as e:
                raise ChEMBLError(f"Database validation failed: {str(e)}") from e

            os.replace(staged_path, self.db_path)
            # Before compacting, not after: the extract is built from the
            # database, so 5.8 GB of archive held on to here would raise the
            # peak by the size of a download nobody needs any more.
            os.remove(archive_path)

            if self.source == "sqlite":
                self._compact_after_download()

            self.logger.info("ChEMBL database download and setup complete")

        except DownloadError as e:
            # The .part file is left where it is: the next call resumes from
            # it rather than fetching 5.8 GB again.
            raise ChEMBLError(f"Download failed: {str(e)}") from e
        except Exception as e:
            # Extraction or validation failed, so the archive is suspect and
            # the half-built database is worthless. Both go; whatever was at
            # db_path before is untouched.
            for path in (staged_path, archive_path):
                if os.path.exists(path):
                    os.remove(path)
            if isinstance(e, ChEMBLError):
                raise
            raise ChEMBLError(f"Database setup failed: {str(e)}") from e

    def _install_from_mysql_dump(self, url: str) -> None:
        """
        Download ChEMBL's MySQL dump and build the extract from it.

        The whole of ``source="mysql"``.  The dump is downloaded to disk rather
        than parsed straight off the socket, because
        [`download_file`][provesid.datasets.download_file] can resume it: an
        interruption at 90% of 2.1 GB should cost the last 10%, not the whole
        transfer again.  That costs 2.1 GB of transient disk, against the 33.4
        GiB the ``sqlite`` route needs.

        The archive is deleted as soon as it has been read, and on any failure
        to build from it.  On success ``db_path`` moves onto the extract.

        Parameters
        ----------
        url : str
            URL of ``chembl_NN_mysql.tar.gz``.

        Raises
        ------
        ChEMBLError
            If the download fails -- leaving a ``.part`` file the next call
            resumes -- or the extract cannot be built, in which case the
            message names ``source="sqlite"`` as the route that does not
            depend on the dump's format.
        """
        archive_name = os.path.basename(urlsplit(url).path) or "chembl_mysql.tar.gz"
        archive_path = os.path.join(self.path, archive_name)
        extract_path = self.compact_path_for(self.db_path)

        try:
            download_file(
                url,
                archive_path,
                description=f"ChEMBL MySQL dump ({archive_name})",
                log=self.logger,
            )
        except DownloadError as e:
            raise ChEMBLError(f"Download failed: {str(e)}") from e

        try:
            self.build_from_mysql_dump(
                archive_path, extract_path, remove_source=True, force=True
            )
        except Exception as e:
            if os.path.exists(archive_path):
                os.remove(archive_path)
            raise ChEMBLError(
                f"Building the ChEMBL extract from the MySQL dump failed: {e}. "
                "CheMBL(source='sqlite') builds the same extract from the "
                "SQLite release instead."
            ) from e
        self.db_path = extract_path

    def _compact_after_download(self) -> None:
        """
        Turn the freshly downloaded release into the extract, and keep only it.

        This is the whole of ``source="sqlite"``: the archive route has to
        write the 27.7 GiB database before anything can be read out of it, so
        the saving is made at the end rather than avoided at the start. The
        extract is verified against the release before the release is deleted
        ([`compact`][provesid.chembl.CheMBL.compact]), so the file that
        survives is the one that was checked.

        ``force=True`` is passed deliberately: a re-download means the caller
        asked for this release again, and an extract left over from a previous
        attempt at the same release would otherwise stop it.

        On success ``db_path`` moves onto the extract, which is the file every
        later query -- and the connection this constructor is about to open --
        will use.

        Notes
        -----
        A failure here is logged, not raised. The release is already in place
        and usable; what is lost is disk, not data, and the alternative would
        be to throw away a 5.8 GB download over a step that can be repeated
        with one call.
        """
        full_path = self.db_path
        extract_path = self.compact_path_for(full_path)
        try:
            self.compact(
                source_path=full_path,
                dest_path=extract_path,
                remove_source=True,
                force=True,
            )
        except Exception as exc:
            self.logger.warning(
                "ChEMBL downloaded and extracted, but building the PROVESID "
                "extract from it failed (%s). The full release is in place at "
                "%s and answers every query; rerun the compaction with "
                "CheMBL.compact(remove_source=True) to reclaim the space.",
                exc, full_path,
            )
            return
        self.db_path = extract_path

    def _extract_database(self, archive_path: str, dest_path: str) -> None:
        """
        Pull the single ``.db`` member out of a ChEMBL archive.

        ChEMBL's archive nests the database at a depth that varies between
        releases --- 37 ships ``chembl_37/chembl_37_sqlite/chembl_37.db`` ---
        so the directories the extraction created are walked back up and
        removed afterwards, stopping at the data directory or at the first one
        that is not empty.

        Args:
            archive_path: The downloaded ``.tar.gz``.
            dest_path: Where the extracted database is moved to. The caller
                passes a staging name, not ``db_path``, so that a failure here
                cannot destroy a release already on disk.

        Raises:
            ChEMBLError: If the archive holds no ``.db`` member.
        """
        with tarfile.open(archive_path, 'r:gz') as tar:
            db_members = [m for m in tar.getmembers() if m.name.endswith('.db')]

            if not db_members:
                raise ChEMBLError("No .db file found in the tar.gz archive")

            db_member = db_members[0]
            self.logger.info(f"Extracting {db_member.name}...")

            with tqdm(total=db_member.size, unit='B', unit_scale=True,
                      desc="Extracting database") as pbar:
                tar.extract(db_member, path=self.path)
                pbar.update(db_member.size)

            extracted_path = os.path.join(self.path, db_member.name)
            os.replace(extracted_path, dest_path)

            extracted_dir = os.path.dirname(extracted_path)
            while os.path.abspath(extracted_dir) != os.path.abspath(self.path):
                try:
                    os.rmdir(extracted_dir)
                except OSError:
                    break
                extracted_dir = os.path.dirname(extracted_dir)

    def _row_to_dict(self, row: sqlite3.Row) -> Dict[str, Any]:
        """Convert sqlite3.Row to dictionary"""
        if row is None:
            return None
        return dict(zip(row.keys(), row))

    def chembl_id_to_molregno(self, chembl_id: str) -> Optional[int]:
        """
        Convert ChEMBL ID to internal molregno identifier.

        Parameters
        ----------
        chembl_id : str
            ChEMBL identifier (e.g., 'CHEMBL25')

        Returns
        -------
        int or None
            Internal molregno ID, or None if not found

        Examples
        --------
        >>> chembl = CheMBL()
        >>> molregno = chembl.chembl_id_to_molregno('CHEMBL25')
        >>> print(molregno)
        1280
        """
        try:
            self.cursor.execute(
                "SELECT entity_id FROM chembl_id_lookup WHERE chembl_id = ? AND entity_type = 'COMPOUND'",
                (chembl_id.upper(),)
            )
            result = self.cursor.fetchone()
            return result[0] if result else None
        except sqlite3.Error as e:
            self.logger.error(f"Database error in chembl_id_to_molregno: {str(e)}")
            return None

    def molregno_to_chembl_id(self, molregno: int) -> Optional[str]:
        """
        Convert internal molregno to ChEMBL ID.

        Parameters
        ----------
        molregno : int
            Internal molecule registry number

        Returns
        -------
        str or None
            ChEMBL identifier, or None if not found

        Examples
        --------
        >>> chembl = CheMBL()
        >>> chembl_id = chembl.molregno_to_chembl_id(1280)
        >>> print(chembl_id)
        CHEMBL25
        """
        try:
            self.cursor.execute(
                "SELECT chembl_id FROM molecule_dictionary WHERE molregno = ?",
                (molregno,)
            )
            result = self.cursor.fetchone()
            return result[0] if result else None
        except sqlite3.Error as e:
            self.logger.error(f"Database error in molregno_to_chembl_id: {str(e)}")
            return None

    def search_by_chembl_id(self, chembl_id: str) -> Optional[Dict[str, Any]]:
        """
        Search for compound by ChEMBL ID.

        Parameters
        ----------
        chembl_id : str
            ChEMBL identifier (e.g., 'CHEMBL25' for aspirin)

        Returns
        -------
        dict or None
            Compound information including structure, or None if not found

        Examples
        --------
        >>> chembl = CheMBL()
        >>> aspirin = chembl.search_by_chembl_id('CHEMBL25')
        >>> print(aspirin['pref_name'])
        ASPIRIN
        """
        molregno = self.chembl_id_to_molregno(chembl_id)
        if molregno is None:
            return None
        return self.get_compound(molregno)

    def search_by_name(
        self, name: str, limit: int = 100, exact: bool = False
    ) -> List[Dict[str, Any]]:
        """
        Search for compounds by name, case-insensitively.

        Searches both preferred names and synonyms.

        Results are ordered by ``molregno``, so the same query returns the same
        compounds in the same order on every call and on every copy of the
        database. That matters more than it sounds: without an ``ORDER BY``,
        ``LIMIT`` silently changes *which* rows come back when a name matches
        more than ``limit`` compounds, and the answer can move after a
        ``VACUUM``, an index change or a SQLite upgrade.

        Parameters
        ----------
        name : str
            Compound name, or partial name when ``exact`` is False.
        limit : int, optional
            Maximum number of results (default: 100). The ``limit`` lowest
            ``molregno`` values among the matches are returned.
        exact : bool, optional
            If True, the name must equal a preferred name or synonym exactly
            (ignoring case). If False (the default), any compound whose
            preferred name or synonym *contains* ``name`` is returned.

            Substring matching is deliberately permissive and will surface
            unrelated compounds for short queries: ``'asprin'`` matches
            PHENYRAMIDOL via its synonym ``'Evasprin'``. Pass ``exact=True``
            when you need the name to actually be the compound's name.

        Returns
        -------
        list of dict
            List of matching compounds with structure information, ordered by
            ``molregno``.

        Notes
        -----
        The two arms of the search — preferred name and synonym — are issued as
        a ``UNION`` of two single-table lookups rather than as an ``OR`` across
        a ``LEFT JOIN``. An ``OR`` whose arms live in different tables cannot
        use an index for either, so the joined form scanned all 2.9 M rows of
        ``molecule_dictionary`` on every call, at a measured 743 ms per exact
        lookup. Each arm of the ``UNION`` is independently indexable, and on a
        database carrying the ``lower(...)`` expression indexes that
        [`compact`][provesid.chembl.CheMBL.compact] builds
        (``ix_md_pref_lower``, ``ix_ms_syn_lower``) an exact lookup costs about
        10 µs — some 77 000× faster. A full ChEMBL release has no such indexes
        and still scans, but it scans two small queries instead of a join.

        ``exact=False`` cannot use those indexes either way: a leading-wildcard
        ``LIKE`` is a scan by construction.

        Examples
        --------
        >>> chembl = CheMBL()
        >>> results = chembl.search_by_name('aspirin', exact=True)
        >>> results[0].get("chembl_id")
        'CHEMBL25'

        A substring search matches more broadly:

        >>> for r in chembl.search_by_name('asprin'):
        ...     print(r["pref_name"])
        PHENYRAMIDOL
        >>> chembl.search_by_name('asprin', exact=True)
        []
        """
        try:
            if exact:
                query = """
                SELECT molregno FROM molecule_dictionary
                WHERE LOWER(pref_name) = LOWER(?)
                UNION
                SELECT molregno FROM molecule_synonyms
                WHERE LOWER(synonyms) = LOWER(?)
                ORDER BY molregno
                LIMIT ?
                """
                search_term = name
            else:
                query = """
                SELECT molregno FROM molecule_dictionary
                WHERE LOWER(pref_name) LIKE LOWER(?)
                UNION
                SELECT molregno FROM molecule_synonyms
                WHERE LOWER(synonyms) LIKE LOWER(?)
                ORDER BY molregno
                LIMIT ?
                """
                search_term = f"%{name}%"

            self.cursor.execute(query, (search_term, search_term, limit))
            results = self.cursor.fetchall()

            compounds = []
            for row in results:
                molregno = row[0]
                compound = self.get_compound(molregno)
                if compound:
                    compounds.append(compound)

            return compounds
        except sqlite3.Error as e:
            self.logger.error(f"Database error in search_by_name: {str(e)}")
            return []

    def search_by_inchi(self, inchi: str) -> Optional[Dict[str, Any]]:
        """
        Search for compound by Standard InChI.

        Parameters
        ----------
        inchi : str
            Standard InChI string

        Returns
        -------
        dict or None
            Compound information, or None if not found

        Examples
        --------
        >>> chembl = CheMBL()
        >>> inchi = 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)'
        >>> compound = chembl.search_by_inchi(inchi)
        >>> print(compound['chembl_id'])
        CHEMBL25
        """
        try:
            self.cursor.execute(
                "SELECT molregno FROM compound_structures WHERE standard_inchi = ?",
                (inchi,)
            )
            result = self.cursor.fetchone()
            if result:
                return self.get_compound(result[0])
            return None
        except sqlite3.Error as e:
            self.logger.error(f"Database error in search_by_inchi: {str(e)}")
            return None

    def search_by_inchikey(self, inchikey: str) -> Optional[Dict[str, Any]]:
        """
        Search for compound by Standard InChI Key.

        Parameters
        ----------
        inchikey : str
            Standard InChI Key (e.g., 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N')

        Returns
        -------
        dict or None
            Compound information, or None if not found

        Examples
        --------
        >>> chembl = CheMBL()
        >>> compound = chembl.search_by_inchikey('BSYNRYMUTXBXSQ-UHFFFAOYSA-N')
        >>> print(compound['pref_name'])
        ASPIRIN
        """
        try:
            self.cursor.execute(
                "SELECT molregno FROM compound_structures WHERE standard_inchi_key = ?",
                (inchikey.upper(),)
            )
            result = self.cursor.fetchone()
            if result:
                return self.get_compound(result[0])
            return None
        except sqlite3.Error as e:
            self.logger.error(f"Database error in search_by_inchikey: {str(e)}")
            return None

    def search_by_smiles(self, smiles: str) -> Optional[Dict[str, Any]]:
        """
        Search for compound by canonical SMILES.

        Note: This performs exact string matching. For similarity searches,
        consider using RDKit or other cheminformatics tools.

        Parameters
        ----------
        smiles : str
            Canonical SMILES string

        Returns
        -------
        dict or None
            Compound information, or None if not found

        Examples
        --------
        >>> chembl = CheMBL()
        >>> compound = chembl.search_by_smiles('CC(=O)Oc1ccccc1C(=O)O')
        >>> print(compound['chembl_id'])
        CHEMBL25
        """
        try:
            self.cursor.execute(
                "SELECT molregno FROM compound_structures WHERE canonical_smiles = ?",
                (smiles,)
            )
            result = self.cursor.fetchone()
            if result:
                return self.get_compound(result[0])
            return None
        except sqlite3.Error as e:
            self.logger.error(f"Database error in search_by_smiles: {str(e)}")
            return None

    def get_compound(self, molregno: int) -> Optional[Dict[str, Any]]:
        """
        Get complete compound information by internal molregno.

        Retrieves data from molecule_dictionary, compound_structures, and molecule_synonyms tables.

        Parameters
        ----------
        molregno : int
            Internal molecule registry number

        Returns
        -------
        dict or None
            Dictionary with compound information including:

            - molregno, chembl_id, pref_name, max_phase, therapeutic_flag,
              molecule_type
            - canonical_smiles, standard_inchi, standard_inchi_key
            - synonyms: list of alternative names

            Returns None if not found.

            No ``molfile`` is returned.  It is a quarter of a full ChEMBL
            database on its own, nothing in PROVESID consumed it, and it is
            absent from the extract [`compact`][provesid.chembl.CheMBL.compact]
            builds.  Build a MOL block from ``canonical_smiles`` with RDKit
            when you need one.

        Examples
        --------
        >>> chembl = CheMBL()
        >>> compound = chembl.get_compound(1280)   # molregno of CHEMBL25
        >>> print(compound['pref_name'])
        ASPIRIN
        >>> print(compound['canonical_smiles'])
        CC(=O)Oc1ccccc1C(=O)O
        >>> print(compound.get("synonyms", [])[:2])
        ['Acetylsalicylic acid', 'Aspirin']
        """
        try:
            query = """
            SELECT
                md.molregno,
                md.chembl_id,
                md.pref_name,
                md.max_phase,
                md.therapeutic_flag,
                md.molecule_type,
                cs.canonical_smiles,
                cs.standard_inchi,
                cs.standard_inchi_key
            FROM molecule_dictionary md
            LEFT JOIN compound_structures cs ON md.molregno = cs.molregno
            WHERE md.molregno = ?
            """
            self.cursor.execute(query, (molregno,))
            result = self.cursor.fetchone()

            if not result:
                return None

            compound = self._row_to_dict(result)

            # Get synonyms
            synonym_query = """
            SELECT synonyms, syn_type
            FROM molecule_synonyms
            WHERE molregno = ?
            ORDER BY syn_type, synonyms
            """
            self.cursor.execute(synonym_query, (molregno,))
            synonym_results = self.cursor.fetchall()

            # Add synonyms as a list to the compound dictionary
            compound['synonyms'] = [row[0] for row in synonym_results] if synonym_results else []

            return compound
        except sqlite3.Error as e:
            self.logger.error(f"Database error in get_compound: {str(e)}")
            return None

    def get_properties(self, molregno: int) -> Optional[Dict[str, Any]]:
        """
        Get physicochemical properties for a compound.

        Parameters
        ----------
        molregno : int
            Internal molecule registry number

        Returns
        -------
        dict or None
            Dictionary with properties including:

            - mw_freebase: Molecular weight
            - alogp: Calculated LogP
            - hba: Hydrogen bond acceptors
            - hbd: Hydrogen bond donors
            - psa: Polar surface area
            - rtb: Rotatable bonds
            - ro3_pass: Rule of 3 compliance
            - num_ro5_violations: Lipinski violations
            - aromatic_rings, heavy_atoms, etc.
            Returns None if not found

        Examples
        --------
        >>> chembl = CheMBL()
        >>> props = chembl.get_properties(1280)  # Aspirin
        >>> print(f"MW: {props['mw_freebase']:.2f}")
        MW: 180.16
        >>> print(f"LogP: {props['alogp']:.2f}")
        LogP: 1.31
        """
        try:
            query = """
            SELECT *
            FROM compound_properties
            WHERE molregno = ?
            """
            self.cursor.execute(query, (molregno,))
            result = self.cursor.fetchone()
            return self._row_to_dict(result) if result else None
        except sqlite3.Error as e:
            self.logger.error(f"Database error in get_properties: {str(e)}")
            return None

    def get_molecule_dictionary(self, molregno: int) -> Optional[Dict[str, Any]]:
        """
        Get complete molecule dictionary information for a compound.

        Retrieves all fields from the molecule_dictionary table including drug classification,
        approval status, administration routes, and other drug-related attributes.

        Parameters
        ----------
        molregno : int
            Internal molecule registry number

        Returns
        -------
        dict or None
            Dictionary with molecule_dictionary fields including:

            - molregno: Internal Primary Key
            - pref_name: Preferred name for the molecule
            - chembl_id: ChEMBL identifier
            - max_phase: Maximum development phase (4=Approved, 3=Phase 3, 2=Phase 2,
                        1=Phase 1, 0.5=Early Phase 1, -1=Clinical Phase unknown, NULL=preclinical)
            - therapeutic_flag: Has therapeutic application (1=yes, 0=no)
            - dosed_ingredient: Drug is dosed in this form (1=yes, 0=no)
            - structure_type: MOL/SEQ/NONE indicating structure availability
            - molecule_type: Small molecule, Protein, Antibody, etc.
            - first_approval: Earliest approval year
            - oral: Administered orally (1=yes, 0=no)
            - parenteral: Administered parenterally (1=yes, 0=no)
            - topical: Administered topically (1=yes, 0=no)
            - black_box_warning: Has black box warning (1=yes, 0=no)
            - natural_product: Is natural product per COCONUT (1=yes, 0=no)
            - first_in_class: First approved drug of its class (1=yes, 0=no, -1=preclinical)
            - chirality: Chirality status (2=achiral, 1=single enantiomer, 0=mixture, -1=unknown)
            - prodrug: Is a pro-drug (1=yes, 0=no, -1=preclinical)
            - inorganic_flag: Is inorganic (1=yes, 0=no, -1=preclinical)
            - usan_year: USAN name application year
            - availability_type: -2=withdrawn, -1=unknown, 0=discontinued, 1=prescription, 2=OTC
            - usan_stem: USAN stem designation
            - polymer_flag: Is small molecule polymer (1=yes, 0=no)
            - usan_substem: USAN substem
            - usan_stem_definition: Definition of USAN stem
            - withdrawn_flag: Withdrawn for toxicity (1=yes, 0=no)
            - chemical_probe: Is chemical probe (1=yes, 0=no)
            - orphan: Orphan designation (1=yes, 0=no, -1=preclinical)
            - veterinary: Has veterinary product (1=yes, 0=no, -1=preclinical)
            Returns None if not found

        Examples
        --------
        >>> chembl = CheMBL()
        >>> mol_dict = chembl.get_molecule_dictionary(1280)  # Aspirin
        >>> print(f"Name: {mol_dict['pref_name']}")
        Name: ASPIRIN
        >>> print(f"Max Phase: {mol_dict['max_phase']}")
        Max Phase: 4
        >>> print(f"First Approval: {mol_dict['first_approval']}")
        First Approval: 1950
        >>> print(f"Oral: {mol_dict['oral']}, Black Box: {mol_dict['black_box_warning']}")
        Oral: 1, Black Box: 0

        Notes
        -----
        This method retrieves all available fields from the molecule_dictionary table.
        For basic compound info with structures, use get_compound() instead.
        """
        try:
            query = """
            SELECT
                molregno,
                pref_name,
                chembl_id,
                max_phase,
                therapeutic_flag,
                dosed_ingredient,
                structure_type,
                molecule_type,
                first_approval,
                oral,
                parenteral,
                topical,
                black_box_warning,
                natural_product,
                first_in_class,
                chirality,
                prodrug,
                inorganic_flag,
                usan_year,
                availability_type,
                usan_stem,
                polymer_flag,
                usan_substem,
                usan_stem_definition,
                withdrawn_flag,
                chemical_probe,
                orphan,
                veterinary
            FROM molecule_dictionary
            WHERE molregno = ?
            """
            self.cursor.execute(query, (molregno,))
            result = self.cursor.fetchone()
            return self._row_to_dict(result) if result else None
        except sqlite3.Error as e:
            self.logger.error(f"Database error in get_molecule_dictionary: {str(e)}")
            return None

    def get_molecule_hierarchy(self, molregno: int) -> Optional[Dict[str, Any]]:
        """
        Get molecule hierarchy information showing parent-salt-metabolite relationships.

        Retrieves data from the molecule_hierarchy table which stores relationships between
        parent compounds, salts, and active metabolites for pro-drugs.

        Parameters
        ----------
        molregno : int
            Internal molecule registry number

        Returns
        -------
        dict or None
            Dictionary with hierarchy information including:

            - molregno: The compound's molregno (has associated data)
            - parent_molregno: Parent compound after removing salts. If same as molregno,
                             no salt component or couldn't be processed
            - active_molregno: For pro-drugs, the active metabolite. If same as
                             parent_molregno, not currently known to be a pro-drug
            - is_parent: Boolean, True if molregno equals parent_molregno
            - is_prodrug: Boolean, True if parent_molregno differs from active_molregno
            - parent_chembl_id: ChEMBL ID of parent compound (if available)
            - active_chembl_id: ChEMBL ID of active metabolite (if available)
            Returns None if not found

        Examples
        --------
        >>> chembl = CheMBL()
        >>> # A salt: valproate sodium (molregno 1991) and its parent
        >>> hierarchy = chembl.get_molecule_hierarchy(1991)
        >>> hierarchy['is_parent'], hierarchy['parent_chembl_id']
        (False, 'CHEMBL109')
        >>> # A pro-drug: enalapril (16847) and the enalaprilat it becomes
        >>> hierarchy = chembl.get_molecule_hierarchy(16847)
        >>> hierarchy['is_prodrug'], hierarchy['active_chembl_id']
        (True, 'CHEMBL577')

        Notes
        -----
        - Parent compounds generated only by removing salts (without their own data)
          appear only in parent_molregno field, not in molregno field
        - When molregno == parent_molregno: compound has no salt or couldn't be processed
        - When parent_molregno == active_molregno: not known to be a pro-drug
        - Compounds with activity data or that are drugs appear in the molregno field
        """
        try:
            query = """
            SELECT
                mh.molregno,
                mh.parent_molregno,
                mh.active_molregno
            FROM molecule_hierarchy mh
            WHERE mh.molregno = ?
            """
            self.cursor.execute(query, (molregno,))
            result = self.cursor.fetchone()

            if not result:
                return None

            hierarchy = self._row_to_dict(result)

            # Add computed flags
            hierarchy['is_parent'] = (hierarchy['molregno'] == hierarchy['parent_molregno'])
            hierarchy['is_prodrug'] = (hierarchy['parent_molregno'] != hierarchy['active_molregno'])

            # Get ChEMBL IDs for parent and active metabolite
            if hierarchy['parent_molregno']:
                hierarchy['parent_chembl_id'] = self.molregno_to_chembl_id(hierarchy['parent_molregno'])

            if hierarchy['active_molregno']:
                hierarchy['active_chembl_id'] = self.molregno_to_chembl_id(hierarchy['active_molregno'])

            return hierarchy
        except sqlite3.Error as e:
            self.logger.error(f"Database error in get_molecule_hierarchy: {str(e)}")
            return None

    def get_pesticide_classifications(self, molregno: int) -> List[Dict[str, Any]]:
        """
        Get pesticide classification information for a compound.

        Retrieves all pesticide classifications (fungicide, herbicide, insecticide)
        associated with a compound according to FRAC, HRAC, and IRAC classification systems.
        A compound may have multiple pesticide classifications.

        Parameters
        ----------
        molregno : int
            Internal molecule registry number

        Returns
        -------
        list of dict
            List of dictionaries containing pesticide classification information:

            - mol_pest_id: Primary key for the mapping
            - pest_class_id: ID of the pesticide classification
            - molregno: Molecule registry number
            - compound_name: Name used in FRAC/HRAC/IRAC classification
            - mec_id: Mechanism of action ID (links to drug_mechanism table)
            - mechanism_comment: Additional mechanism information
            - ref_type: Source of classification (FRAC, HRAC, or IRAC)
            - ref_id: Name of source file
            - ref_url: Full URL to source information
            Returns empty list if no classifications found

        Examples
        --------
        >>> chembl = CheMBL()
        >>> classifications = chembl.get_pesticide_classifications(123456)
        >>> for classification in classifications:
        ...     print(f"Type: {classification['ref_type']}")
        ...     print(f"Name: {classification['compound_name']}")
        ...     print(f"Mechanism: {classification['mechanism_comment']}")
        ...     print(f"Source: {classification['ref_url']}")

        Notes
        -----
        - FRAC: Fungicide Resistance Action Committee
        - HRAC: Herbicide Resistance Action Committee
        - IRAC: Insecticide Resistance Action Committee
        - A single compound may be classified under multiple categories
        - Use compound_name with ref_id and ref_url to locate info in source file
        """
        try:
            query = """
            SELECT
                pcm.mol_pest_id,
                pcm.pest_class_id,
                pcm.molregno,
                pc.compound_name,
                pc.mec_id,
                pc.mechanism_comment,
                pc.ref_type,
                pc.ref_id,
                pc.ref_url
            FROM pesticide_class_mapping pcm
            JOIN pesticide_classification pc ON pcm.pest_class_id = pc.pest_class_id
            WHERE pcm.molregno = ?
            ORDER BY pc.ref_type, pc.compound_name
            """
            self.cursor.execute(query, (molregno,))
            results = self.cursor.fetchall()

            classifications = [self._row_to_dict(row) for row in results]
            return classifications
        except sqlite3.Error as e:
            self.logger.error(f"Database error in get_pesticide_classifications: {str(e)}")
            return []

    def get_pesticide_classification_by_id(self, pest_class_id: int) -> Optional[Dict[str, Any]]:
        """
        Get detailed information for a specific pesticide classification.

        Retrieves pesticide classification details by pest_class_id, including
        mechanism of action information and source references.

        Parameters
        ----------
        pest_class_id : int
            Primary key for the pesticide classification

        Returns
        -------
        dict or None
            Dictionary with pesticide classification details:

            - pest_class_id: Primary key
            - compound_name: Name used in FRAC/HRAC/IRAC classification. Use with
                           ref_id and ref_url to identify row in source file
            - mec_id: Mechanism-of-action identifier (foreign key to drug_mechanism table)
            - mechanism_comment: Additional mechanism information from FRAC/HRAC/IRAC
            - ref_type: Source of classification (FRAC, HRAC, or IRAC)
            - ref_id: Name of file from which classification is derived (see ref_url)
            - ref_url: Full URL for source information
            - associated_molregnos: List of ChEMBL compound molregnos with this classification
            Returns None if not found

        Examples
        --------
        >>> chembl = CheMBL()
        >>> pest_class = chembl.get_pesticide_classification_by_id(42)
        >>> if pest_class:
        ...     print(f"Classification: {pest_class['ref_type']}")
        ...     print(f"Compound: {pest_class['compound_name']}")
        ...     print(f"Mechanism: {pest_class['mechanism_comment']}")
        ...     print(f"Reference: {pest_class['ref_url']}")
        ...     print(f"Applied to {len(pest_class['associated_molregnos'])} compounds")
        Classification: IRAC
        Compound: Nitenpyram
        Mechanism: Nicotinic acetylcholine receptor (nAChR) competitive modulators
        Reference: https://irac-online.org/mode-of-action/
        Applied to 1 compounds

        Notes
        -----
        The compound_name field should be used in conjunction with ref_id and ref_url
        to identify the appropriate row within the source classification file.
        """
        try:
            # Get pesticide classification details
            query = """
            SELECT
                pest_class_id,
                compound_name,
                mec_id,
                mechanism_comment,
                ref_type,
                ref_id,
                ref_url
            FROM pesticide_classification
            WHERE pest_class_id = ?
            """
            self.cursor.execute(query, (pest_class_id,))
            result = self.cursor.fetchone()

            if not result:
                return None

            classification = self._row_to_dict(result)

            # Get all molregnos associated with this classification
            mapping_query = """
            SELECT molregno
            FROM pesticide_class_mapping
            WHERE pest_class_id = ?
            ORDER BY molregno
            """
            self.cursor.execute(mapping_query, (pest_class_id,))
            mapping_results = self.cursor.fetchall()

            classification['associated_molregnos'] = [row[0] for row in mapping_results]

            return classification
        except sqlite3.Error as e:
            self.logger.error(f"Database error in get_pesticide_classification_by_id: {str(e)}")
            return None

    def search_pesticide_by_name(self, name: str, ref_type: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        Search for pesticide classifications by compound name.

        Performs case-insensitive partial matching on compound names in the
        pesticide classification system.

        Parameters
        ----------
        name : str
            Compound name or partial name to search for
        ref_type : str, optional
            Filter by classification source: 'FRAC', 'HRAC', or 'IRAC'
            If None, searches all types

        Returns
        -------
        list of dict
            List of matching pesticide classifications with fields:

            - pest_class_id: Primary key
            - compound_name: Name in classification
            - mec_id: Mechanism ID
            - mechanism_comment: Mechanism details
            - ref_type: FRAC, HRAC, or IRAC
            - ref_id: Source file name
            - ref_url: Source URL
            - molregno_count: Number of compounds with this classification
            Returns empty list if no matches found

        Examples
        --------
        >>> chembl = CheMBL()
        >>> # Search all pesticide types
        >>> results = chembl.search_pesticide_by_name('chloro')
        >>> for result in results[:3]:
        ...     print(f"{result['compound_name']} ({result['ref_type']})")
        chloroneb (FRAC)
        chlorothalonil (FRAC)
        Chlorotoluron (HRAC)
        >>> # Search only fungicides
        >>> fungicides = chembl.search_pesticide_by_name('azole', ref_type='FRAC')
        >>> print(f"Found {len(fungicides)} fungicides")
        Found 25 fungicides

        Notes
        -----
        Valid ref_type values:

        - 'FRAC': Fungicide Resistance Action Committee
        - 'HRAC': Herbicide Resistance Action Committee
        - 'IRAC': Insecticide Resistance Action Committee
        """
        try:
            search_term = f"%{name}%"

            if ref_type:
                query = """
                SELECT
                    pc.pest_class_id,
                    pc.compound_name,
                    pc.mec_id,
                    pc.mechanism_comment,
                    pc.ref_type,
                    pc.ref_id,
                    pc.ref_url,
                    COUNT(pcm.molregno) as molregno_count
                FROM pesticide_classification pc
                LEFT JOIN pesticide_class_mapping pcm ON pc.pest_class_id = pcm.pest_class_id
                WHERE LOWER(pc.compound_name) LIKE LOWER(?)
                  AND pc.ref_type = ?
                GROUP BY pc.pest_class_id, pc.compound_name, pc.mec_id, pc.mechanism_comment,
                         pc.ref_type, pc.ref_id, pc.ref_url
                ORDER BY pc.compound_name
                """
                self.cursor.execute(query, (search_term, ref_type.upper()))
            else:
                query = """
                SELECT
                    pc.pest_class_id,
                    pc.compound_name,
                    pc.mec_id,
                    pc.mechanism_comment,
                    pc.ref_type,
                    pc.ref_id,
                    pc.ref_url,
                    COUNT(pcm.molregno) as molregno_count
                FROM pesticide_classification pc
                LEFT JOIN pesticide_class_mapping pcm ON pc.pest_class_id = pcm.pest_class_id
                WHERE LOWER(pc.compound_name) LIKE LOWER(?)
                GROUP BY pc.pest_class_id, pc.compound_name, pc.mec_id, pc.mechanism_comment,
                         pc.ref_type, pc.ref_id, pc.ref_url
                ORDER BY pc.ref_type, pc.compound_name
                """
                self.cursor.execute(query, (search_term,))

            results = self.cursor.fetchall()
            return [self._row_to_dict(row) for row in results]
        except sqlite3.Error as e:
            self.logger.error(f"Database error in search_pesticide_by_name: {str(e)}")
            return []