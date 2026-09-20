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
never opens, so :meth:`CheMBL.compact` builds an extract holding only the eight
-- about 2.6 GB, with identical results.  See ``CheMBL.compact`` for details.

For detailed table schema information, see src/provesid/data/schema_documentation.txt
"""

import os
import re
import glob
import random
import sqlite3
import tarfile
import logging
import datetime
from typing import Optional, Dict, List, Any, Sequence, Tuple
from urllib.parse import urlsplit
import requests
from tqdm import tqdm
from .datasets import DownloadError, download_file
from .utils import user_dataset_path


class ChEMBLError(Exception):
    """Custom exception for ChEMBL database errors"""
    pass


class CheMBL:
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
    is_compact : bool
        True when the open database is a PROVESID extract built by
        :meth:`compact` rather than a full ChEMBL release.
    provenance : dict or None
        What the extract was built from -- release, source file, build time,
        PROVESID version and per-table row counts.  None for a full release.
    conn : sqlite3.Connection
        SQLite database connection
    cursor : sqlite3.Cursor
        Database cursor for queries

    Examples
    --------
    >>> chembl = CheMBL()
    >>> compound = chembl.search_by_chembl_id('CHEMBL25')  # Aspirin
    >>> print(compound['pref_name'])
    'ASPIRIN'
    >>> props = chembl.get_properties(compound['molregno'])
    >>> print(f"MW: {props['mw_freebase']}")
    MW: 180.16
    
    Notes
    -----
    The database is large and grows with every release: ChEMBL 37 is ~5.8 GB
    compressed and ~30 GB once extracted. Initial setup downloads and extracts it
    from the EMBL-EBI FTP server.

    Most of that is never read.  :meth:`compact` builds a ~2.6 GB extract holding
    only the eight tables this class queries, answering identically::

        CheMBL.compact(remove_source=True)   # 30 GB -> 2.6 GB, then reclaim

    A later ``CheMBL()`` opens the extract in preference to a full release of the
    same number.

    ``latest/`` is a moving directory: it holds only the current release, so the
    archive name changes with every ChEMBL release.  The version is therefore
    resolved from the directory listing rather than pinned, with
    :data:`DEFAULT_DB_URL` as the fallback when the listing cannot be read.
    """

    LATEST_DIR_URL = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/"

    #: Release used when the ``latest/`` listing cannot be read.
    FALLBACK_RELEASE = 37

    DEFAULT_DB_URL = f"{LATEST_DIR_URL}chembl_{FALLBACK_RELEASE}_sqlite.tar.gz"

    #: Matches the SQLite archives advertised in the ``latest/`` listing.
    _ARCHIVE_RE = re.compile(r"chembl_(\d+)_sqlite\.tar\.gz")

    # ── Compaction ────────────────────────────────────────────────────────────

    #: Filename marker distinguishing a PROVESID extract from a full release:
    #: ``chembl_37.db`` is the 30 GB original, ``chembl_37_provesid.db`` the
    #: extract :meth:`compact` builds from it.
    COMPACT_SUFFIX = "_provesid"

    #: Bumped whenever :data:`COMPACT_TABLES` or the kept columns change, so an
    #: extract built by an older PROVESID can be recognised as incomplete
    #: instead of failing later with ``no such table``.
    COMPACT_SCHEMA_VERSION = 1

    #: The only tables this package reads, with the columns kept for each.
    #:
    #: ``columns=None`` keeps every column.  ``where`` restricts the rows: the
    #: ChEMBL id lookup carries an entry for every entity type (assays, targets,
    #: documents), and :meth:`chembl_id_to_molregno` only ever asks about
    #: compounds.  ``molfile`` is deliberately absent from
    #: ``compound_structures``: it is a quarter of the whole database, nothing in
    #: PROVESID consumes it, and a MOL block is reconstructible from the SMILES
    #: with RDKit.
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

    #: Indexes built on the extract, one per lookup the package performs.
    #:
    #: ChEMBL's own indexes are not copied: most of them serve range queries on
    #: ``compound_properties`` (``alogp``, ``psa``, ``rtb``, …) that this package
    #: never issues.  The two ``lower(...)`` expression indexes are new, and are
    #: what makes :meth:`search_by_name` an index lookup rather than a scan of
    #: all 2.9 M rows: 743 ms per exact lookup on a full release, 10 µs here.
    #: They cost a few tens of MB.
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

    #: Table recording where an extract came from; absent in a full release.
    PROVENANCE_TABLE = "provesid_provenance"

    #: How many compounds :meth:`_verify_compact` re-reads from both databases
    #: before an extract is trusted enough to delete its source.
    _VERIFY_SAMPLE = 500

    def __init__(
        self,
        db_name: Optional[str] = None,
        auto_download: bool = True,
        db_url: Optional[str] = None,
        data_dir: Optional[str] = None,
        db_path: Optional[str] = None,
        redownload: bool = False,
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

        Raises
        ------
        FileNotFoundError
            If database not found and auto_download is False
        ChEMBLError
            If database connection or validation fails
        """
        self.logger = logging.getLogger(__name__)
        self.db_url = db_url

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
                self.db_url = self.db_url or self.resolve_latest_db_url()
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

        # Connect to database
        try:
            self.conn = sqlite3.connect(self.db_path)
            self.conn.row_factory = sqlite3.Row  # Enable column access by name
            self.cursor = self.conn.cursor()
            self.logger.info(f"Connected to ChEMBL database at {self.db_path}")
        except sqlite3.Error as e:
            raise ChEMBLError(f"Failed to connect to database: {str(e)}")

        if self.is_compact:
            self._warn_if_extract_is_stale()
    
    def __del__(self):
        """Close database connection when object is destroyed"""
        if hasattr(self, 'conn') and self.conn:
            try:
                self.conn.close()
                self.logger.debug("Closed ChEMBL database connection")
            except Exception as e:
                self.logger.warning(f"Error closing database connection: {str(e)}")
    
    def _warn_if_extract_is_stale(self) -> None:
        """
        Say so when the open extract predates the tables this version needs.

        An extract is a subset of ChEMBL, so a PROVESID that later reads a ninth
        table would meet ``no such table`` with no hint that the database is
        merely old.  Comparing the stored :data:`COMPACT_SCHEMA_VERSION` turns
        that into an instruction.

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
    def resolve_latest_db_url(cls, timeout: float = 30) -> str:
        """
        Resolve the download URL of the current ChEMBL SQLite archive.

        ``latest/`` only ever holds the newest release, so its archive name
        (``chembl_NN_sqlite.tar.gz``) changes with every ChEMBL release.  This
        reads the directory listing and returns the highest ``NN`` it advertises.

        Parameters
        ----------
        timeout : float, optional
            Timeout in seconds for the listing request (default: 30)

        Returns
        -------
        str
            URL of the newest archive, or :data:`DEFAULT_DB_URL` (the pinned
            fallback release) when the listing cannot be read or names no
            archive.  Never raises — a stale pin is better than a hard failure.

        Examples
        --------
        >>> CheMBL.resolve_latest_db_url()
        'https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_37_sqlite.tar.gz'
        """
        logger = logging.getLogger(__name__)
        try:
            response = requests.get(cls.LATEST_DIR_URL, timeout=timeout)
            response.raise_for_status()
            releases = {
                int(m.group(1)) for m in cls._ARCHIVE_RE.finditer(response.text)
            }
        except requests.RequestException as exc:
            logger.warning(
                "Could not read the ChEMBL release listing at %s (%s); "
                "falling back to the pinned release chembl_%d",
                cls.LATEST_DIR_URL, exc, cls.FALLBACK_RELEASE,
            )
            return cls.DEFAULT_DB_URL

        if not releases:
            logger.warning(
                "No chembl_NN_sqlite.tar.gz found in the listing at %s; "
                "falling back to the pinned release chembl_%d",
                cls.LATEST_DIR_URL, cls.FALLBACK_RELEASE,
            )
            return cls.DEFAULT_DB_URL

        release = max(releases)
        logger.info("Resolved current ChEMBL release: chembl_%d", release)
        return f"{cls.LATEST_DIR_URL}chembl_{release}_sqlite.tar.gz"

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
        :meth:`compact` has run, a plain ``CheMBL()`` should open the 2.6 GB
        extract rather than the 30 GB original that may still sit beside it --
        they answer identically, and one of them costs a tenth of the page
        cache.

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

        An extract carries a :data:`PROVENANCE_TABLE`; a full ChEMBL release
        does not.  The check opens the file read-only and never raises for a
        missing or unreadable path — a file that cannot be read is, for this
        purpose, not an extract.

        Parameters
        ----------
        db_path : str
            Path to the SQLite file to inspect.

        Returns
        -------
        bool
            True if the file is an extract built by :meth:`compact`.

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
        Read the provenance record written by :meth:`compact`.

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
        release number that :class:`CheMBL` parses out of the filename.

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
        opens eight of them (:data:`COMPACT_TABLES`) and reads no bioactivity
        data at all, so almost the whole file is dead weight.  The extract
        copies those eight tables, drops the ``molfile`` column -- a quarter of
        the entire database on its own, and consumed nowhere -- keeps only the
        ``COMPOUND`` rows of ``chembl_id_lookup``, and rebuilds just the indexes
        the package's own queries need.

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
            :meth:`compact_path_for` of the source, i.e. ``chembl_37.db``
            produces ``chembl_37_provesid.db`` beside it.
        data_dir : str, optional
            Directory searched for the source and used for the default
            destination.  Defaults to the shared PROVESID dataset directory.
        keep_inchi : bool, optional
            Keep the ``standard_inchi`` column and its index (default: True).
            Dropping it saves a further ~1.0 GB, but :meth:`search_by_inchi`
            then has no column to match and ``Search`` loses the InChI it
            currently reads straight from ChEMBL.  Leave this alone unless disk
            is genuinely short.
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
        Shrink whatever release is already on disk, then reclaim the space::

            from provesid import CheMBL

            path = CheMBL.compact(remove_source=True)
            print(path)   # .../chembl_37_provesid.db

        A later ``CheMBL()`` picks up the extract automatically, because
        :meth:`_find_local_database` prefers it over a full release of the same
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
        Find the full ChEMBL database :meth:`compact` should read.

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
        Return :data:`COMPACT_TABLES` with ``keep_inchi`` applied.

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
            Table specification from :meth:`_compact_table_spec`.

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
            Table specification from :meth:`_compact_table_spec`.
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

            kept_columns = tables["compound_structures"]["columns"]
            for index_name, target in cls.COMPACT_INDEXES:
                if "standard_inchi)" in target and "standard_inchi" not in kept_columns:
                    continue
                conn.execute(f"CREATE INDEX IF NOT EXISTS {index_name} ON {target}")
            logger.info("  built %d indexes", len(cls.COMPACT_INDEXES))

            cls._write_provenance(conn, source_path, tables, row_counts, keep_inchi)

            conn.commit()
            conn.execute("DETACH DATABASE source")
            conn.execute("ANALYZE main")
            conn.commit()
            conn.execute("VACUUM")
        finally:
            conn.close()
        return row_counts

    @classmethod
    def _write_provenance(
        cls,
        conn: sqlite3.Connection,
        source_path: str,
        tables: Dict[str, Dict[str, Any]],
        row_counts: Dict[str, int],
        keep_inchi: bool,
    ) -> None:
        """
        Record what this extract was built from, and by what.

        An extract is a *subset*, so it goes stale differently from a copy: a
        later PROVESID that needs a ninth table has to be able to say so rather
        than fail with ``no such table``.  :data:`COMPACT_SCHEMA_VERSION` is what
        makes that possible, and the rest of the record makes a database on disk
        able to answer where it came from.

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
            Up to :data:`_VERIFY_SAMPLE` molregnos that exist in
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

    def download_database(self, url: Optional[str] = None, force: bool = False):
        """
        Download and extract ChEMBL SQLite database from EMBL-EBI FTP.
        
        Downloads the compressed tar.gz archive (~5.8 GB for release 37), extracts
        the SQLite database (~30 GB), and validates its integrity by querying the
        molecule_dictionary table.

        The download is resumable. An interrupted transfer leaves a ``.part``
        file beside the archive and the next call continues from it, which
        matters more here than anywhere else in the package: this is the
        largest single download PROVESID makes, and it used to start again from
        zero. Consider :meth:`compact` afterwards -- it reduces the extracted
        30 GB to about 2.6 GB.
        
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
        
        Examples
        --------
        >>> chembl = CheMBL(auto_download=False)      # doctest: +SKIP
        >>> chembl.download_database(force=True)      # doctest: +SKIP

        Both lines are marked ``+SKIP`` deliberately: under
        ``pytest --doctest-modules`` this example would otherwise *execute*, and
        fetch ~5.8 GB into whatever ``db_path`` happens to be -- including over a
        compacted extract.
        """
        url = url or self.db_url or self.resolve_latest_db_url()
        self.db_url = url

        # Check if already exists
        if os.path.exists(self.db_path) and not force:
            self.logger.info(f"Database already exists at {self.db_path}")
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
            os.remove(archive_path)

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
        15
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
        >>> chembl_id = chembl.molregno_to_chembl_id(15)
        >>> print(chembl_id)
        'CHEMBL25'
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
        'ASPIRIN'
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
        :meth:`compact` builds (``ix_md_pref_lower``, ``ix_ms_syn_lower``) an
        exact lookup costs about 10 µs — some 77 000× faster. A full ChEMBL
        release has no such indexes and still scans, but it scans two small
        queries instead of a join.

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
        'CHEMBL25'
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
        'ASPIRIN'
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
        'CHEMBL25'
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
            absent from the extract :meth:`compact` builds.  Build a MOL block
            from ``canonical_smiles`` with RDKit when you need one.
        
        Examples
        --------
        >>> chembl = CheMBL()
        >>> compound = chembl.get_compound(15)
        >>> print(compound['pref_name'])
        'ASPIRIN'
        >>> print(compound['canonical_smiles'])
        'CC(=O)Oc1ccccc1C(=O)O'
        >>> print(compound.get("synonyms", [])[:3])
        ['Acetylsalicylic acid', 'Aspirin', '2-Acetoxybenzoic acid']
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
        >>> props = chembl.get_properties(15)  # Aspirin
        >>> print(f"MW: {props['mw_freebase']:.2f}")
        MW: 180.16
        >>> print(f"LogP: {props['alogp']:.2f}")
        LogP: 1.19
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
        >>> mol_dict = chembl.get_molecule_dictionary(15)  # Aspirin
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
        >>> # Example with a salt form
        >>> hierarchy = chembl.get_molecule_hierarchy(1234567)
        >>> if hierarchy['is_parent']:
        ...     print("This is the parent compound (no salt)")
        ... else:
        ...     print(f"Parent: {hierarchy['parent_chembl_id']}")
        >>> 
        >>> # Example with a pro-drug
        >>> hierarchy = chembl.get_molecule_hierarchy(7654321)
        >>> if hierarchy['is_prodrug']:
        ...     print(f"Active metabolite: {hierarchy['active_chembl_id']}")
        
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
        >>> for result in results:
        ...     print(f"{result['compound_name']} ({result['ref_type']})")
        >>> 
        >>> # Search only fungicides
        >>> fungicides = chembl.search_pesticide_by_name('azole', ref_type='FRAC')
        >>> print(f"Found {len(fungicides)} fungicides")
        
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