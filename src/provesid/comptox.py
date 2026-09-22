"""
CompToxID - Interface to CompTox Chemicals Dashboard SQLite database for fast identifier lookup and conversion.

This class provides access to a local SQLite database containing CompTox chemicals
with their identifiers (DTXSID, DTXCID, CASRN, InChIKey, SMILES, PREFERRED_NAME, etc.)
and chemical properties (molecular formula, average mass, monoisotopic mass, etc.).

The database is read from comptox_chemicals.db file.

Attributes:
    db_path (str): Path to the SQLite database file
    conn (sqlite3.Connection): Database connection

Records are dicts keyed by the database's upper-case column names
(``DTXSID``, ``PREFERRED_NAME``, ``CASRN``, ``INCHIKEY``, ``SMILES`` ...), plus
``identifiers``, the ``IDENTIFIER`` column split into a list.

Examples:
    >>> from provesid import CompToxID
    >>> db = CompToxID()
    >>> result = db.get_by_casrn("50-78-2")  # Aspirin
    >>> result['PREFERRED_NAME'], result['DTXSID']
    ('Aspirin', 'DTXSID5020108')
    >>> db.batch_casrn_to_dtxsid(["50-78-2", "50-00-0"])
    {'50-78-2': 'DTXSID5020108', '50-00-0': 'DTXSID7020637'}
"""

import os
import re
import sqlite3
import logging
import threading

from typing import Dict, List, Optional, Any, Union

from .datasets import download_file
from .sqlite_client import SQLiteClient
from .utils import user_dataset_path

#: The table :meth:`CompToxID.build_name_index` adds to the database: one row
#: per distinct name of each chemical, keyed by :func:`name_key`.
NAME_INDEX_TABLE = "chemical_names"

#: The shape of a CAS Registry Number, for inputs that must be one.
_CAS_NUMBER = re.compile(r"\d{2,7}-\d{2}-\d")

#: Where a name came from, in the order an exact lookup ranks its matches: a
#: chemical *called* the query outranks one that merely lists it as a synonym.
NAME_KINDS = {"preferred": 0, "iupac": 1, "identifier": 2}


def name_key(name: str) -> str:
    """The form a name is indexed and looked up under: stripped and lower-cased.

    Python's :meth:`str.lower` rather than SQLite's ``lower()``, which folds
    ASCII only; the same function builds the index and reads it, so the two
    always agree.

    Args:
        name: A chemical name or other identifier.

    Returns:
        The lookup key.

    Examples:
        >>> name_key("  Acetylsalicylic Acid ")
        'acetylsalicylic acid'
    """
    return name.strip().lower()


def _build_name_index(connection: sqlite3.Connection) -> int:
    """Create and fill :data:`NAME_INDEX_TABLE` in one transaction.

    Every chemical contributes its ``PREFERRED_NAME``, its ``IUPAC_NAME`` and
    each ``|``-separated token of ``IDENTIFIER`` --- synonyms, but also former
    CAS numbers, InChIKeys and registry codes --- once per distinct key, tagged
    with the first :data:`NAME_KINDS` it appeared under.  The table is
    ``WITHOUT ROWID`` with the key leading its primary key, so the text is
    stored once and the lookup is a single B-tree search: on the 2025 release,
    5.1 M rows, ~290 MiB and ~20 s to build.  A separate index over an
    ordinary table measured 540 MiB.

    ``BEGIN IMMEDIATE`` takes the write lock before the existence check, so of
    two processes building at once the second waits and then finds the table
    there.

    Args:
        connection: A read-write connection to a CompTox database.

    Returns:
        The number of rows in the table afterwards.

    Raises:
        sqlite3.OperationalError: If the database is read-only or stays locked
            beyond the connection's timeout.
    """
    connection.execute("BEGIN IMMEDIATE")
    try:
        exists = connection.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
            (NAME_INDEX_TABLE,),
        ).fetchone()
        if not exists:
            connection.execute(
                f"""
                CREATE TABLE {NAME_INDEX_TABLE} (
                    name_key TEXT NOT NULL,
                    kind INTEGER NOT NULL,
                    chemical_rowid INTEGER NOT NULL,
                    PRIMARY KEY (name_key, kind, chemical_rowid)
                ) WITHOUT ROWID
                """
            )

            def names():
                chemicals = connection.cursor().execute(
                    "SELECT rowid, PREFERRED_NAME, IUPAC_NAME, IDENTIFIER FROM chemicals"
                )
                for rowid, preferred, iupac, identifier in chemicals:
                    labelled = [
                        (NAME_KINDS["preferred"], preferred),
                        (NAME_KINDS["iupac"], iupac),
                        *((NAME_KINDS["identifier"], token)
                          for token in (identifier or "").split("|")),
                    ]
                    seen = set()
                    for kind, name in labelled:
                        key = name_key(name) if name else ""
                        if key and key not in seen:
                            seen.add(key)
                            yield key, kind, rowid

            connection.executemany(
                f"INSERT INTO {NAME_INDEX_TABLE} VALUES (?, ?, ?)", names()
            )
        connection.commit()
    except BaseException:
        connection.rollback()
        raise
    return connection.execute(f"SELECT COUNT(*) FROM {NAME_INDEX_TABLE}").fetchone()[0]


class CompToxID(SQLiteClient):
    """
    Interface to CompTox Chemicals Dashboard SQLite database.

    The database file is automatically downloaded on first use when missing.

    Inherits its connection handling from
    :class:`~provesid.sqlite_client.SQLiteClient`: use it as a context
    manager, or call :meth:`~provesid.sqlite_client.SQLiteClient.close` when
    finished, and query it from as many threads as you like --- each gets its
    own connection.

    Examples:
        >>> with CompToxID() as db:                     # doctest: +SKIP
        ...     dtxsid = db.casrn_to_dtxsid("50-78-2")
    """

    # Default database filename
    DEFAULT_DB_NAME = "comptox_chemicals.db"
    DEFAULT_DB_URL = (
        "https://zenodo.org/records/18833587/files/comptox_chemicals.db"
    )
    DEFAULT_DB_SIZE_MB = 856

    def __init__(
        self,
        db_path: Optional[str] = None,
        auto_download: bool = True,
        db_url: Optional[str] = None,
        data_dir: Optional[str] = None,
        redownload: bool = False,
    ):
        """
        Initialize CompToxID database connection.

        Args:
            db_path (str, optional): Path to SQLite database. If None, uses default
                                    location in the persistent user dataset directory.
            auto_download (bool, optional): If True, automatically download the
                database when missing (default: True).
            db_url (str, optional): Custom URL for database download. If None,
                uses the default Zenodo URL.
            data_dir (str, optional): Directory to store the database when
                ``db_path`` is not provided. If None, uses platformdirs-based
                user data directory.
            redownload (bool, optional): If True, force a fresh download when
                ``auto_download`` is enabled.

        Raises:
            FileNotFoundError: If database file doesn't exist and auto_download is False.
        """
        self.logger = logging.getLogger(__name__)

        # Whether exact name lookups can use the name index: None until the
        # first one checks, then True, or False for good when it cannot be
        # built (a read-only file).  The lock keeps two threads from building
        # it twice.
        self._name_index_ready: Optional[bool] = None
        self._name_index_lock = threading.Lock()

        if db_path is None:
            base_dir = data_dir or user_dataset_path()
            db_path = os.path.join(base_dir, self.DEFAULT_DB_NAME)

        self.db_path = os.path.abspath(os.path.expanduser(db_path))
        self.db_url = db_url or self.DEFAULT_DB_URL

        needs_download = redownload or not os.path.exists(self.db_path)

        # Check if database exists
        if needs_download:
            if auto_download:
                if redownload and os.path.exists(self.db_path):
                    self.logger.warning(
                        "Forced CompTox redownload requested for: %s", self.db_path
                    )
                else:
                    self.logger.warning(f"CompTox database not found at: {self.db_path}")
                self.logger.warning(
                    "The CompTox database is large (~856 MB). "
                    "Initial setup may take several minutes depending on your connection."
                )
                self.logger.warning(
                    f"Downloading CompTox database from: {self.db_url}"
                )
                self.download_database(url=self.db_url, force=redownload)
            else:
                raise FileNotFoundError(
                    f"CompTox database not found at: {self.db_path}\n"
                    "Database size: ~856 MB\n"
                    f"Run CompToxID.download_database() or set auto_download=True\n"
                    f"Download URL: {self.db_url}"
                )

        # Connect to the database.  One connection per thread, released by
        # close() or by leaving a ``with`` block --- see
        # :class:`~provesid.sqlite_client.SQLiteClient`.
        self._open_database(self.db_path)

        # Verify the database has the expected table
        self._verify_database()

    def download_database(self, url: Optional[str] = None, force: bool = False) -> str:
        """
        Download the CompTox SQLite database from Zenodo.

        The file is approximately 856 MB and is not shipped with the GitHub
        repository due to size limitations.

        The transfer is resumable: an interrupted download leaves a ``.part``
        file beside the destination and the next call continues from it rather
        than fetching the 856 MB again. The file is checked before it is moved
        into place, so a failed download never replaces a working database.

        Args:
            url (str, optional): Download URL. If None, uses `self.db_url`.
            force (bool, optional): If True, overwrite existing database file.

        Returns:
            str: Path to the downloaded database file.

        Raises:
            FileExistsError: If the database already exists and `force` is False.
            provesid.datasets.DownloadError: If the download could not be
                completed.
            RuntimeError: If the file that arrived is not the CompTox database.

        Examples:
            >>> db = CompToxID()
            >>> db.download_database(force=True)            # doctest: +SKIP
            '/home/me/.local/share/provesid/comptox_chemicals.db'
        """
        download_url = url or self.db_url

        if os.path.exists(self.db_path) and not force:
            raise FileExistsError(
                f"Database already exists at: {self.db_path}. Use force=True to overwrite."
            )

        def must_contain_the_chemicals_table(path):
            """Reject a download that is not the CompTox database.

            Run on the ``.part`` file, before it is moved into place. The
            previous implementation renamed first and checked afterwards, so a
            failed check left the broken file where the good one had been.
            """
            connection = sqlite3.connect(path)
            try:
                found = connection.execute(
                    "SELECT name FROM sqlite_master "
                    "WHERE type='table' AND name='chemicals'"
                ).fetchone()
            finally:
                connection.close()
            if not found:
                raise RuntimeError(
                    "Downloaded database does not contain 'chemicals' table"
                )

        self.logger.warning(
            "CompTox database download starting (~856 MB). "
            "Please ensure you have enough disk space and stable internet."
        )
        download_file(
            download_url,
            self.db_path,
            verify=must_contain_the_chemicals_table,
            description="CompTox database",
            log=self.logger,
        )

        # Built now, while the user is already waiting for a download, rather
        # than on the first name lookup.  A failure here costs nothing that
        # the lazy build in search_by_name will not retry.
        self.logger.warning("Building the CompTox name index (~20 s, ~290 MiB).")
        connection = sqlite3.connect(self.db_path)
        try:
            _build_name_index(connection)
        except sqlite3.Error as exc:
            self.logger.warning("CompTox name index not built: %s", exc)
        finally:
            connection.close()
        self._name_index_ready = None
        return self.db_path

    @property
    def has_name_index(self) -> bool:
        """Whether the database holds the name index (see :meth:`build_name_index`).

        Returns:
            True when :data:`NAME_INDEX_TABLE` exists.

        Examples:
            >>> isinstance(CompToxID().has_name_index, bool)
            True
        """
        return self.conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
            (NAME_INDEX_TABLE,),
        ).fetchone() is not None

    def build_name_index(self) -> int:
        """Index every name of every chemical, so that exact lookups find synonyms.

        The downloaded database indexes ``PREFERRED_NAME`` only; the synonyms,
        former CAS numbers and registry codes sit together in the
        ``|``-separated ``IDENTIFIER`` column, which only a full scan can read.
        This adds a table, :data:`NAME_INDEX_TABLE`, with one row per distinct
        name of each chemical (compared by :func:`name_key`), and
        :meth:`search_by_name` with ``exact=True`` uses it from then on.

        It is built automatically after :meth:`download_database`, and by the
        first exact :meth:`search_by_name` on a database downloaded before the
        index existed.  Call it yourself to pay the ~20 s at a time of your
        choosing.  Building it again is a no-op.

        Returns:
            int: The number of rows in the index (5.1 M on the 2025 release).

        Raises:
            sqlite3.OperationalError: If the database file is read-only.

        Examples:
            >>> with CompToxID() as db:                     # doctest: +SKIP
            ...     db.build_name_index()
            5128983
        """
        with self._name_index_lock:
            rows = _build_name_index(self.conn)
            self._name_index_ready = True
        return rows

    def _ensure_name_index(self) -> bool:
        """Make sure the name index exists, building it once if it does not.

        Returns:
            True when exact lookups can use the index; False when it is missing
            and cannot be built, in which case they fall back to
            ``PREFERRED_NAME`` alone.  The failure is logged once.
        """
        if self._name_index_ready is not None:
            return self._name_index_ready
        with self._name_index_lock:
            if self._name_index_ready is None:
                if self.has_name_index:
                    self._name_index_ready = True
                else:
                    self.logger.warning(
                        "Building the CompTox name index, once, so exact name "
                        "lookups find synonyms (~20 s, ~290 MiB added to %s).",
                        self.db_path,
                    )
                    try:
                        _build_name_index(self.conn)
                        self._name_index_ready = True
                    except sqlite3.OperationalError as exc:
                        self.logger.warning(
                            "CompTox name index could not be built (%s); exact "
                            "name lookups will match preferred names only.",
                            exc,
                        )
                        self._name_index_ready = False
        return self._name_index_ready

    def _verify_database(self):
        """Verify the database has the expected table structure."""
        cursor = self.conn.cursor()
        cursor.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='chemicals'"
        )
        if not cursor.fetchone():
            raise RuntimeError("Database does not contain 'chemicals' table")

        # Check for required columns
        cursor.execute("PRAGMA table_info(chemicals)")
        columns = [row[1] for row in cursor.fetchall()]
        required_columns = {
            "DTXSID",
            "PREFERRED_NAME",
            "CASRN",
            "DTXCID",
            "INCHIKEY",
            "SMILES",
            "MOLECULAR_FORMULA",
        }
        missing = required_columns - set(columns)
        if missing:
            raise RuntimeError(f"Database table missing required columns: {missing}")

    # Basic lookup methods

    @staticmethod
    def _parse_identifiers(identifier_string: Optional[str]) -> List[str]:
        """
        Parse pipe-separated identifier string into list of identifiers.

        The IDENTIFIER column contains pipe-separated identifiers and synonyms.
        Format: "identifier1 | identifier2 | identifier3"

        Args:
            identifier_string: Pipe-separated identifier string

        Returns:
            List of identifiers (stripped of whitespace)
        """
        if not identifier_string:
            return []

        # Split by pipe character, strip whitespace
        identifiers = [id.strip() for id in identifier_string.split("|")]
        # Filter out empty strings
        return [id for id in identifiers if id]

    def get_by_dtxsid(self, dtxsid: str) -> Optional[Dict[str, Any]]:
        """
        Get chemical information by DTXSID.

        Args:
            dtxsid (str): DSSTox Substance ID (e.g., "DTXSID7020001")

        Returns:
            dict: Every column of the ``chemicals`` row --- ``DTXSID``,
            ``DTXCID``, ``PREFERRED_NAME``, ``CASRN``, ``INCHIKEY``,
            ``IUPAC_NAME``, ``SMILES``, ``MS_READY_SMILES``,
            ``QSAR_READY_SMILES``, ``MOLECULAR_FORMULA``, ``AVERAGE_MASS``,
            ``MONOISOTOPIC_MASS`` and ``IDENTIFIER`` --- plus ``identifiers``,
            the last split into a list. None if not found.

        Examples:
            >>> record = CompToxID().get_by_dtxsid("DTXSID5020108")
            >>> record["PREFERRED_NAME"], record["CASRN"], record["identifiers"][:2]
            ('Aspirin', '50-78-2', ['50-78-2', '11126-35-5'])
        """
        cursor = self.conn.cursor()
        cursor.execute(
            """
            SELECT * FROM chemicals WHERE DTXSID = ?
        """,
            (dtxsid,),
        )

        row = cursor.fetchone()
        if not row:
            return None

        result = dict(row)

        # Parse identifiers
        result["identifiers"] = self._parse_identifiers(result.get("IDENTIFIER"))

        return result

    def get_by_casrn(self, casrn: str) -> Optional[Dict[str, Any]]:
        """
        Get chemical information by CAS Registry Number.

        Args:
            casrn (str): CAS Registry Number (e.g., "50-78-2")

        Returns:
            dict: The :meth:`get_by_dtxsid` record, or None if not found. Only
            the chemical's own ``CASRN`` matches; for a retired or alternate
            number see :meth:`get_by_alternate_casrn`.

        Examples:
            >>> CompToxID().get_by_casrn("50-78-2")["DTXSID"]
            'DTXSID5020108'
        """
        cursor = self.conn.cursor()
        cursor.execute(
            """
            SELECT * FROM chemicals WHERE CASRN = ?
        """,
            (casrn,),
        )

        row = cursor.fetchone()
        if not row:
            return None

        result = dict(row)
        result["identifiers"] = self._parse_identifiers(result.get("IDENTIFIER"))
        return result

    def get_by_alternate_casrn(self, casrn: str) -> Optional[Dict[str, Any]]:
        """
        Get the chemical that lists a CAS number other than its own ``CASRN``.

        CAS deletes and merges registry numbers, and old datasets still carry
        the numbers it retired. CompTox keeps them, together with alternate
        numbers, among a chemical's ``IDENTIFIER`` tokens, where
        :meth:`get_by_casrn` does not look. For example, atrazine is
        ``1912-24-9`` but also lists ``39400-72-1``. This method reads the name
        index (see :meth:`build_name_index`) for such a number.

        It answers only when **exactly one** chemical lists the number.
        Of the 83,933 numbers CompTox holds only in ``IDENTIFIER``, four are
        listed by two unrelated chemicals, and picking one of the two would
        be a guess. Call :meth:`get_by_casrn` first: a number that is some
        chemical's own ``CASRN`` belongs to that chemical, whatever else
        lists it.

        Args:
            casrn (str): CAS Registry Number (e.g., "39400-72-1")

        Returns:
            dict: Chemical information, or None if ``casrn`` is not shaped
            like a CAS number, no chemical lists it, more than one does, or the
            name index is unavailable (a read-only database that predates it).

        Examples:
            >>> with CompToxID() as db:                     # doctest: +SKIP
            ...     db.get_by_alternate_casrn("39400-72-1")["PREFERRED_NAME"]
            'Atrazine'
        """
        # IDENTIFIER also holds synonyms; without this, a name one chemical
        # lists would come back as its "alternate CAS number".
        if not _CAS_NUMBER.fullmatch(casrn.strip()):
            return None
        if not self._ensure_name_index():
            return None
        cursor = self.conn.cursor()
        cursor.execute(
            f"""
            SELECT c.* FROM {NAME_INDEX_TABLE} n
            JOIN chemicals c ON c.rowid = n.chemical_rowid
            WHERE n.name_key = ? AND n.kind = ?
            LIMIT 2
        """,
            (name_key(casrn), NAME_KINDS["identifier"]),
        )
        rows = cursor.fetchall()
        if len(rows) != 1:
            if rows:
                self.logger.debug(
                    "CAS %s is listed by more than one CompTox chemical; not guessing.",
                    casrn,
                )
            return None

        result = dict(rows[0])
        result["identifiers"] = self._parse_identifiers(result.get("IDENTIFIER"))
        return result

    def get_by_inchikey(self, inchikey: str) -> Optional[Dict[str, Any]]:
        """
        Get chemical information by InChIKey.

        Args:
            inchikey (str): Standard InChIKey (27 characters)

        Returns:
            dict: The :meth:`get_by_dtxsid` record, or None if not found

        Examples:
            >>> CompToxID().get_by_inchikey("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")["DTXSID"]
            'DTXSID5020108'
        """
        cursor = self.conn.cursor()
        cursor.execute(
            """
            SELECT * FROM chemicals WHERE INCHIKEY = ?
        """,
            (inchikey,),
        )

        row = cursor.fetchone()
        if not row:
            return None

        result = dict(row)
        result["identifiers"] = self._parse_identifiers(result.get("IDENTIFIER"))
        return result

    def get_by_smiles(self, smiles: str) -> Optional[Dict[str, Any]]:
        """
        Get chemical information by SMILES string.

        Args:
            smiles (str): SMILES string, matched as a string against CompTox's
                own: another valid SMILES for the same structure finds nothing

        Returns:
            dict: The :meth:`get_by_dtxsid` record, or None if not found

        Examples:
            >>> db = CompToxID()
            >>> db.get_by_smiles("CC(=O)OC1=C(C=CC=C1)C(O)=O")["CASRN"]
            '50-78-2'
            >>> db.get_by_smiles("CC(=O)OC1=CC=CC=C1C(O)=O") is None
            True
        """
        cursor = self.conn.cursor()
        cursor.execute(
            """
            SELECT * FROM chemicals WHERE SMILES = ?
        """,
            (smiles,),
        )

        row = cursor.fetchone()
        if not row:
            return None

        result = dict(row)
        result["identifiers"] = self._parse_identifiers(result.get("IDENTIFIER"))
        return result

    def get_by_name(self, name: str) -> Optional[Dict[str, Any]]:
        """
        Get chemical information by preferred name (exact match).

        Args:
            name (str): Preferred name, case included. :meth:`search_by_name`
                with ``exact=True`` matches any name, ignoring case.

        Returns:
            dict: The :meth:`get_by_dtxsid` record, or None if not found

        Examples:
            >>> db = CompToxID()
            >>> db.get_by_name("Aspirin")["CASRN"], db.get_by_name("aspirin")
            ('50-78-2', None)
        """
        cursor = self.conn.cursor()
        cursor.execute(
            """
            SELECT * FROM chemicals WHERE PREFERRED_NAME = ?
        """,
            (name,),
        )

        row = cursor.fetchone()
        if not row:
            return None

        result = dict(row)
        result["identifiers"] = self._parse_identifiers(result.get("IDENTIFIER"))
        return result

    def get_by_dtxcid(self, dtxcid: str) -> Optional[Dict[str, Any]]:
        """
        Get chemical information by DTXCID.

        Args:
            dtxcid (str): DSSTox Compound ID (e.g., "DTXCID101")

        Returns:
            dict: The :meth:`get_by_dtxsid` record, or None if not found

        Examples:
            >>> CompToxID().get_by_dtxcid("DTXCID50108")["PREFERRED_NAME"]
            'Aspirin'
        """
        cursor = self.conn.cursor()
        cursor.execute(
            """
            SELECT * FROM chemicals WHERE DTXCID = ?
        """,
            (dtxcid,),
        )

        row = cursor.fetchone()
        if not row:
            return None

        result = dict(row)
        result["identifiers"] = self._parse_identifiers(result.get("IDENTIFIER"))
        return result

    def search_by_name(
        self, name: str, exact: bool = False, limit: int = 10
    ) -> List[Dict[str, Any]]:
        """
        Search chemicals by name or synonym.

        With ``exact=True`` the query is compared, case-insensitively, with
        every name the database holds for a chemical: its preferred name, its
        IUPAC name and each synonym or identifier in ``IDENTIFIER``.  Chemicals
        *called* the query come first, then those whose IUPAC name it is, then
        those that list it as a synonym; ties keep the database's order.  This
        reads the name index, which the first exact call builds if the database
        predates it (see :meth:`build_name_index`); on a read-only file without
        one, only preferred names are matched, case-sensitively, as before the
        index existed.

        With ``exact=False`` the query is matched as a substring, first of the
        preferred name and then of ``IDENTIFIER``.  That is a full scan, about
        3 s a call.

        Args:
            name (str): Chemical name or synonym to search for
            exact (bool): If True, exact (case-insensitive) match on any name.
                If False, partial match (case-insensitive)
            limit (int): Maximum number of results to return

        Returns:
            list: List of matching chemicals

        Examples:
            >>> with CompToxID() as db:                     # doctest: +SKIP
            ...     [r["PREFERRED_NAME"] for r in db.search_by_name("Acetaldoxime", exact=True)]
            ['Acetaldehyde oxime']
        """
        cursor = self.conn.cursor()
        results = []

        if exact and self._ensure_name_index():
            cursor.execute(
                f"""
                SELECT c.* FROM {NAME_INDEX_TABLE} n
                JOIN chemicals c ON c.rowid = n.chemical_rowid
                WHERE n.name_key = ?
                ORDER BY n.kind, n.chemical_rowid
                LIMIT ?
            """,
                (name_key(name), limit),
            )
        elif exact:
            cursor.execute(
                """
                SELECT * FROM chemicals WHERE PREFERRED_NAME = ? LIMIT ?
            """,
                (name, limit),
            )
        else:
            # Partial match with LIKE (case-insensitive)
            search_term = f"%{name}%"
            cursor.execute(
                """
                SELECT * FROM chemicals WHERE PREFERRED_NAME LIKE ? LIMIT ?
            """,
                (search_term, limit),
            )

        rows = cursor.fetchall()
        for row in rows:
            result = dict(row)
            result["identifiers"] = self._parse_identifiers(result.get("IDENTIFIER"))
            results.append(result)

        # If not enough results and not exact, also search in identifiers
        if len(results) < limit and not exact:
            search_term = f"%{name}%"
            cursor.execute(
                """
                SELECT * FROM chemicals WHERE IDENTIFIER LIKE ? LIMIT ?
            """,
                (search_term, limit - len(results)),
            )

            # Need to deduplicate by DTXSID
            seen_dtxsids = {r["DTXSID"] for r in results}
            for row in cursor.fetchall():
                if row["DTXSID"] not in seen_dtxsids:
                    result = dict(row)
                    result["identifiers"] = self._parse_identifiers(
                        result.get("IDENTIFIER")
                    )
                    results.append(result)
                    seen_dtxsids.add(row["DTXSID"])
                if len(results) >= limit:
                    break

        return results

    def search_by_formula(self, formula: str, limit: int = 100) -> List[Dict[str, Any]]:
        """
        Search chemicals by molecular formula.

        Args:
            formula (str): Molecular formula (e.g., "C9H8O4")
            limit (int): Maximum number of results to return

        Returns:
            list: :meth:`get_by_dtxsid` records, in database order

        Examples:
            >>> [r["PREFERRED_NAME"] for r in CompToxID().search_by_formula("C9H8O4", limit=2)]
            ['Aspirin', '3,4-Dihydroxycinnamic acid']
        """
        cursor = self.conn.cursor()
        cursor.execute(
            """
            SELECT * FROM chemicals WHERE MOLECULAR_FORMULA = ? LIMIT ?
        """,
            (formula, limit),
        )

        results = []
        for row in cursor.fetchall():
            result = dict(row)
            result["identifiers"] = self._parse_identifiers(result.get("IDENTIFIER"))
            results.append(result)

        return results

    # Conversion methods

    def casrn_to_dtxsid(self, casrn: str) -> Optional[str]:
        """
        Convert CAS Registry Number to DTXSID.

        Args:
            casrn: CAS Registry Number.

        Returns:
            The DTXSID, or None if not found.

        Examples:
            >>> db = CompToxID()
            >>> db.casrn_to_dtxsid("50-78-2")
            'DTXSID5020108'
        """
        result = self.get_by_casrn(casrn)
        return result["DTXSID"] if result else None

    def casrn_to_inchikey(self, casrn: str) -> Optional[str]:
        """
        Convert CAS Registry Number to InChIKey.

        Args:
            casrn: CAS Registry Number.

        Returns:
            The InChIKey, or None if not found.

        Examples:
            >>> db = CompToxID()
            >>> db.casrn_to_inchikey("50-78-2")
            'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'
        """
        result = self.get_by_casrn(casrn)
        return result["INCHIKEY"] if result else None

    def casrn_to_smiles(self, casrn: str) -> Optional[str]:
        """
        Convert CAS Registry Number to SMILES.

        Args:
            casrn: CAS Registry Number.

        Returns:
            CompTox's SMILES, or None if not found.

        Examples:
            >>> db = CompToxID()
            >>> db.casrn_to_smiles("50-78-2")
            'CC(=O)OC1=C(C=CC=C1)C(O)=O'
        """
        result = self.get_by_casrn(casrn)
        return result["SMILES"] if result else None

    def inchikey_to_casrn(self, inchikey: str) -> Optional[str]:
        """
        Convert InChIKey to CAS Registry Number.

        Args:
            inchikey: Standard InChIKey.

        Returns:
            The CAS number, or None if not found.

        Examples:
            >>> db = CompToxID()
            >>> db.inchikey_to_casrn("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
            '50-78-2'
        """
        result = self.get_by_inchikey(inchikey)
        return result["CASRN"] if result else None

    def inchikey_to_dtxsid(self, inchikey: str) -> Optional[str]:
        """
        Convert InChIKey to DTXSID.

        Args:
            inchikey: Standard InChIKey.

        Returns:
            The DTXSID, or None if not found.

        Examples:
            >>> db = CompToxID()
            >>> db.inchikey_to_dtxsid("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
            'DTXSID5020108'
        """
        result = self.get_by_inchikey(inchikey)
        return result["DTXSID"] if result else None

    def dtxsid_to_casrn(self, dtxsid: str) -> Optional[str]:
        """
        Convert DTXSID to CAS Registry Number.

        Args:
            dtxsid: DSSTox Substance ID.

        Returns:
            The CAS number, or None if not found.

        Examples:
            >>> db = CompToxID()
            >>> db.dtxsid_to_casrn("DTXSID5020108")
            '50-78-2'
        """
        result = self.get_by_dtxsid(dtxsid)
        return result["CASRN"] if result else None

    def dtxsid_to_inchikey(self, dtxsid: str) -> Optional[str]:
        """
        Convert DTXSID to InChIKey.

        Args:
            dtxsid: DSSTox Substance ID.

        Returns:
            The InChIKey, or None if not found.

        Examples:
            >>> db = CompToxID()
            >>> db.dtxsid_to_inchikey("DTXSID5020108")
            'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'
        """
        result = self.get_by_dtxsid(dtxsid)
        return result["INCHIKEY"] if result else None

    def dtxsid_to_smiles(self, dtxsid: str) -> Optional[str]:
        """
        Convert DTXSID to SMILES.

        Args:
            dtxsid: DSSTox Substance ID.

        Returns:
            CompTox's SMILES, or None if not found.

        Examples:
            >>> db = CompToxID()
            >>> db.dtxsid_to_smiles("DTXSID5020108")
            'CC(=O)OC1=C(C=CC=C1)C(O)=O'
        """
        result = self.get_by_dtxsid(dtxsid)
        return result["SMILES"] if result else None

    def smiles_to_casrn(self, smiles: str) -> Optional[str]:
        """
        Convert SMILES to CAS Registry Number.

        Args:
            smiles: SMILES, matched as a string; see :meth:`get_by_smiles`.

        Returns:
            The CAS number, or None if not found.

        Examples:
            >>> db = CompToxID()
            >>> db.smiles_to_casrn("CC(=O)OC1=C(C=CC=C1)C(O)=O")
            '50-78-2'
        """
        result = self.get_by_smiles(smiles)
        return result["CASRN"] if result else None

    def smiles_to_dtxsid(self, smiles: str) -> Optional[str]:
        """
        Convert SMILES to DTXSID.

        Args:
            smiles: SMILES, matched as a string; see :meth:`get_by_smiles`.

        Returns:
            The DTXSID, or None if not found.

        Examples:
            >>> db = CompToxID()
            >>> db.smiles_to_dtxsid("CCO")
            'DTXSID9020584'
        """
        result = self.get_by_smiles(smiles)
        return result["DTXSID"] if result else None

    # Batch conversion methods

    def batch_casrn_to_dtxsid(self, casrn_list: List[str]) -> Dict[str, Optional[str]]:
        """
        Convert multiple CAS numbers to DTXSIDs.

        Args:
            casrn_list (list): List of CAS numbers

        Returns:
            dict: Mapping of CAS -> DTXSID (None if not found)

        Examples:
            >>> CompToxID().batch_casrn_to_dtxsid(["50-78-2", "0-00-0"])
            {'50-78-2': 'DTXSID5020108', '0-00-0': None}
        """
        results = {}
        for casrn in casrn_list:
            results[casrn] = self.casrn_to_dtxsid(casrn)
        return results

    def batch_casrn_to_inchikey(
        self, casrn_list: List[str]
    ) -> Dict[str, Optional[str]]:
        """
        Convert multiple CAS numbers to InChIKeys.

        Args:
            casrn_list (list): List of CAS numbers

        Returns:
            dict: Mapping of CAS -> InChIKey (None if not found)

        Examples:
            >>> CompToxID().batch_casrn_to_inchikey(["50-78-2", "0-00-0"])
            {'50-78-2': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N', '0-00-0': None}
        """
        results = {}
        for casrn in casrn_list:
            results[casrn] = self.casrn_to_inchikey(casrn)
        return results

    def batch_inchikey_to_casrn(
        self, inchikey_list: List[str]
    ) -> Dict[str, Optional[str]]:
        """
        Convert multiple InChIKeys to CAS numbers.

        Args:
            inchikey_list (list): List of InChIKeys

        Returns:
            dict: Mapping of InChIKey -> CAS (None if not found)

        Examples:
            >>> CompToxID().batch_inchikey_to_casrn(
            ...     ["BSYNRYMUTXBXSQ-UHFFFAOYSA-N", "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"])
            {'BSYNRYMUTXBXSQ-UHFFFAOYSA-N': '50-78-2', 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N': '64-17-5'}
        """
        results = {}
        for inchikey in inchikey_list:
            results[inchikey] = self.inchikey_to_casrn(inchikey)
        return results
