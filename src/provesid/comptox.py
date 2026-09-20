"""
CompToxID - Interface to CompTox Chemicals Dashboard SQLite database for fast identifier lookup and conversion.

This class provides access to a local SQLite database containing CompTox chemicals
with their identifiers (DTXSID, DTXCID, CASRN, InChIKey, SMILES, PREFERRED_NAME, etc.)
and chemical properties (molecular formula, average mass, monoisotopic mass, etc.).

The database is read from comptox_chemicals.db file.

Attributes:
    db_path (str): Path to the SQLite database file
    conn (sqlite3.Connection): Database connection

Example:
    >>> from provesid import CompToxID
    >>> db = CompToxID()
    >>>
    >>> # Lookup by CASRN
    >>> result = db.get_by_casrn("50-78-2")  # Aspirin
    >>> print(result['preferred_name'])
    >>>
    >>> # Lookup by InChIKey
    >>> result = db.get_by_inchikey("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
    >>> print(result['dtxsid'])
    >>>
    >>> # Convert CASRN to DTXSID
    >>> dtxsid = db.casrn_to_dtxsid("50-78-2")
    >>>
    >>> # Batch conversion
    >>> results = db.batch_casrn_to_dtxsid(["50-78-2", "50-00-0"])
"""

import os
import sqlite3
import logging

from typing import Dict, List, Optional, Any, Union

from .datasets import download_file
from .sqlite_client import SQLiteClient
from .utils import user_dataset_path


class CompToxID(SQLiteClient):
    """
    Interface to CompTox Chemicals Dashboard SQLite database.

    The database file is automatically downloaded on first use when missing.

    Inherits its connection handling from
    :class:`~provesid.sqlite_client.SQLiteClient`: use it as a context
    manager, or call :meth:`~provesid.sqlite_client.SQLiteClient.close` when
    finished, and query it from as many threads as you like --- each gets its
    own connection.

    Example:
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
        return self.db_path

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
            dict: Chemical information including identifiers and properties, or None if not found
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
            dict: Chemical information, or None if not found
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

    def get_by_inchikey(self, inchikey: str) -> Optional[Dict[str, Any]]:
        """
        Get chemical information by InChIKey.

        Args:
            inchikey (str): Standard InChIKey (27 characters)

        Returns:
            dict: Chemical information, or None if not found
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
            smiles (str): SMILES string

        Returns:
            dict: Chemical information, or None if not found
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
            name (str): Preferred name

        Returns:
            dict: Chemical information, or None if not found
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
            dict: Chemical information, or None if not found
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

        Args:
            name (str): Chemical name or synonym to search for
            exact (bool): If True, exact match only. If False, partial match (case-insensitive)
            limit (int): Maximum number of results to return

        Returns:
            list: List of matching chemicals
        """
        cursor = self.conn.cursor()
        results = []

        if exact:
            # Search in preferred name
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
            list: List of matching chemicals
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
        """Convert CAS Registry Number to DTXSID."""
        result = self.get_by_casrn(casrn)
        return result["DTXSID"] if result else None

    def casrn_to_inchikey(self, casrn: str) -> Optional[str]:
        """Convert CAS Registry Number to InChIKey."""
        result = self.get_by_casrn(casrn)
        return result["INCHIKEY"] if result else None

    def casrn_to_smiles(self, casrn: str) -> Optional[str]:
        """Convert CAS Registry Number to SMILES."""
        result = self.get_by_casrn(casrn)
        return result["SMILES"] if result else None

    def inchikey_to_casrn(self, inchikey: str) -> Optional[str]:
        """Convert InChIKey to CAS Registry Number."""
        result = self.get_by_inchikey(inchikey)
        return result["CASRN"] if result else None

    def inchikey_to_dtxsid(self, inchikey: str) -> Optional[str]:
        """Convert InChIKey to DTXSID."""
        result = self.get_by_inchikey(inchikey)
        return result["DTXSID"] if result else None

    def dtxsid_to_casrn(self, dtxsid: str) -> Optional[str]:
        """Convert DTXSID to CAS Registry Number."""
        result = self.get_by_dtxsid(dtxsid)
        return result["CASRN"] if result else None

    def dtxsid_to_inchikey(self, dtxsid: str) -> Optional[str]:
        """Convert DTXSID to InChIKey."""
        result = self.get_by_dtxsid(dtxsid)
        return result["INCHIKEY"] if result else None

    def dtxsid_to_smiles(self, dtxsid: str) -> Optional[str]:
        """Convert DTXSID to SMILES."""
        result = self.get_by_dtxsid(dtxsid)
        return result["SMILES"] if result else None

    def smiles_to_casrn(self, smiles: str) -> Optional[str]:
        """Convert SMILES to CAS Registry Number."""
        result = self.get_by_smiles(smiles)
        return result["CASRN"] if result else None

    def smiles_to_dtxsid(self, smiles: str) -> Optional[str]:
        """Convert SMILES to DTXSID."""
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
        """
        results = {}
        for casrn in casrn_list:
            results[casrn] = self.casrn_to_dtxsid(casrn)
        return results

    def batch_casrn_to_inchikey(
        self, casrn_list: List[str]
    ) -> Dict[str, Optional[str]]:
        """Convert multiple CAS numbers to InChIKeys."""
        results = {}
        for casrn in casrn_list:
            results[casrn] = self.casrn_to_inchikey(casrn)
        return results

    def batch_inchikey_to_casrn(
        self, inchikey_list: List[str]
    ) -> Dict[str, Optional[str]]:
        """Convert multiple InChIKeys to CAS numbers."""
        results = {}
        for inchikey in inchikey_list:
            results[inchikey] = self.inchikey_to_casrn(inchikey)
        return results
