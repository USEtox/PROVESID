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
import requests
from tqdm import tqdm
from .utils import data_path


class CompToxID:
    """
    Interface to CompTox Chemicals Dashboard SQLite database.

    The database file is automatically downloaded on first use when missing.
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
    ):
        """
        Initialize CompToxID database connection.

        Args:
            db_path (str, optional): Path to SQLite database. If None, uses default
                                    location in data directory.
            auto_download (bool, optional): If True, automatically download the
                database when missing (default: True).
            db_url (str, optional): Custom URL for database download. If None,
                uses the default Zenodo URL.

        Raises:
            FileNotFoundError: If database file doesn't exist and auto_download is False.
        """
        self.logger = logging.getLogger(__name__)

        if db_path is None:
            db_path = os.path.join(data_path(), self.DEFAULT_DB_NAME)

        self.db_path = db_path
        self.db_url = db_url or self.DEFAULT_DB_URL

        # Check if database exists
        if not os.path.exists(db_path):
            if auto_download:
                self.logger.warning(f"CompTox database not found at: {db_path}")
                self.logger.warning(
                    "The CompTox database is large (~856 MB). "
                    "Initial setup may take several minutes depending on your connection."
                )
                self.logger.warning(
                    f"Downloading CompTox database from: {self.db_url}"
                )
                self.download_database(url=self.db_url, force=False)
            else:
                raise FileNotFoundError(
                    f"CompTox database not found at: {db_path}\n"
                    "Database size: ~856 MB\n"
                    f"Run CompToxID.download_database() or set auto_download=True\n"
                    f"Download URL: {self.db_url}"
                )

        # Connect to database
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row  # Access columns by name

        # Verify the database has the expected table
        self._verify_database()

    def __del__(self):
        """Close database connection on deletion."""
        if hasattr(self, "conn"):
            self.conn.close()

    def download_database(self, url: Optional[str] = None, force: bool = False) -> str:
        """
        Download the CompTox SQLite database from Zenodo.

        The file is approximately 856 MB and is not shipped with the GitHub
        repository due to size limitations.

        Args:
            url (str, optional): Download URL. If None, uses `self.db_url`.
            force (bool, optional): If True, overwrite existing database file.

        Returns:
            str: Path to the downloaded database file.

        Raises:
            FileExistsError: If the database already exists and `force` is False.
            RuntimeError: If download or validation fails.
        """
        download_url = url or self.db_url

        if os.path.exists(self.db_path) and not force:
            raise FileExistsError(
                f"Database already exists at: {self.db_path}. Use force=True to overwrite."
            )

        os.makedirs(os.path.dirname(self.db_path), exist_ok=True)
        temp_path = f"{self.db_path}.tmp"

        self.logger.warning(
            "CompTox database download starting (~856 MB). "
            "Please ensure you have enough disk space and stable internet."
        )
        self.logger.warning(f"Source: {download_url}")
        self.logger.warning(f"Destination: {self.db_path}")

        try:
            response = requests.get(download_url, stream=True, timeout=60)
            response.raise_for_status()

            total_size = int(response.headers.get("content-length", 0))
            chunk_size = 1024 * 1024

            with open(temp_path, "wb") as f:
                if total_size > 0:
                    with tqdm(
                        total=total_size,
                        unit="B",
                        unit_scale=True,
                        desc="Downloading CompTox database",
                    ) as progress_bar:
                        for chunk in response.iter_content(chunk_size=chunk_size):
                            if chunk:
                                f.write(chunk)
                                progress_bar.update(len(chunk))
                else:
                    for chunk in response.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)

            if os.path.exists(self.db_path):
                os.remove(self.db_path)
            os.replace(temp_path, self.db_path)

            try:
                test_connection = sqlite3.connect(self.db_path)
                cursor = test_connection.cursor()
                cursor.execute(
                    "SELECT name FROM sqlite_master WHERE type='table' AND name='chemicals'"
                )
                if not cursor.fetchone():
                    raise RuntimeError("Downloaded database does not contain 'chemicals' table")
            finally:
                if "test_connection" in locals():
                    test_connection.close()

            self.logger.warning("CompTox database downloaded and validated successfully.")
            return self.db_path

        except requests.exceptions.RequestException as exc:
            if os.path.exists(temp_path):
                os.remove(temp_path)
            raise RuntimeError(f"Failed to download CompTox database: {exc}") from exc
        except Exception:
            if os.path.exists(temp_path):
                os.remove(temp_path)
            raise

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
