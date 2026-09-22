"""
The ZeroPM global chemical inventory, offline: :class:`ZeroPM`.

ZeroPM (https://zeropm.eu) merged 25 national and regional chemical
inventories --- TSCA, the EC Inventory, Japan's CSCL, China's IECSC and
others --- into one SQLite file, resolved every listed CAS number and name to
structures, and scored the structures for persistence and mobility. Version
0.0.4 holds 164 513 CAS numbers, 283 104 names and 359 221 structures.

Its tables, and the names this class uses for them:

``api_ready_query``
    Every query string --- a CAS number or a name --- with its ``query_id``.
    Names are stored as the inventories spelled them, so matching is exact
    and case-sensitive.
``api_results``
    The structures (``inchi_id``) each query resolved to, with a ``rank``.
    Rank 1 is the resolver's best answer. Lower ranks are often a different
    compound entirely: the name ``formaldehyde`` reaches methane at rank 2.
``substances``
    One row per structure: ``inchi_id``, InChI and InChIKey.
``zeropm_chemicals``, ``pm_probabilities``
    The structures ZeroPM assessed (``zeropm_id``) and their probabilities
    of being persistent and mobile.

The methods that follow every rank (``get_cas_from_name``,
``get_cas_from_inchi`` and the ``get_id_table_from_*`` family) are broad
rather than precise; filter a table on ``rank == 1`` for the best answer. This
is why :class:`~provesid.search.Search` leaves ZeroPM out of its default vote.

SMILES are not stored. They are written from the InChI with RDKit on demand.

Example:
    >>> from provesid import ZeroPM
    >>> with ZeroPM() as zpm:
    ...     zpm.get_smiles_from_cas("64-17-5")
    ...     zpm.get_zeropm_id(cas="64-17-5")
    'CCO'
    1452
"""

import sqlite3
import os
from rdkit import Chem
import logging
from rapidfuzz import process, fuzz, utils
import pandas as pd
from typing import Optional

# Try relative import first, fall back to direct import for testing
from .datasets import download_file
from .sqlite_client import SQLiteClient
from .utils import user_dataset_path


class ZeroPM(SQLiteClient):
    """
    Class to extract data from the ZeroPM SQLite database using SQL queries.
    This class provides the same functionality as ZeroPM but uses SQL instead of pandas.
    SMILES are generated on-the-fly from InChI using RDKit.

    The database file will be automatically downloaded if not found locally.

    Connection handling comes from
    :class:`~provesid.sqlite_client.SQLiteClient`: use the class as a context
    manager, or call :meth:`~provesid.sqlite_client.SQLiteClient.close` when
    finished, and query it from as many threads as you like --- each gets its
    own connection.  This is the one client that also *writes*
    (:meth:`create_indexes`, :meth:`create_view`); SQLite serialises those
    against the reading connections.

    Example
    -------
    >>> with ZeroPM() as zpm:
    ...     df = zpm.get_id_table_from_cas("50-00-0")
    ...     df[["cas", "rank", "inchikey", "zeropm_id"]].to_dict("records")
    [{'cas': '50-00-0', 'rank': 1, 'inchikey': 'WSFSSNUMVMOOMR-UHFFFAOYSA-N', 'zeropm_id': 3224}]
    """

    # Default download URL for the ZeroPM database
    DEFAULT_DB_URL = "https://github.com/ZeroPM-H2020/global-chemical-inventory-database/raw/refs/heads/main/zeropm-v0-0-4.sqlite"

    def __init__(
        self,
        db_name: str = 'zeropm-v0-0-4.sqlite',
        auto_download: bool = True,
        db_url: Optional[str] = None,
        data_dir: Optional[str] = None,
        db_path: Optional[str] = None,
        redownload: bool = False,
    ):
        """
        Initialize connection to the ZeroPM SQLite database.

        Parameters
        ----------
        db_name : str, optional
            Name of the SQLite database file (default: 'zeropm-v0-0-4.sqlite')
        auto_download : bool, optional
            If True, automatically download the database if not found (default: True)
        db_url : str, optional
            Custom URL to download the database from. If None, uses the default GitHub URL.
        data_dir : str, optional
            Directory to store the database when ``db_path`` is not provided.
        db_path : str, optional
            Full path to a database file. Overrides ``db_name``/``data_dir``.
        redownload : bool, optional
            If True, force re-download when ``auto_download`` is enabled.

        Raises
        ------
        FileNotFoundError
            If the database is not on disk and ``auto_download`` is False.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> os.path.basename(zpm.db_path)
        'zeropm-v0-0-4.sqlite'
        >>> ZeroPM(db_path="/no/such/zeropm.sqlite", auto_download=False)
        Traceback (most recent call last):
        ...
        FileNotFoundError: Database not found at: /no/such/zeropm.sqlite
        ...
        """
        self.logger = logging.getLogger(__name__)
        if db_path is None:
            self.path = data_dir or user_dataset_path()
            self.db_path = os.path.join(self.path, db_name)
        else:
            self.db_path = os.path.abspath(os.path.expanduser(db_path))
            self.path = os.path.dirname(self.db_path)

        self.db_url = db_url or self.DEFAULT_DB_URL

        needs_download = redownload or not os.path.exists(self.db_path)

        # Check if database exists, download if needed
        if needs_download:
            if auto_download:
                if redownload and os.path.exists(self.db_path):
                    self.logger.info(
                        "Forced ZeroPM redownload requested for: %s", self.db_path
                    )
                else:
                    self.logger.info(f"Database not found at: {self.db_path}")
                self.logger.info("Downloading database automatically...")
                self.download_database(url=self.db_url, force=redownload)
            else:
                raise FileNotFoundError(
                    f"Database not found at: {self.db_path}\n"
                    f"Please run ZeroPM.download_database() or set auto_download=True"
                )

        # Create the connection.  One per thread, reused for every query on
        # that thread and released by close() or by leaving a ``with`` block.
        # row_factory stays unset: this module's queries index rows by
        # position, and sqlite3.Row would be a behaviour change.
        self._open_database(self.db_path, row_factory=None)

        # Cache chemical names for fuzzy matching (lazy loading)
        self._chemical_names_cache = None

    def download_database(self, url=None, force=False):
        """
        Download the ZeroPM SQLite database from a remote URL.

        The transfer is resumable: an interrupted download leaves a ``.part``
        file beside the destination and the next call continues from it rather
        than starting the 100 MB again. Nothing replaces an existing database
        until the new file has downloaded in full and opened successfully.

        Parameters
        ----------
        url : str, optional
            URL to download the database from. If None, uses the default GitHub URL.
        force : bool, optional
            If True, download even if the database already exists (default: False)

        Returns
        -------
        str
            Path to the downloaded database file

        Raises
        ------
        FileExistsError
            If the database already exists and force=False
        provesid.datasets.DownloadError
            If the download could not be completed, or the file that arrived is
            not a readable SQLite database

        Example
        -------
        >>> zpm = ZeroPM()
        >>> zpm.download_database(force=True)      # doctest: +SKIP
        '/home/me/.local/share/provesid/zeropm-v0-0-4.sqlite'
        """
        download_url = url or self.db_url

        # Check if database already exists
        if os.path.exists(self.db_path) and not force:
            raise FileExistsError(
                f"Database already exists at: {self.db_path}\n"
                f"Use force=True to overwrite"
            )

        def must_be_a_database(path):
            """Reject a download that is not a readable SQLite file.

            Run on the ``.part`` file, before it is moved into place, so a
            damaged download leaves any existing database untouched.
            """
            connection = sqlite3.connect(path)
            try:
                connection.execute(
                    "SELECT name FROM sqlite_master WHERE type='table' LIMIT 1"
                ).fetchone()
            except sqlite3.Error as exc:
                raise RuntimeError(f"Downloaded database is corrupted: {exc}") from exc
            finally:
                connection.close()

        download_file(
            download_url,
            self.db_path,
            verify=must_be_a_database,
            description="ZeroPM database",
            log=self.logger,
        )
        return self.db_path

    def _get_chemical_names_cache(self):
        """
        Lazy load and cache all chemical names for fuzzy matching.
        Returns a list of (name, query_id) tuples.
        """
        if self._chemical_names_cache is None:
            self.cursor.execute("""
                SELECT query, query_id
                FROM api_ready_query
                WHERE type = 'chemical name'
            """)
            self._chemical_names_cache = self.cursor.fetchall()
        return self._chemical_names_cache

    def _inchi_to_smiles(self, inchi):
        """
        Convert InChI string to SMILES using RDKit.

        Parameters
        ----------
        inchi : str
            InChI string

        Returns
        -------
        str or None
            SMILES string, or None if conversion fails
        """
        try:
            mol = Chem.MolFromInchi(inchi)
            if mol is None:
                return None
            return Chem.MolToSmiles(mol)
        except Exception as e:
            logging.warning(f"Error converting InChI to SMILES: {e}")
            return None

    def query_cas(self, cas_rn):
        """
        Returns a query id from the query with the CAS RN to be used with the query_results function.

        Parameters
        ----------
        cas_rn : str
            CAS Registry Number

        Returns
        -------
        int or None
            query_id if found, None otherwise

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> zpm.query_cas("50-00-0")
        8671
        >>> zpm.query_cas("0-00-0") is None
        True
        """
        self.cursor.execute("""
            SELECT query_id
            FROM api_ready_query
            WHERE query = ? AND type = 'CAS Registry Number'
        """, (cas_rn,))
        result = self.cursor.fetchone()
        return result[0] if result else None

    def query_name(self, name):
        """
        Returns a query id from the query with the exact chemical name to be used with the query_results function.

        Parameters
        ----------
        name : str
            Exact chemical name, case included: the inventories' spellings
            are separate queries.

        Returns
        -------
        int or None
            query_id if found, None otherwise

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> zpm.query_name("Formaldehyde"), zpm.query_name("formaldehyde")
        (8672, 325578)
        """
        self.cursor.execute("""
            SELECT query_id
            FROM api_ready_query
            WHERE query = ? AND type = 'chemical name'
        """, (name,))
        result = self.cursor.fetchone()
        return result[0] if result else None

    def query_similar_name(self, name, number_of_results=5, score_cutoff=80):
        """
        Returns number_of_results query ids from a query with similar chemical names
        using fuzzy string matching.

        Parameters
        ----------
        name : str
            Chemical name to search for
        number_of_results : int, optional
            Maximum number of results to return (default: 5)
        score_cutoff : int, optional
            Minimum similarity score (0-100) (default: 80)

        Returns
        -------
        list or None
            List of query_ids, or None if no matches above cutoff

        Notes
        -----
        Scores with ``rapidfuzz``'s ``WRatio``, which rates a short name
        highly whenever it appears inside the query. :meth:`match_similar_name`
        uses a stricter scorer and reports the names and scores.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> zpm.query_similar_name("formaldehyd")
        [8672, 8673, 104113, 325578, 367895]
        """
        names_cache = self._get_chemical_names_cache()
        name_list = [n[0] for n in names_cache]

        res = process.extract(
            name,
            name_list,
            scorer=fuzz.WRatio,
            limit=number_of_results,
            processor=utils.default_process,
        )

        if (len(res) < 1) or (res[0][1] < score_cutoff):
            return None
        else:
            # Get query_ids for matching names
            query_ids = []
            for match in res:
                if match[1] >= score_cutoff:
                    matched_name = match[0]
                    # Find the query_id for this name
                    query_id = next((n[1] for n in names_cache if n[0] == matched_name), None)
                    if query_id:
                        query_ids.append(query_id)
            return query_ids if query_ids else None

    def match_similar_name(self, name, number_of_results=5, score_cutoff=80,
                           scorer=None):
        """
        Fuzzy-match a chemical name and return the matches with their scores.

        Same purpose as :meth:`query_similar_name`, but keeps the matched name
        and the similarity score instead of discarding them, so callers can tell
        *what* matched and *how well*.

        Uses ``rapidfuzz.fuzz.ratio`` rather than the ``WRatio`` used by
        :meth:`query_similar_name`. ``WRatio`` includes a partial-ratio term
        that scores a short name highly whenever it appears anywhere inside the
        query, which over a list of millions of chemical names is a reliable
        source of nonsense: ``WRatio("caffiene", "ne")`` is 90 and
        ``WRatio("zzzznotachemical", "Mica")`` is also 90, while ``ratio`` puts
        both at 40 and still scores the genuine typo
        ``ratio("caffiene", "caffeine")`` at 87.5.

        Parameters
        ----------
        name : str
            Chemical name to search for.
        number_of_results : int, optional
            Maximum number of matches to return (default: 5).
        score_cutoff : int, optional
            Minimum similarity score, 0-100 (default: 80).
        scorer : callable, optional
            A ``rapidfuzz.fuzz`` scorer. Defaults to ``fuzz.ratio``. Pass
            ``fuzz.token_sort_ratio`` when word order may differ; avoid
            ``fuzz.WRatio`` for the reason above.

        Returns
        -------
        list of tuple
            ``(matched_name, query_id, score)`` tuples, best first. Empty when
            nothing scores at or above ``score_cutoff``.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> name, query_id, score = zpm.match_similar_name("formaldehyd")[0]
        >>> name, round(score, 1)
        ('Formaldehyde', 95.7)
        >>> zpm.match_similar_name("zzzznotachemical")
        []
        """
        names_cache = self._get_chemical_names_cache()
        query_id_of = {row[0]: row[1] for row in names_cache}

        matches = process.extract(
            name,
            list(query_id_of),
            scorer=scorer or fuzz.ratio,
            limit=number_of_results,
            processor=utils.default_process,
        )

        return [
            (matched_name, query_id_of[matched_name], score)
            for matched_name, score, _ in matches
            if score >= score_cutoff
        ]

    def get_id_table_from_similar_name(self, name, number_of_results=5, score_cutoff=80):
        """
        Returns identifiers for the chemical whose name best fuzzy-matches *name*.

        The fuzzy counterpart of :meth:`get_id_table_from_name`: use it when the
        name may be misspelled or formatted differently from the database entry.
        The table is built for the single best-scoring match.

        Parameters
        ----------
        name : str
            Chemical name, possibly misspelled.
        number_of_results : int, optional
            How many fuzzy candidates to consider (default: 5). Only the best
            one is turned into a table.
        score_cutoff : int, optional
            Minimum ``rapidfuzz.fuzz.ratio`` score, 0-100 (default: 80); see
            :meth:`match_similar_name`.

        Returns
        -------
        pandas.DataFrame or None
            Same columns as :meth:`get_id_table_from_name`, with an extra
            ``matched_name`` column recording what actually matched, and
            ``match_score`` holding its similarity. None when nothing scores at
            or above ``score_cutoff``.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> df = zpm.get_id_table_from_similar_name("formaldehyd")
        >>> row = df.iloc[0]
        >>> row["name"], row["matched_name"], round(float(row["match_score"]), 1), row["inchikey"]
        ('formaldehyd', 'Formaldehyde', 95.7, 'WSFSSNUMVMOOMR-UHFFFAOYSA-N')
        """
        matches = self.match_similar_name(
            name, number_of_results=number_of_results, score_cutoff=score_cutoff
        )
        if not matches:
            self.logger.debug("No fuzzy name match for '%s' at cutoff %s", name, score_cutoff)
            return None

        matched_name, query_id, score = matches[0]
        table = self._id_table_for_query_id(query_id, name)
        if table is None or table.empty:
            return None

        table["matched_name"] = matched_name
        table["match_score"] = score
        return table

    def get_inchi_id(self, query_id):
        """
        Returns all the inchi_id and ranks of a query with a given query_id.

        Parameters
        ----------
        query_id : int
            Query identifier

        Returns
        -------
        tuple of (list, list)
            (inchi_ids, ranks) sorted by rank, with duplicates removed; two
            empty lists when the query has no structures

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> zpm.get_inchi_id(zpm.query_cas("50-00-0"))
        ([32227], [1])
        >>> zpm.get_inchi_id(zpm.query_name("formaldehyde"))
        ([32227, 73275, 27053, 119941], [1, 2, 3, 4])
        """
        self.cursor.execute("""
            SELECT DISTINCT inchi_id, rank
            FROM api_results
            WHERE query_id = ?
            ORDER BY rank
        """, (query_id,))
        results = self.cursor.fetchall()

        if not results:
            return [], []

        # Separate inchi_ids and ranks
        inchi_ids = [r[0] for r in results]
        ranks = [r[1] for r in results]

        return inchi_ids, ranks

    def get_inchi(self, inchi_id):
        """
        Returns the inchi and inchikey string of a given inchi_id.

        Parameters
        ----------
        inchi_id : int
            InChI identifier

        Returns
        -------
        tuple of (str, str) or (None, None)
            (inchi, inchikey) if found, (None, None) otherwise

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> zpm.get_inchi(32227)
        ('InChI=1S/CH2O/c1-2/h1H2', 'WSFSSNUMVMOOMR-UHFFFAOYSA-N')
        """
        self.cursor.execute("""
            SELECT inchi, inchikey
            FROM substances
            WHERE inchi_id = ?
        """, (inchi_id,))
        result = self.cursor.fetchone()
        return (result[0], result[1]) if result else (None, None)

    def get_names(self, cas_rn):
        """
        Returns all the names for a CAS number.

        Parameters
        ----------
        cas_rn : str
            CAS Registry Number

        Returns
        -------
        list
            The distinct names the inventories list under this CAS number,
            excluding the CAS number itself, in no particular order. Empty
            when the CAS number is not in the database.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> sorted(zpm.get_names("64-17-5"))[:4]
        ['Alcohol', 'ETHANOL', 'ETHYL ALCOHOL', 'Ethanol']
        """
        query_id = self.query_cas(cas_rn)
        if query_id is None:
            return []

        # Get inventory_ids from inventory_summary
        self.cursor.execute("""
            SELECT inventory_id
            FROM inventory_summary
            WHERE query_id = ?
        """, (query_id,))
        inventory_ids = [row[0] for row in self.cursor.fetchall()]

        if len(inventory_ids) == 0:
            return []

        # Get identifiers from inventories
        names = set()
        for inv_id in inventory_ids:
            self.cursor.execute("""
                SELECT identifier
                FROM inventories
                WHERE inventory_id = ?
            """, (inv_id,))
            result = self.cursor.fetchone()
            if result:
                # Split by semicolon and add to set
                identifier_string = result[0]
                for name in identifier_string.split(';'):
                    name = name.strip()
                    if name and name != cas_rn:
                        names.add(name)

        return list(names)

    def get_smiles_from_cas(self, cas_rn):
        """
        Returns the SMILES from a CAS number.
        SMILES is generated on-the-fly from InChI using RDKit.

        Parameters
        ----------
        cas_rn : str
            CAS Registry Number

        Returns
        -------
        str or None
            SMILES string of the rank-1 structure, or None if not found

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> zpm.get_smiles_from_cas("50-00-0")
        'C=O'
        """
        query_id = self.query_cas(cas_rn)
        if query_id is None:
            return None

        # Get inchi_id from the query_id
        inchi_ids, _ = self.get_inchi_id(query_id)
        if len(inchi_ids) == 0:
            return None

        # Get InChI and convert to SMILES
        inchi, _ = self.get_inchi(inchi_ids[0])
        if inchi is None:
            return None

        return self._inchi_to_smiles(inchi)

    def get_cas_from_inchi(self, inchi):
        """
        Returns the CAS number(s) from an InChI string.

        Parameters
        ----------
        inchi : str
            InChI string

        Returns
        -------
        str, list, or None
            CAS number, list of CAS numbers, or None if not found. Every CAS
            number whose query reaches this structure at any rank, in no
            particular order: formaldehyde's list includes carbon monoxide's
            ``630-08-0``, which reaches it at rank 2.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> cas = zpm.get_cas_from_inchi("InChI=1S/CH2O/c1-2/h1H2")
        >>> "50-00-0" in cas, "630-08-0" in cas
        (True, True)
        >>> zpm.get_cas_from_inchi("InChI=1S/Xx") is None
        True
        """
        # First, find the inchi_id
        self.cursor.execute("""
            SELECT inchi_id
            FROM substances
            WHERE inchi = ?
        """, (inchi,))
        result = self.cursor.fetchone()
        if not result:
            return None

        inchi_id = result[0]

        # Find all query_ids for this inchi_id that are CAS numbers
        self.cursor.execute("""
            SELECT DISTINCT aq.query
            FROM api_results ar
            JOIN api_ready_query aq ON ar.query_id = aq.query_id
            WHERE ar.inchi_id = ? AND aq.type = 'CAS Registry Number'
        """, (inchi_id,))
        cas_numbers = [row[0] for row in self.cursor.fetchall()]

        if not cas_numbers:
            return None
        elif len(cas_numbers) == 1:
            return cas_numbers[0]
        else:
            return cas_numbers

    def get_cas_from_inchikey(self, inchikey):
        """
        Returns the CAS number(s) from an InChIKey.

        Parameters
        ----------
        inchikey : str
            InChIKey string

        Returns
        -------
        str, list, or None
            CAS number, list of CAS numbers, or None if not found; as broad
            as :meth:`get_cas_from_inchi`

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> "64-17-5" in zpm.get_cas_from_inchikey("LFQSCWFLJHTTHZ-UHFFFAOYSA-N")
        True
        """
        # First, find the inchi_id
        self.cursor.execute("""
            SELECT inchi_id
            FROM substances
            WHERE inchikey = ?
        """, (inchikey,))
        result = self.cursor.fetchone()
        if not result:
            return None

        inchi_id = result[0]

        # Find all query_ids for this inchi_id that are CAS numbers
        self.cursor.execute("""
            SELECT DISTINCT aq.query
            FROM api_results ar
            JOIN api_ready_query aq ON ar.query_id = aq.query_id
            WHERE ar.inchi_id = ? AND aq.type = 'CAS Registry Number'
        """, (inchi_id,))
        cas_numbers = [row[0] for row in self.cursor.fetchall()]

        if not cas_numbers:
            return None
        elif len(cas_numbers) == 1:
            return cas_numbers[0]
        else:
            return cas_numbers

    def get_smiles_from_inchikey(self, inchikey):
        """
        Returns the SMILES from an InChIKey.
        SMILES is generated on-the-fly from InChI using RDKit.

        Parameters
        ----------
        inchikey : str
            InChIKey string

        Returns
        -------
        str or None
            SMILES string, or None if not found

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> zpm.get_smiles_from_inchikey("LFQSCWFLJHTTHZ-UHFFFAOYSA-N")
        'CCO'
        """
        # Get InChI from InChIKey
        self.cursor.execute("""
            SELECT inchi
            FROM substances
            WHERE inchikey = ?
        """, (inchikey,))
        result = self.cursor.fetchone()

        if not result:
            return None

        inchi = result[0]
        return self._inchi_to_smiles(inchi)

    def get_cas_from_smiles(self, smiles):
        """
        Returns the CAS number from a SMILES string.
        This is done by converting the SMILES to InChI and then to CAS number.

        Parameters
        ----------
        smiles : str
            SMILES string

        Returns
        -------
        str, list, or None
            CAS number, list of CAS numbers, or None if not found or the
            SMILES cannot be parsed; as broad as :meth:`get_cas_from_inchi`

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> "64-17-5" in zpm.get_cas_from_smiles("OCC")
        True
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logging.warning(f"Invalid SMILES: {smiles}")
                return None
            inchi = Chem.MolToInchi(mol)
        except Exception as e:
            logging.warning(f"Error converting SMILES to InChI for smiles: {smiles}. Error: {e}")
            return None

        return self.get_cas_from_inchi(inchi)

    def get_cas_from_name(self, name):
        """
        Returns the CAS number(s) associated with a chemical name.

        This method performs an exact match search for the chemical name in the database.
        For fuzzy matching, use query_similar_name() first to get query_ids.

        The answer is broad: it is every CAS number that reaches any of the
        structures the name resolved to, at any rank. For
        ``"formaldehyde"`` that is 31 numbers, methane's and carbon's among
        them. :meth:`get_id_table_from_name` shows where each came from.

        Parameters
        ----------
        name : str
            Chemical name (exact match)

        Returns
        -------
        str, list, or None
            CAS number, list of CAS numbers, or None if not found

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> cas = zpm.get_cas_from_name("formaldehyde")
        >>> len(cas), "50-00-0" in cas, "74-82-8" in cas
        (31, True, True)
        >>> zpm.get_cas_from_name("acetylsalicylic acid") is None
        True
        """
        # Get query_id for this name
        query_id = self.query_name(name)
        if query_id is None:
            return None

        # Get inchi_ids for this query_id
        inchi_ids, _ = self.get_inchi_id(query_id)
        if not inchi_ids:
            return None

        # Collect all CAS numbers for all inchi_ids
        all_cas = set()
        for inchi_id in inchi_ids:
            # Get InChI for this inchi_id
            inchi, _ = self.get_inchi(inchi_id)
            if inchi:
                cas_result = self.get_cas_from_inchi(inchi)
                if cas_result:
                    if isinstance(cas_result, list):
                        all_cas.update(cas_result)
                    else:
                        all_cas.add(cas_result)

        if not all_cas:
            return None
        elif len(all_cas) == 1:
            return list(all_cas)[0]
        else:
            return sorted(list(all_cas))

    def get_cas_from_formula(self, formula):
        """
        Returns CAS numbers for chemicals matching a molecular formula.

        Note: Molecular formulas are not unique identifiers - many different chemicals
        can have the same formula (isomers). This method returns all CAS numbers
        for chemicals matching the given formula.

        Parameters
        ----------
        formula : str
            Molecular formula (e.g., "H2O", "C6H12O6", "CH2O")

        Returns
        -------
        list or None
            List of CAS numbers matching the formula, or None if not found

        Warning
        -------
        This method can be slow as it needs to parse all InChI strings to extract
        molecular formulas. Consider caching results for frequently used formulas.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> zpm.get_cas_from_formula("CH2O")  # Formaldehyde
        ['108-62-3', '1664-98-8', '30525-89-4', '3228-27-1', '50-00-0', '630-08-0', '63101-50-8']
        """
        # Normalize formula (basic normalization - can be improved)
        formula = formula.replace(" ", "")

        # Query all substances and check their formulas
        # InChI format: InChI=1S/CH2O/c1-2/h1H2
        # Formula is between the first two slashes
        self.cursor.execute("""
            SELECT DISTINCT s.inchi_id, s.inchi
            FROM substances s
            WHERE s.inchi IS NOT NULL
        """)

        matching_inchi_ids = []
        for inchi_id, inchi in self.cursor.fetchall():
            try:
                # Extract formula from InChI
                # Format: InChI=1S/FORMULA/...
                parts = inchi.split('/')
                if len(parts) >= 2:
                    inchi_formula = parts[1]
                    if inchi_formula == formula:
                        matching_inchi_ids.append(inchi_id)
            except Exception:
                continue

        if not matching_inchi_ids:
            return None

        # Get all CAS numbers for matching inchi_ids
        all_cas = set()
        for inchi_id in matching_inchi_ids:
            self.cursor.execute("""
                SELECT DISTINCT aq.query
                FROM api_results ar
                JOIN api_ready_query aq ON ar.query_id = aq.query_id
                WHERE ar.inchi_id = ? AND aq.type = 'CAS Registry Number'
            """, (inchi_id,))
            cas_results = [row[0] for row in self.cursor.fetchall()]
            all_cas.update(cas_results)

        return sorted(list(all_cas)) if all_cas else None

    def batch_get_cas_from_smiles(self, smiles_list):
        """
        Get CAS numbers for multiple SMILES strings at once.

        Parameters
        ----------
        smiles_list : list of str
            List of SMILES strings

        Returns
        -------
        dict
            Dictionary mapping SMILES strings to CAS numbers (or None if not found)

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> zpm.batch_get_cas_from_smiles(["CC", "not a smiles"])
        {'CC': ['74-84-0', '9002-88-4'], 'not a smiles': None}
        """
        return {smiles: self.get_cas_from_smiles(smiles) for smiles in smiles_list}

    def batch_get_cas_from_name(self, name_list):
        """
        Get CAS numbers for multiple chemical names at once.

        Parameters
        ----------
        name_list : list of str
            List of chemical names (exact match)

        Returns
        -------
        dict
            Dictionary mapping chemical names to CAS numbers (or None if not found)

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> results = zpm.batch_get_cas_from_name(["Formaldehyde", "xyzzy"])
        >>> "50-00-0" in results["Formaldehyde"], results["xyzzy"]
        (True, None)
        """
        return {name: self.get_cas_from_name(name) for name in name_list}

    def batch_get_cas_from_formula(self, formula_list):
        """
        Get CAS numbers for multiple molecular formulas at once.

        Parameters
        ----------
        formula_list : list of str
            List of molecular formulas

        Returns
        -------
        dict
            Dictionary mapping formulas to lists of CAS numbers

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> results = zpm.batch_get_cas_from_formula(["CH2O", "C2H6O"])
        >>> {formula: len(cas) for formula, cas in results.items()}
        {'CH2O': 7, 'C2H6O': 10}
        """
        return {formula: self.get_cas_from_formula(formula) for formula in formula_list}

    def get_id_table_from_cas(self, cas):
        """
        Returns a pandas DataFrame containing all identifiers for a given CAS number.

        This method retrieves all query_ids associated with the CAS number, then for each query_id,
        it retrieves all associated inchi_ids and their corresponding InChI and InChIKey values.
        Synonyms (chemical names) and data sources are also included.

        Parameters
        ----------
        cas : str
            CAS Registry Number

        Returns
        -------
        pandas.DataFrame
            DataFrame with columns: 'cas', 'query_id', 'inchi_id', 'rank', 'inchi', 'inchikey', 'zeropm_id', 'synonyms', 'sources'
            Returns None if the CAS number is not found in the database.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> df = zpm.get_id_table_from_cas("50-00-0")
        >>> df[["cas", "query_id", "inchi_id", "rank", "inchikey", "zeropm_id"]]
               cas  query_id  inchi_id  rank                     inchikey  zeropm_id
        0  50-00-0      8671     32227     1  WSFSSNUMVMOOMR-UHFFFAOYSA-N       3224
        >>> df.loc[0, "sources"]
        'Chemical Data Reporting Inventory, Industrial ...'
        """
        # Get all query_ids for this CAS (using fetchall in case there are multiple)
        self.cursor.execute("""
            SELECT query_id
            FROM api_ready_query
            WHERE query = ? AND type = 'CAS Registry Number'
        """, (cas,))
        query_ids = [row[0] for row in self.cursor.fetchall()]

        if not query_ids:
            self.logger.debug("CAS number %s not found in database", cas)
            return None

        # Get synonyms for this CAS
        synonyms = self.get_names(cas)
        synonyms_str = "; ".join(synonyms) if synonyms else ""

        # Get sources for this CAS
        self.cursor.execute("""
            SELECT DISTINCT s.source_name
            FROM inventory_summary issum
            JOIN inventories inv ON issum.inventory_id = inv.inventory_id
            JOIN sources s ON inv.source_id = s.source_id
            WHERE issum.query_id IN ({})
        """.format(','.join('?' * len(query_ids))), query_ids)
        sources = [row[0] for row in self.cursor.fetchall()]
        sources_str = "; ".join(sources) if sources else ""

        # Collect all data
        rows = []
        for query_id in query_ids:
            # Get all inchi_ids for this query_id
            inchi_ids, ranks = self.get_inchi_id(query_id)

            if not inchi_ids:
                # If no inchi_ids found, still add a row with the query_id
                rows.append({
                    'cas': cas,
                    'query_id': query_id,
                    'inchi_id': None,
                    'rank': None,
                    'inchi': None,
                    'inchikey': None,
                    'zeropm_id': None,
                    'synonyms': synonyms_str,
                    'sources': sources_str
                })
            else:
                # For each inchi_id, get the inchi and inchikey
                for inchi_id, rank in zip(inchi_ids, ranks):
                    inchi, inchikey = self.get_inchi(inchi_id)
                    # Get zeropm_id for this inchi_id
                    self.cursor.execute("""
                        SELECT zeropm_id
                        FROM zeropm_chemicals
                        WHERE inchi_id = ?
                    """, (inchi_id,))
                    zeropm_result = self.cursor.fetchone()
                    zeropm_id = zeropm_result[0] if zeropm_result else None

                    rows.append({
                        'cas': cas,
                        'query_id': query_id,
                        'inchi_id': inchi_id,
                        'rank': rank,
                        'inchi': inchi,
                        'inchikey': inchikey,
                        'zeropm_id': zeropm_id,
                        'synonyms': synonyms_str,
                        'sources': sources_str
                    })

        # Create DataFrame
        df = pd.DataFrame(rows)
        # Convert zeropm_id to nullable integer type
        if not df.empty and 'zeropm_id' in df.columns:
            df['zeropm_id'] = df['zeropm_id'].astype('Int64')
        return df

    def get_id_table_from_zeropm_id(self, zeropm_id):
        """
        Returns a pandas DataFrame containing all identifiers for a given zeropm_id.

        This method retrieves the inchi_id associated with the zeropm_id, then finds all
        query_ids (CAS numbers) linked to that inchi_id and builds a comprehensive table
        with InChI, InChIKey, synonyms, and data sources.

        Parameters
        ----------
        zeropm_id : int
            ZeroPM identifier

        Returns
        -------
        pandas.DataFrame
            DataFrame with columns: 'cas', 'query_id', 'inchi_id', 'rank', 'inchi', 'inchikey', 'zeropm_id', 'synonyms', 'sources'
            Returns None if the zeropm_id is not found in the database.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> df = zpm.get_id_table_from_zeropm_id(3224)   # formaldehyde
        >>> df[["cas", "rank"]].sort_values(["rank", "cas"]).values.tolist()
        [['30525-89-4', 1], ['50-00-0', 1], ['108-62-3', 2], ['1664-98-8', 2], ['630-08-0', 2], ['63101-50-8', 2]]
        """
        # Get inchi_id for this zeropm_id
        self.cursor.execute("""
            SELECT inchi_id
            FROM zeropm_chemicals
            WHERE zeropm_id = ?
        """, (zeropm_id,))
        result = self.cursor.fetchone()

        if not result:
            self.logger.debug("zeropm_id %s not found in database", zeropm_id)
            return None

        inchi_id = result[0]

        # Get InChI and InChIKey
        inchi, inchikey = self.get_inchi(inchi_id)

        # Get all query_ids (CAS numbers) associated with this inchi_id.
        # DISTINCT because api_results can hold the same (query, structure,
        # rank) more than once, differing only in columns not read here.
        self.cursor.execute("""
            SELECT DISTINCT ar.query_id, ar.rank, aq.query
            FROM api_results ar
            JOIN api_ready_query aq ON ar.query_id = aq.query_id
            WHERE ar.inchi_id = ? AND aq.type = 'CAS Registry Number'
        """, (inchi_id,))
        query_results = self.cursor.fetchall()

        if not query_results:
            logging.warning(f"No CAS numbers found for zeropm_id {zeropm_id}")
            return None

        # Collect all data
        rows = []
        for query_id, rank, cas in query_results:
            # Get synonyms for this CAS
            synonyms = self.get_names(cas)
            synonyms_str = "; ".join(synonyms) if synonyms else ""

            # Get sources for this query_id
            self.cursor.execute("""
                SELECT DISTINCT s.source_name
                FROM inventory_summary issum
                JOIN inventories inv ON issum.inventory_id = inv.inventory_id
                JOIN sources s ON inv.source_id = s.source_id
                WHERE issum.query_id = ?
            """, (query_id,))
            sources = [row[0] for row in self.cursor.fetchall()]
            sources_str = "; ".join(sources) if sources else ""

            rows.append({
                'cas': cas,
                'query_id': query_id,
                'inchi_id': inchi_id,
                'rank': rank,
                'inchi': inchi,
                'inchikey': inchikey,
                'zeropm_id': zeropm_id,
                'synonyms': synonyms_str,
                'sources': sources_str
            })

        # Create DataFrame
        df = pd.DataFrame(rows)
        # Convert zeropm_id to nullable integer type
        if not df.empty and 'zeropm_id' in df.columns:
            df['zeropm_id'] = df['zeropm_id'].astype('Int64')
        return df

    def batch_get_id_table_from_cas(self, cas_list):
        """
        Returns a pandas DataFrame containing all identifiers for a list of CAS numbers.

        This method calls get_id_table_from_cas for each CAS number in the list and
        combines the results into a single DataFrame. CAS numbers not found in the
        database are logged but skipped in the output.

        Parameters
        ----------
        cas_list : list of str
            List of CAS Registry Numbers

        Returns
        -------
        pandas.DataFrame
            Combined DataFrame with columns: 'cas', 'query_id', 'inchi_id', 'rank', 'inchi', 'inchikey', 'zeropm_id', 'synonyms', 'sources'
            Returns an empty DataFrame if no CAS numbers are found in the database.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> cas_numbers = ["50-00-0", "50-78-2", "64-17-5"]  # formaldehyde, aspirin, ethanol
        >>> df = zpm.batch_get_id_table_from_cas(cas_numbers)
        >>> df[["cas", "rank", "inchikey", "zeropm_id"]]
               cas  rank                     inchikey  zeropm_id
        0  50-00-0     1  WSFSSNUMVMOOMR-UHFFFAOYSA-N       3224
        1  50-78-2     1  BSYNRYMUTXBXSQ-UHFFFAOYSA-N       4267
        2  50-78-2     2  BSYNRYMUTXBXSQ-UHFFFAOYSA-M       <NA>
        3  50-78-2     3  XDZMPRGFOOFSBL-UHFFFAOYSA-N       6402
        4  50-78-2     4  BSYNRYMUTXBXSQ-FIBGUPNXSA-N       <NA>
        5  64-17-5     1  LFQSCWFLJHTTHZ-UHFFFAOYSA-N       1452
        """
        if not cas_list:
            logging.warning("Empty CAS list provided")
            return pd.DataFrame(columns=['cas', 'query_id', 'inchi_id', 'rank', 'inchi', 'inchikey', 'zeropm_id', 'synonyms', 'sources'])

        # Collect DataFrames for each CAS
        dataframes = []
        for cas in cas_list:
            df = self.get_id_table_from_cas(cas)
            if df is not None:
                dataframes.append(df)

        # Combine all DataFrames
        if not dataframes:
            logging.warning("None of the provided CAS numbers were found in the database")
            return pd.DataFrame(columns=['cas', 'query_id', 'inchi_id', 'rank', 'inchi', 'inchikey', 'zeropm_id', 'synonyms', 'sources'])

        # Concatenate all DataFrames and reset index
        combined_df = pd.concat(dataframes, ignore_index=True)
        return combined_df

    def batch_get_id_table_from_cas_filtered(self, cas_list, rank=None, have_zeropm_id=None):
        """
        Returns a filtered pandas DataFrame containing identifiers for a list of CAS numbers.

        This method calls batch_get_id_table_from_cas and applies optional filters to the results.

        Parameters
        ----------
        cas_list : list of str
            List of CAS Registry Numbers
        rank : int, optional
            If specified, only include rows with this rank value (e.g., rank=1 for top results)
            If None, no rank filtering is applied (default: None)
        have_zeropm_id : bool, optional
            If True, only include rows where zeropm_id is not None
            If False, only include rows where zeropm_id is None
            If None, no zeropm_id filtering is applied (default: None)

        Returns
        -------
        pandas.DataFrame
            Filtered DataFrame with columns: 'cas', 'query_id', 'inchi_id', 'rank', 'inchi', 'inchikey', 'zeropm_id', 'synonyms', 'sources'
            Returns an empty DataFrame if no CAS numbers match the filters.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> cas_numbers = ["50-00-0", "50-78-2", "64-17-5"]
        >>> # Get only rank=1 results with zeropm_id
        >>> df = zpm.batch_get_id_table_from_cas_filtered(cas_numbers, rank=1, have_zeropm_id=True)
        >>> df[["cas", "rank", "zeropm_id"]]
               cas  rank  zeropm_id
        0  50-00-0     1       3224
        1  50-78-2     1       4267
        2  64-17-5     1       1452

        See Also
        --------
        batch_get_id_table_from_cas : Returns all results without filtering
        """
        # Get the full id table
        df = self.batch_get_id_table_from_cas(cas_list)

        # Return empty if no results
        if df.empty:
            return df

        # Apply rank filter if specified
        if rank is not None:
            df = df[df['rank'] == rank]

        # Apply zeropm_id filter if specified
        if have_zeropm_id is not None:
            if have_zeropm_id:
                df = df[df['zeropm_id'].notna()]
            else:
                df = df[df['zeropm_id'].isna()]

        # Reset index
        df = df.reset_index(drop=True)

        return df

    def get_id_table_from_inchi(self, inchi):
        """
        Returns a pandas DataFrame containing all identifiers for a given InChI.

        This method retrieves the inchi_id for the InChI, then finds all associated
        query_ids and their CAS numbers. It also includes synonyms and sources.

        Parameters
        ----------
        inchi : str
            InChI string

        Returns
        -------
        pandas.DataFrame
            DataFrame with columns: 'inchi', 'inchikey', 'inchi_id', 'query_id', 'rank', 'cas', 'sources', 'synonyms'
            Returns None if the InChI is not found in the database.
            One row per query --- CAS number or name --- that reaches the
            structure, best rank first; ``cas`` is NaN for a name query. The
            synonyms are those of the first CAS number, on every row.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> df = zpm.get_id_table_from_inchi("InChI=1S/CH2O/c1-2/h1H2")
        >>> df.dropna(subset=["cas"])[["query_id", "rank", "cas"]].head(2)
           query_id  rank         cas
        0      8671     1     50-00-0
        3     35725     1  30525-89-4
        """
        # Get inchi_id and inchikey from InChI
        self.cursor.execute("""
            SELECT inchi_id, inchikey
            FROM substances
            WHERE inchi = ?
        """, (inchi,))
        result = self.cursor.fetchone()

        if not result:
            self.logger.debug("InChI %s not found in database", inchi)
            return None

        inchi_id, inchikey = result

        # Get all query_ids and ranks for this inchi_id
        self.cursor.execute("""
            SELECT DISTINCT ar.query_id, ar.rank
            FROM api_results ar
            WHERE ar.inchi_id = ?
            ORDER BY ar.rank
        """, (inchi_id,))
        query_results = self.cursor.fetchall()

        if not query_results:
            # If no query_ids found, still return basic info
            return pd.DataFrame([{
                'inchi': inchi,
                'inchikey': inchikey,
                'inchi_id': inchi_id,
                'query_id': None,
                'rank': None,
                'cas': None,
                'synonyms': '',
                'sources': ''
            }])

        # Get CAS numbers for these query_ids
        rows = []
        primary_cas = None
        query_ids_list = [q[0] for q in query_results]

        # Get sources for all query_ids at once
        if query_ids_list:
            placeholders = ','.join('?' * len(query_ids_list))
            self.cursor.execute(f"""
                SELECT DISTINCT s.source_name
                FROM inventory_summary issum
                JOIN inventories inv ON issum.inventory_id = inv.inventory_id
                JOIN sources s ON inv.source_id = s.source_id
                WHERE issum.query_id IN ({placeholders})
            """, query_ids_list)
            sources = [row[0] for row in self.cursor.fetchall()]
            sources_str = "; ".join(sources) if sources else ""
        else:
            sources_str = ""

        for query_id, rank in query_results:
            # Get CAS number for this query_id
            self.cursor.execute("""
                SELECT query
                FROM api_ready_query
                WHERE query_id = ? AND type = 'CAS Registry Number'
            """, (query_id,))
            cas_result = self.cursor.fetchone()
            cas = cas_result[0] if cas_result else None

            # Use first CAS as primary for synonyms
            if cas and primary_cas is None:
                primary_cas = cas

            rows.append({
                'inchi': inchi,
                'inchikey': inchikey,
                'inchi_id': inchi_id,
                'query_id': query_id,
                'rank': rank,
                'cas': cas,
                'sources': sources_str
            })

        # Get synonyms from primary CAS
        synonyms_str = ''
        if primary_cas:
            synonyms = self.get_names(primary_cas)
            synonyms_str = "; ".join(synonyms) if synonyms else ""

        # Add synonyms to all rows
        for row in rows:
            row['synonyms'] = synonyms_str

        return pd.DataFrame(rows)

    def batch_get_id_table_from_inchi(self, inchi_list):
        """
        Returns a pandas DataFrame containing all identifiers for a list of InChI strings.

        Parameters
        ----------
        inchi_list : list of str
            List of InChI strings

        Returns
        -------
        pandas.DataFrame
            Combined DataFrame with columns: 'inchi', 'inchikey', 'inchi_id', 'query_id', 'rank', 'cas', 'sources', 'synonyms'
            Returns an empty DataFrame if no InChIs are found in the database.
            InChIs not found are left out.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> df = zpm.batch_get_id_table_from_inchi(
        ...     ["InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3", "InChI=1S/Xx"])
        >>> df["inchikey"].unique().tolist(), len(df)
        (['LFQSCWFLJHTTHZ-UHFFFAOYSA-N'], 43)
        """
        if not inchi_list:
            logging.warning("Empty InChI list provided")
            return pd.DataFrame(columns=['inchi', 'inchikey', 'inchi_id', 'query_id', 'rank', 'cas', 'synonyms', 'sources'])

        dataframes = []
        for inchi in inchi_list:
            df = self.get_id_table_from_inchi(inchi)
            if df is not None:
                dataframes.append(df)

        if not dataframes:
            logging.warning("None of the provided InChIs were found in the database")
            return pd.DataFrame(columns=['inchi', 'inchikey', 'inchi_id', 'query_id', 'rank', 'cas', 'synonyms', 'sources'])

        combined_df = pd.concat(dataframes, ignore_index=True)
        return combined_df

    def get_id_table_from_inchikey(self, inchikey):
        """
        Returns a pandas DataFrame containing all identifiers for a given InChIKey.

        This method retrieves the inchi_id for the InChIKey, then finds all associated
        query_ids and their CAS numbers. It also includes synonyms and sources.

        Parameters
        ----------
        inchikey : str
            InChIKey string

        Returns
        -------
        pandas.DataFrame
            DataFrame with columns: 'inchikey', 'inchi', 'inchi_id', 'query_id', 'rank', 'cas', 'sources', 'synonyms'
            Returns None if the InChIKey is not found in the database.
            Shaped as :meth:`get_id_table_from_inchi` describes.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> df = zpm.get_id_table_from_inchikey("WSFSSNUMVMOOMR-UHFFFAOYSA-N")
        >>> df.dropna(subset=["cas"])[["rank", "cas"]].head(2)
           rank         cas
        0     1     50-00-0
        3     1  30525-89-4
        """
        # Get inchi_id and inchi from InChIKey
        self.cursor.execute("""
            SELECT inchi_id, inchi
            FROM substances
            WHERE inchikey = ?
        """, (inchikey,))
        result = self.cursor.fetchone()

        if not result:
            self.logger.debug("InChIKey %s not found in database", inchikey)
            return None

        inchi_id, inchi = result

        # Get all query_ids and ranks for this inchi_id
        self.cursor.execute("""
            SELECT DISTINCT ar.query_id, ar.rank
            FROM api_results ar
            WHERE ar.inchi_id = ?
            ORDER BY ar.rank
        """, (inchi_id,))
        query_results = self.cursor.fetchall()

        if not query_results:
            # If no query_ids found, still return basic info
            return pd.DataFrame([{
                'inchikey': inchikey,
                'inchi': inchi,
                'inchi_id': inchi_id,
                'query_id': None,
                'rank': None,
                'cas': None,
                'synonyms': '',
                'sources': ''
            }])

        # Get CAS numbers for these query_ids
        rows = []
        primary_cas = None
        query_ids_list = [q[0] for q in query_results]

        # Get sources for all query_ids at once
        if query_ids_list:
            placeholders = ','.join('?' * len(query_ids_list))
            self.cursor.execute(f"""
                SELECT DISTINCT s.source_name
                FROM inventory_summary issum
                JOIN inventories inv ON issum.inventory_id = inv.inventory_id
                JOIN sources s ON inv.source_id = s.source_id
                WHERE issum.query_id IN ({placeholders})
            """, query_ids_list)
            sources = [row[0] for row in self.cursor.fetchall()]
            sources_str = "; ".join(sources) if sources else ""
        else:
            sources_str = ""

        for query_id, rank in query_results:
            # Get CAS number for this query_id
            self.cursor.execute("""
                SELECT query
                FROM api_ready_query
                WHERE query_id = ? AND type = 'CAS Registry Number'
            """, (query_id,))
            cas_result = self.cursor.fetchone()
            cas = cas_result[0] if cas_result else None

            # Use first CAS as primary for synonyms
            if cas and primary_cas is None:
                primary_cas = cas

            rows.append({
                'inchikey': inchikey,
                'inchi': inchi,
                'inchi_id': inchi_id,
                'query_id': query_id,
                'rank': rank,
                'cas': cas,
                'sources': sources_str
            })

        # Get synonyms from primary CAS
        synonyms_str = ''
        if primary_cas:
            synonyms = self.get_names(primary_cas)
            synonyms_str = "; ".join(synonyms) if synonyms else ""

        # Add synonyms to all rows
        for row in rows:
            row['synonyms'] = synonyms_str

        return pd.DataFrame(rows)

    def batch_get_id_table_from_inchikey(self, inchikey_list):
        """
        Returns a pandas DataFrame containing all identifiers for a list of InChIKey strings.

        Parameters
        ----------
        inchikey_list : list of str
            List of InChIKey strings

        Returns
        -------
        pandas.DataFrame
            Combined DataFrame with columns: 'inchikey', 'inchi', 'inchi_id', 'query_id', 'rank', 'cas', 'sources', 'synonyms'
            Returns an empty DataFrame if no InChIKeys are found in the database.
            InChIKeys not found are left out.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> df = zpm.batch_get_id_table_from_inchikey(
        ...     ["LFQSCWFLJHTTHZ-UHFFFAOYSA-N", "XXXXXXXXXXXXXX-XXXXXXXXXX-X"])
        >>> df["inchikey"].unique().tolist()
        ['LFQSCWFLJHTTHZ-UHFFFAOYSA-N']
        """
        if not inchikey_list:
            logging.warning("Empty InChIKey list provided")
            return pd.DataFrame(columns=['inchikey', 'inchi', 'inchi_id', 'query_id', 'rank', 'cas', 'synonyms', 'sources'])

        dataframes = []
        for inchikey in inchikey_list:
            df = self.get_id_table_from_inchikey(inchikey)
            if df is not None:
                dataframes.append(df)

        if not dataframes:
            logging.warning("None of the provided InChIKeys were found in the database")
            return pd.DataFrame(columns=['inchikey', 'inchi', 'inchi_id', 'query_id', 'rank', 'cas', 'synonyms', 'sources'])

        combined_df = pd.concat(dataframes, ignore_index=True)
        return combined_df

    def get_id_table_from_name(self, name):
        """
        Returns a pandas DataFrame containing all identifiers for a given chemical name.

        This method searches for an exact match of the chemical name, then retrieves all
        associated inchi_ids and their corresponding InChI, InChIKey, CAS numbers, and sources.

        Parameters
        ----------
        name : str
            Chemical name (exact match)

        Returns
        -------
        pandas.DataFrame
            DataFrame with columns: 'name', 'query_id', 'inchi_id', 'rank', 'inchi', 'inchikey', 'cas', 'sources'
            Returns None if the name is not found in the database.
            One row per (structure, CAS number): each structure the name
            resolved to, at every rank, with every CAS number that reaches
            that structure.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> df = zpm.get_id_table_from_name("Formaldehyde")
        >>> df.groupby("rank")["inchikey"].first().to_dict()
        {1: 'WSFSSNUMVMOOMR-UHFFFAOYSA-N', 2: 'VNWKTOKETHGBQD-UHFFFAOYSA-N', 3: 'MDYZKJNTKZIUSK-UHFFFAOYSA-N', 4: 'SYCNHFWYTQQMNG-UHFFFAOYSA-N'}
        >>> df[df["rank"] == 1]["cas"].tolist()[:2]
        ['50-00-0', '30525-89-4']
        """
        # Get query_id for this name
        query_id = self.query_name(name)

        if query_id is None:
            self.logger.debug("Chemical name '%s' not found in database", name)
            return None

        return self._id_table_for_query_id(query_id, name)

    def _id_table_for_query_id(self, query_id, name):
        """
        Build the identifier table for one already-resolved query_id.

        This is the shared body of :meth:`get_id_table_from_name` and
        :meth:`get_id_table_from_similar_name`; the only difference between them
        is how the ``query_id`` was found.

        Parameters
        ----------
        query_id : int
            An ``api_ready_query`` id.
        name : str
            Name to record in the ``name`` column of the result.

        Returns
        -------
        pandas.DataFrame
            DataFrame with columns: 'name', 'query_id', 'inchi_id', 'rank',
            'inchi', 'inchikey', 'cas', 'sources'.
        """
        # Get sources for this query_id
        self.cursor.execute("""
            SELECT DISTINCT s.source_name
            FROM inventory_summary issum
            JOIN inventories inv ON issum.inventory_id = inv.inventory_id
            JOIN sources s ON inv.source_id = s.source_id
            WHERE issum.query_id = ?
        """, (query_id,))
        sources = [row[0] for row in self.cursor.fetchall()]
        sources_str = "; ".join(sources) if sources else ""

        # Get all inchi_ids and ranks for this query_id
        inchi_ids, ranks = self.get_inchi_id(query_id)

        if not inchi_ids:
            # If no inchi_ids found, still return basic info
            return pd.DataFrame([{
                'name': name,
                'query_id': query_id,
                'inchi_id': None,
                'rank': None,
                'inchi': None,
                'inchikey': None,
                'cas': None,
                'sources': sources_str
            }])

        # Collect all data
        rows = []
        for inchi_id, rank in zip(inchi_ids, ranks):
            # Get InChI and InChIKey
            inchi, inchikey = self.get_inchi(inchi_id)

            # Get CAS number(s) for this inchi_id
            self.cursor.execute("""
                SELECT DISTINCT aq.query
                FROM api_results ar
                JOIN api_ready_query aq ON ar.query_id = aq.query_id
                WHERE ar.inchi_id = ? AND aq.type = 'CAS Registry Number'
            """, (inchi_id,))
            cas_results = [row[0] for row in self.cursor.fetchall()]

            # If multiple CAS numbers, create a row for each
            if cas_results:
                for cas in cas_results:
                    rows.append({
                        'name': name,
                        'query_id': query_id,
                        'inchi_id': inchi_id,
                        'rank': rank,
                        'inchi': inchi,
                        'inchikey': inchikey,
                        'cas': cas,
                        'sources': sources_str
                    })
            else:
                # No CAS found, still add the row
                rows.append({
                    'name': name,
                    'query_id': query_id,
                    'inchi_id': inchi_id,
                    'rank': rank,
                    'inchi': inchi,
                    'inchikey': inchikey,
                    'cas': None,
                    'sources': sources_str
                })

        return pd.DataFrame(rows)

    def batch_get_id_table_from_name(self, name_list):
        """
        Returns a pandas DataFrame containing all identifiers for a list of chemical names.

        Parameters
        ----------
        name_list : list of str
            List of chemical names (exact match)

        Returns
        -------
        pandas.DataFrame
            Combined DataFrame with columns: 'name', 'query_id', 'inchi_id', 'rank', 'inchi', 'inchikey', 'cas', 'sources'
            Returns an empty DataFrame if no names are found in the database.
            Names not found are left out.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> df = zpm.batch_get_id_table_from_name(["Formaldehyde", "ethanol", "xyzzy"])
        >>> df[df["rank"] == 1].groupby("name")["cas"].first().to_dict()
        {'Formaldehyde': '50-00-0', 'ethanol': '64-17-5'}
        """
        if not name_list:
            logging.warning("Empty name list provided")
            return pd.DataFrame(columns=['name', 'query_id', 'inchi_id', 'rank', 'inchi', 'inchikey', 'cas', 'sources'])

        dataframes = []
        for name in name_list:
            df = self.get_id_table_from_name(name)
            if df is not None:
                dataframes.append(df)

        if not dataframes:
            logging.warning("None of the provided names were found in the database")
            return pd.DataFrame(columns=['name', 'query_id', 'inchi_id', 'rank', 'inchi', 'inchikey', 'cas', 'sources'])

        combined_df = pd.concat(dataframes, ignore_index=True)
        return combined_df

    # ==================== Performance Enhancement Methods ====================

    def create_indexes(self, force=False):
        """
        Create indexes on frequently queried columns to improve performance.
        Indexes are created on query, type, query_id, inchi_id, inchi, and inchikey.

        Parameters
        ----------
        force : bool, optional
            If True, drop existing indexes before creating new ones (default: False)

        Returns
        -------
        dict
            Dictionary with index names as keys and status ('created', 'exists', 'error') as values.
            Without ``force`` every index reads ``'exists'``, whether or not
            it was just built: ``CREATE INDEX IF NOT EXISTS`` does not say.

        Notes
        -----
        This writes to the database file. An index that already exists under
        its name is left alone, so a second call does nothing.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> zpm.create_indexes()["idx_query"]     # doctest: +SKIP
        'exists'
        """
        indexes = {
            'idx_query': 'CREATE INDEX IF NOT EXISTS idx_query ON api_ready_query(query)',
            'idx_type': 'CREATE INDEX IF NOT EXISTS idx_type ON api_ready_query(type)',
            'idx_query_id_results': 'CREATE INDEX IF NOT EXISTS idx_query_id_results ON api_results(query_id)',
            'idx_inchi_id_results': 'CREATE INDEX IF NOT EXISTS idx_inchi_id_results ON api_results(inchi_id)',
            'idx_inchi': 'CREATE INDEX IF NOT EXISTS idx_inchi ON substances(inchi)',
            'idx_inchikey': 'CREATE INDEX IF NOT EXISTS idx_inchikey ON substances(inchikey)',
            'idx_inventory_query': 'CREATE INDEX IF NOT EXISTS idx_inventory_query ON inventory_summary(query_id)',
            'idx_inventory_id': 'CREATE INDEX IF NOT EXISTS idx_inventory_id ON inventories(inventory_id)',
        }

        results = {}

        if force:
            # Drop existing indexes
            for idx_name in indexes.keys():
                try:
                    self.cursor.execute(f"DROP INDEX IF EXISTS {idx_name}")
                except Exception as e:
                    logging.warning(f"Could not drop index {idx_name}: {e}")

        # Create indexes
        for idx_name, sql in indexes.items():
            try:
                self.cursor.execute(sql)
                self.conn.commit()
                results[idx_name] = 'created' if force else 'exists'
            except Exception as e:
                logging.error(f"Error creating index {idx_name}: {e}")
                results[idx_name] = 'error'

        return results

    # ==================== Batch Query Methods ====================

    def batch_query_cas(self, cas_list):
        """
        Query multiple CAS numbers at once.

        Parameters
        ----------
        cas_list : list of str
            List of CAS Registry Numbers

        Returns
        -------
        dict
            Dictionary mapping CAS numbers to query_ids (or None if not found)

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> zpm.batch_query_cas(["50-00-0", "64-17-5", "0-00-0"])
        {'50-00-0': 8671, '64-17-5': 3904, '0-00-0': None}
        """
        if not cas_list:
            return {}

        # Use parameterized query with IN clause
        placeholders = ','.join('?' * len(cas_list))
        self.cursor.execute(f"""
            SELECT query, query_id
            FROM api_ready_query
            WHERE query IN ({placeholders}) AND type = 'CAS Registry Number'
        """, cas_list)

        results = {row[0]: row[1] for row in self.cursor.fetchall()}

        # Add None for CAS numbers not found
        return {cas: results.get(cas) for cas in cas_list}

    def batch_get_smiles_from_cas(self, cas_list):
        """
        Get SMILES for multiple CAS numbers at once.

        Parameters
        ----------
        cas_list : list of str
            List of CAS Registry Numbers

        Returns
        -------
        dict
            Dictionary mapping CAS numbers to SMILES strings (or None if not found)

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> zpm.batch_get_smiles_from_cas(["50-00-0", "64-17-5", "0-00-0"])
        {'50-00-0': 'C=O', '64-17-5': 'CCO', '0-00-0': None}
        """
        query_ids = self.batch_query_cas(cas_list)
        results = {}

        for cas, query_id in query_ids.items():
            if query_id is None:
                results[cas] = None
            else:
                results[cas] = self.get_smiles_from_cas(cas)

        return results

    def batch_get_names(self, cas_list):
        """
        Get all names for multiple CAS numbers at once.

        Parameters
        ----------
        cas_list : list of str
            List of CAS Registry Numbers

        Returns
        -------
        dict
            Dictionary mapping CAS numbers to lists of names, as
            :meth:`get_names` returns them (empty when not found)

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> names = zpm.batch_get_names(["64-17-5", "0-00-0"])
        >>> "Ethanol" in names["64-17-5"], names["0-00-0"]
        (True, [])
        """
        return {cas: self.get_names(cas) for cas in cas_list}

    def batch_get_cas_from_inchikey(self, inchikey_list):
        """
        Get CAS numbers for multiple InChIKeys at once.

        Parameters
        ----------
        inchikey_list : list of str
            List of InChIKey strings

        Returns
        -------
        dict
            Dictionary mapping InChIKeys to CAS numbers (or None if not found);
            one number as a string, several as a list, as broad as
            :meth:`get_cas_from_inchi`

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> found = zpm.batch_get_cas_from_inchikey(
        ...     ["WSFSSNUMVMOOMR-UHFFFAOYSA-N", "XXXXXXXXXXXXXX-XXXXXXXXXX-X"])
        >>> "50-00-0" in found["WSFSSNUMVMOOMR-UHFFFAOYSA-N"]
        True
        >>> found["XXXXXXXXXXXXXX-XXXXXXXXXX-X"] is None
        True
        """
        if not inchikey_list:
            return {}

        # First, get inchi_ids for all inchikeys
        placeholders = ','.join('?' * len(inchikey_list))
        self.cursor.execute(f"""
            SELECT inchikey, inchi_id
            FROM substances
            WHERE inchikey IN ({placeholders})
        """, inchikey_list)

        inchikey_to_id = {row[0]: row[1] for row in self.cursor.fetchall()}

        # Get all CAS numbers for these inchi_ids
        if not inchikey_to_id:
            return {key: None for key in inchikey_list}

        inchi_ids = list(inchikey_to_id.values())
        placeholders = ','.join('?' * len(inchi_ids))
        self.cursor.execute(f"""
            SELECT DISTINCT ar.inchi_id, aq.query
            FROM api_results ar
            JOIN api_ready_query aq ON ar.query_id = aq.query_id
            WHERE ar.inchi_id IN ({placeholders}) AND aq.type = 'CAS Registry Number'
        """, inchi_ids)

        # Group CAS numbers by inchi_id
        inchi_to_cas = {}
        for inchi_id, cas in self.cursor.fetchall():
            if inchi_id not in inchi_to_cas:
                inchi_to_cas[inchi_id] = []
            inchi_to_cas[inchi_id].append(cas)

        # Map back to inchikeys
        results = {}
        for inchikey in inchikey_list:
            inchi_id = inchikey_to_id.get(inchikey)
            if inchi_id and inchi_id in inchi_to_cas:
                cas_list = inchi_to_cas[inchi_id]
                results[inchikey] = cas_list[0] if len(cas_list) == 1 else cas_list
            else:
                results[inchikey] = None

        return results

    # ==================== Advanced Search Methods ====================

    def query_name_regex(self, pattern, case_sensitive=False, limit=100):
        """
        Search for chemical names with a simple wildcard pattern.

        Not a full regular expression: ``.*`` matches any run of characters
        and ``.`` any single character, and everything else is literal. The
        pattern is translated to SQL ``LIKE`` (case-insensitive) or ``GLOB``
        (case-sensitive), so it must match the whole name.

        Parameters
        ----------
        pattern : str
            Pattern using ``.*`` and ``.`` as wildcards
        case_sensitive : bool, optional
            Whether the search is case-sensitive (default: False)
        limit : int, optional
            Maximum number of results to return (default: 100)

        Returns
        -------
        list of tuple
            List of (query_id, name) tuples matching the pattern, in database
            order

        Note
        ----
        Use '.*pattern.*' for substring matching. A case-insensitive pattern
        may also use ``%`` and ``_``, which ``LIKE`` reads as wildcards.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> zpm.query_name_regex("formaldehyde.*", limit=2)
        [(8672, 'Formaldehyde'), (8673, 'formaldehyde ... %')]
        >>> zpm.query_name_regex("formaldehyde.*", case_sensitive=True, limit=2)
        [(8673, 'formaldehyde ... %'), (104113, 'formaldehyde ...%')]
        """
        if case_sensitive:
            # LIKE ignores case for ASCII letters whatever the pattern says;
            # GLOB does not, and takes * and ? as its wildcards.
            pattern = pattern.replace('.*', '*').replace('.', '?')
            self.cursor.execute(f"""
                SELECT query_id, query
                FROM api_ready_query
                WHERE type = 'chemical name' AND query GLOB ?
                LIMIT ?
            """, (pattern, limit))
        else:
            # Case-insensitive search
            pattern = pattern.replace('.*', '%').replace('.', '_')
            self.cursor.execute(f"""
                SELECT query_id, query
                FROM api_ready_query
                WHERE type = 'chemical name' AND LOWER(query) LIKE LOWER(?)
                LIMIT ?
            """, (pattern, limit))

        return self.cursor.fetchall()

    def get_cas_by_substructure(self, smarts_pattern, max_results=100):
        """
        Search for chemicals containing a specific substructure.
        This method converts InChIs to molecules and performs substructure
        matching using RDKit, in database order.

        Parameters
        ----------
        smarts_pattern : str
            SMARTS pattern for substructure search
        max_results : int, optional
            Maximum number of results to return (default: 100)

        Returns
        -------
        list of dict
            List of dictionaries with keys: 'cas', 'inchi', 'inchikey', 'smiles'.
            ``cas`` is as :meth:`get_cas_from_inchi` returns it. Empty for an
            invalid SMARTS pattern.

        Warning
        -------
        Only the first 10 000 of the database's ~359 000 structures are
        searched, so a structure beyond them is never found. Converting each
        InChI costs time, and RDKit logs a warning for many of them.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> hits = zpm.get_cas_by_substructure("c1ccccc1C(=O)O", max_results=2)
        >>> [hit["smiles"] for hit in hits]
        ['COc1ccc(C(=O)O)cc1', 'O=C(O)c1ccc(C(=O)O)cc1']
        """
        try:
            pattern_mol = Chem.MolFromSmarts(smarts_pattern)
            if pattern_mol is None:
                logging.error(f"Invalid SMARTS pattern: {smarts_pattern}")
                return []
        except Exception as e:
            logging.error(f"Error parsing SMARTS pattern: {e}")
            return []

        # Get all substances (this could be optimized with pagination)
        self.cursor.execute("""
            SELECT s.inchi_id, s.inchi, s.inchikey
            FROM substances s
            LIMIT 10000
        """)

        results = []
        count = 0

        for inchi_id, inchi, inchikey in self.cursor.fetchall():
            if count >= max_results:
                break

            # Convert InChI to mol
            try:
                mol = Chem.MolFromInchi(inchi)
                if mol is None:
                    continue

                # Check for substructure match
                if mol.HasSubstructMatch(pattern_mol):
                    # Get CAS number
                    cas = self.get_cas_from_inchi(inchi)
                    smiles = Chem.MolToSmiles(mol)

                    results.append({
                        'cas': cas,
                        'inchi': inchi,
                        'inchikey': inchikey,
                        'smiles': smiles
                    })
                    count += 1
            except Exception as e:
                continue

        return results

    # ==================== Export Methods ====================

    def export_to_csv(self, query_results, filename, columns=None):
        """
        Export query results to a CSV file.

        Parameters
        ----------
        query_results : list or dict
            Query results to export (list of tuples or dictionary)
        filename : str
            Output CSV filename. A relative name is written into the
            database's directory, beside the database; pass an absolute path
            to write anywhere else.
        columns : list of str, optional
            Column names for the CSV header. A dict gets ``key,value`` when
            none are given; a list gets no header.

        Returns
        -------
        str
            Path to the created CSV file

        Examples
        --------
        >>> import tempfile
        >>> zpm = ZeroPM()
        >>> path = os.path.join(tempfile.mkdtemp(), "smiles.csv")
        >>> zpm.export_to_csv({"50-00-0": "C=O"}, path, columns=["cas", "smiles"]) == path
        True
        >>> print(open(path).read())
        cas,smiles
        50-00-0,C=O
        <BLANKLINE>
        """
        import csv

        output_path = os.path.join(self.path, filename)

        with open(output_path, 'w', newline='', encoding='utf-8') as f:
            if isinstance(query_results, dict):
                # Handle dictionary results
                writer = csv.writer(f)
                if columns:
                    writer.writerow(columns)
                else:
                    writer.writerow(['key', 'value'])

                for key, value in query_results.items():
                    writer.writerow([key, value])
            else:
                # Handle list of tuples/lists
                writer = csv.writer(f)
                if columns:
                    writer.writerow(columns)

                for row in query_results:
                    writer.writerow(row)

        return output_path

    def create_view(self, view_name, sql_query):
        """
        Create a custom view in the database for frequently used queries.

        Parameters
        ----------
        view_name : str
            Name of the view to create
        sql_query : str
            SQL SELECT statement defining the view

        Returns
        -------
        bool
            True if view was created successfully, False otherwise. A view
            of the same name is replaced.

        Notes
        -----
        This writes to the database file.

        Example
        -------
        >>> zpm = ZeroPM()
        >>> sql = '''
        ...     SELECT aq.query AS cas, s.inchi, s.inchikey
        ...     FROM api_ready_query aq
        ...     JOIN api_results ar ON aq.query_id = ar.query_id
        ...     JOIN substances s ON ar.inchi_id = s.inchi_id
        ...     WHERE aq.type = 'CAS Registry Number' AND ar.rank = 1
        ... '''
        >>> zpm.create_view('cas_to_inchi', sql)            # doctest: +SKIP
        True
        """
        try:
            # Drop view if it exists
            self.cursor.execute(f"DROP VIEW IF EXISTS {view_name}")

            # Create new view
            self.cursor.execute(f"CREATE VIEW {view_name} AS {sql_query}")
            self.conn.commit()

            logging.info(f"View '{view_name}' created successfully")
            return True
        except Exception as e:
            logging.error(f"Error creating view '{view_name}': {e}")
            return False

    def export_query_results(self, sql_query, filename, include_headers=True):
        """
        Execute a custom SQL query and export results to CSV.

        Parameters
        ----------
        sql_query : str
            SQL query to execute
        filename : str
            Output CSV filename
        include_headers : bool, optional
            Include column headers in CSV (default: True)

        Returns
        -------
        str
            Path to the created CSV file; see :meth:`export_to_csv` for
            where a relative ``filename`` goes

        Examples
        --------
        >>> import tempfile
        >>> zpm = ZeroPM()
        >>> path = os.path.join(tempfile.mkdtemp(), "regions.csv")
        >>> _ = zpm.export_query_results(
        ...     "SELECT region_id, region FROM global_regions ORDER BY region_id", path)
        >>> print(open(path).read().splitlines()[:3])
        ['region_id,region', '1,North America', '2,Europe']
        """
        import csv

        self.cursor.execute(sql_query)
        results = self.cursor.fetchall()

        # Get column names from cursor description
        columns = [desc[0] for desc in self.cursor.description] if include_headers else None

        return self.export_to_csv(results, filename, columns)

    def get_database_stats(self):
        """
        Get statistics about the database contents.

        Returns
        -------
        dict
            Row counts of ``api_ready_query``, ``api_results``,
            ``substances``, ``inventories``, ``inventory_summary``,
            ``cleanventory_chemicals``, ``zeropm_chemicals``, ``components``
            and ``multi_components``, plus ``unique_cas_numbers`` and
            ``unique_chemical_names``. A table that cannot be counted holds
            its error message instead.

        Examples
        --------
        >>> stats = ZeroPM().get_database_stats()
        >>> stats["unique_cas_numbers"], stats["zeropm_chemicals"]
        (164513, 126369)
        """
        tables = [
            'api_ready_query', 'api_results', 'substances',
            'inventories', 'inventory_summary', 'cleanventory_chemicals',
            'zeropm_chemicals', 'components', 'multi_components'
        ]

        stats = {}

        for table in tables:
            try:
                self.cursor.execute(f"SELECT COUNT(*) FROM {table}")
                count = self.cursor.fetchone()[0]
                stats[table] = count
            except Exception as e:
                stats[table] = f"Error: {e}"

        # Additional statistics
        self.cursor.execute("""
            SELECT COUNT(DISTINCT query)
            FROM api_ready_query
            WHERE type = 'CAS Registry Number'
        """)
        stats['unique_cas_numbers'] = self.cursor.fetchone()[0]

        self.cursor.execute("""
            SELECT COUNT(DISTINCT query)
            FROM api_ready_query
            WHERE type = 'chemical name'
        """)
        stats['unique_chemical_names'] = self.cursor.fetchone()[0]

        return stats

    # ==================== Inventory, Country, and Region Query Methods ====================

    def get_all_inventories(self):
        """
        Get all available inventory sources.

        Returns
        -------
        list of dict
            List of dictionaries with keys: 'source_id', 'source_name', 'country_scope', 'link', 'type',
            ordered by name. Some names carry stray spaces, as stored.

        Examples
        --------
        >>> inventories = ZeroPM().get_all_inventories()
        >>> len(inventories)
        25
        >>> [(i["source_id"], i["country_scope"]) for i in inventories if "TSCA" in i["source_name"]]
        [(24, 'United States of America')]
        """
        self.cursor.execute("""
            SELECT source_id, source_name, country_scope, link, type
            FROM sources
            ORDER BY source_name
        """)

        inventories = []
        for row in self.cursor.fetchall():
            inventories.append({
                'source_id': row[0],
                'source_name': row[1],
                'country_scope': row[2],
                'link': row[3],
                'type': row[4]
            })

        return inventories

    def get_all_countries(self):
        """
        Get all countries in the database.

        Returns
        -------
        list of dict
            List of dictionaries with keys: 'country_id', 'country', ordered by name

        Examples
        --------
        >>> countries = ZeroPM().get_all_countries()
        >>> len(countries), countries[0]
        (38, {'country_id': 1, 'country': 'Australia'})
        """
        self.cursor.execute("""
            SELECT country_id, country
            FROM countries
            ORDER BY country
        """)

        countries = []
        for row in self.cursor.fetchall():
            countries.append({
                'country_id': row[0],
                'country': row[1]
            })

        return countries

    def get_all_regions(self):
        """
        Get all global regions in the database.

        Returns
        -------
        list of dict
            List of dictionaries with keys: 'region_id', 'region', ordered by name

        Examples
        --------
        >>> [r["region"] for r in ZeroPM().get_all_regions()]
        ['Asia', 'Europe', 'North America', 'Oceania', 'Scandinavia']
        """
        self.cursor.execute("""
            SELECT region_id, region
            FROM global_regions
            ORDER BY region
        """)

        regions = []
        for row in self.cursor.fetchall():
            regions.append({
                'region_id': row[0],
                'region': row[1]
            })

        return regions

    def query_by_inventory(self, source_name=None, source_id=None):
        """
        Query chemicals by inventory source.

        Parameters
        ----------
        source_name : str, optional
            Name of the inventory source (case-insensitive partial match)
        source_id : int, optional
            Source ID (exact match)

        Returns
        -------
        list of dict
            List of chemicals with keys: 'cas', 'query_id', 'inchi_id', 'source_name',
            ordered by CAS number. A CAS number appears once per structure
            it resolves to, at any rank, and once per matching inventory.

        Raises
        ------
        ValueError
            If neither ``source_name`` nor ``source_id`` is given.

        Note
        ----
        Either source_name or source_id must be provided.
        :meth:`count_chemicals_by_inventory` counts distinct CAS numbers
        without building the list.

        Examples
        --------
        >>> rows = ZeroPM().query_by_inventory(source_name="TSCA")
        >>> rows[0]
        {'cas': '100-00-5', 'query_id': 1927, 'inchi_id': 1, 'source_name': 'Toxic Substances Control Act (TSCA) Chemical Substance Inventory'}
        """
        if source_name is None and source_id is None:
            raise ValueError("Either source_name or source_id must be provided")

        if source_id is not None:
            # Query by source_id
            self.cursor.execute("""
                SELECT DISTINCT aq.query, aq.query_id, ar.inchi_id, s.source_name
                FROM api_ready_query aq
                JOIN inventory_summary issum ON aq.query_id = issum.query_id
                JOIN inventories inv ON issum.inventory_id = inv.inventory_id
                JOIN sources s ON inv.source_id = s.source_id
                JOIN api_results ar ON aq.query_id = ar.query_id
                WHERE aq.type = 'CAS Registry Number' AND s.source_id = ?
                ORDER BY aq.query
            """, (source_id,))
        else:
            # Query by source_name (partial, case-insensitive)
            self.cursor.execute("""
                SELECT DISTINCT aq.query, aq.query_id, ar.inchi_id, s.source_name
                FROM api_ready_query aq
                JOIN inventory_summary issum ON aq.query_id = issum.query_id
                JOIN inventories inv ON issum.inventory_id = inv.inventory_id
                JOIN sources s ON inv.source_id = s.source_id
                JOIN api_results ar ON aq.query_id = ar.query_id
                WHERE aq.type = 'CAS Registry Number' AND LOWER(s.source_name) LIKE LOWER(?)
                ORDER BY aq.query
            """, (f'%{source_name}%',))

        results = []
        for row in self.cursor.fetchall():
            results.append({
                'cas': row[0],
                'query_id': row[1],
                'inchi_id': row[2],
                'source_name': row[3]
            })

        return results

    def query_by_country(self, country_name=None, country_id=None):
        """
        Query chemicals by country.

        Parameters
        ----------
        country_name : str, optional
            Name of the country (case-insensitive partial match)
        country_id : int, optional
            Country ID (exact match)

        Returns
        -------
        list of dict
            List of chemicals with keys: 'cas', 'query_id', 'inchi_id', 'country', 'source_name',
            ordered by CAS number, repeated as :meth:`query_by_inventory`
            describes

        Raises
        ------
        ValueError
            If neither ``country_name`` nor ``country_id`` is given.

        Note
        ----
        Either country_name or country_id must be provided.

        Examples
        --------
        >>> rows = ZeroPM().query_by_country("Japan")
        >>> rows[0]["cas"], rows[0]["source_name"]
        ('100-00-5', 'NITE')
        """
        if country_name is None and country_id is None:
            raise ValueError("Either country_name or country_id must be provided")

        if country_id is not None:
            # Query by country_id
            self.cursor.execute("""
                SELECT DISTINCT aq.query, aq.query_id, ar.inchi_id, c.country, s.source_name
                FROM api_ready_query aq
                JOIN inventory_summary issum ON aq.query_id = issum.query_id
                JOIN inventories inv ON issum.inventory_id = inv.inventory_id
                JOIN sources s ON inv.source_id = s.source_id
                JOIN country_sources_index csi ON s.source_id = csi.source_id
                JOIN countries c ON csi.country_id = c.country_id
                JOIN api_results ar ON aq.query_id = ar.query_id
                WHERE aq.type = 'CAS Registry Number' AND c.country_id = ?
                ORDER BY aq.query
            """, (country_id,))
        else:
            # Query by country_name (partial, case-insensitive)
            self.cursor.execute("""
                SELECT DISTINCT aq.query, aq.query_id, ar.inchi_id, c.country, s.source_name
                FROM api_ready_query aq
                JOIN inventory_summary issum ON aq.query_id = issum.query_id
                JOIN inventories inv ON issum.inventory_id = inv.inventory_id
                JOIN sources s ON inv.source_id = s.source_id
                JOIN country_sources_index csi ON s.source_id = csi.source_id
                JOIN countries c ON csi.country_id = c.country_id
                JOIN api_results ar ON aq.query_id = ar.query_id
                WHERE aq.type = 'CAS Registry Number' AND LOWER(c.country) LIKE LOWER(?)
                ORDER BY aq.query
            """, (f'%{country_name}%',))

        results = []
        for row in self.cursor.fetchall():
            results.append({
                'cas': row[0],
                'query_id': row[1],
                'inchi_id': row[2],
                'country': row[3],
                'source_name': row[4]
            })

        return results

    def query_by_region(self, region_name=None, region_id=None):
        """
        Query chemicals by global region.

        Parameters
        ----------
        region_name : str, optional
            Name of the region (case-insensitive partial match)
        region_id : int, optional
            Region ID (exact match)

        Returns
        -------
        list of dict
            List of chemicals with keys: 'cas', 'query_id', 'inchi_id', 'region', 'country', 'source_name',
            ordered by CAS number, repeated as :meth:`query_by_inventory`
            describes

        Raises
        ------
        ValueError
            If neither ``region_name`` nor ``region_id`` is given.

        Note
        ----
        Either region_name or region_id must be provided.

        Examples
        --------
        >>> rows = ZeroPM().query_by_region("Oceania")
        >>> rows[0]["cas"], rows[0]["country"]
        ('100-00-5', 'New Zealand')
        """
        if region_name is None and region_id is None:
            raise ValueError("Either region_name or region_id must be provided")

        if region_id is not None:
            # Query by region_id
            self.cursor.execute("""
                SELECT DISTINCT aq.query, aq.query_id, ar.inchi_id, gr.region, c.country, s.source_name
                FROM api_ready_query aq
                JOIN inventory_summary issum ON aq.query_id = issum.query_id
                JOIN inventories inv ON issum.inventory_id = inv.inventory_id
                JOIN sources s ON inv.source_id = s.source_id
                JOIN country_sources_index csi ON s.source_id = csi.source_id
                JOIN countries c ON csi.country_id = c.country_id
                JOIN region_country_index rci ON c.country_id = rci.country_id
                JOIN global_regions gr ON rci.region_id = gr.region_id
                JOIN api_results ar ON aq.query_id = ar.query_id
                WHERE aq.type = 'CAS Registry Number' AND gr.region_id = ?
                ORDER BY aq.query
            """, (region_id,))
        else:
            # Query by region_name (partial, case-insensitive)
            self.cursor.execute("""
                SELECT DISTINCT aq.query, aq.query_id, ar.inchi_id, gr.region, c.country, s.source_name
                FROM api_ready_query aq
                JOIN inventory_summary issum ON aq.query_id = issum.query_id
                JOIN inventories inv ON issum.inventory_id = inv.inventory_id
                JOIN sources s ON inv.source_id = s.source_id
                JOIN country_sources_index csi ON s.source_id = csi.source_id
                JOIN countries c ON csi.country_id = c.country_id
                JOIN region_country_index rci ON c.country_id = rci.country_id
                JOIN global_regions gr ON rci.region_id = gr.region_id
                JOIN api_results ar ON aq.query_id = ar.query_id
                WHERE aq.type = 'CAS Registry Number' AND LOWER(gr.region) LIKE LOWER(?)
                ORDER BY aq.query
            """, (f'%{region_name}%',))

        results = []
        for row in self.cursor.fetchall():
            results.append({
                'cas': row[0],
                'query_id': row[1],
                'inchi_id': row[2],
                'region': row[3],
                'country': row[4],
                'source_name': row[5]
            })

        return results

    def get_countries_for_region(self, region_name=None, region_id=None):
        """
        Get all countries in a specific region.

        Parameters
        ----------
        region_name : str, optional
            Name of the region (case-insensitive partial match)
        region_id : int, optional
            Region ID (exact match)

        Returns
        -------
        list of dict
            List of dictionaries with keys: 'country_id', 'country', 'region'

        Raises
        ------
        ValueError
            If neither ``region_name`` nor ``region_id`` is given.

        Note
        ----
        Either region_name or region_id must be provided.

        Examples
        --------
        >>> [c["country"] for c in ZeroPM().get_countries_for_region("Scandinavia")]
        ['Denmark', 'Finland', 'Norway', 'Sweden']
        """
        if region_name is None and region_id is None:
            raise ValueError("Either region_name or region_id must be provided")

        if region_id is not None:
            self.cursor.execute("""
                SELECT DISTINCT c.country_id, c.country, gr.region
                FROM countries c
                JOIN region_country_index rci ON c.country_id = rci.country_id
                JOIN global_regions gr ON rci.region_id = gr.region_id
                WHERE gr.region_id = ?
                ORDER BY c.country
            """, (region_id,))
        else:
            self.cursor.execute("""
                SELECT DISTINCT c.country_id, c.country, gr.region
                FROM countries c
                JOIN region_country_index rci ON c.country_id = rci.country_id
                JOIN global_regions gr ON rci.region_id = gr.region_id
                WHERE LOWER(gr.region) LIKE LOWER(?)
                ORDER BY c.country
            """, (f'%{region_name}%',))

        countries = []
        for row in self.cursor.fetchall():
            countries.append({
                'country_id': row[0],
                'country': row[1],
                'region': row[2]
            })

        return countries

    def get_inventories_for_country(self, country_name=None, country_id=None):
        """
        Get all inventory sources for a specific country.

        Parameters
        ----------
        country_name : str, optional
            Name of the country (case-insensitive partial match)
        country_id : int, optional
            Country ID (exact match)

        Returns
        -------
        list of dict
            List of dictionaries with keys: 'source_id', 'source_name', 'country', 'link', 'type'

        Raises
        ------
        ValueError
            If neither ``country_name`` nor ``country_id`` is given.

        Note
        ----
        Either country_name or country_id must be provided.

        Examples
        --------
        >>> [i["source_id"] for i in ZeroPM().get_inventories_for_country("Japan")]
        [8, 9, 10, 11]
        """
        if country_name is None and country_id is None:
            raise ValueError("Either country_name or country_id must be provided")

        if country_id is not None:
            self.cursor.execute("""
                SELECT DISTINCT s.source_id, s.source_name, c.country, s.link, s.type
                FROM sources s
                JOIN country_sources_index csi ON s.source_id = csi.source_id
                JOIN countries c ON csi.country_id = c.country_id
                WHERE c.country_id = ?
                ORDER BY s.source_name
            """, (country_id,))
        else:
            self.cursor.execute("""
                SELECT DISTINCT s.source_id, s.source_name, c.country, s.link, s.type
                FROM sources s
                JOIN country_sources_index csi ON s.source_id = csi.source_id
                JOIN countries c ON csi.country_id = c.country_id
                WHERE LOWER(c.country) LIKE LOWER(?)
                ORDER BY s.source_name
            """, (f'%{country_name}%',))

        inventories = []
        for row in self.cursor.fetchall():
            inventories.append({
                'source_id': row[0],
                'source_name': row[1],
                'country': row[2],
                'link': row[3],
                'type': row[4]
            })

        return inventories

    def count_chemicals_by_inventory(self, source_id):
        """
        Count the number of chemicals in a specific inventory.

        Parameters
        ----------
        source_id : int
            Source ID

        Returns
        -------
        int
            Number of unique CAS numbers in the inventory

        Examples
        --------
        >>> ZeroPM().count_chemicals_by_inventory(12)   # South Korea's
        21580
        """
        self.cursor.execute("""
            SELECT COUNT(DISTINCT aq.query)
            FROM api_ready_query aq
            JOIN inventory_summary issum ON aq.query_id = issum.query_id
            JOIN inventories inv ON issum.inventory_id = inv.inventory_id
            WHERE aq.type = 'CAS Registry Number' AND inv.source_id = ?
        """, (source_id,))

        return self.cursor.fetchone()[0]

    def count_chemicals_by_country(self, country_id):
        """
        Count the number of chemicals registered in a specific country.

        Parameters
        ----------
        country_id : int
            Country ID

        Returns
        -------
        int
            Number of unique CAS numbers in the country, over all its
            inventories

        Examples
        --------
        >>> ZeroPM().count_chemicals_by_country(1)      # Australia
        25183
        """
        self.cursor.execute("""
            SELECT COUNT(DISTINCT aq.query)
            FROM api_ready_query aq
            JOIN inventory_summary issum ON aq.query_id = issum.query_id
            JOIN inventories inv ON issum.inventory_id = inv.inventory_id
            JOIN sources s ON inv.source_id = s.source_id
            JOIN country_sources_index csi ON s.source_id = csi.source_id
            WHERE aq.type = 'CAS Registry Number' AND csi.country_id = ?
        """, (country_id,))

        return self.cursor.fetchone()[0]

    def count_chemicals_by_region(self, region_id):
        """
        Count the number of chemicals registered in a specific region.

        Parameters
        ----------
        region_id : int
            Region ID

        Returns
        -------
        int
            Number of unique CAS numbers in the region

        Examples
        --------
        >>> ZeroPM().count_chemicals_by_region(5)       # Oceania
        34175
        """
        self.cursor.execute("""
            SELECT COUNT(DISTINCT aq.query)
            FROM api_ready_query aq
            JOIN inventory_summary issum ON aq.query_id = issum.query_id
            JOIN inventories inv ON issum.inventory_id = inv.inventory_id
            JOIN sources s ON inv.source_id = s.source_id
            JOIN country_sources_index csi ON s.source_id = csi.source_id
            JOIN countries c ON csi.country_id = c.country_id
            JOIN region_country_index rci ON c.country_id = rci.country_id
            WHERE aq.type = 'CAS Registry Number' AND rci.region_id = ?
        """, (region_id,))

        return self.cursor.fetchone()[0]

    # ==================== ZeroPM Specific Methods (v0-0-4) ====================

    def get_zeropm_id(self, cas=None, inchi_id=None):
        """
        Get the zeropm_id for a chemical from CAS number or inchi_id.

        Parameters
        ----------
        cas : str, optional
            CAS Registry Number
        inchi_id : int, optional
            InChI identifier

        Returns
        -------
        int or None
            zeropm_id if found, None otherwise. A CAS number is resolved to
            its rank-1 structure first.

        Raises
        ------
        ValueError
            If neither ``cas`` nor ``inchi_id`` is given.

        Note
        ----
        Either cas or inchi_id must be provided.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> zpm.get_zeropm_id(cas="50-00-0"), zpm.get_zeropm_id(inchi_id=32227)
        (3224, 3224)
        """
        if cas is None and inchi_id is None:
            raise ValueError("Either cas or inchi_id must be provided")

        if inchi_id is None:
            inchi_id = self._inchi_id_from_cas(cas)
            if inchi_id is None:
                return None

        # Get zeropm_id from inchi_id
        self.cursor.execute("""
            SELECT zeropm_id
            FROM zeropm_chemicals
            WHERE inchi_id = ?
        """, (inchi_id,))
        result = self.cursor.fetchone()
        return result[0] if result else None

    def zeropm_id_to_inchi_id(self, zeropm_id):
        """
        Get the inchi_id for a zeropm_id — the reverse of :meth:`get_zeropm_id`.

        Parameters
        ----------
        zeropm_id : int
            ZeroPM identifier.

        Returns
        -------
        int or None
            The inchi_id, or None when the zeropm_id is not in the database.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> zpm.zeropm_id_to_inchi_id(1)
        6210
        >>> zpm.zeropm_id_to_inchi_id(3224)   # formaldehyde
        32227
        """
        self.cursor.execute("""
            SELECT inchi_id
            FROM zeropm_chemicals
            WHERE zeropm_id = ?
        """, (zeropm_id,))
        result = self.cursor.fetchone()
        return result[0] if result else None

    def _inchi_id_from_cas(self, cas):
        """
        Resolve a CAS number to its first inchi_id.

        Parameters
        ----------
        cas : str
            CAS Registry Number.

        Returns
        -------
        int or None
            The first inchi_id for the CAS, or None when it is not found.
        """
        query_id = self.query_cas(cas)
        if query_id is None:
            return None
        inchi_ids, _ = self.get_inchi_id(query_id)
        return inchi_ids[0] if inchi_ids else None


    def get_pm_probabilities(self, cas=None, inchi_id=None, zeropm_id=None):
        """
        Get P/M (Persistent/Mobile) probability data for a chemical.

        Parameters
        ----------
        cas : str, optional
            CAS Registry Number
        inchi_id : int, optional
            InChI identifier
        zeropm_id : int, optional
            ZeroPM identifier

        Returns
        -------
        dict or None
            Dictionary with probability data:
            - probability_of_not_p: Probability of NOT persistent
            - probability_of_p_or_vp: Probability of persistent OR very persistent
            - probability_of_p: Probability of persistent
            - probability_of_vp: Probability of very persistent
            - probability_of_not_m: Probability of NOT mobile
            - probability_of_m_or_vm: Probability of mobile OR very mobile
            - probability_of_m: Probability of mobile
            - probability_of_vm: Probability of very mobile
            - n: Sample size
            Returns None if not found, or if ZeroPM assessed the chemical
            but published no probabilities for it --- formaldehyde is one.

        Raises
        ------
        ValueError
            If none of ``cas``, ``inchi_id`` or ``zeropm_id`` is provided.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> probs = zpm.get_pm_probabilities(inchi_id=6210)
        >>> round(probs["probability_of_p"], 3)
        0.4
        >>> round(zpm.get_pm_probabilities(cas="64-17-5")["probability_of_m"], 3)
        0.682
        >>> zpm.get_pm_probabilities(cas="50-00-0") is None
        True

        Note
        ----
        ``pm_probabilities`` is keyed on ``inchi_id``, so a ``zeropm_id`` is
        translated first via :meth:`zeropm_id_to_inchi_id`.
        """
        if cas is None and inchi_id is None and zeropm_id is None:
            raise ValueError("One of cas, inchi_id or zeropm_id must be provided")

        if inchi_id is None:
            if zeropm_id is not None:
                inchi_id = self.zeropm_id_to_inchi_id(zeropm_id)
            else:
                inchi_id = self._inchi_id_from_cas(cas)
            if inchi_id is None:
                return None

        self.cursor.execute("""
            SELECT probability_of_not_p, probability_of_p_or_vp, probability_of_p, probability_of_vp,
                   probability_of_not_m, probability_of_m_or_vm, probability_of_m, probability_of_vm, n
            FROM pm_probabilities
            WHERE inchi_id = ?
        """, (inchi_id,))
        result = self.cursor.fetchone()

        if not result:
            return None

        return {
            'probability_of_not_p': result[0],
            'probability_of_p_or_vp': result[1],
            'probability_of_p': result[2],
            'probability_of_vp': result[3],
            'probability_of_not_m': result[4],
            'probability_of_m_or_vm': result[5],
            'probability_of_m': result[6],
            'probability_of_vm': result[7],
            'n': result[8]
        }

    def is_in_zeropm(self, cas=None, inchi_id=None):
        """
        Check if a chemical is in the ZeroPM database.

        Parameters
        ----------
        cas : str, optional
            CAS Registry Number
        inchi_id : int, optional
            InChI identifier

        Returns
        -------
        bool
            True if the chemical has a ``zeropm_id`` (ZeroPM assessed it),
            False otherwise --- including a CAS number that is in the
            inventories but was not assessed

        Raises
        ------
        ValueError
            If neither ``cas`` nor ``inchi_id`` is given.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> zpm.is_in_zeropm(cas="50-00-0"), zpm.is_in_zeropm(cas="0-00-0")
        (True, False)
        """
        return self.get_zeropm_id(cas=cas, inchi_id=inchi_id) is not None

    def is_multicomponent(self, inchi_id):
        """
        Check if a substance is a multi-component substance.

        Parameters
        ----------
        inchi_id : int
            InChI identifier

        Returns
        -------
        bool
            True if substance is multi-component, False otherwise

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> zpm.is_multicomponent(5), zpm.is_multicomponent(32227)
        (True, False)
        """
        self.cursor.execute("""
            SELECT mc_id
            FROM multi_components
            WHERE inchi_id = ?
        """, (inchi_id,))
        return self.cursor.fetchone() is not None

    def get_multicomponent_id(self, inchi_id):
        """
        Get the multi-component ID for a substance.

        Parameters
        ----------
        inchi_id : int
            InChI identifier

        Returns
        -------
        int or None
            mc_id if found, None otherwise

        Examples
        --------
        >>> ZeroPM().get_multicomponent_id(5)
        1
        """
        self.cursor.execute("""
            SELECT mc_id
            FROM multi_components
            WHERE inchi_id = ?
        """, (inchi_id,))
        result = self.cursor.fetchone()
        return result[0] if result else None

    def get_components(self, mc_id):
        """
        Get all components of a multi-component substance.

        Parameters
        ----------
        mc_id : int
            Multi-component identifier

        Returns
        -------
        list of dict
            List of component information with keys:
            - component_id: Component identifier
            - component_frequency: How often the component appears
            - inchi_id: InChI identifier of the component
            - inchi: InChI string of the component
            - inchikey: InChIKey of the component
            Most frequent first. Empty for an unknown ``mc_id``.

        Examples
        --------
        >>> [c["inchi"] for c in ZeroPM().get_components(1)]
        ['InChI=1S/ClH/h1H/p-1', 'InChI=1S/C8H10N3/c1-11(2)8-5-3-7(10-9)4-6-8/h3-6H,1-2H3/q+1']
        """
        self.cursor.execute("""
            SELECT ci.component_id, ci.component_frequency, c.inchi_id, s.inchi, s.inchikey
            FROM component_index ci
            JOIN components c ON ci.component_id = c.component_id
            JOIN substances s ON c.inchi_id = s.inchi_id
            WHERE ci.mc_id = ?
            ORDER BY ci.component_frequency DESC
        """, (mc_id,))

        components = []
        for row in self.cursor.fetchall():
            components.append({
                'component_id': row[0],
                'component_frequency': row[1],
                'inchi_id': row[2],
                'inchi': row[3],
                'inchikey': row[4]
            })

        return components

    def get_multicomponent_info(self, cas=None, inchi_id=None):
        """
        Get complete multi-component information for a substance.

        Parameters
        ----------
        cas : str, optional
            CAS Registry Number
        inchi_id : int, optional
            InChI identifier

        Returns
        -------
        dict or None
            Dictionary with:
            - mc_id: Multi-component identifier
            - inchi_id: InChI identifier of the multi-component
            - inchi: InChI of the multi-component
            - inchikey: InChIKey of the multi-component
            - components: List of component dictionaries, as
              :meth:`get_components` returns them
            Returns None if not a multi-component substance. A CAS number is
            resolved to its rank-1 structure first.

        Raises
        ------
        ValueError
            If neither ``cas`` nor ``inchi_id`` is given.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> info = zpm.get_multicomponent_info(inchi_id=5)
        >>> info["mc_id"], info["inchikey"], len(info["components"])
        (1, 'CCIAVEMREXZXAK-UHFFFAOYSA-M', 2)
        >>> zpm.get_multicomponent_info(cas="50-00-0") is None
        True
        """
        if inchi_id is None:
            if cas is None:
                raise ValueError("Either cas or inchi_id must be provided")
            query_id = self.query_cas(cas)
            if query_id is None:
                return None
            inchi_ids, _ = self.get_inchi_id(query_id)
            if not inchi_ids:
                return None
            inchi_id = inchi_ids[0]

        # Check if it's a multi-component
        mc_id = self.get_multicomponent_id(inchi_id)
        if mc_id is None:
            return None

        # Get multi-component info
        self.cursor.execute("""
            SELECT mc.inchi_id, s.inchi, s.inchikey
            FROM multi_components mc
            JOIN substances s ON mc.inchi_id = s.inchi_id
            WHERE mc.mc_id = ?
        """, (mc_id,))
        result = self.cursor.fetchone()

        if not result:
            return None

        # Get components
        components = self.get_components(mc_id)

        return {
            'mc_id': mc_id,
            'inchi_id': result[0],
            'inchi': result[1],
            'inchikey': result[2],
            'components': components
        }

    def is_in_cleanventory(self, cas=None, inchi_id=None):
        """
        Check if a chemical is in the Cleanventory database.

        Parameters
        ----------
        cas : str, optional
            CAS Registry Number
        inchi_id : int, optional
            InChI identifier

        Returns
        -------
        bool
            True if chemical is in Cleanventory, False otherwise. A CAS
            number is resolved to its rank-1 structure first.

        Raises
        ------
        ValueError
            If neither ``cas`` nor ``inchi_id`` is given.

        Examples
        --------
        >>> zpm = ZeroPM()
        >>> zpm.is_in_cleanventory(cas="50-00-0"), zpm.is_in_cleanventory(cas="0-00-0")
        (True, False)
        """
        if inchi_id is None:
            if cas is None:
                raise ValueError("Either cas or inchi_id must be provided")
            query_id = self.query_cas(cas)
            if query_id is None:
                return False
            inchi_ids, _ = self.get_inchi_id(query_id)
            if not inchi_ids:
                return False
            inchi_id = inchi_ids[0]

        self.cursor.execute("""
            SELECT cleanventory_id
            FROM cleanventory_chemicals
            WHERE inchi_id = ?
        """, (inchi_id,))
        return self.cursor.fetchone() is not None

    def get_consensus_score(self, cas=None, inchi_id=None):
        """
        Get consensus scoring information for a chemical.

        Parameters
        ----------
        cas : str, optional
            CAS Registry Number
        inchi_id : int, optional
            InChI identifier

        Returns
        -------
        list of dict or None
            List of consensus scores from different inventories, each with:
            - inventory_id: Inventory identifier
            - consensus_score: Consensus score value
            - consensus_count: Count of consensus
            Returns None if not found. The values are as stored: in v0.0.4
            ``consensus_score`` is a string of a small integer and
            ``consensus_count`` a fraction between 0 and 1.

        Raises
        ------
        ValueError
            If neither ``cas`` nor ``inchi_id`` is given.

        Examples
        --------
        >>> scores = ZeroPM().get_consensus_score(cas="50-00-0")
        >>> scores[0]
        {'inventory_id': 692, 'consensus_score': '2', 'consensus_count': 0.198675496688742}
        """
        if inchi_id is None:
            if cas is None:
                raise ValueError("Either cas or inchi_id must be provided")
            query_id = self.query_cas(cas)
            if query_id is None:
                return None
            inchi_ids, _ = self.get_inchi_id(query_id)
            if not inchi_ids:
                return None
            inchi_id = inchi_ids[0]

        self.cursor.execute("""
            SELECT inventory_id, consensus_score, consensus_count
            FROM consensus_index
            WHERE inchi_id = ?
        """, (inchi_id,))

        results = self.cursor.fetchall()
        if not results:
            return None

        consensus_data = []
        for row in results:
            consensus_data.append({
                'inventory_id': row[0],
                'consensus_score': row[1],
                'consensus_count': row[2]
            })

        return consensus_data

    def get_all_zeropm_chemicals(self, limit=None, include_pm_probs=False):
        """
        Get all chemicals in the ZeroPM database.

        Parameters
        ----------
        limit : int, optional
            Maximum number of results to return
        include_pm_probs : bool, optional
            If True, include P/M probability data (default: False)

        Returns
        -------
        pandas.DataFrame
            DataFrame with zeropm_id, inchi_id, inchi, inchikey
            If include_pm_probs=True, also includes all probability columns,
            NaN where none were published

        Examples
        --------
        >>> df = ZeroPM().get_all_zeropm_chemicals(limit=2, include_pm_probs=True)
        >>> df[["zeropm_id", "inchi_id", "probability_of_p", "n"]].round(3)
           zeropm_id  inchi_id  probability_of_p  n
        0          1      6210             0.400  1
        1          2    101901             0.438  1
        """
        if include_pm_probs:
            query = """
                SELECT zc.zeropm_id, zc.inchi_id, s.inchi, s.inchikey,
                       pm.probability_of_not_p, pm.probability_of_p_or_vp,
                       pm.probability_of_p, pm.probability_of_vp,
                       pm.probability_of_not_m, pm.probability_of_m_or_vm,
                       pm.probability_of_m, pm.probability_of_vm, pm.n
                FROM zeropm_chemicals zc
                JOIN substances s ON zc.inchi_id = s.inchi_id
                LEFT JOIN pm_probabilities pm ON zc.inchi_id = pm.inchi_id
            """
            columns = ['zeropm_id', 'inchi_id', 'inchi', 'inchikey',
                      'probability_of_not_p', 'probability_of_p_or_vp',
                      'probability_of_p', 'probability_of_vp',
                      'probability_of_not_m', 'probability_of_m_or_vm',
                      'probability_of_m', 'probability_of_vm', 'n']
        else:
            query = """
                SELECT zc.zeropm_id, zc.inchi_id, s.inchi, s.inchikey
                FROM zeropm_chemicals zc
                JOIN substances s ON zc.inchi_id = s.inchi_id
            """
            columns = ['zeropm_id', 'inchi_id', 'inchi', 'inchikey']

        if limit:
            query += f" LIMIT {limit}"

        self.cursor.execute(query)
        results = self.cursor.fetchall()

        return pd.DataFrame(results, columns=columns)

    def get_all_multicomponent_substances(self, limit=None):
        """
        Get all multi-component substances.

        Parameters
        ----------
        limit : int, optional
            Maximum number of results to return

        Returns
        -------
        pandas.DataFrame
            DataFrame with mc_id, inchi_id, inchi, inchikey, component_count

        Examples
        --------
        >>> ZeroPM().get_all_multicomponent_substances(limit=2)[["mc_id", "inchi_id", "component_count"]]
           mc_id  inchi_id  component_count
        0      1         5                2
        1      2         6                2
        """
        query = """
            SELECT mc.mc_id, mc.inchi_id, s.inchi, s.inchikey,
                   COUNT(ci.component_id) as component_count
            FROM multi_components mc
            JOIN substances s ON mc.inchi_id = s.inchi_id
            LEFT JOIN component_index ci ON mc.mc_id = ci.mc_id
            GROUP BY mc.mc_id, mc.inchi_id, s.inchi, s.inchikey
        """

        if limit:
            query += f" LIMIT {limit}"

        self.cursor.execute(query)
        results = self.cursor.fetchall()

        return pd.DataFrame(results, columns=['mc_id', 'inchi_id', 'inchi', 'inchikey', 'component_count'])

    def batch_get_pm_probabilities(self, cas_list=None, inchi_id_list=None):
        """
        Get P/M probabilities for multiple chemicals at once.

        Parameters
        ----------
        cas_list : list of str, optional
            List of CAS Registry Numbers
        inchi_id_list : list of int, optional
            List of InChI identifiers

        Returns
        -------
        pandas.DataFrame
            DataFrame with columns for identifiers and all probability values:
            one row per chemical ZeroPM assessed, with ``cas`` first when
            ``cas_list`` was given. Chemicals it did not assess, and CAS
            numbers not in the database, have no row; one assessed without
            published probabilities has NaN. Empty when nothing is found.

        Examples
        --------
        >>> df = ZeroPM().batch_get_pm_probabilities(cas_list=["50-00-0", "64-17-5", "0-00-0"])
        >>> df[["cas", "probability_of_p", "probability_of_m"]].round(3)
               cas  probability_of_p  probability_of_m
        0  50-00-0               NaN               NaN
        1  64-17-5             0.307             0.682
        """
        if cas_list is not None:
            # Convert CAS to inchi_ids
            inchi_id_list = []
            cas_to_inchi_id = {}
            for cas in cas_list:
                query_id = self.query_cas(cas)
                if query_id:
                    inchi_ids, _ = self.get_inchi_id(query_id)
                    if inchi_ids:
                        inchi_id = inchi_ids[0]
                        inchi_id_list.append(inchi_id)
                        cas_to_inchi_id[inchi_id] = cas

        if not inchi_id_list:
            return pd.DataFrame()

        # Query all at once
        placeholders = ','.join('?' * len(inchi_id_list))
        query = f"""
            SELECT zc.inchi_id, s.inchi, s.inchikey,
                   pm.probability_of_not_p, pm.probability_of_p_or_vp,
                   pm.probability_of_p, pm.probability_of_vp,
                   pm.probability_of_not_m, pm.probability_of_m_or_vm,
                   pm.probability_of_m, pm.probability_of_vm, pm.n
            FROM zeropm_chemicals zc
            JOIN substances s ON zc.inchi_id = s.inchi_id
            LEFT JOIN pm_probabilities pm ON zc.inchi_id = pm.inchi_id
            WHERE zc.inchi_id IN ({placeholders})
        """

        self.cursor.execute(query, inchi_id_list)
        results = self.cursor.fetchall()

        df = pd.DataFrame(results, columns=[
            'inchi_id', 'inchi', 'inchikey',
            'probability_of_not_p', 'probability_of_p_or_vp',
            'probability_of_p', 'probability_of_vp',
            'probability_of_not_m', 'probability_of_m_or_vm',
            'probability_of_m', 'probability_of_vm', 'n'
        ])

        # Add CAS if available
        if cas_list is not None:
            df['cas'] = df['inchi_id'].map(cas_to_inchi_id)
            # Reorder columns to put cas first
            cols = ['cas'] + [col for col in df.columns if col != 'cas']
            df = df[cols]

        return df
