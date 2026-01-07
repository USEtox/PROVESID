import sqlite3
import os
from rdkit import Chem
import logging
from rapidfuzz import process, fuzz, utils
import requests
from tqdm import tqdm
import pandas as pd

# Try relative import first, fall back to direct import for testing
from .utils import data_path as ZeroPM_path


class ZeroPM:
    """
    Class to extract data from the ZeroPM SQLite database using SQL queries.
    This class provides the same functionality as ZeroPM but uses SQL instead of pandas.
    SMILES are generated on-the-fly from InChI using RDKit.
    
    The database file will be automatically downloaded if not found locally.
    """
    
    # Default download URL for the ZeroPM database
    DEFAULT_DB_URL = "https://github.com/ZeroPM-H2020/global-chemical-inventory-database/raw/refs/heads/main/zeropm-v0-0-3.sqlite"
    
    def __init__(self, db_name='zeropm-v0-0-3.sqlite', auto_download=True, db_url=None):
        """
        Initialize connection to the ZeroPM SQLite database.
        
        Parameters
        ----------
        db_name : str, optional
            Name of the SQLite database file (default: 'zeropm-v0-0-3.sqlite')
        auto_download : bool, optional
            If True, automatically download the database if not found (default: True)
        db_url : str, optional
            Custom URL to download the database from. If None, uses the default GitHub URL.
        """
        self.path = ZeroPM_path()
        self.db_path = os.path.join(self.path, db_name)
        self.db_url = db_url or self.DEFAULT_DB_URL
        
        # Check if database exists, download if needed
        if not os.path.exists(self.db_path):
            if auto_download:
                logging.info(f"Database not found at: {self.db_path}")
                logging.info("Downloading database automatically...")
                self.download_database(url=self.db_url, force=False)
            else:
                raise FileNotFoundError(
                    f"Database not found at: {self.db_path}\n"
                    f"Please run ZeroPM.download_database() or set auto_download=True"
                )
        
        # Create connection (will be reused for queries)
        self.conn = sqlite3.connect(self.db_path)
        self.cursor = self.conn.cursor()
        
        # Cache chemical names for fuzzy matching (lazy loading)
        self._chemical_names_cache = None
    
    def __del__(self):
        """Close database connection when object is destroyed."""
        if hasattr(self, 'conn'):
            self.conn.close()
    
    def download_database(self, url=None, force=False):
        """
        Download the ZeroPM SQLite database from a remote URL.
        
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
        requests.exceptions.RequestException
            If the download fails
            
        Example
        -------
        >>> zpm = ZeroPM(auto_download=False)  # Don't auto-download
        >>> zpm.download_database()  # Manually trigger download
        """
        download_url = url or self.db_url
        
        # Check if database already exists
        if os.path.exists(self.db_path) and not force:
            raise FileExistsError(
                f"Database already exists at: {self.db_path}\n"
                f"Use force=True to overwrite"
            )
        
        # Create data directory if it doesn't exist
        os.makedirs(self.path, exist_ok=True)
        
        logging.info(f"Downloading ZeroPM database from: {download_url}")
        logging.info(f"Destination: {self.db_path}")
        
        try:
            # Stream the download with progress bar
            response = requests.get(download_url, stream=True, timeout=30)
            response.raise_for_status()
            
            # Get total file size
            total_size = int(response.headers.get('content-length', 0))
            
            # Download with progress bar
            temp_path = self.db_path + '.tmp'
            with open(temp_path, 'wb') as f:
                if total_size > 0:
                    # Show progress bar
                    with tqdm(total=total_size, unit='B', unit_scale=True, 
                             desc="Downloading ZeroPM database") as pbar:
                        for chunk in response.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                                pbar.update(len(chunk))
                else:
                    # No content-length header, download without progress
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
            
            # Move temp file to final location
            if os.path.exists(self.db_path):
                os.remove(self.db_path)
            os.rename(temp_path, self.db_path)
            
            logging.info(f"✓ Database downloaded successfully to: {self.db_path}")
            
            # Verify the database is valid
            try:
                test_conn = sqlite3.connect(self.db_path)
                test_cursor = test_conn.cursor()
                test_cursor.execute("SELECT name FROM sqlite_master WHERE type='table' LIMIT 1")
                test_cursor.fetchone()
                test_conn.close()
                logging.info("✓ Database verified successfully")
            except sqlite3.Error as e:
                os.remove(self.db_path)
                raise RuntimeError(f"Downloaded database is corrupted: {e}")
            
            return self.db_path
            
        except requests.exceptions.RequestException as e:
            # Clean up temp file if it exists
            if os.path.exists(temp_path):
                os.remove(temp_path)
            logging.error(f"Failed to download database: {e}")
            raise
    
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
            Exact chemical name
            
        Returns
        -------
        int or None
            query_id if found, None otherwise
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
    
    def get_inchi_id(self, query_id):
        """
        Returns all the inchi_id and ranks of a query with a given query_id.
        
        Parameters
        ----------
        query_id : int
            Query identifier
            
        Returns
        -------
        tuple of (numpy.ndarray, numpy.ndarray)
            (inchi_ids, ranks) sorted by rank, with duplicates removed
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
            List of chemical names (excluding the CAS number itself)
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
            SMILES string, or None if not found
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
            CAS number, list of CAS numbers, or None if not found
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
            CAS number, list of CAS numbers, or None if not found
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
            CAS number, list of CAS numbers, or None if not found
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
    
    def get_id_table_from_cas(self, cas):
        """
        Returns a pandas DataFrame containing all identifiers for a given CAS number.
        
        This method retrieves all query_ids associated with the CAS number, then for each query_id,
        it retrieves all associated inchi_ids and their corresponding InChI and InChIKey values.
        Synonyms (chemical names) are also included.
        
        Parameters
        ----------
        cas : str
            CAS Registry Number
            
        Returns
        -------
        pandas.DataFrame
            DataFrame with columns: 'cas', 'query_id', 'inchi_id', 'rank', 'inchi', 'inchikey', 'synonyms'
            Returns None if the CAS number is not found in the database.
            
        Examples
        --------
        >>> zpm = ZeroPM()
        >>> df = zpm.get_id_table_from_cas("50-00-0")
        >>> print(df)
        """
        # Get all query_ids for this CAS (using fetchall in case there are multiple)
        self.cursor.execute("""
            SELECT query_id 
            FROM api_ready_query 
            WHERE query = ? AND type = 'CAS Registry Number'
        """, (cas,))
        query_ids = [row[0] for row in self.cursor.fetchall()]
        
        if not query_ids:
            logging.warning(f"CAS number {cas} not found in database")
            return None
        
        # Get synonyms for this CAS
        synonyms = self.get_names(cas)
        synonyms_str = "; ".join(synonyms) if synonyms else ""
        
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
                    'synonyms': synonyms_str
                })
            else:
                # For each inchi_id, get the inchi and inchikey
                for inchi_id, rank in zip(inchi_ids, ranks):
                    inchi, inchikey = self.get_inchi(inchi_id)
                    rows.append({
                        'cas': cas,
                        'query_id': query_id,
                        'inchi_id': inchi_id,
                        'rank': rank,
                        'inchi': inchi,
                        'inchikey': inchikey,
                        'synonyms': synonyms_str
                    })
        
        # Create DataFrame
        df = pd.DataFrame(rows)
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
            Combined DataFrame with columns: 'cas', 'query_id', 'inchi_id', 'rank', 'inchi', 'inchikey', 'synonyms'
            Returns an empty DataFrame if no CAS numbers are found in the database.
            
        Examples
        --------
        >>> zpm = ZeroPM()
        >>> cas_numbers = ["50-00-0", "50-78-2", "64-17-5"]  # formaldehyde, aspirin, ethanol
        >>> df = zpm.batch_get_id_table_from_cas(cas_numbers)
        >>> print(df)
        >>> # Group by CAS to see counts
        >>> print(df.groupby('cas').size())
        """
        if not cas_list:
            logging.warning("Empty CAS list provided")
            return pd.DataFrame(columns=['cas', 'query_id', 'inchi_id', 'rank', 'inchi', 'inchikey', 'synonyms'])
        
        # Collect DataFrames for each CAS
        dataframes = []
        for cas in cas_list:
            df = self.get_id_table_from_cas(cas)
            if df is not None:
                dataframes.append(df)
        
        # Combine all DataFrames
        if not dataframes:
            logging.warning("None of the provided CAS numbers were found in the database")
            return pd.DataFrame(columns=['cas', 'query_id', 'inchi_id', 'rank', 'inchi', 'inchikey', 'synonyms'])
        
        # Concatenate all DataFrames and reset index
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
            Dictionary with index names as keys and status ('created', 'exists', 'error') as values
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
            Dictionary mapping CAS numbers to lists of names
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
            Dictionary mapping InChIKeys to CAS numbers (or None if not found)
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
        Search for chemical names using regular expressions.
        
        Parameters
        ----------
        pattern : str
            Regular expression pattern (SQLite REGEXP syntax)
        case_sensitive : bool, optional
            Whether the search is case-sensitive (default: False)
        limit : int, optional
            Maximum number of results to return (default: 100)
            
        Returns
        -------
        list of tuple
            List of (query_id, name) tuples matching the pattern
            
        Note
        ----
        SQLite's REGEXP requires the pattern to match the entire string unless
        wildcards are used. Use '.*pattern.*' for substring matching.
        """
        # SQLite REGEXP is case-sensitive by default
        # For case-insensitive, we use LIKE with wildcards
        if case_sensitive:
            # Enable REGEXP (requires loading regexp extension or using LIKE alternative)
            # Using LIKE as fallback with % wildcards
            pattern = pattern.replace('.*', '%').replace('.', '_')
            self.cursor.execute(f"""
                SELECT query_id, query 
                FROM api_ready_query 
                WHERE type = 'chemical name' AND query LIKE ?
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
        This method retrieves all SMILES from the database and performs
        substructure matching using RDKit.
        
        Parameters
        ----------
        smarts_pattern : str
            SMARTS pattern for substructure search
        max_results : int, optional
            Maximum number of results to return (default: 100)
            
        Returns
        -------
        list of dict
            List of dictionaries with keys: 'cas', 'inchi', 'inchikey', 'smiles'
            
        Warning
        -------
        This method can be slow for large databases as it needs to convert
        all InChI to SMILES and perform substructure matching.
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
            Output CSV filename
        columns : list of str, optional
            Column names for the CSV header
            
        Returns
        -------
        str
            Path to the created CSV file
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
            True if view was created successfully, False otherwise
            
        Example
        -------
        >>> sql = '''
        ...     SELECT aq.query AS cas, s.inchi, s.inchikey
        ...     FROM api_ready_query aq
        ...     JOIN api_results ar ON aq.query_id = ar.query_id
        ...     JOIN substances s ON ar.inchi_id = s.inchi_id
        ...     WHERE aq.type = 'CAS Registry Number' AND ar.rank = 1
        ... '''
        >>> zpm.create_view('cas_to_inchi', sql)
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
            Path to the created CSV file
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
            Dictionary with statistics about each table
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
            List of dictionaries with keys: 'source_id', 'source_name', 'country_scope', 'link', 'type'
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
            List of dictionaries with keys: 'country_id', 'country'
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
            List of dictionaries with keys: 'region_id', 'region'
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
            List of chemicals with keys: 'cas', 'query_id', 'inchi_id', 'source_name'
            
        Note
        ----
        Either source_name or source_id must be provided.
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
            List of chemicals with keys: 'cas', 'query_id', 'inchi_id', 'country', 'source_name'
            
        Note
        ----
        Either country_name or country_id must be provided.
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
            List of chemicals with keys: 'cas', 'query_id', 'inchi_id', 'region', 'country', 'source_name'
            
        Note
        ----
        Either region_name or region_id must be provided.
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
            
        Note
        ----
        Either region_name or region_id must be provided.
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
            
        Note
        ----
        Either country_name or country_id must be provided.
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
            Number of unique CAS numbers in the country
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
