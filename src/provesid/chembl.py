"""
ChEMBL Database Interface

This module provides access to the ChEMBL SQLite database for querying chemical compounds,
structures, and properties. ChEMBL is a manually curated database of bioactive molecules
with drug-like properties maintained by EMBL-EBI.

Database tables accessed:
- molecule_dictionary: Primary compound information (ChEMBL ID, names, max_phase, drug 
                       classifications, approval status, administration routes)
- molecule_hierarchy: Parent-salt-metabolite relationships for compounds and pro-drugs
- compound_structures: Chemical structures (SMILES, InChI, InChIKey, molfile)
- compound_properties: Physicochemical properties (MW, ALogP, HBA, HBD, PSA, etc.)
- molecule_synonyms: Alternative names and synonyms
- chembl_id_lookup: ChEMBL ID to internal ID mappings
- pesticide_classification: Pesticide mechanism classifications (FRAC, HRAC, IRAC)
- pesticide_class_mapping: Links compounds to pesticide classifications

For detailed table schema information, see src/provesid/data/schema_documentation.txt
"""

import os
import sqlite3
import tarfile
import logging
from typing import Optional, Dict, List, Any
import requests
from tqdm import tqdm
from .utils import data_path


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
        Name of the SQLite database file (default: 'chembl_36.db')
    auto_download : bool, optional
        If True, automatically download database if not found (default: True)
    db_url : str, optional
        Custom URL for database download (default: ChEMBL FTP URL)
    
    Attributes
    ----------
    path : str
        Path to the data directory
    db_path : str
        Full path to the SQLite database file
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
    The database file is approximately 5GB. Initial setup will download and extract
    the database from the EMBL-EBI FTP server (~1.5GB compressed).
    """
    
    DEFAULT_DB_URL = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_36_sqlite.tar.gz"
    
    def __init__(self, db_name: str = 'chembl_36.db', auto_download: bool = True, db_url: Optional[str] = None):
        """
        Initialize ChEMBL database interface.
        
        Parameters
        ----------
        db_name : str, optional
            Database filename (default: 'chembl_36.db')
        auto_download : bool, optional
            Auto-download if database missing (default: True)
        db_url : str, optional
            Custom download URL (default: ChEMBL FTP)
        
        Raises
        ------
        FileNotFoundError
            If database not found and auto_download is False
        ChEMBLError
            If database connection or validation fails
        """
        self.path = data_path()
        self.db_path = os.path.join(self.path, db_name)
        self.db_url = db_url or self.DEFAULT_DB_URL
        self.logger = logging.getLogger(__name__)
        
        # Ensure data directory exists
        os.makedirs(self.path, exist_ok=True)
        
        # Check if database exists
        if not os.path.exists(self.db_path):
            if auto_download:
                self.logger.info(f"Database not found at {self.db_path}, downloading...")
                self.download_database(url=self.db_url, force=False)
            else:
                raise FileNotFoundError(
                    f"ChEMBL database not found at {self.db_path}. "
                    f"Set auto_download=True or manually download from {self.DEFAULT_DB_URL}"
                )
        
        # Connect to database
        try:
            self.conn = sqlite3.connect(self.db_path)
            self.conn.row_factory = sqlite3.Row  # Enable column access by name
            self.cursor = self.conn.cursor()
            self.logger.info(f"Connected to ChEMBL database at {self.db_path}")
        except sqlite3.Error as e:
            raise ChEMBLError(f"Failed to connect to database: {str(e)}")
    
    def __del__(self):
        """Close database connection when object is destroyed"""
        if hasattr(self, 'conn') and self.conn:
            try:
                self.conn.close()
                self.logger.debug("Closed ChEMBL database connection")
            except Exception as e:
                self.logger.warning(f"Error closing database connection: {str(e)}")
    
    def download_database(self, url: Optional[str] = None, force: bool = False):
        """
        Download and extract ChEMBL SQLite database from EMBL-EBI FTP.
        
        Downloads the compressed tar.gz archive (~1.5GB), extracts the SQLite database
        (~5GB), and validates the database integrity by querying the molecule_dictionary table.
        
        Parameters
        ----------
        url : str, optional
            Download URL (default: DEFAULT_DB_URL)
        force : bool, optional
            If True, re-download even if database exists (default: False)
        
        Raises
        ------
        ChEMBLError
            If download, extraction, or validation fails
        
        Examples
        --------
        >>> chembl = CheMBL(auto_download=False)  # Will raise FileNotFoundError
        >>> chembl.download_database(force=True)  # Explicit download
        """
        url = url or self.db_url
        
        # Check if already exists
        if os.path.exists(self.db_path) and not force:
            self.logger.info(f"Database already exists at {self.db_path}")
            return
        
        self.logger.info(f"Downloading ChEMBL database from {url}")
        
        # Temporary files
        tar_gz_path = self.db_path + ".tar.gz.tmp"
        
        try:
            # Download compressed archive with progress bar
            response = requests.get(url, stream=True, timeout=30)
            response.raise_for_status()
            
            total_size = int(response.headers.get('content-length', 0))
            
            with open(tar_gz_path, 'wb') as f:
                with tqdm(total=total_size, unit='B', unit_scale=True, 
                         desc="Downloading ChEMBL database") as pbar:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))
            
            self.logger.info(f"Download complete. Extracting database...")
            
            # Extract tar.gz archive
            # ChEMBL archive structure: chembl_36_sqlite.tar.gz -> chembl_36/chembl_36.db
            with tarfile.open(tar_gz_path, 'r:gz') as tar:
                # Find the .db file in the archive
                db_members = [m for m in tar.getmembers() if m.name.endswith('.db')]
                
                if not db_members:
                    raise ChEMBLError("No .db file found in the tar.gz archive")
                
                db_member = db_members[0]
                self.logger.info(f"Extracting {db_member.name}...")
                
                # Extract with progress bar
                with tqdm(total=db_member.size, unit='B', unit_scale=True,
                         desc="Extracting database") as pbar:
                    # Extract to temp location first
                    tar.extract(db_member, path=self.path)
                    pbar.update(db_member.size)
                
                # Move extracted file to final location
                extracted_path = os.path.join(self.path, db_member.name)
                if extracted_path != self.db_path:
                    os.rename(extracted_path, self.db_path)
                    # Clean up extracted directory if it exists
                    extracted_dir = os.path.join(self.path, os.path.dirname(db_member.name))
                    if os.path.exists(extracted_dir) and os.path.isdir(extracted_dir):
                        os.rmdir(extracted_dir)
            
            # Validate database integrity
            self.logger.info("Validating database integrity...")
            try:
                test_conn = sqlite3.connect(self.db_path)
                test_cursor = test_conn.cursor()
                test_cursor.execute("SELECT COUNT(*) FROM molecule_dictionary")
                count = test_cursor.fetchone()[0]
                test_conn.close()
                self.logger.info(f"Database validated successfully. Contains {count:,} compounds.")
            except sqlite3.Error as e:
                os.remove(self.db_path)
                raise ChEMBLError(f"Database validation failed: {str(e)}")
            
            # Clean up temporary files
            if os.path.exists(tar_gz_path):
                os.remove(tar_gz_path)
            
            self.logger.info("ChEMBL database download and setup complete")
            
        except requests.RequestException as e:
            if os.path.exists(tar_gz_path):
                os.remove(tar_gz_path)
            raise ChEMBLError(f"Download failed: {str(e)}")
        except Exception as e:
            # Clean up on any error
            if os.path.exists(tar_gz_path):
                os.remove(tar_gz_path)
            if os.path.exists(self.db_path):
                os.remove(self.db_path)
            raise ChEMBLError(f"Database setup failed: {str(e)}")
    
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
    
    def search_by_name(self, name: str, limit: int = 100) -> List[Dict[str, Any]]:
        """
        Search for compounds by name (case-insensitive partial match).
        
        Searches both preferred names and synonyms.
        
        Parameters
        ----------
        name : str
            Compound name or partial name to search
        limit : int, optional
            Maximum number of results (default: 100)
        
        Returns
        -------
        list of dict
            List of matching compounds with structure information
        
        Examples
        --------
        >>> chembl = CheMBL()
        >>> results = chembl.search_by_name('aspirin')
        >>> print(len(results))
        1
        >>> print(results[0]['chembl_id'])
        'CHEMBL25'
        """
        try:
            # Search in preferred names and synonyms
            query = """
            SELECT DISTINCT md.molregno
            FROM molecule_dictionary md
            LEFT JOIN molecule_synonyms ms ON md.molregno = ms.molregno
            WHERE LOWER(md.pref_name) LIKE LOWER(?) 
               OR LOWER(ms.synonyms) LIKE LOWER(?)
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
            - molregno, chembl_id, pref_name, max_phase
            - canonical_smiles, standard_inchi, standard_inchi_key
            - molfile (if available)
            - synonyms: list of alternative names
            Returns None if not found
        
        Examples
        --------
        >>> chembl = CheMBL()
        >>> compound = chembl.get_compound(15)
        >>> print(compound['pref_name'])
        'ASPIRIN'
        >>> print(compound['canonical_smiles'])
        'CC(=O)Oc1ccccc1C(=O)O'
        >>> print(compound['synonyms'][:3])
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
                cs.standard_inchi_key,
                cs.molfile
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