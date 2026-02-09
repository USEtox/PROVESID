"""
ChEBI (Chemical Entities of Biological Interest) API interface.

This module provides a Python interface to the ChEBI REST API for retrieving
chemical compound information from the ChEBI database, as well as a local
SDF file parser for offline access to ChEBI data.

Author: USEtox team
Date: August 2025
"""

import requests
import logging
import time
import os
import pickle
import gzip
import shutil
from typing import Dict, List, Optional, Union, Any
import xml.etree.ElementTree as ET
from rdkit import Chem
import pandas as pd
from tqdm import tqdm
from .utils import data_path


class ChEBIError(Exception):
    """Custom exception for ChEBI API errors."""
    pass


class ChEBI:
    """
    Interface for the ChEBI (Chemical Entities of Biological Interest) REST API.
    
    The ChEBI database is a freely available dictionary of molecular entities 
    focused on 'small' chemical compounds. This class provides methods to search
    for and retrieve compound information from ChEBI.
    
    Attributes:
        base_url (str): Base URL for ChEBI API
        timeout (int): Request timeout in seconds
        session (requests.Session): HTTP session for connection pooling
    
    Example:
        >>> chebi = ChEBI()
        >>> compound = chebi.get_complete_entity(15377)  # ChEBI:15377 (water)
        >>> print(compound['chebiAsciiName'])
        water
    """
    
    def __init__(self, timeout: int = 30):
        """
        Initialize ChEBI API client.
        
        Args:
            timeout (int): Request timeout in seconds (default: 30)
        """
        self.base_url = "https://www.ebi.ac.uk/webservices/chebi/2.0/test"
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'PROVESID-ChEBI-Client/1.0',
            'Accept': 'application/xml, text/xml'
        })
        
        # Setup logging
        self.logger = logging.getLogger(__name__)
    
    def _make_request(self, endpoint: str, params: Optional[Dict] = None) -> requests.Response:
        """
        Make HTTP request to ChEBI API.
        
        Args:
            endpoint (str): API endpoint
            params (dict, optional): Query parameters
            
        Returns:
            requests.Response: HTTP response object
            
        Raises:
            ChEBIError: If request fails or returns error status
        """
        url = f"{self.base_url}/{endpoint}"
        
        try:
            response = self.session.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            return response
            
        except requests.exceptions.Timeout:
            raise ChEBIError(f"Request timeout after {self.timeout} seconds")
        except requests.exceptions.RequestException as e:
            raise ChEBIError(f"Request failed: {str(e)}")
    
    def _parse_xml_response(self, response: requests.Response) -> Any:
        """
        Parse XML response from ChEBI API.
        
        Args:
            response (requests.Response): HTTP response with XML content
            
        Returns:
            Any: Parsed XML data as dictionary, list, or string
            
        Raises:
            ChEBIError: If XML parsing fails
        """
        try:
            root = ET.fromstring(response.content)
            return self._xml_to_dict(root)
        except ET.ParseError as e:
            raise ChEBIError(f"Failed to parse XML response: {str(e)}")
    
    def _xml_to_dict(self, element: ET.Element) -> Union[Dict, List, str]:
        """
        Convert XML element to dictionary.
        
        Args:
            element (ET.Element): XML element to convert
            
        Returns:
            Union[dict, list, str]: Converted data structure
        """
        # Remove namespace from tag
        tag = element.tag.split('}')[-1] if '}' in element.tag else element.tag
        
        # If element has children, process recursively
        if len(element) > 0:
            result = {}
            for child in element:
                child_tag = child.tag.split('}')[-1] if '}' in child.tag else child.tag
                child_data = self._xml_to_dict(child)
                
                # Handle multiple elements with same tag
                if child_tag in result:
                    if not isinstance(result[child_tag], list):
                        result[child_tag] = [result[child_tag]]
                    result[child_tag].append(child_data)
                else:
                    result[child_tag] = child_data
            
            # Add attributes if any
            if element.attrib:
                result.update(element.attrib)
                
            return result
        else:
            # Leaf element - return text content
            return element.text if element.text else ""
    
    def get_complete_entity(self, chebi_id: Union[int, str]) -> Optional[Dict[str, Any]]:
        """
        Get complete entity information for a ChEBI ID.
        
        Args:
            chebi_id (Union[int, str]): ChEBI ID (with or without 'CHEBI:' prefix)
            
        Returns:
            dict: Complete entity information, None if not found
            
        Example:
            >>> chebi = ChEBI()
            >>> water = chebi.get_complete_entity(15377)
            >>> print(water['chebiAsciiName'])
            water
        """
        # Ensure ChEBI ID is properly formatted
        if isinstance(chebi_id, int):
            chebi_id = f"CHEBI:{chebi_id}"
        elif not str(chebi_id).startswith("CHEBI:"):
            chebi_id = f"CHEBI:{chebi_id}"
        
        try:
            response = self._make_request("getCompleteEntity", {"chebiId": chebi_id})
            data = self._parse_xml_response(response)
            
            # Navigate to the actual entity data
            if 'Body' in data and 'getCompleteEntityResponse' in data['Body']:
                response_data = data['Body']['getCompleteEntityResponse']
                # Handle case where response is an error string
                if isinstance(response_data, str):
                    self.logger.warning(f"ChEBI API returned error: {response_data}")
                    return None
                
                entity = response_data.get('return')
                return entity
            
            return None
            
        except ChEBIError as e:
            self.logger.warning(f"Failed to get complete entity for {chebi_id}: {e}")
            return None
    
    def get_lite_entity(self, chebi_id: Union[int, str]) -> Optional[Dict[str, Any]]:
        """
        Get lite entity information for a ChEBI ID (basic information only).
        
        Args:
            chebi_id (Union[int, str]): ChEBI ID (with or without 'CHEBI:' prefix)
            
        Returns:
            dict: Lite entity information, None if not found
        """
        # Ensure ChEBI ID is properly formatted
        if isinstance(chebi_id, int):
            chebi_id = f"CHEBI:{chebi_id}"
        elif not str(chebi_id).startswith("CHEBI:"):
            chebi_id = f"CHEBI:{chebi_id}"
        
        try:
            response = self._make_request("getLiteEntity", {"chebiId": chebi_id})
            data = self._parse_xml_response(response)
            
            # Navigate to the actual entity data
            if 'Body' in data and 'getLiteEntityResponse' in data['Body']:
                entity = data['Body']['getLiteEntityResponse'].get('return')
                return entity
            
            return None
            
        except ChEBIError as e:
            self.logger.warning(f"Failed to get lite entity for {chebi_id}: {e}")
            return None
    
    def search_by_name(self, search_text: str, search_category: str = "ALL", 
                      max_results: int = 50, stars: str = "ALL") -> List[Dict[str, Any]]:
        """
        Search ChEBI by compound name.
        
        Args:
            search_text (str): Text to search for
            search_category (str): Search category ('ALL', 'CHEBI_NAME', 'DEFINITION', etc.)
            max_results (int): Maximum number of results to return
            stars (str): Star category ('ALL', 'TWO_ONLY', 'THREE_ONLY')
            
        Returns:
            list: List of matching entities
            
        Example:
            >>> chebi = ChEBI()
            >>> results = chebi.search_by_name("water")
            >>> for result in results[:3]:
            ...     print(f"{result['chebiId']}: {result['chebiAsciiName']}")
        """
        params = {
            "search": search_text,
            "searchCategory": search_category,
            "maxResults": max_results,
            "stars": stars
        }
        
        try:
            response = self._make_request("getLiteEntity", params)
            data = self._parse_xml_response(response)
            
            # Navigate to the search results
            if 'Body' in data and 'getLiteEntityResponse' in data['Body']:
                results = data['Body']['getLiteEntityResponse'].get('return', [])
                if not isinstance(results, list):
                    results = [results] if results else []
                return results
            
            return []
            
        except ChEBIError as e:
            self.logger.warning(f"Search failed for '{search_text}': {e}")
            return []
    
    def get_ontology_parents(self, chebi_id: Union[int, str]) -> List[Dict[str, Any]]:
        """
        Get ontology parents for a ChEBI ID.
        
        Args:
            chebi_id (Union[int, str]): ChEBI ID (with or without 'CHEBI:' prefix)
            
        Returns:
            list: List of parent entities in the ontology
        """
        # Ensure ChEBI ID is properly formatted
        if isinstance(chebi_id, int):
            chebi_id = f"CHEBI:{chebi_id}"
        elif not str(chebi_id).startswith("CHEBI:"):
            chebi_id = f"CHEBI:{chebi_id}"
        
        try:
            response = self._make_request("getOntologyParents", {"chebiId": chebi_id})
            data = self._parse_xml_response(response)
            
            if 'Body' in data and 'getOntologyParentsResponse' in data['Body']:
                parents = data['Body']['getOntologyParentsResponse'].get('return', [])
                if not isinstance(parents, list):
                    parents = [parents] if parents else []
                return parents
            
            return []
            
        except ChEBIError as e:
            self.logger.warning(f"Failed to get ontology parents for {chebi_id}: {e}")
            return []
    
    def get_ontology_children(self, chebi_id: Union[int, str]) -> List[Dict[str, Any]]:
        """
        Get ontology children for a ChEBI ID.
        
        Args:
            chebi_id (Union[int, str]): ChEBI ID (with or without 'CHEBI:' prefix)
            
        Returns:
            list: List of child entities in the ontology
        """
        # Ensure ChEBI ID is properly formatted
        if isinstance(chebi_id, int):
            chebi_id = f"CHEBI:{chebi_id}"
        elif not str(chebi_id).startswith("CHEBI:"):
            chebi_id = f"CHEBI:{chebi_id}"
        
        try:
            response = self._make_request("getOntologyChildren", {"chebiId": chebi_id})
            data = self._parse_xml_response(response)
            
            if 'Body' in data and 'getOntologyChildrenResponse' in data['Body']:
                children = data['Body']['getOntologyChildrenResponse'].get('return', [])
                if not isinstance(children, list):
                    children = [children] if children else []
                return children
            
            return []
            
        except ChEBIError as e:
            self.logger.warning(f"Failed to get ontology children for {chebi_id}: {e}")
            return []
    
    def get_all_ontology_children_in_path(self, chebi_id: Union[int, str], 
                                        ontology_type: str = "is_a") -> List[Dict[str, Any]]:
        """
        Get all ontology children in path for a ChEBI ID.
        
        Args:
            chebi_id (Union[int, str]): ChEBI ID (with or without 'CHEBI:' prefix)
            ontology_type (str): Type of ontology relationship ('is_a', 'has_part', etc.)
            
        Returns:
            list: List of all children in the ontology path
        """
        # Ensure ChEBI ID is properly formatted
        if isinstance(chebi_id, int):
            chebi_id = f"CHEBI:{chebi_id}"
        elif not str(chebi_id).startswith("CHEBI:"):
            chebi_id = f"CHEBI:{chebi_id}"
        
        params = {
            "chebiId": chebi_id,
            "ontologyType": ontology_type
        }
        
        try:
            response = self._make_request("getAllOntologyChildrenInPath", params)
            data = self._parse_xml_response(response)
            
            if 'Body' in data and 'getAllOntologyChildrenInPathResponse' in data['Body']:
                children = data['Body']['getAllOntologyChildrenInPathResponse'].get('return', [])
                if not isinstance(children, list):
                    children = [children] if children else []
                return children
            
            return []
            
        except ChEBIError as e:
            self.logger.warning(f"Failed to get all ontology children for {chebi_id}: {e}")
            return []
    
    def get_structure(self, chebi_id: Union[int, str], structure_type: str = "mol") -> Optional[str]:
        """
        Get chemical structure for a ChEBI ID.
        
        Args:
            chebi_id (Union[int, str]): ChEBI ID (with or without 'CHEBI:' prefix)
            structure_type (str): Structure format ('mol', 'sdf', 'smiles', 'inchi')
            
        Returns:
            str: Chemical structure in requested format, None if not found
        """
        # Ensure ChEBI ID is properly formatted
        if isinstance(chebi_id, int):
            chebi_id = f"CHEBI:{chebi_id}"
        elif not str(chebi_id).startswith("CHEBI:"):
            chebi_id = f"CHEBI:{chebi_id}"
        
        params = {
            "chebiId": chebi_id,
            "structureType": structure_type.upper()
        }
        
        try:
            response = self._make_request("getStructure", params)
            data = self._parse_xml_response(response)
            
            if 'Body' in data and 'getStructureResponse' in data['Body']:
                response_data = data['Body']['getStructureResponse']
                # Handle case where response is an error string
                if isinstance(response_data, str):
                    self.logger.warning(f"ChEBI API returned error: {response_data}")
                    return None
                
                structure = response_data.get('return')
                return structure
            
            return None
            
        except ChEBIError as e:
            self.logger.warning(f"Failed to get structure for {chebi_id}: {e}")
            return None
    
    def batch_get_entities(self, chebi_ids: List[Union[int, str]], 
                          pause_time: float = 0.1) -> Dict[str, Dict[str, Any]]:
        """
        Get complete entity information for multiple ChEBI IDs.
        
        Args:
            chebi_ids (List[Union[int, str]]): List of ChEBI IDs
            pause_time (float): Pause between requests to be respectful to the API
            
        Returns:
            dict: Dictionary mapping ChEBI IDs to entity information
            
        Example:
            >>> chebi = ChEBI()
            >>> results = chebi.batch_get_entities([15377, 16236])  # water, ethanol
            >>> for chebi_id, data in results.items():
            ...     print(f"{chebi_id}: {data['chebiAsciiName']}")
        """
        results = {}
        
        for chebi_id in chebi_ids:
            # Format the ID for the key
            key = f"CHEBI:{chebi_id}" if not str(chebi_id).startswith("CHEBI:") else str(chebi_id)
            
            entity = self.get_complete_entity(chebi_id)
            if entity:
                results[key] = entity
            
            # Be respectful to the API
            if pause_time > 0:
                time.sleep(pause_time)
        
        return results
    
    def __repr__(self) -> str:
        """String representation of ChEBI client."""
        return f"ChEBI(base_url='{self.base_url}', timeout={self.timeout})"


# Convenience function for quick lookups
def get_chebi_entity(chebi_id: Union[int, str]) -> Optional[Dict[str, Any]]:
    """
    Convenience function to get ChEBI entity information.
    
    Args:
        chebi_id (Union[int, str]): ChEBI ID (with or without 'CHEBI:' prefix)
        
    Returns:
        dict: Entity information, None if not found
        
    Example:
        >>> from provesid import get_chebi_entity
        >>> water = get_chebi_entity(15377)
        >>> print(water['chebiAsciiName'])
        water
    """
    chebi = ChEBI()
    return chebi.get_complete_entity(chebi_id)


def search_chebi(search_text: str, max_results: int = 10) -> List[Dict[str, Any]]:
    """
    Convenience function to search ChEBI by name.
    
    Args:
        search_text (str): Text to search for
        max_results (int): Maximum number of results to return
        
    Returns:
        list: List of matching entities
        
    Example:
        >>> from provesid import search_chebi
        >>> results = search_chebi("aspirin")
        >>> for result in results[:3]:
        ...     print(f"{result['chebiId']}: {result['chebiAsciiName']}")
    """
    chebi = ChEBI()
    return chebi.search_by_name(search_text, max_results=max_results)


class ChebiSDF:
    """
    Parser for ChEBI SDF (Structure-Data File) for offline access to ChEBI data.
    
    This class provides efficient querying of the ChEBI SDF file containing
    ~190,000 compounds with structure data, chemical properties, synonyms,
    and cross-references to 80+ external databases.
    
    The class builds an index on first use for fast lookups. The index is
    cached to disk for faster subsequent initializations.
    
    If the SDF file is not found, it can be automatically downloaded from the
    ChEBI FTP server.
    
    Attributes:
        sdf_path (str): Path to ChEBI SDF file
        index (dict): In-memory index for fast lookups
        
    Example:
        >>> chebi_sdf = ChebiSDF()
        >>> compound = chebi_sdf.get_compound_by_id("CHEBI:15377")
        >>> print(compound['ChEBI NAME'])
        water
    """
    
    # Default download URL for ChEBI SDF file
    DEFAULT_SDF_URL = "https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/chebi.sdf.gz"
    
    def __init__(self, sdf_path: Optional[str] = None, rebuild_index: bool = False, 
                 auto_download: bool = True, sdf_url: Optional[str] = None):
        """
        Initialize ChebiSDF parser.
        
        Args:
            sdf_path (str, optional): Path to ChEBI SDF file. If None, uses default location.
            rebuild_index (bool): If True, rebuild index even if cache exists (default: False)
            auto_download (bool): If True, automatically download SDF file if not found (default: True)
            sdf_url (str, optional): Custom URL to download the SDF from. If None, uses default.
        """
        
        if sdf_path is None:
            data_dir = data_path()
            sdf_path = os.path.join(data_dir, 'chebi.sdf')
        
        self.sdf_path = sdf_path
        self.index_path = sdf_path + '.index.pkl'
        self.sdf_url = sdf_url or self.DEFAULT_SDF_URL
        self.logger = logging.getLogger(__name__)
        
        # Check if SDF file exists, download if needed
        if not os.path.exists(self.sdf_path):
            if auto_download:
                self.logger.info(f"ChEBI SDF file not found at: {self.sdf_path}")
                self.logger.info("Downloading ChEBI SDF file automatically...")
                self.download_sdf(url=self.sdf_url, force=False)
            else:
                raise FileNotFoundError(
                    f"ChEBI SDF file not found at: {self.sdf_path}\n"
                    f"Please run ChebiSDF.download_sdf() or set auto_download=True\n"
                    f"Or download manually from: https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/"
                )
        
        # Load or build index
        if rebuild_index or not os.path.exists(self.index_path):
            self.logger.info("Building index from SDF file...")
            self.index = self._build_index()
            self._save_index()
        else:
            self.logger.info("Loading cached index...")
            self.index = self._load_index()
    
    def download_sdf(self, url: Optional[str] = None, force: bool = False):
        """
        Download the ChEBI SDF file from the ChEBI FTP server.
        
        The file is downloaded as a gzip archive (~250 MB) and automatically
        extracted to the data directory (~868 MB uncompressed).
        
        Args:
            url (str, optional): URL to download from. If None, uses default ChEBI FTP URL.
            force (bool): If True, download even if file already exists (default: False)
            
        Returns:
            str: Path to the downloaded and extracted SDF file
            
        Raises:
            FileExistsError: If the file already exists and force=False
            requests.exceptions.RequestException: If the download fails
            
        Example:
            >>> chebi_sdf = ChebiSDF(auto_download=False)
            >>> chebi_sdf.download_sdf()  # Manually trigger download
        """
        download_url = url or self.sdf_url
        
        # Check if file already exists
        if os.path.exists(self.sdf_path) and not force:
            raise FileExistsError(
                f"ChEBI SDF file already exists at: {self.sdf_path}\n"
                f"Use force=True to overwrite"
            )
        
        # Create data directory if it doesn't exist
        data_dir = os.path.dirname(self.sdf_path)
        os.makedirs(data_dir, exist_ok=True)
        
        self.logger.info(f"Downloading ChEBI SDF from: {download_url}")
        self.logger.info(f"Destination: {self.sdf_path}")
        
        try:
            # Stream the download with progress bar
            response = requests.get(download_url, stream=True, timeout=60)
            response.raise_for_status()
            
            # Get total file size
            total_size = int(response.headers.get('content-length', 0))
            
            # Download gzipped file to temporary location
            gz_path = self.sdf_path + '.gz.tmp'
            with open(gz_path, 'wb') as f:
                if total_size > 0:
                    with tqdm(total=total_size, unit='B', unit_scale=True, 
                             desc="Downloading ChEBI SDF") as pbar:
                        for chunk in response.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                                pbar.update(len(chunk))
                else:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
            
            self.logger.info("Download complete. Extracting gzip archive...")
            
            # Extract gzipped file
            temp_path = self.sdf_path + '.tmp'
            with gzip.open(gz_path, 'rb') as f_in:
                with open(temp_path, 'wb') as f_out:
                    # Copy with progress bar
                    file_size = os.path.getsize(gz_path)
                    with tqdm(total=file_size, unit='B', unit_scale=True, 
                             desc="Extracting ChEBI SDF") as pbar:
                        while True:
                            chunk = f_in.read(8192)
                            if not chunk:
                                break
                            f_out.write(chunk)
                            pbar.update(len(chunk))
            
            # Move temp file to final location
            if os.path.exists(self.sdf_path):
                os.remove(self.sdf_path)
            os.rename(temp_path, self.sdf_path)
            
            # Clean up temporary gzip file
            if os.path.exists(gz_path):
                os.remove(gz_path)
            
            self.logger.info(f"✓ ChEBI SDF file downloaded and extracted successfully")
            self.logger.info(f"✓ File location: {self.sdf_path}")
            
            # Verify the file is valid by checking first few lines
            try:
                with open(self.sdf_path, 'r', encoding='utf-8', errors='ignore') as f:
                    first_lines = [f.readline() for _ in range(5)]
                    if not any(first_lines):
                        raise RuntimeError("Downloaded file appears to be empty")
                self.logger.info("✓ File verified successfully")
            except Exception as e:
                os.remove(self.sdf_path)
                raise RuntimeError(f"Downloaded file is corrupted: {e}")
            
            return self.sdf_path
            
        except requests.exceptions.RequestException as e:
            # Clean up temp files if they exist
            if os.path.exists(gz_path):
                os.remove(gz_path)
            if os.path.exists(temp_path):
                os.remove(temp_path)
            self.logger.error(f"Failed to download ChEBI SDF: {e}")
            raise
    
    def _build_index(self) -> Dict:
        """
        Build index from SDF file for fast lookups.
        
        Returns:
            dict: Index containing mappings for various query types
        """
        index = {
            'id_to_offset': {},           # ChEBI ID -> file offset
            'name_to_ids': {},            # lowercase name -> list of ChEBI IDs
            'inchikey_to_id': {},         # InChIKey -> ChEBI ID
            'inchi_to_id': {},            # InChI -> ChEBI ID
            'formula_to_ids': {},         # formula -> list of ChEBI IDs
            'cas_to_ids': {},             # CAS number -> list of ChEBI IDs
            'synonym_to_ids': {},         # lowercase synonym -> list of ChEBI IDs
        }
        
        # Parse SDF file and build index
        with open(self.sdf_path, 'r', encoding='utf-8', errors='ignore') as f:
            file_offset = 0
            current_mol_offset = 0
            in_mol = False
            current_data = {}
            
            # Use tqdm for progress bar
            file_size = os.path.getsize(self.sdf_path)
            pbar = tqdm(total=file_size, unit='B', unit_scale=True, desc="Indexing ChEBI SDF")
            
            for line in f:
                if not in_mol and line.strip():
                    # Start of a new molecule
                    in_mol = True
                    current_mol_offset = file_offset
                    current_data = {}
                
                # Check for property tags
                if line.startswith('> <'):
                    field_name = line.strip()[3:-1]  # Extract field name
                    value_line = next(f, '')
                    file_offset += len(line.encode('utf-8'))
                    file_offset += len(value_line.encode('utf-8'))
                    current_data[field_name] = value_line.strip()
                    pbar.update(len(line.encode('utf-8')) + len(value_line.encode('utf-8')))
                    continue
                
                # Check for end of molecule
                if line.startswith('$$$$'):
                    in_mol = False
                    
                    # Index this molecule
                    if 'ChEBI ID' in current_data:
                        chebi_id = current_data['ChEBI ID']
                        index['id_to_offset'][chebi_id] = current_mol_offset
                        
                        # Index by name
                        if 'ChEBI NAME' in current_data:
                            name_lower = current_data['ChEBI NAME'].lower()
                            if name_lower not in index['name_to_ids']:
                                index['name_to_ids'][name_lower] = []
                            index['name_to_ids'][name_lower].append(chebi_id)
                        
                        # Index by InChIKey
                        if 'INCHIKEY' in current_data:
                            index['inchikey_to_id'][current_data['INCHIKEY']] = chebi_id
                        
                        # Index by InChI
                        if 'INCHI' in current_data:
                            index['inchi_to_id'][current_data['INCHI']] = chebi_id
                        
                        # Index by formula
                        if 'FORMULA' in current_data:
                            formula = current_data['FORMULA']
                            if formula not in index['formula_to_ids']:
                                index['formula_to_ids'][formula] = []
                            index['formula_to_ids'][formula].append(chebi_id)
                        
                        # Index by CAS
                        if 'CAS Registry Numbers' in current_data:
                            cas_numbers = current_data['CAS Registry Numbers'].split(';')
                            for cas in cas_numbers:
                                cas = cas.strip()
                                if cas:
                                    if cas not in index['cas_to_ids']:
                                        index['cas_to_ids'][cas] = []
                                    index['cas_to_ids'][cas].append(chebi_id)
                        
                        # Index by synonyms
                        if 'SYNONYM' in current_data:
                            synonyms = current_data['SYNONYM'].split(';')
                            for syn in synonyms:
                                syn_lower = syn.strip().lower()
                                if syn_lower:
                                    if syn_lower not in index['synonym_to_ids']:
                                        index['synonym_to_ids'][syn_lower] = []
                                    index['synonym_to_ids'][syn_lower].append(chebi_id)
                
                file_offset += len(line.encode('utf-8'))
                pbar.update(len(line.encode('utf-8')))
            
            pbar.close()
        
        self.logger.info(f"Index built: {len(index['id_to_offset'])} compounds indexed")
        return index
    
    def _save_index(self):
        """Save index to disk for faster subsequent loads."""
        try:
            with open(self.index_path, 'wb') as f:
                pickle.dump(self.index, f, protocol=pickle.HIGHEST_PROTOCOL)
            self.logger.info(f"Index saved to {self.index_path}")
        except Exception as e:
            self.logger.warning(f"Failed to save index: {e}")
    
    def _load_index(self) -> Dict:
        """Load index from disk."""
        try:
            with open(self.index_path, 'rb') as f:
                index = pickle.load(f)
            self.logger.info(f"Index loaded: {len(index['id_to_offset'])} compounds")
            return index
        except Exception as e:
            self.logger.warning(f"Failed to load index: {e}. Rebuilding...")
            return self._build_index()
    
    def _read_mol_at_offset(self, offset: int) -> Dict[str, str]:
        """
        Read a molecule entry from the SDF file at a specific offset.
        
        Args:
            offset (int): File offset where molecule starts
            
        Returns:
            dict: Molecule data including molfile and properties
        """
        data = {'molfile': ''}
        
        with open(self.sdf_path, 'r', encoding='utf-8', errors='ignore') as f:
            f.seek(offset)
            in_molfile = True
            
            for line in f:
                if in_molfile:
                    data['molfile'] += line
                    if line.startswith('M  END'):
                        in_molfile = False
                    continue
                
                # Parse property tags
                if line.startswith('> <'):
                    field_name = line.strip()[3:-1]
                    value_line = next(f, '').strip()
                    data[field_name] = value_line
                    continue
                
                # End of molecule
                if line.startswith('$$$$'):
                    break
        
        return data
    
    def get_compound_by_id(self, chebi_id: str) -> Optional[Dict[str, str]]:
        """
        Get compound data by ChEBI ID.
        
        Args:
            chebi_id (str): ChEBI ID (e.g., "CHEBI:15377" or "15377")
            
        Returns:
            dict: Compound data, or None if not found
            
        Example:
            >>> chebi_sdf = ChebiSDF()
            >>> water = chebi_sdf.get_compound_by_id("CHEBI:15377")
            >>> print(water['ChEBI NAME'])
            water
        """
        # Normalize ID format
        if not chebi_id.startswith('CHEBI:'):
            chebi_id = f'CHEBI:{chebi_id}'
        
        offset = self.index['id_to_offset'].get(chebi_id)
        if offset is None:
            return None
        
        return self._read_mol_at_offset(offset)
    
    def search_by_name(self, name: str, exact: bool = True) -> List[Dict[str, str]]:
        """
        Search compounds by name.
        
        Args:
            name (str): Compound name to search for
            exact (bool): If True, exact match; if False, partial match (default: True)
            
        Returns:
            list: List of matching compound data
            
        Example:
            >>> chebi_sdf = ChebiSDF()
            >>> results = chebi_sdf.search_by_name("water")
            >>> print(len(results))
            1
        """
        name_lower = name.lower()
        results = []
        
        if exact:
            chebi_ids = self.index['name_to_ids'].get(name_lower, [])
        else:
            # Partial match
            chebi_ids = []
            for indexed_name, ids in self.index['name_to_ids'].items():
                if name_lower in indexed_name:
                    chebi_ids.extend(ids)
        
        for chebi_id in chebi_ids:
            compound = self.get_compound_by_id(chebi_id)
            if compound:
                results.append(compound)
        
        return results
    
    def search_by_synonym(self, synonym: str, exact: bool = True) -> List[Dict[str, str]]:
        """
        Search compounds by synonym.
        
        Args:
            synonym (str): Synonym to search for
            exact (bool): If True, exact match; if False, partial match (default: True)
            
        Returns:
            list: List of matching compound data
        """
        synonym_lower = synonym.lower()
        results = []
        
        if exact:
            chebi_ids = self.index['synonym_to_ids'].get(synonym_lower, [])
        else:
            # Partial match
            chebi_ids = []
            for indexed_syn, ids in self.index['synonym_to_ids'].items():
                if synonym_lower in indexed_syn:
                    chebi_ids.extend(ids)
            # Remove duplicates
            chebi_ids = list(set(chebi_ids))
        
        for chebi_id in chebi_ids:
            compound = self.get_compound_by_id(chebi_id)
            if compound:
                results.append(compound)
        
        return results
    
    def search_by_inchikey(self, inchikey: str) -> Optional[Dict[str, str]]:
        """
        Search compound by InChIKey.
        
        Args:
            inchikey (str): InChIKey to search for
            
        Returns:
            dict: Compound data, or None if not found
        """
        chebi_id = self.index['inchikey_to_id'].get(inchikey)
        if chebi_id:
            return self.get_compound_by_id(chebi_id)
        return None
    
    def search_by_inchi(self, inchi: str) -> Optional[Dict[str, str]]:
        """
        Search compound by InChI.
        
        Args:
            inchi (str): InChI string to search for
            
        Returns:
            dict: Compound data, or None if not found
        """
        chebi_id = self.index['inchi_to_id'].get(inchi)
        if chebi_id:
            return self.get_compound_by_id(chebi_id)
        return None
    
    def search_by_cas(self, cas: str) -> List[Dict[str, str]]:
        """
        Search compounds by CAS Registry Number.
        
        Args:
            cas (str): CAS Registry Number
            
        Returns:
            list: List of matching compound data
        """
        chebi_ids = self.index['cas_to_ids'].get(cas, [])
        results = []
        
        for chebi_id in chebi_ids:
            compound = self.get_compound_by_id(chebi_id)
            if compound:
                results.append(compound)
        
        return results
    
    def search_by_formula(self, formula: str) -> List[Dict[str, str]]:
        """
        Search compounds by molecular formula.
        
        Args:
            formula (str): Molecular formula (e.g., "H2O")
            
        Returns:
            list: List of matching compound data
        """
        chebi_ids = self.index['formula_to_ids'].get(formula, [])
        results = []
        
        for chebi_id in chebi_ids:
            compound = self.get_compound_by_id(chebi_id)
            if compound:
                results.append(compound)
        
        return results
    
    def filter_by_star_rating(self, min_stars: int = 3) -> List[str]:
        """
        Get ChEBI IDs of compounds with minimum star rating.
        
        Args:
            min_stars (int): Minimum star rating (1-3, default: 3)
            
        Returns:
            list: List of ChEBI IDs matching criteria
        """
        matching_ids = []
        
        for chebi_id in tqdm(self.index['id_to_offset'].keys(), desc="Filtering by star rating"):
            compound = self.get_compound_by_id(chebi_id)
            if compound and 'STAR' in compound:
                try:
                    star = int(compound['STAR'])
                    if star >= min_stars:
                        matching_ids.append(chebi_id)
                except ValueError:
                    continue
        
        return matching_ids
    
    def get_compounds_by_ids(self, chebi_ids: List[str]) -> List[Dict[str, str]]:
        """
        Get multiple compounds by ChEBI IDs.
        
        Args:
            chebi_ids (list): List of ChEBI IDs
            
        Returns:
            list: List of compound data dictionaries
        """
        results = []
        for chebi_id in chebi_ids:
            compound = self.get_compound_by_id(chebi_id)
            if compound:
                results.append(compound)
        return results
    
    def export_to_dataframe(self, chebi_ids: Optional[List[str]] = None, 
                           fields: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Export compounds to pandas DataFrame.
        
        Args:
            chebi_ids (list, optional): List of ChEBI IDs. If None, exports all compounds.
            fields (list, optional): List of field names to include. If None, includes common fields.
            
        Returns:
            pd.DataFrame: DataFrame with compound data
            
        Example:
            >>> chebi_sdf = ChebiSDF()
            >>> df = chebi_sdf.export_to_dataframe(["CHEBI:15377", "CHEBI:16236"])
            >>> print(df[['ChEBI ID', 'ChEBI NAME', 'FORMULA']])
        """
        if fields is None:
            fields = ['ChEBI ID', 'ChEBI NAME', 'STAR', 'FORMULA', 'MASS', 
                     'SMILES', 'INCHI', 'INCHIKEY', 'CAS Registry Numbers']
        
        if chebi_ids is None:
            chebi_ids = list(self.index['id_to_offset'].keys())
        
        data = []
        for chebi_id in tqdm(chebi_ids, desc="Exporting to DataFrame"):
            compound = self.get_compound_by_id(chebi_id)
            if compound:
                row = {field: compound.get(field, None) for field in fields}
                data.append(row)
        
        return pd.DataFrame(data)
    
    def get_database_stats(self) -> Dict[str, int]:
        """
        Get statistics about the ChEBI SDF database.
        
        Returns:
            dict: Statistics including counts of various indexed fields
        """
        return {
            'total_compounds': len(self.index['id_to_offset']),
            'compounds_with_inchikey': len(self.index['inchikey_to_id']),
            'compounds_with_inchi': len(self.index['inchi_to_id']),
            'compounds_with_cas': len(self.index['cas_to_ids']),
            'unique_formulas': len(self.index['formula_to_ids']),
            'indexed_names': len(self.index['name_to_ids']),
            'indexed_synonyms': len(self.index['synonym_to_ids']),
        }
