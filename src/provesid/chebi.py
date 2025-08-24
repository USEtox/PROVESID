"""
ChEBI (Chemical Entities of Biological Interest) API interface.

This module provides a Python interface to the ChEBI REST API for retrieving
chemical compound information from the ChEBI database.

Author: USEtox team
Date: August 2025
"""

import requests
import logging
import time
from typing import Dict, List, Optional, Union, Any
import xml.etree.ElementTree as ET


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
