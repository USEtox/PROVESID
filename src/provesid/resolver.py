import requests
import time
import logging
from urllib.parse import quote
from typing import Dict, List, Union, Optional, Any

class NCIResolverError(Exception):
    """Custom exception for NCI Chemical Identifier Resolver errors"""
    pass

class NCIResolverNotFoundError(NCIResolverError):
    """Exception raised when chemical identifier is not found"""
    pass

class NCIResolverTimeoutError(NCIResolverError):
    """Exception raised when request times out"""
    pass

class NCIChemicalIdentifierResolver:
    """
    A Python interface to the NCI Chemical Identifier Resolver web service
    
    This class provides methods to interact with the NCI CADD Group's Chemical Identifier
    Resolver service for converting between different chemical structure identifiers.
    
    The service can resolve various types of chemical identifiers and convert them
    into different representations.
    
    URL API scheme: https://cactus.nci.nih.gov/chemical/structure/{identifier}/{representation}
    
    Usage examples:
        resolver = NCIChemicalIdentifierResolver()
        
        # Convert SMILES to InChI
        inchi = resolver.resolve('CCO', 'stdinchi')
        
        # Get all names for a compound
        names = resolver.resolve('aspirin', 'names')
        
        # Get molecular weight
        mw = resolver.resolve('caffeine', 'mw')
        
        # Get comprehensive molecular data
        mol_data = resolver.get_molecular_data('50-00-0')  # formaldehyde CAS
    """
    
    def __init__(self, base_url: str = "https://cactus.nci.nih.gov/chemical/structure", 
                 timeout: int = 30, pause_time: float = 0.1):
        """
        Initialize NCI Chemical Identifier Resolver client
        
        Args:
            base_url: Base URL for the NCI resolver service
            timeout: Request timeout in seconds
            pause_time: Minimum time between API calls in seconds
        """
        self.base_url = base_url.rstrip('/')
        self.timeout = timeout
        self.pause_time = pause_time
        self.last_request_time = 0
        
        # Available representation methods
        self.representations = {
            # Structure identifiers
            'stdinchi': 'Standard InChI',
            'stdinchikey': 'Standard InChIKey', 
            'smiles': 'Unique SMILES',
            'ficts': 'NCI/CADD FICTS identifier',
            'ficus': 'NCI/CADD FICuS identifier',
            'uuuuu': 'NCI/CADD uuuuu identifier',
            'hashisy': 'CACTVS HASHISY hashcode',
            # File formats
            'sdf': 'SD file format',
            # Names and properties
            'names': 'Chemical names list',
            'iupac_name': 'IUPAC name',
            'cas': 'CAS Registry Number',
            'mw': 'Molecular weight',
            'formula': 'Molecular formula',
            # Images
            'image': 'Chemical structure image',
            # Additional properties (may vary by compound)
            'exactmass': 'Exact mass',
            'charge': 'Formal charge',
            'h_bond_acceptor_count': 'Hydrogen bond acceptor count',
            'h_bond_donor_count': 'Hydrogen bond donor count',
            'rotor_count': 'Rotatable bond count',
            'effective_rotor_count': 'Effective rotor count',
            'ring_count': 'Ring count',
            'ringsys_count': 'Ring system count'
        }
    
    def _rate_limit(self):
        """Enforce rate limiting between requests"""
        if self.pause_time > 0:
            current_time = time.time()
            time_since_last = current_time - self.last_request_time
            if time_since_last < self.pause_time:
                time.sleep(self.pause_time - time_since_last)
            self.last_request_time = time.time()
    
    def _make_request(self, url: str) -> str:
        """
        Make HTTP request with error handling and rate limiting
        
        Args:
            url: Request URL
            
        Returns:
            Response text
            
        Raises:
            NCIResolverTimeoutError: If request times out
            NCIResolverNotFoundError: If identifier not found (404)
            NCIResolverError: For other HTTP errors
        """
        self._rate_limit()
        
        try:
            response = requests.get(url, timeout=self.timeout)
            
            if response.status_code == 200:
                return response.text.strip()
            elif response.status_code == 404:
                raise NCIResolverNotFoundError("Chemical identifier not found")
            elif response.status_code == 500:
                raise NCIResolverError("Internal server error")
            else:
                raise NCIResolverError(f"HTTP error {response.status_code}: {response.text}")
                
        except requests.Timeout:
            raise NCIResolverTimeoutError("Request timed out")
        except requests.RequestException as e:
            raise NCIResolverError(f"Request failed: {str(e)}")
    
    def _build_url(self, identifier: str, representation: str, xml_format: bool = False) -> str:
        """
        Build URL for the NCI resolver service
        
        Args:
            identifier: Chemical structure identifier
            representation: Desired representation
            xml_format: Whether to request XML format
            
        Returns:
            Complete URL string
        """
        # URL-encode the identifier for special characters
        encoded_identifier = quote(identifier, safe='')
        
        # Build URL components
        url_parts = [self.base_url, encoded_identifier, representation]
        
        if xml_format:
            url_parts.append('xml')
        
        return '/'.join(url_parts)
    
    def resolve(self, identifier: str, representation: str, xml_format: bool = False) -> str:
        """
        Resolve a chemical identifier to another representation
        
        Args:
            identifier: Input chemical identifier (name, SMILES, InChI, CAS, etc.)
            representation: Target representation (see self.representations for options)
            xml_format: Whether to request XML format response
            
        Returns:
            Resolved representation as string
            
        Raises:
            ValueError: If representation is not supported
            NCIResolverNotFoundError: If identifier cannot be resolved
            NCIResolverError: For other resolver errors
        """
        if not identifier or not identifier.strip():
            raise NCIResolverError("Empty or invalid identifier provided")
            
        if representation not in self.representations:
            available = ', '.join(self.representations.keys())
            raise ValueError(f"Unsupported representation '{representation}'. "
                           f"Available: {available}")
        
        url = self._build_url(identifier, representation, xml_format)
        return self._make_request(url)
    
    def get_available_representations(self) -> List[str]:
        """
        Get list of available representation types
        
        Returns:
            List of available representation keys
        """
        return list(self.representations.keys())
    
    def resolve_multiple(self, identifier: str, representations: List[str]) -> Dict[str, str]:
        """
        Resolve a single identifier to multiple representations
        
        Args:
            identifier: Input chemical identifier
            representations: List of target representations
            
        Returns:
            Dictionary mapping representation to resolved value
        """
        results = {}
        for representation in representations:
            try:
                results[representation] = self.resolve(identifier, representation)
            except NCIResolverError as e:
                results[representation] = None
                logging.warning(f"Failed to resolve {identifier} to {representation}: {e}")
        
        return results
    
    def get_molecular_data(self, identifier: str) -> Dict[str, Any]:
        """
        Get comprehensive molecular data for a chemical identifier
        
        This method attempts to retrieve multiple common properties and identifiers
        for a given chemical, similar to the original nci_cas_to_mol function.
        
        Args:
            identifier: Input chemical identifier
            
        Returns:
            Dictionary with molecular data and metadata
        """
        # Standard representations to retrieve
        standard_reps = [
            'stdinchi', 'stdinchikey', 'smiles', 'names', 'iupac_name', 
            'cas', 'mw', 'formula', 'ficts', 'ficus', 'uuuuu', 'hashisy'
        ]
        
        result = {
            'found_by': identifier,
            'success': True,
            'error': None,
            'available_data': {}
        }
        
        success_count = 0
        
        for rep in standard_reps:
            try:
                value = self.resolve(identifier, rep)
                
                # Process specific data types
                if rep == 'names':
                    # Split names by newline and filter empty strings
                    names_list = [name.strip() for name in value.split('\n') if name.strip()]
                    result['available_data'][rep] = names_list
                elif rep == 'mw':
                    # Try to convert molecular weight to float
                    try:
                        result['available_data'][rep] = float(value)
                    except ValueError:
                        result['available_data'][rep] = value
                else:
                    result['available_data'][rep] = value
                
                success_count += 1
                
            except NCIResolverError as e:
                result['available_data'][rep] = None
                logging.debug(f"Could not resolve {identifier} to {rep}: {e}")
        
        # Set overall success status
        if success_count == 0:
            result['success'] = False
            result['error'] = "No representations could be resolved"
        
        # Add convenience accessors for backwards compatibility
        data = result['available_data']
        result.update({
            'stdinchi': data.get('stdinchi'),
            'stdinchikey': data.get('stdinchikey'), 
            'smiles': data.get('smiles'),
            'names': data.get('names'),
            'iupac_name': data.get('iupac_name'),
            'cas': data.get('cas'),
            'mw': data.get('mw'),
            'formula': data.get('formula'),
            'ficts': data.get('ficts'),
            'ficus': data.get('ficus'), 
            'uuuuu': data.get('uuuuu'),
            'hashisy': data.get('hashisy'),
            'note': 'OK' if result['success'] else 'Error calling the NCI web API'
        })
        
        return result
    
    def get_image_url(self, identifier: str, image_format: str = 'gif', 
                     width: int = 200, height: int = 200) -> str:
        """
        Get URL for chemical structure image
        
        Args:
            identifier: Chemical identifier
            image_format: Image format ('gif' or 'png')
            width: Image width in pixels
            height: Image height in pixels
            
        Returns:
            URL for the structure image
        """
        url = self._build_url(identifier, 'image')
        
        # Add image format and size parameters
        params = []
        if image_format.lower() in ['gif', 'png']:
            params.append(f"format={image_format.lower()}")
        if width != 200 or height != 200:
            params.append(f"width={width}")
            params.append(f"height={height}")
        
        if params:
            url += '?' + '&'.join(params)
        
        return url
    
    def download_image(self, identifier: str, filename: str, 
                      image_format: str = 'gif', width: int = 200, height: int = 200) -> bool:
        """
        Download chemical structure image to file
        
        Args:
            identifier: Chemical identifier
            filename: Output filename
            image_format: Image format ('gif' or 'png')
            width: Image width in pixels
            height: Image height in pixels
            
        Returns:
            True if download successful, False otherwise
        """
        try:
            image_url = self.get_image_url(identifier, image_format, width, height)
            self._rate_limit()
            
            response = requests.get(image_url, timeout=self.timeout)
            response.raise_for_status()
            
            with open(filename, 'wb') as f:
                f.write(response.content)
            
            return True
            
        except Exception as e:
            logging.error(f"Failed to download image for {identifier}: {e}")
            return False
    
    def batch_resolve(self, identifiers: List[str], representation: str) -> Dict[str, str]:
        """
        Resolve multiple identifiers to a single representation
        
        Args:
            identifiers: List of chemical identifiers
            representation: Target representation
            
        Returns:
            Dictionary mapping identifier to resolved value (None if failed)
        """
        results = {}
        
        for identifier in identifiers:
            try:
                results[identifier] = self.resolve(identifier, representation)
            except NCIResolverError as e:
                results[identifier] = None
                logging.warning(f"Failed to resolve {identifier}: {e}")
        
        return results
    
    def is_valid_identifier(self, identifier: str) -> bool:
        """
        Check if an identifier can be resolved by the service
        
        Args:
            identifier: Chemical identifier to test
            
        Returns:
            True if identifier can be resolved, False otherwise
        """
        try:
            # Try to get SMILES as a basic test
            self.resolve(identifier, 'smiles')
            return True
        except NCIResolverError:
            return False
    
    def search_by_partial_name(self, partial_name: str) -> List[str]:
        """
        Search for compounds by partial name match
        Note: This is a basic implementation - the NCI service doesn't have
        a dedicated partial matching endpoint, so this tries the exact name first.
        
        Args:
            partial_name: Partial chemical name
            
        Returns:
            List of matching names (may be empty)
        """
        try:
            names = self.resolve(partial_name, 'names')
            return [name.strip() for name in names.split('\n') if name.strip()]
        except NCIResolverError:
            return []

# Convenience functions for backwards compatibility and ease of use

def nci_cas_to_mol(cas_rn: str) -> Dict[str, Any]:
    """
    Convert a CAS RN to a molecule data structure using the NCI web API
    
    This function maintains compatibility with the original nci_cas_to_mol function
    while using the new NCIChemicalIdentifierResolver class.
    
    Args:
        cas_rn: CAS Registry Number
        
    Returns:
        Dictionary with molecular data
    """
    resolver = NCIChemicalIdentifierResolver()
    return resolver.get_molecular_data(cas_rn)

def nci_id_to_mol(identifier: str) -> Dict[str, Any]:
    """
    Convert any chemical identifier to a molecule data structure
    
    Args:
        identifier: Chemical identifier (CAS, name, SMILES, InChI, etc.)
        
    Returns:
        Dictionary with molecular data
    """
    resolver = NCIChemicalIdentifierResolver()
    return resolver.get_molecular_data(identifier)

def nci_resolver(input_value: str, output_type: str, timeout: int = 30) -> Optional[str]:
    """
    Simple resolver function for converting between identifier types
    
    This function maintains compatibility with the original nci_resolver function.
    
    Args:
        input_value: Input chemical identifier
        output_type: Desired output representation
        timeout: Request timeout in seconds
        
    Returns:
        Resolved representation as string, None if failed
    """
    try:
        resolver = NCIChemicalIdentifierResolver(timeout=timeout)
        return resolver.resolve(input_value, output_type)
    except NCIResolverError:
        return None

def nci_smiles_to_names(smiles: str) -> List[str]:
    """
    Get chemical names for a SMILES string
    
    Args:
        smiles: SMILES string
        
    Returns:
        List of chemical names
    """
    try:
        resolver = NCIChemicalIdentifierResolver()
        names_str = resolver.resolve(smiles, 'names')
        return [name.strip() for name in names_str.split('\n') if name.strip()]
    except NCIResolverError:
        return []

def nci_name_to_smiles(name: str) -> Optional[str]:
    """
    Convert chemical name to SMILES
    
    Args:
        name: Chemical name
        
    Returns:
        SMILES string or None if not found
    """
    try:
        resolver = NCIChemicalIdentifierResolver()
        return resolver.resolve(name, 'smiles')
    except NCIResolverError:
        return None

def nci_inchi_to_smiles(inchi: str) -> Optional[str]:
    """
    Convert InChI to SMILES
    
    Args:
        inchi: InChI string
        
    Returns:
        SMILES string or None if not found
    """
    try:
        resolver = NCIChemicalIdentifierResolver()
        return resolver.resolve(inchi, 'smiles')
    except NCIResolverError:
        return None

def nci_cas_to_inchi(cas_rn: str) -> Optional[str]:
    """
    Convert CAS Registry Number to Standard InChI
    
    Args:
        cas_rn: CAS Registry Number
        
    Returns:
        Standard InChI string or None if not found
    """
    try:
        resolver = NCIChemicalIdentifierResolver()
        return resolver.resolve(cas_rn, 'stdinchi')
    except NCIResolverError:
        return None

def nci_get_molecular_weight(identifier: str) -> Optional[float]:
    """
    Get molecular weight for any chemical identifier
    
    Args:
        identifier: Chemical identifier
        
    Returns:
        Molecular weight as float or None if not found
    """
    try:
        resolver = NCIChemicalIdentifierResolver()
        mw_str = resolver.resolve(identifier, 'mw')
        return float(mw_str)
    except (NCIResolverError, ValueError):
        return None

def nci_get_formula(identifier: str) -> Optional[str]:
    """
    Get molecular formula for any chemical identifier
    
    Args:
        identifier: Chemical identifier
        
    Returns:
        Molecular formula string or None if not found
    """
    try:
        resolver = NCIChemicalIdentifierResolver()
        return resolver.resolve(identifier, 'formula')
    except NCIResolverError:
        return None