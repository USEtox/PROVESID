# a pubchem package with a limited number of functionalities of pubchempy 
# but with a simpler interface that serves our purpose in PROVES

import requests
import time
import json
import logging
import re
import os
import pandas as pd
from functools import lru_cache
from typing import Dict, List, Union, Optional, Any
from urllib.parse import quote
from .cache import cached

pugrest_prolog = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
pause_between_calls = 0.2 # seconds

# create an enumerate class called domain
# <domain> = substance | compound | assay | gene | protein | pathway | taxonomy | cell
class Domain:
    SUBSTANCE = "substance"
    COMPOUND = "compound"
    ASSAY = "assay"
    GENE = "gene"
    PROTEIN = "protein"
    PATHWAY = "pathway"
    TAXONOMY = "taxonomy"
    CELL = "cell"

class CompoundDomainNamespace:
    # compound domain <namespace> = cid | name | smiles | inchi | sdf | inchikey | formula | <structure search> | <xref> | <mass> | listkey | <fast search>
    CID = "cid"
    NAME = "name"
    SMILES = "smiles"
    INCHI = "inchi"
    SDF = "sdf"
    INCHIKEY = "inchikey"
    FORMULA = "formula"
    STRUCTURE_SEARCH = "structure search"
    XREF = "xref"
    MASS = "mass"
    LISTKEY = "listkey"
    FAST_SEARCH = "fast search"

class SubstanceDomainNamespace:
    SID = "sid"
    SOURCEID = "sourceid"
    SOURCEALL = "sourceall"
    NAME = "name"
    XREF = "xref"
    LISTKEY = "listkey"

class AssayDomainNamespace:
    AID = "aid"
    LISTKEY = "listkey"
    TYPE = "type"
    SOURCEALL = "sourceall"
    TARGET = "target"
    ACTIVITY = "activity"

class StructureSearch:
    # <structure search> = { substructure | superstructure | similarity | identity } / { smiles | inchi | sdf | cid}
    SUBSTRUCTURE = "substructure"
    SUPERSTRUCTURE = "superstructure"
    SIMILARITY = "similarity"
    IDENTITY = "identity"

class StructureSearchQueryType:
    SMILES = "smiles"
    INCHI = "inchi"
    SDF = "sdf"
    CID = "cid"

class FastSearch:
    FASTIDENTITY = "fastidentity"
    FASTSIMILARITY_2D = "fastsimilarity_2d"
    FASTSIMILARITY_3D = "fastsimilarity_3d"
    FASTSUBSTRUCTURE = "fastsubstructure"
    FASTSUPERSTRUCTURE = "fastsuperstructure"
    FASTFORMULA = "fastformula"

class Operation:
    # Compound operations
    RECORD = "record"
    PROPERTY = "property"
    SYNONYMS = "synonyms"
    SIDS = "sids"
    CIDS = "cids"
    AIDS = "aids"
    ASSAYSUMMARY = "assaysummary"
    CLASSIFICATION = "classification"
    XREFS = "xrefs"
    DESCRIPTION = "description"
    CONFORMERS = "conformers"
    DATES = "dates"
    # Assay operations
    CONCISE = "concise"
    TARGETS = "targets"
    DOSERESPONSE = "doseresponse"
    SUMMARY = "summary"

class OutputFormat:
    XML = "XML"
    JSON = "JSON"
    JSONP = "JSONP"
    SDF = "SDF"
    CSV = "CSV"
    PNG = "PNG"
    TXT = "TXT"
    ASNT = "ASNT"
    ASNB = "ASNB"

class CompoundProperties:
    """Available compound properties for property tables"""
    MOLECULAR_FORMULA = "MolecularFormula"
    MOLECULAR_WEIGHT = "MolecularWeight"
    SMILES = "SMILES"
    CONNECTIVITY_SMILES = "ConnectivitySMILES"
    INCHI = "InChI"
    INCHIKEY = "InChIKey"
    IUPAC_NAME = "IUPACName"
    TITLE = "Title"
    XLOGP = "XLogP"
    EXACT_MASS = "ExactMass"
    MONOISOTOPIC_MASS = "MonoisotopicMass"
    TPSA = "TPSA"
    COMPLEXITY = "Complexity"
    CHARGE = "Charge"
    HBOND_DONOR_COUNT = "HBondDonorCount"
    HBOND_ACCEPTOR_COUNT = "HBondAcceptorCount"
    ROTATABLE_BOND_COUNT = "RotatableBondCount"
    HEAVY_ATOM_COUNT = "HeavyAtomCount"
    ISOTOPE_ATOM_COUNT = "IsotopeAtomCount"
    ATOM_STEREO_COUNT = "AtomStereoCount"
    DEFINED_ATOM_STEREO_COUNT = "DefinedAtomStereoCount"
    UNDEFINED_ATOM_STEREO_COUNT = "UndefinedAtomStereoCount"
    BOND_STEREO_COUNT = "BondStereoCount"
    DEFINED_BOND_STEREO_COUNT = "DefinedBondStereoCount"
    UNDEFINED_BOND_STEREO_COUNT = "UndefinedBondStereoCount"
    COVALENT_UNIT_COUNT = "CovalentUnitCount"
    PATENT_COUNT = "PatentCount"
    PATENT_FAMILY_COUNT = "PatentFamilyCount"
    ANNOTATION_TYPES = "AnnotationTypes"
    ANNOTATION_TYPE_COUNT = "AnnotationTypeCount"
    SOURCE_CATEGORIES = "SourceCategories"
    LITERATURE_COUNT = "LiteratureCount"
    VOLUME_3D = "Volume3D"
    X_STERIC_QUADRUPOLE_3D = "XStericQuadrupole3D"
    Y_STERIC_QUADRUPOLE_3D = "YStericQuadrupole3D"
    Z_STERIC_QUADRUPOLE_3D = "ZStericQuadrupole3D"
    FEATURE_COUNT_3D = "FeatureCount3D"
    FEATURE_ACCEPTOR_COUNT_3D = "FeatureAcceptorCount3D"
    FEATURE_DONOR_COUNT_3D = "FeatureDonorCount3D"
    FEATURE_ANION_COUNT_3D = "FeatureAnionCount3D"
    FEATURE_CATION_COUNT_3D = "FeatureCationCount3D"
    FEATURE_RING_COUNT_3D = "FeatureRingCount3D"
    FEATURE_HYDROPHOBE_COUNT_3D = "FeatureHydrophobeCount3D"
    CONFORMER_MODEL_RMSD_3D = "ConformerModelRMSD3D"
    EFFECTIVE_ROTOR_COUNT_3D = "EffectiveRotorCount3D"
    CONFORMER_COUNT_3D = "ConformerCount3D"
    FINGERPRINT_2D = "Fingerprint2D"

class PubChemError(Exception):
    """Custom exception for PubChem API errors"""
    pass

class PubChemTimeoutError(PubChemError):
    """Exception raised when request times out"""
    pass

class PubChemNotFoundError(PubChemError):
    """Exception raised when resource is not found"""
    pass

class PubChemServerError(PubChemError):
    """Exception raised when server error occurs"""
    pass

class PubChemAPI:
    """
    A Python interface to the PubChem REST API (PUG-REST)
    
    This class provides methods to interact with PubChem's REST API for retrieving
    chemical compound, substance, and assay information.
    
    Usage examples:
        api = PubChemAPI()
        
        # Get compound by CID
        compound = api.get_compound_by_cid(2244)
        
        # Get compound properties
        props = api.get_compound_properties([2244, 5793], ['MolecularFormula', 'MolecularWeight'])
        
        # Search by name
        compounds = api.get_compounds_by_name('aspirin')
        
        # Structure search
        similar = api.similarity_search('CCO', threshold=90)
    """
    
    def __init__(self, base_url: str = pugrest_prolog, pause_time: float = pause_between_calls, use_cache: bool = True):
        """
        Initialize PubChem API client
        
        Args:
            base_url: Base URL for PubChem REST API
            pause_time: Minimum time between API calls in seconds
            use_cache: Whether to use cache for lookups (default: True). 
                      When False, skips cache lookup but still stores results.
            
        Note:
            Caching is now unlimited by default with persistent storage.
            Use provesid.cache functions for cache management.
        """
        self.base_url = base_url.rstrip('/')
        self.pause_time = pause_time
        self.last_request_time = 0
        self.use_cache = use_cache
        
    def clear_cache(self):
        """Clear all cached results for PubChem API"""
        from .cache import clear_pubchem_cache
        clear_pubchem_cache()
            
    def get_cache_info(self):
        """Get cache statistics for PubChem API cached methods"""
        from .cache import get_pubchem_cache_info
        return get_pubchem_cache_info()
                    
        return cache_info
        
    def _rate_limit(self):
        """Enforce rate limiting between requests"""
        current_time = time.time()
        time_since_last = current_time - self.last_request_time
        if time_since_last < self.pause_time:
            time.sleep(self.pause_time - time_since_last)
        self.last_request_time = time.time()
    
    def _make_request(self, url: str, method: str = 'GET', data: Optional[Dict] = None, 
                     timeout: int = 30, headers: Optional[Dict] = None) -> requests.Response:
        """
        Make HTTP request with error handling and rate limiting
        
        Args:
            url: Request URL
            method: HTTP method (GET or POST)
            data: POST data
            timeout: Request timeout in seconds
            headers: HTTP headers
            
        Returns:
            Response object
            
        Raises:
            PubChemTimeoutError: If request times out
            PubChemNotFoundError: If resource not found (404)
            PubChemServerError: If server error occurs (5xx)
            PubChemError: For other HTTP errors
        """
        self._rate_limit()
        
        try:
            if method.upper() == 'GET':
                response = requests.get(url, timeout=timeout, headers=headers)
            elif method.upper() == 'POST':
                response = requests.post(url, data=data, timeout=timeout, headers=headers)
            else:
                raise ValueError(f"Unsupported HTTP method: {method}")
            
            # Handle different HTTP status codes
            if response.status_code == 200:
                return response
            elif response.status_code == 202:
                # Accepted - asynchronous operation pending
                logging.warning("Asynchronous operation pending - may need to poll for results")
                return response
            elif response.status_code == 400:
                raise PubChemError(f"Bad request: {response.text}")
            elif response.status_code == 404:
                raise PubChemNotFoundError("Resource not found")
            elif response.status_code == 405:
                raise PubChemError("Method not allowed")
            elif response.status_code == 500:
                raise PubChemServerError("Internal server error")
            elif response.status_code == 501:
                raise PubChemError("Not implemented")
            elif response.status_code == 503:
                raise PubChemServerError("Server busy - try again later")
            elif response.status_code == 504:
                raise PubChemTimeoutError("Request timed out")
            else:
                raise PubChemError(f"HTTP error {response.status_code}: {response.text}")
                
        except requests.Timeout:
            raise PubChemTimeoutError("Request timed out")
        except requests.RequestException as e:
            raise PubChemError(f"Request failed: {str(e)}")
    
    def _build_url(self, domain: str, namespace: str, identifiers: Union[str, int, List[Union[str, int]]], 
                   operation: Optional[str] = None, output_format: str = OutputFormat.JSON,
                   **options) -> str:
        """
        Build PubChem REST API URL
        
        Args:
            domain: API domain (compound, substance, assay, etc.)
            namespace: Namespace within domain (cid, name, smiles, etc.)
            identifiers: Single identifier or list of identifiers
            operation: Operation to perform
            output_format: Desired output format
            **options: Additional URL parameters
            
        Returns:
            Complete URL string
        """
        # Handle identifiers
        if isinstance(identifiers, list):
            identifiers_str = ','.join(map(str, identifiers))
        else:
            identifiers_str = str(identifiers)
        
        # URL-encode identifiers for special characters
        identifiers_str = quote(identifiers_str, safe=',')
        
        # Build base URL path
        url_parts = [self.base_url, domain, namespace, identifiers_str]
        
        if operation:
            url_parts.append(operation)
        
        if output_format:
            url_parts.append(output_format)
        
        url = '/'.join(url_parts)
        
        # Add query parameters
        if options:
            params = []
            for key, value in options.items():
                if isinstance(value, bool):
                    value = str(value).lower()
                params.append(f"{key}={quote(str(value))}")
            if params:
                url += '?' + '&'.join(params)
        
        return url
    
    def _parse_response(self, response: requests.Response, output_format: str = OutputFormat.JSON) -> Any:
        """
        Parse API response based on format
        
        Args:
            response: HTTP response object
            output_format: Expected output format
            
        Returns:
            Parsed response data
        """
        if output_format == OutputFormat.JSON:
            try:
                return response.json()
            except json.JSONDecodeError:
                return response.text
        elif output_format in [OutputFormat.XML, OutputFormat.SDF, OutputFormat.CSV, 
                              OutputFormat.TXT, OutputFormat.ASNT, OutputFormat.ASNB]:
            return response.text
        elif output_format == OutputFormat.PNG:
            return response.content
        else:
            return response.text
    
    # Public cached methods
    @cached(service='pubchem')
    def get_compound_by_cid(self, cid: Union[int, str], output_format: str = OutputFormat.JSON) -> Any:
        """
        Get compound record by CID
        
        Args:
            cid: Compound ID
            output_format: Desired output format
            
        Returns:
            Compound data (automatically extracts from PC_Compounds wrapper for JSON format)
        """
        url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.CID, cid, 
                             Operation.RECORD, output_format)
        response = self._make_request(url)
        result = self._parse_response(response, output_format)
        
        # For JSON format, automatically extract the compound data from the wrapper
        if output_format == OutputFormat.JSON and isinstance(result, dict):
            if "PC_Compounds" in result and isinstance(result["PC_Compounds"], list) and len(result["PC_Compounds"]) > 0:
                return result["PC_Compounds"][0]
        
        return result
    
    def _get_compounds_by_name_impl(self, name: str, output_format: str = OutputFormat.JSON,
                                   name_type: str = "word") -> Any:
        """Implementation method for get_compounds_by_name with caching"""
        url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.NAME, name,
                             Operation.RECORD, output_format, name_type=name_type)
        response = self._make_request(url)
        result = self._parse_response(response, output_format)
        
        # For JSON format, automatically extract the compound data from the wrapper
        if output_format == OutputFormat.JSON and isinstance(result, dict):
            if "PC_Compounds" in result and isinstance(result["PC_Compounds"], list) and len(result["PC_Compounds"]) > 0:
                # If there's only one compound, return it directly, otherwise return the list
                if len(result["PC_Compounds"]) == 1:
                    return result["PC_Compounds"][0]
                else:
                    return result["PC_Compounds"]
        
        return result

    @cached(service='pubchem')
    def get_compounds_by_name(self, name: str, output_format: str = OutputFormat.JSON,
                             name_type: str = "word") -> Any:
        """
        Get compounds by name
        
        Args:
            name: Compound name
            output_format: Desired output format
            name_type: Name search type ("word" or "complete")
            
        Returns:
            Compound data (automatically extracts from PC_Compounds wrapper for JSON format)
        """
        url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.NAME, name,
                             Operation.RECORD, output_format, name_type=name_type)
        response = self._make_request(url)
        result = self._parse_response(response, output_format)
        
        # For JSON format, automatically extract the compound data from the wrapper
        if output_format == OutputFormat.JSON and isinstance(result, dict):
            if "PC_Compounds" in result and isinstance(result["PC_Compounds"], list) and len(result["PC_Compounds"]) > 0:
                # If there's only one compound, return it directly, otherwise return the list
                if len(result["PC_Compounds"]) == 1:
                    return result["PC_Compounds"][0]
                else:
                    return result["PC_Compounds"]
        
        return result
    
    def _get_compounds_by_smiles_impl(self, smiles: str, output_format: str = OutputFormat.JSON) -> Any:
        """Implementation method for get_compounds_by_smiles with caching"""
        url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.SMILES, smiles,
                             Operation.RECORD, output_format)
        response = self._make_request(url)
        result = self._parse_response(response, output_format)
        
        # For JSON format, automatically extract the compound data from the wrapper
        if output_format == OutputFormat.JSON and isinstance(result, dict):
            if "PC_Compounds" in result and isinstance(result["PC_Compounds"], list) and len(result["PC_Compounds"]) > 0:
                # If there's only one compound, return it directly, otherwise return the list
                if len(result["PC_Compounds"]) == 1:
                    return result["PC_Compounds"][0]
                else:
                    return result["PC_Compounds"]
        
        return result

    @cached(service='pubchem')
    def get_compounds_by_smiles(self, smiles: str, output_format: str = OutputFormat.JSON) -> Any:
        """
        Get compounds by SMILES
        
        Args:
            smiles: SMILES string
            output_format: Desired output format
            
        Returns:
            Compound data (automatically extracts from PC_Compounds wrapper for JSON format)
        """
        url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.SMILES, smiles,
                             Operation.RECORD, output_format)
        response = self._make_request(url)
        result = self._parse_response(response, output_format)
        
        # For JSON format, automatically extract the compound data from the wrapper
        if output_format == OutputFormat.JSON and isinstance(result, dict):
            if "PC_Compounds" in result and isinstance(result["PC_Compounds"], list) and len(result["PC_Compounds"]) > 0:
                # If there's only one compound, return it directly, otherwise return the list
                if len(result["PC_Compounds"]) == 1:
                    return result["PC_Compounds"][0]
                else:
                    return result["PC_Compounds"]
        
        return result
    
    def _get_compounds_by_inchikey_impl(self, inchikey: str, output_format: str = OutputFormat.JSON) -> Any:
        """Implementation method for get_compounds_by_inchikey with caching"""
        url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.INCHIKEY, inchikey,
                             Operation.RECORD, output_format)
        response = self._make_request(url)
        result = self._parse_response(response, output_format)
        
        # For JSON format, automatically extract the compound data from the wrapper
        if output_format == OutputFormat.JSON and isinstance(result, dict):
            if "PC_Compounds" in result and isinstance(result["PC_Compounds"], list) and len(result["PC_Compounds"]) > 0:
                # If there's only one compound, return it directly, otherwise return the list
                if len(result["PC_Compounds"]) == 1:
                    return result["PC_Compounds"][0]
                else:
                    return result["PC_Compounds"]
        
        return result

    @cached(service='pubchem')
    def get_compounds_by_inchikey(self, inchikey: str, output_format: str = OutputFormat.JSON) -> Any:
        """
        Get compounds by InChIKey
        
        Args:
            inchikey: InChI Key
            output_format: Desired output format
            
        Returns:
            Compound data (automatically extracts from PC_Compounds wrapper for JSON format)
        """
        url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.INCHIKEY, inchikey,
                             Operation.RECORD, output_format)
        response = self._make_request(url)
        result = self._parse_response(response, output_format)
        
        # For JSON format, automatically extract the compound data from the wrapper
        if output_format == OutputFormat.JSON and isinstance(result, dict):
            if "PC_Compounds" in result and isinstance(result["PC_Compounds"], list) and len(result["PC_Compounds"]) > 0:
                # If there's only one compound, return it directly, otherwise return the list
                if len(result["PC_Compounds"]) == 1:
                    return result["PC_Compounds"][0]
                else:
                    return result["PC_Compounds"]
        
        return result
    
    @cached(service='pubchem')
    def _cached_get_compound_properties(self, cid: Union[int, str], 
                                       properties_tuple: tuple, 
                                       include_synonyms: bool = True,
                                       output_format: str = OutputFormat.JSON) -> Dict[str, Any]:
        """Cached implementation of get_compound_properties"""
        properties = list(properties_tuple)
        try:
            props_str = ','.join(properties)
            operation = f"{Operation.PROPERTY}/{props_str}"
            url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.CID, cid,
                                 operation, output_format)
            response = self._make_request(url)
            prop_data = self._parse_response(response, output_format)
            
            # Get synonyms if requested
            synonyms_data = None
            if include_synonyms:
                try:
                    synonyms_data = self.get_compound_synonyms(cid, output_format)
                except Exception as e:
                    # Don't fail the whole request if synonyms fail
                    synonyms_data = []
            
            # Extract properties from nested structure
            if prop_data and 'PropertyTable' in prop_data and 'Properties' in prop_data['PropertyTable']:
                properties_dict = prop_data['PropertyTable']['Properties'][0]
                
                # Add metadata keys
                properties_dict['success'] = True
                properties_dict['cid'] = cid
                properties_dict['error'] = None
                if include_synonyms:
                    properties_dict['synonyms'] = synonyms_data
                
                return properties_dict
            else:
                # Fallback for unexpected structure
                result = {
                    "success": False,
                    "cid": cid,
                    "error": "Unexpected response structure",
                    "raw_response": prop_data
                }
                if include_synonyms:
                    result['synonyms'] = synonyms_data
                return result
                
        except Exception as e:
            result = {
                "success": False,
                "cid": cid,
                "error": str(e)
            }
            if include_synonyms:
                result['synonyms'] = None
            return result

    @cached(service='pubchem')
    def get_compound_properties(self, cid: Union[int, str], 
                               properties: List[str], 
                               include_synonyms: bool = True,
                               output_format: str = OutputFormat.JSON) -> Dict[str, Any]:
        """
        Get compound properties and synonyms by CID with direct property access
        
        Args:
            cid: Single Compound ID
            properties: List of property names
            include_synonyms: Whether to include synonyms in the output
            output_format: Desired output format
            
        Returns:
            Dictionary with properties directly accessible at top level,
            plus 'success', 'cid', 'synonyms', and 'error' metadata keys
        """
        # Convert list to tuple for caching
        properties_tuple = tuple(properties)
        return self._cached_get_compound_properties(cid, properties_tuple, include_synonyms, output_format)

    def get_compound_properties_batch(self, cids: List[Union[int, str]], 
                                     properties: List[str], 
                                     output_format: str = OutputFormat.JSON) -> List[Dict[str, Any]]:
        """
        Get compound properties for multiple CIDs (legacy batch method)
        
        Args:
            cids: List of CIDs
            properties: List of property names
            output_format: Desired output format
            
        Returns:
            List of dictionaries, each with flat structure like get_compound_properties
        """
        results = []
        for cid in cids:
            try:
                result = self.get_compound_properties(cid, properties, include_synonyms=False, output_format=output_format)
                results.append(result)
            except Exception as e:
                results.append({
                    "success": False,
                    "cid": cid,
                    "error": str(e)
                })
        return results
    
    @cached(service='pubchem')
    def get_compound_synonyms(self, cid: Union[int, str], output_format: str = OutputFormat.JSON) -> List[str]:
        """
        Get compound synonyms by CID
        
        Args:
            cid: Compound ID
            output_format: Desired output format
            
        Returns:
            List of synonyms (flattened from nested structure)
        """
        try:
            url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.CID, cid,
                                 Operation.SYNONYMS, output_format)
            response = self._make_request(url)
            raw_data = self._parse_response(response, output_format)
            
            # Extract synonyms from nested structure
            if raw_data and 'InformationList' in raw_data:
                info_list = raw_data['InformationList'].get('Information', [])
                if info_list and len(info_list) > 0:
                    return info_list[0].get('Synonym', [])
            
            # Return empty list if no synonyms found
            return []
            
        except Exception:
            # Return empty list on error to maintain consistent return type
            return []
    
    @cached(service='pubchem')
    def get_cids_by_name(self, name: str, output_format: str = OutputFormat.JSON,
                        name_type: str = "word", domain: str = Domain.COMPOUND) -> Any:
        """
        Get CIDs by name from compound or substance domain
        
        Args:
            name: Compound or substance name
            output_format: Desired output format
            name_type: Name search type ("word" or "complete")
            domain: Search domain (Domain.COMPOUND or Domain.SUBSTANCE)
            
        Returns:
            CID list (extracted from nested response structure)
            
        Note:
            When searching in the substance domain, this can find CIDs for substances
            that may not be directly searchable in the compound domain.
        """
        # Choose appropriate namespace based on domain
        if domain == Domain.COMPOUND:
            namespace = CompoundDomainNamespace.NAME
        elif domain == Domain.SUBSTANCE:
            namespace = SubstanceDomainNamespace.NAME
        else:
            raise ValueError(f"Unsupported domain: {domain}. Use Domain.COMPOUND or Domain.SUBSTANCE")
        
        url = self._build_url(domain, namespace, name,
                             Operation.CIDS, output_format, name_type=name_type)
        response = self._make_request(url)
        parsed_response = self._parse_response(response, output_format)
        
        # Extract CID list from nested structure if JSON format
        if output_format == OutputFormat.JSON and isinstance(parsed_response, dict):
            # Handle compound domain response structure
            if 'IdentifierList' in parsed_response and 'CID' in parsed_response['IdentifierList']:
                return parsed_response['IdentifierList']['CID']
            # Handle substance domain response structure
            elif 'InformationList' in parsed_response and 'Information' in parsed_response['InformationList']:
                cids = []
                for info in parsed_response['InformationList']['Information']:
                    if 'CID' in info:
                        cids.extend(info['CID'])
                return list(set(cids))  # Remove duplicates and return unique CIDs
            elif 'Fault' in parsed_response:
                # Handle API fault response
                raise PubChemNotFoundError(f"No CIDs found for name: {name}")
        
        # Return original response for non-JSON formats or if structure is different
        return parsed_response
    
    @cached(service='pubchem')
    def get_cids_by_smiles(self, smiles: str, output_format: str = OutputFormat.JSON) -> Any:
        """
        Get CIDs by SMILES
        
        Args:
            smiles: SMILES string
            output_format: Desired output format
            
        Returns:
            CID list (extracted from nested response structure)
        """
        url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.SMILES, smiles,
                             Operation.CIDS, output_format)
        response = self._make_request(url)
        parsed_response = self._parse_response(response, output_format)
        
        # Extract CID list from nested structure if JSON format
        if output_format == OutputFormat.JSON and isinstance(parsed_response, dict):
            if 'IdentifierList' in parsed_response and 'CID' in parsed_response['IdentifierList']:
                return parsed_response['IdentifierList']['CID']
            elif 'Fault' in parsed_response:
                # Handle API fault response
                raise PubChemNotFoundError(f"No CIDs found for SMILES: {smiles}")
        
        # Return original response for non-JSON formats or if structure is different
        return parsed_response
    
    @cached(service='pubchem')
    def get_cids_by_inchikey(self, inchikey: str, output_format: str = OutputFormat.JSON) -> Any:
        """
        Get CIDs by InChI Key
        
        Args:
            inchikey: InChI Key string
            output_format: Desired output format
            
        Returns:
            CID list (extracted from nested response structure)
        """
        url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.INCHIKEY, inchikey,
                             Operation.CIDS, output_format)
        response = self._make_request(url)
        parsed_response = self._parse_response(response, output_format)
        
        # Extract CID list from nested structure if JSON format
        if output_format == OutputFormat.JSON and isinstance(parsed_response, dict):
            if 'IdentifierList' in parsed_response and 'CID' in parsed_response['IdentifierList']:
                return parsed_response['IdentifierList']['CID']
            elif 'Fault' in parsed_response:
                # Handle API fault response
                raise PubChemNotFoundError(f"No CIDs found for InChI Key: {inchikey}")
        
        # Return original response for non-JSON formats or if structure is different
        return parsed_response
    
    @cached(service='pubchem')
    def get_cids_by_formula(self, formula: str, output_format: str = OutputFormat.JSON,
                           allow_other_elements: bool = False) -> Any:
        """
        Get CIDs by molecular formula using fast search
        
        Args:
            formula: Molecular formula
            output_format: Desired output format
            allow_other_elements: Allow other elements beyond those specified
            
        Returns:
            CID list
        """
        url = self._build_url(Domain.COMPOUND, FastSearch.FASTFORMULA, formula,
                             Operation.CIDS, output_format, 
                             AllowOtherElements=allow_other_elements)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
    # Structure search methods
    def _make_options_hashable(self, **options: Any) -> tuple:
        """Convert options dict to a hashable tuple for caching"""
        if not options:
            return ()
        # Sort items to ensure consistent hashing
        sorted_items = tuple(sorted(options.items()))
        return sorted_items
    
    @cached(service='pubchem')
    def _cached_substructure_search(self, query: str, query_type: str, output_format: str, options_tuple: tuple) -> Any:
        """Cached implementation of substructure search"""
        options = dict(options_tuple) if options_tuple else {}
        search_type = f"{FastSearch.FASTSUBSTRUCTURE}/{query_type}"
        url = self._build_url(Domain.COMPOUND, search_type, query,
                             Operation.CIDS, output_format, **options)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
    def substructure_search(self, query: str, query_type: str = "smiles", 
                           output_format: str = OutputFormat.JSON, **options: Any) -> Any:
        """
        Perform substructure search
        
        Args:
            query: Query structure (SMILES, CID, etc.)
            query_type: Type of query (smiles, cid, etc.)
            output_format: Desired output format
            **options: Search options (MatchIsotopes, MaxRecords, etc.)
            
        Returns:
            Search results
        """
        options_tuple = self._make_options_hashable(**options)
        return self._cached_substructure_search(query, query_type, output_format, options_tuple)
    
    @cached(service='pubchem')
    def _cached_superstructure_search(self, query: str, query_type: str, output_format: str, options_tuple: tuple) -> Any:
        """Cached implementation of superstructure search"""
        options = dict(options_tuple) if options_tuple else {}
        search_type = f"{FastSearch.FASTSUPERSTRUCTURE}/{query_type}"
        url = self._build_url(Domain.COMPOUND, search_type, query,
                             Operation.CIDS, output_format, **options)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
    def superstructure_search(self, query: str, query_type: str = "smiles",
                             output_format: str = OutputFormat.JSON, **options: Any) -> Any:
        """
        Perform superstructure search
        
        Args:
            query: Query structure (SMILES, CID, etc.)
            query_type: Type of query (smiles, cid, etc.)
            output_format: Desired output format
            **options: Search options
            
        Returns:
            Search results
        """
        options_tuple = self._make_options_hashable(**options)
        return self._cached_superstructure_search(query, query_type, output_format, options_tuple)
    
    @cached(service='pubchem')
    def _cached_similarity_search(self, query: str, query_type: str, threshold: int, output_format: str, options_tuple: tuple) -> Any:
        """Cached implementation of similarity search"""
        options = dict(options_tuple) if options_tuple else {}
        search_type = f"{FastSearch.FASTSIMILARITY_2D}/{query_type}"
        url = self._build_url(Domain.COMPOUND, search_type, query,
                             Operation.CIDS, output_format, 
                             Threshold=threshold, **options)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
    @cached(service='pubchem')
    def similarity_search(self, query: str, query_type: str = "smiles",
                         threshold: int = 90, output_format: str = OutputFormat.JSON,
                         **options: Any) -> Any:
        """
        Perform 2D similarity search
        
        Args:
            query: Query structure (SMILES, CID, etc.)
            query_type: Type of query (smiles, cid, etc.)
            threshold: Similarity threshold (0-100)
            output_format: Desired output format
            **options: Search options
            
        Returns:
            Search results
        """
        options_tuple = self._make_options_hashable(**options)
        return self._cached_similarity_search(query, query_type, threshold, output_format, options_tuple)
    
    @cached(service='pubchem')
    def _cached_identity_search(self, query: str, query_type: str, identity_type: str, output_format: str, options_tuple: tuple) -> Any:
        """Cached implementation of identity search"""
        options = dict(options_tuple) if options_tuple else {}
        search_type = f"{FastSearch.FASTIDENTITY}/{query_type}"
        url = self._build_url(Domain.COMPOUND, search_type, query,
                             Operation.CIDS, output_format,
                             identity_type=identity_type, **options)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
    @cached(service='pubchem')
    def identity_search(self, query: str, query_type: str = "smiles",
                       identity_type: str = "same_stereo_isotope",
                       output_format: str = OutputFormat.JSON, **options: Any) -> Any:
        """
        Perform identity search
        
        Args:
            query: Query structure (SMILES, CID, etc.)
            query_type: Type of query (smiles, cid, etc.)
            identity_type: Type of identity match
            output_format: Desired output format
            **options: Search options
            
        Returns:
            Search results
        """
        options_tuple = self._make_options_hashable(**options)
        return self._cached_identity_search(query, query_type, identity_type, output_format, options_tuple)
    
    # Substance methods
    @cached(service='pubchem')
    def get_substance_by_sid(self, sid: Union[int, str], output_format: str = OutputFormat.JSON) -> Any:
        """
        Get substance by SID
        
        Args:
            sid: Substance ID
            output_format: Desired output format
            
        Returns:
            Substance data (automatically extracts from PC_Substances wrapper for JSON format)
        """
        url = self._build_url(Domain.SUBSTANCE, SubstanceDomainNamespace.SID, sid,
                             Operation.RECORD, output_format)
        response = self._make_request(url)
        result = self._parse_response(response, output_format)
        
        # For JSON format, automatically extract the substance data from the wrapper
        if output_format == OutputFormat.JSON and isinstance(result, dict):
            if "PC_Substances" in result and isinstance(result["PC_Substances"], list) and len(result["PC_Substances"]) > 0:
                return result["PC_Substances"][0]
        
        return result
    
    @cached(service='pubchem')
    def get_substances_by_name(self, name: str, output_format: str = OutputFormat.JSON) -> Any:
        """
        Get substances by name
        
        Args:
            name: Substance name
            output_format: Desired output format
            
        Returns:
            Substance data (automatically extracts from PC_Substances wrapper for JSON format)
        """
        url = self._build_url(Domain.SUBSTANCE, SubstanceDomainNamespace.NAME, name,
                             Operation.RECORD, output_format)
        response = self._make_request(url)
        result = self._parse_response(response, output_format)
        
        # For JSON format, automatically extract the substance data from the wrapper
        if output_format == OutputFormat.JSON and isinstance(result, dict):
            if "PC_Substances" in result and isinstance(result["PC_Substances"], list) and len(result["PC_Substances"]) > 0:
                # If there's only one substance, return it directly, otherwise return the list
                if len(result["PC_Substances"]) == 1:
                    return result["PC_Substances"][0]
                else:
                    return result["PC_Substances"]
        
        return result
    
    @cached(service='pubchem')
    def get_sids_by_name(self, name: str, output_format: str = OutputFormat.JSON,
                        sourcename: Optional[str] = None) -> Any:
        """
        Get SIDs by name
        
        Args:
            name: Substance name
            output_format: Desired output format
            sourcename: Restrict to specific source
            
        Returns:
            SID list (extracted from nested response structure)
        """
        options = {}
        if sourcename:
            options['sourcename'] = sourcename
            
        url = self._build_url(Domain.SUBSTANCE, SubstanceDomainNamespace.NAME, name,
                             Operation.SIDS, output_format, **options)
        response = self._make_request(url)
        parsed_response = self._parse_response(response, output_format)
        
        # Extract SID list from nested structure if JSON format
        if output_format == OutputFormat.JSON and isinstance(parsed_response, dict):
            if 'IdentifierList' in parsed_response and 'SID' in parsed_response['IdentifierList']:
                return parsed_response['IdentifierList']['SID']
            elif 'Fault' in parsed_response:
                # Handle API fault response
                raise PubChemNotFoundError(f"No SIDs found for name: {name}")
        
        # Return original response for non-JSON formats or if structure is different
        return parsed_response
    
    # Assay methods
    @cached(service='pubchem')
    def get_assay_by_aid(self, aid: Union[int, str], output_format: str = OutputFormat.JSON) -> Any:
        """
        Get assay by AID
        
        Args:
            aid: Assay ID
            output_format: Desired output format
            
        Returns:
            Assay data
        """
        url = self._build_url(Domain.ASSAY, AssayDomainNamespace.AID, aid,
                             Operation.RECORD, output_format)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
    @cached(service='pubchem')
    def get_assay_summary(self, cids: Union[int, str, List[Union[int, str]]],
                         output_format: str = OutputFormat.JSON) -> Any:
        """
        Get assay summary for compounds
        
        Args:
            cids: Single CID or list of CIDs
            output_format: Desired output format
            
        Returns:
            Assay summary data
        """
        url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.CID, cids,
                             Operation.ASSAYSUMMARY, output_format)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
    # Convenience methods for common use cases
    @cached(service='pubchem')
    def search_compound(self, query: str, search_type: str = "name") -> Dict[str, Any]:
        """
        Search for compound with automatic format detection
        
        Args:
            query: Search query (name, SMILES, InChIKey, etc.)
            search_type: Type of search ("name", "smiles", "inchikey", "cid")
            
        Returns:
            Dictionary with search results and metadata
        """
        try:
            if search_type == "name":
                result = self.get_compounds_by_name(query)
            elif search_type == "smiles":
                result = self.get_compounds_by_smiles(query)
            elif search_type == "inchikey":
                result = self.get_compounds_by_inchikey(query)
            elif search_type == "cid":
                result = self.get_compound_by_cid(query)
            else:
                raise ValueError(f"Unsupported search type: {search_type}")
            
            return {
                "success": True,
                "query": query,
                "search_type": search_type,
                "data": result,
                "error": None
            }
            
        except Exception as e:
            return {
                "success": False,
                "query": query,
                "search_type": search_type,
                "data": None,
                "error": str(e)
            }

    def format_search_compound_result(self, search_result: Dict[str, Any], 
                                      index: Optional[int] = None) -> Union[Dict[str, Any], List[Dict[str, Any]]]:
        """
        Convert search_compound output to a nicely formatted dictionary with flat structure
        
        Extracts properties from the nested props structure in search_compound results
        and creates a dictionary similar to get_all_compound_info output.
        
        Args:
            search_result: The result dictionary from search_compound()
            index: If search returns multiple results, specify which one to format (0-based).
                   If None and multiple results exist, returns a list of all formatted results.
            
        Returns:
            Dictionary with formatted properties, or list of dictionaries if multiple results.
            If search was unsuccessful or data is missing, returns a dictionary with 
            success=False and error message.
            
        Examples:
            >>> pch = PubChemAPI()
            >>> # Single result (e.g., CAS number)
            >>> res = pch.search_compound("50-00-0")
            >>> formatted = pch.format_search_compound_result(res)
            >>> print(formatted.get("MolecularFormula"))
            'CH2O'
            
            >>> # Multiple results (e.g., common name)
            >>> res = pch.search_compound("aspirin")
            >>> formatted = pch.format_search_compound_result(res, index=0)  # Get first result
            >>> # Or get all results
            >>> all_formatted = pch.format_search_compound_result(res)  # Returns list
        """
        # Check if search was successful
        if not search_result.get("success"):
            return {
                "success": False,
                "error": search_result.get("error", "Search was not successful"),
                "query": search_result.get("query"),
                "search_type": search_result.get("search_type")
            }
        
        # Check if data exists
        data = search_result.get("data")
        if not data:
            return {
                "success": False,
                "error": "No data found in search result",
                "query": search_result.get("query"),
                "search_type": search_result.get("search_type")
            }
        
        # Handle multiple results (list of compounds)
        if isinstance(data, list):
            if index is not None:
                # Format specific index
                if 0 <= index < len(data):
                    return self._format_single_compound(data[index], search_result)
                else:
                    return {
                        "success": False,
                        "error": f"Index {index} out of range (0-{len(data)-1})",
                        "query": search_result.get("query"),
                        "search_type": search_result.get("search_type")
                    }
            else:
                # Format all results
                return [self._format_single_compound(compound, search_result) 
                        for compound in data]
        
        # Handle single result (dict)
        if "props" not in data:
            return {
                "success": False,
                "error": "No properties data found in search result",
                "query": search_result.get("query"),
                "search_type": search_result.get("search_type")
            }
        
        return self._format_single_compound(data, search_result)
    
    def _format_single_compound(self, data: Dict[str, Any], 
                                search_result: Dict[str, Any]) -> Dict[str, Any]:
        """
        Internal helper to format a single compound's data
        
        Args:
            data: Single compound data dictionary with 'props' key
            search_result: Original search result for metadata
            
        Returns:
            Formatted dictionary with flat structure
        """
        # Initialize result dictionary with metadata
        formatted = {
            "success": True,
            "query": search_result.get("query"),
            "search_type": search_result.get("search_type"),
            "error": None
        }
        
        # Add CID if available
        if "id" in data and "id" in data["id"] and "cid" in data["id"]["id"]:
            formatted["CID"] = data["id"]["id"]["cid"]
        
        # Map of property labels/names to standardized keys
        # This maps the PubChem record format to property table format
        property_mapping = {
            ("Molecular Formula", ""): "MolecularFormula",
            ("Molecular Weight", ""): "MolecularWeight",
            ("SMILES", "Absolute"): "SMILES",
            ("SMILES", "Connectivity"): "ConnectivitySMILES",
            ("InChI", "Standard"): "InChI",
            ("InChIKey", "Standard"): "InChIKey",
            ("IUPAC Name", "Preferred"): "IUPACName",
            ("Log P", "XLogP3-AA"): "XLogP",
            ("Mass", "Exact"): "ExactMass",
            ("Weight", "MonoIsotopic"): "MonoisotopicMass",
            ("Topological", "Polar Surface Area"): "TPSA",
            ("Compound Complexity", ""): "Complexity",
            ("Count", "Hydrogen Bond Donor"): "HBondDonorCount",
            ("Count", "Hydrogen Bond Acceptor"): "HBondAcceptorCount",
            ("Count", "Rotatable Bond"): "RotatableBondCount",
            ("Fingerprint", "SubStructure Keys"): "Fingerprint2D",
            ("IUPAC Name", "Allowed"): "IUPACName_Allowed",
            ("IUPAC Name", "CAS-like Style"): "IUPACName_CASStyle",
            ("IUPAC Name", "Markup"): "IUPACName_Markup",
            ("IUPAC Name", "Systematic"): "IUPACName_Systematic",
            ("IUPAC Name", "Traditional"): "IUPACName_Traditional",
        }
        
        # Process each property
        for prop in data["props"]:
            urn = prop.get("urn", {})
            label = urn.get("label", "")
            name = urn.get("name", "")
            
            # Get value based on type
            value_obj = prop.get("value", {})
            if "sval" in value_obj:
                value = value_obj["sval"]
            elif "ival" in value_obj:
                value = value_obj["ival"]
            elif "fval" in value_obj:
                value = value_obj["fval"]
            elif "binary" in value_obj:
                value = value_obj["binary"]
            else:
                continue  # Skip if no recognized value type
            
            # Check if we have a mapping for this property
            key_tuple = (label, name)
            if key_tuple in property_mapping:
                formatted[property_mapping[key_tuple]] = value
            else:
                # For unmapped properties, create a key from label and name
                if name:
                    key = f"{label}_{name}".replace(" ", "").replace("-", "")
                else:
                    key = label.replace(" ", "").replace("-", "")
                formatted[key] = value
        
        return formatted

    @cached(service='pubchem')
    def get_basic_compound_info(self, cid: Union[int, str], 
                                include_synonyms: bool = False) -> Dict[str, Any]:
        """
        Get basic compound information including formula, molecular weight, 
        and structure, and IUPAC name. Synonyms can be included optionally.

        Args:
            cid: Compound ID
            include_synonyms: Whether to include synonyms in the response

        Returns:
            Dictionary with compound properties directly accessible at top level,
            plus 'success', 'cid', 'synonyms', and 'error' metadata keys
        """
        # Get basic properties with synonyms
        properties = [
            CompoundProperties.MOLECULAR_FORMULA,
            CompoundProperties.MOLECULAR_WEIGHT,
            CompoundProperties.SMILES,
            CompoundProperties.INCHI,
            CompoundProperties.INCHIKEY,
            CompoundProperties.IUPAC_NAME
        ]
        
        # Use the new get_compound_properties method which already includes synonyms and metadata
        return self.get_compound_properties(cid, properties, include_synonyms=include_synonyms)

    @cached(service='pubchem')
    def get_all_compound_info(self, cid: Union[int, str]) -> Dict[str, Any]:
        """
        Get all compound properties as listed in CompoundProperties

        Args:
            cid: Compound ID

        Returns:
            Dictionary with compound properties directly accessible at top level,
            plus 'success', 'cid', and 'error' metadata keys
        """
        # Get all property values from CompoundProperties class
        properties = []
        for attr_name in dir(CompoundProperties):
            if not attr_name.startswith("_"):
                prop_value = getattr(CompoundProperties, attr_name)
                if isinstance(prop_value, str):
                    properties.append(prop_value)
        
        # Use the new get_compound_properties method which already returns flat data
        return self.get_compound_properties(cid, properties, include_synonyms=False)

    def extract_identifiers_from_synonyms(self, synonyms: List[str]) -> Dict[str, List[str]]:
        """
        Extract chemical identifiers from a list of synonyms
        
        Args:
            synonyms: List of synonym strings
            
        Returns:
            Dictionary with lists of unique identifiers for each type:
            - casrn: CAS Registry Numbers (format: 2-5 digit-2 digit-single digit)
            - nsc: NSC numbers (begins with NSC)
            - dtxsid: DTXSID identifiers (begins with DTXSID)
            - dtxcid: DTXCID identifiers (begins with DTXCID)
            - ec_number: EC numbers (format: NNN-NNN-N)
            - chebi_id: ChEBI IDs (begins with CHEBI)
            - chembl: ChEMBL numbers (begins with CHEMBL)
        """
        identifiers = {
            'casrn': [],
            'nsc': [],
            'dtxsid': [],
            'dtxcid': [],
            'ec_number': [],
            'chebi_id': [],
            'chembl': []
        }
        
        for synonym in synonyms:
            if not isinstance(synonym, str):
                continue
                
            synonym_upper = synonym.upper().strip()
            
            # CAS Registry Number: 2-5 digits, hyphen, 2 digits, hyphen, 1 digit
            # May or may not begin with "CAS"
            cas_patterns = [
                r'\b(?:CAS\s*[:\-]?\s*)?(\d{2,7}-\d{2}-\d)\b',  # With optional CAS prefix and separators
            ]
            for pattern in cas_patterns:
                matches = re.findall(pattern, synonym_upper)
                for match in matches:
                    # Validate CAS number format more strictly
                    if re.match(r'^\d{2,7}-\d{2}-\d$', match):
                        if match not in identifiers['casrn']:
                            identifiers['casrn'].append(match)
            
            # NSC Number: begins with NSC
            nsc_match = re.search(r'\b(NSC\s*\d+)\b', synonym_upper)
            if nsc_match:
                nsc = nsc_match.group(1).replace(' ', '')
                if nsc not in identifiers['nsc']:
                    identifiers['nsc'].append(nsc)
            
            # DTXSID: begins with DTXSID
            dtxsid_match = re.search(r'\b(DTXSID\d+)\b', synonym_upper)
            if dtxsid_match:
                dtxsid = dtxsid_match.group(1)
                if dtxsid not in identifiers['dtxsid']:
                    identifiers['dtxsid'].append(dtxsid)
            
            # DTXCID: begins with DTXCID
            dtxcid_match = re.search(r'\b(DTXCID\d+)\b', synonym_upper)
            if dtxcid_match:
                dtxcid = dtxcid_match.group(1)
                if dtxcid not in identifiers['dtxcid']:
                    identifiers['dtxcid'].append(dtxcid)
            
            # EC Number: standard format is N.N.N.N (enzyme classification)
            # Only accept the standard dot-separated format
            ec_pattern = r'\b(?:EC\s*[:\-]?\s*)?(\d{1,2}\.\d{1,3}\.\d{1,3}\.(?:\d{1,3}|\-))\b'
            matches = re.findall(ec_pattern, synonym_upper)
            for match in matches:
                if match not in identifiers['ec_number']:
                    identifiers['ec_number'].append(match)
            
            # ChEBI ID: begins with CHEBI
            chebi_match = re.search(r'\b(CHEBI:?\s*\d+)\b', synonym_upper)
            if chebi_match:
                chebi = chebi_match.group(1).replace(' ', '').replace(':', ':')
                # Standardize format to CHEBI:XXXXX
                if not chebi.startswith('CHEBI:'):
                    chebi = chebi.replace('CHEBI', 'CHEBI:')
                if chebi not in identifiers['chebi_id']:
                    identifiers['chebi_id'].append(chebi)
            
            # ChEMBL: begins with CHEMBL
            chembl_match = re.search(r'\b(CHEMBL\d+)\b', synonym_upper)
            if chembl_match:
                chembl = chembl_match.group(1)
                if chembl not in identifiers['chembl']:
                    identifiers['chembl'].append(chembl)
        
        return identifiers

    def get_compound_identifiers(self, cid: Union[int, str]) -> Dict[str, Any]:
        """
        Get compound identifiers extracted from synonyms
        
        Args:
            cid: Compound ID
            
        Returns:
            Dictionary with 'success', 'cid', 'error' metadata and extracted identifiers
        """
        try:
            # Get synonyms
            synonyms_list = self.get_compound_synonyms(cid)
            
            # Extract identifiers
            identifiers = self.extract_identifiers_from_synonyms(synonyms_list)
            
            # Add metadata
            result = {
                'success': True,
                'cid': cid,
                'error': None,
                'total_synonyms': len(synonyms_list)
            }
            result.update(identifiers)
            
            return result
            
        except Exception as e:
            return {
                'success': False,
                'cid': cid,
                'error': str(e),
                'casrn': [],
                'nsc': [],
                'dtxsid': [],
                'dtxcid': [],
                'ec_number': [],
                'chebi_id': [],
                'chembl': [],
                'total_synonyms': 0
            }

    def find_cids_comprehensive(self, name: str, name_type: str = "word") -> Dict[str, Any]:
        """
        Search for CIDs in both compound and substance domains
        
        This method first searches in the compound domain, and if no results are found,
        it searches in the substance domain. This is useful for comprehensive searching
        when you're not sure which domain contains the identifier.
        
        Args:
            name: Compound or substance name (including CAS numbers, trade names, etc.)
            name_type: Name search type ("word" or "complete")
            
        Returns:
            Dictionary with search results from both domains
        """
        results = {
            "query": name,
            "name_type": name_type,
            "compound_domain": {"cids": [], "success": False, "error": None},
            "substance_domain": {"cids": [], "success": False, "error": None},
            "total_unique_cids": [],
            "recommended_domain": None
        }
        
        # Try compound domain first
        try:
            compound_cids = self.get_cids_by_name(name, name_type=name_type, domain=Domain.COMPOUND)
            results["compound_domain"]["cids"] = compound_cids
            results["compound_domain"]["success"] = True
            results["total_unique_cids"].extend(compound_cids)
        except Exception as e:
            results["compound_domain"]["error"] = str(e)
        
        # Try substance domain
        try:
            substance_cids = self.get_cids_by_name(name, name_type=name_type, domain=Domain.SUBSTANCE)
            results["substance_domain"]["cids"] = substance_cids
            results["substance_domain"]["success"] = True
            results["total_unique_cids"].extend(substance_cids)
        except Exception as e:
            results["substance_domain"]["error"] = str(e)
        
        # Remove duplicates and determine recommended domain
        results["total_unique_cids"] = list(set(results["total_unique_cids"]))
        
        if results["compound_domain"]["success"] and results["substance_domain"]["success"]:
            # Both succeeded - recommend the one with more results
            compound_count = len(results["compound_domain"]["cids"])
            substance_count = len(results["substance_domain"]["cids"])
            results["recommended_domain"] = "compound" if compound_count >= substance_count else "substance"
        elif results["compound_domain"]["success"]:
            results["recommended_domain"] = "compound"
        elif results["substance_domain"]["success"]:
            results["recommended_domain"] = "substance"
        else:
            results["recommended_domain"] = None
        
        return results


class PubChemID:
    """
    Interface to PubChem ID SQLite database for fast identifier lookup and conversion.
    
    This class provides access to a local SQLite database containing ~1.6M PubChem compounds
    with their identifiers (CID, CAS, InChI, InChIKey, SMILES) and chemical properties
    (molecular formula, molecular weight, LogP, complexity, etc.).
    
    The database is built from PubChem_CAS_202601.csv using the build_pubchem_id_db.py script.
    
    Attributes:
        db_path (str): Path to the SQLite database file
        conn (sqlite3.Connection): Database connection
    
    Example:
        >>> from provesid import PubChemID
        >>> db = PubChemID()
        >>> 
        >>> # Lookup by CAS
        >>> result = db.get_by_cas("50-78-2")  # Aspirin
        >>> print(result['inchi'])
        >>> 
        >>> # Lookup by InChIKey
        >>> result = db.get_by_inchikey("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        >>> print(result['cid'])
        >>> 
        >>> # Convert CAS to InChI
        >>> inchi = db.cas_to_inchi("50-78-2")
        >>> 
        >>> # Batch conversion
        >>> results = db.batch_cas_to_cid(["50-78-2", "50-00-0"])
    """
    
    def __init__(self, db_path: Optional[str] = None, auto_download: bool = True):
        """
        Initialize PubChemID database connection.
        
        Args:
            db_path (str, optional): Path to SQLite database. If None, uses default
                                    location in data directory.
            auto_download (bool): If True, automatically download database if not found.
                                 Default is True.
        
        Raises:
            FileNotFoundError: If database file doesn't exist and auto_download is False
        """
        import sqlite3
        from .utils import data_path
        
        if db_path is None:
            db_path = os.path.join(data_path(), 'pubchem_id.db')
        
        self.db_path = db_path
        
        if not os.path.exists(db_path):
            if auto_download:
                print(f"Database not found at {db_path}")
                print("Attempting to download from Zenodo...")
                self.download_database()
            else:
                raise FileNotFoundError(
                    f"PubChem ID database not found at {db_path}. "
                    "Set auto_download=True or run scripts/build_pubchem_id_db.py to create it."
                )
        
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row  # Access columns by name
    
    def __del__(self):
        """Close database connection on deletion."""
        if hasattr(self, 'conn'):
            self.conn.close()
    
    @staticmethod
    def download_database(db_path: Optional[str] = None, zenodo_url: Optional[str] = None) -> str:
        """
        Download PubChem ID database from Zenodo.
        
        Args:
            db_path (str, optional): Path where to save the database. If None, uses default
                                    location in data directory.
            zenodo_url (str, optional): URL to download from. If None, uses default Zenodo URL.
                                       Format: https://zenodo.org/record/XXXXXX/files/pubchem_id.db
        
        Returns:
            str: Path to the downloaded database file
        
        Example:
            >>> from provesid import PubChemID
            >>> # Download to default location
            >>> PubChemID.download_database()
            >>> 
            >>> # Or specify custom location
            >>> PubChemID.download_database(db_path='/path/to/pubchem_id.db')
        
        Note:
            After uploading to Zenodo, update the zenodo_url parameter with the actual URL.
            The database file is ~2.2 GB, so download may take several minutes.
        """
        import requests
        from .utils import data_path
        from tqdm import tqdm
        
        if db_path is None:
            db_path = os.path.join(data_path(), 'pubchem_id.db')
        
        if zenodo_url is None:
            zenodo_url = "https://zenodo.org/records/18173204/files/pubchem_id.db"
        
        print(f"Downloading PubChem ID database from Zenodo...")
        print(f"URL: {zenodo_url}")
        print(f"Destination: {db_path}")
        print("This is a large file (~2.2 GB), please be patient...")
        
        # Create temporary file path
        temp_path = db_path + '.tmp'
        
        try:
            # Download with progress bar
            response = requests.get(zenodo_url, stream=True)
            response.raise_for_status()
            
            total_size = int(response.headers.get('content-length', 0))
            
            with open(temp_path, 'wb') as f:
                with tqdm(total=total_size, unit='B', unit_scale=True, desc="Downloading") as pbar:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))
            
            print("Download complete. Verifying...")
            
            # Verify it's a valid SQLite database
            import sqlite3
            try:
                conn = sqlite3.connect(temp_path)
                cursor = conn.cursor()
                cursor.execute("SELECT COUNT(*) FROM compounds")
                count = cursor.fetchone()[0]
                conn.close()
                print(f"✓ Database verified: {count:,} compounds")
            except Exception as e:
                raise RuntimeError(f"Downloaded file is not a valid database: {e}")
            
            # Move to final location
            if os.path.exists(db_path):
                os.remove(db_path)
            os.rename(temp_path, db_path)
            
            print(f"✓ Database ready at {db_path}")
            return db_path
            
        except requests.exceptions.RequestException as e:
            if os.path.exists(temp_path):
                os.remove(temp_path)
            raise RuntimeError(f"Failed to download database: {e}")
        except Exception as e:
            if os.path.exists(temp_path):
                os.remove(temp_path)
            raise RuntimeError(f"Error during download: {e}")
    
    def get_by_cid(self, cid: int) -> Optional[Dict[str, Any]]:
        """
        Get compound information by PubChem CID.
        
        Args:
            cid (int): PubChem Compound ID
        
        Returns:
            dict: Compound information including identifiers and properties, or None if not found
        
        Example:
            >>> db = PubChemID()
            >>> result = db.get_by_cid(2244)  # Aspirin
            >>> print(result['cmpdname'])
            'Aspirin'
        """
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT * FROM compounds WHERE cid = ?
        """, (cid,))
        
        row = cursor.fetchone()
        if not row:
            return None
        
        result = dict(row)
        
        # Add CAS numbers
        cursor.execute("""
            SELECT cas FROM cas_numbers WHERE cid = ?
        """, (cid,))
        result['cas_numbers'] = [r[0] for r in cursor.fetchall()]
        
        # Add synonyms
        cursor.execute("""
            SELECT synonym FROM synonyms WHERE cid = ? LIMIT 100
        """, (cid,))
        result['synonyms'] = [r[0] for r in cursor.fetchall()]
        
        return result
    
    def get_by_cas(self, cas: str) -> Optional[Dict[str, Any]]:
        """
        Get compound information by CAS Registry Number.
        
        Args:
            cas (str): CAS Registry Number (e.g., "50-78-2")
        
        Returns:
            dict: Compound information, or None if not found
        
        Example:
            >>> db = PubChemID()
            >>> result = db.get_by_cas("50-78-2")  # Aspirin
            >>> print(result['inchi'])
        """
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT cid FROM cas_numbers WHERE cas = ? LIMIT 1
        """, (cas,))
        
        row = cursor.fetchone()
        if not row:
            return None
        
        return self.get_by_cid(row[0])
    
    def get_by_inchikey(self, inchikey: str) -> Optional[Dict[str, Any]]:
        """
        Get compound information by InChIKey.
        
        Args:
            inchikey (str): Standard InChIKey (27 characters)
        
        Returns:
            dict: Compound information, or None if not found
        
        Example:
            >>> db = PubChemID()
            >>> result = db.get_by_inchikey("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
            >>> print(result['cmpdname'])
        """
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT * FROM compounds WHERE inchikey = ?
        """, (inchikey,))
        
        row = cursor.fetchone()
        if not row:
            return None
        
        cid = row['cid']
        return self.get_by_cid(cid)
    
    def get_by_inchi(self, inchi: str) -> Optional[Dict[str, Any]]:
        """
        Get compound information by InChI string.
        
        Args:
            inchi (str): Standard InChI string
        
        Returns:
            dict: Compound information, or None if not found
        
        Example:
            >>> db = PubChemID()
            >>> result = db.get_by_inchi("InChI=1S/C9H8O4/c1-6(10)...")
            >>> print(result['cmpdname'])
        """
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT * FROM compounds WHERE inchi = ?
        """, (inchi,))
        
        row = cursor.fetchone()
        if not row:
            return None
        
        cid = row['cid']
        return self.get_by_cid(cid)
    
    def search_by_name(self, name: str, exact: bool = False, limit: int = 10) -> List[Dict[str, Any]]:
        """
        Search compounds by name or synonym.
        
        Args:
            name (str): Compound name or synonym to search for
            exact (bool): If True, exact match only. If False, partial match (case-insensitive)
            limit (int): Maximum number of results to return
        
        Returns:
            list: List of matching compounds
        
        Example:
            >>> db = PubChemID()
            >>> results = db.search_by_name("aspirin", exact=False)
            >>> for r in results:
            ...     print(r['cid'], r['cmpdname'])
        """
        cursor = self.conn.cursor()
        
        results = []
        
        if exact:
            # Search in main compound name
            cursor.execute("""
                SELECT cid FROM compounds WHERE cmpdname = ? LIMIT ?
            """, (name, limit))
            
            cids = [r[0] for r in cursor.fetchall()]
            
            # Also search in synonyms
            if len(cids) < limit:
                cursor.execute("""
                    SELECT DISTINCT cid FROM synonyms WHERE synonym = ? LIMIT ?
                """, (name, limit - len(cids)))
                cids.extend([r[0] for r in cursor.fetchall()])
        else:
            # Partial match with LIKE
            search_term = f"%{name}%"
            
            # Search in main compound name
            cursor.execute("""
                SELECT cid FROM compounds WHERE cmpdname LIKE ? LIMIT ?
            """, (search_term, limit))
            
            cids = [r[0] for r in cursor.fetchall()]
            
            # Also search in synonyms
            if len(cids) < limit:
                cursor.execute("""
                    SELECT DISTINCT cid FROM synonyms WHERE synonym LIKE ? LIMIT ?
                """, (search_term, limit - len(cids)))
                cids.extend([r[0] for r in cursor.fetchall()])
        
        # Get full compound info for each CID
        for cid in cids[:limit]:
            compound = self.get_by_cid(cid)
            if compound:
                results.append(compound)
        
        return results
    
    def search_by_formula(self, formula: str, limit: int = 100) -> List[Dict[str, Any]]:
        """
        Search compounds by molecular formula.
        
        Args:
            formula (str): Molecular formula (e.g., "C9H8O4")
            limit (int): Maximum number of results to return
        
        Returns:
            list: List of matching compounds
        
        Example:
            >>> db = PubChemID()
            >>> results = db.search_by_formula("C9H8O4")
            >>> print(f"Found {len(results)} compounds with formula C9H8O4")
        """
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT cid FROM compounds WHERE mf = ? LIMIT ?
        """, (formula, limit))
        
        results = []
        for row in cursor.fetchall():
            compound = self.get_by_cid(row[0])
            if compound:
                results.append(compound)
        
        return results
    
    # Conversion methods
    
    def cas_to_cid(self, cas: str) -> Optional[int]:
        """Convert CAS number to PubChem CID."""
        result = self.get_by_cas(cas)
        return result['cid'] if result else None
    
    def cas_to_inchi(self, cas: str) -> Optional[str]:
        """Convert CAS number to InChI."""
        result = self.get_by_cas(cas)
        return result['inchi'] if result else None
    
    def cas_to_inchikey(self, cas: str) -> Optional[str]:
        """Convert CAS number to InChIKey."""
        result = self.get_by_cas(cas)
        return result['inchikey'] if result else None
    
    def cas_to_smiles(self, cas: str) -> Optional[str]:
        """Convert CAS number to SMILES."""
        result = self.get_by_cas(cas)
        return result['smiles'] if result else None
    
    def inchikey_to_cid(self, inchikey: str) -> Optional[int]:
        """Convert InChIKey to PubChem CID."""
        result = self.get_by_inchikey(inchikey)
        return result['cid'] if result else None
    
    def inchikey_to_cas(self, inchikey: str) -> Optional[List[str]]:
        """Convert InChIKey to CAS number(s)."""
        result = self.get_by_inchikey(inchikey)
        return result['cas_numbers'] if result else None
    
    def inchi_to_cid(self, inchi: str) -> Optional[int]:
        """Convert InChI to PubChem CID."""
        result = self.get_by_inchi(inchi)
        return result['cid'] if result else None
    
    def inchi_to_cas(self, inchi: str) -> Optional[List[str]]:
        """Convert InChI to CAS number(s)."""
        result = self.get_by_inchi(inchi)
        return result['cas_numbers'] if result else None
    
    def cid_to_cas(self, cid: int) -> Optional[List[str]]:
        """Convert PubChem CID to CAS number(s)."""
        result = self.get_by_cid(cid)
        return result['cas_numbers'] if result else None
    
    def cid_to_inchikey(self, cid: int) -> Optional[str]:
        """Convert PubChem CID to InChIKey."""
        result = self.get_by_cid(cid)
        return result['inchikey'] if result else None
    
    def cid_to_inchi(self, cid: int) -> Optional[str]:
        """Convert PubChem CID to InChI."""
        result = self.get_by_cid(cid)
        return result['inchi'] if result else None
    
    # Batch conversion methods
    
    def batch_cas_to_cid(self, cas_list: List[str]) -> Dict[str, Optional[int]]:
        """
        Convert multiple CAS numbers to CIDs.
        
        Args:
            cas_list (list): List of CAS numbers
        
        Returns:
            dict: Mapping of CAS -> CID (None if not found)
        
        Example:
            >>> db = PubChemID()
            >>> results = db.batch_cas_to_cid(["50-78-2", "50-00-0"])
            >>> print(results)
            {'50-78-2': 2244, '50-00-0': 712}
        """
        results = {}
        for cas in cas_list:
            results[cas] = self.cas_to_cid(cas)
        return results
    
    def batch_cas_to_inchikey(self, cas_list: List[str]) -> Dict[str, Optional[str]]:
        """Convert multiple CAS numbers to InChIKeys."""
        results = {}
        for cas in cas_list:
            results[cas] = self.cas_to_inchikey(cas)
        return results
    
    def batch_cid_to_cas(self, cid_list: List[int]) -> Dict[int, Optional[List[str]]]:
        """Convert multiple CIDs to CAS numbers."""
        results = {}
        for cid in cid_list:
            results[cid] = self.cid_to_cas(cid)
        return results
    
    def get_by_cas_batch(self, cas_list: List[str]) -> 'pd.DataFrame':
        """
        Get complete compound information for multiple CAS numbers as a DataFrame.
        
        This method returns all available data including identifiers, chemical properties,
        and physical properties for each CAS number.
        
        Args:
            cas_list (list): List of CAS Registry Numbers
        
        Returns:
            pandas.DataFrame: DataFrame with columns for all compound properties including:
                             cid, cas, inchi, inchikey, smiles, cmpdname, iupacname, mf, mw,
                             polararea, complexity, xlogp, heavycnt, hbonddonor, hbondacc,
                             rotbonds, exactmass, charge, cidcdate
        
        Example:
            >>> db = PubChemID()
            >>> cas_list = ["50-78-2", "50-00-0", "64-17-5"]
            >>> df = db.get_by_cas_batch(cas_list)
            >>> print(df[['cas', 'cmpdname', 'mf', 'mw']])
        """
        
        rows = []
        for cas in cas_list:
            result = self.get_by_cas(cas)
            if result:
                # Create row with all properties
                row = {
                    'cid': result.get('cid'),
                    'cas': cas,
                    'inchi': result.get('inchi', ''),
                    'inchikey': result.get('inchikey', ''),
                    'smiles': result.get('smiles', ''),
                    'cmpdname': result.get('cmpdname', ''),
                    'iupacname': result.get('iupacname', ''),
                    'mf': result.get('mf', ''),
                    'mw': result.get('mw'),
                    'polararea': result.get('polararea'),
                    'complexity': result.get('complexity'),
                    'xlogp': result.get('xlogp'),
                    'heavycnt': result.get('heavycnt'),
                    'hbonddonor': result.get('hbonddonor'),
                    'hbondacc': result.get('hbondacc'),
                    'rotbonds': result.get('rotbonds'),
                    'exactmass': result.get('exactmass'),
                    'charge': result.get('charge'),
                    'cidcdate': result.get('cidcdate', '')
                }
                rows.append(row)
        
        if not rows:
            # Return empty DataFrame with correct columns
            return pd.DataFrame(columns=[
                'cid', 'cas', 'inchi', 'inchikey', 'smiles', 'cmpdname', 'iupacname',
                'mf', 'mw', 'polararea', 'complexity', 'xlogp', 'heavycnt',
                'hbonddonor', 'hbondacc', 'rotbonds', 'exactmass', 'charge', 'cidcdate'
            ])
        
        return pd.DataFrame(rows)
    
    def get_id_table_from_cas(self, cas: str) -> Optional['pd.DataFrame']:
        """
        Get identifier table for a CAS number (similar to ZeroPM format).
        
        Args:
            cas (str): CAS Registry Number
        
        Returns:
            pandas.DataFrame: Table with columns [cid, cas, inchi, inchikey, smiles, 
                             cmpdname, mf, mw] or None if not found
        
        Example:
            >>> db = PubChemID()
            >>> df = db.get_id_table_from_cas("50-78-2")
            >>> print(df)
        """
        import pandas as pd
        
        result = self.get_by_cas(cas)
        if not result:
            return None
        
        # Create DataFrame with main identifiers and properties
        df = pd.DataFrame([{
            'cid': result['cid'],
            'cas': cas,
            'inchi': result.get('inchi', ''),
            'inchikey': result.get('inchikey', ''),
            'smiles': result.get('smiles', ''),
            'cmpdname': result.get('cmpdname', ''),
            'mf': result.get('mf', ''),
            'mw': result.get('mw', None)
        }])
        
        return df
    
    def batch_get_id_table_from_cas(self, cas_list: List[str]) -> 'pd.DataFrame':
        """
        Get identifier tables for multiple CAS numbers.
        
        Args:
            cas_list (list): List of CAS Registry Numbers
        
        Returns:
            pandas.DataFrame: Combined table for all CAS numbers
        
        Example:
            >>> db = PubChemID()
            >>> df = db.batch_get_id_table_from_cas(["50-78-2", "50-00-0"])
            >>> print(df)
        """
        import pandas as pd
        
        tables = []
        for cas in cas_list:
            df = self.get_id_table_from_cas(cas)
            if df is not None:
                tables.append(df)
        
        if not tables:
            # Return empty DataFrame with correct columns
            return pd.DataFrame(columns=['cid', 'cas', 'inchi', 'inchikey', 
                                        'smiles', 'cmpdname', 'mf', 'mw'])
        
        return pd.concat(tables, ignore_index=True)
    
    def get_stats(self) -> Dict[str, int]:
        """
        Get database statistics.
        
        Returns:
            dict: Statistics about the database
        
        Example:
            >>> db = PubChemID()
            >>> stats = db.get_stats()
            >>> print(f"Total compounds: {stats['total_compounds']:,}")
        """
        cursor = self.conn.cursor()
        
        cursor.execute("SELECT COUNT(*) FROM compounds")
        total_compounds = cursor.fetchone()[0]
        
        cursor.execute("SELECT COUNT(*) FROM cas_numbers")
        total_cas = cursor.fetchone()[0]
        
        cursor.execute("SELECT COUNT(DISTINCT cid) FROM cas_numbers")
        compounds_with_cas = cursor.fetchone()[0]
        
        cursor.execute("SELECT COUNT(*) FROM synonyms")
        total_synonyms = cursor.fetchone()[0]
        
        cursor.execute("SELECT COUNT(*) FROM compounds WHERE inchikey IS NOT NULL AND inchikey != ''")
        compounds_with_inchikey = cursor.fetchone()[0]
        
        return {
            'total_compounds': total_compounds,
            'total_cas_numbers': total_cas,
            'compounds_with_cas': compounds_with_cas,
            'total_synonyms': total_synonyms,
            'compounds_with_inchikey': compounds_with_inchikey,
            'database_path': self.db_path,
            'database_size_mb': os.path.getsize(self.db_path) / (1024**2)
        }
