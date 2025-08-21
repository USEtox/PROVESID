# a pubchem package with a limited number of functionalities of pubchempy 
# but with a simpler interface that serves our purpose in PROVES

import requests
import time
import json
import logging
from typing import Dict, List, Union, Optional, Any
from urllib.parse import quote

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
    SMILES = "ConnectivitySMILES"
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
    
    def __init__(self, base_url: str = pugrest_prolog, pause_time: float = pause_between_calls):
        """
        Initialize PubChem API client
        
        Args:
            base_url: Base URL for PubChem REST API
            pause_time: Minimum time between API calls in seconds
        """
        self.base_url = base_url.rstrip('/')
        self.pause_time = pause_time
        self.last_request_time = 0
        
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
    
    # Compound methods
    def get_compound_by_cid(self, cid: Union[int, str], output_format: str = OutputFormat.JSON) -> Any:
        """
        Get compound record by CID
        
        Args:
            cid: Compound ID
            output_format: Desired output format
            
        Returns:
            Compound data
        """
        url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.CID, cid, 
                             Operation.RECORD, output_format)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
    def get_compounds_by_name(self, name: str, output_format: str = OutputFormat.JSON,
                             name_type: str = "word") -> Any:
        """
        Get compounds by name
        
        Args:
            name: Compound name
            output_format: Desired output format
            name_type: Name search type ("word" or "complete")
            
        Returns:
            Compound data
        """
        url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.NAME, name,
                             Operation.RECORD, output_format, name_type=name_type)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
    def get_compounds_by_smiles(self, smiles: str, output_format: str = OutputFormat.JSON) -> Any:
        """
        Get compounds by SMILES
        
        Args:
            smiles: SMILES string
            output_format: Desired output format
            
        Returns:
            Compound data
        """
        url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.SMILES, smiles,
                             Operation.RECORD, output_format)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
    def get_compounds_by_inchikey(self, inchikey: str, output_format: str = OutputFormat.JSON) -> Any:
        """
        Get compounds by InChIKey
        
        Args:
            inchikey: InChI Key
            output_format: Desired output format
            
        Returns:
            Compound data
        """
        url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.INCHIKEY, inchikey,
                             Operation.RECORD, output_format)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
    def get_compound_properties(self, cids: Union[int, str, List[Union[int, str]]], 
                               properties: List[str], 
                               output_format: str = OutputFormat.JSON) -> Any:
        """
        Get compound properties by CID(s)
        
        Args:
            cids: Single CID or list of CIDs
            properties: List of property names
            output_format: Desired output format
            
        Returns:
            Property data
        """
        props_str = ','.join(properties)
        operation = f"{Operation.PROPERTY}/{props_str}"
        url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.CID, cids,
                             operation, output_format)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
    def get_compound_synonyms(self, cid: Union[int, str], output_format: str = OutputFormat.JSON) -> Any:
        """
        Get compound synonyms by CID
        
        Args:
            cid: Compound ID
            output_format: Desired output format
            
        Returns:
            Synonyms data
        """
        url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.CID, cid,
                             Operation.SYNONYMS, output_format)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
    def get_cids_by_name(self, name: str, output_format: str = OutputFormat.JSON,
                        name_type: str = "word") -> Any:
        """
        Get CIDs by compound name
        
        Args:
            name: Compound name
            output_format: Desired output format
            name_type: Name search type ("word" or "complete")
            
        Returns:
            CID list
        """
        url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.NAME, name,
                             Operation.CIDS, output_format, name_type=name_type)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
    def get_cids_by_smiles(self, smiles: str, output_format: str = OutputFormat.JSON) -> Any:
        """
        Get CIDs by SMILES
        
        Args:
            smiles: SMILES string
            output_format: Desired output format
            
        Returns:
            CID list
        """
        url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.SMILES, smiles,
                             Operation.CIDS, output_format)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
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
    def substructure_search(self, query: str, query_type: str = "smiles", 
                           output_format: str = OutputFormat.JSON, **options) -> Any:
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
        search_type = f"{FastSearch.FASTSUBSTRUCTURE}/{query_type}"
        url = self._build_url(Domain.COMPOUND, search_type, query,
                             Operation.CIDS, output_format, **options)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
    def superstructure_search(self, query: str, query_type: str = "smiles",
                             output_format: str = OutputFormat.JSON, **options) -> Any:
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
        search_type = f"{FastSearch.FASTSUPERSTRUCTURE}/{query_type}"
        url = self._build_url(Domain.COMPOUND, search_type, query,
                             Operation.CIDS, output_format, **options)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
    def similarity_search(self, query: str, query_type: str = "smiles",
                         threshold: int = 90, output_format: str = OutputFormat.JSON,
                         **options) -> Any:
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
        search_type = f"{FastSearch.FASTSIMILARITY_2D}/{query_type}"
        url = self._build_url(Domain.COMPOUND, search_type, query,
                             Operation.CIDS, output_format, 
                             Threshold=threshold, **options)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
    def identity_search(self, query: str, query_type: str = "smiles",
                       identity_type: str = "same_stereo_isotope",
                       output_format: str = OutputFormat.JSON, **options) -> Any:
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
        search_type = f"{FastSearch.FASTIDENTITY}/{query_type}"
        url = self._build_url(Domain.COMPOUND, search_type, query,
                             Operation.CIDS, output_format,
                             identity_type=identity_type, **options)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
    # Substance methods
    def get_substance_by_sid(self, sid: Union[int, str], output_format: str = OutputFormat.JSON) -> Any:
        """
        Get substance by SID
        
        Args:
            sid: Substance ID
            output_format: Desired output format
            
        Returns:
            Substance data
        """
        url = self._build_url(Domain.SUBSTANCE, SubstanceDomainNamespace.SID, sid,
                             Operation.RECORD, output_format)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
    def get_substances_by_name(self, name: str, output_format: str = OutputFormat.JSON) -> Any:
        """
        Get substances by name
        
        Args:
            name: Substance name
            output_format: Desired output format
            
        Returns:
            Substance data
        """
        url = self._build_url(Domain.SUBSTANCE, SubstanceDomainNamespace.NAME, name,
                             Operation.RECORD, output_format)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
    def get_sids_by_name(self, name: str, output_format: str = OutputFormat.JSON,
                        sourcename: Optional[str] = None) -> Any:
        """
        Get SIDs by name
        
        Args:
            name: Substance name
            output_format: Desired output format
            sourcename: Restrict to specific source
            
        Returns:
            SID list
        """
        options = {}
        if sourcename:
            options['sourcename'] = sourcename
            
        url = self._build_url(Domain.SUBSTANCE, SubstanceDomainNamespace.NAME, name,
                             Operation.SIDS, output_format, **options)
        response = self._make_request(url)
        return self._parse_response(response, output_format)
    
    # Assay methods
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
    
    def get_basic_compound_info(self, cid: Union[int, str]) -> Dict[str, Any]:
        """
        Get basic compound information including properties and synonyms
        
        Args:
            cid: Compound ID
            
        Returns:
            Dictionary with compound information
        """
        try:
            # Get basic properties
            properties = [
                CompoundProperties.MOLECULAR_FORMULA,
                CompoundProperties.MOLECULAR_WEIGHT,
                CompoundProperties.SMILES,
                CompoundProperties.INCHI,
                CompoundProperties.INCHIKEY,
                CompoundProperties.IUPAC_NAME
            ]
            
            prop_data = self.get_compound_properties(cid, properties)
            synonyms_data = self.get_compound_synonyms(cid)
            
            return {
                "success": True,
                "cid": cid,
                "properties": prop_data,
                "synonyms": synonyms_data,
                "error": None
            }
            
        except Exception as e:
            return {
                "success": False,
                "cid": cid,
                "properties": None,
                "synonyms": None,
                "error": str(e)
            }

    def get_all_compound_info(self, cid: Union[int, str]) -> Dict[str, Any]:
        """
        Get all compound properties as listed in CompoundProperties

        Args:
            cid: Compound ID

        Returns:
            Dictionary with compound information
        """
        try:
            properties = [prop for prop in dir(CompoundProperties) if not prop.startswith("_")]
            prop_data = self.get_compound_properties(cid, properties)

            return {
                "success": True,
                "cid": cid,
                "properties": prop_data,
                "error": None
            }

        except Exception as e:
            return {
                "success": False,
                "cid": cid,
                "properties": None,
                "error": str(e)
            }
