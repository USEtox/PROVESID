import os
import json
import requests
import logging
from functools import lru_cache
from .cache import cached
from .config import get_cas_api_key
CASCommonChem_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

class CASCommonChem:
    """
    calling the CAS Common Chemistry API v2.0 to get the information for a given CAS RN
    Requires API key authentication via X-API-KEY header
    """
    def __init__(self, swagger_file_name='commonchemistry-swagger.json', use_cache: bool = True, 
                 api_key: str = None, api_key_file: str = None):
        """
        Initialize CAS Common Chemistry API client
        
        Args:
            swagger_file_name: Name of the swagger JSON file
            use_cache: Whether to use cache for lookups
            api_key: Direct API key string (takes precedence over api_key_file)
            api_key_file: Path to file containing API key
        """
        self.data_folder = CASCommonChem_path
        self.swagger_file_path = os.path.join(self.data_folder, swagger_file_name)
        
        # Load the swagger file
        with open(self.swagger_file_path, 'r') as f:
            self.swagger = json.load(f)
        
        # Set up API configuration
        host = self.swagger["host"]
        schemes = self.swagger["schemes"][0]
        base_path = self.swagger["basePath"]
        self.base_url = f"{schemes}://{host}{base_path}"
        self.query_url = ["/detail", "/export", "/search"]
        self.responses = {200: "Success", 400: "Invalid Request", 404: "Invalid Request", 500: "Internal Server Error"}
        self.use_cache = use_cache
        
        # Handle API key - check multiple sources in priority order
        self.api_key = self._load_api_key(api_key, api_key_file)
        if not self.api_key:
            raise ValueError(
                "API key is required for CAS Common Chemistry API v2.0.\n"
                "Options:\n"
                "1. Provide api_key parameter: CASCommonChem(api_key='your-key')\n"
                "2. Provide api_key_file parameter: CASCommonChem(api_key_file='path/to/key.txt')\n"
                "3. Set persistent API key: from provesid.config import set_cas_api_key; set_cas_api_key('your-key')\n"
                "4. Set environment variable: CCC_API_KEY or CAS_API_KEY"
            )
    
    def _load_api_key(self, api_key: str = None, api_key_file: str = None) -> str:
        """
        Load API key from multiple sources in priority order:
        1. Direct api_key parameter
        2. API key file (api_key_file parameter)
        3. Persistent configuration (set via set_cas_api_key())
        4. Environment variables (CCC_API_KEY, CAS_API_KEY)
        """
        # Priority 1: Direct parameter
        if api_key:
            return api_key.strip()
        
        # Priority 2: API key file
        if api_key_file:
            try:
                with open(api_key_file, 'r', encoding='utf-8') as f:
                    key = f.read().strip()
                    if key:
                        return key
                    else:
                        logging.warning(f"API key file {api_key_file} is empty")
            except FileNotFoundError:
                logging.warning(f"API key file not found: {api_key_file}")
            except Exception as e:
                logging.warning(f"Error reading API key file {api_key_file}: {e}")
        
        # Priority 3: Persistent configuration
        try:
            config_key = get_cas_api_key()
            if config_key:
                return config_key.strip()
        except Exception as e:
            logging.debug(f"Could not load API key from config: {e}")
        
        # Priority 4: Environment variables
        env_key = os.environ.get('CCC_API_KEY') or os.environ.get('CAS_API_KEY')
        if env_key:
            return env_key.strip()
        
        return None
    
    def _get_headers(self) -> dict:
        """Get headers with API key for authentication"""
        return {
            'X-API-KEY': self.api_key,
            'Accept': 'application/json',
            'Content-Type': 'application/json'
        }
    
    @cached(service='cas')
    def cas_to_detail(self, cas_rn: str, timeout=30):
        """
        Returns a dictionary with the data for a given CAS RN using API v2.0
        
        Args:
            cas_rn: CAS Registry Number (with or without hyphens)
            timeout: Request timeout in seconds
            
        Returns:
            Dictionary with compound details including:
            - status: Request status
            - canonicalSmile: Canonical SMILES string
            - experimentalProperties: List of experimental properties
            - hasMolfile: Boolean indicating molfile availability
            - images: List of compound images
            - inchi: InChI string
            - inchiKey: InChI Key
            - molecularFormula: Molecular formula
            - molecularMass: Molecular mass
            - name: Compound name
            - synonyms: List of synonyms
            - uri: Compound URI
        """
        url = self.base_url + self.query_url[0] + "?cas_rn=" + cas_rn
        res = self._empty_res()
        
        try:
            response = requests.get(url, headers=self._get_headers(), timeout=timeout)
            res["status"] = self.responses.get(response.status_code, "Unknown Status")
            
            if response.status_code == 200:
                data = response.json()
                for key in data.keys():
                    res[key] = data[key]
                res["found"] = True
            else:
                res["found"] = False
                if response.status_code == 404:
                    res["status"] = "Not Found"
                elif response.status_code == 401:
                    res["status"] = "Unauthorized - Check API Key"
                    logging.error("CAS API authentication failed. Check your API key.")
                
        except requests.exceptions.Timeout:
            res["status"] = "Timeout"
            logging.error(f"Request timeout for CAS RN: {cas_rn}")
        except requests.exceptions.RequestException as e:
            res["status"] = "Network Error"
            logging.error(f"Network error for CAS RN {cas_rn}: {e}")
        except Exception as e:
            res["status"] = "Error"
            logging.error(f"Unexpected error for CAS RN {cas_rn}: {e}")
            
        return res
    
    @cached(service='cas')
    def name_to_detail(self, name: str, timeout=30):
        """
        Returns compound details for a given name or SMILES using API v2.0
        
        Args:
            name: Compound name or SMILES string
            timeout: Request timeout in seconds
            
        Returns:
            Dictionary with compound details (same as cas_to_detail)
        """
        res = self._empty_res()
        url = self.base_url + self.query_url[2] + "?q=" + requests.utils.quote(name)
        
        try:
            response = requests.get(url, headers=self._get_headers(), timeout=timeout)
            
            if response.status_code == 401:
                res["status"] = "Unauthorized - Check API Key"
                res["found"] = False
                logging.error("CAS API authentication failed. Check your API key.")
                return res
            elif response.status_code != 200:
                res["status"] = self.responses.get(response.status_code, "Unknown Status")
                res["found"] = False
                return res
                
            res_call = response.json()
            
            if "count" not in res_call or res_call["count"] == 0:
                res["status"] = "Not found"
                res["found"] = False
                return res
                
            if res_call["count"] > 1:
                logging.warning(f"Multiple compounds found for '{name}', using first result")
                
            # Get CAS RN from first result and fetch details
            cas_rn = res_call["results"][0]["rn"]
            return self.cas_to_detail(cas_rn)
            
        except requests.exceptions.Timeout:
            res["status"] = "Timeout"
            res["found"] = False
            logging.error(f"Request timeout for name: {name}")
        except requests.exceptions.RequestException as e:
            res["status"] = "Network Error"
            res["found"] = False
            logging.error(f"Network error for name '{name}': {e}")
        except Exception as e:
            res["status"] = "Error"
            res["found"] = False
            logging.error(f"Unexpected error for name '{name}': {e}")
            
        return res
    
    def smiles_to_detail(self, smiles: str, timeout=30):
        return self.name_to_detail(smiles, timeout)
    
    def clear_cache(self):
        """Clear all cached results for CAS Common Chemistry"""
        from .cache import clear_cas_cache
        clear_cas_cache()
    
    def get_cache_info(self):
        """Get cache information for CAS Common Chemistry cached methods"""
        from .cache import get_cas_cache_info
        return get_cas_cache_info()
    
    @staticmethod
    def _empty_res():
        return {
                "cas_rn": "",
                "status": "",
                "canonicalSmile": "",
                "experimentalProperties": [],
                "hasMolfile": False,
                "images": [],
                "inchi": "",
                "inchiKey": "",
                "molecularFormula": "",
                "molecularMass": "",
                "name": "",
                "propertyCitations": [],
                "replacedRns": [],
                "rn": "",
                "smile": "",
                "synonyms": [],
                "uri": ""
            }
        
        
        
