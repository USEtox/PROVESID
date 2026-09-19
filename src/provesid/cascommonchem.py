import os
import json
import requests
import logging
from typing import Any
from .cache import cached
from .config import get_cas_api_key
from .http import HTTPClient, NotFoundError, ServiceError, ServiceTimeoutError
CASCommonChem_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

#: Seconds between two requests to CAS Common Chemistry. CAS does not publish a
#: per-second figure alongside the API key, so this is politeness: five
#: requests a second is what the rest of the package uses for a service that
#: has not said otherwise.
CAS_MIN_INTERVAL = 0.2


class CASCommonChemError(ServiceError):
    """
    Base exception for CAS Common Chemistry failures.

    Raised by the transport rather than by the methods below, which catch it
    and report the failure in the ``status`` key of the dict they return ---
    the contract those methods have always had.

    Example:
        >>> issubclass(CASCommonChemNotFoundError, CASCommonChemError)
        True
    """
    pass


class CASCommonChemNotFoundError(CASCommonChemError, NotFoundError):
    """
    CAS answered, and its answer was that there is no such substance.

    Example:
        >>> issubclass(CASCommonChemNotFoundError, NotFoundError)
        True
    """
    pass


class CASCommonChemTimeoutError(CASCommonChemError, ServiceTimeoutError):
    """
    Every attempt timed out or the connection could not be made.

    Example:
        >>> issubclass(CASCommonChemTimeoutError, ServiceTimeoutError)
        True
    """
    pass


def _lookup_failed(result: Any) -> bool:
    """
    Report whether a CAS lookup did not come back with a substance.

    Both methods in this module report absence *and* failure the same way, by
    returning a dict whose ``found`` is False --- so neither can be cached.
    Caching the failure would turn one throttled request into a permanent
    "no such CAS number"; caching the absence would do the same thing for a
    substance CAS adds next month. Re-asking costs one cheap request, and
    :func:`~provesid.cache.is_empty_result` makes the same trade for every
    other client in the package.

    Args:
        result: The value returned by ``cas_to_detail`` or ``name_to_detail``.

    Returns:
        True unless the result is a dict that reports ``found`` as True.

    Example:
        >>> _lookup_failed({"found": False, "status": "Timeout"})
        True
        >>> _lookup_failed({"found": True, "rn": "7732-18-5"})
        False
    """
    return not (isinstance(result, dict) and result.get("found") is True)


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

        self.logger = logging.getLogger(__name__)

        # One shared transport. Before this the module had no pacing and no
        # retry at all: a throttled request became a "Network Error" entry in
        # the returned dict, and the cache remembered it. CAS uses its status
        # codes honestly --- 404 for an unknown CAS number, 401 for a rejected
        # key --- so the default classifier reads them correctly, and the key
        # goes on the client where every request picks it up.
        self._http = HTTPClient(
            min_interval=CAS_MIN_INTERVAL,
            headers=self._get_headers(),
            error_cls=CASCommonChemError,
            not_found_cls=CASCommonChemNotFoundError,
            timeout_cls=CASCommonChemTimeoutError,
            pace_host=self.base_url,
            logger=self.logger,
        )

    def _failure_status(self, exc: CASCommonChemError) -> str:
        """
        Name the failure an exception describes, the way this module always has.

        Args:
            exc: The exception the transport raised.

        Returns:
            The string to put in the ``status`` key of the returned dict.

        Example:
            >>> cas = CASCommonChem.__new__(CASCommonChem)
            >>> cas.responses = {200: "Success", 404: "Invalid Request"}
            >>> cas._failure_status(CASCommonChemError("x", status_code=401))
            'Unauthorized - Check API Key'
            >>> cas._failure_status(CASCommonChemTimeoutError("x"))
            'Timeout'
        """
        if isinstance(exc, CASCommonChemTimeoutError):
            return "Timeout"
        if exc.status_code == 401:
            return "Unauthorized - Check API Key"
        if isinstance(exc, CASCommonChemNotFoundError):
            return "Not Found"
        if exc.status_code is None:
            # No response arrived at all, and it was not a timeout.
            return "Network Error"
        return self.responses.get(exc.status_code, "Unknown Status")

    #: Bumped whenever an entry written by an earlier version must not be
    #: served. Version 2 is this change: the old code cached failures, so a CAS
    #: number that was looked up during one network outage is on disk as a
    #: permanent ``{"found": False, "status": "Network Error"}``. Version 1
    #: entries are made unreachable rather than left to be served as fact.
    CACHE_SCHEMA_VERSION = 2

    def __cache_key__(self) -> tuple:
        """
        Identify this client for cache-key purposes.

        The API key is deliberately not part of it: two keys reach the same
        registry and get the same answer, so keying on it would only mean a new
        key started from an empty cache. The endpoint and the schema version are
        what change what a call returns.

        Returns:
            Tuple of the class path, the configured base URL and
            :attr:`CACHE_SCHEMA_VERSION`.

        Example:
            >>> cas = CASCommonChem.__new__(CASCommonChem)
            >>> cas.base_url = "https://commonchemistry.cas.org/api"
            >>> cas.__cache_key__()
            ('provesid.cascommonchem.CASCommonChem', 'https://commonchemistry.cas.org/api', 2)
        """
        return ("provesid.cascommonchem.CASCommonChem", self.base_url,
                self.CACHE_SCHEMA_VERSION)

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
    
    @cached(service='cas', skip_if=_lookup_failed)
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
            data = self._http.get_json(url, timeout=timeout)
        except CASCommonChemError as e:
            res["status"] = self._failure_status(e)
            res["found"] = False
            if res["status"] == "Unauthorized - Check API Key":
                self.logger.error("CAS API authentication failed. Check your API key.")
            else:
                self.logger.warning(f"CAS lookup failed for CAS RN {cas_rn}: {e}")
            return res
        except Exception as e:
            res["status"] = "Error"
            res["found"] = False
            self.logger.error(f"Unexpected error for CAS RN {cas_rn}: {e}")
            return res

        res["status"] = self.responses[200]
        for key in data.keys():
            res[key] = data[key]
        res["found"] = True
        return res
    
    @cached(service='cas', skip_if=_lookup_failed)
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
            res_call = self._http.get_json(url, timeout=timeout)
        except CASCommonChemError as e:
            res["status"] = self._failure_status(e)
            res["found"] = False
            if res["status"] == "Unauthorized - Check API Key":
                self.logger.error("CAS API authentication failed. Check your API key.")
            else:
                self.logger.warning(f"CAS search failed for name '{name}': {e}")
            return res
        except Exception as e:
            res["status"] = "Error"
            res["found"] = False
            self.logger.error(f"Unexpected error for name '{name}': {e}")
            return res

        if not res_call.get("count"):
            res["status"] = "Not found"
            res["found"] = False
            return res

        if res_call["count"] > 1:
            self.logger.warning(f"Multiple compounds found for '{name}', using first result")

        # Get CAS RN from first result and fetch details
        cas_rn = res_call["results"][0]["rn"]
        return self.cas_to_detail(cas_rn)
    
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
        
        
        
