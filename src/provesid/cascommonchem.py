"""
CAS Common Chemistry, online: [`CASCommonChem`][provesid.cascommonchem.CASCommonChem].

CAS Common Chemistry (https://commonchemistry.cas.org) is CAS's free subset
of the CAS Registry: about half a million substances with their CAS number,
names, structure and some experimental properties. Its API (v2.0) needs a
free API key, sent in the ``X-API-KEY`` header; see
`CASCommonChem._load_api_key` for where it is looked for.

The lookups report failure in the dict they return (``found`` False and a
``status``) rather than raising, and are cached only when they succeed.

Examples:
    >>> from provesid.cascommonchem import CASCommonChem
    >>> cas = CASCommonChem(api_key="your-cas-api-key")          # doctest: +SKIP
    >>> result = cas.cas_to_detail("7732-18-5")                  # doctest: +SKIP
    >>> result["found"], result["rn"]                            # doctest: +SKIP
    (True, '7732-18-5')
"""

import os
import re
import json
import requests
import logging
from typing import Any
from rdkit import Chem, rdBase
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from .cache import cached
from .config import get_cas_api_key
from .http import HTTPClient, NotFoundError, ServiceError, ServiceTimeoutError
CASCommonChem_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

CAS_MIN_INTERVAL = 0.2
"""Seconds between two requests to CAS Common Chemistry. CAS does not publish a
per-second figure alongside the API key, so this is politeness: five
requests a second is what the rest of the package uses for a service that
has not said otherwise.
"""


CAS_SMILES_MAX_HITS = 50
"""Most search hits [`smiles_to_detail`][provesid.cascommonchem.CASCommonChem.smiles_to_detail]
fetches a record for. Every hit costs one request the first time, because a
hit carries no formula. Water's InChI finds 30, most of them water clusters.
"""


class CASCommonChemError(ServiceError):
    """
    Base exception for CAS Common Chemistry failures.

    Raised by the transport rather than by the methods below, which catch it
    and report the failure in the ``status`` key of the dict they return ---
    the contract those methods have always had.

    Examples:
        >>> issubclass(CASCommonChemNotFoundError, CASCommonChemError)
        True
    """
    pass


class CASCommonChemNotFoundError(CASCommonChemError, NotFoundError):
    """
    CAS answered, and its answer was that there is no such substance.

    Examples:
        >>> issubclass(CASCommonChemNotFoundError, NotFoundError)
        True
    """
    pass


class CASCommonChemTimeoutError(CASCommonChemError, ServiceTimeoutError):
    """
    Every attempt timed out or the connection could not be made.

    Examples:
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
    [`is_empty_result`][provesid.cache.is_empty_result] makes the same trade for every
    other client in the package.

    Args:
        result: The value returned by ``cas_to_detail`` or ``name_to_detail``.

    Returns:
        True unless the result is a dict that reports ``found`` as True.

    Examples:
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

    The endpoint comes from the Swagger file shipped in the package data
    directory. Every lookup needs the network.

    Raises:
        ValueError: At construction, when no API key can be found.

    Examples:
        >>> cas = CASCommonChem(api_key="your-cas-api-key")
        >>> cas.base_url
        'https://commonchemistry.cas.org/api'
        >>> cas.name_to_detail("aspirin")["rn"]              # doctest: +SKIP
        '50-78-2'
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

        Raises:
            ValueError: If no API key is found in any of the places
                `_load_api_key` looks.
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

        Examples:
            >>> cas = CASCommonChem.__new__(CASCommonChem)
            >>> cas.responses = {200: "Success", 404: "Invalid Request"}
            >>> cas._failure_status(CASCommonChemError("x", status_code=401))
            'Unauthorized - Check API Key'
            >>> cas._failure_status(CASCommonChemError("x", status_code=403))
            'Unauthorized - Check API Key'
            >>> cas._failure_status(CASCommonChemTimeoutError("x"))
            'Timeout'
        """
        if isinstance(exc, CASCommonChemTimeoutError):
            return "Timeout"
        # CAS answers a missing or rejected key with 403 ("API key required"),
        # checked live on 2026-09-22; 401 is kept in case it ever uses it.
        if exc.status_code in (401, 403):
            return "Unauthorized - Check API Key"
        if isinstance(exc, CASCommonChemNotFoundError):
            return "Not Found"
        if exc.status_code is None:
            # No response arrived at all, and it was not a timeout.
            return "Network Error"
        return self.responses.get(exc.status_code, "Unknown Status")

    CACHE_SCHEMA_VERSION = 2
    """Bumped whenever an entry written by an earlier version must not be
    served. Version 2 is this change: the old code cached failures, so a CAS
    number that was looked up during one network outage is on disk as a
    permanent ``{"found": False, "status": "Network Error"}``. Version 1
    entries are made unreachable rather than left to be served as fact.
    """

    def __cache_key__(self) -> tuple:
        """
        Identify this client for cache-key purposes.

        The API key is deliberately not part of it: two keys reach the same
        registry and get the same answer, so keying on it would only mean a new
        key started from an empty cache. The endpoint and the schema version are
        what change what a call returns.

        Returns:
            Tuple of the class path, the configured base URL and
            [`CACHE_SCHEMA_VERSION`][provesid.cascommonchem.CASCommonChem.CACHE_SCHEMA_VERSION].

        Examples:
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

            - found: True when CAS returned the substance
            - status: Request status: ``"Success"``, or the failure, e.g.
              ``"Not Found"``, ``"Timeout"`` or ``"Unauthorized - Check API Key"``
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
            The keys are CAS's own, so the set can grow with the API. Nothing
            is raised: every failure is reported in ``found`` and ``status``,
            and is not cached.

        Examples:
            >>> cas = CASCommonChem(api_key="your-cas-api-key")
            >>> result = cas.cas_to_detail("7732-18-5")          # doctest: +SKIP
            >>> result["found"], result["rn"]                    # doctest: +SKIP
            (True, '7732-18-5')
            >>> CASCommonChem(api_key="not-a-key").cas_to_detail("7732-18-5")["status"]  # doctest: +SKIP
            'Unauthorized - Check API Key'
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
        Returns compound details for a given name using API v2.0

        CAS's search matches a name only if CAS lists it as a synonym
        ("acetone" is found, "propan-2-one" is not). For a structure, use
        [`smiles_to_detail`][provesid.cascommonchem.CASCommonChem.smiles_to_detail].

        Args:
            name: Compound name
            timeout: Request timeout in seconds

        Returns:
            Dictionary with compound details (same as cas_to_detail), for the
            first search hit when there are several (logged). ``status`` is
            ``"Not found"`` when the search finds nothing. Two requests: a
            search, then a detail lookup.

        Examples:
            >>> cas = CASCommonChem(api_key="your-cas-api-key")
            >>> cas.name_to_detail("aspirin")["rn"]              # doctest: +SKIP
            '50-78-2'
        """
        hits, failure = self._search(name, timeout)
        if failure is not None:
            return failure
        if not hits:
            return self._not_found()

        if len(hits) > 1:
            self.logger.warning(f"Multiple compounds found for '{name}', using first result")

        return self.cas_to_detail(hits[0]["rn"])

    @cached(service='cas', skip_if=_lookup_failed)
    def smiles_to_detail(self, smiles: str, timeout=30):
        """
        Look a substance up by structure: search CAS by InChI, keep the exact match.

        CAS's search matches a SMILES only as the exact string CAS stores
        (``OCC`` finds ethanol, ``CCO`` does not), and even then its first hit
        can be an isotopologue. So the SMILES is not sent. RDKit writes its
        standard InChI, which CAS matches whatever the SMILES spelling, and
        which keeps isotopes apart. The hits are then narrowed down:

        1. The record's InChI must equal the query's. CAS gives oligomers and
           polymers the InChI of their repeat unit, so for ethanol it also
           finds "Ethanol, dimer".
        2. Records whose formula is RDKit's formula for the SMILES are
           preferred. That drops the dimer, whose formula is
           ``(C2H6O)2``. When no record passes, all of step 1's are kept, since
           CAS writes some salts differently (``C2H4O2.Na``).
        3. Of those left, the record with the most synonyms wins: sodium
           chloride (65) over rock salt (6), both ``ClNa``. When more than
           one was left, the others are logged at WARNING.

        Each hit costs one detail request the first time it is seen, up to
        ``CAS_SMILES_MAX_HITS``. The detail answers are cached, and so is the
        result. A convenience layer over CAS's ``/search`` and ``/detail``
        calls.

        Args:
            smiles: SMILES string.
            timeout: Request timeout in seconds, per request.

        Returns:
            Dictionary with compound details, as
            [`cas_to_detail`][provesid.cascommonchem.CASCommonChem.cas_to_detail]
            returns them. ``status`` is ``"Invalid SMILES"`` when RDKit cannot
            read the SMILES (no request is made), and ``"Not found"`` when no
            record has its InChI.

        Examples:
            >>> cas = CASCommonChem(api_key="your-cas-api-key")
            >>> cas.smiles_to_detail("CCO")["rn"]                # doctest: +SKIP
            '64-17-5'
            >>> cas.smiles_to_detail("C1=CC=CC=C1")["name"]      # doctest: +SKIP
            'Benzene'
            >>> cas.smiles_to_detail("not-a-smiles")["status"]
            'Invalid SMILES'
        """
        with rdBase.BlockLogs():
            mol = Chem.MolFromSmiles(smiles) if smiles and smiles.strip() else None
            inchi = Chem.MolToInchi(mol) if mol is not None else ""
        if not inchi:
            res = self._empty_res()
            res["status"] = "Invalid SMILES"
            res["found"] = False
            return res

        hits, failure = self._search(inchi, timeout)
        if failure is not None:
            return failure
        if len(hits) > CAS_SMILES_MAX_HITS:
            self.logger.warning(
                f"{len(hits)} CAS hits for '{smiles}'; examining the first {CAS_SMILES_MAX_HITS}")

        records = [self.cas_to_detail(hit["rn"], timeout)
                   for hit in hits[:CAS_SMILES_MAX_HITS]]
        same_inchi = [r for r in records if r["found"] and r.get("inchi") == inchi]
        formula = CalcMolFormula(mol)
        same_formula = [r for r in same_inchi
                        if self._plain_formula(r.get("molecularFormula")) == formula]
        candidates = same_formula or same_inchi
        if not candidates:
            self.logger.debug(f"No CAS record has the InChI of '{smiles}': {inchi}")
            return self._not_found()

        best = max(candidates, key=lambda r: len(r.get("synonyms") or []))
        if len(candidates) > 1:
            others = ", ".join(r["rn"] for r in candidates if r is not best)
            self.logger.warning(
                f"Several CAS records match '{smiles}'; using {best['rn']}, not {others}")
        return best

    def _search(self, query: str, timeout=30):
        """
        Run CAS's ``/search`` and return ``(hits, failure)``.

        ``hits`` is the list of ``{"rn", "name", "images"}`` dicts CAS
        returned (empty when it found nothing), and ``failure`` is None. When
        the request fails, ``hits`` is None and ``failure`` is the result dict
        to hand back, with its ``status``, as
        [`cas_to_detail`][provesid.cascommonchem.CASCommonChem.cas_to_detail]
        reports failures.
        """
        url = self.base_url + self.query_url[2] + "?q=" + requests.utils.quote(query)
        try:
            data = self._http.get_json(url, timeout=timeout)
        except CASCommonChemError as e:
            res = self._empty_res()
            res["status"] = self._failure_status(e)
            res["found"] = False
            if res["status"] == "Unauthorized - Check API Key":
                self.logger.error("CAS API authentication failed. Check your API key.")
            else:
                self.logger.warning(f"CAS search failed for '{query}': {e}")
            return None, res
        except Exception as e:
            res = self._empty_res()
            res["status"] = "Error"
            res["found"] = False
            self.logger.error(f"Unexpected error searching CAS for '{query}': {e}")
            return None, res
        return (data.get("results") or [] if data.get("count") else []), None

    def _not_found(self):
        """The result dict for a search that CAS answered with no match."""
        res = self._empty_res()
        res["status"] = "Not found"
        res["found"] = False
        return res

    @staticmethod
    def _plain_formula(formula) -> str:
        """
        CAS's formula without its HTML markup.

        Examples:
            >>> CASCommonChem._plain_formula("(C<sub>2</sub>H<sub>6</sub>O)<sub>2</sub>")
            '(C2H6O)2'
        """
        return re.sub(r"<[^>]+>", "", formula or "")

    def clear_cache(self):
        """
        Delete every cached CAS Common Chemistry answer, in memory and on disk.

        The same as ``provesid.clear_cache(service='cas')``.

        Examples:
            >>> cas = CASCommonChem(api_key="your-cas-api-key")
            >>> cas.clear_cache()
            >>> cas.get_cache_info()['file_count']
            0
        """
        from .cache import clear_cache
        clear_cache(service='cas')

    def get_cache_info(self):
        """
        Size and location of the CAS Common Chemistry cache.

        Returns:
            The statistics
            [`provesid.cache.get_cache_info`][provesid.cache.get_cache_info]
            reports for
                ``service='cas'``.

        Examples:
            >>> CASCommonChem(api_key="your-cas-api-key").get_cache_info()['cache_directory'].endswith('cas')
            True
        """
        from .cache import get_cache_info
        return get_cache_info(service='cas')

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



