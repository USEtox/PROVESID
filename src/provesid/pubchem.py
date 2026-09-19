# a pubchem package with a limited number of functionalities of pubchempy 
# but with a simpler interface that serves our purpose in PROVES

import requests
import json
import logging
import re
import os
import pandas as pd
from typing import Dict, List, Union, Optional, Any
from urllib.parse import quote
from .cache import cached, is_empty_result
from .http import (
    HTTPClient,
    Outcome,
    ServiceError,
    NotFoundError,
    ServiceTimeoutError,
)
from .utils import user_dataset_path

pugrest_prolog = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
pause_between_calls = 0.2 # seconds

#: Longest identifier list PROVESID will put in a URL path before switching to
#: POST. PUG-REST documents a ceiling of about 2000 characters for the whole
#: URL; the margin left here covers the prolog, the operation and the property
#: list that share the path with the identifiers.
URL_IDENTIFIER_LIMIT = 1600

#: Longest total time a PubChem client will spend waiting between retries, in
#: seconds: "do not make the caller wait longer than this".
#:
#: Ten seconds leaves the whole cheap back-off curve intact --- 1 + 2 + 4 for a
#: transient 500 or a timeout, which is where retrying earns its keep --- while
#: declining to sit out PubChem's ``Retry-After: 30``. That is deliberate.
#: PubChem sends that header when it has throttled or blacklisted an IP, and a
#: block like that does not lift in thirty seconds: measured on 2026-09-19,
#: waiting the full thirty and asking again returned the same 503. So the wait
#: buys nothing while every call pays it, which for a caller resolving a
#: thousand names is hours instead of an immediate "you are blocked".
#:
#: A caller who does want to wait a throttle out raises it
#: (``api._http.max_elapsed = 180``), which is the right setting for an
#: unattended bulk job.
RETRY_WAIT_BUDGET = 10.0

#: Default number of CIDs per bulk property request. PubChem answers several
#: hundred at a time without complaint, but a smaller chunk costs less to redo
#: when one request has to be retried.
PROPERTY_CHUNK_SIZE = 200

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

class PubChemError(ServiceError):
    """Custom exception for PubChem API errors"""
    pass

class PubChemTimeoutError(PubChemError, ServiceTimeoutError):
    """Exception raised when request times out"""
    pass

class PubChemNotFoundError(PubChemError, NotFoundError):
    """Exception raised when resource is not found"""
    pass

class PubChemServerError(PubChemError):
    """Exception raised when server error occurs"""
    pass


#: Fault codes PubChem returns when it is shedding load rather than reporting
#: absence. PubChem does not always pair them with a 5xx status, so the code in
#: the body — not the HTTP status alone — decides whether the answer is "no such
#: data" or "ask again later". Both services' codes live here because both
#: :mod:`provesid.pubchem` and :mod:`provesid.pubchemview` read this one set;
#: PUG-REST never emits a ``PUGVIEW.*`` code, so carrying them is inert there.
TRANSIENT_FAULT_CODES = frozenset({
    "PUGREST.ServerBusy",
    "PUGREST.ServerError",
    "PUGREST.Timeout",
    "PUGVIEW.ServerBusy",
    "PUGVIEW.ServerError",
    "PUGVIEW.Timeout",
})


#: Fault codes that report genuine absence: no such compound, no such heading.
#: Permanent answers, so a request carrying one is never retried.
ABSENCE_FAULT_CODES = frozenset({
    "PUGREST.NotFound",
    "PUGVIEW.NotFound",
    "PUGVIEW.BadRequest",
})


def fault_code(response: requests.Response) -> Optional[str]:
    """
    Read PubChem's ``Fault.Code`` out of an error response.

    PUG-REST describes every error in the body, for example
    ``{"Fault": {"Code": "PUGREST.NotFound", "Message": "No synonyms found ..."}}``.
    That code is the only reliable way to tell a compound that genuinely has no
    such data from a service that is momentarily refusing work.

    Args:
        response: The HTTP response to inspect.

    Returns:
        The fault code, or None when the body is not a PubChem fault.
    """
    try:
        return response.json().get("Fault", {}).get("Code")
    except Exception:
        return None


def _classify_fault(response: requests.Response, bare_400: Outcome) -> Outcome:
    """
    Classify a PubChem response, reading the fault code before the status.

    Shared by :func:`pugrest_classify` and :func:`pugview_classify`, which
    differ in one place only --- see ``bare_400``.

    Args:
        response: The response to classify.
        bare_400: What a 400 that carries no fault code means for this service.

    Returns:
        The :class:`~provesid.http.Outcome` for this response.
    """
    status = response.status_code
    if 200 <= status < 300:
        return Outcome.OK

    code = fault_code(response)
    if code in TRANSIENT_FAULT_CODES:
        return Outcome.RETRY
    if code in ABSENCE_FAULT_CODES:
        return Outcome.ABSENT
    if status == 429 or status >= 500:
        return Outcome.RETRY
    if status == 404:
        # An unknown compound, on either service. Permanent even with no fault.
        return Outcome.ABSENT
    if status == 400:
        return bare_400
    if 400 <= status < 500:
        return Outcome.FATAL
    return Outcome.RETRY


def pugrest_classify(response: requests.Response) -> Outcome:
    """
    Decide what a PUG-REST response means, reading the fault code first.

    PubChem's status codes are not a reliable guide on their own: under load it
    answers with ``PUGREST.ServerBusy`` behind a 404 as readily as behind a
    503, and treating that as absence caches a momentary outage as fact (see
    §17.6 of the refactor plan). The fault code in the body is what the service
    always sends, so it is what decides.

    A 400 is where this differs from :func:`pugview_classify`: PUG-REST takes
    its query in the URL path, so a 400 means the path was wrong --- an
    unknown property name, a malformed identifier. That is a permanent error in
    the request, not a statement that the compound has no such data, and the
    caller needs the body to learn which property it misspelled.

    Args:
        response: The response to classify.

    Returns:
        :attr:`~provesid.http.Outcome.OK` for a 2xx, ``RETRY`` for any
        transient fault code and for 429/5xx, ``ABSENT`` for an absence fault
        code and for a bare 404, ``FATAL`` for a 400 and any other client
        error.

    Example:
        >>> class R:
        ...     status_code = 404
        ...     def json(self): return {"Fault": {"Code": "PUGREST.ServerBusy"}}
        >>> pugrest_classify(R()).name
        'RETRY'
        >>> class R:
        ...     status_code = 400
        ...     def json(self): return {"Fault": {"Code": "PUGREST.BadRequest"}}
        >>> pugrest_classify(R()).name
        'FATAL'
    """
    return _classify_fault(response, Outcome.FATAL)


def pugview_classify(response: requests.Response) -> Outcome:
    """
    Decide what a PUG-View response means, reading the fault code first.

    Identical to :func:`pugrest_classify` but for a bare 400, which PUG-View
    uses to say "no such heading" --- the heading is a query parameter, and
    asking for one a compound does not have is absence, not a bad request.
    ``PUGVIEW.BadRequest`` is already in :data:`ABSENCE_FAULT_CODES` for the
    same reason; this covers the case where the body carries no fault at all.

    Args:
        response: The response to classify.

    Returns:
        As :func:`pugrest_classify`, except that a bare 400 is ``ABSENT``.

    Example:
        >>> class R:
        ...     status_code = 400
        ...     def json(self): return None
        >>> pugview_classify(R()).name
        'ABSENT'
    """
    return _classify_fault(response, Outcome.ABSENT)


def _synonyms_incomplete(result: Any) -> bool:
    """
    Report whether a property result carries a failed synonym fetch.

    ``get_compound_properties`` returns the properties it did retrieve even when
    the follow-up synonym request failed, so the dict looks successful while its
    ``synonyms`` entry is missing. Caching that would make the gap permanent, so
    it is passed to :func:`cached` as a ``skip_if`` predicate.

    Args:
        result: Return value of ``get_compound_properties``.

    Returns:
        True when the result records a synonym fetch error.
    """
    return isinstance(result, dict) and result.get('synonyms_error') is not None


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
        self.use_cache = use_cache
        self.logger = logging.getLogger(__name__)

        # One shared transport. PUG-REST describes every error in the body, so
        # it supplies its own classifier rather than trusting the status code:
        # see :func:`pugrest_classify`. ``pace_host`` is what keeps this client
        # and any PubChemView in the same process inside PubChem's five
        # requests per second, which is a per-IP budget and not a per-object
        # one.
        self._http = HTTPClient(
            min_interval=pause_time,
            max_elapsed=RETRY_WAIT_BUDGET,
            classify=pugrest_classify,
            error_cls=PubChemError,
            not_found_cls=PubChemNotFoundError,
            timeout_cls=PubChemTimeoutError,
            rate_limit_cls=PubChemServerError,
            retry_exhausted_cls=PubChemServerError,
            pace_host=self.base_url,
            logger=self.logger,
        )

    @property
    def pause_time(self) -> float:
        """
        Minimum seconds between two requests to PubChem.

        Held by the transport, and settable: raising it slows a long batch
        down, and the new value takes effect on the next request, retries
        included.

        Example:
            >>> api = PubChemAPI()
            >>> api.pause_time
            0.2
            >>> api.pause_time = 1.0      # gentler, for a long batch
        """
        return self._http.min_interval

    @pause_time.setter
    def pause_time(self, seconds: float) -> None:
        self._http.min_interval = seconds

    @property
    def last_request_time(self) -> float:
        """
        When this client last made a request, as a Unix timestamp.

        Kept on the transport; exposed here because it describes this client's
        own pacing. The clock the pacing is *measured* against is shared with
        every other client aimed at PubChem --- see
        :class:`provesid.http.RateLimiter`.

        Returns:
            Seconds since the epoch, or 0.0 before the first request.

        Example:
            >>> PubChemAPI().last_request_time
            0.0
        """
        return self._http.last_request_time

    def __cache_key__(self) -> tuple:
        """
        Identify this client for cache-key purposes.

        Only the endpoint distinguishes two clients' results; the pause time and
        the ``use_cache`` flag change how a call is made, not what it returns.
        Returning a stable value instead of the default object ``repr`` is what
        lets a cache entry written in one process be found in the next.

        Returns:
            Tuple of the class path and the configured base URL.
        """
        return ("provesid.pubchem.PubChemAPI", self.base_url)

    def clear_cache(self):
        """Clear all cached results for PubChem API"""
        from .cache import clear_pubchem_cache
        clear_pubchem_cache()
            
    def get_cache_info(self):
        """Get cache statistics for PubChem API cached methods"""
        from .cache import get_pubchem_cache_info
        return get_pubchem_cache_info()
        
    def _rate_limit(self):
        """
        Sleep, if needed, so requests to PubChem stay :attr:`pause_time` apart.

        Delegates to the shared transport, which paces every request it makes
        including retries.

        Example:
            >>> PubChemAPI(pause_time=0)._rate_limit()
        """
        self._http.rate_limit()

    def _make_request(self, url: str, method: str = 'GET', data: Optional[Dict] = None,
                     timeout: int = 30, headers: Optional[Dict] = None) -> requests.Response:
        """
        Make one request to PUG-REST and return the raw response.

        The transport handles the pacing, the retries and the back-off; what
        stays here is the shape of a PUG-REST request and the one status code
        that is neither success nor failure. Classification is
        :func:`pugrest_classify`, which reads the fault code in the body
        because PubChem's status codes alone do not distinguish a compound that
        has no such data from a service shedding load.

        Args:
            url: Request URL
            method: HTTP method (GET or POST)
            data: POST data
            timeout: Request timeout in seconds
            headers: HTTP headers

        Returns:
            Response object. Returned unparsed because the caller knows which
            of JSON, text, SDF or PNG it asked for --- see
            :meth:`_parse_response`.

        Raises:
            PubChemNotFoundError: The compound, substance or assay does not
                exist. Never retried; asking again cannot change it.
            PubChemServerError: PubChem stayed busy, or kept throttling, for
                every attempt. Transient fault codes, 429 and 5xx are retried
                first.
            PubChemTimeoutError: Every attempt timed out or could not connect.
            PubChemError: A permanent error in the request itself --- a 400
                naming a property that does not exist, a 403. The message
                carries the start of PubChem's own explanation.
        """
        response = self._http.request(method, url, data=data, timeout=timeout,
                                      headers=headers)
        if response.status_code == 202:
            # A list-key operation that PubChem has accepted but not finished.
            # It is a real answer --- the body holds the key to poll with ---
            # so it is returned rather than retried.
            self.logger.warning(
                "Asynchronous operation pending - may need to poll for results"
            )
        return response

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
            synonyms_error = None
            if include_synonyms:
                try:
                    synonyms_data = self.get_compound_synonyms(cid, output_format)
                except Exception as e:
                    # Don't fail the whole request if synonyms fail, but record
                    # the failure so the partial result is not cached.
                    logging.warning(f"Synonym lookup failed for CID {cid}: {e}")
                    synonyms_data = []
                    synonyms_error = str(e)
            
            # Extract properties from nested structure
            if prop_data and 'PropertyTable' in prop_data and 'Properties' in prop_data['PropertyTable']:
                properties_dict = prop_data['PropertyTable']['Properties'][0]
                
                # Add metadata keys
                properties_dict['success'] = True
                properties_dict['cid'] = cid
                properties_dict['error'] = None
                if include_synonyms:
                    properties_dict['synonyms'] = synonyms_data
                    if synonyms_error is not None:
                        properties_dict['synonyms_error'] = synonyms_error
                
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
                    if synonyms_error is not None:
                        result['synonyms_error'] = synonyms_error
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

    @cached(service='pubchem', skip_if=_synonyms_incomplete)
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

    def _build_post_url(self, domain: str, namespace: str,
                        operation: Optional[str] = None,
                        output_format: str = OutputFormat.JSON) -> str:
        """
        Build a PUG-REST URL whose identifiers travel in the request body.

        A GET request carries its identifiers in the path, which PubChem caps at
        roughly 2000 characters. The POST form of the same call leaves the
        identifier segment out of the path entirely and takes it as a form
        field instead, so the only limit is the body size.

        Args:
            domain: API domain (compound, substance, assay, ...).
            namespace: Namespace within the domain (cid, name, smiles, ...).
            operation: Operation to perform, e.g. ``"property/MolecularWeight"``.
            output_format: Desired output format.

        Returns:
            URL string with no identifier segment, for use with
            ``self._make_request(url, method='POST', data={namespace: ...})``.

        Example:
            >>> api = PubChemAPI()
            >>> api._build_post_url('compound', 'cid', 'property/MolecularWeight')
            'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/MolecularWeight/JSON'
        """
        url_parts = [self.base_url, domain, namespace]
        if operation:
            url_parts.append(operation)
        if output_format:
            url_parts.append(output_format)
        return '/'.join(url_parts)

    @cached(service='pubchem', skip_if=is_empty_result)
    def _get_property_rows(self, cids_tuple: tuple, properties_tuple: tuple) -> List[Dict[str, Any]]:
        """
        Fetch one ``PropertyTable`` covering several CIDs in a single request.

        The identifiers go in the URL path while they fit within
        :data:`URL_IDENTIFIER_LIMIT`, and in a POST body when they do not. Both
        forms hit the same endpoint and return the same payload; the choice is
        purely about the URL length ceiling.

        Args:
            cids_tuple: CIDs to look up, as a tuple so the result is cacheable.
            properties_tuple: Property names, as a tuple for the same reason.

        Returns:
            The rows of ``PropertyTable.Properties``, one per CID that PubChem
            recognised. A CID it does not know still yields a row, but that row
            carries only the ``CID`` key. An empty list means the request
            reported no data at all.

        Raises:
            PubChemError: If the request could not be completed. A failed fetch
                is never reported as an empty list, so a transient error cannot
                be mistaken for — or cached as — "these compounds have no data".
        """
        cids = list(cids_tuple)
        properties = list(properties_tuple)
        operation = f"{Operation.PROPERTY}/{','.join(properties)}"
        identifiers = ','.join(str(cid) for cid in cids)

        try:
            if len(identifiers) > URL_IDENTIFIER_LIMIT:
                url = self._build_post_url(Domain.COMPOUND, CompoundDomainNamespace.CID,
                                           operation, OutputFormat.JSON)
                logging.debug(f"POSTing {len(cids)} CIDs to {url}")
                response = self._make_request(url, method='POST',
                                              data={CompoundDomainNamespace.CID: identifiers})
            else:
                url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.CID, cids,
                                      operation, OutputFormat.JSON)
                response = self._make_request(url)
            data = self._parse_response(response, OutputFormat.JSON)
        except PubChemNotFoundError:
            logging.debug(f"No property record for any of {len(cids)} CIDs")
            return []

        if isinstance(data, dict):
            return data.get('PropertyTable', {}).get('Properties', [])
        return []

    def get_properties_for_cids(self, cids: List[Union[int, str]],
                                properties: List[str],
                                chunk_size: int = PROPERTY_CHUNK_SIZE) -> List[Dict[str, Any]]:
        """
        Get a property table for many CIDs, a few hundred compounds per request.

        This is the bulk counterpart to :meth:`get_compound_properties`, which
        asks about one compound at a time. PubChem's property endpoint accepts a
        comma-separated CID list and answers the whole set in one round trip, so
        a thousand compounds cost five requests here instead of a thousand.

        Args:
            cids: CIDs to look up. Duplicates are collapsed, and the original
                order of first appearance is preserved.
            properties: Property names, e.g.
                ``['MolecularWeight', 'SMILES']``. See
                :class:`CompoundProperties` for the full list.
            chunk_size: How many CIDs to put in one request. The default is
                deliberately below what PubChem tolerates: a smaller chunk
                wastes less work when a request has to be retried.

        Returns:
            One row per distinct CID PubChem answered for, each a dict with a
            ``CID`` key plus the properties it holds. A CID PubChem does not
            know yields a row carrying only ``CID``; a property the compound has
            no value for is absent from the row rather than ``None``, which is
            how PubChem itself reports it.

        Raises:
            ValueError: If ``properties`` is empty or ``chunk_size`` is not positive.
            PubChemError: If a request could not be completed. Partial results
                are not returned: either every chunk succeeded or the failure
                surfaces, so a short table never has to be second-guessed.

        Example:
            >>> api = PubChemAPI()
            >>> rows = api.get_properties_for_cids([2244, 702],
            ...                                    ['MolecularFormula', 'MolecularWeight'])
            >>> {row['CID']: row['MolecularFormula'] for row in rows}
            {2244: 'C9H8O4', 702: 'C2H6O'}
        """
        if not properties:
            raise ValueError("properties must name at least one property")
        if chunk_size <= 0:
            raise ValueError(f"chunk_size must be positive, got {chunk_size}")

        unique_cids = list(dict.fromkeys(cids))
        if not unique_cids:
            return []

        rows: List[Dict[str, Any]] = []
        for start in range(0, len(unique_cids), chunk_size):
            chunk = unique_cids[start:start + chunk_size]
            rows.extend(self._get_property_rows(tuple(chunk), tuple(properties)))
        return rows

    def get_compound_properties_batch(self, cids: List[Union[int, str]],
                                      properties: List[str],
                                      chunk_size: int = PROPERTY_CHUNK_SIZE) -> List[Dict[str, Any]]:
        """
        Get compound properties for multiple CIDs, one row per CID.

        A convenience layer over :meth:`get_properties_for_cids` that reshapes
        the raw property table into the flat, self-describing dicts
        :meth:`get_compound_properties` returns, so a caller can iterate over
        the result without checking which CIDs came back.

        Args:
            cids: List of CIDs. Duplicates are answered once each, in the order
                they first appear.
            properties: List of property names.
            chunk_size: How many CIDs to put in one request.

        Returns:
            One dict per distinct CID, in the order requested, each carrying the
            retrieved properties plus ``success``, ``cid`` and ``error`` keys. A
            CID PubChem does not know is reported with ``success=False`` and an
            explanatory ``error`` rather than being dropped, so the result can
            be zipped against the input.

        Raises:
            ValueError: If ``properties`` is empty or ``chunk_size`` is not positive.
            PubChemError: If a request could not be completed. A transport
                failure is not reported as a row of missing properties.

        Note:
            This asks PubChem about ``chunk_size`` compounds per request rather
            than one compound per request, so it does not share cache entries
            with :meth:`get_compound_properties`. Synonyms are not included;
            they need one request per compound, which defeats the point of
            batching. Use :meth:`get_compound_synonyms` where they are needed.

        Example:
            >>> api = PubChemAPI()
            >>> rows = api.get_compound_properties_batch([2244, 702], ['MolecularFormula'])
            >>> for row in rows:
            ...     print(row['cid'], row.get('MolecularFormula'))
            2244 C9H8O4
            702 C2H6O
        """
        unique_cids = list(dict.fromkeys(cids))
        rows = self.get_properties_for_cids(unique_cids, properties, chunk_size=chunk_size)

        # PubChem answers with whatever CIDs it recognised, in its own order, so
        # index the table before walking the caller's list.
        by_cid = {str(row.get('CID')): row for row in rows}

        results: List[Dict[str, Any]] = []
        for cid in unique_cids:
            row = by_cid.get(str(cid))
            if row is None:
                results.append({
                    "success": False,
                    "cid": cid,
                    "error": "CID not present in the property table returned by PubChem",
                })
                continue

            record = {key: value for key, value in row.items() if key != 'CID'}
            if not record:
                # A bare CID row is how PubChem reports a compound it has no
                # record for; it is an answer, not a failed request.
                results.append({
                    "success": False,
                    "cid": cid,
                    "error": "No such compound",
                })
                continue

            record['success'] = True
            record['cid'] = cid
            record['error'] = None
            results.append(record)

        return results

    @cached(service='pubchem', skip_if=is_empty_result)
    def get_compound_synonyms(self, cid: Union[int, str], output_format: str = OutputFormat.JSON) -> List[str]:
        """
        Get compound synonyms by CID
        
        Args:
            cid: Compound ID
            output_format: Desired output format
            
        Returns:
            List of synonyms (flattened from nested structure). An empty list
            means PubChem lists no synonyms for this CID.

        Raises:
            PubChemError: If the request could not be completed. A failed fetch
                is never reported as an empty list, so that a transient error
                does not get cached as "this compound has no synonyms".
        """
        try:
            url = self._build_url(Domain.COMPOUND, CompoundDomainNamespace.CID, cid,
                                 Operation.SYNONYMS, output_format)
            response = self._make_request(url)
            raw_data = self._parse_response(response, output_format)
        except PubChemNotFoundError:
            logging.debug(f"No synonym record for CID {cid}")
            return []
        
        # Extract synonyms from nested structure
        if raw_data and 'InformationList' in raw_data:
            info_list = raw_data['InformationList'].get('Information', [])
            if info_list and len(info_list) > 0:
                return info_list[0].get('Synonym', [])
        
        # Return empty list if no synonyms found
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
    
    DEFAULT_DB_NAME = "pubchem_id.db"
    DEFAULT_DB_URL = "https://zenodo.org/records/18173204/files/pubchem_id.db"

    #: PubChem property names the local database can answer, mapped to their
    #: column in the ``compounds`` table. The database is built from PubChem's
    #: CAS export, which carries the identifiers and the cheap computed
    #: descriptors but not the full property set, so anything outside this
    #: mapping — ``MonoisotopicMass``, ``ConnectivitySMILES``, the 3D
    #: descriptors, the patent and literature counts — needs the online API.
    #: Note that ``smiles`` holds the isomeric SMILES, which is what PubChem now
    #: calls ``SMILES``; the stereochemistry-free ``ConnectivitySMILES`` is not
    #: stored locally.
    OFFLINE_PROPERTIES = {
        'MolecularFormula': 'mf',
        'MolecularWeight': 'mw',
        'SMILES': 'smiles',
        'InChI': 'inchi',
        'InChIKey': 'inchikey',
        'IUPACName': 'iupacname',
        'Title': 'cmpdname',
        'XLogP': 'xlogp',
        'TPSA': 'polararea',
        'Complexity': 'complexity',
        'Charge': 'charge',
        'HBondDonorCount': 'hbonddonor',
        'HBondAcceptorCount': 'hbondacc',
        'RotatableBondCount': 'rotbonds',
        'HeavyAtomCount': 'heavycnt',
        'ExactMass': 'exactmass',
    }

    #: What :meth:`properties` retrieves when the caller names no properties:
    #: everything available without touching the network.
    DEFAULT_PROPERTIES = tuple(OFFLINE_PROPERTIES)

    #: Type each property is normalised to, so that a table assembled from both
    #: sources is usable as one table. PUG-REST reports ``MolecularWeight`` and
    #: ``ExactMass`` as strings while the local database holds floats.
    _PROPERTY_CASTS = {
        'MolecularWeight': float,
        'ExactMass': float,
        'MonoisotopicMass': float,
        'XLogP': float,
        'TPSA': float,
        'Complexity': float,
        'Charge': int,
        'HBondDonorCount': int,
        'HBondAcceptorCount': int,
        'RotatableBondCount': int,
        'HeavyAtomCount': int,
    }

    #: Bound parameters per ``IN`` clause. SQLite's own default ceiling is 999.
    _SQL_PARAMETER_LIMIT = 500


    def __init__(
        self,
        db_path: Optional[str] = None,
        auto_download: bool = True,
        data_dir: Optional[str] = None,
        db_url: Optional[str] = None,
        redownload: bool = False,
        api: Optional['PubChemAPI'] = None,
    ):
        """
        Initialize PubChemID database connection.
        
        Args:
            db_path (str, optional): Path to SQLite database. If None, uses default
                                    location in the persistent user dataset directory.
            auto_download (bool): If True, automatically download database if not found.
                                 Default is True.
            data_dir (str, optional): Directory to store the database when
                ``db_path`` is not provided.
            db_url (str, optional): Download URL for the database. If None,
                uses the package default URL.
            redownload (bool): If True, force re-download when
                ``auto_download`` is enabled.
            api (PubChemAPI, optional): Online client used by
                :meth:`properties` when the local database cannot answer a
                request. One is created on first use if none is given, so
                passing this is only needed to share a client or to configure
                its pause time.
        
        Raises:
            FileNotFoundError: If database file doesn't exist and auto_download is False
        """
        import sqlite3

        self.logger = logging.getLogger(__name__)
        self.db_url = db_url or self.DEFAULT_DB_URL
        self._api = api
        
        if db_path is None:
            base_dir = data_dir or user_dataset_path()
            db_path = os.path.join(base_dir, self.DEFAULT_DB_NAME)
        
        self.db_path = os.path.abspath(os.path.expanduser(db_path))

        needs_download = redownload or not os.path.exists(self.db_path)

        if needs_download:
            if auto_download:
                if redownload and os.path.exists(self.db_path):
                    self.logger.info(
                        "Forced PubChemID redownload requested for: %s", self.db_path
                    )
                else:
                    self.logger.info("Database not found at %s", self.db_path)
                self.logger.info("Attempting to download from configured source...")
                self.download_database(
                    db_path=self.db_path,
                    zenodo_url=self.db_url,
                    force=redownload,
                )
            else:
                raise FileNotFoundError(
                    f"PubChem ID database not found at {self.db_path}. "
                    "Set auto_download=True or run PubChemID.download_database()."
                )
        
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row  # Access columns by name
    
    def __del__(self):
        """Close database connection on deletion."""
        if hasattr(self, 'conn'):
            self.conn.close()
    
    @staticmethod
    def download_database(
        db_path: Optional[str] = None,
        zenodo_url: Optional[str] = None,
        force: bool = False,
    ) -> str:
        """
        Download PubChem ID database from Zenodo.
        
        Args:
            db_path (str, optional): Path where to save the database. If None, uses default
                                    location in the persistent user dataset directory.
            zenodo_url (str, optional): URL to download from. If None, uses default Zenodo URL.
                                       Format: https://zenodo.org/record/XXXXXX/files/pubchem_id.db
            force (bool): If True, overwrite an existing local database file.
        
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
        from tqdm import tqdm

        logger = logging.getLogger(__name__)
        
        if db_path is None:
            db_path = os.path.join(
                user_dataset_path(),
                PubChemID.DEFAULT_DB_NAME,
            )
        else:
            db_path = os.path.abspath(os.path.expanduser(db_path))

        os.makedirs(os.path.dirname(db_path), exist_ok=True)

        if os.path.exists(db_path) and not force:
            raise FileExistsError(
                f"Database already exists at: {db_path}. Use force=True to overwrite."
            )
        
        if zenodo_url is None:
            zenodo_url = PubChemID.DEFAULT_DB_URL
        
        logger.info("Downloading PubChem ID database from: %s", zenodo_url)
        logger.info("Destination: %s", db_path)
        logger.info("This is a large file (~2.2 GB), please be patient.")
        
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
            
            logger.info("Download complete. Verifying...")
            
            # Verify it's a valid SQLite database
            import sqlite3
            try:
                conn = sqlite3.connect(temp_path)
                cursor = conn.cursor()
                cursor.execute("SELECT COUNT(*) FROM compounds")
                count = cursor.fetchone()[0]
                conn.close()
                logger.info("Database verified: %s compounds", f"{count:,}")
            except Exception as e:
                raise RuntimeError(f"Downloaded file is not a valid database: {e}")
            
            # Move to final location
            if os.path.exists(db_path):
                os.remove(db_path)
            os.rename(temp_path, db_path)
            
            logger.info("Database ready at %s", db_path)
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

    def get_by_smiles(self, smiles: str) -> Optional[Dict[str, Any]]:
        """
        Get compound information by SMILES string.
        
        Args:
            smiles (str): SMILES string
        
        Returns:
            dict: Compound information, or None if not found
        
        Example:
            >>> db = PubChemID()
            >>> result = db.get_by_smiles("CC(=O)OC1=CC=CC=C1C(=O)O")  # Aspirin
            >>> print(result['cmpdname'])
        """
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT * FROM compounds WHERE smiles = ?
        """, (smiles,))
        
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
    
    def cid_to_smiles(self, cid: int) -> Optional[str]:
        """Convert PubChem CID to SMILES."""
        result = self.get_by_cid(cid)
        return result['smiles'] if result else None
    
    def smiles_to_cid(self, smiles: str) -> Optional[int]:
        """Convert SMILES string to PubChem CID."""
        result = self.get_by_smiles(smiles)
        return result['cid'] if result else None
    
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
    
    def batch_smiles_to_cid(self, smiles_list: List[str]) -> Dict[str, Optional[int]]:
        """
        Convert multiple SMILES strings to CIDs.
        
        Args:
            smiles_list (list): List of SMILES strings
        
        Returns:
            dict: Mapping of SMILES -> CID (None if not found)
        
        Example:
            >>> db = PubChemID()
            >>> results = db.batch_smiles_to_cid(["CC(=O)OC1=CC=CC=C1C(=O)O", "C"])
            >>> print(results)
        """
        results = {}
        for smiles in smiles_list:
            results[smiles] = self.smiles_to_cid(smiles)
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
    
    def get_by_smiles_batch(self, smiles_list: List[str]) -> 'pd.DataFrame':
        """
        Get complete compound information for multiple SMILES strings as a DataFrame.
        
        This method returns all available data including identifiers, chemical properties,
        and physical properties for each SMILES string.
        
        Args:
            smiles_list (list): List of SMILES strings
        
        Returns:
            pandas.DataFrame: DataFrame with columns for all compound properties including:
                             cid, cas, inchi, inchikey, smiles, cmpdname, iupacname, mf, mw,
                             polararea, complexity, xlogp, heavycnt, hbonddonor, hbondacc,
                             rotbonds, exactmass, charge, cidcdate
        
        Example:
            >>> db = PubChemID()
            >>> smiles_list = ["CC(=O)OC1=CC=CC=C1C(=O)O", "C", "CCO"]
            >>> df = db.get_by_smiles_batch(smiles_list)
            >>> print(df[['smiles', 'cmpdname', 'mf', 'mw']])
        """
        import pandas as pd
        
        rows = []
        for smiles in smiles_list:
            result = self.get_by_smiles(smiles)
            if result:
                # Create row with all properties
                row = {
                    'cid': result.get('cid'),
                    'cas': result.get('cas_numbers', [])[0] if result.get('cas_numbers') else None,
                    'inchi': result.get('inchi', ''),
                    'inchikey': result.get('inchikey', ''),
                    'smiles': smiles,
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
    
    # Additional CAS conversion methods
    
    def smiles_to_cas(self, smiles: str) -> Optional[List[str]]:
        """
        Convert SMILES string to CAS number(s).
        
        Args:
            smiles (str): SMILES string
        
        Returns:
            list: List of CAS numbers, or None if not found
        
        Example:
            >>> db = PubChemID()
            >>> cas_list = db.smiles_to_cas("CC(=O)OC1=CC=CC=C1C(=O)O")  # Aspirin
            >>> print(cas_list)
        """
        # First convert SMILES to InChI using RDKit
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            inchi = Chem.MolToInchi(mol)
        except Exception:
            return None
        
        # Then look up by InChI
        return self.inchi_to_cas(inchi)
    
    def name_to_cas(self, name: str, exact: bool = True) -> Optional[List[str]]:
        """
        Convert chemical name to CAS number(s).
        
        Args:
            name (str): Chemical name or synonym
            exact (bool): If True, exact match only. If False, returns first match from search.
        
        Returns:
            list: List of CAS numbers, or None if not found
        
        Example:
            >>> db = PubChemID()
            >>> cas_list = db.name_to_cas("aspirin")
            >>> print(cas_list)
        
        Note:
            For exact=False, only the first match from the search is returned.
            Use search_by_name() for more control over multiple matches.
        """
        results = self.search_by_name(name, exact=exact, limit=1)
        if not results:
            return None
        return results[0].get('cas_numbers')
    
    def formula_to_cas(self, formula: str, limit: int = 100) -> Optional[List[str]]:
        """
        Convert molecular formula to CAS numbers.
        
        Note: Molecular formulas are not unique - many isomers can share the same formula.
        This method returns CAS numbers for all compounds matching the formula.
        
        Args:
            formula (str): Molecular formula (e.g., "C9H8O4", "CH2O")
            limit (int): Maximum number of compounds to retrieve
        
        Returns:
            list: List of CAS numbers for all compounds with this formula, or None if not found
        
        Example:
            >>> db = PubChemID()
            >>> cas_list = db.formula_to_cas("C9H8O4")
            >>> print(f"Found {len(cas_list)} CAS numbers for C9H8O4")
        
        Warning:
            Can return many results for common formulas. Use limit parameter to control.
        """
        results = self.search_by_formula(formula, limit=limit)
        if not results:
            return None
        
        # Collect all unique CAS numbers from all matching compounds
        all_cas = []
        for compound in results:
            cas_numbers = compound.get('cas_numbers', [])
            if cas_numbers:
                all_cas.extend(cas_numbers)
        
        # Remove duplicates and sort
        unique_cas = sorted(set(all_cas))
        return unique_cas if unique_cas else None
    
    def batch_smiles_to_cas(self, smiles_list: List[str]) -> Dict[str, Optional[List[str]]]:
        """
        Convert multiple SMILES strings to CAS numbers.
        
        Args:
            smiles_list (list): List of SMILES strings
        
        Returns:
            dict: Mapping of SMILES -> list of CAS numbers (None if not found)
        
        Example:
            >>> db = PubChemID()
            >>> smiles = ["C", "CO", "CCO"]
            >>> results = db.batch_smiles_to_cas(smiles)
            >>> for smi, cas in results.items():
            ...     print(f"{smi}: {cas}")
        """
        return {smiles: self.smiles_to_cas(smiles) for smiles in smiles_list}
    
    def batch_name_to_cas(self, name_list: List[str], exact: bool = True) -> Dict[str, Optional[List[str]]]:
        """
        Convert multiple chemical names to CAS numbers.
        
        Args:
            name_list (list): List of chemical names
            exact (bool): If True, exact match only
        
        Returns:
            dict: Mapping of name -> list of CAS numbers (None if not found)
        
        Example:
            >>> db = PubChemID()
            >>> names = ["aspirin", "caffeine", "glucose"]
            >>> results = db.batch_name_to_cas(names)
            >>> for name, cas in results.items():
            ...     print(f"{name}: {cas}")
        """
        return {name: self.name_to_cas(name, exact=exact) for name in name_list}
    
    def batch_formula_to_cas(self, formula_list: List[str], limit: int = 100) -> Dict[str, Optional[List[str]]]:
        """
        Convert multiple molecular formulas to CAS numbers.
        
        Args:
            formula_list (list): List of molecular formulas
            limit (int): Maximum number of compounds per formula
        
        Returns:
            dict: Mapping of formula -> list of CAS numbers (None if not found)
        
        Example:
            >>> db = PubChemID()
            >>> formulas = ["H2O", "CH4", "C9H8O4"]
            >>> results = db.batch_formula_to_cas(formulas)
            >>> for formula, cas_list in results.items():
            ...     if cas_list:
            ...         print(f"{formula}: {len(cas_list)} compounds")
        """
        return {formula: self.formula_to_cas(formula, limit=limit) for formula in formula_list}
    
    @property
    def api(self) -> 'PubChemAPI':
        """
        The online client used to answer what the local database cannot.

        Created on first use rather than in ``__init__``, so a strictly offline
        session never builds one.

        Returns:
            The :class:`PubChemAPI` instance passed to ``__init__``, or one
            created with default settings.
        """
        if self._api is None:
            self._api = PubChemAPI()
        return self._api

    def properties(self, cid: Union[int, str],
                   properties: Optional[List[str]] = None,
                   use_online_fallback: bool = True) -> Optional[Dict[str, Any]]:
        """
        Look up computed properties for one compound, offline first.

        The local database answers from disk in microseconds; the online API is
        consulted only when the local database cannot serve the request, either
        because it holds no row for this CID or because a requested property is
        not one of the columns it carries (see :attr:`OFFLINE_PROPERTIES`).

        Args:
            cid: PubChem Compound ID.
            properties: Property names to retrieve, e.g.
                ``['MolecularWeight', 'XLogP']``. Defaults to every property the
                local database can answer, :attr:`DEFAULT_PROPERTIES`.
            use_online_fallback: When True (default), fall back to PUG-REST for
                anything the local database cannot answer. When False, the
                lookup is strictly offline and an unavailable property is simply
                absent from the result.

        Returns:
            A dict carrying ``CID``, a ``Source`` of ``'offline'`` or
            ``'online'``, and one key per property that has a value. A property
            the compound has no value for is omitted rather than set to None,
            which is how PubChem itself reports it — so ``'XLogP' not in
            result`` means PubChem computes no logP for this compound, not that
            the lookup fell short. Returns None when neither source knows the
            CID.

        Raises:
            ValueError: If ``cid`` is not an integer, or ``properties`` is an
                empty list.
            PubChemError: If the online fallback was needed and its request
                could not be completed. An incomplete answer is never passed off
                as a complete one.

        Example:
            >>> db = PubChemID()
            >>> db.properties(2244, ['MolecularFormula', 'MolecularWeight'])
            {'CID': 2244, 'Source': 'offline', 'MolecularFormula': 'C9H8O4', 'MolecularWeight': 180.16}
            >>> # MonoisotopicMass is not in the local database, so this one goes online
            >>> db.properties(2244, ['MonoisotopicMass'])['Source']
            'online'
        """
        rows = self.properties_for_cids([cid], properties,
                                        use_online_fallback=use_online_fallback)
        return rows[0] if rows else None

    def properties_for_cids(self, cids: List[Union[int, str]],
                            properties: Optional[List[str]] = None,
                            use_online_fallback: bool = True,
                            chunk_size: int = PROPERTY_CHUNK_SIZE) -> List[Dict[str, Any]]:
        """
        Look up computed properties for many compounds, offline first.

        Everything the local database can answer is read in a handful of SQL
        statements; only the remainder is requested from PubChem, in bulk, a few
        hundred compounds per request. A list of ten thousand CIDs that the
        local database covers therefore costs no network traffic at all.

        Args:
            cids: PubChem Compound IDs. Duplicates are collapsed and the order
                of first appearance is preserved.
            properties: Property names to retrieve. Defaults to
                :attr:`DEFAULT_PROPERTIES`.
            use_online_fallback: When True (default), CIDs the local database
                does not cover are requested from PUG-REST.
            chunk_size: How many CIDs to put in one online request.

        Returns:
            One dict per CID that could be answered, in the order requested,
            each carrying ``CID``, a ``Source`` of ``'offline'`` or
            ``'online'``, and one key per property that has a value. CIDs
            neither source knows are omitted; use :meth:`properties_table` to
            get a row for every CID asked about.

        Raises:
            ValueError: If a CID is not an integer, ``properties`` is an empty
                list, or ``chunk_size`` is not positive.
            PubChemError: If an online request could not be completed.

        Note:
            If *any* requested property lies outside
            :attr:`OFFLINE_PROPERTIES`, the whole request goes online: the
            missing property would need a request per compound anyway, so
            splitting the property list between the two sources would cost the
            same traffic and return rows assembled from two different PubChem
            snapshots.

        Example:
            >>> db = PubChemID()
            >>> rows = db.properties_for_cids([2244, 702], ['MolecularFormula'])
            >>> for row in rows:
            ...     print(row['CID'], row['MolecularFormula'], row['Source'])
            2244 C9H8O4 offline
            702 C2H6O offline
        """
        if properties is not None and not properties:
            raise ValueError("properties must name at least one property, or be None")
        if chunk_size <= 0:
            raise ValueError(f"chunk_size must be positive, got {chunk_size}")

        requested_properties = list(properties) if properties else list(self.DEFAULT_PROPERTIES)
        wanted_cids = [self._coerce_cid(cid) for cid in cids]
        wanted_cids = list(dict.fromkeys(wanted_cids))
        if not wanted_cids:
            return []

        online_only = [name for name in requested_properties
                       if name not in self.OFFLINE_PROPERTIES]

        found: Dict[int, Dict[str, Any]] = {}
        if online_only:
            self.logger.debug(
                "Going straight online for %d CIDs: %s not in the local database",
                len(wanted_cids), ', '.join(online_only))
            missing = wanted_cids
        else:
            found = self._offline_properties(wanted_cids, requested_properties)
            missing = [cid for cid in wanted_cids if cid not in found]
            self.logger.debug("Served %d/%d CIDs offline", len(found), len(wanted_cids))

        if missing and use_online_fallback:
            self.logger.debug("Falling back online for %d CIDs", len(missing))
            found.update(self._online_properties(missing, requested_properties, chunk_size))

        return [found[cid] for cid in wanted_cids if cid in found]

    def properties_table(self, cids: List[Union[int, str]],
                         properties: Optional[List[str]] = None,
                         use_online_fallback: bool = True,
                         chunk_size: int = PROPERTY_CHUNK_SIZE) -> 'pd.DataFrame':
        """
        Offline-first property lookup for many compounds, as a DataFrame.

        Same lookup as :meth:`properties_for_cids`, reshaped so that every CID
        asked about has a row whether or not it could be answered. That makes
        the frame safe to concatenate or join against the caller's own table.

        Args:
            cids: PubChem Compound IDs. Duplicates are collapsed.
            properties: Property names to retrieve. Defaults to
                :attr:`DEFAULT_PROPERTIES`.
            use_online_fallback: When True (default), consult PUG-REST for CIDs
                the local database does not cover.
            chunk_size: How many CIDs to put in one online request.

        Returns:
            A DataFrame with one row per distinct CID in the order requested.
            Columns are ``CID``, ``Source`` and the requested properties.
            ``Source`` reads ``'offline'``, ``'online'``, or ``'missing'`` for a
            CID neither source knows; a property with no value is NaN/None.

        Raises:
            ValueError: If a CID is not an integer, ``properties`` is an empty
                list, or ``chunk_size`` is not positive.
            PubChemError: If an online request could not be completed.

        Example:
            >>> db = PubChemID()
            >>> table = db.properties_table([2244, 702], ['MolecularWeight'])
            >>> table[['CID', 'MolecularWeight', 'Source']].to_dict('records')
            [{'CID': 2244, 'MolecularWeight': 180.16, 'Source': 'offline'},
             {'CID': 702, 'MolecularWeight': 46.07, 'Source': 'offline'}]
        """
        requested_properties = list(properties) if properties else list(self.DEFAULT_PROPERTIES)
        rows = self.properties_for_cids(cids, requested_properties,
                                        use_online_fallback=use_online_fallback,
                                        chunk_size=chunk_size)
        by_cid = {row['CID']: row for row in rows}

        records = []
        for cid in dict.fromkeys(self._coerce_cid(cid) for cid in cids):
            row = by_cid.get(cid, {'CID': cid, 'Source': 'missing'})
            records.append({'CID': cid, 'Source': row['Source'],
                            **{name: row.get(name) for name in requested_properties}})

        return pd.DataFrame(records, columns=['CID', 'Source'] + requested_properties)

    @staticmethod
    def _coerce_cid(cid: Union[int, str]) -> int:
        """
        Normalise a CID to an int so that ``2244`` and ``"2244"`` share a row.

        Args:
            cid: CID as an int or a string of digits.

        Returns:
            The CID as an int.

        Raises:
            ValueError: If ``cid`` is not an integer. PubChem answers a
                malformed CID with a blanket ``PUGREST.BadRequest`` that fails
                the whole batch, so it is worth catching here, where the
                offending value can be named.
        """
        try:
            return int(cid)
        except (TypeError, ValueError):
            raise ValueError(f"CID must be an integer, got {cid!r}")

    def _offline_properties(self, cids: List[int],
                            properties: List[str]) -> Dict[int, Dict[str, Any]]:
        """
        Read properties for the given CIDs from the local database.

        Args:
            cids: CIDs to look up, already coerced to int.
            properties: Property names, all of which must be keys of
                :attr:`OFFLINE_PROPERTIES`.

        Returns:
            A dict keyed by CID, holding one record per CID present in the
            database. A record carries ``CID``, ``Source='offline'`` and the
            properties that have a value; a NULL column is left out, matching
            PubChem, which omits a property rather than reporting it as null.
        """
        columns = [self.OFFLINE_PROPERTIES[name] for name in properties]
        cursor = self.conn.cursor()
        found: Dict[int, Dict[str, Any]] = {}

        # SQLite allows a limited number of bound parameters per statement
        # (999 by default), so the IN list is filled in batches.
        for start in range(0, len(cids), self._SQL_PARAMETER_LIMIT):
            batch = cids[start:start + self._SQL_PARAMETER_LIMIT]
            placeholders = ','.join('?' * len(batch))
            cursor.execute(
                f"SELECT cid, {', '.join(columns)} FROM compounds "
                f"WHERE cid IN ({placeholders})", batch)
            for row in cursor.fetchall():
                record = {'CID': row['cid'], 'Source': 'offline'}
                for name, column in zip(properties, columns):
                    value = row[column]
                    if value is not None and value != '':
                        record[name] = self._cast_property(name, value)
                found[row['cid']] = record

        return found

    def _online_properties(self, cids: List[int], properties: List[str],
                           chunk_size: int) -> Dict[int, Dict[str, Any]]:
        """
        Fetch properties for the given CIDs from PUG-REST.

        Args:
            cids: CIDs the local database could not answer.
            properties: Property names to request.
            chunk_size: How many CIDs to put in one request.

        Returns:
            A dict keyed by CID, holding one record per CID PubChem answered
            for, shaped like the offline records but with
            ``Source='online'``. PubChem returns a bare CID for a compound it
            has no record of; such a row is dropped, so the CID is reported as
            unknown rather than as a compound with no properties.

        Raises:
            PubChemError: If a request could not be completed.
        """
        rows = self.api.get_properties_for_cids(cids, properties, chunk_size=chunk_size)

        found: Dict[int, Dict[str, Any]] = {}
        for row in rows:
            values = {name: self._cast_property(name, row[name])
                      for name in properties
                      if row.get(name) is not None and row.get(name) != ''}
            if not values:
                continue
            cid = self._coerce_cid(row['CID'])
            found[cid] = {'CID': cid, 'Source': 'online', **values}

        return found

    @classmethod
    def _cast_property(cls, name: str, value: Any) -> Any:
        """
        Coerce a property value to one consistent type across both sources.

        The two sources disagree on types for the same property: PUG-REST
        returns ``MolecularWeight`` as the string ``"180.16"`` while the local
        database holds it as a float, and a table assembled from both sources
        has to be usable as one table.

        Args:
            name: Property name.
            value: Raw value from either source.

        Returns:
            The value cast to the type recorded in ``_PROPERTY_CASTS``, or
            unchanged when no cast is recorded or the cast does not apply. An
            uncastable value is returned as-is rather than discarded: a
            surprising value is more useful to the caller than a silent hole.
        """
        cast = cls._PROPERTY_CASTS.get(name)
        if cast is None:
            return value
        try:
            return cast(value)
        except (TypeError, ValueError):
            logging.debug("Could not cast %s=%r with %s", name, value, cast.__name__)
            return value

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
