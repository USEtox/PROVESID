"""
The NCI/CADD Chemical Identifier Resolver (CACTUS), online.

CACTUS (https://cactus.nci.nih.gov/chemical/structure) turns any identifier it
recognises --- a name, a CAS number, SMILES, InChI or InChIKey --- into
another representation, through one URL shape:
``/chemical/structure/{identifier}/{representation}``.
[`NCIChemicalIdentifierResolver`][provesid.resolver.NCIChemicalIdentifierResolver]
wraps it with pacing, retries and caching, and the ``nci_*`` functions are
one-line shortcuts that return None instead of raising.

Every call needs the network. [`Search`][provesid.search.Search] asks CACTUS
only when its offline sources cannot answer and ``online_fallback`` allows.

Examples:
    >>> from provesid.resolver import NCIChemicalIdentifierResolver, nci_get_formula
    >>> resolver = NCIChemicalIdentifierResolver()
    >>> resolver.resolve("CCO", "stdinchikey")                 # doctest: +SKIP
    'InChIKey=LFQSCWFLJHTTHZ-UHFFFAOYSA-N'
    >>> nci_get_formula("caffeine")                            # doctest: +SKIP
    'C8H10N4O2'
"""

import logging
import re
from urllib.parse import quote
from typing import Dict, List, Optional, Any
import requests
from .cache import cached
from .http import (
    HTTPClient,
    Outcome,
    ServiceError,
    NotFoundError,
    ServiceTimeoutError,
    default_classify,
)

logger = logging.getLogger(__name__)

class NCIResolverError(ServiceError):
    """
    Custom exception for NCI Chemical Identifier Resolver errors.

    The base of the two below, and raised on its own for an empty identifier
    or a failure that outlived the retries.

    Examples:
        >>> NCIChemicalIdentifierResolver().resolve("  ", "smiles")
        Traceback (most recent call last):
        ...
        provesid.resolver.NCIResolverError: Empty or invalid identifier provided
    """
    pass

class NCIResolverNotFoundError(NCIResolverError, NotFoundError):
    """
    Exception raised when chemical identifier is not found.

    CACTUS reports this as HTTP 500 with a "Page not found" body;
    [`nci_classify`][provesid.resolver.nci_classify] reads the body so that it
    is not retried as a server fault.

    Examples:
        >>> NCIChemicalIdentifierResolver().resolve(
        ...     "this_is_definitely_not_a_chemical_12345", "smiles")  # doctest: +SKIP
        Traceback (most recent call last):
        ...
        provesid.resolver.NCIResolverNotFoundError: No data for https://cactus.nci.nih.gov/chemical/structure/this_is_definitely_not_a_chemical_12345/smiles (HTTP 500)
    """
    pass

class NCIResolverTimeoutError(NCIResolverError, ServiceTimeoutError):
    """
    Exception raised when request times out.

    Raised after every attempt timed out or failed to connect. Also a
    [`ServiceTimeoutError`][provesid.http.ServiceTimeoutError].

    Examples:
        >>> issubclass(NCIResolverTimeoutError, ServiceTimeoutError)
        True
    """
    pass

_NOT_FOUND_BODY = re.compile(r"Page not found", re.IGNORECASE)
r"""CACTUS reports an identifier it cannot resolve with HTTP 500 and a body of
``<h1>Page not found (404)</h1>``. Its status code is not a reliable guide,
so the body decides --- verified live on 2026-09-19 against
``this_is_definitely_not_a_chemical_12345``, which answers 500/404-body while
``\u03b1-glucose`` answers 200.
"""


def nci_classify(response: requests.Response) -> Outcome:
    """
    Decide what a CACTUS response means, reading the body behind a 500.

    The resolver signals "I cannot resolve this" with a 500 whose body is a
    Django 404 page. Treating that as a server fault costs four requests and
    several seconds of back-off to learn a permanent answer, so the body is
    checked before the status is believed. Everything else is classified by
    status in the usual way.

    Args:
        response: The response to classify.

    Returns:
        [`ABSENT`][provesid.http.Outcome] for a 404, and for a 5xx whose
        body is CACTUS's not-found page; otherwise whatever
        [`default_classify`][provesid.http.default_classify] says.

    Examples:
        >>> class R:
        ...     status_code = 500
        ...     text = '<h1>Page not found (404)</h1>'
        >>> nci_classify(R()).name
        'ABSENT'
    """
    if response.status_code >= 500:
        try:
            body = response.text
        except Exception:
            body = ""
        if _NOT_FOUND_BODY.search(body or ""):
            return Outcome.ABSENT
    return default_classify(response)


def _lines(text: str) -> List[str]:
    """
    Split a CACTUS answer that holds one value per line.

    Args:
        text: The answer, as ``resolve`` returns it.

    Returns:
        The non-empty lines, stripped, in CACTUS's order.

    Examples:
        >>> _lines("64-17-5\\n8024-45-1\\n")
        ['64-17-5', '8024-45-1']
    """
    return [line.strip() for line in text.splitlines() if line.strip()]


def _strip_inchikey_prefix(text: str) -> str:
    """
    Remove the ``InChIKey=`` CACTUS writes before a key.

    Args:
        text: The ``stdinchikey`` answer.

    Returns:
        The bare 27-character key.

    Examples:
        >>> _strip_inchikey_prefix("InChIKey=LFQSCWFLJHTTHZ-UHFFFAOYSA-N")
        'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'
    """
    text = text.strip()
    return text[len("InChIKey="):] if text.startswith("InChIKey=") else text


class NCIChemicalIdentifierResolver:
    """
    A Python interface to the NCI Chemical Identifier Resolver web service

    This class provides methods to interact with the NCI CADD Group's Chemical Identifier
    Resolver service for converting between different chemical structure identifiers.

    The service can resolve various types of chemical identifiers and convert them
    into different representations.

    URL API scheme: https://cactus.nci.nih.gov/chemical/structure/{identifier}/{representation}

    [`resolve`][provesid.resolver.NCIChemicalIdentifierResolver.resolve]
    answers as CACTUS writes: text, several values one per line (``names``,
    ``cas``), InChIKeys with their ``InChIKey=`` prefix.
    [`get_molecular_data`][provesid.resolver.NCIChemicalIdentifierResolver.get_molecular_data]
    and the ``nci_*_to_mol`` functions parse it into lists, a bare key and a
    float.

    Examples:
        >>> resolver = NCIChemicalIdentifierResolver()
        >>> resolver.resolve('CCO', 'stdinchi')                  # doctest: +SKIP
        'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'
        >>> resolver.resolve('aspirin', 'names').splitlines()[:2]  # doctest: +SKIP
        ['2-acetyloxybenzoic acid', '2-Acetoxybenzoic acid']
        >>> resolver.resolve('caffeine', 'mw')                   # doctest: +SKIP
        '194.1926'
    """

    def __init__(self, base_url: str = "https://cactus.nci.nih.gov/chemical/structure",
                 timeout: int = 30, pause_time: float = 0.1, use_cache: bool = True):
        """
        Initialize NCI Chemical Identifier Resolver client

        Args:
            base_url: Base URL for the NCI resolver service
            timeout: Request timeout in seconds
            pause_time: Minimum time between API calls in seconds
            use_cache: Whether to use cache for lookups (default: True).
                      When False, skips cache lookup but still stores results.

        Examples:
            >>> resolver = NCIChemicalIdentifierResolver(timeout=60, pause_time=0.5)
            >>> resolver.timeout, resolver.pause_time
            (60, 0.5)
        """
        self.base_url = base_url.rstrip('/')
        self.timeout = timeout
        self.use_cache = use_cache
        self.logger = logger

        # One shared transport, configured with this service's exceptions so
        # callers keep catching NCIResolverError and friends. CACTUS hides a
        # not-found behind a 500, hence its own classifier. ``pace_host`` gives
        # every resolver aimed at CACTUS one clock, so two instances in a
        # process cannot together ask twice as fast as either promises.
        self._http = HTTPClient(
            min_interval=pause_time,
            timeout=timeout,
            pace_host=self.base_url,
            classify=nci_classify,
            error_cls=NCIResolverError,
            not_found_cls=NCIResolverNotFoundError,
            timeout_cls=NCIResolverTimeoutError,
            logger=self.logger,
        )

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

    def __cache_key__(self) -> tuple:
        """
        Identify this client for cache-key purposes.

        Only the endpoint distinguishes two clients' results; timeout, pause
        time and ``use_cache`` change how a call is made, not what it returns.

        Returns:
            Tuple of the class path and the configured base URL.
        """
        return ("provesid.resolver.NCIChemicalIdentifierResolver", self.base_url)

    def clear_cache(self):
        """
        Delete every cached CACTUS answer, in memory and on disk.

        The same as ``provesid.clear_cache(service='nci')``.

        Examples:
            >>> resolver = NCIChemicalIdentifierResolver()
            >>> resolver.clear_cache()
            >>> resolver.get_cache_info()['file_count']
            0
        """
        from .cache import clear_cache
        clear_cache(service='nci')

    def get_cache_info(self) -> Dict[str, Any]:
        """
        Size and location of the CACTUS cache.

        Returns:
            The statistics
            [`provesid.cache.get_cache_info`][provesid.cache.get_cache_info]
            reports for
                ``service='nci'``.

        Examples:
            >>> NCIChemicalIdentifierResolver().get_cache_info()['cache_directory'].endswith('nci')
            True
        """
        from .cache import get_cache_info
        return get_cache_info(service='nci')

    def _rate_limit(self):
        """
        Sleep, if needed, so this client's requests stay ``pause_time`` apart.

        Delegates to the shared transport, which paces every request it makes
        including retries.

        Examples:
            >>> NCIChemicalIdentifierResolver(pause_time=0)._rate_limit()
        """
        self._http.rate_limit()

    @property
    def pause_time(self) -> float:
        """
        Minimum seconds between two requests from this client.

        Held by the transport, and settable: raising it slows a long batch
        down, and the new value takes effect on the next request, retries
        included.

        Examples:
            >>> r = NCIChemicalIdentifierResolver()
            >>> r.pause_time
            0.1
            >>> r.pause_time = 1.0   # gentler, for a long batch
        """
        return self._http.min_interval

    @pause_time.setter
    def pause_time(self, seconds: float) -> None:
        self._http.min_interval = seconds

    @property
    def last_request_time(self) -> float:
        """
        When this client last made a request, as a Unix timestamp.

        Kept on the transport; exposed here because it describes the client's
        own pacing.

        Returns:
            Seconds since the epoch, or 0.0 before the first request.

        Examples:
            >>> NCIChemicalIdentifierResolver().last_request_time
            0.0
        """
        return self._http.last_request_time

    def _make_request(self, url: str) -> str:
        """
        Make one request to the resolver and return its body as text.

        The transport handles the pacing, the retries and the mapping of
        status codes onto this module's exceptions; this method exists so the
        service's one request shape has a name.

        Args:
            url: Request URL

        Returns:
            Response text, stripped.

        Raises:
            NCIResolverTimeoutError: Every attempt timed out or failed to
                connect.
            NCIResolverNotFoundError: The identifier could not be resolved
                (404).
            NCIResolverError: Any other error, including a 429 or 5xx that
                outlived the retry budget.
        """
        return self._http.get_text(url)

    def _build_url(self, identifier: str, representation: str, xml_format: bool = False) -> str:
        """
        Build URL for the NCI resolver service

        Args:
            identifier: Chemical structure identifier
            representation: Desired representation
            xml_format: Whether to request XML format

        Returns:
            Complete URL string

        Examples:
            >>> NCIChemicalIdentifierResolver()._build_url("C/C=C/C", "stdinchikey")
            'https://cactus.nci.nih.gov/chemical/structure/C/C%3DC/C/stdinchikey'
        """
        # URL-encode the identifier for special characters, except "/".
        # CACTUS's front end answers 404 to a path holding %2F, which made
        # every InChI --- and every SMILES with a stereo bond --- look
        # unresolvable. A literal slash reaches the resolver intact.
        encoded_identifier = quote(identifier, safe='/')

        # Build URL components
        url_parts = [self.base_url, encoded_identifier, representation]

        if xml_format:
            url_parts.append('xml')

        return '/'.join(url_parts)

    @cached(service='nci')
    def resolve(self, identifier: str, representation: str, xml_format: bool = False) -> str:
        """
        Resolve a chemical identifier to another representation

        Args:
            identifier: Input chemical identifier (name, SMILES, InChI, CAS, etc.)
            representation: Target representation (see self.representations for options)
            xml_format: Whether to request XML format response

        Returns:
            Resolved representation as string, as CACTUS wrote it: several
            values are one per line.

        Raises:
            ValueError: If representation is not supported
            NCIResolverNotFoundError: If identifier cannot be resolved
            NCIResolverError: For other resolver errors, and for an empty
                identifier

        Examples:
            >>> resolver = NCIChemicalIdentifierResolver()
            >>> resolver.resolve("50-00-0", "smiles")              # doctest: +SKIP
            'C=O'
            >>> resolver.resolve("aspirin", "cas").splitlines()[:2]  # doctest: +SKIP
            ['50-78-2', '11126-35-5']
            >>> resolver.resolve("CCO", "nope")
            Traceback (most recent call last):
            ...
            ValueError: Unsupported representation 'nope'. Available: stdinchi, ...
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

        Examples:
            >>> NCIChemicalIdentifierResolver().get_available_representations()[:4]
            ['stdinchi', 'stdinchikey', 'smiles', 'ficts']
        """
        return list(self.representations.keys())

    @cached(service='nci')
    def resolve_multiple(self, identifier: str, representations: List[str]) -> Dict[str, str]:
        """
        Resolve a single identifier to multiple representations

        Args:
            identifier: Input chemical identifier
            representations: List of target representations

        Returns:
            Dictionary mapping representation to resolved value; None for a
            representation that failed, with a warning

        Examples:
            >>> NCIChemicalIdentifierResolver().resolve_multiple(
            ...     "ethanol", ["formula", "mw"])                     # doctest: +SKIP
            {'formula': 'C2H6O', 'mw': '46.0688'}
        """
        results = {}
        for representation in representations:
            try:
                results[representation] = self.resolve(identifier, representation)
            except NCIResolverError as e:
                results[representation] = None
                self.logger.warning(f"Failed to resolve {identifier} to {representation}: {e}")

        return results

    # version=2: cas became a list and stdinchikey lost its InChIKey= prefix.
    @cached(service='nci', version=2)
    def get_molecular_data(self, identifier: str) -> Dict[str, Any]:
        """
        Get comprehensive molecular data for a chemical identifier

        This method attempts to retrieve multiple common properties and identifiers
        for a given chemical, similar to the original nci_cas_to_mol function.

        Unlike [`resolve`][provesid.resolver.NCIChemicalIdentifierResolver.resolve],
        this parses CACTUS's text: ``names`` and ``cas`` become lists,
        ``stdinchikey`` loses its ``InChIKey=`` prefix and ``mw`` becomes a
        float.

        Args:
            identifier: Input chemical identifier

        Returns:
            Dictionary with molecular data and metadata: ``found_by``,
            ``success`` (False only when nothing resolved), ``error``,
            ``available_data`` and, at the top level for convenience, each
            representation (``stdinchi``, ``stdinchikey``, ``smiles``,
            ``names`` as a list, ``iupac_name``, ``cas`` as a list, ``mw`` as
            a float, ``formula``, ``ficts``, ``ficus``, ``uuuuu``,
            ``hashisy``), None where it failed. Twelve requests.

            ``cas`` is every CAS number CACTUS associates with the structure,
            in its order, which is not a ranking: ethanol's list starts with
            121182-78-3, not 64-17-5. Treat it as a set of candidates.

        Examples:
            >>> data = NCIChemicalIdentifierResolver().get_molecular_data("64-17-5")  # doctest: +SKIP
            >>> data["formula"], data["mw"], data["stdinchikey"]                      # doctest: +SKIP
            ('C2H6O', 46.0688, 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N')
            >>> "64-17-5" in data["cas"]                                              # doctest: +SKIP
            True
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
                if rep in ('names', 'cas'):
                    result['available_data'][rep] = _lines(value)
                elif rep == 'stdinchikey':
                    result['available_data'][rep] = _strip_inchikey_prefix(value)
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
                self.logger.debug(f"Could not resolve {identifier} to {rep}: {e}")

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
            URL for the structure image. Builds the URL only; nothing is
            requested.

        Examples:
            >>> NCIChemicalIdentifierResolver().get_image_url("aspirin", "png", 300, 300)
            'https://cactus.nci.nih.gov/chemical/structure/aspirin/image?format=png&width=300&height=300'
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

        Examples:
            >>> r = NCIChemicalIdentifierResolver()
            >>> r.download_image('aspirin', 'aspirin.png', 'png')   # doctest: +SKIP
            True
        """
        try:
            image_url = self.get_image_url(identifier, image_format, width, height)
            response = self._http.get(image_url)

            with open(filename, 'wb') as f:
                f.write(response.content)

            return True

        except Exception as e:
            self.logger.error(f"Failed to download image for {identifier}: {e}")
            return False

    @cached(service='nci')
    def batch_resolve(self, identifiers: List[str], representation: str) -> Dict[str, str]:
        """
        Resolve multiple identifiers to a single representation

        Args:
            identifiers: List of chemical identifiers
            representation: Target representation

        Returns:
            Dictionary mapping identifier to resolved value (None if failed)

        Examples:
            >>> NCIChemicalIdentifierResolver().batch_resolve(
            ...     ["ethanol", "methanol"], "formula")              # doctest: +SKIP
            {'ethanol': 'C2H6O', 'methanol': 'CH4O'}
        """
        results = {}

        for identifier in identifiers:
            try:
                results[identifier] = self.resolve(identifier, representation)
            except NCIResolverError as e:
                results[identifier] = None
                self.logger.warning(f"Failed to resolve {identifier}: {e}")

        return results

    @cached(service='nci')
    def is_valid_identifier(self, identifier: str) -> bool:
        """
        Check if an identifier can be resolved by the service

        Args:
            identifier: Chemical identifier to test

        Returns:
            True if identifier can be resolved, False otherwise. A failed
            request also reads as False.

        Examples:
            >>> resolver = NCIChemicalIdentifierResolver()
            >>> resolver.is_valid_identifier("aspirin")                 # doctest: +SKIP
            True
            >>> resolver.is_valid_identifier("this_is_definitely_not_a_chemical_12345")  # doctest: +SKIP
            False
        """
        try:
            # Try to get SMILES as a basic test
            self.resolve(identifier, 'smiles')
            return True
        except NCIResolverError:
            return False

    @cached(service='nci')
    def search_by_partial_name(self, partial_name: str) -> List[str]:
        """
        Search for compounds by partial name match
        Note: This is a basic implementation - the NCI service doesn't have
        a dedicated partial matching endpoint, so this tries the exact name first.

        Args:
            partial_name: Partial chemical name

        Returns:
            List of matching names (may be empty): the ``names`` of the
            compound the text resolves to as a whole

        Examples:
            >>> NCIChemicalIdentifierResolver().search_by_partial_name("aspirin")[:2]  # doctest: +SKIP
            ['2-acetyloxybenzoic acid', '2-Acetoxybenzoic acid']
        """
        try:
            names = self.resolve(partial_name, 'names')
            return [name.strip() for name in names.split('\n') if name.strip()]
        except NCIResolverError:
            return []

# Convenience functions for backwards compatibility and ease of use

@cached(service='nci', version=2)
def nci_cas_to_mol(cas_rn: str) -> Dict[str, Any]:
    """
    Convert a CAS RN to a molecule data structure using the NCI web API

    This function maintains compatibility with the original nci_cas_to_mol function
    while using the new NCIChemicalIdentifierResolver class.

    Args:
        cas_rn: CAS Registry Number

    Returns:
        Dictionary with molecular data, as
        [`NCIChemicalIdentifierResolver.get_molecular_data`][provesid.resolver.NCIChemicalIdentifierResolver.get_molecular_data]
        returns it

    Examples:
        >>> nci_cas_to_mol("64-17-5")["formula"]                 # doctest: +SKIP
        'C2H6O'
    """
    resolver = NCIChemicalIdentifierResolver()
    return resolver.get_molecular_data(cas_rn)

@cached(service='nci', version=2)
def nci_id_to_mol(identifier: str) -> Dict[str, Any]:
    """
    Convert any chemical identifier to a molecule data structure

    Args:
        identifier: Chemical identifier (CAS, name, SMILES, InChI, etc.)

    Returns:
        Dictionary with molecular data, as
        [`NCIChemicalIdentifierResolver.get_molecular_data`][provesid.resolver.NCIChemicalIdentifierResolver.get_molecular_data]
        returns it

    Examples:
        >>> nci_id_to_mol("CCO")["stdinchikey"]                  # doctest: +SKIP
        'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'
    """
    resolver = NCIChemicalIdentifierResolver()
    return resolver.get_molecular_data(identifier)

@cached(service='nci')
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

    Raises:
        ValueError: If ``output_type`` is not a supported representation.

    Examples:
        >>> nci_resolver("ethanol", "stdinchikey")                # doctest: +SKIP
        'InChIKey=LFQSCWFLJHTTHZ-UHFFFAOYSA-N'
    """
    try:
        resolver = NCIChemicalIdentifierResolver(timeout=timeout)
        return resolver.resolve(input_value, output_type)
    except NCIResolverError:
        return None

@cached(service='nci')
def nci_smiles_to_names(smiles: str) -> List[str]:
    """
    Get chemical names for a SMILES string

    Args:
        smiles: SMILES string

    Returns:
        List of chemical names, CAS numbers among them; empty if not found

    Examples:
        >>> nci_smiles_to_names("CCO")[:3]                        # doctest: +SKIP
        ['ethanol', '121182-78-3', '64-17-5']
    """
    try:
        resolver = NCIChemicalIdentifierResolver()
        names_str = resolver.resolve(smiles, 'names')
        return [name.strip() for name in names_str.split('\n') if name.strip()]
    except NCIResolverError:
        return []

@cached(service='nci')
def nci_name_to_smiles(name: str) -> Optional[str]:
    """
    Convert chemical name to SMILES

    Args:
        name: Chemical name

    Returns:
        SMILES string or None if not found

    Examples:
        >>> nci_name_to_smiles("caffeine")                        # doctest: +SKIP
        'Cn1cnc2N(C)C(=O)N(C)C(=O)c12'
    """
    try:
        resolver = NCIChemicalIdentifierResolver()
        return resolver.resolve(name, 'smiles')
    except NCIResolverError:
        return None

@cached(service='nci')
def nci_inchi_to_smiles(inchi: str) -> Optional[str]:
    """
    Convert InChI to SMILES

    Args:
        inchi: InChI string

    Returns:
        SMILES string or None if not found

    Examples:
        >>> nci_inchi_to_smiles("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")  # doctest: +SKIP
        'CCO'
    """
    try:
        resolver = NCIChemicalIdentifierResolver()
        return resolver.resolve(inchi, 'smiles')
    except NCIResolverError:
        return None

@cached(service='nci')
def nci_cas_to_inchi(cas_rn: str) -> Optional[str]:
    """
    Convert CAS Registry Number to Standard InChI

    Args:
        cas_rn: CAS Registry Number

    Returns:
        Standard InChI string or None if not found

    Examples:
        >>> nci_cas_to_inchi("50-00-0")                           # doctest: +SKIP
        'InChI=1S/CH2O/c1-2/h1H2'
    """
    try:
        resolver = NCIChemicalIdentifierResolver()
        return resolver.resolve(cas_rn, 'stdinchi')
    except NCIResolverError:
        return None

@cached(service='nci')
def nci_get_molecular_weight(identifier: str) -> Optional[float]:
    """
    Get molecular weight for any chemical identifier

    Args:
        identifier: Chemical identifier

    Returns:
        Molecular weight as float or None if not found

    Examples:
        >>> nci_get_molecular_weight("caffeine")                  # doctest: +SKIP
        194.1926
    """
    try:
        resolver = NCIChemicalIdentifierResolver()
        mw_str = resolver.resolve(identifier, 'mw')
        return float(mw_str)
    except (NCIResolverError, ValueError):
        return None

@cached(service='nci')
def nci_get_formula(identifier: str) -> Optional[str]:
    """
    Get molecular formula for any chemical identifier

    Args:
        identifier: Chemical identifier

    Returns:
        Molecular formula string or None if not found

    Examples:
        >>> nci_get_formula("caffeine")                           # doctest: +SKIP
        'C8H10N4O2'
    """
    try:
        resolver = NCIChemicalIdentifierResolver()
        return resolver.resolve(identifier, 'formula')
    except NCIResolverError:
        return None
