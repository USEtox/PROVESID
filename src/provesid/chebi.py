"""
ChEBI (Chemical Entities of Biological Interest) API interface.

This module provides a Python interface to the ChEBI 2.0 REST API for retrieving
chemical compound information from the ChEBI database. Offline access through
ChEBI's SDF release is :class:`~provesid.chebi_sdf.ChebiSDF`, in
:mod:`provesid.chebi_sdf`.

API documentation: https://www.ebi.ac.uk/chebi/backend/api/docs/

Author: USEtox team
Date: August 2025
"""

import requests
import logging
import time
from typing import Dict, List, Optional, Union, Any
from .http import HTTPClient, NotFoundError, ServiceError, ServiceTimeoutError


#: Seconds between two requests to ChEBI from this process. EBI publishes no
#: per-IP figure for the ChEBI 2.0 API, so this is politeness rather than a
#: quoted limit: ten requests a second is well inside what a walk over an
#: ontology subtree needs.
CHEBI_MIN_INTERVAL = 0.1


class ChEBIError(ServiceError):
    """
    Custom exception for ChEBI API errors.

    Also the base for the two more specific failures below, so a caller that
    catches this one keeps catching everything --- which is what every method
    in this module does internally before returning None.

    Examples:
        >>> issubclass(ChEBINotFoundError, ChEBIError)
        True
    """
    pass


class ChEBINotFoundError(ChEBIError, NotFoundError):
    """
    ChEBI answered, and its answer was that there is no such record.

    A statement about the data rather than the request, so it is never retried.

    Examples:
        >>> issubclass(ChEBINotFoundError, NotFoundError)
        True
    """
    pass


class ChEBITimeoutError(ChEBIError, ServiceTimeoutError):
    """
    Every attempt timed out or the connection could not be made.

    Examples:
        >>> issubclass(ChEBITimeoutError, ServiceTimeoutError)
        True
    """
    pass


class ChEBI:
    """
    Interface for the ChEBI 2.0 (Chemical Entities of Biological Interest) REST API.

    The ChEBI database is a freely available dictionary of molecular entities
    focused on 'small' chemical compounds. This class provides methods to search
    for and retrieve compound information from the ChEBI 2.0 API.

    Attributes:
        base_url (str): Base URL for ChEBI 2.0 API
        timeout (int): Request timeout in seconds
        session (requests.Session): HTTP session for connection pooling

    Every method needs the network. Lookups that fail return None (or an
    empty collection) after logging a warning, rather than raising; the
    offline alternative for records, names and structures is
    :class:`~provesid.chebi_sdf.ChebiSDF`.

    Examples:
        >>> chebi = ChEBI()
        >>> compound = chebi.get_compound(15377)  # doctest: +SKIP
        >>> print(compound['name'])               # doctest: +SKIP
        water
    """

    # Valid ontology relation types for the ChEBI 2.0 API
    VALID_RELATIONS = [
        "has_functional_parent", "has_parent_hydride", "has_part", "has_role",
        "is_a", "is_conjugate_acid_of", "is_conjugate_base_of",
        "is_enantiomer_of", "is_part_of", "is_substituent_group_from",
        "is_tautomer_of",
    ]

    # Valid structure search types
    VALID_SEARCH_TYPES = ["connectivity", "similarity", "substructure"]

    def __init__(self, timeout: int = 30):
        """
        Initialize ChEBI 2.0 API client.

        Args:
            timeout (int): Request timeout in seconds (default: 30)

        Examples:
            >>> ChEBI(timeout=60)
            ChEBI(base_url='https://www.ebi.ac.uk/chebi/backend/api/public', timeout=60)
        """
        self.base_url = "https://www.ebi.ac.uk/chebi/backend/api/public"
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'PROVESID-ChEBI-Client/2.0',
            'Accept': 'application/json',
        })

        # Setup logging
        self.logger = logging.getLogger(__name__)

        # One shared transport, making its calls through the session above so
        # that the pooled connection and the persistent headers survive the
        # move. ChEBI uses its status codes honestly --- unlike PubChem and
        # CACTUS, it has no fault code in the body to read --- so the default
        # classifier is the right one.
        #
        # Two retries with a half-second base, rather than the transport's
        # three and one: an ontology walk makes many small requests to a fast
        # service, where a long back-off curve costs more than the request it
        # is protecting.
        self._http = HTTPClient(
            session=self.session,
            min_interval=CHEBI_MIN_INTERVAL,
            timeout=timeout,
            max_retries=2,
            backoff=0.5,
            error_cls=ChEBIError,
            not_found_cls=ChEBINotFoundError,
            timeout_cls=ChEBITimeoutError,
            pace_host=self.base_url,
            logger=self.logger,
        )

    @staticmethod
    def _format_chebi_id(chebi_id: Union[int, str]) -> str:
        """
        Format a ChEBI ID to the canonical ``CHEBI:<number>`` form.

        The ChEBI 2.0 API accepts IDs with or without the prefix, but this
        helper ensures consistency.

        Args:
            chebi_id: ChEBI ID (int, bare number string, or ``CHEBI:…`` string)

        Returns:
            str: ID in ``CHEBI:<number>`` form
        """
        chebi_id_str = str(chebi_id).strip()
        if not chebi_id_str.upper().startswith("CHEBI:"):
            chebi_id_str = f"CHEBI:{chebi_id_str}"
        return chebi_id_str

    # ------------------------------------------------------------------
    # Low-level HTTP helpers
    # ------------------------------------------------------------------

    def _json_or_text(self, response: requests.Response) -> Any:
        """
        Return a response body as JSON when ChEBI says it is JSON, else as text.

        ChEBI serves molfiles and SVG from the same API as its records, and
        labels them honestly in ``Content-Type``, so the header is what decides.

        Args:
            response: The response to read.

        Returns:
            The decoded JSON, or the body as text.

        Raises:
            ChEBIError: The body was labelled JSON and is not.
        """
        if response.headers.get("Content-Type", "").startswith("application/json"):
            return self._http.decode_json(response)
        return response.text

    def _get(self, endpoint: str, params: Optional[Dict] = None) -> Any:
        """
        Perform a GET request and return the parsed JSON body.

        Pacing, retries and back-off belong to the shared transport; what stays
        here is ChEBI's own contract --- where its endpoints live and how it
        labels what it returns.

        Args:
            endpoint (str): Path relative to *base_url* (e.g. ``compound/15377/``).
            params (dict, optional): Query-string parameters.

        Returns:
            Parsed JSON response (dict / list / str).

        Raises:
            ChEBINotFoundError: ChEBI reported no such record (404).
            ChEBITimeoutError: Every attempt timed out or could not connect.
            ChEBIError: Any other HTTP or JSON failure, including a transient
                one that outlived the retries.
        """
        url = f"{self.base_url}/{endpoint}"
        return self._json_or_text(self._http.get(url, params=params))

    def _get_raw(self, endpoint: str, params: Optional[Dict] = None) -> requests.Response:
        """
        Perform a GET request and return the raw :class:`requests.Response`.

        Useful for endpoints that return non-JSON content (SVG, molfile, images).

        Args:
            endpoint (str): Path relative to *base_url*.
            params (dict, optional): Query-string parameters.

        Returns:
            The response, already classified as a success.

        Raises:
            ChEBINotFoundError: ChEBI reported no such record (404).
            ChEBITimeoutError: Every attempt timed out or could not connect.
            ChEBIError: Any other HTTP failure.
        """
        url = f"{self.base_url}/{endpoint}"
        return self._http.get(url, params=params)

    def _post_json(self, endpoint: str, json_body: Any = None,
                   params: Optional[Dict] = None) -> Any:
        """
        Perform a POST request with a JSON body and return parsed JSON.

        Args:
            endpoint (str): Path relative to *base_url*.
            json_body: Object to serialise as JSON request body.
            params (dict, optional): Query-string parameters.

        Returns:
            Parsed JSON response.

        Raises:
            ChEBINotFoundError: ChEBI reported no such record (404).
            ChEBITimeoutError: Every attempt timed out or could not connect.
            ChEBIError: Any other HTTP or JSON failure.
        """
        url = f"{self.base_url}/{endpoint}"
        return self._json_or_text(
            self._http.post(url, json=json_body, params=params)
        )

    def _post_text(self, endpoint: str, text_body: str,
                   params: Optional[Dict] = None) -> str:
        """
        Perform a POST request with a ``text/plain`` body and return the
        response text.  Used by the structure-calculation endpoints.

        Args:
            endpoint (str): Path relative to *base_url*.
            text_body (str): The body to send, a SMILES string or a formula.
            params (dict, optional): Query-string parameters.

        Returns:
            The response body as text.

        Raises:
            ChEBINotFoundError: ChEBI reported no such record (404).
            ChEBITimeoutError: Every attempt timed out or could not connect.
            ChEBIError: Any other HTTP failure.
        """
        return self._post_raw(endpoint, text_body, params=params).text

    def _post_raw(self, endpoint: str, text_body: str,
                  params: Optional[Dict] = None) -> requests.Response:
        """
        POST a ``text/plain`` body and return the raw response.

        The structure-calculation endpoints want the text; ``depict_structure``
        wants the bytes of a PNG. Two headers are set per request, and
        ``requests`` merges them over the session's own: ``Content-Type``, and
        ``Accept``. These endpoints answer ``text/plain`` or ``image/png``, and
        refuse the session's ``Accept: application/json`` with HTTP 406, which
        until 2026-09-22 made every one of them return None.

        Args:
            endpoint (str): Path relative to *base_url*.
            text_body (str): The body to send.
            params (dict, optional): Query-string parameters.

        Returns:
            The response, already classified as a success.

        Raises:
            ChEBINotFoundError: ChEBI reported no such record (404).
            ChEBITimeoutError: Every attempt timed out or could not connect.
            ChEBIError: Any other HTTP failure.
        """
        url = f"{self.base_url}/{endpoint}"
        return self._http.post(
            url, data=text_body, params=params,
            headers={"Content-Type": "text/plain;charset=UTF-8", "Accept": "*/*"},
        )

    # ------------------------------------------------------------------
    # Compound retrieval
    # ------------------------------------------------------------------

    def get_compound(self, chebi_id: Union[int, str], *,
                     only_ontology_parents: bool = False,
                     only_ontology_children: bool = False) -> Optional[Dict[str, Any]]:
        """
        Retrieve information about a single compound.

        Endpoint: ``GET /compound/{chebi_id}/``

        Args:
            chebi_id: ChEBI ID (with or without ``CHEBI:`` prefix).
            only_ontology_parents: If *True*, return only ontology parents.
            only_ontology_children: If *True*, return only ontology children.

        Returns:
            dict with compound data, or *None* if not found.

        Examples:
            >>> chebi = ChEBI()
            >>> water = chebi.get_compound(15377)
        """
        chebi_id_str = self._format_chebi_id(chebi_id)
        params: Dict[str, Any] = {}
        if only_ontology_parents:
            params["only_ontology_parents"] = "true"
        if only_ontology_children:
            params["only_ontology_children"] = "true"

        try:
            return self._get(f"compound/{chebi_id_str}/", params=params or None)
        except ChEBIError as e:
            self.logger.warning(f"Failed to get compound {chebi_id_str}: {e}")
            return None

    def get_complete_entity(
        self,
        chebi_id: Union[int, str],
        *,
        only_ontology_parents: bool = False,
        only_ontology_children: bool = False,
    ) -> Optional[Dict[str, Any]]:
        """Backward-compatible alias for :meth:`get_compound`.

        Older examples and user code refer to ``get_complete_entity``. The
        ChEBI 2.0 client uses ``get_compound`` as the canonical method name.

        Args:
            chebi_id: ChEBI ID (with or without ``CHEBI:`` prefix).
            only_ontology_parents: As for :meth:`get_compound`.
            only_ontology_children: As for :meth:`get_compound`.

        Returns:
            The :meth:`get_compound` record, or None on error.

        Examples:
            >>> ChEBI().get_complete_entity("CHEBI:15377")["name"]  # doctest: +SKIP
            'water'
        """
        return self.get_compound(
            chebi_id,
            only_ontology_parents=only_ontology_parents,
            only_ontology_children=only_ontology_children,
        )

    def get_compounds(self, chebi_ids: List[Union[int, str]]) -> Optional[Any]:
        """
        Retrieve information about one or more compounds in a single call.

        Endpoint: ``POST /compounds/``

        Args:
            chebi_ids: List of ChEBI IDs.

        Returns:
            API response, a dict keyed by ``CHEBI:<id>`` holding each compound's
            record, or *None* on error.

        Examples:
            >>> chebi = ChEBI()
            >>> results = chebi.get_compounds(["CHEBI:15377", "CHEBI:15365"])  # doctest: +SKIP
            >>> {key: record["name"] for key, record in results.items()}      # doctest: +SKIP
            {'CHEBI:15377': 'water', 'CHEBI:15365': 'acetylsalicylic acid'}
        """
        ids_formatted = [self._format_chebi_id(cid) for cid in chebi_ids]
        try:
            return self._post_json("compounds/", json_body={"chebi_ids": ids_formatted})
        except ChEBIError as e:
            self.logger.warning(f"Failed to get compounds: {e}")
            return None

    def batch_get_entities(
        self,
        chebi_ids: List[Union[int, str]],
        pause_time: float = 0.0,
    ) -> Dict[str, Dict[str, Any]]:
        """Backward-compatible batch helper for tutorial workflows.

        Args:
            chebi_ids: List of ChEBI IDs (with or without CHEBI: prefix).
            pause_time: Optional sleep time between requests in seconds.

        Returns:
            Mapping of canonical ``CHEBI:<id>`` string to compound payload.
            IDs that could not be fetched are left out.

        Examples:
            >>> found = ChEBI().batch_get_entities([15377, "CHEBI:15365"])  # doctest: +SKIP
            >>> {key: record["name"] for key, record in found.items()}     # doctest: +SKIP
            {'CHEBI:15377': 'water', 'CHEBI:15365': 'acetylsalicylic acid'}
        """
        results: Dict[str, Dict[str, Any]] = {}
        for chebi_id in chebi_ids:
            canonical_id = self._format_chebi_id(chebi_id)
            entity = self.get_compound(chebi_id)
            if entity:
                results[canonical_id] = entity
            if pause_time > 0:
                time.sleep(pause_time)
        return results

    # ------------------------------------------------------------------
    # Search
    # ------------------------------------------------------------------

    def search(self, term: str, *, page: int = 1, size: int = 15) -> Optional[Any]:
        """
        General text search (Elasticsearch-backed).

        Endpoint: ``GET /es_search/``

        You can search by ChEBI name, brand name, IUPAC name, synonym,
        InChIKey, formula, SMILES, InChI, CAS number, database cross-reference
        IDs, PubMed IDs, and more.

        Args:
            term: Search term (name, SMILES, InChIKey, CAS, formula, …).
            page: Page number for pagination (default 1).
            size: Page size (default 15).

        Returns:
            Search results dict, or *None* on error: ``results`` (the hits,
            best first, each with its record under ``_source``), ``total`` and
            ``number_pages``.

        Examples:
            >>> chebi = ChEBI()
            >>> results = chebi.search("paracetamol")               # doctest: +SKIP
            >>> results["results"][0]["_source"]["chebi_accession"]  # doctest: +SKIP
            'CHEBI:46195'
        """
        params: Dict[str, Any] = {"term": term, "page": page, "size": size}
        try:
            return self._get("es_search/", params=params)
        except ChEBIError as e:
            self.logger.warning(f"Search failed for '{term}': {e}")
            return None

    def search_by_name(self, search_text: str, *, page: int = 1,
                       size: int = 15) -> List[Dict[str, Any]]:
        """
        Search ChEBI by compound name (convenience wrapper around :meth:`search`).

        Args:
            search_text: Text to search for.
            page: Page number for pagination (default 1).
            size: Page size (default 15).

        Returns:
            list of matching entity dicts (may be empty): the ``results`` of
            :meth:`search`, each carrying ``_score`` and the record under
            ``_source``.

        Examples:
            >>> hits = ChEBI().search_by_name("paracetamol", size=2)  # doctest: +SKIP
            >>> [(hit["_source"]["chebi_accession"], hit["_source"]["name"]) for hit in hits]  # doctest: +SKIP
            [('CHEBI:46195', 'paracetamol'), ('CHEBI:74529', 'antidote to paracetamol poisoning')]
        """
        result = self.search(search_text, page=page, size=size)
        if result is None:
            return []
        # The es_search endpoint typically returns a dict with a results list
        if isinstance(result, dict):
            return result.get("results", result.get("items", [result]))
        if isinstance(result, list):
            return result
        return []

    # ------------------------------------------------------------------
    # Advanced search
    # ------------------------------------------------------------------

    def advanced_search(
        self,
        specification: Dict[str, Any],
        *,
        three_star_only: bool = True,
        has_structure: Optional[bool] = None,
        page: int = 1,
        size: int = 15,
        download: bool = False,
    ) -> Optional[Any]:
        """
        Perform an advanced compound search using specification objects.

        Endpoint: ``POST /advanced_search/``

        Seven specification types are supported (combine with ``and_specification``,
        ``or_specification``, ``but_not_specification``):

        - ``ontology_specification``
        - ``formula_specification``
        - ``mass_specification``
        - ``monoisotopicmass_specification``
        - ``charge_specification``
        - ``database_name_specification``
        - ``text_search_specification``

        Args:
            specification: Request body following the ``FullSpecification`` schema.
            three_star_only: Only include 3-star compounds (default *True*).
            has_structure: Filter by structure availability (*None* = no filter).
            page: Page number (default 1).
            size: Page size (default 15).
            download: Return results in download format (default *False*).

        Returns:
            API response dict/list, or *None* on error.

        Examples:
            >>> chebi = ChEBI()
            >>> results = chebi.advanced_search({
            ...     "formula_specification": {
            ...         "and_specification": [{"term": "C6H12O7"}]
            ...     }
            ... }, three_star_only=False)                         # doctest: +SKIP
            >>> results["total"]                                  # doctest: +SKIP
            15
        """
        params: Dict[str, Any] = {
            "three_star_only": str(three_star_only).lower(),
            "page": page,
            "size": size,
            "download": str(download).lower(),
        }
        if has_structure is not None:
            params["has_structure"] = str(has_structure).lower()

        try:
            return self._post_json("advanced_search/", json_body=specification,
                                   params=params)
        except ChEBIError as e:
            self.logger.warning(f"Advanced search failed: {e}")
            return None

    def get_sources_list(self) -> Optional[Any]:
        """
        Retrieve the list of available database sources for advanced search.

        Endpoint: ``GET /advanced_search/sources_list``

        Returns:
            Sources information, or *None* on error: a list of
            ``{"prefix", "name"}`` dicts, the names usable in a
            ``database_name_specification``.

        Examples:
            >>> ChEBI().get_sources_list()[:2]                    # doctest: +SKIP
            [{'prefix': 'agr', 'name': 'Agricola'}, {'prefix': 'pesticides', 'name': "Alan Wood's Pesticides"}]
        """
        try:
            return self._get("advanced_search/sources_list")
        except ChEBIError as e:
            self.logger.warning(f"Failed to get sources list: {e}")
            return None

    # ------------------------------------------------------------------
    # Ontology
    # ------------------------------------------------------------------

    def get_ontology_parents(self, chebi_id: Union[int, str]) -> Optional[Any]:
        """
        Get the ontology parents of a compound.

        Endpoint: ``GET /ontology/parents/{chebi_id}/``

        Args:
            chebi_id: ChEBI ID (with or without ``CHEBI:`` prefix).

        Returns:
            Ontology parent data, or *None* on error: the compound's ``id`` and
            ``chebi_accession``, with its relations under
            ``ontology_relations["outgoing_relations"]``, each naming the
            ``relation_type`` and the parent (``final_id``, ``final_name``).

        Examples:
            >>> parents = ChEBI().get_ontology_parents(15377)     # doctest: +SKIP
            >>> first = parents["ontology_relations"]["outgoing_relations"][0]  # doctest: +SKIP
            >>> first["relation_type"], first["final_id"]         # doctest: +SKIP
            ('has role', 75772)
        """
        chebi_id_str = self._format_chebi_id(chebi_id)
        try:
            return self._get(f"ontology/parents/{chebi_id_str}/")
        except ChEBIError as e:
            self.logger.warning(f"Failed to get ontology parents for {chebi_id_str}: {e}")
            return None

    def get_ontology_children(self, chebi_id: Union[int, str]) -> Optional[Any]:
        """
        Get the ontology children of a compound.

        Endpoint: ``GET /ontology/children/{chebi_id}/``

        Args:
            chebi_id: ChEBI ID (with or without ``CHEBI:`` prefix).

        Returns:
            Ontology children data, or *None* on error: as for
            :meth:`get_ontology_parents`, with the relations under
            ``ontology_relations["incoming_relations"]`` and the child named by
            ``init_id`` and ``init_name``.

        Examples:
            >>> children = ChEBI().get_ontology_children(15377)   # doctest: +SKIP
            >>> children["ontology_relations"]["incoming_relations"][0]["init_name"]  # doctest: +SKIP
            'methane clathrate'
        """
        chebi_id_str = self._format_chebi_id(chebi_id)
        try:
            return self._get(f"ontology/children/{chebi_id_str}/")
        except ChEBIError as e:
            self.logger.warning(f"Failed to get ontology children for {chebi_id_str}: {e}")
            return None

    def get_all_ontology_children_in_path(
        self,
        relation: str,
        entity: Union[int, str],
        *,
        three_star_only: bool = True,
        has_structure: Optional[bool] = None,
        page: int = 1,
        size: int = 15,
        download: bool = False,
    ) -> Optional[Any]:
        """
        Search all compounds in the ontology matching a relation and entity.

        Endpoint: ``GET /ontology/all_children_in_path/``

        Args:
            relation: Ontology relation type (e.g. ``is_a``, ``has_role``).
            entity: ChEBI ID of the entity to find children of.
            three_star_only: Only include 3-star compounds (default *True*).
            has_structure: Filter by structure availability.
            page: Page number (default 1).
            size: Page size (default 15).
            download: Return results in download format.

        Returns:
            API response, or *None* on error.

        Examples:
            >>> chebi = ChEBI()
            >>> # Get all compounds that are alcohols
            >>> results = chebi.get_all_ontology_children_in_path(
            ...     relation="is_a", entity="CHEBI:30879"
            ... )                                                 # doctest: +SKIP
            >>> results["total"]                                  # doctest: +SKIP
            4029
        """
        entity_str = self._format_chebi_id(entity)
        params: Dict[str, Any] = {
            "relation": relation,
            "entity": entity_str,
            "three_star_only": str(three_star_only).lower(),
            "page": page,
            "size": size,
            "download": str(download).lower(),
        }
        if has_structure is not None:
            params["has_structure"] = str(has_structure).lower()

        try:
            return self._get("ontology/all_children_in_path/", params=params)
        except ChEBIError as e:
            self.logger.warning(f"Failed to get all ontology children: {e}")
            return None

    # ------------------------------------------------------------------
    # Structures
    # ------------------------------------------------------------------

    def get_compound_structure(self, chebi_id: Union[int, str], *,
                               width: int = 300, height: int = 300) -> Optional[str]:
        """
        Get the default SVG structure for a compound.

        Endpoint: ``GET /compound/{id}/structure/``

        Args:
            chebi_id: ChEBI ID (numeric, the primary key of the compound).
            width: Width of the SVG (default 300).
            height: Height of the SVG (default 300).

        Returns:
            Raw SVG string, or *None* on error.

        Examples:
            >>> ChEBI().get_compound_structure(15377)[:21]        # doctest: +SKIP
            "<?xml version='1.0' e"
        """
        # This endpoint expects the numeric compound PK
        numeric_id = str(chebi_id)
        if numeric_id.upper().startswith("CHEBI:"):
            numeric_id = numeric_id.split(":")[-1]

        params: Dict[str, Any] = {"width": width, "height": height}
        try:
            resp = self._get_raw(f"compound/{numeric_id}/structure/", params=params)
            return resp.text
        except ChEBIError as e:
            self.logger.warning(f"Failed to get compound structure for {chebi_id}: {e}")
            return None

    def get_structure(self, structure_id: int, *, width: int = 300,
                      height: int = 300) -> Optional[str]:
        """
        Get raw SVG contents of a structure by its primary key.

        Endpoint: ``GET /structure/{id}/``

        Args:
            structure_id: Primary key of the structure.
            width: Width of the SVG (default 300).
            height: Height of the SVG (default 300).

        Returns:
            Raw SVG string, or *None* on error.

        Examples:
            >>> chebi = ChEBI()
            >>> structure_id = chebi.get_compound(15377)["default_structure"]["id"]  # doctest: +SKIP
            >>> structure_id, chebi.get_structure(structure_id)[:5]  # doctest: +SKIP
            (2018, '<?xml')
        """
        params: Dict[str, Any] = {"width": width, "height": height}
        try:
            resp = self._get_raw(f"structure/{structure_id}/", params=params)
            return resp.text
        except ChEBIError as e:
            self.logger.warning(f"Failed to get structure {structure_id}: {e}")
            return None

    def get_molfile(self, compound_id: int) -> Optional[str]:
        """
        Download the Mol file for a compound's default structure.

        Endpoint: ``GET /molfile/{id}/``

        Args:
            compound_id: Primary key of the compound.

        Returns:
            Mol file contents as string, or *None* on error.

        Examples:
            >>> molfile = ChEBI().get_molfile(15377)              # doctest: +SKIP
            >>> molfile.splitlines()[3][:6]                       # doctest: +SKIP
            '  3  2'
        """
        try:
            resp = self._get_raw(f"molfile/{compound_id}/")
            return resp.text
        except ChEBIError as e:
            self.logger.warning(f"Failed to get molfile for {compound_id}: {e}")
            return None

    # ------------------------------------------------------------------
    # Structure search
    # ------------------------------------------------------------------

    def structure_search(
        self,
        smiles: str,
        search_type: str = "connectivity",
        *,
        similarity: Optional[float] = None,
        three_star_only: bool = True,
        page: int = 1,
        size: int = 15,
        download: bool = False,
    ) -> Optional[Any]:
        """
        Search compounds by chemical structure.

        Endpoint: ``GET /structure_search/``

        Args:
            smiles: Molecule structure in SMILES representation.
            search_type: One of ``connectivity``, ``similarity``, ``substructure``.
            similarity: Similarity threshold for similarity search (0.4–1.0).
            three_star_only: Only include 3-star compounds (default *True*).
            page: Page number (default 1).
            size: Page size (default 15).
            download: Return in download format.

        Returns:
            Search results, or *None* on error.

        Examples:
            >>> chebi = ChEBI()
            >>> results = chebi.structure_search("CCO", "connectivity")  # doctest: +SKIP
            >>> results["total"]                                         # doctest: +SKIP
            2
        """
        params: Dict[str, Any] = {
            "smiles": smiles,
            "search_type": search_type,
            "three_star_only": str(three_star_only).lower(),
            "page": page,
            "size": size,
            "download": str(download).lower(),
        }
        if similarity is not None:
            params["similarity"] = similarity

        try:
            return self._get("structure_search/", params=params)
        except ChEBIError as e:
            self.logger.warning(f"Structure search failed: {e}")
            return None

    # ------------------------------------------------------------------
    # Structure calculations
    # ------------------------------------------------------------------

    def calculate_avg_mass(self, structure: str) -> Optional[str]:
        """
        Calculate the average mass from a structure (a molfile).

        Endpoint: ``POST /structure-calculations/avg-mass/``

        Args:
            structure: A molfile. ChEBI does not accept SMILES here.

        Returns:
            The average mass as text, e.g. ``'46.069'``, or *None* on error (a structure ChEBI cannot parse is
            HTTP 400, logged).

        Examples:
            >>> chebi = ChEBI()
            >>> from rdkit import Chem
            >>> molfile = Chem.MolToMolBlock(Chem.MolFromSmiles("CCO"))
            >>> chebi.calculate_avg_mass(molfile)  # doctest: +SKIP
            '46.069'
        """
        try:
            return self._post_text("structure-calculations/avg-mass/", structure)
        except ChEBIError as e:
            self.logger.warning(f"avg-mass calculation failed: {e}")
            return None

    def calculate_avg_mass_from_formula(self, formula: str) -> Optional[str]:
        """
        Calculate the average mass from a molecular formula.

        Endpoint: ``POST /structure-calculations/avg-mass/from-formula/``

        Args:
            formula: A molecular formula, e.g. ``'C2H6O'``.

        Returns:
            The average mass as text, or *None* on error (a structure ChEBI cannot parse is
            HTTP 400, logged).

        Examples:
            >>> chebi = ChEBI()
            >>> chebi.calculate_avg_mass_from_formula('C2H6O')  # doctest: +SKIP
            '46.069'
        """
        try:
            return self._post_text("structure-calculations/avg-mass/from-formula/", formula)
        except ChEBIError as e:
            self.logger.warning(f"avg-mass-from-formula calculation failed: {e}")
            return None

    def calculate_mol_formula(self, structure: str) -> Optional[str]:
        """
        Calculate the molecular formula from a structure (a molfile).

        Endpoint: ``POST /structure-calculations/mol-formula/``

        Args:
            structure: A molfile. ChEBI does not accept SMILES here.

        Returns:
            The formula as text, or *None* on error (a structure ChEBI cannot parse is
            HTTP 400, logged).

        Examples:
            >>> chebi = ChEBI()
            >>> from rdkit import Chem
            >>> molfile = Chem.MolToMolBlock(Chem.MolFromSmiles("CCO"))
            >>> chebi.calculate_mol_formula(molfile)  # doctest: +SKIP
            'C2H6O'
        """
        try:
            return self._post_text("structure-calculations/mol-formula/", structure)
        except ChEBIError as e:
            self.logger.warning(f"mol-formula calculation failed: {e}")
            return None

    def calculate_monoisotopic_mass(self, structure: str) -> Optional[str]:
        """
        Calculate the monoisotopic mass from a structure (a molfile).

        Endpoint: ``POST /structure-calculations/monoisotopic-mass/``

        Args:
            structure: A molfile. ChEBI does not accept SMILES here.

        Returns:
            The monoisotopic mass as text, or *None* on error (a structure ChEBI cannot parse is
            HTTP 400, logged).

        Examples:
            >>> chebi = ChEBI()
            >>> from rdkit import Chem
            >>> molfile = Chem.MolToMolBlock(Chem.MolFromSmiles("CCO"))
            >>> chebi.calculate_monoisotopic_mass(molfile)  # doctest: +SKIP
            '46.041864812'
        """
        try:
            return self._post_text("structure-calculations/monoisotopic-mass/", structure)
        except ChEBIError as e:
            self.logger.warning(f"monoisotopic-mass calculation failed: {e}")
            return None

    def calculate_monoisotopic_mass_from_formula(self, formula: str) -> Optional[str]:
        """
        Calculate the monoisotopic mass from a molecular formula.

        Endpoint: ``POST /structure-calculations/monoisotopic-mass/from-formula/``

        Args:
            formula: A molecular formula, e.g. ``'C2H6O'``.

        Returns:
            The monoisotopic mass as text, or *None* on error (a structure ChEBI cannot parse is
            HTTP 400, logged).

        Examples:
            >>> chebi = ChEBI()
            >>> chebi.calculate_monoisotopic_mass_from_formula('C2H6O')  # doctest: +SKIP
            '46.041864812'
        """
        try:
            return self._post_text(
                "structure-calculations/monoisotopic-mass/from-formula/", formula,
            )
        except ChEBIError as e:
            self.logger.warning(f"monoisotopic-mass-from-formula calculation failed: {e}")
            return None

    def calculate_net_charge(self, structure: str) -> Optional[str]:
        """
        Calculate the net charge from a structure (a molfile).

        Endpoint: ``POST /structure-calculations/net-charge/``

        Args:
            structure: A molfile. ChEBI does not accept SMILES here.

        Returns:
            The net charge as text, e.g. ``'-1'``, or *None* on error (a structure ChEBI cannot parse is
            HTTP 400, logged).

        Examples:
            >>> chebi = ChEBI()
            >>> from rdkit import Chem
            >>> molfile = Chem.MolToMolBlock(Chem.MolFromSmiles("CCO"))
            >>> chebi.calculate_net_charge(molfile)  # doctest: +SKIP
            '0'
        """
        try:
            return self._post_text("structure-calculations/net-charge/", structure)
        except ChEBIError as e:
            self.logger.warning(f"net-charge calculation failed: {e}")
            return None

    def depict_structure(self, structure: str, *, width: int = 300,
                         height: int = 300,
                         transparent_bg: bool = False) -> Optional[bytes]:
        """
        Generate a PNG depiction of a structure using Indigo.

        Endpoint: ``POST /structure-calculations/depict-indigo/``

        Args:
            structure: Structure in SMILES or molfile format.
            width: Image width in pixels (default 300).
            height: Image height in pixels (default 300).
            transparent_bg: Use transparent background (default *False*).

        Returns:
            PNG image data as bytes, or *None* on error.

        Examples:
            >>> png = ChEBI().depict_structure("CCO")             # doctest: +SKIP
            >>> png[:4]                                           # doctest: +SKIP
            b'\x89PNG'
        """
        params: Dict[str, Any] = {
            "width": width,
            "height": height,
            "transbg": str(transparent_bg).lower(),
        }
        try:
            response = self._post_raw(
                "structure-calculations/depict-indigo/", structure, params=params,
            )
            return response.content
        except ChEBIError as e:
            self.logger.warning(f"depict-indigo failed: {e}")
            return None

    # ------------------------------------------------------------------
    # Batch helpers
    # ------------------------------------------------------------------

    def batch_get_compounds(self, chebi_ids: List[Union[int, str]],
                            pause_time: float = 0.1) -> Dict[str, Dict[str, Any]]:
        """
        Retrieve compound info for multiple ChEBI IDs one-by-one.

        For small batches prefer :meth:`get_compounds` which uses the bulk
        endpoint.  This method calls :meth:`get_compound` in a loop with an
        optional pause between requests.

        Args:
            chebi_ids: List of ChEBI IDs.
            pause_time: Seconds to pause between requests (default 0.1).

        Returns:
            dict mapping ``CHEBI:<id>`` to compound data. IDs that could not
            be fetched are left out.

        Examples:
            >>> found = ChEBI().batch_get_compounds([15377, 15365])  # doctest: +SKIP
            >>> {key: record["name"] for key, record in found.items()}  # doctest: +SKIP
            {'CHEBI:15377': 'water', 'CHEBI:15365': 'acetylsalicylic acid'}
        """
        results: Dict[str, Dict[str, Any]] = {}
        for chebi_id in chebi_ids:
            key = self._format_chebi_id(chebi_id)
            entity = self.get_compound(chebi_id)
            if entity:
                results[key] = entity
            if pause_time > 0:
                time.sleep(pause_time)
        return results

    # ------------------------------------------------------------------
    # Dunder
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        """String representation of ChEBI client."""
        return f"ChEBI(base_url='{self.base_url}', timeout={self.timeout})"


# ------------------------------------------------------------------
# Convenience functions
# ------------------------------------------------------------------

def get_chebi_entity(chebi_id: Union[int, str]) -> Optional[Dict[str, Any]]:
    """
    Convenience function to get ChEBI compound information.

    Args:
        chebi_id: ChEBI ID (with or without ``CHEBI:`` prefix).

    Returns:
        dict with compound data, or *None* if not found.

    Examples:
        >>> from provesid import get_chebi_entity
        >>> water = get_chebi_entity(15377)                   # doctest: +SKIP
        >>> water["default_structure"]["smiles"]              # doctest: +SKIP
        '[H]O[H]'
    """
    chebi = ChEBI()
    return chebi.get_compound(chebi_id)


def search_chebi(search_text: str, max_results: int = 10) -> List[Dict[str, Any]]:
    """
    Convenience function to search ChEBI by name.

    Args:
        search_text: Text to search for.
        max_results: Maximum number of results to return.

    Returns:
        list of matching entity dicts.

    Examples:
        >>> from provesid import search_chebi
        >>> results = search_chebi("paracetamol", max_results=1)  # doctest: +SKIP
        >>> results[0]["_source"]["name"]                         # doctest: +SKIP
        'paracetamol'
    """
    chebi = ChEBI()
    return chebi.search_by_name(search_text, size=max_results)
