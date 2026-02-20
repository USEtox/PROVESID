"""
ChEBI (Chemical Entities of Biological Interest) API interface.

This module provides a Python interface to the ChEBI 2.0 REST API for retrieving
chemical compound information from the ChEBI database, as well as a local
SDF file parser for offline access to ChEBI data.

API documentation: https://www.ebi.ac.uk/chebi/backend/api/docs/

Author: USEtox team
Date: August 2025
"""

import requests
import logging
import time
import os
import pickle
import gzip
import shutil
from typing import Dict, List, Optional, Union, Any
from rdkit import Chem
import pandas as pd
from tqdm import tqdm
from .utils import data_path


class ChEBIError(Exception):
    """Custom exception for ChEBI API errors."""
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

    Example:
        >>> chebi = ChEBI()
        >>> compound = chebi.get_compound(15377)  # CHEBI:15377 (water)
        >>> print(compound['name'])
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

    def _get(self, endpoint: str, params: Optional[Dict] = None) -> Any:
        """
        Perform a GET request and return the parsed JSON body.

        Args:
            endpoint (str): Path relative to *base_url* (e.g. ``compound/15377/``).
            params (dict, optional): Query-string parameters.

        Returns:
            Parsed JSON response (dict / list / str).

        Raises:
            ChEBIError: On network / HTTP / JSON errors.
        """
        url = f"{self.base_url}/{endpoint}"
        try:
            response = self.session.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            if response.headers.get("Content-Type", "").startswith("application/json"):
                return response.json()
            return response.text
        except requests.exceptions.Timeout:
            raise ChEBIError(f"Request timeout after {self.timeout} seconds")
        except requests.exceptions.RequestException as e:
            raise ChEBIError(f"Request failed: {str(e)}")

    def _get_raw(self, endpoint: str, params: Optional[Dict] = None) -> requests.Response:
        """
        Perform a GET request and return the raw :class:`requests.Response`.

        Useful for endpoints that return non-JSON content (SVG, molfile, images).
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
            ChEBIError: On network / HTTP / JSON errors.
        """
        url = f"{self.base_url}/{endpoint}"
        try:
            response = self.session.post(
                url, json=json_body, params=params, timeout=self.timeout,
            )
            response.raise_for_status()
            if response.headers.get("Content-Type", "").startswith("application/json"):
                return response.json()
            return response.text
        except requests.exceptions.Timeout:
            raise ChEBIError(f"Request timeout after {self.timeout} seconds")
        except requests.exceptions.RequestException as e:
            raise ChEBIError(f"Request failed: {str(e)}")

    def _post_text(self, endpoint: str, text_body: str,
                   params: Optional[Dict] = None) -> str:
        """
        Perform a POST request with a ``text/plain`` body and return the
        response text.  Used by the structure-calculation endpoints.
        """
        url = f"{self.base_url}/{endpoint}"
        try:
            response = self.session.post(
                url, data=text_body, params=params, timeout=self.timeout,
                headers={**self.session.headers, "Content-Type": "text/plain;charset=UTF-8"},
            )
            response.raise_for_status()
            return response.text
        except requests.exceptions.Timeout:
            raise ChEBIError(f"Request timeout after {self.timeout} seconds")
        except requests.exceptions.RequestException as e:
            raise ChEBIError(f"Request failed: {str(e)}")

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

        Example:
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

    def get_compounds(self, chebi_ids: List[Union[int, str]]) -> Optional[Any]:
        """
        Retrieve information about one or more compounds in a single call.

        Endpoint: ``POST /compounds/``

        Args:
            chebi_ids: List of ChEBI IDs.

        Returns:
            API response (typically a list of compound dicts), or *None* on error.

        Example:
            >>> chebi = ChEBI()
            >>> results = chebi.get_compounds(["CHEBI:15377", "CHEBI:16236"])
        """
        ids_formatted = [self._format_chebi_id(cid) for cid in chebi_ids]
        try:
            return self._post_json("compounds/", json_body={"chebi_ids": ids_formatted})
        except ChEBIError as e:
            self.logger.warning(f"Failed to get compounds: {e}")
            return None

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
            Search results dict, or *None* on error.

        Example:
            >>> chebi = ChEBI()
            >>> results = chebi.search("paracetamol")
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
            list of matching entity dicts (may be empty).
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

        Example:
            >>> chebi = ChEBI()
            >>> results = chebi.advanced_search({
            ...     "formula_specification": {
            ...         "and_specification": [{"term": "C6H12O7"}]
            ...     }
            ... }, three_star_only=False)
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
            Sources information, or *None* on error.
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
            Ontology parent data, or *None* on error.
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
            Ontology children data, or *None* on error.
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

        Example:
            >>> chebi = ChEBI()
            >>> # Get all compounds that are alcohols
            >>> results = chebi.get_all_ontology_children_in_path(
            ...     relation="is_a", entity="CHEBI:30879"
            ... )
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

        Example:
            >>> chebi = ChEBI()
            >>> results = chebi.structure_search("c1ccccc1", "substructure")
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
        Calculate the average mass from a structure (SMILES or molfile).

        Endpoint: ``POST /structure-calculations/avg-mass/``
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
        """
        try:
            return self._post_text("structure-calculations/avg-mass/from-formula/", formula)
        except ChEBIError as e:
            self.logger.warning(f"avg-mass-from-formula calculation failed: {e}")
            return None

    def calculate_mol_formula(self, structure: str) -> Optional[str]:
        """
        Calculate the molecular formula from a structure (SMILES or molfile).

        Endpoint: ``POST /structure-calculations/mol-formula/``
        """
        try:
            return self._post_text("structure-calculations/mol-formula/", structure)
        except ChEBIError as e:
            self.logger.warning(f"mol-formula calculation failed: {e}")
            return None

    def calculate_monoisotopic_mass(self, structure: str) -> Optional[str]:
        """
        Calculate the monoisotopic mass from a structure (SMILES or molfile).

        Endpoint: ``POST /structure-calculations/monoisotopic-mass/``
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
        Calculate the net charge from a structure (SMILES or molfile).

        Endpoint: ``POST /structure-calculations/net-charge/``
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
        """
        url = f"{self.base_url}/structure-calculations/depict-indigo/"
        params: Dict[str, Any] = {
            "width": width,
            "height": height,
            "transbg": str(transparent_bg).lower(),
        }
        try:
            response = self.session.post(
                url, data=structure, params=params, timeout=self.timeout,
                headers={**self.session.headers, "Content-Type": "text/plain;charset=UTF-8"},
            )
            response.raise_for_status()
            return response.content
        except requests.exceptions.Timeout:
            raise ChEBIError(f"Request timeout after {self.timeout} seconds")
        except requests.exceptions.RequestException as e:
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
            dict mapping ``CHEBI:<id>`` to compound data.
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

    Example:
        >>> from provesid import get_chebi_entity
        >>> water = get_chebi_entity(15377)
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

    Example:
        >>> from provesid import search_chebi
        >>> results = search_chebi("aspirin")
    """
    chebi = ChEBI()
    return chebi.search_by_name(search_text, size=max_results)


class ChebiSDF:
    """
    Parser for ChEBI SDF (Structure-Data File) for offline access to ChEBI data.
    
    This class provides efficient querying of the ChEBI SDF file containing
    ~190,000 compounds with structure data, chemical properties, synonyms,
    and cross-references to 80+ external databases.
    
    The class builds an index on first use for fast lookups. The index is
    cached to disk for faster subsequent initializations.
    
    If the SDF file is not found, it can be automatically downloaded from the
    ChEBI FTP server.
    
    Attributes:
        sdf_path (str): Path to ChEBI SDF file
        index (dict): In-memory index for fast lookups
        
    Example:
        >>> chebi_sdf = ChebiSDF()
        >>> compound = chebi_sdf.get_compound_by_id("CHEBI:15377")
        >>> print(compound['ChEBI NAME'])
        water
    """
    
    # Default download URL for ChEBI SDF file
    DEFAULT_SDF_URL = "https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/chebi.sdf.gz"
    
    def __init__(self, sdf_path: Optional[str] = None, rebuild_index: bool = False, 
                 auto_download: bool = True, sdf_url: Optional[str] = None):
        """
        Initialize ChebiSDF parser.
        
        Args:
            sdf_path (str, optional): Path to ChEBI SDF file. If None, uses default location.
            rebuild_index (bool): If True, rebuild index even if cache exists (default: False)
            auto_download (bool): If True, automatically download SDF file if not found (default: True)
            sdf_url (str, optional): Custom URL to download the SDF from. If None, uses default.
        """
        
        if sdf_path is None:
            data_dir = data_path()
            sdf_path = os.path.join(data_dir, 'chebi.sdf')
        
        self.sdf_path = sdf_path
        self.index_path = sdf_path + '.index.pkl'
        self.sdf_url = sdf_url or self.DEFAULT_SDF_URL
        self.logger = logging.getLogger(__name__)
        
        # Check if SDF file exists, download if needed
        if not os.path.exists(self.sdf_path):
            if auto_download:
                self.logger.info(f"ChEBI SDF file not found at: {self.sdf_path}")
                self.logger.info("Downloading ChEBI SDF file automatically...")
                self.download_sdf(url=self.sdf_url, force=False)
            else:
                raise FileNotFoundError(
                    f"ChEBI SDF file not found at: {self.sdf_path}\n"
                    f"Please run ChebiSDF.download_sdf() or set auto_download=True\n"
                    f"Or download manually from: https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/"
                )
        
        # Load or build index
        if rebuild_index or not os.path.exists(self.index_path):
            self.logger.info("Building index from SDF file...")
            self.index = self._build_index()
            self._save_index()
        else:
            self.logger.info("Loading cached index...")
            self.index = self._load_index()
    
    def download_sdf(self, url: Optional[str] = None, force: bool = False):
        """
        Download the ChEBI SDF file from the ChEBI FTP server.
        
        The file is downloaded as a gzip archive (~250 MB) and automatically
        extracted to the data directory (~868 MB uncompressed).
        
        Args:
            url (str, optional): URL to download from. If None, uses default ChEBI FTP URL.
            force (bool): If True, download even if file already exists (default: False)
            
        Returns:
            str: Path to the downloaded and extracted SDF file
            
        Raises:
            FileExistsError: If the file already exists and force=False
            requests.exceptions.RequestException: If the download fails
            
        Example:
            >>> chebi_sdf = ChebiSDF(auto_download=False)
            >>> chebi_sdf.download_sdf()  # Manually trigger download
        """
        download_url = url or self.sdf_url
        
        # Check if file already exists
        if os.path.exists(self.sdf_path) and not force:
            raise FileExistsError(
                f"ChEBI SDF file already exists at: {self.sdf_path}\n"
                f"Use force=True to overwrite"
            )
        
        # Create data directory if it doesn't exist
        data_dir = os.path.dirname(self.sdf_path)
        os.makedirs(data_dir, exist_ok=True)
        
        self.logger.info(f"Downloading ChEBI SDF from: {download_url}")
        self.logger.info(f"Destination: {self.sdf_path}")
        
        try:
            # Stream the download with progress bar
            response = requests.get(download_url, stream=True, timeout=60)
            response.raise_for_status()
            
            # Get total file size
            total_size = int(response.headers.get('content-length', 0))
            
            # Download gzipped file to temporary location
            gz_path = self.sdf_path + '.gz.tmp'
            with open(gz_path, 'wb') as f:
                if total_size > 0:
                    with tqdm(total=total_size, unit='B', unit_scale=True, 
                             desc="Downloading ChEBI SDF") as pbar:
                        for chunk in response.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                                pbar.update(len(chunk))
                else:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
            
            self.logger.info("Download complete. Extracting gzip archive...")
            
            # Extract gzipped file
            temp_path = self.sdf_path + '.tmp'
            with gzip.open(gz_path, 'rb') as f_in:
                with open(temp_path, 'wb') as f_out:
                    # Copy with progress bar
                    file_size = os.path.getsize(gz_path)
                    with tqdm(total=file_size, unit='B', unit_scale=True, 
                             desc="Extracting ChEBI SDF") as pbar:
                        while True:
                            chunk = f_in.read(8192)
                            if not chunk:
                                break
                            f_out.write(chunk)
                            pbar.update(len(chunk))
            
            # Move temp file to final location
            if os.path.exists(self.sdf_path):
                os.remove(self.sdf_path)
            os.rename(temp_path, self.sdf_path)
            
            # Clean up temporary gzip file
            if os.path.exists(gz_path):
                os.remove(gz_path)
            
            self.logger.info(f"✓ ChEBI SDF file downloaded and extracted successfully")
            self.logger.info(f"✓ File location: {self.sdf_path}")
            
            # Verify the file is valid by checking first few lines
            try:
                with open(self.sdf_path, 'r', encoding='utf-8', errors='ignore') as f:
                    first_lines = [f.readline() for _ in range(5)]
                    if not any(first_lines):
                        raise RuntimeError("Downloaded file appears to be empty")
                self.logger.info("✓ File verified successfully")
            except Exception as e:
                os.remove(self.sdf_path)
                raise RuntimeError(f"Downloaded file is corrupted: {e}")
            
            return self.sdf_path
            
        except requests.exceptions.RequestException as e:
            # Clean up temp files if they exist
            if os.path.exists(gz_path):
                os.remove(gz_path)
            if os.path.exists(temp_path):
                os.remove(temp_path)
            self.logger.error(f"Failed to download ChEBI SDF: {e}")
            raise
    
    def _build_index(self) -> Dict:
        """
        Build index from SDF file for fast lookups.
        
        Returns:
            dict: Index containing mappings for various query types
        """
        index = {
            'id_to_offset': {},           # ChEBI ID -> file offset
            'name_to_ids': {},            # lowercase name -> list of ChEBI IDs
            'inchikey_to_id': {},         # InChIKey -> ChEBI ID
            'inchi_to_id': {},            # InChI -> ChEBI ID
            'formula_to_ids': {},         # formula -> list of ChEBI IDs
            'cas_to_ids': {},             # CAS number -> list of ChEBI IDs
            'synonym_to_ids': {},         # lowercase synonym -> list of ChEBI IDs
        }
        
        # Parse SDF file and build index
        with open(self.sdf_path, 'r', encoding='utf-8', errors='ignore') as f:
            file_offset = 0
            current_mol_offset = 0
            in_mol = False
            current_data = {}
            
            # Use tqdm for progress bar
            file_size = os.path.getsize(self.sdf_path)
            pbar = tqdm(total=file_size, unit='B', unit_scale=True, desc="Indexing ChEBI SDF")
            
            for line in f:
                if not in_mol and line.strip():
                    # Start of a new molecule
                    in_mol = True
                    current_mol_offset = file_offset
                    current_data = {}
                
                # Check for property tags
                if line.startswith('> <'):
                    field_name = line.strip()[3:-1]  # Extract field name
                    value_line = next(f, '')
                    file_offset += len(line.encode('utf-8'))
                    file_offset += len(value_line.encode('utf-8'))
                    current_data[field_name] = value_line.strip()
                    pbar.update(len(line.encode('utf-8')) + len(value_line.encode('utf-8')))
                    continue
                
                # Check for end of molecule
                if line.startswith('$$$$'):
                    in_mol = False
                    
                    # Index this molecule
                    if 'ChEBI ID' in current_data:
                        chebi_id = current_data['ChEBI ID']
                        index['id_to_offset'][chebi_id] = current_mol_offset
                        
                        # Index by name
                        if 'ChEBI NAME' in current_data:
                            name_lower = current_data['ChEBI NAME'].lower()
                            if name_lower not in index['name_to_ids']:
                                index['name_to_ids'][name_lower] = []
                            index['name_to_ids'][name_lower].append(chebi_id)
                        
                        # Index by InChIKey
                        if 'INCHIKEY' in current_data:
                            index['inchikey_to_id'][current_data['INCHIKEY']] = chebi_id
                        
                        # Index by InChI
                        if 'INCHI' in current_data:
                            index['inchi_to_id'][current_data['INCHI']] = chebi_id
                        
                        # Index by formula
                        if 'FORMULA' in current_data:
                            formula = current_data['FORMULA']
                            if formula not in index['formula_to_ids']:
                                index['formula_to_ids'][formula] = []
                            index['formula_to_ids'][formula].append(chebi_id)
                        
                        # Index by CAS
                        if 'CAS Registry Numbers' in current_data:
                            cas_numbers = current_data['CAS Registry Numbers'].split(';')
                            for cas in cas_numbers:
                                cas = cas.strip()
                                if cas:
                                    if cas not in index['cas_to_ids']:
                                        index['cas_to_ids'][cas] = []
                                    index['cas_to_ids'][cas].append(chebi_id)
                        
                        # Index by synonyms
                        if 'SYNONYM' in current_data:
                            synonyms = current_data['SYNONYM'].split(';')
                            for syn in synonyms:
                                syn_lower = syn.strip().lower()
                                if syn_lower:
                                    if syn_lower not in index['synonym_to_ids']:
                                        index['synonym_to_ids'][syn_lower] = []
                                    index['synonym_to_ids'][syn_lower].append(chebi_id)
                
                file_offset += len(line.encode('utf-8'))
                pbar.update(len(line.encode('utf-8')))
            
            pbar.close()
        
        self.logger.info(f"Index built: {len(index['id_to_offset'])} compounds indexed")
        return index
    
    def _save_index(self):
        """Save index to disk for faster subsequent loads."""
        try:
            with open(self.index_path, 'wb') as f:
                pickle.dump(self.index, f, protocol=pickle.HIGHEST_PROTOCOL)
            self.logger.info(f"Index saved to {self.index_path}")
        except Exception as e:
            self.logger.warning(f"Failed to save index: {e}")
    
    def _load_index(self) -> Dict:
        """Load index from disk."""
        try:
            with open(self.index_path, 'rb') as f:
                index = pickle.load(f)
            self.logger.info(f"Index loaded: {len(index['id_to_offset'])} compounds")
            return index
        except Exception as e:
            self.logger.warning(f"Failed to load index: {e}. Rebuilding...")
            return self._build_index()
    
    def _read_mol_at_offset(self, offset: int) -> Dict[str, str]:
        """
        Read a molecule entry from the SDF file at a specific offset.
        
        Args:
            offset (int): File offset where molecule starts
            
        Returns:
            dict: Molecule data including molfile and properties
        """
        data = {'molfile': ''}
        
        with open(self.sdf_path, 'r', encoding='utf-8', errors='ignore') as f:
            f.seek(offset)
            in_molfile = True
            
            for line in f:
                if in_molfile:
                    data['molfile'] += line
                    if line.startswith('M  END'):
                        in_molfile = False
                    continue
                
                # Parse property tags
                if line.startswith('> <'):
                    field_name = line.strip()[3:-1]
                    value_line = next(f, '').strip()
                    data[field_name] = value_line
                    continue
                
                # End of molecule
                if line.startswith('$$$$'):
                    break
        
        return data
    
    def get_compound_by_id(self, chebi_id: str) -> Optional[Dict[str, str]]:
        """
        Get compound data by ChEBI ID.
        
        Args:
            chebi_id (str): ChEBI ID (e.g., "CHEBI:15377" or "15377")
            
        Returns:
            dict: Compound data, or None if not found
            
        Example:
            >>> chebi_sdf = ChebiSDF()
            >>> water = chebi_sdf.get_compound_by_id("CHEBI:15377")
            >>> print(water['ChEBI NAME'])
            water
        """
        # Normalize ID format
        if not chebi_id.startswith('CHEBI:'):
            chebi_id = f'CHEBI:{chebi_id}'
        
        offset = self.index['id_to_offset'].get(chebi_id)
        if offset is None:
            return None
        
        return self._read_mol_at_offset(offset)
    
    def search_by_name(self, name: str, exact: bool = True) -> List[Dict[str, str]]:
        """
        Search compounds by name.
        
        Args:
            name (str): Compound name to search for
            exact (bool): If True, exact match; if False, partial match (default: True)
            
        Returns:
            list: List of matching compound data
            
        Example:
            >>> chebi_sdf = ChebiSDF()
            >>> results = chebi_sdf.search_by_name("water")
            >>> print(len(results))
            1
        """
        name_lower = name.lower()
        results = []
        
        if exact:
            chebi_ids = self.index['name_to_ids'].get(name_lower, [])
        else:
            # Partial match
            chebi_ids = []
            for indexed_name, ids in self.index['name_to_ids'].items():
                if name_lower in indexed_name:
                    chebi_ids.extend(ids)
        
        for chebi_id in chebi_ids:
            compound = self.get_compound_by_id(chebi_id)
            if compound:
                results.append(compound)
        
        return results
    
    def search_by_synonym(self, synonym: str, exact: bool = True) -> List[Dict[str, str]]:
        """
        Search compounds by synonym.
        
        Args:
            synonym (str): Synonym to search for
            exact (bool): If True, exact match; if False, partial match (default: True)
            
        Returns:
            list: List of matching compound data
        """
        synonym_lower = synonym.lower()
        results = []
        
        if exact:
            chebi_ids = self.index['synonym_to_ids'].get(synonym_lower, [])
        else:
            # Partial match
            chebi_ids = []
            for indexed_syn, ids in self.index['synonym_to_ids'].items():
                if synonym_lower in indexed_syn:
                    chebi_ids.extend(ids)
            # Remove duplicates
            chebi_ids = list(set(chebi_ids))
        
        for chebi_id in chebi_ids:
            compound = self.get_compound_by_id(chebi_id)
            if compound:
                results.append(compound)
        
        return results
    
    def search_by_inchikey(self, inchikey: str) -> Optional[Dict[str, str]]:
        """
        Search compound by InChIKey.
        
        Args:
            inchikey (str): InChIKey to search for
            
        Returns:
            dict: Compound data, or None if not found
        """
        chebi_id = self.index['inchikey_to_id'].get(inchikey)
        if chebi_id:
            return self.get_compound_by_id(chebi_id)
        return None
    
    def search_by_inchi(self, inchi: str) -> Optional[Dict[str, str]]:
        """
        Search compound by InChI.
        
        Args:
            inchi (str): InChI string to search for
            
        Returns:
            dict: Compound data, or None if not found
        """
        chebi_id = self.index['inchi_to_id'].get(inchi)
        if chebi_id:
            return self.get_compound_by_id(chebi_id)
        return None
    
    def search_by_cas(self, cas: str) -> List[Dict[str, str]]:
        """
        Search compounds by CAS Registry Number.
        
        Args:
            cas (str): CAS Registry Number
            
        Returns:
            list: List of matching compound data
        """
        chebi_ids = self.index['cas_to_ids'].get(cas, [])
        results = []
        
        for chebi_id in chebi_ids:
            compound = self.get_compound_by_id(chebi_id)
            if compound:
                results.append(compound)
        
        return results
    
    def search_by_formula(self, formula: str) -> List[Dict[str, str]]:
        """
        Search compounds by molecular formula.
        
        Args:
            formula (str): Molecular formula (e.g., "H2O")
            
        Returns:
            list: List of matching compound data
        """
        chebi_ids = self.index['formula_to_ids'].get(formula, [])
        results = []
        
        for chebi_id in chebi_ids:
            compound = self.get_compound_by_id(chebi_id)
            if compound:
                results.append(compound)
        
        return results
    
    def filter_by_star_rating(self, min_stars: int = 3) -> List[str]:
        """
        Get ChEBI IDs of compounds with minimum star rating.
        
        Args:
            min_stars (int): Minimum star rating (1-3, default: 3)
            
        Returns:
            list: List of ChEBI IDs matching criteria
        """
        matching_ids = []
        
        for chebi_id in tqdm(self.index['id_to_offset'].keys(), desc="Filtering by star rating"):
            compound = self.get_compound_by_id(chebi_id)
            if compound and 'STAR' in compound:
                try:
                    star = int(compound['STAR'])
                    if star >= min_stars:
                        matching_ids.append(chebi_id)
                except ValueError:
                    continue
        
        return matching_ids
    
    def get_compounds_by_ids(self, chebi_ids: List[str]) -> List[Dict[str, str]]:
        """
        Get multiple compounds by ChEBI IDs.
        
        Args:
            chebi_ids (list): List of ChEBI IDs
            
        Returns:
            list: List of compound data dictionaries
        """
        results = []
        for chebi_id in chebi_ids:
            compound = self.get_compound_by_id(chebi_id)
            if compound:
                results.append(compound)
        return results
    
    def export_to_dataframe(self, chebi_ids: Optional[List[str]] = None, 
                           fields: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Export compounds to pandas DataFrame.
        
        Args:
            chebi_ids (list, optional): List of ChEBI IDs. If None, exports all compounds.
            fields (list, optional): List of field names to include. If None, includes common fields.
            
        Returns:
            pd.DataFrame: DataFrame with compound data
            
        Example:
            >>> chebi_sdf = ChebiSDF()
            >>> df = chebi_sdf.export_to_dataframe(["CHEBI:15377", "CHEBI:16236"])
            >>> print(df[['ChEBI ID', 'ChEBI NAME', 'FORMULA']])
        """
        if fields is None:
            fields = ['ChEBI ID', 'ChEBI NAME', 'STAR', 'FORMULA', 'MASS', 
                     'SMILES', 'INCHI', 'INCHIKEY', 'CAS Registry Numbers']
        
        if chebi_ids is None:
            chebi_ids = list(self.index['id_to_offset'].keys())
        
        data = []
        for chebi_id in tqdm(chebi_ids, desc="Exporting to DataFrame"):
            compound = self.get_compound_by_id(chebi_id)
            if compound:
                row = {field: compound.get(field, None) for field in fields}
                data.append(row)
        
        return pd.DataFrame(data)
    
    def get_database_stats(self) -> Dict[str, int]:
        """
        Get statistics about the ChEBI SDF database.
        
        Returns:
            dict: Statistics including counts of various indexed fields
        """
        return {
            'total_compounds': len(self.index['id_to_offset']),
            'compounds_with_inchikey': len(self.index['inchikey_to_id']),
            'compounds_with_inchi': len(self.index['inchi_to_id']),
            'compounds_with_cas': len(self.index['cas_to_ids']),
            'unique_formulas': len(self.index['formula_to_ids']),
            'indexed_names': len(self.index['name_to_ids']),
            'indexed_synonyms': len(self.index['synonym_to_ids']),
        }
