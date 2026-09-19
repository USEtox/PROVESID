import requests
import logging
import time
from urllib.parse import quote
from typing import Dict, Iterator, List, Optional, Any, Tuple, Union
from dataclasses import dataclass
import pandas as pd
import re
from .cache import cached
from .pubchemview_parse import ParsedValue, parse_value


@dataclass
class PropertyData:
    """
    One value PubChem holds for one property, with its provenance.

    Attributes:
        value: The value exactly as PUG-View reported it, e.g. ``"138-140 °C"``.
            Always populated; nothing the parser cannot read is ever lost.
        unit: The unit, normalised — a convenience copy of ``parsed.unit``.
        conditions: The measurement conditions, a copy of ``parsed.conditions``.
        reference: The first reference string attached to the value.
        reference_number: PubChem's reference number, which
            :meth:`PubChemView.get_property_table` resolves to a full citation.
        description: PUG-View's own description of the value, when it gives one.
        name: PUG-View's ``Name`` field for the value, when it gives one.
        heading: The section heading the value was found under, e.g.
            ``"Melting Point"``. This is what tells the parser that a bare
            number is a temperature.
        parsed: The value turned into numbers: see :class:`ParsedValue` for the
            range, the SI conversion, the measurement temperature and any
            comparison operator or qualitative term. None only when the value
            string was empty.
    """
    value: str
    unit: Optional[str] = None
    conditions: Optional[str] = None
    reference: Optional[str] = None
    reference_number: Optional[int] = None
    description: Optional[str] = None
    name: Optional[str] = None
    heading: Optional[str] = None
    parsed: Optional[ParsedValue] = None


class PubChemViewError(Exception):
    """Base exception class for PubChem View API errors"""
    pass


#: Fault codes PubChem returns when it is shedding load rather than reporting
#: absence. PubChem does not always pair them with a 5xx status, so the code in
#: the body — not the HTTP status alone — decides whether a retry makes sense.
TRANSIENT_FAULT_CODES = frozenset({
    "PUGVIEW.ServerBusy",
    "PUGVIEW.ServerError",
    "PUGVIEW.Timeout",
    "PUGREST.ServerBusy",
    "PUGREST.ServerError",
    "PUGREST.Timeout",
})


def fault_code(response: "requests.Response") -> Optional[str]:
    """
    Read PubChem's ``Fault.Code`` out of an error response.

    Both PUG-REST and PUG-View describe every error in the body, for example
    ``{"Fault": {"Code": "PUGVIEW.NotFound", "Message": "No data found"}}``.
    That code is the only reliable way to tell a compound that genuinely lacks
    data from a service that is momentarily refusing work.

    Args:
        response: The HTTP response to inspect.

    Returns:
        The fault code, or None when the body is not a PubChem fault.

    Example:
        >>> fault_code(response)   # doctest: +SKIP
        'PUGVIEW.NotFound'
    """
    try:
        return response.json().get("Fault", {}).get("Code")
    except Exception:
        return None


def _is_empty_lookup(result: Any) -> bool:
    """
    Report whether an extraction returned nothing.

    Absence is never cached. PubChem answers a genuinely empty property with a
    ``NotFound`` fault, but it has also been observed returning one while under
    load, and a wrongly cached "no data" is permanent whereas re-fetching an
    empty result costs one cheap request. Positive results stay cached.

    Args:
        result: Return value of an extraction method --- a list, a dict or a
            DataFrame.

    Returns:
        True when the result carries no extracted data.

    Example:
        >>> _is_empty_lookup([])
        True
        >>> _is_empty_lookup(['3.47'])
        False
    """
    if result is None:
        return True
    try:
        return len(result) == 0
    except TypeError:
        # Not a sized value; treat it as a real result.
        return False


class PubChemViewNotFoundError(PubChemViewError):
    """Exception raised when compound or property is not found"""
    pass


class PubChemView:
    """
    A class that uses PUG View for extracting properties reported for each substance in PubChem but are not
    included in the standard API response for substance, compound, assay, etc.
    The response to these queries is a large JSON object that requires some post-processing to extract the
    relevant information.
    """
    
    def __init__(self, base_url: str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view",
                 timeout: int = 30, max_retries: int = 3, backoff_factor: float = 1.0, use_cache: bool = True):
        """
        Initialize PubChemView API client
        
        Args:
            base_url: Base URL for PubChem PUG View API
            timeout: Request timeout in seconds
            max_retries: Maximum number of retry attempts
            backoff_factor: Backoff factor for retries
            use_cache: Whether to use cache for lookups (default: True). 
                      When False, skips cache lookup but still stores results.
        """
        self.base_url = base_url.rstrip('/')
        self.timeout = timeout
        self.max_retries = max_retries
        self.backoff_factor = backoff_factor
        self.logger = logging.getLogger(__name__)
        self.use_cache = use_cache
        
        # Rate limiting
        self.last_request_time = 0
        self.min_request_interval = 0.2  # 5 requests per second max
        
        # Standard experimental property headings
        self.experimental_properties = {
            "Accelerating Rate Calorimetry (ARC)": "Accelerating+Rate+Calorimetry+(ARC)",
            "Acid Value": "Acid+Value",
            "Autoignition Temperature": "Autoignition+Temperature", 
            "Boiling Point": "Boiling+Point",
            "Caco2 Permeability": "Caco2+Permeability",
            "Collision Cross Section": "Collision+Cross+Section",
            "Color/Form": "Color/Form",
            "Corrosivity": "Corrosivity",
            "Decomposition": "Decomposition",
            "Density": "Density",
            "Dielectric Constant": "Dielectric+Constant",
            "Differential Scanning Calorimetry (DSC)": "Differential+Scanning+Calorimetry+(DSC)",
            "Dispersion": "Dispersion",
            "Dissociation Constants": "Dissociation+Constants",
            "Enthalpy of Sublimation": "Enthalpy+of+Sublimation",
            "Flash Point": "Flash+Point",
            "Heat of Combustion": "Heat+of+Combustion",
            "Heat of Vaporization": "Heat+of+Vaporization",
            "Henry's Law Constant": "Henry's+Law+Constant",
            "Hydrophobicity": "Hydrophobicity",
            "Ionization Efficiency": "Ionization+Efficiency",
            "Ionization Potential": "Ionization+Potential",
            "Isoelectric Point": "Isoelectric+Point",
            "Kovats Retention Index": "Kovats+Retention+Index",
            "LogP": "LogP",
            "LogS": "LogS",
            "Melting Point": "Melting+Point",
            "Odor": "Odor",
            "Odor Threshold": "Odor+Threshold",
            "Optical Rotation": "Optical+Rotation",
            "Other Experimental Properties": "Other+Experimental+Properties",
            "pH": "pH",
            "Physical Description": "Physical+Description",
            "Polymerization": "Polymerization",
            "Refractive Index": "Refractive+Index",
            "Relative Evaporation Rate": "Relative+Evaporation+Rate",
            "Self-Accelerating Decomposition Temperature (SADT)": "Self-Accelerating+Decomposition+Temperature+(SADT)",
            "Solubility": "Solubility",
            "Stability/Shelf Life": "Stability/Shelf+Life",
            "Surface Tension": "Surface+Tension",
            "Taste": "Taste",
            "Vapor Density": "Vapor+Density",
            "Vapor Pressure": "Vapor+Pressure",
            "Viscosity": "Viscosity"
        }
    
    #: Bumped whenever the *shape* of what a cached method returns changes, so
    #: that entries written by an earlier version become unreachable instead of
    #: being deserialised into the wrong structure. Version 2 introduced
    #: ``PropertyData.parsed``: a version-1 entry unpickles into an object with
    #: no such attribute, and every caller that reads it would raise.
    CACHE_SCHEMA_VERSION = 2

    def __cache_key__(self) -> tuple:
        """
        Identify this client for cache-key purposes.

        Only the endpoint distinguishes two clients' results; timeout, retry and
        ``use_cache`` settings change how a call is made, not what it returns.
        The schema version is part of the key so that an upgrade cannot serve a
        stale entry of the previous shape.

        Returns:
            Tuple of the class path, the configured base URL and
            :attr:`CACHE_SCHEMA_VERSION`.
        """
        return ("provesid.pubchemview.PubChemView", self.base_url,
                self.CACHE_SCHEMA_VERSION)

    def clear_cache(self):
        """Clear all cached results for PubChem View"""
        from .cache import clear_pubchemview_cache
        clear_pubchemview_cache()
    
    def get_cache_info(self):
        """Get cache statistics for PubChem View cached methods"""
        from .cache import get_pubchemview_cache_info
        return get_pubchemview_cache_info()
    
    def _rate_limit(self):
        """Implement rate limiting to respect PubChem's usage policy"""
        current_time = time.time()
        time_since_last = current_time - self.last_request_time
        if time_since_last < self.min_request_interval:
            time.sleep(self.min_request_interval - time_since_last)
        self.last_request_time = time.time()
    
    def _make_request(self, url: str) -> Dict[str, Any]:
        """
        Make HTTP request with retries and error handling
        
        Args:
            url: Request URL
            
        Returns:
            JSON response as dictionary

        Raises:
            PubChemViewNotFoundError: When the compound or heading does not
                exist. PUG-View answers 404 for an unknown compound and 400
                (``PUGVIEW.BadRequest``) for an unknown heading; both are
                permanent, so both are reported as absence rather than retried.
            PubChemViewError: For any other API error, including a transient
                ``ServerBusy`` that survived every retry. Only 429 and 5xx
                responses, timeouts and connection errors are retried — a 4xx
                cannot be fixed by asking again.
        """
        self._rate_limit()
        
        for attempt in range(self.max_retries + 1):
            try:
                self.logger.debug(f"Making request to: {url}")
                response = requests.get(url, timeout=self.timeout)
                
                if response.status_code == 200:
                    return response.json()

                code = fault_code(response)

                if (code in TRANSIENT_FAULT_CODES
                        or response.status_code == 429
                        or response.status_code >= 500):
                    # Transient. Routed through the retry handler below, which
                    # is the single place that decides when to give up.
                    raise requests.exceptions.HTTPError(
                        f"HTTP {response.status_code} {code or 'error'} for {url}"
                    )

                if response.status_code in (400, 404):
                    # PUGVIEW.NotFound (no data for this compound) or
                    # PUGVIEW.BadRequest (no such heading). Both permanent.
                    self.logger.debug(f"{code or response.status_code} for {url}")
                    raise PubChemViewNotFoundError(
                        f"No data for this request ({code or response.status_code}): {url}"
                    )

                # Any other client error is permanent; retrying cannot help.
                raise PubChemViewError(
                    f"HTTP {response.status_code} {code or ''} for {url}: "
                    f"{response.text[:200]}"
                )
                    
            except requests.exceptions.RequestException as e:
                if attempt == self.max_retries:
                    raise PubChemViewError(f"Request failed after {self.max_retries + 1} attempts: {e}")
                
                wait_time = self.backoff_factor * (2 ** attempt)
                self.logger.warning(f"Request failed, retrying in {wait_time:.1f}s: {e}")
                time.sleep(wait_time)
        
        # This should never be reached due to exceptions above
        raise PubChemViewError("Unexpected error in request handling")
        
    @cached(service='pubchemview')
    def get_experimental_properties(self, cid: Union[int, str]) -> Dict[str, Any]:
        """
        Get all experimental properties for a compound
        
        Args:
            cid: PubChem Compound ID
            
        Returns:
            Raw JSON response containing all experimental properties
        """
        url = f"{self.base_url}/data/compound/{cid}/JSON?heading=Experimental+Properties"
        return self._make_request(url)
    
    @cached(service='pubchemview')
    def get_property(self, cid: Union[int, str], property_name: str) -> Dict[str, Any]:
        """
        Get a specific property for a compound
        
        Args:
            cid: PubChem Compound ID
            property_name: Name of the property (can be with spaces or plus signs)
            
        Returns:
            Raw JSON response for the specific property
        """
        # Convert property name to URL-safe format
        if property_name in self.experimental_properties:
            url_property = self.experimental_properties[property_name]
        else:
            url_property = property_name.replace(' ', '+')
        
        url = f"{self.base_url}/data/compound/{cid}/JSON?heading={url_property}"
        return self._make_request(url)
    

    def _extract_value_info(self, info_item: Dict[str, Any],
                            heading: Optional[str] = None) -> PropertyData:
        """
        Turn one PUG-View ``Information`` item into a :class:`PropertyData`.

        Args:
            info_item: A single ``Information`` dictionary from a PUG-View
                response.
            heading: The section heading the item was found under. Passed
                through to :func:`~provesid.pubchemview_parse.parse_value`,
                which needs it to read a bare number: 138 under "Melting Point"
                is 138 °C, while 1.19 under "LogP" is dimensionless.

        Returns:
            A PropertyData carrying the original string, the parsed numbers and
            the reference it came from.
        """
        value_str = ""
        if "Value" in info_item and "StringWithMarkup" in info_item["Value"]:
            markup = info_item["Value"]["StringWithMarkup"]
            if markup:
                value_str = markup[0].get("String", "")

        # Some headings report a plain number instead of prose, with the unit
        # in its own field; the computed-property sections do this throughout.
        value_block = info_item.get("Value", {})
        if not value_str and value_block.get("Number"):
            number = value_block["Number"][0]
            unit = value_block.get("Unit", "")
            value_str = f"{number} {unit}".strip()

        reference = None
        reference_number = info_item.get("ReferenceNumber")
        if "Reference" in info_item and info_item["Reference"]:
            reference = info_item["Reference"][0]

        name = info_item.get("Name", "")
        description = info_item.get("Description", "")

        parsed = parse_value(value_str, heading)

        return PropertyData(
            value=value_str,
            unit=parsed.unit,
            conditions=parsed.conditions,
            reference=reference,
            reference_number=reference_number,
            description=description if description else None,
            name=name if name else None,
            heading=heading,
            parsed=parsed,
        )

    def _iter_information(self, response: Dict[str, Any]
                          ) -> Iterator[Tuple[Optional[str], Dict[str, Any]]]:
        """
        Walk a PUG-View record and yield every value with its heading.

        PUG-View nests a requested heading wherever it sits in the compound's
        table of contents, and that position differs by heading: "Melting Point"
        arrives under *Chemical and Physical Properties → Experimental
        Properties*, "GHS Classification" under *Safety and Hazards → Hazards
        Identification*, "Drug Indication" under *Drug and Medication
        Information*. Walking the tree instead of following one hard-coded path
        is what lets every heading work, rather than only the experimental ones.

        Args:
            response: Raw JSON response from PUG-View.

        Yields:
            ``(heading, information_item)`` pairs, the heading being the
            innermost ``TOCHeading`` that owns the item.
        """
        def walk(node: Dict[str, Any], inherited: Optional[str]):
            """Yield this node's values, then recurse into its subsections."""
            # A subsection without its own heading belongs to its parent's.
            heading = node.get("TOCHeading", inherited)
            for item in node.get("Information", []) or []:
                yield heading, item
            for child in node.get("Section", []) or []:
                yield from walk(child, heading)

        if not isinstance(response, dict):
            return
        record = response.get("Record", {})
        for section in record.get("Section", []) or []:
            yield from walk(section, None)

    @staticmethod
    def _carries_a_value(info_item: Dict[str, Any]) -> bool:
        """
        Whether an ``Information`` item states a value at all.

        Some items are pointers rather than measurements — PubChem records a
        compound's IUPAC pKa data as ``{"Value": {"ExternalTableName":
        "iupacpka"}}``, which carries a reference and nothing else. Including
        one would put a blank row in a property table.

        Args:
            info_item: A single ``Information`` dictionary.

        Returns:
            True when the item has a string or a number to parse.
        """
        value = info_item.get("Value", {})
        markup = value.get("StringWithMarkup") or []
        return bool((markup and markup[0].get("String")) or value.get("Number"))

    @cached(service='pubchemview', skip_if=_is_empty_lookup)
    def extract_property_data(self, cid: Union[int, str], property_name: str) -> List[PropertyData]:
        """
        Extract structured property data for a specific property.

        Args:
            cid: PubChem Compound ID
            property_name: The PUG-View heading, e.g. ``"Melting Point"``. Any
                heading works, not only the experimental ones.

        Returns:
            List of PropertyData objects, each carrying the value as PubChem
            wrote it plus a ``parsed`` :class:`ParsedValue` holding the numbers
            recovered from it. An empty list means PubChem holds no such
            property for this compound.

        Raises:
            PubChemViewError: If the request could not be completed, for example
                a PUG-View ``ServerBusy`` response that survived every retry.
                This is deliberately *not* reported as an empty list: a
                transient failure must stay distinguishable from real absence,
                or it gets cached as "no data" and never retried.

        Example:
            >>> view = PubChemView()
            >>> data = view.extract_property_data(2244, "Dissociation Constants")
            >>> data[0].value
            '3.47'
            >>> data[0].parsed.value
            3.47
        """
        try:
            response = self.get_property(cid, property_name)
        except PubChemViewNotFoundError:
            self.logger.debug(f"Property '{property_name}' not present for CID {cid}")
            return []
        return self._parse_property_response(response)
    
    def _parse_property_response(self, response: Dict[str, Any]) -> List[PropertyData]:
        """
        Parse a single-heading response into structured values.

        Args:
            response: Raw JSON response from PUG-View, as returned by
                :meth:`get_property`.

        Returns:
            One PropertyData per value in the response, in document order.
            Items that carry no value — PubChem's pointers to its own external
            tables — are left out rather than returned blank. An empty list
            means the response carried no values, not that parsing failed.
        """
        return [self._extract_value_info(item, heading)
                for heading, item in self._iter_information(response)
                if self._carries_a_value(item)]

    @cached(service='pubchemview', skip_if=_is_empty_lookup)
    def extract_all_experimental_properties(self, cid: Union[int, str]) -> Dict[str, List[PropertyData]]:
        """
        Extract all experimental properties for a compound in structured format
        
        Args:
            cid: PubChem Compound ID
            
        Returns:
            Dictionary mapping property names to lists of PropertyData objects
        """
        try:
            response = self.get_experimental_properties(cid)
            return self._parse_all_properties_response(response)
        except PubChemViewNotFoundError:
            self.logger.warning(f"No experimental properties found for CID {cid}")
            return {}
    
    def _parse_all_properties_response(self, response: Dict[str, Any]
                                       ) -> Dict[str, List[PropertyData]]:
        """
        Parse a multi-property response, grouped by heading.

        Args:
            response: Raw JSON response from PUG-View, as returned by
                :meth:`get_experimental_properties`.

        Returns:
            A dict mapping each heading that carried values to its
            PropertyData list, in the order PubChem listed them.
        """
        grouped: Dict[str, List[PropertyData]] = {}
        for heading, item in self._iter_information(response):
            if not self._carries_a_value(item):
                continue
            grouped.setdefault(heading or "Unknown", []).append(
                self._extract_value_info(item, heading))
        return grouped

    @cached(service='pubchemview', skip_if=_is_empty_lookup)
    def get_available_properties(self, cid: Union[int, str]) -> List[str]:
        """
        Get list of available experimental properties for a compound
        
        Args:
            cid: PubChem Compound ID
            
        Returns:
            List of available property names
        """
        try:
            response = self.get_experimental_properties(cid)
            return list(self._parse_all_properties_response(response).keys())
        except PubChemViewNotFoundError:
            return []
    
    @cached(service='pubchemview', skip_if=lambda result: not result.get('values'))
    def get_property_summary(self, cid: Union[int, str], property_name: str) -> Dict[str, Any]:
        """
        Get a summary of a property including all values, units, and references
        
        Args:
            cid: PubChem Compound ID
            property_name: Name of the property
            
        Returns:
            Dictionary with the raw strings under ``values``, the numbers
            recovered from them under ``numeric_values`` and
            ``numeric_values_si``, plus the references, units and conditions
            seen across the entries and a ``count``.
        """
        property_data = self.extract_property_data(cid, property_name)
        
        if not property_data:
            return {"property": property_name, "values": [], "numeric_values": [],
                    "numeric_values_si": [], "references": [], "units": [],
                    "units_si": [], "conditions": [], "count": 0}
        
        parsed_values = [data.parsed or parse_value(data.value, data.heading)
                         for data in property_data]

        summary = {
            "property": property_name,
            "values": [data.value for data in property_data],
            # The numbers behind those strings, for the common case of wanting
            # a range or a mean rather than the prose.
            "numeric_values": [parsed.value for parsed in parsed_values
                               if parsed.value is not None],
            "numeric_values_si": [parsed.value_si for parsed in parsed_values
                                  if parsed.value_si is not None],
            "references": [data.reference for data in property_data if data.reference],
            "units": list(set([data.unit for data in property_data if data.unit])),
            "units_si": list({parsed.unit_si for parsed in parsed_values
                              if parsed.unit_si}),
            "conditions": list(set([data.conditions for data in property_data if data.conditions])),
            "count": len(property_data)
        }
        
        return summary
    
    # Convenience methods for common properties
    @cached(service='pubchemview', skip_if=_is_empty_lookup)
    def get_melting_point(self, cid: Union[int, str]) -> List[PropertyData]:
        """Get melting point data for a compound"""
        return self.extract_property_data(cid, "Melting Point")
    
    @cached(service='pubchemview', skip_if=_is_empty_lookup)
    def get_boiling_point(self, cid: Union[int, str]) -> List[PropertyData]:
        """Get boiling point data for a compound"""
        return self.extract_property_data(cid, "Boiling Point")
    
    @cached(service='pubchemview', skip_if=_is_empty_lookup)
    def get_density(self, cid: Union[int, str]) -> List[PropertyData]:
        """Get density data for a compound"""
        return self.extract_property_data(cid, "Density")
    
    @cached(service='pubchemview', skip_if=_is_empty_lookup)
    def get_solubility(self, cid: Union[int, str]) -> List[PropertyData]:
        """Get solubility data for a compound"""
        return self.extract_property_data(cid, "Solubility")
    
    @cached(service='pubchemview', skip_if=_is_empty_lookup)
    def get_flash_point(self, cid: Union[int, str]) -> List[PropertyData]:
        """Get flash point data for a compound"""
        return self.extract_property_data(cid, "Flash Point")
    
    @cached(service='pubchemview', skip_if=_is_empty_lookup)
    def get_vapor_pressure(self, cid: Union[int, str]) -> List[PropertyData]:
        """Get vapor pressure data for a compound"""
        return self.extract_property_data(cid, "Vapor Pressure")
    
    @cached(service='pubchemview', skip_if=_is_empty_lookup)
    def get_viscosity(self, cid: Union[int, str]) -> List[PropertyData]:
        """Get viscosity data for a compound"""
        return self.extract_property_data(cid, "Viscosity")
    
    @cached(service='pubchemview', skip_if=_is_empty_lookup)
    def get_logp(self, cid: Union[int, str]) -> List[PropertyData]:
        """Get LogP data for a compound"""
        return self.extract_property_data(cid, "LogP")
    
    @cached(service='pubchemview', skip_if=_is_empty_lookup)
    def get_refractive_index(self, cid: Union[int, str]) -> List[PropertyData]:
        """Get refractive index data for a compound"""
        return self.extract_property_data(cid, "Refractive Index")
    
    def batch_extract_properties(self, cid: Union[int, str], 
                                property_names: List[str]) -> Dict[str, List[PropertyData]]:
        """
        Extract multiple properties for a compound
        
        Args:
            cid: PubChem Compound ID
            property_names: List of property names to extract

        Returns:
            Dictionary mapping property names to PropertyData lists.

        Note:
            One property failing does not abort the batch: that entry is logged
            at WARNING and comes back as an empty list, so here — unlike in
            :meth:`extract_property_data` — an empty list does not prove the
            property is absent. Call ``extract_property_data`` directly when the
            difference matters.
        """
        results = {}
        for prop_name in property_names:
            try:
                results[prop_name] = self.extract_property_data(cid, prop_name)
            except Exception as e:
                self.logger.warning(f"Failed to extract {prop_name} for CID {cid}: {e}")
                results[prop_name] = []
        
        return results
    
    def export_properties_to_dict(self, property_data_list: List[PropertyData]) -> List[Dict[str, Any]]:
        """
        Convert PropertyData objects to plain dictionaries for serialization.

        Args:
            property_data_list: List of PropertyData objects, as returned by
                :meth:`extract_property_data` or any of the convenience getters.

        Returns:
            One dict per value, carrying the original string, the reference, and
            the parsed numbers flattened into top-level keys
            (``numeric_value``, ``value_min``, ``value_si``, ``unit_si``,
            ``operator``, ``qualitative``, ``temperature_c``) so the result can
            go straight into ``json.dumps`` or ``pd.DataFrame``.
        """
        rows = []
        for data in property_data_list:
            parsed = data.parsed or parse_value(data.value, data.heading)
            rows.append({
                "value": data.value,
                "unit": data.unit,
                "conditions": data.conditions,
                "reference": data.reference,
                "reference_number": data.reference_number,
                "description": data.description,
                "name": data.name,
                "heading": data.heading,
                # The parsed numbers are flattened rather than nested, so the
                # result stays directly serialisable to JSON or a DataFrame.
                "numeric_value": parsed.value,
                "value_min": parsed.value_min,
                "value_max": parsed.value_max,
                "value_si": parsed.value_si,
                "value_min_si": parsed.value_min_si,
                "value_max_si": parsed.value_max_si,
                "unit_si": parsed.unit_si,
                "operator": parsed.operator,
                "qualitative": parsed.qualitative,
                "temperature_c": parsed.temperature_c,
            })
        return rows


    #: Columns of the frame :meth:`get_property_table` returns, in order.
    PROPERTY_TABLE_COLUMNS = [
        "CID", "Heading", "StringWithMarkup", "ExperimentalValue",
        "ValueMin", "ValueMax", "Unit", "ValueSI", "ValueMinSI", "ValueMaxSI",
        "UnitSI", "Operator", "Qualitative", "Temperature", "Conditions",
        "FullReference",
    ]

    @cached(service='pubchemview', skip_if=_is_empty_lookup)
    def get_property_table(self, cid: Union[int, str], property_name: str) -> pd.DataFrame:
        """
        Get a table of property values, parsed into numbers, with references.

        Args:
            cid: PubChem Compound ID.
            property_name: The PUG-View heading, e.g. ``"Melting Point"``. Any
                heading works, not only the experimental ones.

        Returns:
            A DataFrame with the columns in :attr:`PROPERTY_TABLE_COLUMNS`:

            - ``StringWithMarkup`` — the value exactly as PubChem wrote it.
            - ``ExperimentalValue`` — the single number, as a float; NaN when
              the entry reports a range or no number at all.
            - ``ValueMin`` / ``ValueMax`` — the bounds, equal to
              ``ExperimentalValue`` for a single value, so a numeric filter
              needs no special case for ranges.
            - ``Unit`` — normalised, so ``torr`` and ``mm Hg`` agree.
            - ``ValueSI`` / ``ValueMinSI`` / ``ValueMaxSI`` / ``UnitSI`` — the
              same quantity in SI units, or NaN where the unit has no
              unambiguous SI equivalent (a percentage, a ppm).
            - ``Operator`` — ``>``, ``<``, ``>=``, ``<=`` or ``~`` when the
              entry bounds the value rather than stating it.
            - ``Qualitative`` — the word that replaced the number, such as
              ``insoluble``.
            - ``Temperature`` — the temperature the measurement was made at, in
              °C. A *condition*, not the value.
            - ``FullReference`` — the resolved citation.

            The frame is empty, with these columns, when PubChem holds no such
            property for this compound.

        Raises:
            PubChemViewError: If the request could not be completed. An empty
                frame always means "no such data", never "the fetch failed".

        Example:
            >>> view = PubChemView()
            >>> table = view.get_property_table(2244, "Melting Point")
            >>> table[["ExperimentalValue", "Unit", "ValueSI", "UnitSI"]].head(1)
               ExperimentalValue Unit  ValueSI UnitSI
            0              135.0   °C   408.15      K
        """
        try:
            response = self.get_property(cid, property_name)
            reference_map = self._extract_reference_map(response)
            property_data = self._parse_property_response(response)
        except PubChemViewNotFoundError:
            self.logger.debug(f"Property '{property_name}' not present for CID {cid}")
            return pd.DataFrame(columns=self.PROPERTY_TABLE_COLUMNS)
        except PubChemViewError:
            # Transport failure: let it surface rather than pass an empty table
            # off as "this compound has no data".
            raise

        rows = []
        for data in property_data:
            full_reference = ""
            if data.reference_number and data.reference_number in reference_map:
                full_reference = reference_map[data.reference_number]
            elif data.reference:
                full_reference = data.reference

            parsed = data.parsed or parse_value(data.value, data.heading)
            rows.append({
                "CID": cid,
                "Heading": data.heading,
                "StringWithMarkup": data.value,
                "ExperimentalValue": parsed.value,
                "ValueMin": parsed.value_min,
                "ValueMax": parsed.value_max,
                "Unit": parsed.unit,
                "ValueSI": parsed.value_si,
                "ValueMinSI": parsed.value_min_si,
                "ValueMaxSI": parsed.value_max_si,
                "UnitSI": parsed.unit_si,
                "Operator": parsed.operator,
                "Qualitative": parsed.qualitative,
                "Temperature": parsed.temperature_c,
                "Conditions": parsed.conditions,
                "FullReference": full_reference,
            })

        return pd.DataFrame(rows, columns=self.PROPERTY_TABLE_COLUMNS)

    def _extract_reference_map(self, response: Dict[str, Any]) -> Dict[int, str]:
        """
        Extract mapping of reference numbers to full reference strings
        
        Args:
            response: Raw JSON response from PUG View
            
        Returns:
            Dictionary mapping reference numbers to full reference strings
        """
        reference_map = {}
        
        try:
            record = response.get("Record", {})
            references = record.get("Reference", [])
            
            for ref in references:
                ref_num = ref.get("ReferenceNumber")
                if ref_num:
                    # Build full reference string from available fields
                    ref_parts = []
                    
                    # Add source name
                    if "SourceName" in ref:
                        ref_parts.append(ref["SourceName"])
                    
                    # Add name/title
                    if "Name" in ref:
                        ref_parts.append(ref["Name"])
                    
                    # Add description
                    if "Description" in ref:
                        description = ref["Description"]
                        # Truncate very long descriptions
                        if len(description) > 200:
                            description = description[:200] + "..."
                        ref_parts.append(description)
                    
                    # Add URL if available
                    if "URL" in ref:
                        ref_parts.append(f"URL: {ref['URL']}")
                    
                    # Combine parts with proper separation
                    full_ref = " | ".join(ref_parts) if ref_parts else f"Reference #{ref_num}"
                    reference_map[ref_num] = full_ref
                    
        except Exception as e:
            self.logger.warning(f"Error extracting reference map: {e}")
            
        return reference_map
    

# Convenience functions for easy access.
#
# These are not cached themselves: each one delegates to a cached PubChemView
# method, so a second cache here would keep a duplicate copy of the same payload
# on disk — and one whose key does not carry PubChemView.CACHE_SCHEMA_VERSION,
# which is how an upgrade would end up serving an entry of the previous shape.
def get_experimental_property(cid: Union[int, str], property_name: str) -> List[PropertyData]:
    """
    Convenience function to get experimental property data
    
    Args:
        cid: PubChem Compound ID
        property_name: Name of the property
        
    Returns:
        List of PropertyData objects
    """
    pugview = PubChemView()
    return pugview.extract_property_data(cid, property_name)


def get_all_experimental_properties(cid: Union[int, str]) -> Dict[str, List[PropertyData]]:
    """
    Convenience function to get all experimental properties
    
    Args:
        cid: PubChem Compound ID
        
    Returns:
        Dictionary mapping property names to PropertyData lists
    """
    pugview = PubChemView()
    return pugview.extract_all_experimental_properties(cid)


def get_property_values_only(cid: Union[int, str], property_name: str) -> List[str]:
    """
    Convenience function to get just the property values as strings
    
    Args:
        cid: PubChem Compound ID
        property_name: Name of the property
        
    Returns:
        List of property value strings
    """
    pugview = PubChemView()
    property_data = pugview.extract_property_data(cid, property_name)
    return [data.value for data in property_data if data.value]


def get_property_table(cid: Union[int, str], property_name: str) -> pd.DataFrame:
    """
    Convenience function to get a comprehensive property table with full references
    
    Args:
        cid: PubChem Compound ID
        property_name: Name of the experimental property
        
    Returns:
        pandas DataFrame with the columns listed in
        :attr:`PubChemView.PROPERTY_TABLE_COLUMNS`.
    """
    pugview = PubChemView()
    return pugview.get_property_table(cid, property_name)
