import requests
import logging
import time
from urllib.parse import quote
from typing import Dict, List, Optional, Any, Union
from dataclasses import dataclass
import pandas as pd
import re


@dataclass
class PropertyData:
    """Data structure for holding extracted property information"""
    value: str
    unit: Optional[str] = None
    conditions: Optional[str] = None
    reference: Optional[str] = None
    reference_number: Optional[int] = None
    description: Optional[str] = None
    name: Optional[str] = None


class PubChemViewError(Exception):
    """Base exception class for PubChem View API errors"""
    pass


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
                 timeout: int = 30, max_retries: int = 3, backoff_factor: float = 1.0):
        """
        Initialize PubChemView API client
        
        Args:
            base_url: Base URL for PubChem PUG View API
            timeout: Request timeout in seconds
            max_retries: Maximum number of retry attempts
            backoff_factor: Backoff factor for retries
        """
        self.base_url = base_url.rstrip('/')
        self.timeout = timeout
        self.max_retries = max_retries
        self.backoff_factor = backoff_factor
        self.logger = logging.getLogger(__name__)
        
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
            PubChemViewError: For API errors
            PubChemViewNotFoundError: When resource not found
        """
        self._rate_limit()
        
        for attempt in range(self.max_retries + 1):
            try:
                self.logger.debug(f"Making request to: {url}")
                response = requests.get(url, timeout=self.timeout)
                
                if response.status_code == 200:
                    return response.json()
                elif response.status_code == 404:
                    raise PubChemViewNotFoundError(f"Resource not found: {url}")
                else:
                    response.raise_for_status()
                    
            except requests.exceptions.RequestException as e:
                if attempt == self.max_retries:
                    raise PubChemViewError(f"Request failed after {self.max_retries + 1} attempts: {e}")
                
                wait_time = self.backoff_factor * (2 ** attempt)
                self.logger.warning(f"Request failed, retrying in {wait_time:.1f}s: {e}")
                time.sleep(wait_time)
        
        # This should never be reached due to exceptions above
        raise PubChemViewError("Unexpected error in request handling")
        
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
    

    def _extract_value_info(self, info_item: Dict[str, Any]) -> PropertyData:
        """
        Extract structured information from an Information item
        
        Args:
            info_item: Single Information dictionary from PUG View response
            
        Returns:
            PropertyData object with extracted information
        """
        # Extract value
        value_str = ""
        if "Value" in info_item and "StringWithMarkup" in info_item["Value"]:
            value_str = info_item["Value"]["StringWithMarkup"][0]["String"]
        
        # Extract reference
        reference = None
        reference_number = info_item.get("ReferenceNumber")
        if "Reference" in info_item and info_item["Reference"]:
            reference = info_item["Reference"][0]
        
        # Extract name/description if available
        name = info_item.get("Name", "")
        description = info_item.get("Description", "")
        
        # Attempt to parse units and conditions from value string
        unit, conditions = self._parse_value_string(value_str)
        
        return PropertyData(
            value=value_str,
            unit=unit,
            conditions=conditions,
            reference=reference,
            reference_number=reference_number,
            description=description if description else None,
            name=name if name else None
        )
    
    def _parse_value_string(self, value_str: str) -> tuple[Optional[str], Optional[str]]:
        """
        Attempt to parse units and conditions from a value string
        
        Args:
            value_str: Value string from PUG View
            
        Returns:
            Tuple of (unit, conditions)
        """
        import re
        
        # Common temperature patterns
        temp_pattern = r'at\s+(-?\d+(?:\.\d+)?)\s*°?C'
        temp_match = re.search(temp_pattern, value_str)
        conditions = None
        if temp_match:
            conditions = f"at {temp_match.group(1)}°C"
        
        # Common unit patterns
        unit_patterns = [
            r'(\d+(?:\.\d+)?)\s*(cP|mPa·s|Pa·s)', # viscosity
            r'(\d+(?:\.\d+)?)\s*(°C|K)', # temperature  
            r'(\d+(?:\.\d+)?)\s*(g/cm³|g/mL|kg/m³)', # density
            r'(\d+(?:\.\d+)?)\s*(mmHg|kPa|Pa|atm|bar)', # pressure
            r'(\d+(?:\.\d+)?)\s*(mN/m|N/m|dyn/cm)', # surface tension
            r'(\d+(?:\.\d+)?)\s*(g/L|mg/L|%|ppm)', # solubility/concentration
        ]
        
        unit = None
        for pattern in unit_patterns:
            match = re.search(pattern, value_str)
            if match:
                unit = match.group(2)
                break
        
        return unit, conditions
    
    def extract_property_data(self, cid: Union[int, str], property_name: str) -> List[PropertyData]:
        """
        Extract structured property data for a specific property
        
        Args:
            cid: PubChem Compound ID
            property_name: Name of the property to extract
            
        Returns:
            List of PropertyData objects with extracted information
        """
        try:
            response = self.get_property(cid, property_name)
            return self._parse_property_response(response)
        except (PubChemViewNotFoundError, PubChemViewError):
            self.logger.warning(f"Property '{property_name}' not found for CID {cid}")
            return []
    
    def _parse_property_response(self, response: Dict[str, Any]) -> List[PropertyData]:
        """
        Parse a property response and extract structured data
        
        Args:
            response: Raw JSON response from PUG View
            
        Returns:
            List of PropertyData objects
        """
        property_data = []
        
        try:
            # Navigate to the Information section
            record = response.get("Record", {})
            sections = record.get("Section", [])
            
            # Find the experimental properties section
            for section in sections:
                if section.get("TOCHeading") == "Chemical and Physical Properties":
                    exp_sections = section.get("Section", [])
                    for exp_section in exp_sections:
                        if exp_section.get("TOCHeading") == "Experimental Properties":
                            prop_sections = exp_section.get("Section", [])
                            
                            # Find the specific property section
                            for prop_section in prop_sections:
                                information_items = prop_section.get("Information", [])
                                for info_item in information_items:
                                    property_data.append(self._extract_value_info(info_item))
                                        
        except Exception as e:
            self.logger.error(f"Error parsing property response: {e}")
            
        return property_data
    
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
    
    def _parse_all_properties_response(self, response: Dict[str, Any]) -> Dict[str, List[PropertyData]]:
        """
        Parse a full experimental properties response
        
        Args:
            response: Raw JSON response from PUG View
            
        Returns:
            Dictionary mapping property names to PropertyData lists
        """
        all_properties = {}
        
        try:
            # Navigate to the experimental properties section
            record = response.get("Record", {})
            sections = record.get("Section", [])
            
            for section in sections:
                if section.get("TOCHeading") == "Chemical and Physical Properties":
                    exp_sections = section.get("Section", [])
                    for exp_section in exp_sections:
                        if exp_section.get("TOCHeading") == "Experimental Properties":
                            prop_sections = exp_section.get("Section", [])
                            
                            # Extract each property
                            for prop_section in prop_sections:
                                prop_name = prop_section.get("TOCHeading", "Unknown")
                                property_data = []
                                
                                information_items = prop_section.get("Information", [])
                                for info_item in information_items:
                                    property_data.append(self._extract_value_info(info_item))
                                
                                if property_data:
                                    all_properties[prop_name] = property_data
                                        
        except Exception as e:
            self.logger.error(f"Error parsing all properties response: {e}")
            
        return all_properties
    
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
    
    def get_property_summary(self, cid: Union[int, str], property_name: str) -> Dict[str, Any]:
        """
        Get a summary of a property including all values, units, and references
        
        Args:
            cid: PubChem Compound ID
            property_name: Name of the property
            
        Returns:
            Dictionary with property summary
        """
        property_data = self.extract_property_data(cid, property_name)
        
        if not property_data:
            return {"property": property_name, "values": [], "references": [], "units": set()}
        
        summary = {
            "property": property_name,
            "values": [data.value for data in property_data],
            "references": [data.reference for data in property_data if data.reference],
            "units": list(set([data.unit for data in property_data if data.unit])),
            "conditions": list(set([data.conditions for data in property_data if data.conditions])),
            "count": len(property_data)
        }
        
        return summary
    
    # Convenience methods for common properties
    def get_melting_point(self, cid: Union[int, str]) -> List[PropertyData]:
        """Get melting point data for a compound"""
        return self.extract_property_data(cid, "Melting Point")
    
    def get_boiling_point(self, cid: Union[int, str]) -> List[PropertyData]:
        """Get boiling point data for a compound"""
        return self.extract_property_data(cid, "Boiling Point")
    
    def get_density(self, cid: Union[int, str]) -> List[PropertyData]:
        """Get density data for a compound"""
        return self.extract_property_data(cid, "Density")
    
    def get_solubility(self, cid: Union[int, str]) -> List[PropertyData]:
        """Get solubility data for a compound"""
        return self.extract_property_data(cid, "Solubility")
    
    def get_flash_point(self, cid: Union[int, str]) -> List[PropertyData]:
        """Get flash point data for a compound"""
        return self.extract_property_data(cid, "Flash Point")
    
    def get_vapor_pressure(self, cid: Union[int, str]) -> List[PropertyData]:
        """Get vapor pressure data for a compound"""
        return self.extract_property_data(cid, "Vapor Pressure")
    
    def get_viscosity(self, cid: Union[int, str]) -> List[PropertyData]:
        """Get viscosity data for a compound"""
        return self.extract_property_data(cid, "Viscosity")
    
    def get_logp(self, cid: Union[int, str]) -> List[PropertyData]:
        """Get LogP data for a compound"""
        return self.extract_property_data(cid, "LogP")
    
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
            Dictionary mapping property names to PropertyData lists
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
        Convert PropertyData objects to dictionaries for easy serialization
        
        Args:
            property_data_list: List of PropertyData objects
            
        Returns:
            List of dictionaries
        """
        return [
            {
                "value": data.value,
                "unit": data.unit,
                "conditions": data.conditions,
                "reference": data.reference,
                "reference_number": data.reference_number,
                "description": data.description,
                "name": data.name
            }
            for data in property_data_list
        ]


    def get_property_table(self, cid: Union[int, str], property_name: str) -> pd.DataFrame:
        """
        Get a comprehensive table of property data with full reference information
        
        Args:
            cid: PubChem Compound ID
            property_name: Name of the experimental property
            
        Returns:
            pandas DataFrame with columns: CID, StringWithMarkup, ExperimentalValue, Unit, Temperature, Conditions, FullReference
        """
        try:
            # Get the raw response to extract full reference information
            response = self.get_property(cid, property_name)
            
            # Extract reference mapping from the response
            reference_map = self._extract_reference_map(response)
            
            # Get structured property data
            property_data = self._parse_property_response(response)
            
            # Build table data
            table_data = []
            for data in property_data:
                # Get full reference string
                full_reference = ""
                if data.reference_number and data.reference_number in reference_map:
                    full_reference = reference_map[data.reference_number]
                elif data.reference:
                    full_reference = data.reference
                
                # Parse experimental value, unit, temperature, and conditions from the StringWithMarkup
                exp_value, unit, temperature, conditions = self._extract_experimental_value_and_unit(data.value, property_name)
                
                table_data.append({
                    "CID": cid,
                    "StringWithMarkup": data.value,
                    "ExperimentalValue": exp_value,
                    "Unit": unit,  # Use only the parsed unit from improved extraction
                    "Temperature": temperature,
                    "Conditions": conditions,
                    "FullReference": full_reference
                })
            
            return pd.DataFrame(table_data)
            
        except Exception as e:
            self.logger.error(f"Error creating property table for CID {cid}, property {property_name}: {e}")
            # Return empty DataFrame with expected columns
            return pd.DataFrame(columns=["CID", "StringWithMarkup", "ExperimentalValue", "Unit", "Temperature", "Conditions", "FullReference"])
    
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
    
    def _extract_experimental_value_and_unit(self, value_str: str, property_name: str = None) -> tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
        """
        Extract numerical experimental value, unit, temperature, and conditions from StringWithMarkup
        
        Args:
            value_str: Full StringWithMarkup value
            property_name: Name of the property being extracted (for property-specific patterns)
            
        Returns:
            Tuple of (experimental_value, unit, temperature, conditions)
        """
        if not value_str:
            return None, None, None, None
        
        # Helper function to extract temperature and conditions
        def extract_temperature_and_conditions(text: str) -> tuple[Optional[str], Optional[str]]:
            """Extract temperature and conditions from text"""
            # Pattern for "at 25°C", "@ 25°C", "at 20 °C", etc.
            temp_match = re.search(r'(?:at|@)\s+(-?\d+(?:\.\d+)?)\s*[°]?\s*([CF]|K)\b', text, re.IGNORECASE)
            if temp_match:
                temp_value = temp_match.group(1)
                temp_unit = temp_match.group(2)
                temperature = f"{temp_value}°{temp_unit}"
                
                # Extract any additional conditions
                conditions_parts = []
                # Look for pressure conditions
                pressure_cond = re.search(r'(\d+(?:\.\d+)?)\s*(mmHg|kPa|Pa|atm|bar)', text, re.IGNORECASE)
                if pressure_cond:
                    conditions_parts.append(f"{pressure_cond.group(1)} {pressure_cond.group(2)}")
                
                # Look for other descriptive conditions
                descriptors = re.findall(r'/([^/]+)/', text)
                conditions_parts.extend(descriptors)
                
                conditions = "; ".join(conditions_parts) if conditions_parts else None
                return temperature, conditions
            
            return None, None
        
        # Property-specific patterns (highest priority)
        if property_name:
            prop_name_lower = property_name.lower()
            
            # Vapor Pressure specific patterns
            if 'vapor pressure' in prop_name_lower:
                # Pattern: "Vapor pressure at 20°C: negligible" - return None for both value and unit
                negligible_match = re.search(r'negligible', value_str, re.IGNORECASE)
                if negligible_match:
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return None, None, temperature, conditions
                
                # Pattern: "8.5X10-5 mm Hg at 25 °C" (X represents multiplication)
                scientific_x_match = re.search(r'^(\d+(?:\.\d+)?)X10([+-]?\d+)\s*(mmHg|mm\s+Hg|kPa|Pa|atm|bar|torr)', value_str, re.IGNORECASE)
                if scientific_x_match:
                    # Convert XNotation to E notation
                    mantissa = scientific_x_match.group(1)
                    exponent = scientific_x_match.group(2)
                    unit = scientific_x_match.group(3)
                    scientific_value = f"{mantissa}e{exponent}"
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return scientific_value, unit, temperature, conditions
                
                # Pattern: "2.7X10+0 at 25 °C /Estimated/" (X notation without unit)
                scientific_x_no_unit_match = re.search(r'^(\d+(?:\.\d+)?)X10([+-]?\d+)\s', value_str, re.IGNORECASE)
                if scientific_x_no_unit_match:
                    # Convert XNotation to E notation
                    mantissa = scientific_x_no_unit_match.group(1)
                    exponent = scientific_x_no_unit_match.group(2)
                    scientific_value = f"{mantissa}e{exponent}"
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return scientific_value, None, temperature, conditions
                
                # Pattern: "0.05 [mmHg]" (brackets around unit)
                bracketed_unit_match = re.search(r'^(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s*\[(mmHg|mm\s+Hg|kPa|Pa|atm|bar|torr)\]', value_str, re.IGNORECASE)
                if bracketed_unit_match:
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return bracketed_unit_match.group(1), bracketed_unit_match.group(2), temperature, conditions
                
                # Pattern: "Vapor pressure, kPa at 20°C: 24"
                vp_colon_match = re.search(r'vapor\s+pressure[^:]*:\s*(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)', value_str, re.IGNORECASE)
                if vp_colon_match:
                    # Extract unit from before the colon
                    unit_match = re.search(r'vapor\s+pressure[,\s]*([a-zA-Z]+)', value_str, re.IGNORECASE)
                    unit = unit_match.group(1) if unit_match else None
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return vp_colon_match.group(1), unit, temperature, conditions
                
                # Pattern: "kPa at 20°C: 24" or similar
                unit_colon_match = re.search(r'(mmHg|mm\s+Hg|kPa|Pa|atm|bar|torr)[^:]*:\s*(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)', value_str, re.IGNORECASE)
                if unit_colon_match:
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return unit_colon_match.group(2), unit_colon_match.group(1), temperature, conditions
                
                # Standard pressure patterns at start
                pressure_match = re.search(r'^(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s*(mmHg|mm\s+Hg|kPa|Pa|atm|bar|torr)\b', value_str, re.IGNORECASE)
                if pressure_match:
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return pressure_match.group(1), pressure_match.group(2), temperature, conditions
                
                # If no vapor pressure pattern matched, return None (no fallback for property-specific extraction)
                return None, None, None, None
            
            # LogP specific patterns
            elif 'logp' in prop_name_lower:
                # Pattern: "LogP: -2.3" or "log P = 1.5" or "log Kow = 1.19"
                logp_colon_match = re.search(r'log\s*(?:p|kow|k[ow]{1,2})\s*[=:]\s*(-?\d+(?:\.\d+)?)', value_str, re.IGNORECASE)
                if logp_colon_match:
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return logp_colon_match.group(1), None, temperature, conditions
                
                # Simple number at start for LogP (usually unitless)
                logp_start_match = re.search(r'^(-?\d+(?:\.\d+)?)\s*$', value_str.strip())
                if logp_start_match:
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return logp_start_match.group(1), None, temperature, conditions
                
                # If no LogP pattern matched, return None (no fallback for property-specific extraction)
                return None, None, None, None
            
            # Dissociation Constants specific patterns
            elif 'dissociation' in prop_name_lower:
                # Pattern: "pKa = 14.31 @ 25 °C" or "pKa: 3.6" or "pKb = 10.2"
                pka_colon_match = re.search(r'p[Kk][abAB]\s*[=:@]\s*(-?\d+(?:\.\d+)?)', value_str, re.IGNORECASE)
                if pka_colon_match:
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return pka_colon_match.group(1), None, temperature, conditions
                
                # Pattern: "Ka: 2.5e-4" or "Kb = 1.0e-10" or "K1=3.3X10-5"
                ka_colon_match = re.search(r'[Kk][abAB\d]*[=:]\s*([^\s;]+)', value_str, re.IGNORECASE)
                if ka_colon_match:
                    # Convert X notation to e notation if present (X10 -> e)
                    value = re.sub(r'([Xx])10', r'e', ka_colon_match.group(1))
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return value, None, temperature, conditions
                
                # Pattern: "2.91 (at 25 °C)" - just a number with temperature in parentheses
                number_with_temp_match = re.search(r'^([^\s;]+)\s*\(.*?°.*?\)', value_str.strip())
                if number_with_temp_match:
                    # Convert X notation to e notation if present (X10 -> e)
                    value = re.sub(r'([Xx])10', r'e', number_with_temp_match.group(1))
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return value, None, temperature, conditions
                
                # Simple number at start for dissociation constants (usually unitless pKa/pKb values)
                dc_start_match = re.search(r'^([^\s;]+)\s*$', value_str.strip())
                if dc_start_match:
                    # Convert X notation to e notation if present (X10 -> e)
                    value = re.sub(r'([Xx])10', r'e', dc_start_match.group(1))
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return value, None, temperature, conditions
                
                # If no dissociation constants pattern matched, return None (no fallback for property-specific extraction)
                return None, None, None, None
            
            # Melting Point / Boiling Point specific patterns
            elif any(temp_prop in prop_name_lower for temp_prop in ['melting point', 'boiling point', 'temperature']):
                # Range patterns with temperature units
                temp_range_match = re.search(r'^(-?\d+(?:\.\d+)?-\d+(?:\.\d+)?)\s*[°]?\s*([CF]|K)\b', value_str)
                if temp_range_match:
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return temp_range_match.group(1), f"°{temp_range_match.group(2)}", temperature, conditions
                
                # Single temperature values
                temp_match = re.search(r'^(-?\d+(?:\.\d+)?)\s*[°]?\s*([CF]|K)\b', value_str)
                if temp_match:
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return temp_match.group(1), f"°{temp_match.group(2)}", temperature, conditions
                
                # If no temperature pattern matched, return None (no fallback for property-specific extraction)
                return None, None, None, None
            
            # Density specific patterns
            elif 'density' in prop_name_lower:
                density_match = re.search(r'^(-?\d+(?:\.\d+)?)\s*(g/cm³|g/mL|kg/m³|g/L)\b', value_str)
                if density_match:
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return density_match.group(1), density_match.group(2), temperature, conditions
                
                # If no density pattern matched, return None (no fallback for property-specific extraction)
                return None, None, None, None
            
            # Viscosity specific patterns
            elif 'viscosity' in prop_name_lower:
                viscosity_match = re.search(r'^(-?\d+(?:\.\d+)?)\s*(cP|mPa·s|Pa·s|cSt)\b', value_str)
                if viscosity_match:
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return viscosity_match.group(1), viscosity_match.group(2), temperature, conditions
                
                # If no viscosity pattern matched, return None (no fallback for property-specific extraction)
                return None, None, None, None
            
            # Solubility specific patterns
            elif 'solubility' in prop_name_lower:
                # Pattern: "greater than or equal to 100 mg/mL" - comparison operators
                # Split into two separate searches for clarity
                sol_comparison_match = re.search(r'(?:greater than or equal to|greater than|less than or equal to|less than|≥|≤|>|<|>=|<=)\s*(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s*(μg|ug|ng|pg|g|mg|kg)/(mL|L|l|100mL|100ml|dl|dL)', value_str, re.IGNORECASE)
                if not sol_comparison_match:
                    sol_comparison_match = re.search(r'(?:greater than or equal to|greater than|less than or equal to|less than|≥|≤|>|<|>=|<=)\s*(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s*(g/L|g/l|mg/L|mg/l|g/100mL|mg/mL|mol/L|M|%|ppm|ppb)', value_str, re.IGNORECASE)
                
                if sol_comparison_match:
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    value = sol_comparison_match.group(1)
                    unit = sol_comparison_match.group(2)
                    
                    # Handle the compound unit format (e.g., μg/mL) vs simple unit format (e.g., g/L)
                    if '/' in unit:
                        # Simple unit format like g/L, mg/L
                        return value, unit, temperature, conditions
                    else:
                        # Compound unit format - need to get the volume part from the original match
                        if sol_comparison_match.lastindex >= 3:
                            volume_part = sol_comparison_match.group(3) if sol_comparison_match.group(3) else ''
                            if unit.lower() in ['μg', 'ug']:
                                unit = 'μg'
                            final_unit = f"{unit}/{volume_part}" if volume_part else unit
                            return value, final_unit, temperature, conditions
                        else:
                            return value, unit, temperature, conditions
                
                # Pattern: "1.2 [ug/mL] (additional info)" - bracketed units with microgram notation
                sol_bracketed_unit_match = re.search(r'(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s*\[(μg|ug|ng|pg|g|mg|kg)/(?:mL|L|100mL|100ml)\]', value_str, re.IGNORECASE)
                if sol_bracketed_unit_match:
                    value = sol_bracketed_unit_match.group(1)
                    unit_part = sol_bracketed_unit_match.group(2)
                    # Normalize microgram notation
                    if unit_part.lower() in ['μg', 'ug']:
                        unit_part = 'μg'
                    elif unit_part.lower() == 'ng':
                        unit_part = 'ng'
                    elif unit_part.lower() == 'pg':
                        unit_part = 'pg'
                    
                    # Extract volume part from the match
                    if '/mL' in sol_bracketed_unit_match.group(0):
                        unit = f"{unit_part}/mL"
                    elif '/L' in sol_bracketed_unit_match.group(0):
                        unit = f"{unit_part}/L"
                    elif '/100mL' in sol_bracketed_unit_match.group(0) or '/100ml' in sol_bracketed_unit_match.group(0):
                        unit = f"{unit_part}/100mL"
                    else:
                        unit = sol_bracketed_unit_match.group(2)
                    
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return value, unit, temperature, conditions
                
                # Pattern: "Soluble in water: 5.6 g/L at 20°C"
                sol_colon_match = re.search(r'[sS]olub[a-z]*[^:]*:\s*(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s*(μg|ug|ng|pg|g|mg|kg)/(mL|L|100mL|100ml|dl|dL)|[sS]olub[a-z]*[^:]*:\s*(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s*(g/L|mg/L|g/100mL|mg/mL|mol/L|M|%|ppm|ppb)\b', value_str, re.IGNORECASE)
                if sol_colon_match:
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    if sol_colon_match.group(1):  # First pattern (microgram notation)
                        unit_part = sol_colon_match.group(2)
                        if unit_part.lower() in ['μg', 'ug']:
                            unit_part = 'μg'
                        volume_part = sol_colon_match.group(3)
                        return sol_colon_match.group(1), f"{unit_part}/{volume_part}", temperature, conditions
                    else:  # Second pattern (standard units)
                        return sol_colon_match.group(4), sol_colon_match.group(5), temperature, conditions
                
                # Pattern: "5.6 g/L at 20°C" (starts with value and unit) - updated with microgram support
                sol_start_match = re.search(r'^(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s*(μg|ug|ng|pg|g|mg|kg)/(mL|L|100mL|100ml|dl|dL)|^(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s*(g/L|mg/L|g/100mL|mg/mL|mol/L|M|%|ppm|ppb)\b', value_str, re.IGNORECASE)
                if sol_start_match:
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    if sol_start_match.group(1):  # First pattern (microgram notation)
                        unit_part = sol_start_match.group(2)
                        if unit_part.lower() in ['μg', 'ug']:
                            unit_part = 'μg'
                        volume_part = sol_start_match.group(3)
                        return sol_start_match.group(1), f"{unit_part}/{volume_part}", temperature, conditions
                    else:  # Second pattern (standard units)
                        return sol_start_match.group(4), sol_start_match.group(5), temperature, conditions
                
                # Pattern: "2.5X10-3 g/L" (X notation)
                sol_scientific_x_match = re.search(r'^(-?\d+(?:\.\d+)?)X10([+-]?\d+)\s*(g/L|mg/L|g/100mL|mg/mL|mol/L|M|%|ppm|ppb)\b', value_str, re.IGNORECASE)
                if sol_scientific_x_match:
                    # Convert XNotation to E notation
                    mantissa = sol_scientific_x_match.group(1)
                    exponent = sol_scientific_x_match.group(2)
                    unit = sol_scientific_x_match.group(3)
                    scientific_value = f"{mantissa}e{exponent}"
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return scientific_value, unit, temperature, conditions
                
                # Pattern: "[2.5] g/L" (brackets around value)
                sol_bracketed_val_match = re.search(r'^\[(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\]\s*(g/L|mg/L|g/100mL|mg/mL|mol/L|M|%|ppm|ppb)\b', value_str, re.IGNORECASE)
                if sol_bracketed_val_match:
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return sol_bracketed_val_match.group(1), sol_bracketed_val_match.group(2), temperature, conditions
                
                # Pattern: "1.35X10+5 mg/l" appearing anywhere in text (X notation anywhere)
                sol_scientific_x_anywhere_match = re.search(r'(-?\d+(?:\.\d+)?)X10([+-]?\d+)\s*(mg/l|g/L|mg/L|g/100mL|mg/mL|mol/L|M|%|ppm|ppb)\b', value_str, re.IGNORECASE)
                if sol_scientific_x_anywhere_match:
                    # Convert XNotation to E notation
                    mantissa = sol_scientific_x_anywhere_match.group(1)
                    exponent = sol_scientific_x_anywhere_match.group(2)
                    unit = sol_scientific_x_anywhere_match.group(3)
                    scientific_value = f"{mantissa}e{exponent}"
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return scientific_value, unit, temperature, conditions
                
                # Pattern: "0.9%" or "0.86% wt" (standalone percentage values)
                sol_percentage_match = re.search(r'^(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)%(?:\s+wt|w/w|v/v)?$', value_str, re.IGNORECASE)
                if sol_percentage_match:
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return sol_percentage_match.group(1), "%", temperature, conditions
                
                # Pattern: "Solubility in water: 0.86% wt" (colon-based percentage)
                sol_colon_percentage_match = re.search(r'[sS]olub[a-z]*[^:]*:\s*(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)%(?:\s+wt|w/w|v/v)?$', value_str, re.IGNORECASE)
                if sol_colon_percentage_match:
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return sol_colon_percentage_match.group(1), "%", temperature, conditions
                
                # Pattern: "Insoluble" or "slightly soluble" - return None for both value and unit
                insoluble_match = re.search(r'(insoluble|practically insoluble|very slightly soluble)', value_str, re.IGNORECASE)
                if insoluble_match:
                    temperature, conditions = extract_temperature_and_conditions(value_str)
                    return None, None, temperature, conditions
                
                # If no solubility pattern matched, return None (no fallback for property-specific extraction)
                return None, None, None, None
            
            # For other property types, return None to avoid fallback to generic patterns
            return None, None, None, None
        
        # General patterns (fallback only when no property name is specified)
        
        # 1. Range patterns at start (e.g., "138-140", "135-140 °C")
        range_match = re.search(r'^(-?\d+(?:\.\d+)?-\d+(?:\.\d+)?)\s*([°]?[CF]|K)?\b', value_str)
        if range_match:
            temperature, conditions = extract_temperature_and_conditions(value_str)
            return range_match.group(1), range_match.group(2), temperature, conditions
        
        # 2. Temperature patterns at start (main value, not conditions)
        temp_match = re.search(r'^(-?\d+(?:\.\d+)?)\s*[°]?([CF]|K)\b', value_str)
        if temp_match:
            temperature, conditions = extract_temperature_and_conditions(value_str)
            return temp_match.group(1), temp_match.group(2), temperature, conditions
        
        # 3. Density patterns
        density_match = re.search(r'^(-?\d+(?:\.\d+)?)\s*(g/cm³|g/mL|kg/m³|g/L)\b', value_str)
        if density_match:
            temperature, conditions = extract_temperature_and_conditions(value_str)
            return density_match.group(1), density_match.group(2), temperature, conditions
        
        # 4. Pressure patterns (not conditions)
        pressure_match = re.search(r'^(-?\d+(?:\.\d+)?)\s*(mmHg|mm Hg|kPa|Pa|atm|bar|torr)\b', value_str)
        if pressure_match:
            temperature, conditions = extract_temperature_and_conditions(value_str)
            return pressure_match.group(1), pressure_match.group(2), temperature, conditions
        
        # 5. Viscosity patterns
        viscosity_match = re.search(r'^(\d+(?:\.\d+)?)\s*(cP|mPa·s|Pa·s|cSt)\b', value_str)
        if viscosity_match:
            temperature, conditions = extract_temperature_and_conditions(value_str)
            return viscosity_match.group(1), viscosity_match.group(2), temperature, conditions
        
        # 6. General numeric value at the start with various units
        general_match = re.search(r'^(-?\d+(?:\.\d+)?)\s*([a-zA-Z/²³·°%]+)?', value_str)
        if general_match:
            value = general_match.group(1)
            unit = general_match.group(2) if general_match.group(2) else None
            
            # Clean up unit
            if unit:
                unit = unit.strip()
                # Filter out common non-unit words that might be captured
                if unit.lower() in ['at', 'in', 'on', 'to', 'from', 'with', 'and', 'or']:
                    unit = None
                elif not unit or unit in ['', ' ']:
                    unit = None
            
            temperature, conditions = extract_temperature_and_conditions(value_str)
            return value, unit, temperature, conditions
        
        # 7. Scientific notation
        sci_match = re.search(r'^(-?\d+(?:\.\d+)?[Ee][+-]?\d+)\s*([a-zA-Z/²³·°%]+)?', value_str)
        if sci_match:
            temperature, conditions = extract_temperature_and_conditions(value_str)
            return sci_match.group(1), sci_match.group(2), temperature, conditions
        
        # 8. Fallback: any number in the string
        number_match = re.search(r'(-?\d+(?:\.\d+)?)', value_str)
        if number_match:
            temperature, conditions = extract_temperature_and_conditions(value_str)
            return number_match.group(1), None, temperature, conditions
        
        return None, None, None, None


# Convenience functions for easy access
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
        pandas DataFrame with columns: CID, StringWithMarkup, ExperimentalValue, Unit, FullReference
    """
    pugview = PubChemView()
    return pugview.get_property_table(cid, property_name)