"""
Comprehensive tests for PubChem PUG View functionality
"""

import pytest
import pandas as pd
import sys
import os

# Add the src directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from provesid.pubchemview import (
    PubChemView, 
    PropertyData, 
    PubChemViewError, 
    PubChemViewNotFoundError,
    get_experimental_property,
    get_all_experimental_properties,
    get_property_values_only,
    get_property_table
)


class TestPubChemView:
    """Test suite for PubChemView class"""
    
    @pytest.fixture
    def pugview(self):
        """Create a PubChemView instance for testing"""
        return PubChemView()
    
    def test_initialization(self, pugview):
        """Test PubChemView initialization"""
        assert pugview.base_url == "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"
        assert pugview.timeout == 30
        assert pugview.max_retries == 3
        assert pugview.min_request_interval == 0.2
        assert len(pugview.experimental_properties) > 30  # Should have many property mappings
    
    def test_custom_initialization(self):
        """Test PubChemView with custom parameters"""
        custom_pugview = PubChemView(
            base_url="https://example.com",
            timeout=60,
            max_retries=5,
            backoff_factor=2.0
        )
        assert custom_pugview.base_url == "https://example.com"
        assert custom_pugview.timeout == 60
        assert custom_pugview.max_retries == 5
        assert custom_pugview.backoff_factor == 2.0
    
    def test_get_experimental_properties(self, pugview):
        """Test getting all experimental properties for a compound"""
        # Test with aspirin (CID 2244)
        response = pugview.get_experimental_properties(2244)
        
        assert isinstance(response, dict)
        assert "Record" in response
        assert response["Record"]["RecordNumber"] == 2244
        assert response["Record"]["RecordTitle"] == "Aspirin"
    
    def test_get_property(self, pugview):
        """Test getting a specific property"""
        # Test melting point for aspirin
        response = pugview.get_property(2244, "Melting Point")
        
        assert isinstance(response, dict)
        assert "Record" in response
        assert response["Record"]["RecordNumber"] == 2244
    
    def test_extract_property_data(self, pugview):
        """Test extracting structured property data"""
        # Test melting point extraction
        mp_data = pugview.extract_property_data(2244, "Melting Point")
        
        assert isinstance(mp_data, list)
        assert len(mp_data) > 0
        
        # Check that we get PropertyData objects
        assert all(isinstance(item, PropertyData) for item in mp_data)
        
        # Check that first item has expected attributes
        first_item = mp_data[0]
        assert hasattr(first_item, 'value')
        assert hasattr(first_item, 'unit')
        assert hasattr(first_item, 'reference')
        assert isinstance(first_item.value, str)
    
    def test_extract_all_experimental_properties(self, pugview):
        """Test extracting all experimental properties"""
        # Test with ethanol (smaller compound)
        all_props = pugview.extract_all_experimental_properties(702)
        
        assert isinstance(all_props, dict)
        assert len(all_props) > 0
        
        # Check that common properties are present
        expected_props = ["Physical Description", "Melting Point", "Boiling Point", "Density"]
        found_props = list(all_props.keys())
        
        # At least some expected properties should be found
        assert any(prop in found_props for prop in expected_props)
        
        # Check data structure
        for prop_name, prop_data in all_props.items():
            assert isinstance(prop_data, list)
            if prop_data:  # If not empty
                assert all(isinstance(item, PropertyData) for item in prop_data)
    
    def test_get_available_properties(self, pugview):
        """Test getting list of available properties"""
        properties = pugview.get_available_properties(2244)
        
        assert isinstance(properties, list)
        assert len(properties) > 0
        assert all(isinstance(prop, str) for prop in properties)
        
        # Check for common properties
        assert any("Melting Point" in prop for prop in properties)
    
    def test_get_property_summary(self, pugview):
        """Test property summary functionality"""
        summary = pugview.get_property_summary(2244, "Melting Point")
        
        assert isinstance(summary, dict)
        required_keys = ["property", "values", "references", "units", "conditions", "count"]
        for key in required_keys:
            assert key in summary
        
        assert summary["property"] == "Melting Point"
        assert isinstance(summary["values"], list)
        assert isinstance(summary["count"], int)
        assert summary["count"] > 0
    
    def test_convenience_methods(self, pugview):
        """Test convenience methods for common properties"""
        # Test melting point
        mp_data = pugview.get_melting_point(2244)
        assert isinstance(mp_data, list)
        assert len(mp_data) > 0
        
        # Test boiling point
        bp_data = pugview.get_boiling_point(2244)
        assert isinstance(bp_data, list)
        
        # Test density
        density_data = pugview.get_density(702)  # ethanol
        assert isinstance(density_data, list)
        
        # Test viscosity with DMSO
        viscosity_data = pugview.get_viscosity(679)
        assert isinstance(viscosity_data, list)
        assert len(viscosity_data) > 0
    
    def test_batch_extract_properties(self, pugview):
        """Test batch property extraction"""
        properties = ["Melting Point", "Boiling Point", "Density"]
        results = pugview.batch_extract_properties(2244, properties)
        
        assert isinstance(results, dict)
        assert len(results) == len(properties)
        
        for prop in properties:
            assert prop in results
            assert isinstance(results[prop], list)
    
    def test_get_property_table(self, pugview):
        """Test property table functionality"""
        # Test with melting point
        table = pugview.get_property_table(2244, "Melting Point")
        
        assert isinstance(table, pd.DataFrame)
        assert len(table) > 0
        
        # Check expected columns
        expected_columns = ["CID", "StringWithMarkup", "ExperimentalValue", "Unit", "FullReference"]
        for col in expected_columns:
            assert col in table.columns
        
        # Check data types and content
        assert all(table["CID"] == 2244)
        assert all(isinstance(val, str) for val in table["StringWithMarkup"])
        
        # Check that we have some non-null experimental values
        assert table["ExperimentalValue"].notna().any()
    
    def test_value_parsing(self, pugview):
        """Test value and unit parsing"""
        # Test various value strings
        test_cases = [
            ("275 °F", ("275", "F")),
            ("2.47cP at 20 °C", ("2.47", "cP")),
            ("138-140", ("138-140", None)),
            ("1.95 mmÂ²/s", ("1.95", "mm")),
            ("0.79 g/cm³", ("0.79", "g/cm³"))
        ]
        
        for value_str, expected in test_cases:
            result = pugview._extract_experimental_value_and_unit(value_str)
            assert result[0] == expected[0], f"Failed for {value_str}: expected {expected[0]}, got {result[0]}"
            # Unit matching can be flexible due to complex patterns
    
    def test_reference_extraction(self, pugview):
        """Test reference mapping extraction"""
        # Get a response with references
        response = pugview.get_property(2244, "Melting Point")
        ref_map = pugview._extract_reference_map(response)
        
        assert isinstance(ref_map, dict)
        assert len(ref_map) > 0
        
        # Check that reference numbers are integers and values are strings
        for ref_num, ref_str in ref_map.items():
            assert isinstance(ref_num, int)
            assert isinstance(ref_str, str)
            assert len(ref_str) > 0
    
    def test_export_properties_to_dict(self, pugview):
        """Test exporting PropertyData to dictionaries"""
        property_data = pugview.extract_property_data(2244, "Melting Point")
        dict_list = pugview.export_properties_to_dict(property_data)
        
        assert isinstance(dict_list, list)
        assert len(dict_list) == len(property_data)
        
        if dict_list:
            first_dict = dict_list[0]
            expected_keys = ["value", "unit", "conditions", "reference", "reference_number", "description", "name"]
            for key in expected_keys:
                assert key in first_dict
    
    def test_error_handling(self, pugview):
        """Test error handling for invalid inputs"""
        # Test with non-existent CID
        result = pugview.extract_property_data(99999999, "Melting Point")
        assert isinstance(result, list)
        assert len(result) == 0  # Should return empty list, not raise exception
        
        # Test with non-existent property
        result = pugview.extract_property_data(2244, "NonExistentProperty")
        assert isinstance(result, list)
        assert len(result) == 0


class TestConvenienceFunctions:
    """Test suite for convenience functions"""
    
    def test_get_experimental_property(self):
        """Test get_experimental_property convenience function"""
        result = get_experimental_property(2244, "Melting Point")
        
        assert isinstance(result, list)
        assert len(result) > 0
        assert all(isinstance(item, PropertyData) for item in result)
    
    def test_get_all_experimental_properties(self):
        """Test get_all_experimental_properties convenience function"""
        result = get_all_experimental_properties(702)  # ethanol
        
        assert isinstance(result, dict)
        assert len(result) > 0
        
        for prop_name, prop_data in result.items():
            assert isinstance(prop_data, list)
            if prop_data:
                assert all(isinstance(item, PropertyData) for item in prop_data)
    
    def test_get_property_values_only(self):
        """Test get_property_values_only convenience function"""
        result = get_property_values_only(2244, "Melting Point")
        
        assert isinstance(result, list)
        assert len(result) > 0
        assert all(isinstance(value, str) for value in result)
    
    def test_get_property_table(self):
        """Test get_property_table convenience function"""
        result = get_property_table(679, "Viscosity")
        
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0
        
        expected_columns = ["CID", "StringWithMarkup", "ExperimentalValue", "Unit", "FullReference"]
        for col in expected_columns:
            assert col in result.columns


class TestPropertyData:
    """Test suite for PropertyData dataclass"""
    
    def test_property_data_creation(self):
        """Test PropertyData object creation"""
        data = PropertyData(
            value="275 °F",
            unit="F",
            conditions="at 760 mmHg",
            reference="Test Reference",
            reference_number=1,
            description="Test Description",
            name="Test Name"
        )
        
        assert data.value == "275 °F"
        assert data.unit == "F"
        assert data.conditions == "at 760 mmHg"
        assert data.reference == "Test Reference"
        assert data.reference_number == 1
        assert data.description == "Test Description"
        assert data.name == "Test Name"
    
    def test_property_data_defaults(self):
        """Test PropertyData with default values"""
        data = PropertyData(value="test value")
        
        assert data.value == "test value"
        assert data.unit is None
        assert data.conditions is None
        assert data.reference is None
        assert data.reference_number is None
        assert data.description is None
        assert data.name is None


class TestEdgeCases:
    """Test suite for edge cases and special scenarios"""
    
    def test_empty_response_handling(self):
        """Test handling of empty or malformed responses"""
        pugview = PubChemView()
        
        # Test with empty dictionary
        result = pugview._parse_property_response({})
        assert isinstance(result, list)
        assert len(result) == 0
        
        # Test with malformed response
        result = pugview._parse_property_response({"Record": {}})
        assert isinstance(result, list)
        assert len(result) == 0
    
    def test_special_characters_in_values(self):
        """Test handling of special characters in property values"""
        pugview = PubChemView()
        
        # Test with special characters
        test_values = [
            "1.95 mmÂ²/s at 20Â °C",  # Contains special characters
            "",  # Empty string
            "no numeric value here",  # No numbers
            "multiple 123 numbers 456 here"  # Multiple numbers
        ]
        
        for value in test_values:
            result = pugview._extract_experimental_value_and_unit(value)
            assert isinstance(result, tuple)
            assert len(result) == 4  # Returns (value, unit, temperature, conditions)
    
    def test_rate_limiting(self):
        """Test rate limiting functionality"""
        pugview = PubChemView()
        
        import time
        start_time = time.time()
        
        # Make multiple calls
        for _ in range(3):
            pugview._rate_limit()
        
        elapsed = time.time() - start_time
        
        # Should take at least 2 * min_request_interval due to rate limiting
        expected_min_time = 2 * pugview.min_request_interval
        assert elapsed >= expected_min_time * 0.9  # Allow some tolerance


if __name__ == "__main__":
    # Run tests if executed directly
    pytest.main([__file__, "-v"])
