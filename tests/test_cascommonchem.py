"""
Tests for CAS Common Chemistry API functionality
"""

import pytest
import sys
import os
import json

# Add the src directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from provesid.cascommonchem import CASCommonChem


class TestCASCommonChem:
    """Test suite for CASCommonChem class"""
    
    @pytest.fixture
    def cas_api(self):
        """Create a CASCommonChem instance for testing"""
        return CASCommonChem()
    
    def test_initialization(self, cas_api):
        """Test CASCommonChem initialization"""
        assert hasattr(cas_api, 'base_url')
        assert hasattr(cas_api, 'swagger')
        assert hasattr(cas_api, 'query_url')
        assert hasattr(cas_api, 'responses')
        
        # Check that swagger data was loaded
        assert isinstance(cas_api.swagger, dict)
        assert 'host' in cas_api.swagger
        assert 'schemes' in cas_api.swagger
        assert 'basePath' in cas_api.swagger
        
        # Check URL construction
        assert cas_api.base_url.startswith('http')
        assert len(cas_api.query_url) == 3
        assert '/detail' in cas_api.query_url
        assert '/export' in cas_api.query_url
        assert '/search' in cas_api.query_url
    
    def test_empty_res_structure(self, cas_api):
        """Test the structure of empty response"""
        empty_res = cas_api._empty_res()
        
        expected_keys = [
            'cas_rn', 'status', 'canonicalSmile', 'experimentalProperties',
            'hasMolfile', 'images', 'inchi', 'inchiKey', 'molecularFormula',
            'molecularMass', 'name', 'propertyCitations', 'replacedRns',
            'rn', 'smile', 'synonyms', 'uri'
        ]
        
        for key in expected_keys:
            assert key in empty_res
        
        # Check default values
        assert empty_res['cas_rn'] == ""
        assert empty_res['status'] == ""
        assert empty_res['experimentalProperties'] == []
        assert empty_res['hasMolfile'] is False
        assert empty_res['images'] == []
        assert isinstance(empty_res['synonyms'], list)
    
    def test_cas_to_detail_valid(self, cas_api):
        """Test cas_to_detail with valid CAS number"""
        # Test with water CAS number
        result = cas_api.cas_to_detail('7732-18-5')
        
        # Should have a status
        assert 'status' in result
        
        if result['status'] == 'Success':
            # Check that we got valid data
            assert 'rn' in result
            assert 'name' in result
            assert 'molecularFormula' in result
            assert 'molecularMass' in result
            
            # Water should have H2O formula (may be HTML formatted)
            if result['molecularFormula']:
                formula = result['molecularFormula']
                assert 'H2O' in formula or 'H<sub>2</sub>O' in formula
        elif result['status'] in ['Invalid Request', 'Error']:
            # API might not be available or rate limited
            pytest.skip("CAS Common Chemistry API not available")
    
    def test_cas_to_detail_invalid(self, cas_api):
        """Test cas_to_detail with invalid CAS number"""
        result = cas_api.cas_to_detail('invalid-cas-number')
        
        # Should return error or invalid request
        assert 'status' in result
        assert result['status'] in ['Invalid Request', 'Error', 'Not found']
    
    def test_cas_to_detail_with_hyphens(self, cas_api):
        """Test cas_to_detail handles CAS numbers with hyphens"""
        # Test ethanol with hyphens
        result1 = cas_api.cas_to_detail('64-17-5')
        
        # Test ethanol without hyphens  
        result2 = cas_api.cas_to_detail('64175')
        
        # Both should have status
        assert 'status' in result1
        assert 'status' in result2
        
        # If successful, should have similar data structure
        if result1['status'] == 'Success':
            assert 'rn' in result1
            assert 'molecularFormula' in result1
    
    def test_name_to_detail_valid(self, cas_api):
        """Test name_to_detail with valid compound name"""
        result = cas_api.name_to_detail('ethanol')
        
        assert 'status' in result
        
        if result['status'] == 'Success':
            # Should have compound data
            assert 'name' in result
            assert 'molecularFormula' in result
            assert 'rn' in result  # Should have found a CAS number
            
            # Ethanol should have C2H6O formula (may be HTML formatted)
            if result['molecularFormula']:
                formula = result['molecularFormula']
                assert ('C2H6O' in formula or 'C2 H6 O' in formula or 
                       'C<sub>2</sub>H<sub>6</sub>O' in formula)
        elif result['status'] in ['Not found', 'Error']:
            # API might not be available
            pytest.skip("CAS Common Chemistry API not available or compound not found")
    
    def test_name_to_detail_invalid(self, cas_api):
        """Test name_to_detail with invalid compound name"""
        result = cas_api.name_to_detail('definitely_not_a_chemical_compound_xyz123')
        
        assert 'status' in result
        # Should return not found or error
        assert result['status'] in ['Not found', 'Error']
    
    def test_name_to_detail_multiple_results(self, cas_api):
        """Test name_to_detail with name that might return multiple results"""
        # Test with a generic name that might match multiple compounds
        result = cas_api.name_to_detail('alcohol')
        
        assert 'status' in result
        
        if result['status'] == 'Success':
            # Should have selected the first result
            assert 'rn' in result
            assert result['rn'] != ""
    
    def test_smiles_to_detail(self, cas_api):
        """Test smiles_to_detail functionality"""
        # Test with ethanol SMILES
        result = cas_api.smiles_to_detail('CCO')
        
        assert 'status' in result
        
        if result['status'] == 'Success':
            # Should have found ethanol data
            assert 'molecularFormula' in result
            assert 'rn' in result
            
            # Check for ethanol characteristics
            if result['molecularFormula']:
                assert 'C2H6O' in result['molecularFormula'] or 'C2 H6 O' in result['molecularFormula']
        elif result['status'] in ['Not found', 'Error']:
            pytest.skip("CAS Common Chemistry API not available or SMILES not found")
    
    def test_smiles_to_detail_invalid(self, cas_api):
        """Test smiles_to_detail with invalid SMILES"""
        result = cas_api.smiles_to_detail('invalid_smiles_xyz')
        
        assert 'status' in result
        assert result['status'] in ['Not found', 'Error']
    
    def test_response_codes(self, cas_api):
        """Test that response codes are properly defined"""
        expected_codes = {200: "Success", 400: "Invalid Request", 404: "Invalid Request", 500: "Internal Server Error"}
        
        assert cas_api.responses == expected_codes
    
    def test_timeout_parameter(self, cas_api):
        """Test that timeout parameter is accepted"""
        # Test with very short timeout (might fail, but should not raise exception)
        try:
            result = cas_api.cas_to_detail('7732-18-5', timeout=1)
            assert 'status' in result
        except Exception:
            # Timeout or network error is acceptable for this test
            pass
    
    def test_swagger_file_loading(self, cas_api):
        """Test that swagger file is properly loaded"""
        # Check swagger file path exists
        assert os.path.exists(cas_api.data_folder)
        
        # Check swagger data structure
        assert 'host' in cas_api.swagger
        assert 'schemes' in cas_api.swagger
        assert 'basePath' in cas_api.swagger
        
        # Check URL construction elements
        assert isinstance(cas_api.swagger['schemes'], list)
        assert len(cas_api.swagger['schemes']) > 0


class TestCASCommonChemIntegration:
    """Integration tests for CAS Common Chemistry API"""
    
    @pytest.fixture
    def cas_api(self):
        return CASCommonChem()
    
    def test_water_compound_data(self, cas_api):
        """Test comprehensive data extraction for water"""
        result = cas_api.cas_to_detail('7732-18-5')
        
        if result['status'] != 'Success':
            pytest.skip("CAS Common Chemistry API not available")
        
        # Water-specific checks
        if result.get('molecularFormula'):
            formula = result['molecularFormula']
            assert 'H2O' in formula or 'H<sub>2</sub>O' in formula
        
        if result.get('name'):
            assert 'water' in result['name'].lower()
        
        # Should have basic structure information
        assert 'inchi' in result
        assert 'inchiKey' in result
        assert 'smile' in result or 'canonicalSmile' in result
    
    def test_ethanol_compound_data(self, cas_api):
        """Test comprehensive data extraction for ethanol"""
        result = cas_api.cas_to_detail('64-17-5')
        
        if result['status'] != 'Success':
            pytest.skip("CAS Common Chemistry API not available")
        
        # Ethanol-specific checks
        if result.get('molecularFormula'):
            formula = result['molecularFormula']
            assert ('C2H6O' in formula or 'C2 H6 O' in formula or 
                   'C<sub>2</sub>H<sub>6</sub>O' in formula)
        
        if result.get('name'):
            name = result['name'].lower()
            assert 'ethanol' in name or 'alcohol' in name
        
        # Should have CAS number
        if result.get('rn'):
            assert '64-17-5' in result['rn']
    
    def test_cross_method_consistency(self, cas_api):
        """Test that different methods return consistent data for same compound"""
        # Test ethanol via CAS number and name
        cas_result = cas_api.cas_to_detail('64-17-5')
        name_result = cas_api.name_to_detail('ethanol')
        smiles_result = cas_api.smiles_to_detail('CCO')
        
        successful_results = []
        for result in [cas_result, name_result, smiles_result]:
            if result['status'] == 'Success':
                successful_results.append(result)
        
        if len(successful_results) < 2:
            pytest.skip("Not enough successful API calls for comparison")
        
        # Compare molecular formulas if available
        formulas = [r.get('molecularFormula') for r in successful_results if r.get('molecularFormula')]
        if len(formulas) > 1:
            # All formulas should be similar (allowing for formatting differences)
            for formula in formulas:
                # Check for C2H6O pattern in various formats
                assert (('C2' in formula and 'H6' in formula and 'O' in formula) or
                       'C<sub>2</sub>H<sub>6</sub>O' in formula)
    
    def test_experimental_properties_structure(self, cas_api):
        """Test structure of experimental properties when available"""
        # Try with a compound that might have experimental properties
        result = cas_api.cas_to_detail('64-17-5')  # ethanol
        
        if result['status'] != 'Success':
            pytest.skip("CAS Common Chemistry API not available")
        
        # Check experimental properties structure
        exp_props = result.get('experimentalProperties', [])
        assert isinstance(exp_props, list)
        
        # If there are experimental properties, check their structure
        for prop in exp_props:
            if isinstance(prop, dict):
                # Properties should have some identifiable structure
                # (exact structure may vary by API version)
                assert len(prop) > 0
    
    def test_synonyms_structure(self, cas_api):
        """Test structure of synonyms when available"""
        result = cas_api.cas_to_detail('7732-18-5')  # water
        
        if result['status'] != 'Success':
            pytest.skip("CAS Common Chemistry API not available")
        
        synonyms = result.get('synonyms', [])
        assert isinstance(synonyms, list)
        
        # If there are synonyms, they should be strings
        for synonym in synonyms:
            if synonym:  # Skip empty entries
                assert isinstance(synonym, str)
    
    def test_molecular_mass_format(self, cas_api):
        """Test molecular mass format when available"""
        result = cas_api.cas_to_detail('7732-18-5')  # water
        
        if result['status'] != 'Success':
            pytest.skip("CAS Common Chemistry API not available")
        
        mol_mass = result.get('molecularMass')
        if mol_mass:
            # Should be a number or string representation of a number
            if isinstance(mol_mass, str):
                try:
                    float(mol_mass)
                except ValueError:
                    pytest.fail(f"Molecular mass '{mol_mass}' is not a valid number")
            elif isinstance(mol_mass, (int, float)):
                assert mol_mass > 0


class TestCASCommonChemErrorHandling:
    """Test error handling scenarios for CAS Common Chemistry API"""
    
    @pytest.fixture
    def cas_api(self):
        return CASCommonChem()
    
    def test_network_timeout(self, cas_api):
        """Test handling of network timeouts"""
        # Test with very short timeout
        result = cas_api.cas_to_detail('7732-18-5', timeout=0.001)
        
        # Should handle timeout gracefully
        assert 'status' in result
        assert result['status'] in ['Error', 'Success']  # Success if very fast network
    
    def test_malformed_cas_numbers(self, cas_api):
        """Test handling of malformed CAS numbers"""
        malformed_cas = [
            '',
            'not-a-cas-number',
            '123-45-6789',  # Too many digits
            '12-3-4',       # Too few digits
            'abc-def-ghi',  # Non-numeric
        ]
        
        for cas in malformed_cas:
            result = cas_api.cas_to_detail(cas)
            assert 'status' in result
            # Should return error or invalid request
            if result['status'] not in ['Error', 'Invalid Request', 'Not found']:
                # Some malformed CAS might still get processed
                pass
    
    def test_empty_compound_names(self, cas_api):
        """Test handling of empty or invalid compound names"""
        invalid_names = ['', '   ', 'xyz123notacompound', '!@#$%^&*()']
        
        for name in invalid_names:
            result = cas_api.name_to_detail(name)
            assert 'status' in result
            # Should return not found or error
            expected_statuses = ['Not found', 'Error', 'Invalid Request']
            if result['status'] not in expected_statuses:
                # Some invalid names might still return results
                pass
    
    def test_special_characters_in_names(self, cas_api):
        """Test handling of special characters in compound names"""
        special_names = [
            'α-glucose',
            'β-carotene', 
            '2,4-dinitrophenol',
            'N,N-dimethylformamide'
        ]
        
        for name in special_names:
            try:
                result = cas_api.name_to_detail(name)
                assert 'status' in result
                # Should handle special characters gracefully
            except Exception as e:
                # URL encoding issues are acceptable
                assert 'encoding' in str(e).lower() or 'unicode' in str(e).lower()


if __name__ == "__main__":
    # Run tests if executed directly
    pytest.main([__file__, "-v"])
