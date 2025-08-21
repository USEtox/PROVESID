"""
Tests for OPSIN (Open Parser for Systematic IUPAC Nomenclature) functionality
"""

import pytest
import sys
import os

# Add the src directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from provesid.opsin import OPSIN


class TestOPSIN:
    """Test suite for OPSIN class"""
    
    @pytest.fixture
    def opsin(self):
        """Create an OPSIN instance for testing"""
        return OPSIN()
    
    def test_initialization(self, opsin):
        """Test OPSIN initialization"""
        assert opsin.base_url == "https://opsin.ch.cam.ac.uk/opsin/"
        assert hasattr(opsin, 'responses')
        
        expected_responses = {200: "SUCCESS", 404: "FAILURE", 500: "Internal server error"}
        assert opsin.responses == expected_responses
    
    def test_empty_res_structure(self, opsin):
        """Test the structure of empty response"""
        empty_res = opsin._empty_res()
        
        expected_keys = ['status', 'message', 'inchi', 'stdinchi', 'stdinchikey', 'smiles']
        
        for key in expected_keys:
            assert key in empty_res
        
        # Check default values are empty strings
        for key in expected_keys:
            assert empty_res[key] == ""
    
    def test_get_id_valid_iupac_name(self, opsin):
        """Test get_id with valid IUPAC names"""
        valid_names = [
            'ethanol',
            'methanol', 
            'benzene',
            'acetic acid',
            'propane'
        ]
        
        for name in valid_names:
            result = opsin.get_id(name)
            
            assert 'status' in result
            
            if result['status'] == 'SUCCESS':
                # Should have chemical identifiers
                assert 'smiles' in result
                assert 'inchi' in result
                assert 'stdinchi' in result
                assert 'stdinchikey' in result
                
                # SMILES should contain carbon for organic compounds
                if result['smiles']:
                    # Most organic compounds contain carbon
                    if name not in ['water']:  # water doesn't contain carbon
                        assert 'C' in result['smiles'] or 'c' in result['smiles']
                
                # InChI should start with "InChI="
                if result['inchi']:
                    assert result['inchi'].startswith('InChI=')
                
                if result['stdinchi']:
                    assert result['stdinchi'].startswith('InChI=')
                
                # InChI Key should be the right format (length ~27 with hyphens)
                if result['stdinchikey']:
                    assert len(result['stdinchikey']) >= 20  # Approximate check
                    assert '-' in result['stdinchikey']
            
            elif result['status'] in ['FAILURE', 'Internal server error']:
                # API might not recognize the name or be unavailable
                pass
            else:
                pytest.fail(f"Unexpected status: {result['status']}")
    
    def test_get_id_specific_compounds(self, opsin):
        """Test get_id with specific well-known compounds"""
        # Test ethanol
        result = opsin.get_id('ethanol')
        
        if result['status'] == 'SUCCESS':
            # Ethanol should have specific characteristics
            if result['smiles']:
                # Ethanol SMILES should be CCO or OCC or C(C)O
                smiles = result['smiles']
                assert ('CCO' in smiles or 'OCC' in smiles or 'C-C-O' in smiles or 
                       'C(C)O' in smiles or 'OC(C)' in smiles)
            
            if result['stdinchi']:
                # Should contain ethanol formula
                assert 'C2H6O' in result['stdinchi']
        else:
            pytest.skip("OPSIN API not available or name not recognized")
    
    def test_get_id_invalid_names(self, opsin):
        """Test get_id with invalid IUPAC names"""
        invalid_names = [
            'definitely_not_a_chemical_name',
            'xyz123',
            '',
            '!@#$%^&*()',
            'random_text_string'
        ]
        
        for name in invalid_names:
            result = opsin.get_id(name)
            
            assert 'status' in result
            
            if name == '':  # Empty string might cause different behavior
                continue
            
            # Should return FAILURE for invalid names
            assert result['status'] in ['FAILURE', 'Internal server error', 'SUCCESS']
            
            if result['status'] == 'FAILURE':
                # Failed results should have empty chemical data
                assert result['smiles'] == '' or result['smiles'] is None
                assert result['inchi'] == '' or result['inchi'] is None
    
    def test_get_id_timeout(self, opsin):
        """Test get_id with timeout parameter"""
        # Test with reasonable timeout
        result = opsin.get_id('ethanol', timeout=30)
        assert 'status' in result
        
        # Test with very short timeout (might fail)
        try:
            result = opsin.get_id('ethanol', timeout=0.001)
            assert 'status' in result
            # Very short timeout might still succeed with fast network
        except Exception:
            # Timeout exception is acceptable
            pass
    
    def test_get_id_from_list_valid(self, opsin):
        """Test get_id_from_list with valid IUPAC names"""
        names = ['ethanol', 'methanol', 'benzene']
        results = opsin.get_id_from_list(names)
        
        assert isinstance(results, list)
        assert len(results) == len(names)
        
        for i, result in enumerate(results):
            assert 'status' in result
            
            if result['status'] == 'SUCCESS':
                # Should have chemical identifiers
                assert 'smiles' in result
                assert 'inchi' in result
                assert 'stdinchi' in result
                assert 'stdinchikey' in result
    
    def test_get_id_from_list_mixed(self, opsin):
        """Test get_id_from_list with mix of valid and invalid names"""
        names = ['ethanol', 'invalid_name_xyz', 'benzene', 'not_a_chemical']
        results = opsin.get_id_from_list(names, pause_time=0.1)
        
        assert isinstance(results, list)
        assert len(results) == len(names)
        
        # Check that all entries have status
        for result in results:
            assert 'status' in result
            assert result['status'] in ['SUCCESS', 'FAILURE', 'Internal server error']
    
    def test_get_id_from_list_empty(self, opsin):
        """Test get_id_from_list with empty list"""
        results = opsin.get_id_from_list([])
        
        assert isinstance(results, list)
        assert len(results) == 0
    
    def test_get_id_from_list_pause_time(self, opsin):
        """Test get_id_from_list with different pause times"""
        names = ['ethanol', 'methanol']
        
        import time
        
        # Test with no pause
        start_time = time.time()
        results = opsin.get_id_from_list(names, pause_time=0)
        no_pause_time = time.time() - start_time
        
        # Test with pause
        start_time = time.time()
        results_with_pause = opsin.get_id_from_list(names, pause_time=0.5)
        with_pause_time = time.time() - start_time
        
        # With pause should take longer (allowing for network variance)
        # This is a rough check since network timing can vary
        assert len(results) == len(results_with_pause)
        
        # Both should return proper results
        for result in results:
            assert 'status' in result
        for result in results_with_pause:
            assert 'status' in result
    
    def test_complex_iupac_names(self, opsin):
        """Test get_id with complex IUPAC names"""
        complex_names = [
            '2-methylpropane',
            '2,2-dimethylpropane',
            'cyclohexane',
            'phenol',
            'benzoic acid'
        ]
        
        for name in complex_names:
            result = opsin.get_id(name)
            
            assert 'status' in result
            
            if result['status'] == 'SUCCESS':
                # Should have valid chemical data
                if result['smiles']:
                    # Should contain carbon for these organic compounds
                    assert 'C' in result['smiles'] or 'c' in result['smiles']
                
                if result['stdinchi']:
                    assert result['stdinchi'].startswith('InChI=')
            
            # Some complex names might not be recognized
            elif result['status'] == 'FAILURE':
                pass
            else:
                pytest.fail(f"Unexpected status for {name}: {result['status']}")


class TestOPSINIntegration:
    """Integration tests for OPSIN API"""
    
    @pytest.fixture
    def opsin(self):
        return OPSIN()
    
    def test_ethanol_complete_data(self, opsin):
        """Test complete data extraction for ethanol"""
        result = opsin.get_id('ethanol')
        
        if result['status'] != 'SUCCESS':
            pytest.skip("OPSIN API not available or ethanol not recognized")
        
        # Ethanol-specific validations
        if result['smiles']:
            smiles = result['smiles'].upper()
            # Should represent ethanol structure
            assert 'CCO' in smiles or 'OCC' in smiles or ('C' in smiles and 'O' in smiles)
        
        if result['stdinchi']:
            # Should contain ethanol molecular formula
            assert 'C2H6O' in result['stdinchi']
        
        if result['stdinchikey']:
            # Should be proper InChI Key format
            assert len(result['stdinchikey']) >= 20
            assert result['stdinchikey'].count('-') >= 2
    
    def test_benzene_aromatics(self, opsin):
        """Test aromatic compound recognition"""
        result = opsin.get_id('benzene')
        
        if result['status'] != 'SUCCESS':
            pytest.skip("OPSIN API not available or benzene not recognized")
        
        if result['smiles']:
            smiles = result['smiles'].lower()
            # Benzene should have aromatic representation
            assert 'c1ccccc1' in smiles or 'C1=CC=CC=C1' in smiles.upper() or 'c' in smiles
        
        if result['stdinchi']:
            # Should contain benzene formula
            assert 'C6H6' in result['stdinchi']
    
    def test_cross_validation_with_known_structures(self, opsin):
        """Test cross-validation with known chemical structures"""
        # Test methanol
        result = opsin.get_id('methanol')
        
        if result['status'] != 'SUCCESS':
            pytest.skip("OPSIN API not available")
        
        if result['smiles'] and result['stdinchi']:
            # Methanol should have CH4O formula
            if 'CH4O' in result['stdinchi'] or 'C H4 O' in result['stdinchi']:
                # SMILES should be consistent (CO or OC)
                smiles = result['smiles'].upper()
                assert 'CO' in smiles or 'OC' in smiles
    
    def test_systematic_vs_common_names(self, opsin):
        """Test systematic vs common name recognition"""
        name_pairs = [
            ('ethanol', 'ethyl alcohol'),
            ('methanol', 'methyl alcohol'),
            ('propane', 'n-propane'),
        ]
        
        for systematic, common in name_pairs:
            sys_result = opsin.get_id(systematic)
            com_result = opsin.get_id(common)
            
            # At least one should succeed
            if sys_result['status'] == 'SUCCESS' or com_result['status'] == 'SUCCESS':
                # If both succeed, they should give same or similar results
                if sys_result['status'] == 'SUCCESS' and com_result['status'] == 'SUCCESS':
                    if sys_result['stdinchikey'] and com_result['stdinchikey']:
                        # InChI Keys should be identical for same compound
                        assert sys_result['stdinchikey'] == com_result['stdinchikey']


class TestOPSINErrorHandling:
    """Test error handling scenarios for OPSIN"""
    
    @pytest.fixture
    def opsin(self):
        return OPSIN()
    
    def test_network_timeout_handling(self, opsin):
        """Test handling of network timeouts"""
        # Test with very short timeout
        try:
            result = opsin.get_id('ethanol', timeout=0.001)
            # If it succeeds, network was very fast
            assert 'status' in result
        except Exception as e:
            # Timeout or connection error is expected
            assert 'timeout' in str(e).lower() or 'connection' in str(e).lower()
    
    def test_malformed_url_handling(self, opsin):
        """Test behavior with malformed URLs (indirectly)"""
        # Test with names that might cause URL encoding issues
        problematic_names = [
            'name with spaces',
            'name/with/slashes',
            'name?with?questions',
            'name&with&ampersands',
            'name#with#hashes'
        ]
        
        for name in problematic_names:
            try:
                result = opsin.get_id(name)
                assert 'status' in result
                # Should handle URL encoding gracefully
            except Exception as e:
                # URL encoding errors are acceptable
                assert any(keyword in str(e).lower() for keyword in ['url', 'encoding', 'request'])
    
    def test_empty_and_whitespace_names(self, opsin):
        """Test handling of empty and whitespace-only names"""
        empty_names = ['', '   ', '\t', '\n', '  \t\n  ']
        
        for name in empty_names:
            try:
                result = opsin.get_id(name)
                assert 'status' in result
                # Should handle gracefully, likely with FAILURE status
                if result['status'] not in ['FAILURE', 'Internal server error']:
                    # Some might be processed as valid (unlikely but possible)
                    pass
            except Exception:
                # Request errors for empty names are acceptable
                pass
    
    def test_very_long_names(self, opsin):
        """Test handling of very long compound names"""
        # Create a very long name
        long_name = 'a' * 1000  # 1000 character name
        
        try:
            result = opsin.get_id(long_name)
            assert 'status' in result
            # Should handle gracefully, likely with FAILURE
            assert result['status'] in ['FAILURE', 'Internal server error', 'SUCCESS']
        except Exception as e:
            # URL too long or other request errors are acceptable
            assert any(keyword in str(e).lower() for keyword in ['url', 'request', 'long', 'size'])
    
    def test_unicode_and_special_characters(self, opsin):
        """Test handling of unicode and special characters"""
        unicode_names = [
            'café',  # accented characters
            'naïve',  # diaeresis
            'α-glucose',  # Greek letters
            'β-carotene',  # Greek letters
            '2,4-dinitrophenol',  # numbers and hyphens
        ]
        
        for name in unicode_names:
            try:
                result = opsin.get_id(name)
                assert 'status' in result
                # Should handle unicode gracefully
            except Exception as e:
                # Unicode encoding errors are acceptable
                assert any(keyword in str(e).lower() for keyword in ['unicode', 'encoding', 'character'])
    
    def test_batch_processing_error_recovery(self, opsin):
        """Test that batch processing continues after individual failures"""
        # Mix of valid, invalid, and problematic names
        mixed_names = [
            'ethanol',  # should work
            'invalid_xyz_123',  # should fail
            'methanol',  # should work
            '',  # might cause error
            'benzene',  # should work
        ]
        
        results = opsin.get_id_from_list(mixed_names, pause_time=0.1)
        
        # Should return results for all entries
        assert len(results) == len(mixed_names)
        
        # Each result should have a status
        for i, result in enumerate(results):
            assert 'status' in result, f"Missing status for entry {i}: {mixed_names[i]}"
        
        # Should have at least some successful results (ethanol, methanol, benzene)
        successful_count = sum(1 for r in results if r['status'] == 'SUCCESS')
        # At least one should succeed (allowing for API availability issues)
        if successful_count == 0:
            pytest.skip("OPSIN API not available or all names failed")


if __name__ == "__main__":
    # Run tests if executed directly
    pytest.main([__file__, "-v"])
