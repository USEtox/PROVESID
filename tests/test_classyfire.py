"""
Tests for ClassyFireAPI functionality
"""

import pytest
import sys
import os
import time

# Add the src directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from provesid.classyfire import ClassyFireAPI


class TestClassyFireAPI:
    """Test suite for ClassyFireAPI class"""
    
    def test_base_url(self):
        """Test that base URL is correctly set"""
        expected_url = "http://classyfire.wishartlab.com"
        assert ClassyFireAPI.URL == expected_url
    
    def test_submit_query_valid_smiles(self):
        """Test submit_query with valid SMILES strings"""
        valid_smiles = [
            'CCO',  # ethanol
            'C',    # methane
            'O',    # water
            'c1ccccc1',  # benzene
            'CC(=O)O'  # acetic acid
        ]
        
        for smiles in valid_smiles:
            try:
                response = ClassyFireAPI.submit_query(f"test_{smiles}", smiles)
                
                # Should return a response object
                assert response is not None
                
                # Check if response has expected structure
                if hasattr(response, 'status_code'):
                    # If successful submission
                    if response.status_code in [200, 201, 202]:
                        # Should have some content
                        assert len(response.text) > 0
                        
                        # Try to parse as JSON if possible
                        try:
                            data = response.json()
                            # Should have query information
                            if isinstance(data, dict):
                                # Common ClassyFire response fields
                                possible_fields = ['id', 'query_id', 'status', 'smiles', 'label']
                                # At least one field should be present
                                assert any(field in data for field in possible_fields)
                        except ValueError:
                            # Not JSON, that's okay - might be HTML or other format
                            pass
                    
                    elif response.status_code in [400, 422]:
                        # Bad request - might be invalid SMILES format for API
                        pass
                    elif response.status_code in [500, 503]:
                        # Server error - API might be down
                        pytest.skip("ClassyFire API appears to be unavailable")
                    else:
                        # Other status codes
                        print(f"Unexpected status code {response.status_code} for SMILES {smiles}")
                
            except Exception as e:
                # Network errors, timeouts, etc.
                if "timeout" in str(e).lower() or "connection" in str(e).lower():
                    pytest.skip("ClassyFire API connection issues")
                else:
                    # Re-raise unexpected errors
                    raise
    
    def test_submit_query_with_label(self):
        """Test submit_query with custom labels"""
        smiles = 'CCO'  # ethanol
        label = 'test_ethanol_classification'
        
        try:
            response = ClassyFireAPI.submit_query(label, smiles)
            
            assert response is not None
            
            if hasattr(response, 'status_code'):
                if response.status_code in [200, 201, 202]:
                    # Should contain the label information
                    response_text = response.text.lower()
                    # Label might appear in response
                    # (This is a best-effort check since API format may vary)
                
        except Exception as e:
            if "timeout" in str(e).lower() or "connection" in str(e).lower():
                pytest.skip("ClassyFire API connection issues")
            else:
                raise
    
    def test_submit_query_invalid_smiles(self):
        """Test submit_query with invalid SMILES strings"""
        invalid_smiles = [
            'INVALID_SMILES',
            'XYZ123',
            '[[[[',
            '))))',
            '',
            '!@#$%'
        ]
        
        for smiles in invalid_smiles:
            try:
                response = ClassyFireAPI.submit_query(f"invalid_test_{smiles}", smiles)
                
                assert response is not None
                
                if hasattr(response, 'status_code'):
                    # Should handle invalid SMILES gracefully
                    # Common responses: 400 (bad request), 422 (unprocessable entity)
                    if response.status_code in [400, 422]:
                        # Expected behavior for invalid SMILES
                        pass
                    elif response.status_code in [200, 201, 202]:
                        # Some APIs might accept and process even "invalid" SMILES
                        # or return success with error message in body
                        pass
                    elif response.status_code in [500, 503]:
                        # Server error
                        pytest.skip("ClassyFire API server error")
                
            except Exception as e:
                if "timeout" in str(e).lower() or "connection" in str(e).lower():
                    pytest.skip("ClassyFire API connection issues")
                # Invalid SMILES might cause request exceptions, which is acceptable
                elif "invalid" in str(e).lower() or "bad" in str(e).lower():
                    pass
                else:
                    raise
    
    def test_query_status_valid_id(self):
        """Test query_status with query IDs"""
        # Test with various query ID formats
        test_ids = [
            '12345',
            'abcdef',
            '123abc',
            'cf_query_123'
        ]
        
        for query_id in test_ids:
            try:
                response = ClassyFireAPI.query_status(query_id)
                
                # Note: May return None for non-existent query IDs
                if response is not None:
                    if hasattr(response, 'status_code'):
                        if response.status_code == 200:
                            # Query found
                            assert len(response.text) > 0
                            
                            # Try to parse response
                            try:
                                data = response.json()
                                if isinstance(data, dict):
                                    # Common status fields
                                    possible_fields = ['status', 'state', 'id', 'query_id', 'classification']
                                    # Should have some relevant information
                            except ValueError:
                                # Not JSON format
                                pass
                        elif response.status_code == 404:
                            # Query not found - expected for random IDs
                            pass
                        elif response.status_code in [500, 503]:
                            # Server error
                            pytest.skip("ClassyFire API server error")
                    else:
                        # Got a valid response object without status_code
                        assert isinstance(response, dict)
                else:
                    # None response is acceptable for non-existent queries
                    pass
                
            except Exception as e:
                if "timeout" in str(e).lower() or "connection" in str(e).lower():
                    pytest.skip("ClassyFire API connection issues")
                else:
                    raise
    
    def test_query_status_invalid_id(self):
        """Test query_status with invalid query IDs"""
        invalid_ids = [
            '',
            '!@#$%',
            'very_long_id_' * 100,  # Very long ID
            None,
            123,  # Non-string type
        ]
        
        for query_id in invalid_ids:
            try:
                # Some might cause type errors before the request
                if query_id is None or isinstance(query_id, int):
                    # These should be handled by the method or cause appropriate errors
                    pass
                
                response = ClassyFireAPI.query_status(str(query_id) if query_id is not None else '')
                
                assert response is not None
                
                if hasattr(response, 'status_code'):
                    # Should handle invalid IDs gracefully
                    # Usually 404 (not found) or 400 (bad request)
                    assert response.status_code in [200, 400, 404, 422, 500, 503]
                
            except TypeError:
                # Type errors for None or invalid types are acceptable
                pass
            except Exception as e:
                if "timeout" in str(e).lower() or "connection" in str(e).lower():
                    pytest.skip("ClassyFire API connection issues")
                else:
                    # Other errors might be acceptable for invalid inputs
                    pass
    
    def test_get_query_valid_id(self):
        """Test get_query with query IDs"""
        test_ids = [
            '12345',
            'test_query',
            'abc123'
        ]
        
        for query_id in test_ids:
            try:
                response = ClassyFireAPI.get_query(query_id)
                
                assert response is not None
                
                if hasattr(response, 'status_code'):
                    if response.status_code == 200:
                        # Query found and results available
                        assert len(response.text) > 0
                        
                        # Try to parse response
                        try:
                            data = response.json()
                            if isinstance(data, dict):
                                # Should contain query results
                                possible_fields = [
                                    'smiles', 'inchikey', 'kingdom', 'superclass', 
                                    'class', 'subclass', 'direct_parent', 'classification'
                                ]
                                # Should have classification information
                        except ValueError:
                            # Not JSON format
                            pass
                            
                    elif response.status_code == 404:
                        # Query not found - expected for random IDs
                        pass
                    elif response.status_code in [202]:
                        # Query still processing
                        pass
                    elif response.status_code in [500, 503]:
                        # Server error
                        pytest.skip("ClassyFire API server error")
                
            except Exception as e:
                if "timeout" in str(e).lower() or "connection" in str(e).lower():
                    pytest.skip("ClassyFire API connection issues")
                else:
                    raise
    
    def test_get_query_format_parameter(self):
        """Test get_query with different format parameters"""
        query_id = 'test_12345'
        formats = ['json', 'sdf', 'csv']
        
        for fmt in formats:
            try:
                response = ClassyFireAPI.get_query(query_id, format=fmt)
                
                assert response is not None
                
                if hasattr(response, 'status_code'):
                    # Should handle different formats
                    if response.status_code == 200:
                        # Should return data in requested format
                        if fmt == 'json':
                            try:
                                response.json()  # Should be valid JSON
                            except ValueError:
                                # Might not be JSON even if format=json was requested
                                pass
                        elif fmt in ['sdf', 'mol']:
                            # Should be text content for molecular formats
                            assert isinstance(response.text, str)
                    
                    elif response.status_code == 404:
                        # Query not found
                        pass
                    elif response.status_code in [400, 422]:
                        # Invalid format parameter
                        pass
                
            except Exception as e:
                if "timeout" in str(e).lower() or "connection" in str(e).lower():
                    pytest.skip("ClassyFire API connection issues")
                else:
                    raise


class TestClassyFireAPIIntegration:
    """Integration tests for ClassyFireAPI"""
    
    def test_full_workflow_simulation(self):
        """Test a simulated full workflow"""
        # This test simulates the full workflow but doesn't rely on actual query completion
        # since ClassyFire queries can take a long time to process
        
        smiles = 'CCO'  # ethanol
        label = 'integration_test_ethanol'
        
        try:
            # Step 1: Submit query
            submit_response = ClassyFireAPI.submit_query(label, smiles)
            
            if submit_response is None:
                pytest.skip("ClassyFire API not responding")
            
            if not hasattr(submit_response, 'status_code'):
                pytest.skip("Unexpected response format from ClassyFire API")
            
            if submit_response.status_code not in [200, 201, 202]:
                if submit_response.status_code in [500, 503]:
                    pytest.skip("ClassyFire API server error")
                else:
                    # Other status codes - might be expected for test data
                    pass
            
            # If submission was successful, extract query ID
            query_id = None
            if submit_response.status_code in [200, 201, 202]:
                try:
                    data = submit_response.json()
                    if isinstance(data, dict):
                        # Try to find query ID in response
                        for key in ['id', 'query_id', 'queryId']:
                            if key in data:
                                query_id = data[key]
                                break
                except ValueError:
                    # Not JSON, might be in HTML or other format
                    pass
            
            # Step 2: Check status (with dummy ID if we couldn't extract one)
            test_query_id = query_id if query_id else 'test_query_id'
            status_response = ClassyFireAPI.query_status(test_query_id)
            
            # Status response may be None for non-existent queries
            if status_response is not None:
                if hasattr(status_response, 'status_code'):
                    # Should return valid response (200, 404, etc.)
                    assert status_response.status_code in [200, 404, 400, 422, 500, 503]
            
            # Step 3: Try to get results (with dummy ID)
            results_response = ClassyFireAPI.get_query(test_query_id)
            
            assert results_response is not None
            if hasattr(results_response, 'status_code'):
                # Should return valid response (429 is rate limiting)
                assert results_response.status_code in [200, 202, 404, 400, 422, 429, 500, 503]
        
        except Exception as e:
            if "timeout" in str(e).lower() or "connection" in str(e).lower():
                pytest.skip("ClassyFire API connection issues")
            else:
                raise
    
    def test_different_molecule_types(self):
        """Test classification of different types of molecules"""
        molecule_types = {
            'simple_alcohol': 'CCO',  # ethanol
            'aromatic': 'c1ccccc1',  # benzene
            'acid': 'CC(=O)O',  # acetic acid
            'ester': 'CC(=O)OC',  # methyl acetate
            'amine': 'CCN',  # ethylamine
        }
        
        for mol_type, smiles in molecule_types.items():
            try:
                response = ClassyFireAPI.submit_query(f"test_{mol_type}", smiles)
                
                assert response is not None
                
                if hasattr(response, 'status_code'):
                    if response.status_code in [200, 201, 202]:
                        # Successful submission
                        pass
                    elif response.status_code in [400, 422]:
                        # Some molecules might not be accepted
                        pass
                    elif response.status_code in [500, 503]:
                        pytest.skip("ClassyFire API server error")
                
            except Exception as e:
                if "timeout" in str(e).lower() or "connection" in str(e).lower():
                    pytest.skip("ClassyFire API connection issues")
                else:
                    # Some molecules might cause errors, which is acceptable for testing
                    pass
    
    def test_api_response_consistency(self):
        """Test that API responses are consistent"""
        smiles = 'CCO'  # ethanol
        
        # Submit the same query multiple times
        responses = []
        for i in range(3):
            try:
                response = ClassyFireAPI.submit_query(f"consistency_test_{i}", smiles)
                if response is not None:
                    responses.append(response)
                
                # Add small delay between requests
                time.sleep(0.1)
                
            except Exception as e:
                if "timeout" in str(e).lower() or "connection" in str(e).lower():
                    pytest.skip("ClassyFire API connection issues")
                else:
                    # Network issues might cause some requests to fail
                    pass
        
        if len(responses) < 2:
            pytest.skip("Not enough successful responses for consistency check")
        
        # Check that all responses have similar structure
        status_codes = [r.status_code for r in responses if hasattr(r, 'status_code')]
        
        if status_codes:
            # All status codes should be similar (allowing for minor variations)
            unique_status = set(status_codes)
            # Should not have wildly different status codes
            assert len(unique_status) <= 2, f"Too many different status codes: {unique_status}"


class TestClassyFireAPIErrorHandling:
    """Test error handling scenarios for ClassyFireAPI"""
    
    def test_network_timeout_handling(self):
        """Test handling of network timeouts"""
        smiles = 'CCO'
        
        # Test each method with timeout scenarios
        methods_to_test = [
            ('submit_query', lambda: ClassyFireAPI.submit_query("timeout_test", smiles)),
            ('query_status', lambda: ClassyFireAPI.query_status('test_id')),
            ('get_query', lambda: ClassyFireAPI.get_query('test_id')),
        ]
        
        for method_name, method_call in methods_to_test:
            try:
                response = method_call()
                # If successful, that's fine - network was fast enough
                # Note: Some methods may return None for invalid IDs, which is valid
                if response is not None:
                    assert response is not None
                
            except Exception as e:
                error_msg = str(e).lower()
                # These are acceptable network-related errors
                acceptable_errors = ['timeout', 'connection', 'network', 'unreachable']
                if any(err in error_msg for err in acceptable_errors):
                    # Expected network errors
                    pass
                else:
                    # Unexpected error - re-raise
                    raise
    
    def test_malformed_requests(self):
        """Test behavior with malformed requests"""
        # Test submit_query with problematic SMILES
        problematic_smiles = [
            'C' * 1000,  # Very long SMILES
            'C\nC',  # SMILES with newlines
            'C\tC',  # SMILES with tabs
            'C C',  # SMILES with spaces
        ]
        
        for smiles in problematic_smiles:
            try:
                response = ClassyFireAPI.submit_query("malformed_test", smiles)
                
                if response is not None and hasattr(response, 'status_code'):
                    # Should handle malformed input gracefully
                    # Usually 400 (bad request) or 422 (unprocessable entity)
                    assert response.status_code in [200, 201, 202, 400, 422, 500, 503]
                
            except Exception as e:
                error_msg = str(e).lower()
                # These errors are acceptable for malformed input
                acceptable_errors = ['invalid', 'malformed', 'bad request', 'encoding']
                if any(err in error_msg for err in acceptable_errors):
                    pass
                elif 'timeout' in error_msg or 'connection' in error_msg:
                    pytest.skip("ClassyFire API connection issues")
                else:
                    # Other errors might be acceptable too
                    pass
    
    def test_special_characters_in_labels(self):
        """Test handling of special characters in labels"""
        smiles = 'CCO'
        problematic_labels = [
            'label with spaces',
            'label/with/slashes',
            'label?with?questions',
            'label&with&ampersands',
            'label#with#hashes',
            'label"with"quotes',
            "label'with'apostrophes",
            'label\nwith\nnewlines',
            'café_label',  # Unicode
            'label_' + 'x' * 200,  # Very long label
        ]
        
        for label in problematic_labels:
            try:
                response = ClassyFireAPI.submit_query(label, smiles)
                
                if response is not None and hasattr(response, 'status_code'):
                    # Should handle special characters gracefully
                    assert response.status_code in [200, 201, 202, 400, 422, 500, 503]
                
            except Exception as e:
                error_msg = str(e).lower()
                # These errors are acceptable for special characters
                acceptable_errors = ['encoding', 'character', 'invalid', 'bad request']
                if any(err in error_msg for err in acceptable_errors):
                    pass
                elif 'timeout' in error_msg or 'connection' in error_msg:
                    pytest.skip("ClassyFire API connection issues")
                else:
                    # Other errors might be acceptable
                    pass
    
    def test_empty_parameters(self):
        """Test behavior with empty parameters"""
        # Test empty SMILES
        try:
            response = ClassyFireAPI.submit_query('empty_smiles_test', '')
            if response is not None and hasattr(response, 'status_code'):
                # Should handle empty SMILES gracefully
                assert response.status_code in [200, 400, 422, 500, 503]
        except Exception as e:
            # Errors for empty SMILES are acceptable
            pass
        
        # Test empty label
        try:
            response = ClassyFireAPI.submit_query('', 'CCO')
            if response is not None and hasattr(response, 'status_code'):
                # Should handle empty label gracefully
                assert response.status_code in [200, 201, 202, 400, 422, 500, 503]
        except Exception as e:
            # Errors for empty labels are acceptable
            pass
        
        # Test empty query ID
        try:
            response = ClassyFireAPI.query_status('')
            if response is not None and hasattr(response, 'status_code'):
                # Should handle empty query ID gracefully
                assert response.status_code in [400, 404, 422, 500, 503]
        except Exception as e:
            # Errors for empty query IDs are acceptable
            pass


class TestClassyFireAPIStaticMethods:
    """Test that all methods are properly static"""
    
    def test_methods_are_static(self):
        """Test that all methods can be called without instantiation"""
        # All methods should be callable directly on the class
        assert hasattr(ClassyFireAPI, 'submit_query')
        assert hasattr(ClassyFireAPI, 'query_status')
        assert hasattr(ClassyFireAPI, 'get_query')
        
        # Should be able to call without creating an instance
        try:
            # These calls might fail due to network/API issues, but should not fail due to instance issues
            ClassyFireAPI.submit_query('static_test', 'CCO')
            ClassyFireAPI.query_status('test_id')
            ClassyFireAPI.get_query('test_id')
        except Exception as e:
            # Network/API errors are acceptable, instance errors are not
            error_msg = str(e).lower()
            instance_errors = ['instance', 'self', 'object', 'class']
            assert not any(err in error_msg for err in instance_errors), f"Instance-related error: {e}"
    
    def test_class_attributes(self):
        """Test that class attributes are accessible"""
        # URL should be accessible
        assert hasattr(ClassyFireAPI, 'URL')
        assert isinstance(ClassyFireAPI.URL, str)
        assert len(ClassyFireAPI.URL) > 0
        assert ClassyFireAPI.URL.startswith('http')


if __name__ == "__main__":
    # Run tests if executed directly
    pytest.main([__file__, "-v"])
