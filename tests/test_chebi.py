"""
Tests for ChEBI API interface.

This module contains comprehensive tests for the ChEBI class and related functionality.
"""

import pytest
import requests
from unittest.mock import patch, Mock
import xml.etree.ElementTree as ET

from provesid.chebi import ChEBI, ChEBIError, get_chebi_entity, search_chebi


class TestChEBI:
    """Test cases for ChEBI class."""
    
    def test_initialization(self):
        """Test ChEBI class initialization."""
        chebi = ChEBI()
        assert chebi.base_url == "https://www.ebi.ac.uk/webservices/chebi/2.0/test"
        assert chebi.timeout == 30
        assert hasattr(chebi, 'session')
        user_agent = chebi.session.headers.get('User-Agent', '')
        assert 'PROVESID-ChEBI-Client' in str(user_agent)
    
    def test_initialization_custom_timeout(self):
        """Test ChEBI initialization with custom timeout."""
        chebi = ChEBI(timeout=60)
        assert chebi.timeout == 60
    
    def test_xml_to_dict_simple(self):
        """Test XML to dictionary conversion with simple element."""
        chebi = ChEBI()
        xml_str = "<name>water</name>"
        element = ET.fromstring(xml_str)
        result = chebi._xml_to_dict(element)
        assert result == "water"
    
    def test_xml_to_dict_complex(self):
        """Test XML to dictionary conversion with complex structure."""
        chebi = ChEBI()
        xml_str = """
        <entity>
            <chebiId>CHEBI:15377</chebiId>
            <chebiAsciiName>water</chebiAsciiName>
            <formulaConnectivity>H2O</formulaConnectivity>
        </entity>
        """
        element = ET.fromstring(xml_str)
        result = chebi._xml_to_dict(element)
        
        assert isinstance(result, dict)
        assert result['chebiId'] == 'CHEBI:15377'
        assert result['chebiAsciiName'] == 'water'
        assert result['formulaConnectivity'] == 'H2O'
    
    def test_xml_to_dict_with_list(self):
        """Test XML to dictionary conversion with repeated elements."""
        chebi = ChEBI()
        xml_str = """
        <synonyms>
            <synonym>water</synonym>
            <synonym>dihydrogen oxide</synonym>
            <synonym>H2O</synonym>
        </synonyms>
        """
        element = ET.fromstring(xml_str)
        result = chebi._xml_to_dict(element)
        
        assert isinstance(result, dict)
        assert isinstance(result['synonym'], list)
        assert len(result['synonym']) == 3
        assert 'water' in result['synonym']
        assert 'dihydrogen oxide' in result['synonym']
    
    @patch('requests.Session.get')
    def test_make_request_success(self, mock_get):
        """Test successful API request."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b'<response>success</response>'
        mock_get.return_value = mock_response
        
        chebi = ChEBI()
        response = chebi._make_request("test_endpoint")
        
        assert response == mock_response
        mock_get.assert_called_once()
    
    @patch('requests.Session.get')
    def test_make_request_timeout(self, mock_get):
        """Test request timeout handling."""
        mock_get.side_effect = requests.exceptions.Timeout()
        
        chebi = ChEBI()
        with pytest.raises(ChEBIError, match="Request timeout"):
            chebi._make_request("test_endpoint")
    
    @patch('requests.Session.get')
    def test_make_request_network_error(self, mock_get):
        """Test network error handling."""
        mock_get.side_effect = requests.exceptions.ConnectionError("Network error")
        
        chebi = ChEBI()
        with pytest.raises(ChEBIError, match="Request failed"):
            chebi._make_request("test_endpoint")
    
    @patch('requests.Session.get')
    def test_get_complete_entity_success(self, mock_get):
        """Test successful get_complete_entity call."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b'''
        <soap:Envelope xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/">
            <soap:Body>
                <getCompleteEntityResponse>
                    <return>
                        <chebiId>CHEBI:15377</chebiId>
                        <chebiAsciiName>water</chebiAsciiName>
                        <definition>An oxygen hydride consisting of an oxygen atom that is covalently bonded to two hydrogen atoms.</definition>
                        <formulaConnectivity>H2O</formulaConnectivity>
                        <mass>18.01056</mass>
                        <monoisotopicMass>18.01056</monoisotopicMass>
                        <charge>0</charge>
                    </return>
                </getCompleteEntityResponse>
            </soap:Body>
        </soap:Envelope>
        '''
        mock_get.return_value = mock_response
        
        chebi = ChEBI()
        result = chebi.get_complete_entity(15377)
        
        assert result is not None
        assert result['chebiId'] == 'CHEBI:15377'
        assert result['chebiAsciiName'] == 'water'
        assert result['formulaConnectivity'] == 'H2O'
        assert result['mass'] == '18.01056'
    
    @patch('requests.Session.get')
    def test_get_complete_entity_with_chebi_prefix(self, mock_get):
        """Test get_complete_entity with CHEBI: prefix."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b'''
        <soap:Envelope xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/">
            <soap:Body>
                <getCompleteEntityResponse>
                    <return>
                        <chebiId>CHEBI:15377</chebiId>
                        <chebiAsciiName>water</chebiAsciiName>
                    </return>
                </getCompleteEntityResponse>
            </soap:Body>
        </soap:Envelope>
        '''
        mock_get.return_value = mock_response
        
        chebi = ChEBI()
        result = chebi.get_complete_entity("CHEBI:15377")
        
        assert result is not None
        assert result['chebiId'] == 'CHEBI:15377'
    
    @patch('requests.Session.get')
    def test_get_complete_entity_not_found(self, mock_get):
        """Test get_complete_entity when entity not found."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b'''
        <soap:Envelope xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/">
            <soap:Body>
                <getCompleteEntityResponse>
                </getCompleteEntityResponse>
            </soap:Body>
        </soap:Envelope>
        '''
        mock_get.return_value = mock_response
        
        chebi = ChEBI()
        result = chebi.get_complete_entity(999999)
        
        assert result is None
    
    @patch('requests.Session.get')
    def test_get_lite_entity_success(self, mock_get):
        """Test successful get_lite_entity call."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b'''
        <soap:Envelope xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/">
            <soap:Body>
                <getLiteEntityResponse>
                    <return>
                        <chebiId>CHEBI:15377</chebiId>
                        <chebiAsciiName>water</chebiAsciiName>
                        <searchScore>1.0</searchScore>
                        <entityStar>3</entityStar>
                    </return>
                </getLiteEntityResponse>
            </soap:Body>
        </soap:Envelope>
        '''
        mock_get.return_value = mock_response
        
        chebi = ChEBI()
        result = chebi.get_lite_entity(15377)
        
        assert result is not None
        assert result['chebiId'] == 'CHEBI:15377'
        assert result['chebiAsciiName'] == 'water'
        assert result['entityStar'] == '3'
    
    @patch('requests.Session.get')
    def test_get_structure_success(self, mock_get):
        """Test successful get_structure call."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b'''
        <soap:Envelope xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/">
            <soap:Body>
                <getStructureResponse>
                    <return>O</return>
                </getStructureResponse>
            </soap:Body>
        </soap:Envelope>
        '''
        mock_get.return_value = mock_response
        
        chebi = ChEBI()
        result = chebi.get_structure(15377, "smiles")
        
        assert result == "O"
    
    @patch('requests.Session.get')
    def test_get_ontology_parents_success(self, mock_get):
        """Test successful get_ontology_parents call."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b'''
        <soap:Envelope xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/">
            <soap:Body>
                <getOntologyParentsResponse>
                    <return>
                        <chebiId>CHEBI:24431</chebiId>
                        <chebiName>chemical entity</chebiName>
                        <type>is_a</type>
                        <status>C</status>
                    </return>
                    <return>
                        <chebiId>CHEBI:33259</chebiId>
                        <chebiName>elemental molecule</chebiName>
                        <type>is_a</type>
                        <status>C</status>
                    </return>
                </getOntologyParentsResponse>
            </soap:Body>
        </soap:Envelope>
        '''
        mock_get.return_value = mock_response
        
        chebi = ChEBI()
        result = chebi.get_ontology_parents(15377)
        
        assert isinstance(result, list)
        assert len(result) == 2
        assert result[0]['chebiId'] == 'CHEBI:24431'
        assert result[1]['chebiId'] == 'CHEBI:33259'
    
    @patch('requests.Session.get')
    def test_get_ontology_children_success(self, mock_get):
        """Test successful get_ontology_children call."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b'''
        <soap:Envelope xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/">
            <soap:Body>
                <getOntologyChildrenResponse>
                    <return>
                        <chebiId>CHEBI:27313</chebiId>
                        <chebiName>water-18O</chebiName>
                        <type>is_a</type>
                        <status>C</status>
                    </return>
                </getOntologyChildrenResponse>
            </soap:Body>
        </soap:Envelope>
        '''
        mock_get.return_value = mock_response
        
        chebi = ChEBI()
        result = chebi.get_ontology_children(15377)
        
        assert isinstance(result, list)
        assert len(result) == 1
        assert result[0]['chebiId'] == 'CHEBI:27313'
    
    @patch('provesid.chebi.time.sleep')
    @patch('requests.Session.get')
    def test_batch_get_entities(self, mock_get, mock_sleep):
        """Test batch entity retrieval."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b'''
        <soap:Envelope xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/">
            <soap:Body>
                <getCompleteEntityResponse>
                    <return>
                        <chebiId>CHEBI:15377</chebiId>
                        <chebiAsciiName>water</chebiAsciiName>
                    </return>
                </getCompleteEntityResponse>
            </soap:Body>
        </soap:Envelope>
        '''
        mock_get.return_value = mock_response
        
        chebi = ChEBI()
        result = chebi.batch_get_entities([15377, "CHEBI:16236"], pause_time=0.1)
        
        assert len(result) == 2
        assert "CHEBI:15377" in result
        assert "CHEBI:16236" in result
        assert mock_sleep.call_count == 2  # Called after each request
    
    def test_repr(self):
        """Test string representation of ChEBI object."""
        chebi = ChEBI(timeout=45)
        repr_str = repr(chebi)
        assert "ChEBI" in repr_str
        assert "timeout=45" in repr_str
        assert chebi.base_url in repr_str


class TestChEBIIntegration:
    """Integration tests for ChEBI API (requires internet connection)."""
    
    @pytest.mark.skip(reason="Integration test - requires internet connection")
    def test_get_water_complete_entity(self):
        """Test getting complete entity for water (CHEBI:15377)."""
        chebi = ChEBI()
        result = chebi.get_complete_entity(15377)
        
        assert result is not None
        assert result['chebiId'] == 'CHEBI:15377'
        assert 'water' in result['chebiAsciiName'].lower()
        assert 'H2O' in result.get('formulaConnectivity', '')
    
    @pytest.mark.skip(reason="Integration test - requires internet connection")
    def test_get_ethanol_lite_entity(self):
        """Test getting lite entity for ethanol (CHEBI:16236)."""
        chebi = ChEBI()
        result = chebi.get_lite_entity("CHEBI:16236")
        
        assert result is not None
        assert result['chebiId'] == 'CHEBI:16236'
        assert 'ethanol' in result['chebiAsciiName'].lower()
    
    @pytest.mark.skip(reason="Integration test - requires internet connection")
    def test_search_by_name_water(self):
        """Test searching for water by name."""
        chebi = ChEBI()
        results = chebi.search_by_name("water", max_results=5)
        
        assert isinstance(results, list)
        assert len(results) > 0
        
        # Check if water (CHEBI:15377) is in results
        water_found = False
        for result in results:
            if result.get('chebiId') == 'CHEBI:15377':
                water_found = True
                break
        
        assert water_found, "Water (CHEBI:15377) should be found in search results"
    
    @pytest.mark.skip(reason="Integration test - requires internet connection")
    def test_get_ontology_structure(self):
        """Test getting ontology parents and children."""
        chebi = ChEBI()
        
        # Test parents
        parents = chebi.get_ontology_parents(15377)  # water
        assert isinstance(parents, list)
        
        # Test children
        children = chebi.get_ontology_children(15377)  # water
        assert isinstance(children, list)


class TestChEBIErrorHandling:
    """Test error handling in ChEBI class."""
    
    @patch('requests.Session.get')
    def test_network_timeout_handling(self, mock_get):
        """Test handling of network timeouts."""
        mock_get.side_effect = requests.exceptions.Timeout()
        
        chebi = ChEBI()
        result = chebi.get_complete_entity(15377)
        
        assert result is None  # Should return None on error, not raise
    
    @patch('requests.Session.get')
    def test_malformed_xml_handling(self, mock_get):
        """Test handling of malformed XML responses."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b'<invalid>xml<content>'  # Malformed XML
        mock_get.return_value = mock_response
        
        chebi = ChEBI()
        result = chebi.get_complete_entity(15377)
        
        assert result is None  # Should return None on XML parse error
    
    def test_invalid_chebi_id_formats(self):
        """Test handling of various ChEBI ID formats."""
        chebi = ChEBI()
        
        # Test with different ID formats - should not raise exceptions
        test_ids = [
            15377,
            "15377",
            "CHEBI:15377",
            "chebi:15377",  # lowercase
            "invalid_id",
            "",
            None
        ]
        
        for test_id in test_ids:
            try:
                # These calls might fail due to network/API issues
                # but should not raise exceptions due to ID format
                result = chebi.get_complete_entity(test_id)
                # Just ensure no exception is raised
            except ChEBIError:
                # ChEBIError is acceptable (network issues, etc.)
                pass
            except Exception as e:
                pytest.fail(f"Unexpected exception for ID {test_id}: {e}")
    
    def test_empty_search_text(self):
        """Test handling of empty search text."""
        chebi = ChEBI()
        
        # Empty search should return empty list, not raise exception
        result = chebi.search_by_name("", max_results=5)
        assert isinstance(result, list)
    
    def test_very_long_search_text(self):
        """Test handling of very long search text."""
        chebi = ChEBI()
        
        # Very long search text should not crash
        long_text = "a" * 1000
        result = chebi.search_by_name(long_text, max_results=5)
        assert isinstance(result, list)


class TestChEBIConvenienceFunctions:
    """Test convenience functions for ChEBI."""
    
    @patch('provesid.chebi.ChEBI.get_complete_entity')
    def test_get_chebi_entity(self, mock_get_entity):
        """Test get_chebi_entity convenience function."""
        mock_get_entity.return_value = {
            'chebiId': 'CHEBI:15377',
            'chebiAsciiName': 'water'
        }
        
        result = get_chebi_entity(15377)
        
        assert result is not None
        assert result['chebiId'] == 'CHEBI:15377'
        mock_get_entity.assert_called_once_with(15377)
    
    @patch('provesid.chebi.ChEBI.search_by_name')
    def test_search_chebi(self, mock_search):
        """Test search_chebi convenience function."""
        mock_search.return_value = [
            {'chebiId': 'CHEBI:15377', 'chebiAsciiName': 'water'},
            {'chebiId': 'CHEBI:27313', 'chebiAsciiName': 'water-18O'}
        ]
        
        result = search_chebi("water", max_results=5)
        
        assert isinstance(result, list)
        assert len(result) == 2
        mock_search.assert_called_once_with("water", max_results=5)


class TestChEBIStructureTypes:
    """Test different structure format requests."""
    
    @patch('requests.Session.get')
    def test_get_structure_different_formats(self, mock_get):
        """Test getting structures in different formats."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b'''
        <soap:Envelope xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/">
            <soap:Body>
                <getStructureResponse>
                    <return>O</return>
                </getStructureResponse>
            </soap:Body>
        </soap:Envelope>
        '''
        mock_get.return_value = mock_response
        
        chebi = ChEBI()
        
        # Test different structure formats
        formats = ["mol", "sdf", "smiles", "inchi"]
        for fmt in formats:
            result = chebi.get_structure(15377, fmt)
            assert result == "O"
    
    @patch('requests.Session.get')
    def test_get_structure_not_available(self, mock_get):
        """Test getting structure when not available."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b'''
        <soap:Envelope xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/">
            <soap:Body>
                <getStructureResponse>
                </getStructureResponse>
            </soap:Body>
        </soap:Envelope>
        '''
        mock_get.return_value = mock_response
        
        chebi = ChEBI()
        result = chebi.get_structure(15377, "mol")
        
        assert result is None
