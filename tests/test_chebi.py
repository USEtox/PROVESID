"""
Tests for ChEBI 2.0 API interface.

This module contains comprehensive tests for the ChEBI class and related functionality.
"""

import pytest
import requests
from unittest.mock import patch, Mock, MagicMock

from provesid.chebi import ChEBI, ChEBIError, get_chebi_entity, search_chebi


class TestChEBI:
    """Test cases for ChEBI class."""

    def test_initialization(self):
        """Test ChEBI class initialization."""
        chebi = ChEBI()
        assert chebi.base_url == "https://www.ebi.ac.uk/chebi/backend/api/public"
        assert chebi.timeout == 30
        assert hasattr(chebi, 'session')
        user_agent = chebi.session.headers.get('User-Agent', '')
        assert 'PROVESID-ChEBI-Client' in str(user_agent)

    def test_initialization_custom_timeout(self):
        """Test ChEBI initialization with custom timeout."""
        chebi = ChEBI(timeout=60)
        assert chebi.timeout == 60

    def test_format_chebi_id_integer(self):
        """Test _format_chebi_id with integer input."""
        assert ChEBI._format_chebi_id(15377) == "CHEBI:15377"

    def test_format_chebi_id_bare_string(self):
        """Test _format_chebi_id with bare number string."""
        assert ChEBI._format_chebi_id("15377") == "CHEBI:15377"

    def test_format_chebi_id_with_prefix(self):
        """Test _format_chebi_id with existing CHEBI: prefix."""
        assert ChEBI._format_chebi_id("CHEBI:15377") == "CHEBI:15377"

    def test_format_chebi_id_lowercase_prefix(self):
        """Test _format_chebi_id with lowercase prefix."""
        assert ChEBI._format_chebi_id("chebi:15377") == "chebi:15377"

    @patch('requests.Session.get')
    def test_get_success(self, mock_get):
        """Test successful _get call."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.headers = {"Content-Type": "application/json"}
        mock_response.json.return_value = {"chebi_id": "CHEBI:15377", "name": "water"}
        mock_get.return_value = mock_response

        chebi = ChEBI()
        result = chebi._get("compound/CHEBI:15377/")

        assert result == {"chebi_id": "CHEBI:15377", "name": "water"}
        mock_get.assert_called_once()

    @patch('requests.Session.get')
    def test_get_timeout(self, mock_get):
        """Test request timeout handling."""
        mock_get.side_effect = requests.exceptions.Timeout()

        chebi = ChEBI()
        with pytest.raises(ChEBIError, match="Request timeout"):
            chebi._get("compound/CHEBI:15377/")

    @patch('requests.Session.get')
    def test_get_network_error(self, mock_get):
        """Test network error handling."""
        mock_get.side_effect = requests.exceptions.ConnectionError("Network error")

        chebi = ChEBI()
        with pytest.raises(ChEBIError, match="Request failed"):
            chebi._get("compound/CHEBI:15377/")

    @patch('requests.Session.get')
    def test_get_compound_success(self, mock_get):
        """Test successful get_compound call."""
        compound_data = {
            "chebi_id": "CHEBI:15377",
            "name": "water",
            "definition": "An oxygen hydride.",
            "formula": "H2O",
            "mass": 18.01056,
            "monoisotopicmass": 18.01056,
            "charge": 0,
        }
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.headers = {"Content-Type": "application/json"}
        mock_response.json.return_value = compound_data
        mock_get.return_value = mock_response

        chebi = ChEBI()
        result = chebi.get_compound(15377)

        assert result is not None
        assert result['chebi_id'] == 'CHEBI:15377'
        assert result['name'] == 'water'
        assert result['formula'] == 'H2O'

    @patch('requests.Session.get')
    def test_get_compound_with_chebi_prefix(self, mock_get):
        """Test get_compound with CHEBI: prefix."""
        compound_data = {"chebi_id": "CHEBI:15377", "name": "water"}
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.headers = {"Content-Type": "application/json"}
        mock_response.json.return_value = compound_data
        mock_get.return_value = mock_response

        chebi = ChEBI()
        result = chebi.get_compound("CHEBI:15377")

        assert result is not None
        assert result['chebi_id'] == 'CHEBI:15377'

    @patch('requests.Session.get')
    def test_get_compound_not_found(self, mock_get):
        """Test get_compound when entity not found."""
        mock_get.side_effect = requests.exceptions.HTTPError(
            response=Mock(status_code=404)
        )

        chebi = ChEBI()
        result = chebi.get_compound(999999)

        assert result is None

    @patch('requests.Session.post')
    def test_get_compounds_success(self, mock_post):
        """Test successful get_compounds (batch) call."""
        compounds_data = [
            {"chebi_id": "CHEBI:15377", "name": "water"},
            {"chebi_id": "CHEBI:16236", "name": "ethanol"},
        ]
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.headers = {"Content-Type": "application/json"}
        mock_response.json.return_value = compounds_data
        mock_post.return_value = mock_response

        chebi = ChEBI()
        result = chebi.get_compounds(["CHEBI:15377", "CHEBI:16236"])

        assert result is not None
        assert len(result) == 2
        assert result[0]['chebi_id'] == 'CHEBI:15377'
        assert result[1]['chebi_id'] == 'CHEBI:16236'

    @patch('requests.Session.get')
    def test_search_success(self, mock_get):
        """Test successful search call."""
        search_data = {
            "results": [
                {"chebi_id": "CHEBI:15377", "name": "water"},
                {"chebi_id": "CHEBI:27313", "name": "water-18O"},
            ]
        }
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.headers = {"Content-Type": "application/json"}
        mock_response.json.return_value = search_data
        mock_get.return_value = mock_response

        chebi = ChEBI()
        result = chebi.search("water")

        assert result is not None
        assert 'results' in result
        assert len(result['results']) == 2

    @patch('requests.Session.get')
    def test_search_by_name_success(self, mock_get):
        """Test search_by_name convenience wrapper."""
        search_data = {
            "results": [
                {"chebi_id": "CHEBI:15377", "name": "water"},
            ]
        }
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.headers = {"Content-Type": "application/json"}
        mock_response.json.return_value = search_data
        mock_get.return_value = mock_response

        chebi = ChEBI()
        result = chebi.search_by_name("water")

        assert isinstance(result, list)
        assert len(result) >= 1

    @patch('requests.Session.get')
    def test_get_ontology_parents_success(self, mock_get):
        """Test successful get_ontology_parents call."""
        parents_data = [
            {"chebi_id": "CHEBI:24431", "name": "chemical entity", "type": "is_a"},
            {"chebi_id": "CHEBI:33259", "name": "elemental molecule", "type": "is_a"},
        ]
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.headers = {"Content-Type": "application/json"}
        mock_response.json.return_value = parents_data
        mock_get.return_value = mock_response

        chebi = ChEBI()
        result = chebi.get_ontology_parents(15377)

        assert isinstance(result, list)
        assert len(result) == 2
        assert result[0]['chebi_id'] == 'CHEBI:24431'
        assert result[1]['chebi_id'] == 'CHEBI:33259'

    @patch('requests.Session.get')
    def test_get_ontology_children_success(self, mock_get):
        """Test successful get_ontology_children call."""
        children_data = [
            {"chebi_id": "CHEBI:27313", "name": "water-18O", "type": "is_a"},
        ]
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.headers = {"Content-Type": "application/json"}
        mock_response.json.return_value = children_data
        mock_get.return_value = mock_response

        chebi = ChEBI()
        result = chebi.get_ontology_children(15377)

        assert isinstance(result, list)
        assert len(result) == 1
        assert result[0]['chebi_id'] == 'CHEBI:27313'

    @patch('requests.Session.get')
    def test_get_all_ontology_children_in_path(self, mock_get):
        """Test get_all_ontology_children_in_path call."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.headers = {"Content-Type": "application/json"}
        mock_response.json.return_value = {"results": []}
        mock_get.return_value = mock_response

        chebi = ChEBI()
        result = chebi.get_all_ontology_children_in_path(
            relation="is_a", entity="CHEBI:30879"
        )

        assert result is not None
        mock_get.assert_called_once()
        call_args = mock_get.call_args
        assert "relation" in str(call_args)

    @patch('provesid.chebi.time.sleep')
    @patch('requests.Session.get')
    def test_batch_get_compounds(self, mock_get, mock_sleep):
        """Test batch compound retrieval."""
        compound_data = {"chebi_id": "CHEBI:15377", "name": "water"}
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.headers = {"Content-Type": "application/json"}
        mock_response.json.return_value = compound_data
        mock_get.return_value = mock_response

        chebi = ChEBI()
        result = chebi.batch_get_compounds([15377, "CHEBI:16236"], pause_time=0.1)

        assert len(result) == 2
        assert "CHEBI:15377" in result
        assert "CHEBI:16236" in result
        assert mock_sleep.call_count == 2

    @patch('requests.Session.get')
    def test_get_compound_structure_success(self, mock_get):
        """Test successful compound structure retrieval (SVG)."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.text = '<svg>...</svg>'
        mock_get.return_value = mock_response

        chebi = ChEBI()
        result = chebi.get_compound_structure(15377)

        assert result == '<svg>...</svg>'

    @patch('requests.Session.get')
    def test_get_molfile_success(self, mock_get):
        """Test successful molfile retrieval."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.text = 'molfile contents'
        mock_get.return_value = mock_response

        chebi = ChEBI()
        result = chebi.get_molfile(15377)

        assert result == 'molfile contents'

    @patch('requests.Session.get')
    def test_structure_search_success(self, mock_get):
        """Test successful structure search."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.headers = {"Content-Type": "application/json"}
        mock_response.json.return_value = {"results": []}
        mock_get.return_value = mock_response

        chebi = ChEBI()
        result = chebi.structure_search("c1ccccc1", "substructure")

        assert result is not None

    @patch('requests.Session.post')
    def test_advanced_search_success(self, mock_post):
        """Test successful advanced search."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.headers = {"Content-Type": "application/json"}
        mock_response.json.return_value = {"results": []}
        mock_post.return_value = mock_response

        chebi = ChEBI()
        result = chebi.advanced_search({
            "formula_specification": {
                "and_specification": [{"term": "C6H12O7"}]
            }
        }, three_star_only=False)

        assert result is not None

    @patch('requests.Session.post')
    def test_calculate_avg_mass(self, mock_post):
        """Test average mass calculation."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.text = "18.015"
        mock_post.return_value = mock_response

        chebi = ChEBI()
        result = chebi.calculate_avg_mass("O")

        assert result == "18.015"

    @patch('requests.Session.post')
    def test_calculate_mol_formula(self, mock_post):
        """Test molecular formula calculation."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.text = "H2O"
        mock_post.return_value = mock_response

        chebi = ChEBI()
        result = chebi.calculate_mol_formula("O")

        assert result == "H2O"

    @patch('requests.Session.post')
    def test_calculate_net_charge(self, mock_post):
        """Test net charge calculation."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.text = "0"
        mock_post.return_value = mock_response

        chebi = ChEBI()
        result = chebi.calculate_net_charge("O")

        assert result == "0"

    def test_repr(self):
        """Test string representation of ChEBI object."""
        chebi = ChEBI(timeout=45)
        repr_str = repr(chebi)
        assert "ChEBI" in repr_str
        assert "timeout=45" in repr_str
        assert chebi.base_url in repr_str


class TestChEBIIntegration:
    """Integration tests for ChEBI 2.0 API (requires internet connection)."""

    @pytest.mark.skip(reason="Integration test - requires internet connection")
    def test_get_water_compound(self):
        """Test getting compound data for water (CHEBI:15377)."""
        chebi = ChEBI()
        result = chebi.get_compound(15377)

        assert result is not None
        assert 'CHEBI:15377' in str(result.get('chebi_id', result.get('chebi_accession', '')))

    @pytest.mark.skip(reason="Integration test - requires internet connection")
    def test_search_water(self):
        """Test searching for water by name."""
        chebi = ChEBI()
        results = chebi.search("water")

        assert results is not None

    @pytest.mark.skip(reason="Integration test - requires internet connection")
    def test_get_ontology_structure(self):
        """Test getting ontology parents and children."""
        chebi = ChEBI()

        parents = chebi.get_ontology_parents(15377)
        assert parents is not None

        children = chebi.get_ontology_children(15377)
        assert children is not None

    @pytest.mark.skip(reason="Integration test - requires internet connection")
    def test_structure_search(self):
        """Test structure search via SMILES."""
        chebi = ChEBI()
        results = chebi.structure_search("c1ccccc1", "connectivity")
        assert results is not None


class TestChEBIErrorHandling:
    """Test error handling in ChEBI class."""

    @patch('requests.Session.get')
    def test_network_timeout_handling(self, mock_get):
        """Test handling of network timeouts."""
        mock_get.side_effect = requests.exceptions.Timeout()

        chebi = ChEBI()
        result = chebi.get_compound(15377)

        assert result is None  # Should return None on error, not raise

    @patch('requests.Session.get')
    def test_http_error_handling(self, mock_get):
        """Test handling of HTTP errors."""
        mock_get.side_effect = requests.exceptions.HTTPError(
            response=Mock(status_code=500)
        )

        chebi = ChEBI()
        result = chebi.get_compound(15377)

        assert result is None

    def test_invalid_chebi_id_formats(self):
        """Test handling of various ChEBI ID formats."""
        chebi = ChEBI()

        test_ids = [
            15377,
            "15377",
            "CHEBI:15377",
            "chebi:15377",
            "invalid_id",
            "",
        ]

        for test_id in test_ids:
            try:
                result = chebi.get_compound(test_id)
            except ChEBIError:
                pass
            except Exception as e:
                pytest.fail(f"Unexpected exception for ID {test_id}: {e}")

    def test_empty_search_text(self):
        """Test handling of empty search text."""
        chebi = ChEBI()

        result = chebi.search_by_name("", size=5)
        assert isinstance(result, list)

    def test_very_long_search_text(self):
        """Test handling of very long search text."""
        chebi = ChEBI()

        long_text = "a" * 1000
        result = chebi.search_by_name(long_text, size=5)
        assert isinstance(result, list)


class TestChEBIConvenienceFunctions:
    """Test convenience functions for ChEBI."""

    @patch('provesid.chebi.ChEBI.get_compound')
    def test_get_chebi_entity(self, mock_get_compound):
        """Test get_chebi_entity convenience function."""
        mock_get_compound.return_value = {
            'chebi_id': 'CHEBI:15377',
            'name': 'water',
        }

        result = get_chebi_entity(15377)

        assert result is not None
        assert result['chebi_id'] == 'CHEBI:15377'
        mock_get_compound.assert_called_once_with(15377)

    @patch('provesid.chebi.ChEBI.search_by_name')
    def test_search_chebi(self, mock_search):
        """Test search_chebi convenience function."""
        mock_search.return_value = [
            {'chebi_id': 'CHEBI:15377', 'name': 'water'},
            {'chebi_id': 'CHEBI:27313', 'name': 'water-18O'},
        ]

        result = search_chebi("water", max_results=5)

        assert isinstance(result, list)
        assert len(result) == 2
        mock_search.assert_called_once_with("water", size=5)


class TestChEBIStructureEndpoints:
    """Test structure-related endpoints."""

    @patch('requests.Session.get')
    def test_get_compound_structure(self, mock_get):
        """Test getting compound SVG structure."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.text = '<svg xmlns="http://www.w3.org/2000/svg">...</svg>'
        mock_get.return_value = mock_response

        chebi = ChEBI()
        result = chebi.get_compound_structure(15377, width=400, height=400)

        assert result is not None
        assert '<svg' in result

    @patch('requests.Session.get')
    def test_get_structure_by_pk(self, mock_get):
        """Test getting structure SVG by primary key."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.text = '<svg>structure</svg>'
        mock_get.return_value = mock_response

        chebi = ChEBI()
        result = chebi.get_structure(12345, width=200, height=200)

        assert result == '<svg>structure</svg>'

    @patch('requests.Session.get')
    def test_get_structure_not_available(self, mock_get):
        """Test getting structure when not available."""
        mock_get.side_effect = requests.exceptions.HTTPError(
            response=Mock(status_code=404)
        )

        chebi = ChEBI()
        result = chebi.get_compound_structure(999999)

        assert result is None

    @patch('requests.Session.get')
    def test_structure_search_with_similarity(self, mock_get):
        """Test structure search with similarity threshold."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.headers = {"Content-Type": "application/json"}
        mock_response.json.return_value = {"results": []}
        mock_get.return_value = mock_response

        chebi = ChEBI()
        result = chebi.structure_search(
            "c1ccccc1", "similarity", similarity=0.8
        )

        assert result is not None
        call_args = mock_get.call_args
        params = call_args[1].get('params', {})
        assert params.get('similarity') == 0.8

    @patch('requests.Session.post')
    def test_depict_structure(self, mock_post):
        """Test structure depiction (PNG)."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b'\x89PNG\x0d\x0a\x1a\x0a'  # PNG header
        mock_post.return_value = mock_response

        chebi = ChEBI()
        result = chebi.depict_structure("CCO", width=400, height=400)

        assert result is not None
        assert isinstance(result, bytes)
