"""Tests for PubChem API functionality"""

import pytest
from provesid.pubchem import PubChemAPI, PubChemError, PubChemNotFoundError

@pytest.mark.unit  
class TestPubChemAPI:
    
    @pytest.fixture
    def api(self):
        """Create a PubChemAPI instance for testing"""
        return PubChemAPI()
    
    def test_initialization(self, api):
        """Test PubChemAPI initialization"""
        assert api.base_url == "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        assert api.pause_time == 0.2
        
    def test_custom_initialization(self):
        """Test PubChemAPI with custom parameters"""
        custom_api = PubChemAPI(
            base_url="https://example.com",
            pause_time=0.5
        )
        assert custom_api.base_url == "https://example.com"
        assert custom_api.pause_time == 0.5

    @pytest.mark.api
    def test_compound_by_cid(self, api):
        """Test retrieving compound by CID - updated for new flat API"""
        result = api.get_compound_by_cid(2244, output_format='JSON')
        
        # New API returns compound object directly, not wrapped in PC_Compounds
        assert isinstance(result, dict)
        assert result['id']['id']['cid'] == 2244

    @pytest.mark.api
    def test_cids_by_name(self, api):
        """Test getting CIDs by compound name - updated for new flat API"""
        result = api.get_cids_by_name('aspirin')
        
        # New API returns CIDs directly as a list, not wrapped in IdentifierList.CID
        assert isinstance(result, list)
        assert 2244 in result  # Known aspirin CID

    @pytest.mark.api
    def test_properties(self, api):
        """Test getting compound properties - updated for new flat API"""
        # Single property - now uses single CID, not list
        mw = api.get_compound_properties(2244, ['MolecularWeight'])
        
        # New API returns flat structure with properties directly accessible
        assert isinstance(mw, dict)
        assert mw['success'] == True
        assert mw['CID'] == 2244
        assert 'MolecularWeight' in mw
        assert float(mw['MolecularWeight']) > 0

    @pytest.mark.api
    def test_synonyms(self, api):
        """Test getting compound synonyms - updated for new flat API"""
        synonyms = api.get_compound_synonyms(2244)
        
        # New API returns flat list of synonyms, not nested structure
        assert isinstance(synonyms, list)
        assert len(synonyms) > 0
        
        # Aspirin should have "aspirin" as a synonym
        synonyms_list = [s.lower() for s in synonyms]
        assert 'aspirin' in synonyms_list

    def test_error_handling_invalid_cid(self, api):
        """Test error handling for invalid CID"""
        with pytest.raises(PubChemNotFoundError):
            api.get_compound_by_cid(999999999)

    def test_error_handling_invalid_name(self, api):
        """Test error handling for invalid name"""
        with pytest.raises(PubChemNotFoundError):
            api.get_cids_by_name('this_is_definitely_not_a_chemical_name_12345')

    def test_malformed_property_names(self, api):
        """Test with malformed property names - updated for new error handling"""
        result = api.get_compound_properties(2244, ['InvalidPropertyName'])
        
        # New API returns error result instead of raising exception
        assert isinstance(result, dict)
        assert result['success'] == False
        assert 'error' in result
        assert 'Invalid property' in result['error']
