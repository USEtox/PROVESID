"""
Tests for PubChem API functionality
"""

import pytest
import sys
import os

# Add the src directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from provesid.pubchem import (
    PubChemAPI,
    PubChemError,
    PubChemNotFoundError,
    PubChemTimeoutError,
    PubChemServerError,
    CompoundProperties,
    OutputFormat,
    Domain,
    CompoundDomainNamespace,
    Operation
)


class TestPubChemAPI:
    """Test suite for PubChemAPI class"""
    
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
    
    def test_compound_by_name(self, api):
        """Test retrieving compound by name"""
        # Test with a well-known compound
        result = api.get_compounds_by_name('aspirin', output_format='JSON')
        
        assert isinstance(result, dict)
        assert 'PC_Compounds' in result
        assert len(result['PC_Compounds']) > 0
        
        compound = result['PC_Compounds'][0]
        assert 'id' in compound
        assert 'atoms' in compound
        assert 'bonds' in compound
    
    def test_compound_by_cid(self, api):
        """Test retrieving compound by CID"""
        # Test with aspirin CID (2244)
        result = api.get_compound_by_cid(2244, output_format='JSON')
        
        assert isinstance(result, dict)
        assert 'PC_Compounds' in result
        assert len(result['PC_Compounds']) > 0
        
        compound = result['PC_Compounds'][0]
        assert compound['id']['id']['cid'] == 2244
    
    def test_compound_by_smiles(self, api):
        """Test retrieving compound by SMILES"""
        # Test with simple molecule (ethanol)
        result = api.get_compounds_by_smiles('CCO', output_format='JSON')
        
        assert isinstance(result, dict)
        assert 'PC_Compounds' in result
        assert len(result['PC_Compounds']) > 0
    
    def test_cids_by_name(self, api):
        """Test getting CIDs by compound name"""
        result = api.get_cids_by_name('aspirin')
        
        assert isinstance(result, dict)
        assert 'IdentifierList' in result
        assert 'CID' in result['IdentifierList']
        assert 2244 in result['IdentifierList']['CID']  # Known aspirin CID
    
    def test_cids_by_smiles(self, api):
        """Test getting CIDs by SMILES"""
        result = api.get_cids_by_smiles('CCO')  # ethanol
        
        assert isinstance(result, dict)
        assert 'IdentifierList' in result
        assert 'CID' in result['IdentifierList']
        assert len(result['IdentifierList']['CID']) > 0
    
    def test_properties(self, api):
        """Test getting compound properties"""
        # Test single property
        mw = api.get_compound_properties([2244], ['MolecularWeight'])
        
        assert isinstance(mw, dict)
        assert 'PropertyTable' in mw
        assert len(mw['PropertyTable']['Properties']) > 0
        
        prop = mw['PropertyTable']['Properties'][0]
        assert prop['CID'] == 2244
        assert 'MolecularWeight' in prop
        assert float(prop['MolecularWeight']) > 0
    
    def test_multiple_properties(self, api):
        """Test getting multiple properties"""
        props = api.get_compound_properties([2244], ['MolecularWeight', 'MolecularFormula', 'ConnectivitySMILES'])
        
        assert isinstance(props, dict)
        assert 'PropertyTable' in props
        
        prop = props['PropertyTable']['Properties'][0]
        assert prop['CID'] == 2244
        assert 'MolecularWeight' in prop
        assert 'MolecularFormula' in prop
        assert 'ConnectivitySMILES' in prop
    
    def test_synonyms(self, api):
        """Test getting compound synonyms"""
        synonyms = api.get_compound_synonyms(2244)
        
        assert isinstance(synonyms, dict)
        assert 'InformationList' in synonyms
        assert len(synonyms['InformationList']['Information']) > 0
        
        info = synonyms['InformationList']['Information'][0]
        assert info['CID'] == 2244
        assert 'Synonym' in info
        assert len(info['Synonym']) > 0
        
        # Aspirin should have "aspirin" as a synonym
        synonyms_list = [s.lower() for s in info['Synonym']]
        assert 'aspirin' in synonyms_list
    
    def test_search_compound(self, api):
        """Test compound search functionality"""
        result = api.search_compound('aspirin', search_type='name')
        
        assert isinstance(result, dict)
        assert result['success'] == True
        assert result['query'] == 'aspirin'
        assert result['search_type'] == 'name'
        assert result['data'] is not None
    
    def test_basic_compound_info(self, api):
        """Test getting basic compound information"""
        info = api.get_basic_compound_info(2244)  # aspirin
        
        assert isinstance(info, dict)
        assert info['success'] == True
        assert info['cid'] == 2244
        assert info['properties'] is not None
        assert info['synonyms'] is not None
    
    def test_similarity_search(self, api):
        """Test similarity search"""
        # Search for compounds similar to aspirin
        results = api.similarity_search(2244, threshold=95)
        
        assert isinstance(results, dict)
        assert 'IdentifierList' in results
        assert 'CID' in results['IdentifierList']
        assert len(results['IdentifierList']['CID']) > 0
        
        # Should include aspirin itself
        assert 2244 in results['IdentifierList']['CID']
    
    def test_substructure_search(self, api):
        """Test substructure search"""
        # Search for compounds containing benzene ring
        results = api.substructure_search('c1ccccc1')  # benzene SMILES
        
        assert isinstance(results, dict)
        assert 'IdentifierList' in results
        assert 'CID' in results['IdentifierList']
        assert len(results['IdentifierList']['CID']) > 0
    
    def test_batch_operations(self, api):
        """Test batch operations with multiple CIDs"""
        cids = [2244, 702]  # aspirin and ethanol
        
        # Test batch properties
        props = api.get_compound_properties(cids, ['MolecularWeight', 'MolecularFormula'])
        
        assert isinstance(props, dict)
        assert 'PropertyTable' in props
        assert len(props['PropertyTable']['Properties']) == 2
        
        # Check both compounds are present
        found_cids = [prop['CID'] for prop in props['PropertyTable']['Properties']]
        assert 2244 in found_cids
        assert 702 in found_cids
    
    def test_error_handling(self, api):
        """Test error handling for invalid inputs"""
        # Test with non-existent CID
        with pytest.raises(PubChemNotFoundError):
            api.get_compound_by_cid(999999999)  # Very unlikely to exist
        
        # Test with invalid name
        with pytest.raises(PubChemNotFoundError):
            api.get_cids_by_name('this_is_definitely_not_a_chemical_name_12345')
    
    def test_rate_limiting(self, api):
        """Test rate limiting functionality"""
        import time
        
        start_time = time.time()
        
        # Make multiple API calls
        for _ in range(3):
            api._rate_limit()
        
        elapsed = time.time() - start_time
        
        # Should take at least 2 * pause_time
        expected_min_time = 2 * api.pause_time
        assert elapsed >= expected_min_time * 0.9  # Allow some tolerance


class TestSpecialCases:
    """Test suite for special cases and edge conditions"""
    
    def test_common_drugs(self):
        """Test with common pharmaceutical compounds"""
        drugs = {
            'aspirin': 2244,
            'caffeine': 2519,
            'ibuprofen': 3672
        }
        
        api = PubChemAPI()
        
        for drug_name, expected_cid in drugs.items():
            result = api.get_cids_by_name(drug_name)
            cids = result['IdentifierList']['CID']
            assert expected_cid in cids
    
    def test_simple_molecules(self):
        """Test with simple molecules"""
        molecules = {
            'water': 962,
            'ethanol': 702,
            'methanol': 887
        }
        
        api = PubChemAPI()
        
        for mol_name, expected_cid in molecules.items():
            result = api.get_cids_by_name(mol_name)
            cids = result['IdentifierList']['CID']
            assert expected_cid in cids
    
    def test_large_batch_operations(self):
        """Test batch operations with larger lists"""
        api = PubChemAPI()
        
        # Test with multiple CIDs
        cids = [2244, 702, 887, 2519, 3672]  # Various compounds

        props = api.get_compound_properties(cids, ['MolecularWeight']) # type: ignore

        assert isinstance(props, dict)
        assert 'PropertyTable' in props
        assert len(props['PropertyTable']['Properties']) == len(cids)
        
        # Check all CIDs are present
        found_cids = [prop['CID'] for prop in props['PropertyTable']['Properties']]
        for cid in cids:
            assert cid in found_cids
    
    def test_different_output_formats(self):
        """Test different output formats"""
        api = PubChemAPI()
        
        # Test JSON format (default)
        json_result = api.get_compound_by_cid(2244, output_format='JSON')
        assert isinstance(json_result, dict)
        
        # Test SDF format
        sdf_result = api.get_compound_by_cid(2244, output_format='SDF')
        assert isinstance(sdf_result, str)
        assert 'M  END' in sdf_result
    
    def test_cas_number_lookup(self):
        """Test looking up compounds by CAS numbers"""
        api = PubChemAPI()
        
        # Test aspirin CAS number
        try:
            result = api.get_cids_by_name('50-78-2')  # aspirin CAS
            assert isinstance(result, dict)
            assert 'IdentifierList' in result
            assert len(result['IdentifierList']['CID']) > 0
        except PubChemNotFoundError:
            # CAS lookup might not always work
            pass


class TestErrorHandling:
    """Test suite for error handling scenarios"""
    
    def test_invalid_cid(self):
        """Test with invalid CID"""
        api = PubChemAPI()
        
        with pytest.raises(PubChemNotFoundError):
            api.get_compound_by_cid(999999999)
    
    def test_invalid_smiles(self):
        """Test with invalid SMILES"""
        api = PubChemAPI()
        
        with pytest.raises((PubChemError, PubChemNotFoundError)):
            api.get_compounds_by_smiles('invalid_smiles_string')
    
    def test_empty_compound_name(self):
        """Test with empty compound name"""
        api = PubChemAPI()
        
        with pytest.raises((PubChemError, PubChemNotFoundError)):
            api.get_cids_by_name('')
    
    def test_malformed_property_names(self):
        """Test with malformed property names"""
        api = PubChemAPI()
        
        with pytest.raises(PubChemError):
            api.get_compound_properties([2244], ['InvalidPropertyName'])
    
    def test_empty_cid_list(self):
        """Test with empty CID list"""
        api = PubChemAPI()
        
        with pytest.raises((PubChemError, ValueError)):
            api.get_compound_properties([], ['MolecularWeight'])


if __name__ == "__main__":
    # Run tests if executed directly
    pytest.main([__file__, "-v"])
