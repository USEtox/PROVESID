"""
Tests for ChebiSDF class

Tests the ChEBI SDF file parser and query methods.
"""

import pytest
import sys
import os

# Add the src directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from provesid.chebi import ChebiSDF


class TestChebiSDFInitialization:
    """Test suite for ChebiSDF initialization"""
    
    @pytest.fixture
    def chebi_sdf(self):
        """Create a ChebiSDF instance for testing"""
        return ChebiSDF()
    
    def test_initialization_default(self, chebi_sdf):
        """Test ChebiSDF initialization with default SDF file"""
        assert chebi_sdf.sdf_path is not None
        assert os.path.exists(chebi_sdf.sdf_path)
        assert chebi_sdf.index is not None
        assert len(chebi_sdf.index['id_to_offset']) > 0
    
    def test_index_structure(self, chebi_sdf):
        """Test that index has the expected structure"""
        required_keys = ['id_to_offset', 'name_to_ids', 'inchikey_to_id', 
                        'inchi_to_id', 'formula_to_ids', 'cas_to_ids', 'synonym_to_ids']
        for key in required_keys:
            assert key in chebi_sdf.index
    
    def test_index_cache_created(self, chebi_sdf):
        """Test that index cache file is created"""
        assert os.path.exists(chebi_sdf.index_path)
    
    def test_initialization_with_auto_download_false(self):
        """Test that FileNotFoundError is raised when auto_download=False and file missing"""
        # This test assumes the file exists, so we test with a fake path
        with pytest.raises(FileNotFoundError):
            ChebiSDF(sdf_path="/nonexistent/path/chebi.sdf", auto_download=False)


class TestChebiSDFQueries:
    """Test suite for ChebiSDF query methods"""
    
    @pytest.fixture
    def chebi_sdf(self):
        """Create a ChebiSDF instance for testing"""
        return ChebiSDF()
    
    def test_get_compound_by_id(self, chebi_sdf):
        """Test getting compound by ChEBI ID"""
        # Water: CHEBI:15377
        water = chebi_sdf.get_compound_by_id("CHEBI:15377")
        
        assert water is not None
        assert 'ChEBI ID' in water
        assert water['ChEBI ID'] == 'CHEBI:15377'
        assert 'ChEBI NAME' in water
        assert water['ChEBI NAME'].lower() == 'water'
    
    def test_get_compound_by_id_without_prefix(self, chebi_sdf):
        """Test getting compound by ChEBI ID without CHEBI: prefix"""
        water = chebi_sdf.get_compound_by_id("15377")
        
        assert water is not None
        assert water['ChEBI ID'] == 'CHEBI:15377'
    
    def test_get_compound_by_id_nonexistent(self, chebi_sdf):
        """Test getting nonexistent compound"""
        result = chebi_sdf.get_compound_by_id("CHEBI:9999999")
        assert result is None
    
    def test_search_by_name_exact(self, chebi_sdf):
        """Test exact name search"""
        results = chebi_sdf.search_by_name("water", exact=True)
        
        assert len(results) > 0
        assert any(r['ChEBI NAME'].lower() == 'water' for r in results)
    
    def test_search_by_name_partial(self, chebi_sdf):
        """Test partial name search"""
        results = chebi_sdf.search_by_name("glucose", exact=False)
        
        assert len(results) > 0
        # Should find compounds with "glucose" in the name
        assert any('glucose' in r['ChEBI NAME'].lower() for r in results)
    
    def test_search_by_inchikey(self, chebi_sdf):
        """Test search by InChIKey"""
        # Water InChIKey
        inchikey = "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
        result = chebi_sdf.search_by_inchikey(inchikey)
        
        assert result is not None
        assert result['INCHIKEY'] == inchikey
    
    def test_search_by_cas(self, chebi_sdf):
        """Test search by CAS number"""
        # Aspirin CAS: 50-78-2
        results = chebi_sdf.search_by_cas("50-78-2")
        
        assert len(results) > 0
        # Check that CAS is in the results
        for result in results:
            if 'CAS Registry Numbers' in result:
                assert '50-78-2' in result['CAS Registry Numbers']
    
    def test_search_by_formula(self, chebi_sdf):
        """Test search by molecular formula"""
        # Search for water: H2O
        results = chebi_sdf.search_by_formula("H2O")
        
        assert len(results) > 0
        assert any(r.get('FORMULA') == 'H2O' for r in results)
    
    def test_search_by_synonym(self, chebi_sdf):
        """Test search by synonym"""
        results = chebi_sdf.search_by_synonym("acetylsalicylic acid", exact=False)
        
        # Should find aspirin
        assert len(results) > 0


class TestChebiSDFBatchOperations:
    """Test suite for batch operations"""
    
    @pytest.fixture
    def chebi_sdf(self):
        """Create a ChebiSDF instance for testing"""
        return ChebiSDF()
    
    def test_get_compounds_by_ids(self, chebi_sdf):
        """Test getting multiple compounds at once"""
        chebi_ids = ["CHEBI:15377", "CHEBI:16236"]  # water, ethanol
        results = chebi_sdf.get_compounds_by_ids(chebi_ids)
        
        assert len(results) == 2
        assert all('ChEBI ID' in r for r in results)
    
    def test_export_to_dataframe(self, chebi_sdf):
        """Test exporting to pandas DataFrame"""
        import pandas as pd
        
        chebi_ids = ["CHEBI:15377", "CHEBI:16236"]
        df = chebi_sdf.export_to_dataframe(chebi_ids)
        
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 2
        assert 'ChEBI ID' in df.columns
        assert 'ChEBI NAME' in df.columns
    
    def test_get_database_stats(self, chebi_sdf):
        """Test getting database statistics"""
        stats = chebi_sdf.get_database_stats()
        
        assert 'total_compounds' in stats
        assert stats['total_compounds'] > 0
        assert 'compounds_with_inchikey' in stats
        assert 'compounds_with_cas' in stats
