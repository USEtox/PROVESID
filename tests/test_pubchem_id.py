"""
Tests for PubChemID class - local SQLite database identifier lookup and conversion.

These tests require the PubChem ID database to be built first using:
    python scripts/build_pubchem_id_db.py
"""

import pytest
import os
import pandas as pd
from provesid.utils import data_path


# Check if database exists
db_path = os.path.join(data_path(), 'pubchem_id.db')
DB_EXISTS = os.path.exists(db_path)

pytestmark = pytest.mark.skipif(
    not DB_EXISTS,
    reason="PubChem ID database not found. Run scripts/build_pubchem_id_db.py first."
)


@pytest.fixture
def pubchem_id():
    """Fixture to create PubChemID instance."""
    from provesid import PubChemID
    return PubChemID()


class TestPubChemIDInitialization:
    """Test database initialization and connection."""
    
    def test_init_default_path(self, pubchem_id):
        """Test initialization with default database path."""
        assert pubchem_id is not None
        assert pubchem_id.conn is not None
        assert os.path.exists(pubchem_id.db_path)
    
    def test_init_custom_path(self):
        """Test initialization with custom database path."""
        from provesid import PubChemID
        db = PubChemID(db_path=db_path)
        assert db.db_path == db_path
    
    def test_init_nonexistent_path(self):
        """Test initialization with non-existent database path."""
        from provesid import PubChemID
        with pytest.raises(FileNotFoundError):
            PubChemID(db_path='nonexistent.db')
    
    def test_get_stats(self, pubchem_id):
        """Test database statistics."""
        stats = pubchem_id.get_stats()
        assert 'total_compounds' in stats
        assert 'total_cas_numbers' in stats
        assert 'compounds_with_cas' in stats
        assert 'total_synonyms' in stats
        assert stats['total_compounds'] > 0
        assert stats['total_cas_numbers'] > 0


class TestPubChemIDLookups:
    """Test basic lookup methods."""
    
    def test_get_by_cas_aspirin(self, pubchem_id):
        """Test lookup by CAS number (Aspirin)."""
        result = pubchem_id.get_by_cas("50-78-2")
        if result:  # May not be in every database version
            assert result['cid'] is not None
            assert '50-78-2' in result['cas_numbers']
            assert result['inchikey'] is not None
    
    def test_get_by_cas_formaldehyde(self, pubchem_id):
        """Test lookup by CAS number (Formaldehyde)."""
        result = pubchem_id.get_by_cas("50-00-0")
        if result:
            assert result['cid'] is not None
            assert '50-00-0' in result['cas_numbers']
    
    def test_get_by_cas_nonexistent(self, pubchem_id):
        """Test lookup with non-existent CAS number."""
        result = pubchem_id.get_by_cas("99999-99-9")
        assert result is None
    
    def test_get_by_cid(self, pubchem_id):
        """Test lookup by PubChem CID."""
        # Try a common CID
        result = pubchem_id.get_by_cid(2244)  # Aspirin
        if result:
            assert result['cid'] == 2244
            assert result['inchikey'] is not None
    
    def test_get_by_inchikey(self, pubchem_id):
        """Test lookup by InChIKey."""
        # Aspirin InChIKey
        result = pubchem_id.get_by_inchikey("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        if result:
            assert result['inchikey'] == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
            assert result['cid'] is not None
    
    def test_search_by_name_exact(self, pubchem_id):
        """Test exact name search."""
        results = pubchem_id.search_by_name("Aspirin", exact=True)
        # Results may vary, just check structure
        assert isinstance(results, list)
    
    def test_search_by_name_partial(self, pubchem_id):
        """Test partial name search."""
        results = pubchem_id.search_by_name("aspirin", exact=False, limit=5)
        assert isinstance(results, list)
        assert len(results) <= 5
    
    def test_search_by_formula(self, pubchem_id):
        """Test search by molecular formula."""
        results = pubchem_id.search_by_formula("C9H8O4", limit=10)
        assert isinstance(results, list)
        assert len(results) <= 10


class TestPubChemIDConversions:
    """Test identifier conversion methods."""
    
    def test_cas_to_cid(self, pubchem_id):
        """Test CAS to CID conversion."""
        cid = pubchem_id.cas_to_cid("50-78-2")
        if cid:
            assert isinstance(cid, int)
            assert cid > 0
    
    def test_cas_to_inchi(self, pubchem_id):
        """Test CAS to InChI conversion."""
        inchi = pubchem_id.cas_to_inchi("50-78-2")
        if inchi:
            assert isinstance(inchi, str)
            assert inchi.startswith("InChI=")
    
    def test_cas_to_inchikey(self, pubchem_id):
        """Test CAS to InChIKey conversion."""
        inchikey = pubchem_id.cas_to_inchikey("50-78-2")
        if inchikey:
            assert isinstance(inchikey, str)
            assert len(inchikey) == 27
    
    def test_cas_to_smiles(self, pubchem_id):
        """Test CAS to SMILES conversion."""
        smiles = pubchem_id.cas_to_smiles("50-78-2")
        if smiles:
            assert isinstance(smiles, str)
    
    def test_inchikey_to_cid(self, pubchem_id):
        """Test InChIKey to CID conversion."""
        cid = pubchem_id.inchikey_to_cid("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        if cid:
            assert isinstance(cid, int)
    
    def test_cid_to_cas(self, pubchem_id):
        """Test CID to CAS conversion."""
        cas_list = pubchem_id.cid_to_cas(2244)  # Aspirin
        if cas_list:
            assert isinstance(cas_list, list)
            assert len(cas_list) > 0
    
    def test_conversion_nonexistent(self, pubchem_id):
        """Test conversion with non-existent identifier."""
        result = pubchem_id.cas_to_cid("99999-99-9")
        assert result is None


class TestPubChemIDBatchOperations:
    """Test batch conversion methods."""
    
    def test_batch_cas_to_cid(self, pubchem_id):
        """Test batch CAS to CID conversion."""
        cas_list = ["50-78-2", "50-00-0", "99999-99-9"]
        results = pubchem_id.batch_cas_to_cid(cas_list)
        
        assert isinstance(results, dict)
        assert len(results) == 3
        assert "50-78-2" in results
        assert "99999-99-9" in results
        assert results["99999-99-9"] is None  # Non-existent
    
    def test_batch_cas_to_inchikey(self, pubchem_id):
        """Test batch CAS to InChIKey conversion."""
        cas_list = ["50-78-2", "50-00-0"]
        results = pubchem_id.batch_cas_to_inchikey(cas_list)
        
        assert isinstance(results, dict)
        assert len(results) == 2
    
    def test_batch_cid_to_cas(self, pubchem_id):
        """Test batch CID to CAS conversion."""
        cid_list = [2244, 712]  # Aspirin, Formaldehyde
        results = pubchem_id.batch_cid_to_cas(cid_list)
        
        assert isinstance(results, dict)
        assert len(results) == 2
    
    def test_get_id_table_from_cas(self, pubchem_id):
        """Test getting identifier table for single CAS."""
        df = pubchem_id.get_id_table_from_cas("50-78-2")
        
        if df is not None:
            assert 'cid' in df.columns
            assert 'cas' in df.columns
            assert 'inchi' in df.columns
            assert 'inchikey' in df.columns
            assert len(df) == 1
            assert df['cas'].iloc[0] == "50-78-2"
    
    def test_get_id_table_from_cas_nonexistent(self, pubchem_id):
        """Test getting identifier table for non-existent CAS."""
        df = pubchem_id.get_id_table_from_cas("99999-99-9")
        assert df is None
    
    def test_batch_get_id_table_from_cas(self, pubchem_id):
        """Test getting identifier table for multiple CAS numbers."""
        cas_list = ["50-78-2", "50-00-0"]
        df = pubchem_id.batch_get_id_table_from_cas(cas_list)
        
        assert 'cid' in df.columns
        assert 'cas' in df.columns
        assert 'inchi' in df.columns
        # Should have results for valid CAS numbers
        if len(df) > 0:
            assert df['cas'].isin(cas_list).any()
    
    def test_batch_get_id_table_empty_list(self, pubchem_id):
        """Test batch get with empty list."""
        df = pubchem_id.batch_get_id_table_from_cas([])
        
        assert len(df) == 0
        assert 'cid' in df.columns
        assert 'cas' in df.columns
    
    def test_batch_get_id_table_nonexistent_only(self, pubchem_id):
        """Test batch get with only non-existent CAS numbers."""
        df = pubchem_id.batch_get_id_table_from_cas(["99999-99-9", "88888-88-8"])
        
        assert len(df) == 0
        assert 'cid' in df.columns
    
    def test_get_by_cas_batch(self, pubchem_id):
        """Test get_by_cas_batch with multiple CAS numbers."""
        cas_list = ["50-78-2", "50-00-0", "64-17-5"]  # Aspirin, Formaldehyde, Ethanol
        df = pubchem_id.get_by_cas_batch(cas_list)
        
        # Check DataFrame structure
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0
        
        # Check all expected columns exist
        expected_cols = ['cid', 'cas', 'inchi', 'inchikey', 'smiles', 'cmpdname', 
                        'iupacname', 'mf', 'mw', 'polararea', 'complexity', 'xlogp',
                        'heavycnt', 'hbonddonor', 'hbondacc', 'rotbonds', 'exactmass',
                        'charge', 'cidcdate']
        for col in expected_cols:
            assert col in df.columns
        
        # Check that CAS numbers are in the result
        assert set(cas_list).issubset(set(df['cas'].tolist()))
    
    def test_get_by_cas_batch_empty_list(self, pubchem_id):
        """Test get_by_cas_batch with empty list."""
        df = pubchem_id.get_by_cas_batch([])
        
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0
        assert 'cid' in df.columns
        assert 'cas' in df.columns
    
    def test_get_by_cas_batch_nonexistent(self, pubchem_id):
        """Test get_by_cas_batch with non-existent CAS numbers."""
        df = pubchem_id.get_by_cas_batch(["99999-99-9", "88888-88-8"])
        
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0
    
    def test_get_by_cas_batch_mixed(self, pubchem_id):
        """Test get_by_cas_batch with mix of valid and invalid CAS numbers."""
        cas_list = ["50-78-2", "99999-99-9", "50-00-0"]
        df = pubchem_id.get_by_cas_batch(cas_list)
        
        # Should only have valid CAS numbers
        assert len(df) == 2
        assert "50-78-2" in df['cas'].tolist()
        assert "50-00-0" in df['cas'].tolist()
        assert "99999-99-9" not in df['cas'].tolist()


class TestPubChemIDCASConversionMethods:
    """Test new CAS conversion methods (smiles_to_cas, name_to_cas, formula_to_cas)."""
    
    def test_smiles_to_cas(self, pubchem_id):
        """Test converting SMILES to CAS."""
        # Test with simple molecules
        # Methane
        cas_list = pubchem_id.smiles_to_cas("C")
        if cas_list:
            assert isinstance(cas_list, list)
            assert len(cas_list) > 0
            assert all(isinstance(cas, str) for cas in cas_list)
    
    def test_smiles_to_cas_aspirin(self, pubchem_id):
        """Test converting aspirin SMILES to CAS."""
        # Aspirin SMILES
        cas_list = pubchem_id.smiles_to_cas("CC(=O)OC1=CC=CC=C1C(=O)O")
        if cas_list:
            assert isinstance(cas_list, list)
            # Should include aspirin CAS
            assert "50-78-2" in cas_list
    
    def test_smiles_to_cas_invalid(self, pubchem_id):
        """Test converting invalid SMILES."""
        cas_list = pubchem_id.smiles_to_cas("INVALID_SMILES_123")
        assert cas_list is None
    
    def test_name_to_cas_exact(self, pubchem_id):
        """Test converting chemical name to CAS (exact match)."""
        # Search for aspirin (case sensitive exact match may not work)
        # Try with proper capitalization
        cas_list = pubchem_id.name_to_cas("Aspirin", exact=False)
        if cas_list:
            assert isinstance(cas_list, list)
            assert len(cas_list) > 0
            # Should include aspirin CAS
            assert "50-78-2" in cas_list
    
    def test_name_to_cas_nonexistent(self, pubchem_id):
        """Test converting non-existent name to CAS."""
        cas_list = pubchem_id.name_to_cas("NonexistentChemical12345XYZ")
        assert cas_list is None
    
    def test_formula_to_cas(self, pubchem_id):
        """Test converting molecular formula to CAS."""
        # Water
        cas_list = pubchem_id.formula_to_cas("H2O")
        if cas_list:
            assert isinstance(cas_list, list)
            assert len(cas_list) > 0
            # Should be sorted
            assert cas_list == sorted(cas_list)
    
    def test_formula_to_cas_common_formula(self, pubchem_id):
        """Test converting common formula (should return multiple)."""
        # Aspirin formula
        cas_list = pubchem_id.formula_to_cas("C9H8O4")
        if cas_list:
            assert isinstance(cas_list, list)
            assert len(cas_list) > 0
            # Should include aspirin
            assert "50-78-2" in cas_list
    
    def test_formula_to_cas_nonexistent(self, pubchem_id):
        """Test converting non-existent formula."""
        cas_list = pubchem_id.formula_to_cas("Zr999Xe999")
        assert cas_list is None
    
    def test_formula_to_cas_limit(self, pubchem_id):
        """Test formula conversion with limit parameter."""
        # Use a common formula
        cas_list = pubchem_id.formula_to_cas("CH4O", limit=5)
        if cas_list:
            # Due to aggregation, may have more than 5 CAS numbers
            # but we requested at most 5 compounds
            assert isinstance(cas_list, list)
    
    def test_batch_smiles_to_cas(self, pubchem_id):
        """Test batch SMILES to CAS conversion."""
        smiles_list = ["C", "CO", "CCO"]
        results = pubchem_id.batch_smiles_to_cas(smiles_list)
        
        assert isinstance(results, dict)
        assert len(results) == 3
        
        # Check all SMILES are in results
        for smiles in smiles_list:
            assert smiles in results
            # Results can be None or list
            assert results[smiles] is None or isinstance(results[smiles], list)
    
    def test_batch_smiles_to_cas_empty(self, pubchem_id):
        """Test batch SMILES conversion with empty list."""
        results = pubchem_id.batch_smiles_to_cas([])
        assert isinstance(results, dict)
        assert len(results) == 0
    
    def test_batch_name_to_cas(self, pubchem_id):
        """Test batch name to CAS conversion."""
        names = ["aspirin", "caffeine", "glucose"]
        results = pubchem_id.batch_name_to_cas(names, exact=False)
        
        assert isinstance(results, dict)
        assert len(results) == 3
        
        # Check all names are in results
        for name in names:
            assert name in results
            # Results can be None or list
            assert results[name] is None or isinstance(results[name], list)
    
    def test_batch_name_to_cas_empty(self, pubchem_id):
        """Test batch name conversion with empty list."""
        results = pubchem_id.batch_name_to_cas([])
        assert isinstance(results, dict)
        assert len(results) == 0
    
    def test_batch_formula_to_cas(self, pubchem_id):
        """Test batch formula to CAS conversion."""
        formulas = ["H2O", "CH4", "C9H8O4"]
        results = pubchem_id.batch_formula_to_cas(formulas)
        
        assert isinstance(results, dict)
        assert len(results) == 3
        
        # Check all formulas are in results
        for formula in formulas:
            assert formula in results
            # Results can be None or list
            assert results[formula] is None or isinstance(results[formula], list)
    
    def test_batch_formula_to_cas_empty(self, pubchem_id):
        """Test batch formula conversion with empty list."""
        results = pubchem_id.batch_formula_to_cas([])
        assert isinstance(results, dict)
        assert len(results) == 0
    
    def test_smiles_to_cas_integration(self, pubchem_id):
        """Integration test: SMILES -> CAS -> back to SMILES."""
        # Start with a SMILES
        original_smiles = "CCO"  # Ethanol
        
        # Convert to CAS
        cas_list = pubchem_id.smiles_to_cas(original_smiles)
        if cas_list:
            # Convert first CAS back to SMILES
            smiles_result = pubchem_id.cas_to_smiles(cas_list[0])
            
            # SMILES might not be identical due to canonicalization
            # but should be valid
            if smiles_result:
                assert isinstance(smiles_result, str)
                assert len(smiles_result) > 0
