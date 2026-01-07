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
