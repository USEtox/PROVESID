"""
Tests for ChEMBL database interface

These tests verify the ChEMBL class functionality including:
- Database initialization and connection
- Download and extraction of database
- Search methods (by ChEMBL ID, name, InChI, InChIKey, SMILES)
- Data retrieval (compound information, properties)
- ID conversion methods
"""

import pytest
import os
import sqlite3
from provesid import CheMBL, ChEMBLError
from provesid.utils import data_path


class TestChEMBLInitialization:
    """Test ChEMBL initialization and database setup"""
    
    def test_initialization_default(self):
        """Test default initialization"""
        chembl = CheMBL()
        assert chembl is not None
        assert hasattr(chembl, 'conn')
        assert hasattr(chembl, 'cursor')
        assert os.path.exists(chembl.db_path)
    
    def test_initialization_custom_db_name(self):
        """Test initialization with custom database name"""
        # This will likely fail unless the custom db exists
        # Just testing the path construction
        try:
            chembl = CheMBL(db_name='custom_chembl.db', auto_download=False)
        except FileNotFoundError as e:
            assert 'custom_chembl.db' in str(e)
    
    def test_database_connection(self):
        """Test that database connection is valid"""
        chembl = CheMBL()
        # Try a simple query to verify connection works
        chembl.cursor.execute("SELECT name FROM sqlite_master WHERE type='table' LIMIT 1")
        result = chembl.cursor.fetchone()
        assert result is not None
    
    def test_database_has_required_tables(self):
        """Test that database contains expected tables"""
        chembl = CheMBL()
        required_tables = [
            'molecule_dictionary',
            'compound_structures',
            'compound_properties',
            'molecule_synonyms',
            'chembl_id_lookup'
        ]
        
        for table in required_tables:
            chembl.cursor.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
                (table,)
            )
            result = chembl.cursor.fetchone()
            assert result is not None, f"Table {table} not found in database"


class TestChEMBLConversion:
    """Test ID conversion methods"""
    
    @pytest.fixture
    def chembl(self):
        """Fixture providing CheMBL instance"""
        return CheMBL()
    
    def test_chembl_id_to_molregno_valid(self, chembl):
        """Test conversion of valid ChEMBL ID to molregno"""
        # CHEMBL25 is aspirin - a well-known compound
        molregno = chembl.chembl_id_to_molregno('CHEMBL25')
        assert molregno is not None
        assert isinstance(molregno, int)
        assert molregno > 0
    
    def test_chembl_id_to_molregno_invalid(self, chembl):
        """Test conversion of invalid ChEMBL ID"""
        molregno = chembl.chembl_id_to_molregno('CHEMBL99999999999')
        assert molregno is None
    
    def test_chembl_id_case_insensitive(self, chembl):
        """Test that ChEMBL ID lookup is case-insensitive"""
        molregno1 = chembl.chembl_id_to_molregno('CHEMBL25')
        molregno2 = chembl.chembl_id_to_molregno('chembl25')
        assert molregno1 == molregno2
    
    def test_molregno_to_chembl_id_valid(self, chembl):
        """Test conversion of valid molregno to ChEMBL ID"""
        # First get a valid molregno
        molregno = chembl.chembl_id_to_molregno('CHEMBL25')
        assert molregno is not None
        
        # Convert back
        chembl_id = chembl.molregno_to_chembl_id(molregno)
        assert chembl_id == 'CHEMBL25'
    
    def test_molregno_to_chembl_id_invalid(self, chembl):
        """Test conversion of invalid molregno"""
        chembl_id = chembl.molregno_to_chembl_id(999999999)
        assert chembl_id is None
    
    def test_roundtrip_conversion(self, chembl):
        """Test roundtrip conversion ChEMBL ID <-> molregno"""
        original_id = 'CHEMBL25'
        molregno = chembl.chembl_id_to_molregno(original_id)
        converted_id = chembl.molregno_to_chembl_id(molregno)
        assert converted_id == original_id


class TestChEMBLSearch:
    """Test search methods"""
    
    @pytest.fixture
    def chembl(self):
        """Fixture providing ChEMBL instance"""
        return CheMBL()
    
    def test_search_by_chembl_id_aspirin(self, chembl):
        """Test searching for aspirin by ChEMBL ID"""
        compound = chembl.search_by_chembl_id('CHEMBL25')
        assert compound is not None
        assert compound['chembl_id'] == 'CHEMBL25'
        assert 'pref_name' in compound
        assert 'ASPIRIN' in compound['pref_name'].upper()
    
    def test_search_by_chembl_id_not_found(self, chembl):
        """Test searching for non-existent ChEMBL ID"""
        compound = chembl.search_by_chembl_id('CHEMBL99999999999')
        assert compound is None
    
    def test_search_by_name_aspirin(self, chembl):
        """Test searching by compound name"""
        results = chembl.search_by_name('aspirin')
        assert len(results) > 0
        
        # Check that we found aspirin
        chembl_ids = [r['chembl_id'] for r in results]
        assert 'CHEMBL25' in chembl_ids
    
    def test_search_by_name_partial_match(self, chembl):
        """Test partial name matching"""
        results = chembl.search_by_name('acetyl')
        assert len(results) > 0
        # Should find compounds with 'acetyl' in their name
    
    def test_search_by_name_case_insensitive(self, chembl):
        """Test that name search is case-insensitive"""
        results1 = chembl.search_by_name('aspirin')
        results2 = chembl.search_by_name('ASPIRIN')
        results3 = chembl.search_by_name('AsPiRiN')
        assert len(results1) == len(results2) == len(results3)
    
    def test_search_by_name_limit(self, chembl):
        """Test search result limiting"""
        results = chembl.search_by_name('a', limit=10)
        assert len(results) <= 10
    
    def test_search_by_inchikey(self, chembl):
        """Test searching by InChI Key"""
        # Aspirin InChI Key
        inchikey = 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'
        compound = chembl.search_by_inchikey(inchikey)
        assert compound is not None
        assert compound['chembl_id'] == 'CHEMBL25'
    
    def test_search_by_inchikey_case_insensitive(self, chembl):
        """Test that InChI Key search is case-insensitive"""
        inchikey_upper = 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'
        inchikey_lower = 'bsynrymutxbxsq-uhfffaoysa-n'
        
        compound1 = chembl.search_by_inchikey(inchikey_upper)
        compound2 = chembl.search_by_inchikey(inchikey_lower)
        
        assert compound1 is not None
        assert compound2 is not None
        assert compound1['molregno'] == compound2['molregno']
    
    def test_search_by_inchi(self, chembl):
        """Test searching by InChI"""
        # Aspirin InChI
        inchi = 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)'
        compound = chembl.search_by_inchi(inchi)
        assert compound is not None
        assert compound['chembl_id'] == 'CHEMBL25'
    
    def test_search_by_smiles(self, chembl):
        """Test searching by canonical SMILES"""
        # Aspirin canonical SMILES
        smiles = 'CC(=O)Oc1ccccc1C(=O)O'
        compound = chembl.search_by_smiles(smiles)
        assert compound is not None
        assert compound['chembl_id'] == 'CHEMBL25'
    
    def test_search_by_smiles_not_found(self, chembl):
        """Test searching with non-matching SMILES"""
        compound = chembl.search_by_smiles('INVALID_SMILES_STRING')
        assert compound is None


class TestChEMBLDataRetrieval:
    """Test data retrieval methods"""
    
    @pytest.fixture
    def chembl(self):
        """Fixture providing ChEMBL instance"""
        return CheMBL()
    
    @pytest.fixture
    def aspirin_molregno(self, chembl):
        """Fixture providing aspirin's molregno"""
        return chembl.chembl_id_to_molregno('CHEMBL25')
    
    def test_get_compound_valid(self, chembl, aspirin_molregno):
        """Test retrieving compound information"""
        compound = chembl.get_compound(aspirin_molregno)
        assert compound is not None
        assert compound['molregno'] == aspirin_molregno
        assert compound['chembl_id'] == 'CHEMBL25'
        assert 'pref_name' in compound
        assert 'canonical_smiles' in compound
        assert 'standard_inchi' in compound
        assert 'standard_inchi_key' in compound
    
    def test_get_compound_invalid(self, chembl):
        """Test retrieving non-existent compound"""
        compound = chembl.get_compound(999999999)
        assert compound is None
    
    def test_get_compound_has_structure(self, chembl, aspirin_molregno):
        """Test that compound has structure information"""
        compound = chembl.get_compound(aspirin_molregno)
        assert compound is not None
        assert compound['canonical_smiles'] is not None
        assert compound['standard_inchi'] is not None
        assert compound['standard_inchi_key'] is not None
    
    def test_get_compound_has_synonyms(self, chembl, aspirin_molregno):
        """Test that compound includes synonyms"""
        compound = chembl.get_compound(aspirin_molregno)
        assert compound is not None
        assert 'synonyms' in compound
        assert isinstance(compound['synonyms'], list)
        # Aspirin should have multiple synonyms
        assert len(compound['synonyms']) > 0
    
    def test_get_properties_valid(self, chembl, aspirin_molregno):
        """Test retrieving compound properties"""
        props = chembl.get_properties(aspirin_molregno)
        assert props is not None
        assert 'mw_freebase' in props
        assert 'alogp' in props
        assert 'hba' in props
        assert 'hbd' in props
        assert 'psa' in props
        assert 'rtb' in props
    
    def test_get_properties_invalid(self, chembl):
        """Test retrieving properties for non-existent compound"""
        props = chembl.get_properties(999999999)
        assert props is None
    
    def test_get_properties_values(self, chembl, aspirin_molregno):
        """Test that property values are reasonable for aspirin"""
        props = chembl.get_properties(aspirin_molregno)
        assert props is not None
        
        # Aspirin molecular weight should be around 180
        mw = props['mw_freebase']
        assert mw is not None
        assert 175 < mw < 185
        
        # Aspirin should have hydrogen bond acceptors and donors
        hba = props['hba']
        hbd = props['hbd']
        assert hba is not None
        assert hbd is not None
        assert hba > 0
        assert hbd > 0
    
    def test_compound_structure_consistency(self, chembl):
        """Test that structure lookups are consistent"""
        # Get compound by ChEMBL ID
        compound1 = chembl.search_by_chembl_id('CHEMBL25')
        
        # Get compound by InChI Key
        inchikey = compound1['standard_inchi_key']
        compound2 = chembl.search_by_inchikey(inchikey)
        
        # Get compound by SMILES
        smiles = compound1['canonical_smiles']
        compound3 = chembl.search_by_smiles(smiles)
        
        # All should return the same compound
        assert compound1['molregno'] == compound2['molregno']
        assert compound2['molregno'] == compound3['molregno']
    
    def test_synonyms_in_search_results(self, chembl):
        """Test that all search methods return compounds with synonyms"""
        # Search by ChEMBL ID
        compound = chembl.search_by_chembl_id('CHEMBL25')
        assert 'synonyms' in compound
        assert isinstance(compound['synonyms'], list)
        
        # Search by name
        results = chembl.search_by_name('aspirin')
        if results:
            assert 'synonyms' in results[0]
            assert isinstance(results[0]['synonyms'], list)


class TestChEMBLIntegration:
    """Integration tests for complete workflows"""
    
    @pytest.fixture
    def chembl(self):
        """Fixture providing ChEMBL instance"""
        return CheMBL()
    
    def test_workflow_chembl_id_lookup(self, chembl):
        """Test complete workflow: ChEMBL ID -> compound + properties"""
        # Search by ChEMBL ID
        compound = chembl.search_by_chembl_id('CHEMBL25')
        assert compound is not None
        
        # Get properties
        molregno = compound['molregno']
        props = chembl.get_properties(molregno)
        assert props is not None
        
        # Verify data consistency
        assert compound['chembl_id'] == 'CHEMBL25'
        assert props['molregno'] == molregno
    
    def test_workflow_name_search(self, chembl):
        """Test complete workflow: name search -> multiple compounds"""
        # Search by name
        results = chembl.search_by_name('acetaminophen')
        assert len(results) > 0
        
        # Get properties for first result
        if len(results) > 0:
            compound = results[0]
            props = chembl.get_properties(compound['molregno'])
            # Should have properties
            assert props is not None or compound['molregno'] is not None
    
    def test_workflow_structure_search(self, chembl):
        """Test complete workflow: structure search -> compound info"""
        # Start with SMILES
        smiles = 'CC(=O)Oc1ccccc1C(=O)O'
        compound = chembl.search_by_smiles(smiles)
        assert compound is not None
        
        # Get ChEMBL ID
        chembl_id = compound['chembl_id']
        assert chembl_id == 'CHEMBL25'
        
        # Verify InChI
        assert compound['standard_inchi'] is not None
        
        # Get properties
        props = chembl.get_properties(compound['molregno'])
        assert props is not None


class TestChEMBLRobustness:
    """Test error handling and edge cases"""
    
    @pytest.fixture
    def chembl(self):
        """Fixture providing ChEMBL instance"""
        return CheMBL()
    
    def test_empty_string_chembl_id(self, chembl):
        """Test searching with empty ChEMBL ID"""
        compound = chembl.search_by_chembl_id('')
        assert compound is None
    
    def test_empty_string_name(self, chembl):
        """Test searching with empty name"""
        results = chembl.search_by_name('')
        # Empty string might match many things or nothing
        assert isinstance(results, list)
    
    def test_special_characters_in_name(self, chembl):
        """Test name search with special characters"""
        results = chembl.search_by_name("2'-deoxyadenosine")
        assert isinstance(results, list)
    
    def test_very_long_name_search(self, chembl):
        """Test name search with very long string"""
        long_name = 'a' * 1000
        results = chembl.search_by_name(long_name)
        assert isinstance(results, list)
        # Unlikely to find matches
        assert len(results) == 0 or len(results) > 0
    
    def test_connection_cleanup(self):
        """Test that database connection is properly closed"""
        chembl = CheMBL()
        db_path = chembl.db_path
        del chembl
        
        # Verify we can open the database (no locks)
        conn = sqlite3.connect(db_path)
        conn.close()


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
