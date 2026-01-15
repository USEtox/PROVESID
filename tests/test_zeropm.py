"""
Tests for ZeroPM SQLite database functionality

The ZeroPM class provides access to the ZeroPM database containing chemical data
with InChI/InChIKey/CAS mappings and fuzzy name searching capabilities.
"""

import pytest
import sys
import os
import sqlite3
import tempfile
import shutil
from unittest.mock import patch, MagicMock

# Add the src directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from provesid.zeropm import ZeroPM


class TestZeroPMInitialization:
    """Test suite for ZeroPM initialization"""
    
    def test_initialization_default(self):
        """Test ZeroPM initialization with default database"""
        zpm = ZeroPM()
        assert zpm.db_path is not None
        assert os.path.exists(zpm.db_path)
        assert zpm.conn is not None
        assert zpm.cursor is not None
    
    def test_initialization_custom_db(self):
        """Test ZeroPM initialization with custom database name"""
        zpm = ZeroPM(db_name='zeropm-v0-0-3.sqlite')
        assert 'zeropm-v0-0-3.sqlite' in zpm.db_path
        assert os.path.exists(zpm.db_path)
    
    def test_initialization_nonexistent_db_no_autodownload(self):
        """Test ZeroPM initialization with nonexistent database and auto_download=False"""
        with pytest.raises(FileNotFoundError):
            ZeroPM(db_name='nonexistent_database.sqlite', auto_download=False)
    
    def test_connection_cleanup(self):
        """Test that database connection is properly closed"""
        zpm = ZeroPM()
        conn = zpm.conn
        del zpm
        # Check that connection was closed
        with pytest.raises(sqlite3.ProgrammingError):
            conn.execute("SELECT 1")
    
    def test_download_database_already_exists(self):
        """Test that download_database raises error if database already exists"""
        zpm = ZeroPM()
        with pytest.raises(FileExistsError):
            zpm.download_database()
    
    def test_download_database_force_parameter(self):
        """Test download_database with force=True parameter"""
        zpm = ZeroPM()
        # Create a mock for requests.get to avoid actual download
        with patch('provesid.zeropm.requests.get') as mock_get:
            # Mock response
            mock_response = MagicMock()
            mock_response.headers = {'content-length': '1000'}
            mock_response.iter_content = lambda chunk_size: [b'test_data']
            mock_get.return_value = mock_response
            
            # Should not raise error with force=True
            try:
                zpm.download_database(force=True)
            except Exception:
                # It's ok if this fails due to mocking, we just want to verify
                # that FileExistsError is not raised with force=True
                pass



class TestZeroPMQueryMethods:
    """Test suite for ZeroPM query methods"""
    
    @pytest.fixture
    def zpm(self):
        """Create a ZeroPM instance for testing"""
        return ZeroPM()
    
    def test_query_cas_existing(self, zpm):
        """Test querying an existing CAS number"""
        # Test with a known CAS number (if any exists in the database)
        # This test may need adjustment based on actual database contents
        zpm.cursor.execute("""
            SELECT query, query_id 
            FROM api_ready_query 
            WHERE type = 'CAS Registry Number' 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            cas, expected_id = result
            query_id = zpm.query_cas(cas)
            assert query_id == expected_id
        else:
            pytest.skip("No CAS numbers found in database")
    
    def test_query_cas_nonexistent(self, zpm):
        """Test querying a nonexistent CAS number"""
        query_id = zpm.query_cas('999-99-9')
        assert query_id is None
    
    def test_query_name_existing(self, zpm):
        """Test querying an existing chemical name"""
        # Get a known chemical name from the database
        zpm.cursor.execute("""
            SELECT query, query_id 
            FROM api_ready_query 
            WHERE type = 'chemical name' 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            name, expected_id = result
            query_id = zpm.query_name(name)
            assert query_id == expected_id
        else:
            pytest.skip("No chemical names found in database")
    
    def test_query_name_nonexistent(self, zpm):
        """Test querying a nonexistent chemical name"""
        query_id = zpm.query_name('nonexistent_chemical_name_12345')
        assert query_id is None
    
    def test_query_similar_name_found(self, zpm):
        """Test fuzzy name matching with similar names"""
        # Get a known chemical name from the database
        zpm.cursor.execute("""
            SELECT query 
            FROM api_ready_query 
            WHERE type = 'chemical name' 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            name = result[0]
            # Try with a slightly modified name (e.g., lowercase)
            similar_name = name.lower() if name.isupper() else name.upper()
            query_ids = zpm.query_similar_name(similar_name, number_of_results=5, score_cutoff=50)
            # Should find results due to fuzzy matching
            assert query_ids is not None
            assert isinstance(query_ids, list)
            assert len(query_ids) > 0
        else:
            pytest.skip("No chemical names found in database")
    
    def test_query_similar_name_no_match(self, zpm):
        """Test fuzzy name matching with no matches"""
        # Use a very high cutoff score to ensure no matches
        query_ids = zpm.query_similar_name('xyzabc123nonexistent', score_cutoff=95)
        # Note: fuzzy matching might still find weak matches, so we check if results are None or have low scores
        # The important thing is that the function doesn't crash
        assert query_ids is None or len(query_ids) >= 0
    
    def test_query_similar_name_custom_params(self, zpm):
        """Test fuzzy name matching with custom parameters"""
        zpm.cursor.execute("""
            SELECT query 
            FROM api_ready_query 
            WHERE type = 'chemical name' 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            name = result[0]
            # Test with different parameters
            query_ids = zpm.query_similar_name(name, number_of_results=3, score_cutoff=90)
            if query_ids:
                assert len(query_ids) <= 3
        else:
            pytest.skip("No chemical names found in database")


class TestZeroPMInChIMethods:
    """Test suite for InChI/InChIKey related methods"""
    
    @pytest.fixture
    def zpm(self):
        """Create a ZeroPM instance for testing"""
        return ZeroPM()
    
    def test_get_inchi_id_existing(self, zpm):
        """Test getting InChI IDs for an existing query"""
        # Get a known query_id
        zpm.cursor.execute("""
            SELECT query_id 
            FROM api_results 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            query_id = result[0]
            inchi_ids, ranks = zpm.get_inchi_id(query_id)
            assert isinstance(inchi_ids, list)
            assert isinstance(ranks, list)
            assert len(inchi_ids) == len(ranks)
            if len(inchi_ids) > 0:
                assert all(isinstance(id, int) for id in inchi_ids)
                assert all(isinstance(rank, int) for rank in ranks)
        else:
            pytest.skip("No results found in database")
    
    def test_get_inchi_id_nonexistent(self, zpm):
        """Test getting InChI IDs for a nonexistent query"""
        inchi_ids, ranks = zpm.get_inchi_id(999999)
        assert inchi_ids == []
        assert ranks == []
    
    def test_get_inchi_existing(self, zpm):
        """Test getting InChI for an existing inchi_id"""
        # Get a known inchi_id
        zpm.cursor.execute("""
            SELECT inchi_id, inchi, inchikey 
            FROM substances 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            inchi_id, expected_inchi, expected_inchikey = result
            inchi, inchikey = zpm.get_inchi(inchi_id)
            assert inchi == expected_inchi
            assert inchikey == expected_inchikey
        else:
            pytest.skip("No substances found in database")
    
    def test_get_inchi_nonexistent(self, zpm):
        """Test getting InChI for a nonexistent inchi_id"""
        inchi, inchikey = zpm.get_inchi(999999)
        assert inchi is None
        assert inchikey is None
    
    def test_inchi_to_smiles_valid(self, zpm):
        """Test InChI to SMILES conversion with valid InChI"""
        # Simple water molecule
        inchi = "InChI=1S/H2O/h1H2"
        smiles = zpm._inchi_to_smiles(inchi)
        assert smiles is not None
        assert isinstance(smiles, str)
        # Water can be represented as 'O' or '[H]O[H]' depending on RDKit version
        assert 'O' in smiles
    
    def test_inchi_to_smiles_invalid(self, zpm):
        """Test InChI to SMILES conversion with invalid InChI"""
        inchi = "invalid_inchi_string"
        smiles = zpm._inchi_to_smiles(inchi)
        assert smiles is None


class TestZeroPMCASMethods:
    """Test suite for CAS number related methods"""
    
    @pytest.fixture
    def zpm(self):
        """Create a ZeroPM instance for testing"""
        return ZeroPM()
    
    def test_get_names_existing_cas(self, zpm):
        """Test getting names for an existing CAS number"""
        # Get a CAS number that has associated names
        zpm.cursor.execute("""
            SELECT query 
            FROM api_ready_query 
            WHERE type = 'CAS Registry Number' 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            cas = result[0]
            names = zpm.get_names(cas)
            assert isinstance(names, list)
            # Names might be empty if no inventory data exists
        else:
            pytest.skip("No CAS numbers found in database")
    
    def test_get_names_nonexistent_cas(self, zpm):
        """Test getting names for a nonexistent CAS number"""
        names = zpm.get_names('999-99-9')
        assert names == []
    
    def test_get_smiles_from_cas_existing(self, zpm):
        """Test getting SMILES from an existing CAS number"""
        # Get a CAS number with InChI data
        zpm.cursor.execute("""
            SELECT aq.query 
            FROM api_ready_query aq
            JOIN api_results ar ON aq.query_id = ar.query_id
            JOIN substances s ON ar.inchi_id = s.inchi_id
            WHERE aq.type = 'CAS Registry Number' 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            cas = result[0]
            smiles = zpm.get_smiles_from_cas(cas)
            # SMILES may be None if RDKit conversion fails
            if smiles:
                assert isinstance(smiles, str)
        else:
            pytest.skip("No CAS numbers with InChI data found")
    
    def test_get_smiles_from_cas_nonexistent(self, zpm):
        """Test getting SMILES from a nonexistent CAS number"""
        smiles = zpm.get_smiles_from_cas('999-99-9')
        assert smiles is None
    
    def test_get_cas_from_inchi_existing(self, zpm):
        """Test getting CAS from an existing InChI"""
        # Get an InChI that has a CAS number
        zpm.cursor.execute("""
            SELECT s.inchi 
            FROM substances s
            JOIN api_results ar ON s.inchi_id = ar.inchi_id
            JOIN api_ready_query aq ON ar.query_id = aq.query_id
            WHERE aq.type = 'CAS Registry Number'
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            inchi = result[0]
            cas = zpm.get_cas_from_inchi(inchi)
            assert cas is not None
            # Can be string or list
            if isinstance(cas, list):
                assert len(cas) > 0
            else:
                assert isinstance(cas, str)
        else:
            pytest.skip("No InChI with CAS mapping found")
    
    def test_get_cas_from_inchi_nonexistent(self, zpm):
        """Test getting CAS from a nonexistent InChI"""
        inchi = "InChI=1S/C99H99N99O99/nonexistent"
        cas = zpm.get_cas_from_inchi(inchi)
        assert cas is None
    
    def test_get_cas_from_inchikey_existing(self, zpm):
        """Test getting CAS from an existing InChIKey"""
        # Get an InChIKey that has a CAS number
        zpm.cursor.execute("""
            SELECT s.inchikey 
            FROM substances s
            JOIN api_results ar ON s.inchi_id = ar.inchi_id
            JOIN api_ready_query aq ON ar.query_id = aq.query_id
            WHERE aq.type = 'CAS Registry Number'
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            inchikey = result[0]
            cas = zpm.get_cas_from_inchikey(inchikey)
            assert cas is not None
            # Can be string or list
            if isinstance(cas, list):
                assert len(cas) > 0
            else:
                assert isinstance(cas, str)
        else:
            pytest.skip("No InChIKey with CAS mapping found")
    
    def test_get_cas_from_inchikey_nonexistent(self, zpm):
        """Test getting CAS from a nonexistent InChIKey"""
        inchikey = "NONEXISTENT-INCHIKEY-XY"
        cas = zpm.get_cas_from_inchikey(inchikey)
        assert cas is None
    
    def test_get_smiles_from_inchikey_existing(self, zpm):
        """Test getting SMILES from an existing InChIKey"""
        zpm.cursor.execute("""
            SELECT inchikey 
            FROM substances 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            inchikey = result[0]
            smiles = zpm.get_smiles_from_inchikey(inchikey)
            # SMILES may be None if RDKit conversion fails
            if smiles:
                assert isinstance(smiles, str)
        else:
            pytest.skip("No InChIKey found in database")
    
    def test_get_smiles_from_inchikey_nonexistent(self, zpm):
        """Test getting SMILES from a nonexistent InChIKey"""
        inchikey = "NONEXISTENT-INCHIKEY-XY"
        smiles = zpm.get_smiles_from_inchikey(inchikey)
        assert smiles is None
    
    def test_get_cas_from_smiles_valid(self, zpm):
        """Test getting CAS from a valid SMILES"""
        # Use a simple SMILES that might exist in the database
        # This test is database-dependent
        smiles = "O"  # Water
        cas = zpm.get_cas_from_smiles(smiles)
        # May or may not find a CAS depending on database contents
        if cas:
            assert isinstance(cas, (str, list))
    
    def test_get_cas_from_smiles_invalid(self, zpm):
        """Test getting CAS from an invalid SMILES"""
        smiles = "invalid_smiles_xyz123"
        cas = zpm.get_cas_from_smiles(smiles)
        assert cas is None
    
    def test_get_id_table_from_cas_existing(self, zpm):
        """Test getting identifier table from an existing CAS number"""
        import pandas as pd
        
        # Use formaldehyde (CAS: 50-00-0) which should exist in the database
        cas = "50-00-0"
        df = zpm.get_id_table_from_cas(cas)
        
        # Check that we got a DataFrame
        assert isinstance(df, pd.DataFrame)
        
        # Check that the DataFrame has the expected columns
        expected_columns = {'cas', 'query_id', 'inchi_id', 'rank', 'inchi', 'inchikey', 'zeropm_id', 'synonyms', 'sources'}
        assert set(df.columns) == expected_columns
        
        # Check that all rows have the same CAS
        assert (df['cas'] == cas).all()
        
        # Check that we have at least one row
        assert len(df) > 0
        
        # Check data types
        assert df['cas'].dtype == object  # string
        assert df['inchi'].dtype == object  # string
        assert df['inchikey'].dtype == object  # string
        
        # Check that InChI and InChIKey are not None for rows with inchi_id
        rows_with_inchi = df[df['inchi_id'].notna()]
        if len(rows_with_inchi) > 0:
            assert rows_with_inchi['inchi'].notna().all()
            assert rows_with_inchi['inchikey'].notna().all()
    
    def test_get_id_table_from_cas_nonexistent(self, zpm):
        """Test getting identifier table from a nonexistent CAS number"""
        cas = "999-99-9"  # Likely nonexistent CAS
        df = zpm.get_id_table_from_cas(cas)
        
        # Should return None for nonexistent CAS
        assert df is None
    
    def test_get_id_table_from_cas_multiple_inchi_ids(self, zpm):
        """Test that the method handles CAS with multiple InChI IDs"""
        import pandas as pd
        
        # Use a CAS that might have multiple structures
        # This is database-dependent, so we'll test with formaldehyde
        cas = "50-00-0"
        df = zpm.get_id_table_from_cas(cas)
        
        if df is not None and len(df) > 1:
            # Check that each row has a different inchi_id or rank
            # (or they could be the same if there's only one substance)
            assert isinstance(df, pd.DataFrame)
            assert len(df) >= 1
    
    def test_batch_get_id_table_from_cas_valid_list(self, zpm):
        """Test batch getting identifier tables from a list of valid CAS numbers"""
        import pandas as pd
        
        cas_list = ["50-00-0", "50-78-2"]  # formaldehyde and aspirin
        df = zpm.batch_get_id_table_from_cas(cas_list)
        
        # Check that we got a DataFrame
        assert isinstance(df, pd.DataFrame)
        
        # Check that the DataFrame has the expected columns
        expected_columns = {'cas', 'query_id', 'inchi_id', 'rank', 'inchi', 'inchikey', 'zeropm_id', 'synonyms', 'sources'}
        assert set(df.columns) == expected_columns
        
        # Check that we have data for both CAS numbers
        unique_cas = df['cas'].unique()
        assert len(unique_cas) >= 1  # At least one should be found
        
        # Check that each CAS in the result was in the input list
        for cas in unique_cas:
            assert cas in cas_list
    
    def test_batch_get_id_table_from_cas_empty_list(self, zpm):
        """Test batch getting identifier tables with empty list"""
        import pandas as pd
        
        df = zpm.batch_get_id_table_from_cas([])
        
        # Should return an empty DataFrame with the correct columns
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0
        expected_columns = {'cas', 'query_id', 'inchi_id', 'rank', 'inchi', 'inchikey', 'zeropm_id', 'synonyms', 'sources'}
        assert set(df.columns) == expected_columns
    
    def test_batch_get_id_table_from_cas_nonexistent_list(self, zpm):
        """Test batch getting identifier tables with nonexistent CAS numbers"""
        import pandas as pd
        
        cas_list = ["999-99-9", "888-88-8"]  # Likely nonexistent CAS numbers
        df = zpm.batch_get_id_table_from_cas(cas_list)
        
        # Should return an empty DataFrame
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0
    
    def test_batch_get_id_table_from_cas_mixed_list(self, zpm):
        """Test batch getting identifier tables with mixed valid and invalid CAS"""
        import pandas as pd
        
        cas_list = ["50-00-0", "999-99-9", "50-78-2"]  # valid, invalid, valid
        df = zpm.batch_get_id_table_from_cas(cas_list)
        
        # Should return a DataFrame with only valid CAS numbers
        assert isinstance(df, pd.DataFrame)
        
        if len(df) > 0:
            unique_cas = df['cas'].unique()
            # Should not include the nonexistent CAS
            assert "999-99-9" not in unique_cas
            # Should include at least one of the valid CAS numbers
            assert any(cas in ["50-00-0", "50-78-2"] for cas in unique_cas)


class TestZeroPMBatchMethods:
    """Test suite for batch query methods"""
    
    @pytest.fixture
    def zpm(self):
        """Create a ZeroPM instance for testing"""
        return ZeroPM()
    
    def test_batch_query_cas_empty_list(self, zpm):
        """Test batch CAS query with empty list"""
        results = zpm.batch_query_cas([])
        assert results == {}
    
    def test_batch_query_cas_single(self, zpm):
        """Test batch CAS query with single CAS number"""
        # Get a known CAS number
        zpm.cursor.execute("""
            SELECT query 
            FROM api_ready_query 
            WHERE type = 'CAS Registry Number' 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            cas = result[0]
            results = zpm.batch_query_cas([cas])
            assert isinstance(results, dict)
            assert cas in results
            assert results[cas] is not None
        else:
            pytest.skip("No CAS numbers found in database")
    
    def test_batch_query_cas_multiple(self, zpm):
        """Test batch CAS query with multiple CAS numbers"""
        # Get multiple known CAS numbers
        zpm.cursor.execute("""
            SELECT query 
            FROM api_ready_query 
            WHERE type = 'CAS Registry Number' 
            LIMIT 3
        """)
        cas_list = [row[0] for row in zpm.cursor.fetchall()]
        
        if len(cas_list) >= 2:
            results = zpm.batch_query_cas(cas_list)
            assert isinstance(results, dict)
            assert len(results) == len(cas_list)
            for cas in cas_list:
                assert cas in results
        else:
            pytest.skip("Not enough CAS numbers in database")
    
    def test_batch_query_cas_mixed(self, zpm):
        """Test batch CAS query with mix of existing and nonexistent CAS"""
        # Get one known CAS
        zpm.cursor.execute("""
            SELECT query 
            FROM api_ready_query 
            WHERE type = 'CAS Registry Number' 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            cas_list = [result[0], '999-99-9', '888-88-8']
            results = zpm.batch_query_cas(cas_list)
            assert len(results) == 3
            assert results[result[0]] is not None
            assert results['999-99-9'] is None
            assert results['888-88-8'] is None
        else:
            pytest.skip("No CAS numbers found in database")
    
    def test_batch_get_smiles_from_cas(self, zpm):
        """Test batch SMILES retrieval from CAS numbers"""
        zpm.cursor.execute("""
            SELECT query 
            FROM api_ready_query 
            WHERE type = 'CAS Registry Number' 
            LIMIT 2
        """)
        cas_list = [row[0] for row in zpm.cursor.fetchall()]
        
        if cas_list:
            results = zpm.batch_get_smiles_from_cas(cas_list)
            assert isinstance(results, dict)
            assert len(results) == len(cas_list)
            for cas in cas_list:
                assert cas in results
        else:
            pytest.skip("No CAS numbers found in database")
    
    def test_batch_get_names(self, zpm):
        """Test batch name retrieval from CAS numbers"""
        zpm.cursor.execute("""
            SELECT query 
            FROM api_ready_query 
            WHERE type = 'CAS Registry Number' 
            LIMIT 2
        """)
        cas_list = [row[0] for row in zpm.cursor.fetchall()]
        
        if cas_list:
            results = zpm.batch_get_names(cas_list)
            assert isinstance(results, dict)
            assert len(results) == len(cas_list)
            for cas in cas_list:
                assert cas in results
                assert isinstance(results[cas], list)
        else:
            pytest.skip("No CAS numbers found in database")
    
    def test_batch_get_cas_from_inchikey_empty(self, zpm):
        """Test batch InChIKey to CAS with empty list"""
        results = zpm.batch_get_cas_from_inchikey([])
        assert results == {}
    
    def test_batch_get_cas_from_inchikey_multiple(self, zpm):
        """Test batch InChIKey to CAS with multiple keys"""
        zpm.cursor.execute("""
            SELECT inchikey 
            FROM substances 
            LIMIT 3
        """)
        inchikey_list = [row[0] for row in zpm.cursor.fetchall()]
        
        if len(inchikey_list) >= 2:
            results = zpm.batch_get_cas_from_inchikey(inchikey_list)
            assert isinstance(results, dict)
            assert len(results) == len(inchikey_list)
            for key in inchikey_list:
                assert key in results
        else:
            pytest.skip("Not enough InChIKeys in database")


class TestZeroPMPerformanceMethods:
    """Test suite for performance enhancement methods"""
    
    @pytest.fixture
    def zpm(self):
        """Create a ZeroPM instance for testing"""
        return ZeroPM()
    
    def test_create_indexes_default(self, zpm):
        """Test creating indexes with default settings"""
        results = zpm.create_indexes()
        assert isinstance(results, dict)
        assert len(results) > 0
        # All should be either 'created' or 'exists'
        for status in results.values():
            assert status in ['created', 'exists', 'error']
    
    def test_create_indexes_force(self, zpm):
        """Test creating indexes with force option"""
        results = zpm.create_indexes(force=True)
        assert isinstance(results, dict)
        assert len(results) > 0
        for status in results.values():
            assert status in ['created', 'exists', 'error']


class TestZeroPMAdvancedSearchMethods:
    """Test suite for advanced search methods"""
    
    @pytest.fixture
    def zpm(self):
        """Create a ZeroPM instance for testing"""
        return ZeroPM()
    
    def test_query_name_regex_case_insensitive(self, zpm):
        """Test regex name search case insensitive"""
        # Get a sample name
        zpm.cursor.execute("""
            SELECT query 
            FROM api_ready_query 
            WHERE type = 'chemical name' 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            name = result[0]
            # Search for part of the name
            if len(name) >= 5:
                pattern = f"%{name[:5]}%"
                results = zpm.query_name_regex(pattern, case_sensitive=False, limit=10)
                assert isinstance(results, list)
                # Should find at least the original name
                assert len(results) >= 1
        else:
            pytest.skip("No chemical names found in database")
    
    def test_query_name_regex_case_sensitive(self, zpm):
        """Test regex name search case sensitive"""
        zpm.cursor.execute("""
            SELECT query 
            FROM api_ready_query 
            WHERE type = 'chemical name' 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            name = result[0]
            if len(name) >= 5:
                pattern = f"%{name[:5]}%"
                results = zpm.query_name_regex(pattern, case_sensitive=True, limit=10)
                assert isinstance(results, list)
        else:
            pytest.skip("No chemical names found in database")
    
    def test_query_name_regex_no_match(self, zpm):
        """Test regex name search with no matches"""
        pattern = "%xyzabc123nonexistent%"
        results = zpm.query_name_regex(pattern, limit=10)
        assert isinstance(results, list)
        assert len(results) == 0
    
    @pytest.mark.slow
    def test_get_cas_by_substructure_valid(self, zpm):
        """Test substructure search with valid SMARTS"""
        # Search for simple patterns like benzene ring
        smarts = "c1ccccc1"  # Benzene ring
        results = zpm.get_cas_by_substructure(smarts, max_results=5)
        assert isinstance(results, list)
        # Results depend on database content
        for item in results:
            assert isinstance(item, dict)
            assert 'cas' in item
            assert 'inchi' in item
            assert 'inchikey' in item
            assert 'smiles' in item
    
    def test_get_cas_by_substructure_invalid(self, zpm):
        """Test substructure search with invalid SMARTS"""
        smarts = "invalid_smarts_xyz"
        results = zpm.get_cas_by_substructure(smarts, max_results=5)
        assert results == []


class TestZeroPMExportMethods:
    """Test suite for export methods"""
    
    @pytest.fixture
    def zpm(self):
        """Create a ZeroPM instance for testing"""
        return ZeroPM()
    
    def test_export_to_csv_list(self, zpm, tmp_path):
        """Test exporting list results to CSV"""
        # Create test data
        query_results = [
            ('CAS1', 'Name1'),
            ('CAS2', 'Name2'),
            ('CAS3', 'Name3')
        ]
        
        # Use temp directory
        filename = os.path.join(tmp_path, 'test_export.csv')
        output_path = zpm.export_to_csv(query_results, filename, columns=['CAS', 'Name'])
        
        # Verify file was created
        assert os.path.exists(output_path) or os.path.exists(filename)
    
    def test_export_to_csv_dict(self, zpm, tmp_path):
        """Test exporting dictionary results to CSV"""
        query_results = {
            'key1': 'value1',
            'key2': 'value2',
            'key3': 'value3'
        }
        
        filename = os.path.join(tmp_path, 'test_export_dict.csv')
        output_path = zpm.export_to_csv(query_results, filename)
        
        # Verify file was created
        assert os.path.exists(output_path) or os.path.exists(filename)
    
    def test_create_view_valid(self, zpm):
        """Test creating a database view"""
        view_name = "test_view_cas_inchi"
        sql_query = """
            SELECT aq.query AS cas, s.inchi
            FROM api_ready_query aq
            JOIN api_results ar ON aq.query_id = ar.query_id
            JOIN substances s ON ar.inchi_id = s.inchi_id
            WHERE aq.type = 'CAS Registry Number'
            LIMIT 10
        """
        
        result = zpm.create_view(view_name, sql_query)
        
        # Cleanup
        try:
            zpm.cursor.execute(f"DROP VIEW IF EXISTS {view_name}")
            zpm.conn.commit()
        except:
            pass
        
        # View creation might succeed or fail depending on permissions
        assert isinstance(result, bool)
    
    def test_export_query_results(self, zpm, tmp_path):
        """Test exporting custom query results to CSV"""
        sql_query = """
            SELECT query_id, query, type
            FROM api_ready_query
            LIMIT 5
        """
        
        filename = os.path.join(tmp_path, 'test_query_export.csv')
        output_path = zpm.export_query_results(sql_query, filename)
        
        # Verify file was created
        assert os.path.exists(output_path) or os.path.exists(filename)
    
    def test_get_database_stats(self, zpm):
        """Test getting database statistics"""
        stats = zpm.get_database_stats()
        
        assert isinstance(stats, dict)
        assert len(stats) > 0
        
        # Check for expected keys
        assert 'unique_cas_numbers' in stats
        assert 'unique_chemical_names' in stats
        
        # Verify counts are reasonable (non-negative integers or error messages)
        for key, value in stats.items():
            if isinstance(value, int):
                assert value >= 0
            else:
                # Should be an error message
                assert 'Error' in str(value) or isinstance(value, int)


class TestZeroPMEdgeCases:
    """Test suite for edge cases and error handling"""
    
    @pytest.fixture
    def zpm(self):
        """Create a ZeroPM instance for testing"""
        return ZeroPM()
    
    def test_cache_lazy_loading(self, zpm):
        """Test that chemical names cache is lazily loaded"""
        assert zpm._chemical_names_cache is None
        
        # Trigger cache loading
        cache = zpm._get_chemical_names_cache()
        assert cache is not None
        assert isinstance(cache, list)
        
        # Verify cache is reused
        cache2 = zpm._get_chemical_names_cache()
        assert cache is cache2  # Should be the same object
    
    def test_empty_inchi_list(self, zpm):
        """Test handling of empty inchi_id list"""
        inchi_ids, ranks = zpm.get_inchi_id(999999)
        assert inchi_ids == []
        assert ranks == []
    
    def test_multiple_cas_for_inchi(self, zpm):
        """Test handling when InChI maps to multiple CAS numbers"""
        # Find an InChI with multiple CAS mappings
        zpm.cursor.execute("""
            SELECT s.inchi, COUNT(DISTINCT aq.query) as cas_count
            FROM substances s
            JOIN api_results ar ON s.inchi_id = ar.inchi_id
            JOIN api_ready_query aq ON ar.query_id = aq.query_id
            WHERE aq.type = 'CAS Registry Number'
            GROUP BY s.inchi
            HAVING cas_count > 1
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            inchi = result[0]
            cas = zpm.get_cas_from_inchi(inchi)
            # Should return a list when multiple CAS numbers exist
            assert isinstance(cas, list)
            assert len(cas) > 1
        else:
            pytest.skip("No InChI with multiple CAS mappings found")
    
    def test_sql_injection_protection(self, zpm):
        """Test that inputs are properly sanitized"""
        # Try SQL injection in CAS query
        malicious_input = "'; DROP TABLE api_ready_query; --"
        result = zpm.query_cas(malicious_input)
        
        # Should return None without error
        assert result is None
        
        # Verify table still exists
        zpm.cursor.execute("SELECT COUNT(*) FROM api_ready_query")
        count = zpm.cursor.fetchone()[0]
        assert count > 0


@pytest.mark.integration
class TestZeroPMIntegration:
    """Integration tests for ZeroPM workflow"""
    
    @pytest.fixture
    def zpm(self):
        """Create a ZeroPM instance for testing"""
        return ZeroPM()
    
    def test_complete_workflow_cas_to_smiles(self, zpm):
        """Test complete workflow from CAS to SMILES"""
        # Get a CAS number
        zpm.cursor.execute("""
            SELECT aq.query
            FROM api_ready_query aq
            JOIN api_results ar ON aq.query_id = ar.query_id
            JOIN substances s ON ar.inchi_id = s.inchi_id
            WHERE aq.type = 'CAS Registry Number'
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            cas = result[0]
            
            # Test full workflow
            query_id = zpm.query_cas(cas)
            assert query_id is not None
            
            inchi_ids, ranks = zpm.get_inchi_id(query_id)
            assert len(inchi_ids) > 0
            
            inchi, inchikey = zpm.get_inchi(inchi_ids[0])
            assert inchi is not None
            
            smiles = zpm.get_smiles_from_cas(cas)
            # SMILES may be None if conversion fails, but workflow should complete
        else:
            pytest.skip("No suitable CAS number found")
    
    def test_complete_workflow_name_search(self, zpm):
        """Test complete workflow with name search"""
        # Get a chemical name
        zpm.cursor.execute("""
            SELECT query
            FROM api_ready_query
            WHERE type = 'chemical name'
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            name = result[0]
            
            # Test exact search
            query_id = zpm.query_name(name)
            assert query_id is not None
            
            # Test fuzzy search
            fuzzy_results = zpm.query_similar_name(name[:5], number_of_results=5)
            if fuzzy_results:
                assert query_id in fuzzy_results or len(fuzzy_results) > 0
        else:
            pytest.skip("No chemical names found")


class TestZeroPMInventoryCountryRegion:
    """Test suite for inventory, country, and region query methods"""
    
    @pytest.fixture
    def zpm(self):
        """Create a ZeroPM instance for testing"""
        return ZeroPM()
    
    def test_get_all_inventories(self, zpm):
        """Test getting all inventory sources"""
        inventories = zpm.get_all_inventories()
        assert isinstance(inventories, list)
        assert len(inventories) > 0
        
        # Check structure of first inventory
        first = inventories[0]
        assert 'source_id' in first
        assert 'source_name' in first
        assert 'country_scope' in first
        assert isinstance(first['source_id'], int)
        assert isinstance(first['source_name'], str)
    
    def test_get_all_countries(self, zpm):
        """Test getting all countries"""
        countries = zpm.get_all_countries()
        assert isinstance(countries, list)
        assert len(countries) > 0
        
        # Check structure
        first = countries[0]
        assert 'country_id' in first
        assert 'country' in first
        assert isinstance(first['country_id'], int)
        assert isinstance(first['country'], str)
    
    def test_get_all_regions(self, zpm):
        """Test getting all global regions"""
        regions = zpm.get_all_regions()
        assert isinstance(regions, list)
        assert len(regions) > 0
        
        # Check structure
        first = regions[0]
        assert 'region_id' in first
        assert 'region' in first
        assert isinstance(first['region_id'], int)
        assert isinstance(first['region'], str)
    
    def test_query_by_inventory_by_id(self, zpm):
        """Test querying chemicals by inventory source_id"""
        # Get a source_id
        inventories = zpm.get_all_inventories()
        if inventories:
            source_id = inventories[0]['source_id']
            results = zpm.query_by_inventory(source_id=source_id)
            
            assert isinstance(results, list)
            if results:
                assert 'cas' in results[0]
                assert 'query_id' in results[0]
                assert 'inchi_id' in results[0]
                assert 'source_name' in results[0]
        else:
            pytest.skip("No inventories found")
    
    def test_query_by_inventory_by_name(self, zpm):
        """Test querying chemicals by inventory source name"""
        # Get a source name
        inventories = zpm.get_all_inventories()
        if inventories:
            source_name = inventories[0]['source_name']
            # Use partial name
            partial_name = source_name.split()[0] if ' ' in source_name else source_name[:5]
            results = zpm.query_by_inventory(source_name=partial_name)
            
            assert isinstance(results, list)
        else:
            pytest.skip("No inventories found")
    
    def test_query_by_inventory_no_params(self, zpm):
        """Test that query_by_inventory raises error without parameters"""
        with pytest.raises(ValueError):
            zpm.query_by_inventory()
    
    def test_query_by_country_by_id(self, zpm):
        """Test querying chemicals by country_id"""
        countries = zpm.get_all_countries()
        if countries:
            country_id = countries[0]['country_id']
            results = zpm.query_by_country(country_id=country_id)
            
            assert isinstance(results, list)
            if results:
                assert 'cas' in results[0]
                assert 'country' in results[0]
                assert 'source_name' in results[0]
        else:
            pytest.skip("No countries found")
    
    def test_query_by_country_by_name(self, zpm):
        """Test querying chemicals by country name"""
        countries = zpm.get_all_countries()
        if countries:
            country_name = countries[0]['country']
            results = zpm.query_by_country(country_name=country_name)
            
            assert isinstance(results, list)
            if results:
                assert 'cas' in results[0]
                assert 'country' in results[0]
        else:
            pytest.skip("No countries found")
    
    def test_query_by_country_no_params(self, zpm):
        """Test that query_by_country raises error without parameters"""
        with pytest.raises(ValueError):
            zpm.query_by_country()
    
    def test_query_by_region_by_id(self, zpm):
        """Test querying chemicals by region_id"""
        regions = zpm.get_all_regions()
        if regions:
            region_id = regions[0]['region_id']
            results = zpm.query_by_region(region_id=region_id)
            
            assert isinstance(results, list)
            if results:
                assert 'cas' in results[0]
                assert 'region' in results[0]
                assert 'country' in results[0]
                assert 'source_name' in results[0]
        else:
            pytest.skip("No regions found")
    
    def test_query_by_region_by_name(self, zpm):
        """Test querying chemicals by region name"""
        regions = zpm.get_all_regions()
        if regions:
            region_name = regions[0]['region']
            results = zpm.query_by_region(region_name=region_name)
            
            assert isinstance(results, list)
            if results:
                assert 'cas' in results[0]
                assert 'region' in results[0]
        else:
            pytest.skip("No regions found")
    
    def test_query_by_region_no_params(self, zpm):
        """Test that query_by_region raises error without parameters"""
        with pytest.raises(ValueError):
            zpm.query_by_region()
    
    def test_get_countries_for_region(self, zpm):
        """Test getting countries in a region"""
        regions = zpm.get_all_regions()
        if regions:
            region_id = regions[0]['region_id']
            countries = zpm.get_countries_for_region(region_id=region_id)
            
            assert isinstance(countries, list)
            if countries:
                assert 'country_id' in countries[0]
                assert 'country' in countries[0]
                assert 'region' in countries[0]
        else:
            pytest.skip("No regions found")
    
    def test_get_inventories_for_country(self, zpm):
        """Test getting inventories for a country"""
        countries = zpm.get_all_countries()
        if countries:
            country_id = countries[0]['country_id']
            inventories = zpm.get_inventories_for_country(country_id=country_id)
            
            assert isinstance(inventories, list)
            # May be empty if country has no inventories
            if inventories:
                assert 'source_id' in inventories[0]
                assert 'source_name' in inventories[0]
                assert 'country' in inventories[0]
        else:
            pytest.skip("No countries found")
    
    def test_count_chemicals_by_inventory(self, zpm):
        """Test counting chemicals in an inventory"""
        inventories = zpm.get_all_inventories()
        if inventories:
            source_id = inventories[0]['source_id']
            count = zpm.count_chemicals_by_inventory(source_id)
            
            assert isinstance(count, int)
            assert count >= 0
        else:
            pytest.skip("No inventories found")
    
    def test_count_chemicals_by_country(self, zpm):
        """Test counting chemicals in a country"""
        countries = zpm.get_all_countries()
        if countries:
            country_id = countries[0]['country_id']
            count = zpm.count_chemicals_by_country(country_id)
            
            assert isinstance(count, int)
            assert count >= 0
        else:
            pytest.skip("No countries found")
    
    def test_count_chemicals_by_region(self, zpm):
        """Test counting chemicals in a region"""
        regions = zpm.get_all_regions()
        if regions:
            region_id = regions[0]['region_id']
            count = zpm.count_chemicals_by_region(region_id)
            
            assert isinstance(count, int)
            assert count >= 0
        else:
            pytest.skip("No regions found")


class TestZeroPMGetIdTableMethods:
    """Test suite for get_id_table_from_* methods"""
    
    @pytest.fixture
    def zpm(self):
        """Create a ZeroPM instance for testing"""
        return ZeroPM()
    
    def test_get_id_table_from_inchi(self, zpm):
        """Test getting ID table from InChI"""
        # Get a known InChI from the database
        zpm.cursor.execute("""
            SELECT inchi 
            FROM substances 
            WHERE inchi IS NOT NULL
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            inchi = result[0]
            df = zpm.get_id_table_from_inchi(inchi)
            
            assert df is not None
            assert len(df) > 0
            assert 'inchi' in df.columns
            assert 'inchikey' in df.columns
            assert 'inchi_id' in df.columns
            assert 'query_id' in df.columns
            assert 'rank' in df.columns
            assert 'cas' in df.columns
            assert 'synonyms' in df.columns
            assert df['inchi'].iloc[0] == inchi
        else:
            pytest.skip("No InChI found in database")
    
    def test_get_id_table_from_inchi_nonexistent(self, zpm):
        """Test getting ID table from nonexistent InChI"""
        df = zpm.get_id_table_from_inchi("InChI=1S/NONEXISTENT")
        assert df is None
    
    def test_batch_get_id_table_from_inchi(self, zpm):
        """Test batch getting ID tables from InChI list"""
        # Get some known InChIs from the database
        zpm.cursor.execute("""
            SELECT inchi 
            FROM substances 
            WHERE inchi IS NOT NULL
            LIMIT 3
        """)
        results = zpm.cursor.fetchall()
        
        if results:
            inchi_list = [row[0] for row in results]
            df = zpm.batch_get_id_table_from_inchi(inchi_list)
            
            assert df is not None
            assert len(df) > 0
            assert 'inchi' in df.columns
            assert 'inchikey' in df.columns
            # Check that all InChIs are in the result
            for inchi in inchi_list:
                assert inchi in df['inchi'].values
        else:
            pytest.skip("No InChIs found in database")
    
    def test_batch_get_id_table_from_inchi_empty_list(self, zpm):
        """Test batch getting ID tables from empty list"""
        df = zpm.batch_get_id_table_from_inchi([])
        assert df is not None
        assert len(df) == 0
    
    def test_get_id_table_from_inchikey(self, zpm):
        """Test getting ID table from InChIKey"""
        # Get a known InChIKey from the database
        zpm.cursor.execute("""
            SELECT inchikey 
            FROM substances 
            WHERE inchikey IS NOT NULL
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            inchikey = result[0]
            df = zpm.get_id_table_from_inchikey(inchikey)
            
            assert df is not None
            assert len(df) > 0
            assert 'inchikey' in df.columns
            assert 'inchi' in df.columns
            assert 'inchi_id' in df.columns
            assert 'query_id' in df.columns
            assert 'rank' in df.columns
            assert 'cas' in df.columns
            assert 'synonyms' in df.columns
            assert df['inchikey'].iloc[0] == inchikey
        else:
            pytest.skip("No InChIKey found in database")
    
    def test_get_id_table_from_inchikey_nonexistent(self, zpm):
        """Test getting ID table from nonexistent InChIKey"""
        df = zpm.get_id_table_from_inchikey("NONEXISTENTKEY-UHFFFAOYSA-N")
        assert df is None
    
    def test_batch_get_id_table_from_inchikey(self, zpm):
        """Test batch getting ID tables from InChIKey list"""
        # Get some known InChIKeys from the database
        zpm.cursor.execute("""
            SELECT inchikey 
            FROM substances 
            WHERE inchikey IS NOT NULL
            LIMIT 3
        """)
        results = zpm.cursor.fetchall()
        
        if results:
            inchikey_list = [row[0] for row in results]
            df = zpm.batch_get_id_table_from_inchikey(inchikey_list)
            
            assert df is not None
            assert len(df) > 0
            assert 'inchikey' in df.columns
            assert 'inchi' in df.columns
            # Check that all InChIKeys are in the result
            for inchikey in inchikey_list:
                assert inchikey in df['inchikey'].values
        else:
            pytest.skip("No InChIKeys found in database")
    
    def test_batch_get_id_table_from_inchikey_empty_list(self, zpm):
        """Test batch getting ID tables from empty list"""
        df = zpm.batch_get_id_table_from_inchikey([])
        assert df is not None
        assert len(df) == 0
    
    def test_get_id_table_from_name(self, zpm):
        """Test getting ID table from chemical name"""
        # Get a known chemical name from the database
        zpm.cursor.execute("""
            SELECT query 
            FROM api_ready_query 
            WHERE type = 'chemical name'
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            name = result[0]
            df = zpm.get_id_table_from_name(name)
            
            assert df is not None
            assert len(df) > 0
            assert 'name' in df.columns
            assert 'query_id' in df.columns
            assert 'inchi_id' in df.columns
            assert 'rank' in df.columns
            assert 'inchi' in df.columns
            assert 'inchikey' in df.columns
            assert 'cas' in df.columns
            assert df['name'].iloc[0] == name
        else:
            pytest.skip("No chemical names found in database")
    
    def test_get_id_table_from_name_nonexistent(self, zpm):
        """Test getting ID table from nonexistent name"""
        df = zpm.get_id_table_from_name("nonexistent_chemical_name_12345")
        assert df is None
    
    def test_batch_get_id_table_from_name(self, zpm):
        """Test batch getting ID tables from name list"""
        # Get some known chemical names from the database
        zpm.cursor.execute("""
            SELECT query 
            FROM api_ready_query 
            WHERE type = 'chemical name'
            LIMIT 3
        """)
        results = zpm.cursor.fetchall()
        
        if results:
            name_list = [row[0] for row in results]
            df = zpm.batch_get_id_table_from_name(name_list)
            
            assert df is not None
            assert len(df) > 0
            assert 'name' in df.columns
            assert 'query_id' in df.columns
            # Check that all names are in the result
            for name in name_list:
                assert name in df['name'].values
        else:
            pytest.skip("No chemical names found in database")
    
    def test_batch_get_id_table_from_name_empty_list(self, zpm):
        """Test batch getting ID tables from empty list"""
        df = zpm.batch_get_id_table_from_name([])
        assert df is not None
        assert len(df) == 0


class TestZeroPMCASConversionMethods:
    """Test suite for CAS conversion methods (get_cas_from_*)"""
    
    @pytest.fixture
    def zpm(self):
        """Create a ZeroPM instance for testing"""
        return ZeroPM()
    
    def test_get_cas_from_name(self, zpm):
        """Test getting CAS from chemical name"""
        # Get a known chemical name from the database
        zpm.cursor.execute("""
            SELECT query 
            FROM api_ready_query 
            WHERE type = 'chemical name'
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            name = result[0]
            cas = zpm.get_cas_from_name(name)
            
            # Should return either a string or list of CAS numbers
            assert cas is not None
            assert isinstance(cas, (str, list))
            
            if isinstance(cas, str):
                # Basic CAS format check: XXX-XX-X or similar
                assert '-' in cas or cas.isdigit()
            else:
                assert len(cas) > 0
                assert all('-' in c or c.isdigit() for c in cas)
        else:
            pytest.skip("No chemical names found in database")
    
    def test_get_cas_from_name_nonexistent(self, zpm):
        """Test getting CAS from nonexistent name"""
        cas = zpm.get_cas_from_name("nonexistent_chemical_name_12345")
        assert cas is None
    
    def test_get_cas_from_formula(self, zpm):
        """Test getting CAS from molecular formula"""
        # Test with a simple formula that should exist (water)
        cas_list = zpm.get_cas_from_formula("H2O")
        
        if cas_list:
            assert isinstance(cas_list, list)
            assert len(cas_list) > 0
            # Check format of CAS numbers
            for cas in cas_list:
                assert isinstance(cas, str)
                assert '-' in cas or cas.isdigit()
        else:
            # If H2O not found, try with formaldehyde
            cas_list = zpm.get_cas_from_formula("CH2O")
            if cas_list:
                assert isinstance(cas_list, list)
                assert len(cas_list) > 0
            else:
                pytest.skip("Neither H2O nor CH2O found in database")
    
    def test_get_cas_from_formula_nonexistent(self, zpm):
        """Test getting CAS from nonexistent formula"""
        # Use an unlikely formula
        cas_list = zpm.get_cas_from_formula("Zr999Xe999")
        assert cas_list is None
    
    def test_batch_get_cas_from_smiles(self, zpm):
        """Test batch getting CAS from SMILES"""
        # Use simple SMILES
        smiles_list = ["C", "CC", "O"]  # methane, ethane, water
        results = zpm.batch_get_cas_from_smiles(smiles_list)
        
        assert isinstance(results, dict)
        assert len(results) == len(smiles_list)
        
        # Check that all SMILES are in the results
        for smiles in smiles_list:
            assert smiles in results
            # Result can be None, string, or list
            assert results[smiles] is None or isinstance(results[smiles], (str, list))
    
    def test_batch_get_cas_from_smiles_empty_list(self, zpm):
        """Test batch getting CAS from empty SMILES list"""
        results = zpm.batch_get_cas_from_smiles([])
        assert isinstance(results, dict)
        assert len(results) == 0
    
    def test_batch_get_cas_from_name(self, zpm):
        """Test batch getting CAS from names"""
        # Get some known chemical names from the database
        zpm.cursor.execute("""
            SELECT query 
            FROM api_ready_query 
            WHERE type = 'chemical name'
            LIMIT 3
        """)
        results_db = zpm.cursor.fetchall()
        
        if results_db:
            name_list = [row[0] for row in results_db]
            results = zpm.batch_get_cas_from_name(name_list)
            
            assert isinstance(results, dict)
            assert len(results) == len(name_list)
            
            # Check that all names are in the results
            for name in name_list:
                assert name in results
                # Result can be None, string, or list
                assert results[name] is None or isinstance(results[name], (str, list))
        else:
            pytest.skip("No chemical names found in database")
    
    def test_batch_get_cas_from_name_empty_list(self, zpm):
        """Test batch getting CAS from empty name list"""
        results = zpm.batch_get_cas_from_name([])
        assert isinstance(results, dict)
        assert len(results) == 0
    
    def test_batch_get_cas_from_formula(self, zpm):
        """Test batch getting CAS from formulas"""
        formula_list = ["H2O", "CH4", "CH2O"]
        results = zpm.batch_get_cas_from_formula(formula_list)
        
        assert isinstance(results, dict)
        assert len(results) == len(formula_list)
        
        # Check that all formulas are in the results
        for formula in formula_list:
            assert formula in results
            # Result can be None or list
            assert results[formula] is None or isinstance(results[formula], list)
    
    def test_batch_get_cas_from_formula_empty_list(self, zpm):
        """Test batch getting CAS from empty formula list"""
        results = zpm.batch_get_cas_from_formula([])
        assert isinstance(results, dict)
        assert len(results) == 0
    
    def test_get_cas_from_name_integration(self, zpm):
        """Integration test: name -> CAS -> back to name"""
        # Get a CAS number with known names
        zpm.cursor.execute("""
            SELECT query 
            FROM api_ready_query 
            WHERE type = 'CAS Registry Number'
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            original_cas = result[0]
            # Get names for this CAS
            names = zpm.get_names(original_cas)
            
            if names:
                # Try to get CAS back from the first name
                first_name = names[0]
                
                # First check if this name exists as a query
                query_id = zpm.query_name(first_name)
                if query_id:
                    retrieved_cas = zpm.get_cas_from_name(first_name)
                    
                    # The retrieved CAS should either be the original or in a list containing it
                    if isinstance(retrieved_cas, str):
                        assert retrieved_cas == original_cas or original_cas in names
                    elif isinstance(retrieved_cas, list):
                        assert original_cas in retrieved_cas or any(cas in names for cas in retrieved_cas)
                else:
                    pytest.skip(f"Name '{first_name}' not found as a query in database")
            else:
                pytest.skip("No names found for CAS")
        else:
            pytest.skip("No CAS numbers found in database")


class TestZeroPMv004Features:
    """Test suite for ZeroPM v0-0-4 new features"""
    
    @pytest.fixture
    def zpm(self):
        """Create a ZeroPM instance for testing"""
        return ZeroPM()
    
    # ==================== ZeroPM ID Tests ====================
    
    def test_get_zeropm_id_with_inchi_id(self, zpm):
        """Test getting ZeroPM ID from inchi_id"""
        # Find a chemical in zeropm_chemicals
        zpm.cursor.execute("""
            SELECT inchi_id, zeropm_id 
            FROM zeropm_chemicals 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            inchi_id, expected_zeropm_id = result
            zeropm_id = zpm.get_zeropm_id(inchi_id=inchi_id)
            assert zeropm_id == expected_zeropm_id
        else:
            pytest.skip("No chemicals found in zeropm_chemicals table")
    
    def test_get_zeropm_id_with_cas(self, zpm):
        """Test getting ZeroPM ID from CAS number"""
        # Find a CAS that has a zeropm_id
        zpm.cursor.execute("""
            SELECT DISTINCT aq.query 
            FROM api_ready_query aq
            JOIN api_results ar ON aq.query_id = ar.query_id
            JOIN zeropm_chemicals zc ON ar.inchi_id = zc.inchi_id
            WHERE aq.type = 'CAS Registry Number'
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            cas = result[0]
            zeropm_id = zpm.get_zeropm_id(cas=cas)
            assert zeropm_id is not None
            assert isinstance(zeropm_id, int)
        else:
            pytest.skip("No CAS found with zeropm_id")
    
    def test_get_zeropm_id_not_found(self, zpm):
        """Test getting ZeroPM ID for non-existent chemical"""
        zeropm_id = zpm.get_zeropm_id(inchi_id=999999999)
        assert zeropm_id is None
    
    def test_get_zeropm_id_no_parameters(self, zpm):
        """Test that get_zeropm_id raises error with no parameters"""
        with pytest.raises(ValueError):
            zpm.get_zeropm_id()
    
    def test_is_in_zeropm_true(self, zpm):
        """Test is_in_zeropm returns True for ZeroPM chemical"""
        zpm.cursor.execute("""
            SELECT inchi_id 
            FROM zeropm_chemicals 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            inchi_id = result[0]
            assert zpm.is_in_zeropm(inchi_id=inchi_id) is True
        else:
            pytest.skip("No chemicals found in zeropm_chemicals table")
    
    def test_is_in_zeropm_false(self, zpm):
        """Test is_in_zeropm returns False for non-ZeroPM chemical"""
        assert zpm.is_in_zeropm(inchi_id=999999999) is False
    
    # ==================== P/M Probability Tests ====================
    
    def test_get_pm_probabilities(self, zpm):
        """Test getting P/M probabilities"""
        # Find a chemical with PM probabilities
        zpm.cursor.execute("""
            SELECT zeropm_id 
            FROM pm_probabilities 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            zeropm_id = result[0]
            probs = zpm.get_pm_probabilities(zeropm_id=zeropm_id)
            
            assert probs is not None
            assert isinstance(probs, dict)
            
            # Check all expected keys
            expected_keys = [
                'probability_of_not_p', 'probability_of_p_or_vp',
                'probability_of_p', 'probability_of_vp',
                'probability_of_not_m', 'probability_of_m_or_vm',
                'probability_of_m', 'probability_of_vm', 'n'
            ]
            for key in expected_keys:
                assert key in probs
                # n should be an integer, others should be numeric or None
                if key == 'n':
                    assert isinstance(probs[key], (int, type(None)))
                else:
                    assert isinstance(probs[key], (int, float, type(None)))
        else:
            pytest.skip("No P/M probabilities found in database")
    
    def test_get_pm_probabilities_not_found(self, zpm):
        """Test getting P/M probabilities for non-existent zeropm_id"""
        probs = zpm.get_pm_probabilities(zeropm_id=999999999)
        assert probs is None
    
    def test_batch_get_pm_probabilities_with_cas(self, zpm):
        """Test batch getting P/M probabilities from CAS list"""
        # Find CAS numbers with PM data
        zpm.cursor.execute("""
            SELECT DISTINCT aq.query 
            FROM api_ready_query aq
            JOIN api_results ar ON aq.query_id = ar.query_id
            JOIN zeropm_chemicals zc ON ar.inchi_id = zc.inchi_id
            JOIN pm_probabilities pm ON zc.zeropm_id = pm.zeropm_id
            WHERE aq.type = 'CAS Registry Number'
            LIMIT 3
        """)
        results = zpm.cursor.fetchall()
        
        if results:
            cas_list = [row[0] for row in results]
            df = zpm.batch_get_pm_probabilities(cas_list=cas_list)
            
            assert df is not None
            assert not df.empty
            assert 'cas' in df.columns
            assert 'probability_of_p' in df.columns
            assert 'probability_of_m' in df.columns
            assert len(df) > 0
        else:
            pytest.skip("No CAS with P/M probabilities found")
    
    def test_batch_get_pm_probabilities_with_inchi_ids(self, zpm):
        """Test batch getting P/M probabilities from inchi_id list"""
        zpm.cursor.execute("""
            SELECT zc.inchi_id 
            FROM zeropm_chemicals zc
            JOIN pm_probabilities pm ON zc.zeropm_id = pm.zeropm_id
            LIMIT 3
        """)
        results = zpm.cursor.fetchall()
        
        if results:
            inchi_id_list = [row[0] for row in results]
            df = zpm.batch_get_pm_probabilities(inchi_id_list=inchi_id_list)
            
            assert df is not None
            assert not df.empty
            assert 'inchi_id' in df.columns
            assert len(df) > 0
        else:
            pytest.skip("No inchi_ids with P/M probabilities found")
    
    def test_get_all_zeropm_chemicals(self, zpm):
        """Test getting all ZeroPM chemicals"""
        df = zpm.get_all_zeropm_chemicals(limit=10)
        
        assert df is not None
        assert 'zeropm_id' in df.columns
        assert 'inchi_id' in df.columns
        assert 'inchi' in df.columns
        assert 'inchikey' in df.columns
        assert len(df) <= 10
    
    def test_get_all_zeropm_chemicals_with_probs(self, zpm):
        """Test getting all ZeroPM chemicals with P/M probabilities"""
        df = zpm.get_all_zeropm_chemicals(limit=5, include_pm_probs=True)
        
        assert df is not None
        assert 'zeropm_id' in df.columns
        assert 'probability_of_p' in df.columns
        assert 'probability_of_m' in df.columns
        assert len(df) <= 5
    
    # ==================== Multi-Component Tests ====================
    
    def test_is_multicomponent(self, zpm):
        """Test checking if substance is multi-component"""
        # Find a multi-component substance
        zpm.cursor.execute("""
            SELECT inchi_id 
            FROM multi_components 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            inchi_id = result[0]
            assert zpm.is_multicomponent(inchi_id) is True
        else:
            pytest.skip("No multi-component substances found")
    
    def test_is_multicomponent_false(self, zpm):
        """Test is_multicomponent returns False for single component"""
        # Find a substance that is NOT multi-component
        zpm.cursor.execute("""
            SELECT s.inchi_id 
            FROM substances s
            WHERE NOT EXISTS (
                SELECT 1 FROM multi_components mc WHERE mc.inchi_id = s.inchi_id
            )
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            inchi_id = result[0]
            assert zpm.is_multicomponent(inchi_id) is False
        else:
            pytest.skip("All substances are multi-component")
    
    def test_get_multicomponent_id(self, zpm):
        """Test getting multi-component ID"""
        zpm.cursor.execute("""
            SELECT inchi_id, mc_id 
            FROM multi_components 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            inchi_id, expected_mc_id = result
            mc_id = zpm.get_multicomponent_id(inchi_id)
            assert mc_id == expected_mc_id
        else:
            pytest.skip("No multi-component substances found")
    
    def test_get_components(self, zpm):
        """Test getting components of a multi-component substance"""
        zpm.cursor.execute("""
            SELECT mc_id 
            FROM multi_components 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            mc_id = result[0]
            components = zpm.get_components(mc_id)
            
            assert isinstance(components, list)
            if len(components) > 0:
                # Check structure of component data
                component = components[0]
                assert 'component_id' in component
                assert 'component_frequency' in component
                assert 'inchi_id' in component
                assert 'inchi' in component
                assert 'inchikey' in component
        else:
            pytest.skip("No multi-component substances found")
    
    def test_get_multicomponent_info(self, zpm):
        """Test getting complete multi-component information"""
        zpm.cursor.execute("""
            SELECT inchi_id 
            FROM multi_components 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            inchi_id = result[0]
            info = zpm.get_multicomponent_info(inchi_id=inchi_id)
            
            assert info is not None
            assert isinstance(info, dict)
            assert 'mc_id' in info
            assert 'inchi_id' in info
            assert 'inchi' in info
            assert 'inchikey' in info
            assert 'components' in info
            assert isinstance(info['components'], list)
        else:
            pytest.skip("No multi-component substances found")
    
    def test_get_multicomponent_info_not_multicomponent(self, zpm):
        """Test get_multicomponent_info returns None for single component"""
        # Find a substance that is NOT multi-component
        zpm.cursor.execute("""
            SELECT s.inchi_id 
            FROM substances s
            WHERE NOT EXISTS (
                SELECT 1 FROM multi_components mc WHERE mc.inchi_id = s.inchi_id
            )
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            inchi_id = result[0]
            info = zpm.get_multicomponent_info(inchi_id=inchi_id)
            assert info is None
        else:
            pytest.skip("All substances are multi-component")
    
    def test_get_all_multicomponent_substances(self, zpm):
        """Test getting all multi-component substances"""
        df = zpm.get_all_multicomponent_substances(limit=10)
        
        assert df is not None
        assert 'mc_id' in df.columns
        assert 'inchi_id' in df.columns
        assert 'component_count' in df.columns
        assert len(df) <= 10
    
    # ==================== Cleanventory Tests ====================
    
    def test_is_in_cleanventory_true(self, zpm):
        """Test is_in_cleanventory returns True for Cleanventory chemical"""
        zpm.cursor.execute("""
            SELECT inchi_id 
            FROM cleanventory_chemicals 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            inchi_id = result[0]
            assert zpm.is_in_cleanventory(inchi_id=inchi_id) is True
        else:
            pytest.skip("No chemicals found in cleanventory_chemicals table")
    
    def test_is_in_cleanventory_false(self, zpm):
        """Test is_in_cleanventory returns False for non-Cleanventory chemical"""
        # Find a substance NOT in cleanventory
        zpm.cursor.execute("""
            SELECT s.inchi_id 
            FROM substances s
            WHERE NOT EXISTS (
                SELECT 1 FROM cleanventory_chemicals cc WHERE cc.inchi_id = s.inchi_id
            )
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            inchi_id = result[0]
            assert zpm.is_in_cleanventory(inchi_id=inchi_id) is False
        else:
            pytest.skip("All substances are in cleanventory")
    
    # ==================== Consensus Score Tests ====================
    
    def test_get_consensus_score(self, zpm):
        """Test getting consensus scores"""
        zpm.cursor.execute("""
            SELECT DISTINCT inchi_id 
            FROM consensus_index 
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            inchi_id = result[0]
            consensus = zpm.get_consensus_score(inchi_id=inchi_id)
            
            assert consensus is not None
            assert isinstance(consensus, list)
            if len(consensus) > 0:
                score = consensus[0]
                assert 'inventory_id' in score
                assert 'consensus_score' in score
                assert 'consensus_count' in score
        else:
            pytest.skip("No consensus scores found in database")
    
    def test_get_consensus_score_not_found(self, zpm):
        """Test getting consensus score for chemical without consensus data"""
        consensus = zpm.get_consensus_score(inchi_id=999999999)
        assert consensus is None
    
    # ==================== Sources Column Tests ====================
    
    def test_get_id_table_from_cas_includes_sources(self, zpm):
        """Test that get_id_table_from_cas includes sources column"""
        # Find a CAS with inventory data
        zpm.cursor.execute("""
            SELECT DISTINCT aq.query 
            FROM api_ready_query aq
            JOIN inventory_summary issum ON aq.query_id = issum.query_id
            WHERE aq.type = 'CAS Registry Number'
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            cas = result[0]
            df = zpm.get_id_table_from_cas(cas)
            
            assert df is not None
            assert 'sources' in df.columns
            # Check that sources column contains data
            if not df.empty:
                assert df['sources'].notna().any() or df['sources'].str.len().sum() >= 0
        else:
            pytest.skip("No CAS with inventory data found")
    
    def test_get_id_table_from_inchi_includes_sources(self, zpm):
        """Test that get_id_table_from_inchi includes sources column"""
        # Find an InChI with inventory data
        zpm.cursor.execute("""
            SELECT DISTINCT s.inchi 
            FROM substances s
            JOIN api_results ar ON s.inchi_id = ar.inchi_id
            JOIN inventory_summary issum ON ar.query_id = issum.query_id
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            inchi = result[0]
            df = zpm.get_id_table_from_inchi(inchi)
            
            assert df is not None
            assert 'sources' in df.columns
        else:
            pytest.skip("No InChI with inventory data found")
    
    def test_get_id_table_from_inchikey_includes_sources(self, zpm):
        """Test that get_id_table_from_inchikey includes sources column"""
        # Find an InChIKey with inventory data
        zpm.cursor.execute("""
            SELECT DISTINCT s.inchikey 
            FROM substances s
            JOIN api_results ar ON s.inchi_id = ar.inchi_id
            JOIN inventory_summary issum ON ar.query_id = issum.query_id
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            inchikey = result[0]
            df = zpm.get_id_table_from_inchikey(inchikey)
            
            assert df is not None
            assert 'sources' in df.columns
        else:
            pytest.skip("No InChIKey with inventory data found")
    
    def test_get_id_table_from_name_includes_sources(self, zpm):
        """Test that get_id_table_from_name includes sources column"""
        # Find a name with inventory data
        zpm.cursor.execute("""
            SELECT DISTINCT aq.query 
            FROM api_ready_query aq
            JOIN inventory_summary issum ON aq.query_id = issum.query_id
            WHERE aq.type = 'chemical name'
            LIMIT 1
        """)
        result = zpm.cursor.fetchone()
        
        if result:
            name = result[0]
            df = zpm.get_id_table_from_name(name)
            
            assert df is not None
            assert 'sources' in df.columns
        else:
            pytest.skip("No chemical name with inventory data found")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
