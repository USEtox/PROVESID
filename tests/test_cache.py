"""
Test the new advanced cache system for PROVESID package.
"""

import pytest
import tempfile
import os
from pathlib import Path
import warnings
import json

import provesid
from provesid.cache import CacheManager, cached


class TestCacheManager:
    """Test the CacheManager class functionality."""
    
    def test_initialization(self):
        """Test cache manager initialization."""
        # Test with default directory
        cache_mgr = CacheManager()
        assert cache_mgr.cache_dir.exists()
        assert cache_mgr.max_size_gb == 5.0
        assert cache_mgr.enable_warnings == True
        
        # Test with custom directory
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_mgr = CacheManager(cache_dir=temp_dir, max_size_gb=10.0, enable_warnings=False)
            assert str(cache_mgr.cache_dir) == temp_dir
            assert cache_mgr.max_size_gb == 10.0
            assert cache_mgr.enable_warnings == False
    
    def test_cache_key_generation(self):
        """Test cache key generation."""
        cache_mgr = CacheManager()
        
        # Same arguments should produce same key
        key1 = cache_mgr._get_cache_key("test_func", (1, 2), {"a": 3})
        key2 = cache_mgr._get_cache_key("test_func", (1, 2), {"a": 3})
        assert key1 == key2
        
        # Different arguments should produce different keys
        key3 = cache_mgr._get_cache_key("test_func", (1, 3), {"a": 3})
        assert key1 != key3
        
        # Different function names should produce different keys
        key4 = cache_mgr._get_cache_key("other_func", (1, 2), {"a": 3})
        assert key1 != key4
    
    def test_cache_operations(self):
        """Test basic cache get/set operations."""
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_mgr = CacheManager(cache_dir=temp_dir)
            
            # Test cache miss
            found, value = cache_mgr.get("test_func", (1, 2), {})
            assert not found
            assert value is None
            
            # Test cache set and hit
            test_value = {"result": "cached_data"}
            cache_mgr.set("test_func", (1, 2), {}, test_value)
            
            found, value = cache_mgr.get("test_func", (1, 2), {})
            assert found
            assert value == test_value
    
    def test_persistent_storage(self):
        """Test that cache persists across manager instances."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create first manager and cache some data
            cache_mgr1 = CacheManager(cache_dir=temp_dir)
            test_value = {"result": "persistent_data"}
            cache_mgr1.set("test_func", ("persistent",), {}, test_value)
            
            # Create second manager with same directory
            cache_mgr2 = CacheManager(cache_dir=temp_dir)
            found, value = cache_mgr2.get("test_func", ("persistent",), {})
            assert found
            assert value == test_value
    
    def test_cache_info(self):
        """Test cache information retrieval."""
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_mgr = CacheManager(cache_dir=temp_dir)
            
            # Initially empty
            info = cache_mgr.get_cache_info()
            assert info['memory_entries'] == 0
            assert info['disk_entries'] == 0
            assert info['total_size_bytes'] >= 0
            
            # Add some data
            cache_mgr.set("test_func", (1,), {}, "data1")
            cache_mgr.set("test_func", (2,), {}, "data2")
            
            info = cache_mgr.get_cache_info()
            assert info['memory_entries'] == 2
            assert info['disk_entries'] == 2
            assert info['total_size_bytes'] > 0
    
    def test_cache_size_monitoring(self):
        """Test cache size monitoring and warnings."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Set very low threshold to trigger warning
            cache_mgr = CacheManager(cache_dir=temp_dir, max_size_gb=0.000001, enable_warnings=True)
            
            # This should trigger a warning when size is checked
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                
                # Add some data that will exceed the tiny threshold
                large_data = "x" * 10000  # 10KB of data to ensure we exceed threshold
                cache_mgr.set("test_func", (1,), {}, large_data)
                cache_mgr.set("test_func", (2,), {}, large_data)  # Add more data
                
                # Save metadata to update size tracking
                cache_mgr._save_metadata()
                
                # Force a size check
                cache_mgr._check_cache_size()
                
                # Should have a warning
                assert len(w) > 0
                assert "Cache size" in str(w[0].message)
    
    def test_clear_cache(self):
        """Test cache clearing."""
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_mgr = CacheManager(cache_dir=temp_dir)
            
            # Add some data
            cache_mgr.set("test_func", (1,), {}, "data1")
            cache_mgr.set("test_func", (2,), {}, "data2")
            
            # Verify data exists
            info = cache_mgr.get_cache_info()
            assert info['memory_entries'] == 2
            assert info['disk_entries'] == 2
            
            # Clear cache
            cache_mgr.clear()
            
            # Verify cache is empty
            info = cache_mgr.get_cache_info()
            assert info['memory_entries'] == 0
            assert info['disk_entries'] == 0
    
    def test_export_import_cache(self):
        """Test cache export and import functionality."""
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir1 = os.path.join(temp_dir, "cache1")
            cache_dir2 = os.path.join(temp_dir, "cache2")
            export_file = os.path.join(temp_dir, "cache_export.pkl")
            
            # Create cache with some data
            cache_mgr1 = CacheManager(cache_dir=cache_dir1)
            test_data = {
                "data1": [1, 2, 3],
                "data2": {"nested": "value"}
            }
            cache_mgr1.set("func1", (1,), {}, test_data["data1"])
            cache_mgr1.set("func2", ("arg",), {"kwarg": "value"}, test_data["data2"])
            
            # Export cache
            success = cache_mgr1.export_cache(export_file)
            assert success
            assert os.path.exists(export_file)
            
            # Create new cache and import
            cache_mgr2 = CacheManager(cache_dir=cache_dir2)
            success = cache_mgr2.import_cache(export_file)
            assert success
            
            # Verify imported data
            found1, value1 = cache_mgr2.get("func1", (1,), {})
            assert found1
            assert value1 == test_data["data1"]
            
            found2, value2 = cache_mgr2.get("func2", ("arg",), {"kwarg": "value"})
            assert found2
            assert value2 == test_data["data2"]
    
    def test_json_export_import(self):
        """Test JSON format export/import."""
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = os.path.join(temp_dir, "cache")
            export_file = os.path.join(temp_dir, "cache_export.json")
            
            # Create cache with JSON-serializable data
            cache_mgr = CacheManager(cache_dir=cache_dir)
            test_data = {"simple": "value", "number": 42}
            cache_mgr.set("func", ("arg",), {}, test_data)
            
            # Export as JSON
            success = cache_mgr.export_cache(export_file, format='json')
            assert success
            assert os.path.exists(export_file)
            
            # Verify it's valid JSON
            with open(export_file, 'r') as f:
                exported_data = json.load(f)
                assert 'metadata' in exported_data
                assert 'cache_data' in exported_data
                assert 'export_info' in exported_data
    
    
class TestCachedDecorator:
    """Test the @cached decorator functionality."""
    
    def test_basic_caching(self):
        """Test basic function caching."""
        call_count = 0
        
        @cached
        def test_function(x, y):
            nonlocal call_count
            call_count += 1
            return x + y
        
        # First call
        result1 = test_function(1, 2)
        assert result1 == 3
        assert call_count == 1
        
        # Second call with same arguments - should be cached
        result2 = test_function(1, 2)
        assert result2 == 3
        assert call_count == 1  # Should not increment
        
        # Call with different arguments
        result3 = test_function(2, 3)
        assert result3 == 5
        assert call_count == 2  # Should increment
    
    def test_kwargs_caching(self):
        """Test caching with keyword arguments."""
        call_count = 0
        
        @cached
        def test_function(x, y=10, z=None):
            nonlocal call_count
            call_count += 1
            return x + y + (z or 0)
        
        # Test different call patterns
        result1 = test_function(1, y=2, z=3)
        assert result1 == 6
        assert call_count == 1
        
        # Same arguments in different order - should be cached
        result2 = test_function(1, z=3, y=2)
        assert result2 == 6
        assert call_count == 1
        
        # Different arguments - default z=None becomes 0
        result3 = test_function(1, y=5)
        assert result3 == 6  # 1 + 5 + 0 (z defaults to None, becomes 0)
        assert call_count == 2


class TestProvesidCacheIntegration:
    """Test integration with PROVESID modules."""
    
    def test_cache_functions_available(self):
        """Test that cache management functions are available."""
        # Test functions are imported
        assert hasattr(provesid, 'clear_cache')
        assert hasattr(provesid, 'get_cache_info')
        assert hasattr(provesid, 'export_cache')
        assert hasattr(provesid, 'import_cache')
        assert hasattr(provesid, 'get_cache_size')
        assert hasattr(provesid, 'set_cache_warning_threshold')
        assert hasattr(provesid, 'enable_cache_warnings')
    
    def test_global_cache_operations(self):
        """Test global cache operations."""
        # Clear any existing cache
        provesid.clear_cache()
        
        # Get initial cache info
        info = provesid.get_cache_info()
        assert 'total_size_bytes' in info
        assert 'memory_entries' in info
        
        # Get cache size
        size_info = provesid.get_cache_size()
        assert 'bytes' in size_info
        assert 'mb' in size_info
        assert 'gb' in size_info
    
    def test_cache_warning_threshold(self):
        """Test cache warning threshold management."""
        # Set warning threshold
        provesid.set_cache_warning_threshold(10.0)
        
        # Enable/disable warnings
        provesid.enable_cache_warnings(False)
        provesid.enable_cache_warnings(True)
    
    @pytest.mark.skipif(True, reason="Requires network access - integration test")
    def test_pubchem_caching(self):
        """Test that PubChem API uses caching (requires network)."""
        # This test would require network access and API calls
        # Skip for now but shows how to test caching with real API calls
        api = provesid.PubChemAPI()
        
        # Clear cache first
        api.clear_cache()
        
        # Make API call (would be cached)
        # result = api.get_compound_by_cid(2244)  # Aspirin
        
        # Get cache info
        info = api.get_cache_info()
        # assert info['memory_entries'] > 0
    
    def test_nci_resolver_caching(self):
        """Test that NCI resolver has cache methods."""
        resolver = provesid.NCIChemicalIdentifierResolver()
        
        # Test cache methods are available
        assert hasattr(resolver, 'clear_cache')
        assert hasattr(resolver, 'get_cache_info')
        
        # Test cache info works
        info = resolver.get_cache_info()
        assert 'memory_entries' in info


def test_cache_export_import_integration():
    """Integration test for cache export/import."""
    with tempfile.TemporaryDirectory() as temp_dir:
        export_file = os.path.join(temp_dir, "test_export.pkl")
        
        # Clear cache
        provesid.clear_cache()
        
        # Create some test cached functions
        @cached
        def test_func1(x):
            return x * 2
        
        @cached  
        def test_func2(x, y):
            return x + y
        
        # Call functions to populate cache
        result1 = test_func1(5)
        result2 = test_func2(3, 4)
        
        assert result1 == 10
        assert result2 == 7
        
        # Export cache
        success = provesid.export_cache(export_file)
        assert success
        
        # Clear cache
        provesid.clear_cache()
        
        # Import cache
        success = provesid.import_cache(export_file)
        assert success
        
        # Test that cached results are available (functions should not be called again)
        # Note: This is tricky to test because the functions are defined in this test
        # In a real scenario, you'd import the cache after restarting the application