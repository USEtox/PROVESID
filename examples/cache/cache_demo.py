#!/usr/bin/env python3
"""
Example demonstrating the advanced cache features of PROVESID.

This script shows:
1. Unlimited caching with persistent storage
2. Cache size monitoring and warnings
3. Cache export/import functionality
4. Cache management utilities

Run this script to see the cache system in action.

The demo clears caches and lowers the size-warning threshold, so it redirects
PROVESID_CACHE_DIR into a throwaway directory first --- the real cache is
persistent now (``~/.cache/provesid`` and friends, no longer the system temp
directory), and a demo has no business deleting responses you paid network
time for. See ``cache_layout_and_versioning_demo.py`` for the layout itself.
"""

import os
import tempfile

_SANDBOX = tempfile.TemporaryDirectory(prefix="provesid-cache-demo-")
os.environ["PROVESID_CACHE_DIR"] = _SANDBOX.name

import provesid  # noqa: E402
import time  # noqa: E402

def main():
    print("🚀 PROVESID Advanced Cache System Demo")
    print("=" * 50)
    
    # 1. Demonstrate unlimited caching
    print("\n1. Creating PubChem API instance with unlimited cache...")
    api = provesid.PubChemAPI()
    
    # Clear any existing cache for clean demo
    print("   Clearing existing cache...")
    provesid.clear_cache()
    
    # Show initial cache state
    cache_info = provesid.get_cache_info()
    print(f"   Initial cache: {cache_info['memory_entries']} entries, {cache_info['total_size_mb']:.2f} MB")
    
    # 2. Demonstrate caching with mock functions (since we might not have network)
    print("\n2. Demonstrating caching behavior...")
    
    from provesid.cache import cached
    
    # Use a mutable container to track call count (works with nested functions)
    call_stats = {'count': 0}
    
    @cached
    def expensive_calculation(n):
        """Simulate an expensive API call or calculation."""
        call_stats['count'] += 1
        print(f"      🔄 Actually computing for n={n} (call #{call_stats['count']})")
        time.sleep(0.1)  # Simulate work
        return n * n + n
    
    # First calls - will be cached
    print("   Making function calls...")
    result1 = expensive_calculation(5)
    result2 = expensive_calculation(10)
    result3 = expensive_calculation(5)  # This should be cached!
    
    print(f"   Results: {result1}, {result2}, {result3}")
    print(f"   Total actual function calls: {call_stats['count']} (should be 2, not 3)")
    
    # 3. Show cache statistics
    print("\n3. Cache statistics...")
    print("   (these are the *global* cache: @cached with no service= writes here)")
    for name, info in provesid.get_all_cache_info().items():
        print(f"   {name:12s} {info['disk_entries']:3d} entries")
    cache_info = provesid.get_cache_info()
    print(f"   Cache entries: {cache_info['memory_entries']} in memory, {cache_info['disk_entries']} on disk")
    print(f"   Cache size: {cache_info['total_size_bytes']} bytes ({cache_info['total_size_mb']:.6f} MB)")
    print(f"   Cache directory: {cache_info['cache_directory']}")
    
    # 4. Demonstrate cache size monitoring
    print("\n4. Cache size monitoring...")
    print(f"   Warning threshold: {cache_info['warning_threshold_gb']} GB")
    
    # Set a very low threshold to demonstrate warning
    print("   Setting warning threshold to 0.001 GB to trigger warning...")
    provesid.set_cache_warning_threshold(0.001)
    
    # Add more data to trigger warning
    for i in range(10):
        expensive_calculation(i + 20)
    
    # 5. Demonstrate cache export
    print("\n5. Cache export/import...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        export_file = os.path.join(temp_dir, "cache_backup.pkl")
        
        # Export cache
        print(f"   Exporting cache to: {export_file}")
        success = provesid.export_cache(export_file)
        if success:
            file_size = os.path.getsize(export_file)
            print(f"   ✅ Export successful! File size: {file_size} bytes")
        else:
            print("   ❌ Export failed!")
        
        # Clear cache
        print("   Clearing cache...")
        provesid.clear_cache()
        cache_info = provesid.get_cache_info()
        print(f"   Cache after clear: {cache_info['memory_entries']} entries")
        
        # Import cache
        print("   Importing cache from backup...")
        success = provesid.import_cache(export_file)
        if success:
            cache_info = provesid.get_cache_info()
            print(f"   ✅ Import successful! Restored {cache_info['disk_entries']} entries")
        else:
            print("   ❌ Import failed!")
    
    # 6. Test that imported cache works
    print("\n6. Testing restored cache...")
    call_count_before = call_stats['count']
    result = expensive_calculation(5)  # Should be cached from import
    if call_stats['count'] == call_count_before:
        print(f"   ✅ Cache working! Got {result} without function call")
    else:
        print(f"   ⚠️  Cache miss - function was called again")
    
    # 7. Cache management utilities
    print("\n7. Cache management...")
    
    # Get cache size breakdown
    size_info = provesid.get_cache_size()
    print(f"   Detailed size info:")
    print(f"   - {size_info['files']} files")
    print(f"   - {size_info['bytes']} bytes")
    print(f"   - {size_info['mb']:.6f} MB")
    print(f"   - {size_info['gb']:.9f} GB")
    
    # Final cache info
    print("\n📊 Final Cache Statistics:")
    final_info = provesid.get_cache_info()
    for key, value in final_info.items():
        if isinstance(value, float):
            print(f"   {key}: {value:.6f}")
        else:
            print(f"   {key}: {value}")
    
    print("\n🎉 Demo completed!")
    print("\n💡 Key Benefits:")
    print("   • Unlimited cache size (no more 512-entry limits)")
    print("   • Persistent storage across sessions")
    print("   • Automatic size monitoring with 5GB warning")
    print("   • Export/import for backup and sharing")
    print("   • Zero-configuration - just import and use!")
    
    print("\n📖 Usage in your code:")
    print("   import provesid")
    print("   api = provesid.PubChemAPI()  # Unlimited cache automatically enabled")
    print("   provesid.export_cache('my_cache.pkl', service='pubchem')")
    print("   provesid.import_cache('shared_cache.pkl', service='pubchem')")
    print("   provesid.clear_cache(service='pubchem')  # or all_services=True")


if __name__ == "__main__":
    try:
        main()
    finally:
        _SANDBOX.cleanup()