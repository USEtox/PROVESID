# Advanced Caching in PROVESID

PROVESID now features an advanced caching system with unlimited storage, persistent caching across sessions, size monitoring, and import/export functionality.

## Key Features

### 🚀 Unlimited Caching
- **No more 512-entry limits**: Cache as many API calls as you need
- **Persistent storage**: Cache survives restarts and reinstalls
- **Automatic**: Zero configuration required - just import and use

### 📊 Size Monitoring
- **Smart warnings**: Get notified when cache exceeds 5GB (configurable)
- **Size tracking**: Monitor cache size in bytes, MB, and GB
- **File counting**: Track number of cached entries

### 💾 Export/Import
- **Backup your cache**: Export to pickle or JSON files
- **Share cache files**: Import cache data from shared files
- **Team collaboration**: Share expensive API results with team members
- **Offline mode**: Use cached data when APIs are down

### 🛠️ Cache Management
- **Clear when needed**: Remove all cached data
- **Get statistics**: Detailed cache information
- **Configure warnings**: Adjust size thresholds
- **Enable/disable monitoring**: Control warning behavior

## Quick Start

```python
import provesid

# All APIs now use unlimited caching automatically
pubchem_api = provesid.PubChemAPI()
nci_resolver = provesid.NCIChemicalIdentifierResolver()
cas_api = provesid.CASCommonChem()

# All API calls are cached forever
result1 = pubchem_api.get_compound_by_cid(2244)  # Cached
result2 = nci_resolver.resolve('aspirin', 'smiles')  # Cached  
result3 = cas_api.cas_to_detail('50-00-0')  # Cached

# Cache management works across all APIs
info = provesid.get_cache_info()
print(f"Cache size: {info['total_size_mb']:.2f} MB")

# Export your valuable cache (includes all API data)
provesid.export_cache('my_research_cache.pkl')

# Import shared cache (works for all APIs)
provesid.import_cache('shared_cache.pkl')

# Clear when needed (clears all API caches)
provesid.clear_cache()
```

## Migration from Previous Versions

**Before (limited cache):**
```python
api = PubChemAPI(cache_size=512)  # Limited to 512 entries
```

**Now (unlimited cache):**
```python
api = PubChemAPI()  # Unlimited cache automatically
# No cache_size parameter needed!
```

Your existing code will work unchanged, but now with unlimited caching!

## Cache Functions Reference

### `provesid.get_cache_info() -> dict`
Get comprehensive cache statistics:
```python
info = provesid.get_cache_info()
print(info)
# {
#     'cache_directory': '/path/to/cache',
#     'memory_entries': 42,
#     'disk_entries': 42, 
#     'total_size_bytes': 1048576,
#     'total_size_mb': 1.0,
#     'total_size_gb': 0.001,
#     'file_count': 42,
#     'warning_threshold_gb': 5.0,
#     'warnings_enabled': True
# }
```

### `provesid.get_cache_size() -> dict`
Get detailed size information:
```python
size = provesid.get_cache_size()
print(f"Cache: {size['mb']:.2f} MB ({size['files']} files)")
```

### `provesid.export_cache(path, format='pickle') -> bool`
Export cache to file:
```python
# Export as pickle (recommended)
success = provesid.export_cache('cache_backup.pkl')

# Export as JSON (human-readable, but limited data types)
success = provesid.export_cache('cache_backup.json', format='json')
```

### `provesid.import_cache(path, merge=True) -> bool`
Import cache from file:
```python
# Merge with existing cache
success = provesid.import_cache('cache_backup.pkl')

# Replace existing cache
success = provesid.import_cache('cache_backup.pkl', merge=False)
```

### `provesid.clear_cache()`
Clear all cached data:
```python
provesid.clear_cache()
```

### `provesid.set_cache_warning_threshold(size_gb)`
Set size warning threshold:
```python
# Warn when cache exceeds 10 GB
provesid.set_cache_warning_threshold(10.0)
```

### `provesid.enable_cache_warnings(enabled=True)`
Enable/disable size warnings:
```python
# Disable warnings
provesid.enable_cache_warnings(False)

# Re-enable warnings  
provesid.enable_cache_warnings(True)
```

## Use Cases

### 1. Long Research Projects
```python
# Start your research project
api = provesid.PubChemAPI()

# Make expensive API calls - all cached automatically
for compound in my_compound_list:
    data = api.get_compound_by_cid(compound)
    properties = api.get_compound_properties(compound, ['MolecularWeight', 'LogP'])

# Export cache at end of day
provesid.export_cache('research_day1.pkl')

# Next day: import and continue
provesid.import_cache('research_day1.pkl')
# All previous calls are cached!
```

### 2. Team Collaboration
```python
# Team member 1: Gather data from multiple APIs
pubchem_api = provesid.PubChemAPI()
nci_resolver = provesid.NCIChemicalIdentifierResolver()

for cid in expensive_compound_list:
    pubchem_api.get_compound_by_cid(cid)
    
for cas in cas_number_list:
    nci_resolver.get_molecular_data(cas)

# Share the cache (includes data from all APIs)
provesid.export_cache('team_shared_cache.pkl')

# Team member 2: Use shared data
provesid.import_cache('team_shared_cache.pkl')
# Instant access to all the data without API calls!
```

### 3. Offline Development
```python
# When online: gather data
api = provesid.PubChemAPI()
test_data = [api.get_compound_by_cid(cid) for cid in test_compounds]
provesid.export_cache('offline_cache.pkl')

# When offline: use cached data
provesid.import_cache('offline_cache.pkl')
api = provesid.PubChemAPI()
# All test compounds available from cache
result = api.get_compound_by_cid(2244)  # Works offline!
```

### 4. Cache Monitoring
```python
import provesid

# Monitor cache size during processing
def process_compounds(compound_list):
    api = provesid.PubChemAPI()
    
    for i, compound in enumerate(compound_list):
        result = api.get_compound_by_cid(compound)
        
        # Check cache size every 100 compounds
        if i % 100 == 0:
            size = provesid.get_cache_size()
            print(f"Processed {i} compounds, cache: {size['mb']:.1f} MB")
            
            # Export backup every 1000 compounds
            if i % 1000 == 0 and i > 0:
                provesid.export_cache(f'backup_{i}.pkl')
```

## Performance Benefits

### Before (Limited Cache)
```
First 512 calls: Fast (cached)
Call 513+: Slow (cache full, LRU eviction)
After restart: Slow (cache lost)
```

### Now (Unlimited Cache)
```
All calls: Fast after first time
After restart: Fast (persistent storage)
Across sessions: Fast (cache preserved)
Team sharing: Instant (import cache)
```

## Cache Storage Location

Cache files are stored in:
- **Windows**: `%TEMP%\provesid_cache\`
- **macOS/Linux**: `/tmp/provesid_cache/`

Each cached API call is stored as a separate file with metadata tracking.

## Best Practices

### 1. Regular Backups
```python
# Export cache regularly during long-running processes
if batch_count % 10 == 0:
    provesid.export_cache(f'backup_batch_{batch_count}.pkl')
```

### 2. Share Team Caches
```python
# At end of data collection phase
provesid.export_cache('project_phase1_cache.pkl')
# Share this file with team members
```

### 3. Monitor Size
```python
# Check cache size for large projects
size = provesid.get_cache_size()
if size['gb'] > 2.0:
    print(f"Large cache: {size['gb']:.2f} GB - consider archiving")
```

### 4. Clean Up When Done
```python
# Clear cache when switching projects
provesid.clear_cache()
```

## Troubleshooting

### Cache Warnings
If you see cache size warnings:
```python
# Option 1: Increase threshold
provesid.set_cache_warning_threshold(10.0)  # 10 GB

# Option 2: Export and clear
provesid.export_cache('archive.pkl')
provesid.clear_cache()

# Option 3: Disable warnings
provesid.enable_cache_warnings(False)
```

### Import/Export Failures
```python
# Always check return values
success = provesid.export_cache('backup.pkl')
if not success:
    print("Export failed - check disk space and permissions")

success = provesid.import_cache('backup.pkl') 
if not success:
    print("Import failed - check file exists and is valid")
```

### Cache Location Issues
```python
# Check cache location
info = provesid.get_cache_info()
print(f"Cache directory: {info['cache_directory']}")

# Verify directory is writable
import os
cache_dir = info['cache_directory']
print(f"Directory writable: {os.access(cache_dir, os.W_OK)}")
```

## Technical Details

### Cache Implementation
- **Storage**: Pickle serialization for Python objects
- **Indexing**: SHA256 hashes of function calls + arguments
- **Metadata**: JSON tracking for size and timestamps
- **Memory**: LRU memory cache backed by persistent disk storage

### Security Considerations
- Cache files contain API response data
- Exported cache files are not encrypted
- Consider security when sharing cache files
- Cache directory permissions follow system defaults

### Performance
- **Memory**: Fast lookup for recently accessed items
- **Disk**: Automatic persistence with minimal overhead
- **Network**: Eliminates repeated API calls entirely
- **Startup**: Quick loading of existing cache metadata