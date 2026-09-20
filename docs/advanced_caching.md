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

# Each service has its own cache; address one with service=
for name, info in provesid.get_all_cache_info().items():
    print(f"{name:12s} {info['total_size_mb']:8.2f} MB")

# Export your valuable PubChem results
provesid.export_cache('my_research_cache.pkl', service='pubchem')

# Import a colleague's
provesid.import_cache('shared_cache.pkl', service='pubchem')

# Clear when needed
provesid.clear_cache(service='pubchem')   # one service
provesid.clear_cache(all_services=True)   # all of them
```

## One Cache per Service

Each online service gets its own directory, so a PubChem cache can be cleared,
sized, exported or shared without touching the others:

```python
import provesid

provesid.CACHE_SERVICES
# ('pubchem', 'cas', 'nci', 'pubchemview', 'classyfire', 'opsin', 'chebifier')
```

Every cache function takes the same optional `service=` argument. Omit it and
you address the *global* cache, which holds entries from `@cached` functions
that name no service:

```python
provesid.get_cache_info(service='pubchem')   # one service
provesid.get_cache_info()                    # the global cache
provesid.get_all_cache_info()                # all of them, keyed by name
provesid.clear_cache(all_services=True)      # everything, global included
```

A service name that is not in `CACHE_SERVICES` raises `ValueError` rather than
falling back to the global cache — a typo would otherwise write entries to a
place the matching `clear_cache` call never looks.

## Cache Functions Reference

### `provesid.get_cache_info(service=None) -> dict`
Get comprehensive cache statistics:
```python
info = provesid.get_cache_info(service='pubchem')
print(info)
# {
#     'cache_directory': '/home/you/.cache/provesid/pubchem',
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

### `provesid.get_all_cache_info() -> dict`
The same, for the global cache and every service at once:
```python
for name, info in provesid.get_all_cache_info().items():
    print(f"{name:12s} {info['total_size_mb']:8.2f} MB")
```

### `provesid.get_cache_size(service=None) -> dict`
Get detailed size information:
```python
size = provesid.get_cache_size(service='pubchem')
print(f"Cache: {size['mb']:.2f} MB ({size['files']} files)")
```

### `provesid.export_cache(path, format='pickle', service=None) -> bool`
Export cache to file:
```python
# Export as pickle (recommended)
success = provesid.export_cache('cache_backup.pkl', service='pubchem')

# Export as JSON (human-readable, but limited data types)
success = provesid.export_cache('cache_backup.json', format='json')
```

### `provesid.import_cache(path, merge=True, service=None) -> bool`
Import cache from file:
```python
# Merge with existing cache
success = provesid.import_cache('cache_backup.pkl', service='pubchem')

# Replace existing cache
success = provesid.import_cache('cache_backup.pkl', merge=False)
```

### `provesid.clear_cache(service=None, all_services=False)`
Clear cached data:
```python
provesid.clear_cache(service='pubchem')   # one service
provesid.clear_cache()                    # the global cache only
provesid.clear_cache(all_services=True)   # every cache
```

### `provesid.set_cache_warning_threshold(size_gb, service=None)`
Set size warning threshold:
```python
# Warn when the PubChem cache exceeds 10 GB
provesid.set_cache_warning_threshold(10.0, service='pubchem')
```

### `provesid.enable_cache_warnings(enabled=True, service=None)`
Enable/disable size warnings:
```python
# Disable warnings
provesid.enable_cache_warnings(False, service='pubchem')

# Re-enable warnings  
provesid.enable_cache_warnings(True, service='pubchem')
```

### `provesid.get_service_cache(service=None) -> CacheManager`
The underlying manager, for code that needs `get`/`set` directly (as
`ChebifierClassifier` does with its InChIKey-keyed entries).

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
provesid.export_cache('research_day1.pkl', service='pubchem')

# Next day: import and continue
provesid.import_cache('research_day1.pkl', service='pubchem')
# All previous calls are cached!
```

Nothing is lost if you skip the export — the cache is already persistent. The
export is for moving entries to another machine or archiving them.

### 2. Team Collaboration
```python
# Team member 1: gather data from two services
pubchem_api = provesid.PubChemAPI()
nci_resolver = provesid.NCIChemicalIdentifierResolver()

for cid in expensive_compound_list:
    pubchem_api.get_compound_by_cid(cid)
    
for cas in cas_number_list:
    nci_resolver.get_molecular_data(cas)

# One file per service: an export covers one cache, not all of them
provesid.export_cache('team_pubchem.pkl', service='pubchem')
provesid.export_cache('team_nci.pkl', service='nci')

# Team member 2: use shared data
provesid.import_cache('team_pubchem.pkl', service='pubchem')
provesid.import_cache('team_nci.pkl', service='nci')
# Instant access to all the data without API calls!
```

### 3. Offline Development
```python
# When online: gather data
api = provesid.PubChemAPI()
test_data = [api.get_compound_by_cid(cid) for cid in test_compounds]
provesid.export_cache('offline_cache.pkl', service='pubchem')

# On another machine: import and work offline
provesid.import_cache('offline_cache.pkl', service='pubchem')
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
            size = provesid.get_cache_size(service='pubchem')
            print(f"Processed {i} compounds, cache: {size['mb']:.1f} MB")
            
            # Export backup every 1000 compounds
            if i % 1000 == 0 and i > 0:
                provesid.export_cache(f'backup_{i}.pkl', service='pubchem')
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
After reboot: Fast (a real cache directory, not /tmp)
Team sharing: Instant (import cache)
```

## Cache Storage Location

Cache files live in the per-user cache directory that `platformdirs` picks for
the platform, with one subdirectory per service:

- **Linux**: `~/.cache/provesid/<service>/`
- **macOS**: `~/Library/Caches/provesid/<service>/`
- **Windows**: `%LOCALAPPDATA%\USEtox\provesid\Cache\<service>\`

Set `PROVESID_CACHE_DIR` to put them somewhere else — a scratch disk, or a
throwaway directory for a test run:

```bash
export PROVESID_CACHE_DIR=/scratch/provesid-cache
```

This is *not* the system temp directory, and that is the point: `/tmp` is
cleared on boot on most Linux systems, so the earlier default quietly discarded
the entire cache between sessions while this page promised it persisted.

Cached responses are disposable — everything here can be re-fetched — which is
why they sit under a different root from the bulk datasets
(`PROVESID_DATA_DIR`, see `provesid.datasets`) and can be deleted independently.

Each cached API call is stored as a separate file with metadata tracking.

## Best Practices

### 1. Regular Backups
```python
# Export cache regularly during long-running processes
if batch_count % 10 == 0:
    provesid.export_cache(f'backup_batch_{batch_count}.pkl', service='pubchem')
```

### 2. Share Team Caches
```python
# At end of data collection phase
provesid.export_cache('project_phase1_cache.pkl', service='pubchem')
# Share this file with team members
```

### 3. Monitor Size
```python
# Check cache size for large projects
size = provesid.get_cache_size(service='pubchem')
if size['gb'] > 2.0:
    print(f"Large cache: {size['gb']:.2f} GB - consider archiving")
```

### 4. Clean Up When Done
```python
# Clear cache when switching projects
provesid.clear_cache(all_services=True)
```

## Troubleshooting

### Cache Warnings
If you see cache size warnings:
```python
# Option 1: Increase threshold
provesid.set_cache_warning_threshold(10.0, service='pubchem')  # 10 GB

# Option 2: Export and clear
provesid.export_cache('archive.pkl', service='pubchem')
provesid.clear_cache(service='pubchem')

# Option 3: Disable warnings
provesid.enable_cache_warnings(False, service='pubchem')
```

### Import/Export Failures
```python
# Always check return values
success = provesid.export_cache('backup.pkl', service='pubchem')
if not success:
    print("Export failed - check disk space and permissions")

success = provesid.import_cache('backup.pkl', service='pubchem') 
if not success:
    print("Import failed - check file exists and is valid")
```

### Cache Location Issues
```python
# Check cache location
info = provesid.get_cache_info(service='pubchem')
print(f"Cache directory: {info['cache_directory']}")

# Verify directory is writable
import os
cache_dir = info['cache_directory']
print(f"Directory writable: {os.access(cache_dir, os.W_OK)}")
```

## What Is and Is Not Cached

### Failed lookups are never stored

A cache that remembers errors is worse than no cache: one HTTP 429 or one
PUG-View `ServerBusy` would turn into a permanent "this compound has no data".
So `@cached` writes nothing when:

- the function raises — the exception propagates uncached; or
- the return value is an in-band failure report, i.e. a dict whose `success`
  key is `False` (`provesid.cache.is_failure_result` recognises these); or
- an optional `skip_if` predicate on the decorator accepts the value.

```python
from provesid import PubChemView

view = PubChemView()

# A transient PUG-View failure raises instead of returning an empty list,
# and nothing is written to the cache, so the next call retries.
data = view.extract_property_data(2244, "Melting Point")

# An empty list therefore has one meaning only: PubChem holds no such
# property for this compound.
assert isinstance(data, list)
```

`use_cache=False` means "do not *read* the cache". A successful result is still
written, so later calls benefit.

### Absence is not cached either

A lookup that legitimately finds nothing is also not stored. PubChem reports a
genuinely empty result with a `NotFound` fault code, but it has been observed
sending one while merely under load, and a wrongly cached "this compound has no
melting point" would be permanent. Re-checking an empty result costs one cheap
request, so that is the trade taken:

- positive results are cached (chemical data rarely changes);
- empty results are re-fetched every time;
- the raw response fetchers (`get_property`, `get_experimental_properties`)
  raise on absence, so they cache unconditionally.

### Cache keys are stable across processes

The key is a SHA-256 digest of the function's qualified name plus its
normalised arguments. Normalisation matters for methods: the first argument is
the client instance, and its default `repr` contains a memory address, which
would make every key unique to one object in one process. Clients therefore
declare their identity:

```python
class PubChemView:
    def __cache_key__(self):
        return ("provesid.pubchemview.PubChemView", self.base_url)
```

Two clients pointing at the same endpoint share cache entries; two different
endpoints keep separate ones. Any class whose cached methods should behave this
way can declare `__cache_key__`; without it, the class's fully qualified name is
used.

### Keys carry a version, so an upgrade cannot serve the wrong shape

Pickle stores an instance's `__dict__`. When a cached return value gains a
field, an entry written by the previous release unpickles into an object that
is *missing that attribute*, and every caller reading it raises — on a machine
where the only thing that changed was the package version. The cache cannot
detect this by itself: the key is the same, the file is readable, the value is
simply the wrong shape.

So three version numbers are folded into the key, and the right one to bump is
the narrowest that covers the change:

| Bump | Retires | Where |
|---|---|---|
| `@cached(version=N)` | one function's entries | the decorator |
| `Client.CACHE_SCHEMA_VERSION` | one client's entries | returned from `__cache_key__` |
| `provesid.cache.CACHE_KEY_VERSION` | every entry, everywhere | module constant |

```python
from provesid.cache import cached

@cached(service='pubchem', version=2)   # v1 entries are now unreachable
def fetch_summary(cid):
    ...
```

Retired entries are not deleted, just never read again; `clear_cache` reclaims
the space. There is no migration path and none is wanted — the value of a
cache entry is one avoided HTTP request, which is not worth the risk of
deserialising a stale shape.

## Technical Details

### Cache Implementation
- **Storage**: Pickle serialization for Python objects
- **Indexing**: SHA256 hashes of function calls + normalised arguments (see above)
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