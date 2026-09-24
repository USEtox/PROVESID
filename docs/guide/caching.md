# Caching

The online clients cache their successful answers on disk, so a question asked
once is not asked again: not later in the session, not after a restart, and
not by a colleague who imports your cache. The offline databases are not part
of this; they are installed data, covered in
[Installing the offline databases](datasets.md).

Caching is on by default and needs no set-up. It covers `PubChemAPI`,
`PubChemView`, `NCIChemicalIdentifierResolver`, `CASCommonChem`, `OPSIN` and
`ClassyFireAPI`, and the Chebifier classifier; the online `ChEBI` client is
not cached. `PubChemAPI`, `PubChemView`, `CASCommonChem` and `OPSIN` take
`use_cache=False`, and `ClassyFireAPI`'s calls take it per call. It means "do
not *read* the cache": a successful result is still written.

## One cache per service

Each service has its own directory, so one can be sized, cleared, exported or
shared without touching the others:

```python
import provesid

provesid.CACHE_SERVICES
# ('pubchem', 'cas', 'nci', 'pubchemview', 'classyfire', 'opsin', 'chebifier')

for name, info in provesid.get_all_cache_info().items():
    print(f"{name:12s} {info['total_size_mb']:8.2f} MB")

provesid.get_cache_info(service="pubchem")
# {'cache_directory': '/home/you/.cache/provesid/pubchem', 'memory_entries': 0,
#  'disk_entries': 6, 'total_size_bytes': 7493, 'total_size_mb': 0.0071, ...,
#  'file_count': 6, 'warning_threshold_gb': 5.0, 'warnings_enabled': True}
```

Every cache function takes the same optional `service=`. Omit it and you
address the *global* cache, which holds entries from `@cached` functions that
name no service. A name that is not in `CACHE_SERVICES` raises `ValueError`
rather than falling back to the global cache, since a typo would otherwise
write entries to a place the matching `clear_cache` never looks.

## Managing it

```python
# Share a service's cache: export it, and import it on another machine
provesid.export_cache("pubchem_cache.pkl", service="pubchem")
provesid.import_cache("pubchem_cache.pkl", service="pubchem")              # merge
provesid.import_cache("pubchem_cache.pkl", service="pubchem", merge=False) # replace

# Reclaim the space
provesid.clear_cache(service="pubchem")      # one service
provesid.clear_cache(all_services=True)      # every service and the global cache

# A UserWarning is issued when a cache passes 5 GB
provesid.set_cache_warning_threshold(10.0, service="pubchem")
provesid.enable_cache_warnings(False, service="pubchem")
```

`export_cache` and `import_cache` return `True` or `False` rather than raising.
An export covers one cache; `format="json"` writes a readable file, but only for
values JSON can hold. Exported files are plain pickles: import only files you
trust.

## Where it lives

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

This is deliberately not the system temp directory, which many Linux systems
clear on every boot.

Cached responses are disposable — everything here can be re-fetched — which is
why they sit under a different root from the bulk datasets
(`PROVESID_DATA_DIR`, see [Installing the offline databases](datasets.md)) and
can be deleted independently.

Each cached API call is stored as a separate file with metadata tracking.

## What is and is not cached

### Failures are never stored

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
