# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.0] - 2026-04-16

### Added
- **`Search` class** — unified cross-database chemical identifier resolver (`provesid.Search`)
  - Accepts CAS, name, SMILES, InChI, InChIKey, DTXSID, or molecular formula as input
  - Queries all five offline databases (ChEBI, CompTox, PubChemID, ZeroPM, ChEMBL) in a single call
  - Returns a `pandas.DataFrame` with 23 standardised columns including `confidence`, `match_method`, `source_details`, and `source_match_scores`
  - Optional **salt stripping** (`strip_salts=True`): uses RDKit `SaltRemover` + largest-fragment picker; parent SMILES and InChIKey stored in dedicated columns
  - Optional **fuzzy name matching** (`fuzzy=True`): RapidFuzz-based candidate search with configurable similarity threshold
  - Optional **Tanimoto similarity** search (`similarity_threshold>0`): Morgan fingerprint fallback when no exact match is found
  - Optional **InChIKey skeleton** search (`inchikey_skeleton=True`): matches stereoisomers via 14-character skeleton prefix
  - Accepts `str`, `list[str]`, `pd.DataFrame` (with `column=` kwarg), or a CSV/Parquet file path as input
  - Confidence scoring model: base per match method × (0.5 + 0.5 × cross-source consensus score)
  - Full per-source traceability in `source_details` column
- **`normalize_structure(smiles)`** — RDKit helper returning canonical SMILES, Kekulized SMILES, InChI, InChIKey, and molecular weight
- **`strip_salts(smiles)`** — standalone salt-stripping utility exported from `provesid`
- **`OUTPUT_COLUMNS`** — list of all 23 column names in the `Search` result schema, exported from `provesid`
- Example scripts: `examples/search/search_by_cas_demo.py`, `search_by_name_demo.py`, `salt_stripping_demo.py`, `similarity_search_demo.py`
- API documentation page `docs/api/search.md`

### Fixed
- `strip_salts`: when `SaltRemover` strips all fragments, fall back to the largest fragment of the original molecule instead of returning an empty string

## [0.2.0] - 2025-09-29

### Added
- **🚀 Unlimited Caching System**: Complete overhaul of caching infrastructure
  - Unlimited cache by default (no size limits)
  - Persistent cache storage across sessions
  - 5GB warning threshold with configurable monitoring
  - Memory + disk hybrid caching for optimal performance
  - SHA256 cache keys for security and uniqueness
  - Import/Export functionality for team collaboration (pickle and JSON formats)
  - Global cache management functions: `clear_cache()`, `get_cache_info()`, `export_cache()`, `import_cache()`

- **📊 Comprehensive API Caching**: All major APIs now support unlimited caching
  - **PubChemAPI**: 19 cached methods including `get_compounds()`, `get_properties()`, `get_synonyms()`, etc.
  - **CASCommonChem**: 2 cached methods (`cas_to_detail()`, `name_to_detail()`)
  - **NCIChemicalIdentifierResolver**: 15 cached methods including all convenience functions
  - **PubChemView**: 15+ cached methods for experimental property extraction
  - **ClassyFireAPI**: 3 cached methods (`submit_query()`, `query_status()`, `get_query()`)
  - **OPSIN**: 2 cached methods (`get_id()`, `get_id_from_list()`)

- **🔧 Cache Management Methods**: Each API class now includes:
  - `clear_cache()`: Clear cached data for that specific API
  - `get_cache_info()`: Get detailed cache statistics and information

- **📈 Performance Improvements**:
  - Significant speed improvements for repeated API calls
  - Reduced API rate limiting issues
  - Offline capability when APIs are unavailable
  - Cross-session data persistence

### Changed
- **Breaking**: Removed `cache_size` parameter from PubChemAPI constructor (now unlimited by default)
- **Breaking**: Replaced all `@lru_cache(maxsize=X)` decorators with unlimited `@cached` decorator
- Cache behavior is now consistent across all APIs
- Cache storage moved from memory-only to persistent disk storage

### Enhanced
- **Test Coverage**: Added comprehensive caching tests
  - 168 tests passing with new caching system
  - Cache persistence, import/export, and size monitoring tests
- **Documentation**: Updated API documentation to reflect caching capabilities
- **Error Handling**: Improved cache error handling and recovery

### Technical Details
- New `cache.py` module with `CacheManager` class
- Automatic cache directory creation in system temp folder
- Thread-safe cache operations
- Configurable warning thresholds and cache policies
- Full backward compatibility (existing code continues to work)

### Performance Metrics
- Cache hit rates: Near 100% for repeated identical requests
- Memory usage: Efficient hybrid memory/disk storage
- Disk usage: Automatic monitoring with configurable warnings
- Speed improvement: 10-100x faster for cached requests

## [0.1.0] - Initial Release

### Added
- Initial implementation of PROVESID package
- PubChemAPI for PubChem REST API access
- CASCommonChem for CAS Common Chemistry API
- NCIChemicalIdentifierResolver for NCI resolver
- PubChemView for experimental properties
- ClassyFireAPI for chemical classification
- OPSIN for IUPAC name to structure conversion
- ChEBI API integration
- Basic caching with lru_cache (limited size)
- Comprehensive test suite
- Documentation and examples

[0.2.0]: https://github.com/USEtox/PROVESID/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/USEtox/PROVESID/releases/tag/v0.1.0