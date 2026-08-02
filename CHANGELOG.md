# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **`Search.enrich(df, column)`** — attach resolved identifier columns to a
  DataFrame, searching each *distinct* value once and broadcasting the result to
  every row that carries it. Row order and index are preserved; added columns are
  namespaced with a configurable `prefix` (default `provesid_`).
- **`resolve_cascade(df, stages, ...)`** — resolve rows through an ordered list of
  `Search` stages, passing each stage only the rows still unresolved, so every row
  is resolved by the most reliable identifier it actually has. Records
  `resolved_by` and `validated_by` per row, with an optional RDKit fallback that
  derives identifiers from the row's own structure.
- **`mw_within(tolerance, reference_column=...)`** — validator factory for
  `resolve_cascade`'s `accept` argument: accepts a hit only when its molecular
  weight agrees with the structure the dataset already carried, and reports any
  additional SMILES/name agreement. This is what stops a confident-but-wrong
  identifier match from ending a cascade.
- `CheMBL.search_by_name(..., exact=True)` for exact (case-insensitive) name and
  synonym matching.
- `ZeroPM.match_similar_name()` and `ZeroPM.get_id_table_from_similar_name()` —
  fuzzy name matching that reports *what* matched and *how well*, instead of
  discarding it.

### Fixed
- **`Search` resolved misspelled names to unrelated compounds and labelled them
  exact matches.** `Search("name", fuzzy=True).search("asprin")` returned
  PHENYRAMIDOL with `match_method="exact_name"`. Three defects combined:
  - `CheMBL.search_by_name` had no exact mode, so `Search`'s exact pass received a
    substring match — PHENYRAMIDOL carries the synonym `"Evasprin"`, which
    contains `"asprin"` — and tagged it `exact_name`.
  - The "did the exact pass find a strong match?" test used a fuzzy score, and
    `WRatio("asprin", "Evasprin")` is 85.7, clearing the cut-off. That suppressed
    the fuzzy widening which would have found aspirin. The test now requires an
    actual (case- and whitespace-insensitive) name equality.
  - A fuzzy match's confidence base was the raw similarity (up to 1.0) while an
    exact name match was pinned at 0.80, so **a typo could score higher than the
    correct spelling**. The fuzzy base is now scaled by the exact-name base, so an
    approximate match can never outrank an exact one.
- `Search`'s ZeroPM fuzzy branch was dead code: it acted only on a `DataFrame`,
  but `query_similar_name` returns a list of ids, so ZeroPM — the only source
  doing true fuzzy *retrieval* — never contributed to fuzzy name search.
- ChEMBL was queried in `Search`'s exact name pass but omitted from the fuzzy
  widening pass; it now participates in both.
- **`ZeroPM.get_pm_probabilities`, `batch_get_pm_probabilities` and
  `get_all_zeropm_chemicals` could never return P/M probability data.** All three
  keyed `pm_probabilities` on `zeropm_id`, but that table is keyed on `inchi_id`,
  so every call raised `sqlite3.OperationalError: no such column: zeropm_id`.
  `get_pm_probabilities` now translates a `zeropm_id` via the new
  `zeropm_id_to_inchi_id`, and the two joins use `inchi_id`.
- **`ChebifierClassifier` could not build the default ensemble at all.**
  `chemlog_extra` reads its element-class mapping files from a *working-directory
  relative* path and rebuilds them from the ChEBI graph when absent — but that
  rebuild crashes, because 288 of the graph's 205k nodes carry `name: None` and
  the builder does `" molecular entity" in properties["name"]`. Every
  `classify()` call failed with `TypeError: argument of type 'NoneType' is not
  iterable`. The new `ensure_element_class_mappings` writes both files into the
  PROVESID chebifier data directory using upstream's own derivation rules
  (skipping unnamed nodes), and the ensemble is now constructed with that
  directory as the working directory.

### Added (continued)
- `ZeroPM.zeropm_id_to_inchi_id()` — the reverse of `get_zeropm_id`.
- `taxonomy.ensure_element_class_mappings()`, `default_ensemble_available()` and
  `missing_ensemble_modules()`. The last two report on the *whole* default
  ensemble (transformer, graph, rule-based and c3p models each live in a separate
  package), so a partial install can be detected up front instead of failing with
  a bare `ModuleNotFoundError` at predict time.

### Fixed (tests)
- `test_pubchem_id.py::test_init_nonexistent_path` omitted `auto_download=False`,
  so instead of asserting `FileNotFoundError` it **downloaded the ~2.3 GB
  database into the repository root** on every run. It now passes the flag and
  uses `tmp_path`.
- `test_zeropm.py::test_get_cas_from_name_integration` asserted a
  name → CAS → same-CAS round trip that the data model does not support (names
  are many-to-many with CAS numbers), and its fallback clause compared a CAS
  against a list of *names*, so it could never hold. It passed only because its
  `SELECT ... LIMIT 1` had no `ORDER BY`: creating indexes — which another test in
  the file does to the shared database — changed which row came back and broke
  it. Replaced with a deterministic test of the guarantee that does hold.
- The three P/M probability tests repeated the same wrong join key in their own
  SQL, and their `else: pytest.skip(...)` branches turned the resulting error
  into a silent skip — which is how the production bug survived. They now assert
  that the fixture query found data, and check the returned values against the
  database so a wrong join cannot pass again.
- `test_zeropm.py` asserted `dtype == object` for string columns; pandas 3 gives
  `StringDtype`. Now uses `pandas.api.types.is_string_dtype`.
- `test_search.py::test_exact_inchikey_with_low_consensus` expected 0.5 for a
  zero consensus score. Zero consensus only occurs when no source matched at all
  (one source scores 1.0; two fully disagreeing sources score 0.5), so 0.0 is
  correct and the test encoded the un-special-cased formula. Replaced with tests
  for both the zero and the partial-agreement cases.
### Removed (tests)
- The end-to-end chebifier classification tests (`TestLiveClassification`).
  chebifier stays an optional extra, and its full model stack (transformer,
  graph/GNN, rule-based and c3p models — separate packages, some git-only) is
  awkward enough to install that the test suite should not depend on it. The
  remaining `test_taxonomy.py` tests all pass with the extra absent, verified by
  running them with the stack made unimportable. The now-unused `chebifier`
  pytest marker was dropped from `pyproject.toml`.

### Changed
- **`Search`'s default `fuzzy_scorer` is now `"ratio"`, was `"WRatio"`.** `WRatio`
  includes a partial-ratio term that scores a short name highly whenever it
  appears anywhere inside the query, which makes `fuzzy_score_cutoff` stop
  discriminating: `WRatio("caffiene", "ne")` and
  `WRatio("zzzznotachemical", "Mica")` are both 90, while `ratio` puts both at 40
  and still scores the genuine typo `ratio("caffiene", "caffeine")` at 87.5. Pass
  `fuzzy_scorer="WRatio"` to restore the old behaviour.
  Confidence values for `fuzzy_name` matches change as a result.
- `ZeroPM` "not found in database" messages moved from `WARNING` to `DEBUG` and
  onto the instance logger. They fired on the root logger during successful
  `Search` runs.

### Removed
- GitHub Actions workflows for running tests (`test.yml`,
  `test-with-api-keys.yml`) and for releasing (`release.yml`). Tests are run
  locally; releases are made manually with `twine`. Only the documentation
  deploy workflow remains.

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