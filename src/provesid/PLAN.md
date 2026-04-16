# PROVESID Search Class — Implementation Plan

## Overview

Extend the offline identifier resolution in `tools.py` with a unified `Search` class
that accepts any chemical identifier (CAS, name, SMILES, InChI, InChIKey, DTXSID,
molecular formula) and returns a standardised result DataFrame. The class adds
structure-aware matching, confidence scoring, fuzzy name search, Tanimoto similarity
search, InChIKey-skeleton matching, and salt/solvent stripping.

The existing standalone functions (`ids_from_CAS`, `ids_from_name`,
`ids_from_SMILES`, `casrn_to_compounds`, `iupac_name_to_id`) remain available as
thin wrappers around `Search` for backward compatibility.

---

## Phase 1 — `Search` class core

### 1.1 New module: `src/provesid/search.py`

Create a new `Search` class that:

- Accepts a **single identifier type** at init time via an enum-like parameter
  (`identifier_type`): `"cas"`, `"name"`, `"smiles"`, `"inchi"`, `"inchikey"`,
  `"dtxsid"`, `"formula"`.
- Holds references to lazy-initialized source clients (ChEBI, CompTox, PubChemID,
  ZeroPM, ChEMBL). Clients are initialized once and reused across queries.
- Exposes a single entry point: `search(queries, ...)` where `queries` can be:
  - A single string
  - A `list[str]`
  - A `pandas.DataFrame` (with a column name matching the identifier type)
  - A file path (CSV or Parquet) — read into a DataFrame internally
- Returns a `pandas.DataFrame` with the standardized output schema (see §1.2).

```python
class Search:
    """Unified chemical identifier resolver using offline databases.

    Attributes:
        identifier_type: The type of input identifier.
        strip_salts: Whether to strip salts/solvents and search parent molecule.
        fuzzy: Whether to enable fuzzy name matching.
        similarity_threshold: Minimum Tanimoto similarity for structure search.
        inchikey_skeleton: Whether to enable InChIKey first-block matching.
    """

    def __init__(
        self,
        identifier_type: str = "cas",
        *,
        strip_salts: bool = False,
        fuzzy: bool = False,
        similarity_threshold: float = 0.0,
        inchikey_skeleton: bool = False,
        show_progress: bool = True,
        chebi: ChebiSDF | None = None,
        comptox: CompToxID | None = None,
        pubchem: PubChemID | None = None,
        zeropm: ZeroPM | None = None,
        chembl: CheMBL | None = None,
    ): ...

    def search(
        self,
        queries: str | list[str] | pd.DataFrame | Path,
        *,
        column: str | None = None,
    ) -> pd.DataFrame: ...
```

### 1.2 Output schema

Every result row includes:

| Column              | Type            | Description                                         |
|---------------------|-----------------|-----------------------------------------------------|
| `query`             | `str`           | Original input value                                |
| `CASRN`             | `str \| None`   | CAS Registry Number                                 |
| `name`              | `str \| None`   | Preferred common name                               |
| `IUPAC_name`        | `str \| None`   | IUPAC systematic name                               |
| `molecular_formula` | `str \| None`   | Molecular formula                                   |
| `SMILES`            | `str \| None`   | Original SMILES from source                         |
| `canonical_smiles`  | `str \| None`   | RDKit-canonical SMILES                              |
| `kekulized_smiles`  | `str \| None`   | Kekulized SMILES from RDKit                         |
| `InChI`             | `str \| None`   | InChI string                                        |
| `InChIKey`          | `str \| None`   | Full InChIKey (always reported)                     |
| `DTXSID`            | `str \| None`   | CompTox substance identifier                        |
| `molecular_mass`    | `float \| None` | Molecular weight                                    |
| `Synonyms`          | `str \| None`   | Semicolon-separated synonyms                        |
| `parent_smiles`     | `str \| None`   | SMILES of parent (salt-stripped) molecule, if opt-in |
| `parent_inchikey`   | `str \| None`   | InChIKey of parent molecule, if opt-in               |
| `foundby`           | `str`           | How the match was found (e.g. `CASRN`, `name`, `inchikey_skeleton`, `tanimoto`) |
| `source`            | `str \| None`   | Data source that provided the primary SMILES         |
| `source_details`    | `dict`          | Per-source flags: which sources contributed data     |
| `confidence`        | `float`         | Overall confidence score in [0, 1]                   |
| `match_method`      | `str`           | Matching method used (see §2)                        |
| `match_score`       | `float`         | Cross-source consensus score in [0, 1]               |
| `consensus_source`  | `str \| None`   | Source chosen by consensus algorithm                 |
| `source_match_scores` | `dict`        | Per-source agreement scores                          |

### 1.3 Identifier-specific lookup dispatch

Internally `search()` dispatches to a private method per identifier type:

- `_resolve_cas(cas)` — wraps existing `ids_from_CAS` logic
- `_resolve_name(name)` — wraps existing `ids_from_name` logic + new fuzzy path
- `_resolve_smiles(smiles)` — wraps existing `ids_from_SMILES` logic
- `_resolve_inchi(inchi)` — new: convert to InChIKey via RDKit, then look up by InChIKey
- `_resolve_inchikey(inchikey)` — new: direct InChIKey lookup across all sources
- `_resolve_dtxsid(dtxsid)` — new: CompTox primary, cross-reference other sources
- `_resolve_formula(formula)` — new: search CompTox/PubChem by formula, rank by completeness

Each resolver returns the standardized result dict (matching the output schema) with
confidence and match_method populated.

### 1.4 Client initialization

Move the repeated client-init pattern from `casrn_to_compounds` into a shared
`_init_clients()` method on `Search`. Allow users to inject pre-initialized clients
(for testing and advanced use) or let the class create them lazily.

---

## Phase 2 — Confidence scoring

### 2.1 Confidence model

Assign confidence based on **match method** and **cross-source agreement**:

| Match method                  | Base confidence |
|-------------------------------|-----------------|
| Exact InChIKey match          | 1.0             |
| Exact canonical SMILES match  | 0.95            |
| Exact CAS lookup              | 0.90            |
| InChIKey skeleton (14-char)   | 0.75            |
| Exact name match              | 0.80            |
| Fuzzy name match              | `rapidfuzz.ratio / 100` (typically 0.2–0.9) |
| Tanimoto similarity           | Tanimoto score × 0.85 |
| Formula-only match            | 0.30            |

The base confidence is then modulated by the cross-source consensus score
(from existing `_compute_consensus`):

```
final_confidence = base_confidence × (0.5 + 0.5 × consensus_score)
```

This ensures that a high-confidence match method with zero cross-source agreement
still gets at least half credit, while full agreement brings it to the base value.

### 2.2 `match_method` field

Populated as a string: `"exact_inchikey"`, `"exact_smiles"`, `"exact_cas"`,
`"inchikey_skeleton"`, `"exact_name"`, `"fuzzy_name"`, `"tanimoto"`, `"formula"`.

---

## Phase 3 — Structure-aware matching

### 3.1 SMILES normalization

Add a `_normalize_structure(smiles)` utility that returns a dict:

```python
{
    "canonical_smiles": str | None,
    "kekulized_smiles": str | None,
    "inchi": str | None,
    "inchikey": str | None,
    "mol": Chem.Mol | None,   # kept in memory only, not serialized
}
```

- Canonicalization via `Chem.MolToSmiles(mol, canonical=True)`
- Kekulization via `Chem.Kekulize(mol)` then `Chem.MolToSmiles(mol, kekuleSmiles=True)`
- InChI/InChIKey generation via `Chem.MolToInchi` / `Chem.InchiToInchiKey`

This replaces the current `_smiles_to_canonical_and_mass` with a richer return value.
The existing function remains as a thin wrapper for backward compatibility within
`tools.py`.

### 3.2 InChIKey always reported

Ensure every result row has `InChIKey` populated whenever a SMILES or InChI is
available, by running the normalization pipeline. If the source provides an InChIKey,
validate it against the RDKit-derived one and log a warning on mismatch.

### 3.3 Tanimoto similarity search

New method `_similarity_search(query_smiles, threshold)`:

1. Compute Morgan fingerprint (radius=2, 2048 bits) for the query molecule.
2. For each source that exposes structure data (ChEBI, CompTox, PubChem, ChEMBL),
   compute fingerprints for candidate molecules and filter by
   `DataStructs.TanimotoSimilarity >= threshold`.
3. Return ranked list of candidates above threshold.

**Performance note:** Pre-computed fingerprint indices per source will be needed for
large-scale use. Initially, iterate over source results. Later (Phase 8) move to
Parquet + vectorized fingerprint columns.

### 3.4 InChIKey skeleton matching

New method `_inchikey_skeleton_search(inchikey)`:

- Extract the first 14 characters of the query InChIKey (the connectivity layer).
- Search each source for InChIKeys starting with the same 14-char prefix.
- This provides a "fuzzy structural match" — same skeleton, possibly different
  stereochemistry or charge.
- Matched results get `match_method="inchikey_skeleton"` and appropriate confidence.

---

## Phase 4 — Salt and solvent stripping

### 4.1 Parent molecule resolution

When `strip_salts=True`:

1. Parse SMILES with RDKit.
2. Use `rdkit.Chem.SaltRemover.SaltRemover` (with default salt definitions, plus
   optional user-supplied SMARTS) to strip salts.
3. Also handle multi-component SMILES (`.`-separated) by picking the largest fragment
   via `Chem.rdmolops.GetMolFrags(mol, asMols=True)` sorted by heavy atom count.
4. Generate `parent_smiles` and `parent_inchikey` from the stripped molecule.
5. Optionally re-run the full lookup pipeline on the parent structure to enrich
   results.

### 4.2 Configuration

```python
Search(
    strip_salts=True,
    salt_smarts: list[str] | None = None,  # additional SMARTS to strip
)
```

---

## Phase 5 — Fuzzy name matching

### 5.1 Integration with rapidfuzz

`rapidfuzz` is already a dependency. Extend name matching in `_resolve_name`:

1. First try exact match across all sources (existing behavior).
2. If no exact match and `fuzzy=True`:
   a. Normalize query: lowercase, strip whitespace, remove common prefixes
      (e.g., "d-", "l-", "dl-", "(±)-", "rac-").
   b. Run `rapidfuzz.process.extract` against the synonym lists of each source
      with `scorer=rapidfuzz.fuzz.WRatio` and `score_cutoff=60`.
   c. Rank results by fuzzy score and cross-source agreement.
   d. Set `match_method="fuzzy_name"` and `confidence` based on fuzzy ratio.

### 5.2 Name normalization utility

Add `_normalize_name(name: str) -> str`:

- Lowercase
- Strip leading/trailing whitespace
- Remove stereochemistry prefixes for matching purposes (keep original in output)
- Collapse multiple spaces
- Handle common abbreviations (e.g., "MEK" → "methyl ethyl ketone") — optional,
  via a small built-in alias dict

---

## Phase 6 — New identifier search methods

### 6.1 `_resolve_inchi(inchi)`

1. Validate InChI format (`starts with "InChI="`).
2. Convert to InChIKey via RDKit.
3. Delegate to `_resolve_inchikey(inchikey)`.
4. Also try direct InChI matching in sources that store it (CompTox, ChEBI).

### 6.2 `_resolve_inchikey(inchikey)`

1. Validate InChIKey format (27-char, `XXXXXXXXXXXXXX-XXXXXXXXXX-X`).
2. Query each source:
   - ChEBI: `search_by_inchikey(inchikey)`
   - CompTox: `get_by_inchikey(inchikey)`
   - PubChemID: needs new `get_by_inchikey()` method (add to `pubchem.py`)
   - ChEMBL: lookup via `standard_inchi_key` column
3. Run consensus and return.

### 6.3 `_resolve_dtxsid(dtxsid)`

1. Query CompTox by DTXSID (add `get_by_dtxsid()` to `comptox.py` if missing).
2. Extract SMILES/InChIKey from CompTox result.
3. Cross-reference other sources using InChIKey.

### 6.4 `_resolve_formula(formula)`

1. Query CompTox and PubChem by molecular formula.
2. Return all matches ranked by completeness (number of non-null identifier fields).
3. Low base confidence (0.30) since formula is not unique.

---

## Phase 7 — Traceability and source reporting

### 7.1 `source_details` field

Each result row includes a `source_details` dict showing exactly what each source
contributed:

```python
{
    "ChEBI": {"found": True, "fields": ["name", "SMILES", "InChIKey"]},
    "CompTox": {"found": True, "fields": ["DTXSID", "molecular_formula"]},
    "PubChemID": {"found": False, "fields": []},
    "ZeroPM": {"found": False, "fields": []},
    "ChEMBL": {"found": True, "fields": ["molecular_mass"]},
}
```

### 7.2 Tracking in `_apply_candidate_to_result`

Modify the existing `_apply_candidate_to_result` helper to also record which fields
were populated from each source, building the `source_details` dict incrementally.

---

## Phase 8 — Performance: Parquet and Polars (later stage)

> **NOTE:** This phase is planned for a later stage of development. The initial
> implementation uses the existing pandas + SQLite/CSV/SDF stack. The Parquet/Polars
> migration is designed to be a drop-in replacement for the data layer without
> changing the public API.

### 8.1 Data format migration: CSV/SQLite → Parquet

- Convert the offline databases to Parquet files:
  - `comptox_chemicals.db` → `comptox_chemicals.parquet`
  - `pubchem_id.db` → `pubchem_id.parquet`
  - `zeropm-v0-0-4.sqlite` → `zeropm.parquet`
  - `chembl_36.db` → `chembl.parquet` (or partitioned Parquet for the ~5 GB DB)
  - `chebi.sdf` → `chebi.parquet` (pre-parsed SDF fields)
- Add a one-time migration script in `scripts/convert_to_parquet.py`.
- Keep the download scripts producing the original formats; add a post-download
  conversion step.
- Include pre-computed columns in Parquet files:
  - `canonical_smiles` (RDKit canonical)
  - `inchikey` (derived from SMILES if missing)
  - `inchikey_skeleton` (first 14 chars of InChIKey, for fast prefix matching)
  - `morgan_fp_2048` (Morgan fingerprint as bit vector bytes, for Tanimoto search)

### 8.2 Internal switch from pandas to polars

- Replace `pd.read_sql` / `pd.read_csv` with `pl.scan_parquet` (lazy evaluation).
- Use Polars expressions for filtering, joining, and aggregation.
- Keep the public API returning `pandas.DataFrame` — call `.to_pandas()` at the
  boundary. Users who want Polars can access `.to_polars()` on the result.
- Add `polars>=0.20.0` as a dependency.

### 8.3 Performance improvements from Parquet + Polars

- **Columnar reads:** Only read the columns needed for each query type.
- **Predicate pushdown:** Filter at the Parquet level (e.g., InChIKey prefix scan).
- **Vectorized fingerprint search:** Store Morgan fingerprints as binary columns;
  compute Tanimoto in batch using numpy/polars bit operations instead of per-row
  RDKit calls.
- **Lazy evaluation:** Chain filters before materializing results.

### 8.4 Migration strategy

1. Add Parquet read/write alongside existing SQLite readers (feature flag).
2. Run benchmarks comparing SQLite vs Parquet for each query type.
3. Once validated, make Parquet the default and deprecate SQLite readers.
4. Update download scripts to produce Parquet directly.

---

## Phase 9 — Backward compatibility wrappers

### 9.1 Keep existing functions

The existing public functions in `tools.py` (`ids_from_CAS`, `ids_from_name`,
`ids_from_SMILES`, `casrn_to_compounds`, `iupac_name_to_id`) remain as thin wrappers
that instantiate a `Search` object internally and delegate:

```python
def ids_from_CAS(cas, chebi=None, comptox=None, pubchem=None, zeropm=None, chembl=None):
    s = Search(
        identifier_type="cas",
        chebi=chebi, comptox=comptox, pubchem=pubchem,
        zeropm=zeropm, chembl=chembl,
        show_progress=False,
    )
    result = s.search(cas)
    # Convert back to legacy dict format (drop new columns)
    return result.iloc[0].to_dict() if not result.empty else _empty_cas_result(cas)
```

### 9.2 Deprecation path

No deprecation warnings needed (per dev-principles: no backward compatibility
burden). But keeping the wrappers avoids breaking existing notebooks and examples
until they are updated.

---

## Phase 10 — Add missing source methods

Some identifier lookups require new methods on the source classes:

| Source class | New method needed          | Purpose                         |
|-------------|----------------------------|---------------------------------|
| `PubChemID` | `get_by_inchikey(key)`     | Lookup by InChIKey in SQLite    |
| `PubChemID` | `get_by_inchi(inchi)`      | Lookup by InChI in SQLite       |
| `PubChemID` | `search_by_formula(f)`     | Lookup by molecular formula     |
| `CompToxID` | `get_by_dtxsid(dtxsid)`   | Lookup by DTXSID (if missing)   |
| `CompToxID` | `search_by_formula(f)`     | Lookup by molecular formula     |
| `CheMBL`    | `search_by_inchikey(key)`  | Lookup by standard_inchi_key    |
| `ChebiSDF`  | `search_by_formula(f)`     | Filter parsed SDF by formula    |

Check which methods already exist before implementing. Only add what is missing.

---

## Phase 11 — Tests

### 11.1 Test file: `tests/test_search.py`

Follow the existing mock pattern from `tests/test_tools.py` (mock source classes
that return canned data without needing real databases).

#### Unit tests for `Search` class:

- **Initialization:**
  - `test_search_init_default` — default params
  - `test_search_init_custom_clients` — injected clients
  - `test_search_init_lazy_clients` — clients created on first search call

- **Dispatch by identifier type:**
  - `test_search_cas` — routes to CAS resolver
  - `test_search_name` — routes to name resolver
  - `test_search_smiles` — routes to SMILES resolver
  - `test_search_inchi` — routes to InChI resolver
  - `test_search_inchikey` — routes to InChIKey resolver
  - `test_search_dtxsid` — routes to DTXSID resolver
  - `test_search_formula` — routes to formula resolver
  - `test_search_invalid_identifier_type` — raises ValueError

- **Input formats:**
  - `test_search_single_string` — single string input
  - `test_search_list_of_strings` — list input
  - `test_search_dataframe_input` — DataFrame input with column
  - `test_search_csv_file_input` — file path input (use tmp_path fixture)

- **Output schema:**
  - `test_output_has_all_columns` — all expected columns present
  - `test_inchikey_always_populated` — InChIKey derived from SMILES when source omits it
  - `test_kekulized_smiles_populated` — kekulized SMILES present when SMILES available

#### Unit tests for confidence scoring:

- `test_confidence_exact_inchikey` — score = 1.0
- `test_confidence_exact_cas` — score ≈ 0.90
- `test_confidence_fuzzy_name` — score proportional to fuzzy ratio
- `test_confidence_tanimoto` — score proportional to similarity
- `test_confidence_modulated_by_consensus` — consensus multiplier applied

#### Unit tests for structure-aware matching:

- `test_normalize_structure` — canonical + kekulized + InChI + InChIKey
- `test_normalize_structure_invalid_smiles` — returns None fields gracefully
- `test_tanimoto_similarity_above_threshold` — returns match
- `test_tanimoto_similarity_below_threshold` — no match
- `test_inchikey_skeleton_match` — 14-char prefix match found
- `test_inchikey_skeleton_no_match` — prefix not in any source

#### Unit tests for salt stripping:

- `test_strip_salts_nacl` — e.g., `[Na+].[Cl-].CC(=O)O` → acetic acid
- `test_strip_salts_largest_fragment` — picks largest fragment
- `test_parent_smiles_populated` — `parent_smiles` and `parent_inchikey` filled
- `test_strip_salts_off_by_default` — no salt stripping when `strip_salts=False`

#### Unit tests for fuzzy name matching:

- `test_fuzzy_name_typo` — "asprin" → "aspirin"
- `test_fuzzy_name_case_insensitive` — "ASPIRIN" matches "aspirin"
- `test_fuzzy_name_off_by_default` — exact match only when `fuzzy=False`
- `test_fuzzy_name_score_cutoff` — low-similarity names rejected

#### Unit tests for traceability:

- `test_source_details_populated` — `source_details` dict present and correct
- `test_source_details_no_match` — all sources show `found=False`

### 11.2 Test file: `tests/test_search_new_methods.py`

Tests for new methods added to source classes (§Phase 10):

- `test_pubchem_get_by_inchikey`
- `test_pubchem_get_by_inchi`
- `test_pubchem_search_by_formula`
- `test_comptox_get_by_dtxsid`
- `test_comptox_search_by_formula`
- `test_chembl_search_by_inchikey`
- `test_chebi_search_by_formula`

### 11.3 Existing tests

`tests/test_tools.py` remains unchanged. The backward-compatible wrappers ensure
existing tests pass without modification.

### 11.4 Running tests

```bash
uv run pytest tests/test_search.py -v
uv run pytest tests/test_search_new_methods.py -v
uv run pytest tests/test_tools.py -v  # regression
```

---

## Phase 12 — Docstrings and documentation

### 12.1 Module-level docstring

`search.py` gets a module docstring explaining the purpose, quick usage example, and
listing the supported identifier types.

### 12.2 Class and method docstrings

All public methods on `Search` follow Google-style docstrings (per dev-principles §6)
with: summary, extended description, Args, Returns, Raises, and Example sections.

### 12.3 Private method docstrings

All `_resolve_*`, `_normalize_*`, `_similarity_search`, and helper methods get
docstrings explaining their logic when non-obvious.

### 12.4 Documentation pages

- Add `docs/api/search.md` with mkdocstrings autodoc for the `Search` class.
- Update `docs/quickstart.md` with a section showing `Search` usage.
- Add example scripts:
  - `examples/search/search_by_cas_demo.py`
  - `examples/search/search_by_name_demo.py`
  - `examples/search/fuzzy_name_demo.py`
  - `examples/search/similarity_search_demo.py`
  - `examples/search/salt_stripping_demo.py`

---

## Phase 13 — Package integration

### 13.1 Exports

Add to `src/provesid/__init__.py`:

```python
from .search import Search
```

### 13.2 Dependencies

Verify these are already in `pyproject.toml` (all should be):

- `rdkit` — structure normalization, fingerprints, salt removal
- `rapidfuzz` — fuzzy name matching
- `pandas` — DataFrame I/O
- `tqdm` — progress bars

No new dependencies needed for the initial implementation.

For Phase 8 (later), add:

- `polars>=0.20.0`
- `pyarrow>=10.0.0` (for Parquet read/write)

---

## Implementation order

| Step | Phase | Description                                   | Depends on |
|------|-------|-----------------------------------------------|-----------|
| 1    | 10    | Add missing source methods                    | —         |
| 2    | 1     | `Search` class skeleton + client init         | —         |
| 3    | 3.1   | `_normalize_structure` utility                | 2         |
| 4    | 2     | Confidence scoring model                      | 2         |
| 5    | 1.3   | CAS / name / SMILES resolvers (wrap existing) | 2, 3      |
| 6    | 6     | InChI / InChIKey / DTXSID / formula resolvers | 1, 2, 3   |
| 7    | 4     | Salt and solvent stripping                    | 3         |
| 8    | 5     | Fuzzy name matching                           | 5         |
| 9    | 3.3   | Tanimoto similarity search                    | 3         |
| 10   | 3.4   | InChIKey skeleton matching                    | 6         |
| 11   | 7     | Traceability / source_details                 | 5, 6      |
| 12   | 9     | Backward compatibility wrappers               | 5         |
| 13   | 11    | Tests                                         | all above |
| 14   | 12    | Docstrings and documentation                  | all above |
| 15   | 13    | Package integration + exports                 | all above |
| 16   | 8     | Parquet + Polars migration (later stage)      | all above |
