# Plan: package-wide refactor — code quality, docs, examples, Search helpers

**Date:** 2026-08-02
**Author:** Ali A Eftekhari + Claude Code
**Scope:** all of `src/provesid/`, `docs/`, `examples/`, `README.md`
**Status:** Proposed. Three items already landed — see §16.

---

## 1. Method

Read every module in `src/provesid/` (36k lines total), ran the web-API test files
against the live services, ran `Search` against the local databases, and measured
docstring coverage with an AST pass. Also read the `Search`-based enrichment
scripts in the private `EXProves` repo (`src/exproves/data/harmonize/`) to extract
the recurring usage patterns for workstream G.

Live service check (2026-08-02):

| Service | Status |
|---|---|
| PubChem PUG REST | 200, working |
| PubChem PUG View | 200, working |
| NCI CACTUS resolver | 200, working |
| ChEBI backend API | 200, working |
| OPSIN web | 301 → `https://www.ebi.ac.uk/opsin/ws/` (works only because `requests` follows redirects) |
| CAS Common Chemistry | 403 without a key; all 23 tests skip |
| **ClassyFire** | **defunct — see below** |

`ClassyFireAPI.submit_query` returns HTTP 201, but the record it returns is stamped
`created_at: 2023-02-10` and the job stays `"In Queue"` indefinitely. The static
entity endpoint (`/entities/InChIKey=BSYNRYMUTXBXSQ-UHFFFAOYSA-N.json`) returns 404.
The service has not classified anything since February 2023. The existing 17 tests
pass because they assert on HTTP status codes and never on a completed
classification.

---

## 2. Current state

### 2.1 Docstring coverage

415 public classes/functions/methods: **32 have no docstring at all**, and **305
have a docstring with no usage example**. Only 78 (19%) carry an example.

| Module | Public | No docstring | No example |
|---|---:|---:|---:|
| pubchem.py | 83 | 9 | 50 |
| zeropm.py | 63 | 0 | 49 |
| chebi.py | 42 | 0 | 27 |
| cache.py | 37 | 2 | 35 |
| pubchemview.py | 30 | 0 | 30 |
| resolver.py | 25 | 0 | 24 |
| comptox.py | 24 | 0 | 24 |
| reach.py | 20 | 2 | 18 |
| chembl.py | 18 | 0 | 2 |
| opsin.py | 16 | **10** | 6 |
| config.py | 14 | 0 | 13 |
| taxonomy.py | 11 | 1 | 8 |
| search.py | 10 | 5 | 1 |
| cascommonchem.py | 7 | 1 | 6 |
| classyfire.py | 6 | 2 | 3 |
| tools.py | 6 | 0 | 6 |
| utils.py | 3 | 0 | 3 |

### 2.2 Per-module assessment

| Module | Lines | Grade | Main problems |
|---|---:|---|---|
| `search.py` | 2561 | good | ~350 lines of near-identical 5-source `try/except` ladders repeated across 10 resolvers; 68 `except Exception`; 2 dead methods (~120 lines); 6 unused imports; fuzzy match mislabels its method |
| `tools.py` | 1048 | poor | 5 legacy resolvers (~470 lines) that duplicate `Search` and are not exported; 23 `except Exception` |
| `zeropm.py` | 2918 | fair | 63 public methods, 0 examples; logs to the **root logger** via `logging.warning`; emits misleading "not found in database" warnings during successful `Search` runs; no module docstring |
| `pubchem.py` | 2458 | poor | online `PubChemAPI` (~1500 lines) and offline `PubChemID` (~950) in one file; 9 missing docstrings; 11 `except Exception`; own copy of the transport layer |
| `chebi.py` | 1469 | fair | online `ChEBI` + offline `ChebiSDF` in one file; own copy of the transport layer |
| `pubchemview.py` | 1088 | poor | `_extract_experimental_value_and_unit` is a **390-line, 111-branch** regex cascade (lines 637–1027); own copy of the transport layer |
| `chembl.py` | 1091 | good | — |
| `comptox.py` | 605 | good | 24 public methods, 0 examples |
| `resolver.py` | 587 | good | own copy of the transport layer |
| `reach.py` | 560 | fair | 106 lines of hand-rolled XLSX parsing (`_read_xlsx_with_stdlib`) plus a 3-branch fallback, solely to avoid depending on `openpyxl` — which is not declared anywhere |
| `cache.py` | 516 | fair | 14 near-identical `clear_<svc>_cache` / `get_<svc>_cache_info` functions; a cached `None` is indistinguishable from a cache miss, so functions returning `None` never cache |
| `taxonomy.py` | 496 | good | — |
| `cascommonchem.py` | 252 | fair | no rate limiting or retry at all |
| `opsin.py` | 188 | poor | 10 of 16 public members undocumented; stale base URL; `if reqdata.status_code != list(self.responses.keys())[0]` instead of `!= 200`; a `message` key that is always `""` but is read back in a log line |
| `config.py` | 155 | fair | `print()` with emoji in library code (4 call sites) |
| `utils.py` | 67 | good | `check_CASRN` lacks an example |
| `classyfire.py` | 116 | dead | see §1 |

### 2.3 Duplicated transport layer

`_rate_limit` + `_make_request` are reimplemented, differently, in four modules
(`pubchem.py:235`, `pubchemview.py:124`, `resolver.py:107`, and `chebi.py:105` as
`_get`/`_get_raw`/`_post_json`/`_post_text`), and are **absent** from
`cascommonchem.py`, `opsin.py`, and `classyfire.py`. **No module handles HTTP 429.**
Only `pubchemview.py` retries at all. This is the direct opposite of what
`.claude/skills/dev-principles` §8 asks for.

### 2.4 Documentation

- `docs/api/` is 3145 lines of which only 8 lines are `:::` mkdocstrings
  directives. `chebi.md`, `classyfire.md`, `opsin.md`, `cascommonchem.md` (692
  lines) and `index.md` contain **zero** directives — they are hand-written prose
  that has drifted.
- `docs/api/index.md` documents `PubChemPUGViewAPI`, a class that does not exist
  (it is `PubChemView`).
- No API page at all for `comptox`, `zeropm`, `reach`, `taxonomy`, `tools`,
  `config`, `cache`, `utils`.
- `docs/plans/` (6 files, incl. a CSV) is published as a "Modernization" nav
  section on the public site.
- `site/` — 85 built HTML/JS artifacts — is committed to git.
- README references three `.ipynb` tutorials that no longer exist, and shows
  `zpm.query_similar_name("formaldehyde", threshold=80)` and
  `zpm.query_by_inventory(inventory_name="REACH")`; the real signatures are
  `query_similar_name(name, number_of_results=5, score_cutoff=80)` and
  `query_by_inventory(source_name=None, source_id=None)`. Both examples raise.

### 2.5 Examples

No notebooks anywhere. What exists instead: 11 MyST-markdown tutorials converted
away from notebooks, 8 loose `.py` demo scripts, one marimo notebook
(`examples/notebooks/notebooks.py`) driving the legacy `tools.ids_from_*`
functions, 4 `demo_*.py` / `debug_props.py` files living inside `tests/`, and two
`__pycache__` directories committed under `examples/chebifier/`. `search/` has
four one-off `.py` demos of 34–45 lines each. `comptox` and `reach` have no
examples at all.

---

## 3. Decisions taken

Confirmed with the user before writing this plan:

1. **ClassyFire** — keep the module, but make it **raise on use** with a clear
   message pointing at `ChebifierClassifier`.
2. **`tools.py`** — delete the 5 legacy resolvers; keep `tools.py` as the
   documented home for the candidate/consensus helpers that `Search` uses.
3. **Examples** — real `.ipynb` with **saved outputs**, committed, rendered by
   `mkdocs-jupyter` with `execute: false`.
4. **Module split** — split both: `pubchem.py` → `pubchem.py` + `pubchem_id.py`,
   `chebi.py` → `chebi.py` + `chebi_sdf.py`. Public imports from `provesid` stay
   identical.

Standing constraints from `.claude/skills/dev-principles`: no backward-compat
shims; offline first; do not alter raw API call logic or response shapes — extend
with wrappers; readability over micro-optimisation; tests run locally only.

---

## 4. Workstream A — one shared HTTP layer

**New file: `src/provesid/http.py`** (~150 lines, one class + one exception base).

```python
class ServiceError(Exception): ...
class NotFoundError(ServiceError): ...
class RateLimitError(ServiceError): ...

class HTTPClient:
    """Rate-limited HTTP client with retry and back-off, shared by every
    PROVESID web-API module."""
    def __init__(self, base_url, *, min_interval=0.2, timeout=30,
                 max_retries=3, backoff=0.5, headers=None,
                 error_cls=ServiceError): ...
    def get(self, path, *, params=None) -> requests.Response: ...
    def get_json(self, path, *, params=None) -> Any: ...
    def post_json(self, path, *, json=None, params=None) -> Any: ...
```

Behaviour, in one place:

- Sleep to honour `min_interval` between requests.
- Retry on 429 (respecting `Retry-After`), 5xx, timeouts and connection errors,
  with exponential back-off up to `max_retries`.
- Map 404 → `NotFoundError`, exhausted retries → `error_cls`. No raw
  `requests` exception ever escapes.
- Log retries at `WARNING` and each request at `DEBUG`, on a module logger. No
  `print`.

**Migration, module by module.** In each case the URL construction and the
response parsing — the parts that encode the upstream contract — are moved
verbatim; only the transport is swapped.

| Module | Change |
|---|---|
| `pubchem.py` | delete `_rate_limit` + `_make_request` (L235–300); `PubChemError`/`PubChemNotFoundError`/`PubChemTimeoutError`/`PubChemServerError` become subclasses of the shared bases so existing `except` clauses keep working |
| `pubchemview.py` | delete `_rate_limit` + `_make_request` (L124–168) |
| `resolver.py` | delete `_rate_limit` + `_make_request` (L107–149) |
| `chebi.py` | `_get`/`_get_raw`/`_post_json`/`_post_text` (L105–199) become 1–3 line calls |
| `cascommonchem.py` | **gains** rate limiting + retry; the per-method `try/except requests...` blocks (~40 lines across 2 methods) collapse |
| `opsin.py` | **gains** rate limiting + retry; fix base URL to `https://www.ebi.ac.uk/opsin/ws/`; replace `list(self.responses.keys())[0]` with `200`; drop the always-empty `message` key |

Net: roughly −250 lines, and 429 handling everywhere instead of nowhere.

### A.2 Separate raw API calls from convenience wrappers

Per dev-principle §3 the raw calls are untouched; what changes is that each web
module gets two clearly labelled sections and the convenience layer is documented
as such.

- **`pubchem.py`** — raw endpoint methods (`get_compound_by_cid`,
  `get_cids_by_name`, `get_compound_properties`, the four structure searches, the
  substance and assay methods) in a `# --- Raw PUG-REST endpoints ---` section;
  the derived layer (`get_basic_compound_info`, `get_all_compound_info`,
  `search_compound`, `format_search_compound_result`, `_format_single_compound`,
  `extract_identifiers_from_synonyms`, `get_compound_identifiers`,
  `find_cids_comprehensive`) under `# --- Convenience layer ---`, each docstring
  naming the raw method it builds on. Also collapse the four
  `_cached_*_search` / `*_search` pairs (L816–931) — the split exists only to make
  `**options` hashable, which the shared `@cached` key builder already handles.
- **`pubchemview.py`** — move the parsing out of the client entirely. New file
  `src/provesid/pubchemview_parse.py` holding pure, testable functions with no
  HTTP: `parse_value(value_str, property_name=None) -> ParsedValue`, plus
  `extract_temperature`, `extract_conditions`. The 390-line `if/elif` tree becomes
  a module-level table

  ```python
  PROPERTY_PATTERNS: dict[str, list[Pattern]] = {
      "vapor pressure": [...],
      "logp": [...],
      ...
  }
  ```

  read by one ~40-line matcher. `PubChemView` keeps a thin
  `parse_value` re-export so `tests/test_pubchemview.py:201,356` and the tutorial
  keep working against a **public** name instead of a private method. Target:
  390 lines → ~140 plus a data table.
- **`classyfire.py`** — every method raises `ServiceUnavailableError` with the
  Feb-2023 evidence and a pointer to `ChebifierClassifier`. The URL-building code
  stays in place, commented, so the module can be revived if the service returns.
  Tests rewritten to assert the raise.

**Size:** L. **Verification:** the web-API test files still pass; add tests for
`http.py` (429 + `Retry-After`, 5xx retry, 404 mapping) against a local stub.

---

## 5. Workstream B — module cleanup

Worst-first. Each item is independently committable.

### B.1 `tools.py` — delete the legacy resolvers (S)

Delete `ids_from_CAS`, `ids_from_name`, `ids_from_SMILES`, `casrn_to_compounds`,
`iupac_name_to_id`, `smiles_to_canonical`, `_best_candidate_by_name` — ~570 lines,
none exported from `__init__.py`, all superseded by `Search`. `tools.py` keeps the
candidate/consensus helpers, renamed without the leading underscore and given
proper docstrings, since `search.py` imports them across module boundaries:

`is_missing`, `pick_first`, `normalize_synonyms`, `to_float`, `text_similarity`,
`extract_cas_values`, `inchi_to_smiles`, `inchikey_from_smiles`, `first_cas`,
`make_candidate`, `candidate_similarity`, `candidate_compatible_with_consensus`,
`apply_candidate_to_result`, `compute_consensus`, `smiles_to_canonical_and_mass`,
and the six `candidate_from_*_row` adapters.

Also delete `tests/test_tools.py`'s legacy-resolver tests and rewrite
`examples/notebooks/notebooks.py` as a `Search` notebook (workstream E).

### B.2 `search.py` — collapse the 5-source ladders (M)

Ten resolvers each contain the same shape: `if self._chebi is not None: try: rows
= ...; except Exception as exc: log.warning(...)`, five times over. Replace the
five `self._chebi` … `self._chembl` attributes with one `self._clients: dict[str,
Any]`, and move the per-source call sites into a plain lookup table in a new
`src/provesid/sources.py`:

```python
# sources.py — no classes, just small functions in a table.
LOOKUPS: dict[str, dict[str, Callable[[Any, str], list[dict]]]] = {
    "cas": {
        "chebi":   lambda c, q: c.search_by_cas(q),
        "comptox": lambda c, q: one(c.get_by_casrn(q)),
        "pubchem": lambda c, q: one(c.get_by_cas(q)),
        "zeropm":  lambda c, q: rows(c.get_id_table_from_cas(q)),
    },
    "inchikey": {...}, "inchi": {...}, "smiles": {...},
    "name": {...}, "formula": {...},
}
TO_CANDIDATE = {"chebi": candidate_from_chebi_row, ...}
```

and one driver on `Search`:

```python
def _collect(self, id_type, query, method, k=1, score=1.0):
    """Query every available source for `query`; return a tagged candidate pool."""
```

`_resolve_cas` drops from 65 lines to ~8; the same for `_resolve_inchikey`,
`_resolve_inchi`, `_resolve_dtxsid`, `_resolve_smiles`, `_resolve_formula`,
`_inchikey_pool`, and both passes of `_candidate_pool_from_name`. Net ≈ −250
lines, and adding a sixth source becomes a one-file edit.

Also in `search.py`:

- Delete `_fuzzy_name_candidates` (L1745–1805) and `_candidates_from_name`
  (L1302–1360) — dead since the multi-hit refactor, ~120 lines. Delete
  `_most_complete_row` (L2461) and the 4 tests that exercise it.
- Delete 6 unused imports from `tools` (`normalize_synonyms`,
  `extract_cas_values`, `first_cas`, `inchi_to_smiles`, `to_float`,
  `smiles_to_canonical_and_mass`).
- Give the 5 undocumented nested helpers (`find`, `union`, `add`, `add_fuzzy`) a
  one-line docstring each, or inline them.
- **Correctness:** `Search("name", fuzzy=True).search("asprin")` currently returns
  PHENYRAMIDOL labelled `match_method="exact_name"` with confidence 0.74. A typo
  fed to a fuzzy search should either resolve to aspirin or return no hit — and
  must never be labelled `exact_name`. Fix the labelling and add a case to
  `tests/test_search_precision_regression.py`. This is behaviour-affecting; flag
  it in `CHANGELOG.md`.

### B.3 `pubchem.py` → `pubchem.py` + `pubchem_id.py` (M)

Move `PubChemID` (L1508–2458, ~950 lines) into `pubchem_id.py` verbatim; keep the
`Domain` / `CompoundProperties` / namespace constant classes with the online
client. Give the 9 undocumented constant classes a docstring each.
`from provesid import PubChemAPI, PubChemID, ...` is unchanged.

In `pubchem_id.py`, the 20 one-line `cas_to_*` / `*_to_*` converters (L1958–2023)
and their 6 `batch_*` twins are fine as they are — they are the readable public
surface — but they need examples.

### B.4 `chebi.py` → `chebi.py` + `chebi_sdf.py` (S)

Move `ChebiSDF` (L850–1469) and the two module functions `get_chebi_entity` /
`search_chebi` stay with the online client. No logic changes.

### B.5 `zeropm.py` (M)

- Add a module docstring and a class docstring that says what the database is and
  where it comes from.
- Replace all 24 `logging.<level>` root-logger calls with `self.logger` (already
  created in `__init__`).
- Downgrade the five `"... not found in database"` messages from `WARNING` to
  `DEBUG`. They currently fire during perfectly successful `Search` runs — a
  `Search("name").search(["aspirin"])` prints `WARNING:root:Chemical name
  'aspirin' not found in database` and then returns aspirin.
- 63 public methods, 0 examples: covered by workstream C.
- Do **not** split this file in this pass; it is one coherent SQLite interface and
  splitting it has no payoff. Revisit if it keeps growing.

### B.6 `reach.py` (S)

Add `openpyxl` to `[project.dependencies]`, delete `_read_xlsx_with_stdlib` (106
lines) and reduce `_load_dataframe` to a single `pd.read_excel` call with one
error wrap. Give `column_index` and `read_cell` docstrings or delete them with
their parent.

### B.7 `cache.py` (S)

- Replace the 14 hand-written `clear_<service>_cache` / `get_<service>_cache_info`
  functions with two parameterised ones, `clear_cache(service=None)` and
  `get_cache_info(service=None)`, and update the seven call sites in
  `classyfire.py`, `opsin.py`, `cascommonchem.py`, `pubchem.py`,
  `pubchemview.py`, `resolver.py`. Net ≈ −80 lines and a smaller public surface.
- Fix the sentinel bug: `_load_from_disk` returns `None` both for "absent" and for
  "cached value is `None`", so `get()` reports a miss and the call is repeated
  forever. Use a `_MISS` sentinel.
- `@cached` currently caches the result even when `use_cache=False` (L388–390).
  Decide and document: `use_cache=False` should mean "don't read", which is what
  it does — say so in the docstring, or make it mean "bypass entirely".
- Document `decorator` and `wrapper`.

### B.8 `config.py` (S)

Replace the 4 `print()` calls (`"✅ CAS API key saved to: ..."` etc.) with
`logger.info`. A library should not print, and should not print emoji. `show_config()`
is a deliberate user-facing reporter — have it **return** the info dict and let
the caller print, or keep the print and document it as interactive-only.

### B.9 `opsin.py` (S)

10 of 16 public members have no docstring at all, including the `OPSIN` class
itself. Document all of them; the six `PYOPSIN.get_*` one-liners are fine as code
and only need docstrings with an example. Fold `_empty_res` into a module-level
constant.

**Workstream B total:** roughly −1400 lines with no loss of capability.

---

## 6. Workstream C — docstrings with examples everywhere

Target: **every public class, function and method has a Google-style docstring
with a runnable example.** 32 to write from scratch, 305 to extend with an
example. This is the largest mechanical task in the plan.

House style, one block, no exceptions:

```python
def get_by_cas(self, cas: str) -> Optional[Dict[str, Any]]:
    """Look up one compound by CAS Registry Number.

    Args:
        cas: CAS Registry Number, with or without hyphens.

    Returns:
        Row dict with ``cid``, ``smiles``, ``inchi``, ``inchikey``,
        ``mf``, ``mw`` and ``cmpdname``, or ``None`` when the CAS is
        not in the database.

    Example:
        >>> db = PubChemID()
        >>> db.get_by_cas("50-78-2")["inchikey"]
        'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'
    """
```

Rules:

- `>>>` doctest form, because it renders well in mkdocstrings and can be checked.
- Offline-database examples use real values from the shipped databases and give
  real output. Online examples show a plausible truncated result.
- Examples that need network or a 30 GB database are **not** collected by default;
  add `--doctest-modules` to a separate opt-in pytest marker so `pytest -m doctest`
  can verify the offline ones locally.

Order of work, largest deficit first: `pubchem.py`/`pubchem_id.py` (50) →
`zeropm.py` (49) → `cache.py` (35) → `pubchemview.py` (30) → `chebi.py`/`chebi_sdf.py`
(27) → `resolver.py` (24) → `comptox.py` (24) → `reach.py` (18) → `config.py` (13) →
`taxonomy.py` (8) → the rest.

**Size:** L. This is the prerequisite for workstream D — the docs are only as good
as these docstrings.

---

## 7. Workstream D — rebuild the documentation from docstrings

Delete `docs/` and rebuild. The rule: **anything that can rot lives in a
docstring, not in a `.md` file.**

### D.1 Delete

`docs/api/*.md` (all 10, 3145 lines), `docs/quickstart.md`, `docs/data_methods.md`,
`docs/advanced_caching.md`, `docs/plans/` (6 files — move the two real plans into
`plans/`, delete the rest), and `site/` (85 committed build artifacts; add `site/`
to `.gitignore` and let the deploy workflow build it).

### D.2 New structure

```
docs/
  index.md              hand-written, ~60 lines: what it is, install, where to go
  api-keys.md           CAS key setup (moved from root API_KEY_GUIDE.md)
  install-chebifier.md  the chebifier/torch install steps (from docs/chebifier.md)
  api/
    index.md            a table: module -> one-line purpose -> link
    search.md           ::: provesid.search
    tools.md            ::: provesid.tools
    sources.md          ::: provesid.sources
    pubchem.md          ::: provesid.pubchem
    pubchem_id.md       ::: provesid.pubchem_id
    pubchemview.md      ::: provesid.pubchemview
    pubchemview_parse.md
    chebi.md            ::: provesid.chebi
    chebi_sdf.md        ::: provesid.chebi_sdf
    chembl.md  comptox.md  zeropm.md  reach.md
    resolver.md  opsin.md  cascommonchem.md  classyfire.md
    taxonomy.md  http.md  cache.md  config.md  utils.md
```

Every `api/*.md` is exactly a `# Title` line plus one `:::` directive. Nothing
hand-written, so nothing can drift. Tutorials in the nav are the notebooks from
workstream E, rendered by `mkdocs-jupyter`.

### D.3 `mkdocs.yml`

- Nav: Home / Getting started (api-keys, install-chebifier) / Tutorials (the
  notebooks) / API reference (the generated pages). Drop the "Modernization" section.
- `mkdocs-jupyter`: `execute: false`, point at `examples/**/*.ipynb`. Drop the
  `jupytext` docs dependency and `scripts/convert_notebooks_to_myst.sh` /
  `scripts/generate_notebooks_from_myst.sh` / `scripts/validate_docs_local.sh`.
- `mkdocstrings`: `docstring_style: google`, `show_source: true`,
  `members_order: source`, `show_if_no_docstring: false`, and
  `show_signature_annotations: true`.
- Fix `extra.social` — it points at `github.com/provesid/provesid`, which is not
  the repo.

**Verification:** `mkdocs build --strict` must pass. Strict mode fails the build on
any broken reference, which is what caught `PubChemPUGViewAPI` in the first place.

**Size:** M (small once workstream C is done; it is mostly deletion).

---

## 8. Workstream E — examples as Jupyter notebooks

Delete every MyST `.md` tutorial, every loose `.py` demo, `examples/notebooks/`
(marimo), the committed `__pycache__` directories, and the four `demo_*.py` /
`debug_props.py` files sitting in `tests/`. Replace with 17 notebooks, real
`.ipynb`, outputs saved, one per feature area:

| Notebook | Covers |
|---|---|
| `search/01_search_basics.ipynb` | all 7 identifier types, output schema, confidence |
| `search/02_tuning_and_multi_hit.ipynb` | `n_hits`, `fuzzy`, `use_opsin`, `strip_salts`, `similarity_threshold`, `min_confidence`, `return_alternatives` |
| `search/03_enriching_datasets.ipynb` | workstream G — the three EXProves recipes |
| `pubchem/pubchem_api.ipynb` | `PubChemAPI`: raw endpoints then the convenience layer |
| `pubchem/pubchem_id_offline.ipynb` | `PubChemID` SQLite lookups and batch converters |
| `pubchemview/experimental_properties.ipynb` | `PubChemView`, `get_property_table`, value parsing |
| `chebi/chebi_api.ipynb` | `ChEBI` web API, ontology walks, structure search |
| `chebi/chebi_sdf_offline.ipynb` | `ChebiSDF` index build and lookups |
| `chembl/chembl.ipynb` | `CheMBL` |
| `comptox/comptox.ipynb` | `CompToxID` — **new, no example exists today** |
| `zeropm/zeropm.ipynb` | `ZeroPM` inventories, countries, PM probabilities |
| `reach/reach.ipynb` | `REACHDossierID` — **new, no example exists today** |
| `resolver/nci_resolver.ipynb` | `NCIChemicalIdentifierResolver` + the `nci_*` functions |
| `opsin/opsin.ipynb` | `OPSIN` (web) and `PYOPSIN` (local) |
| `cascommonchem/cas_common_chemistry.ipynb` | `CASCommonChem` incl. key setup |
| `taxonomy/chebifier.ipynb` | `ChebifierClassifier`, replacing the two `.py` demos |
| `cache_and_config/caching.ipynb` | `@cached`, cache sizing, export/import, config |

Conventions: first cell is a markdown intro stating what the reader will learn and
what it costs (network? which database? how large?); imports in one cell; every
code cell has a short markdown lead-in; outputs saved so GitHub renders the
notebook without execution; no cell takes longer than ~30 s where avoidable.

`examples/README.md` gets a table of the 17 notebooks with a one-line description
and the prerequisites of each.

**Size:** L.

---

## 9. Workstream F — README

Rewrite to one minimal, correct example per public component, in the existing
voice: prose paragraphs, no emoji, no "✨ Enhanced" banners, no feature-matrix
tables.

Structure:

1. What PROVESID is (keep the current PROVES-family paragraph, it is good).
2. Installation — keep the `uv` recommendation and the disk-space warning.
3. **`Search` first**, because it is the entry point most users want:
   ```python
   from provesid import Search
   df = Search("cas").search(["50-00-0", "64-17-5", "50-78-2"])
   df[["CASRN", "name", "canonical_smiles", "InChIKey", "confidence"]]
   ```
4. Then 3–6 lines per component, verified against the code: `PubChemAPI`,
   `PubChemView`, `PubChemID`, `NCIChemicalIdentifierResolver`, `OPSIN`/`PYOPSIN`,
   `CASCommonChem`, `ChEBI`/`ChebiSDF`, `CheMBL`, `CompToxID`, `ZeroPM`,
   `REACHDossierID`, `ChebifierClassifier`. Nothing more — the notebooks carry
   the depth.
5. Offline data directory and the `PROVESID_DATA_DIR` / `data_dir` / `redownload`
   controls (this section is already accurate; keep it).
6. Note that `ClassyFireAPI` is retained but non-functional, and point at
   `ChebifierClassifier`.
7. Related tools and TODO — keep as-is.

Fixes required regardless: the three dead `.ipynb` links, the
`query_similar_name(threshold=...)` example, and the
`query_by_inventory(inventory_name=...)` example. Every snippet in the new README
gets executed once before commit.

**Size:** S.

---

## 10. Workstream G — `Search` helpers for data enrichment

The `EXProves` harmonization scripts (22 files, ~3900 lines) all reimplement the
same three patterns on top of `Search`. Lifting them into PROVESID removes that
duplication and makes `Search` genuinely useful to anyone with a half-annotated
dataset. Three additions, deliberately small — two functions and one method.

Note: `EXProves` stays private. The example notebook ships a small sample built
from public identifiers (CAS numbers, names, SMILES) shaped like the real inputs;
no EXProves data file is copied.

### G.1 `Search.enrich(df, column, prefix="provesid_")` (method)

The `harmonize_chemicals.py` pattern: search only the **unique** values of a
column, then left-merge the results back so repeated identifiers cost nothing.

```python
def enrich(self, df, column, *, prefix="provesid_", n_hits=None):
    """Search the unique values of `column` and merge the results onto `df`.

    Duplicated identifiers are searched once. The output has one row per
    input row, with every Search column added under `prefix`.

    Example:
        >>> df = pd.read_csv("yaws_boiling_points.csv")   # 8000 rows, 3000 unique CAS
        >>> out = Search("cas").enrich(df, "CAS")
        >>> out["provesid_InChIKey"].notna().mean()
        0.83
    """
```

### G.2 `resolve_cascade(...)` (module-level function in `search.py`)

The ONS / eawag / sangster-logP / pKa / VEGA pattern: try identifier types in
order, stop at the first hit that passes validation, fall back to RDKit, record
how each row was resolved.

```python
def resolve_cascade(df, stages, *, accept=None, fallback="rdkit",
                    prefix="provesid_"):
    """Resolve each row through a cascade of Search stages, first success wins.

    Args:
        df: Input frame.
        stages: List of ``(label, Search, column)``. Tried in order; only
            rows still unresolved are passed to the next stage.
        accept: Optional ``(hit_row, input_row) -> bool`` predicate. A hit
            that fails it is discarded and the row stays pending. Use
            `mw_within` for the usual molecular-weight check.
        fallback: ``"rdkit"`` derives SMILES/InChI/InChIKey/mass from the
            row's own structure when no stage validated; ``None`` disables.
        prefix: Prefix for the added columns.

    Returns:
        `df` with the resolved columns added, plus ``<prefix>resolved_by``
        and ``<prefix>validated_by``.

    Example:
        >>> out = resolve_cascade(
        ...     df,
        ...     stages=[("cas",    Search("cas"),                  "CASRN"),
        ...             ("name",   Search("name", use_opsin=True), "name"),
        ...             ("smiles", Search("smiles"),               "SMILES")],
        ...     accept=mw_within(0.5, reference_column="SMILES"),
        ... )
        >>> out["provesid_resolved_by"].value_counts()
        cas       412
        name       98
        smiles     31
        rdkit      14
    """
```

### G.3 `mw_within(tolerance, reference_column)` (validator factory)

The validation step every EXProves script writes by hand: compare the RDKit
molecular weight of the hit's structure against the weight computed from the
structure the dataset already had. Returns the predicate `resolve_cascade` wants,
and records which checks passed (`mw`, `mw+smiles`, `mw+smiles+name`) in
`validated_by`. Ships alongside two smaller ones, `smiles_matches` and
`name_matches`, so validators compose.

### G.4 Notebook `search/03_enriching_datasets.ipynb`

Three worked recipes, taken straight from the real usage:

1. **Unique-value lookup and merge** — a measurement table keyed by CAS, enriched
   with `Search("cas").enrich(df, "CAS")`. Shows the hit rate and what to do with
   the misses.
2. **Cascade for a missing identifier** — an AqSolDB-shaped frame that has
   InChIKey/InChI/name/SMILES but no CAS: cascade
   `inchikey → inchi → name → name+opsin` and report the CAS recovery rate at
   each stage.
3. **Validated cascade with RDKit fallback** — the ONS recipe:
   `cas → name(+opsin) → smiles`, every hit checked with `mw_within(0.5)`,
   RDKit fallback for the remainder, and a final table of
   `resolved_by` × `validated_by` so the reader can see exactly how much of the
   result is trustworthy and why.

Each recipe ends with the one-liner it replaces, so a reader can lift it directly.

**Size:** M. **Verification:** unit tests for `enrich` (dedup correctness, no row
fan-out on duplicate keys, column collision handling), `resolve_cascade` (stage
ordering, pending-set bookkeeping, fallback), and `mw_within` (tolerance edges,
unparseable SMILES).

---

## 11. Repo hygiene

- `git rm -r site/` and add `site/` to `.gitignore`; the deploy workflow builds it.
- Remove `examples/chebifier/__pycache__/` from git; `.gitignore` covers it going
  forward.
- Move `tests/demo_*.py` and `tests/debug_props.py` out of `tests/` — their content
  belongs in notebooks.
- Root `.md` files: `API_KEY_GUIDE.md` → `docs/api-keys.md`; keep `TESTING.md`,
  `CHANGELOG.md`, `LICENSE`, `README.md`.
- `pyproject.toml`: `requires-python = ">=3.12"` but the classifiers advertise
  3.8–3.11 and `[tool.black] target-version = ['py38']`. Fix both to 3.12.
- `pyproject.toml` declares `[build-system] hatchling` **and** a
  `[tool.setuptools]` section — dead configuration; delete the setuptools block.
- Add `openpyxl` (workstream B.6).
- **Needs your call:** `.github/workflows/` contains `test.yml` and
  `test-with-api-keys.yml`, and the README carries a Tests badge — but
  `dev-principles` §5 says tests run locally only. Remove the two test workflows
  and the badge, keeping only `mkdocs-deploy.yml` and `release.yml`? I have not
  planned this either way.

---

## 12. Sequencing

Each step is a commit that leaves the package working and the tests green.

| # | Commit | Workstream | Size |
|---:|---|---|---|
| 1 | repo hygiene: drop `site/`, `__pycache__`, fix `pyproject.toml` metadata | 11 | S |
| 2 | delete legacy resolvers from `tools.py`, rename helpers public | B.1 | S |
| 3 | add `http.py`; migrate `resolver.py` and `pubchemview.py` onto it | A | M |
| 4 | migrate `pubchem.py`, `chebi.py`, `cascommonchem.py`, `opsin.py` onto `http.py` | A | M |
| 5 | `classyfire.py` raises; rewrite its tests | A.2 | S |
| 6 | extract `pubchemview_parse.py`; shrink the 390-line regex tree | A.2 | M |
| 7 | split `pubchem.py` → `+ pubchem_id.py`; split `chebi.py` → `+ chebi_sdf.py` | B.3, B.4 | M |
| 8 | add `sources.py`; collapse the `Search` source ladders; delete dead code | B.2 | M |
| 9 | fix the fuzzy-name mislabelling + regression test | B.2 | S |
| 10 | `cache.py` parameterisation and `_MISS` sentinel | B.7 | S |
| 11 | `zeropm.py` logging, `reach.py` xlsx, `config.py` printing, `opsin.py` docs | B.5–B.9 | M |
| 12 | workstream G: `enrich`, `resolve_cascade`, `mw_within` + tests | G | M |
| 13 | docstrings with examples, module by module (several commits) | C | L |
| 14 | rebuild `docs/` from docstrings; `mkdocs build --strict` clean | D | M |
| 15 | the 17 notebooks (several commits) | E | L |
| 16 | rewrite `README.md`; verify every snippet | F | S |
| 17 | `CHANGELOG.md` entry for the release | — | S |

Steps 1–12 are code and can proceed independently of 13–16. Step 13 must land
before 14, and 12 before the enrichment notebook in 15.

---

## 13. Verification

- `pytest` fully green locally at every commit; no new GitHub Actions test jobs.
- `mkdocs build --strict` clean — this is the gate that keeps the docs honest.
- `pytest -m doctest` for the offline docstring examples.
- All 17 notebooks execute top to bottom on a machine with the databases present,
  before their outputs are committed.
- Every README snippet executed once before commit.
- `Search` behaviour: `tests/test_search_precision_regression.py` extended with
  the `asprin` case and re-run before/after step 8 to prove the ladder collapse is
  behaviour-preserving.

## 14. Explicitly not doing

- Not splitting `zeropm.py` (2918 lines). It is one coherent SQLite interface;
  splitting buys nothing.
- Not touching the confidence/consensus **algorithm** in `Search` beyond the
  mislabelling fix. That is a separate piece of work with its own evaluation set.
- Not adding new data sources.
- Not adding backward-compat shims for the deleted `tools.py` resolvers, per
  dev-principle §1.
- Not reviving ClassyFire. The module stays, raising, in case the service returns.

## 15. Estimated effect

| | Before | After |
|---|---:|---:|
| `src/provesid/` lines | ~16 700 | ~15 300 |
| Public objects with a docstring | 383 / 415 | 415 / 415 |
| Public objects with an example | 78 / 415 | 415 / 415 |
| Copies of the HTTP transport layer | 4 (+3 modules with none) | 1 |
| Modules handling HTTP 429 | 0 | all |
| `docs/` hand-written lines that can rot | ~3 800 | ~120 |
| Notebooks | 0 | 17 |
| Longest function | 390 lines | ~40 |

---

## 16. Landed on 2026-08-02

Three items from this plan are done. The rest of the plan is unchanged.

### 16.1 The `Search` name-resolution bug (§5 B.2, step 9)

`Search("name", fuzzy=True).search("asprin")` returned PHENYRAMIDOL labelled
`exact_name`. Root-caused to **four** defects, all fixed:

1. `CheMBL.search_by_name` had no exact mode — it is
   `LIKE '%name%'` over preferred names *and* synonyms. PHENYRAMIDOL carries the
   synonym `"Evasprin"`, which contains `"asprin"`. `Search`'s exact pass called
   it and tagged the result `exact_name`. Added `exact=False` (default preserves
   the old behaviour for direct callers); `Search`'s exact pass passes `exact=True`.
2. The `strong` test that gates fuzzy widening used a fuzzy score, and
   `WRatio("asprin", "Evasprin") == 85.7` cleared the 80 cut-off — so the widening
   that would have found aspirin never ran. It now requires actual name equality
   via a new `_matches_name_exactly` helper.
3. **Confidence inversion:** a fuzzy match's base was the raw similarity (→1.0)
   while `exact_name` was pinned at 0.80, so a typo scored *higher* than the
   correct spelling (0.9025 vs 0.7806). The fuzzy base is now scaled by the
   exact-name base.
4. The ZeroPM fuzzy branch was **dead code** — it acted only on a `DataFrame`
   while `query_similar_name` returns a list of ids. ZeroPM is the only source
   doing true fuzzy *retrieval*, so `fuzzy=True` could only ever find typos that
   happened to be registered synonyms. Added
   `ZeroPM.match_similar_name` / `get_id_table_from_similar_name` (which report
   what matched and how well) and wired them in.

Enabling ZeroPM retrieval initially made things *worse* — `caffiene` matched a
compound named `"ne"` and `zzzznotachemical` matched `"Mica"` — because `WRatio`'s
partial-ratio term scores a short candidate highly whenever it appears inside the
query. Hence one behaviour change: **the default `fuzzy_scorer` is now `"ratio"`**.

| query | before | after |
|---|---|---|
| `asprin` | PHENYRAMIDOL, `exact_name`, 0.7429 | acetylsalicylic acid, `fuzzy_name`, 0.7199 |
| `aspirin` | acetylsalicylic acid, 0.7806 | unchanged |
| `caffiene` | not found | caffeine, `fuzzy_name` |
| `tolune` | not found | Toluene, `fuzzy_name` |
| `zzzznotachemical` | not found | not found |
| `asprin` (fuzzy off) | not found | not found |

Tests: 8 new cases in `tests/test_search_precision_regression.py` (right compound,
never `exact_name`, typo scores below correct spelling, nonsense matches nothing,
typo needs `fuzzy=True`), 4 in `tests/test_chembl.py`, 1 in `tests/test_search.py`.
Two existing confidence tests were updated — they encoded the inverted semantics.
The stubs in `tests/test_search.py` were corrected to match the real signatures;
they had drifted, so those code paths were silently failing into `except` blocks.

### 16.2 The enrichment helpers (§10, step 12)

`Search.enrich`, `resolve_cascade` and `mw_within` are implemented and exported
from `provesid`, with 41 unit tests in `tests/test_search_enrichment.py` covering
dedup, no row fan-out, index/column preservation, stage bookkeeping, both
`accept` return forms, the RDKit fallback and the error guards.

Still outstanding from §10: the `search/03_enriching_datasets.ipynb` notebook
(depends on workstream E) and `smiles_matches` / `name_matches` as standalone
validators — `mw_within` currently reports SMILES and name agreement itself,
which covers the observed EXProves usage.

### 16.3 CI removal (§11)

Deleted `test.yml`, `test-with-api-keys.yml` and `release.yml`; only
`mkdocs-deploy.yml` remains. Releases are made manually with `twine`. Removed the
Tests badge from `README.md` and the CI sections from `TESTING.md`.

### 16.4 The 10 pre-existing test failures — all resolved

`pytest tests/` is now **653 passed, 32 skipped, 0 failed** (was 10 failed, 635
passed; 4 chebifier tests were removed rather than kept, see below). Verified green twice: once with the ZeroPM `idx_*` indexes absent and once
with them present, because those change SQLite's query plans (see below).

Investigating them turned up **two more pieces of dead production code**, neither of
which had ever worked:

1. **ZeroPM P/M probabilities (5 tests).** `get_pm_probabilities`,
   `batch_get_pm_probabilities` and `get_all_zeropm_chemicals(include_pm_probs=True)`
   all keyed `pm_probabilities` on `zeropm_id`; that table is keyed on `inchi_id`.
   Every call raised `OperationalError`. Fixed by translating through the new
   `zeropm_id_to_inchi_id` and correcting the two joins — the corrected join yields
   97,491 rows, and `zeropm_id` ↔ `inchi_id` is 1:1 in `zeropm_chemicals`.
   The tests repeated the same wrong key in their own SQL, and their
   `else: pytest.skip(...)` branches converted the resulting error into a silent
   skip — which is exactly how the production bug survived. They now assert the
   fixture query found data and compare returned values against the database.
2. **chebifier live classification (2 tests).** `chemlog_extra` reads
   `data/chebi_v244/<Classifier>_element_class_mapping.csv` **relative to the
   working directory** and rebuilds it from the ChEBI graph when missing — but that
   rebuild crashes, because 288 of the graph's 205,592 nodes carry `name: None` and
   the builder evaluates `" molecular entity" in properties["name"]`. So
   `ChebifierClassifier.classify()` failed for everyone with
   `TypeError: argument of type 'NoneType' is not iterable`. (It works in the private
   EXProves repo only because that repo happens to have the two CSVs on disk.) Fixed
   with `ensure_element_class_mappings`, which derives both files using upstream's
   own rules — None-safe — into the PROVESID chebifier data dir, plus a scoped
   `contextlib.chdir` around ensemble construction. This reproduces upstream's files
   exactly (117 and 36 entries). Unrelated to [[chebifier-gnn-index-drift]], which is
   a separate torch/index problem.

   **The two live tests were then removed** at the user's direction (2026-08-02):
   chebifier remains an optional extra, and its model stack — transformer, graph/GNN,
   rule-based and c3p models, each a separate package and some git-only — is awkward
   enough to install that the suite must not depend on it. The production fix stays,
   since it is what makes the extra work at all. `test_taxonomy.py` keeps its 18
   tests, all of which pass with the stack absent (verified by running them with the
   modules made unimportable), and the unused `chebifier` pytest marker was dropped
   from `pyproject.toml`. `default_ensemble_available()` /
   `missing_ensemble_modules()` were kept and are now the supported way to check
   whether a partial install can actually classify — they matter more, not less,
   without live tests. `ensure_element_class_mappings` is consequently untested;
   exercising it requires the extra.

The remaining three were test-only defects:

| Test | Verdict |
|---|---|
| `test_pubchem_id.py::test_init_nonexistent_path` | Omitted `auto_download=False`, so instead of asserting `FileNotFoundError` it **downloaded the ~2.3 GB database into the repo root**. Fixed; now uses `tmp_path` and the file runs in 1.2 s |
| `test_search.py::test_exact_inchikey_with_low_consensus` | The code was right. `_compute_consensus` returns 0.0 **only** when there are no candidates — one source scores 1.0, two fully disagreeing sources score 0.5 — so zero consensus is the no-match case and must not report half confidence. Replaced with tests for the zero and the partial-agreement cases, and documented the short-circuit in `_compute_confidence` |
| `test_zeropm.py::test_get_id_table_from_cas_existing` | Asserted `dtype == object`; pandas 3 gives `StringDtype`. Now uses `pandas.api.types.is_string_dtype`, correct on both |

One further failure surfaced once the others were fixed.
`test_zeropm.py::test_get_cas_from_name_integration` asserted a name → CAS →
same-CAS round trip that the data model does not support: names are many-to-many
with CAS numbers, so `get_names("121-20-0")` includes `"Jasmolin II"` and
`get_cas_from_name("Jasmolin II")` correctly returns `"1172-63-0"`. Across 25
sampled CAS numbers even the weaker "at least one of its names maps back" holds for
only 23. Its fallback clause compared a CAS against a list of *names*, so that could
never hold either. It passed only because its `SELECT ... LIMIT 1` had no
`ORDER BY`, and `test_create_indexes_*` — same file, same shared database file —
creates `idx_type ON api_ready_query(type)`, which changes which row comes back.
Replaced with a deterministic test of the guarantee that does hold.

**Two hygiene items this exposed, not fixed:**

- `tests/test_zeropm.py` has **52 `LIMIT` clauses and only 3 `ORDER BY`**. Every
  unordered `LIMIT` is a latent order-dependent flake of the kind above, because
  `test_create_indexes_*` mutates the shared database mid-run. Worth a sweep.
- `ZeroPM.create_indexes()` reports `'exists'` for indexes it has just created; it
  reports `'created'` only when `force=True`. Cosmetic, but misleading.

---

## 17. Landed on 2026-09-19 — cache correctness (partly §5 B.7, partly new)

Four defects in the caching layer, found while auditing the PubChem modules.
Three of them were **not** in this plan and are the reason the rest of it matters
less than it looked: the persistent cache was not merely untidy, it was inert.

### 17.1 The cache key embedded a memory address (new)

`CacheManager._get_cache_key` hashed `json.dumps(key_data, sort_keys=True,
default=str)` over the call's arguments. For a bound method `args[0]` is `self`,
and `str()` of a client object yields `<provesid.pubchemview.PubChemView object
at 0x784d4cb8b080>`. Every instance, and every interpreter, therefore produced a
different key. **No entry written by one run was ever reachable from the next**,
for any `@cached` method in `pubchem`, `pubchemview`, `resolver`,
`cascommonchem`, `classyfire` or `opsin` — while `/tmp/provesid_cache` filled
with one unreachable pickle per call.

Fixed with `cache.stable_key_part`, which normalises arguments into a
JSON-serialisable, process-stable form: an object reduces to the value it
declares through `__cache_key__()`, else to its fully qualified class name.
`PubChemAPI`, `PubChemView` and `NCIChemicalIdentifierResolver` declare
`__cache_key__` as `(class path, base_url)` — two clients on one endpoint share
entries, two endpoints stay apart. `default=str` is gone, so an un-normalised
argument can no longer slip through.

Verified end to end: a live `extract_property_data(2244, "Dissociation
Constants")` in one process, then the same call in a second process with
`_make_request` replaced by a raising stub — served from disk.

Existing entries are orphaned by the new scheme. Inert, not wrong;
`clear_all_service_caches()` reclaims the space.

### 17.2 Failed lookups were stored as answers (new)

`@cached` wrote whatever the function returned, including the
`{'success': False, 'error': ...}` dicts and empty lists the clients produce
after a 429, a PUG-View `ServerBusy` or a timeout — and it wrote them **even
when `use_cache=False`**. One transient error became a permanent "no data".

`@cached` now skips storage for any result `cache.is_failure_result` accepts,
and takes a `skip_if` predicate for other failure shapes. A raised exception was
already uncached; that is now stated in the docstring.

### 17.3 `PubChemView` reported a failed fetch as an absent property (new)

`extract_property_data` caught `PubChemViewError` — the base class, so
"failed after N attempts" included — logged `"Property not found"` and returned
`[]`. With 17.2 this cached a 503 as fact. Observed live: `get_boiling_point(2244)`
returned `[]` while `get_available_properties(2244)` listed Boiling Point as
present; on retry it returned four values.

Transport failures now propagate; only genuine absence yields `[]`.
`get_property_table` makes the same distinction. Same fix in
`PubChemAPI.get_compound_synonyms`, which swallowed every error into `[]`;
`get_compound_properties` now records a failed synonym fetch under
`synonyms_error` and is not cached while incomplete.

This exposed a second bug: PUG-View answers an unknown heading with **400
`PUGVIEW.BadRequest`**, which `_make_request` sent through `raise_for_status()`
into the retry loop — four requests to learn a permanent answer. 4xx is now
classified: 400 and 404 are absence, other non-429 4xx are non-retryable
errors, and only 429/5xx/timeouts/connection errors retry. Tuning the retry
budget and honouring `Retry-After` stay with §4 Workstream A.

### 17.4 A cached `None` was a cache miss (§5 B.7, as planned)

`_load_from_disk` returned `None` for both "no entry" and "the value is
`None`", so any function returning `None` re-ran every time. Now a `_MISS`
sentinel.

### 17.5 Tests and docs

New offline test files, 29 tests, no network: `tests/test_cache_correctness.py`
(key stability across instances and across a subprocess, `__cache_key__`
separation, the `None` sentinel, and that failures/exceptions/`skip_if` results
are not stored) and `tests/test_pubchem_failure_reporting.py` (absence vs
failure for both clients, and the 4xx/5xx classification). `docs/advanced_caching.md`
gained a "What Is and Is Not Cached" section; `CHANGELOG.md` records all four.

Also removed an unreachable `return cache_info` in `PubChemAPI.get_cache_info`.

### 17.6 What the two failing tests taught (and the fix that followed)

Running the full suite after 17.1–17.4 left two failures:
`test_pubchem.py::test_synonyms` and
`test_pubchemview.py::test_get_property_table`, both asserting a non-empty
result and getting an empty one. Cache inspection showed the stored entries were
**correct** (698 synonyms), so these were live empty responses, not stale reads.

That falsified an assumption in 17.3: that PubChem's 404 and 400 reliably mean
absence. Under load — the suite was competing with concurrent manual checks —
they do not. What *is* reliable is the fault code in the body, which both
services always send:

| Response | Meaning |
|---|---|
| `PUGREST.NotFound` / `PUGVIEW.NotFound` | genuinely no such data |
| `PUGVIEW.BadRequest` | no such heading |
| `...ServerBusy` / `...ServerError` / `...Timeout` | ask again later |

So classification now reads `Fault.Code` (`pubchem.fault_code`,
`pubchemview.fault_code`): a transient code is retried whatever status carries
it, absence codes are not retried, other 4xx are non-retryable errors.

And because no classification is worth trusting absolutely, **absence is no
longer cached at all**. A wrongly cached "no data" is permanent; re-fetching a
genuinely empty result costs one cheap request. Every cached extraction method
in `pubchemview` (17 of them) plus `get_compound_synonyms` carries a `skip_if`
predicate. The two raw fetchers keep caching unconditionally — they raise on
absence, so they never had an empty to store.

This is the reason 17.1 mattered more than it first appeared: with unstable keys
the poisoning was invisible, because no poisoned entry was ever read back. Fixing
the key made the pre-existing defect reproducible, which is how it got found.

### 17.7 Still open in this area

- The default cache directory is `tempfile.gettempdir()/provesid_cache`.
  `docs/advanced_caching.md` promises "cache survives restarts", which `/tmp`
  does not honour on most Linux systems. Moving it under
  `XDG_CACHE_HOME`/`%LOCALAPPDATA%` is a behaviour change, left for §5 B.7.
- `batch_extract_properties` still degrades a per-property failure to `[]`, by
  design, so one bad property does not abort a batch. Documented as such; a
  proper `errors` mapping fits better once Workstream A lands.
- The 14 `clear_<svc>_cache` / `get_<svc>_cache_info` pairs are untouched (§5 B.7).

## 18. Landed on 2026-09-19 — bulk and offline-first properties

Item 3 of the PubChem assessment: the property layer asked PubChem about one
compound at a time and never consulted the local database it ships with.

### 18.1 One request per compound (new)

`get_compound_properties_batch` looped over `get_compound_properties`, so N
compounds cost N requests plus the 0.2 s pause between each. PubChem's property
endpoint has always accepted a comma-separated CID list and answered the whole
set in one payload.

`PubChemAPI.get_properties_for_cids(cids, properties)` does that, in chunks of
`PROPERTY_CHUNK_SIZE` (200). Measured against the live API: 450 CIDs in 3
requests and 2.1 s, versus 450 requests and roughly 95 s for the loop. The chunk
is deliberately well below what PubChem tolerates — a 500-CID POST returns fine,
throttling green — because a smaller chunk wastes less work on a retry.

`get_compound_properties_batch` is now a reshaping layer over it and keeps its
return contract, with one improvement: a CID PubChem has no record of is
reported with `success=False` instead of being dropped, so the result can be
zipped against the input. Its `output_format` parameter is gone; only JSON was
ever reshapeable into per-CID dicts.

### 18.2 Long identifier lists need POST (new)

PUG-REST caps a URL at about 2000 characters. Two hundred CIDs is roughly 1400,
so the chunk size alone nearly reaches it, and any larger `chunk_size` a caller
passes would break. `_build_post_url` builds the endpoint without the identifier
segment and `_get_property_rows` switches to POST above
`URL_IDENTIFIER_LIMIT` (1600), identifiers in the body. Below the limit it stays
with GET, which PubChem prefers and which caches at intermediaries.

### 18.3 The local database was never consulted for properties (§2/§9, new)

`pubchem_id.db` holds 1.59M compounds and its `compounds` table carries 16
columns that are PubChem properties under another name — `mf`, `mw`,
`polararea`, `xlogp`, `exactmass` and so on. Nothing in the package read them:
every property lookup went to the network, including for the 1.59M compounds
already on disk.

`PubChemID.properties()`, `.properties_for_cids()` and `.properties_table()`
implement the two-stage lookup from §9. Routing rules:

- A CID with no local row goes online.
- A requested property with no local column sends the *whole* request online.
  Splitting the property list would cost the same traffic — that property needs
  a request either way — and would return a row assembled from two different
  PubChem snapshots.
- `use_online_fallback=False` is a hard guarantee; the tests assert no request
  is made.

Measured: 473 of 1000 low CIDs served from disk, the remaining 527 in three
batched requests, 1.9 s total. Strictly offline, 5000 CIDs resolved in 0.05 s.

Three details the schema forced:

- The `smiles` column is the **isomeric** SMILES, verified against CID 5793
  (D-glucose): the column carries the stereocentres, matching what PubChem now
  calls `SMILES`. `ConnectivitySMILES` is not stored and stays online-only.
- A NULL column means the compound genuinely has no such value, not a gap in the
  export. Checked against the live API for five CIDs whose `xlogp` is NULL
  (233, 234, 271, 533, 544): PubChem omits `XLogP` from its own answer for all
  five. So a NULL is reported as an absent key — matching PubChem, which omits
  a property rather than reporting it as null — and does not trigger a fallback.
- The two sources disagree on types: PUG-REST returns `MolecularWeight` and
  `ExactMass` as strings, SQLite as floats. `_cast_property` normalises both, or
  a table assembled from the two sources would not be usable as one table. An
  uncastable value passes through unchanged rather than being dropped.

### 18.4 Tests, docs, example

`tests/test_pubchem_properties.py` — 27 tests, fully offline. The transport is a
recorder, so a test asserts both the answer and *how many requests it took*,
which is the only way to test a batching or an offline-first decision. The local
database is a temporary SQLite file with the real schema, holding a complete
compound, one with a NULL `xlogp`, and deliberately omitting a CID so it can
stand for an online-only one.

`docs/api/pubchem.md` gained "Properties Without the Network";
`examples/pubchem_properties_demo.py` runs all of it against the live API and
prints which source answered.

### 18.5 Still open

- Item 4 of the assessment is untouched: `pubchemview` still has two divergent
  parsers and no numeric/unit normalisation, so an experimental property comes
  back as the string `"140 °C"`. That is the remaining sense in which "the API
  for retrieving chemical properties is incomplete".
- `PubChemID.properties()` takes CIDs only. Resolving a CAS or an InChIKey to a
  CID first is already offline (`cas_to_cid` and friends), but a caller has to
  chain the two calls.
- The offline route does not serve synonyms, though the database has a
  `synonyms` table; `get_compound_properties(include_synonyms=True)` is still
  online-only.

## 19. Landed on 2026-09-19 — one parser, typed values

Item 4 of the PubChem assessment, and the last part of "the api for retrieving
chemical properties is incomplete": the values came back as prose.

### 19.1 Two parsers that disagreed (new)

`pubchemview` had two:

- `_parse_value_string` (25 lines), behind `_extract_value_info` and therefore
  behind `PropertyData.unit` and every convenience getter. It knew six unit
  families and returned `(unit, conditions)` — no value at all.
- `_extract_experimental_value_and_unit` (390 lines), behind
  `get_property_table` only. A property-specific cascade — a branch each for
  vapor pressure, logP, dissociation constants, melting/boiling point, density,
  viscosity and solubility, each with four to ten regexes — returning four
  strings.

The same string parsed differently depending on which method a caller used, and
neither returned a number. `get_melting_point(2244)` handed back the string
`"138-140"`; comparing two compounds meant writing a third parser.

`src/provesid/pubchemview_parse.py` replaces both with `parse_value(text,
heading)` → `ParsedValue`. Four ordered patterns (range, number-then-unit,
labelled value, any number) do the work the cascade did, because the notations
that drove most of its branches are normalised first: `8.5X10-5` → `8.5e-5`,
`4,600` → `4600`, NBSP → space, U+2212 → hyphen. 415 lines of cascade became a
115-line parser plus two lookup tables.

### 19.2 What the typed output carries

`value` / `value_min` / `value_max`, `unit`, `value_si` / `value_min_si` /
`value_max_si` / `unit_si`, `temperature_c`, `operator`, `qualitative`,
`conditions`, and always `text`. Decisions worth recording:

- **A single value fills both bounds.** A caller filtering numerically should
  not have to ask which shape it got. A range leaves `value` None, so "is this
  one number" stays answerable.
- **Range ends convert separately.** A temperature conversion is affine, so
  138-140 °C is 411.15-413.15 K; scaling the range by the °C factor alone would
  be nonsense. That is why there are three `*_si` fields rather than one.
- **`temperature_c` is a condition, not the value.** `"2.47 cP at 20 °C"`
  reports a viscosity. Reading the condition's number as the value is the
  easiest way to get this wrong, so the clause is blanked out of the string
  before a value is looked for. Two spellings count as a condition: `at 20 °C`
  and the parenthesised label `(77 °F):`. A bare `(135 °C)` does not — the
  trailing colon is what distinguishes a label from a parenthesised value.
- **Nothing is invented.** An unrecognised unit is reported as written with
  `unit_si=None`. `%`, `ppm` and `ppb` are never converted: a composition needs
  a density to become a concentration, and a guess would be indistinguishable
  from a measurement. `M` for molar is not recognised at all — indistinguishable
  from metres, from the M of a molecular weight, and from a stray capital.

Two bugs the work surfaced, both from combining a permissive unit pattern with
short aliases, both now covered by a test:

- a `'c': '°C'` alias turned the leading letter of `cP` into a temperature, so
  viscosity in centipoise came back as 275.62 K. Single letters are now matched
  only as a whole token;
- case-insensitive matching of the alias `pa` read the `pa` of `parts` as
  pascals, so `"5 parts water"` became 5 Pa. A known spelling must now end
  where the token ends.

### 19.3 The heading argument (new)

Two things a string alone cannot settle, both from real data: a bare `138` under
"Melting Point" is °C (depositors routinely omit the unit), and a bare `1.19`
under "LogP" is dimensionless and complete rather than missing its unit.
`CELSIUS_HEADINGS` and `DIMENSIONLESS_HEADINGS` encode those; passing no heading
declines both hints, so the parser is usable on a string of unknown provenance.

### 19.4 Non-experimental subtrees were reported as absent (new)

Both response parsers walked a hard-coded path: `Record → Section[TOCHeading ==
"Chemical and Physical Properties"] → Section[TOCHeading == "Experimental
Properties"] → Section[*] → Information`. PUG-View nests a requested heading
wherever it sits in that compound's table of contents, which is not the same
place for every heading. Verified against the live service:

| heading | path |
| --- | --- |
| `Melting Point` | Chemical and Physical Properties → Experimental Properties |
| `GHS Classification` | Safety and Hazards → Hazards Identification |
| `Drug Indication` | Drug and Medication Information |
| `Absorption, Distribution and Excretion` | Pharmacology and Biochemistry |
| `Computed Properties` | Chemical and Physical Properties → Computed Properties → *one section per property* |

So `get_property_table(2244, "GHS Classification")` fetched a full response,
found nothing at the hard-coded path, and returned an empty frame — the same
answer it gives for a compound that genuinely has no such data. `_iter_information`
now walks the record and yields `(heading, item)` pairs, which also means the
heading reaching the parser is the innermost one that owns the value rather than
a value the caller passed in.

While walking, items that state no value are skipped: PubChem records some
sections as a pointer, `{"Value": {"ExternalTableName": "iupacpka"}}`, which
previously produced a blank row in the table. All seven "unparseable" strings in
the coverage measurement below were these.

### 19.5 Coverage measured on real data

Ten compounds × ten properties = 365 real value strings: **90% yield a number,
8% a qualitative term**, and the remainder stated no value at all and are now
filtered out. So every non-empty string PubChem returned for that sample parsed
into either a number or a term.

Of the numbers, 76% also carry an SI conversion. The other 24% are the
dimensionless quantities (logP, pKa, refractive index — where the number is the
whole answer) and the compositions (`%`, `ppm`), which is why the docs tell
callers to filter on `UnitSI` rather than assume `ValueSI` is populated.

The cross-check that the conversions are right: aspirin's melting point is
deposited variously as `135 °C`, `135 °C (rapid heating)`, `275 °F` and
`275 °F (NTP, 1992)`, and all of them now read 408.15 K. Ethanol comes out at
−114.0 °C, caffeine at 236.7 °C, ibuprofen at 76.0 °C — all correct.

### 19.6 Changing a cached return shape needs a cache version (new)

Noticed while re-measuring coverage: the first run still reported blank rows
that fresh fetches did not produce. They were cache entries written minutes
earlier, before the value-less pointers were filtered out.

That is the benign version of a real hazard. `PropertyData` gained two fields,
and pickle stores an instance's `__dict__`, so an entry written by the previous
version restores as an object with **no `parsed` attribute** — and every caller
that reads `data.parsed` would raise, on a machine where the only thing that
changed was the package version. Fixing the cache key in §17.1 is what made
these entries reachable in the first place, so this is a hazard that §17
created and §19 had to answer.

`PubChemView.CACHE_SCHEMA_VERSION` (now 2) is part of `__cache_key__`, so every
entry of the previous shape is unreachable. It is a constant to bump whenever a
cached return shape changes; a schema change is not a code change the cache can
detect by itself.

The four module-level convenience functions were a loophole in this: they were
`@cached` *and* delegated to cached methods, so they kept a duplicate payload
under a key with no `self` in it and therefore no version. Dropping their
decorators fixes the staleness and halves the disk they used.

### 19.7 Tests, docs, example

`tests/test_pubchemview_parse.py` — 69 tests, fully offline, each case a string
shape PubChem actually returns, so the file doubles as a record of what the data
looks like. Two tests in `test_pubchemview.py` that exercised the deleted
cascade were rewritten against the parser, and one column assertion now reads
`PROPERTY_TABLE_COLUMNS` rather than repeating the list.

`docs/api/pubchemview.md` — the "Data Structures" and "Advanced Usage" sections
documented a `PropertyData` with fields that never existed (`cid`,
`string_with_markup`, `reference_doi`, …) and a method
`experimental_properties_to_dataframe` that does not exist. Replaced with the
real dataclass, a `ParsedValue` field table, and worked examples.

`examples/pubchemview_parsed_values_demo.py` runs all of it live.

### 19.8 Still open

- `batch_extract_properties` still degrades a per-property failure to `[]` (see
  §17.7).
- The convenience getters still return `List[PropertyData]`; a caller wanting a
  number reaches through `.parsed`. A `melting_point_si(cid)`-style accessor
  would be the obvious next convenience layer, but it needs a policy for which
  of several entries to believe, which is a data-quality question rather than a
  parsing one.
- `PubChemID.properties()` covers computed properties offline; experimental
  ones have no offline source at all, so PUG-View is still always a network
  call.

---

## 20. Landed on 2026-09-19 — `Search` stopped targeting ZeroPM

Not planned as such: it came out of §5 B.2. Collapsing the five source ladders
first required knowing whether all five rungs deserve a vote. One does not.

### 20.1 An inventory harvest voting as a curated source

Since 0.6.0 `confidence` and `min_source_support` are driven by corroboration:
each source that carries a structure is one independent vote. ZeroPM was one of
the five, but it is not the same kind of thing as the other four — it aggregates
regulatory inventories rather than curating compounds, so its name→structure
rows are materially noisier while weighing exactly as much as ChEBI's. A wrong
ZeroPM row could out-vote a right one and, at `min_source_support ≥ 2`, admit a
structure no curated source backed.

`Search` now targets four sources — ChEBI, CompTox, PubChemID, ChEMBL — and
`use_zeropm=True` restores the previous five. The `ZeroPM` class is untouched
and stays fully available for direct use; only the resolver stopped consulting
it, and with it off, ZeroPM's database is never even opened.

### 20.2 "Disabled" must not depend on how the caller built the instance

A `zeropm=` client handed to the constructor is now ignored, with a warning,
unless `use_zeropm=True` is set too. Otherwise `Search(zeropm=client)` would
have quietly queried a source the default says is off.

`_SOURCE_KEYS` became a per-instance attribute holding the sources that instance
targets; the catalogue is `_ALL_SOURCE_KEYS` and the default set
`_DEFAULT_SOURCE_KEYS`. `_ensure_clients` iterates `_SOURCE_KEYS` instead of a
hard-coded list of five, and the consensus pass in `_finalize` does the same
rather than naming four keys inline.

### 20.3 What it costs: recall on typos, not precision

ZeroPM is the only source doing true fuzzy *retrieval* — the others are
substring-matched with `exact=False` — so it was the one that could reach a typo
sharing no usable substring with the real name. `"caffiene"` now returns no
match by default instead of caffeine.

That is a recall loss, and the regression suite now says so explicitly rather
than leaving it to be rediscovered. `test_search_precision_regression.py` runs
the misspelling table twice:

- `test_misspelled_name_resolves_to_the_right_compound` — with `use_zeropm=True`,
  asserting the typo still resolves (unchanged behaviour, now opt-in).
- `test_misspelling_never_resolves_to_a_wrong_compound` — on the default four,
  asserting a typo resolves correctly *or not at all*. Empty is an accepted
  outcome; a different compound never is.

Precision is what the four-source default protects, and that is the test that
protects it.

### 20.4 Consequences for callers

- `source_details` no longer carries a `"ZeroPM"` entry.
- `sources_available` lists four keys; the degraded-run warning says "of 4
  sources".
- Confidence values shift slightly wherever ZeroPM used to vote.

### 20.5 Tests and docs

`tests/test_search.py` gained `TestZeroPMOptIn` (6 tests: not targeted by
default, targeted when opted in, client ignored / kept, absent from reported
availability, contributes only when opted in) and two `source_details` cases.
The two source-availability fixtures in `test_search_precision_regression.py`
and `test_search_scoring_truth.py` stopped hard-coding the five client
attributes and now read `sources_unavailable`, so they follow the target list
instead of restating it.

208 tests pass across the three search files. `docs/api/search.md` gained a
"Targeted sources" section, and the CAS example's comment was corrected.

### 20.6 Still open

- §5 B.2 itself — the source ladders are still ten copies of the same shape, now
  four rungs deep instead of five, and `sources.py` does not exist.
- Fuzzy retrieval has no home among the default sources. If typo recall matters
  later, it wants a real fuzzy index over the curated names, not ZeroPM's rows
  as a proxy for one.
