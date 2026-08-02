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
