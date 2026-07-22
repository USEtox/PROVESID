# Plan: Chemical taxonomy classification (ClassyFire + chebifier)

**Date:** 2026-07-02
**Author:** Ali A Eftekhari + Claude Code
**Target module:** `src/provesid/taxonomy.py` (new), building on `classyfire.py`, `chebi.py`, `cache.py`
**Status:** Partially implemented — **chebifier backend DONE** (`ChebifierClassifier`, §9 phase 2);
**ClassyFire backend + `ChemicalTaxonomy` facade to be added later** (§9 phases 1 & 3, not yet started).
**Consumer:** PROVES (`proves.taxonomy`, `proves.compare.compare_models_by_group`) — see §6

---

## 1. Problem statement

PROVES needs to attach a **chemical-class taxonomy label** to every chemical it evaluates so it
can report **per-chemical-class QSAR/QSPR model performance** (per-class RMSE) and, later, do
class-aware best-model selection. This must cover **the whole EXProves chemical space**
(tens of thousands of substances), not just a curated subset, and the labels must be **cached
once and reused** across runs.

Two classification sources are wanted, and both belong in PROVESID (alongside identifier
resolution, `ClassyFireAPI`, and `ChEBI`):

1. **ClassyFire** — a well-established, hierarchical structural taxonomy
   (kingdom → superclass → class → subclass). Accurate, but its public server is **slow and
   frequently down**, so batch classification of a large set is painful and must be resumable
   and aggressively cached.
2. **chebifier** — an open-source, **AI-based ChEBI-ontology classifier** that runs
   **offline**. It is the practical backend for classifying at the scale of all of EXProves;
   ClassyFire results are reused where they already exist.

PROVESID today exposes only the **low-level** `ClassyFireAPI` (`submit_query` / `query_status`
/ `get_query`) and a `ChEBI` web-API client. There is **no** high-level "give me class labels
for these N chemicals" capability, no chebifier backend, and no unified taxonomy output. This
plan adds them.

> Scope note: this chemistry taxonomy is **separate** from PROVES's organic/PFAS/pesticide
> filtering (`proves.classify`). It is a finer taxonomy used only for analysis/grouping.

---

## 2. Design decisions

- **One tidy output schema, regardless of backend.** `classify(...)` returns a
  `pandas.DataFrame`, one row per chemical:

  | column | meaning |
  |--------|---------|
  | `inchikey` | InChIKey (the cache/merge key) |
  | `smiles` | the structure classified (as submitted / canonicalised) |
  | `kingdom`, `superclass`, `class`, `subclass` | ClassyFire levels (names, may be null) |
  | `chebi_ids` | `|`-joined ChEBI class IDs (chebifier) |
  | `chebi_names` | `|`-joined ChEBI class names (chebifier, resolved via `ChEBI`) |
  | `source` | `"classyfire"` or `"chebifier"` |
  | `confidence` | per-label probability where the backend provides one (chebifier) |

  ClassyFire fills the four level columns; chebifier fills `chebi_ids` / `chebi_names`. A
  `to_labels(df, level=...)` helper collapses this to the `{inchikey: label}` dict PROVES wants.

- **InChIKey-keyed, on-disk, resumable cache.** A chemical is classified **once ever** and
  reused. Reuse the existing `provesid` cache (`cache.py`, the `@cached` decorator ClassyFire
  already uses) keyed by **InChIKey** (structure-stable), not by the transient ClassyFire
  `query_id`. Partial progress survives a crash / server outage: re-running skips already-cached
  chemicals.

- **Import existing saved results.** ClassyFire results already produced for USEtox3 (and any
  other saved JSON) can be **loaded straight into the cache** so they are never re-queried.

- **Backends are pluggable and independent.** `ClassyFireClassifier` and `ChebifierClassifier`
  implement the same `.classify(chemicals) -> DataFrame`. chebifier is an **optional extra**
  (`provesid[chebifier]`) so the core install stays light and offline-capable without it; its
  heavy PyTorch stack and model weights are only pulled in when a user opts in (§10).

- **chebifier model weights live in the shared systemwide dataset dir.** Like every other large
  PROVESID dataset (ChEBI SDF, ChEMBL/CompTox SQLite, …), chebifier's downloaded model weights
  and any auxiliary data are stored **once per machine** under `user_dataset_path()` /
  `PROVESID_DATA_DIR`, not re-downloaded into each virtualenv's site-packages or the default
  `~/.cache/huggingface`. See §10.2. This keeps a single copy shared across every Python env on
  the host.

- **Native taxonomies are kept native.** ClassyFire and ChEBI are different ontologies; we do
  **not** force a mapping between them. Each backend returns its own taxonomy columns.
  (A coarse harmonised label is a possible later addition — see §7.)

- **Batch + retry + backoff for ClassyFire.** Submit in batches, poll `query_status` with
  backoff, tolerate `In Queue` / server-down by re-queuing, and cache each entity by InChIKey as
  results arrive.

---

## 3. Proposed API (`provesid.taxonomy`)

```python
from provesid.taxonomy import ChemicalTaxonomy, classify_classyfire, classify_chebifier

# Facade — pick a backend, classify a list/DataFrame of chemicals, get the tidy table.
tax = ChemicalTaxonomy(backend="chebifier")          # or "classyfire"
labels_df = tax.classify(smiles_list, inchikeys=...)  # uses + fills the on-disk cache

# Convenience functions (thin wrappers):
df = classify_classyfire(smiles_list, inchikeys=None, levels=("class",), use_cache=True)
df = classify_chebifier(smiles_list, inchikeys=None, use_cache=True)

# Collapse to the {inchikey: label} map PROVES consumes:
labels = ChemicalTaxonomy.to_labels(df, level="class")     # ClassyFire level, or
labels = ChemicalTaxonomy.to_labels(df, level="chebi_names")

# Seed the cache from previously saved ClassyFire JSON (e.g. USEtox3):
ChemicalTaxonomy.import_classyfire_json("usetox3_classyfire.json")
```

New helpers/classes:

- `ChemicalTaxonomy(backend, cache=True, level=..., data_dir=None)` — facade with `.classify()`,
  `.to_labels()`, `.import_classyfire_json()`, `.cache_info()`. `data_dir` defaults to
  `user_dataset_path("chebifier")` for the chebifier backend (§10.2), matching the `data_dir`
  convention of `ChebiSDF`/`ChEMBL`/`CompTox`.
- `ClassyFireClassifier` — wraps `ClassyFireAPI`: batch `submit_query` (SMILES, one per line) →
  poll `query_status` until `Done` → `get_query(format="json")` → parse entities → per-level
  labels; cache each entity by InChIKey; resumable.
- `ChebifierClassifier` — wraps chebifier's `predict(smiles) -> {chebi_id: probability}`; resolve
  ChEBI IDs to names with the existing `ChEBI` client; cache by InChIKey.
- `parse_classyfire_entities(results, levels)` — pure function (mirrors the parser PROVES
  already has, moved/duplicated here so PROVESID owns it).

Export the facade + convenience functions from `provesid/__init__.py`.

---

## 4. Backends

### 4.1 ClassyFire (`ClassyFireClassifier`)
- Input: SMILES (preferred) or InChI; `submit_query(label, input=<newline-joined SMILES>,
  type="STRUCTURE")`. Respect ClassyFire's batch-size limit (chunk large sets).
- Poll `query_status(query_id)` with exponential backoff; treat `In Queue`/`In Progress`/HTTP
  errors as "retry later" (do not fail the whole batch). Persist `query_id` so a resumed run can
  re-poll instead of re-submitting.
- On `Done`: `get_query(query_id, "json")`, parse `entities[*]` → `{inchikey, kingdom,
  superclass, class, subclass}` (each level is a `{"name": ...}` node, possibly null), write each
  to the InChIKey cache.
- Uptime reality: this is for **incremental** classification + reuse, not a fast bulk pass.

### 4.2 chebifier (`ChebifierClassifier`)
- The offline, AI-based ChEBI classifier: the ChEB-AI **`chebifier`** package
  (`pip install chebifier`, MIT, [github.com/ChEB-AI/python-chebifier](https://github.com/ChEB-AI/python-chebifier)).
  It is an **ensemble** (deep-learning + rule-based + generative models) covering **1,742 ChEBI
  classes** and runs fully offline once its model weights are present.
- **Real interface** (verified against the repo):

  ```python
  from chebifier import BaseEnsemble
  ensemble = BaseEnsemble()                          # loads/downloads model weights
  preds = ensemble.predict_smiles_list(smiles_list)  # list, one entry per SMILES
  # each entry is a prediction of ChEBI classes for that structure (None if it failed to parse)
  ```

  `ChebifierClassifier` wraps this: build one `BaseEnsemble` per process (expensive — construct
  lazily and reuse), call `predict_smiles_list` in batches, normalise each entry into
  `chebi_ids` (+ `confidence` where the ensemble exposes a per-class score), keep predictions
  above a configurable probability threshold, and resolve ChEBI IDs → names via the existing
  `ChEBI` client (cached). Cache each result by **InChIKey** (§5).
- **Heavy, torch-based, GPU-optional.** The package pulls in a deep-learning stack (PyTorch) and
  downloads model weights from Hugging Face on first use. This is why it is an **optional extra**
  (§10) and why its weights/caches must be redirected to the shared PROVESID dataset directory
  (§10.2) rather than duplicated per virtualenv.
- This is the **scalable** backend: classify all of EXProves once, cache, reuse.

---

## 5. Caching

- Reuse `provesid.cache` (the `@cached(service=...)` mechanism ClassyFire already uses); add
  `service='taxonomy'` (or split `'chebifier'`) to the `_service_caches` registry in `cache.py`,
  plus the matching `clear_*_cache()` / `get_*_cache_info()` helpers (mirrors every other
  service). Cache key = **InChIKey + backend + level set** for ClassyFire, and **InChIKey +
  backend + chebifier model-version + threshold** for chebifier — so an upstream model bump or a
  threshold change produces a fresh entry instead of silently reusing stale labels.
- **This label cache is distinct from the model-weights store (§10.2).** The `@cached` pickle
  cache holds small per-chemical *results*; the multi-GB chebifier *weights* live in the shared
  dataset dir. Keep them separate so `clear_taxonomy_cache()` never nukes the expensive weights.
- `ChemicalTaxonomy.import_classyfire_json(path)` seeds the cache from saved results.
- `cache_info()` / `clear_cache()` exposed for visibility (mirrors the existing cache helpers).

---

## 6. How PROVES will consume it

PROVES already has the downstream plumbing (`proves/taxonomy.py`,
`proves.compare.compare_models_by_group`); it currently has **stubs** for the live backends:

```python
# proves/taxonomy.py (today) — these NotImplementedError stubs get wired to PROVESID:
def classyfire_classifier(level="class"): ...   # -> provesid.classify_classyfire
def chebifier_classifier(level="class"): ...    # -> provesid.classify_chebifier
```

Flow once this feature lands:

1. PROVES builds its **evaluation anchor** (every chemical with experimental data, with SMILES).
2. PROVES calls `provesid` to classify them (chebifier for the full set; ClassyFire reused where
   saved) → `{inchikey: label}` via `ChemicalTaxonomy.to_labels`.
3. `proves.taxonomy.attach_class_labels` adds a `chem_class` column.
4. `proves.compare.compare_models(..., group_column="chem_class")` reports **per-class RMSE**,
   stored in the decision JSON. (Class-aware *selection* stays deferred — PROVES decision QG.)

So PROVESID owns *classification*; PROVES owns *grouping + per-class metrics*.

---

## 7. Open questions

- ~~**chebifier packaging** — is it pip-installable with bundled weights, or does it need a
  separate model download / GPU?~~ **Resolved (see §10):** `pip install chebifier` works;
  weights are **not** bundled but auto-downloaded from Hugging Face on first `BaseEnsemble()`
  construction; PyTorch-based, **CPU is sufficient** (GPU optional). Shipped as the
  `provesid[chebifier]` extra with weights redirected to the shared dataset dir.
- ~~**chebifier version pinning**~~ **Decided:** pin **`chebifier==1.2.1`** (2026-02-17, latest)
  for reproducible cached labels. Three model components (chebai-graph, chemlog-extra, c3p) are
  installed manually from GitHub via `uv pip` (§10.1). Revisit only if a future release ships
  breaking changes; a deliberate bump invalidates the `backend + model-version` cache namespace (§5).
- **Taxonomy harmonisation** — do we ever need a single label space across ClassyFire and ChEBI
  (e.g. to compare per-class RMSE from mixed sources)? Default: no; report per backend. A coarse
  harmonised `class` mapping could be added later if needed.
- **Grouping level** — which level PROVES groups on (`class` vs `superclass`) is a PROVES
  concern (QG); PROVESID just returns all levels.

---

## 8. Testing

- `parse_classyfire_entities` — pure unit tests on a saved JSON fixture (levels, null levels,
  `InChIKey=` prefix stripping).
- Cache round-trip — classify → cache → re-classify hits cache (no network); `import_classyfire_json`
  seeds and is reused.
- chebifier — a few known chemicals (e.g. benzene, glucose) to a stable ChEBI class;
  `@pytest.mark.skipif(not chebifier_available())`. Add a `chebifier` pytest marker (mirrors the
  existing `classyfire`/`chebi` markers in `pyproject.toml`).
- chebifier **missing-extra** path — with `chebifier` uninstalled, `ChemicalTaxonomy(backend=
  "chebifier")` raises `ChebifierMissingError` with the install hint (no raw `ModuleNotFoundError`);
  this test runs even in the core environment.
- chebifier **storage redirect** — `_configure_chebifier_storage(tmp_dir)` sets `HF_HOME`/
  `HF_HUB_CACHE`/`TORCH_HOME` under the dataset dir and honors a pre-set env var (no override) and
  `PROVESID_DATA_DIR`; a pure env-var test, no model download required.
- ClassyFire live path — a tiny submit→poll→parse smoke test, **skipped** if the server is
  unreachable (it usually is).

---

## 9. Implementation phases

1. `parse_classyfire_entities` + `ClassyFireClassifier` (batch, cached, resumable) + tidy schema
   + `import_classyfire_json`. Seed cache from existing USEtox3 results. **(not yet implemented)**
2. ✅ **DONE — `ChebifierClassifier` (optional extra).** Implemented in `src/provesid/taxonomy.py`:
   `provesid[chebifier]` extra pinned to `chebifier==1.2.1` in `pyproject.toml`; the full
   Linux/CPU model stack (incl. graph models) installed via `scripts/install_chebifier.sh`;
   guarded lazy import (`_load_chebifier`, `ChebifierMissingError`, `chebifier_available()`);
   `_configure_chebifier_storage` redirecting HF/torch caches to `user_dataset_path("chebifier")`;
   `ensure_v244_indices()` self-healing the checkpoint/index drift (§10.4); lazy single
   `BaseEnsemble` instance + batched `predict_smiles_list`; InChIKey-keyed resumable cache
   (`chebifier` cache service + `clear_chebifier_cache`/`get_chebifier_cache_info`); tidy
   `DataFrame` schema (`TAXONOMY_COLUMNS`) + `to_labels`; optional ChEBI-name resolution; exports
   from `provesid/__init__`; `chebifier` pytest marker; tests in `tests/test_taxonomy.py`; docs at
   `docs/chebifier.md`; example at `examples/chebifier/chebifier_example.py`. **Verified end-to-end**
   on CPU (all models incl. GNN), with systemwide storage + cache round-trip.
3. `ChemicalTaxonomy` facade (unifying ClassyFire + chebifier behind one `backend=` selector) +
   shared `to_labels` + `import_classyfire_json`. **(deferred — chebifier's own `to_labels` and
   `classify_chebifier` cover the chebifier path today; the facade lands with phase 1.)**
4. (In PROVES) wire `proves.taxonomy.classyfire_classifier` / `chebifier_classifier` to these
   PROVESID entry points; enable `compare_models_by_group` per-class RMSE end to end.

---

## 10. chebifier as an optional dependency with systemwide model storage

This section is the detailed design for the two hard constraints on chebifier: **it must not
bloat the core install**, and **its large model weights must be stored once per machine** (shared
across virtualenvs), exactly like the other big PROVESID datasets.

### 10.1 Optional-dependency packaging

- **Extra, not a core dependency.** Pin the current latest release, **`chebifier==1.2.1`**
  (released 2026-02-17). Add to `pyproject.toml`:

  ```toml
  [project.optional-dependencies]
  chebifier = [
      "chebifier==1.2.1",   # pulls in torch + the ChEB-AI ensemble
  ]
  ```

  Core `pip install provesid` stays light and fully offline-capable; users who want AI
  classification run `pip install "provesid[chebifier]"`. Everything else in PROVESID keeps
  working without torch installed. We stay on 1.2.1 for now; if a later release brings breaking
  changes we decide then whether to bump (a bump invalidates the `model-version` cache namespace,
  §5).

- **The full model stack (incl. graph/GNN models) needs extra deps + a checkpoint fix — handled
  by an install script.** Per the python-chebifier README, "not all models can be installed
  automatically." Base `chebifier==1.2.1` runs the transformer + rule-based models, but the graph
  models (`gat-aug`, `resgated-aug`) need `chebai`, `chebai-graph`, `chebi-utils`, and the PyG
  compiled stack. This was **verified by actually running the full ensemble on CPU in a throwaway
  env** (benzene/aspirin/glucose → sensible ChEBI classes). We target **Linux/CPU** and ship a
  one-shot installer, **`scripts/install_chebifier.sh`**, that reproduces the working recipe and
  verifies both imports and a prediction. Two non-obvious constraints it encodes (see §10.4):

  1. **torch is pinned to 2.11.** The graph stack needs `torch_cluster`, whose newest prebuilt CPU
     wheel is for **torch 2.11** (there is **no** torch 2.12 wheel). All PyG extensions
     (`pyg_lib`/`torch_scatter`/`torch_sparse`/`torch_cluster`) come from the PyG wheel index
     matching the exact torch version + platform, and must be installed **before** `chebai-graph`.
  2. **chebai-graph's property index files are patched** to the v244-checkpoint-matching state
     (§10.4) — without this the GNN models fail to load with a tensor-shape error.

  ```bash
  bash scripts/install_chebifier.sh   # torch 2.11 cpu + full PyG + chebai(-graph) + chebi-utils
                                       #   + chemlog-extra/c3p (git) + provesid[chebifier] + index patch
  ```

  The git deps (`chemlog-extra`, `c3p`) are **not** placed in the `pyproject.toml` extra: PyPI
  rejects direct-URL (`git+…`) references in published dependency metadata. `chebi-utils` is an
  **undeclared** import of `chebai-graph`, so the script installs it explicitly. The script (plus a
  note in the mkdocs taxonomy page + `data/README.md`) is the documented setup step. The ensemble
  degrades gracefully when a member is absent, so the feature is usable in a reduced form (Electra
  + ChemLog + C3P) even without the PyG stack.

- **Lazy, guarded import.** `chebifier` is imported **inside** `ChebifierClassifier`, never at
  module top level, so importing `provesid.taxonomy` (or `provesid`) never requires torch. If the
  extra is missing, raise a clear, actionable error — do **not** let a raw `ModuleNotFoundError`
  escape (dev-principle #8):

  ```python
  class ChebifierMissingError(ProvesidError):
      """Raised when the optional 'chebifier' extra is not installed."""

  def _load_chebifier():
      try:
          from chebifier import BaseEnsemble
      except ImportError as exc:
          raise ChebifierMissingError(
              "The chebifier backend requires the optional extra. Install with:\n"
              "    pip install 'provesid[chebifier]'"
          ) from exc
      return BaseEnsemble
  ```

- **Feature detection.** Expose `provesid.taxonomy.chebifier_available() -> bool` (a cheap
  `importlib.util.find_spec("chebifier") is not None`) so callers (and the PROVES side) can
  branch, and so tests can `pytest.mark.skipif(not chebifier_available())`.

- **`ChemicalTaxonomy(backend="chebifier")`** fails fast with `ChebifierMissingError` at
  construction when the extra is absent, rather than deep inside a batch run.

### 10.2 Systemwide model-weights storage (the core requirement)

chebifier auto-downloads its ensemble weights from Hugging Face on first `BaseEnsemble()`. By
default those land in `~/.cache/huggingface` — which is per-user but **not** under PROVESID's
control, and some setups scatter caches per-env. We redirect them into the same shared dataset
root every other PROVESID dataset uses, so there is exactly **one copy per machine**, overridable
with `PROVESID_DATA_DIR`.

- **Single dataset subdir.** Reserve `user_dataset_path("chebifier")` (resolves under
  `platformdirs` `user_data_dir("provesid", "USEtox")`, or `$PROVESID_DATA_DIR/chebifier`). This
  mirrors how `ChebiSDF`, `ChEMBL`, `CompTox`, `ZeroPM`, and `PubChem` already default their
  `data_dir` to `user_dataset_path()`.

- **Redirect Hugging Face's cache into it.** Before the guarded import constructs `BaseEnsemble`,
  point the HF env vars at the shared dir **if the user has not already set them** (respect an
  explicit user override):

  ```python
  def _configure_chebifier_storage(data_dir: str | None = None) -> str:
      base = data_dir or user_dataset_path("chebifier")     # honors PROVESID_DATA_DIR
      os.environ.setdefault("HF_HOME", os.path.join(base, "huggingface"))
      os.environ.setdefault("HF_HUB_CACHE", os.path.join(base, "huggingface", "hub"))
      # torch.hub weights, if any component uses them:
      os.environ.setdefault("TORCH_HOME", os.path.join(base, "torch"))
      return base
  ```

  Set these **before** `chebifier`/`torch`/`huggingface_hub` are imported (they read the env at
  import time), i.e. inside `_load_chebifier()` prior to `from chebifier import BaseEnsemble`.

- **`data_dir` passthrough for parity with other classes.** `ChebifierClassifier(data_dir=...)`
  and `ChemicalTaxonomy(..., data_dir=...)` accept an explicit directory, defaulting to
  `user_dataset_path("chebifier")` — same signature convention as `ChebiSDF`, `ChEMBL`, etc. If
  chebifier's own API exposes a weights/config path argument, thread `data_dir` through to it as
  well as (or instead of) the env-var redirect.

- **First-run download UX.** The first `BaseEnsemble()` triggers a multi-hundred-MB download.
  Log it at INFO ("Downloading chebifier model weights to <path> (first run only)…") via the
  standard `logging` module (never `print`), consistent with the ChEBI SDF downloader. Weights
  persist and are reused by every env thereafter.

- **Docs + `.gitignore`.** Extend `src/provesid/data/README.md` with a "chebifier" entry stating
  the weights live in the shared dataset dir (`PROVESID_DATA_DIR` to override) and are downloaded
  on first use; confirm the shared dataset dir is already git-ignored (it lives outside the repo).

### 10.3 Runtime cost management

- **One ensemble per process.** `BaseEnsemble()` construction (weight load) is expensive; build
  it lazily on the first `classify()` and cache the instance on the classifier. Never construct
  per-chemical or per-batch.
- **Batch `predict_smiles_list`.** Feed structures in chunks; combine with the InChIKey label
  cache (§5) so already-classified chemicals never reach the model. For a re-run over EXProves,
  the ensemble ideally loads only if there are genuinely-new structures.
- **CPU default.** Do not assume CUDA; chebifier runs on CPU. If a GPU is present the underlying
  torch stack will use it, but nothing in PROVESID should require it.

### 10.4 Getting the graph (GNN) models to actually run — the checkpoint/index fix

The graph models did **not** work out of the box from a fresh install, and diagnosing why took
running the ensemble for real. Documenting it here so it isn't rediscovered painfully.

**Symptom.** With a current install, `BaseEnsemble()` loads every model but prediction crashes in
the graph models with `RuntimeError: mat1 and mat2 shapes cannot be multiplied (…x12 and 11x256)`
(and analogous for the plain `gat`).

**Root cause.** chebifier 1.2.1 ships GNN checkpoints ("v244") whose edge/node feature widths are
baked in (gat edge_dim=7, gat-aug edge_dim=11). Those widths are produced by `chebai-graph`'s
one-hot **property index vocabularies** (`.../chebai_graph/preprocessing/bin/<Property>/indices_one_hot.txt`).
chebai-graph commit **`ea77f36` (2026-03-02, *after* chebifier 1.2.1 shipped 2026-02-17)**
*appended* tokens to three of them — `BondType` (+`UNSPECIFIED`), `AtomNumHs` (+1), `NumAtomBonds`
(+2) — widening the graph models' feature vectors by one and breaking checkpoint compatibility.
Because the tokens were **appended** (not inserted), earlier one-hot positions stay aligned, so it
is a clean dimension mismatch rather than scrambled predictions. Node-feature drift is absorbed by
chebifier's padding (`max_len_node_properties`); **edge features are not padded**, so `BondType` is
what actually crashes the GNNs. This is also why the deployed chebifier **web app works** (its
frozen environment predates the drift) while a fresh unpinned install does not — the web repo's
`requirements.txt` is unpinned, so the "missing link" is not a version pin but the index files.

**Fix (baked into the installer).** After install, revert the three index files inside the
installed `chebai-graph` to their pre-drift (commit `677d44b`) contents, which match the v244
checkpoints:

| index file | pre-drift (v244) | drifted |
|---|---|---|
| `BondType/indices_one_hot.txt` | `DATIVE, SINGLE, AROMATIC, TRIPLE, DOUBLE` (5) | +`UNSPECIFIED` (6) |
| `AtomNumHs/indices_one_hot.txt` | `0,3,2,4,1,5,6` (7) | +`8` (8) |
| `NumAtomBonds/indices_one_hot.txt` | `0,1,2,4,5,3,6,8,7,10,12` (11) | +`11,9` (13) |

The script embeds these contents (tiny + immutable given the `chebifier==1.2.1` pin), so the patch
needs no network. **Verified:** with the patch, a pristine `install_chebifier.sh` env runs the
full ensemble on CPU — `gat-aug` and `resgated-aug` execute at 100% and the ensemble returns
sensible ChEBI classes for benzene/aspirin/glucose.

**Known minor limitation.** The *non-augmented* `gat` model logs "failed to parse a SMILES string"
and abstains (a pre-existing bug in chebai-graph's non-aug reader, unrelated to the index fix). It
degrades gracefully — the two augmented graph models plus Electra/ChemLog/C3P carry the ensemble.

**Upstream follow-up (optional).** The clean long-term fix is for ChEB-AI to ship
checkpoint-matched index files with the model repos (or retrain), removing the need for this patch.
Worth an upstream issue; until then the pin (`chebifier==1.2.1`) + index patch is the reproducible
combination.
