# Chebifier taxonomy backend

PROVESID can classify molecules into **ChEBI ontology classes** using the
[ChEB-AI `chebifier`](https://github.com/ChEB-AI/python-chebifier) ensemble — an
offline, AI-based classifier that combines a transformer model, graph neural
networks, and rule-based models. In PROVESID this is exposed through
`provesid.taxonomy.ChebifierClassifier`.

`chebifier` is **heavy** (it pulls in PyTorch and, for the graph models, part of
the PyG stack) and its model weights are large, so it is an **optional**
dependency. The core `pip install provesid` stays light and does not require any
of it.

---

## TL;DR — install

Linux / CPU, all models including the graph (GNN) models:

```bash
bash scripts/install_chebifier.sh
```

That is a thin wrapper around two pip commands, which you can also run by hand:

```bash
uv pip install "chebifier[models]"
uv pip install torch==2.12.0 torch_scatter torch_geometric \
    -f https://data.pyg.org/whl/torch-2.12.0+cpu.html
uv pip install "provesid[chebifier]"
```

Installing torch from the **CPU** index instead is worth it — measured **1.6 GB**
of site-packages against **5.4 GB**, because plain PyPI `torch` adds 2.7 GB of
CUDA wheels plus a 700 MB triton. Do it first, so nothing pulls a CUDA torch that
is then thrown away (this is what the script does by default):

```bash
uv pip install torch==2.12.0 --index-url https://download.pytorch.org/whl/cpu
uv pip install "chebifier[models]"
uv pip install torch_scatter torch_geometric \
    -f https://data.pyg.org/whl/torch-2.12.0+cpu.html
uv pip install "provesid[chebifier]"
```

Then:

```python
from provesid.taxonomy import ChebifierClassifier

clf = ChebifierClassifier()                 # model weights download on first use
df = clf.classify(["c1ccccc1", "OCC1OC(O)C(O)C(O)C1O"])   # benzene, glucose
print(df[["smiles", "chebi_ids"]])
```

The first `classify()` call downloads the model weights once into the shared
PROVESID dataset directory (see [Model storage](#model-storage)); subsequent runs
reuse them and hit the on-disk cache.

---

## Why an install script?

Much less than it used to be. **chebifier 1.2.2** ships a `models` extra that
pins the entire model stack, and every piece of it is now on PyPI — so the git
installs, the `chebi-utils` workaround and the index patch that earlier PROVESID
releases needed are all gone. What the script still does:

### 1. `torch_scatter` must come from the PyG wheel index

The graph models depend on `chebai-graph`, which needs the compiled PyG extension
**`torch_scatter`** (plus the pure-python `torch_geometric`). `torch_scatter` has
**no source install** — it must come from the PyG wheel index matching the
*exact* torch version and platform:

```
https://data.pyg.org/whl/torch-<version>+cpu.html
```

This cannot be expressed in `pyproject.toml`, hence a script (or the two manual
commands above).

`chebai-graph` 1.0.0 needs **only** `torch_scatter` — not `torch_sparse`,
`torch_cluster` or `pyg_lib`. `torch_cluster` was the package with no wheel past
torch 2.11 and the sole reason the stack used to be pinned to torch 2.11; without
it, **torch 2.12** works (verified end-to-end on CPU).

### 2. It defaults to the CPU torch build and verifies the result

Plain `pip install torch` pulls the CUDA wheels (5.4 GB of site-packages); the
script installs torch from `https://download.pytorch.org/whl/cpu` first (1.6 GB)
and then imports the whole stack with the *target* interpreter, to fail loudly on
a partial install. Set `TORCH_INDEX_URL=""` to opt out.

### What gets installed

| Dependency | Version | Comes from | Needed for |
|---|---|---|---|
| `chebifier[models]` | 1.2.2 | PyPI | the ensemble + all model deps below |
| `chebai` | 1.2.0 | `models` extra | electra transformer + model readers |
| `chebai-graph` | 1.0.0 | `models` extra | graph (GNN) models |
| `chemlog-extra` | 1.0.1 | `models` extra | `chemlog_*` models |
| `c3p` | 0.5.0 | `models` extra | `c3p` model |
| `torch` | 2.12.0 (+cpu) | PyTorch CPU index | everything |
| `torch_scatter`, `torch_geometric` | 2.1.2+pt212cpu, 2.8.0 | PyG wheel index | graph models |

Because the whole model stack is now expressible as PyPI requirements, the
`pyproject.toml` extra carries it directly:

```toml
chebifier = ["chebifier[models]==1.2.2"]
```

`pip install 'provesid[chebifier]'` therefore gets the transformer and
rule-based models working on its own; only the graph models need the extra PyG
step.

---

## The graph-model checkpoint issue (historical — handled by the pin)

Applies to `chebai-graph` **1.1.0+**. With `chebifier[models]==1.2.2`, which pins
`chebai-graph==1.0.0`, this does not occur — `ensure_v244_indices()` reports
every index as `"ok"` on a clean install. Kept here because the failure is
otherwise baffling and reappears the moment `chebai-graph` is upgraded by hand.

With a drifted `chebai-graph`, the ensemble *loads* every model but prediction
crashes in the graph models with:

```
RuntimeError: mat1 and mat2 shapes cannot be multiplied (…x12 and 11x256)
```

**Root cause.** chebifier 1.2.x ships GNN checkpoints ("v244") whose edge/node
feature widths are fixed (e.g. `gat-aug` edge dimension = 11). Those widths are
produced by `chebai-graph`'s one-hot **property index vocabularies**, stored as
text files inside the package:

```
.../chebai_graph/preprocessing/bin/<Property>/indices_one_hot.txt
```

A `chebai-graph` commit (`ea77f36`, **2026-03-02** — *after* chebifier 1.2.1
shipped on 2026-02-17) **appended** tokens to three of them, widening the graph
feature vectors by one and breaking the checkpoints:

| index file | v244 (pre-drift) | drifted |
|---|---|---|
| `BondType` | `DATIVE, SINGLE, AROMATIC, TRIPLE, DOUBLE` (5) | + `UNSPECIFIED` (6) |
| `AtomNumHs` | `0,3,2,4,1,5,6` (7) | + `8` (8) |
| `NumAtomBonds` | `0,1,2,4,5,3,6,8,7,10,12` (11) | + `11,9` (13) |

Because the tokens were *appended* (not inserted), the earlier one-hot positions
stay aligned, so it is a clean **dimension** mismatch rather than scrambled
predictions. Node-feature drift is absorbed by chebifier's internal padding;
edge features are **not** padded, which is why `BondType` is what actually
crashes the GNNs.

This also explained why the hosted [chebifier web
app](https://chebifier.hastingslab.org/) kept working while a fresh unpinned
install did not: the web deployment's environment predates the drift.

**The fix.** Two layers, both idempotent:

* **The pin** — `chebifier[models]==1.2.2` resolves `chebai-graph==1.0.0`
  (released 2025-12-08, before the 2026-03-02 drift), whose index files already
  match the v244 checkpoints. This is the primary fix and needs no patching.
* **At runtime, as a safety net** — `ChebifierClassifier` calls
  `provesid.taxonomy.ensure_v244_indices()` before loading the ensemble, which
  restores the three index files if a drifted `chebai-graph` is ever installed
  over the pin. You can also call it directly:

  ```python
  from provesid.taxonomy import ensure_v244_indices
  ensure_v244_indices()   # {'BondType': 'ok', 'AtomNumHs': 'ok', 'NumAtomBonds': 'ok'}
  ```

**Upstream follow-up.** The clean long-term fix is for ChEB-AI to ship
checkpoint-matched index files with the model repos (or retrain), so newer
`chebai-graph` releases can be used with these checkpoints.

---

## Model storage

Like every other large PROVESID dataset (ChEBI SDF, ChEMBL/CompTox SQLite, …),
chebifier's model weights are stored **once per machine**, shared across virtual
environments, rather than re-downloaded into each one. `ChebifierClassifier`
redirects the Hugging Face / torch caches into the shared PROVESID dataset
directory:

- Default location: `user_dataset_path("chebifier")` (platform-specific via
  `platformdirs`).
- Override with the `PROVESID_DATA_DIR` environment variable, or the
  `data_dir=` constructor argument.
- Pre-existing `HF_HOME` / `HF_HUB_CACHE` / `TORCH_HOME` environment variables
  set by the user are respected and never overridden.

---

## Usage

```python
from provesid.taxonomy import ChebifierClassifier, chebifier_available

# Feature-detect without importing torch:
if not chebifier_available():
    raise SystemExit("run scripts/install_chebifier.sh first")

clf = ChebifierClassifier(
    data_dir=None,        # default: shared PROVESID dataset dir
    use_cache=True,       # cache predictions on disk by InChIKey
    resolve_names=False,  # set True to also resolve ChEBI ids -> names
)

df = clf.classify(["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", "OCC1OC(O)C(O)C(O)C1O"])
# columns: inchikey, smiles, kingdom, superclass, class, subclass,
#          chebi_ids, chebi_names, source, confidence
# (kingdom/superclass/class/subclass are None for this backend — they are the
#  ClassyFire levels of the shared taxonomy schema.)

# Collapse to a {inchikey: label} mapping for downstream grouping:
labels = ChebifierClassifier.to_labels(df, level="chebi_ids")
```

For one-off calls there is a convenience wrapper (constructs and discards a
classifier, so prefer the class for repeated use):

```python
from provesid.taxonomy import classify_chebifier
df = classify_chebifier(["c1ccccc1"])
```

### Caching

Predictions are cached on disk under the `chebifier` cache service, keyed by
**InChIKey + chebifier version + configuration**. A re-run over the same
structures is served from cache without loading the model. Manage it with:

```python
from provesid.cache import clear_chebifier_cache, get_chebifier_cache_info
get_chebifier_cache_info()
clear_chebifier_cache()
```

---

## Known limitations

- **Non-augmented `gat` model.** With the older stack, the plain
  `gat_chebi50_v244` model logged `failed to parse a SMILES string` for some
  inputs and abstained (a `chebai-graph` reader bug). Not reproduced with
  `chebifier[models]==1.2.2` / `chebai-graph==1.0.0`, but it degrades gracefully
  in any case: the augmented graph models plus the transformer and rule-based
  models carry the ensemble.
- **Per-label confidence.** chebifier's default `predict_smiles_list` returns the
  set of predicted ChEBI classes (not per-class probabilities), so the
  `confidence` column is populated only if a configuration returns scores.
- **Platform.** The install script targets **Linux/CPU**. On other platforms the
  PyG wheel tags differ — check what exists under
  <https://data.pyg.org/whl/> for your torch version.

---

## See also

- `scripts/install_chebifier.sh` — the installer.
- `scripts/README.md` — installer notes.
- `plans/2026-07-02-chemical-taxonomy-classyfire-chebifier.md` — design plan
  (§10 covers the optional-dependency + storage + checkpoint-fix design in full).
