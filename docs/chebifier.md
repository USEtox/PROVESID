# Chebifier taxonomy backend

PROVESID can classify molecules into **ChEBI ontology classes** using the
[ChEB-AI `chebifier`](https://github.com/ChEB-AI/python-chebifier) ensemble — an
offline, AI-based classifier that combines a transformer model, graph neural
networks, and rule-based models. In PROVESID this is exposed through
`provesid.taxonomy.ChebifierClassifier`.

`chebifier` is **heavy** (it pulls in PyTorch and, for the graph models, the PyG
stack) and its model weights are large, so it is an **optional** dependency. The
core `pip install provesid` stays light and does not require any of it.

---

## TL;DR — install

Linux / CPU, all models including the graph (GNN) models:

```bash
bash scripts/install_chebifier.sh
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

`chebifier` cannot be fully installed with a single `pip install`. Two problems
had to be solved, both verified by actually running the ensemble on CPU:

### 1. The graph stack pins torch to 2.11

The graph models depend on `chebai-graph`, which needs the PyG compiled
extensions `pyg_lib`, `torch_scatter`, `torch_sparse`, and **`torch_cluster`**.
These have **no source install** — they must come from the PyG wheel index that
matches the *exact* installed torch version and platform:

```
https://data.pyg.org/whl/torch-<version>+cpu.html
```

`torch_cluster`'s newest prebuilt CPU wheel is for **torch 2.11** — there is no
torch 2.12 wheel. So the whole stack is pinned to **torch 2.11.0** on CPU, and
the PyG extensions are installed *before* `chebai-graph`.

### 2. Some dependencies are git-only or undeclared

| Dependency | How it's installed | Needed for |
|---|---|---|
| `chebifier==1.2.1` | PyPI | the ensemble |
| `chebai==1.3.0` | PyPI | all model readers |
| `chebai-graph==1.1.0` | PyPI | graph (GNN) models |
| `chebi-utils` | PyPI (**undeclared** import of `chebai-graph`) | graph models |
| `chemlog-extra` | `git+https://github.com/ChEB-AI/chemlog-extra.git` | `chemlog_*` models |
| `c3p` | `git+https://github.com/sfluegel05/c3p.git` | `c3p` model |

The git dependencies are **not** placed in the `pyproject.toml` `chebifier`
extra, because PyPI rejects direct-URL (`git+…`) references in published package
metadata. The extra therefore contains only `chebifier==1.2.1`, and the script is
the documented step for everything else.

---

## The graph-model checkpoint issue (and its fix)

This is the subtle one. On a fresh install the ensemble *loads* every model, but
prediction crashes in the graph models with:

```
RuntimeError: mat1 and mat2 shapes cannot be multiplied (…x12 and 11x256)
```

**Root cause.** chebifier 1.2.1 ships GNN checkpoints ("v244") whose edge/node
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

This also explains why the hosted [chebifier web
app](https://chebifier.hastingslab.org/) works while a fresh install does not:
the web deployment's environment predates the drift, and the web repo's
`requirements.txt` is unpinned — so the "missing link" is not a version pin, it's
these index files.

**The fix.** Revert the three index files to their pre-drift (`677d44b`)
contents. This is done two ways, both idempotent:

* **Install script** — `install_chebifier.sh` patches them after install.
* **At runtime** — `ChebifierClassifier` calls
  `provesid.taxonomy.ensure_v244_indices()` before loading the ensemble, so it
  self-heals even if `chebai-graph` is later reinstalled. You can also call it
  directly:

  ```python
  from provesid.taxonomy import ensure_v244_indices
  ensure_v244_indices()   # {'BondType': 'patched', 'AtomNumHs': 'ok', ...}
  ```

**Upstream follow-up.** The clean long-term fix is for ChEB-AI to ship
checkpoint-matched index files with the model repos (or retrain). Until then, the
pin (`chebifier==1.2.1`) plus this patch is the reproducible combination.

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

- **Non-augmented `gat` model.** The plain `gat_chebi50_v244` model logs
  `failed to parse a SMILES string` for some inputs and abstains — a pre-existing
  bug in `chebai-graph`'s non-augmented reader, unrelated to the index fix. It
  degrades gracefully: the two augmented graph models (`gat-aug`,
  `resgated-aug`) plus the transformer and rule-based models carry the ensemble.
- **Per-label confidence.** chebifier's default `predict_smiles_list` returns the
  set of predicted ChEBI classes (not per-class probabilities), so the
  `confidence` column is populated only if a configuration returns scores.
- **Platform.** The install script targets **Linux/CPU**. On other platforms the
  PyG wheel tags and `torch_cluster` availability differ.

---

## See also

- `scripts/install_chebifier.sh` — the installer.
- `scripts/README.md` — installer notes.
- `plans/2026-07-02-chemical-taxonomy-classyfire-chebifier.md` — design plan
  (§10 covers the optional-dependency + storage + checkpoint-fix design in full).
