# Database Build Scripts

This directory contains scripts for building local databases used by PROVESID.

## PubChem ID Database

### build_pubchem_id_db.py

Builds `pubchem_id.db` from a monthly snapshot of PubChem's FTP site. A thin
command-line wrapper over `provesid.pubchem_ftp.build_pubchem_id_db`, which is
also what `PubChemID()` runs when the database is missing; the script exists
for refreshing the copy on Zenodo that `PubChemID(source="zenodo")` downloads.

**Usage:**
```bash
python scripts/build_pubchem_id_db.py --list                  # available releases
python scripts/build_pubchem_id_db.py                         # newest snapshot
python scripts/build_pubchem_id_db.py --release 2026-09-01 --out ./pubchem_id.db
```

**Options:** `--no-inchi` (skip the 7.4 GB InChI file and compute InChI with
RDKit), `--no-synonyms`, `--keep-downloads`, `--force`.

**Cost:** 15.4 GB transferred, one file at a time; a 2.5 GB database, plus
7.4 GB free at the worst moment; about 12 minutes of processing on top of the
download. See `docs/api/pubchem.md#the-local-database` for what
goes in and how, and `PubChemID(db_path=...).provenance()` for what a finished
database records about itself.

**Refreshing Zenodo:** build, check `provenance()["release"]`, upload the file
to a new version of the Zenodo record, and point `PubChemID.DEFAULT_DB_URL` at
it.

## chebifier backend installer

### install_chebifier.sh

Installs the optional `chebifier` AI-classification backend for PROVESID on
**Linux/CPU**, with **all models working incl. the graph/GNN models**. Verified
end-to-end (benzene/aspirin/glucose → sensible ChEBI classes).

**Why a script** — since chebifier **1.2.2** upstream ships a `models` extra that
pins the whole model stack (all on PyPI), so this is now essentially two pip
commands:

```bash
uv pip install "chebifier[models]"
uv pip install torch==2.12.0 torch_scatter torch_geometric \
    -f https://data.pyg.org/whl/torch-2.12.0+cpu.html
```

The script wraps those because `torch_scatter` is a compiled extension with no
source install — it must come from the PyG wheel index matching the exact torch
version, which cannot be expressed in `pyproject.toml`. It also installs torch
from the PyTorch **CPU** index (1.6 GB of site-packages, vs 5.4 GB when plain PyPI
torch adds 2.7 GB of CUDA wheels plus triton) and verifies every model module
imports.

Install order: torch 2.12 (CPU) → `chebifier[models]` → `torch_scatter` +
`torch_geometric` from the PyG index → `provesid[chebifier]` → verify.

**Usage:**
```bash
bash scripts/install_chebifier.sh
```

**Env vars:** `TORCH_VERSION` (2.12.0), `CHEBIFIER_VERSION` (1.2.2),
`TORCH_INDEX_URL` (PyTorch CPU index; set to `""` for plain PyPI), `PIP`
(defaults to `uv pip`, else `python -m pip`).

**Notes:**
- `chebai-graph` 1.0.0 needs only `torch_scatter` — no `torch_sparse`,
  `torch_cluster` or `pyg_lib`. `torch_cluster` (no wheel past torch 2.11) was
  what previously pinned the stack to torch 2.11.
- No index patching needed anymore: `chebai-graph==1.0.0` predates the property
  index drift that broke the `v244` GNN checkpoints.
  `provesid.taxonomy.ensure_v244_indices()` still runs as a safety net and
  reports `ok` on a clean install. Root-cause writeup: `docs/chebifier.md`.
- Model weights are **not** installed here; they download on first
  `BaseEnsemble()` use into the shared PROVESID dataset dir
  (`PROVESID_DATA_DIR` to override). See the taxonomy plan §10.2.

## Using the Database

```python
from provesid import PubChemID

db = PubChemID()                     # builds from FTP if pubchem_id.db is missing
db.cas_to_cid("50-78-2")             # 2244
db.get_by_cas_batch(["50-78-2", "50-00-0"])
db.provenance()["release"]           # which PubChem snapshot answered
db.xrefs(2244)                       # DSSTox, ChEBI, ChEMBL, EC and UNII links
```
