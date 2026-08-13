# Database Build Scripts

This directory contains scripts for building local databases used by PROVESID.

## PubChem ID Database

### build_pubchem_id_db.py

Builds a SQLite database from PubChem CAS CSV file for fast local identifier lookup.

**Input:** `src/provesid/data/PubChem_CAS_202601.csv` (~2 GB, 1.6M compounds)

**Output:** `src/provesid/data/pubchem_id.db` (SQLite database)

**What it does:**
1. Extracts CAS numbers, InChI, and InChIKey from the `cmpdsynonym` column
2. Creates separate tables for compounds, CAS numbers, and synonyms
3. Includes chemical properties: molecular formula, molecular weight, LogP, complexity, etc.
4. Builds indexes for fast lookups by CAS, InChIKey, InChI, formula, and synonym

**Usage:**
```bash
cd c:\projects\git\PROVESID
python scripts/build_pubchem_id_db.py
```

**Processing time:** ~10-15 minutes (depends on system)

**Database structure:**
- `compounds` table: Main compound data with identifiers and properties
- `cas_numbers` table: CAS Registry Numbers (one-to-many relationship)
- `synonyms` table: Chemical synonyms (one-to-many relationship)

**Indexes created:**
- CAS number lookup
- InChIKey lookup
- InChI lookup
- Molecular formula lookup
- Synonym search

**Expected output:**
```
Processing: c:\projects\git\PROVESID\src\provesid\data\PubChem_CAS_202601.csv
Output database: c:\projects\git\PROVESID\src\provesid\data\pubchem_id.db
Reading CSV file...
Counting rows...
Processing 1,589,912 compounds...
Processing compounds: 100%|████████████| 1589912/1589912
Creating indexes...

✓ Database created successfully!
  - 1,589,912 compounds
  - XXX,XXX CAS numbers
  - XXX,XXX synonyms
  - Database size: X.XX GB
```

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

After building the database, use the `PubChemID` class:

```python
from provesid import PubChemID

# Initialize
db = PubChemID()

# Lookup by CAS
result = db.get_by_cas("50-78-2")  # Aspirin
print(result['inchi'])

# Convert identifiers
cid = db.cas_to_cid("50-78-2")
inchikey = db.cas_to_inchikey("50-78-2")

# Batch operations
results = db.batch_cas_to_cid(["50-78-2", "50-00-0"])

# Get identifier table
df = db.get_id_table_from_cas("50-78-2")
```

## Notes

- The CSV file and SQLite database are excluded from git (see `.gitignore`)
- The database is ~1-2 GB depending on content
- First-time build required before using `PubChemID` class
- Database includes only identifiers and chemical properties (not annotations or bioassay data)
