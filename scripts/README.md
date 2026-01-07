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
