# PubChem ID Database

## Overview

The `pubchem_id.db` is a SQLite database containing ~1.6M PubChem compounds with their identifiers and chemical properties. This database is built from the PubChem_CAS_202601.csv file and provides fast local lookup for identifier conversion.

**The database can be automatically downloaded from Zenodo on first use, or built locally from the CSV file.** This version of the database is based on the `csv` file downloaded from [pubchem](https://pubchem.ncbi.nlm.nih.gov/classification/#hid=72) for compounds that have CAS entries in the Pubchem datasets. You can download it yourself by going to "Names and Identifiers -> Other Identifiers -> CAS".  

## Quick Start

The database will auto-download on first use:

```python
from provesid import PubChemID

# First time: automatically downloads from Zenodo
db = PubChemID()  # Downloads database if not found

# Or disable auto-download
db = PubChemID(auto_download=False)  # Raises error if not found
```

## Hosting on Zenodo

### Uploading to Zenodo (One-time setup)

1. **Create Zenodo account**: Go to [zenodo.org](https://zenodo.org) and sign up
2. **Create new upload**: Click "New upload" button
3. **Upload database file**: 
   - File: `pubchem_id.db` (~2.2 GB)
   - Location: `src/provesid/data/pubchem_id.db`
4. **Add metadata**:
   - **Title**: "PubChem ID Database - Chemical Identifier Lookup Database"
   - **Description**: 
     ```
     SQLite database containing 1.6M PubChem compounds with identifiers 
     (CID, CAS, InChI, InChIKey, SMILES) and chemical properties 
     (molecular formula, molecular weight, LogP, complexity, etc.).
     
     Built from PubChem CAS dataset (January 2026).
     Used by PROVESID package for fast local identifier conversion.
     ```
   - **Keywords**: PubChem, chemical identifiers, CAS numbers, InChI, InChIKey, SMILES
   - **License**: CC BY 4.0 or CC0 (Public Domain)
   - **Upload type**: Dataset
5. **Publish**: Click "Publish" to make it public
6. **Get download URL**: Copy the direct download URL (format: `https://zenodo.org/records/18173204/files/pubchem_id.db`)
7. **Update code**: Edit [src/provesid/pubchem.py](../pubchem.py) and replace the zenodo_url in `PubChemID.download_database()` method

### Manual Download

If needed, you can download manually:

```python
from provesid import PubChemID

# Download to default location
PubChemID.download_database()

# Or specify custom location
PubChemID.download_database(db_path='/path/to/pubchem_id.db')

# Or with custom Zenodo URL
PubChemID.download_database(zenodo_url='https://zenodo.org/record/XXXXXX/files/pubchem_id.db')
```

## Database Structure

### Tables

#### 1. compounds (1,589,910 rows)
Main compound data with identifiers and chemical properties.

**Columns:**
- `cid` (INTEGER, PRIMARY KEY) - PubChem Compound ID
- `cmpdname` (TEXT) - Common compound name
- `iupacname` (TEXT) - IUPAC systematic name
- `inchi` (TEXT) - International Chemical Identifier
- `smiles` (TEXT) - Simplified Molecular Input Line Entry System
- `inchikey` (TEXT) - Hashed InChI (27 characters)
- `mf` (TEXT) - Molecular formula
- `mw` (REAL) - Molecular weight
- `exactmass` (REAL) - Exact mass
- `polararea` (REAL) - Topological polar surface area
- `complexity` (REAL) - Molecular complexity score
- `xlogp` (REAL) - Partition coefficient (log P)
- `heavycnt` (INTEGER) - Heavy atom count
- `hbonddonor` (INTEGER) - Hydrogen bond donor count
- `hbondacc` (INTEGER) - Hydrogen bond acceptor count
- `rotbonds` (INTEGER) - Rotatable bond count
- `charge` (INTEGER) - Formal charge
- `cidcdate` (TEXT) - Creation date (YYYYMMDD)

**Indexes:**
- `idx_compounds_inchikey` - Fast lookup by InChIKey
- `idx_compounds_inchi` - Fast lookup by InChI
- `idx_compounds_mf` - Fast lookup by molecular formula

#### 2. cas_numbers (1,400,544 rows)
CAS Registry Numbers with one-to-many relationship to compounds.

**Columns:**
- `id` (INTEGER, PRIMARY KEY) - Auto-increment ID
- `cid` (INTEGER) - PubChem Compound ID (foreign key to compounds.cid)
- `cas` (TEXT) - CAS Registry Number (e.g., "50-78-2")

**Indexes:**
- `idx_cas_cas` - Fast lookup by CAS number
- `idx_cas_cid` - Fast lookup by compound ID

**Notes:**
- One compound can have multiple CAS numbers
- 1,323,167 compounds have at least one CAS number

#### 3. synonyms (13,987,682 rows)
Chemical synonyms with one-to-many relationship to compounds.

**Columns:**
- `id` (INTEGER, PRIMARY KEY) - Auto-increment ID
- `cid` (INTEGER) - PubChem Compound ID (foreign key to compounds.cid)
- `synonym` (TEXT) - Chemical synonym or alternative name

**Indexes:**
- `idx_synonyms_synonym` - Fast text search by synonym
- `idx_synonyms_cid` - Fast lookup by compound ID

**Notes:**
- Excludes CAS numbers, InChI, and InChIKey (stored in separate tables)
- Average ~8.8 synonyms per compound

## Building the Database

The database is built from `PubChem_CAS_202601.csv` using the build script:

```bash
python scripts/build_pubchem_id_db.py
```

**Process:**
1. Reads the CSV file (~2 GB, 1.6M compounds)
2. Extracts CAS numbers, InChI, and InChIKey from the `cmpdsynonym` column
3. Creates separate tables for compounds, CAS numbers, and synonyms
4. Builds indexes for fast lookups
5. Outputs `pubchem_id.db` (~2.2 GB)

**Build time:** ~1.5 minutes

## Usage

Use the `PubChemID` class to query the database:

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
cas_list = ["50-78-2", "50-00-0", "64-17-5"]
df = db.get_by_cas_batch(cas_list)
print(df[['cas', 'cmpdname', 'mf', 'mw']])

# Search by name
results = db.search_by_name("aspirin", exact=False)

# Search by formula
results = db.search_by_formula("C9H8O4")

# Get statistics
stats = db.get_stats()
print(f"Total compounds: {stats['total_compounds']:,}")
```

## Database Statistics

- **Total compounds:** 1,589,910
- **Total CAS numbers:** 1,400,544
- **Compounds with CAS:** 1,323,167 (83%)
- **Total synonyms:** 13,987,682
- **Compounds with InChIKey:** 1,589,910 (100%)
- **Database size:** 2.2 GB

## Notes

- The CSV file and database are excluded from version control (see `.gitignore`)
- Database must be built locally before using the `PubChemID` class
- Only identifiers and chemical properties are included (not annotations or bioassay data)
- CAS numbers are extracted from the `cmpdsynonym` field using regex pattern matching
