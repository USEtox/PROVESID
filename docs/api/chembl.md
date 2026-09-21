# ChEMBL Database Interface

::: provesid.chembl.CheMBL
    options:
      show_root_heading: true
      show_source: false
      heading_level: 2
      members:
        - __init__
        - compact
        - download_database
        - search_by_chembl_id
        - search_by_name
        - search_by_inchi
        - search_by_inchikey
        - search_by_smiles
        - get_compound
        - get_properties
        - chembl_id_to_molregno
        - molregno_to_chembl_id

::: provesid.chembl.ChEMBLError
    options:
      show_root_heading: true
      show_source: false
      heading_level: 2

## Overview

The ChEMBL module provides access to the ChEMBL SQLite database, a manually curated database of bioactive molecules with drug-like properties maintained by EMBL-EBI. The database contains over 2.9 million compounds (release 37) with chemical structures, properties, and bioactivity data.

## Database Information

- **Database**: current ChEMBL release, resolved from `latest/` (v37 as of 2026-08)
- **Format**: SQLite (~27.7 GiB uncompressed for release 37)
- **Source**: [EMBL-EBI FTP](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/)
- **Auto-download**: Yes (on first use)
- **Compressed size**: ~5.8 GB (release 37)
- **Installed size**: ~2.4 GiB. `CheMBL(source="sqlite")`, the default, compacts
  the release into the eight tables PROVESID reads and deletes the rest;
  `CheMBL(source="full")` keeps all 74 tables at ~27.7 GiB.
- **Free disk needed**: ~33.4 GiB during installation by either of those routes —
  the archive and the full release both exist before either is removed.
  `CheMBL(source="mysql")` needs ~4.5 GiB instead: it downloads ChEMBL's
  2.1 GB MySQL dump and builds the same extract from it directly, so the full
  release is never written (see [Installing from the MySQL dump](#installing-from-the-mysql-dump)).

## Key Features

- **Multiple search methods**: Search by ChEMBL ID, name, InChI, InChI Key, or SMILES
- **Local database**: Fast queries with no API rate limits
- **Offline access**: Works offline after initial download
- **Comprehensive data**: Structures, properties, identifiers, and metadata
- **Easy integration**: Simple Python API consistent with other PROVESID modules

## Database Schema

The ChEMBL class queries the following main tables:

### Core Tables

- **molecule_dictionary**: Primary compound information
    - `molregno`: Internal molecule registry number (primary key)
    - `chembl_id`: ChEMBL identifier (e.g., CHEMBL25)
    - `pref_name`: Preferred compound name
    - `max_phase`: Maximum clinical trial phase (0-4)
    - `therapeutic_flag`: Drug/therapeutic indicator
    - `molecule_type`: Type classification

- **compound_structures**: Chemical structure representations
    - `molregno`: Foreign key to molecule_dictionary
    - `canonical_smiles`: Canonical SMILES string
    - `standard_inchi`: Standard InChI representation
    - `standard_inchi_key`: Standard InChI Key
    - `molfile`: Molfile structure data

- **compound_properties**: Physicochemical properties
    - `molregno`: Foreign key to molecule_dictionary
    - `mw_freebase`: Molecular weight
    - `alogp`: Calculated LogP (lipophilicity)
    - `hba`: Hydrogen bond acceptors
    - `hbd`: Hydrogen bond donors
    - `psa`: Polar surface area
    - `rtb`: Rotatable bonds
    - `aromatic_rings`: Number of aromatic rings
    - `heavy_atoms`: Heavy atom count
    - `num_ro5_violations`: Lipinski Rule of Five violations

- **molecule_synonyms**: Alternative compound names
    - `molregno`: Foreign key to molecule_dictionary
    - `synonyms`: Synonym/alternative name
    - `syn_type`: Type of synonym
    - **Note**: All compound lookups automatically include a list of synonyms

- **chembl_id_lookup**: ChEMBL ID mapping table
    - `chembl_id`: ChEMBL identifier
    - `entity_type`: Type of entity (COMPOUND, TARGET, ASSAY, etc.)
    - `entity_id`: Internal ID (e.g., molregno for compounds)

For complete schema documentation, see `src/provesid/data/schema_documentation.txt`.

## Quick Start

```python
from provesid import CheMBL

chembl = CheMBL()

# Search by ChEMBL ID
compound = chembl.search_by_chembl_id('CHEMBL25')  # Aspirin
print(compound['pref_name'])  # 'ASPIRIN'

# Get molecular properties
props = chembl.get_properties(compound['molregno'])
print(f"MW: {props['mw_freebase']:.2f}")
print(f"LogP: {props['alogp']:.2f}")

# View synonyms
print(f"Synonyms: {compound['synonyms'][:5]}")  # First 5 synonyms
```

## Usage Examples

### Example 1: Search by Name

```python
from provesid import CheMBL

chembl = CheMBL()

# Search for compounds by name
results = chembl.search_by_name('caffeine')

for compound in results:
    print(f"{compound['chembl_id']}: {compound['pref_name']}")
    print(f"  SMILES: {compound['canonical_smiles']}")
    if compound['synonyms']:
        print(f"  Synonyms: {', '.join(compound['synonyms'][:3])}")
```

### Example 2: Structure-Based Search

```python
from provesid import CheMBL

chembl = CheMBL()

# Search by SMILES
smiles = 'CC(=O)Oc1ccccc1C(=O)O'
compound = chembl.search_by_smiles(smiles)

# Search by InChI Key
inchikey = 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'
compound = chembl.search_by_inchikey(inchikey)

# Both return the same compound (aspirin)
print(compound['chembl_id'])  # 'CHEMBL25'
```

### Example 3: Property Analysis

```python
from provesid import CheMBL

chembl = CheMBL()

# Get compound and properties
compound = chembl.search_by_chembl_id('CHEMBL25')
props = chembl.get_properties(compound['molregno'])

# Check Lipinski's Rule of Five
print("Lipinski's Rule of Five:")
print(f"  MW < 500: {props['mw_freebase'] < 500}")
print(f"  LogP < 5: {props['alogp'] < 5}")
print(f"  HBA < 10: {props['hba'] < 10}")
print(f"  HBD < 5: {props['hbd'] < 5}")
print(f"  Total violations: {props['num_ro5_violations']}")
```

### Example 4: Batch Processing

```python
from provesid import CheMBL
import pandas as pd

chembl = CheMBL()

# Process multiple compounds
drug_ids = ['CHEMBL25', 'CHEMBL521', 'CHEMBL112']
data = []

for chembl_id in drug_ids:
    compound = chembl.search_by_chembl_id(chembl_id)
    if compound:
        props = chembl.get_properties(compound['molregno'])
        data.append({
            'ChEMBL ID': chembl_id,
            'Name': compound['pref_name'],
            'MW': props['mw_freebase'],
            'LogP': props['alogp']
        })

df = pd.DataFrame(data)
print(df)
```

### Example 5: ID Conversion

```python
from provesid import CheMBL

chembl = CheMBL()

# Convert ChEMBL ID to internal molregno
molregno = chembl.chembl_id_to_molregno('CHEMBL25')
print(f"CHEMBL25 -> molregno: {molregno}")

# Convert back
chembl_id = chembl.molregno_to_chembl_id(molregno)
print(f"molregno {molregno} -> {chembl_id}")
```

## Manual Database Download

If you prefer to download the database manually:

```python
from provesid import CheMBL

# Initialize without auto-download
chembl = CheMBL(auto_download=False)

# Or download explicitly
chembl.download_database(force=True)
```

Download from command line. `latest/` holds only the current release, so check
which archive it advertises rather than assuming a release number:

```bash
cd src/provesid/data
curl -s https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/ \
  | grep -o 'chembl_[0-9]*_sqlite.tar.gz' | sort -u
# e.g. chembl_37_sqlite.tar.gz
wget https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_37_sqlite.tar.gz
tar -xzf chembl_37_sqlite.tar.gz
# the .db sits inside the extracted tree; its depth varies between releases
find chembl_37 -name 'chembl_37.db' -exec mv {} . \;
rm -r chembl_37
```

## Installing from the MySQL dump

ChEMBL also publishes each release as a plain-text MySQL dump,
`chembl_NN_mysql.tar.gz` (2.1 GB). `CheMBL(source="mysql")` downloads that
instead of the SQLite archive and reads the eight tables PROVESID uses straight
out of it as it is decompressed; the other 66 tables stream past unparsed. No
MySQL server is involved.

```python
from provesid import CheMBL

chembl = CheMBL(source="mysql")        # 2.1 GB down, ~4.5 GiB free needed
chembl.provenance["source_format"]     # 'mysql'

# ...or from a dump you downloaded yourself
path = CheMBL.build_from_mysql_dump("chembl_37_mysql.tar.gz", remove_source=True)
```

| route | download | free disk at peak | installed |
|---|---:|---:|---:|
| `source="sqlite"` (default) | 5.8 GB | 33.4 GiB | 2.4 GiB |
| `source="mysql"` | 2.1 GB | ~4.5 GiB | 2.4 GiB |
| `source="full"` | 5.8 GB | 33.4 GiB | 27.7 GiB |

The result is the same extract as the `sqlite` route: same tables, columns,
column types and indexes. `CheMBL.extract_digest()` fingerprints each table
(row count plus a content hash that includes each value's SQLite type), which
is how the two routes are compared:

```python
a = CheMBL.extract_digest("chembl_37_provesid.db")      # built from the SQLite release
b = CheMBL.extract_digest("from_mysql/chembl_37_provesid.db")
assert a == b
```

The reader understands exactly what `mysqldump` writes and refuses anything
else, naming the table and position, rather than guessing a value. If a future
ChEMBL dump changes shape, the error says so and `source="sqlite"` still works.

## Release Handling

`CheMBL` does not pin a release number. On first use it reads the `latest/`
directory listing, downloads the archive it advertises, and names the local
database after it (`chembl_37_sqlite.tar.gz` → `chembl_37.db`):

```python
from provesid import CheMBL

CheMBL.resolve_latest_db_url()
# 'https://.../latest/chembl_37_sqlite.tar.gz'

chembl = CheMBL()
chembl.release   # 37
chembl.db_path   # .../chembl_37.db
```

- Any `chembl_*.db` already in the data directory is reused — a new ChEMBL
  release does not trigger a multi-gigabyte re-download. Pass `redownload=True` to move to
  the newest release.
- If the listing cannot be read, `CheMBL.DEFAULT_DB_URL` (pinned to release
  `CheMBL.FALLBACK_RELEASE`) is used and a warning is logged.
- `db_url=`, `db_name=`, and `db_path=` still override everything.

## Performance Notes

- **First query**: May take a few seconds as SQLite loads indexes
- **Subsequent queries**: Very fast (local SQLite, no network overhead)
- **Name searches**: May be slower due to LIKE queries and synonym matching
- **Exact matches**: InChI Key and SMILES searches are fast (indexed)

## Error Handling

```python
from provesid import CheMBL, ChEMBLError

try:
    chembl = CheMBL()
    compound = chembl.search_by_chembl_id('CHEMBL25')
    
    if compound is None:
        print("Compound not found")
    else:
        print(f"Found: {compound['pref_name']}")
        
except ChEMBLError as e:
    print(f"ChEMBL error: {e}")
except FileNotFoundError as e:
    print(f"Database not found: {e}")
```

## Comparison with Other Data Sources

| Feature | ChEMBL | PubChem | ChEBI |
|---------|--------|---------|-------|
| Database Size | 2.3M compounds | 110M+ compounds | 190K+ entities |
| Focus | Bioactive drugs | All chemistry | Biology-focused |
| API | Local SQLite | REST API | REST API + SDF |
| Offline | Yes | No | Partial (SDF) |
| Speed | Very fast | Rate limited | Moderate |
| Bioactivity | Yes | Yes | Limited |

## References

- ChEMBL Database: [https://www.ebi.ac.uk/chembl/](https://www.ebi.ac.uk/chembl/)
- ChEMBL Documentation: [https://chembl.gitbook.io/chembl-interface-documentation/](https://chembl.gitbook.io/chembl-interface-documentation/)
- Database Downloads: [https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/)
- Schema Documentation: [https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/schema_documentation.html](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/schema_documentation.html)

## See Also

- [PubChem API](pubchem.md) - For broader compound coverage
- [ChEBI](chebi.md) - For biological entities and ontology
- [NCI Resolver](nci_resolver.md) - For identifier resolution
