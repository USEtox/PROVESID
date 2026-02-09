# ChEMBL Database Interface

::: provesid.chembl.CheMBL
    options:
      show_root_heading: true
      show_source: false
      heading_level: 2
      members:
        - __init__
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

The ChEMBL module provides access to the ChEMBL SQLite database, a manually curated database of bioactive molecules with drug-like properties maintained by EMBL-EBI. The database contains over 2.3 million compounds with chemical structures, properties, and bioactivity data.

## Database Information

- **Database**: ChEMBL v36
- **Format**: SQLite (~5GB uncompressed)
- **Source**: [EMBL-EBI FTP](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/)
- **Auto-download**: Yes (on first use)
- **Compressed size**: ~1.5GB

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

Download from command line:

```bash
cd src/provesid/data
wget https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_36_sqlite.tar.gz
tar -xzf chembl_36_sqlite.tar.gz
mv chembl_36/chembl_36.db .
```

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
