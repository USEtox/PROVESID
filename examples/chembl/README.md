# ChEMBL Examples

This directory contains examples and tutorials for using the ChEMBL interface in PROVESID.

## What is ChEMBL?

ChEMBL is a manually curated database of bioactive molecules with drug-like properties, maintained by the European Bioinformatics Institute (EMBL-EBI). It contains chemical, bioactivity, and genomic data to aid drug discovery.

## Files

- `chembl_tutorial.md` - Comprehensive tutorial covering all ChEMBL functionality
- `compact_demo.py` - Shrink a ChEMBL release from ~30 GB to ~2.6 GB
- `name_search_demo.py` - Name and synonym lookup, exact and substring
- `demo_implementation.py` - Minimal end-to-end usage

## Quick Start

```python
from provesid import CheMBL

# Initialize (auto-downloads database if needed)
chembl = CheMBL()

# Search for aspirin by ChEMBL ID
compound = chembl.search_by_chembl_id('CHEMBL25')
print(compound['pref_name'])  # 'ASPIRIN'
print(f"Synonyms: {', '.join(compound['synonyms'][:3])}")

# Get molecular properties
props = chembl.get_properties(compound['molregno'])
print(f"Molecular Weight: {props['mw_freebase']}")
print(f"LogP: {props['alogp']}")
```

## Database Information

- **Database**: current ChEMBL release SQLite, resolved from `latest/` (v37 as of 2026-08)
- **Size**: ~30 GB uncompressed, ~5.8 GB compressed (release 37; grows each release).
  `CheMBL.compact()` reduces an installed release to ~2.6 GB — see
  `compact_demo.py`.
- **Source**: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/
- **Auto-download**: Yes (on first use)

## Key Features

1. **Multiple Search Methods**:
   - ChEMBL ID
   - Compound name or synonym, exact or partial (`search_by_name`)
   - InChI / InChI Key
   - SMILES

2. **Data Retrieval**:
   - Complete compound information
   - Physicochemical properties (MW, LogP, HBA, HBD, PSA, etc.)
   - Chemical structures (SMILES, InChI, InChIKey)

3. **Local Database**:
   - Fast queries (local SQLite, no API calls)
   - Offline access after initial download
   - No rate limits

## Example Workflows

### Workflow 1: Search by Name and Get Properties

```python
from provesid import CheMBL

chembl = CheMBL()

# Search by name
results = chembl.search_by_name('caffeine')
compound = results[0]

# Get properties
props = chembl.get_properties(compound['molregno'])
print(f"MW: {props['mw_freebase']}")
print(f"HBA: {props['hba']}")
print(f"HBD: {props['hbd']}")
```

### Workflow 2: Structure-Based Lookup

```python
from provesid import CheMBL

chembl = CheMBL()

# Search by SMILES
smiles = 'CC(=O)Oc1ccccc1C(=O)O'
compound = chembl.search_by_smiles(smiles)

print(f"ChEMBL ID: {compound['chembl_id']}")
print(f"Name: {compound['pref_name']}")
print(f"InChI Key: {compound['standard_inchi_key']}")
```

### Workflow 3: ID Conversion

```python
from provesid import CheMBL

chembl = CheMBL()

# Convert ChEMBL ID to internal molregno
molregno = chembl.chembl_id_to_molregno('CHEMBL25')

# Convert back
chembl_id = chembl.molregno_to_chembl_id(molregno)
print(chembl_id)  # 'CHEMBL25'
```

## Database Schema

The ChEMBL class queries the following main tables:

- **molecule_dictionary**: Basic compound information (ChEMBL ID, name, max_phase)
- **compound_structures**: Chemical structures (SMILES, InChI, InChIKey). The
  `molfile` column is a quarter of the whole database and nothing here reads
  it; build a MOL block from the SMILES with RDKit when you need one.
- **compound_properties**: Physicochemical properties (MW, LogP, descriptors)
- **molecule_synonyms**: Alternative compound names
- **chembl_id_lookup**: ID mapping table

For detailed schema documentation, see: `src/provesid/data/schema_documentation.txt`

## References

- ChEMBL Database: https://www.ebi.ac.uk/chembl/
- ChEMBL Documentation: https://chembl.gitbook.io/chembl-interface-documentation/
- Database Schema: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/schema_documentation.html
