# ChEBI

This directory contains examples for using PROVESID's ChEBI interfaces:

## Files

- **ChEBI_tutorial.ipynb**: Tutorial for the ChEBI REST API (`ChEBI` class)
  - Online access to ChEBI database
  - Search by name, ID, or formula
  - Retrieve complete entity information

- **chebi_sdf_tutorial.ipynb**: Tutorial for the ChEBI SDF parser (`ChebiSDF` class)
  - Offline access to ChEBI SDF file
  - Fast queries with pre-built index (~190,000 compounds)
  - Search by ID, name, CAS, InChI, InChIKey, formula, synonyms
  - Export to pandas DataFrame
  - Access to 80+ database cross-references

## ChEBI SDF Setup

The `ChebiSDF` class will automatically download the ChEBI SDF file on first use if not found locally.

### Automatic Download (Recommended)

```python
from provesid import ChebiSDF

# On first use, will automatically download ~250 MB gzip file
# and extract to ~868 MB SDF file (takes ~2-3 minutes)
chebi_sdf = ChebiSDF()  # auto_download=True by default
```

### Manual Download (Optional)

If you prefer to download manually:

1. Download from: https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/chebi.sdf.gz
2. Extract and place in: `src/provesid/data/chebi.sdf`

Or use the download method:

```python
chebi_sdf = ChebiSDF(auto_download=False)
chebi_sdf.download_sdf()  # Downloads and extracts automatically
```

The file contains ~190,807 compounds. On first use, an index will be built (~15 seconds) and cached for faster subsequent loads.

## Quick Start

### ChEBI API (Online)
```python
from provesid import ChEBI

chebi = ChEBI()
results = chebi.search_by_name("aspirin")
water = chebi.get_complete_entity(15377)  # CHEBI:15377
```

### ChEBI SDF (Offline)
```python
from provesid import ChebiSDF

chebi_sdf = ChebiSDF()
water = chebi_sdf.get_compound_by_id("CHEBI:15377")
results = chebi_sdf.search_by_cas("50-78-2")  # aspirin
df = chebi_sdf.export_to_dataframe(["CHEBI:15377", "CHEBI:16236"])
```