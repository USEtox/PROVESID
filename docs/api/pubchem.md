# PubChem API

The PubChem API module provides access to PubChem's REST services for retrieving compound and substance data. This module has been recently enhanced with improved data access patterns and new search methods.

::: provesid.pubchem

## Quick Start

```python
from provesid import PubChemAPI, Domain, CompoundProperties

# Initialize the API client
pc = PubChemAPI()

# Search for compounds by name
cids = pc.get_cids_by_name('aspirin')
print(f"Found CIDs: {cids}")

# Get basic compound information
basic_info = pc.get_basic_compound_info(cids[0])
print(f"Molecular Formula: {basic_info['MolecularFormula']}")
print(f"Molecular Weight: {basic_info['MolecularWeight']}")

# Get compound by CID (improved - no wrapper needed!)
compound = pc.get_compound_by_cid(cids[0])
print(f"Compound keys: {list(compound.keys())}")
```

## Key Improvements

### Elegant Data Access ✨

**Before (redundant wrapper access):**
```python
# Old way required nested access
substance = pc.get_substance_by_sid(sid)
data = substance["PC_Substances"][0]  # Redundant wrapper

compound = pc.get_compound_by_cid(cid)
data = compound["PC_Compounds"][0]    # Redundant wrapper
```

**Now (direct access):**
```python
# New way provides direct access
substance = pc.get_substance_by_sid(sid)  # Direct access!
compound = pc.get_compound_by_cid(cid)    # Direct access!
```

### Enhanced Search Methods

#### Multiple Search Domains
```python
# Search in compound domain (default)
cids = pc.get_cids_by_name('aspirin', domain=Domain.COMPOUND)

# Search in substance domain  
cids = pc.get_cids_by_name('8000-78-0', domain=Domain.SUBSTANCE)

# Comprehensive search across both domains
results = pc.find_cids_comprehensive('8000-78-0')
```

#### Structure-Based Searching
```python
# Search by SMILES (now returns clean list)
smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # aspirin
cids = pc.get_cids_by_smiles(smiles)

# Search by InChI Key (newly implemented)
inchi_key = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"  # aspirin
cids = pc.get_cids_by_inchikey(inchi_key)

# Get compound records directly
compound = pc.get_compounds_by_smiles(smiles)
compound = pc.get_compounds_by_inchikey(inchi_key)
```

## Available Methods

### Compound Search Methods
- `get_cids_by_name()` - Search by compound name
- `get_cids_by_smiles()` - Search by SMILES string  
- `get_cids_by_inchikey()` - Search by InChI Key ✨ *New*
- `find_cids_comprehensive()` - Multi-domain search

### Compound Data Methods
- `get_compound_by_cid()` - Get compound record ✨ *Improved*
- `get_compounds_by_name()` - Get compounds by name ✨ *Improved*
- `get_compounds_by_smiles()` - Get compounds by SMILES ✨ *Improved*
- `get_compounds_by_inchikey()` - Get compounds by InChI Key ✨ *Improved*

### Substance Methods
- `get_substance_by_sid()` - Get substance record ✨ *Improved*
- `get_substances_by_name()` - Get substances by name ✨ *Improved*
- `get_sids_by_name()` - Search for substance IDs

### Property Methods
- `get_basic_compound_info()` - Essential compound properties
- `get_compound_properties()` - Selected properties
- `get_all_compound_info()` - All available properties
- `get_compound_properties_batch()` - Batch processing

### Utility Methods
- `get_compound_synonyms()` - Get compound synonyms
- `get_compound_identifiers()` - Extract specific identifiers

## Batch Processing

```python
# Process multiple compounds efficiently
compound_names = ["aspirin", "caffeine", "acetaminophen", "ibuprofen"]
all_cids = []

for name in compound_names:
    cids = pc.get_cids_by_name(name)
    if cids:
        all_cids.append(cids[0])

# Batch property retrieval
properties = [CompoundProperties.MOLECULAR_WEIGHT, 
              CompoundProperties.MOLECULAR_FORMULA,
              CompoundProperties.SMILES]

batch_results = pc.get_compound_properties_batch(all_cids, properties)
```

## Error Handling

```python
from provesid import PubChemNotFoundError, PubChemError

try:
    cids = pc.get_cids_by_name('invalid_compound_name')
    if not cids:
        print("No compounds found")
except PubChemNotFoundError:
    print("Compound not found in PubChem")
except PubChemError as e:
    print(f"PubChem API error: {e}")
```

## Best Practices

- Use batch methods for multiple compounds
- Always handle potential errors with try/except blocks
- Use domain-specific searches when appropriate
- Check data availability before processing
- Respect PubChem's rate limits (built into the API)

## Tutorial

For a comprehensive tutorial with examples, see: [PubChem Tutorial](../examples/pubchem/pubchem_tutorial.ipynb)
