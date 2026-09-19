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
- `get_compound_properties()` - Selected properties for one compound, with synonyms
- `get_all_compound_info()` - All available properties
- `get_properties_for_cids()` - Raw property table for many CIDs, in bulk
- `get_compound_properties_batch()` - Bulk retrieval reshaped to one dict per CID

### Utility Methods
- `get_compound_synonyms()` - Get compound synonyms
- `get_compound_identifiers()` - Extract specific identifiers

## Batch Processing

PubChem's property endpoint answers a whole list of CIDs in one round trip, so
asking about a thousand compounds costs a handful of requests rather than a
thousand.

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

`get_compound_properties_batch()` returns one dict per CID, in the order asked,
each carrying the retrieved properties plus `success`, `cid` and `error`. A CID
PubChem has no record of comes back with `success=False` rather than being
dropped, so the result can be zipped against the input.

`get_properties_for_cids()` is the same request without the reshaping — it
returns PubChem's own `PropertyTable` rows, where a compound with no record
yields a row holding only `CID`:

```python
rows = pc.get_properties_for_cids(range(2000, 3000), ["MolecularWeight"])
# 1000 compounds, 5 requests (PROPERTY_CHUNK_SIZE = 200 CIDs each)
```

Both methods switch from a URL path to a POST body once the identifier list
outgrows `URL_IDENTIFIER_LIMIT`, which is how PubChem asks for long lists to be
sent. Neither includes synonyms: those need one request per compound, which
would defeat the batching. Use `get_compound_synonyms()` where they are needed.

## Properties Without the Network

`PubChemID` — the local SQLite database of ~1.6M compounds — carries the
identifiers and the cheap computed descriptors, so most property lookups need no
request at all. `properties()` reads it first and consults PUG-REST only for
what it cannot answer:

```python
from provesid import PubChemID

db = PubChemID()

db.properties(2244, ["MolecularFormula", "MolecularWeight"])
# {'CID': 2244, 'Source': 'offline', 'MolecularFormula': 'C9H8O4',
#  'MolecularWeight': 180.16}

db.properties(2244, ["MonoisotopicMass"])
# {'CID': 2244, 'Source': 'online', 'MonoisotopicMass': 180.04225873}
```

The `Source` key records which one answered. Two rules decide it:

- a CID with no row in the local database goes online;
- a property with no local column sends the *whole* request online, because that
  property would need a request anyway and a row assembled from two PubChem
  snapshots is worse than a row from one.

`PubChemID.OFFLINE_PROPERTIES` lists what can be served locally: the formula,
weight, exact mass, isomeric SMILES, InChI, InChIKey, IUPAC name, title, XLogP,
TPSA, complexity, charge and the H-bond, rotatable-bond and heavy-atom counts.
Everything else — `MonoisotopicMass`, `ConnectivitySMILES`, the 3D descriptors,
the patent and literature counts — is online-only.

Bulk lookups split themselves between the two sources automatically, and the
online remainder is fetched in one batched request:

```python
rows = db.properties_for_cids(my_ten_thousand_cids, ["MolecularWeight", "InChIKey"])
table = db.properties_table(my_ten_thousand_cids, ["MolecularWeight"])
```

`properties_table()` returns a DataFrame with a row for every CID asked about,
`Source` reading `offline`, `online` or `missing`, so it can be joined against
your own table. Pass `use_online_fallback=False` to keep a lookup strictly
local — a hard guarantee, not a preference.

A property the compound has no value for is **absent** from the result rather
than `None`, which is how PubChem itself reports it: `'XLogP' not in result`
means PubChem computes no logP for that compound, not that the lookup fell
short. Values are normalised to one type across both sources, since PUG-REST
reports `MolecularWeight` as a string where the local database holds a float.

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

For a comprehensive tutorial with examples, see: [PubChem Tutorial](../examples/pubchem/pubchem_tutorial.md)
