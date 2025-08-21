# PubChem REST API Python Interface

This document describes the Python interface for the PubChem REST API (PUG-REST) implemented in `src/provesid/pubchem.py`.

## Overview

The `PubChemAPI` class provides a comprehensive interface to interact with PubChem's REST API for retrieving chemical compound, substance, and assay information. It includes support for:

- Compound and substance retrieval by various identifiers
- Property tables and synonyms
- Structure-based searches (substructure, similarity, identity)
- Assay information
- Error handling and rate limiting
- Convenience methods for common use cases

## Quick Start

```python
from provesid.pubchem import PubChemAPI, CompoundProperties

# Initialize the API
api = PubChemAPI()

# Get compound by CID
compound = api.get_compound_by_cid(2244)  # Aspirin

# Get compound properties
properties = [
    CompoundProperties.MOLECULAR_FORMULA,
    CompoundProperties.MOLECULAR_WEIGHT,
    CompoundProperties.SMILES,
    CompoundProperties.INCHIKEY
]
props = api.get_compound_properties(2244, properties)

# Search by name
cids = api.get_cids_by_name('aspirin')

# Structure similarity search
similar = api.similarity_search('CCO', threshold=90)  # Ethanol similarity
```

## API Classes and Constants

### Domain Classes
- `Domain`: API domains (COMPOUND, SUBSTANCE, ASSAY, etc.)
- `CompoundDomainNamespace`: Compound namespaces (CID, NAME, SMILES, etc.)
- `SubstanceDomainNamespace`: Substance namespaces (SID, SOURCEID, etc.)
- `AssayDomainNamespace`: Assay namespaces (AID, TYPE, etc.)

### Operation and Output Constants
- `Operation`: Available operations (RECORD, PROPERTY, SYNONYMS, etc.)
- `OutputFormat`: Output formats (JSON, XML, SDF, CSV, etc.)
- `CompoundProperties`: Available compound properties for property tables

### Structure Search Constants
- `StructureSearch`: Structure search types (SUBSTRUCTURE, SIMILARITY, etc.)
- `FastSearch`: Fast search methods (FASTIDENTITY, FASTSIMILARITY_2D, etc.)

## Main API Methods

### Compound Methods

#### `get_compound_by_cid(cid, output_format=OutputFormat.JSON)`
Get complete compound record by CID.

**Parameters:**
- `cid`: Compound ID (int or str)
- `output_format`: Desired output format

**Returns:** Compound data

#### `get_compounds_by_name(name, output_format=OutputFormat.JSON, name_type="word")`
Get compounds by name search.

**Parameters:**
- `name`: Compound name
- `output_format`: Desired output format
- `name_type`: Search type ("word" or "complete")

#### `get_compound_properties(cids, properties, output_format=OutputFormat.JSON)`
Get compound properties by CID(s).

**Parameters:**
- `cids`: Single CID or list of CIDs
- `properties`: List of property names (use `CompoundProperties` constants)
- `output_format`: Desired output format

**Example:**
```python
props = api.get_compound_properties(
    [2244, 5793], 
    [CompoundProperties.MOLECULAR_FORMULA, CompoundProperties.MOLECULAR_WEIGHT]
)
```

#### `get_compound_synonyms(cid, output_format=OutputFormat.JSON)`
Get compound synonyms by CID.

#### `get_cids_by_name(name, output_format=OutputFormat.JSON, name_type="word")`
Get CIDs by compound name.

#### `get_cids_by_smiles(smiles, output_format=OutputFormat.JSON)`
Get CIDs by SMILES string.

#### `get_cids_by_formula(formula, output_format=OutputFormat.JSON, allow_other_elements=False)`
Get CIDs by molecular formula using fast search.

### Structure Search Methods

#### `substructure_search(query, query_type="smiles", output_format=OutputFormat.JSON, **options)`
Perform substructure search.

**Parameters:**
- `query`: Query structure (SMILES, CID, etc.)
- `query_type`: Type of query ("smiles", "cid", etc.)
- `output_format`: Desired output format
- `**options`: Search options (MatchIsotopes, MaxRecords, etc.)

#### `similarity_search(query, query_type="smiles", threshold=90, output_format=OutputFormat.JSON, **options)`
Perform 2D similarity search.

**Parameters:**
- `threshold`: Similarity threshold (0-100)

#### `identity_search(query, query_type="smiles", identity_type="same_stereo_isotope", output_format=OutputFormat.JSON, **options)`
Perform identity search.

**Parameters:**
- `identity_type`: Type of identity match

### Substance Methods

#### `get_substance_by_sid(sid, output_format=OutputFormat.JSON)`
Get substance by SID.

#### `get_substances_by_name(name, output_format=OutputFormat.JSON)`
Get substances by name.

### Assay Methods

#### `get_assay_by_aid(aid, output_format=OutputFormat.JSON)`
Get assay by AID.

#### `get_assay_summary(cids, output_format=OutputFormat.JSON)`
Get assay summary for compounds.

### Convenience Methods

#### `search_compound(query, search_type="name")`
Search for compound with automatic format detection.

**Parameters:**
- `query`: Search query
- `search_type`: Type of search ("name", "smiles", "inchikey", "cid")

**Returns:** Dictionary with search results and metadata

#### `get_basic_compound_info(cid)`
Get basic compound information including properties and synonyms.

## Error Handling

The API includes custom exception classes:

- `PubChemError`: Base exception for PubChem API errors
- `PubChemTimeoutError`: Request timeout errors
- `PubChemNotFoundError`: Resource not found (404)
- `PubChemServerError`: Server errors (5xx)

Example:
```python
try:
    result = api.get_compound_by_cid(invalid_cid)
except PubChemNotFoundError:
    print("Compound not found")
except PubChemServerError:
    print("Server error - try again later")
```

## Rate Limiting

The API automatically enforces rate limiting with a default pause of 0.2 seconds between requests. This can be configured:

```python
api = PubChemAPI(pause_time=0.5)  # 0.5 seconds between requests
```

## Available Compound Properties

The `CompoundProperties` class defines constants for all available properties:

- `MOLECULAR_FORMULA`: Molecular formula
- `MOLECULAR_WEIGHT`: Molecular weight
- `SMILES`: SMILES string (connectivity)
- `INCHI`: InChI string
- `INCHIKEY`: InChI Key
- `IUPAC_NAME`: IUPAC name
- `EXACT_MASS`: Exact mass
- `TPSA`: Topological polar surface area
- `COMPLEXITY`: Molecular complexity
- `CHARGE`: Net charge
- `HBOND_DONOR_COUNT`: Hydrogen bond donors
- `HBOND_ACCEPTOR_COUNT`: Hydrogen bond acceptors
- `ROTATABLE_BOND_COUNT`: Rotatable bonds
- And many more...

## URL Building and Request Handling

The API automatically:
- Builds proper URLs according to PubChem REST API specification
- Handles URL encoding for special characters
- Manages HTTP headers and content types
- Parses responses according to requested format
- Implements retry logic and error handling

## Integration with PROVESID

The PubChem API is integrated into the PROVESID package and can be imported alongside other modules:

```python
from provesid import PubChemAPI, CompoundProperties
from provesid.utils import pubchem_cas_to_mol  # Existing utility functions
```

## Performance Considerations

- Use batch operations when possible (e.g., `get_compound_properties` with multiple CIDs)
- Consider output format: JSON for programmatic use, SDF for chemical data, CSV for tabular data
- Implement appropriate error handling for network issues
- Be mindful of PubChem's usage policies and rate limits

## Examples

### Basic Property Retrieval
```python
api = PubChemAPI()

# Get multiple properties for aspirin
properties = [
    CompoundProperties.MOLECULAR_FORMULA,
    CompoundProperties.MOLECULAR_WEIGHT,
    CompoundProperties.SMILES,
    CompoundProperties.IUPAC_NAME
]

result = api.get_compound_properties(2244, properties)
props = result['PropertyTable']['Properties'][0]

print(f"Formula: {props['MolecularFormula']}")
print(f"Weight: {props['MolecularWeight']}")
print(f"SMILES: {props['ConnectivitySMILES']}")
```

### Structure Search
```python
# Find compounds similar to caffeine
caffeine_smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
results = api.similarity_search(caffeine_smiles, threshold=85, MaxRecords=50)

if 'IdentifierList' in results:
    similar_cids = results['IdentifierList']['CID']
    print(f"Found {len(similar_cids)} similar compounds")
```

### Batch Processing
```python
# Get properties for multiple compounds
compound_list = [2244, 5793, 702, 5461]  # aspirin, caffeine, ethanol, ibuprofen
properties = [CompoundProperties.MOLECULAR_FORMULA, CompoundProperties.MOLECULAR_WEIGHT]

result = api.get_compound_properties(compound_list, properties)
for prop in result['PropertyTable']['Properties']:
    print(f"CID {prop['CID']}: {prop['MolecularFormula']} - {prop['MolecularWeight']}")
```

This implementation provides a comprehensive and user-friendly interface to the PubChem REST API, following the patterns established in the PROVESID project while adding robust error handling and extensive functionality.
