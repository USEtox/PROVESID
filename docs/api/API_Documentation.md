# PROVESID API Documentation

This document provides comprehensive documentation for the PROVESID Python package, which includes interfaces to the PubChem REST API and the NCI Chemical Identifier Resolver.

## Table of Contents

1. [Installation](#installation)
2. [PubChem API](#pubchem-api)
3. [PubChem PUG View](#pubchem-pug-view)
4. [NCI Chemical Identifier Resolver](#nci-chemical-identifier-resolver)
5. [Quick Start Examples](#quick-start-examples)
6. [Error Handling](#error-handling)
7. [Advanced Usage](#advanced-usage)

## Installation

```python
# Import the PROVESID package
from provesid import (
    PubChemAPI, 
    PubChemView,
    NCIChemicalIdentifierResolver,
    get_experimental_property,
    get_property_table,
    nci_cas_to_mol,
    nci_name_to_smiles,
    nci_get_molecular_weight
)
```

## PubChem API

The `PubChemAPI` class provides a comprehensive interface to the PubChem REST API.

### Basic Usage

```python
from provesid import PubChemAPI

# Initialize the API
pubchem = PubChemAPI()

# Get compound information by name
compound_data = pubchem.get_compound_by_name("aspirin")
print(f"CID: {compound_data['PC_Compounds'][0]['id']['id']['cid']}")

# Get properties
properties = pubchem.get_properties("aspirin", ["MolecularFormula", "MolecularWeight"])
print(f"Formula: {properties[0]['MolecularFormula']}")
print(f"MW: {properties[0]['MolecularWeight']}")

# Get synonyms
synonyms = pubchem.get_synonyms("aspirin")
print(f"Synonyms: {synonyms[:5]}")  # First 5 synonyms
```

### Available Methods

#### Compound Retrieval
- `get_compound_by_name(name, namespace='name', domain='compound')`
- `get_compound_by_cid(cid, record_type='2d', response_type='JSON')`
- `get_compound_by_smiles(smiles, record_type='2d', response_type='JSON')`
- `get_compound_by_inchi(inchi, record_type='2d', response_type='JSON')`

#### Properties and Data
- `get_properties(identifier, properties, namespace='name')`
- `get_synonyms(identifier, namespace='name')`
- `get_sdf(identifier, namespace='name')`
- `get_png_image(identifier, namespace='name', image_size='large')`

#### Searching
- `search_compounds(query, max_records=20)`
- `search_by_formula(formula, max_records=20)`
- `search_by_structure(smiles, search_type='similarity', threshold=90)`

#### Batch Operations
- `batch_get_properties(identifiers, properties, namespace='name')`

### PubChem Properties

Available properties for `get_properties()`:
- `MolecularFormula`
- `MolecularWeight`
- `CanonicalSMILES`
- `InChI`
- `InChIKey`
- `IUPACName`
- `XLogP`
- `TPSA`
- `Complexity`
- `Charge`
- `HBondDonorCount`
- `HBondAcceptorCount`
- `RotatableBondCount`

## PubChem PUG View

The `PubChemView` class provides access to experimental properties from PubChem that are not available through the standard PubChem API. This includes data like melting points, boiling points, density, solubility, and many other experimentally determined properties.

### Basic Usage

```python
from provesid import PubChemView

# Initialize the PUG View client
pugview = PubChemView()

# Get experimental properties for aspirin (CID 2244)
melting_points = pugview.get_melting_point(2244)
print(f"Melting point: {melting_points[0].value}")
print(f"Reference: {melting_points[0].reference}")

# Get all available experimental properties
all_properties = pugview.extract_all_experimental_properties(2244)
print(f"Found {len(all_properties)} property types")

# Get property summary
summary = pugview.get_property_summary(2244, "Boiling Point")
print(f"Boiling point values: {summary['values']}")

# Get comprehensive property table
table = pugview.get_property_table(2244, "Melting Point")
print(f"Table columns: {table.columns.tolist()}")
print(f"First entry: {table.iloc[0]['ExperimentalValue']} {table.iloc[0]['Unit']}")
```

### Available Methods

#### Property Extraction
- `extract_property_data(cid, property_name)` - Extract structured data for a specific property
- `extract_all_experimental_properties(cid)` - Get all experimental properties for a compound
- `get_property_summary(cid, property_name)` - Get summary statistics for a property
- `get_available_properties(cid)` - List all available properties for a compound
- `get_property_table(cid, property_name)` - Get comprehensive table with parsed values and full references

#### Convenience Methods for Common Properties
- `get_melting_point(cid)` - Melting point data
- `get_boiling_point(cid)` - Boiling point data
- `get_density(cid)` - Density data
- `get_solubility(cid)` - Solubility data
- `get_flash_point(cid)` - Flash point data
- `get_vapor_pressure(cid)` - Vapor pressure data
- `get_viscosity(cid)` - Viscosity data
- `get_logp(cid)` - LogP data
- `get_refractive_index(cid)` - Refractive index data

#### Batch Operations
- `batch_extract_properties(cid, property_names)` - Extract multiple properties at once

### PropertyData Object

Each experimental property is returned as a `PropertyData` object with the following attributes:
- `value` - The property value as a string
- `unit` - Extracted unit (if identifiable)
- `conditions` - Experimental conditions (e.g., temperature)
- `reference` - Literature reference
- `reference_number` - Reference number in PubChem
- `description` - Additional description
- `name` - Property name or identifier

### Supported Experimental Properties

The PUG View interface supports all experimental properties available in PubChem, including:

- **Physical Properties**: Melting Point, Boiling Point, Density, Refractive Index
- **Thermal Properties**: Flash Point, Autoignition Temperature, Heat of Combustion
- **Chemical Properties**: pH, Dissociation Constants, LogP, LogS
- **Transport Properties**: Viscosity, Surface Tension, Vapor Pressure
- **Sensory Properties**: Odor, Taste, Color/Form, Physical Description
- **Safety Properties**: Corrosivity, Stability, Decomposition
- **Spectroscopic**: Collision Cross Section, Kovats Retention Index
- **Biological**: Caco2 Permeability, Hydrophobicity

### Convenience Functions

For quick access without creating a PubChemView instance:

```python
from provesid import get_experimental_property, get_property_values_only, get_property_table

# Get full property data
density_data = get_experimental_property(2244, "Density")
print(f"Density: {density_data[0].value}")

# Get just the values as strings
bp_values = get_property_values_only(702, "Boiling Point")
print(f"Ethanol boiling points: {bp_values}")

# Get comprehensive property table
table = get_property_table(2244, "Melting Point")
print(f"Found {len(table)} experimental entries")
print(table[['ExperimentalValue', 'Unit', 'FullReference']].head())
```

### Property Table Format

The `get_property_table()` function returns a pandas DataFrame with the following columns:

- **CID**: PubChem Compound ID
- **StringWithMarkup**: Original experimental text from PubChem
- **ExperimentalValue**: Parsed numerical value
- **Unit**: Extracted unit of measurement
- **FullReference**: Complete citation with source, description, and URL

## NCI Chemical Identifier Resolver

The `NCIChemicalIdentifierResolver` class provides access to the NCI Chemical Identifier Resolver service.

### Basic Usage

```python
from provesid import NCIChemicalIdentifierResolver

# Initialize the resolver
resolver = NCIChemicalIdentifierResolver()

# Convert name to SMILES
smiles = resolver.resolve("aspirin", "smiles")
print(f"Aspirin SMILES: {smiles}")

# Get molecular weight
mw = resolver.resolve("caffeine", "mw")
print(f"Caffeine MW: {mw}")

# Get comprehensive molecular data
mol_data = resolver.get_molecular_data("ethanol")
print(f"Ethanol data: {mol_data}")
```

### Available Methods

#### Basic Resolution
- `resolve(identifier, representation, output_format='text')`
- `resolve_multiple(identifier, representations)`
- `get_molecular_data(identifier)`

#### Batch Operations
- `batch_resolve(identifiers, representation, output_format='text')`

#### Utility Methods
- `get_image_url(identifier, format='gif')`
- `get_available_representations()`

### Supported Representations

The resolver supports many chemical representations:
- `smiles` - Canonical SMILES
- `stdinchi` - Standard InChI
- `stdinchikey` - Standard InChI Key
- `formula` - Molecular formula
- `mw` - Molecular weight
- `names` - Chemical names and synonyms
- `iupac_name` - IUPAC name
- `cas` - CAS Registry Number
- `sdf` - Structure Data File
- `mol` - MDL Molfile
- `image` - Structure image

### Convenience Functions

For common operations, use the convenience functions:

```python
from provesid import (
    nci_cas_to_mol,
    nci_name_to_smiles,
    nci_get_molecular_weight,
    nci_get_formula
)

# Convert CAS number to molecular data
mol_data = nci_cas_to_mol("64-17-5")  # ethanol
print(f"SMILES: {mol_data['smiles']}")

# Convert name to SMILES
smiles = nci_name_to_smiles("glucose")
print(f"Glucose SMILES: {smiles}")

# Get molecular weight
mw = nci_get_molecular_weight("water")
print(f"Water MW: {mw}")

# Get formula
formula = nci_get_formula("benzene")
print(f"Benzene formula: {formula}")
```

## Quick Start Examples

### Example 1: Compare Data from All Three APIs

```python
from provesid import PubChemAPI, PubChemView, NCIChemicalIdentifierResolver

pubchem = PubChemAPI()
pugview = PubChemView()
nci = NCIChemicalIdentifierResolver()

compound = "aspirin"
cid = 2244  # aspirin CID

# Get basic data from PubChem
pc_props = pubchem.get_properties(compound, ["MolecularFormula", "MolecularWeight", "CanonicalSMILES"])
print(f"PubChem - Formula: {pc_props[0]['MolecularFormula']}, MW: {pc_props[0]['MolecularWeight']}")

# Get experimental properties from PUG View
melting_point = pugview.get_melting_point(cid)
boiling_point = pugview.get_boiling_point(cid)
print(f"PUG View - MP: {melting_point[0].value}, BP: {boiling_point[0].value}")

# Get data from NCI
nci_data = nci.get_molecular_data(compound)
print(f"NCI - Formula: {nci_data['formula']}, MW: {nci_data['mw']}")
```

### Example 2: Comprehensive Property Analysis

```python
from provesid import PubChemView

pugview = PubChemView()

# Get all experimental properties for ethanol
all_props = pugview.extract_all_experimental_properties(702)  # ethanol CID

# Focus on physical properties
physical_props = ["Melting Point", "Boiling Point", "Density", "Viscosity", "Refractive Index"]
for prop in physical_props:
    if prop in all_props:
        data_list = all_props[prop]
        print(f"{prop}: {len(data_list)} values")
        print(f"  First value: {data_list[0].value}")
        if data_list[0].unit:
            print(f"  Unit: {data_list[0].unit}")
```

### Example 3: Batch Processing with PUG View

```python
from provesid import PubChemView

pugview = PubChemView()

# Analyze multiple compounds
compounds = {"aspirin": 2244, "caffeine": 2519, "ibuprofen": 3672}
properties = ["Melting Point", "Boiling Point", "Solubility"]

for name, cid in compounds.items():
    print(f"\n{name.title()}:")
    batch_results = pugview.batch_extract_properties(cid, properties)
    
    for prop_name, prop_data in batch_results.items():
        if prop_data:
            print(f"  {prop_name}: {prop_data[0].value}")
        else:
            print(f"  {prop_name}: No data available")
```

## Error Handling

All three APIs provide specific exception classes for error handling:

### PubChem Errors

```python
from provesid import PubChemAPI, PubChemError, PubChemNotFoundError

pubchem = PubChemAPI()

try:
    data = pubchem.get_compound_by_name("nonexistent_compound_12345")
except PubChemNotFoundError:
    print("Compound not found in PubChem")
except PubChemError as e:
    print(f"PubChem API error: {e}")
```

### PubChem PUG View Errors

```python
from provesid import PubChemView, PubChemViewError, PubChemViewNotFoundError

pugview = PubChemView()

try:
    data = pugview.extract_property_data(99999999, "Melting Point")
except PubChemViewNotFoundError:
    print("Property or compound not found in PUG View")
except PubChemViewError as e:
    print(f"PUG View error: {e}")
```

### NCI Resolver Errors

```python
from provesid import NCIChemicalIdentifierResolver, NCIResolverError, NCIResolverNotFoundError

resolver = NCIChemicalIdentifierResolver()

try:
    smiles = resolver.resolve("invalid_compound_xyz", "smiles")
except NCIResolverNotFoundError:
    print("Compound not found in NCI resolver")
except NCIResolverError as e:
    print(f"NCI resolver error: {e}")
```

## Advanced Usage

### Custom Request Configuration

```python
from provesid import PubChemAPI, PubChemView

# Initialize with custom settings
pubchem = PubChemAPI(
    base_url="https://pubchem.ncbi.nlm.nih.gov/rest/pug",
    timeout=30,
    max_retries=3,
    backoff_factor=1.0
)

pugview = PubChemView(
    base_url="https://pubchem.ncbi.nlm.nih.gov/rest/pug_view",
    timeout=30,
    max_retries=3
)
```

### Rate Limiting

All APIs implement automatic rate limiting to respect server limits:
- PubChem: Maximum 5 requests per second
- PubChem PUG View: Maximum 5 requests per second
- NCI Resolver: Maximum 3 requests per second

### Logging

Enable logging to monitor API calls:

```python
import logging
from provesid import PubChemAPI, PubChemView, NCIChemicalIdentifierResolver

# Enable debug logging
logging.basicConfig(level=logging.DEBUG)

pubchem = PubChemAPI()
pugview = PubChemView()
resolver = NCIChemicalIdentifierResolver()
```

## Notes and Limitations

1. **PubChem API**: 
   - Rate limited to 5 requests per second
   - Large batch requests may take time
   - Some properties may not be available for all compounds

2. **PubChem PUG View**:
   - Rate limited to 5 requests per second
   - Provides experimental data not available in standard API
   - Property values may have complex formatting requiring parsing
   - Not all compounds have experimental property data

3. **NCI Resolver**:
   - Rate limited to 3 requests per second  
   - Not all identifiers may be recognized
   - Service availability may vary

4. **General**:
   - Network connectivity required
   - Response times depend on server load
   - Some chemical names may have ambiguous mappings

For more examples and updates, see the test files and demo scripts in the repository.
