# OPSIN API

The OPSIN (Open Parser for Systematic IUPAC Nomenclature) class provides an interface to convert IUPAC chemical names to molecular structures and identifiers.

## Overview

OPSIN is a web service that converts systematic IUPAC chemical names into chemical structures. This class provides a Python interface to the OPSIN web API hosted at Cambridge University.

## Class: OPSIN

### Initialization

```python
from provesid.opsin import OPSIN

opsin = OPSIN()
```

The OPSIN class initializes with:
- `base_url`: "https://opsin.ch.cam.ac.uk/opsin/"
- `responses`: Status code mappings

### Methods

#### `get_id(name, timeout=30)`

Convert a single IUPAC name to chemical identifiers.

**Parameters:**
- `name` (str): IUPAC chemical name
- `timeout` (int, optional): Request timeout in seconds (default: 30)

**Returns:**
- `dict`: Dictionary containing:
  - `status`: "SUCCESS", "FAILURE", or "Internal server error"
  - `message`: Status message
  - `inchi`: InChI string
  - `stdinchi`: Standard InChI string
  - `stdinchikey`: Standard InChI Key
  - `smiles`: SMILES string

**Example:**
```python
opsin = OPSIN()
result = opsin.get_id('ethanol')

if result['status'] == 'SUCCESS':
    print(f"SMILES: {result['smiles']}")
    print(f"InChI: {result['stdinchi']}")
    print(f"InChI Key: {result['stdinchikey']}")
```

#### `get_id_from_list(name_list, pause_time=1, timeout=30)`

Convert multiple IUPAC names to chemical identifiers.

**Parameters:**
- `name_list` (list): List of IUPAC chemical names
- `pause_time` (float, optional): Delay between requests in seconds (default: 1)
- `timeout` (int, optional): Request timeout in seconds (default: 30)

**Returns:**
- `list`: List of dictionaries with same structure as `get_id()`

**Example:**
```python
opsin = OPSIN()
names = ['ethanol', 'methanol', 'benzene']
results = opsin.get_id_from_list(names, pause_time=0.5)

for i, result in enumerate(results):
    if result['status'] == 'SUCCESS':
        print(f"{names[i]}: {result['smiles']}")
```

### Supported Name Types

OPSIN supports various types of systematic chemical names:

- **Simple organic compounds**: ethanol, methanol, propane
- **Substituted compounds**: 2-methylpropane, 2,2-dimethylpropane
- **Cyclic compounds**: cyclohexane, benzene
- **Functional groups**: phenol, benzoic acid, acetic acid
- **Complex systematic names**: Following IUPAC nomenclature rules

### Error Handling

The OPSIN class handles various error conditions:

- **Network timeouts**: Configurable timeout parameter
- **Invalid names**: Returns status "FAILURE" for unrecognized names
- **API unavailability**: Returns appropriate error status
- **URL encoding**: Automatically handles special characters in names

### Best Practices

1. **Rate limiting**: Use `pause_time` parameter when processing multiple names
2. **Timeout handling**: Adjust timeout for slow network connections
3. **Validation**: Always check the `status` field before using results
4. **Batch processing**: Use `get_id_from_list()` for multiple names

### Common Use Cases

#### Single Name Conversion
```python
opsin = OPSIN()
result = opsin.get_id('2-methylbutane')

if result['status'] == 'SUCCESS':
    smiles = result['smiles']
    inchi_key = result['stdinchikey']
```

#### Batch Processing with Error Handling
```python
opsin = OPSIN()
names = ['ethanol', 'invalid_name', 'benzene']
results = opsin.get_id_from_list(names, pause_time=0.5)

successful_results = []
failed_names = []

for i, result in enumerate(results):
    if result['status'] == 'SUCCESS':
        successful_results.append({
            'name': names[i],
            'smiles': result['smiles'],
            'inchi_key': result['stdinchikey']
        })
    else:
        failed_names.append(names[i])

print(f"Successfully converted: {len(successful_results)}")
print(f"Failed conversions: {failed_names}")
```

#### Cross-validation with Other Services
```python
from provesid.opsin import OPSIN
from provesid.pubchem import PubChemAPI

opsin = OPSIN()
pubchem = PubChemAPI()

# Get structure from IUPAC name
name = "ethanol"
opsin_result = opsin.get_id(name)

if opsin_result['status'] == 'SUCCESS':
    smiles = opsin_result['smiles']
    
    # Validate with PubChem
    pubchem_result = pubchem.get_compound_by_smiles(smiles)
    
    if pubchem_result:
        print(f"OPSIN SMILES: {smiles}")
        print(f"PubChem validation: Success")
```

### Limitations

- **Network dependency**: Requires internet connection to Cambridge OPSIN server
- **Name recognition**: Limited to names that OPSIN can parse
- **Processing time**: Some complex names may take longer to process
- **Rate limits**: Cambridge server may have rate limiting policies

### Troubleshooting

#### Common Issues

1. **Timeout errors**: Increase timeout parameter or check network connection
2. **Name not recognized**: Verify IUPAC naming conventions
3. **Empty results**: Check if the name is supported by OPSIN
4. **Network errors**: Verify internet connectivity and server availability

#### Debug Information

```python
opsin = OPSIN()
result = opsin.get_id('problematic_name')

print(f"Status: {result['status']}")
print(f"Message: {result['message']}")

if result['status'] == 'FAILURE':
    print("Name not recognized by OPSIN")
elif result['status'] == 'Internal server error':
    print("OPSIN server error - try again later")
```
