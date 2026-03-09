---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.1
kernelspec:
  display_name: physchem
  language: python
  name: python3
---

# OPSIN (Open Parser for Systematic IUPAC Nomenclature) Tutorial

OPSIN is a chemical name-to-structure service that converts IUPAC chemical names into chemical structures. This tutorial demonstrates how to use the `OPSIN` class from the `provesid` package to convert IUPAC names to SMILES, InChI, and InChI keys.

OPSIN was developed at the University of Cambridge and provides a reliable service for converting systematic chemical names to molecular representations.

```{code-cell} ipython3
from provesid import OPSIN
opsin = OPSIN()
print("OPSIN initialized successfully!")
print(f"Base URL: {opsin.base_url}")
```

## 1. Basic Name-to-Structure Conversion

The primary function of OPSIN is to convert IUPAC chemical names into molecular representations. Let's start with some simple examples:

```{code-cell} ipython3
# Convert a simple IUPAC name to structure
methane_result = opsin.get_id("methane")
print("Methane:")
print(f"  Status: {methane_result['status']}")
print(f"  SMILES: {methane_result['smiles']}")
print(f"  InChI: {methane_result['inchi']}")
print(f"  Standard InChI: {methane_result['stdinchi']}")
print(f"  InChI Key: {methane_result['stdinchikey']}")
```

```{code-cell} ipython3
methane_result
```

```{code-cell} ipython3
# Convert ethanol
ethanol_result = opsin.get_id("ethanol")
print("Ethanol:")
print(f"  Status: {ethanol_result['status']}")
print(f"  SMILES: {ethanol_result['smiles']}")
print(f"  InChI: {ethanol_result['inchi']}")
print(f"  InChI Key: {ethanol_result['stdinchikey']}")
```

## 2. Complex IUPAC Names

OPSIN excels at parsing complex systematic IUPAC names. Let's try some more challenging examples:

```{code-cell} ipython3
# Convert a more complex name
aspirin_name = "2-acetoxybenzoic acid"
aspirin_result = opsin.get_id(aspirin_name)
print(f"'{aspirin_name}':")
print(f"  Status: {aspirin_result['status']}")
print(f"  SMILES: {aspirin_result['smiles']}")
print(f"  InChI Key: {aspirin_result['stdinchikey']}")
```

```{code-cell} ipython3
# Convert a systematic name for caffeine
caffeine_name = "1,3,7-trimethylpurine-2,6-dione"
caffeine_result = opsin.get_id(caffeine_name)
print(f"'{caffeine_name}':")
print(f"  Status: {caffeine_result['status']}")
print(f"  SMILES: {caffeine_result['smiles']}")
print(f"  InChI Key: {caffeine_result['stdinchikey']}")
```

```{code-cell} ipython3
# Convert a complex organic molecule
complex_name = "2-phenylethyl acetate"
complex_result = opsin.get_id(complex_name)
print(f"'{complex_name}':")
print(f"  Status: {complex_result['status']}")
print(f"  SMILES: {complex_result['smiles']}")
print(f"  InChI: {complex_result['inchi']}")
```

## 3. Batch Processing

OPSIN provides a convenient method to process multiple chemical names at once:

```{code-cell} ipython3
# Process a list of IUPAC names
iupac_names = [
    "benzene",
    "toluene", 
    "phenol",
    "aniline",
    "benzoic acid"
]

print("Processing multiple IUPAC names:")
results = opsin.get_id_from_list(iupac_names)

for i, result in enumerate(results):
    name = iupac_names[i]
    print(f"\n{i+1}. {name}:")
    print(f"   Status: {result['status']}")
    if result['status'] == 'SUCCESS':
        print(f"   SMILES: {result['smiles']}")
        print(f"   InChI Key: {result['stdinchikey']}")
    else:
        print(f"   Error: Could not parse '{name}'")
```

## 4. Error Handling

OPSIN handles various types of input errors gracefully. Let's see how it responds to invalid or ambiguous names:

```{code-cell} ipython3
# Try an invalid chemical name
invalid_result = opsin.get_id("notarealchemicalname")
print("Invalid name 'notarealchemicalname':")
print(f"  Status: {invalid_result['status']}")
print(f"  SMILES: {invalid_result['smiles']}")

# Try an empty string
empty_result = opsin.get_id("")
print("\nEmpty string:")
print(f"  Status: {empty_result['status']}")

# Try a common name that might not be recognized
common_name_result = opsin.get_id("table salt")
print("\nCommon name 'table salt' (not IUPAC):")
print(f"  Status: {common_name_result['status']}")
```

## 5. Comparing Systematic vs Common Names

OPSIN works best with systematic IUPAC names. Let's compare results for systematic vs common names:

```{code-cell} ipython3
# Compare systematic vs common names
test_cases = [
    ("systematic", "propan-2-ol", "acetone (systematic)"),
    ("common", "isopropanol", "acetone (common name)"),
    ("systematic", "propanone", "acetone (systematic)"),
    ("common", "acetone", "acetone (common name)")
]

print("Comparing systematic vs common names:")
for name_type, name, description in test_cases:
    result = opsin.get_id(name)
    print(f"\n{description}:")
    print(f"  Name: '{name}' ({name_type})")
    print(f"  Status: {result['status']}")
    if result['status'] == 'SUCCESS':
        print(f"  SMILES: {result['smiles']}")
```

## 6. Practical Applications

Here are some practical use cases for the OPSIN class:

```{code-cell} ipython3
# Use case 1: Convert IUPAC names to SMILES for database storage
def name_to_smiles_converter(iupac_names):
    """Convert a list of IUPAC names to SMILES"""
    results = []
    for name in iupac_names:
        result = opsin.get_id(name)
        if result['status'] == 'SUCCESS':
            results.append({
                'name': name,
                'smiles': result['smiles'],
                'inchi_key': result['stdinchikey']
            })
        else:
            results.append({
                'name': name,
                'smiles': None,
                'inchi_key': None,
                'error': 'Failed to parse'
            })
    return results

# Test the converter
test_names = ["hexane", "cyclohexane", "benzene", "invalid_name"]
converted = name_to_smiles_converter(test_names)

print("IUPAC to SMILES conversion results:")
for item in converted:
    print(f"  {item['name']}: {item['smiles']}")
```

```{code-cell} ipython3
# Use case 2: Validate IUPAC names
def validate_iupac_names(names_list):
    """Check which names in a list are valid IUPAC names"""
    valid_names = []
    invalid_names = []
    
    for name in names_list:
        result = opsin.get_id(name)
        if result['status'] == 'SUCCESS':
            valid_names.append(name)
        else:
            invalid_names.append(name)
    
    return valid_names, invalid_names

# Test validation
test_names = [
    "methane",
    "ethanol", 
    "propanoic acid",
    "water",  # Common name, might not work
    "H2O",    # Formula, not IUPAC
    "benzene",
    "toluene"
]

valid, invalid = validate_iupac_names(test_names)
print("IUPAC Name Validation:")
print(f"✅ Valid names ({len(valid)}): {valid}")
print(f"❌ Invalid names ({len(invalid)}): {invalid}")
```

```{code-cell} ipython3
# Use case 3: Generate molecular identifiers for a compound database
def generate_molecular_identifiers(compound_name):
    """Generate complete molecular identifiers for a compound"""
    result = opsin.get_id(compound_name)
    
    if result['status'] == 'SUCCESS':
        return {
            'iupac_name': compound_name,
            'smiles': result['smiles'],
            'inchi': result['inchi'],
            'standard_inchi': result['stdinchi'],
            'inchi_key': result['stdinchikey'],
            'status': 'success'
        }
    else:
        return {
            'iupac_name': compound_name,
            'status': 'failed',
            'error': 'Could not parse IUPAC name'
        }

# Example: Create identifiers for pharmaceutical compounds
pharmaceutical_names = [
    "2-acetoxybenzoic acid",  # Aspirin
    "N-(4-hydroxyphenyl)acetamide",  # Paracetamol/Acetaminophen
    "2-phenylpropionic acid"  # Ibuprofen (simplified name)
]

print("Pharmaceutical compound identifiers:")
for name in pharmaceutical_names:
    identifiers = generate_molecular_identifiers(name)
    print(f"\n{name}:")
    if identifiers['status'] == 'success':
        print(f"  SMILES: {identifiers['smiles']}")
        print(f"  InChI Key: {identifiers['inchi_key']}")
    else:
        print(f"  Error: {identifiers['error']}")
```

## 7. Performance Considerations

When processing multiple compounds, OPSIN includes built-in rate limiting to be respectful to the service:

```{code-cell} ipython3
import time

# Demonstrate batch processing with timing
large_compound_list = [
    "methane", "ethane", "propane", "butane", "pentane",
    "hexane", "heptane", "octane", "nonane", "decane"
]

print("Processing 10 compounds with default pause (0.5s between requests):")
start_time = time.time()
results = opsin.get_id_from_list(large_compound_list)
end_time = time.time()

successful = sum(1 for r in results if r['status'] == 'SUCCESS')
print(f"\nResults:")
print(f"  Total compounds: {len(large_compound_list)}")
print(f"  Successful conversions: {successful}")
print(f"  Failed conversions: {len(results) - successful}")
print(f"  Total time: {end_time - start_time:.2f} seconds")
print(f"  Average time per compound: {(end_time - start_time)/len(large_compound_list):.2f} seconds")
```

## Summary

The `OPSIN` class provides two main methods:

1. **`get_id(iupac_name)`**: Convert a single IUPAC name to molecular identifiers
2. **`get_id_from_list(iupac_names)`**: Convert multiple IUPAC names with built-in rate limiting

### Key Features:
- ✅ **Systematic IUPAC name parsing**: Converts complex chemical names to structures
- ✅ **Multiple output formats**: SMILES, InChI, Standard InChI, InChI Key
- ✅ **Batch processing**: Handle multiple compounds efficiently
- ✅ **Error handling**: Graceful handling of invalid names
- ✅ **Rate limiting**: Built-in delays for batch processing
- ✅ **Free service**: No API key required

### Returned Data:
- **SMILES**: Simplified molecular-input line-entry system
- **InChI**: International Chemical Identifier
- **Standard InChI**: Standardized version of InChI
- **InChI Key**: Fixed-length hash of InChI
- **Status**: SUCCESS/FAILURE indication

### Best Practices:
1. Use systematic IUPAC names for best results
2. Handle failures gracefully (some names may not be recognized)
3. Use batch processing for multiple compounds
4. Respect rate limits when making many requests
5. Validate results before using in downstream applications

OPSIN is particularly valuable for:
- Converting literature chemical names to machine-readable formats
- Validating IUPAC nomenclature
- Building chemical databases from text sources
- Educational applications for learning chemical nomenclature

The service is provided by the University of Cambridge and is free to use for academic and commercial applications.
