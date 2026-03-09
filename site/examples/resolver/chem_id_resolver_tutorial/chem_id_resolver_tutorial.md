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

# NCI Chemical Identifier Resolver Tutorial

The NCI Chemical Identifier Resolver is a powerful web service provided by the National Cancer Institute (NCI) that can convert between different types of chemical structure identifiers. This tutorial demonstrates how to use the `NCIChemicalIdentifierResolver` class and convenience functions from the `provesid` package to access this service.

The NCI resolver can handle various types of chemical identifiers including:
- Chemical names (common and IUPAC)
- CAS Registry Numbers
- SMILES notation
- InChI and InChIKey
- Chemical structure files (SDF)
- Chemical structure images

The service supports conversion between these formats and can provide additional molecular properties such as molecular weight, formula, and various structural descriptors.

**Base URL**: https://cactus.nci.nih.gov/chemical/structure/

**Service Pattern**: `{base_url}/{identifier}/{representation}`

```{code-cell} ipython3
from provesid import (
    NCIChemicalIdentifierResolver, 
    NCIResolverError, 
    NCIResolverNotFoundError,
    nci_cas_to_mol,
    nci_id_to_mol,
    nci_resolver,
    nci_smiles_to_names,
    nci_name_to_smiles,
    nci_inchi_to_smiles,
    nci_cas_to_inchi,
    nci_get_molecular_weight,
    nci_get_formula
)

# Initialize the NCI resolver
resolver = NCIChemicalIdentifierResolver()
print("NCI Chemical Identifier Resolver initialized successfully!")
print(f"Base URL: {resolver.base_url}")
print(f"Timeout: {resolver.timeout} seconds")
print(f"Rate limiting: {resolver.pause_time} seconds between requests")
```

## Available Representations

The NCI resolver supports many different chemical representations. Let's explore what's available:

```{code-cell} ipython3
# Display available representations
print("Available chemical representations:")
print("=" * 50)

for key, description in resolver.representations.items():
    print(f"  {key:<25} : {description}")

print(f"\nTotal representations available: {len(resolver.representations)}")
```

## 1. Basic Usage - Converting Between Identifiers

The primary method for converting chemical identifiers is `resolve()`. Let's start with some basic examples:

```{code-cell} ipython3
# Convert chemical names to SMILES
print("Converting chemical names to SMILES:")
compounds = ["aspirin", "caffeine", "water", "ethanol"]

for compound in compounds:
    try:
        smiles = resolver.resolve(compound, 'smiles')
        print(f"  {compound:<10} → {smiles}")
    except NCIResolverError as e:
        print(f"  {compound:<10} → Error: {e}")

print("\n" + "="*50)

# Convert SMILES to InChI
print("Converting SMILES to InChI:")
smiles_examples = ["CCO", "CC(=O)OC1=CC=CC=C1C(=O)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"]
names = ["ethanol", "aspirin", "caffeine"]

for smiles, name in zip(smiles_examples, names):
    try:
        inchi = resolver.resolve(smiles, 'stdinchi')
        print(f"  {name} ({smiles}):")
        print(f"    InChI: {inchi[:60]}..." if len(inchi) > 60 else f"    InChI: {inchi}")
    except NCIResolverError as e:
        print(f"  {name} → Error: {e}")
    print()
```

```{code-cell} ipython3
# Working with CAS Registry Numbers
print("Converting CAS Registry Numbers:")
cas_numbers = [
    ("50-78-2", "aspirin"),
    ("58-08-2", "caffeine"), 
    ("64-17-5", "ethanol"),
    ("7732-18-5", "water")
]

for cas, expected_name in cas_numbers:
    try:
        # Get IUPAC name
        iupac_name = resolver.resolve(cas, 'iupac_name')
        # Get SMILES
        smiles = resolver.resolve(cas, 'smiles')
        print(f"  CAS {cas} ({expected_name}):")
        print(f"    IUPAC Name: {iupac_name}")
        print(f"    SMILES: {smiles}")
    except NCIResolverError as e:
        print(f"  CAS {cas} → Error: {e}")
    print()
```

## 2. Getting Comprehensive Molecular Data

The `get_molecular_data()` method retrieves multiple properties and identifiers for a compound in a single call:

```{code-cell} ipython3
# Get comprehensive data for caffeine
caffeine_data = resolver.get_molecular_data("caffeine")

print("Comprehensive molecular data for caffeine:")
print("=" * 50)
print(f"Found by: {caffeine_data['found_by']}")
print(f"Success: {caffeine_data['success']}")
print(f"Note: {caffeine_data['note']}")
print()

# Display basic identifiers
print("Basic Identifiers:")
print(f"  SMILES: {caffeine_data.get('smiles')}")
print(f"  InChI: {caffeine_data.get('stdinchi')}")
print(f"  InChI Key: {caffeine_data.get('stdinchikey')}")
print(f"  CAS Number: {caffeine_data.get('cas')}")
print(f"  IUPAC Name: {caffeine_data.get('iupac_name')}")
print()

# Display molecular properties
print("Molecular Properties:")
print(f"  Formula: {caffeine_data.get('formula')}")
print(f"  Molecular Weight: {caffeine_data.get('mw')}")
print()

# Display NCI identifiers
print("NCI/CADD Identifiers:")
print(f"  FICTS: {caffeine_data.get('ficts')}")
print(f"  FICuS: {caffeine_data.get('ficus')}")
print(f"  uuuuu: {caffeine_data.get('uuuuu')}")
print(f"  HASHISY: {caffeine_data.get('hashisy')}")
print()

# Display names
names = caffeine_data.get('names', [])
if names:
    print(f"Chemical Names ({len(names)} found):")
    for i, name in enumerate(names[:5], 1):  # Show first 5 names
        print(f"  {i}. {name}")
    if len(names) > 5:
        print(f"  ... and {len(names) - 5} more names")
```

## 3. Using Convenience Functions

The package provides convenient functions for common operations without creating a resolver instance:

```{code-cell} ipython3
# Using convenience functions for quick conversions
print("Convenience function examples:")
print("=" * 40)

# Name to SMILES
aspirin_smiles = nci_name_to_smiles("aspirin")
print(f"aspirin → SMILES: {aspirin_smiles}")

# SMILES to names
if aspirin_smiles:
    names = nci_smiles_to_names(aspirin_smiles)
    print(f"SMILES → Names: {names[:3]}...")  # Show first 3 names

# CAS to InChI
cas_aspirin = "50-78-2"
inchi = nci_cas_to_inchi(cas_aspirin)
print(f"CAS {cas_aspirin} → InChI: {inchi[:50]}..." if inchi and len(inchi) > 50 else f"CAS {cas_aspirin} → InChI: {inchi}")

# Get molecular weight and formula
mw = nci_get_molecular_weight("caffeine")
formula = nci_get_formula("caffeine")
print(f"caffeine → MW: {mw}, Formula: {formula}")

print()
print("Using the legacy nci_cas_to_mol function:")
aspirin_mol = nci_cas_to_mol("50-78-2")
if aspirin_mol['success']:
    print(f"  Success: {aspirin_mol['success']}")
    print(f"  SMILES: {aspirin_mol.get('smiles')}")
    print(f"  Formula: {aspirin_mol.get('formula')}")
    print(f"  MW: {aspirin_mol.get('mw')}")
else:
    print(f"  Error: {aspirin_mol.get('error')}")

print()
print("Using the general nci_id_to_mol function:")
ethanol_mol = nci_id_to_mol("ethanol")
if ethanol_mol['success']:
    print(f"  Compound: ethanol")
    print(f"  SMILES: {ethanol_mol.get('smiles')}")
    print(f"  InChI Key: {ethanol_mol.get('stdinchikey')}")
    print(f"  Formula: {ethanol_mol.get('formula')}")
    print(f"  MW: {ethanol_mol.get('mw')}")
```

## 4. Batch Processing

For processing multiple compounds, use batch methods with built-in rate limiting:

```{code-cell} ipython3
# Batch resolve multiple compounds to SMILES
compounds = ["aspirin", "caffeine", "ibuprofen", "acetaminophen", "water"]

print("Batch conversion of compound names to SMILES:")
print("=" * 50)

smiles_results = resolver.batch_resolve(compounds, 'smiles')
for compound, smiles in smiles_results.items():
    status = "✓" if smiles else "✗"
    print(f"  {status} {compound:<15} → {smiles if smiles else 'Not found'}")

print()
print("Batch conversion to molecular weights:")
mw_results = resolver.batch_resolve(compounds, 'mw')
for compound, mw in mw_results.items():
    status = "✓" if mw else "✗"
    mw_value = f"{float(mw):.2f}" if mw and mw.replace('.', '').isdigit() else mw
    print(f"  {status} {compound:<15} → {mw_value if mw else 'Not found'}")

print()
print("Getting multiple representations for a single compound:")
representations = ['smiles', 'stdinchi', 'formula', 'mw', 'cas']
multi_results = resolver.resolve_multiple("aspirin", representations)

print("Aspirin in multiple formats:")
for rep, value in multi_results.items():
    status = "✓" if value else "✗"
    display_value = value if value else "Not available"
    if rep == 'stdinchi' and value and len(value) > 50:
        display_value = value[:50] + "..."
    print(f"  {status} {rep:<12} → {display_value}")
```

## 5. Error Handling

The NCI resolver provides robust error handling for various scenarios:

```{code-cell} ipython3
# Test error handling with invalid inputs
print("Testing error handling:")
print("=" * 30)

# Test invalid compound name
print("1. Invalid compound name:")
try:
    result = resolver.resolve("nonexistentcompound12345", 'smiles')
    print(f"   Result: {result}")
except NCIResolverNotFoundError as e:
    print(f"   ✓ Caught NCIResolverNotFoundError: {e}")
except NCIResolverError as e:
    print(f"   ✓ Caught NCIResolverError: {e}")

print()

# Test invalid representation
print("2. Invalid representation:")
try:
    result = resolver.resolve("aspirin", 'invalid_representation')
    print(f"   Result: {result}")
except ValueError as e:
    print(f"   ✓ Caught ValueError: {e}")

print()

# Test identifier validation
print("3. Identifier validation:")
valid_compounds = ["aspirin", "caffeine", "nonexistentcompound"]
for compound in valid_compounds:
    is_valid = resolver.is_valid_identifier(compound)
    status = "✓" if is_valid else "✗"
    print(f"   {status} {compound}: {'Valid' if is_valid else 'Invalid'}")

print()

# Test with timeout (using a very short timeout to demonstrate)
print("4. Timeout handling:")
try:
    timeout_resolver = NCIChemicalIdentifierResolver(timeout=0.001)  # Very short timeout
    result = timeout_resolver.resolve("aspirin", 'smiles')
    print(f"   Unexpected success: {result}")
except Exception as e:
    print(f"   ✓ Caught timeout/error: {type(e).__name__}: {e}")

print()
print("5. Safe function calls (return None on error):")
safe_results = [
    ("Valid compound", nci_name_to_smiles("aspirin")),
    ("Invalid compound", nci_name_to_smiles("nonexistentcompound12345")),
    ("Valid MW", nci_get_molecular_weight("caffeine")),
    ("Invalid MW", nci_get_molecular_weight("nonexistentcompound"))
]

for description, result in safe_results:
    status = "✓" if result is not None else "✗"
    print(f"   {status} {description}: {result}")
```

## 6. Working with Chemical Structure Images

The NCI resolver can generate chemical structure images:

```{code-cell} ipython3
# Generate image URLs for chemical structures
compounds = ["aspirin", "caffeine", "ibuprofen"]

print("Chemical structure image URLs:")
print("=" * 40)

for compound in compounds:
    try:
        # Get standard image URL
        image_url = resolver.get_image_url(compound)
        print(f"{compound}:")
        print(f"  Standard image: {image_url}")
        
        # Get larger PNG image URL
        large_image_url = resolver.get_image_url(compound, image_format='png', width=400, height=400)
        print(f"  Large PNG image: {large_image_url}")
        print()
        
    except NCIResolverError as e:
        print(f"{compound}: Error generating image URL - {e}")

print("Note: You can copy these URLs into a web browser to view the chemical structures.")
print()

# Demonstrate downloading an image (commented out to avoid file creation)
print("Example of downloading a structure image:")
print("# To download an image file:")
print("# success = resolver.download_image('aspirin', 'aspirin_structure.gif')")
print("# if success:")
print("#     print('Image downloaded successfully!')")
print("# else:")
print("#     print('Failed to download image')")
```

## 7. Practical Applications

Here are some practical use cases for the NCI Chemical Identifier Resolver:

```{code-cell} ipython3
# Use case 1: Build a compound database
def build_compound_database(compound_list):
    """Build a comprehensive database of chemical compounds"""
    database = {}
    
    print(f"Building database for {len(compound_list)} compounds...")
    
    for compound in compound_list:
        print(f"  Processing: {compound}")
        mol_data = nci_id_to_mol(compound)
        
        if mol_data['success']:
            database[compound] = {
                'smiles': mol_data.get('smiles'),
                'inchi_key': mol_data.get('stdinchikey'),
                'formula': mol_data.get('formula'),
                'molecular_weight': mol_data.get('mw'),
                'cas_number': mol_data.get('cas'),
                'iupac_name': mol_data.get('iupac_name'),
                'names': mol_data.get('names', [])[:5],  # First 5 names
                'identifiers': {
                    'ficts': mol_data.get('ficts'),
                    'ficus': mol_data.get('ficus'),
                    'uuuuu': mol_data.get('uuuuu')
                }
            }
        else:
            database[compound] = {'error': mol_data.get('error')}
    
    return database

# Example: Create database for common pharmaceuticals
pharmaceuticals = ["aspirin", "ibuprofen", "acetaminophen", "caffeine"]
pharma_db = build_compound_database(pharmaceuticals)

print("\\nPharmaceutical Database:")
print("=" * 50)
for compound, data in pharma_db.items():
    if 'error' not in data:
        print(f"{compound.upper()}:")
        print(f"  Formula: {data['formula']}")
        print(f"  MW: {data['molecular_weight']}")
        print(f"  SMILES: {data['smiles']}")
        print(f"  CAS: {data['cas_number']}")
        print(f"  Names: {', '.join(data['names'][:3])}...")
    else:
        print(f"{compound.upper()}: {data['error']}")
    print()
```

```{code-cell} ipython3
# Use case 2: Identifier conversion and standardization
def standardize_identifiers(mixed_identifiers):
    """Convert mixed chemical identifiers to standardized format"""
    standardized = {}
    
    for identifier in mixed_identifiers:
        print(f"Standardizing: {identifier}")
        
        # Get comprehensive data
        mol_data = nci_id_to_mol(identifier)
        
        if mol_data['success']:
            # Create standardized entry
            standardized[identifier] = {
                'input_identifier': identifier,
                'canonical_smiles': mol_data.get('smiles'),
                'standard_inchi': mol_data.get('stdinchi'),
                'inchi_key': mol_data.get('stdinchikey'),
                'molecular_formula': mol_data.get('formula'),
                'exact_mass': mol_data.get('mw'),
                'preferred_name': mol_data.get('iupac_name'),
                'cas_registry': mol_data.get('cas'),
                'alternative_names': mol_data.get('names', [])
            }
        else:
            standardized[identifier] = {
                'input_identifier': identifier,
                'error': 'Could not resolve identifier'
            }
    
    return standardized

# Example with mixed identifier types
mixed_ids = [
    "50-78-2",           # CAS number for aspirin
    "aspirin",           # Common name
    "acetylsalicylic acid", # Chemical name
    "CCO",               # SMILES for ethanol
    "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"  # InChI for ethanol
]

print("\\nIdentifier Standardization Example:")
print("=" * 45)
standardized_data = standardize_identifiers(mixed_ids)

for original_id, data in standardized_data.items():
    print(f"\\nOriginal ID: {original_id}")
    if 'error' not in data:
        print(f"  Canonical SMILES: {data['canonical_smiles']}")
        print(f"  InChI Key: {data['inchi_key']}")
        print(f"  Formula: {data['molecular_formula']}")
        print(f"  Preferred Name: {data['preferred_name']}")
        print(f"  CAS Number: {data['cas_registry']}")
    else:
        print(f"  Error: {data['error']}")
```

```{code-cell} ipython3
# Use case 3: Molecular property analysis
def analyze_molecular_properties(compound_list):
    """Analyze and compare molecular properties of compounds"""
    properties = []
    
    for compound in compound_list:
        mol_data = nci_id_to_mol(compound)
        if mol_data['success']:
            mw = mol_data.get('mw')
            formula = mol_data.get('formula')
            
            # Calculate some basic descriptors from molecular weight
            if mw and isinstance(mw, (int, float)):
                properties.append({
                    'name': compound,
                    'formula': formula,
                    'molecular_weight': float(mw),
                    'smiles': mol_data.get('smiles'),
                    'heavy_atom_estimate': len([c for c in formula if c.isupper()]) if formula else 0
                })
    
    return properties

# Analyze a series of related compounds
nsaids = ["aspirin", "ibuprofen", "naproxen", "diclofenac"]
print("\\nMolecular Property Analysis - NSAIDs:")
print("=" * 50)

nsaid_properties = analyze_molecular_properties(nsaids)

# Sort by molecular weight
nsaid_properties.sort(key=lambda x: x['molecular_weight'])

print("Compounds sorted by molecular weight:")
for prop in nsaid_properties:
    print(f"  {prop['name']:<12} | MW: {prop['molecular_weight']:>6.1f} | Formula: {prop['formula']:<12} | SMILES: {prop['smiles']}")

print()
print("Summary statistics:")
if nsaid_properties:
    mw_values = [p['molecular_weight'] for p in nsaid_properties]
    print(f"  Average MW: {sum(mw_values)/len(mw_values):.1f}")
    print(f"  MW Range: {min(mw_values):.1f} - {max(mw_values):.1f}")
    print(f"  Total compounds analyzed: {len(nsaid_properties)}")
```

## Summary

The `NCIChemicalIdentifierResolver` class and convenience functions provide comprehensive access to the NCI Chemical Identifier Resolver service:

### Main NCIChemicalIdentifierResolver Class Methods:
1. **`resolve(identifier, representation)`**: Convert any identifier to any representation
2. **`get_molecular_data(identifier)`**: Get comprehensive molecular data
3. **`resolve_multiple(identifier, representations)`**: Get multiple representations for one compound
4. **`batch_resolve(identifiers, representation)`**: Process multiple compounds efficiently
5. **`get_image_url(identifier)`** / **`download_image(identifier, filename)`**: Chemical structure images
6. **`is_valid_identifier(identifier)`**: Test if an identifier can be resolved

### Convenience Functions:
- **`nci_cas_to_mol(cas_rn)`**: Legacy function for CAS number conversion
- **`nci_id_to_mol(identifier)`**: General identifier conversion
- **`nci_resolver(input_value, output_type)`**: Simple conversion function
- **`nci_name_to_smiles(name)`** / **`nci_smiles_to_names(smiles)`**: Name-SMILES conversion
- **`nci_inchi_to_smiles(inchi)`** / **`nci_cas_to_inchi(cas_rn)`**: Structure format conversion
- **`nci_get_molecular_weight(identifier)`** / **`nci_get_formula(identifier)`**: Property extraction

### Supported Input Identifiers:
- Chemical names (common, trade, systematic)
- CAS Registry Numbers
- SMILES notation
- InChI and InChIKey
- Various database identifiers

### Supported Output Representations:
- **Structure identifiers**: SMILES, InChI, InChIKey, FICTS, FICuS, uuuuu, HASHISY
- **Properties**: Molecular weight, formula, exact mass, charge
- **Names**: IUPAC names, chemical names list, CAS numbers
- **Files**: SDF format
- **Images**: GIF, PNG structure images
- **Descriptors**: H-bond counts, rotatable bonds, ring counts

### Key Features:
- ✅ **Free Service**: No API key required
- ✅ **Rate Limiting**: Built-in delays for respectful API usage  
- ✅ **Error Handling**: Comprehensive exception handling with custom error types
- ✅ **Batch Processing**: Efficient handling of multiple compounds
- ✅ **Flexible Input**: Accepts various identifier types
- ✅ **Multiple Formats**: Convert between many representation types
- ✅ **Structure Images**: Generate and download chemical structure images
- ✅ **Legacy Support**: Maintains compatibility with older function interfaces

### Best Use Cases:
- Chemical identifier conversion and standardization
- Building chemical compound databases
- Molecular property analysis
- Chemical structure visualization
- Data integration from multiple chemical sources
- Chemical informatics research

### Service Information:
- **Provider**: NCI CADD Group (National Cancer Institute)
- **Base URL**: https://cactus.nci.nih.gov/chemical/structure/
- **API Pattern**: `{base_url}/{identifier}/{representation}`
- **Rate Limiting**: Recommended 0.1-1 second between requests
- **Availability**: Generally high, but occasional service interruptions possible

The NCI Chemical Identifier Resolver is an excellent choice for chemical identifier conversion tasks, offering broad coverage of chemical space and reliable performance for most chemical informatics applications.
