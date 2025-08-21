# Quick Start Guide

This guide will help you get up and running with PROVESID quickly. We'll cover the most common use cases and provide practical examples.

## Basic Compound Lookup

### Get Compound Information from PubChem

```python
from provesid import PubChemAPI

# Initialize the API
api = PubChemAPI()

# Get compound by name
aspirin = api.get_compounds_by_name('aspirin')
print(f"Found aspirin compound: {aspirin}")

# Get compound by CID (Compound ID)
compound = api.get_compound_by_cid(2244)  # Aspirin CID
print(f"Aspirin structure: {compound}")

# Get multiple properties at once
properties = api.get_compound_properties(
    [2244], 
    ['MolecularWeight', 'MolecularFormula', 'ConnectivitySMILES']
)
print(f"Aspirin properties: {properties}")
```

### Chemical Identifier Conversion

```python
from provesid import NCIChemicalIdentifierResolver

# Initialize the resolver
resolver = NCIChemicalIdentifierResolver()

# Convert name to SMILES
smiles = resolver.resolve('caffeine', 'smiles')
print(f"Caffeine SMILES: {smiles}")

# Convert SMILES to InChI
inchi = resolver.resolve('CCO', 'stdinchi')  # Ethanol
print(f"Ethanol InChI: {inchi}")

# Get comprehensive molecular data
mol_data = resolver.get_molecular_data('aspirin')
print(f"Aspirin data: {mol_data}")
```

## Experimental Property Extraction

### Extract Properties from PubChem PUG View

```python
from provesid import PubChemView
from provesid.pubchemview import get_experimental_properties_table

# Initialize PubChem View
view = PubChemView()

# Get experimental melting points for aspirin
melting_points = view.get_experimental_properties(2244, 'Melting Point')
for prop in melting_points:
    print(f"Value: {prop.value}, Unit: {prop.unit}, Reference: {prop.reference_title}")

# Get structured table format
table = get_experimental_properties_table(2244, 'Melting Point')
print(table)
```

### Working with DataFrames

```python
# Convert to pandas DataFrame for analysis
df = view.experimental_properties_to_dataframe(2244, 'Boiling Point')
print(df.head())

# Analyze the data
print(f"Number of experimental values: {len(df)}")
print(f"Average boiling point: {df['Value'].mean():.1f} °C")
print(f"Reference sources: {df['Reference'].nunique()}")
```

## Batch Processing

### Process Multiple Compounds

```python
# List of compound CIDs to process
compound_cids = [2244, 2519, 3672]  # Aspirin, Caffeine, Ibuprofen

# Get properties for all compounds
all_properties = api.get_compound_properties(
    compound_cids,
    ['MolecularWeight', 'MolecularFormula']
)

for prop in all_properties['PropertyTable']['Properties']:
    cid = prop['CID']
    mw = prop['MolecularWeight']
    formula = prop['MolecularFormula']
    print(f"CID {cid}: {formula}, MW = {mw}")
```

### Batch Identifier Resolution

```python
# Convert multiple names to SMILES
compound_names = ['aspirin', 'caffeine', 'ibuprofen']
smiles_results = resolver.batch_resolve(compound_names, 'smiles')

for name, smiles in smiles_results.items():
    if smiles:
        print(f"{name}: {smiles}")
    else:
        print(f"{name}: Not found")
```

## Error Handling

PROVESID includes comprehensive error handling:

```python
from provesid.pubchem import PubChemNotFoundError, PubChemError
from provesid.resolver import NCIResolverNotFoundError

try:
    # This will raise an error for non-existent compound
    result = api.get_compound_by_cid(999999999)
except PubChemNotFoundError:
    print("Compound not found in PubChem")
except PubChemError as e:
    print(f"PubChem API error: {e}")

try:
    # This will raise an error for invalid identifier
    smiles = resolver.resolve('invalid_compound_name_xyz', 'smiles')
except NCIResolverNotFoundError:
    print("Identifier not found in NCI resolver")
```

## Working with Search Results

### Similarity and Substructure Searches

```python
# Find compounds similar to aspirin
similar_compounds = api.similarity_search('CCO', threshold=90)
print(f"Found {len(similar_compounds['IdentifierList']['CID'])} similar compounds")

# Search for compounds containing benzene ring
benzene_compounds = api.substructure_search('c1ccccc1')  # Benzene SMILES
print(f"Found {len(benzene_compounds['IdentifierList']['CID'])} compounds with benzene ring")
```

### Get Compound Synonyms

```python
# Get all synonyms for a compound
synonyms = api.get_compound_synonyms(2244)  # Aspirin
synonym_list = synonyms['InformationList']['Information'][0]['Synonym']
print(f"Aspirin synonyms: {synonym_list[:5]}")  # First 5 synonyms
```

## Configuration and Performance

### Adjust API Call Frequency

```python
# Slower API calls for rate-limited scenarios
slow_api = PubChemAPI(pause_time=1.0)  # 1 second between calls

# Faster calls for internal/unlimited access
fast_resolver = NCIChemicalIdentifierResolver(pause_time=0.1)
```

### Logging and Debugging

```python
import logging

# Enable debug logging to see API calls
logging.basicConfig(level=logging.DEBUG)

# Now API calls will be logged
result = api.get_compound_by_cid(2244)
```

## Next Steps

Now that you've learned the basics, explore more advanced features:

- [Property Extraction Examples](examples/property_extraction.md) - Advanced property extraction techniques
- [Batch Processing Guide](examples/batch_processing.md) - Efficient processing of large datasets
- [API Reference](api/pubchem.md) - Complete API documentation

## Common Patterns

### Property Extraction Pipeline

```python
def extract_property_data(cid, property_name):
    """Extract and summarize experimental property data"""
    view = PubChemView()
    
    # Get the data
    properties = view.get_experimental_properties(cid, property_name)
    if not properties:
        return None
    
    # Convert to DataFrame
    df = view.experimental_properties_to_dataframe(cid, property_name)
    
    # Basic statistics
    summary = {
        'count': len(df),
        'mean_value': df['Value'].mean() if len(df) > 0 else None,
        'references': df['Reference'].nunique(),
        'units': df['Unit'].unique().tolist()
    }
    
    return summary

# Use the pipeline
summary = extract_property_data(2244, 'Melting Point')
print(summary)
```

### Identifier Validation

```python
def validate_and_convert(identifier, target_format='smiles'):
    """Validate identifier and convert to target format"""
    resolver = NCIChemicalIdentifierResolver()
    
    try:
        result = resolver.resolve(identifier, target_format)
        return {'success': True, 'result': result, 'error': None}
    except Exception as e:
        return {'success': False, 'result': None, 'error': str(e)}

# Validate multiple identifiers
identifiers = ['aspirin', 'CCO', '50-78-2']  # Name, SMILES, CAS
for ident in identifiers:
    validation = validate_and_convert(ident)
    print(f"{ident}: {validation}")
```
