# NCI Chemical Identifier Resolver API

The NCI Chemical Identifier Resolver provides a comprehensive interface to the NCI CADD Group's Chemical Identifier Resolver web service for converting between different chemical structure identifiers.

::: provesid.resolver

## Quick Start

```python
from provesid import NCIChemicalIdentifierResolver

# Initialize the resolver
resolver = NCIChemicalIdentifierResolver()

# Convert name to SMILES
smiles = resolver.resolve('aspirin', 'smiles')
print(f"Aspirin SMILES: {smiles}")

# Convert SMILES to InChI
inchi = resolver.resolve('CCO', 'stdinchi')  # Ethanol
print(f"Ethanol InChI: {inchi}")

# Get comprehensive molecular data
mol_data = resolver.get_molecular_data('caffeine')
print(f"Caffeine data: {mol_data}")
```

## Supported Representations

The NCI Resolver supports conversion between numerous chemical identifier formats:

### Structure Identifiers
- **smiles** - Unique SMILES strings
- **stdinchi** - Standard InChI identifiers
- **stdinchikey** - Standard InChI Keys
- **ficts** - NCI/CADD FICTS identifiers
- **ficus** - NCI/CADD FICuS identifiers
- **uuuuu** - NCI/CADD uuuuu identifiers
- **hashisy** - CACTVS HASHISY hashcodes

### Chemical Names and Properties
- **names** - Chemical names list
- **iupac_name** - IUPAC systematic names
- **cas** - CAS Registry Numbers
- **mw** - Molecular weight
- **formula** - Molecular formula
- **exactmass** - Exact molecular mass

### Physical Properties
- **charge** - Formal charge
- **h_bond_acceptor_count** - Hydrogen bond acceptor count
- **h_bond_donor_count** - Hydrogen bond donor count
- **rotor_count** - Rotatable bond count
- **ring_count** - Ring count

### File Formats
- **sdf** - Structure Data File format
- **image** - Chemical structure images

## Core Methods

### Basic Resolution

```python
# Single identifier conversion
smiles = resolver.resolve('aspirin', 'smiles')
inchi = resolver.resolve('50-78-2', 'stdinchi')  # CAS to InChI
formula = resolver.resolve('CCO', 'formula')     # SMILES to formula
```

### Multiple Representations

```python
# Get multiple representations at once
representations = ['smiles', 'stdinchi', 'mw', 'formula']
results = resolver.resolve_multiple('caffeine', representations)

print(results)
# {
#     'smiles': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
#     'stdinchi': 'InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3',
#     'mw': '194.1906',
#     'formula': 'C8H10N4O2'
# }
```

### Comprehensive Molecular Data

```python
# Get extensive molecular information
mol_data = resolver.get_molecular_data('aspirin')

# Access individual properties
print(f"SMILES: {mol_data['smiles']}")
print(f"Formula: {mol_data['formula']}")
print(f"Molecular Weight: {mol_data['mw']}")
print(f"Names: {mol_data['names']}")
print(f"CAS Number: {mol_data['cas']}")
```

## Convenience Functions

The module provides simple functions for common conversions:

### CAS Number Conversions

```python
from provesid.resolver import nci_cas_to_mol, nci_cas_to_inchi

# Convert CAS to comprehensive molecular data
mol_data = nci_cas_to_mol('50-78-2')  # Aspirin CAS
print(mol_data['formula'])  # C9H8O4

# Convert CAS to InChI
inchi = nci_cas_to_inchi('64-17-5')  # Ethanol CAS
```

### Name-Based Conversions

```python
from provesid.resolver import nci_name_to_smiles, nci_get_molecular_weight

# Convert name to SMILES
smiles = nci_name_to_smiles('caffeine')

# Get molecular weight from name
mw = nci_get_molecular_weight('water')
print(f"Water MW: {mw}")  # ~18.015
```

### SMILES Conversions

```python
from provesid.resolver import nci_smiles_to_names, nci_inchi_to_smiles

# Get names for a SMILES string
names = nci_smiles_to_names('CCO')  # Ethanol
print(names)  # ['ethanol', 'ethyl alcohol', ...]

# Convert InChI to SMILES
smiles = nci_inchi_to_smiles('InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3')
```

## Batch Processing

### Process Multiple Identifiers

```python
# Batch resolve multiple compounds to SMILES
compounds = ['aspirin', 'caffeine', 'ibuprofen', 'water']
smiles_results = resolver.batch_resolve(compounds, 'smiles')

for compound, smiles in smiles_results.items():
    if smiles:
        print(f"{compound}: {smiles}")
    else:
        print(f"{compound}: Not found")
```

### Validation and Filtering

```python
# Check which identifiers are valid
test_compounds = ['aspirin', 'invalid_name_xyz', 'caffeine', 'fake_compound']
valid_compounds = []

for compound in test_compounds:
    if resolver.is_valid_identifier(compound):
        valid_compounds.append(compound)
        print(f"✓ {compound}")
    else:
        print(f"✗ {compound}")

print(f"Valid compounds: {valid_compounds}")
```

## Working with Images

### Generate Structure Images

```python
# Get image URL for a compound
image_url = resolver.get_image_url('aspirin', image_format='png', width=400, height=400)
print(f"Aspirin structure image: {image_url}")

# Download image to file
success = resolver.download_image('caffeine', 'caffeine_structure.png')
if success:
    print("Image downloaded successfully")
```

## Error Handling

The module provides specific exception types for different error conditions:

```python
from provesid.resolver import (
    NCIResolverError, 
    NCIResolverNotFoundError, 
    NCIResolverTimeoutError
)

try:
    smiles = resolver.resolve('definitely_not_a_chemical', 'smiles')
except NCIResolverNotFoundError:
    print("Chemical identifier not found")
except NCIResolverTimeoutError:
    print("Request timed out")
except NCIResolverError as e:
    print(f"General resolver error: {e}")
```

## Advanced Usage

### Custom Configuration

```python
# Configure with custom settings
custom_resolver = NCIChemicalIdentifierResolver(
    base_url="https://cactus.nci.nih.gov/chemical/structure",
    timeout=60,  # Longer timeout
    pause_time=0.5  # Slower requests
)

# Use custom resolver
result = custom_resolver.resolve('aspirin', 'smiles')
```

### Partial Name Searching

```python
# Search for compounds by partial name
matches = resolver.search_by_partial_name('acetyl')
print(f"Found compounds with 'acetyl': {matches}")
```

### Available Representations

```python
# Get list of all available representations
representations = resolver.get_available_representations()
print(f"Available formats: {representations}")
```

## Integration Examples

### Combine with PubChem Data

```python
from provesid import PubChemAPI, NCIChemicalIdentifierResolver

def cross_validate_identifiers(compound_name):
    """Validate identifier across multiple services"""
    resolver = NCIChemicalIdentifierResolver()
    api = PubChemAPI()
    
    # Get SMILES from NCI
    nci_smiles = resolver.resolve(compound_name, 'smiles')
    
    # Get compound from PubChem using the same name
    pubchem_cids = api.get_cids_by_name(compound_name)
    
    if pubchem_cids and nci_smiles:
        # Get PubChem SMILES for comparison
        cid = pubchem_cids['IdentifierList']['CID'][0]
        pubchem_props = api.get_compound_properties([cid], ['ConnectivitySMILES'])
        pubchem_smiles = pubchem_props['PropertyTable']['Properties'][0]['ConnectivitySMILES']
        
        return {
            'name': compound_name,
            'nci_smiles': nci_smiles,
            'pubchem_smiles': pubchem_smiles,
            'pubchem_cid': cid,
            'match': nci_smiles == pubchem_smiles
        }
    
    return None

# Cross-validate aspirin
validation = cross_validate_identifiers('aspirin')
print(validation)
```

### Data Pipeline Integration

```python
def chemical_identifier_pipeline(identifiers, target_format='smiles'):
    """Process multiple identifiers through resolution pipeline"""
    resolver = NCIChemicalIdentifierResolver()
    results = []
    
    for identifier in identifiers:
        try:
            # Attempt resolution
            result = resolver.resolve(identifier, target_format)
            
            # Get additional data if successful
            mol_data = resolver.get_molecular_data(identifier)
            
            results.append({
                'input': identifier,
                'output': result,
                'formula': mol_data.get('formula'),
                'mw': mol_data.get('mw'),
                'status': 'success'
            })
            
        except Exception as e:
            results.append({
                'input': identifier,
                'output': None,
                'formula': None,
                'mw': None,
                'status': f'error: {str(e)}'
            })
    
    return results

# Process a mixed list of identifiers
identifiers = ['aspirin', 'CCO', '50-78-2', 'invalid_name']
pipeline_results = chemical_identifier_pipeline(identifiers)
```

## Performance Considerations

### Rate Limiting

The resolver includes automatic rate limiting to respect server limits:

```python
# Default rate limiting (3 requests per second)
resolver = NCIChemicalIdentifierResolver()

# Slower rate for large batch jobs
slow_resolver = NCIChemicalIdentifierResolver(pause_time=1.0)

# Faster rate for development (use with caution)
fast_resolver = NCIChemicalIdentifierResolver(pause_time=0.05)
```

### Caching Strategies

Implement caching for frequently accessed data:

```python
import json
from pathlib import Path

class CachedResolver:
    def __init__(self, cache_file='resolver_cache.json'):
        self.resolver = NCIChemicalIdentifierResolver()
        self.cache_file = Path(cache_file)
        self.cache = self._load_cache()
    
    def _load_cache(self):
        if self.cache_file.exists():
            with open(self.cache_file, 'r') as f:
                return json.load(f)
        return {}
    
    def _save_cache(self):
        with open(self.cache_file, 'w') as f:
            json.dump(self.cache, f, indent=2)
    
    def resolve(self, identifier, representation):
        cache_key = f"{identifier}_{representation}"
        
        if cache_key in self.cache:
            return self.cache[cache_key]
        
        # Resolve and cache
        result = self.resolver.resolve(identifier, representation)
        self.cache[cache_key] = result
        self._save_cache()
        
        return result

# Use cached resolver
cached = CachedResolver()
smiles = cached.resolve('aspirin', 'smiles')  # First call - API request
smiles2 = cached.resolve('aspirin', 'smiles')  # Second call - cached
```

## See Also

- [PubChem API](pubchem.md) - For complementary compound data
- [PubChem View](pubchemview.md) - For experimental property data
- [Basic Usage Examples](../examples/basic_usage.md) - Common use cases
