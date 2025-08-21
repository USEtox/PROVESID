# PubChem View API

The PubChem View module provides access to experimental property data from PubChem's PUG View service. This includes extraction of experimental values, units, and comprehensive reference information.

::: provesid.pubchemview

## Quick Start

```python
from provesid import PubChemView
from provesid.pubchemview import get_experimental_properties_table

# Initialize the view client
view = PubChemView()

# Get experimental melting points for aspirin (CID 2244)
properties = view.get_experimental_properties(2244, 'Melting Point')
for prop in properties:
    print(f"Value: {prop.value} {prop.unit}")
    print(f"Reference: {prop.reference_title}")

# Get a structured DataFrame
df = view.experimental_properties_to_dataframe(2244, 'Melting Point')
print(df.head())

# Use the convenience function for a complete table
table = get_experimental_properties_table(2244, 'Boiling Point')
print(table)
```

## Available Property Types

The PubChem View service provides access to various experimental properties:

### Physical Properties
- **Melting Point** - Melting point temperatures
- **Boiling Point** - Boiling point temperatures  
- **Density** - Density measurements
- **Vapor Pressure** - Vapor pressure data
- **Solubility** - Solubility in various solvents
- **LogP** - Partition coefficient data
- **Viscosity** - Viscosity measurements
- **Refractive Index** - Refractive index values

### Spectroscopic Properties
- **UV/Vis Spectrum** - UV-Visible spectroscopy data
- **IR Spectrum** - Infrared spectroscopy data
- **NMR Spectrum** - Nuclear magnetic resonance data
- **Mass Spectrum** - Mass spectrometry data

### Safety and Toxicity
- **Flash Point** - Flash point temperatures
- **Auto-Ignition Temperature** - Auto-ignition data
- **LD50** - Lethal dose data
- **LC50** - Lethal concentration data

## Data Structures

### PropertyData Class

The `PropertyData` dataclass represents a single experimental property measurement:

```python
@dataclass
class PropertyData:
    cid: int
    heading: str
    string_with_markup: str
    value: Optional[str]
    unit: Optional[str] 
    reference_number: Optional[str]
    reference_title: Optional[str]
    reference_authors: Optional[str]
    reference_journal: Optional[str]
    reference_year: Optional[str]
    reference_doi: Optional[str]
    reference_pmid: Optional[str]
    full_reference: Optional[str]
```

**Fields:**
- `cid`: PubChem Compound ID
- `heading`: Property type (e.g., "Melting Point")
- `string_with_markup`: Original text with markup
- `value`: Extracted numeric/text value
- `unit`: Unit of measurement
- `reference_*`: Citation information
- `full_reference`: Complete formatted reference

## Advanced Usage

### Custom Value Parsing

The module includes sophisticated value parsing that handles:

```python
# Complex value strings
examples = [
    "139-140 °C",           # Range values
    "25.5 ± 0.2 °C",       # Values with uncertainty
    "< 100 °C",            # Comparison operators
    "decomp. at 180 °C",   # Qualitative descriptions
    "760 mmHg at 20 °C"    # Conditional values
]

# The parser extracts the main numeric value
for example in examples:
    # Internal parsing would extract the primary value
    pass
```

### Reference Information Extraction

Complete bibliographic information is extracted and structured:

```python
# Get properties with full reference details
properties = view.get_experimental_properties(2244, 'Melting Point')
for prop in properties:
    if prop.reference_doi:
        print(f"DOI: {prop.reference_doi}")
    if prop.reference_pmid:
        print(f"PubMed ID: {prop.reference_pmid}")
    if prop.full_reference:
        print(f"Full citation: {prop.full_reference}")
```

### DataFrame Operations

Convert to pandas DataFrame for data analysis:

```python
import pandas as pd

# Get DataFrame with all experimental data
df = view.experimental_properties_to_dataframe(2244, 'Solubility')

# Filter by specific units
water_solubility = df[df['Unit'].str.contains('g/L', na=False)]

# Group by reference source
by_source = df.groupby('Reference').agg({
    'Value': ['count', 'mean'],
    'Unit': lambda x: list(x.unique())
})

# Export to CSV
df.to_csv('solubility_data.csv', index=False)
```

## Error Handling

The module provides specific exception classes:

```python
from provesid.pubchemview import PubChemViewError, PubChemViewNotFoundError

try:
    properties = view.get_experimental_properties(999999, 'Melting Point')
except PubChemViewNotFoundError:
    print("Compound or property not found")
except PubChemViewError as e:
    print(f"API error: {e}")
```

## Batch Processing

Process multiple compounds efficiently:

```python
def batch_extract_properties(cids, property_type):
    """Extract properties for multiple compounds"""
    results = {}
    view = PubChemView()
    
    for cid in cids:
        try:
            df = view.experimental_properties_to_dataframe(cid, property_type)
            if not df.empty:
                results[cid] = {
                    'count': len(df),
                    'mean_value': df['Value'].mean() if df['Value'].notna().any() else None,
                    'units': df['Unit'].unique().tolist()
                }
        except Exception as e:
            results[cid] = {'error': str(e)}
    
    return results

# Process multiple compounds
cids = [2244, 2519, 3672]  # Aspirin, Caffeine, Ibuprofen
melting_data = batch_extract_properties(cids, 'Melting Point')
```

## Integration with Other APIs

Combine with PubChem compound data:

```python
from provesid import PubChemAPI, PubChemView

def comprehensive_compound_analysis(cid):
    """Get both computed and experimental data"""
    api = PubChemAPI()
    view = PubChemView()
    
    # Get computed properties
    computed = api.get_compound_properties(
        [cid], 
        ['MolecularWeight', 'MolecularFormula', 'ConnectivitySMILES']
    )
    
    # Get experimental properties
    experimental = {}
    for prop_type in ['Melting Point', 'Boiling Point', 'Density']:
        try:
            df = view.experimental_properties_to_dataframe(cid, prop_type)
            if not df.empty:
                experimental[prop_type] = df
        except:
            pass
    
    return {
        'computed': computed,
        'experimental': experimental
    }

# Analyze aspirin
analysis = comprehensive_compound_analysis(2244)
```

## Performance Considerations

### Rate Limiting

The PubChemView client includes automatic rate limiting:

```python
# Adjust request frequency for large batch jobs
view = PubChemView(pause_time=1.0)  # 1 second between requests

# For development/testing with faster requests
view_fast = PubChemView(pause_time=0.1)  # 100ms between requests
```

### Caching

Consider implementing caching for frequently accessed data:

```python
import pickle
from pathlib import Path

def cached_property_extraction(cid, property_type, cache_dir='cache'):
    """Extract properties with file-based caching"""
    cache_path = Path(cache_dir) / f"{cid}_{property_type}.pkl"
    
    if cache_path.exists():
        with open(cache_path, 'rb') as f:
            return pickle.load(f)
    
    # Extract fresh data
    view = PubChemView()
    df = view.experimental_properties_to_dataframe(cid, property_type)
    
    # Cache the results
    cache_path.parent.mkdir(exist_ok=True)
    with open(cache_path, 'wb') as f:
        pickle.dump(df, f)
    
    return df
```

## See Also

- [PubChem API](pubchem.md) - For computed compound properties
- [Property Extraction Examples](../examples/property_extraction.md) - Detailed usage examples
- [Batch Processing Guide](../examples/batch_processing.md) - Efficient processing strategies
