# CAS Common Chemistry API

The CASCommonChem class provides an interface to the CAS (Chemical Abstracts Service) Common Chemistry database using API v2.0 with authentication.

## Overview

CAS Common Chemistry provides reliable chemical information for over 500,000 chemical substances from CAS REGISTRY. This class enables searches by CAS Registry Number, chemical name, and SMILES notation using the authenticated v2.0 API gateway.

## Class: CASCommonChem

### Initialization

```python
from provesid import CASCommonChem, set_cas_api_key

# Recommended: Set API key once for persistent storage
set_cas_api_key("your-api-key-here")
cas_api = CASCommonChem()  # Automatically uses stored API key

# Alternative: Initialize with API key directly  
cas_api = CASCommonChem(api_key="your-api-key-here")

# Or initialize with API key from file
cas_api = CASCommonChem(api_key_file="path/to/api/key/file.txt")
```

**Parameters:**
- `api_key` (str, optional): CAS API key for authentication
- `api_key_file` (str, optional): Path to file containing API key
- `use_cache` (bool, default=True): Enable/disable result caching

**API Key Priority Order:**
1. Direct `api_key` parameter
2. `api_key_file` parameter  
3. Persistent configuration (via `set_cas_api_key()`)
4. Environment variables (`CCC_API_KEY`, `CAS_API_KEY`)

The CASCommonChem class initializes with:
- `base_url`: "https://commonchemistry.cas.org/api" (v2.0 gateway)
- API key authentication using X-API-KEY header
- Enhanced error handling for authentication failures
- Service-specific caching system

## Persistent API Key Configuration

The recommended approach is to configure your API key once using the persistent configuration system:

```python
# One-time setup
from provesid import set_cas_api_key, show_config

# Set your API key (stored securely in user config directory)
set_cas_api_key("your-cas-api-key")

# Verify configuration
show_config()

# Now use anywhere without specifying API key
from provesid import CASCommonChem
cas_api = CASCommonChem()  # Automatically loads stored API key
```

**Configuration Management:**
```python
from provesid import get_cas_api_key, remove_cas_api_key, show_config

# Check current configuration
show_config()

# Get stored API key
api_key = get_cas_api_key()

# Remove stored API key
remove_cas_api_key()
```

### Methods

#### `cas_to_detail(cas_number)`

Retrieve detailed information for a compound using its CAS Registry Number.

**Parameters:**
- `cas_number` (str): CAS Registry Number (e.g., "64-17-5")

**Returns:**
- `dict` or `None`: Detailed compound information if found, None if not found or error

**Example:**
```python
cas_api = CASCommonChem(api_key="your-api-key")

# Get details for ethanol (CAS: 64-17-5)
ethanol_data = cas_api.cas_to_detail("64-17-5")

if ethanol_data and ethanol_data['found']:
    print(f"Status: {ethanol_data['status']}")
    print(f"Name: {ethanol_data['name']}")
    print(f"Molecular Formula: {ethanol_data['molecularFormula']}")
    print(f"SMILES: {ethanol_data['smile']}")
    print(f"InChI: {ethanol_data['inchi']}")
    print(f"InChI Key: {ethanol_data['inchiKey']}")
else:
    print(f"Error: {ethanol_data['status'] if ethanol_data else 'No response'}")
```

#### `name_to_detail(name)`

Search for compound information using chemical name.

**Parameters:**
- `name` (str): Chemical name (common or systematic)

**Returns:**
- `dict` or `None`: Detailed compound information if found, None if not found or error

**Example:**
```python
cas_api = CASCommonChem()

# Search by common name
ethanol_data = cas_api.name_to_detail("ethanol")

# Search by systematic name
ethanol_data = cas_api.name_to_detail("ethyl alcohol")

if ethanol_data:
    print(f"CAS Number: {ethanol_data['rn']}")
    print(f"Molecular Formula: {ethanol_data['molecularFormula']}")
```

#### `smiles_to_detail(smiles)`

Search for compound information using SMILES notation.

**Parameters:**
- `smiles` (str): SMILES string representation

**Returns:**
- `dict` or `None`: Detailed compound information if found, None if not found or error

**Example:**
```python
cas_api = CASCommonChem()

# Search by SMILES
compound_data = cas_api.smiles_to_detail("CCO")  # Ethanol

if compound_data:
    print(f"Name: {compound_data['name']}")
    print(f"CAS Number: {compound_data['rn']}")
```

### Data Structure

The API returns detailed compound information in the following structure:

```python
{
    "uri": "substance/pt/64175",
    "rn": "64-17-5",                    # CAS Registry Number
    "name": "Ethanol",                  # Primary name
    "image": "<svg>...</svg>",          # Chemical structure image (SVG)
    "inchi": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
    "inchiKey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
    "smile": "CCO",                     # SMILES notation
    "canonicalSmile": "CCO",            # Canonical SMILES
    "molecularFormula": "C<sub>2</sub>H<sub>6</sub>O",
    "molecularMass": "46.07",           # Molecular mass in g/mol
    "experimentalProperties": [         # Experimental data
        {
            "name": "Boiling Point",
            "property": "78.37 °C",
            "sourceNumber": 1
        }
    ],
    "propertyCitations": [              # Literature references
        {
            "docId": 1,
            "title": "CRC Handbook of Chemistry and Physics"
        }
    ],
    "synonyms": [                       # Alternative names
        "Ethyl alcohol",
        "Grain alcohol",
        "EtOH"
    ],
    "replacedRns": [],                  # Historical CAS numbers
    "hasMolfile": true                  # Molfile availability
}
```

### Comprehensive Usage Examples

#### Single Compound Lookup

```python
from provesid.cascommonchem import CASCommonChem

cas_api = CASCommonChem()

def get_compound_info(identifier, search_type="cas"):
    """Get comprehensive compound information"""
    
    if search_type == "cas":
        data = cas_api.cas_to_detail(identifier)
    elif search_type == "name":
        data = cas_api.name_to_detail(identifier)
    elif search_type == "smiles":
        data = cas_api.smiles_to_detail(identifier)
    else:
        print("Invalid search type")
        return None
    
    if data:
        print(f"=== {data['name']} ===")
        print(f"CAS Number: {data['rn']}")
        print(f"Molecular Formula: {data['molecularFormula']}")
        print(f"Molecular Mass: {data['molecularMass']} g/mol")
        print(f"SMILES: {data['smile']}")
        print(f"InChI Key: {data['inchiKey']}")
        
        # Synonyms
        if 'synonyms' in data and data['synonyms']:
            print(f"Synonyms: {', '.join(data['synonyms'][:5])}")  # First 5
        
        # Experimental properties
        if 'experimentalProperties' in data:
            print("Experimental Properties:")
            for prop in data['experimentalProperties']:
                print(f"  - {prop['name']}: {prop['property']}")
        
        return data
    else:
        print(f"No data found for {identifier}")
        return None

# Usage examples
ethanol_by_cas = get_compound_info("64-17-5", "cas")
ethanol_by_name = get_compound_info("ethanol", "name")
ethanol_by_smiles = get_compound_info("CCO", "smiles")
```

#### Batch Processing

```python
def batch_cas_lookup(cas_numbers):
    """Process multiple CAS numbers with error handling"""
    cas_api = CASCommonChem()
    results = {}
    
    for cas_num in cas_numbers:
        print(f"Processing CAS: {cas_num}")
        
        try:
            data = cas_api.cas_to_detail(cas_num)
            
            if data:
                results[cas_num] = {
                    'name': data['name'],
                    'formula': data['molecularFormula'],
                    'smiles': data['smile'],
                    'mass': data['molecularMass']
                }
                print(f"  ✓ Found: {data['name']}")
            else:
                results[cas_num] = None
                print(f"  ✗ Not found")
                
        except Exception as e:
            print(f"  ✗ Error: {e}")
            results[cas_num] = None
            
        # Rate limiting - be respectful to the API
        time.sleep(0.5)
    
    return results

# Example usage
cas_list = ["64-17-5", "67-56-1", "108-88-3", "71-43-2"]
batch_results = batch_cas_lookup(cas_list)

for cas, data in batch_results.items():
    if data:
        print(f"{cas}: {data['name']} ({data['formula']})")
    else:
        print(f"{cas}: Not found")
```

#### Cross-Reference with Other APIs

```python
from provesid.cascommonchem import CASCommonChem
from provesid.pubchem import PubChemAPI
from provesid.opsin import OPSIN

def comprehensive_compound_lookup(compound_name):
    """Cross-reference compound across multiple databases"""
    
    results = {
        'input_name': compound_name,
        'cas_data': None,
        'pubchem_data': None,
        'opsin_data': None,
        'cross_validated': False
    }
    
    # 1. Get structure from OPSIN
    opsin = OPSIN()
    opsin_result = opsin.get_id(compound_name)
    
    if opsin_result['status'] == 'SUCCESS':
        results['opsin_data'] = opsin_result
        smiles = opsin_result['smiles']
        
        # 2. Search CAS by SMILES
        cas_api = CASCommonChem()
        cas_data = cas_api.smiles_to_detail(smiles)
        
        if cas_data:
            results['cas_data'] = cas_data
            
            # 3. Search PubChem by CAS
            pubchem = PubChemAPI()
            pubchem_data = pubchem.get_compound_by_name(cas_data['rn'])
            
            if pubchem_data:
                results['pubchem_data'] = pubchem_data
                results['cross_validated'] = True
    
    # Alternative: try direct name search in CAS
    if not results['cas_data']:
        cas_data = cas_api.name_to_detail(compound_name)
        if cas_data:
            results['cas_data'] = cas_data
    
    return results

# Usage
compound_info = comprehensive_compound_lookup("caffeine")

if compound_info['cross_validated']:
    print("Successfully cross-validated across all databases!")
    print(f"CAS: {compound_info['cas_data']['rn']}")
    print(f"SMILES: {compound_info['opsin_data']['smiles']}")
    print(f"PubChem CID: {compound_info['pubchem_data'].get('CID')}")
```

#### Property Analysis

```python
def analyze_experimental_properties(cas_number):
    """Analyze experimental properties from CAS data"""
    cas_api = CASCommonChem()
    data = cas_api.cas_to_detail(cas_number)
    
    if not data or 'experimentalProperties' not in data:
        print("No experimental properties available")
        return None
    
    properties = data['experimentalProperties']
    
    print(f"=== Properties for {data['name']} ===")
    
    # Categorize properties
    physical_props = []
    thermal_props = []
    other_props = []
    
    for prop in properties:
        prop_name = prop['name'].lower()
        
        if any(term in prop_name for term in ['boiling', 'melting', 'flash']):
            thermal_props.append(prop)
        elif any(term in prop_name for term in ['density', 'refractive', 'viscosity']):
            physical_props.append(prop)
        else:
            other_props.append(prop)
    
    # Display categorized properties
    categories = [
        ("Thermal Properties", thermal_props),
        ("Physical Properties", physical_props),
        ("Other Properties", other_props)
    ]
    
    for category_name, props in categories:
        if props:
            print(f"\n{category_name}:")
            for prop in props:
                print(f"  {prop['name']}: {prop['property']}")
    
    return {
        'thermal': thermal_props,
        'physical': physical_props,
        'other': other_props
    }

# Example
analyze_experimental_properties("64-17-5")  # Ethanol
```

### Best Practices

#### 1. Error Handling and Validation

```python
def safe_cas_lookup(cas_number):
    """Safely lookup CAS number with validation"""
    import re
    
    # Validate CAS number format
    cas_pattern = r'^\d{1,7}-\d{2}-\d$'
    if not re.match(cas_pattern, cas_number):
        print(f"Invalid CAS format: {cas_number}")
        return None
    
    cas_api = CASCommonChem()
    
    try:
        data = cas_api.cas_to_detail(cas_number)
        
        if data:
            # Validate essential fields
            required_fields = ['name', 'molecularFormula', 'smile']
            missing_fields = [field for field in required_fields if field not in data]
            
            if missing_fields:
                print(f"Warning: Missing fields: {missing_fields}")
            
            return data
        else:
            print(f"No data found for CAS {cas_number}")
            return None
            
    except Exception as e:
        print(f"Error looking up CAS {cas_number}: {e}")
        return None

# Usage
data = safe_cas_lookup("64-17-5")
```

#### 2. Rate Limiting and Caching

```python
import time
from functools import lru_cache

class CachedCASCommonChem:
    """CAS API with caching and rate limiting"""
    
    def __init__(self, delay=0.5):
        self.api = CASCommonChem()
        self.delay = delay
        self.last_request_time = 0
    
    def _rate_limit(self):
        """Implement rate limiting"""
        elapsed = time.time() - self.last_request_time
        if elapsed < self.delay:
            time.sleep(self.delay - elapsed)
        self.last_request_time = time.time()
    
    @lru_cache(maxsize=1000)
    def cached_cas_lookup(self, cas_number):
        """Cached CAS lookup"""
        self._rate_limit()
        return self.api.cas_to_detail(cas_number)
    
    @lru_cache(maxsize=1000)
    def cached_name_lookup(self, name):
        """Cached name lookup"""
        self._rate_limit()
        return self.api.name_to_detail(name)

# Usage
cached_api = CachedCASCommonChem(delay=1.0)

# First call - hits API
data1 = cached_api.cached_cas_lookup("64-17-5")

# Second call - uses cache
data2 = cached_api.cached_cas_lookup("64-17-5")
```

#### 3. Data Export and Analysis

```python
import pandas as pd
import json

def export_cas_data(cas_numbers, output_format='csv'):
    """Export CAS data to various formats"""
    cas_api = CASCommonChem()
    results = []
    
    for cas_num in cas_numbers:
        data = cas_api.cas_to_detail(cas_num)
        
        if data:
            # Flatten data for tabular export
            flat_data = {
                'cas_number': data['rn'],
                'name': data['name'],
                'molecular_formula': data['molecularFormula'],
                'molecular_mass': data['molecularMass'],
                'smiles': data['smile'],
                'inchi_key': data['inchiKey'],
                'synonyms_count': len(data.get('synonyms', [])),
                'properties_count': len(data.get('experimentalProperties', []))
            }
            
            # Add first few synonyms
            synonyms = data.get('synonyms', [])
            for i in range(min(3, len(synonyms))):
                flat_data[f'synonym_{i+1}'] = synonyms[i]
            
            results.append(flat_data)
    
    if output_format == 'csv':
        df = pd.DataFrame(results)
        df.to_csv('cas_data.csv', index=False)
        print(f"Exported {len(results)} compounds to cas_data.csv")
        
    elif output_format == 'json':
        with open('cas_data.json', 'w') as f:
            json.dump(results, f, indent=2)
        print(f"Exported {len(results)} compounds to cas_data.json")
    
    return results

# Usage
cas_list = ["64-17-5", "67-56-1", "108-88-3"]
exported_data = export_cas_data(cas_list, 'csv')
```

### Integration Patterns

#### Database Integration

```python
import sqlite3

def create_compound_database():
    """Create SQLite database for compound data"""
    conn = sqlite3.connect('compounds.db')
    
    conn.execute('''
        CREATE TABLE IF NOT EXISTS compounds (
            cas_number TEXT PRIMARY KEY,
            name TEXT,
            molecular_formula TEXT,
            molecular_mass REAL,
            smiles TEXT,
            inchi_key TEXT,
            data_json TEXT,
            last_updated TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    conn.commit()
    return conn

def store_cas_data(cas_number, force_update=False):
    """Store or update CAS data in database"""
    conn = create_compound_database()
    
    # Check if already exists
    cursor = conn.execute(
        'SELECT cas_number FROM compounds WHERE cas_number = ?',
        (cas_number,)
    )
    
    exists = cursor.fetchone() is not None
    
    if exists and not force_update:
        print(f"CAS {cas_number} already in database")
        return
    
    # Fetch from API
    cas_api = CASCommonChem()
    data = cas_api.cas_to_detail(cas_number)
    
    if data:
        conn.execute('''
            INSERT OR REPLACE INTO compounds 
            (cas_number, name, molecular_formula, molecular_mass, 
             smiles, inchi_key, data_json)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        ''', (
            data['rn'],
            data['name'],
            data['molecularFormula'],
            float(data['molecularMass']),
            data['smile'],
            data['inchiKey'],
            json.dumps(data)
        ))
        
        conn.commit()
        print(f"Stored {data['name']} (CAS: {cas_number})")
    
    conn.close()

# Usage
store_cas_data("64-17-5")
```

### Limitations and Considerations

#### Data Coverage
- Free tier covers ~500,000 substances
- Focus on commonly used chemicals
- May not include very specialized or proprietary compounds

#### Rate Limiting
- No official documented limits
- Recommended to implement delays between requests
- Monitor for HTTP 429 (Too Many Requests) responses

#### Data Quality
- High-quality data from CAS REGISTRY
- Experimental properties may vary in precision
- Cross-validation with other sources recommended

#### API Stability
- Production-grade API with good uptime
- Data structure may evolve over time
- Always check for API updates and changes

## Cache Management

The CASCommonChem class includes service-specific caching to improve performance and reduce API calls.

### Cache Control

```python
# Enable/disable caching during initialization
cas_api = CASCommonChem(use_cache=True)  # Default

# Force fresh data retrieval (bypasses cache but still stores results)
cas_api.use_cache = False
fresh_data = cas_api.cas_to_detail("64-17-5")

# Re-enable caching
cas_api.use_cache = True
```

### Cache Management Methods

```python
# Get cache information
cache_info = cas_api.get_cache_info()
print(f"Cache directory: {cache_info['cache_directory']}")
print(f"Memory entries: {cache_info['memory_entries']}")
print(f"Disk entries: {cache_info['disk_entries']}")
print(f"Total size: {cache_info['total_size_mb']:.2f} MB")

# Clear CAS-specific cache only
cas_api.clear_cache()
```

## Authentication Requirements

The CAS Common Chemistry API v2.0 requires authentication:

1. **API Key**: Obtain from CAS SciFinder or institutional access
2. **Header Authentication**: Uses `X-API-KEY` header
3. **Error Handling**: Comprehensive authentication error management

### API Key Setup

```python
# Method 1: Direct API key
cas_api = CASCommonChem(api_key="your-cas-api-key")

# Method 2: File-based key
cas_api = CASCommonChem(api_key_file="/path/to/api-key.txt")

# Method 3: Default DTU location (if applicable)
cas_api = CASCommonChem()  # Uses DTU OneDrive path
```

## Migration from v1.0

If migrating from the previous v1.0 implementation:

1. **Add API Key**: v2.0 requires authentication
2. **Error Handling**: New authentication error responses
3. **Status Fields**: Enhanced status reporting with `found` boolean
4. **Cache System**: Now uses service-specific caching
