# ClassyFire API

The ClassyFireAPI class provides an interface to the ClassyFire chemical classification system, which automatically assigns chemical taxonomies to molecular structures.

## Overview

ClassyFire is a web-based application for the automated structural classification of chemical entities. This class provides a Python interface to submit queries, check status, and retrieve classification results.

## Class: ClassyFireAPI

### Class Attributes

- `URL`: Base URL for the ClassyFire API ("http://classyfire.wishartlab.com")

### Static Methods

All methods are static and can be called directly on the class without instantiation.

#### `submit_query(label, input, type='STRUCTURE')`

Submit a chemical structure for classification.

**Parameters:**
- `label` (str): User-defined label for the query
- `input` (str): Chemical structure (SMILES, InChI, etc.)
- `type` (str, optional): Query type, default is 'STRUCTURE'

**Returns:**
- `requests.Response`: HTTP response object from the API

**Example:**
```python
from provesid.classyfire import ClassyFireAPI

# Submit a query
response = ClassyFireAPI.submit_query("Ethanol Example", "CCO")

if response.status_code == 200:
    query_data = response.json()
    query_id = query_data['id']
    print(f"Query submitted successfully. ID: {query_id}")
else:
    print(f"Submission failed: {response.status_code}")
```

#### `query_status(query_id)`

Check the status of a submitted query.

**Parameters:**
- `query_id` (str): The ID of the query to check

**Returns:**
- `requests.Response`: HTTP response object, or `None` if request fails

**Example:**
```python
# Check query status
status_response = ClassyFireAPI.query_status(query_id)

if status_response and status_response.status_code == 200:
    status_data = status_response.json()
    print(f"Query status: {status_data.get('status', 'Unknown')}")
```

#### `get_query(query_id, format="json")`

Retrieve the results of a completed query.

**Parameters:**
- `query_id` (str): The ID of the query to retrieve
- `format` (str, optional): Output format - "json", "sdf", or "csv" (default: "json")

**Returns:**
- `requests.Response`: HTTP response object with the results

**Example:**
```python
# Get query results
results_response = ClassyFireAPI.get_query(query_id, format="json")

if results_response.status_code == 200:
    classification = results_response.json()
    print("Classification Results:")
    print(f"Kingdom: {classification.get('kingdom', {}).get('name')}")
    print(f"Superclass: {classification.get('superclass', {}).get('name')}")
    print(f"Class: {classification.get('class', {}).get('name')}")
    print(f"Subclass: {classification.get('subclass', {}).get('name')}")
```

### Complete Workflow Example

```python
import time
import json
from provesid.classyfire import ClassyFireAPI

# SMILES string for para-chloronitrobenzene
smiles = "C1=CC(=CC=C1[N+](=O)[O-])Cl"
label = "Para-chloronitrobenzene Classification"

# Step 1: Submit the query
print("Submitting query...")
response = ClassyFireAPI.submit_query(label, smiles)

if response.status_code == 200:
    query_result = response.json()
    query_id = query_result['id']
    print(f"Query submitted successfully. Query ID: {query_id}")
    
    # Step 2: Wait for processing
    print("Waiting for processing...")
    time.sleep(10)  # Initial wait
    
    # Step 3: Check status periodically
    max_attempts = 30
    for attempt in range(max_attempts):
        status_response = ClassyFireAPI.query_status(query_id)
        
        if status_response and status_response.status_code == 200:
            status_data = status_response.json()
            status = status_data.get('status', 'Unknown')
            
            print(f"Attempt {attempt + 1}: Status = {status}")
            
            if status == 'Done':
                break
            elif status == 'In Queue' or status == 'In Progress':
                time.sleep(30)  # Wait 30 seconds before next check
            else:
                print(f"Unexpected status: {status}")
                break
        else:
            print("Failed to check status")
            break
    
    # Step 4: Retrieve results
    if status == 'Done':
        print("Retrieving classification results...")
        result_response = ClassyFireAPI.get_query(query_id, format="json")
        
        if result_response.status_code == 200:
            classification_results = result_response.json()
            print("\nClassification Results:")
            print(json.dumps(classification_results, indent=2))
            
            # Extract key classification levels
            print("\nKey Classifications:")
            levels = ['kingdom', 'superclass', 'class', 'subclass', 'direct_parent']
            for level in levels:
                if level in classification_results:
                    name = classification_results[level].get('name', 'N/A')
                    print(f"{level.title()}: {name}")
        else:
            print(f"Failed to retrieve results. Status: {result_response.status_code}")
    else:
        print("Query did not complete successfully")
else:
    print(f"Failed to submit query. Status: {response.status_code}")
```

### Output Formats

#### JSON Format (Default)
The JSON format provides the most comprehensive classification information:

```json
{
  "smiles": "CCO",
  "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
  "kingdom": {
    "name": "Organic compounds",
    "description": "Compounds that contain at least one carbon atom...",
    "chemont_id": "CHEMONTID:0000000"
  },
  "superclass": {
    "name": "Organic oxygen compounds",
    "description": "Organic compounds containing oxygen...",
    "chemont_id": "CHEMONTID:0000000"
  },
  "class": {
    "name": "Organooxygen compounds",
    "description": "Organic compounds containing oxygen...",
    "chemont_id": "CHEMONTID:0000000"
  },
  "subclass": {
    "name": "Alcohols and polyols",
    "description": "Compounds in which a hydroxy group...",
    "chemont_id": "CHEMONTID:0000000"
  },
  "direct_parent": {
    "name": "Primary alcohols",
    "description": "Compounds comprising the primary alcohol...",
    "chemont_id": "CHEMONTID:0000000"
  }
}
```

#### SDF Format
Returns molecular structure data in SDF format with embedded classification annotations.

#### CSV Format
Returns tabular data suitable for spreadsheet applications.

### Classification Hierarchy

ClassyFire organizes compounds into a hierarchical classification system:

1. **Kingdom**: Broadest level (e.g., "Organic compounds")
2. **Superclass**: Major structural categories (e.g., "Organic oxygen compounds")
3. **Class**: Functional group families (e.g., "Organooxygen compounds")
4. **Subclass**: Specific functional groups (e.g., "Alcohols and polyols")
5. **Direct Parent**: Most specific classification (e.g., "Primary alcohols")

### Supported Input Formats

ClassyFire accepts various molecular structure formats:

- **SMILES**: Simplified molecular input line entry system
- **InChI**: International Chemical Identifier
- **MOL**: Molfile format
- **SDF**: Structure data file format

### Best Practices

#### 1. Query Management
```python
# Use descriptive labels
response = ClassyFireAPI.submit_query(
    "Aspirin_analysis_2024", 
    "CC(=O)OC1=CC=CC=C1C(=O)O"
)

# Store query IDs for later retrieval
query_ids = []
if response.status_code == 200:
    query_id = response.json()['id']
    query_ids.append(query_id)
```

#### 2. Status Monitoring
```python
def wait_for_completion(query_id, max_wait_time=1800):
    """Wait for query completion with timeout"""
    start_time = time.time()
    
    while time.time() - start_time < max_wait_time:
        status_response = ClassyFireAPI.query_status(query_id)
        
        if status_response and status_response.status_code == 200:
            status = status_response.json().get('status')
            
            if status == 'Done':
                return True
            elif status in ['In Queue', 'In Progress']:
                time.sleep(30)
            else:
                print(f"Unexpected status: {status}")
                return False
        else:
            print("Failed to check status")
            return False
    
    print("Timeout waiting for completion")
    return False
```

#### 3. Batch Processing
```python
def classify_compounds(compounds_dict, delay=60):
    """
    Classify multiple compounds with proper delays
    
    Args:
        compounds_dict: {"label": "smiles", ...}
        delay: seconds between submissions
    """
    query_ids = {}
    
    # Submit all queries
    for label, smiles in compounds_dict.items():
        response = ClassyFireAPI.submit_query(label, smiles)
        
        if response.status_code == 200:
            query_id = response.json()['id']
            query_ids[label] = query_id
            print(f"Submitted {label}: {query_id}")
        else:
            print(f"Failed to submit {label}")
        
        time.sleep(delay)  # Rate limiting
    
    return query_ids
```

### Error Handling

#### Common HTTP Status Codes

- **200**: Success
- **201**: Created (query submitted)
- **202**: Accepted (query in progress)
- **400**: Bad request (invalid input)
- **404**: Not found (invalid query ID)
- **422**: Unprocessable entity (invalid structure)
- **500**: Internal server error
- **503**: Service unavailable

#### Error Handling Example

```python
def safe_submit_query(label, structure):
    """Submit query with comprehensive error handling"""
    try:
        response = ClassyFireAPI.submit_query(label, structure)
        
        if response.status_code == 200:
            return response.json()['id']
        elif response.status_code == 400:
            print(f"Bad request - check input format: {structure}")
        elif response.status_code == 422:
            print(f"Invalid structure: {structure}")
        elif response.status_code >= 500:
            print("Server error - try again later")
        else:
            print(f"Unexpected status: {response.status_code}")
            
    except requests.exceptions.Timeout:
        print("Request timed out")
    except requests.exceptions.ConnectionError:
        print("Connection error")
    except Exception as e:
        print(f"Unexpected error: {e}")
    
    return None
```

### Integration with Other APIs

```python
from provesid.classyfire import ClassyFireAPI
from provesid.pubchem import PubChemAPI
from provesid.opsin import OPSIN

def complete_compound_analysis(compound_name):
    """Complete analysis using multiple APIs"""
    
    # Step 1: Get structure from IUPAC name
    opsin = OPSIN()
    structure_result = opsin.get_id(compound_name)
    
    if structure_result['status'] != 'SUCCESS':
        print(f"Could not convert {compound_name} to structure")
        return None
    
    smiles = structure_result['smiles']
    
    # Step 2: Get PubChem data
    pubchem = PubChemAPI()
    pubchem_data = pubchem.get_compound_by_smiles(smiles)
    
    # Step 3: Get ClassyFire classification
    response = ClassyFireAPI.submit_query(
        f"{compound_name}_analysis", 
        smiles
    )
    
    if response.status_code == 200:
        query_id = response.json()['id']
        
        # Wait and retrieve classification
        if wait_for_completion(query_id):
            classification = ClassyFireAPI.get_query(query_id)
            
            return {
                'name': compound_name,
                'smiles': smiles,
                'pubchem_data': pubchem_data,
                'classification': classification.json()
            }
    
    return None
```

### Limitations and Considerations

#### Processing Time
- Queries can take several minutes to hours to complete
- Complex structures may require longer processing
- Server load affects processing time

#### Rate Limits
- No official rate limits documented
- Recommended to space submissions by 30-60 seconds
- Monitor server response for rate limiting

#### Structure Requirements
- Must be valid chemical structures
- Some exotic or very large structures may not be supported
- Inorganic compounds may have limited classification

#### Service Availability
- Web-based service with potential downtime
- No guaranteed uptime or SLA
- Always implement error handling and retries
