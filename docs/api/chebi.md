# ChEBI API

The ChEBI (Chemical Entities of Biological Interest) class provides access to the ChEBI database, which is a freely available dictionary of molecular entities focused on 'small' chemical compounds.

## Overview

ChEBI is a database and ontology of chemical entities of biological interest. It includes molecular entities that are either products of nature or synthetic products used to intervene in the processes of living organisms.

## Basic Usage

```python
from provesid import ChEBI, get_chebi_entity, search_chebi

# Initialize the ChEBI client
chebi = ChEBI()

# Get complete entity information
water = chebi.get_complete_entity(15377)  # ChEBI:15377 is water
print(f"Name: {water['chebiAsciiName']}")
print(f"Formula: {water['formulaConnectivity']}")
print(f"Mass: {water['mass']}")

# Search by name
results = chebi.search_by_name("aspirin")
for result in results[:3]:
    print(f"{result['chebiId']}: {result['chebiAsciiName']}")
```

## ChEBI Class

### Constructor

```python
ChEBI(timeout=30)
```

**Parameters:**
- `timeout` (int): Request timeout in seconds (default: 30)

### Methods

#### get_complete_entity(chebi_id)

Get complete entity information for a ChEBI ID.

**Parameters:**
- `chebi_id` (Union[int, str]): ChEBI ID (with or without 'CHEBI:' prefix)

**Returns:**
- `dict`: Complete entity information, None if not found

**Example:**
```python
chebi = ChEBI()
entity = chebi.get_complete_entity(15377)
if entity:
    print(f"Name: {entity['chebiAsciiName']}")
    print(f"Definition: {entity['definition']}")
    print(f"Formula: {entity['formulaConnectivity']}")
    print(f"Mass: {entity['mass']}")
```

#### get_lite_entity(chebi_id)

Get basic entity information for a ChEBI ID (lightweight version).

**Parameters:**
- `chebi_id` (Union[int, str]): ChEBI ID (with or without 'CHEBI:' prefix)

**Returns:**
- `dict`: Basic entity information, None if not found

**Example:**
```python
chebi = ChEBI()
entity = chebi.get_lite_entity("CHEBI:16236")  # ethanol
if entity:
    print(f"Name: {entity['chebiAsciiName']}")
    print(f"Search Score: {entity['searchScore']}")
    print(f"Entity Star: {entity['entityStar']}")
```

#### search_by_name(search_text, search_category="ALL", max_results=50, stars="ALL")

Search ChEBI database by compound name.

**Parameters:**
- `search_text` (str): Text to search for
- `search_category` (str): Search category ('ALL', 'CHEBI_NAME', 'DEFINITION', etc.)
- `max_results` (int): Maximum number of results to return
- `stars` (str): Star category ('ALL', 'TWO_ONLY', 'THREE_ONLY')

**Returns:**
- `list`: List of matching entities

**Example:**
```python
chebi = ChEBI()
results = chebi.search_by_name("caffeine", max_results=10)
for result in results:
    print(f"{result['chebiId']}: {result['chebiAsciiName']}")
```

#### get_structure(chebi_id, structure_type="mol")

Get chemical structure for a ChEBI ID.

**Parameters:**
- `chebi_id` (Union[int, str]): ChEBI ID (with or without 'CHEBI:' prefix)
- `structure_type` (str): Structure format ('mol', 'sdf', 'smiles', 'inchi')

**Returns:**
- `str`: Chemical structure in requested format, None if not found

**Example:**
```python
chebi = ChEBI()
smiles = chebi.get_structure(15377, "smiles")  # water
mol_file = chebi.get_structure(15377, "mol")
inchi = chebi.get_structure(15377, "inchi")

print(f"SMILES: {smiles}")
print(f"InChI: {inchi}")
```

#### get_ontology_parents(chebi_id)

Get ontology parents for a ChEBI ID.

**Parameters:**
- `chebi_id` (Union[int, str]): ChEBI ID (with or without 'CHEBI:' prefix)

**Returns:**
- `list`: List of parent entities in the ontology

**Example:**
```python
chebi = ChEBI()
parents = chebi.get_ontology_parents(15377)  # water
for parent in parents:
    print(f"Parent: {parent['chebiId']} - {parent['chebiName']}")
    print(f"Relationship: {parent['type']}")
```

#### get_ontology_children(chebi_id)

Get ontology children for a ChEBI ID.

**Parameters:**
- `chebi_id` (Union[int, str]): ChEBI ID (with or without 'CHEBI:' prefix)

**Returns:**
- `list`: List of child entities in the ontology

**Example:**
```python
chebi = ChEBI()
children = chebi.get_ontology_children(24431)  # chemical entity
for child in children[:5]:  # Show first 5
    print(f"Child: {child['chebiId']} - {child['chebiName']}")
```

#### batch_get_entities(chebi_ids, pause_time=0.1)

Get complete entity information for multiple ChEBI IDs.

**Parameters:**
- `chebi_ids` (List[Union[int, str]]): List of ChEBI IDs
- `pause_time` (float): Pause between requests to be respectful to the API

**Returns:**
- `dict`: Dictionary mapping ChEBI IDs to entity information

**Example:**
```python
chebi = ChEBI()
ids = [15377, 16236, 27732]  # water, ethanol, caffeine
results = chebi.batch_get_entities(ids)

for chebi_id, data in results.items():
    print(f"{chebi_id}: {data['chebiAsciiName']}")
```

## Convenience Functions

### get_chebi_entity(chebi_id)

Quick function to get ChEBI entity information.

**Parameters:**
- `chebi_id` (Union[int, str]): ChEBI ID (with or without 'CHEBI:' prefix)

**Returns:**
- `dict`: Entity information, None if not found

**Example:**
```python
from provesid import get_chebi_entity

water = get_chebi_entity(15377)
if water:
    print(f"Name: {water['chebiAsciiName']}")
    print(f"Formula: {water['formulaConnectivity']}")
```

### search_chebi(search_text, max_results=10)

Quick function to search ChEBI by name.

**Parameters:**
- `search_text` (str): Text to search for
- `max_results` (int): Maximum number of results to return

**Returns:**
- `list`: List of matching entities

**Example:**
```python
from provesid import search_chebi

results = search_chebi("glucose")
for result in results[:3]:
    print(f"{result['chebiId']}: {result['chebiAsciiName']}")
```

## Common ChEBI IDs

Here are some commonly used ChEBI IDs:

- **CHEBI:15377** - water
- **CHEBI:16236** - ethanol  
- **CHEBI:27732** - caffeine
- **CHEBI:15996** - GTP
- **CHEBI:15422** - ATP
- **CHEBI:17234** - glucose
- **CHEBI:16467** - cholesterol
- **CHEBI:15365** - aspirin

## Error Handling

The ChEBI class includes comprehensive error handling:

```python
from provesid import ChEBI, ChEBIError

chebi = ChEBI()

try:
    entity = chebi.get_complete_entity(15377)
    if entity is None:
        print("Entity not found")
    else:
        print(f"Found: {entity['chebiAsciiName']}")
except ChEBIError as e:
    print(f"ChEBI API error: {e}")
```

## Data Structures

### Complete Entity Response

A complete entity response includes:

```python
{
    'chebiId': 'CHEBI:15377',
    'chebiAsciiName': 'water',
    'definition': 'An oxygen hydride consisting of an oxygen atom...',
    'formulaConnectivity': 'H2O',
    'mass': '18.01056',
    'monoisotopicMass': '18.01056',
    'charge': '0',
    'synonyms': [...],  # List of synonyms
    'iupacNames': [...],  # List of IUPAC names
    'databaseLinks': [...],  # External database links
    'chemicalStructures': [...],  # Chemical structures
    'registryNumbers': [...]  # Registry numbers
}
```

### Lite Entity Response

A lite entity response includes basic information:

```python
{
    'chebiId': 'CHEBI:15377',
    'chebiAsciiName': 'water',
    'searchScore': '1.0',
    'entityStar': '3'
}
```

### Ontology Relationship

Ontology relationships include:

```python
{
    'chebiId': 'CHEBI:24431',
    'chebiName': 'chemical entity',
    'type': 'is_a',
    'status': 'C'
}
```

## Performance Tips

1. **Use lite entities** when you only need basic information
2. **Batch requests** when getting multiple entities
3. **Add pause time** in batch operations to be respectful to the API
4. **Cache results** for frequently accessed entities
5. **Use specific search categories** to improve search accuracy

## Integration Example

```python
from provesid import ChEBI
import pandas as pd

# Initialize ChEBI client
chebi = ChEBI()

# List of compounds to look up
compound_names = ["water", "ethanol", "glucose", "caffeine"]

# Search and collect data
results = []
for name in compound_names:
    search_results = chebi.search_by_name(name, max_results=1)
    if search_results:
        chebi_id = search_results[0]['chebiId']
        complete_entity = chebi.get_complete_entity(chebi_id)
        if complete_entity:
            results.append({
                'search_name': name,
                'chebi_id': chebi_id,
                'name': complete_entity['chebiAsciiName'],
                'formula': complete_entity.get('formulaConnectivity', ''),
                'mass': complete_entity.get('mass', ''),
                'definition': complete_entity.get('definition', '')[:100] + '...'
            })

# Create DataFrame
df = pd.DataFrame(results)
print(df)
```

## Links

- [ChEBI Database](https://www.ebi.ac.uk/chebi/)
- [ChEBI Web Services](https://www.ebi.ac.uk/chebi/webServices.do)
- [ChEBI Ontology](https://www.ebi.ac.uk/ols/ontologies/chebi)
