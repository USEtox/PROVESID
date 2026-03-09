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

# ChEBI (Chemical Entities of Biological Interest) Tutorial

ChEBI is a freely available dictionary of molecular entities focused on 'small' chemical compounds. This tutorial demonstrates how to use the `ChEBI` class and convenience functions from the `provesid` package to access chemical information from the ChEBI database.

**Important Note**: The ChEBI search API is currently experiencing intermittent issues. This tutorial shows both the intended search functionality and alternative approaches using direct ChEBI ID lookups when search is unavailable.

ChEBI provides comprehensive chemical information including:
- Chemical structures (SMILES, InChI, MOL)
- Ontological relationships
- Biological roles and activities
- Cross-references to other databases
- Detailed chemical properties

The database is maintained by the European Bioinformatics Institute (EBI) and contains over 185,000 chemical entities.

```{code-cell} ipython3
from provesid import ChEBI, ChEBIError, get_chebi_entity, search_chebi

chebi = ChEBI()
print("ChEBI initialized successfully!")
print(f"Base URL: {chebi.base_url}")
print(f"Timeout: {chebi.timeout} seconds")
```

## 1. Getting Complete Entity Information

The primary method for retrieving detailed chemical information is `get_complete_entity()`. Let's look up some common compounds:

```{code-cell} ipython3
# Get complete information for water (ChEBI:15377)
water = chebi.get_complete_entity(15377)
if water:
    print("Water (CHEBI:15377):")
    print(f"  ASCII Name: {water.get('chebiAsciiName')}")
    print(f"  IUPAC Name: {water.get('iupacName')}")
    print(f"  Definition: {water.get('definition')}")
    # Access chemical formula from the Formulae dict
    formulae_data = water.get('Formulae', {})
    formula = formulae_data.get('data', 'N/A') if formulae_data else 'N/A'
    print(f"  Molecular Formula: {formula}")
    print(f"  SMILES: {water.get('smiles')}")
    print(f"  InChI: {water.get('inchi')}")
    print(f"  InChI Key: {water.get('inchiKey')}")
    print(f"  Mass: {water.get('mass')}")
    print(f"  Charge: {water.get('charge')}")
else:
    print("Water information not found")
```

```{code-cell} ipython3
# Get information for aspirin (ChEBI:15365)
aspirin = chebi.get_complete_entity(15365)
if aspirin:
    print("Aspirin (CHEBI:15365):")
    print(f"  ASCII Name: {aspirin.get('chebiAsciiName')}")
    print(f"  Definition: {aspirin.get('definition')}")
    print(f"  SMILES: {aspirin.get('smiles')}")
    print(f"  InChI Key: {aspirin.get('inchiKey')}")
    print(f"  Mass: {aspirin.get('mass')}")
    
    # Show synonyms if available
    synonyms = aspirin.get('Synonyms', [])
    if synonyms:
        print(f"  Number of synonyms: {len(synonyms)}")
        print(f"  First 3 synonyms: {[syn.get('data') for syn in synonyms[:3]]}")
else:
    print("Aspirin information not found")
```

## 2. Using Convenience Functions

The package provides convenient functions for quick lookups without creating a ChEBI instance:

```{code-cell} ipython3
# Using the convenience function get_chebi_entity
caffeine = get_chebi_entity(27732)  # Caffeine
if caffeine:
    print("Caffeine (using convenience function):")
    print(f"  Name: {caffeine.get('chebiAsciiName')}")
    print(f"  Definition: {caffeine.get('definition')}")
    print(f"  SMILES: {caffeine.get('smiles')}")
    # Access chemical formula correctly
    formulae_data = caffeine.get('Formulae', {})
    formula = formulae_data.get('data', 'N/A') if formulae_data else 'N/A'
    print(f"  Molecular formula: {formula}")
else:
    print("Caffeine not found")
```

## Important: Accessing Chemical Formula Data

ChEBI returns chemical formula data in a specific structure. The formula is stored in the `Formulae` key (with capital F) as a dictionary with `data` and `source` fields. Here's the correct way to access it:

```{code-cell} ipython3
# Demonstrate the correct way to access chemical formula data
aspirin = get_chebi_entity(15365)  # Aspirin
if aspirin:
    print("Correct way to access chemical formula:")
    print(f"  Compound: {aspirin.get('chebiAsciiName')}")
    
    # Method 1: Safe access with get()
    formulae_data = aspirin.get('Formulae', {})
    if formulae_data:
        formula = formulae_data.get('data', 'N/A')
        source = formulae_data.get('source', 'N/A')
        print(f"  Formula: {formula} (source: {source})")
    else:
        print("  Formula: Not available")
    
    # Method 2: One-liner (also safe)
    formula_oneliner = aspirin.get('Formulae', {}).get('data', 'N/A')
    print(f"  Formula (one-liner): {formula_oneliner}")
    
    print()
    print("❌ WRONG way (this would cause errors):")
    print("  # aspirin['formulae'][0]['data']  # Wrong - lowercase 'formulae' doesn't exist")
    print("  # aspirin['Formulae']['data']     # Wrong - no error checking")
    
    print()
    print("✅ CORRECT way:")
    print("  formulae_data = aspirin.get('Formulae', {})")
    print("  formula = formulae_data.get('data', 'N/A') if formulae_data else 'N/A'")
```

```{code-cell} ipython3
# Using the search convenience function - Note: ChEBI search API may have intermittent issues
print("Searching for 'glucose':")
print("(Note: ChEBI search API is currently experiencing issues)")
glucose_results = search_chebi("glucose", max_results=5)
if glucose_results:
    for i, result in enumerate(glucose_results[:3], 1):
        print(f"  {i}. {result.get('chebiId')}: {result.get('chebiAsciiName')} - {result.get('definition', 'No definition')[:100]}...")
else:
    print("  Search returned no results (API may be temporarily unavailable)")
    print("  Alternative: Use get_chebi_entity() with known ChEBI IDs")
    print("  For example, glucose is ChEBI:17234")
    
    # Show alternative approach
    glucose = get_chebi_entity(17234)  # D-glucose
    if glucose:
        print(f"  Direct lookup - ChEBI:17234: {glucose.get('chebiAsciiName')}")
```

## 3. Searching by Name

ChEBI provides powerful search capabilities to find compounds by name, synonym, or definition:

```{code-cell} ipython3
# Search for compounds containing "ethanol" - Note: API may have issues
ethanol_results = chebi.search_by_name("ethanol", size=10)
print(f"Found {len(ethanol_results)} results for 'ethanol':")

if ethanol_results:
    for i, result in enumerate(ethanol_results[:5], 1):
        print(f"  {i}. {result.get('chebiId')}: {result.get('chebiAsciiName')}")
        print(f"     Definition: {result.get('definition', 'No definition')[:80]}...")
        print()
else:
    print("  Search API temporarily unavailable. Using direct lookup instead:")
    # Alternative: direct lookup for ethanol (ChEBI:16236)
    ethanol = get_chebi_entity(16236)
    if ethanol:
        print(f"  Direct lookup - ChEBI:16236: {ethanol.get('chebiAsciiName')}")
        print(f"     Definition: {ethanol.get('definition', 'No definition')[:80]}...")
```

```{code-cell} ipython3
# Search for vitamin compounds - Note: API may have issues
vitamin_results = chebi.search_by_name("vitamin", size=8)
print(f"Found {len(vitamin_results)} results for 'vitamin':")

if vitamin_results:
    for result in vitamin_results[:5]:
        print(f"  • {result.get('chebiId')}: {result.get('chebiAsciiName')}")
else:
    print("  Search API temporarily unavailable. Using known vitamin ChEBI IDs:")
    vitamin_ids = [29073, 18405, 17015]  # Vitamin C, Vitamin E, Vitamin D
    for vid in vitamin_ids:
        vitamin = get_chebi_entity(vid)
        if vitamin:
            print(f"  • CHEBI:{vid}: {vitamin.get('chebiAsciiName')}")
```

## 4. Getting Chemical Structures

ChEBI can provide chemical structures in various formats:

```{code-cell} ipython3
# Get structure-related fields for ethanol (ChEBI:16236)
ethanol_id = 16236

print("Ethanol structure in different formats:")

ethanol = chebi.get_complete_entity(ethanol_id)
if ethanol:
    print(f"  SMILES: {ethanol.get('smiles')}")
    print(f"  InChI: {ethanol.get('inchi')}")

# Retrieve MOL file content from dedicated endpoint
mol_structure = chebi.get_molfile(ethanol_id)
if mol_structure:
    mol_lines = mol_structure.split('\n')[:5]
    print("  MOL format (first 5 lines):")
    for line in mol_lines:
        print(f"    {line}")
```

## 5. Ontological Relationships

ChEBI organizes compounds in an ontological hierarchy. You can explore parent-child relationships:

```{code-cell} ipython3
# Get ontology parents for ethanol
ethanol_parents = chebi.get_ontology_parents(16236)
print("Ethanol ontology parents:")

parent_items = ethanol_parents if isinstance(ethanol_parents, list) else []
if isinstance(ethanol_parents, dict):
    for key in ("items", "results", "data"):
        if isinstance(ethanol_parents.get(key), list):
            parent_items = ethanol_parents[key]
            break

for parent in parent_items[:5]:
    print(f"  • {parent.get('chebiId')}: {parent.get('chebiName')} ({parent.get('type')})")
```

```{code-cell} ipython3
# Get ontology children for alcohols (ChEBI:30879)
alcohol_children = chebi.get_ontology_children(30879)
child_items = alcohol_children if isinstance(alcohol_children, list) else []
if isinstance(alcohol_children, dict):
    for key in ("items", "results", "data"):
        if isinstance(alcohol_children.get(key), list):
            child_items = alcohol_children[key]
            break

print(f"Found {len(child_items)} children for 'alcohol' (first 5):")
for child in child_items[:5]:
    print(f"  • {child.get('chebiId')}: {child.get('chebiName')} ({child.get('type')})")
```

## 6. Batch Processing

For multiple compounds, use batch processing with built-in rate limiting:

```{code-cell} ipython3
# Process multiple ChEBI IDs at once
compound_ids = [15377, 16236, 15365, 27732, 17234]  # water, ethanol, aspirin, caffeine, glucose
compound_names = ["water", "ethanol", "aspirin", "caffeine", "glucose"]

print("Batch processing multiple compounds:")
batch_results = chebi.batch_get_entities(compound_ids, pause_time=0.2)

for i, (chebi_id, name) in enumerate(zip([f"CHEBI:{id}" for id in compound_ids], compound_names)):
    if chebi_id in batch_results:
        compound = batch_results[chebi_id]
        print(f"  {i+1}. {name} ({chebi_id}):")
        print(f"     Name: {compound.get('chebiAsciiName')}")
        # Access formula correctly from Formulae dict
        formulae_data = compound.get('Formulae', {})
        formula = formulae_data.get('data', 'N/A') if formulae_data else 'N/A'
        print(f"     Formula: {formula}")
        print(f"     Mass: {compound.get('mass')}")
    else:
        print(f"  {i+1}. {name} ({chebi_id}): Not found")
    print()
```

## 7. Error Handling

ChEBI provides robust error handling for various scenarios:

```{code-cell} ipython3
# Try to get information for invalid ChEBI IDs
print("Testing error handling:")

# Invalid ChEBI ID
invalid_result = chebi.get_complete_entity(999999999)
print(f"Invalid ID (999999999): {invalid_result}")

# Non-existent compound search
empty_search = chebi.search_by_name("thiscompounddoesnotexist12345")
print(f"Empty search results: {len(empty_search)} results")

# Handle ChEBIError exceptions
try:
    # This might cause a timeout or network error
    chebi_timeout = ChEBI(timeout=0.001)  # Very short timeout
    result = chebi_timeout.get_complete_entity(15377)
except ChEBIError as e:
    print(f"ChEBIError caught: {e}")
except Exception as e:
    print(f"Other error: {e}")
```

## 8. Exploring Compound Details

Let's explore the comprehensive information available for a complex biological molecule:

```{code-cell} ipython3
# Remove the problematic line that tries to access vitamin_c["Formulae"]["data"] directly

# This was causing an error since we should use .get() for safe access
vitamin_c = chebi.get_complete_entity(29073)

if vitamin_c:
    print("Vitamin C (CHEBI:29073) - Detailed Information:")
    print(f"  ASCII Name: {vitamin_c.get('chebiAsciiName')}")
    print(f"  IUPAC Name: {vitamin_c.get('iupacName')}")
    print(f"  Definition: {vitamin_c.get('definition')}")
    print(f"  SMILES: {vitamin_c.get('smiles')}")
    print(f"  Mass: {vitamin_c.get('mass')}")
    print(f"  Charge: {vitamin_c.get('charge')}")
    
    # Explore formulas
    formulae = vitamin_c.get('Formulae', {})
    if formulae:
        print(f"  Chemical Formula:")
        print(f"    • {formulae.get('data')} (source: {formulae.get('source')})")
    
    # Explore synonyms
    synonyms = vitamin_c.get('Synonyms', [])
    if synonyms:
        print(f"  Synonyms ({len(synonyms)} total, showing first 5):")
        for syn in synonyms[:5]:
            print(f"    • {syn.get('data')} ({syn.get('source')})")
    
    # Explore database links
    db_links = vitamin_c.get('DatabaseLinks', [])
    if db_links:
        print(f"  Database Links ({len(db_links)} total, showing first 5):")
        for link in db_links[:5]:
            print(f"    • {link.get('type')}: {link.get('data')}")
```

```{code-cell} ipython3
print(vitamin_c.get("Formulae", {}).get("data", "N/A"))
```

## 9. Practical Applications

Here are some practical use cases for the ChEBI API:

```{code-cell} ipython3
# Use case 1: Get basic chemical identifiers for a list of compounds
def get_chemical_identifiers(chebi_ids):
    """Get basic chemical identifiers for multiple compounds"""
    results = []
    for chebi_id in chebi_ids:
        compound = get_chebi_entity(chebi_id)
        if compound:
            results.append({
                'chebi_id': f"CHEBI:{chebi_id}",
                'name': compound.get('chebiAsciiName'),
                'smiles': compound.get('smiles'),
                'inchi_key': compound.get('inchiKey'),
                'mass': compound.get('mass'),
                # Access formula correctly from Formulae dict
                'formula': compound.get('Formulae', {}).get('data') if compound.get('Formulae') else None
            })
        else:
            results.append({
                'chebi_id': f"CHEBI:{chebi_id}",
                'error': 'Not found'
            })
    return results

# Test with common metabolites
metabolite_ids = [15377, 16236, 17234, 15422, 16526]  # water, ethanol, glucose, adenosine triphosphate, carbon dioxide
metabolite_data = get_chemical_identifiers(metabolite_ids)

print("Chemical identifiers for common metabolites:")
for data in metabolite_data:
    if 'error' not in data:
        print(f"  {data['chebi_id']}: {data['name']}")
        print(f"    Formula: {data['formula']}, Mass: {data['mass']}")
        print(f"    SMILES: {data['smiles']}")
    else:
        print(f"  {data['chebi_id']}: {data['error']}")
    print()
```

```{code-cell} ipython3
# Use case 2: Find compounds by biological role
def find_compounds_by_role(search_term, max_results=10):
    """Find compounds related to a biological role or function"""
    # Note: Search API may have issues, so we'll show an alternative approach
    print(f"Note: ChEBI search API is currently having issues.")
    print(f"For demonstration, showing some known compounds related to '{search_term}':")
    
    # Example compounds for different search terms
    known_compounds = {
        'antioxidant': [29073, 16236, 27732],  # Vitamin C, ethanol, caffeine
        'vitamin': [29073, 18405, 17015],      # Vitamin C, E, D
        'hormone': [15365, 27732],             # Example compounds
    }
    
    compound_ids = known_compounds.get(search_term.lower(), [29073, 16236])  # Default examples
    compounds = []
    
    for compound_id in compound_ids[:max_results]:
        detailed = get_chebi_entity(compound_id)
        if detailed:
            compounds.append({
                'chebi_id': f"CHEBI:{compound_id}",
                'name': detailed.get('chebiAsciiName'),
                'definition': detailed.get('definition'),
                'smiles': detailed.get('smiles'),
                'mass': detailed.get('mass')
            })
    
    return compounds

# Search for antioxidants
antioxidants = find_compounds_by_role("antioxidant", max_results=5)
print("Compounds related to 'antioxidant':")
for compound in antioxidants[:3]:
    print(f"  • {compound['chebi_id']}: {compound['name']}")
    print(f"    Definition: {compound['definition'][:100]}...")
    print(f"    SMILES: {compound['smiles']}")
    print()
```

```{code-cell} ipython3
# Use case 3: Build a compound database with cross-references
def build_compound_database(chebi_ids):
    """Build a comprehensive compound database with cross-references"""
    database = {}
    
    for chebi_id in chebi_ids:
        compound = get_chebi_entity(chebi_id)
        if compound:
            # Extract cross-references
            db_links = compound.get('DatabaseLinks', [])
            cross_refs = {}
            for link in db_links:
                db_type = link.get('type', 'Unknown')
                if db_type not in cross_refs:
                    cross_refs[db_type] = []
                cross_refs[db_type].append(link.get('data'))
            
            database[f"CHEBI:{chebi_id}"] = {
                'name': compound.get('chebiAsciiName'),
                'iupac_name': compound.get('iupacName'),
                'definition': compound.get('definition'),
                'smiles': compound.get('smiles'),
                'inchi_key': compound.get('inchiKey'),
                'mass': compound.get('mass'),
                'charge': compound.get('charge'),
                'cross_references': cross_refs
            }
    
    return database

# Build database for some pharmaceutical compounds
pharma_ids = [15365, 27732, 3002]  # aspirin, caffeine, morphine
pharma_db = build_compound_database(pharma_ids)

print("Pharmaceutical compound database:")
for chebi_id, data in pharma_db.items():
    print(f"\n{chebi_id}: {data['name']}")
    print(f"  IUPAC: {data['iupac_name']}")
    print(f"  Mass: {data['mass']}")
    print(f"  Cross-references:")
    for db_name, refs in data['cross_references'].items():
        if refs:  # Only show non-empty references
            print(f"    {db_name}: {refs[:2]}...")  # Show first 2 references
```

## Summary

The `ChEBI` class and convenience functions provide comprehensive access to the ChEBI database:

### Main ChEBI Class Methods:
1. **`get_complete_entity(chebi_id)`**: Get detailed information for a ChEBI ID
2. **`get_lite_entity(chebi_id)`**: Get basic information only
3. **`search_by_name(search_text)`**: Search compounds by name
4. **`get_complete_entity(chebi_id)` + `get_molfile(compound_id)`**: Get structure fields and MOL data
5. **`get_ontology_parents(chebi_id)`**: Get parent entities in ontology
6. **`get_ontology_children(chebi_id)`**: Get child entities in ontology
7. **`batch_get_entities(chebi_ids)`**: Process multiple IDs efficiently

### Convenience Functions:
- **`get_chebi_entity(chebi_id)`**: Quick entity lookup
- **`search_chebi(search_text)`**: Quick search functionality

### Key Features:
- ✅ **Comprehensive Data**: Names, structures, properties, ontology, cross-references
- ✅ **Multiple Formats**: SMILES, InChI, MOL files for structures
- ✅ **Ontological Browsing**: Navigate parent-child relationships
- ✅ **Batch Processing**: Efficient handling of multiple compounds
- ✅ **Error Handling**: Robust error management with custom exceptions
- ✅ **Rate Limiting**: Built-in delays for respectful API usage
- ✅ **Free Access**: No API key required

### Returned Data Includes:
- Chemical identifiers (name, IUPAC name, synonyms)
- Molecular structures (SMILES, InChI, InChI Key)
- Physical properties (mass, charge, molecular formula)
- Biological information (definition, role, function)
- Database cross-references (PubChem, UniProt, KEGG, etc.)
- Ontological relationships (parents, children, classifications)

### Best Use Cases:
- Chemical database integration
- Biological pathway analysis
- Drug discovery research
- Metabolomics studies
- Chemical ontology exploration
- Cross-database linking

ChEBI is particularly valuable for researchers working with biological systems, as it focuses on chemical entities relevant to biological processes and provides rich ontological context for understanding chemical relationships.
