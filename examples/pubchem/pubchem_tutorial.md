---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.1
kernelspec:
  display_name: mychem
  language: python
  name: python3
---

# PubChem API tutorial
This tutorial explains some use cases of the PubChem API. Note that we also have another API in this package called `pubchemview` that has its separate tutorial.

```{code-cell} ipython3
# Import the required modules for PubChem API
from provesid.pubchem import PubChemAPI, Domain, CompoundProperties
import json # mostly for nicer printing :-)

# Initialize the PubChem API client
pc = PubChemAPI()
```

# CID, SID, and AID
CID (Compound ID), SID (Substance ID), and AID (Assay ID) are unique identifiers used by PubChem:

- **CID (Compound ID):** Identifies a unique chemical structure in the PubChem Compound database. Each distinct molecule has a single CID, regardless of how it was submitted or by whom. Example: formaldehyde has CID 712.

- **SID (Substance ID):** Identifies a record in the PubChem Substance database, which represents a substance as submitted by a depositor. Multiple SIDs can map to the same CID if different sources submit the same compound. Example: formaldehyde may have many SIDs from different submitters.

- **AID (Assay ID):** Identifies a bioassay record in the PubChem BioAssay database. Each AID corresponds to a specific biological test or experiment, which may reference one or more CIDs or SIDs.

In summary:  
- **CID** = unique chemical structure  
- **SID** = depositor-submitted sample/record  
- **AID** = bioassay/experiment

Whatever we need to retrieve from PubChem, we first need to look for these IDs.

+++

# How to search for IDs?
The next code cell demonstrates how to use the `PubChemAPI` to search for CIDs (Compound IDs) and SIDs (Substance IDs) in PubChem by querying with different types of identifiers such as chemical names, SMILES strings, or CAS numbers.

- `pc.get_cids_by_name('aspirin')` looks up CIDs by the compound name "aspirin".
- `pc.get_cids_by_name('water', domain=Domain.COMPOUND)` searches for CIDs in the compound domain using the name "water".
- `pc.get_cids_by_name('8000-78-0', domain=Domain.SUBSTANCE)` searches for CIDs in the substance domain using a CAS number.
- `pc.find_cids_comprehensive('8000-78-0')` performs a comprehensive search across both compound and substance domains for the given CAS number, returning the found cid numbers and how they are found.
- `pc.get_sids_by_name('8000-78-0')` retrieves SIDs by searching with a CAS number.

You can use these methods to retrieve PubChem IDs by providing any identifier (name, SMILES, CAS, etc.). The API will return the corresponding CIDs or SIDs, making it easy to map between different chemical identifiers and PubChem records. This is useful for integrating chemical data from various sources or for further property lookups in PubChem.

```{code-cell} ipython3
# 1. Default behavior (backward compatible)
cids_aspirin = pc.get_cids_by_name('aspirin')  # Returns clean list of CIDs
print(f"cids found for the name aspirin: {cids_aspirin}")

# 2. Explicit compound domain
cids_water = pc.get_cids_by_name('water', domain=Domain.COMPOUND)
print(f"cids found for the name water: {cids_water}")

# 3. Search in substance domain (new capability)
cids_garlic_oil = pc.get_cids_by_name('8000-78-0', domain=Domain.SUBSTANCE)  # [6850738]
print(f"cid found for the CAS number 8000-78-0: {cids_garlic_oil}")

# 4. Comprehensive search across both domains
results = pc.find_cids_comprehensive('8000-78-0')
# Returns detailed results with recommendations
print(f"comprehensive search results for the CAS number 8000-78-0: {results}")

# 5. Enhanced SID search
sids = pc.get_sids_by_name('8000-78-0')  # Returns clean list of SIDs
print(f"sids found for the CAS number 8000-78-0: {sids}")
```

```{code-cell} ipython3
pc.find_cids_comprehensive('8000-78-0')
```

One of the main use cases for me was to look up a compound by its CAS number and if nothing is found look up a substance by the same CAS number, especially for those that are not found in CAS Common Registry.

+++

# Synonyms and specific IDs
After finding the `cid` using "one" identifier, we can obtain a list of synonyms (e.g. chemical name, CAS number, etc.) and also extract certain identifiers from the list: 

```{code-cell} ipython3
synonyms = pc.get_compound_synonyms(cids_aspirin[0])
ids = pc.get_compound_identifiers(cids_aspirin[0])
print(synonyms[:5])
print(ids)
```

# Compound property
After obtaining the `cid` of a compound, we can obtain [compound properties](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest#section=Compound-Property-Tables) by calling one of the following functions, that can give basic, selected, or all available properties.  

```{code-cell} ipython3
res_basic = pc.get_basic_compound_info(cids_aspirin[0])
# Pretty print the result
print(json.dumps(res_basic, indent=2))
```

```{code-cell} ipython3
res_selected = pc.get_compound_properties(cids_aspirin[0],
                                          [CompoundProperties.SMILES,
                                           CompoundProperties.INCHI,
                                           CompoundProperties.INCHIKEY],
                                          include_synonyms=False)
print(json.dumps(res_selected, indent=2))
```

```{code-cell} ipython3
res_all = pc.get_all_compound_info(cids_aspirin[0])
print(json.dumps(res_all, indent=2))
```

# Substance properties
The substance properties can be found by providing a `sid` to the 

```{code-cell} ipython3
sids = pc.get_sids_by_name("garlic oil")
res = pc.get_substance_by_sid(sids[1])
res
```

# PubChem View

```{code-cell} ipython3
from provesid import PubChemView, get_property_table
```

```{code-cell} ipython3
logp_table = get_property_table(cids_aspirin[0], "LogP")
logp_table
```

```{code-cell} ipython3
pcv = PubChemView()
res_logP = pcv.get_property_summary(cids_aspirin[0], "LogP")
print(json.dumps(res_logP, indent=2))
```

# Advanced PubChem API Features

The PubChem API has been improved to provide more elegant data access. Previously, methods like `get_substance_by_sid()` and `get_compound_by_cid()` returned data wrapped in redundant structures requiring access like `result["PC_Substances"][0]` or `result["PC_Compounds"][0]`. Now these methods automatically extract the relevant data for easier access.

## Batch Processing and Multiple Compounds

Let's explore how to work with multiple compounds and batch processing:

```{code-cell} ipython3
# Batch processing for multiple compounds
compound_names = ["aspirin", "caffeine", "acetaminophen", "ibuprofen"]
all_cids = []

for name in compound_names:
    cids = pc.get_cids_by_name(name)
    if cids:
        all_cids.append(cids[0])  # Take the first CID for each compound
        print(f"{name}: CID {cids[0]}")

print(f"\nCollected CIDs: {all_cids}")

# Batch property retrieval
properties = [CompoundProperties.MOLECULAR_WEIGHT, 
              CompoundProperties.MOLECULAR_FORMULA,
              CompoundProperties.SMILES]

batch_results = pc.get_compound_properties_batch(all_cids, properties)
print("\nBatch property results:")
print(json.dumps(batch_results, indent=2))
```

## Chemical Structure Searching

PubChem API supports searching by various chemical identifiers including SMILES and InChI keys. Both `get_cids_by_smiles()` and `get_cids_by_inchikey()` methods return clean lists of CIDs, and the corresponding `get_compounds_by_*()` methods return the compound data directly without wrapper structures:

```{code-cell} ipython3
# Search by SMILES string (caffeine)
caffeine_smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
cids_by_smiles = pc.get_cids_by_smiles(caffeine_smiles)
print(f"CIDs found by SMILES: {cids_by_smiles}")

# Get compound record by SMILES (new improved method - no wrapper needed!)
compound_by_smiles = pc.get_compounds_by_smiles(caffeine_smiles)
print(f"Compound data type: {type(compound_by_smiles)}")
if isinstance(compound_by_smiles, dict):
    print(f"Direct access to compound keys: {list(compound_by_smiles.keys())}")

# Search by InChI Key
inchikey = "RYYVLZVUVIJVGH-UHFFFAOYSA-N"  # caffeine InChI key
cids_by_inchikey = pc.get_cids_by_inchikey(inchikey)
print(f"CIDs found by InChI Key: {cids_by_inchikey}")

# Get compound record by InChI Key (new improved method)
compound_by_inchikey = pc.get_compounds_by_inchikey(inchikey)
print(f"Compound by InChI Key - type: {type(compound_by_inchikey)}")
```

# Comprehensive PubChem View Tutorial

PubChemView provides access to experimental properties that are not available through the standard PubChem API. These include physical and chemical properties like melting point, boiling point, solubility, vapor pressure, and many others.

## Available Experimental Properties

Let's first explore what experimental properties are available for a compound:

```{code-cell} ipython3
# Check what experimental properties are available for aspirin
available_props = pcv.get_available_properties(cids_aspirin[0])
print(f"Available experimental properties for aspirin ({len(available_props)} total):")
for prop in available_props:
    print(f"  - {prop}")

# Show the standard experimental properties mapping
print(f"\nTotal standard experimental properties supported: {len(pcv.experimental_properties)}")
print("Some examples:")
for i, (key, value) in enumerate(list(pcv.experimental_properties.items())[:10]):
    print(f"  {key} -> {value}")
print("  ...")
```

## Common Physical Properties

Let's extract some common physical and chemical properties using the convenience methods:

```{code-cell} ipython3
# Melting Point
melting_point = pcv.get_melting_point(cids_aspirin[0])
print("Melting Point data:")
for i, mp in enumerate(melting_point[:3]):  # Show first 3 entries
    print(f"  {i+1}: {mp.value} (Ref: {mp.reference_number})")

# Boiling Point  
boiling_point = pcv.get_boiling_point(cids_aspirin[0])
print(f"\nBoiling Point data ({len(boiling_point)} entries):")
for i, bp in enumerate(boiling_point[:2]):
    print(f"  {i+1}: {bp.value}")

# Solubility
solubility = pcv.get_solubility(cids_aspirin[0])
print(f"\nSolubility data ({len(solubility)} entries):")
for i, sol in enumerate(solubility[:3]):
    print(f"  {i+1}: {sol.value}")

# Density
density = pcv.get_density(cids_aspirin[0])
print(f"\nDensity data ({len(density)} entries):")
for i, dens in enumerate(density[:2]):
    print(f"  {i+1}: {dens.value}")
```

Note that the experimental data are not reported homogeneously and therefore it becomes difficult to come up with a single method to extract values, units, and experimental conditions from the reported data that are always in `string` format. We will gradually improve this feature by adding more formats to our `regex` code as we encounted them.

+++

## Property Tables with Full References

The `get_property_table()` function provides comprehensive property data in a pandas DataFrame format with full reference information and parsed experimental values:

```{code-cell} ipython3
# Get comprehensive LogP data with references
logp_table = get_property_table(cids_aspirin[0], "LogP")
print("LogP Property Table:")
print(logp_table)
print(f"\nColumns: {list(logp_table.columns)}")

# Show some specific data
if len(logp_table) > 0:
    print(f"\nExample extracted values:")
    for i, row in logp_table.iterrows():
        if row['ExperimentalValue'] is not None:
            print(f"  Original: '{row['StringWithMarkup']}'")
            print(f"  Extracted: {row['ExperimentalValue']} {row['Unit'] if row['Unit'] else '(unitless)'}")
            break
```

```{code-cell} ipython3
# Compare different properties for aspirin
properties_to_check = ["Vapor Pressure", "Melting Point", "Boiling Point", "Solubility"]

for prop in properties_to_check:
    table = get_property_table(cids_aspirin[0], prop)
    if len(table) > 0:
        valid_values = table[table['ExperimentalValue'].notna()]
        print(f"{prop}: {len(valid_values)} experimental values extracted from {len(table)} total entries")
        if len(valid_values) > 0:
            # Show one example
            example = valid_values.iloc[0]
            print(f"  Example: {example['ExperimentalValue']} {example['Unit'] if example['Unit'] else ''}")
    else:
        print(f"{prop}: No data available")
    print()
```

## Advanced Pattern Recognition

PubChemView includes sophisticated pattern recognition for extracting experimental values from various text formats. The recent improvements include support for formats like "log Kow = 1.19" for LogP data:

```{code-cell} ipython3
# Demonstrate the improved LogP pattern recognition
logp_data = pcv.extract_property_data(cids_aspirin[0], "LogP")
print("LogP pattern recognition examples:")
for i, data in enumerate(logp_data):
    if data.value:  # Only show non-empty values
        # Test the extraction function directly
        exp_value, unit, temp, cond = pcv._extract_experimental_value_and_unit(data.value, "LogP")
        print(f"  {i+1}: '{data.value}' -> {exp_value} {unit if unit else '(unitless)'}")

# Test with a compound that has vapor pressure data
caffeine_cid = pc.get_cids_by_name("caffeine")[0]
vp_data = pcv.extract_property_data(caffeine_cid, "Vapor Pressure")
print(f"\nVapor Pressure pattern recognition examples (CID {caffeine_cid}):")
for i, data in enumerate(vp_data[:3]):  # Show first 3
    if data.value:
        exp_value, unit, temp, cond = pcv._extract_experimental_value_and_unit(data.value, "Vapor Pressure")
        print(f"  {i+1}: '{data.value}' -> {exp_value} {unit if unit else 'no unit'}")
```

## Batch Property Extraction

For multiple compounds, you can extract properties in batch:

```{code-cell} ipython3
# Batch extraction for multiple properties of aspirin
target_properties = ["LogP", "Melting Point", "Boiling Point", "Solubility", "Vapor Pressure"]
batch_results = pcv.batch_extract_properties(cids_aspirin[0], target_properties)

print("Batch extraction results for aspirin:")
for prop_name, prop_data in batch_results.items():
    print(f"\n{prop_name}: {len(prop_data)} entries")
    if prop_data:
        # Show first non-empty value
        for data in prop_data:
            if data.value and data.value.strip():
                print(f"  Example: {data.value}")
                break

# Extract all experimental properties for a compound
print(f"\n" + "="*50)
print("ALL EXPERIMENTAL PROPERTIES")
print("="*50)

all_properties = pcv.extract_all_experimental_properties(cids_aspirin[0])
print(f"Total experimental property categories found: {len(all_properties)}")
for prop_name, data_list in list(all_properties.items())[:5]:  # Show first 5
    print(f"{prop_name}: {len(data_list)} entries")
```

# Practical Use Cases

## Building a Property Database

Here's how you might build a property database for multiple compounds:

```{code-cell} ipython3
import pandas as pd

# Example: Build a small database of pharmaceutical compounds
pharma_compounds = {
    "aspirin": 2244,
    "ibuprofen": 3672, 
    "acetaminophen": 1983,
    "caffeine": 2519
}

# Create a comprehensive database
database_records = []

for name, cid in pharma_compounds.items():
    print(f"Processing {name} (CID: {cid})...")
    
    # Get basic compound info
    basic_info = pc.get_basic_compound_info(cid)
    
    # Get experimental properties
    logp_table = get_property_table(cid, "LogP")
    mp_table = get_property_table(cid, "Melting Point")
    sol_table = get_property_table(cid, "Solubility")
    
    # Extract first valid experimental value for each property
    logp_exp = logp_table[logp_table['ExperimentalValue'].notna()]['ExperimentalValue'].iloc[0] if len(logp_table[logp_table['ExperimentalValue'].notna()]) > 0 else None
    mp_exp = mp_table[mp_table['ExperimentalValue'].notna()]['ExperimentalValue'].iloc[0] if len(mp_table[mp_table['ExperimentalValue'].notna()]) > 0 else None
    sol_count = len(sol_table[sol_table['ExperimentalValue'].notna()])
    
    record = {
        'name': name,
        'cid': cid,
        'molecular_formula': basic_info.get('MolecularFormula'),
        'molecular_weight': basic_info.get('MolecularWeight'),
        'smiles': basic_info.get('CanonicalSMILES'),
        'logp_experimental': logp_exp,
        'melting_point_experimental': mp_exp,
        'solubility_data_points': sol_count
    }
    database_records.append(record)

# Create DataFrame
pharma_db = pd.DataFrame(database_records)
print("\nPharmaceutical Compounds Database:")
print(pharma_db.to_string(index=False))
```

## Error Handling and Best Practices

When working with PubChem APIs, it's important to handle errors gracefully:

```{code-cell} ipython3
from provesid import PubChemNotFoundError, PubChemError, PubChemViewNotFoundError

# Example of robust compound lookup
def safe_compound_lookup(identifier, search_type="name"):
    """Safely look up a compound with error handling"""
    try:
        if search_type == "name":
            cids = pc.get_cids_by_name(identifier)
        elif search_type == "smiles":
            cids = pc.get_cids_by_smiles(identifier)
        else:
            raise ValueError(f"Unsupported search type: {search_type}")
        
        if not cids:
            print(f"No compounds found for '{identifier}'")
            return None
            
        print(f"Found {len(cids)} compound(s) for '{identifier}': {cids[:5]}...")  # Show first 5
        return cids[0]  # Return first CID
        
    except PubChemNotFoundError:
        print(f"Compound '{identifier}' not found in PubChem")
        return None
    except PubChemError as e:
        print(f"PubChem API error: {e}")
        return None
    except Exception as e:
        print(f"Unexpected error: {e}")
        return None

# Test with various inputs
test_compounds = [
    ("aspirin", "name"),
    ("invalid_compound_name_xyz", "name"),
    ("CC(=O)OC1=CC=CC=C1C(=O)O", "smiles"),  # aspirin SMILES
    ("invalid_smiles", "smiles")
]

for compound, search_type in test_compounds:
    print(f"\nTesting: {compound} (search type: {search_type})")
    cid = safe_compound_lookup(compound, search_type)
    if cid:
        print(f"  Success! CID: {cid}")

# Safe property extraction
def safe_property_extraction(cid, property_name):
    """Safely extract property data with error handling"""
    try:
        data = pcv.extract_property_data(cid, property_name)
        return data
    except PubChemViewNotFoundError:
        print(f"Property '{property_name}' not found for CID {cid}")
        return []
    except Exception as e:
        print(f"Error extracting {property_name} for CID {cid}: {e}")
        return []

print(f"\n" + "="*40)
print("Safe property extraction example:")
logp_safe = safe_property_extraction(cids_aspirin[0], "LogP")
print(f"LogP data extracted safely: {len(logp_safe)} entries")
```

# Summary

This tutorial covered the comprehensive functionality of both PubChem APIs in the PROVESID package:

## PubChemAPI (Standard API)
- **Improved Data Access**: Methods like `get_compound_by_cid()` and `get_substance_by_sid()` now return data directly without redundant wrapper structures
- **Multiple Search Methods**: Search by name, SMILES, InChI key, CAS number
- **Comprehensive ID Resolution**: Find CIDs across both compound and substance domains
- **Batch Processing**: Handle multiple compounds efficiently
- **Property Extraction**: Get basic, selected, or all compound properties

## PubChemView (Experimental Properties API)
- **Experimental Properties**: Access to 40+ experimental properties not available in the standard API
- **Advanced Pattern Recognition**: Sophisticated text parsing for extracting numerical values from diverse formats
- **Property Tables**: Comprehensive DataFrames with full reference information
- **Batch Extraction**: Extract multiple properties for compounds efficiently
- **Temperature and Conditions**: Automatic extraction of experimental conditions

## Key Improvements
1. **Elegant Data Access**: No more `["PC_Compounds"][0]` or `["PC_Substances"][0]` needed
2. **Enhanced LogP Recognition**: Now supports "log Kow = 1.19" and similar formats
3. **Robust Error Handling**: Proper exception handling for network and data issues
4. **Comprehensive Pattern Support**: Handles scientific notation, comparison operators, and diverse units

## Best Practices
- Always use error handling for production code
- Use batch methods for multiple compounds
- Check data availability before processing
- Respect PubChem's rate limits (built into the APIs)

The PROVESID package provides a powerful and user-friendly interface to PubChem's vast chemical database, making it easy to integrate chemical data into your research workflows.
