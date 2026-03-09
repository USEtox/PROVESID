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

```{code-cell} ipython3
# Reload the module to pick up the new method
import importlib
import provesid.zeropm
importlib.reload(provesid.zeropm)
from provesid.zeropm import ZeroPM

# Initialize the ZeroPM database connection
zpm = ZeroPM()
print("✓ ZeroPM reloaded and initialized successfully!")
```

# ZeroPM Tutorial

This notebook demonstrates the functionality of the ZeroPM class, which provides efficient access to the ZeroPM SQLite database containing chemical identifiers and properties.

## Overview

ZeroPM allows you to:
- Query chemicals by CAS number or name
- Query chemicals by regulatory inventory, country, or region
- Convert between different chemical identifiers (CAS, InChI, InChIKey, SMILES)
- Perform fuzzy name searches
- Batch process multiple chemicals
- Search by substructure
- Export results to CSV

## Database Auto-Download

The ZeroPM database (~400MB) is automatically downloaded from [GitHub](https://github.com/ZeroPM-H2020/global-chemical-inventory-database) on first use. The download happens once and the database is cached locally for future use.

Let's get started!

+++

## 1. Installation and Setup

First, make sure you have the required dependencies installed.

```{code-cell} ipython3
# Import the ZeroPM class
from provesid.zeropm import ZeroPM

# Initialize the ZeroPM database connection
zpm = ZeroPM()

print("✓ ZeroPM initialized successfully!")
print(f"Database path: {zpm.db_path}")
```

```{code-cell} ipython3
# playing around with the ZeroPM instance
query_id = zpm.query_cas("50-00-0")
inchi_ids, ranks = zpm.get_inchi_id(query_id)
print(f"InChI IDs for CAS 50-00-0: {inchi_ids} with ranks {ranks}")
```

## Testing the new get_id_table_from_cas method

This method returns a pandas DataFrame with all identifiers (CAS, InChI, InChIKey, synonyms) for a given CAS number.

```{code-cell} ipython3
# Get identifier table for formaldehyde (CAS: 50-00-0)
df = zpm.get_id_table_from_cas("50-00-0")
print(df)
print(f"\nNumber of rows: {len(df)}")
```

```{code-cell} ipython3
# Display the full table with better formatting
print("\nFull table with all columns:")
print(df.to_string())

# Show specific columns
print("\n\nKey columns:")
print(df[['cas', 'inchi', 'inchikey']].to_string())
```

```{code-cell} ipython3
# Test with another compound that might have multiple entries
# Let's try with aspirin
df2 = zpm.get_id_table_from_cas("50-78-2")  # Aspirin
if df2 is not None:
    print("\nAspirin (CAS: 50-78-2):")
    print(df2[['cas', 'inchi_id', 'rank', 'inchikey']].to_string())
else:
    print("\nAspirin not found in database")
```

### Benefits of get_id_table_from_cas

The `get_id_table_from_cas()` method provides several advantages:

1. **Comprehensive identifier retrieval**: Returns CAS, InChI, InChIKey, and synonyms in one call
2. **Handles multiple structures**: Some CAS numbers map to multiple InChI structures (different forms, salts, etc.)
3. **Rank information**: Includes rank to indicate the relevance/confidence of each match
4. **Easy data manipulation**: Returns a pandas DataFrame for easy filtering, sorting, and export

```{code-cell} ipython3
df2
```

### Batch Processing Multiple CAS Numbers

For processing multiple CAS numbers at once, use `batch_get_id_table_from_cas()`:

```{code-cell} ipython3
# Process multiple CAS numbers at once
cas_list = ["50-00-0", "50-78-2", "64-17-5"]  # formaldehyde, aspirin, ethanol

batch_df = zpm.batch_get_id_table_from_cas(cas_list)

print(f"Total rows: {len(batch_df)}")
print(f"\nNumber of structures per CAS:")
print(batch_df.groupby('cas').size())
print(f"\nFirst few rows:")
print(batch_df[['cas', 'inchi_id', 'rank', 'inchikey']].head(10))
```

```{code-cell} ipython3
# The batch method handles missing CAS numbers gracefully
mixed_cas_list = ["50-00-0", "999-99-9", "50-78-2"]  # valid, invalid, valid

batch_df_mixed = zpm.batch_get_id_table_from_cas(mixed_cas_list)
print(f"Requested {len(mixed_cas_list)} CAS numbers")
print(f"Found {len(batch_df_mixed['cas'].unique())} in database")
print(f"CAS numbers found: {list(batch_df_mixed['cas'].unique())}")
```

## 2. Database Statistics

Let's start by exploring what's in the database.

```{code-cell} ipython3
# Get database statistics
stats = zpm.get_database_stats()

print("Database Statistics:")
print("=" * 50)
for key, value in stats.items():
    print(f"{key:30s}: {value:>15,}" if isinstance(value, int) else f"{key:30s}: {value}")
```

## 3. Querying by CAS Number

The most common way to query the database is using a CAS Registry Number.

```{code-cell} ipython3
# Get a sample CAS number from the database
zpm.cursor.execute("""
    SELECT query 
    FROM api_ready_query 
    WHERE type = 'CAS Registry Number' 
    LIMIT 1
""")
sample_cas = zpm.cursor.fetchone()[0]

print(f"Sample CAS number: {sample_cas}")

# Query the CAS number to get a query_id
query_id = zpm.query_cas(sample_cas)
print(f"Query ID: {query_id}")

# Get InChI information
inchi_ids, ranks = zpm.get_inchi_id(query_id)
if inchi_ids:
    inchi, inchikey = zpm.get_inchi(inchi_ids[0])
    print(f"\nInChI: {inchi[:50]}..." if len(inchi) > 50 else f"\nInChI: {inchi}")
    print(f"InChIKey: {inchikey}")
```

## 4. Converting CAS to SMILES

ZeroPM can convert CAS numbers to SMILES using RDKit.

```{code-cell} ipython3
# Get SMILES from CAS number
smiles = zpm.get_smiles_from_cas(sample_cas)
print(f"CAS: {sample_cas}")
print(f"SMILES: {smiles}")

# Get chemical names
names = zpm.get_names(sample_cas)
if names:
    print(f"\nAlternative names ({len(names)}):")
    for i, name in enumerate(names[:5], 1):  # Show first 5 names
        print(f"  {i}. {name}")
    if len(names) > 5:
        print(f"  ... and {len(names) - 5} more")
```

## 5. Querying by Chemical Name

You can search for chemicals by their exact name or use fuzzy matching.

```{code-cell} ipython3
# Get a sample chemical name
zpm.cursor.execute("""
    SELECT query 
    FROM api_ready_query 
    WHERE type = 'chemical name' 
    LIMIT 1
""")
sample_name = zpm.cursor.fetchone()[0]

# Exact match
query_id = zpm.query_name(sample_name)
print(f"Exact search for '{sample_name}'")
print(f"Query ID: {query_id}")

# Fuzzy match (with partial name)
if len(sample_name) >= 5:
    partial_name = sample_name[:5]
    print(f"\nFuzzy search for '{partial_name}':")
    similar_ids = zpm.query_similar_name(partial_name, number_of_results=5, score_cutoff=70)
    if similar_ids:
        print(f"Found {len(similar_ids)} similar matches")
        for qid in similar_ids[:3]:
            zpm.cursor.execute("SELECT query FROM api_ready_query WHERE query_id = ?", (qid,))
            result = zpm.cursor.fetchone()
            if result:
                print(f"  - {result[0]}")
```

## 6. Converting Between Different Identifiers

ZeroPM supports conversion between various chemical identifiers.

```{code-cell} ipython3
# Get a sample InChI and InChIKey
zpm.cursor.execute("""
    SELECT inchi, inchikey 
    FROM substances 
    LIMIT 1
""")
test_inchi, test_inchikey = zpm.cursor.fetchone()

print("Identifier Conversions:")
print("=" * 60)

# InChI to CAS
cas_from_inchi = zpm.get_cas_from_inchi(test_inchi)
print(f"InChI → CAS: {cas_from_inchi}")

# InChIKey to CAS
cas_from_key = zpm.get_cas_from_inchikey(test_inchikey)
print(f"InChIKey → CAS: {cas_from_key}")

# InChIKey to SMILES
smiles_from_key = zpm.get_smiles_from_inchikey(test_inchikey)
print(f"InChIKey → SMILES: {smiles_from_key}")

# SMILES to CAS (if we have a valid SMILES)
if smiles_from_key:
    cas_from_smiles = zpm.get_cas_from_smiles(smiles_from_key)
    print(f"SMILES → CAS: {cas_from_smiles}")
```

## 7. Batch Processing

For efficiency, ZeroPM provides batch methods to process multiple chemicals at once.

```{code-cell} ipython3
# Get multiple CAS numbers
zpm.cursor.execute("""
    SELECT query 
    FROM api_ready_query 
    WHERE type = 'CAS Registry Number' 
    LIMIT 5
""")
cas_list = [row[0] for row in zpm.cursor.fetchall()]

print(f"Batch processing {len(cas_list)} CAS numbers:")
print("=" * 60)

# Batch query CAS numbers
query_ids = zpm.batch_query_cas(cas_list)
for cas, qid in query_ids.items():
    print(f"{cas}: Query ID = {qid}")

print("\n" + "=" * 60)

# Batch get SMILES
smiles_dict = zpm.batch_get_smiles_from_cas(cas_list)
for cas, smiles in smiles_dict.items():
    print(f"{cas}: {smiles if smiles else 'N/A'}")
```

## 8. Batch InChIKey to CAS Conversion

```{code-cell} ipython3
# Get multiple InChIKeys
zpm.cursor.execute("""
    SELECT inchikey 
    FROM substances 
    LIMIT 5
""")
inchikey_list = [row[0] for row in zpm.cursor.fetchall()]

print(f"Batch converting {len(inchikey_list)} InChIKeys to CAS:")
print("=" * 60)

cas_dict = zpm.batch_get_cas_from_inchikey(inchikey_list)
for key, cas in cas_dict.items():
    print(f"{key}: {cas if cas else 'N/A'}")
```

## 9. Advanced Search: Regex Pattern Matching

Search for chemicals using pattern matching.

```{code-cell} ipython3
# Search for chemicals with names containing a pattern
# For example, search for names containing "acid"
pattern = "%acid%"

results = zpm.query_name_regex(pattern, case_sensitive=False, limit=10)

print(f"Chemical names matching pattern '{pattern}':")
print("=" * 60)
for query_id, name in results[:5]:
    print(f"  {name}")
print(f"\nTotal matches: {len(results)}")
```

## 10. Advanced Search: Substructure Search

Search for chemicals containing a specific substructure (SMARTS pattern).

⚠️ **Note:** This operation can be slow for large searches as it needs to check each molecule.

```{code-cell} ipython3
# Search for molecules containing a benzene ring
smarts_pattern = "c1ccccc1"  # Benzene ring

print(f"Searching for molecules with benzene ring (max 5 results)...")
results = zpm.get_cas_by_substructure(smarts_pattern, max_results=5)

print(f"Found {len(results)} molecules with benzene ring:")
print("=" * 60)
for i, compound in enumerate(results, 1):
    print(f"\n{i}. CAS: {compound['cas']}")
    print(f"   SMILES: {compound['smiles']}")
    print(f"   InChIKey: {compound['inchikey'][:27]}...")
```

## 11. Exporting Data to CSV

You can export query results to CSV files for further analysis.

```{code-cell} ipython3
import tempfile
import os

# Create a temporary directory for exports
temp_dir = tempfile.mkdtemp()

# Export batch results
output_file = os.path.join(temp_dir, 'cas_smiles_export.csv')
zpm.export_to_csv(
    list(smiles_dict.items()), 
    output_file,
    columns=['CAS', 'SMILES']
)

print(f"✓ Exported data to: {output_file}")

# Export custom query results
sql_query = """
    SELECT aq.query AS CAS, s.inchikey
    FROM api_ready_query aq
    JOIN api_results ar ON aq.query_id = ar.query_id
    JOIN substances s ON ar.inchi_id = s.inchi_id
    WHERE aq.type = 'CAS Registry Number'
    LIMIT 10
"""
output_file2 = os.path.join(temp_dir, 'cas_inchikey_export.csv')
zpm.export_query_results(sql_query, output_file2, include_headers=True)

print(f"✓ Exported custom query to: {output_file2}")

# List exported files
print(f"\nExported files in {temp_dir}:")
for file in os.listdir(temp_dir):
    filepath = os.path.join(temp_dir, file)
    size = os.path.getsize(filepath)
    print(f"  - {file} ({size:,} bytes)")
```

## 12. Performance Optimization: Creating Indexes

Create database indexes to speed up queries.

```{code-cell} ipython3
# Create indexes for better query performance
print("Creating database indexes...")
index_results = zpm.create_indexes()

print("\nIndex Status:")
print("=" * 60)
for index_name, status in index_results.items():
    print(f"{index_name:30s}: {status}")

print("\n✓ Indexes created successfully!")
print("Note: Subsequent queries will be faster with these indexes.")
```

## 13. Creating Custom Views

Create database views for frequently used queries.

```{code-cell} ipython3
# Create a view for CAS to InChI mapping
view_name = "cas_to_inchi_view"
sql_query = """
    SELECT aq.query AS cas, s.inchi, s.inchikey
    FROM api_ready_query aq
    JOIN api_results ar ON aq.query_id = ar.query_id
    JOIN substances s ON ar.inchi_id = s.inchi_id
    WHERE aq.type = 'CAS Registry Number' AND ar.rank = 1
"""

success = zpm.create_view(view_name, sql_query)

if success:
    print(f"✓ View '{view_name}' created successfully!")
    
    # Query the view
    zpm.cursor.execute(f"SELECT * FROM {view_name} LIMIT 5")
    print(f"\nSample data from view:")
    print("=" * 60)
    for row in zpm.cursor.fetchall():
        print(f"CAS: {row[0]}, InChIKey: {row[2][:27]}...")
    
    # Clean up - drop the view
    zpm.cursor.execute(f"DROP VIEW IF EXISTS {view_name}")
    zpm.conn.commit()
    print(f"\n✓ View dropped for cleanup")
else:
    print(f"✗ Failed to create view")
```

## 14. Complete Example: Workflow for Multiple Chemicals

Here's a complete workflow demonstrating how to process multiple chemicals efficiently.

```{code-cell} ipython3
import pandas as pd

# Get a sample of CAS numbers to process
zpm.cursor.execute("""
    SELECT query 
    FROM api_ready_query 
    WHERE type = 'CAS Registry Number' 
    LIMIT 10
""")
cas_numbers = [row[0] for row in zpm.cursor.fetchall()]

print(f"Processing {len(cas_numbers)} chemicals...")
print("=" * 80)

# Batch get all the data we need
query_ids = zpm.batch_query_cas(cas_numbers)
smiles_data = zpm.batch_get_smiles_from_cas(cas_numbers)
names_data = zpm.batch_get_names(cas_numbers)

# Create a pandas DataFrame
data = []
for cas in cas_numbers:
    data.append({
        'CAS': cas,
        'Query_ID': query_ids.get(cas),
        'SMILES': smiles_data.get(cas),
        'Names_Count': len(names_data.get(cas, [])),
        'First_Name': names_data.get(cas, [''])[0] if names_data.get(cas) else ''
    })

df = pd.DataFrame(data)

print("\nResults Summary:")
print(df.to_string(index=False, max_colwidth=50))

print(f"\n✓ Processed {len(cas_numbers)} chemicals successfully!")
print(f"  - {df['SMILES'].notna().sum()} have SMILES")
print(f"  - {df[df['Names_Count'] > 0].shape[0]} have alternative names")
```

## 15. Summary and Best Practices

### Key Features:
1. **Simple Queries**: `query_cas()`, `query_name()`
2. **Fuzzy Matching**: `query_similar_name()` with configurable score cutoff
3. **Identifier Conversion**: Convert between CAS, InChI, InChIKey, and SMILES
4. **Batch Operations**: Process multiple chemicals efficiently
5. **Advanced Search**: Regex patterns and substructure matching
6. **Export**: Save results to CSV for further analysis
7. **Performance**: Create indexes for faster queries

### Best Practices:
- Use batch methods when processing multiple chemicals
- Create indexes before running many queries
- Use fuzzy matching for user input with potential typos
- Set appropriate score cutoffs for fuzzy matching (70-90 is typical)
- Export results to CSV for sharing or further analysis

### Performance Tips:
- Batch operations are much faster than individual queries
- Create indexes once at the start if doing many queries
- Use `score_cutoff` parameter to limit fuzzy search results
- Limit substructure searches with `max_results` parameter

+++

## 14. New ID Table Methods - Get Complete Identifier Tables

The ZeroPM class now provides four new methods to retrieve complete identifier tables from different starting points: InChI, InChIKey, and chemical names. These complement the existing `get_id_table_from_cas()` method.

+++

### 14.1 Get ID Table from InChI

The `get_id_table_from_inchi()` method returns all identifiers for a given InChI string.

```{code-cell} ipython3
# First, let's get an InChI from a known CAS to use as an example
formaldehyde_table = zpm.get_id_table_from_cas("50-00-0")
if formaldehyde_table is not None and len(formaldehyde_table) > 0:
    example_inchi = formaldehyde_table['inchi'].iloc[0]
    print(f"Example InChI: {example_inchi}\n")
    
    # Now get the ID table from InChI
    df_from_inchi = zpm.get_id_table_from_inchi(example_inchi)
    print("ID Table from InChI:")
    print(df_from_inchi)
    print(f"\nColumns: {list(df_from_inchi.columns)}")
else:
    print("Could not find formaldehyde in database")
```

### 14.2 Batch Get ID Tables from InChI List

Process multiple InChI strings at once:

```{code-cell} ipython3
# Get multiple InChIs from the database to use as examples
example_cas_list = ["50-00-0", "50-78-2", "64-17-5"]  # formaldehyde, aspirin, ethanol
inchi_list = []

for cas in example_cas_list:
    table = zpm.get_id_table_from_cas(cas)
    if table is not None and len(table) > 0:
        inchi_list.append(table['inchi'].iloc[0])

print(f"Testing with {len(inchi_list)} InChI strings\n")

# Batch process
batch_df = zpm.batch_get_id_table_from_inchi(inchi_list)
print("Batch ID Table from InChI list:")
print(batch_df[['inchi', 'inchikey', 'cas', 'rank']].head(10))
print(f"\nTotal rows: {len(batch_df)}")
print(f"Unique InChIs processed: {batch_df['inchi'].nunique()}")
```

### 14.3 Get ID Table from InChIKey

The `get_id_table_from_inchikey()` method returns all identifiers for a given InChIKey string.

```{code-cell} ipython3
# Get an InChIKey example from formaldehyde
if formaldehyde_table is not None and len(formaldehyde_table) > 0:
    example_inchikey = formaldehyde_table['inchikey'].iloc[0]
    print(f"Example InChIKey: {example_inchikey}\n")
    
    # Now get the ID table from InChIKey
    df_from_inchikey = zpm.get_id_table_from_inchikey(example_inchikey)
    print("ID Table from InChIKey:")
    print(df_from_inchikey)
    print(f"\nColumns: {list(df_from_inchikey.columns)}")
else:
    print("Could not find formaldehyde in database")
```

### 14.4 Batch Get ID Tables from InChIKey List

Process multiple InChIKey strings at once:

```{code-cell} ipython3
# Get InChIKeys from our example CAS list
inchikey_list = []
for cas in example_cas_list:
    table = zpm.get_id_table_from_cas(cas)
    if table is not None and len(table) > 0:
        inchikey_list.append(table['inchikey'].iloc[0])

print(f"Testing with {len(inchikey_list)} InChIKey strings:")
for key in inchikey_list:
    print(f"  - {key}")

# Batch process
batch_df_keys = zpm.batch_get_id_table_from_inchikey(inchikey_list)
print("\nBatch ID Table from InChIKey list:")
print(batch_df_keys[['inchikey', 'cas', 'rank']].head(10))
print(f"\nTotal rows: {len(batch_df_keys)}")
print(f"Unique InChIKeys processed: {batch_df_keys['inchikey'].nunique()}")
```

### 14.5 Get ID Table from Chemical Name

The `get_id_table_from_name()` method returns all identifiers for a given chemical name (exact match).

```{code-cell} ipython3
# Try to find a chemical name in the database
# First, let's see what names are available for formaldehyde
formaldehyde_names = zpm.get_names("50-00-0")
print(f"Available names for formaldehyde (CAS 50-00-0): {formaldehyde_names[:5]}\n")

if formaldehyde_names:
    # Use the first name
    example_name = formaldehyde_names[0]
    print(f"Using name: '{example_name}'\n")
    
    # Get ID table from name
    df_from_name = zpm.get_id_table_from_name(example_name)
    if df_from_name is not None:
        print("ID Table from Name:")
        print(df_from_name)
        print(f"\nColumns: {list(df_from_name.columns)}")
    else:
        print(f"Name '{example_name}' not found in database as a query")
else:
    print("No names found for formaldehyde")
```

### 14.6 Batch Get ID Tables from Name List

Process multiple chemical names at once:

```{code-cell} ipython3
# Get a few chemical names from the database
# Query some names directly from the database
zpm.cursor.execute("""
    SELECT query 
    FROM api_ready_query 
    WHERE type = 'chemical name'
    LIMIT 5
""")
name_results = zpm.cursor.fetchall()
name_list = [row[0] for row in name_results]

print(f"Testing with {len(name_list)} chemical names:")
for name in name_list:
    print(f"  - {name}")

# Batch process
batch_df_names = zpm.batch_get_id_table_from_name(name_list)
if batch_df_names is not None and len(batch_df_names) > 0:
    print("\nBatch ID Table from Name list:")
    print(batch_df_names[['name', 'cas', 'inchikey', 'rank']].head(10))
    print(f"\nTotal rows: {len(batch_df_names)}")
    print(f"Unique names processed: {batch_df_names['name'].nunique()}")
else:
    print("\nNo results found for the name list")
```

### 14.7 Summary of ID Table Methods

All six ID table methods provide:
- **Complete identifier mapping**: Returns all related identifiers in a DataFrame
- **Rank information**: Shows the relevance/confidence of each match
- **Batch processing**: Efficient handling of multiple queries
- **Consistent output format**: Easy to combine and analyze results

**Use cases:**
- `get_id_table_from_cas()` - When you have CAS numbers and need all associated identifiers
- `get_id_table_from_inchi()` - When you have InChI strings from calculations or other sources
- `get_id_table_from_inchikey()` - When you have InChIKeys from databases or publications
- `get_id_table_from_name()` - When you have exact chemical names (not fuzzy matching)

The batch versions are recommended when processing multiple chemicals for better performance.

+++

## 15. CAS Conversion Methods - Convert to CAS from Various Identifiers

The ZeroPM class provides methods to convert from different types of chemical identifiers to CAS numbers.

+++

### 15.1 Get CAS from Chemical Name

Convert exact chemical names to CAS numbers:

```{code-cell} ipython3
# Example: Get CAS from a chemical name
# Use a name from our previous examples
if formaldehyde_names:
    test_name = formaldehyde_names[0]  # "Formalin"
    print(f"Chemical name: {test_name}")
    
    cas_from_name = zpm.get_cas_from_name(test_name)
    print(f"CAS number(s): {cas_from_name}")
    print(f"Type: {type(cas_from_name)}")
    
    # If multiple CAS numbers are returned
    if isinstance(cas_from_name, list):
        print(f"\nFound {len(cas_from_name)} CAS numbers for this name")
        for i, cas in enumerate(cas_from_name, 1):
            print(f"  {i}. {cas}")
else:
    # Try with another name
    cas_from_name = zpm.get_cas_from_name("methanol")
    print(f"Methanol CAS: {cas_from_name}")
```

### 15.2 Get CAS from SMILES

Convert SMILES strings to CAS numbers (already existed, shown for completeness):

```{code-cell} ipython3
# Example: Get CAS from SMILES
smiles_examples = {
    "C": "Methane",
    "CO": "Methanol", 
    "CCO": "Ethanol",
    "C=O": "Formaldehyde"
}

print("Converting SMILES to CAS:")
for smiles, name in smiles_examples.items():
    cas = zpm.get_cas_from_smiles(smiles)
    print(f"  {smiles:6s} ({name:15s}): {cas}")
```

### 15.3 Get CAS from Molecular Formula

Convert molecular formulas to CAS numbers. Note that formulas are not unique - many isomers can share the same formula:

```{code-cell} ipython3
# Example: Get CAS from molecular formula
# Warning: This can be slow for the first run as it scans the entire database
formulas_to_test = ["CH2O", "H2O", "CH4O"]

print("Converting molecular formulas to CAS numbers:\n")
for formula in formulas_to_test:
    print(f"Formula: {formula}")
    cas_list = zpm.get_cas_from_formula(formula)
    
    if cas_list:
        print(f"  Found {len(cas_list)} chemicals with this formula")
        # Show first 5 CAS numbers
        for i, cas in enumerate(cas_list[:5], 1):
            print(f"    {i}. {cas}")
        if len(cas_list) > 5:
            print(f"    ... and {len(cas_list) - 5} more")
    else:
        print(f"  No chemicals found with formula {formula}")
    print()
```

### 15.4 Batch CAS Conversion Methods

Process multiple identifiers efficiently:

```{code-cell} ipython3
# Batch convert SMILES to CAS
smiles_list = ["C", "CO", "CCO", "C=O"]
print("Batch SMILES to CAS conversion:")
results_smiles = zpm.batch_get_cas_from_smiles(smiles_list)
for smiles, cas in results_smiles.items():
    print(f"  {smiles:6s} -> {cas}")

print("\n" + "="*50 + "\n")

# Batch convert names to CAS
if formaldehyde_names:
    # Use a few names from formaldehyde
    test_names = formaldehyde_names[:3]
    print(f"Batch names to CAS conversion:")
    results_names = zpm.batch_get_cas_from_name(test_names)
    for name, cas in results_names.items():
        print(f"  {name:20s} -> {cas}")

print("\n" + "="*50 + "\n")

# Batch convert formulas to CAS
formula_batch = ["CH4", "CH4O", "C2H6O"]
print("Batch formulas to CAS conversion:")
results_formulas = zpm.batch_get_cas_from_formula(formula_batch)
for formula, cas_list in results_formulas.items():
    if cas_list:
        print(f"  {formula:8s} -> {len(cas_list)} chemicals found")
    else:
        print(f"  {formula:8s} -> Not found")
```

### 15.5 Summary of CAS Conversion Methods

**Available methods:**
- `get_cas_from_name(name)` - Convert exact chemical name to CAS
- `get_cas_from_smiles(smiles)` - Convert SMILES to CAS (via InChI)
- `get_cas_from_formula(formula)` - Convert molecular formula to CAS (returns list, non-unique)
- `get_cas_from_inchi(inchi)` - Convert InChI to CAS
- `get_cas_from_inchikey(inchikey)` - Convert InChIKey to CAS

**Batch versions:**
- `batch_get_cas_from_name(name_list)`
- `batch_get_cas_from_smiles(smiles_list)`
- `batch_get_cas_from_formula(formula_list)`
- `batch_get_cas_from_inchikey(inchikey_list)`

**Important notes:**
- Name conversion requires exact matches. Use `query_similar_name()` for fuzzy matching
- Formula conversion is not unique - multiple chemicals can have the same formula
- SMILES conversion works by first converting to InChI
- All methods return `None` if no match is found
- Some methods may return a list of CAS numbers if multiple matches exist

+++

## Conclusion

This notebook has demonstrated the main features of the ZeroPM class:

✓ Querying by CAS number and chemical name  
✓ Fuzzy name matching for handling variations  
✓ Converting between chemical identifiers  
✓ Batch processing for efficiency  
✓ Advanced searches (regex, substructure)  
✓ Exporting results to CSV  
✓ Performance optimization with indexes  

The ZeroPM class provides a convenient Python interface to the ZeroPM database, making it easy to work with chemical identifiers in your research and applications.

For more information, see the [PROVESID documentation](https://usetox.github.io/PROVESID/).
