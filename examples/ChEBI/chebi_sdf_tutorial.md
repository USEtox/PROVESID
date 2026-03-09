---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.1
---

## 1. Initialize ChebiSDF

On first initialization, an index will be built (takes ~15 seconds). Subsequent loads use the cached index.

```{code-cell}
from provesid.chebi import ChebiSDF

# Initialize (will use cached index if available)
chebi_sdf = ChebiSDF()

# Get database statistics
stats = chebi_sdf.get_database_stats()
print("ChEBI SDF Database Statistics:")
print("=" * 50)
for key, value in stats.items():
    print(f"{key:30s}: {value:>15,}")
```

## 2. Query by ChEBI ID

```{code-cell}
# Get water (CHEBI:15377)
water = chebi_sdf.get_compound_by_id("CHEBI:15377")

print(f"ChEBI ID: {water['ChEBI ID']}")
print(f"Name: {water['ChEBI NAME']}")
print(f"Formula: {water.get('FORMULA', 'N/A')}")
print(f"Mass: {water.get('MASS', 'N/A')}")
print(f"SMILES: {water.get('SMILES', 'N/A')}")
print(f"InChI: {water.get('INCHI', 'N/A')}")
print(f"InChIKey: {water.get('INCHIKEY', 'N/A')}")
print(f"Star Rating: {water.get('STAR', 'N/A')}")
```

## 3. Search by Name

```{code-cell}
# Exact name search
results = chebi_sdf.search_by_name("glucose", exact=True)
print(f"Found {len(results)} compound(s) with exact name 'glucose'")
for result in results[:3]:
    print(f"  - {result['ChEBI ID']}: {result['ChEBI NAME']}")

print()

# Partial name search
results = chebi_sdf.search_by_name("glucose", exact=False)
print(f"Found {len(results)} compound(s) with 'glucose' in the name")
for result in results[:5]:
    print(f"  - {result['ChEBI ID']}: {result['ChEBI NAME']}")
```

## 4. Search by CAS Number

```{code-cell}
# Search for aspirin (CAS: 50-78-2)
results = chebi_sdf.search_by_cas("50-78-2")

if results:
    aspirin = results[0]
    print(f"ChEBI ID: {aspirin['ChEBI ID']}")
    print(f"Name: {aspirin['ChEBI NAME']}")
    print(f"Formula: {aspirin.get('FORMULA', 'N/A')}")
    print(f"\nSynonyms:")
    if 'SYNONYM' in aspirin:
        synonyms = aspirin['SYNONYM'].split(';')
        for syn in synonyms[:5]:
            print(f"  - {syn.strip()}")
```

## 5. Search by Structure Identifiers

```{code-cell}
# Search by InChIKey (water)
inchikey = "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
result = chebi_sdf.search_by_inchikey(inchikey)

if result:
    print(f"Found compound: {result['ChEBI NAME']} ({result['ChEBI ID']})")

# Search by molecular formula
results = chebi_sdf.search_by_formula("H2O")
print(f"\nFound {len(results)} compound(s) with formula H2O")
for r in results[:3]:
    print(f"  - {r['ChEBI ID']}: {r['ChEBI NAME']}")
```

## 6. Search by Synonym

```{code-cell}
# Search for compounds with "acetylsalicylic acid" as a synonym
results = chebi_sdf.search_by_synonym("acetylsalicylic acid", exact=False)

print(f"Found {len(results)} compound(s)")
for result in results:
    print(f"\nChEBI ID: {result['ChEBI ID']}")
    print(f"Name: {result['ChEBI NAME']}")
    print(f"Formula: {result.get('FORMULA', 'N/A')}")
```

## 7. Batch Operations and DataFrame Export

```{code-cell}
# Get multiple compounds at once
chebi_ids = ["CHEBI:15377", "CHEBI:16236", "CHEBI:17234"]  # water, ethanol, glucose

# Export to DataFrame
df = chebi_sdf.export_to_dataframe(chebi_ids)

print("DataFrame with selected compounds:")
print(df[['ChEBI ID', 'ChEBI NAME', 'FORMULA', 'MASS', 'STAR']])
```

## 8. Working with High-Quality Compounds

ChEBI assigns star ratings (1-3) to compounds, where 3 stars indicates the highest quality/completeness.

```{code-cell}
# Note: This will take a few minutes as it scans all compounds
# Uncomment to run:

# three_star_ids = chebi_sdf.filter_by_star_rating(min_stars=3)
# print(f"Found {len(three_star_ids)} compounds with 3-star rating")

# # Export first 100 to DataFrame
# df_high_quality = chebi_sdf.export_to_dataframe(three_star_ids[:100])
# print(df_high_quality.head())
```

## 9. Accessing External Database Links

```{code-cell}
# Get water and show all available database links
water = chebi_sdf.get_compound_by_id("CHEBI:15377")

print("External Database Links for Water:")
print("=" * 50)
for key, value in water.items():
    if 'Database Links' in key or 'Registry Number' in key:
        print(f"{key}: {value}")
```

## Summary

The `ChebiSDF` class provides:

1. **Fast queries** using pre-built index (~190,000 compounds in <1 second after indexing)
2. **Multiple search methods**: by ID, name, synonym, CAS, InChI, InChIKey, formula
3. **Batch operations**: get multiple compounds at once
4. **DataFrame export**: easy integration with pandas for analysis
5. **Offline access**: no internet required after SDF file is downloaded
6. **Complete data**: structure, properties, synonyms, and 80+ database cross-references
