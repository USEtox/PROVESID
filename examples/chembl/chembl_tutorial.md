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

# ChEMBL Database Tutorial

This tutorial demonstrates how to use the ChEMBL interface in PROVESID to query chemical compounds, structures, and properties.

## What is ChEMBL?

ChEMBL is a manually curated database of bioactive molecules with drug-like properties maintained by EMBL-EBI. It contains:
- Over 2.3 million distinct compounds
- Chemical structures and properties
- Bioactivity data from scientific literature
- Drug/clinical candidate information

## Setup

First, import the ChEMBL class. On first run, the database will be automatically downloaded (~5GB compressed, ~29GB uncompressed).

```{code-cell}
from provesid import CheMBL

# Initialize ChEMBL (auto-downloads database if needed)
chembl = CheMBL()
print("ChEMBL database loaded successfully!")
```

## 1. Search by ChEMBL ID

The most direct way to retrieve a compound is by its ChEMBL ID.

```{code-cell}
# Search for aspirin (CHEMBL25)
aspirin = chembl.search_by_chembl_id('CHEMBL25')

print(f"ChEMBL ID: {aspirin['chembl_id']}")
print(f"Name: {aspirin['pref_name']}")
print(f"SMILES: {aspirin['canonical_smiles']}")
print(f"InChI: {aspirin['standard_inchi']}")
print(f"InChI Key: {aspirin['standard_inchi_key']}")
print(f"Max Phase: {aspirin['max_phase']}")
print(f"\\nSynonyms ({len(aspirin['synonyms'])} total):")
for syn in aspirin['synonyms'][:5]:  # Show first 5 synonyms
    print(f"  - {syn}")
```

## 2. Search by Compound Name

Search for compounds by name (partial matching supported).

```{code-cell}
# Search for caffeine
results = chembl.search_by_name('caffeine')

print(f"Found {len(results)} compound(s) matching 'caffeine'\n")

for compound in results[:3]:  # Show first 3 results
    print(f"ChEMBL ID: {compound['chembl_id']}")
    print(f"Name: {compound['pref_name']}")
    print(f"SMILES: {compound['canonical_smiles']}")
    print("-" * 60)
```

## 3. Search by Chemical Structure

Search using InChI, InChI Key, or SMILES.

```{code-cell}
# Search by SMILES (aspirin)
smiles = 'CC(=O)Oc1ccccc1C(=O)O'
compound = chembl.search_by_smiles(smiles)

if compound:
    print(f"Found compound: {compound['pref_name']} ({compound['chembl_id']})")
else:
    print("Compound not found")
```

```{code-cell}
# Search by InChI Key (aspirin)
inchikey = 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'
compound = chembl.search_by_inchikey(inchikey)

if compound:
    print(f"Found compound: {compound['pref_name']} ({compound['chembl_id']})")
    print(f"SMILES: {compound['canonical_smiles']}")
else:
    print("Compound not found")
```

```{code-cell}
# Search by InChI (aspirin)
inchi = 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)'
compound = chembl.search_by_inchi(inchi)

if compound:
    print(f"Found compound: {compound['pref_name']} ({compound['chembl_id']})")
else:
    print("Compound not found")
```

## 4. Retrieve Physicochemical Properties

Get calculated molecular properties for a compound.

```{code-cell}
# Get aspirin's properties
aspirin = chembl.search_by_chembl_id('CHEMBL25')
props = chembl.get_properties(aspirin['molregno'])

if props:
    print(f"Molecular Properties for {aspirin['pref_name']}:")
    print(f"  Molecular Weight: {props['mw_freebase']:.2f}")
    print(f"  ALogP: {props['alogp']:.2f}")
    print(f"  Hydrogen Bond Acceptors: {props['hba']}")
    print(f"  Hydrogen Bond Donors: {props['hbd']}")
    print(f"  Polar Surface Area: {props['psa']:.2f}")
    print(f"  Rotatable Bonds: {props['rtb']}")
    print(f"  Aromatic Rings: {props['aromatic_rings']}")
    print(f"  Heavy Atoms: {props['heavy_atoms']}")
    print(f"  Lipinski Violations: {props['num_ro5_violations']}")
```

## 5. ID Conversion

Convert between ChEMBL IDs and internal molregno identifiers.

```{code-cell}
# ChEMBL ID to molregno
chembl_id = 'CHEMBL25'
molregno = chembl.chembl_id_to_molregno(chembl_id)
print(f"{chembl_id} -> molregno: {molregno}")

# molregno to ChEMBL ID
converted_id = chembl.molregno_to_chembl_id(molregno)
print(f"molregno {molregno} -> {converted_id}")
```

## 6. Complete Workflow Example

Let's search for a drug, retrieve its properties, and analyze them.

```{code-cell}
# Search for ibuprofen
results = chembl.search_by_name('ibuprofen')

if results:
    ibuprofen = results[0]
    print(f"===== {ibuprofen['pref_name']} =====")
    print(f"\nIdentifiers:")
    print(f"  ChEMBL ID: {ibuprofen['chembl_id']}")
    print(f"  Molregno: {ibuprofen['molregno']}")
    
    print(f"\nStructure:")
    print(f"  SMILES: {ibuprofen['canonical_smiles']}")
    print(f"  InChI Key: {ibuprofen['standard_inchi_key']}")
    
    # Get properties
    props = chembl.get_properties(ibuprofen['molregno'])
    if props:
        print(f"\nProperties:")
        print(f"  MW: {props['mw_freebase']:.2f}")
        print(f"  LogP: {props['alogp']:.2f}")
        print(f"  HBA: {props['hba']}")
        print(f"  HBD: {props['hbd']}")
        print(f"  PSA: {props['psa']:.2f}")
        
        # Lipinski's Rule of Five
        print(f"\nLipinski's Rule of Five:")
        print(f"  MW < 500: {props['mw_freebase'] < 500}")
        print(f"  LogP < 5: {props['alogp'] < 5}")
        print(f"  HBA < 10: {props['hba'] < 10}")
        print(f"  HBD < 5: {props['hbd'] < 5}")
        print(f"  Violations: {props['num_ro5_violations']}")
```

## 7. Batch Processing Example

Process multiple compounds at once.

```{code-cell}
import pandas as pd

# List of common drugs by ChEMBL ID
drug_ids = [
    'CHEMBL25',      # Aspirin
    'CHEMBL521',     # Ibuprofen
    'CHEMBL112',     # Acetaminophen
    'CHEMBL113',     # Caffeine
]

# Collect data
drug_data = []
for chembl_id in drug_ids:
    compound = chembl.search_by_chembl_id(chembl_id)
    if compound:
        props = chembl.get_properties(compound['molregno'])
        if props:
            drug_data.append({
                'ChEMBL ID': chembl_id,
                'Name': compound['pref_name'],
                'MW': props['mw_freebase'],
                'LogP': props['alogp'],
                'HBA': props['hba'],
                'HBD': props['hbd'],
                'PSA': props['psa'],
                'Ro5 Violations': props['num_ro5_violations']
            })

# Display as table
df = pd.DataFrame(drug_data)
print(df.to_string(index=False))
```

## 8. Advanced: Exploring Multiple Synonyms

Compounds can have multiple names and synonyms.

```{code-cell}
# Search for compounds with 'acetyl' in the name
results = chembl.search_by_name('acetyl', limit=5)

print(f"Found {len(results)} compounds with 'acetyl'\n")

for i, compound in enumerate(results, 1):
    print(f"{i}. {compound['pref_name']} ({compound['chembl_id']})")
    props = chembl.get_properties(compound['molregno'])
    if props:
        print(f"   MW: {props['mw_freebase']:.1f}, LogP: {props['alogp']:.2f}")
    print()
```

## Summary

In this tutorial, we covered:

1. ✓ Initializing the ChEMBL database
2. ✓ Searching by ChEMBL ID, name, InChI, InChI Key, and SMILES
3. ✓ Retrieving molecular properties
4. ✓ Converting between ID formats
5. ✓ Complete workflow examples
6. ✓ Batch processing multiple compounds

## Next Steps

- Explore the database schema in `src/provesid/data/schema_documentation.txt`
- Check out bioactivity data tables (activities, assays, targets)
- Combine with other PROVESID tools (PubChem, ChEBI, etc.)

## Resources

- ChEMBL Database: https://www.ebi.ac.uk/chembl/
- ChEMBL Documentation: https://chembl.gitbook.io/chembl-interface-documentation/
- PROVESID Documentation: [Link to your docs]
