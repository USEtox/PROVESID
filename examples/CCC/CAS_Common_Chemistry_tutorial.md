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

# CAS Common Chemistry API Tutorial

This tutorial demonstrates how to use the `CASCommonChem` class from the `provesid` package to access chemical information from the CAS Common Chemistry database. The CAS Common Chemistry API provides access to chemical information for more than 500,000 chemical substances from CAS REGISTRY®.

+++

## 1. Import and Initialize

First, import the class and create an instance:

```{code-cell} ipython3
from provesid import CASCommonChem
ccc = CASCommonChem()
print("CASCommonChem initialized successfully!")
print(f"Base URL: {ccc.base_url}")
```

## 2. Lookup by CAS Registry Number

The most direct way to get chemical information is by using a CAS Registry Number (CAS RN). Let's look up some common compounds:

```{code-cell} ipython3
# Lookup water by CAS RN
water_info = ccc.cas_to_detail("7732-18-5")
print("Water (7732-18-5):")
print(f"  Name: {water_info.get('name')}")
print(f"  Molecular Formula: {water_info.get('molecularFormula')}")
print(f"  Molecular Mass: {water_info.get('molecularMass')}")
print(f"  SMILES: {water_info.get('smile')}")
print(f"  InChI: {water_info.get('inchi')}")
print(f"  Status: {water_info.get('status')}")
```

```{code-cell} ipython3
water_info
```

```{code-cell} ipython3
# Lookup aspirin by CAS RN
aspirin_info = ccc.cas_to_detail("50-78-2")
print("Aspirin (50-78-2):")
print(f"  Name: {aspirin_info.get('name')}")
print(f"  Molecular Formula: {aspirin_info.get('molecularFormula')}")
print(f"  Molecular Mass: {aspirin_info.get('molecularMass')}")
print(f"  SMILES: {aspirin_info.get('smile')}")
print(f"  Number of synonyms: {len(aspirin_info.get('synonyms', []))}")
print(f"  First 5 synonyms: {aspirin_info.get('synonyms', [])[:5]}")
```

## 3. Search by Name

You can search for compounds using their common names or IUPAC names:

```{code-cell} ipython3
# Search by common name
caffeine_info = ccc.name_to_detail("caffeine")
print("Caffeine (by name search):")
print(f"  CAS RN: {caffeine_info.get('rn')}")
print(f"  Name: {caffeine_info.get('name')}")
print(f"  Molecular Formula: {caffeine_info.get('molecularFormula')}")
print(f"  Molecular Mass: {caffeine_info.get('molecularMass')}")
print(f"  Status: {caffeine_info.get('status')}")
```

```{code-cell} ipython3
# Search by IUPAC name
acetone_info = ccc.name_to_detail("propan-2-one")
print("Acetone (by IUPAC name 'propan-2-one'):")
print(f"  CAS RN: {acetone_info.get('rn')}")
print(f"  Name: {acetone_info.get('name')}")
print(f"  Molecular Formula: {acetone_info.get('molecularFormula')}")
print(f"  SMILES: {acetone_info.get('smile')}")
print(f"  Status: {acetone_info.get('status')}")
```

## 4. Search by SMILES

You can also search using SMILES notation:

```{code-cell} ipython3
# Search by SMILES string
ethanol_smiles = "CCO"
ethanol_info = ccc.smiles_to_detail(ethanol_smiles)
print(f"Compound with SMILES '{ethanol_smiles}':")
print(f"  CAS RN: {ethanol_info.get('rn')}")
print(f"  Name: {ethanol_info.get('name')}")
print(f"  Molecular Formula: {ethanol_info.get('molecularFormula')}")
print(f"  Canonical SMILES: {ethanol_info.get('canonicalSmile')}")
print(f"  Status: {ethanol_info.get('status')}")
```

## 5. Exploring Detailed Information

The API returns comprehensive information about each compound. Let's explore what's available:

```{code-cell} ipython3
# Get detailed information for formaldehyde
formaldehyde_info = ccc.cas_to_detail("50-00-0")

print("Formaldehyde - Complete Information:")
print(f"  CAS RN: {formaldehyde_info.get('rn')}")
print(f"  Name: {formaldehyde_info.get('name')}")
print(f"  Molecular Formula: {formaldehyde_info.get('molecularFormula')}")
print(f"  Molecular Mass: {formaldehyde_info.get('molecularMass')}")
print(f"  SMILES: {formaldehyde_info.get('smile')}")
print(f"  Canonical SMILES: {formaldehyde_info.get('canonicalSmile')}")
print(f"  InChI: {formaldehyde_info.get('inchi')}")
print(f"  InChI Key: {formaldehyde_info.get('inchiKey')}")
print(f"  Has Molfile: {formaldehyde_info.get('hasMolfile')}")
print(f"  Number of images: {len(formaldehyde_info.get('images', []))}")
print(f"  Number of experimental properties: {len(formaldehyde_info.get('experimentalProperties', []))}")
print(f"  URI: {formaldehyde_info.get('uri')}")
```

```{code-cell} ipython3
# Explore synonyms
synonyms = formaldehyde_info.get('synonyms', [])
print(f"\nFormaldehyde has {len(synonyms)} synonyms:")
print("First 10 synonyms:")
for i, synonym in enumerate(synonyms[:10], 1):
    print(f"  {i}. {synonym}")
```

```{code-cell} ipython3
# Explore experimental properties
exp_props = formaldehyde_info.get('experimentalProperties', [])
print(f"\nFormaldehyde has {len(exp_props)} experimental properties:")
if exp_props:
    print("First few experimental properties:")
    for i, prop in enumerate(exp_props[:3], 1):
        print(f"  {i}. Property: {prop.get('property', 'N/A')}")
        print(f"     Value: {prop.get('value', 'N/A')}")
        print(f"     Units: {prop.get('units', 'N/A')}")
        print()
```

## 6. Error Handling

The API handles various error conditions gracefully. Let's see what happens with invalid inputs:

```{code-cell} ipython3
# Try an invalid CAS RN
invalid_cas = ccc.cas_to_detail("0000-00-0")
print("Invalid CAS RN (0000-00-0):")
print(f"  Status: {invalid_cas.get('status')}")
print(f"  Name: {invalid_cas.get('name')}")

# Try a non-existent compound name
non_existent = ccc.name_to_detail("thiscompounddoesnotexist12345")
print(f"\nNon-existent compound name:")
print(f"  Status: {non_existent.get('status')}")

# Try an invalid SMILES
invalid_smiles = ccc.smiles_to_detail("INVALID_SMILES")
print(f"\nInvalid SMILES:")
print(f"  Status: {invalid_smiles.get('status')}")
```

## 7. Common Use Cases

Here are some practical examples of how to use the CASCommonChem class:

```{code-cell} ipython3
# Use case 1: Get basic identifiers for a compound
def get_basic_identifiers(cas_rn):
    """Get basic chemical identifiers for a compound"""
    info = ccc.cas_to_detail(cas_rn)
    if info.get('status') == 'Success':
        return {
            'cas_rn': info.get('rn'),
            'name': info.get('name'),
            'formula': info.get('molecularFormula'),
            'mass': info.get('molecularMass'),
            'smiles': info.get('smile'),
            'inchi_key': info.get('inchiKey')
        }
    return None

# Test with benzene
benzene_ids = get_basic_identifiers("71-43-2")
print("Benzene identifiers:")
for key, value in benzene_ids.items():
    print(f"  {key}: {value}")
```

```{code-cell} ipython3
# Use case 2: Find all synonyms for a compound
def get_all_synonyms(cas_rn):
    """Get all synonyms for a compound"""
    info = ccc.cas_to_detail(cas_rn)
    if info.get('status') == 'Success':
        return {
            'name': info.get('name'),
            'cas_rn': info.get('rn'),
            'synonyms': info.get('synonyms', [])
        }
    return None

# Test with glucose
glucose_synonyms = get_all_synonyms("50-99-7")
if glucose_synonyms:
    print(f"Glucose ({glucose_synonyms['cas_rn']}) has {len(glucose_synonyms['synonyms'])} synonyms:")
    print("Sample synonyms:")
    for synonym in glucose_synonyms['synonyms'][:8]:
        print(f"  - {synonym}")
```

```{code-cell} ipython3
# Use case 3: Compare multiple compounds
compounds_to_compare = ["64-17-5", "67-56-1", "78-93-3"]  # Ethanol, Methanol, Butanone

print("Comparison of three compounds:")
print("-" * 70)
for cas_rn in compounds_to_compare:
    info = ccc.cas_to_detail(cas_rn)
    if info.get('status') == 'Success':
        print(f"CAS RN: {cas_rn}")
        print(f"  Name: {info.get('name')}")
        print(f"  Formula: {info.get('molecularFormula')}")
        print(f"  Mass: {info.get('molecularMass')}")
        print(f"  SMILES: {info.get('smile')}")
        print("-" * 70)
```

## Summary

The `CASCommonChem` class provides three main methods:

1. **`cas_to_detail(cas_rn)`**: Look up by CAS Registry Number
2. **`name_to_detail(name)`**: Search by compound name or IUPAC name  
3. **`smiles_to_detail(smiles)`**: Search by SMILES notation

### Key Features:
- ✅ Access to 500,000+ chemical substances
- ✅ Comprehensive chemical data (names, formulas, structures, properties)
- ✅ Multiple search methods (CAS RN, name, SMILES)
- ✅ Robust error handling
- ✅ Rich metadata including synonyms and experimental properties

### Returned Data Includes:
- Basic identifiers (name, CAS RN, molecular formula, mass)
- Structure information (SMILES, InChI, InChI Key)
- Synonyms and alternative names
- Experimental properties
- Images and molecular files (when available)
- Citations and references
