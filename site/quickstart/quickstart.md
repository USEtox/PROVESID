---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.1
kernelspec:
  display_name: python3
  language: python
  name: python3
---

# Quick Start

This page is written as a MyST notebook so code snippets can be executed and re-validated.

For local validation, run:

```bash
./scripts/validate_docs_local.sh
./scripts/validate_docs_local.sh --execute
```

## 1. Imports

```{code-cell} ipython3
from provesid import (
    PubChemAPI,
    PubChemView,
    NCIChemicalIdentifierResolver,
    CASCommonChem,
    CheMBL,
    ZeroPM,
)
from provesid.pubchem import CompoundProperties

print("PROVESID imports succeeded")
```

## 2. Online lookup with PubChem

```{code-cell} ipython3
pc = PubChemAPI()

cids = pc.get_cids_by_name("aspirin")
print("Top aspirin CID candidates:", cids[:5])

aspirin_cid = cids[0]
basic = pc.get_basic_compound_info(aspirin_cid)
print("CID:", aspirin_cid)
print("Title:", basic.get("Title"))
print("MolecularFormula:", basic.get("MolecularFormula"))
```

## 3. Targeted properties from PubChem

```{code-cell} ipython3
props = pc.get_compound_properties(
    aspirin_cid,
    [
        CompoundProperties.MOLECULAR_WEIGHT,
        CompoundProperties.MOLECULAR_FORMULA,
        CompoundProperties.INCHIKEY,
    ],
    include_synonyms=False,
)

print("MolecularWeight:", props.get("MolecularWeight"))
print("MolecularFormula:", props.get("MolecularFormula"))
print("InChIKey:", props.get("InChIKey"))
```

## 4. Experimental property table with PubChemView

```{code-cell} ipython3
pv = PubChemView()
boiling_table = pv.get_property_table(aspirin_cid, "Boiling Point")

print("Rows:", len(boiling_table))
print(boiling_table.head(3))
```

## 5. Identifier conversion with NCI resolver

```{code-cell} ipython3
resolver = NCIChemicalIdentifierResolver()

smiles = resolver.resolve("aspirin", "smiles")
inchi = resolver.resolve(smiles, "stdinchi")

print("SMILES:", smiles)
print("InChI prefix:", inchi[:20])
```

## 6. CAS Common Chemistry lookup

```{code-cell} ipython3
try:
    ccc = CASCommonChem()
    water = ccc.cas_to_detail("7732-18-5")
    print("Name:", water.get("name"))
    print("Formula:", water.get("molecularFormula"))
    print("CAS:", water.get("rn"))
except Exception as exc:
    print("CAS Common Chemistry example skipped:", exc)
```

## 7. Offline-first classes (local database interfaces)

```{code-cell} ipython3
# These classes are local/offline interfaces that can auto-download datasets.
# auto_download=False lets this quickstart remain lightweight when datasets
# are not yet present on disk.

for cls in (CheMBL, ZeroPM):
    try:
        _ = cls(auto_download=False)
        print(f"{cls.__name__}: local dataset is available")
    except Exception as exc:
        print(f"{cls.__name__}: dataset not yet available ({exc.__class__.__name__})")
```

## Next

- [Online and Offline Data Methods](data_methods.md)
- [PubChem Tutorial](examples/pubchem/pubchem_tutorial.md)
- [API Overview](api/index.md)
