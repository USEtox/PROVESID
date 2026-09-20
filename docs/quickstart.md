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

Those are one request per compound. For many compounds at once use
`pc.get_compound_properties_batch(cids, properties)`, which asks PubChem about
200 CIDs per request, and for most properties you can skip the network
altogether: `PubChemID().properties(cid, properties)` reads the local database
first and only falls back online for what it cannot answer. See
[the PubChem API page](api/pubchem.md#properties-without-the-network).

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

## 8. Which datasets are installed, and what they cost

The offline sources read five bulk datasets, together about 6.7 GiB installed
— ChEMBL is compacted to the eight tables PROVESID reads as the last step of
its download, so it costs 2.4 GiB rather than 27.7 GiB. None of them is
downloaded on your behalf: `Search` uses whatever is on disk and reports the
rest.

```{code-cell} ipython3
from provesid import datasets

status = datasets.status()
print(status[["dataset", "present", "size", "release"]].to_string(index=False))
print("in", status.attrs["data_dir"])
```

`datasets.plan()` says what a download would transfer before it starts, and
`datasets.fetch()` installs one by name:

```{code-cell} ipython3
todo = datasets.plan(["pubchem", "chebi"])
print(todo[["dataset", "action", "download", "installed"]].to_string(index=False))
print("transfer:", datasets.human_bytes(todo.attrs["total_download_bytes"]))
# attrs["peak_bytes"] is the free disk needed at the worst moment, which for
# ChEMBL is far more than it installs: 33.4 GiB to leave 2.4 GiB behind.

# datasets.fetch(["pubchem", "chebi"])   # resumable, verified, skips what is present
# datasets.remove("chembl")              # reclaim the space, by name
```

`Search(datasets=...)` chooses what happens when one is absent: `"present"`
(the default) runs on the installed sources, `"auto"` downloads what is
missing, and `"required"` raises and names the `fetch` call.

## 9. Closing a database, and querying one from several threads

`PubChemID`, `CompToxID`, `ZeroPM` and `CheMBL` hold a local SQLite file open.
Use them as context managers, or call `close()`, so the file is released when
you are done with it — which is what lets a later download replace it:

```{code-cell} ipython3
from provesid import PubChemID

try:
    with PubChemID(auto_download=False) as db:
        print(db.cas_to_inchi("50-78-2"))
    print("closed:", db.closed)
except FileNotFoundError:
    print("pubchem_id.db is not installed")
```

Each thread gets its own connection, so a pool over a list of identifiers is
the ordinary thing to write:

```python
from concurrent.futures import ThreadPoolExecutor

with PubChemID() as db:
    with ThreadPoolExecutor(8) as pool:
        rows = list(pool.map(db.get_by_cas, cas_numbers))
```

`Search` is a context manager too, and closes the clients it constructed —
though not one you passed to it yourself, which stays yours:

```python
with Search("cas") as s:
    df = s.search(["50-00-0", "64-17-5"])
```

## Next

- [Online and Offline Data Methods](data_methods.md)
- [PubChem Tutorial](examples/pubchem/pubchem_tutorial.md)
- [API Overview](api/index.md)
