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

# Quick start

A tour of the package, offline first. This page is a MyST notebook: run
`./scripts/validate_docs_local.sh --execute` to execute it.

## 1. Which databases are installed

The offline sources read five bulk datasets. None is downloaded on your
behalf: `Search` uses whatever is on disk and reports the rest.

```{code-cell} ipython3
from provesid import datasets

status = datasets.status()
print(status[["dataset", "present", "size", "release"]].to_string(index=False))
print("in", status.attrs["data_dir"])
```

`datasets.plan()` says what a download would transfer before it starts, and
`datasets.fetch()` installs by name:

```{code-cell} ipython3
todo = datasets.plan(["pubchem", "chebi"])
print(todo[["dataset", "action", "download", "installed"]].to_string(index=False))
print("transfer:", datasets.human_bytes(todo.attrs["total_download_bytes"]))

# datasets.fetch(["pubchem", "chebi"])   # resumable, verified, skips what is present
# datasets.remove("chembl")              # reclaim the space, by name
```

See [Installing the offline databases](guide/datasets.md) for what each one
holds and the different ways of building PubChem and ChEMBL.

## 2. Resolve identifiers with `Search`

```{code-cell} ipython3
from provesid import Search

with Search("cas", show_progress=False) as s:
    df = s.search(["50-00-0", "64-17-5", "1912-24-9"])

print(df[["query", "name", "canonical_smiles", "n_source_support", "confidence"]].to_string())
print("sources used:", df.attrs["sources_available"])
```

The same class resolves names, SMILES, InChIs, InChIKeys, DTXSIDs and
formulas. Presets trade precision for recall:

```{code-cell} ipython3
with Search("name", preset="recall", n_hits=1, show_progress=False) as s:
    print(s.search(["asprin", "caffiene"])[["query", "name", "confidence"]].to_string())
```

With `online_fallback=True`, a query no installed database answers is asked of
PubChem and the NCI resolver. See
[Resolving identifiers with Search](guide/search.md).

## 3. Use a database directly

Each database has its own client. They are context managers, and each thread
that queries one gets its own connection:

```{code-cell} ipython3
from provesid import PubChemID

try:
    with PubChemID(auto_download=False) as db:
        print(db.cas_to_inchi("50-78-2"))
        print(db.properties(2244, ["MolecularFormula", "InChIKey"],
                            use_online_fallback=False))
        print(db.descriptors(2244, ["TPSA", "MolLogP"]))
except FileNotFoundError:
    print("pubchem_id.db is not installed")
```

`auto_download=False` matters here: a client constructed directly downloads
its database if it is missing. See
[Using the local databases directly](guide/local-databases.md).

## 4. Online services

The online clients share one transport that paces and retries requests, and
cache their answers on disk.

```{code-cell} ipython3
from provesid import PubChemAPI, PubChemView, NCIChemicalIdentifierResolver

pc = PubChemAPI()
cid = pc.get_cids_by_name("aspirin")[0]
print("CID:", cid, pc.get_basic_compound_info(cid).get("MolecularFormula"))

melting = PubChemView().get_property_table(cid, "Melting Point")
print(melting[["StringWithMarkup", "ValueSI", "UnitSI"]].head(3))

resolver = NCIChemicalIdentifierResolver()
print("SMILES:", resolver.resolve("aspirin", "smiles"))
```

`CASCommonChem` needs an API key; see [API keys](guide/api-keys.md).

## Next

- [Installing the offline databases](guide/datasets.md)
- [Resolving identifiers with Search](guide/search.md)
- [PubChem tutorial](examples/pubchem/pubchem_tutorial.md)
- [API reference](api/index.md)
