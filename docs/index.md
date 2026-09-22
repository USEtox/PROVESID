# PROVESID

PROVESID resolves chemical identifiers and retrieves chemical data, **offline
first**. It keeps local copies of PubChem, EPA CompTox, ChEBI, ChEMBL and
ZeroPM, answers from them, and asks the online services only when you let it.

```python
from provesid import Search

df = Search("cas").search(["50-00-0", "64-17-5", "1912-24-9"])
df[["query", "name", "canonical_smiles", "InChIKey", "confidence"]]
```

```text
       query          name         canonical_smiles                     InChIKey  confidence
0    50-00-0  formaldehyde                      C=O  WSFSSNUMVMOOMR-UHFFFAOYSA-N      0.9000
1    64-17-5       ethanol                      CCO  LFQSCWFLJHTTHZ-UHFFFAOYSA-N      0.8906
2  1912-24-9      atrazine  CCNc1nc(Cl)nc(NC(C)C)n1  MXWJVTOOROXGIU-UHFFFAOYSA-N      0.9000
```

`Search` asks every installed database about each identifier, keeps the
structure they agree on, and reports how many agreed. No request leaves the
machine unless you pass `online_fallback=True`.

## Install

```bash
uv pip install provesid
```

or, from source, `uv pip install -e .` in a clone of
[the repository](https://github.com/USEtox/PROVESID). The package itself is
small. The databases are not, and none is downloaded until you ask:

```python
from provesid import datasets

datasets.plan()                                  # what each would cost
datasets.fetch(["pubchem", "comptox", "chebi"])  # install by name
```

They go into one per-user directory shared by every virtual environment on the
machine. See [Installing the offline databases](guide/datasets.md).

## What is in it

**Offline.** `Search` over `PubChemID`, `CompToxID`, `ChebiSDF` and `CheMBL`,
with `ZeroPM` on request; each client can also be used directly. PubChem
properties and RDKit descriptors without the network. `REACHDossierID` for
REACH dossiers, and optional ChEBI classification with Chebifier.

**Online.** `PubChemAPI` (PUG-REST), `PubChemView` (experimental properties,
parsed into numbers with SI units), `NCIChemicalIdentifierResolver`, `ChEBI`,
`CASCommonChem`, `OPSIN` and `ClassyFireAPI`. They share one transport that
paces requests per host, retries what is worth retrying, and stops asking a
host that has said to wait; their answers are cached on disk.

## Where to go next

- [Quick start](quickstart.md) — a tour, offline first.
- Guides — [datasets](guide/datasets.md), [`Search`](guide/search.md),
  [the local databases](guide/local-databases.md),
  [experimental properties](guide/experimental-properties.md),
  [network behaviour](guide/network.md), [caching](guide/caching.md),
  [API keys](guide/api-keys.md), [Chebifier](guide/chebifier.md).
- Tutorials — one per service, from the repository's `examples/` folder.
- [API reference](api/index.md) — generated from the docstrings.
