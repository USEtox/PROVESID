# PROVESID

[![Documentation](https://github.com/USEtox/PROVESID/actions/workflows/mkdocs-deploy.yml/badge.svg)](https://usetox.github.io/PROVESID/)
[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

PROVESID resolves chemical identifiers and retrieves chemical data, **offline
first**. It keeps local copies of PubChem, EPA CompTox, ChEBI, ChEMBL and
ZeroPM, answers from them, and asks the online services only when you let it.
It is part of **PROVES**, a family of packages for pre**PRO**cessing and
**VE**rification of **S**ubstance data.

```python
from provesid import Search

with Search("cas") as s:
    df = s.search(["50-00-0", "64-17-5", "1912-24-9"])
df[["query", "name", "canonical_smiles", "InChIKey", "n_source_support", "confidence"]]
```

```text
       query          name         canonical_smiles                     InChIKey  n_source_support  confidence
0    50-00-0  formaldehyde                      C=O  WSFSSNUMVMOOMR-UHFFFAOYSA-N                 3      0.9000
1    64-17-5       ethanol                      CCO  LFQSCWFLJHTTHZ-UHFFFAOYSA-N                 4      0.8906
2  1912-24-9      atrazine  CCNc1nc(Cl)nc(NC(C)C)n1  MXWJVTOOROXGIU-UHFFFAOYSA-N                 4      0.9000
```

`Search` asks every installed database about each identifier, keeps the
structure they agree on, and reports how many sources agreed. The same class
takes names, SMILES, InChIs, InChIKeys, DTXSIDs and formulas. No request
leaves the machine unless you pass `online_fallback=True`.

## Install

```bash
uv pip install provesid                                # from PyPI
uv pip install git+https://github.com/USEtox/PROVESID  # the development version
```

`pip` works as well. Python 3.12 or later is required.

## The offline databases

The package is small, but the databases are large. None of them ships with
the package, and none is downloaded until you ask for it by name:

```python
from provesid import datasets

datasets.status()                                # what is installed, and where
datasets.plan(["pubchem", "comptox", "chebi"])   # what a download would cost
datasets.fetch(["pubchem", "comptox", "chebi"])  # install them
datasets.remove("chembl")                        # reclaim the space
```

| name | client | role | download | on disk |
|---|---|---|---:|---:|
| `pubchem` | `PubChemID` | CAS, name, InChIKey and formula lookups; the broadest source | 14.3 GiB | 2.3 GiB |
| `comptox` | `CompToxID` | DTXSID lookups, and curated CAS–name pairs | 817 MiB | 1.1 GiB |
| `chebi` | `ChebiSDF` | curated structures, synonyms and ChEBI IDs | 250 MiB | 954 MiB |
| `chembl` | `CheMBL` | adds ChEMBL IDs to structures already found | 5.7 GiB | 2.4 GiB |
| `zeropm` | `ZeroPM` | regulatory inventories, persistence and mobility; off in `Search` unless `use_zeropm=True` | 439 MiB | 439 MiB |

All five are about 21.5 GiB to download and 7.2 GiB to keep. While ChEMBL
unpacks, the install needs up to about 38 GiB of free disk space.
Downloads resume after an interruption, and each file is checked before it
replaces an existing one.

The databases go into one per-user directory shared by every virtual
environment on the machine (`~/.local/share/provesid` on Linux). Set
`PROVESID_DATA_DIR`, or pass `data_dir=` to any client, to put them elsewhere.
`Search` uses whatever is installed and reports which sources it used in
`df.attrs["sources_available"]`.

Each database can also be used directly:

```python
from provesid import PubChemID

with PubChemID(auto_download=False) as db:
    db.cas_to_inchi("50-78-2")
    db.properties(2244, ["MolecularFormula", "InChIKey"], use_online_fallback=False)
    db.descriptors(2244, ["TPSA", "MolLogP"])  # RDKit descriptors, computed locally
```

## Online services

| client | service |
|---|---|
| `PubChemAPI` | [PubChem PUG-REST](https://pubchem.ncbi.nlm.nih.gov/) |
| `PubChemView` | PubChem PUG-View: experimental properties, with values parsed into numbers and SI units |
| `NCIChemicalIdentifierResolver` | [NCI/CADD Chemical Identifier Resolver](https://cactus.nci.nih.gov/chemical/structure) |
| `ChEBI` | [ChEBI](https://www.ebi.ac.uk/chebi/) web service |
| `CASCommonChem` | [CAS Common Chemistry](https://commonchemistry.cas.org/); needs an API key |
| `OPSIN` | [OPSIN](https://www.ebi.ac.uk/opsin/) name-to-structure; `PYOPSIN` runs it locally, with Java |

```python
from provesid import PubChemAPI, PubChemView, NCIChemicalIdentifierResolver

pc = PubChemAPI()
cid = pc.get_cids_by_name("aspirin")[0]                       # 2244
melting = PubChemView().get_property_table(cid, "Melting Point")
smiles = NCIChemicalIdentifierResolver().resolve("50-00-0", "smiles")  # "C=O"
```

The online clients share one transport. It paces requests per host, retries
what is worth retrying, and stops asking a host that has said to wait. Their
answers are cached on disk under `~/.cache/provesid/`, or `PROVESID_CACHE_DIR`.

The CAS Common Chemistry key is stored once and picked up by every later
`CASCommonChem()`:

```python
from provesid import set_cas_api_key, CASCommonChem

set_cas_api_key("your-cas-api-key")
CASCommonChem().cas_to_detail("7732-18-5")["name"]  # "Water"
```

`ClassyFireAPI` is still in the package, but the ClassyFire service has not
classified a new structure since February 2023. For ChEBI chemical classes
computed offline, install the `chebifier` extra; see the
[Chebifier guide](https://usetox.github.io/PROVESID/guide/chebifier/).

## Documentation and tutorials

The [documentation](https://usetox.github.io/PROVESID/) has a quick start,
guides and an API reference generated from the docstrings. The tutorials are
executed notebooks in [`examples/`](./examples/):

- [Resolving a dataset with `Search`](./examples/search/search_tutorial.ipynb), the place to start
- [PubChem](./examples/pubchem/pubchem_tutorial.ipynb) and [PubChem View](./examples/pubchemview/pubchem_view_tutorial.ipynb)
- [ChEMBL](./examples/chembl/chembl_tutorial.ipynb)
- [ChEBI](./examples/ChEBI/ChEBI_tutorial.ipynb) and [ChEBI SDF](./examples/ChEBI/chebi_sdf_tutorial.ipynb)
- [ZeroPM](./examples/zeropm/zeropm-example.ipynb)
- [Chemical Identifier Resolver](./examples/resolver/chem_id_resolver_tutorial.ipynb)
- [CAS Common Chemistry](./examples/CCC/CAS_Common_Chemistry_tutorial.ipynb)
- [OPSIN](./examples/OPSIN/opsin_tutorial.ipynb)

Shorter scripts sit beside them, one folder per feature.

## Related tools

PROVESID learned from these packages and resources:

- [PubChemPy](https://github.com/mcs07/PubChemPy) ([docs](https://docs.pubchempy.org/en/latest/))
- [CIRpy](https://github.com/mcs07/CIRpy) ([docs](https://cirpy.readthedocs.io/en/latest/))
- the [IUPAC FAIR Chemistry Cookbook](https://iupac.github.io/WFChemCookbook/intro.html), for tutorials on chemistry web APIs

## Planned

- [UniChem](https://www.ebi.ac.uk/unichem/api/docs) cross-references.
- The [ChEBI ontology](https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/), read with [pronto](https://github.com/althonos/pronto).
- Structure standardisation with the [ChEMBL Structure Pipeline](https://github.com/chembl/ChEMBL_Structure_Pipeline); this may go to `IMPROVES` instead.

Please [open an issue](https://github.com/USEtox/PROVESID/issues) to suggest
another source or to report a problem.
