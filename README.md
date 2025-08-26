# PROVESID

[![Documentation Status](https://github.com/USEtox/PROVESID/actions/workflows/mkdocs-deploy.yml/badge.svg)](https://usetox.github.io/PROVESID/)
[![Tests](https://github.com/USEtox/PROVESID/actions/workflows/test.yml/badge.svg)](https://github.com/USEtox/PROVESID/actions/workflows/test.yml)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

PROVESID is a member of the family of PROVES packages that provides Pythonic access to online services of chemical identifiers and data. The goal is to have a clean interface to the most important online databases with a simple, intuitive (and documented), up-to-date, and extendable interface. We offer interfaces to [PubChem](https://pubchem.ncbi.nlm.nih.gov/), [NCI chemical identifier resolver](https://cactus.nci.nih.gov/chemical/structure), [CAS Common Chemistry](https://commonchemistry.cas.org/), [IUPAC OPSIN](https://www.ebi.ac.uk/opsin/), [ChEBI](https://www.ebi.ac.uk/chebi/beta/), and [ClassyFire](http://classyfire.wishartlab.com/). We highly recommend the new users to jump head-first into [examples folder](./examples/) and get started by playing with the code. We also keep documenting the old and new functionalities [here]().

# Examples

**PubChem**

```python
from provesid.pubchem import PubChemAPI
pc = PubChemAPI()
cids_aspirin = pc.get_cids_by_name('aspirin')
res_basic = pc.get_basic_compound_info(cids_aspirin[0])
```

which returns

```python
{
  "CID": 2244,
  "MolecularFormula": "C9H8O4",
  "MolecularWeight": "180.16",
  "SMILES": "CC(=O)OC1=CC=CC=C1C(=O)O",
  "InChI": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
  "InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
  "IUPACName": "2-acetyloxybenzoic acid",
  "success": true,
  "cid": 2244,
  "error": null
}
```

**PubChem View for data**

```python
from provesid import PubChemView, get_property_table
logp_table = get_property_table(cids_aspirin[0], "LogP")
logp_table
```

which returns a table with the reported values of `logP` for aspirin (including the references for each data point).

**Chemical Identifier Resolver**

```python
from provesid import NCIChemicalIdentifierResolver
resolver = NCIChemicalIdentifierResolver()
smiles = resolver.resolve(compound, 'smiles')
```

**OPSIN**

```python
from provesid import OPSIN
opsin = OPSIN()
methane_result = opsin.get_id("methane")
```

which returns:

```python
{'status': 'SUCCESS',
 'message': '',
 'inchi': 'InChI=1/CH4/h1H4',
 'stdinchi': 'InChI=1S/CH4/h1H4',
 'stdinchikey': 'VNWKTOKETHGBQD-UHFFFAOYSA-N',
 'smiles': 'C'}
 ```

**CAS Common Chemistry**

```python
from provesid import CASCommonChem
ccc = CASCommonChem()
water_info = ccc.cas_to_detail("7732-18-5")
print("Water (7732-18-5):")
print(f"  Name: {water_info.get('name')}")
print(f"  Molecular Formula: {water_info.get('molecularFormula')}")
print(f"  Molecular Mass: {water_info.get('molecularMass')}")
print(f"  SMILES: {water_info.get('smile')}")
print(f"  InChI: {water_info.get('inchi')}")
print(f"  Status: {water_info.get('status')}")
```

which returns

```
Water (7732-18-5):
  Name: Water
  Molecular Formula: H<sub>2</sub>O
  Molecular Mass: 18.02
  SMILES: O
  InChI: InChI=1S/H2O/h1H2
  Status: Success
```

**ClassyFire**

See the [tutorial notebook](./examples/ClassyFire/classyfire_tutorial.ipynb).

# Other tools

Several other Python (and other) packages and sample codes are available. We are inspired by them and tried to improve upon them based on our personal experiences working with chemical identifiers and data.  

  - [PubChemPy](https://github.com/mcs07/PubChemPy) and [docs](https://docs.pubchempy.org/en/latest/)  
  - [CIRpy](https://github.com/mcs07/CIRpy) and [docs](https://cirpy.readthedocs.io/en/latest/)  
  - [IUPAC cookbook](https://iupac.github.io/WFChemCookbook/intro.html) for a tutorial on using various web APIs.  
  - more?

# TODO list

We will provide Python interfaces to more online services, including:  

  - [ZeroPM](https://database.zeropm.eu/) even though there is no web API, the data is available on GitHub. I have written an interface that is not shared here since it can make this codebase too large, and I aim to keep it lean. We will find a way to share it.  
  - More? Please [open an issue](https://github.com/USEtox/PROVESID/issues) and let us know what else you would like to have included.