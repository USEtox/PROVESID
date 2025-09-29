# PROVESID

[![Documentation Status](https://github.com/USEtox/PROVESID/actions/workflows/mkdocs-deploy.yml/badge.svg)](https://usetox.github.io/PROVESID/)
[![Tests](https://github.com/USEtox/PROVESID/actions/workflows/test.yml/badge.svg)](https://github.com/USEtox/PROVESID/actions/workflows/test.yml)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

PROVESID is a member of the family of PROVES packages that provides Pythonic access to online services of chemical identifiers and data. The goal is to have a clean interface to the most important online databases with a simple, intuitive (and documented), up-to-date, and extendable interface. We offer interfaces to [PubChem](https://pubchem.ncbi.nlm.nih.gov/), [NCI chemical identifier resolver](https://cactus.nci.nih.gov/chemical/structure), [CAS Common Chemistry](https://commonchemistry.cas.org/), [IUPAC OPSIN](https://www.ebi.ac.uk/opsin/), [ChEBI](https://www.ebi.ac.uk/chebi/beta/), and [ClassyFire](http://classyfire.wishartlab.com/). We highly recommend the new users to jump head-first into [examples folder](./examples/) and get started by playing with the code. We also keep documenting the old and new functionalities [here](https://usetox.github.io/PROVESID/).

# Installation

The package can be installed from PyPi by running

```
pip install provesid
```

To install the latest development version (for developers and enthusiasts), clone or download this repository, for to the root folder and install it by

```
pip install -e .
```

# ✨ New in v0.2.0: Advanced Unlimited Caching

PROVESID v0.2.0 introduces **unlimited caching** across ALL APIs with persistent storage, automatic size monitoring, and import/export capabilities:

```python
import provesid

# ALL APIs now use unlimited caching automatically - no more 512-entry limits!
pubchem_api = provesid.PubChemAPI()          # 19 cached methods
nci_resolver = provesid.NCIChemicalIdentifierResolver()  # 15 cached methods
cas_api = provesid.CASCommonChem()           # 2 cached methods
pugview = provesid.PubChemView()             # 15+ cached methods
classyfire = provesid.ClassyFireAPI()        # 3 cached methods
opsin = provesid.OPSIN()                     # 2 cached methods

# All API calls cached forever, survive restarts
result1 = pubchem_api.get_compound_by_cid(2244)     # PubChem - cached
result2 = nci_resolver.resolve('aspirin', 'smiles') # NCI - cached  
result3 = cas_api.cas_to_detail('50-00-0')          # CAS - cached
result4 = pugview.get_melting_point(2244)           # PubChemView - cached
result5 = opsin.get_id('ethanol')                   # OPSIN - cached

# Unified cache management across all APIs
provesid.export_cache('my_research_cache.pkl')  # Backup all API cache
provesid.import_cache('shared_cache.pkl')       # Load shared team cache
info = provesid.get_cache_info()               # Monitor cache size (296 entries!)
provesid.clear_cache()                         # Clear all caches

# Individual API cache management
pubchem_api.clear_cache()                      # Clear only PubChem cache
cache_stats = pugview.get_cache_info()         # Get detailed cache statistics
```

**Key benefits:**
- 🚀 **Unlimited cache** - No more entry limits for ANY API (was 512 max)
- 💾 **Persistent storage** - Cache survives restarts and reinstalls
- 📊 **Size monitoring** - Warns at 5GB with detailed statistics
- 🔄 **Import/Export** - Share cache files with team (pickle + JSON)
- ⚡ **Zero config** - Just import and use - caching is automatic!
- 🎯 **Complete coverage** - All 6 major APIs fully cached

See [Advanced Caching Guide](./docs/advanced_caching.md) for details.

# Examples

**PubChem**

```python
from provesid.pubchem import PubChemAPI
pc = PubChemAPI()  # Now with unlimited caching!
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

# For Developers

If you're interested in contributing to PROVESID or need to understand the release workflow, please see our comprehensive [Development Guide](./DEVELOPMENT.md) which includes:

- 🛠️ Development setup and environment configuration
- 🚀 Step-by-step release workflow and version management
- 🧪 Testing procedures and guidelines
- 📚 Documentation building and contribution guidelines
- 🔍 Code quality standards and tools
- 🤝 Contribution workflow and pull request guidelines

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