# PubChem API

The PubChem API module provides access to PubChem's REST services for retrieving compound and substance data. This module has been recently enhanced with improved data access patterns and new search methods.

::: provesid.pubchem

## Quick Start

```python
from provesid import PubChemAPI, Domain, CompoundProperties

# Initialize the API client
pc = PubChemAPI()

# Search for compounds by name
cids = pc.get_cids_by_name('aspirin')
print(f"Found CIDs: {cids}")

# Get basic compound information
basic_info = pc.get_basic_compound_info(cids[0])
print(f"Molecular Formula: {basic_info['MolecularFormula']}")
print(f"Molecular Weight: {basic_info['MolecularWeight']}")

# Get compound by CID (improved - no wrapper needed!)
compound = pc.get_compound_by_cid(cids[0])
print(f"Compound keys: {list(compound.keys())}")
```

## Key Improvements

### Elegant Data Access ✨

**Before (redundant wrapper access):**
```python
# Old way required nested access
substance = pc.get_substance_by_sid(sid)
data = substance["PC_Substances"][0]  # Redundant wrapper

compound = pc.get_compound_by_cid(cid)
data = compound["PC_Compounds"][0]    # Redundant wrapper
```

**Now (direct access):**
```python
# New way provides direct access
substance = pc.get_substance_by_sid(sid)  # Direct access!
compound = pc.get_compound_by_cid(cid)    # Direct access!
```

### Enhanced Search Methods

#### Multiple Search Domains
```python
# Search in compound domain (default)
cids = pc.get_cids_by_name('aspirin', domain=Domain.COMPOUND)

# Search in substance domain  
cids = pc.get_cids_by_name('8000-78-0', domain=Domain.SUBSTANCE)

# Comprehensive search across both domains
results = pc.find_cids_comprehensive('8000-78-0')
```

#### Structure-Based Searching
```python
# Search by SMILES (now returns clean list)
smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # aspirin
cids = pc.get_cids_by_smiles(smiles)

# Search by InChI Key (newly implemented)
inchi_key = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"  # aspirin
cids = pc.get_cids_by_inchikey(inchi_key)

# Get compound records directly
compound = pc.get_compounds_by_smiles(smiles)
compound = pc.get_compounds_by_inchikey(inchi_key)
```

## Available Methods

### Compound Search Methods
- `get_cids_by_name()` - Search by compound name
- `get_cids_by_smiles()` - Search by SMILES string  
- `get_cids_by_inchikey()` - Search by InChI Key ✨ *New*
- `find_cids_comprehensive()` - Multi-domain search

### Compound Data Methods
- `get_compound_by_cid()` - Get compound record ✨ *Improved*
- `get_compounds_by_name()` - Get compounds by name ✨ *Improved*
- `get_compounds_by_smiles()` - Get compounds by SMILES ✨ *Improved*
- `get_compounds_by_inchikey()` - Get compounds by InChI Key ✨ *Improved*

### Substance Methods
- `get_substance_by_sid()` - Get substance record ✨ *Improved*
- `get_substances_by_name()` - Get substances by name ✨ *Improved*
- `get_sids_by_name()` - Search for substance IDs

### Property Methods
- `get_basic_compound_info()` - Essential compound properties
- `get_compound_properties()` - Selected properties for one compound, with synonyms
- `get_all_compound_info()` - All available properties
- `get_properties_for_cids()` - Raw property table for many CIDs, in bulk
- `get_compound_properties_batch()` - Bulk retrieval reshaped to one dict per CID

### Utility Methods
- `get_compound_synonyms()` - Get compound synonyms
- `get_compound_identifiers()` - Extract specific identifiers

## Batch Processing

PubChem's property endpoint answers a whole list of CIDs in one round trip, so
asking about a thousand compounds costs a handful of requests rather than a
thousand.

```python
# Process multiple compounds efficiently
compound_names = ["aspirin", "caffeine", "acetaminophen", "ibuprofen"]
all_cids = []

for name in compound_names:
    cids = pc.get_cids_by_name(name)
    if cids:
        all_cids.append(cids[0])

# Batch property retrieval
properties = [CompoundProperties.MOLECULAR_WEIGHT,
              CompoundProperties.MOLECULAR_FORMULA,
              CompoundProperties.SMILES]

batch_results = pc.get_compound_properties_batch(all_cids, properties)
```

`get_compound_properties_batch()` returns one dict per CID, in the order asked,
each carrying the retrieved properties plus `success`, `cid` and `error`. A CID
PubChem has no record of comes back with `success=False` rather than being
dropped, so the result can be zipped against the input.

`get_properties_for_cids()` is the same request without the reshaping — it
returns PubChem's own `PropertyTable` rows, where a compound with no record
yields a row holding only `CID`:

```python
rows = pc.get_properties_for_cids(range(2000, 3000), ["MolecularWeight"])
# 1000 compounds, 5 requests (PROPERTY_CHUNK_SIZE = 200 CIDs each)
```

Both methods switch from a URL path to a POST body once the identifier list
outgrows `URL_IDENTIFIER_LIMIT`, which is how PubChem asks for long lists to be
sent. Neither includes synonyms: those need one request per compound, which
would defeat the batching. Use `get_compound_synonyms()` where they are needed.

## The Local Database

`PubChemID` reads `pubchem_id.db`, a SQLite file holding the ~1.43 M PubChem
compounds that carry a CAS number: their CAS numbers, title, IUPAC name,
formula, molecular weight, exact and monoisotopic mass, isomeric SMILES, InChI,
InChIKey, creation date and synonyms. When the file is missing, `source` decides
where it comes from:

```python
from provesid import PubChemID

db = PubChemID()                    # source="ftp": build it from PubChem's FTP site
db = PubChemID(source="zenodo")     # download a 2.2 GiB prebuilt copy instead
```

`source` describes how a missing database is acquired; a database already on
disk is opened whichever way it was made.

`PubChemID` lives in `provesid.pubchem_id`, apart from the online client, and
is imported from `provesid` like everything else.

::: provesid.pubchem_id

### Building from PubChem's FTP site

`provesid.pubchem_ftp.build_pubchem_id_db()` — what `PubChemID()` calls — builds
the database from a dated monthly snapshot of `Compound/Extras/`:

```python
from provesid.pubchem_ftp import build_pubchem_id_db, list_releases

list_releases()
# ['2026-09-01', '2026-08-01', '2026-07-01', '2026-06-01', 'current']

build_pubchem_id_db(release="2026-09-01")      # or "latest" (default), "current"
```

- **The compounds** are those PubChem's own identifier file,
  `CID-Identifiers.tsv.gz`, maps to a CAS number, each checked against the CAS
  check digit. The old Zenodo database found its CAS numbers by running a
  regular expression over free-text synonyms, and 1.25% of them are not valid
  CAS numbers; this route rejects about 120 rows in 1.46 M.
- **Every file** is downloaded resumably and checked against the MD5 PubChem
  publishes beside it, read once to keep the compounds in scope, and deleted —
  so the free disk needed is the 2.5 GB database plus the largest file
  (7.4 GB), not the 15.4 GB total. `keep_downloads=True` keeps them for a later
  rebuild. Processing takes about 12 minutes; the rest is download time.
- **`include_inchi=False`** skips the 7.4 GB InChI file and computes InChI and
  InChIKey from the SMILES with RDKit: less to download, a longer build, and
  RDKit's InChI rather than PubChem's.
- **The molecular weight** is not in any FTP file, so it is computed from the
  formula (and from the SMILES for isotopically labelled compounds). It agrees
  with PubChem's own value to the second decimal for most compounds; PubChem
  rounds some weights more coarsely, and where it does this is the more precise
  number.

The database records where it came from, and carries the cross-references
PubChem publishes in the same identifier file:

```python
db = PubChemID()
db.provenance()["release"]          # '2026-09-01', plus every file's URL and MD5
db.xrefs(2244)
# {'chebi': ['CHEBI:15365'], 'chembl': ['CHEMBL25'], 'dtxsid': ['DTXSID5020108'],
#  'ec': ['200-064-1'], 'unii': ['R16CO5Y76E']}
```

For a shell, `python scripts/build_pubchem_id_db.py --help`.

::: provesid.pubchem_ftp

## Properties Without the Network

`PubChemID` answers the properties that are *data* about a compound — its
identifiers, names, formula and masses — from disk, so most property lookups
need no request at all. `properties()` reads it first and consults PUG-REST only
for what it cannot answer:

```python
from provesid import PubChemID

db = PubChemID()

db.properties(2244, ["MolecularFormula", "MonoisotopicMass"])
# {'CID': 2244, 'Source': 'offline', 'MolecularFormula': 'C9H8O4',
#  'MonoisotopicMass': 180.04225873}

db.properties(2244, ["XLogP"])
# {'CID': 2244, 'Source': 'online', 'XLogP': 1.2}
```

The `Source` key records which one answered. Two rules decide it:

- a CID with no row in the local database goes online;
- a property with no local column sends the *whole* request online, because that
  property would need a request anyway and a row assembled from two PubChem
  snapshots is worse than a row from one.

`db.offline_properties` lists what the open database can serve: the formula,
molecular weight, exact and monoisotopic mass, isomeric SMILES, InChI,
InChIKey, IUPAC name and title (a Zenodo copy lacks the monoisotopic mass).
Everything else is online-only — including the computed descriptors `XLogP`,
`TPSA`, `Complexity`, `Charge` and the atom and bond counts, which are
PubChem's model outputs rather than data about the compound, and are served
from PubChem even by a Zenodo copy that still stores them.

Bulk lookups split themselves between the two sources automatically, and the
online remainder is fetched in one batched request:

```python
rows = db.properties_for_cids(my_ten_thousand_cids, ["MolecularWeight", "InChIKey"])
table = db.properties_table(my_ten_thousand_cids, ["MolecularWeight"])
```

`properties_table()` returns a DataFrame with a row for every CID asked about,
`Source` reading `offline`, `online` or `missing`, so it can be joined against
your own table. Pass `use_online_fallback=False` to keep a lookup strictly
local — a hard guarantee, not a preference.

A property the compound has no value for is **absent** from the result rather
than `None`, which is how PubChem itself reports it: `'XLogP' not in result`
means PubChem computes no logP for that compound, not that the lookup fell
short. Values are normalised to one type across both sources, since PUG-REST
reports `MolecularWeight` as a string where the local database holds a float.

## Descriptors on Demand

The local database stores no computed descriptors: XLogP, TPSA and the atom
and bond counts are the output of a model run over the structure, and there is
more than one model. `descriptors()` runs one and names it:

```python
from provesid import PubChemID

db = PubChemID()

db.descriptors(2244)                        # RDKit, from the stored SMILES
# {'CID': 2244, 'Source': 'rdkit', 'MolLogP': 1.3101, 'TPSA': 63.6,
#  'HBondDonorCount': 1, 'HBondAcceptorCount': 3, 'RotatableBondCount': 2,
#  'HeavyAtomCount': 13, 'Charge': 0}

db.descriptors(2244, ["XLogP", "Complexity"], source="pubchem")
# {'CID': 2244, 'Source': 'online', 'XLogP': 1.2, 'Complexity': 212.0}
```

- **`source="rdkit"`** (default) needs no network for a compound the database
  holds, and costs about half a millisecond per compound. For one it does not
  hold, only the SMILES is fetched from PubChem; `use_online_fallback=False`
  prevents even that.
- **`source="pubchem"`** is PubChem's own values over PUG-REST, the same path
  as `properties()`. It is the only way to `XLogP` and `Complexity`.

The names are PubChem's wherever the quantity is the same one, so a table can
switch source without renaming its columns. The logP is the exception: RDKit's
is Crippen's model, not XLogP3, so it is called `MolLogP`, and asking RDKit for
`XLogP` raises an error that says so. `Complexity` has no RDKit counterpart.

The numbers differ even where the names agree, because PubChem computes its
descriptors with Cactvs. Against PubChem's values for 20 000 random compounds
in the database:

| descriptor | RDKit equals PubChem |
|---|---:|
| `HeavyAtomCount`, `Charge` | all |
| `HBondDonorCount` | 94% |
| `RotatableBondCount` | 74% |
| `TPSA` | 70% |
| `HBondAcceptorCount` | 63% |
| `MolLogP` vs `XLogP` | 62% within 0.5 log units |

Use one source throughout an analysis. `descriptors_for_cids()` and
`descriptors_table()` are the bulk forms, shaped like their `properties_`
counterparts. For a structure that is not in PubChem at all:

```python
from provesid.pubchem_id import rdkit_descriptors

rdkit_descriptors("CCO", ["TPSA", "MolLogP"])
# {'TPSA': 20.23, 'MolLogP': -0.0014}
```

## Error Handling

```python
from provesid import PubChemNotFoundError, PubChemError

try:
    cids = pc.get_cids_by_name('invalid_compound_name')
    if not cids:
        print("No compounds found")
except PubChemNotFoundError:
    print("Compound not found in PubChem")
except PubChemError as e:
    print(f"PubChem API error: {e}")
```

## Best Practices

- Use batch methods for multiple compounds
- Always handle potential errors with try/except blocks
- Use domain-specific searches when appropriate
- Check data availability before processing
- Respect PubChem's rate limits (built into the API)

## Tutorial

For a comprehensive tutorial with examples, see: [PubChem Tutorial](../examples/pubchem/pubchem_tutorial.md)
