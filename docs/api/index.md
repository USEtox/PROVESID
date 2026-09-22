# API reference

Every page here is generated from the docstrings in `src/provesid/`, and every
example on it is run as a doctest. Everything listed is importable from
`provesid` itself unless the page says otherwise.

## Resolving identifiers

| | |
|---|---|
| [`Search`](search.md) | resolve CAS numbers, names, SMILES, InChIs, InChIKeys, DTXSIDs and formulas across the offline databases |
| [`datasets`](datasets.md) | install, inspect and remove the offline databases |

## Offline databases

| Client | Reads | Dataset |
|---|---|---|
| [`PubChemID`](pubchem_id.md) | PubChem compounds with a CAS number | `pubchem` |
| [`CompToxID`](comptox.md) | EPA CompTox chemicals | `comptox` |
| [`ChebiSDF`](chebi_sdf.md) | ChEBI's SDF release | `chebi` |
| [`CheMBL`](chembl.md) | ChEMBL | `chembl` |
| [`ZeroPM`](zeropm.md) | the ZeroPM global inventory | `zeropm` |
| [`REACHDossierID`](reach.md) | REACH dossiers | ships with the package |
| [`SQLiteClient`](sqlite_clients.md) | what the four SQLite clients share | |

## Online services

| Client | Service |
|---|---|
| [`PubChemAPI`](pubchem.md) | PubChem PUG-REST |
| [`PubChemView`](pubchemview.md) | PubChem PUG-View experimental properties |
| [`NCIChemicalIdentifierResolver`](nci_resolver.md) | NCI/CADD resolver (CACTUS) |
| [`ChEBI`](chebi.md) | ChEBI 2.0 REST |
| [`CASCommonChem`](cascommonchem.md) | CAS Common Chemistry (needs a key) |
| [`OPSIN`, `PYOPSIN`](opsin.md) | OPSIN name-to-structure, online or local |
| [`ClassyFireAPI`](classyfire.md) | ClassyFire (frozen since 2023) |

## Classification

| | |
|---|---|
| [`provesid.taxonomy`](taxonomy.md) | offline ChEBI classification with Chebifier (optional) |

## Infrastructure

| | |
|---|---|
| [`provesid.http`](http.md) | the transport every online client shares |
| [`provesid.cache`](cache.md) | the on-disk cache of online answers |
| [`provesid.config`](config.md) | stored API keys |
| [`provesid.tools`, `provesid.utils`](utilities.md) | the consensus vote behind `Search`, and small helpers |
