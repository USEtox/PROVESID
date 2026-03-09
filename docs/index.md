# PROVESID

PROVESID is a Python package for chemical identifier resolution, property retrieval, and local dataset lookup across multiple chemistry data systems.

## What PROVESID Covers

PROVESID includes interfaces for online APIs and local/offline databases.

### Online interfaces

- `PubChemAPI` for PUG-REST queries (compound, substance, assay, and identifier workflows)
- `PubChemView` for experimental property extraction with references and tabular outputs
- `NCIChemicalIdentifierResolver` for identifier conversion through NCI resolver endpoints
- `CASCommonChem` for CAS Common Chemistry lookups
- `ChEBI` for ChEBI API access
- `OPSIN` and `PYOPSIN` for chemical name-to-structure conversion
- `ClassyFireAPI` for chemical taxonomy/classification

### Offline and local dataset interfaces

- `CheMBL` local SQLite interface (auto-download supported)
- `PubChemID` local SQLite identifier database (CID/CAS/InChI/InChIKey and metadata)
- `CompToxID` local SQLite interface (auto-download supported)
- `ZeroPM` local SQLite interface (auto-download supported)
- `REACHDossierID` local REACH dossier dataset lookup
- `ChebiSDF` local ChEBI SDF parser

## Recommended Installation Method

Use `uv` as the primary installation workflow:

```bash
uv pip install provesid
```

For development from source:

```bash
git clone https://github.com/USEtox/PROVESID.git
cd PROVESID
uv pip install -e .
```

`uv` is recommended because PROVESID can work with large local data files and database assets. Using `uv` helps avoid repeated package/data copies across many environments.

## Start Here

- [Quick Start](quickstart.md)
- [Online and Offline Data Methods](data_methods.md)
- [Advanced Caching](advanced_caching.md)

## Tutorials

- [PubChem Tutorial](examples/pubchem/pubchem_tutorial.md)
- [CAS Common Chemistry Tutorial](examples/CCC/CAS_Common_Chemistry_tutorial.md)
- [ChEBI Tutorial](examples/ChEBI/ChEBI_tutorial.md)
- [ChEBI SDF Tutorial](examples/ChEBI/chebi_sdf_tutorial.md)
- [ClassyFire Tutorial](examples/ClassyFire/classyfire_tutorial.md)
- [OPSIN Tutorial](examples/OPSIN/opsin_tutorial.md)
- [Chemical ID Resolver Tutorial](examples/resolver/chem_id_resolver_tutorial.md)
- [ChEMBL Tutorial](examples/chembl/chembl_tutorial.md)
- [PubChem View Tutorial](examples/pubchemview/pubchem_view_tutorial.md)
- [ZeroPM Tutorial](examples/zeropm/zeropm-example.md)

## API Reference

- [API Overview](api/index.md)
- [PubChem API](api/pubchem.md)
- [PubChem View](api/pubchemview.md)
- [NCI Resolver](api/nci_resolver.md)
- [CAS Common Chemistry](api/cascommonchem.md)
- [ChEBI](api/chebi.md)
- [ClassyFire](api/classyfire.md)
- [OPSIN](api/opsin.md)
- [ChEMBL](api/chembl.md)
