# PubChem ID Database

`pubchem_id.db` is the SQLite database `PubChemID` reads: the ~1.43 M PubChem
compounds that carry a CAS number, with their CAS numbers, names, synonyms,
formula, masses, SMILES, InChI and InChIKey, and cross-references to DSSTox,
ChEBI, ChEMBL, EC and UNII.

## Where it comes from

```python
from provesid import PubChemID

db = PubChemID()                  # builds it from PubChem's FTP site if missing
db = PubChemID(source="zenodo")   # or downloads a prebuilt copy
```

The default route, `provesid.pubchem_ftp.build_pubchem_id_db()`, builds it from
a dated monthly snapshot of `Compound/Extras/` on PubChem's FTP site, and
records the release and every source file's MD5 in the database itself
(`db.provenance()`). `docs/guide/datasets.md` explains the build;
`scripts/build_pubchem_id_db.py` runs it from a shell.

The Zenodo copy (`PubChemID.DEFAULT_DB_URL`) is refreshed by hand from an FTP
build every few months. Copies made before 2026-09 were built from a CSV export
of PubChem's classification browser, with CAS numbers found by a regular
expression over the synonyms; they have the eight computed descriptor columns
listed below and no `provenance` or `xrefs` tables.

## Tables

| table | rows (2026-09-01) | |
|---|---:|---|
| `compounds` | see `db.get_stats()` | one row per CID |
| `cas_numbers` | | `(cid, cas)`, one-to-many |
| `synonyms` | | `(cid, synonym)`, in PubChem's order, best first |
| `xrefs` | | `(cid, source, identifier)`; source is `dtxsid`, `chebi`, `chembl`, `ec` or `unii` |
| `provenance` | | `(key, value)`: release, build time, row counts |
| `provenance_files` | | one row per source file: URL, MD5, bytes, lines read and kept |

`compounds` columns: `cid`, `cmpdname` (PubChem's title), `mf`, `inchi`,
`smiles` (isomeric), `inchikey`, `iupacname`, `mw`, `exactmass`,
`monoisotopicmass`, `cidcdate` (`YYYY-MM-DD`).

`mw` is computed from the formula with `provesid.pubchem_ftp.ATOMIC_WEIGHTS`,
and from the SMILES for isotopically labelled compounds, because no FTP file
carries PubChem's own molecular weight.

A Zenodo copy made before 2026-09 has instead `polararea`, `complexity`,
`xlogp`, `heavycnt`, `hbonddonor`, `hbondacc`, `rotbonds` and `charge`, no
`monoisotopicmass`, and `cidcdate` as `YYYYMMDD`. `PubChemID.properties()`
never serves those eight descriptors from disk; it asks PubChem.

Indexes: `compounds(inchikey)`, `compounds(inchi)`, `compounds(mf)`,
`cas_numbers(cas)`, `cas_numbers(cid)`, `synonyms(synonym)`, `synonyms(cid)`,
`xrefs(cid)`, `xrefs(identifier)`.

## Notes

- The database is excluded from version control (see `.gitignore`).
- Only identifiers, names and structures are included, not annotations or
  bioassay data.
