# Installing the offline databases

PROVESID answers most questions from five local databases. They are large, so
none of them ships with the package and none is downloaded behind your back:
`Search` runs on whatever is installed and says what is missing, and
[`provesid.datasets`](../api/datasets.md) installs them by name.

| name | client | role | download | on disk |
|---|---|---|---:|---:|
| `pubchem` | [`PubChemID`](../api/pubchem_id.md) | CAS, name, InChIKey and formula lookups; the broadest source | 14.3 GiB | 2.3 GiB |
| `comptox` | [`CompToxID`](../api/comptox.md) | DTXSID lookups, and curated CAS–name pairs | 817 MiB | 1.1 GiB |
| `chebi` | [`ChebiSDF`](../api/chebi_sdf.md) | curated structures, synonyms and ChEBI IDs | 250 MiB | 954 MiB |
| `chembl` | [`CheMBL`](../api/chembl.md) | enrichment: adds ChEMBL IDs to structures already found | 5.7 GiB | 2.4 GiB |
| `zeropm` | [`ZeroPM`](../api/zeropm.md) | regulatory inventories; off unless `Search(use_zeropm=True)` | 439 MiB | 439 MiB |

The first four are what `Search` queries by default. All five together are
about 21.5 GiB to download and 7.2 GiB to keep; see below for why the free
space needed on the way in is larger than either.

## Status, plan, fetch, remove

```python
from provesid import datasets

datasets.status()                       # what is installed, where, and which release
datasets.plan(["pubchem", "chebi"])     # what a fetch would transfer and leave behind
datasets.fetch(["pubchem", "chebi"])    # install them
datasets.remove("chembl")               # reclaim the space
```

On an empty data directory, `plan()` reads:

```text
dataset   action  download installed
pubchem download  14.3 GiB   2.3 GiB
comptox download 816.6 MiB   1.1 GiB
  chebi download 250.0 MiB 954.2 MiB
 chembl download   5.7 GiB   2.4 GiB
 zeropm download 438.7 MiB 438.7 MiB
```

and its `attrs` carry the totals: `total_download_bytes`,
`total_resident_bytes` and `peak_bytes`, the free space needed at the worst
moment of the install. `datasets.human_bytes` formats any of them.

- **`fetch` has no "all" default.** Everything is ~21.5 GiB of transfer and,
  while ChEMBL unpacks, up to ~38 GiB of free disk. That has to be asked for
  by name.
- **Every transfer resumes.** An interrupted download leaves a `.part` file
  that the next `fetch` continues with an HTTP `Range` request. A file is
  checked against the MD5 its server publishes, where it publishes one, and
  against the size the server declared, before it is moved into place, so a
  failed download never replaces a good database.
- **`fetch` skips what is present.** Calling it again on the same list is
  cheap. `force=True` downloads again; for ChEMBL that means the current
  release, which may be newer than the one installed.
- **`remove(..., dry_run=True)`** lists the files it would delete, including
  the partial downloads and build leftovers each dataset can leave behind.

## Where they go

The data directory is the per-user one `platformdirs` picks —
`~/.local/share/provesid` on Linux, `~/Library/Application Support/provesid`
on macOS, `%LOCALAPPDATA%\USEtox\provesid` on Windows — and is shared by every
virtual environment on the machine. Set `PROVESID_DATA_DIR` to put it
elsewhere, such as a larger disk:

```bash
export PROVESID_DATA_DIR=/data/provesid
```

Every client and every `datasets` function also takes `data_dir=`. The data
directory holds data you would have to download again; the
[response cache](caching.md) is a separate directory of things that can be
re-fetched at any time.

## What `Search` does when one is missing

`Search(datasets=...)` decides:

| `datasets=` | when a database is missing |
|---|---|
| `"present"` (default) | run on the ones installed, and log once which are missing and what installing them would cost; `df.attrs["sources_unavailable"]` names them |
| `"auto"` | download what is missing, then run |
| `"required"` | raise `MissingDatasetError` in the constructor, naming the `fetch` call |

```text
MissingDatasetError: 1 dataset(s) missing from /home/me/.local/share/provesid:
  PubChem identifiers (pubchem): 14.3 GiB to download, 2.3 GiB on disk --- CAS, name, InChIKey and formula lookups; the broadest source

Install them with:
  provesid.datasets.fetch('pubchem')
```

`datasets.missing()` answers the same question without constructing anything,
and `datasets.require(names)` raises that error for your own code.

!!! warning "The clients themselves do download"
    Constructing a client directly — `PubChemID()`, `CheMBL()` and the others —
    downloads its database if it is absent, because `auto_download` defaults to
    `True`. Pass `auto_download=False` to get a `FileNotFoundError` instead.
    Only `Search` defaults to not downloading.

## PubChem: built from the FTP site, or a prebuilt copy

`PubChemID` reads `pubchem_id.db`: the ~1.43 million PubChem compounds that
carry a CAS number, with their CAS numbers, title, IUPAC name, formula, masses,
SMILES, InChI, InChIKey and synonyms. By default it is **built** from a dated
monthly snapshot of PubChem's FTP site:

```python
from provesid import PubChemID
from provesid.pubchem_ftp import build_pubchem_id_db, list_releases

PubChemID()                     # source="ftp": build from the newest snapshot
PubChemID(source="zenodo")      # download a 2.2 GiB prebuilt copy instead

list_releases()                 # ['2026-09-01', '2026-08-01', ..., 'current']
build_pubchem_id_db(release="2026-09-01")
```

- The compounds are those PubChem's own `CID-Identifiers.tsv.gz` maps to a CAS
  number, each checked against the CAS check digit.
- Eight files are downloaded, read and deleted one at a time, so the free space
  needed is the finished database plus the largest file (7.4 GB), not the
  whole 15.4 GB. Processing takes about twelve minutes on top of the download.
  `keep_downloads=True` keeps them for a later rebuild.
- `include_inchi=False` skips the 7.4 GB InChI file and computes InChI and
  InChIKey with RDKit instead: less to download, a longer build, and RDKit's
  InChI rather than PubChem's.
- The molecular weight is in no FTP file, so it is computed from the formula.
- The built database records its release and every source file's URL and MD5
  (`PubChemID().provenance()`), and carries PubChem's cross-references to
  ChEBI, ChEMBL, CompTox, EC and UNII (`PubChemID().xrefs(cid)`).

The Zenodo copy is older, finds its CAS numbers by pattern-matching synonyms
(1.25% of them are not valid CAS numbers), and has no provenance or
cross-reference tables. `source=` only decides how a *missing* database is
obtained; one already on disk is opened whichever way it was made.

From a shell: `python scripts/build_pubchem_id_db.py --help`.

## ChEMBL: three routes to the same extract

A ChEMBL release is a 27.7 GiB SQLite database with 74 tables, of which
PROVESID reads eight. `CheMBL` installs an **extract** of those eight tables,
2.4 GiB, and the route decides only what is downloaded to build it:

| `CheMBL(source=...)` | download | free disk at peak | installed |
|---|---:|---:|---:|
| `"sqlite"` (default) | 5.8 GB | 33.4 GiB | 2.4 GiB |
| `"mysql"` | 2.1 GB | ~4.5 GiB | 2.4 GiB |
| `"full"` | 5.8 GB | 33.4 GiB | 27.7 GiB |

- **`"sqlite"`** downloads the SQLite archive, unpacks it, compacts it and
  deletes the full release. If the compaction fails, the full database is
  kept and opened instead, so the download is not wasted.
- **`"mysql"`** reads the eight tables straight out of ChEMBL's plain-text
  MySQL dump as it is decompressed. No MySQL server is involved. Its result is
  the same as the SQLite route's, table for table; `CheMBL.extract_digest()`
  fingerprints an extract so that two can be compared.
- **`"full"`** keeps the whole release, for queries outside PROVESID.

A full release already on disk can be compacted in place:

```python
from provesid import CheMBL

CheMBL.compact(remove_source=True)        # 27.7 GiB -> 2.4 GiB
```

`CheMBL` does not pin a release. On first use it reads the `latest/` listing on
EBI's FTP site and downloads what it advertises; after that, any installed
ChEMBL database is reused, so a new release does not trigger a download.
`datasets.fetch("chembl", force=True)` moves to the newest one.

## CompTox, ChEBI and ZeroPM

- **CompTox** is downloaded as one SQLite file from Zenodo. A name index is
  added to it straight after the download (about 20 s, 290 MiB), so an exact
  name lookup finds synonyms as well as the preferred name. A database
  downloaded before the index existed builds it on its first exact name
  lookup; `CompToxID().build_name_index()` does it at a time of your choosing.
- **ChEBI** is downloaded as a gzipped SDF, expanded, and indexed on first
  use, which takes a few minutes once.
- **ZeroPM** is one SQLite file from the ZeroPM project's repository.
  `Search` leaves it out unless asked: its records are harvested from
  regulatory inventories rather than curated compound by compound, so it is
  noisier than the other four. It is the one source that retrieves misspelled
  names, which is why `Search(preset="recall")` turns it on.
