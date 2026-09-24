# Dataset examples

`provesid.datasets` is two things: the one downloader every bulk dataset in
PROVESID goes through — ChEMBL's 5.8 GB archive, PubChem's 2.2 GB identifier
database, CompTox, ZeroPM and the ChEBI SDF — and the manager that decides
which of them are on this machine.

## Files

- `dataset_manager_demo.py` - what is installed, what a download would cost,
  and the three `Search(datasets=...)` policies. Downloads nothing.
- `resumable_download_demo.py` - resumption, checksums and rejection, shown
  against a local server that drops the connection on purpose. Downloads
  nothing large and needs no network.

## Deciding what to install

```python
from provesid import datasets

datasets.status()                    # what is on disk, and what it occupies
datasets.plan(["pubchem", "chebi"])  # what a download would transfer
datasets.fetch("pubchem")            # install one, resumable and verified
datasets.remove("chembl")            # reclaim the space, by name
```

`Search("cas").search("50-00-0")` on a clean machine used to fetch about 32 GB
without asking, because each source client defaults to `auto_download=True`.
It no longer does: `Search(datasets="present")` is the default and uses
whatever is installed, `"auto"` restores the old behaviour, and `"required"`
raises in the constructor naming the exact `fetch` call.

## Why it exists

Each of those five modules used to carry its own copy of "stream the response
into a temporary file with a progress bar". None retried, none resumed, and
none checked a checksum, so an interrupted 5.8 GB download started again from
zero. Three of them also renamed the file into place *before* checking it, so a
failed download replaced a working database with a broken one.

## What `download_file` does

```python
from provesid.datasets import download_file

download_file(
    "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz",
    "/data/CID-SMILES.gz",
    checksum_url="https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz.md5",
)
```

It streams into `<dest>.part`, resumes that file with an HTTP `Range` request
if the transfer is interrupted, and then applies four checks in order — each of
which leaves the destination untouched if it fails:

1. the byte count against the size the server declared;
2. the MD5, when one is available;
3. your own `verify` callback, handed the finished file;
4. an atomic rename.

A file rejected by its checksum or by `verify` is deleted: it is already
complete, so resuming it would fail the same check again. A transfer that
merely stopped is kept, and the next call continues from it — but only if it
came from the same URL, which is recorded in a `.part.source` marker beside it.
A partial left by some other download is discarded rather than spliced onto
this one.

## Rejecting a bad file before it lands

```python
import sqlite3

def must_be_a_database(path):
    sqlite3.connect(path).execute("SELECT COUNT(*) FROM compounds")

download_file(url, "/data/pubchem_id.db", verify=must_be_a_database)
```

This is what each client now passes: the query it needs the file to be able to
answer. A file that cannot answer it never reaches the destination.
