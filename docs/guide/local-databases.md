# Using the local databases directly

`Search` is the usual way in, but each database has its own client, with
lookups `Search` does not expose. This page covers what they share and the
few things particular to each. Installing them is covered in
[Installing the offline databases](datasets.md).

## Closing, and `with`

`PubChemID`, `CompToxID`, `ZeroPM` and `CheMBL` hold a SQLite file open. They
share their connection handling ([`SQLiteClient`](../api/sqlite_clients.md)),
so all four close the same way:

```python
from provesid import PubChemID

with PubChemID() as db:
    inchi = db.cas_to_inchi("50-78-2")
# the file is no longer held open here

db = PubChemID()
db.close()          # idempotent
db.closed           # True
```

Closing matters on Windows in particular, where an open file cannot be
replaced: a download cannot overwrite a database a live object still holds.
Querying a closed client raises `DatabaseClosedError`, which names the client
and the file.

## Threads

Each thread gets its own connection, opened on that thread's first query and
closed with the client, so a thread pool over a list of identifiers is the
ordinary thing to write:

```python
from concurrent.futures import ThreadPoolExecutor

with PubChemID() as db:
    with ThreadPoolExecutor(8) as pool:
        rows = list(pool.map(db.get_by_cas, cas_numbers))
```

`close()` closes every thread's connection, from whichever thread calls it. A
query still running on another thread when it does will fail: closing a client
mid-query is a caller error.

A pool is worth it when each item does something *else* as well — a web
lookup, an RDKit canonicalisation, a file read. Measured on 400 CAS numbers,
each a `PubChemID` lookup followed by 20 ms of waiting: 8.52 s serially,
1.17 s on eight threads. A tight loop of nothing but local lookups is the
opposite: these queries take tens of microseconds, far less than the GIL
handoff around each one, so 5 000 lookups that take 0.29 s serially take
4.80 s on four threads. That is `sqlite3` under CPython, and the remedy is not
to use a pool for that shape of work.

## PubChem properties without the network

`PubChemID` answers the properties that are *data* about a compound — its
identifiers, names, formula and masses — from disk. `properties()` reads it
first and asks PUG-REST only for what it cannot answer:

```python
from provesid import PubChemID

db = PubChemID()

db.properties(2244, ["MolecularFormula", "MonoisotopicMass"])   # an FTP-built database
# {'CID': 2244, 'Source': 'offline', 'MolecularFormula': 'C9H8O4',
#  'MonoisotopicMass': 180.04225873}

db.properties(2244, ["XLogP"])
# {'CID': 2244, 'Source': 'online', 'XLogP': 1.2}
```

The `Source` key records which one answered. Two rules decide it:

- a CID with no row in the local database goes online;
- a property with no local column sends the *whole* request online, because
  that property needs a request anyway, and a row assembled from two PubChem
  snapshots is worse than a row from one.

`db.offline_properties` lists what the open database can serve: the formula,
molecular weight, exact and monoisotopic mass, isomeric SMILES, InChI,
InChIKey, IUPAC name and title (a Zenodo copy lacks the monoisotopic mass).
Everything else is online-only, including PubChem's computed descriptors
`XLogP`, `TPSA`, `Complexity`, `Charge` and the atom and bond counts.

The bulk forms split themselves between the two sources, and fetch the online
remainder in batched requests:

```python
rows = db.properties_for_cids(cids, ["MolecularWeight", "InChIKey"])
table = db.properties_table(cids, ["MolecularWeight"])
```

`properties_table()` returns a DataFrame with a row for every CID asked about,
`Source` reading `offline`, `online` or `missing`. Pass
`use_online_fallback=False` to keep a lookup strictly local — a guarantee, not
a preference.

A property the compound has no value for is **absent** from the result rather
than `None`, as PubChem itself reports it: `'XLogP' not in result` means
PubChem computes no logP for that compound.

## Descriptors on demand

The local database stores no computed descriptors: XLogP, TPSA and the atom
and bond counts are the output of a model run over the structure, and there is
more than one model. `descriptors()` runs one and says which:

```python
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
- **`source="pubchem"`** is PubChem's own values, over PUG-REST. It is the only
  way to `XLogP` and `Complexity`.

The names are PubChem's wherever the quantity is the same one, so a table can
switch source without renaming its columns. The logP is the exception: RDKit's
is Crippen's model, not XLogP3, so it is called `MolLogP`, and asking RDKit for
`XLogP` raises an error that says so.

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
`descriptors_table()` are the bulk forms. For a structure that is not in
PubChem at all, `provesid.pubchem_id.rdkit_descriptors("CCO")` computes the
same set from a SMILES.

## CompTox: synonyms and retired CAS numbers

CompTox stores a chemical's synonyms, former CAS numbers and registry codes
together in one `|`-separated column, which only a full scan can read.
`CompToxID` adds a table holding one row per distinct name of each chemical,
so `search_by_name(name, exact=True)` matches any of them, case-insensitively,
in about 40 µs:

```python
from provesid import CompToxID

with CompToxID() as db:
    db.search_by_name("Acetaldoxime", exact=True)[0]["PREFERRED_NAME"]  # 'Acetaldehyde oxime'
    db.get_by_casrn("39400-72-1")                                        # None: not a current number
    db.get_by_alternate_casrn("39400-72-1")["PREFERRED_NAME"]            # 'Atrazine'
```

`get_by_alternate_casrn` answers only when exactly one chemical lists the
number. Four of the 83 933 such numbers are listed by two unrelated chemicals,
and those four go unanswered. `Search("cas")` asks it whenever `get_by_casrn`
misses, which is how a retired number in an old dataset still resolves
offline. On a read-only database file with no name index, exact lookups fall
back to preferred names alone, with one warning.

## ChEMBL

`CheMBL` is used by `Search` for enrichment — it adds ChEMBL IDs to a
structure another source found — but it answers lookups by ChEMBL ID, name,
SMILES, InChI and InChIKey itself, and reports molecular properties, the
parent–salt hierarchy and pesticide classifications. A
compacted extract answers every `CheMBL` method exactly as the full release
does. See the [ChEMBL tutorial](../examples/chembl/chembl_tutorial.ipynb).
