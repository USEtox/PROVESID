# SQLite Clients

Four PROVESID clients read a local SQLite database: `PubChemID` (2.2 GB),
`CompToxID` (856 MB, 1.1 GiB once its name index is built), `ZeroPM` (460 MB) and `CheMBL` (2.4 GiB as a PROVESID
extract, 27.7 GiB as a full release). They inherit their connection handling
from one mixin, so all four close the same way, work in a `with` block, and
can be queried from several threads at once.

## What a caller sees

```python
from provesid import PubChemID

# A context manager
with PubChemID() as db:
    inchi = db.cas_to_inchi("50-78-2")
# the 2.2 GB file is no longer held open here

# Or close it by hand; idempotent
db = PubChemID()
db.close()
db.closed        # True
```

Querying a closed client raises
[`DatabaseClosedError`][provesid.sqlite_client.DatabaseClosedError], which
names the client and the file it was reading. It is a `RuntimeError` subclass,
so existing `except RuntimeError` handlers still catch it.

`Search` is a context manager too. It closes the source clients it
constructed, and leaves alone any you passed to its constructor — those belong
to you, and you may still be using them:

```python
from provesid import Search, PubChemID

with Search("cas") as s:
    df = s.search(["50-00-0", "64-17-5"])
# CompTox, PubChemID and ChEMBL are closed here

db = PubChemID()
with Search("cas", pubchem=db) as s:
    df = s.search("50-00-0")
db.get_by_cas("64-17-5")        # still open
```

## Threads

Each thread gets its own connection and cursor, created on that thread's first
query and closed with the client:

```python
from concurrent.futures import ThreadPoolExecutor

with PubChemID() as db:
    with ThreadPoolExecutor(8) as pool:
        rows = list(pool.map(db.get_by_cas, cas_numbers))
```

`close()` closes every one of them, from whichever thread calls it — so `with`
still means what it says when workers are involved. A query executing on
another thread while `close()` runs will fail; closing a client mid-query is a
caller error, and the alternative — blocking until the workers finish — makes
`with` unable to guarantee anything.

### Threads make this possible, not fast

A pool is worth it when each item does something *else* as well — a web
lookup, an RDKit canonicalisation, a file read. Measured on 400 CAS numbers,
each a `PubChemID` lookup followed by 20 ms of waiting: 8.52 s serially,
1.17 s on eight threads, **7.3×**.

A tight loop of nothing but local lookups is the opposite. These queries take
tens of microseconds against a warm page cache, which is far less than the
GIL handoff around each one costs, so 5 000 lookups that take 0.29 s serially
take 4.80 s on four threads and 14.65 s on eight. That is `sqlite3` under
CPython, not something this package can fix — plain `sqlite3.connect` per
thread measures the same — and the remedy is simply not to use a pool for
that shape of work. What changed here is that the pool no longer *fails*; how
fast it runs is still the caller's to measure.

## Why it exists

Each client opened its connection in `__init__` and closed it in `__del__`,
and nowhere else. Three ordinary consequences:

- **Nothing to call when you are finished.** A notebook cell re-running
  `db = PubChemID()` leaked the previous connection until the collector ran,
  and on Windows the file stayed locked while it did — which is why a
  re-download could not replace a database a live object was holding.
- **Nothing to put in a `with` block.** Every other resource in the package is
  scoped; these four were not.
- **One connection, shared by every thread.** `sqlite3` refuses by default to
  use a connection from a thread other than the one that created it, so the
  obvious way to resolve ten thousand CAS numbers against a local 2.2 GB
  database died on the first worker with
  `sqlite3.ProgrammingError: SQLite objects created in a thread can only be
  used in that same thread`.

## Why a connection per thread, and not `check_same_thread=False`

That flag removes the *check*, not the problem. Every thread would still share
one connection and, worse, one `self.cursor` — and these clients are written as
`cursor.execute(...)` followed by `cursor.fetchone()`, so two threads
interleaving those two statements would read each other's rows. Silently wrong
answers are a far worse failure than the exception they replace.

One connection per thread makes the existing query code correct as written.
SQLite serialises the writes `ZeroPM` performs when it builds an index or a
view. `check_same_thread=False` is still passed, but only so that `close()`,
called from whichever thread owns the object, may close the connections the
workers opened.

## Adding a client

A class joins in by inheriting the mixin and calling `_open_database` once the
file it wants exists:

```python
from provesid import SQLiteClient

class MyClient(SQLiteClient):
    def __init__(self, db_path):
        self.db_path = db_path
        self._open_database(db_path)          # row_factory=None for tuples

    def lookup(self, key):
        self.cursor.execute("SELECT v FROM t WHERE k = ?", (key,))
        row = self.cursor.fetchone()
        return None if row is None else row["v"]
```

The connection is opened eagerly, so an unreadable file still fails inside
`__init__`. `close()`, `with`, the per-thread handles and the `__del__`
backstop come with that one call.

`conn` has no setter: assigning to it would leave the replaced connection open
and unreachable, with nothing left to close it. A client that should run over
a connection it did not open — an in-memory database, or a stub in a test that
builds the instance with `object.__new__` — calls `_adopt_connection` instead,
which registers the connection so `close()` still reaches it.

::: provesid.sqlite_client

## CompTox's name index

The CompTox download indexes `PREFERRED_NAME` only. Its synonyms, former CAS
numbers and registry codes are stored together in one `|`-separated
`IDENTIFIER` column, and only a full scan can read that. `CompToxID` therefore
adds a table, `chemical_names`, holding one row per distinct name of each
chemical, so that `search_by_name(name, exact=True)` matches any of them,
case-insensitively, in about 40 µs:

```python
from provesid import CompToxID

with CompToxID() as db:
    db.search_by_name("Acetaldoxime", exact=True)[0]["PREFERRED_NAME"]  # 'Acetaldehyde oxime'
    db.search_by_name("39400-72-1", exact=True)[0]["PREFERRED_NAME"]    # 'Atrazine' (a retired CAS)
```

The table is built straight after a download. A database downloaded before
the index existed gets it on its first exact name lookup, which takes about
20 s once and adds 290 MiB; call `db.build_name_index()` to pay that at a time
of your choosing. On a read-only file, exact lookups fall back to preferred
names alone, with one warning.

The same table answers **retired and alternate CAS numbers**, the ones
CompTox lists among a chemical's identifiers but not as its `CASRN`:

```python
with CompToxID() as db:
    db.get_by_casrn("39400-72-1")                                # None: not a current number
    db.get_by_alternate_casrn("39400-72-1")["PREFERRED_NAME"]    # 'Atrazine'
```

`get_by_alternate_casrn` answers only when exactly one chemical lists the
number. Four of the 83,933 such numbers are listed by two unrelated
chemicals, and those four go unanswered. `Search("cas")` asks it whenever
`get_by_casrn` misses, so a retired number in an old dataset resolves
offline.
