# SQLite client examples

Four PROVESID clients read a local SQLite database: `PubChemID`, `CompToxID`,
`ZeroPM` and `CheMBL`. They all inherit their connection handling from
`provesid.SQLiteClient`, so they all behave the same way about closing,
`with`, and threads.

## Files

- `connection_lifetime_demo.py` — `close()`, `with`, per-thread connections
  and what a closed client does, first against a ten-row database the script
  builds in a temporary directory and then against whichever real datasets are
  installed. Downloads nothing and needs no network.

## The three things it gives you

```python
from provesid import PubChemID

# 1. A context manager
with PubChemID() as db:
    inchi = db.cas_to_inchi("50-78-2")
# the 2.2 GB file is no longer held open here

# 2. close(), for when a with block does not fit
db = PubChemID()
db.close()
db.closed          # True; querying it now raises DatabaseClosedError

# 3. A connection per thread
from concurrent.futures import ThreadPoolExecutor

with PubChemID() as db:
    with ThreadPoolExecutor(8) as pool:
        rows = list(pool.map(db.get_by_cas, cas_numbers))
```

`Search` closes the clients it constructed, and leaves alone any you passed to
it yourself:

```python
from provesid import Search

with Search("cas") as s:
    df = s.search(["50-00-0", "64-17-5"])
# CompTox, PubChemID and ChEMBL are closed here
```

## Why it exists

Each client opened its connection in `__init__` and closed it in `__del__`,
and nowhere else. So there was nothing to call when you were finished — a
notebook cell re-running `db = PubChemID()` leaked the previous connection, and
on Windows the file stayed locked while it did, which is why a re-download
could not replace a database a live object was holding. There was nothing to
put in a `with` block either.

The third consequence was the one users met first. A single connection is
bound to the thread that created it, so the obvious way to resolve ten
thousand CAS numbers against a local 2.2 GB database — a `ThreadPoolExecutor`
— died on the first worker with:

```
sqlite3.ProgrammingError: SQLite objects created in a thread can only be
used in that same thread.
```

## Why a connection per thread, and not `check_same_thread=False`

That flag removes the check, not the problem. Every thread would still share
one connection and, worse, one `self.cursor` — and these clients are written as
`cursor.execute(...)` followed by `cursor.fetchone()`, so two threads
interleaving those two statements would read each other's rows. Silently wrong
answers are a much worse failure than the exception they replace.

One connection per thread makes the existing query code correct as written.
SQLite serialises the writes `ZeroPM` performs when it builds an index or a
view, and `close()` closes every connection the client opened, on every thread,
so `with` still means what it says.

## Threads make this possible, not fast

A pool pays when each item does something else as well — a web lookup, an
RDKit canonicalisation, a file read. Measured on 400 CAS numbers, each a
`PubChemID` lookup plus 20 ms of waiting: 8.52 s serially against 1.17 s on
eight threads, 7.3×.

A tight loop of nothing but local lookups is the opposite: those queries take
tens of microseconds against a warm cache, less than the GIL handoff around
each one, so 5 000 of them take 0.29 s serially and 14.65 s on eight threads.
That is `sqlite3` under CPython — plain `sqlite3.connect` per thread measures
the same — and the answer is to not use a pool for that shape of work. What
changed here is that the pool no longer *fails*.
