"""
Closing a PROVESID database, and querying one from several threads.

Four clients read a local SQLite file: ``PubChemID`` (2.2 GB), ``CompToxID``
(856 MB), ``ZeroPM`` (460 MB) and ``CheMBL`` (2.4 GiB as an extract, 27.7 GiB
as a full release).  Each of them used to open its connection in ``__init__``
and close it only in ``__del__``, which left three ordinary problems:

* nothing to call when you are finished -- a notebook cell re-running
  ``db = PubChemID()`` leaked the previous connection, and on Windows the file
  stayed locked while it did, so the next download could not replace it;
* nothing to put in a ``with`` statement;
* one connection shared by every thread, which ``sqlite3`` refuses to allow --
  so the obvious way to resolve ten thousand CAS numbers against a local
  2.2 GB database, a ``ThreadPoolExecutor``, failed on the first worker with
  ``ProgrammingError: SQLite objects created in a thread can only be used in
  that same thread``.

All four now inherit :class:`provesid.SQLiteClient`, so all four have
``close()``, work as context managers, and hand each thread its own connection.
``Search`` closes the clients it constructed -- and only those.

Sections 1-4 run against a small database this script builds in a temporary
directory, so they need no installed dataset and no network.  Section 5 uses
the real clients if they happen to be installed, and says so if they are not.

Run with::

    python examples/sqlite_clients/connection_lifetime_demo.py
"""

import os
import sqlite3
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor

from provesid import DatabaseClosedError, SQLiteClient
from provesid.datasets import status


def rule(title):
    print(f"\n{'─' * 78}\n{title}\n")


class ElementTable(SQLiteClient):
    """A stand-in for the real clients: same mixin, 118 rows instead of 1.6 M.

    The only thing a client has to do is call :meth:`_open_database` once the
    file exists.  Everything below --- ``close``, ``with``, the per-thread
    connection --- comes with that one call.
    """

    def __init__(self, db_path):
        self.db_path = str(db_path)
        self._open_database(self.db_path)

    def symbol_of(self, number):
        """Return the symbol for an atomic number, through this thread's cursor."""
        self.cursor.execute(
            "SELECT symbol FROM elements WHERE number = ?", (number,)
        )
        row = self.cursor.fetchone()
        return None if row is None else row["symbol"]


def build_demo_database(path):
    """Write a one-table SQLite file with a handful of elements in it."""
    elements = [
        (1, "H"), (6, "C"), (7, "N"), (8, "O"), (9, "F"),
        (15, "P"), (16, "S"), (17, "Cl"), (35, "Br"), (53, "I"),
    ]
    conn = sqlite3.connect(path)
    conn.execute("CREATE TABLE elements (number INTEGER PRIMARY KEY, symbol TEXT)")
    conn.executemany("INSERT INTO elements VALUES (?, ?)", elements)
    conn.commit()
    conn.close()
    return [number for number, _ in elements]


def main():
    with tempfile.TemporaryDirectory() as tmp:
        db_path = os.path.join(tmp, "elements.db")
        numbers = build_demo_database(db_path)

        # ── 1. with ──────────────────────────────────────────────────────────
        rule("1. A client is a context manager, so the file is released at the "
             "end of\n   the block whatever happens inside it.")

        with ElementTable(db_path) as db:
            print(f"   element 6 is {db.symbol_of(6)}")
            print(f"   closed inside the block: {db.closed}")
        print(f"   closed after the block:  {db.closed}")

        # ── 2. close() ───────────────────────────────────────────────────────
        rule("2. Or close it by hand.  Idempotent, and a closed client says so "
             "rather\n   than failing somewhere deeper.")

        db = ElementTable(db_path)
        db.symbol_of(8)
        db.close()
        db.close()          # no complaint: closing twice is not an error
        try:
            db.symbol_of(8)
        except DatabaseClosedError as exc:
            print(f"   {exc}")

        print("\n   That is a RuntimeError subclass, so `except RuntimeError` "
              "still catches it,\n   and the message names the client and the "
              "file rather than leaving you with\n   sqlite3's 'Cannot operate "
              "on a closed database'.")

        # ── 3. Threads ───────────────────────────────────────────────────────
        rule("3. One connection per thread, opened on that thread's first "
             "query.  This is\n   the loop that used to raise "
             "sqlite3.ProgrammingError on the first worker.")

        with ElementTable(db_path) as db:
            queries = numbers * 200          # 2 000 lookups
            started = time.perf_counter()
            with ThreadPoolExecutor(8) as pool:
                symbols = list(pool.map(db.symbol_of, queries))
            elapsed = time.perf_counter() - started

            print(f"   {len(symbols)} lookups on 8 threads in {elapsed * 1000:.0f} ms")
            print(f"   distinct answers: {sorted(set(symbols))}")
            print(f"   connections open: {len(db._open_connections)} "
                  "(the main thread's, plus one per worker)")

        print("\n   Per-thread connections rather than check_same_thread=False: "
              "that flag only\n   removes the *check*, leaving every thread "
              "sharing one cursor.  These clients\n   are written as "
              "`cursor.execute(...)` then `cursor.fetchone()`, so two threads\n"
              "   interleaving those two statements would read each other's "
              "rows -- silently\n   wrong answers instead of a loud exception.")

        print("\n   Note what this does and does not buy.  The pool now *works*; "
              "for a tight\n   loop of nothing but local lookups it is slower "
              "than a plain loop, because\n   each query takes tens of "
              "microseconds and the GIL handoff around it costs\n   more.  A "
              "pool pays when each item also waits on something -- a request, "
              "a file,\n   an RDKit call.")

        # ── 4. Closing reaches every thread ──────────────────────────────────
        rule("4. close() closes the connections the worker threads opened, not "
             "just the\n   one on the thread that calls it.")

        db = ElementTable(db_path)
        barrier = threading.Barrier(4, timeout=10)

        def query(number):
            barrier.wait()                   # keep all four threads alive at once
            return db.symbol_of(number)

        with ThreadPoolExecutor(4) as pool:
            list(pool.map(query, [1, 6, 7, 8]))

        opened = list(db._open_connections)
        print(f"   connections before close: {len(opened)}")
        db.close()
        still_open = sum(1 for c in opened if _is_usable(c))
        print(f"   still open after close:   {still_open}")

    # ── 5. The real clients ──────────────────────────────────────────────────
    rule("5. The same three properties, on the databases you actually query.")

    installed = {
        row.dataset: row.present
        for row in status().itertuples()
    }
    print("   installed here: "
          + ", ".join(name for name, present in installed.items() if present)
          + "\n")

    demos = [
        ("pubchem", "PubChemID", _demo_pubchem),
        ("comptox", "CompToxID", _demo_comptox),
        ("chembl", "CheMBL", _demo_chembl),
    ]
    for key, label, demo in demos:
        if installed.get(key):
            demo()
        else:
            print(f"   {label:<10} not installed -- "
                  f"provesid.datasets.fetch({key!r}) would install it")

    print("""
   And Search closes what it opened:

       with Search("cas") as s:
           df = s.search(["50-00-0", "64-17-5"])
       # CompTox, PubChemID and ChEMBL are closed here

   A client you passed in yourself is left alone -- it is yours, and you may
   still be using it:

       db = PubChemID()
       with Search("cas", pubchem=db) as s:
           df = s.search("50-00-0")
       db.get_by_cas("64-17-5")        # still open
""")


def _is_usable(connection):
    """Whether a sqlite3 connection still answers a query."""
    try:
        connection.execute("SELECT 1")
        return True
    except sqlite3.ProgrammingError:
        return False


def _demo_pubchem():
    """Query PubChemID from a pool of threads and close it afterwards."""
    from provesid import PubChemID

    cas_numbers = ["50-00-0", "64-17-5", "71-43-2", "108-88-3", "50-78-2"] * 4
    with PubChemID(auto_download=False) as db:
        with ThreadPoolExecutor(8) as pool:
            rows = list(pool.map(db.get_by_cas, cas_numbers))
        found = [row for row in rows if row]
        print(f"   PubChemID  {len(found)}/{len(rows)} CAS numbers resolved on "
              f"8 threads; {len(db._open_connections)} connections")
    print(f"              closed: {db.closed}")


def _demo_comptox():
    """Query CompToxID inside a with block."""
    from provesid import CompToxID

    with CompToxID(auto_download=False) as db:
        print(f"   CompToxID  50-78-2 -> {db.casrn_to_dtxsid('50-78-2')}")
    print(f"              closed: {db.closed}")


def _demo_chembl():
    """Query CheMBL inside a with block."""
    from provesid import CheMBL

    with CheMBL(auto_download=False) as db:
        compound = db.search_by_chembl_id("CHEMBL25")
        print(f"   CheMBL     CHEMBL25 -> {compound['pref_name']}")
    print(f"              closed: {db.closed}")


if __name__ == "__main__":
    main()
