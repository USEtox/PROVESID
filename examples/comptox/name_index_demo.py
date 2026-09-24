"""
Look chemicals up in CompTox by any of their names, and time it.

The CompTox database indexes each chemical's preferred name and nothing else.
Its synonyms, retired CAS numbers and registry codes are stored in a single
``|``-separated ``IDENTIFIER`` column. An exact lookup therefore used to miss
"Acetaldoxime" (a synonym of acetaldehyde oxime), and the only way to reach it
was a substring scan at about 3 s a query.

``CompToxID`` now keeps a name index, the ``chemical_names`` table, with one
row per distinct name of each chemical. ``search_by_name(..., exact=True)``
matches any of them, case-insensitively, in tens of microseconds. A chemical
*called* the query is returned before one that only lists it as a synonym.

The index is built after a download. A database downloaded earlier gets it
on its first exact name lookup: about 20 s, once, adding 290 MiB to the file.
It also answers retired CAS numbers: ``get_by_alternate_casrn``.

Run with::

    uv run python examples/comptox/name_index_demo.py
"""

import logging
import time

from provesid import CompToxID

# The one-time build announces itself at WARNING.
logging.basicConfig(level=logging.WARNING)

queries = [
    "Aspirin",                     # the preferred name
    "ACETYLSALICYLIC ACID",        # a synonym, in any case
    "Acetaldoxime",                # a synonym the old lookup missed
    "2-acetyloxybenzoic acid",     # an IUPAC name
    "39400-72-1",                  # a retired CAS number of atrazine
]

with CompToxID() as db:
    if not db.has_name_index:
        started = time.perf_counter()
        rows = db.build_name_index()
        print(f"Built the name index: {rows:,} names in {time.perf_counter() - started:.0f} s\n")

    for query in queries:
        started = time.perf_counter()
        hits = db.search_by_name(query, exact=True, limit=3)
        elapsed_us = (time.perf_counter() - started) * 1e6
        found = ", ".join(f"{h['PREFERRED_NAME']} ({h['DTXSID']})" for h in hits) or "-"
        print(f"{query:26s} {elapsed_us:7.0f} µs  {found}")

    # get_by_name still means the preferred name, and only that.
    print("\nget_by_name('Acetaldoxime'):", db.get_by_name("Acetaldoxime"))

    # A retired CAS number is no chemical's CASRN, but CompTox still lists it.
    # Search("cas") falls back to this whenever get_by_casrn misses.
    for cas in ["39400-72-1", "11126-35-5"]:
        hit = db.get_by_alternate_casrn(cas)
        print(f"{cas}: get_by_casrn={db.get_by_casrn(cas)}, "
              f"alternate -> {hit['PREFERRED_NAME']} (current CAS {hit['CASRN']})")
