"""
Look compounds up in ChEMBL by name, exactly or loosely, and time it.

``CheMBL.search_by_name`` answers from two places at once: the preferred name in
``molecule_dictionary`` and every trade name and alias in ``molecule_synonyms``.
It used to ask for both in a single ``OR`` across a ``LEFT JOIN``, which no
index can serve -- every lookup scanned all 2.9 million compounds, at about
740 ms a call.  It now asks the two questions separately and merges the answers,
which lets each one use an index.

On an extract built by ``CheMBL.compact()`` -- which carries the two
``lower(...)`` indexes this needs -- an exact lookup costs roughly **10 us**
instead of 740 ms.  A full ChEMBL release has no such indexes and still scans,
but it scans two small queries rather than a join, and lands around 220 ms.

Results are also ordered by ``molregno`` now.  That sounds cosmetic and is not:
with ``limit`` smaller than the number of matches, an unordered query changes
*which* compounds it returns from run to run.

Run with::

    python examples/chembl/name_search_demo.py

Needs a ChEMBL database on disk; see ``examples/chembl/compact_demo.py`` for how
to shrink one to 2.6 GB.
"""

import time

from provesid import CheMBL


def show(title, compounds):
    """Print a search result as one line per compound."""
    print(f"\n{title}")
    if not compounds:
        print("    (no matches)")
        return
    for compound in compounds:
        name = compound.get("pref_name") or "(no preferred name)"
        print(f"    {compound['chembl_id']:<12} {name}")


def main():
    chembl = CheMBL()
    kind = "extract" if chembl.is_compact else "full release"
    print(f"ChEMBL {chembl.release} ({kind}) at {chembl.db_path}")

    # ── 1. Exact and substring are different questions ───────────────────────
    # exact=True means the name *is* the compound's name or one of its
    # synonyms.  The default is a substring match, which is permissive on
    # purpose and will surface compounds you did not mean.
    show("aspirin, exactly:", chembl.search_by_name("aspirin", exact=True))
    show("aspirin, as a substring:", chembl.search_by_name("aspirin", limit=5))

    # A misspelling that still matches something, because 'Evasprin' is a
    # synonym of PHENYRAMIDOL and contains 'asprin'.
    show("asprin (misspelt), as a substring:", chembl.search_by_name("asprin"))
    show("asprin (misspelt), exactly:", chembl.search_by_name("asprin", exact=True))

    # ── 2. Synonyms count as names ───────────────────────────────────────────
    # A trade name finds the compound just as a preferred name does.
    show("Tylenol (a trade name):", chembl.search_by_name("Tylenol", exact=True))

    # ── 3. Truncated results are reproducible ────────────────────────────────
    # With a limit below the match count, the lowest molregnos come back -- the
    # same ones on every call, on every copy of the database.
    first = [c["chembl_id"] for c in chembl.search_by_name("acid", limit=5)]
    again = [c["chembl_id"] for c in chembl.search_by_name("acid", limit=5)]
    print(f"\n'acid', limit 5, twice:\n    {first}\n    {again}")
    print(f"    identical: {first == again}")

    # ── 4. What the rewrite costs now ────────────────────────────────────────
    # Time the SQL alone, without the per-result get_compound calls, so the
    # number is about the lookup rather than about how many compounds matched.
    query = """
    SELECT molregno FROM molecule_dictionary WHERE LOWER(pref_name) = LOWER(?)
    UNION
    SELECT molregno FROM molecule_synonyms WHERE LOWER(synonyms) = LOWER(?)
    ORDER BY molregno LIMIT 100
    """
    names = [row[0] for row in chembl.conn.execute(
        "SELECT pref_name FROM molecule_dictionary "
        "WHERE pref_name IS NOT NULL LIMIT 200"
    )]
    repeats = 200 if chembl.is_compact else 5
    start = time.perf_counter()
    for i in range(repeats):
        name = names[i % len(names)]
        chembl.conn.execute(query, (name, name)).fetchall()
    per_call = (time.perf_counter() - start) / repeats
    print(f"\nexact name lookup on this database: {per_call * 1e6:.1f} us per call")
    if not chembl.is_compact:
        print("    (a compact() extract carries the indexes that make this ~10 us)")

    # ── 5. The plan, if you want to see why ──────────────────────────────────
    print("\nquery plan:")
    for row in chembl.conn.execute("EXPLAIN QUERY PLAN " + query, ("x", "x")):
        print(f"    {row[-1]}")


if __name__ == "__main__":
    main()
