"""Validate name-based resolution against CASRN ground truth (CompTox).

Strategy (suggested by the maintainer):

  CompTox carries both CASRN and PREFERRED_NAME for every chemical.  A CASRN
  lookup is (almost) always correct, so it serves as ground truth.  We:

    1. Sample N CompTox rows that have CASRN, PREFERRED_NAME and INCHIKEY.
    2. Treat each row's (DTXSID, INCHIKEY) as the truth for that CASRN.
    3. Resolve the *same* chemical by its PREFERRED_NAME through a CompTox-only
       ``Search`` and check whether we recover the same compound.
    4. Report accuracy and dump the mismatches — those are the real
       "wrong name hit" cases, ready to become regression tests.

The Search is restricted to the CompTox source so the experiment isolates the
name-matching / ranking logic against the very table we are trying to
reproduce ("can you reproduce the CASRN table from the names?").

Usage::

    uv run python scripts/validate_name_search.py --n 300
    uv run python scripts/validate_name_search.py --n 300 --tuned
    uv run python scripts/validate_name_search.py --n 300 --out mismatches.csv
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd

from provesid.comptox import CompToxID
from provesid.search import Search


def _skeleton(ik):
    return ik[:14] if isinstance(ik, str) and len(ik) >= 14 else ik


def sample_rows(comptox: CompToxID, n: int) -> pd.DataFrame:
    """Deterministically sample n rows spread across the CompTox table.

    Uses a ``rowid``-stride so the sample is reproducible without RNG and is
    spread across the whole table rather than clustered alphabetically.
    """
    total = comptox.conn.execute("SELECT COUNT(*) FROM chemicals").fetchone()[0]
    step = max(1, total // (n * 4))  # over-sample; we filter nulls/dupes below
    cur = comptox.conn.execute(
        f"""
        SELECT DTXSID, CASRN, PREFERRED_NAME, IUPAC_NAME, INCHIKEY, IDENTIFIER
        FROM chemicals
        WHERE CASRN IS NOT NULL AND CASRN != ''
          AND PREFERRED_NAME IS NOT NULL AND PREFERRED_NAME != ''
          AND INCHIKEY IS NOT NULL AND INCHIKEY != ''
          AND (rowid % {step}) = 1
        LIMIT {n * 4}
        """
    )
    rows = [dict(r) for r in cur.fetchall()]
    return pd.DataFrame(rows)


def pick_query_name(row, name_source: str):
    """Choose the query name for a row given the requested name source.

    Returns ``None`` when no suitable name exists (row is then skipped).
    """
    if name_source == "preferred":
        return row["PREFERRED_NAME"]
    if name_source == "iupac":
        return row.get("IUPAC_NAME") or None
    if name_source == "synonym":
        ident = row.get("IDENTIFIER") or ""
        parts = [p.strip() for p in ident.split("|") if p.strip()]
        pref = (row["PREFERRED_NAME"] or "").strip().lower()
        cas = (row["CASRN"] or "").strip()
        # Skip the CAS token and the preferred name; take a different synonym.
        for p in parts:
            if p == cas or p.lower() == pref:
                continue
            if any(ch.isalpha() for ch in p):  # looks like a name, not a code
                return p
        return None
    raise ValueError(name_source)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=300, help="number of chemicals to test")
    ap.add_argument("--tuned", action="store_true",
                    help="enable fuzzy + OPSIN + larger top_k tuning")
    ap.add_argument("--all-sources", action="store_true",
                    help="resolve through all 5 sources (real scenario where "
                         "cross-source disagreement produces wrong hits)")
    ap.add_argument("--name-source", choices=["preferred", "iupac", "synonym"],
                    default="preferred",
                    help="which name to query by (default: preferred)")
    ap.add_argument("--out", default=None, help="write mismatches to this CSV")
    args = ap.parse_args(argv)

    comptox = CompToxID()
    df = sample_rows(comptox, args.n)

    # Keep only CASRNs that map to exactly one chemical (the documented caveat).
    cas_counts = df.groupby("CASRN")["DTXSID"].nunique()
    unique_cas = set(cas_counts[cas_counts == 1].index)
    df = df[df["CASRN"].isin(unique_cas)].drop_duplicates("CASRN")

    # Resolve the query name per the chosen source; drop rows without one.
    df["_query"] = df.apply(lambda r: pick_query_name(r, args.name_source), axis=1)
    df = df[df["_query"].notna()].drop_duplicates("_query").head(args.n)
    df = df.reset_index(drop=True)
    print(f"Testing {len(df)} chemicals | name-source={args.name_source} | "
          f"sources={'all' if args.all_sources else 'comptox-only'} | "
          f"{'tuned' if args.tuned else 'default'} config\n")

    kwargs = dict(identifier_type="name", show_progress=True)
    if args.all_sources:
        # Pass no clients so Search lazily builds all 5 from local data.
        pass
    else:
        # CompTox-only: pass comptox, force the rest off.
        kwargs.update(comptox=comptox, chebi=None, pubchem=None,
                      zeropm=None, chembl=None)
    if args.tuned:
        kwargs.update(fuzzy=True, use_opsin=True, top_k_per_source=10,
                      fuzzy_score_cutoff=70.0)
    s = Search(**kwargs)

    res = s.search(df["_query"].tolist())

    truth_dtxsid = df["DTXSID"].tolist()
    truth_ik = df["INCHIKEY"].tolist()

    # Correctness is structure identity (skeleton), which is source-independent.
    # Three outcomes per query:
    #   correct   – returned a structure whose skeleton matches the truth
    #   wrong_hit – returned a *different* structure (the precision failure)
    #   not_found – returned no structure at all (a recall gap, not a wrong hit)
    correct = wrong = not_found = 0
    wrong_hits = []
    for i, row in res.iterrows():
        got_ik = row.get("InChIKey")
        truth_skel = _skeleton(truth_ik[i])
        if got_ik is None or (isinstance(got_ik, float)):
            not_found += 1
            continue
        if _skeleton(got_ik) == truth_skel:
            correct += 1
        else:
            wrong += 1
            wrong_hits.append({
                "query_name": df["_query"][i],
                "preferred_name": df["PREFERRED_NAME"][i],
                "truth_CASRN": df["CASRN"][i],
                "truth_DTXSID": truth_dtxsid[i],
                "truth_InChIKey": truth_ik[i],
                "got_name": row.get("name"),
                "got_InChIKey": got_ik,
                "got_DTXSID": row.get("DTXSID"),
                "confidence": row.get("confidence"),
                "match_method": row.get("match_method"),
            })

    n = len(df)
    print(f"\n  correct  (structure reproduced): {correct}/{n}  ({correct/n:.1%})")
    print(f"  WRONG HIT (different structure):  {wrong}/{n}  ({wrong/n:.1%})")
    print(f"  not found (no structure):         {not_found}/{n}  ({not_found/n:.1%})")
    found = correct + wrong
    if found:
        print(f"  precision when a hit is returned: {correct}/{found}  ({correct/found:.1%})")

    if wrong_hits:
        wdf = pd.DataFrame(wrong_hits)
        print(f"\nWRONG-HIT cases (returned the wrong compound) — regression candidates:\n")
        with pd.option_context("display.max_colwidth", 34, "display.width", 220):
            print(wdf[["query_name", "preferred_name", "got_name",
                       "match_method", "confidence"]].head(25).to_string(index=False))
        if args.out:
            wdf.to_csv(args.out, index=False)
            print(f"\nWrote {len(wdf)} wrong-hit cases to {args.out}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
