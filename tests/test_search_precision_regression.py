"""Precision regression test for name-based resolution.

Uses CompTox as ground truth: every chemical carries a CASRN and several
names/synonyms.  We resolve a sample of chemicals *by a synonym* through a
CompTox-only :class:`~provesid.Search` (precision-first default config) and
assert that the resolver never returns a **structurally different** compound
(a "wrong hit").

Correctness is measured by InChIKey *skeleton* (the connectivity layer), which
is source-independent and tolerant of stereo/charge/tautomer naming variants.
Three outcomes are possible per query:

- **correct**   — returned a structure whose skeleton matches the truth,
- **wrong_hit** — returned a *different* structure (the precision failure we
  guard against), or
- **not_found** — returned no structure (a recall gap, explicitly allowed).

The guard asserts ``wrong_hit == 0``: the resolver may fail to find a compound,
but it must not confidently return the wrong one.  This locks in the behaviour
established by the candidate-pool / clustering / query-aware-ranking rework.

Marked ``integration`` + ``slow``; skipped automatically when the offline
CompTox database is unavailable (e.g. a fresh checkout with no downloaded data).

Run with::

    uv run pytest tests/test_search_precision_regression.py -v -m integration
"""

from __future__ import annotations

import pytest

from provesid.search import Search

# Sample size — kept modest so the test runs in a couple of seconds with the
# CompTox-only resolver while still exercising a representative spread.
SAMPLE_SIZE = 80


def _skeleton(inchikey):
    """Return the 14-char InChIKey connectivity skeleton, or the input as-is."""
    if isinstance(inchikey, str) and len(inchikey) >= 14:
        return inchikey[:14]
    return inchikey


def _pick_synonym(row):
    """Choose a query synonym for a CompTox row that is not its CAS or name.

    Args:
        row: A ``sqlite3.Row`` (mapping) with ``CASRN``, ``PREFERRED_NAME`` and
            the pipe-delimited ``IDENTIFIER`` field.

    Returns:
        A synonym string distinct from the CASRN and preferred name, or
        ``None`` when no suitable name-like token exists.
    """
    ident = row["IDENTIFIER"] or ""
    pref = (row["PREFERRED_NAME"] or "").strip().lower()
    cas = (row["CASRN"] or "").strip()
    for part in (p.strip() for p in ident.split("|")):
        if not part or part == cas or part.lower() == pref:
            continue
        if any(ch.isalpha() for ch in part):  # name-like, not a bare code
            return part
    return None


@pytest.fixture(scope="module")
def comptox():
    """Yield an offline CompToxID client, skipping if its database is absent."""
    from provesid.comptox import CompToxID

    try:
        client = CompToxID()
        # Touch the connection so a missing/empty DB skips rather than errors.
        client.conn.execute("SELECT 1 FROM chemicals LIMIT 1").fetchone()
    except Exception as exc:  # pragma: no cover - environment dependent
        pytest.skip(f"CompTox offline database unavailable: {exc}")
    return client


@pytest.fixture(scope="module")
def synonym_sample(comptox):
    """Build a deterministic sample of (query synonym, truth skeleton) pairs.

    Samples rows spread across the CompTox table via a ``rowid`` stride (no RNG,
    so the sample is reproducible), keeps only CASRNs that map to a single
    structure, and selects a synonym distinct from the preferred name.

    Returns:
        List of ``(synonym, truth_inchikey)`` tuples of length up to
        :data:`SAMPLE_SIZE`.
    """
    total = comptox.conn.execute("SELECT COUNT(*) FROM chemicals").fetchone()[0]
    step = max(1, total // (SAMPLE_SIZE * 6))
    cur = comptox.conn.execute(
        f"""
        SELECT CASRN, PREFERRED_NAME, INCHIKEY, IDENTIFIER
        FROM chemicals
        WHERE CASRN IS NOT NULL AND CASRN != ''
          AND PREFERRED_NAME IS NOT NULL AND PREFERRED_NAME != ''
          AND INCHIKEY IS NOT NULL AND INCHIKEY != ''
          AND IDENTIFIER IS NOT NULL AND IDENTIFIER != ''
          AND (rowid % {step}) = 1
        LIMIT {SAMPLE_SIZE * 6}
        """
    )
    rows = [dict(r) for r in cur.fetchall()]

    # Drop CASRNs that map to more than one structure (the documented caveat).
    cas_structs = {}
    for r in rows:
        cas_structs.setdefault(r["CASRN"], set()).add(r["INCHIKEY"])
    ambiguous = {c for c, ss in cas_structs.items() if len(ss) > 1}

    sample = []
    seen = set()
    for r in rows:
        if r["CASRN"] in ambiguous:
            continue
        syn = _pick_synonym(r)
        if syn is None or syn.lower() in seen:
            continue
        seen.add(syn.lower())
        sample.append((syn, r["INCHIKEY"]))
        if len(sample) >= SAMPLE_SIZE:
            break

    if len(sample) < 10:  # pragma: no cover - environment dependent
        pytest.skip("Not enough CompTox synonym samples to run the regression.")
    return sample


@pytest.mark.integration
@pytest.mark.slow
def test_name_resolution_has_zero_wrong_hits(comptox, synonym_sample):
    """The resolver must never return a structurally different compound.

    Resolves each sampled synonym through a CompTox-only, precision-first
    ``Search`` and asserts that no query yields a structure whose InChIKey
    skeleton differs from the CASRN-derived ground truth.  Not-found results
    are permitted (recall is not asserted here).
    """
    s = Search(
        "name",
        show_progress=False,
        comptox=comptox,
        chebi=None,
        pubchem=None,
        zeropm=None,
        chembl=None,
    )

    queries = [syn for syn, _ in synonym_sample]
    res = s.search(queries)

    wrong_hits = []
    correct = not_found = 0
    for (syn, truth_ik), (_, row) in zip(synonym_sample, res.iterrows()):
        got_ik = row.get("InChIKey")
        if not isinstance(got_ik, str):  # None / NaN
            not_found += 1
            continue
        if _skeleton(got_ik) == _skeleton(truth_ik):
            correct += 1
        else:
            wrong_hits.append((syn, row.get("name"), truth_ik, got_ik))

    # Recall (correct vs not_found) is informational; precision is the contract.
    assert not wrong_hits, (
        f"{len(wrong_hits)} wrong hit(s) out of {len(synonym_sample)} "
        f"(correct={correct}, not_found={not_found}). Examples: {wrong_hits[:5]}"
    )
