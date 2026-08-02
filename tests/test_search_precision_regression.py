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


# ─────────────────────────────────────────────────────────────────────────────
# Misspelled-name regression
#
# A typo used to resolve to an unrelated compound, labelled as an *exact* match:
#
#   Search("name", fuzzy=True).search("asprin")
#     -> PHENYRAMIDOL, match_method="exact_name", confidence=0.74
#
# Three defects combined to produce that:
#   1. CheMBL.search_by_name had no exact mode, so Search's exact pass got a
#      substring match -- PHENYRAMIDOL carries the synonym "Evasprin", which
#      contains "asprin" -- and tagged it "exact_name".
#   2. The "is the exact pass strong?" test used a fuzzy score, and
#      WRatio("asprin", "Evasprin") == 85.7 cleared the 80 cut-off, so the fuzzy
#      widening that would have found aspirin never ran.
#   3. A fuzzy match's confidence base was the raw similarity (up to 1.0) while
#      an exact name match was pinned at 0.80, so a typo could outrank the
#      correctly spelled name.
# ─────────────────────────────────────────────────────────────────────────────

ASPIRIN_SKELETON = "BSYNRYMUTXBXSQ"
CAFFEINE_SKELETON = "RYYVLZVUVIJVGH"

# (query, expected skeleton) -- common misspellings of well-known compounds.
MISSPELLINGS = [
    ("asprin", ASPIRIN_SKELETON),
    ("caffiene", CAFFEINE_SKELETON),
]


@pytest.fixture(scope="module")
def all_sources_available():
    """Skip unless every offline source Search needs is present."""
    s = Search("name", fuzzy=True, show_progress=False)
    s._ensure_clients()
    missing = [
        key
        for key, client in [
            ("chebi", s._chebi),
            ("comptox", s._comptox),
            ("pubchem", s._pubchem),
            ("zeropm", s._zeropm),
            ("chembl", s._chembl),
        ]
        if client is None
    ]
    if missing:  # pragma: no cover - environment dependent
        pytest.skip(f"Offline sources unavailable: {', '.join(missing)}")
    return True


@pytest.mark.integration
@pytest.mark.slow
@pytest.mark.parametrize("query,expected_skeleton", MISSPELLINGS)
def test_misspelled_name_resolves_to_the_right_compound(
    all_sources_available, query, expected_skeleton
):
    """A misspelled name must resolve to the intended compound, or to nothing."""
    row = Search("name", fuzzy=True, show_progress=False).search(query).iloc[0]

    got = row["InChIKey"]
    assert isinstance(got, str), f"{query!r} resolved to nothing"
    assert _skeleton(got) == expected_skeleton, (
        f"{query!r} resolved to {row['name']!r} ({got}), "
        f"expected skeleton {expected_skeleton}"
    )


@pytest.mark.integration
@pytest.mark.slow
@pytest.mark.parametrize("query,_expected", MISSPELLINGS)
def test_misspelling_is_not_labelled_an_exact_match(
    all_sources_available, query, _expected
):
    """A typo must never be reported as ``exact_name``."""
    row = Search("name", fuzzy=True, show_progress=False).search(query).iloc[0]
    assert row["match_method"] != "exact_name", (
        f"{query!r} was labelled an exact name match "
        f"(resolved to {row['name']!r})"
    )


@pytest.mark.integration
@pytest.mark.slow
def test_misspelling_scores_below_correct_spelling(all_sources_available):
    """The correctly spelled name must be the more confident of the two."""
    s = Search("name", fuzzy=True, show_progress=False)
    res = s.search(["asprin", "aspirin"]).set_index("query")

    typo, correct = res.loc["asprin"], res.loc["aspirin"]
    assert _skeleton(typo["InChIKey"]) == _skeleton(correct["InChIKey"])
    assert typo["confidence"] < correct["confidence"], (
        f"typo confidence {typo['confidence']} >= "
        f"correct-spelling confidence {correct['confidence']}"
    )


@pytest.mark.integration
@pytest.mark.slow
def test_nonsense_name_resolves_to_nothing(all_sources_available):
    """A string that is not a chemical name must not match anything.

    Guards the fuzzy *retrieval* scorer: rapidfuzz ``WRatio`` scores a short
    candidate name highly whenever it appears inside the query, so
    ``"zzzznotachemical"`` matched "Mica" and ``"caffiene"`` matched a compound
    named "ne". ``ratio`` puts both at 40.
    """
    row = Search("name", fuzzy=True, show_progress=False).search(
        "zzzznotachemical"
    ).iloc[0]
    assert not isinstance(row["InChIKey"], str), (
        f"nonsense query matched {row['name']!r} ({row['InChIKey']})"
    )


@pytest.mark.integration
@pytest.mark.slow
def test_typo_needs_fuzzy_enabled(all_sources_available):
    """Without ``fuzzy=True`` a typo must return nothing, not a guess."""
    row = Search("name", show_progress=False).search("asprin").iloc[0]
    assert not isinstance(row["InChIKey"], str), (
        f"non-fuzzy search matched {row['name']!r} for a misspelling"
    )
