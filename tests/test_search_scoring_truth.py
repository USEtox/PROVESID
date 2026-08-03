"""Truth-based regression tests for the :class:`~provesid.Search` scoring system.

Two independent defects made ``Search("cas")`` return the wrong structure for
18 of 65 pesticide CAS numbers (reported 2026-08-03):

1. **ChEBI record offsets drifted** (``provesid.chebi.ChebiSDF``).  The SDF index
   recorded *text-mode* byte offsets, so every CRLF line in the ChEBI SDF
   under-counted by one byte.  The drift grew to ~59 kB by the end of the file
   and ``get_compound_by_id`` seeked into a **neighbouring record**, so ChEBI
   answered "Mefluidide" for metaldehyde's CAS.  The ChEBI *data* was correct all
   along; only the offsets were wrong.
2. **Corroboration did not count in the ranking** (``provesid.search``).  A hit's
   confidence was driven by which database answered, not by how many agreed: a
   lone ChEBI hit scored a flat 0.90 while a structure CompTox, PubChem and
   ZeroPM all agreed on scored 0.8777, so the uncorroborated hit won.

The tests below lock in both fixes.  Ground truth is the bundled **CompTox**
database, which carries a curated CASRN → structure mapping for ~1.2 M
substances; correctness is measured by InChIKey *skeleton* (the connectivity
layer) so stereo/charge/tautomer representation differences do not count as
disagreements.

Unit tests for the scoring rules run without any database.  The truth-sample
tests are marked ``integration`` + ``slow`` and skip when the offline databases
are not present.

Run with::

    uv run pytest tests/test_search_scoring_truth.py -v
    uv run pytest tests/test_search_scoring_truth.py -v -m integration
"""

from __future__ import annotations

import pytest

from provesid.search import (
    _BASE_CONFIDENCE,
    _SUPPORT_FACTOR,
    Search,
    _has_attachment_point,
    resolve_cascade,
)

# Sample sizes — large enough to be representative, small enough to keep the
# offline run to a few seconds per test.
CAS_SAMPLE_SIZE = 250
NAME_SAMPLE_SIZE = 120
CHEBI_SAMPLE_SIZE = 400

# CompTox and PubChem genuinely disagree about the structure behind a handful of
# CASRNs (mostly organometallics and polymeric salts registered under one number
# with different stoichiometry).  Those are data disagreements, not ranking
# defects, so a small rate is tolerated; the *scoring* contract is asserted
# separately by :func:`test_corroborated_cas_hits_are_never_wrong`.
MAX_WRONG_RATE = 0.01


def _skeleton(inchikey):
    """Return the 14-char InChIKey connectivity skeleton, or the input as-is."""
    if isinstance(inchikey, str) and len(inchikey) >= 14:
        return inchikey[:14]
    return inchikey


# ─────────────────────────────────────────────────────────────────────────────
# Scoring rules — pure unit tests, no database needed
# ─────────────────────────────────────────────────────────────────────────────


@pytest.fixture
def resolver():
    """A Search instance used only for its scoring methods (no client init)."""
    return Search("cas", show_progress=False)


def test_corroboration_beats_provenance(resolver):
    """Three agreeing databases must outrank one uncorroborated database.

    This is the bug itself, in numbers: a lone ChEBI CAS hit used to score
    ``0.90`` while a CompTox/PubChem/ZeroPM consensus scored ``0.8777``, so the
    resolver returned the compound *nothing* corroborated.
    """
    solo = resolver._compute_confidence("exact_cas", 1.0, n_source_support=1)
    trio = resolver._compute_confidence("exact_cas", 0.9506, n_source_support=3)

    assert trio > solo, (
        f"a three-source consensus ({trio}) must outrank a single-source hit "
        f"({solo}) for the same match method"
    )


def test_confidence_rises_with_source_support(resolver):
    """Confidence must be non-decreasing in the number of corroborating sources."""
    scores = [
        resolver._compute_confidence("exact_cas", 0.95, n_source_support=n)
        for n in (1, 2, 3, 4, 5)
    ]
    assert scores == sorted(scores), scores
    assert scores[0] < scores[-1]


def test_single_source_hit_is_penalised(resolver):
    """An uncorroborated hit must not report the full method base confidence."""
    solo = resolver._compute_confidence("exact_cas", 1.0, n_source_support=1)
    assert solo < _BASE_CONFIDENCE["exact_cas"]
    assert solo == pytest.approx(_BASE_CONFIDENCE["exact_cas"] * _SUPPORT_FACTOR[1])


def test_opsin_only_cluster_is_not_penalised_for_support(resolver):
    """OPSIN parses a name into a structure; that is not a corroboration question."""
    opsin = resolver._compute_confidence("opsin", 1.0, n_source_support=0)
    assert opsin == pytest.approx(_BASE_CONFIDENCE["opsin"])


def test_zero_consensus_still_short_circuits(resolver):
    """No candidates at all must stay at zero confidence, not a support-scaled floor."""
    assert resolver._compute_confidence("exact_cas", 0.0, n_source_support=0) == 0.0


@pytest.mark.parametrize("method", sorted(_BASE_CONFIDENCE))
@pytest.mark.parametrize("n", [0, 1, 2, 3, 9])
def test_confidence_stays_in_unit_interval(resolver, method, n):
    """Every method/support combination must produce a score in [0, 1]."""
    score = resolver._compute_confidence(
        method, 1.0, fuzzy_score=1.0, tanimoto=1.0, query_score=1.0, n_source_support=n
    )
    assert 0.0 <= score <= 1.0


# ─────────────────────────────────────────────────────────────────────────────
# Hit filtering — group structures and min_source_support
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("*C(=O)CCCC=CCC=CCCCCC", True),   # ChEBI "…dienoyl group"
        ("[*]CC", True),
        ("CC1OC(C)OC(C)OC(C)O1", False),   # metaldehyde
        (None, False),
        ("", False),
    ],
)
def test_attachment_point_detection(smiles, expected):
    """A dummy atom marks a substituent, which is never a queryable substance."""
    assert _has_attachment_point(smiles) is expected


def _pool_entry(resolver, source_key, smiles, inchikey, name):
    """Build one tagged candidate for ``_finalise_hits`` without touching a database."""
    from provesid.tools import _make_candidate

    cand = _make_candidate(source_key.title(), name=name, smiles=smiles, inchikey=inchikey)
    return resolver._tag_candidate(cand, source_key, 0, "exact_cas", 1.0)


def test_group_structures_never_win_a_lookup(resolver):
    """A group record must not be returned as the structure for an identifier."""
    template = resolver._empty_result("108-62-3", "CAS")
    pool = [
        _pool_entry(
            resolver, "chebi", "*C(=O)CCCC=CCC=CCCCCC", None, "dienoyl group"
        ),
        _pool_entry(
            resolver,
            "comptox",
            "CC1OC(C)OC(C)OC(C)O1",
            "GKKDCARASOJPNG-UHFFFAOYSA-N",
            "Metaldehyde",
        ),
    ]
    hits = resolver._finalise_hits(template, pool, 1, 0.0)

    assert len(hits) == 1
    assert _skeleton(hits[0]["InChIKey"]) == "GKKDCARASOJPNG"
    assert "*" not in str(hits[0]["SMILES"])


def test_group_query_still_resolves_to_the_group(resolver):
    """Querying a group SMILES must not filter the group away."""
    template = resolver._empty_result("*C(=O)CCCC=CCC=CCCCCC", "SMILES")
    pool = [_pool_entry(resolver, "chebi", "*C(=O)CCCC=CCC=CCCCCC", None, "dienoyl group")]

    hits = resolver._finalise_hits(template, pool, 1, 0.0)
    assert hits[0]["name"] == "dienoyl group"


def test_min_source_support_drops_uncorroborated_hits(resolver):
    """``min_source_support`` must remove hits carried by too few databases."""
    template = resolver._empty_result("108-62-3", "CAS")
    pool = [
        _pool_entry(
            resolver, "chebi", "CCO", "LFQSCWFLJHTTHZ-UHFFFAOYSA-N", "ethanol"
        )
    ]

    kept = resolver._finalise_hits(template, list(pool), 1, 0.0, min_source_support=1)
    assert kept[0]["name"] == "ethanol"

    dropped = resolver._finalise_hits(template, list(pool), 1, 0.0, min_source_support=2)
    assert dropped[0]["name"] is None, "a single-source hit must not survive min_source_support=2"


# ─────────────────────────────────────────────────────────────────────────────
# Offline-database fixtures
# ─────────────────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def comptox():
    """Yield an offline CompToxID client, skipping if its database is absent."""
    from provesid.comptox import CompToxID

    try:
        client = CompToxID()
        client.conn.execute("SELECT 1 FROM chemicals LIMIT 1").fetchone()
    except Exception as exc:  # pragma: no cover - environment dependent
        pytest.skip(f"CompTox offline database unavailable: {exc}")
    return client


@pytest.fixture(scope="module")
def chebi():
    """Yield an offline ChebiSDF client, skipping if the SDF is absent."""
    from provesid.chebi import ChebiSDF

    try:
        return ChebiSDF(auto_download=False)
    except Exception as exc:  # pragma: no cover - environment dependent
        pytest.skip(f"ChEBI SDF unavailable: {exc}")


@pytest.fixture(scope="module")
def all_sources():
    """Skip unless every offline source Search resolves against is present."""
    s = Search("cas", show_progress=False)
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
    return s


@pytest.fixture(scope="module")
def cas_truth_sample(comptox):
    """Sample CompTox rows whose CASRN maps to exactly one structure.

    Rows are taken on a ``rowid`` stride (no RNG, so the sample is reproducible)
    and spread across the whole table.  CASRNs that CompTox itself maps to more
    than one structure are dropped, since there is then no single truth to
    compare against.

    Returns:
        List of dicts with ``CASRN``, ``PREFERRED_NAME`` and ``INCHIKEY``.
    """
    total = comptox.conn.execute("SELECT COUNT(*) FROM chemicals").fetchone()[0]
    step = max(1, total // (CAS_SAMPLE_SIZE * 4))
    rows = [
        dict(r)
        for r in comptox.conn.execute(
            f"""
            SELECT CASRN, PREFERRED_NAME, INCHIKEY, SMILES
            FROM chemicals
            WHERE CASRN IS NOT NULL AND CASRN != ''
              AND INCHIKEY IS NOT NULL AND INCHIKEY != ''
              AND (rowid % {step}) = 1
            LIMIT {CAS_SAMPLE_SIZE * 4}
            """
        )
    ]

    structures = {}
    for row in rows:
        structures.setdefault(row["CASRN"], set()).add(_skeleton(row["INCHIKEY"]))

    sample, seen = [], set()
    for row in rows:
        cas = row["CASRN"]
        if len(structures[cas]) > 1 or cas in seen:
            continue
        seen.add(cas)
        sample.append(row)
        if len(sample) >= CAS_SAMPLE_SIZE:
            break

    if len(sample) < 20:  # pragma: no cover - environment dependent
        pytest.skip("Not enough unambiguous CompTox CAS rows to run the regression.")
    return sample


def _classify(sample, result_df):
    """Split resolved rows into correct / wrong / not-found against the truth.

    Args:
        sample: List of CompTox truth dicts (``INCHIKEY`` is the truth).
        result_df: ``Search.search`` output, one row per sample entry, in order.

    Returns:
        Tuple ``(correct, not_found, wrong)`` where ``wrong`` is a list of
        detail tuples for the failure message.
    """
    correct = not_found = 0
    wrong = []
    for truth, (_, got) in zip(sample, result_df.iterrows()):
        found = got.get("InChIKey")
        if not isinstance(found, str):
            not_found += 1
        elif _skeleton(found) == _skeleton(truth["INCHIKEY"]):
            correct += 1
        else:
            wrong.append(
                (
                    truth.get("CASRN") or truth.get("query"),
                    truth["PREFERRED_NAME"],
                    truth["INCHIKEY"],
                    got.get("name"),
                    found,
                    got.get("source"),
                    got.get("n_source_support"),
                    got.get("confidence"),
                )
            )
    return correct, not_found, wrong


# ─────────────────────────────────────────────────────────────────────────────
# ChEBI record integrity — guards the offset drift
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.integration
def test_chebi_records_round_trip_by_id(chebi):
    """Every indexed ChEBI ID must read back the record it points at.

    With text-mode offsets this failed for every record past the first CRLF line
    in the SDF — silently, because a neighbouring record parses perfectly well.
    """
    ids = list(chebi.index["id_to_offset"])
    step = max(1, len(ids) // CHEBI_SAMPLE_SIZE)
    sampled = ids[::step][:CHEBI_SAMPLE_SIZE]

    mismatched = []
    for chebi_id in sampled:
        record = chebi.get_compound_by_id(chebi_id)
        if not record or record.get("ChEBI ID") != chebi_id:
            mismatched.append((chebi_id, (record or {}).get("ChEBI ID")))

    assert not mismatched, (
        f"{len(mismatched)} of {len(sampled)} ChEBI records read back a different "
        f"compound (offset drift). Examples: {mismatched[:5]}"
    )


@pytest.mark.integration
def test_chebi_cas_lookup_returns_a_record_carrying_that_cas(chebi):
    """A CAS lookup must return records that actually list that CAS."""
    cas_numbers = list(chebi.index["cas_to_ids"])
    step = max(1, len(cas_numbers) // 200)
    sampled = cas_numbers[::step][:200]

    mismatched = []
    for cas in sampled:
        for record in chebi.search_by_cas(cas):
            listed = {c.strip() for c in (record.get("CAS Registry Numbers") or "").split(";")}
            if cas not in listed:
                mismatched.append((cas, record.get("ChEBI NAME"), sorted(listed)[:3]))

    assert not mismatched, (
        f"{len(mismatched)} ChEBI CAS lookups returned a record that does not carry "
        f"the queried CAS. Examples: {mismatched[:5]}"
    )


@pytest.mark.integration
def test_chebi_metaldehyde_regression(chebi):
    """CAS 108-62-3 is metaldehyde, not Mefluidide (the reported symptom)."""
    hits = chebi.search_by_cas("108-62-3")
    assert hits, "108-62-3 not found in ChEBI"
    assert hits[0]["ChEBI NAME"].lower() == "metaldehyde"
    assert hits[0]["INCHIKEY"].startswith("GKKDCARASOJPNG")


# ─────────────────────────────────────────────────────────────────────────────
# Reported failures — the 18 CAS numbers that used to resolve to the wrong compound
# ─────────────────────────────────────────────────────────────────────────────

# CAS -> (truth InChIKey skeleton, CompTox preferred name).  Independently
# confirmed for all 18 by the OECD QSAR Toolbox 4.8.2 WebAPI in the study that
# found the bug.
REPORTED_WRONG = {
    "10380-28-6": ("YXLXNENXOJSQEI", "Copper-8-hydroxyquinoline"),
    "100784-20-1": ("FMGZEUWROYGLAY", "Halosulfuron-methyl"),
    "108-62-3": ("GKKDCARASOJPNG", "Metaldehyde"),
    "110488-70-5": ("QNBTYORWCCMPQP", "Dimethomorph"),
    "119446-68-3": ("BQYJATMQXGBDHF", "Difenoconazole"),
    "123312-89-0": ("QHMTXANCGGJZRX", "Pymetrozine"),
    "125401-92-5": ("FUHMZYWBSHTEDZ", "Bispyribac-sodium"),
    "23947-60-6": ("BBXXLROWFHWFQY", "Ethirimol"),
    "25606-41-1": ("MKIMSXGUTQTKJU", "Propamocarb hydrochloride"),
    "33629-47-9": ("SPNQRCTZKIBOAX", "Butralin"),
    "34643-46-4": ("FITIWKDOCAUBQD", "Prothiofos"),
    "38641-94-0": ("ZEKANFGSDXODPD", "Glyphosate isopropylamine"),
    "41814-78-2": ("DQJCHOQLCLEDLL", "Tricyclazole"),
    "57966-95-7": ("XERJKGMBORTKEO", "Cymoxanil"),
    "66230-04-4": ("NYPJDWWKZLNGGM", "Esfenvalerate"),
    "70630-17-0": ("ZQEIXNIJLIKNTD", "Metalaxyl-M"),
    "95266-40-3": ("RVKCCVTVZORVGD", "Trinexapac-ethyl"),
    "99129-21-2": ("SILSDTWXNBZOGF", "Clethodim"),
}


@pytest.mark.integration
@pytest.mark.slow
def test_reported_wrong_cas_numbers_all_resolve_correctly(all_sources):
    """All 18 CAS numbers from the bug report must return the right structure."""
    result = Search("cas", show_progress=False).search(list(REPORTED_WRONG))

    wrong = []
    for query, row in zip(result["query"], result.to_dict("records")):
        expected, label = REPORTED_WRONG[query]
        if _skeleton(row.get("InChIKey")) != expected:
            wrong.append((query, label, expected, row.get("name"), row.get("InChIKey")))

    assert not wrong, f"{len(wrong)} of {len(REPORTED_WRONG)} still wrong: {wrong}"


# ─────────────────────────────────────────────────────────────────────────────
# CompTox-truth sampling — CAS, name and InChIKey resolution
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.integration
@pytest.mark.slow
def test_cas_resolution_matches_comptox_truth(all_sources, cas_truth_sample):
    """Resolving a CAS across all sources must reproduce CompTox's structure."""
    result = Search("cas", show_progress=False).search(
        [row["CASRN"] for row in cas_truth_sample]
    )
    correct, not_found, wrong = _classify(cas_truth_sample, result)

    rate = len(wrong) / len(cas_truth_sample)
    assert rate <= MAX_WRONG_RATE, (
        f"{len(wrong)}/{len(cas_truth_sample)} CAS numbers resolved to a different "
        f"structure than CompTox (correct={correct}, not_found={not_found}). "
        f"Examples: {wrong[:5]}"
    )


@pytest.mark.integration
@pytest.mark.slow
def test_corroborated_cas_hits_are_never_wrong(all_sources, cas_truth_sample):
    """A hit three or more databases agree on must match the truth.

    This is the contract the confidence formula now encodes: corroboration is
    evidence.  If a structure that three independent databases carry can still
    be wrong, the ranking is not trustworthy at any confidence level.
    """
    result = Search("cas", show_progress=False).search(
        [row["CASRN"] for row in cas_truth_sample]
    )

    wrong = []
    for truth, (_, got) in zip(cas_truth_sample, result.iterrows()):
        found = got.get("InChIKey")
        if not isinstance(found, str) or (got.get("n_source_support") or 0) < 3:
            continue
        if _skeleton(found) != _skeleton(truth["INCHIKEY"]):
            wrong.append(
                (truth["CASRN"], truth["PREFERRED_NAME"], truth["INCHIKEY"],
                 got.get("name"), found, got.get("n_source_support"))
            )

    assert not wrong, (
        f"{len(wrong)} hit(s) corroborated by 3+ databases disagree with CompTox: "
        f"{wrong[:5]}"
    )


@pytest.mark.integration
@pytest.mark.slow
def test_returned_hits_are_the_most_corroborated_available(all_sources, cas_truth_sample):
    """The top hit must never be less corroborated than an equally-matched rival.

    The reported bug in one line: a ``n_source_support == 1`` hit was returned
    while a ``n_source_support == 3`` hit sat below it in the same result set.

    Compared within one ``match_method`` — a stronger method legitimately
    outranks a weaker one however many databases carry it (an exact InChIKey hit
    should beat a formula match), so only like-for-like rivals are a violation.
    """
    queries = [row["CASRN"] for row in cas_truth_sample[:60]]
    result = Search("cas", show_progress=False).search(queries, n_hits="all")

    offenders = []
    for query, group in result.groupby("query", sort=False):
        ranked = group.sort_values("hit_rank")
        if len(ranked) < 2:
            continue
        top = ranked.iloc[0]
        rivals = ranked[ranked["match_method"] == top["match_method"]]
        best_support = rivals["n_source_support"].max()
        if (top["n_source_support"] or 0) < best_support:
            offenders.append(
                (query, top["name"], int(top["n_source_support"] or 0), int(best_support))
            )

    assert not offenders, (
        "top hit is less corroborated than a lower-ranked hit for "
        f"{len(offenders)} quer(ies): {offenders[:5]}"
    )


@pytest.mark.integration
@pytest.mark.slow
def test_name_resolution_matches_comptox_truth(all_sources, comptox, cas_truth_sample):
    """Resolving CompTox preferred names must reproduce CompTox's structure."""
    sample = [
        row
        for row in cas_truth_sample
        if row["PREFERRED_NAME"] and any(ch.isalpha() for ch in row["PREFERRED_NAME"])
    ][:NAME_SAMPLE_SIZE]

    result = Search("name", show_progress=False).search(
        [row["PREFERRED_NAME"] for row in sample]
    )
    correct, not_found, wrong = _classify(sample, result)

    rate = len(wrong) / len(sample)
    assert rate <= MAX_WRONG_RATE, (
        f"{len(wrong)}/{len(sample)} names resolved to a different structure than "
        f"CompTox (correct={correct}, not_found={not_found}). Examples: {wrong[:5]}"
    )


@pytest.mark.integration
@pytest.mark.slow
def test_inchikey_resolution_is_structure_preserving(all_sources, cas_truth_sample):
    """An InChIKey query must come back as the very same structure.

    The strongest available identity check: unlike a CAS or a name, an InChIKey
    *is* the structure, so any disagreement here is a resolver defect rather than
    a cross-database annotation difference.
    """
    sample = cas_truth_sample[:NAME_SAMPLE_SIZE]
    result = Search("inchikey", show_progress=False).search(
        [row["INCHIKEY"] for row in sample]
    )
    correct, not_found, wrong = _classify(sample, result)

    assert not wrong, (
        f"{len(wrong)}/{len(sample)} InChIKey queries returned a different structure "
        f"(correct={correct}, not_found={not_found}). Examples: {wrong[:5]}"
    )


@pytest.mark.integration
@pytest.mark.slow
def test_dtxsid_resolution_matches_comptox_truth(all_sources, comptox):
    """A DTXSID must resolve to the structure CompTox files under it."""
    total = comptox.conn.execute("SELECT COUNT(*) FROM chemicals").fetchone()[0]
    step = max(1, total // (NAME_SAMPLE_SIZE * 2))
    sample = [
        dict(r)
        for r in comptox.conn.execute(
            f"""
            SELECT DTXSID, PREFERRED_NAME, CASRN, INCHIKEY
            FROM chemicals
            WHERE DTXSID IS NOT NULL AND DTXSID != ''
              AND INCHIKEY IS NOT NULL AND INCHIKEY != ''
              AND (rowid % {step}) = 1
            LIMIT {NAME_SAMPLE_SIZE}
            """
        )
    ]

    result = Search("dtxsid", show_progress=False).search(
        [row["DTXSID"] for row in sample]
    )
    correct, not_found, wrong = _classify(sample, result)

    assert not wrong, (
        f"{len(wrong)}/{len(sample)} DTXSIDs resolved to a different structure "
        f"(correct={correct}, not_found={not_found}). Examples: {wrong[:5]}"
    )


@pytest.mark.integration
@pytest.mark.slow
def test_smiles_resolution_is_structure_preserving(all_sources, comptox, cas_truth_sample):
    """A SMILES query must return the compound that SMILES describes."""
    sample = [row for row in cas_truth_sample if row.get("SMILES")][:NAME_SAMPLE_SIZE]
    if len(sample) < 10:  # pragma: no cover - environment dependent
        pytest.skip("Not enough CompTox rows with a SMILES in the sample.")

    result = Search("smiles", show_progress=False).search(
        [row["SMILES"] for row in sample]
    )
    correct, not_found, wrong = _classify(sample, result)

    assert not wrong, (
        f"{len(wrong)}/{len(sample)} SMILES queries returned a different structure "
        f"(correct={correct}, not_found={not_found}). Examples: {wrong[:5]}"
    )


@pytest.mark.integration
@pytest.mark.slow
def test_inchi_resolution_is_structure_preserving(all_sources, cas_truth_sample):
    """An InChI query must return the compound that InChI describes."""
    from provesid.search import normalize_structure

    sample, inchis = [], []
    for row in cas_truth_sample:
        norm = normalize_structure(row.get("SMILES"))
        if norm["inchi"] and _skeleton(norm["inchikey"]) == _skeleton(row["INCHIKEY"]):
            sample.append(row)
            inchis.append(norm["inchi"])
        if len(sample) >= 40:
            break
    if len(sample) < 10:  # pragma: no cover - environment dependent
        pytest.skip("RDKit could not derive enough InChIs from the sample.")

    result = Search("inchi", show_progress=False).search(inchis)
    correct, not_found, wrong = _classify(sample, result)

    assert not wrong, (
        f"{len(wrong)}/{len(sample)} InChI queries returned a different structure "
        f"(correct={correct}, not_found={not_found}). Examples: {wrong[:5]}"
    )


# ─────────────────────────────────────────────────────────────────────────────
# Search-based helpers — enrich() and resolve_cascade()
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.integration
@pytest.mark.slow
def test_enrich_matches_direct_search(all_sources, cas_truth_sample):
    """``enrich`` must attach exactly what ``search`` resolves, row for row."""
    import pandas as pd

    sample = cas_truth_sample[:40]
    frame = pd.DataFrame(
        {
            "CAS": [row["CASRN"] for row in sample] * 2,  # duplicated rows
            "measurement": list(range(len(sample) * 2)),
        }
    )

    enriched = Search("cas", show_progress=False).enrich(frame, "CAS")
    assert len(enriched) == len(frame)
    assert list(enriched["measurement"]) == list(frame["measurement"])

    truth = {row["CASRN"]: _skeleton(row["INCHIKEY"]) for row in sample}
    wrong = [
        (cas, got)
        for cas, got in zip(enriched["CAS"], enriched["provesid_InChIKey"])
        if isinstance(got, str) and _skeleton(got) != truth[cas]
    ]
    rate = len(wrong) / len(enriched)
    assert rate <= MAX_WRONG_RATE, f"enrich attached wrong structures: {wrong[:5]}"


@pytest.mark.integration
@pytest.mark.slow
def test_resolve_cascade_prefers_the_reliable_identifier(all_sources, cas_truth_sample):
    """A CAS→name cascade must resolve rows to CompTox's structure.

    Half the rows carry only a name, so this also exercises the fall-through
    path: a row the CAS stage cannot resolve must still be resolved by name.
    """
    import pandas as pd

    sample = cas_truth_sample[:40]
    frame = pd.DataFrame(
        {
            "CASRN": [r["CASRN"] if i % 2 == 0 else "" for i, r in enumerate(sample)],
            "name": [r["PREFERRED_NAME"] for r in sample],
        }
    )

    out = resolve_cascade(
        frame,
        stages=[
            ("cas", Search("cas", show_progress=False), "CASRN"),
            ("name", Search("name", show_progress=False), "name"),
        ],
    )

    wrong = [
        (sample[i]["CASRN"], sample[i]["PREFERRED_NAME"], got)
        for i, got in enumerate(out["provesid_InChIKey"])
        if isinstance(got, str) and _skeleton(got) != _skeleton(sample[i]["INCHIKEY"])
    ]
    rate = len(wrong) / len(out)
    assert rate <= MAX_WRONG_RATE, f"cascade resolved wrong structures: {wrong[:5]}"

    resolved_by = out["provesid_resolved_by"].value_counts().to_dict()
    assert resolved_by.get("name", 0) > 0, (
        f"no row fell through to the name stage: {resolved_by}"
    )
