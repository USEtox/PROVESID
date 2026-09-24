"""
Tests for ``CheMBL.search_by_name`` — the ``UNION`` rewrite and its ordering.

The method used to issue ``SELECT DISTINCT … FROM molecule_dictionary LEFT JOIN
molecule_synonyms … WHERE LOWER(pref_name) = ? OR LOWER(synonyms) = ?`` with no
``ORDER BY``.  That shape has two defects:

1. An ``OR`` whose arms live in different tables across a ``LEFT JOIN`` cannot
   use an index for either arm, so every call scanned all 2.9 M rows of
   ``molecule_dictionary`` — 743 ms per exact lookup on a real release.
2. ``LIMIT`` with no ``ORDER BY`` returns whatever the query plan produces, so
   the *content* of a truncated result set is not reproducible.

These tests pin the rewrite: the same compounds come back, in ``molregno``
order, ``limit`` keeps the lowest ``molregno`` values, and on a database with
the ``lower(...)`` expression indexes the plan is a pair of index searches
rather than a scan.  They build a miniature ChEMBL, so they need neither the
30 GB database nor network access.
"""

import sqlite3

import pytest

from provesid import CheMBL


# ── A miniature ChEMBL built for name search ─────────────────────────────────

#: The query ``search_by_name`` used to issue, kept here so the tests can assert
#: that the rewrite returns the same compounds as the shape it replaced.
LEGACY_EXACT_QUERY = """
SELECT DISTINCT md.molregno
FROM molecule_dictionary md
LEFT JOIN molecule_synonyms ms ON md.molregno = ms.molregno
WHERE LOWER(md.pref_name) = LOWER(?)
   OR LOWER(ms.synonyms) = LOWER(?)
LIMIT ?
"""

LEGACY_SUBSTRING_QUERY = LEGACY_EXACT_QUERY.replace("= LOWER(?)", "LIKE LOWER(?)")


def _make_chembl(path, *, with_indexes=True):
    """Build a small ChEMBL carrying the tables ``search_by_name`` touches.

    Names are laid out to exercise the cases that matter:

    * ``SHARED NAME`` is the preferred name of 40 compounds, so a ``limit``
      smaller than the match count has to choose *which* ones — the point of
      the ``ORDER BY``.
    * ``molregno`` is deliberately not insertion-ordered, so a result that
      happens to be sorted cannot be an accident of physical row order.
    * ``AMBIGUOUS`` is one compound's preferred name and a *different*
      compound's synonym, which is the case the two ``UNION`` arms must merge.
    """
    conn = sqlite3.connect(path)
    conn.executescript(
        """
        CREATE TABLE molecule_dictionary (
            molregno INTEGER PRIMARY KEY, pref_name TEXT, chembl_id TEXT,
            max_phase NUMERIC, therapeutic_flag INTEGER, molecule_type TEXT);
        CREATE TABLE compound_structures (
            molregno INTEGER PRIMARY KEY, standard_inchi TEXT,
            standard_inchi_key TEXT, canonical_smiles TEXT);
        CREATE TABLE compound_properties (
            molregno INTEGER PRIMARY KEY, mw_freebase REAL);
        CREATE TABLE molecule_synonyms (
            molregno INTEGER, syn_type TEXT, molsyn_id INTEGER, synonyms TEXT);
        """
    )

    # Insert in descending molregno so that "sorted" is never the natural order.
    for molregno in range(200, 0, -1):
        if molregno <= 40:
            pref_name = "SHARED NAME"
        elif molregno == 77:
            pref_name = "AMBIGUOUS"
        else:
            pref_name = f"COMPOUND {molregno}"
        conn.execute(
            "INSERT INTO molecule_dictionary VALUES (?,?,?,?,?,?)",
            (molregno, pref_name, f"CHEMBL{molregno}", 4, 1, "Small molecule"),
        )
        conn.execute(
            "INSERT INTO compound_structures VALUES (?,?,?,?)",
            (molregno, f"InChI=1S/C{molregno}", f"KEY{molregno:016d}-A-N", "C" * 3),
        )
        conn.execute(
            "INSERT INTO compound_properties VALUES (?,?)", (molregno, 100.0 + molregno)
        )
        conn.execute(
            "INSERT INTO molecule_synonyms VALUES (?,?,?,?)",
            (molregno, "TRADE_NAME", molregno, f"synonym-of-{molregno}"),
        )

    # A synonym-only match, and a name that is a preferred name on one compound
    # and a synonym on another — the two arms of the UNION have to merge.
    conn.execute(
        "INSERT INTO molecule_synonyms VALUES (150,'TRADE_NAME',9001,'Panacea')"
    )
    conn.execute(
        "INSERT INTO molecule_synonyms VALUES (12,'TRADE_NAME',9002,'AMBIGUOUS')"
    )
    if with_indexes:
        conn.executescript(
            """
            CREATE INDEX ix_md_pref_lower ON molecule_dictionary(lower(pref_name));
            CREATE INDEX ix_ms_syn_lower ON molecule_synonyms(lower(synonyms));
            """
        )
    conn.commit()
    conn.close()
    return str(path)


@pytest.fixture
def chembl(tmp_path):
    """A ``CheMBL`` bound to the miniature database, indexes and all."""
    path = _make_chembl(tmp_path / "chembl_36.db")
    return CheMBL(db_path=path, auto_download=False)


def _molregnos(results):
    return [r["molregno"] for r in results]


# ── The compounds that come back ─────────────────────────────────────────────

class TestWhatIsFound:

    def test_exact_matches_a_preferred_name(self, chembl):
        assert _molregnos(chembl.search_by_name("COMPOUND 123", exact=True)) == [123]

    def test_exact_matches_a_synonym(self, chembl):
        assert _molregnos(chembl.search_by_name("Panacea", exact=True)) == [150]

    def test_both_arms_are_merged(self, chembl):
        """A name that is one compound's pref_name and another's synonym."""
        assert _molregnos(chembl.search_by_name("AMBIGUOUS", exact=True)) == [12, 77]

    def test_a_compound_matching_in_both_arms_appears_once(self, chembl):
        """molregno 150 matches its own synonym twice over; UNION dedupes."""
        assert _molregnos(chembl.search_by_name("synonym-of-150", exact=True)) == [150]

    def test_exact_is_case_insensitive(self, chembl):
        for spelling in ("panacea", "PANACEA", "PaNaCeA"):
            assert _molregnos(chembl.search_by_name(spelling, exact=True)) == [150]

    def test_exact_rejects_a_substring(self, chembl):
        assert chembl.search_by_name("anace", exact=True) == []

    def test_substring_matches_inside_a_name(self, chembl):
        assert _molregnos(chembl.search_by_name("anace")) == [150]

    def test_substring_matches_both_arms(self, chembl):
        assert _molregnos(chembl.search_by_name("AMBIGU")) == [12, 77]

    def test_no_match_is_an_empty_list(self, chembl):
        assert chembl.search_by_name("NOT A COMPOUND", exact=True) == []

    def test_results_carry_structures_and_synonyms(self, chembl):
        (compound,) = chembl.search_by_name("Panacea", exact=True)
        assert compound["chembl_id"] == "CHEMBL150"
        assert compound["canonical_smiles"] == "CCC"
        assert "Panacea" in compound["synonyms"]


# ── Ordering and limit ───────────────────────────────────────────────────────

class TestOrderingIsDeterministic:

    def test_results_are_ordered_by_molregno(self, chembl):
        found = _molregnos(chembl.search_by_name("SHARED NAME", exact=True))
        assert found == sorted(found)
        assert found == list(range(1, 41))

    def test_limit_keeps_the_lowest_molregnos(self, chembl):
        """Not merely *some* five rows — the same five on every call."""
        assert _molregnos(
            chembl.search_by_name("SHARED NAME", limit=5, exact=True)
        ) == [1, 2, 3, 4, 5]

    def test_a_truncated_result_is_stable_across_calls(self, chembl):
        first = _molregnos(chembl.search_by_name("SHARED NAME", limit=7, exact=True))
        for _ in range(3):
            assert (
                _molregnos(chembl.search_by_name("SHARED NAME", limit=7, exact=True))
                == first
            )

    def test_a_truncated_result_survives_a_vacuum(self, chembl):
        """The defect the ORDER BY fixes: physical row order used to decide."""
        before = _molregnos(chembl.search_by_name("SHARED NAME", limit=7, exact=True))
        chembl.conn.execute("VACUUM")
        after = _molregnos(chembl.search_by_name("SHARED NAME", limit=7, exact=True))
        assert before == after

    def test_substring_results_are_ordered_too(self, chembl):
        found = _molregnos(chembl.search_by_name("SHARED"))
        assert found == sorted(found)


# ── Equivalence with the shape it replaced ───────────────────────────────────

class TestSameCompoundsAsTheLegacyQuery:
    """The rewrite may reorder, but it must not change *which* compounds match.

    Compared as sets and with a limit high enough that neither query truncates,
    because the legacy query's order — and so its truncation — is exactly what
    is not reproducible.
    """

    NO_TRUNCATION = 10_000

    def _legacy(self, chembl, query, term):
        return {
            row[0]
            for row in chembl.conn.execute(query, (term, term, self.NO_TRUNCATION))
        }

    @pytest.mark.parametrize(
        "name",
        ["COMPOUND 1", "COMPOUND 199", "SHARED NAME", "AMBIGUOUS", "Panacea",
         "synonym-of-42", "panacea", "nothing at all"],
    )
    def test_exact_finds_the_same_compounds(self, chembl, name):
        rewritten = set(_molregnos(chembl.search_by_name(name, exact=True,
                                                         limit=self.NO_TRUNCATION)))
        assert rewritten == self._legacy(chembl, LEGACY_EXACT_QUERY, name)

    @pytest.mark.parametrize(
        "fragment", ["COMPOUND 1", "SHARED", "AMBIGU", "anace", "synonym-of-4", "zzz"]
    )
    def test_substring_finds_the_same_compounds(self, chembl, fragment):
        rewritten = set(_molregnos(chembl.search_by_name(fragment,
                                                         limit=self.NO_TRUNCATION)))
        assert rewritten == self._legacy(
            chembl, LEGACY_SUBSTRING_QUERY, f"%{fragment}%"
        )


# ── The plan, which is the whole point of the rewrite ────────────────────────

class TestTheQueryPlan:

    def _plan(self, chembl, name, exact):
        """The plan SQLite chooses for the statement the method actually issues."""
        query = _issued_query(chembl, name, exact)
        return "\n".join(
            row[-1] for row in chembl.conn.execute("EXPLAIN QUERY PLAN " + query)
        )

    def test_exact_lookup_uses_the_expression_indexes(self, chembl):
        plan = self._plan(chembl, "COMPOUND 123", exact=True)
        assert "ix_md_pref_lower" in plan
        assert "ix_ms_syn_lower" in plan
        assert "SCAN molecule_dictionary" not in plan

    def test_exact_lookup_without_the_indexes_still_answers(self, tmp_path):
        """A full ChEMBL release has no expression indexes; it scans, correctly."""
        path = _make_chembl(tmp_path / "chembl_36.db", with_indexes=False)
        unindexed = CheMBL(db_path=path, auto_download=False)
        assert _molregnos(unindexed.search_by_name("Panacea", exact=True)) == [150]

    def test_no_arm_of_the_query_joins_the_two_tables(self, chembl):
        """The LEFT JOIN is what defeated indexing; it must be gone."""
        plan = self._plan(chembl, "COMPOUND 123", exact=True)
        assert "LEFT-JOIN" not in plan


def _issued_query(chembl, name, exact):
    """Recover the SQL ``search_by_name`` executes, via SQLite's trace hook.

    The hook reports the statement with its parameters already substituted, so
    the result is runnable — and, more to the point, plannable — as it stands.
    The first statement of the call is the name lookup; the rest are the
    per-result ``get_compound`` queries.
    """
    statements = []
    chembl.conn.set_trace_callback(statements.append)
    try:
        chembl.search_by_name(name, exact=exact)
    finally:
        chembl.conn.set_trace_callback(None)
    return statements[0]
