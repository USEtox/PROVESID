"""Unit tests for the dataset-enrichment helpers in provesid.search.

Covers :meth:`Search.enrich`, :func:`resolve_cascade` and :func:`mw_within`.
Uses a fake resolver rather than real databases, so the tests are fast and
assert the *bookkeeping* — which rows go to which stage, what gets accepted,
how columns and the index come back — independently of any dataset content.

Run with::

    uv run pytest tests/test_search_enrichment.py -v
"""

from __future__ import annotations

import pandas as pd
import pytest

from provesid.search import OUTPUT_COLUMNS, Search, mw_within, resolve_cascade

# Real structures, so mw_within has something honest to compute on.
ETHANOL = {
    "name": "ethanol",
    "SMILES": "CCO",
    "InChIKey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
}
BENZENE = {
    "name": "benzene",
    "SMILES": "c1ccccc1",
    "InChIKey": "UHOVQNZJYSORNB-UHFFFAOYSA-N",
}
FORMALDEHYDE = {
    "name": "formaldehyde",
    "SMILES": "C=O",
    "InChIKey": "WSFSSNUMVMOOMR-UHFFFAOYSA-N",
}


class FakeSearch:
    """Stand-in for Search that answers from a dict and records its queries.

    Args:
        answers: Maps a query string to the compound dict to return. Anything
            not in the mapping resolves to an empty row.
    """

    def __init__(self, answers):
        self.answers = answers
        self.calls = []

    def search(self, queries, *, n_hits=None, **_kwargs):
        if isinstance(queries, str):
            queries = [queries]
        self.calls.append(list(queries))
        rows = []
        for query in queries:
            row = {column: None for column in OUTPUT_COLUMNS}
            row["query"] = query
            row.update(self.answers.get(query, {}))
            rows.append(row)
        return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)


# ─────────────────────────────────────────────────────────────────────────────
# Search.enrich
# ─────────────────────────────────────────────────────────────────────────────

class TestEnrich:
    def _search(self, answers):
        """A real Search whose search() is swapped for the fake resolver."""
        searcher = Search("cas", show_progress=False)
        fake = FakeSearch(answers)
        searcher.search = fake.search
        return searcher, fake

    def test_each_distinct_value_is_searched_once(self):
        df = pd.DataFrame({"CAS": ["64-17-5", "64-17-5", "71-43-2", "64-17-5"]})
        searcher, fake = self._search({"64-17-5": ETHANOL, "71-43-2": BENZENE})

        searcher.enrich(df, "CAS")

        assert len(fake.calls) == 1
        assert sorted(fake.calls[0]) == ["64-17-5", "71-43-2"]

    def test_result_is_broadcast_to_every_matching_row(self):
        df = pd.DataFrame({"CAS": ["64-17-5", "71-43-2", "64-17-5"]})
        searcher, _ = self._search({"64-17-5": ETHANOL, "71-43-2": BENZENE})

        out = searcher.enrich(df, "CAS")

        assert out["provesid_name"].tolist() == ["ethanol", "benzene", "ethanol"]

    def test_row_count_and_order_are_unchanged(self):
        df = pd.DataFrame({"CAS": ["64-17-5", "71-43-2", "64-17-5", "71-43-2"]})
        searcher, _ = self._search({"64-17-5": ETHANOL, "71-43-2": BENZENE})

        out = searcher.enrich(df, "CAS")

        assert len(out) == len(df)
        assert out["CAS"].tolist() == df["CAS"].tolist()

    def test_original_index_is_preserved(self):
        df = pd.DataFrame({"CAS": ["64-17-5", "71-43-2"]}, index=["x", "y"])
        searcher, _ = self._search({"64-17-5": ETHANOL, "71-43-2": BENZENE})

        out = searcher.enrich(df, "CAS")

        assert out.index.tolist() == ["x", "y"]

    def test_original_columns_are_untouched(self):
        df = pd.DataFrame({"CAS": ["64-17-5"], "name": ["EtOH"], "bp": [78.4]})
        searcher, _ = self._search({"64-17-5": ETHANOL})

        out = searcher.enrich(df, "CAS")

        assert out["name"].tolist() == ["EtOH"]      # not the resolved name
        assert out["bp"].tolist() == [78.4]
        assert out["provesid_name"].tolist() == ["ethanol"]

    def test_every_output_column_is_added_with_the_prefix(self):
        df = pd.DataFrame({"CAS": ["64-17-5"]})
        searcher, _ = self._search({"64-17-5": ETHANOL})

        out = searcher.enrich(df, "CAS")

        for column in OUTPUT_COLUMNS:
            assert f"provesid_{column}" in out.columns

    def test_custom_prefix(self):
        df = pd.DataFrame({"CAS": ["64-17-5"]})
        searcher, _ = self._search({"64-17-5": ETHANOL})

        out = searcher.enrich(df, "CAS", prefix="id_")

        assert out["id_name"].tolist() == ["ethanol"]
        assert not any(c.startswith("provesid_") for c in out.columns)

    def test_empty_values_are_not_searched(self):
        df = pd.DataFrame({"CAS": ["64-17-5", "", None, "nan", "  "]})
        searcher, fake = self._search({"64-17-5": ETHANOL})

        out = searcher.enrich(df, "CAS")

        assert fake.calls[0] == ["64-17-5"]
        assert out["provesid_name"].tolist()[0] == "ethanol"
        assert out["provesid_name"].isna().tolist()[1:] == [True, True, True, True]

    def test_unresolved_values_give_empty_columns_not_dropped_rows(self):
        df = pd.DataFrame({"CAS": ["64-17-5", "nonsense"]})
        searcher, _ = self._search({"64-17-5": ETHANOL})

        out = searcher.enrich(df, "CAS")

        assert len(out) == 2
        assert pd.isna(out["provesid_InChIKey"].iloc[1])

    def test_column_with_no_usable_values_still_adds_columns(self):
        df = pd.DataFrame({"CAS": ["", None]})
        searcher, fake = self._search({})

        out = searcher.enrich(df, "CAS")

        assert fake.calls == []                       # nothing searched
        assert len(out) == 2
        assert "provesid_InChIKey" in out.columns
        assert out["provesid_InChIKey"].isna().all()

    def test_missing_column_raises_keyerror(self):
        searcher, _ = self._search({})
        with pytest.raises(KeyError, match="absent"):
            searcher.enrich(pd.DataFrame({"CAS": ["64-17-5"]}), "absent")

    def test_prefix_collision_raises_valueerror(self):
        df = pd.DataFrame({"CAS": ["64-17-5"], "provesid_name": ["stale"]})
        searcher, _ = self._search({"64-17-5": ETHANOL})
        with pytest.raises(ValueError, match="prefix"):
            searcher.enrich(df, "CAS")

    def test_enriching_twice_works_with_distinct_prefixes(self):
        df = pd.DataFrame({"CAS": ["64-17-5"]})
        searcher, _ = self._search({"64-17-5": ETHANOL})

        once = searcher.enrich(df, "CAS")
        twice = searcher.enrich(once, "CAS", prefix="second_")

        assert twice["provesid_name"].tolist() == ["ethanol"]
        assert twice["second_name"].tolist() == ["ethanol"]

    def test_values_are_matched_after_stripping_whitespace(self):
        df = pd.DataFrame({"CAS": [" 64-17-5 ", "64-17-5"]})
        searcher, fake = self._search({"64-17-5": ETHANOL})

        out = searcher.enrich(df, "CAS")

        assert fake.calls[0] == ["64-17-5"]
        assert out["provesid_name"].tolist() == ["ethanol", "ethanol"]


# ─────────────────────────────────────────────────────────────────────────────
# mw_within
# ─────────────────────────────────────────────────────────────────────────────

class TestMwWithin:
    def test_accepts_a_structurally_identical_hit(self):
        accept = mw_within(0.5)
        assert accept({"SMILES": "CCO"}, {"SMILES": "CCO"}) == ["mw", "smiles"]

    def test_rejects_a_different_compound(self):
        accept = mw_within(0.5)
        # formaldehyde (30.03) vs benzene (78.11)
        assert accept({"SMILES": "C=O"}, {"SMILES": "c1ccccc1"}) == []

    def test_reports_mw_only_when_smiles_differ_within_tolerance(self):
        """Same formula, different structure: MW agrees, canonical SMILES do not."""
        accept = mw_within(0.5)
        # ethanol vs dimethyl ether, both C2H6O (46.07)
        assert accept({"SMILES": "COC"}, {"SMILES": "CCO"}) == ["mw"]

    def test_tolerance_is_honoured(self):
        # benzene 78.11 vs pyridine 79.10 -> ~0.99 Da apart
        hit, row = {"SMILES": "c1ccncc1"}, {"SMILES": "c1ccccc1"}
        assert mw_within(0.5)(hit, row) == []
        assert mw_within(1.5)(hit, row) == ["mw"]

    def test_rejects_when_the_reference_structure_is_missing(self):
        accept = mw_within(0.5)
        assert accept({"SMILES": "CCO"}, {"SMILES": None}) == []
        assert accept({"SMILES": "CCO"}, {}) == []

    def test_rejects_when_the_hit_has_no_structure(self):
        accept = mw_within(0.5)
        assert accept({"SMILES": None}, {"SMILES": "CCO"}) == []

    def test_rejects_unparseable_smiles(self):
        accept = mw_within(0.5)
        assert accept({"SMILES": "not-a-smiles"}, {"SMILES": "CCO"}) == []
        assert accept({"SMILES": "CCO"}, {"SMILES": "not-a-smiles"}) == []

    def test_custom_reference_column(self):
        accept = mw_within(0.5, reference_column="canonical_SMILES")
        assert accept({"SMILES": "CCO"}, {"canonical_SMILES": "CCO"}) == ["mw", "smiles"]

    def test_name_agreement_is_reported_when_requested(self):
        accept = mw_within(0.5, name_column="name")
        hit = {"SMILES": "CCO", "name": "Ethanol"}
        assert accept(hit, {"SMILES": "CCO", "name": "ethanol"}) == ["mw", "smiles", "name"]

    def test_name_disagreement_does_not_reject(self):
        """A wrong name is recorded by omission, not treated as a failure."""
        accept = mw_within(0.5, name_column="name")
        hit = {"SMILES": "CCO", "name": "ethyl alcohol"}
        assert accept(hit, {"SMILES": "CCO", "name": "ethanol"}) == ["mw", "smiles"]


# ─────────────────────────────────────────────────────────────────────────────
# resolve_cascade
# ─────────────────────────────────────────────────────────────────────────────

class TestResolveCascade:
    def test_first_stage_wins_and_later_stages_are_skipped(self):
        df = pd.DataFrame({"CAS": ["64-17-5"], "name": ["ethanol"]})
        by_cas = FakeSearch({"64-17-5": ETHANOL})
        by_name = FakeSearch({"ethanol": ETHANOL})

        out = resolve_cascade(
            df, [("cas", by_cas, "CAS"), ("name", by_name, "name")]
        )

        assert out["provesid_resolved_by"].tolist() == ["cas"]
        assert by_name.calls == []               # never consulted

    def test_unresolved_rows_fall_through_to_the_next_stage(self):
        df = pd.DataFrame({"CAS": ["nonsense"], "name": ["ethanol"]})
        by_cas = FakeSearch({})
        by_name = FakeSearch({"ethanol": ETHANOL})

        out = resolve_cascade(
            df, [("cas", by_cas, "CAS"), ("name", by_name, "name")]
        )

        assert out["provesid_resolved_by"].tolist() == ["name"]
        assert out["provesid_name"].tolist() == ["ethanol"]

    def test_only_pending_rows_are_passed_to_later_stages(self):
        df = pd.DataFrame({
            "CAS": ["64-17-5", "nonsense", "71-43-2"],
            "name": ["ethanol", "benzene", "benzene"],
        })
        by_cas = FakeSearch({"64-17-5": ETHANOL, "71-43-2": BENZENE})
        by_name = FakeSearch({"benzene": BENZENE})

        resolve_cascade(df, [("cas", by_cas, "CAS"), ("name", by_name, "name")])

        assert by_cas.calls[0] == ["64-17-5", "nonsense", "71-43-2"]
        assert by_name.calls[0] == ["benzene"]   # only the row that failed

    def test_rows_with_an_empty_stage_column_skip_that_stage(self):
        df = pd.DataFrame({"CAS": ["", "64-17-5"], "name": ["benzene", "ethanol"]})
        by_cas = FakeSearch({"64-17-5": ETHANOL})
        by_name = FakeSearch({"benzene": BENZENE})

        out = resolve_cascade(df, [("cas", by_cas, "CAS"), ("name", by_name, "name")])

        assert by_cas.calls[0] == ["64-17-5"]
        assert out["provesid_resolved_by"].tolist() == ["name", "cas"]

    def test_a_stage_with_no_eligible_rows_is_not_called(self):
        df = pd.DataFrame({"CAS": [""], "name": ["ethanol"]})
        by_cas = FakeSearch({})
        by_name = FakeSearch({"ethanol": ETHANOL})

        resolve_cascade(df, [("cas", by_cas, "CAS"), ("name", by_name, "name")])

        assert by_cas.calls == []

    def test_row_order_index_and_original_columns_survive(self):
        df = pd.DataFrame(
            {"CAS": ["64-17-5", "71-43-2"], "bp": [78.4, 80.1]}, index=["a", "b"]
        )
        by_cas = FakeSearch({"64-17-5": ETHANOL, "71-43-2": BENZENE})

        out = resolve_cascade(df, [("cas", by_cas, "CAS")])

        assert out.index.tolist() == ["a", "b"]
        assert out["bp"].tolist() == [78.4, 80.1]
        assert out["provesid_name"].tolist() == ["ethanol", "benzene"]

    def test_default_accept_requires_an_inchikey(self):
        df = pd.DataFrame({"CAS": ["64-17-5", "no-structure"]})
        by_cas = FakeSearch({
            "64-17-5": ETHANOL,
            "no-structure": {"name": "mystery"},     # a name but no InChIKey
        })

        out = resolve_cascade(df, [("cas", by_cas, "CAS")])

        assert out["provesid_resolved_by"].tolist() == ["cas", "none"]
        assert out["provesid_validated_by"].tolist() == ["inchikey", ""]

    def test_accept_rejecting_a_hit_keeps_the_row_pending(self):
        """A contradictory identifier must not end the cascade.

        The CAS points at formaldehyde while the row's own structure is benzene,
        so the MW check rejects the CAS hit and the name stage resolves the row.
        """
        df = pd.DataFrame({
            "CAS": ["50-00-0"],          # formaldehyde
            "name": ["benzene"],
            "SMILES": ["c1ccccc1"],      # ... but the structure is benzene
        })
        by_cas = FakeSearch({"50-00-0": FORMALDEHYDE})
        by_name = FakeSearch({"benzene": BENZENE})

        out = resolve_cascade(
            df,
            [("cas", by_cas, "CAS"), ("name", by_name, "name")],
            accept=mw_within(0.5, reference_column="SMILES"),
        )

        assert out["provesid_resolved_by"].tolist() == ["name"]
        assert out["provesid_name"].tolist() == ["benzene"]

    def test_accept_may_return_a_bool(self):
        df = pd.DataFrame({"CAS": ["64-17-5", "71-43-2"]})
        by_cas = FakeSearch({"64-17-5": ETHANOL, "71-43-2": BENZENE})

        out = resolve_cascade(
            df,
            [("cas", by_cas, "CAS")],
            accept=lambda hit, row: hit.get("name") == "ethanol",
        )

        assert out["provesid_resolved_by"].tolist() == ["cas", "none"]
        assert out["provesid_validated_by"].tolist() == ["accept", ""]

    def test_accept_check_names_are_joined_into_validated_by(self):
        df = pd.DataFrame({"CAS": ["64-17-5"], "SMILES": ["CCO"], "name": ["ethanol"]})
        by_cas = FakeSearch({"64-17-5": ETHANOL})

        out = resolve_cascade(
            df,
            [("cas", by_cas, "CAS")],
            accept=mw_within(0.5, reference_column="SMILES", name_column="name"),
        )

        assert out["provesid_validated_by"].tolist() == ["mw+smiles+name"]

    def test_rdkit_fallback_derives_identifiers_from_the_given_structure(self):
        df = pd.DataFrame({"CAS": ["nonsense"], "SMILES": ["CCO"]})
        by_cas = FakeSearch({})

        out = resolve_cascade(
            df, [("cas", by_cas, "CAS")], fallback_column="SMILES"
        )

        row = out.iloc[0]
        assert row["provesid_resolved_by"] == "rdkit"
        assert row["provesid_InChIKey"] == ETHANOL["InChIKey"]
        assert row["provesid_source"] == "RDKit"

    def test_rdkit_fallback_gives_up_on_an_unparseable_structure(self):
        df = pd.DataFrame({"CAS": ["nonsense"], "SMILES": ["not-a-smiles"]})
        by_cas = FakeSearch({})

        out = resolve_cascade(
            df, [("cas", by_cas, "CAS")], fallback_column="SMILES"
        )

        assert out["provesid_resolved_by"].tolist() == ["none"]
        assert pd.isna(out["provesid_InChIKey"].iloc[0])

    def test_without_a_fallback_column_unresolved_rows_are_empty(self):
        df = pd.DataFrame({"CAS": ["nonsense"], "SMILES": ["CCO"]})
        by_cas = FakeSearch({})

        out = resolve_cascade(df, [("cas", by_cas, "CAS")])

        assert out["provesid_resolved_by"].tolist() == ["none"]
        assert pd.isna(out["provesid_InChIKey"].iloc[0])

    def test_every_output_column_is_present_with_the_prefix(self):
        df = pd.DataFrame({"CAS": ["64-17-5"]})
        by_cas = FakeSearch({"64-17-5": ETHANOL})

        out = resolve_cascade(df, [("cas", by_cas, "CAS")], prefix="p_")

        for column in OUTPUT_COLUMNS:
            assert f"p_{column}" in out.columns
        assert "p_resolved_by" in out.columns
        assert "p_validated_by" in out.columns

    def test_empty_stages_raises_valueerror(self):
        with pytest.raises(ValueError, match="at least one"):
            resolve_cascade(pd.DataFrame({"CAS": ["64-17-5"]}), [])

    def test_stage_naming_a_missing_column_raises_keyerror(self):
        with pytest.raises(KeyError, match="absent"):
            resolve_cascade(
                pd.DataFrame({"CAS": ["64-17-5"]}),
                [("cas", FakeSearch({}), "absent")],
            )

    def test_empty_dataframe_is_handled(self):
        out = resolve_cascade(
            pd.DataFrame({"CAS": []}), [("cas", FakeSearch({}), "CAS")]
        )
        assert len(out) == 0
