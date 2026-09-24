"""Tests for provesid.sources: the lookup table Search queries, and its driver.

No database is opened.  The table is exercised with small stubs, and
``Search._collect`` with stubs that record what they were asked.
"""

import logging

import pandas as pd
import pytest

from provesid.search import Search
from provesid.sources import (
    LOOKUPS,
    ONLINE_SOURCE_KEYS,
    SOURCE_DISPLAY,
    SOURCE_KEYS,
    Query,
    rank_rows_by_completeness,
)

_ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
_ASPIRIN_IK = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"


class _PubChem:
    """Answers every PubChemID call with aspirin, and records the calls."""

    def __init__(self, row=None):
        self.row = row if row is not None else {
            "cid": 2244, "smiles": _ASPIRIN_SMILES, "inchikey": _ASPIRIN_IK,
            "iupac_name": "2-acetyloxybenzoic acid", "synonyms": ["aspirin"],
            "cas_numbers": ["50-78-2"],
        }
        self.calls = []

    def __getattr__(self, method):
        def answer(*args, **kwargs):
            self.calls.append((method, args, kwargs))
            return [self.row] if method.startswith("search_") else self.row
        return answer


class _Broken:
    """A source whose every call raises."""

    def __getattr__(self, method):
        def fail(*args, **kwargs):
            raise RuntimeError("database is corrupt")
        return fail


def _search(**clients):
    """A Search with exactly these clients and no others."""
    return Search(
        "cas",
        show_progress=False,
        sources=list(clients),
        **clients,
    )


# ─────────────────────────────────────────────────────────────────────────────
# The table
# ─────────────────────────────────────────────────────────────────────────────

class TestTable:
    def test_every_source_named_is_a_known_source(self):
        for kind, row in LOOKUPS.items():
            assert set(row) <= set(SOURCE_KEYS + ONLINE_SOURCE_KEYS), kind

    def test_every_source_has_a_display_name(self):
        assert set(SOURCE_DISPLAY) == set(SOURCE_KEYS + ONLINE_SOURCE_KEYS)

    def test_online_and_offline_keys_are_disjoint(self):
        assert not set(SOURCE_KEYS) & set(ONLINE_SOURCE_KEYS)

    def test_every_resolver_kind_is_in_the_table(self):
        assert set(LOOKUPS) == {
            "cas", "inchikey", "inchikey_skeleton", "inchi", "smiles",
            "dtxsid", "name", "fuzzy_name", "formula",
        }

    def test_chembl_has_no_cas_row(self):
        """ChEMBL records no CAS numbers; Search reaches it by SMILES."""
        assert "chembl" not in LOOKUPS["cas"]

    def test_a_miss_is_an_empty_list(self):
        assert LOOKUPS["cas"]["pubchem"](_PubChem(row={}), Query("0-00-0")) == []

    def test_a_hit_is_an_adapted_candidate(self):
        (cand,) = LOOKUPS["cas"]["pubchem"](_PubChem(), Query("50-78-2"))
        assert cand["source"] == "PubChemID"
        assert cand["InChIKey"] == _ASPIRIN_IK

    def test_name_lookups_pass_k_as_the_limit(self):
        client = _PubChem()
        LOOKUPS["name"]["pubchem"](client, Query("aspirin", k=3))
        assert client.calls == [("search_by_name", ("aspirin",), {"exact": True, "limit": 3})]

    def test_zeropm_candidate_is_named_after_the_label(self):
        table = pd.DataFrame([{"inchi": None, "inchikey": _ASPIRIN_IK, "cas": "50-78-2", "name": "aspirin"}])

        class ZeroPM:
            def get_id_table_from_inchikey(self, ik):
                return table

        (cand,) = LOOKUPS["inchikey"]["zeropm"](
            ZeroPM(), Query(_ASPIRIN_IK, label="DTXSID5020108")
        )
        assert cand["name"] == "DTXSID5020108"

    def test_zeropm_fuzzy_keeps_its_own_score(self):
        table = pd.DataFrame([{
            "inchi": None, "inchikey": _ASPIRIN_IK, "cas": "50-78-2",
            "name": "aspirin", "matched_name": "aspirin", "match_score": 92.0,
        }])

        class ZeroPM:
            def get_id_table_from_similar_name(self, name, number_of_results, score_cutoff):
                assert (number_of_results, score_cutoff) == (4, 85.0)
                return table

        (cand,) = LOOKUPS["fuzzy_name"]["zeropm"](
            ZeroPM(), Query("asprin", k=4, fuzzy_cutoff=85.0)
        )
        assert cand["name"] == "aspirin"
        assert cand["query_match_score"] == pytest.approx(0.92)


class TestQuery:
    def test_label_defaults_to_value(self):
        assert Query("50-78-2").label == "50-78-2"

    def test_label_can_differ(self):
        assert Query(_ASPIRIN_IK, label="DTXSID5020108").label == "DTXSID5020108"


class TestRankRowsByCompleteness:
    def test_most_complete_first(self):
        rows = [{"a": 1, "b": None}, {"a": 1, "b": 2}]
        assert rank_rows_by_completeness(rows) == [{"a": 1, "b": 2}, {"a": 1, "b": None}]

    def test_stable_for_ties(self):
        rows = [{"a": 1}, {"a": 2}]
        assert rank_rows_by_completeness(rows) == rows

    def test_none_is_empty(self):
        assert rank_rows_by_completeness(None) == []


# ─────────────────────────────────────────────────────────────────────────────
# Search._collect, the one caller
# ─────────────────────────────────────────────────────────────────────────────

class TestCollect:
    def test_asks_only_sources_with_a_client(self):
        s = _search(pubchem=_PubChem())
        assert set(s._collect("cas", "50-78-2")) == {"pubchem"}

    def test_skips_sources_with_no_row_for_the_kind(self):
        s = _search(pubchem=_PubChem(), chembl=_PubChem())
        assert set(s._collect("cas", "50-78-2")) == {"pubchem"}

    def test_sources_restricts_the_question(self):
        pubchem, chembl = _PubChem(), _PubChem()
        s = _search(pubchem=pubchem, chembl=chembl)
        s._collect("smiles", _ASPIRIN_SMILES, sources=["chembl"])
        assert pubchem.calls == []
        assert [c[0] for c in chembl.calls] == ["search_by_smiles"]

    def test_a_failing_source_is_logged_and_the_rest_still_answer(self, caplog):
        s = _search(pubchem=_PubChem(), comptox=_Broken())
        with caplog.at_level(logging.WARNING, logger="provesid.search"):
            hits = s._collect("cas", "50-78-2")
        assert set(hits) == {"pubchem"}
        assert any(
            "CompTox cas lookup failed for '50-78-2'" in r.getMessage()
            for r in caplog.records
        )

    def test_pool_is_in_source_order_and_rank_order(self):
        s = _search(chebi=_PubChem(), pubchem=_PubChem())
        hits = {
            "pubchem": [{"n": 1}, {"n": 2}],
            "chebi": [{"n": 0}],
        }
        pool = s._pool(hits, "exact_cas")
        assert [(c["_source_key"], c["_origin_rank"]) for c in pool] == [
            ("chebi", 0), ("pubchem", 0), ("pubchem", 1),
        ]

    def test_pool_scores_with_a_function(self):
        s = _search(pubchem=_PubChem())
        pool = s._pool({"pubchem": [{"n": 3}]}, "formula", lambda c: c["n"] / 10)
        assert pool[0]["query_match_score"] == pytest.approx(0.3)

    def test_cas_reaches_chembl_through_the_smiles_found(self):
        class ChEMBL:
            asked = []

            def search_by_smiles(self, smiles):
                self.asked.append(smiles)
                return {
                    "molregno": 1, "chembl_id": "CHEMBL25", "pref_name": "ASPIRIN",
                    "canonical_smiles": _ASPIRIN_SMILES,
                    "standard_inchi_key": _ASPIRIN_IK,
                }

            def get_properties(self, molregno):
                return {}

        chembl = ChEMBL()
        s = _search(pubchem=_PubChem(), chembl=chembl)
        s._ensure_clients()
        _, pool, _ = s._resolve_cas("50-78-2")
        assert {c["_source_key"] for c in pool} == {"pubchem", "chembl"}
        assert chembl.asked == [_ASPIRIN_SMILES]
