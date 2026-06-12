"""Tests for the multi-hit / tunable Search features.

Covers the candidate-pool refactor, structure clustering, ranked multi-hit
output, query-aware ranking (the wrong-name-hit fix), tunable attributes, the
DataFrame broadcast for multi-hit rows, and OPSIN anchoring (skipped when no
Java/py2opsin runtime is available).

Run with::

    uv run pytest tests/test_search_multihit.py -v
"""

from __future__ import annotations

import pandas as pd
import pytest

from provesid.search import (
    OUTPUT_COLUMNS,
    Search,
    _candidate_cluster_keys,
    _cluster_candidates,
)
from provesid.tools import _make_candidate


# ─────────────────────────────────────────────────────────────────────────────
# Helpers / stubs
# ─────────────────────────────────────────────────────────────────────────────

def _tag(cand, source_key, rank=0, method="exact_name", score=1.0):
    cand["_source_key"] = source_key
    cand["_origin_rank"] = rank
    cand["_match_method"] = method
    cand["query_match_score"] = score
    return cand


class _PubChemNameStub:
    """Returns a fixed list of rows for any name search; no InChIKey hits."""

    def __init__(self, rows):
        self._rows = rows

    def search_by_name(self, name, exact=False, limit=10):
        return self._rows

    def get_by_inchikey(self, ik):
        return None


_ETHANOL = {
    "cmpdname": "Wrong Hit",
    "iupacname": "ethanol",
    "mf": "C2H6O",
    "smiles": "CCO",
    "inchi": None,
    "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
    "mw": 46.07,
    "synonyms": ["ethanol"],
    "cas_numbers": ["64-17-5"],
}

_XYLENE = {
    "cmpdname": "xylene",
    "iupacname": "1,2-dimethylbenzene",
    "mf": "C8H10",
    "smiles": "Cc1ccccc1C",
    "inchi": None,
    "inchikey": "CTQNGGLPUBDAKN-UHFFFAOYSA-N",
    "mw": 106.16,
    "synonyms": ["o-xylene", "xylene"],
    "cas_numbers": ["95-47-6"],
}


def _name_search(rows, **kwargs):
    return Search(
        "name",
        show_progress=False,
        pubchem=_PubChemNameStub(rows),
        chebi=None,
        comptox=None,
        zeropm=None,
        chembl=None,
        **kwargs,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Clustering
# ─────────────────────────────────────────────────────────────────────────────

class TestClustering:
    def test_same_inchikey_merges(self):
        a = _tag(_make_candidate("ChEBI", name="aspirin", inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N"), "chebi")
        b = _tag(_make_candidate("CompTox", name="Aspirin", inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N"), "comptox")
        c = _tag(_make_candidate("PubChemID", name="caffeine", inchikey="RYYVLZVUVIJVGH-UHFFFAOYSA-N"), "pubchem")
        clusters = _cluster_candidates([a, b, c], by_skeleton=True)
        assert len(clusters) == 2
        assert sorted(len(cl["members"]) for cl in clusters) == [1, 2]

    def test_skeleton_merge_toggle(self):
        d = _tag(_make_candidate("ChEBI", name="x", inchikey="ABCDEFGHIJKLMN-AAAAAAAAAA-N"), "chebi")
        e = _tag(_make_candidate("CompTox", name="y", inchikey="ABCDEFGHIJKLMN-BBBBBBBBBB-N"), "comptox")
        assert len(_cluster_candidates([d, e], by_skeleton=True)) == 1
        assert len(_cluster_candidates([d, e], by_skeleton=False)) == 2

    def test_name_only_fallback_singletons(self):
        a = _tag(_make_candidate("ChEBI", name="foo"), "chebi")
        b = _tag(_make_candidate("CompTox", name="bar"), "comptox")
        assert len(_cluster_candidates([a, b], by_skeleton=True)) == 2

    def test_cluster_keys_priority(self):
        ik = _make_candidate("X", inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N", smiles="CCO")
        keys = _candidate_cluster_keys(ik, by_skeleton=True)
        # InChIKey present -> ik + skel keys, no smiles key
        assert ("ik", "BSYNRYMUTXBXSQ-UHFFFAOYSA-N") in keys
        assert ("skel", "BSYNRYMUTXBXSQ") in keys
        assert not any(k[0] == "smi" for k in keys)


# ─────────────────────────────────────────────────────────────────────────────
# n_hits validation
# ─────────────────────────────────────────────────────────────────────────────

class TestNHitsValidation:
    @pytest.mark.parametrize("value", [1, 3, "all", "ALL"])
    def test_valid(self, value):
        Search("name", n_hits=value, show_progress=False)

    @pytest.mark.parametrize("value", [0, -1, "two", 1.5, True])
    def test_invalid(self, value):
        with pytest.raises(ValueError):
            Search("name", n_hits=value, show_progress=False)

    def test_invalid_fuzzy_scorer(self):
        with pytest.raises(ValueError):
            Search("name", fuzzy_scorer="does_not_exist", show_progress=False)


# ─────────────────────────────────────────────────────────────────────────────
# Multi-hit output + ranking
# ─────────────────────────────────────────────────────────────────────────────

class TestMultiHit:
    def test_default_single_row(self):
        s = _name_search([_ETHANOL, _XYLENE])
        df = s.search("xylene")
        assert len(df) == 1
        assert df.iloc[0]["hit_rank"] == 0

    def test_query_aware_ranking_fixes_wrong_hit(self):
        # PubChem returns ethanol as hit #0, xylene as hit #1.  The query
        # "xylene" must still resolve to xylene at rank 0.
        s = _name_search([_ETHANOL, _XYLENE])
        df = s.search("xylene")
        assert df.iloc[0]["name"] == "xylene"

    def test_n_hits_all_returns_distinct_compounds(self):
        s = _name_search([_ETHANOL, _XYLENE])
        df = s.search("xylene", n_hits="all")
        assert len(df) == 2
        assert list(df["hit_rank"]) == [0, 1]
        # confidence is non-increasing with rank
        assert df.iloc[0]["confidence"] >= df.iloc[1]["confidence"]

    def test_n_hits_int_truncates(self):
        s = _name_search([_ETHANOL, _XYLENE])
        assert len(s.search("xylene", n_hits=1)) == 1
        assert len(s.search("xylene", n_hits=5)) == 2  # only 2 distinct compounds

    def test_min_confidence_filters(self):
        s = _name_search([_ETHANOL, _XYLENE])
        # high floor keeps only the strong (exact-name) hit
        df = s.search("xylene", n_hits="all", min_confidence=0.6)
        assert len(df) == 1
        assert df.iloc[0]["name"] == "xylene"

    def test_new_columns_present(self):
        s = _name_search([_XYLENE])
        df = s.search("xylene")
        for col in ("hit_rank", "n_source_support", "opsin_smiles"):
            assert col in df.columns
        assert col in OUTPUT_COLUMNS

    def test_dataframe_broadcast(self):
        s = _name_search([_ETHANOL, _XYLENE])
        inp = pd.DataFrame({"q": ["xylene"], "batch": ["b1"]})
        out = s.search(inp, column="q", n_hits="all")
        assert len(out) == 2
        assert set(out["batch"]) == {"b1"}

    def test_return_alternatives(self):
        s = _name_search([_ETHANOL, _XYLENE], return_alternatives=True)
        df = s.search("xylene")  # n_hits defaults to 1
        assert len(df) == 1
        alts = df.iloc[0]["alternatives"]
        assert isinstance(alts, list) and len(alts) >= 1
        assert alts[0]["name"] == "Wrong Hit"

    def test_no_match_still_one_row(self):
        s = _name_search([])
        df = s.search("nonexistent-compound")
        assert len(df) == 1
        assert df.iloc[0]["confidence"] == 0.0
        assert df.iloc[0]["match_method"] == "exact_name"


# ─────────────────────────────────────────────────────────────────────────────
# OPSIN anchoring (requires Java + py2opsin)
# ─────────────────────────────────────────────────────────────────────────────

@pytest.mark.opsin
class TestOpsinAnchor:
    def test_opsin_disabled_by_default(self):
        s = _name_search([_XYLENE])
        assert s.use_opsin is False
        # OPSIN column present but empty when disabled
        df = s.search("xylene")
        assert df.iloc[0]["opsin_smiles"] is None

    def test_opsin_anchor_when_available(self):
        s = _name_search([_XYLENE], use_opsin=True)
        anchor = s._opsin_anchor("1,2-dimethylbenzene")
        if anchor is None:
            pytest.skip("PYOPSIN/Java unavailable in this environment")
        assert anchor["smiles"]
        df = s.search("1,2-dimethylbenzene")
        assert df.iloc[0]["opsin_smiles"] is not None
