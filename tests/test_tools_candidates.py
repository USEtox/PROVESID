"""
Tests for the candidate adapters in :mod:`provesid.tools`, and for how a
hit's CAS number is chosen from them.

Everything here is offline: ZeroPM is a stub answering the two calls the
adapter makes, with the rows the real v0.0.4 database holds for ethanol, and
``Search`` is given stub clients with rows copied from the real databases.
"""

import pandas as pd
import pytest

from provesid.search import Search
from provesid.tools import (
    candidate_from_cactus,
    candidate_from_chebi_row,
    candidate_from_comptox_row,
    candidate_from_pubchem_row,
    candidate_from_zeropm_smiles,
    extract_cas_values,
    make_candidate,
    pick_casrn,
    sort_cas_by_number,
)

ETHANOL = ("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3", "LFQSCWFLJHTTHZ-UHFFFAOYSA-N")
ETHANOL_13C = ("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3/i2+1", "LFQSCWFLJHTTHZ-VQEHIDDOSA-N")


class _ZeroPMStub:
    """What ZeroPM answers for "CCO", with 13C-ethanol's CAS first."""

    TABLES = {
        "14742-23-5": [ETHANOL_13C, ETHANOL],
        "64-17-5": [ETHANOL],
    }

    def get_cas_from_smiles(self, smiles):
        return list(self.TABLES)

    def get_id_table_from_cas(self, cas):
        return pd.DataFrame(
            [{"cas": cas, "rank": rank, "inchi": inchi, "inchikey": key, "synonyms": ""}
             for rank, (inchi, key) in enumerate(self.TABLES[cas], start=1)]
        )


@pytest.mark.unit
def test_zeropm_structure_candidate_takes_the_query_structure_not_a_relative():
    candidate = candidate_from_zeropm_smiles("CCO", _ZeroPMStub())
    assert candidate["InChIKey"] == ETHANOL[1]
    assert candidate["InChI"] == ETHANOL[0]
    assert candidate["CAS_candidates"] == ["14742-23-5", "64-17-5"]


# ─────────────────────────────────────────────────────────────────────────────
# CAS order (plan §23.5, §26.5)
# ─────────────────────────────────────────────────────────────────────────────

# Rows as the installed databases hold them, cut to the fields that matter.
ATRAZINE_SMILES = "CCNC1=NC(Cl)=NC(NC(C)C)=N1"
COMPTOX_ATRAZINE = {
    "DTXSID": "DTXSID9020112", "PREFERRED_NAME": "Atrazine", "CASRN": "1912-24-9",
    "SMILES": ATRAZINE_SMILES, "INCHIKEY": "MXWJVTOOROXGIU-UHFFFAOYSA-N",
    "identifiers": ["1912-24-9", "11121-31-6", "12040-45-8", "39400-72-1", "Atrazine"],
}
CAMPHOR_R_SMILES = "CC1(C)[C@@H]2CC[C@@]1(C)C(=O)C2"
CAMPHOR_R_KEY = "DSSYKIVIOFKYAU-XCBNKYQSSA-N"
CHEBI_R_CAMPHOR = {
    "ChEBI ID": "CHEBI:15396", "ChEBI NAME": "(R)-camphor", "SMILES": CAMPHOR_R_SMILES,
    "INCHIKEY": CAMPHOR_R_KEY, "CAS Registry Numbers": "464-49-3;76-22-2",
}
COMPTOX_R_CAMPHOR = {
    "DTXSID": "DTXSID1025512", "PREFERRED_NAME": "(+)-Camphor", "CASRN": "464-49-3",
    "SMILES": CAMPHOR_R_SMILES, "INCHIKEY": CAMPHOR_R_KEY,
    "identifiers": ["464-49-3", "(+)-Camphor"],
}


class TestExtractionKeepsTheSourceOrder:
    def test_first_occurrence_order_not_text_order(self):
        assert extract_cas_values(["50-78-2", "11126-35-5 | 50-78-2"]) == ["50-78-2", "11126-35-5"]

    def test_a_set_is_read_sorted(self):
        assert extract_cas_values({"64-17-5", "121182-78-3"}) == ["121182-78-3", "64-17-5"]

    def test_make_candidate_keeps_the_order(self):
        cand = make_candidate("CompTox", cas_candidates=["50-78-2", "11126-35-5", "50-78-2"])
        assert cand["CAS_candidates"] == ["50-78-2", "11126-35-5"]

    def test_comptox_casrn_column_comes_first(self):
        assert candidate_from_comptox_row(COMPTOX_ATRAZINE)["CAS_candidates"][0] == "1912-24-9"

    def test_pubchem_keeps_its_order(self):
        row = {"cmpdname": "caffeine", "cas_numbers": ["58-08-2", "95789-13-2", "114303-55-8"]}
        assert candidate_from_pubchem_row(row)["CAS_candidates"] == ["58-08-2", "95789-13-2", "114303-55-8"]


class TestUnrankedSourcesAreOrderedByNumber:
    def test_sort_cas_by_number(self):
        assert sort_cas_by_number(["11126-35-5", "2349-94-2", "50-78-2"]) == [
            "50-78-2", "2349-94-2", "11126-35-5"]

    def test_chebi_reads_only_its_cas_field(self):
        # ChEBI's InChI for this thioester contains "(12)14-10-6-8".
        row = {"ChEBI NAME": "S-(2,5-dimethyl-3-furyl) 3-methylbutanethioate",
               "INCHI": "InChI=1S/C11H16O2S/c1-7(2)5-11(12)14-10-6-8(3)13-9(10)4/h6-7H,5H2,1-4H3",
               "CAS Registry Numbers": "55764-28-8"}
        assert candidate_from_chebi_row(row)["CAS_candidates"] == ["55764-28-8"]

    def test_chebi_list_is_ordered_by_number(self):
        row = {"CAS Registry Numbers": "11126-35-5;50-78-2"}
        assert candidate_from_chebi_row(row)["CAS_candidates"] == ["50-78-2", "11126-35-5"]

    def test_cactus_list_is_ordered_by_number(self):
        # CACTUS's names for ethanol, 2026-09-24: 121182-78-3 comes first.
        cand = candidate_from_cactus("CCO", ["121182-78-3", "64-17-5", "Ethanol"])
        assert cand["CAS_candidates"] == ["64-17-5", "121182-78-3"]


class TestPickCasrn:
    def test_a_ranked_source_wins_over_an_unranked_one_listed_earlier(self):
        chebi = candidate_from_chebi_row(CHEBI_R_CAMPHOR)
        comptox = candidate_from_comptox_row(COMPTOX_R_CAMPHOR)
        assert chebi["CAS_candidates"][0] == "76-22-2"
        assert pick_casrn([chebi, comptox]) == "464-49-3"

    def test_an_unranked_source_is_used_when_it_is_the_only_one(self):
        cactus = candidate_from_cactus("CCO", ["121182-78-3", "64-17-5"])
        assert pick_casrn([None, cactus]) == "64-17-5"

    def test_an_unranked_source_with_one_number_keeps_its_place(self):
        chebi = make_candidate("ChEBI", cas_candidates=["1345-25-1"])
        pubchem = make_candidate("PubChemID", cas_candidates=["17125-56-3", "1345-25-1"])
        assert pick_casrn([chebi, pubchem]) == "1345-25-1"

    def test_a_ranked_source_without_a_cas_is_passed_over(self):
        chebi = make_candidate("ChEBI", cas_candidates=["50-78-2"])
        pubchem = make_candidate("PubChemID", cas_candidates=[])
        assert pick_casrn([chebi, pubchem]) == "50-78-2"

    def test_no_cas_at_all(self):
        assert pick_casrn([make_candidate("ChEBI"), None]) is None


class _Client:
    """Answers the named methods with fixed values, and anything else with None."""

    def __init__(self, **answers):
        self.answers = answers

    def __getattr__(self, method):
        if method.startswith("__"):
            raise AttributeError(method)
        return lambda *args, **kwargs: self.answers.get(method)


class TestSearchReportsTheCurrentNumber:
    def test_a_search_by_a_retired_number_reports_the_current_one(self):
        comptox = _Client(get_by_alternate_casrn=COMPTOX_ATRAZINE)
        with Search("cas", show_progress=False, sources="comptox", comptox=comptox) as s:
            row = s.search("39400-72-1").iloc[0]
        assert row["query"] == "39400-72-1"
        assert row["CASRN"] == "1912-24-9"

    def test_a_search_by_a_stereoisomer_number_keeps_it(self):
        chebi = _Client(search_by_cas=[CHEBI_R_CAMPHOR])
        comptox = _Client(get_by_casrn=COMPTOX_R_CAMPHOR)
        with Search("cas", show_progress=False, sources=["chebi", "comptox"],
                    chebi=chebi, comptox=comptox) as s:
            row = s.search("464-49-3").iloc[0]
        assert row["source"] == "ChEBI"
        assert row["CASRN"] == "464-49-3"

    def test_a_miss_has_no_casrn(self):
        with Search("cas", show_progress=False, sources="comptox", comptox=_Client()) as s:
            row = s.search("0-00-0").iloc[0]
        assert row["query"] == "0-00-0"
        assert row["CASRN"] is None
