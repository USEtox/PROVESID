"""Structure quality control in Search (plan 2026-10-13).

The unit tests open no database: candidates are built directly, and Search
runs on stub clients.  The regression set at the end resolves the CAS numbers
of the Sanofi structure audit against the installed offline databases, and
is skipped when ChEBI, CompTox or PubChem is missing.

Run with::

    uv run pytest tests/test_structure_qc.py -v
"""

from __future__ import annotations

import logging

import pandas as pd
import pytest
from rdkit import Chem

from provesid import datasets
from provesid.chebi_sdf import _restore_molfile_header
from provesid.search import Search, _cluster_candidates
from provesid.tools import (
    apply_candidate_to_result,
    candidate_from_chebi_row,
    candidate_similarity,
    compute_consensus,
    formula_element_counts,
    identity_consistent,
    make_candidate,
    repair_structure,
)

# ChEBI's SMILES for these drops the ring N-H; the InChIKeys are ChEBI's.
FLUOROURACIL_BROKEN = "O=c1ncc(F)c(=O)n1"
FLUOROURACIL = "O=c1[nH]cc(F)c(=O)[nH]1"
FLUOROURACIL_IK = "GHASVSINZRGABV-UHFFFAOYSA-N"
VALSARTAN_BROKEN = "CCCCC(=O)N(Cc1ccc(-c2ccccc2-c2nnnn2)cc1)[C@H](C(=O)O)C(C)C"
VALSARTAN = "CCCCC(=O)N(Cc1ccc(-c2ccccc2-c2nnn[nH]2)cc1)[C@H](C(=O)O)C(C)C"
VALSARTAN_IK = "ACWBQPMHZXGDFX-QFIPXVFZSA-N"
# Parses, but lacks two hydrogens: C5H3N5O, not guanine's C5H5N5O.
GUANINE_WRONG = "Nc1nc(=O)c2ncnc2n1"
GUANINE = "Nc1nc2[nH]cnc2c(=O)[nH]1"
GUANINE_IK = "UYTPUPDQBNUYGX-UHFFFAOYSA-N"
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
ASPIRIN_IK = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"


def _molfile(smiles: str) -> str:
    return Chem.MolToMolBlock(Chem.MolFromSmiles(smiles))


# ─────────────────────────────────────────────────────────────────────────────
# Defect A: unparsable source SMILES
# ─────────────────────────────────────────────────────────────────────────────


class TestRepairStructure:
    def test_a_parsable_smiles_is_kept(self):
        assert repair_structure(ASPIRIN, inchikey=ASPIRIN_IK) == (ASPIRIN, "smiles", "valid")

    def test_the_molfile_rebuilds_an_unparsable_smiles(self):
        smiles, origin, status = repair_structure(
            VALSARTAN_BROKEN, molfile=_molfile(VALSARTAN), inchikey=VALSARTAN_IK
        )
        assert (origin, status) == ("molfile", "repaired:molfile")
        assert Chem.MolToInchiKey(Chem.MolFromSmiles(smiles)) == VALSARTAN_IK

    def test_the_inchi_rebuilds_when_there_is_no_molfile(self):
        inchi = Chem.MolToInchi(Chem.MolFromSmiles(FLUOROURACIL))
        _, origin, status = repair_structure(FLUOROURACIL_BROKEN, inchi=inchi)
        assert (origin, status) == ("inchi", "repaired:inchi")

    def test_nh_repair_needs_the_stated_key(self):
        smiles, origin, _ = repair_structure(FLUOROURACIL_BROKEN, inchikey=FLUOROURACIL_IK)
        assert origin == "nh_repair"
        assert Chem.MolToInchiKey(Chem.MolFromSmiles(smiles)) == FLUOROURACIL_IK

    def test_a_rebuild_with_another_key_is_refused(self):
        assert repair_structure(
            VALSARTAN_BROKEN, molfile=_molfile(VALSARTAN), inchikey=ASPIRIN_IK
        ) == (None, None, "invalid")

    def test_nothing_to_check_against_is_invalid(self):
        assert repair_structure(VALSARTAN_BROKEN, molfile=_molfile(VALSARTAN)) == (
            None, None, "invalid",
        )

    def test_a_parsable_smiles_of_another_compound_yields_to_the_molfile(self):
        smiles, origin, status = repair_structure(
            GUANINE_WRONG, molfile=_molfile(GUANINE), inchikey=GUANINE_IK
        )
        assert (origin, status) == ("molfile", "repaired:molfile")
        assert Chem.MolToInchiKey(Chem.MolFromSmiles(smiles)) == GUANINE_IK

    def test_without_a_matching_molfile_a_parsable_smiles_stays(self):
        assert repair_structure(GUANINE_WRONG, inchikey=GUANINE_IK) == (
            GUANINE_WRONG, "smiles", "valid",
        )

    def test_no_smiles_is_none(self):
        assert repair_structure(None, inchikey=ASPIRIN_IK) == (None, None, "none")


class TestCandidates:
    def test_an_invalid_smiles_never_reaches_the_candidate(self):
        cand = make_candidate("ChEBI", smiles=VALSARTAN_BROKEN, inchikey="AAAAAAAAAAAAAA-UHFFFAOYSA-N")
        assert cand["SMILES"] is None and cand["canonical_smiles"] is None
        assert cand["structure_status"] == "invalid"
        assert cand["source_smiles_raw"] == VALSARTAN_BROKEN
        assert cand["InChIKey"] == "AAAAAAAAAAAAAA-UHFFFAOYSA-N"

    def test_a_chebi_row_passes_its_molfile(self):
        # The index drops a blank title line; the reader puts it back.
        stripped = _molfile(VALSARTAN).split("\n", 1)[1]
        row = {
            "ChEBI NAME": "valsartan", "SMILES": VALSARTAN_BROKEN,
            "INCHIKEY": VALSARTAN_IK, "molfile": _restore_molfile_header(stripped),
        }
        cand = candidate_from_chebi_row(row)
        assert cand["structure_status"] == "repaired:molfile"
        assert cand["canonical_smiles"] is not None

    def test_restore_molfile_header_leaves_a_whole_molfile_alone(self):
        whole = _molfile(ASPIRIN)
        assert _restore_molfile_header(whole) == whole
        assert Chem.MolFromMolBlock(_restore_molfile_header(whole.split("\n", 1)[1])) is not None

    def test_an_invalid_smiles_does_not_hold_the_result_field(self):
        result = {"SMILES": VALSARTAN_BROKEN, "source": None}
        apply_candidate_to_result(result, make_candidate("CompTox", smiles=VALSARTAN))
        assert result["SMILES"] == VALSARTAN
        assert result["source"] == "CompTox"

    def test_tautomers_with_one_inchikey_agree(self):
        chebi = make_candidate("ChEBI", smiles="Oc1ncnc2[nH]ncc12", inchikey="OFCNXPDARWKPPY-UHFFFAOYSA-N")
        pubchem = make_candidate("PubChemID", smiles="O=c1[nH]cnc2[nH]ncc12",
                                 inchikey="OFCNXPDARWKPPY-UHFFFAOYSA-N")
        assert chebi["canonical_smiles"] != pubchem["canonical_smiles"]
        assert candidate_similarity(chebi, pubchem) == 1.0


# ─────────────────────────────────────────────────────────────────────────────
# Defect E: identity consistency
# ─────────────────────────────────────────────────────────────────────────────


class TestIdentity:
    @pytest.mark.parametrize("left, right", [
        ("C17H21NO.HCl", "C17H22ClNO"),
        ("C72H104O48S8.8Na", "C72H104Na8O48S8"),
        ("C10H15O3", "C10H15O3-"),
        ("C7H5NO5", "C7H5NO5-2"),
        ("C2H3NaO2.3H2O", "C2H9NaO5"),
    ])
    def test_formula_notations_compare_by_element_counts(self, left, right):
        assert formula_element_counts(left) == formula_element_counts(right)

    @pytest.mark.parametrize("formula", ["(C2H4O)n", "C4H6O4R2", "", None, "C6H12O6x"])
    def test_formulas_that_are_not_plain_do_not_compare(self, formula):
        assert formula_element_counts(formula) is None

    def test_a_key_of_another_compound_is_inconsistent(self):
        assert identity_consistent(ASPIRIN, GUANINE_IK) is False

    def test_stereo_alone_is_not_an_inconsistency(self):
        assert identity_consistent("C/C=C/C", "IAQRGUVFOMOMEM-UHFFFAOYSA-N") is True

    def test_a_salt_formula_matches_its_salt(self):
        assert identity_consistent("CNCC[C@@H](Oc1ccccc1C)c1ccccc1.Cl", formula="C17H21NO.HCl")


# ─────────────────────────────────────────────────────────────────────────────
# Search on stub clients
# ─────────────────────────────────────────────────────────────────────────────


class _Source:
    """A source client that answers each named method with fixed rows."""

    def __init__(self, **answers):
        self.answers = answers

    def __getattr__(self, method):
        answers = self.__dict__["answers"]
        if method not in answers:
            raise AttributeError(method)
        return lambda *args, **kwargs: answers[method]


def _chebi(row):
    return _Source(search_by_cas=[row])


def _comptox(row):
    return _Source(get_by_casrn=row, get_by_alternate_casrn=None)


def _pubchem(row):
    return _Source(get_by_cas=row)


def _cas_search(**kwargs):
    clients = {k: v for k, v in kwargs.items() if k in ("chebi", "comptox", "pubchem", "chembl")}
    options = {k: v for k, v in kwargs.items() if k not in clients}
    return Search("cas", show_progress=False, sources=list(clients), **clients, **options)


_CEFQUINOME_CAS = "118443-89-3"
_COMPTOX_WRONG = {
    "DTXSID": "DTXSID00000001", "PREFERRED_NAME": "Cefquinome sulfate", "CASRN": _CEFQUINOME_CAS,
    "SMILES": "CC(=O)Oc1ccccc1C(=O)O", "INCHIKEY": ASPIRIN_IK, "MOLECULAR_FORMULA": "C9H8O4",
}
_PUBCHEM_RIGHT = {
    "cid": 1, "cmpdname": "Cefquinome sulfate", "smiles": "CCO",
    "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N", "mf": "C2H6O", "cas_numbers": [_CEFQUINOME_CAS],
}


class TestSearchStructureStatus:
    def test_a_repaired_chebi_structure_fills_the_row(self):
        row = {"ChEBI NAME": "valsartan", "SMILES": VALSARTAN_BROKEN, "INCHIKEY": VALSARTAN_IK,
               "molfile": _molfile(VALSARTAN), "CAS Registry Numbers": "137862-53-4"}
        hit = _cas_search(chebi=_chebi(row)).search("137862-53-4").iloc[0]
        assert hit["structure_status"] == "repaired:molfile"
        assert hit["canonical_smiles"] is not None
        assert hit["InChIKey"] == VALSARTAN_IK
        assert hit["identity_consistent"]

    def test_an_unrepairable_structure_is_flagged_and_logged(self, caplog):
        row = {"ChEBI NAME": "valsartan", "SMILES": VALSARTAN_BROKEN,
               "INCHIKEY": "AAAAAAAAAAAAAA-UHFFFAOYSA-N", "CAS Registry Numbers": "137862-53-4"}
        with caplog.at_level(logging.WARNING, logger="provesid.search"):
            hit = _cas_search(chebi=_chebi(row)).search("137862-53-4").iloc[0]
        assert hit["structure_status"] == "invalid"
        assert hit["SMILES"] is None and hit["canonical_smiles"] is None
        assert hit["InChIKey"] == "AAAAAAAAAAAAAA-UHFFFAOYSA-N"
        assert any("cannot be read or repaired" in r.getMessage() for r in caplog.records)

    def test_a_compatible_source_supplies_the_structure_chebi_could_not(self):
        chebi = {"ChEBI NAME": "valsartan", "SMILES": VALSARTAN_BROKEN, "INCHIKEY": VALSARTAN_IK,
                 "CAS Registry Numbers": "137862-53-4"}
        comptox = {"DTXSID": "DTXSID1", "PREFERRED_NAME": "Valsartan", "CASRN": "137862-53-4",
                   "SMILES": VALSARTAN, "INCHIKEY": VALSARTAN_IK}
        # No molfile, and nh_repair is the last resort: it puts the H back.
        hit = _cas_search(chebi=_chebi(chebi), comptox=_comptox(comptox)).search("137862-53-4").iloc[0]
        assert hit["canonical_smiles"] is not None
        assert hit["structure_status"] in ("valid", "repaired:nh_repair")

    def test_the_structure_decides_inchikey_and_formula(self):
        comptox = dict(_COMPTOX_WRONG, MOLECULAR_FORMULA="C9H8O4.HCl")
        hit = _cas_search(comptox=_comptox(comptox)).search(_CEFQUINOME_CAS).iloc[0]
        assert hit["molecular_formula"] == "C9H8O4"
        assert hit["source_formula"] == "C9H8O4.HCl"
        assert hit["identity_consistent"] is not None and not hit["identity_consistent"]


class TestDerivedChembl:
    def test_chembl_echo_does_not_vote(self):
        chembl = _Source(search_by_inchikey={
            "molregno": 1, "pref_name": "ASPIRIN", "canonical_smiles": ASPIRIN,
            "standard_inchi_key": ASPIRIN_IK,
        }, get_properties={})
        hit = _cas_search(comptox=_comptox(_COMPTOX_WRONG), chembl=chembl).search(_CEFQUINOME_CAS).iloc[0]
        assert hit["n_source_support"] == 1
        assert set(hit["source_match_scores"]) == {"CompTox"}
        assert hit["source_details"]["ChEMBL"]["found"]

    def test_compute_consensus_ignores_derived_candidates(self):
        a = make_candidate("ChEBI", smiles=ASPIRIN)
        echo = make_candidate("ChEMBL", smiles=ASPIRIN, derived_from="ChEBI")
        b = make_candidate("CompTox", smiles="CCO")
        source, scores, _ = compute_consensus({"chebi": a, "chembl": echo, "comptox": b})
        assert "chembl" not in scores


class TestConflicts:
    def test_disagreeing_sources_are_flagged_and_penalised(self):
        s = _cas_search(comptox=_comptox(_COMPTOX_WRONG), pubchem=_pubchem(_PUBCHEM_RIGHT))
        hit = s.search(_CEFQUINOME_CAS).iloc[0]
        assert hit["conflict"]
        assert hit["n_clusters"] == 2
        assert hit["confidence"] == pytest.approx(0.765 * 0.5)
        assert hit["runner_up_sources"] == ["PubChemID"]

    def test_strict_refuses_a_disputed_query(self):
        s = _cas_search(comptox=_comptox(_COMPTOX_WRONG), pubchem=_pubchem(_PUBCHEM_RIGHT),
                        preset="strict", min_source_support=0)
        hit = s.search(_CEFQUINOME_CAS).iloc[0]
        assert hit["canonical_smiles"] is None
        assert hit["conflict"]

    def test_a_salt_of_the_same_compound_is_no_conflict(self):
        pubchem_salt = dict(_PUBCHEM_RIGHT, smiles="CC(=O)Oc1ccccc1C(=O)[O-].[Na+]",
                            inchikey="BIDNLKIUORFRQP-UHFFFAOYSA-M")
        s = _cas_search(comptox=_comptox(_COMPTOX_WRONG), pubchem=_pubchem(pubchem_salt))
        assert not s.search(_CEFQUINOME_CAS)["conflict"].iloc[0]

    def test_a_tautomer_of_the_same_compound_is_no_conflict(self):
        keto = dict(_COMPTOX_WRONG, SMILES="CC(=O)C", INCHIKEY="CSCPPACGZOOCGX-UHFFFAOYSA-N")
        enol = dict(_PUBCHEM_RIGHT, smiles="CC(O)=C", inchikey="NARVIWMVBMUEOG-UHFFFAOYSA-N")
        s = _cas_search(comptox=_comptox(keto), pubchem=_pubchem(enol))
        hit = s.search(_CEFQUINOME_CAS).iloc[0]
        assert hit["n_clusters"] == 2 and not hit["conflict"]

    def test_agreeing_sources_have_no_conflict(self):
        s = _cas_search(comptox=_comptox(_COMPTOX_WRONG),
                        pubchem=_pubchem(dict(_PUBCHEM_RIGHT, smiles=ASPIRIN, inchikey=ASPIRIN_IK)))
        hit = s.search(_CEFQUINOME_CAS).iloc[0]
        assert not hit["conflict"] and hit["n_clusters"] == 1


class TestDtxsidQueries:
    """A DTXSID query is cross-checked through CompTox's CAS number."""

    _DTXSID = _COMPTOX_WRONG["DTXSID"]

    def _search(self, pubchem, **options):
        comptox = _Source(get_by_dtxsid=_COMPTOX_WRONG)
        return Search("dtxsid", show_progress=False, sources=["comptox", "pubchem"],
                      comptox=comptox, pubchem=pubchem, **options)

    def test_a_source_found_by_comptox_inchikey_does_not_vote(self):
        pubchem = _Source(get_by_cas=None, get_by_inchikey=dict(_PUBCHEM_RIGHT, smiles=ASPIRIN,
                                                                inchikey=ASPIRIN_IK))
        hit = self._search(pubchem).search(self._DTXSID).iloc[0]
        assert hit["n_source_support"] == 1
        assert set(hit["source_match_scores"]) == {"CompTox"}
        assert hit["source_details"]["PubChemID"]["found"]

    def test_a_source_found_by_the_cas_number_votes(self):
        pubchem = _Source(get_by_cas=dict(_PUBCHEM_RIGHT, smiles=ASPIRIN, inchikey=ASPIRIN_IK))
        hit = self._search(pubchem).search(self._DTXSID).iloc[0]
        assert hit["n_source_support"] == 2 and not hit["conflict"]

    def test_a_disagreeing_cas_record_is_a_conflict_the_dtxsid_owner_wins(self):
        hit = self._search(_pubchem(_PUBCHEM_RIGHT)).search(self._DTXSID).iloc[0]
        assert hit["conflict"]
        assert hit["InChIKey"] == ASPIRIN_IK
        assert hit["DTXSID"] == self._DTXSID
        assert hit["runner_up_sources"] == ["PubChemID"]

    def test_strict_refuses_a_disputed_dtxsid(self):
        hit = self._search(_pubchem(_PUBCHEM_RIGHT), preset="strict",
                           min_source_support=0).search(self._DTXSID).iloc[0]
        assert hit["canonical_smiles"] is None and hit["conflict"]


# ─────────────────────────────────────────────────────────────────────────────
# Defect D: a structure-less record outranks a structure
# ─────────────────────────────────────────────────────────────────────────────


_EPRINOMECTIN_CAS = "123997-26-2"
_COMPTOX_MIXTURE = {"DTXSID": "DTXSID8037542", "PREFERRED_NAME": "Eprinomectin", "CASRN": _EPRINOMECTIN_CAS}
_PUBCHEM_EPRINOMECTIN = {
    "cid": 6450531, "cmpdname": "Eprinomectin", "smiles": ASPIRIN, "inchikey": ASPIRIN_IK,
    "cas_numbers": [_EPRINOMECTIN_CAS],
}


class TestStructureless:
    def test_a_structureless_record_joins_the_structure_of_its_substance(self):
        s = _cas_search(comptox=_comptox(_COMPTOX_MIXTURE), pubchem=_pubchem(_PUBCHEM_EPRINOMECTIN))
        hit = s.search(_EPRINOMECTIN_CAS).iloc[0]
        assert hit["canonical_smiles"] == ASPIRIN
        assert hit["DTXSID"] == "DTXSID8037542"
        assert hit["n_clusters"] == 1 and hit["n_source_support"] == 2

    def test_a_structureless_answer_scores_lower(self):
        hit = _cas_search(comptox=_comptox(_COMPTOX_MIXTURE)).search(_EPRINOMECTIN_CAS).iloc[0]
        assert hit["structure_status"] == "none"
        assert hit["confidence"] == pytest.approx(0.765 * 0.8)

    def test_two_structures_leave_the_structureless_record_alone(self):
        def tag(cand, key):
            cand.update(_source_key=key, _origin_rank=0, _match_method="exact_cas", query_match_score=1.0)
            return cand

        pool = [
            tag(make_candidate("CompTox", name="x", cas_candidates=["1-1-1"]), "comptox"),
            tag(make_candidate("PubChemID", name="y", smiles=ASPIRIN, cas_candidates=["1-1-1"]), "pubchem"),
            tag(make_candidate("ChEBI", name="z", smiles="CCO", cas_candidates=["1-1-1"]), "chebi"),
        ]
        assert len(_cluster_candidates(pool)) == 3


# ─────────────────────────────────────────────────────────────────────────────
# Defect F: validation and provenance
# ─────────────────────────────────────────────────────────────────────────────


class TestValidate:
    def _search(self, **validate):
        return _cas_search(comptox=_comptox(_COMPTOX_WRONG), pubchem=_pubchem(_PUBCHEM_RIGHT),
                           validate=validate)

    def test_a_failing_hit_falls_through_to_the_next(self):
        table = pd.DataFrame({"CAS": [_CEFQUINOME_CAS], "MW": [46.07]})
        hit = self._search(mw="MW").search(table, column="CAS").iloc[0]
        assert hit["InChIKey"] == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
        assert hit["validated_by"] == "mw"

    def test_formula_and_reported_name(self):
        table = pd.DataFrame({"CAS": [_CEFQUINOME_CAS], "F": ["C9H8O4"], "N": ["cefquinome sulfate"]})
        hit = self._search(formula="F", name="N").search(table, column="CAS").iloc[0]
        assert hit["InChIKey"] == ASPIRIN_IK
        assert hit["validated_by"] == "formula+name"

    def test_nothing_passing_is_a_no_match_row(self):
        table = pd.DataFrame({"CAS": [_CEFQUINOME_CAS], "MW": [500.0]})
        hit = self._search(mw="MW").search(table, column="CAS").iloc[0]
        assert hit["canonical_smiles"] is None and hit["n_clusters"] == 2

    def test_an_empty_reference_is_not_checked(self):
        table = pd.DataFrame({"CAS": [_CEFQUINOME_CAS], "MW": [None]})
        assert self._search(mw="MW").search(table, column="CAS").iloc[0]["validated_by"] == ""

    def test_validate_needs_a_table(self):
        with pytest.raises(ValueError, match="validate"):
            self._search(mw="MW").search(_CEFQUINOME_CAS)

    def test_unknown_keys_are_refused(self):
        with pytest.raises(ValueError, match="validate keys"):
            Search("cas", validate={"weight": "MW"})


def test_dataset_versions_are_recorded():
    df = _cas_search(comptox=_comptox(_COMPTOX_WRONG)).search(_CEFQUINOME_CAS)
    versions = df.attrs["dataset_versions"]
    assert "provesid" in versions and "comptox" in versions


# ─────────────────────────────────────────────────────────────────────────────
# The regression set of the Sanofi structure audit, against the installed data
# ─────────────────────────────────────────────────────────────────────────────


_NEEDED = ("chebi", "comptox", "pubchem")

# (CAS, accepted parent formulas, defect guarded against).  Formulas compare by
# element counts after salts are stripped; see _parent_formula.
REGRESSION_SET = [
    ("137862-53-4", {"C24H29N5O3"}, "A: valsartan, tetrazole N-H"),
    ("315-30-0", {"C5H4N4O"}, "A: allopurinol"),
    ("446-86-6", {"C9H7N7O2S"}, "A: azathioprine"),
    ("51-21-8", {"C4H3FN2O2"}, "A: 5-fluorouracil"),
    ("103577-45-3", {"C16H14F3N3O2S"}, "A: lansoprazole"),
    ("27314-97-2", {"C7H6N4O2"}, "C: tirapazamine, wrong CompTox formula"),
    ("171500-79-1", {"C88H100Cl2N10O28"}, "B: dalbavancin"),
    ("864070-44-0", {"C23H27ClO7"}, "B: empagliflozin"),
    ("69739-16-8", {"C20H20N6O7S4"}, "B: cefodizime"),
    ("82248-59-7", {"C17H21NO"}, "B: atomoxetine HCl"),
    ("343306-79-6", {"C72H104O48S8"}, "B: sugammadex (the octa-anion), not lumacaftor"),
    # PubChem holds eprinomectin B1b (C49); B1a is C50.
    ("123997-26-2", {"C50H75NO14", "C49H73NO14"}, "D: eprinomectin"),
    # Gentamicin C1, C1a or C2.
    ("1403-66-3", {"C21H43N5O7", "C19H39N5O7", "C20H41N5O7"}, "D: gentamicin"),
]


def _parent_formula(smiles):
    """Element counts of a structure's largest fragment once salts are stripped."""
    from provesid.search import strip_salts
    from provesid.tools import structure_formula

    return formula_element_counts(structure_formula(strip_salts(smiles)))


@pytest.fixture(scope="module")
def audit():
    missing = [name for name in _NEEDED if not datasets.is_present(name)]
    if missing:
        pytest.skip(f"offline datasets not installed: {missing}")
    with Search("cas", show_progress=False) as search:
        yield search.search([cas for cas, _, _ in REGRESSION_SET + [("118443-89-3", None, None)]])


@pytest.mark.integration
@pytest.mark.parametrize("cas, formulas, why", REGRESSION_SET, ids=[r[0] for r in REGRESSION_SET])
def test_audit_regression_set(audit, cas, formulas, why):
    hit = audit.set_index("query").loc[cas]
    assert hit["canonical_smiles"] is not None, why
    assert hit["structure_status"] in ("valid", "repaired:molfile"), why
    assert _parent_formula(hit["canonical_smiles"]) in [formula_element_counts(f) for f in formulas], why


@pytest.mark.integration
def test_audit_no_inchikey_without_a_structure(audit):
    unflagged = audit[audit["InChIKey"].notna() & audit["canonical_smiles"].isna()
                      & ~audit["structure_status"].isin(["invalid", "none"])]
    assert unflagged.empty


@pytest.mark.integration
def test_audit_cefquinome_is_flagged(audit):
    # CompTox's record is a 1,2,4-thiadiazole; PubChem has cefquinome.
    hit = audit.set_index("query").loc["118443-89-3"]
    assert hit["conflict"] or _parent_formula(hit["canonical_smiles"]) == formula_element_counts("C23H24N6O5S2")


@pytest.mark.integration
def test_audit_dtxsid_queries_are_cross_checked():
    missing = [name for name in _NEEDED if not datasets.is_present(name)]
    if missing:
        pytest.skip(f"offline datasets not installed: {missing}")
    with Search("dtxsid", show_progress=False) as search:
        hits = search.search(["DTXSID9020112", "DTXSID30865361", "DTXSID6048593"]).set_index("query")
    # Atrazine: ChEBI and PubChem hold its CAS number and agree with CompTox.
    assert hits.loc["DTXSID9020112", "n_source_support"] >= 3
    # Tirapazamine: CompTox's C8H7N3O2 loses to ChEBI and PubChem.
    tirapazamine = hits.loc["DTXSID30865361"]
    assert tirapazamine["conflict"]
    assert _parent_formula(tirapazamine["canonical_smiles"]) == formula_element_counts("C7H6N4O2")
    # Cefquinome sulfate: no longer confirmed by an InChIKey echo.
    assert hits.loc["DTXSID6048593", "conflict"]
