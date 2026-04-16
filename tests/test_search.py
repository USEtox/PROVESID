"""Unit tests for provesid.search — Search class and related utilities.

Follows the mock-based pattern from test_tools.py: offline source classes are
replaced with lightweight stubs so tests run without real database files.

Run with::

    uv run pytest tests/test_search.py -v
"""

from __future__ import annotations

import pandas as pd
import pytest

import provesid.search as search_module
from provesid.search import (
    OUTPUT_COLUMNS,
    Search,
    _any_candidate,
    _most_complete_row,
    normalize_structure,
    strip_salts,
)


# ─────────────────────────────────────────────────────────────────────────────
# Stub source classes (no real database access)
# ─────────────────────────────────────────────────────────────────────────────

class _ChebiStub:
    """Minimal ChebiSDF stub.  Returns canned data or empty results."""

    def __init__(self, row=None, rows=None):
        self._row = row
        self._rows = rows or ([] if row is None else [row])

    def search_by_cas(self, cas):
        return self._rows

    def search_by_name(self, name, exact=True):
        return self._rows

    def search_by_synonym(self, name, exact=True):
        return self._rows

    def search_by_inchikey(self, inchikey):
        return self._row

    def search_by_inchi(self, inchi):
        return self._row

    def search_by_formula(self, formula):
        return self._rows


class _CompToxStub:
    """Minimal CompToxID stub."""

    def __init__(self, row=None, rows=None):
        self._row = row
        self._rows = rows or ([] if row is None else [row])

    def get_by_casrn(self, cas):
        return self._row

    def get_by_inchikey(self, ik):
        return self._row

    def get_by_smiles(self, smiles):
        return self._row

    def get_by_name(self, name):
        return self._row

    def search_by_name(self, name, exact=False, limit=10):
        return self._rows

    def get_by_dtxsid(self, dtxsid):
        return self._row

    def search_by_formula(self, formula, limit=100):
        return self._rows


class _PubChemStub:
    """Minimal PubChemID stub."""

    def __init__(self, row=None, rows=None):
        self._row = row
        self._rows = rows or ([] if row is None else [row])

    def get_by_cas(self, cas):
        return self._row

    def get_by_inchikey(self, ik):
        return self._row

    def get_by_smiles(self, smiles):
        return self._row

    def get_by_inchi(self, inchi):
        return self._row

    def search_by_name(self, name, exact=False, limit=10):
        return self._rows

    def search_by_formula(self, formula, limit=100):
        return self._rows


class _ZeroPMStub:
    """Minimal ZeroPM stub."""

    def __init__(self, table=None, cas_list=None):
        self._table = table if table is not None else pd.DataFrame()
        self._cas_list = cas_list or []

    def get_id_table_from_cas(self, cas):
        return self._table

    def get_id_table_from_name(self, name):
        return self._table

    def get_id_table_from_inchikey(self, ik):
        return self._table

    def get_id_table_from_inchi(self, inchi):
        return self._table

    def get_cas_from_smiles(self, smiles):
        return self._cas_list

    def get_cas_from_formula(self, formula):
        return self._cas_list

    def query_similar_name(self, name):
        return self._table


class _ChEMBLStub:
    """Minimal CheMBL stub."""

    def __init__(self, compound=None, props=None):
        self._compound = compound
        self._props = props or {}

    def search_by_smiles(self, smiles):
        return self._compound

    def search_by_name(self, name, limit=100):
        return [self._compound] if self._compound else []

    def search_by_inchikey(self, ik):
        return self._compound

    def get_properties(self, molregno):
        return self._props


# ─────────────────────────────────────────────────────────────────────────────
# Canned compound records
# ─────────────────────────────────────────────────────────────────────────────

_ASPIRIN_CHEBI = {
    "ChEBI NAME": "aspirin",
    "FORMULA": "C9H8O4",
    "SMILES": "CC(=O)Oc1ccccc1C(=O)O",
    "INCHI": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
    "INCHIKEY": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
    "SYNONYM": ["Acetylsalicylic acid", "2-acetoxybenzoic acid"],
}

_ASPIRIN_COMPTOX = {
    "PREFERRED_NAME": "Aspirin",
    "IUPAC_NAME": "2-(acetyloxy)benzoic acid",
    "MOLECULAR_FORMULA": "C9H8O4",
    "SMILES": "CC(=O)Oc1ccccc1C(=O)O",
    "INCHI": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
    "INCHIKEY": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
    "DTXSID": "DTXSID2021735",
    "CASRN": "50-78-2",
    "identifiers": "50-78-2; Acetylsalicylic acid",
    "AVERAGE_MASS": 180.16,
}

_ASPIRIN_PUBCHEM = {
    "cmpdname": "Aspirin",
    "iupacname": "2-(acetyloxy)benzoic acid",
    "mf": "C9H8O4",
    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "inchi": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
    "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
    "synonyms": ["Aspirin", "Acetylsalicylic acid"],
    "cas_numbers": ["50-78-2"],
    "mw": 180.16,
}

_ASPIRIN_ZEROPM = pd.DataFrame(
    {
        "cas": ["50-78-2"],
        "inchi": ["InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"],
        "inchikey": ["BSYNRYMUTXBXSQ-UHFFFAOYSA-N"],
        "rank": [1],
        "name": ["aspirin"],
    }
)

_ASPIRIN_CHEMBL = {
    "molregno": 1,
    "pref_name": "ASPIRIN",
    "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "standard_inchi": "InChI=1S/C9H8O4",
    "standard_inchi_key": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
    "synonyms": [],
}


def _make_search(**clients) -> Search:
    """Build a Search("cas") instance injecting the given stub clients."""
    return Search(
        "cas",
        show_progress=False,
        chebi=clients.get("chebi", _ChebiStub()),
        comptox=clients.get("comptox", _CompToxStub()),
        pubchem=clients.get("pubchem", _PubChemStub()),
        zeropm=clients.get("zeropm", _ZeroPMStub()),
        chembl=clients.get("chembl", _ChEMBLStub()),
    )


# ─────────────────────────────────────────────────────────────────────────────
# normalize_structure
# ─────────────────────────────────────────────────────────────────────────────

class TestNormalizeStructure:
    def test_canonical_smiles_produced(self):
        rec = normalize_structure("c1ccccc1")
        assert rec["canonical_smiles"] is not None
        assert "c" in rec["canonical_smiles"] or "C" in rec["canonical_smiles"]

    def test_kekulized_smiles_produced(self):
        rec = normalize_structure("c1ccccc1")
        assert rec["kekulized_smiles"] is not None
        # Kekulized benzene uses uppercase C=C
        assert "C" in rec["kekulized_smiles"]

    def test_inchikey_produced(self):
        rec = normalize_structure("CC(=O)Oc1ccccc1C(=O)O")
        assert rec["inchikey"] == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

    def test_inchi_produced(self):
        rec = normalize_structure("CC(=O)Oc1ccccc1C(=O)O")
        assert rec["inchi"] is not None
        assert rec["inchi"].startswith("InChI=")

    def test_mol_weight_produced(self):
        rec = normalize_structure("CC(=O)Oc1ccccc1C(=O)O")
        assert rec["mol_weight"] is not None
        assert 179 < rec["mol_weight"] < 182

    def test_invalid_smiles_returns_none_fields(self):
        rec = normalize_structure("this-is-not-a-smiles!!!")
        assert rec["canonical_smiles"] is None
        assert rec["inchikey"] is None

    def test_none_input_returns_empty(self):
        rec = normalize_structure(None)
        assert rec["canonical_smiles"] is None

    def test_empty_string_returns_empty(self):
        rec = normalize_structure("")
        assert rec["canonical_smiles"] is None


# ─────────────────────────────────────────────────────────────────────────────
# strip_salts
# ─────────────────────────────────────────────────────────────────────────────

class TestStripSalts:
    def test_nacl_salt_removed(self):
        result = strip_salts("[Na+].[Cl-].CC(=O)O")
        # Parent should be acetic acid
        assert result is not None
        norm = normalize_structure(result)
        # Acetic acid canonical SMILES
        assert norm["inchikey"] is not None

    def test_largest_fragment_selected(self):
        # Two fragments: ethanol (larger) and water (smaller)
        result = strip_salts("CCO.O")
        assert result is not None
        norm = normalize_structure(result)
        assert norm["mol_weight"] is not None
        # Should keep ethanol (~46 Da), not water (~18 Da)
        assert norm["mol_weight"] > 30

    def test_single_fragment_unchanged(self):
        smiles = "CC(=O)Oc1ccccc1C(=O)O"
        result = strip_salts(smiles)
        assert result is not None

    def test_none_input_returns_none(self):
        result = strip_salts(None)
        assert result is None

    def test_extra_smarts_applied(self):
        # Extra SMARTS "[Na]" should strip sodium
        result = strip_salts("[Na]CC(=O)O", extra_smarts=["[Na]"])
        assert result is not None


# ─────────────────────────────────────────────────────────────────────────────
# Search initialisation
# ─────────────────────────────────────────────────────────────────────────────

class TestSearchInit:
    def test_default_identifier_type(self):
        s = Search(show_progress=False)
        assert s.identifier_type == "cas"

    def test_custom_identifier_type(self):
        s = Search("smiles", show_progress=False)
        assert s.identifier_type == "smiles"

    def test_invalid_identifier_type_raises(self):
        with pytest.raises(ValueError, match="identifier_type must be one of"):
            Search("invalid_type")

    def test_injected_clients_not_reinitialised(self):
        chebi = _ChebiStub()
        s = Search("cas", chebi=chebi, show_progress=False)
        assert s._chebi is chebi
        assert s._clients_initialized

    def test_strip_salts_default_false(self):
        s = Search(show_progress=False)
        assert s.strip_salts is False

    def test_fuzzy_default_false(self):
        s = Search(show_progress=False)
        assert s.fuzzy is False

    def test_similarity_threshold_default_zero(self):
        s = Search(show_progress=False)
        assert s.similarity_threshold == 0.0

    def test_inchikey_skeleton_default_false(self):
        s = Search(show_progress=False)
        assert s.inchikey_skeleton is False


# ─────────────────────────────────────────────────────────────────────────────
# Output schema
# ─────────────────────────────────────────────────────────────────────────────

class TestOutputSchema:
    def test_all_output_columns_present(self):
        s = _make_search()
        df = s.search("50-78-2")
        for col in OUTPUT_COLUMNS:
            assert col in df.columns, f"Missing column: {col}"

    def test_single_row_for_single_query(self):
        s = _make_search()
        df = s.search("50-78-2")
        assert len(df) == 1

    def test_one_row_per_query_in_list(self):
        s = _make_search()
        df = s.search(["50-78-2", "64-17-5"])
        assert len(df) == 2

    def test_query_column_preserves_input(self):
        s = _make_search()
        df = s.search("50-78-2")
        assert df.iloc[0]["query"] == "50-78-2"

    def test_dataframe_input_with_column(self):
        s = _make_search()
        input_df = pd.DataFrame({"cas": ["50-78-2"], "extra": ["foo"]})
        df = s.search(input_df, column="cas")
        assert len(df) == 1
        assert "query" in df.columns

    def test_dataframe_input_without_column_raises(self):
        s = _make_search()
        with pytest.raises(ValueError, match="column="):
            s.search(pd.DataFrame({"cas": ["50-78-2"]}))

    def test_csv_file_input(self, tmp_path):
        csv_path = tmp_path / "compounds.csv"
        pd.DataFrame({"CAS": ["50-78-2", "64-17-5"]}).to_csv(csv_path, index=False)
        s = _make_search()
        df = s.search(csv_path, column="CAS")
        assert len(df) == 2

    def test_csv_file_input_without_column_raises(self, tmp_path):
        csv_path = tmp_path / "compounds.csv"
        pd.DataFrame({"CAS": ["50-78-2"]}).to_csv(csv_path, index=False)
        s = _make_search()
        with pytest.raises(ValueError, match="column="):
            s.search(csv_path)


# ─────────────────────────────────────────────────────────────────────────────
# CAS resolver
# ─────────────────────────────────────────────────────────────────────────────

class TestResolveCas:
    def test_chebi_hit_populates_fields(self):
        s = _make_search(chebi=_ChebiStub(_ASPIRIN_CHEBI))
        df = s.search("50-78-2")
        row = df.iloc[0]
        assert row["name"] == "aspirin"
        assert row["molecular_formula"] == "C9H8O4"
        assert row["InChIKey"] == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

    def test_source_set_to_chebi(self):
        s = _make_search(chebi=_ChebiStub(_ASPIRIN_CHEBI))
        df = s.search("50-78-2")
        assert df.iloc[0]["source"] == "ChEBI"

    def test_comptox_fallback_when_chebi_empty(self):
        s = _make_search(
            chebi=_ChebiStub(),
            comptox=_CompToxStub(_ASPIRIN_COMPTOX),
        )
        df = s.search("50-78-2")
        row = df.iloc[0]
        assert row["DTXSID"] == "DTXSID2021735"
        assert row["source"] == "CompTox"

    def test_pubchem_fallback(self):
        s = _make_search(
            chebi=_ChebiStub(),
            comptox=_CompToxStub(),
            pubchem=_PubChemStub(_ASPIRIN_PUBCHEM),
        )
        df = s.search("50-78-2")
        row = df.iloc[0]
        assert row["name"] == "Aspirin"
        assert row["source"] == "PubChemID"

    def test_inchikey_always_derived_from_smiles(self):
        """InChIKey should be populated from RDKit even if source omits it."""
        rec = dict(_ASPIRIN_CHEBI)
        rec.pop("INCHIKEY", None)
        s = _make_search(chebi=_ChebiStub(rec))
        df = s.search("50-78-2")
        assert df.iloc[0]["InChIKey"] == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

    def test_kekulized_smiles_populated(self):
        s = _make_search(chebi=_ChebiStub(_ASPIRIN_CHEBI))
        df = s.search("50-78-2")
        assert df.iloc[0]["kekulized_smiles"] is not None

    def test_canonical_smiles_populated(self):
        s = _make_search(chebi=_ChebiStub(_ASPIRIN_CHEBI))
        df = s.search("50-78-2")
        assert df.iloc[0]["canonical_smiles"] is not None

    def test_confidence_positive(self):
        s = _make_search(chebi=_ChebiStub(_ASPIRIN_CHEBI))
        df = s.search("50-78-2")
        assert df.iloc[0]["confidence"] > 0

    def test_match_method_exact_cas(self):
        s = _make_search()
        df = s.search("50-78-2")
        assert df.iloc[0]["match_method"] == "exact_cas"

    def test_source_details_present(self):
        s = _make_search(chebi=_ChebiStub(_ASPIRIN_CHEBI))
        df = s.search("50-78-2")
        sd = df.iloc[0]["source_details"]
        assert "ChEBI" in sd
        assert sd["ChEBI"]["found"] is True

    def test_source_details_no_match_sources_found_false(self):
        s = _make_search()
        df = s.search("50-78-2")
        sd = df.iloc[0]["source_details"]
        for src_info in sd.values():
            assert src_info["found"] is False

    def test_foundby_casrn(self):
        s = _make_search()
        df = s.search("50-78-2")
        assert df.iloc[0]["foundby"] == "CASRN"


# ─────────────────────────────────────────────────────────────────────────────
# Name resolver
# ─────────────────────────────────────────────────────────────────────────────

class TestResolveName:
    def test_exact_name_hit(self):
        s = Search(
            "name",
            show_progress=False,
            chebi=_ChebiStub(_ASPIRIN_CHEBI),
            comptox=_CompToxStub(),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search("aspirin")
        row = df.iloc[0]
        assert row["name"] == "aspirin"
        assert row["match_method"] == "exact_name"

    def test_foundby_name(self):
        s = Search(
            "name",
            show_progress=False,
            chebi=_ChebiStub(_ASPIRIN_CHEBI),
            comptox=_CompToxStub(),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search("aspirin")
        assert df.iloc[0]["foundby"] == "name"

    def test_fuzzy_not_triggered_when_exact_match_found(self):
        s = Search(
            "name",
            fuzzy=True,
            show_progress=False,
            chebi=_ChebiStub(_ASPIRIN_CHEBI),
            comptox=_CompToxStub(),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search("aspirin")
        assert df.iloc[0]["match_method"] == "exact_name"


# ─────────────────────────────────────────────────────────────────────────────
# SMILES resolver
# ─────────────────────────────────────────────────────────────────────────────

class TestResolveSmiles:
    def test_foundby_smiles(self):
        s = Search(
            "smiles",
            show_progress=False,
            chebi=_ChebiStub(_ASPIRIN_CHEBI),
            comptox=_CompToxStub(),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search("CC(=O)Oc1ccccc1C(=O)O")
        assert df.iloc[0]["foundby"] == "SMILES"

    def test_match_method_exact_smiles(self):
        s = Search(
            "smiles",
            show_progress=False,
            chebi=_ChebiStub(_ASPIRIN_CHEBI),
            comptox=_CompToxStub(),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search("CC(=O)Oc1ccccc1C(=O)O")
        assert df.iloc[0]["match_method"] == "exact_smiles"

    def test_salt_stripping_populates_parent_fields(self):
        s = Search(
            "smiles",
            strip_salts=True,
            show_progress=False,
            chebi=_ChebiStub(),
            comptox=_CompToxStub(),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        # NaCl salt + acetic acid
        df = s.search("[Na+].[Cl-].CC(=O)O")
        row = df.iloc[0]
        # parent_smiles should be acetic acid
        assert row["parent_smiles"] is not None
        assert row["parent_inchikey"] is not None

    def test_strip_salts_off_by_default(self):
        s = Search(
            "smiles",
            show_progress=False,
            chebi=_ChebiStub(),
            comptox=_CompToxStub(),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search("[Na+].[Cl-].CC(=O)O")
        row = df.iloc[0]
        assert row["parent_smiles"] is None


# ─────────────────────────────────────────────────────────────────────────────
# InChI resolver
# ─────────────────────────────────────────────────────────────────────────────

class TestResolveInchi:
    _ASPIRIN_INCHI = (
        "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
    )

    def test_foundby_inchi(self):
        s = Search(
            "inchi",
            show_progress=False,
            chebi=_ChebiStub(_ASPIRIN_CHEBI),
            comptox=_CompToxStub(),
            pubchem=_PubChemStub(_ASPIRIN_PUBCHEM),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search(self._ASPIRIN_INCHI)
        assert df.iloc[0]["foundby"] == "InChI"

    def test_match_method_inchi(self):
        s = Search(
            "inchi",
            show_progress=False,
            chebi=_ChebiStub(_ASPIRIN_CHEBI),
            comptox=_CompToxStub(),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search(self._ASPIRIN_INCHI)
        assert df.iloc[0]["match_method"] == "inchi"

    def test_inchikey_derived(self):
        s = Search(
            "inchi",
            show_progress=False,
            chebi=_ChebiStub(),
            comptox=_CompToxStub(),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search(self._ASPIRIN_INCHI)
        # Even with no source match, InChIKey should be derived from InChI
        assert df.iloc[0]["InChIKey"] == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"


# ─────────────────────────────────────────────────────────────────────────────
# InChIKey resolver
# ─────────────────────────────────────────────────────────────────────────────

class TestResolveInchikey:
    _ASPIRIN_IK = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

    def test_foundby_inchikey(self):
        s = Search(
            "inchikey",
            show_progress=False,
            chebi=_ChebiStub(_ASPIRIN_CHEBI),
            comptox=_CompToxStub(),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search(self._ASPIRIN_IK)
        assert df.iloc[0]["foundby"] == "InChIKey"

    def test_match_method_exact_inchikey(self):
        s = Search(
            "inchikey",
            show_progress=False,
            chebi=_ChebiStub(_ASPIRIN_CHEBI),
            comptox=_CompToxStub(),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search(self._ASPIRIN_IK)
        assert df.iloc[0]["match_method"] == "exact_inchikey"

    def test_confidence_near_1_for_inchikey_hit(self):
        s = Search(
            "inchikey",
            show_progress=False,
            chebi=_ChebiStub(_ASPIRIN_CHEBI),
            comptox=_CompToxStub(),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search(self._ASPIRIN_IK)
        # Base = 1.0, modulated by consensus but should be high
        assert df.iloc[0]["confidence"] >= 0.5

    def test_inchikey_skeleton_fallback(self):
        """When exact match fails and inchikey_skeleton=True, skeleton is tried."""
        s = Search(
            "inchikey",
            inchikey_skeleton=True,
            show_progress=False,
            chebi=_ChebiStub(),  # no exact match
            comptox=_CompToxStub(),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search(self._ASPIRIN_IK)
        # With no match even via skeleton (all stubs empty), match_method still
        # starts as "exact_inchikey" then stays (no skeleton match found either)
        assert df.iloc[0]["match_method"] in ("exact_inchikey", "inchikey_skeleton")

    def test_no_match_confidence_zero(self):
        s = Search(
            "inchikey",
            show_progress=False,
            chebi=_ChebiStub(),
            comptox=_CompToxStub(),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search(self._ASPIRIN_IK)
        # No sources matched; confidence should be at the floor
        assert df.iloc[0]["confidence"] <= 0.5


# ─────────────────────────────────────────────────────────────────────────────
# DTXSID resolver
# ─────────────────────────────────────────────────────────────────────────────

class TestResolveDtxsid:
    def test_foundby_dtxsid(self):
        s = Search(
            "dtxsid",
            show_progress=False,
            chebi=_ChebiStub(_ASPIRIN_CHEBI),
            comptox=_CompToxStub(_ASPIRIN_COMPTOX),
            pubchem=_PubChemStub(_ASPIRIN_PUBCHEM),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search("DTXSID2021735")
        assert df.iloc[0]["foundby"] == "DTXSID"

    def test_dtxsid_preserved_in_result(self):
        s = Search(
            "dtxsid",
            show_progress=False,
            chebi=_ChebiStub(),
            comptox=_CompToxStub(_ASPIRIN_COMPTOX),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search("DTXSID2021735")
        assert df.iloc[0]["DTXSID"] == "DTXSID2021735"

    def test_match_method_dtxsid(self):
        s = Search(
            "dtxsid",
            show_progress=False,
            chebi=_ChebiStub(),
            comptox=_CompToxStub(_ASPIRIN_COMPTOX),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search("DTXSID2021735")
        assert df.iloc[0]["match_method"] == "dtxsid"


# ─────────────────────────────────────────────────────────────────────────────
# Formula resolver
# ─────────────────────────────────────────────────────────────────────────────

class TestResolveFormula:
    def test_foundby_formula(self):
        s = Search(
            "formula",
            show_progress=False,
            chebi=_ChebiStub(_ASPIRIN_CHEBI),
            comptox=_CompToxStub(_ASPIRIN_COMPTOX),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search("C9H8O4")
        assert df.iloc[0]["foundby"] == "formula"

    def test_confidence_capped_at_formula_base(self):
        """Formula matches should have low base confidence (0.30)."""
        s = Search(
            "formula",
            show_progress=False,
            chebi=_ChebiStub(_ASPIRIN_CHEBI),
            comptox=_CompToxStub(_ASPIRIN_COMPTOX),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search("C9H8O4")
        # Even with high cross-source agreement, base is 0.30
        assert df.iloc[0]["confidence"] <= 0.30 * 1.0 + 0.001  # slight float tolerance

    def test_match_method_formula(self):
        s = Search(
            "formula",
            show_progress=False,
            chebi=_ChebiStub(_ASPIRIN_CHEBI),
            comptox=_CompToxStub(),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search("C9H8O4")
        assert df.iloc[0]["match_method"] == "formula"


# ─────────────────────────────────────────────────────────────────────────────
# Confidence scoring
# ─────────────────────────────────────────────────────────────────────────────

class TestConfidenceScoring:
    def _search_with_method(self, method, consensus=1.0, fuzzy_score=None, tanimoto=None):
        s = Search(show_progress=False)
        return s._compute_confidence(
            method, consensus, fuzzy_score=fuzzy_score, tanimoto=tanimoto
        )

    def test_exact_inchikey_high_confidence(self):
        c = self._search_with_method("exact_inchikey", consensus=1.0)
        assert c == 1.0

    def test_exact_inchikey_with_low_consensus(self):
        c = self._search_with_method("exact_inchikey", consensus=0.0)
        assert c == pytest.approx(0.5)

    def test_exact_cas_confidence(self):
        c = self._search_with_method("exact_cas", consensus=1.0)
        assert c == pytest.approx(0.9)

    def test_formula_confidence_capped(self):
        c = self._search_with_method("formula", consensus=1.0)
        assert c == pytest.approx(0.3)

    def test_fuzzy_name_uses_fuzzy_score(self):
        c = self._search_with_method("fuzzy_name", consensus=1.0, fuzzy_score=0.8)
        assert c == pytest.approx(0.8)

    def test_fuzzy_name_low_score(self):
        c = self._search_with_method("fuzzy_name", consensus=1.0, fuzzy_score=0.2)
        assert c == pytest.approx(0.2)

    def test_tanimoto_scaled(self):
        c = self._search_with_method("tanimoto", consensus=1.0, tanimoto=1.0)
        assert c == pytest.approx(0.85)

    def test_consensus_modulation(self):
        """Same match method with different consensus yields different confidence."""
        c_high = self._search_with_method("exact_name", consensus=1.0)
        c_low = self._search_with_method("exact_name", consensus=0.0)
        assert c_high > c_low

    def test_confidence_clamped_to_0_1(self):
        s = Search(show_progress=False)
        assert s._compute_confidence("unknown", 99.0) <= 1.0
        assert s._compute_confidence("unknown", -5.0) >= 0.0


# ─────────────────────────────────────────────────────────────────────────────
# Name normalisation
# ─────────────────────────────────────────────────────────────────────────────

class TestNormalizeName:
    def _norm(self, name):
        return Search(show_progress=False)._normalize_name(name)

    def test_lowercase(self):
        assert self._norm("ASPIRIN") == "aspirin"

    def test_strip_whitespace(self):
        assert self._norm("  aspirin  ") == "aspirin"

    def test_remove_d_prefix(self):
        assert self._norm("D-Aspirin") == "aspirin"

    def test_remove_l_prefix(self):
        assert self._norm("L-alanine") == "alanine"

    def test_remove_rac_prefix(self):
        assert self._norm("rac-ibuprofen") == "ibuprofen"

    def test_abbreviation_expansion(self):
        assert self._norm("MEK") == "methyl ethyl ketone"

    def test_collapse_multiple_spaces(self):
        assert self._norm("methyl  ethyl  ketone") == "methyl ethyl ketone"


# ─────────────────────────────────────────────────────────────────────────────
# Source details
# ─────────────────────────────────────────────────────────────────────────────

class TestSourceDetails:
    def test_found_true_for_matched_source(self):
        s = _make_search(chebi=_ChebiStub(_ASPIRIN_CHEBI))
        df = s.search("50-78-2")
        sd = df.iloc[0]["source_details"]
        assert sd["ChEBI"]["found"] is True

    def test_found_false_for_unmatched_sources(self):
        s = _make_search(
            chebi=_ChebiStub(_ASPIRIN_CHEBI),
            comptox=_CompToxStub(),
            pubchem=_PubChemStub(),
        )
        df = s.search("50-78-2")
        sd = df.iloc[0]["source_details"]
        assert sd["CompTox"]["found"] is False
        assert sd["PubChemID"]["found"] is False

    def test_fields_list_non_empty_for_matched_source(self):
        s = _make_search(chebi=_ChebiStub(_ASPIRIN_CHEBI))
        df = s.search("50-78-2")
        sd = df.iloc[0]["source_details"]
        assert len(sd["ChEBI"]["fields"]) > 0

    def test_all_five_sources_present_in_source_details(self):
        s = _make_search()
        df = s.search("50-78-2")
        sd = df.iloc[0]["source_details"]
        for src in ["ChEBI", "CompTox", "PubChemID", "ZeroPM", "ChEMBL"]:
            assert src in sd

    def test_smiles_field_recorded(self):
        s = _make_search(chebi=_ChebiStub(_ASPIRIN_CHEBI))
        df = s.search("50-78-2")
        sd = df.iloc[0]["source_details"]
        assert "SMILES" in sd["ChEBI"]["fields"]


# ─────────────────────────────────────────────────────────────────────────────
# Salt stripping integration
# ─────────────────────────────────────────────────────────────────────────────

class TestSaltStrippingIntegration:
    def test_parent_smiles_present_when_strip_salts_true(self):
        s = Search(
            "smiles",
            strip_salts=True,
            show_progress=False,
            chebi=_ChebiStub(),
            comptox=_CompToxStub(),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search("[Na+].[Cl-].CC(=O)O")
        assert df.iloc[0]["parent_smiles"] is not None

    def test_parent_inchikey_present_when_strip_salts_true(self):
        s = Search(
            "smiles",
            strip_salts=True,
            show_progress=False,
            chebi=_ChebiStub(),
            comptox=_CompToxStub(),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search("[Na+].[Cl-].CC(=O)O")
        assert df.iloc[0]["parent_inchikey"] is not None

    def test_no_parent_for_single_component_smiles(self):
        s = Search(
            "smiles",
            strip_salts=True,
            show_progress=False,
            chebi=_ChebiStub(),
            comptox=_CompToxStub(),
            pubchem=_PubChemStub(),
            zeropm=_ZeroPMStub(),
            chembl=_ChEMBLStub(),
        )
        df = s.search("CC(=O)Oc1ccccc1C(=O)O")
        # No salt to strip → parent_smiles should be None (same as canonical)
        assert df.iloc[0]["parent_smiles"] is None


# ─────────────────────────────────────────────────────────────────────────────
# Helper function unit tests
# ─────────────────────────────────────────────────────────────────────────────

class TestHelpers:
    def test_any_candidate_true(self):
        assert _any_candidate({"a": {"x": 1}, "b": None}) is True

    def test_any_candidate_false(self):
        assert _any_candidate({"a": None, "b": None}) is False

    def test_any_candidate_empty(self):
        assert _any_candidate({}) is False

    def test_most_complete_row_picks_best(self):
        rows = [
            {"a": 1, "b": None, "c": None},
            {"a": 1, "b": 2, "c": 3},
            {"a": 1, "b": 2, "c": None},
        ]
        best = _most_complete_row(rows)
        assert best["c"] == 3

    def test_most_complete_row_single_item(self):
        rows = [{"a": 1}]
        assert _most_complete_row(rows) == {"a": 1}

    def test_most_complete_row_empty(self):
        assert _most_complete_row([]) == {}


# ─────────────────────────────────────────────────────────────────────────────
# Multi-source consensus
# ─────────────────────────────────────────────────────────────────────────────

class TestMultiSourceConsensus:
    def test_two_agreeing_sources_boost_confidence(self):
        """Two sources with matching SMILES → higher confidence than one source."""
        s_single = _make_search(chebi=_ChebiStub(_ASPIRIN_CHEBI))
        df_single = s_single.search("50-78-2")
        conf_single = df_single.iloc[0]["confidence"]

        s_two = _make_search(
            chebi=_ChebiStub(_ASPIRIN_CHEBI),
            comptox=_CompToxStub(_ASPIRIN_COMPTOX),
        )
        df_two = s_two.search("50-78-2")
        conf_two = df_two.iloc[0]["confidence"]

        assert conf_two >= conf_single

    def test_match_score_positive_with_two_sources(self):
        s = _make_search(
            chebi=_ChebiStub(_ASPIRIN_CHEBI),
            comptox=_CompToxStub(_ASPIRIN_COMPTOX),
        )
        df = s.search("50-78-2")
        assert df.iloc[0]["match_score"] > 0

    def test_source_match_scores_dict_present(self):
        s = _make_search(chebi=_ChebiStub(_ASPIRIN_CHEBI))
        df = s.search("50-78-2")
        assert isinstance(df.iloc[0]["source_match_scores"], dict)
