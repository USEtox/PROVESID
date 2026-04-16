"""Tests for source class methods used by the Search class.

Verifies that the methods required by the Search resolvers are present and
return the expected types.  Uses the real source class interfaces but stubs
out the underlying database connections so no real files are needed.

Run with::

    uv run pytest tests/test_search_new_methods.py -v
"""

from __future__ import annotations

import sqlite3
import tempfile
from pathlib import Path
from typing import Any, Dict, Optional
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest


# ─────────────────────────────────────────────────────────────────────────────
# PubChemID method interface tests
# ─────────────────────────────────────────────────────────────────────────────

class TestPubChemIDMethods:
    """Verify PubChemID exposes the methods used by _resolve_inchikey / _resolve_inchi."""

    def test_get_by_inchikey_method_exists(self):
        from provesid.pubchem import PubChemID
        assert hasattr(PubChemID, "get_by_inchikey"), (
            "PubChemID must have get_by_inchikey() for _resolve_inchikey"
        )
        assert callable(PubChemID.get_by_inchikey)

    def test_get_by_inchi_method_exists(self):
        from provesid.pubchem import PubChemID
        assert hasattr(PubChemID, "get_by_inchi"), (
            "PubChemID must have get_by_inchi() for _resolve_inchi"
        )
        assert callable(PubChemID.get_by_inchi)

    def test_search_by_formula_method_exists(self):
        from provesid.pubchem import PubChemID
        assert hasattr(PubChemID, "search_by_formula"), (
            "PubChemID must have search_by_formula() for _resolve_formula"
        )
        assert callable(PubChemID.search_by_formula)

    def test_get_by_smiles_method_exists(self):
        from provesid.pubchem import PubChemID
        assert hasattr(PubChemID, "get_by_smiles"), (
            "PubChemID must have get_by_smiles() for _resolve_smiles"
        )
        assert callable(PubChemID.get_by_smiles)

    def test_search_by_name_method_exists(self):
        from provesid.pubchem import PubChemID
        assert hasattr(PubChemID, "search_by_name")
        assert callable(PubChemID.search_by_name)

    def test_get_by_cas_method_exists(self):
        from provesid.pubchem import PubChemID
        assert hasattr(PubChemID, "get_by_cas")
        assert callable(PubChemID.get_by_cas)


# ─────────────────────────────────────────────────────────────────────────────
# CompToxID method interface tests
# ─────────────────────────────────────────────────────────────────────────────

class TestCompToxIDMethods:
    """Verify CompToxID exposes the methods used by the Search resolvers."""

    def test_get_by_dtxsid_method_exists(self):
        from provesid.comptox import CompToxID
        assert hasattr(CompToxID, "get_by_dtxsid"), (
            "CompToxID must have get_by_dtxsid() for _resolve_dtxsid"
        )
        assert callable(CompToxID.get_by_dtxsid)

    def test_get_by_inchikey_method_exists(self):
        from provesid.comptox import CompToxID
        assert hasattr(CompToxID, "get_by_inchikey")
        assert callable(CompToxID.get_by_inchikey)

    def test_search_by_formula_method_exists(self):
        from provesid.comptox import CompToxID
        assert hasattr(CompToxID, "search_by_formula"), (
            "CompToxID must have search_by_formula() for _resolve_formula"
        )
        assert callable(CompToxID.search_by_formula)

    def test_get_by_casrn_method_exists(self):
        from provesid.comptox import CompToxID
        assert hasattr(CompToxID, "get_by_casrn")
        assert callable(CompToxID.get_by_casrn)

    def test_get_by_smiles_method_exists(self):
        from provesid.comptox import CompToxID
        assert hasattr(CompToxID, "get_by_smiles")
        assert callable(CompToxID.get_by_smiles)

    def test_get_by_name_method_exists(self):
        from provesid.comptox import CompToxID
        assert hasattr(CompToxID, "get_by_name")
        assert callable(CompToxID.get_by_name)

    def test_search_by_name_method_exists(self):
        from provesid.comptox import CompToxID
        assert hasattr(CompToxID, "search_by_name")
        assert callable(CompToxID.search_by_name)


# ─────────────────────────────────────────────────────────────────────────────
# CheMBL method interface tests
# ─────────────────────────────────────────────────────────────────────────────

class TestCheMBLMethods:
    """Verify CheMBL exposes search_by_inchikey used by _resolve_inchikey."""

    def test_search_by_inchikey_method_exists(self):
        from provesid.chembl import CheMBL
        assert hasattr(CheMBL, "search_by_inchikey"), (
            "CheMBL must have search_by_inchikey() for _resolve_inchikey"
        )
        assert callable(CheMBL.search_by_inchikey)

    def test_search_by_smiles_method_exists(self):
        from provesid.chembl import CheMBL
        assert hasattr(CheMBL, "search_by_smiles")
        assert callable(CheMBL.search_by_smiles)

    def test_search_by_name_method_exists(self):
        from provesid.chembl import CheMBL
        assert hasattr(CheMBL, "search_by_name")
        assert callable(CheMBL.search_by_name)

    def test_get_properties_method_exists(self):
        from provesid.chembl import CheMBL
        assert hasattr(CheMBL, "get_properties")
        assert callable(CheMBL.get_properties)


# ─────────────────────────────────────────────────────────────────────────────
# ChebiSDF method interface tests
# ─────────────────────────────────────────────────────────────────────────────

class TestChebiSDFMethods:
    """Verify ChebiSDF exposes the methods used by the Search resolvers."""

    def test_search_by_formula_method_exists(self):
        from provesid.chebi import ChebiSDF
        assert hasattr(ChebiSDF, "search_by_formula"), (
            "ChebiSDF must have search_by_formula() for _resolve_formula"
        )
        assert callable(ChebiSDF.search_by_formula)

    def test_search_by_inchikey_method_exists(self):
        from provesid.chebi import ChebiSDF
        assert hasattr(ChebiSDF, "search_by_inchikey")
        assert callable(ChebiSDF.search_by_inchikey)

    def test_search_by_inchi_method_exists(self):
        from provesid.chebi import ChebiSDF
        assert hasattr(ChebiSDF, "search_by_inchi")
        assert callable(ChebiSDF.search_by_inchi)

    def test_search_by_cas_method_exists(self):
        from provesid.chebi import ChebiSDF
        assert hasattr(ChebiSDF, "search_by_cas")
        assert callable(ChebiSDF.search_by_cas)

    def test_search_by_name_method_exists(self):
        from provesid.chebi import ChebiSDF
        assert hasattr(ChebiSDF, "search_by_name")
        assert callable(ChebiSDF.search_by_name)

    def test_search_by_synonym_method_exists(self):
        from provesid.chebi import ChebiSDF
        assert hasattr(ChebiSDF, "search_by_synonym")
        assert callable(ChebiSDF.search_by_synonym)


# ─────────────────────────────────────────────────────────────────────────────
# ZeroPM method interface tests
# ─────────────────────────────────────────────────────────────────────────────

class TestZeroPMMethods:
    """Verify ZeroPM exposes the methods used by the Search resolvers."""

    def test_get_id_table_from_inchikey_method_exists(self):
        from provesid.zeropm import ZeroPM
        assert hasattr(ZeroPM, "get_id_table_from_inchikey"), (
            "ZeroPM must have get_id_table_from_inchikey() for _resolve_inchikey"
        )
        assert callable(ZeroPM.get_id_table_from_inchikey)

    def test_get_id_table_from_inchi_method_exists(self):
        from provesid.zeropm import ZeroPM
        assert hasattr(ZeroPM, "get_id_table_from_inchi"), (
            "ZeroPM must have get_id_table_from_inchi() for _resolve_inchi"
        )
        assert callable(ZeroPM.get_id_table_from_inchi)

    def test_get_id_table_from_cas_method_exists(self):
        from provesid.zeropm import ZeroPM
        assert hasattr(ZeroPM, "get_id_table_from_cas")
        assert callable(ZeroPM.get_id_table_from_cas)

    def test_get_id_table_from_name_method_exists(self):
        from provesid.zeropm import ZeroPM
        assert hasattr(ZeroPM, "get_id_table_from_name")
        assert callable(ZeroPM.get_id_table_from_name)

    def test_get_cas_from_smiles_method_exists(self):
        from provesid.zeropm import ZeroPM
        assert hasattr(ZeroPM, "get_cas_from_smiles")
        assert callable(ZeroPM.get_cas_from_smiles)

    def test_query_similar_name_method_exists(self):
        from provesid.zeropm import ZeroPM
        assert hasattr(ZeroPM, "query_similar_name"), (
            "ZeroPM must have query_similar_name() for fuzzy name search"
        )
        assert callable(ZeroPM.query_similar_name)

    def test_get_cas_from_formula_method_exists(self):
        from provesid.zeropm import ZeroPM
        assert hasattr(ZeroPM, "get_cas_from_formula")
        assert callable(ZeroPM.get_cas_from_formula)


# ─────────────────────────────────────────────────────────────────────────────
# Mock-based return-type tests (ensure methods return expected shapes)
# ─────────────────────────────────────────────────────────────────────────────

class TestPubChemIDReturnTypes:
    """Test PubChemID method return types with a mocked SQLite connection."""

    def _make_mock_pubchem(self, rows, columns):
        """Create a PubChemID with mocked _conn."""
        from provesid.pubchem import PubChemID

        mock_conn = MagicMock(spec=sqlite3.Connection)
        mock_cursor = MagicMock()
        mock_cursor.description = [(col,) for col in columns]
        mock_cursor.fetchone.return_value = rows[0] if rows else None
        mock_cursor.fetchall.return_value = rows
        mock_conn.execute.return_value = mock_cursor

        obj = object.__new__(PubChemID)
        obj.conn = mock_conn
        return obj

    def test_get_by_inchikey_returns_dict_or_none(self):
        from provesid.pubchem import PubChemID

        cols = ["cid", "inchikey", "inchi", "smiles", "mf", "mw", "cmpdname", "iupacname",
                "cas_numbers", "synonyms"]
        row = (1, "BSYNRYMUTXBXSQ-UHFFFAOYSA-N", "InChI=1S/...", "CC(=O)O",
               "C9H8O4", 180.16, "Aspirin", "2-(acetyloxy)benzoic acid",
               "50-78-2", "Aspirin")
        obj = self._make_mock_pubchem([row], cols)
        result = obj.get_by_inchikey("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        # Accept dict or None — just check it's callable and doesn't crash
        assert result is None or isinstance(result, dict)

    def test_get_by_inchi_returns_dict_or_none(self):
        from provesid.pubchem import PubChemID

        cols = ["cid", "inchi"]
        obj = self._make_mock_pubchem([], cols)
        result = obj.get_by_inchi("InChI=1S/C9H8O4")
        assert result is None or isinstance(result, dict)

    def test_search_by_formula_returns_list(self):
        from provesid.pubchem import PubChemID

        cols = ["cid", "mf", "smiles"]
        obj = self._make_mock_pubchem([], cols)
        result = obj.search_by_formula("C9H8O4")
        assert isinstance(result, list)


class TestCompToxIDReturnTypes:
    """Test CompToxID method return types with a mocked SQLite connection."""

    def _make_mock_comptox(self, rows, columns):
        from provesid.comptox import CompToxID

        mock_conn = MagicMock(spec=sqlite3.Connection)
        mock_cursor = MagicMock()
        mock_cursor.description = [(col,) for col in columns]
        mock_cursor.fetchone.return_value = rows[0] if rows else None
        mock_cursor.fetchall.return_value = rows
        mock_conn.execute.return_value = mock_cursor

        obj = object.__new__(CompToxID)
        obj.conn = mock_conn
        return obj

    def test_get_by_dtxsid_returns_dict_or_none(self):
        from provesid.comptox import CompToxID

        cols = ["DTXSID", "CASRN", "PREFERRED_NAME", "SMILES", "INCHIKEY"]
        obj = self._make_mock_comptox([], cols)
        result = obj.get_by_dtxsid("DTXSID2021735")
        assert result is None or isinstance(result, dict)

    def test_get_by_inchikey_returns_dict_or_none(self):
        from provesid.comptox import CompToxID

        cols = ["DTXSID", "INCHIKEY", "SMILES"]
        obj = self._make_mock_comptox([], cols)
        result = obj.get_by_inchikey("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        assert result is None or isinstance(result, dict)

    def test_search_by_formula_returns_list(self):
        from provesid.comptox import CompToxID

        cols = ["DTXSID", "MOLECULAR_FORMULA", "SMILES"]
        obj = self._make_mock_comptox([], cols)
        result = obj.search_by_formula("C9H8O4")
        assert isinstance(result, list)


# ─────────────────────────────────────────────────────────────────────────────
# Search package-level export
# ─────────────────────────────────────────────────────────────────────────────

class TestPackageExports:
    """Verify the Search class and utilities are exported from provesid."""

    def test_search_exported(self):
        import provesid
        assert hasattr(provesid, "Search")
        assert provesid.Search is not None

    def test_normalize_structure_exported(self):
        import provesid
        assert hasattr(provesid, "normalize_structure")
        assert callable(provesid.normalize_structure)

    def test_strip_salts_exported(self):
        import provesid
        assert hasattr(provesid, "strip_salts")
        assert callable(provesid.strip_salts)

    def test_output_columns_exported(self):
        import provesid
        assert hasattr(provesid, "OUTPUT_COLUMNS")
        assert isinstance(provesid.OUTPUT_COLUMNS, list)
        assert "query" in provesid.OUTPUT_COLUMNS
        assert "confidence" in provesid.OUTPUT_COLUMNS

    def test_search_class_instantiable(self):
        from provesid import Search
        s = Search("cas", show_progress=False)
        assert s.identifier_type == "cas"

    def test_all_identifier_types_accepted(self):
        from provesid import Search
        for id_type in ["cas", "name", "smiles", "inchi", "inchikey", "dtxsid", "formula"]:
            s = Search(id_type, show_progress=False)
            assert s.identifier_type == id_type
