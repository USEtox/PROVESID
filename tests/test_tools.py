"""Unit tests for offline CAS resolution helpers in provesid.tools."""

import provesid.tools as tools
import pandas as pd


class _ChebiWithResult:
    def __init__(self, rows):
        self.rows = rows

    def search_by_cas(self, cas):
        return self.rows

    def search_by_name(self, name, exact=True):
        return self.rows

    def search_by_synonym(self, name, exact=True):
        return self.rows

    def search_by_inchikey(self, inchikey):
        return self.rows[0] if self.rows else None


class _CompToxWithResult:
    def __init__(self, row):
        self.row = row

    def get_by_casrn(self, cas):
        return self.row

    def get_by_name(self, name):
        return self.row

    def search_by_name(self, name, exact=False, limit=10):
        return [self.row] if self.row else []

    def get_by_smiles(self, smiles):
        return self.row


class _PubChemWithResult:
    def __init__(self, row):
        self.row = row

    def get_by_cas(self, cas):
        return self.row

    def search_by_name(self, name, exact=False, limit=10):
        return [self.row] if self.row else []

    def get_by_smiles(self, smiles):
        return self.row


class _ZeroPMWithResult:
    def __init__(self, table):
        self.table = table

    def get_id_table_from_cas(self, cas):
        return self.table

    def get_id_table_from_name(self, name):
        return self.table

    def get_cas_from_smiles(self, smiles):
        if self.table is None or self.table.empty or "cas" not in self.table.columns:
            return None
        return self.table["cas"].dropna().tolist()


class _ChEMBLWithResult:
    def __init__(self, compound=None, props=None):
        self.compound = compound
        self.props = props or {}

    def search_by_smiles(self, smiles):
        return self.compound

    def search_by_name(self, name, limit=100):
        return [self.compound] if self.compound else []

    def get_properties(self, molregno):
        return self.props


def _minimal_row(cas):
    return {
        "CASRN": cas,
        "name": None,
        "IUPAC_name": None,
        "molecular_formula": None,
        "SMILES": None,
        "canonical_smiles": None,
        "InChI": None,
        "InChIKey": None,
        "DTXSID": None,
        "molecular_mass": None,
        "Synonyms": None,
        "foundby": "CASRN",
        "source": None,
    }


def test_ids_from_cas_prefers_chebi_when_smiles_present(monkeypatch):
    """ChEBI should short-circuit lower-priority sources when it provides SMILES."""
    monkeypatch.setattr(tools, "_smiles_to_canonical_and_mass", lambda smiles: ("CHEBI_CANON", 180.1))

    chebi = _ChebiWithResult([
        {
            "ChEBI NAME": "Aspirin",
            "FORMULA": "C9H8O4",
            "SMILES": "CC(=O)Oc1ccccc1C(=O)O",
            "INCHI": "InChI=1S/C9H8O4/...",
            "INCHIKEY": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
            "SYNONYM": "Acetylsalicylic acid",
        }
    ])

    class _CompToxMustNotRun:
        def get_by_casrn(self, cas):
            raise AssertionError("CompTox should not be queried when ChEBI already has SMILES")

    result = tools.ids_from_CAS("50-78-2", chebi=chebi, comptox=_CompToxMustNotRun())

    assert result["source"] == "ChEBI"
    assert result["name"] == "Aspirin"
    assert result["canonical_smiles"] == "CHEBI_CANON"
    assert result["molecular_mass"] == 180.1


def test_ids_from_cas_falls_back_to_comptox(monkeypatch):
    """CompTox should be used when ChEBI has no match."""
    monkeypatch.setattr(tools, "_smiles_to_canonical_and_mass", lambda smiles: ("COMPTOX_CANON", 46.07))

    chebi = _ChebiWithResult([])
    comptox = _CompToxWithResult(
        {
            "PREFERRED_NAME": "Ethanol",
            "IUPAC_NAME": "ethanol",
            "MOLECULAR_FORMULA": "C2H6O",
            "SMILES": "CCO",
            "INCHIKEY": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            "DTXSID": "DTXSID7020001",
            "identifiers": ["Ethyl alcohol", "EtOH"],
            "AVERAGE_MASS": 46.07,
        }
    )

    result = tools.ids_from_CAS("64-17-5", chebi=chebi, comptox=comptox)

    assert result["source"] == "CompTox"
    assert result["DTXSID"] == "DTXSID7020001"
    assert result["Synonyms"] == "Ethyl alcohol; EtOH"
    assert result["canonical_smiles"] == "COMPTOX_CANON"


def test_ids_from_cas_uses_pubchem_then_chembl_enrichment(monkeypatch):
    """PubChem can provide base identifiers and ChEMBL can enrich missing mass."""
    monkeypatch.setattr(tools, "_smiles_to_canonical_and_mass", lambda smiles: ("PUBCHEM_CANON", 999.0))

    pubchem = _PubChemWithResult(
        {
            "cmpdname": "Acetophenone",
            "iupacname": "1-phenylethan-1-one",
            "mf": "C8H8O",
            "smiles": "CC(=O)c1ccccc1",
            "inchi": "InChI=1S/C8H8O/...",
            "inchikey": "KWYUFKZDYYNOTN-UHFFFAOYSA-N",
            "synonyms": ["Acetophenone"],
            "mw": None,
        }
    )
    chembl = _ChEMBLWithResult(
        compound={"molregno": 123, "pref_name": "ACETOPHENONE"},
        props={"mw_freebase": 120.15},
    )

    result = tools.ids_from_CAS(
        "98-86-2",
        chebi=_ChebiWithResult([]),
        comptox=_CompToxWithResult(None),
        pubchem=pubchem,
        chembl=chembl,
    )

    assert result["source"] == "PubChemID"
    assert result["name"] == "Acetophenone"
    assert result["canonical_smiles"] == "PUBCHEM_CANON"
    assert result["molecular_mass"] == 120.15


def test_casrn_to_compounds_delegates_to_ids_from_cas(monkeypatch):
    """Batch wrapper should initialize clients once and call ids_from_CAS per CAS."""
    created = {"chebi": 0, "comptox": 0, "pubchem": 0, "zeropm": 0, "chembl": 0}

    class _StubChebi:
        def __init__(self):
            created["chebi"] += 1

    class _StubCompTox:
        def __init__(self):
            created["comptox"] += 1

    class _StubPubChem:
        def __init__(self):
            created["pubchem"] += 1

    class _StubZeroPM:
        def __init__(self):
            created["zeropm"] += 1

    class _StubChembl:
        def __init__(self):
            created["chembl"] += 1

    monkeypatch.setattr(tools, "ChebiSDF", _StubChebi)
    monkeypatch.setattr(tools, "CompToxID", _StubCompTox)
    monkeypatch.setattr(tools, "PubChemID", _StubPubChem)
    monkeypatch.setattr(tools, "ZeroPM", _StubZeroPM)
    monkeypatch.setattr(tools, "CheMBL", _StubChembl)

    calls = []

    def _fake_ids_from_cas(cas, **kwargs):
        calls.append((cas, kwargs))
        row = _minimal_row(cas)
        row["name"] = f"name-{cas}"
        row["source"] = "stub"
        return row

    monkeypatch.setattr(tools, "ids_from_CAS", _fake_ids_from_cas)

    df = tools.casrn_to_compounds(["50-78-2", "64-17-5"], show_progress=False)

    assert list(df["CASRN"]) == ["50-78-2", "64-17-5"]
    assert list(df["name"]) == ["name-50-78-2", "name-64-17-5"]

    assert created == {"chebi": 1, "comptox": 1, "pubchem": 1, "zeropm": 1, "chembl": 1}
    assert len(calls) == 2
    assert calls[0][1]["chebi"] is not None
    assert calls[0][1]["comptox"] is not None
    assert calls[0][1]["pubchem"] is not None
    assert calls[0][1]["zeropm"] is not None
    assert calls[0][1]["chembl"] is not None


def test_ids_from_name_uses_consensus_then_priority(monkeypatch):
    """Name lookup should build consensus but still populate using CAS priority."""
    monkeypatch.setattr(tools, "_smiles_to_canonical_and_mass", lambda smiles: ("CANON_ETHANOL", 46.07))

    chebi = _ChebiWithResult(
        [
            {
                "ChEBI NAME": "ethanol",
                "FORMULA": "C2H6O",
                "SMILES": "CCO",
                "INCHI": "InChI=1S/C2H6O/...",
                "INCHIKEY": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                "SYNONYM": "ethyl alcohol",
            }
        ]
    )
    comptox = _CompToxWithResult(
        {
            "PREFERRED_NAME": "Ethanol",
            "IUPAC_NAME": "ethanol",
            "CASRN": "64-17-5",
            "MOLECULAR_FORMULA": "C2H6O",
            "SMILES": "CCO",
            "INCHIKEY": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            "DTXSID": "DTXSID7020001",
            "identifiers": ["ethyl alcohol"],
            "AVERAGE_MASS": 46.07,
        }
    )
    pubchem = _PubChemWithResult(
        {
            "cmpdname": "Ethanol",
            "iupacname": "ethanol",
            "cas_numbers": ["64-17-5"],
            "mf": "C2H6O",
            "smiles": "CCO",
            "inchi": "InChI=1S/C2H6O/...",
            "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            "synonyms": ["ethyl alcohol"],
            "mw": 46.07,
        }
    )
    zeropm = _ZeroPMWithResult(
        pd.DataFrame(
            [
                {
                    "name": "ethanol",
                    "query_id": 1,
                    "inchi_id": 1,
                    "rank": 1,
                    "inchi": "InChI=1S/C2H6O/...",
                    "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                    "cas": "64-17-5",
                    "sources": "ZeroPM",
                }
            ]
        )
    )

    result = tools.ids_from_name("ethanol", chebi=chebi, comptox=comptox, pubchem=pubchem, zeropm=zeropm)

    assert result["foundby"] == "name"
    assert result["canonical_smiles"] == "CANON_ETHANOL"
    assert result["source"] == "ChEBI"
    assert result["CASRN"] == "64-17-5"
    assert result["match_score"] > 0.5
    assert "ChEBI" in result["source_match_scores"]


def test_ids_from_smiles_consensus_includes_score(monkeypatch):
    """SMILES lookup should compare records and return a consensus score."""
    monkeypatch.setattr(tools, "_smiles_to_canonical_and_mass", lambda smiles: ("CANON_INPUT", 46.07))
    monkeypatch.setattr(tools, "_inchikey_from_smiles", lambda smiles: "LFQSCWFLJHTTHZ-UHFFFAOYSA-N")
    monkeypatch.setattr(tools, "_inchi_to_smiles", lambda inchi: "CCO")

    chebi = _ChebiWithResult(
        [
            {
                "ChEBI NAME": "ethanol",
                "FORMULA": "C2H6O",
                "SMILES": "CCO",
                "INCHI": "InChI=1S/C2H6O/...",
                "INCHIKEY": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                "SYNONYM": "ethyl alcohol",
            }
        ]
    )
    comptox = _CompToxWithResult(
        {
            "PREFERRED_NAME": "Ethanol",
            "IUPAC_NAME": "ethanol",
            "CASRN": "64-17-5",
            "MOLECULAR_FORMULA": "C2H6O",
            "SMILES": "CCO",
            "INCHIKEY": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            "DTXSID": "DTXSID7020001",
            "identifiers": ["ethyl alcohol"],
            "AVERAGE_MASS": 46.07,
        }
    )
    pubchem = _PubChemWithResult(
        {
            "cmpdname": "Ethanol",
            "iupacname": "ethanol",
            "cas_numbers": ["64-17-5"],
            "mf": "C2H6O",
            "smiles": "CCO",
            "inchi": "InChI=1S/C2H6O/...",
            "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            "synonyms": ["ethyl alcohol"],
            "mw": 46.07,
        }
    )
    zeropm = _ZeroPMWithResult(
        pd.DataFrame(
            [
                {
                    "cas": "64-17-5",
                    "query_id": 1,
                    "inchi_id": 1,
                    "rank": 1,
                    "inchi": "InChI=1S/C2H6O/...",
                    "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                    "synonyms": "ethanol",
                    "sources": "ZeroPM",
                }
            ]
        )
    )

    result = tools.ids_from_SMILES("CCO", chebi=chebi, comptox=comptox, pubchem=pubchem, zeropm=zeropm)

    assert result["foundby"] == "SMILES"
    assert result["source"] == "ChEBI"
    assert result["CASRN"] == "64-17-5"
    assert result["canonical_smiles"] == "CANON_INPUT"
    assert result["match_score"] > 0.5
    assert "CompTox" in result["source_match_scores"]
