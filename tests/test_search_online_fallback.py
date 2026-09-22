"""Tests for ``Search(online_fallback=...)``: the online rows of the source table.

No socket is opened.  ``provesid.search.PubChemAPI`` and
``provesid.search.NCIChemicalIdentifierResolver`` are patched with stubs that
answer like the real services and record what they were asked, so the tests
can check both what came back and how much network a run would have used.
"""

import logging

import pandas as pd
import pytest

from provesid.http import NotFoundError, ServiceError
from provesid.pubchem import PubChemAPI, PubChemNotFoundError
from provesid.search import Search
from provesid.sources import LOOKUPS, Query
from provesid.tools import candidate_from_cactus, candidate_from_pubchem_online

_ASPIRIN_SMILES = "CC(=O)OC1=CC=CC=C1C(=O)O"
_ASPIRIN_IK = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
_ASPIRIN_INCHI = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"

_ASPIRIN_ROW = {
    "CID": 2244, "Title": "Aspirin", "IUPACName": "2-acetyloxybenzoic acid",
    "MolecularFormula": "C9H8O4", "SMILES": _ASPIRIN_SMILES, "InChI": _ASPIRIN_INCHI,
    "InChIKey": _ASPIRIN_IK, "MolecularWeight": "180.16",
}


class _Empty:
    """An offline source that holds nothing."""

    def __getattr__(self, method):
        return lambda *args, **kwargs: None


class _OfflineAspirin:
    """An offline PubChemID that knows aspirin by every identifier."""

    def __getattr__(self, method):
        row = {"cid": 2244, "smiles": _ASPIRIN_SMILES, "inchikey": _ASPIRIN_IK,
               "cmpdname": "aspirin", "cas_numbers": ["50-78-2"]}
        return lambda *args, **kwargs: [row] if method.startswith("search_") else row


class _PubChemOnline:
    """Answers like PUG-REST: aspirin for aspirin's identifiers, 404 otherwise."""

    instances = []

    def __init__(self, fail=False, held=False):
        self.calls = []
        self.fail = fail
        self.held = held
        _PubChemOnline.instances.append(self)

    def _cids(self, method, value, known):
        self.calls.append((method, value))
        if self.held:
            # What the circuit breaker raises for a host a Retry-After holds.
            exc = ServiceError("pubchem.ncbi.nlm.nih.gov asked for no requests until ...")
            exc.held_until = 1e12
            raise exc
        if self.fail:
            raise ServiceError("PUGREST.ServerBusy")
        if value in known:
            return [2244]
        raise PubChemNotFoundError(f"no CID for {value}")

    def get_cids_by_name(self, value, name_type="word"):
        assert name_type == "complete"
        return self._cids("name", value, {"50-78-2", "aspirin", "DTXSID5020108"})

    def get_cids_by_inchikey(self, value):
        return self._cids("inchikey", value, {_ASPIRIN_IK})

    def get_cids_by_smiles(self, value):
        return self._cids("smiles", value, {_ASPIRIN_SMILES})

    def get_cids_by_inchi(self, value):
        return self._cids("inchi", value, {_ASPIRIN_INCHI})

    def get_properties_for_cids(self, cids, properties):
        self.calls.append(("properties", tuple(cids)))
        return [dict(_ASPIRIN_ROW) for cid in cids if cid == 2244]

    def get_compound_synonyms(self, cid):
        self.calls.append(("synonyms", cid))
        return ["aspirin", "Acetylsalicylic acid", "50-78-2", "DTXSID5020108"]


class _Cactus:
    """Answers like CACTUS: aspirin for aspirin's identifiers, 404 otherwise."""

    instances = []

    def __init__(self):
        self.calls = []
        _Cactus.instances.append(self)

    def resolve(self, identifier, representation):
        self.calls.append((identifier, representation))
        if identifier not in {"50-78-2", "aspirin", _ASPIRIN_IK, _ASPIRIN_SMILES, _ASPIRIN_INCHI}:
            raise NotFoundError(f"{identifier} not found")
        if representation == "smiles":
            return "CC(=O)Oc1ccccc1C(O)=O"
        return "Aspirin\n2-acetyloxybenzoic acid\n50-78-2"


@pytest.fixture
def online(monkeypatch):
    """Patch both online clients; yields the classes so tests can read calls."""
    _PubChemOnline.instances.clear()
    _Cactus.instances.clear()
    monkeypatch.setattr("provesid.search.PubChemAPI", _PubChemOnline)
    monkeypatch.setattr("provesid.search.NCIChemicalIdentifierResolver", _Cactus)
    return _PubChemOnline, _Cactus


def _search(identifier_type="cas", offline=None, **kwargs):
    """A Search on one offline source, which misses unless told otherwise."""
    return Search(
        identifier_type,
        show_progress=False,
        pubchem=offline if offline is not None else _Empty(),
        **kwargs,
    )


# ─────────────────────────────────────────────────────────────────────────────
# When the network is used
# ─────────────────────────────────────────────────────────────────────────────

class TestWhenOnline:
    def test_off_by_default_builds_no_online_client(self, online):
        df = _search().search("50-78-2")
        assert online[0].instances == [] and online[1].instances == []
        assert df.iloc[0]["confidence"] == 0.0
        assert df.attrs["online_fallbacks"] == 0
        assert df.attrs["online_resolved"] == 0

    def test_off_by_default_leaves_source_details_offline(self, online):
        df = _search().search("50-78-2")
        assert "PubChem (online)" not in df.iloc[0]["source_details"]

    def test_offline_answer_is_not_retried(self, online):
        df = _search(offline=_OfflineAspirin(), online_fallback=True).search("50-78-2")
        assert df.iloc[0]["InChIKey"] == _ASPIRIN_IK
        assert online[0].instances == []
        assert df.attrs["online_fallbacks"] == 0

    def test_offline_miss_is_asked_online(self, online):
        df = _search(online_fallback=True).search("50-78-2")
        row = df.iloc[0]
        assert row["InChIKey"] == _ASPIRIN_IK
        assert row["CASRN"] == "50-78-2"
        assert row["source"] == "PubChem (online)"
        assert row["source_details"]["PubChem (online)"]["found"] is True
        assert row["source_details"]["CACTUS"]["found"] is True
        assert row["source_details"]["PubChemID"]["found"] is False
        assert row["n_source_support"] == 2
        assert row["match_method"] == "exact_cas"
        assert df.attrs["online_fallbacks"] == 1
        assert df.attrs["online_resolved"] == 1

    def test_only_the_misses_go_online(self, online):
        class _KnowsOne(_Empty):
            def get_by_cas(self, cas):
                if cas == "64-17-5":
                    return {"cid": 702, "smiles": "CCO", "cmpdname": "ethanol"}
                return None

        df = _search(offline=_KnowsOne(), online_fallback=True).search(
            ["64-17-5", "50-78-2", "0000-00-0"]
        )
        assert list(df["source"].iloc[:2]) == ["PubChemID", "PubChem (online)"]
        assert pd.isna(df["source"].iloc[2])
        asked = [value for method, value in online[0].instances[0].calls if method == "name"]
        assert asked == ["50-78-2", "0000-00-0"]
        assert df.attrs["online_fallbacks"] == 2
        assert df.attrs["online_resolved"] == 1

    def test_clients_are_built_once_per_search_instance(self, online):
        s = _search(online_fallback=True)
        s.search(["50-78-2", "aspirin-cas-miss"])
        s.search("50-78-2")
        assert len(online[0].instances) == 1 and len(online[1].instances) == 1

    def test_counts_are_per_call(self, online):
        s = _search(online_fallback=True)
        s.search(["50-78-2", "58-08-2"])
        df = s.search("50-78-2")
        assert df.attrs["online_fallbacks"] == 1

    def test_formula_is_never_asked_online(self, online):
        df = _search("formula", online_fallback=True).search("C9H8O4")
        assert online[0].instances == []
        assert df.attrs["online_fallbacks"] == 0

    def test_fallback_is_logged_at_debug(self, online, caplog):
        with caplog.at_level(logging.DEBUG, logger="provesid.search"):
            _search(online_fallback=True).search("50-78-2")
        debug = [r.getMessage() for r in caplog.records if r.levelno == logging.DEBUG]
        assert any("No offline source answered cas '50-78-2'" in m for m in debug)
        assert any("PubChem (online) 1, CACTUS 1" in m for m in debug)


# ─────────────────────────────────────────────────────────────────────────────
# Every identifier type
# ─────────────────────────────────────────────────────────────────────────────

class TestIdentifierTypes:
    @pytest.mark.parametrize("identifier_type, query, expected_pubchem, cactus_asked", [
        ("cas", "50-78-2", ("name", "50-78-2"), True),
        ("name", "aspirin", ("name", "aspirin"), True),
        ("smiles", _ASPIRIN_SMILES, ("smiles", _ASPIRIN_SMILES), True),
        ("inchi", _ASPIRIN_INCHI, ("inchi", _ASPIRIN_INCHI), True),
        ("inchikey", _ASPIRIN_IK, ("inchikey", _ASPIRIN_IK), True),
        ("dtxsid", "DTXSID5020108", ("name", "DTXSID5020108"), False),
    ])
    def test_query_is_asked_as_its_own_type(
        self, online, identifier_type, query, expected_pubchem, cactus_asked
    ):
        df = _search(identifier_type, online_fallback=True).search(query)
        assert df.iloc[0]["InChIKey"] == _ASPIRIN_IK
        assert online[0].instances[0].calls[0] == expected_pubchem
        assert bool(online[1].instances[0].calls) is cactus_asked

    def test_name_hits_are_scored_against_the_query(self, online):
        df = _search("name", online_fallback=True).search("aspirin")
        row = df.iloc[0]
        assert row["match_method"] == "exact_name"
        assert row["source_match_scores"]
        assert row["confidence"] > 0.5

    def test_name_takes_top_k_cids(self, online):
        class _ManyCids(_PubChemOnline):
            def get_cids_by_name(self, value, name_type="word"):
                self.calls.append(("name", value))
                return [2244, 1, 2, 3, 4, 5, 6]

        online_api = _ManyCids()
        s = _search("name", online_fallback=True, top_k_per_source=3)
        s._online_clients_built = True
        s._clients["pubchem_online"] = online_api
        s.search("aspirin")
        assert ("properties", (2244, 1, 2)) in online_api.calls


# ─────────────────────────────────────────────────────────────────────────────
# Failures
# ─────────────────────────────────────────────────────────────────────────────

class TestFailures:
    def test_not_found_is_a_miss_not_a_failure(self, online, caplog):
        with caplog.at_level(logging.WARNING, logger="provesid.search"):
            df = _search(online_fallback=True).search("0000-00-0")
        assert df.iloc[0]["confidence"] == 0.0
        assert not [r for r in caplog.records if "lookup failed" in r.getMessage()]
        assert df.attrs["online_fallbacks"] == 1
        assert df.attrs["online_resolved"] == 0

    def test_a_failing_service_costs_its_own_vote(self, online, monkeypatch, caplog):
        monkeypatch.setattr(
            "provesid.search.PubChemAPI", lambda: _PubChemOnline(fail=True)
        )
        with caplog.at_level(logging.WARNING, logger="provesid.search"):
            df = _search(online_fallback=True).search("50-78-2")
        row = df.iloc[0]
        assert row["source"] == "CACTUS"
        assert row["n_source_support"] == 1
        assert any(
            "PubChem (online) cas lookup failed for '50-78-2'" in r.getMessage()
            for r in caplog.records
        )


    def test_a_held_host_is_not_warned_about_per_query(self, online, monkeypatch, caplog):
        """
        The transport warned once, when PubChem's Retry-After was recorded; a
        batch refused by the circuit breaker must not repeat that per query.
        """
        monkeypatch.setattr(
            "provesid.search.PubChemAPI", lambda: _PubChemOnline(held=True)
        )
        with caplog.at_level(logging.DEBUG, logger="provesid.search"):
            df = _search(online_fallback=True).search("50-78-2")
        assert df.iloc[0]["source"] == "CACTUS"
        failed = [r for r in caplog.records if "lookup failed" in r.getMessage()]
        assert failed, "the refusal was not logged at all"
        assert all(r.levelno == logging.DEBUG for r in failed)

# ─────────────────────────────────────────────────────────────────────────────
# The lookups and adapters, on their own
# ─────────────────────────────────────────────────────────────────────────────

class TestLookups:
    def test_cid_zero_is_a_miss(self):
        class _Zero(_PubChemOnline):
            def get_cids_by_inchikey(self, value):
                return [0]

        assert LOOKUPS["inchikey"]["pubchem_online"](_Zero(), Query(_ASPIRIN_IK)) == []

    def test_cid_without_properties_is_dropped(self):
        class _Bare(_PubChemOnline):
            def get_properties_for_cids(self, cids, properties):
                return [{"CID": cid} for cid in cids]

        assert LOOKUPS["cas"]["pubchem_online"](_Bare(), Query("50-78-2")) == []

    def test_cactus_takes_the_first_of_several_structures(self):
        class _Ambiguous(_Cactus):
            def resolve(self, identifier, representation):
                if representation == "smiles":
                    return "CCO\nCC(=O)O\n"
                return "ethanol"

        [cand] = LOOKUPS["name"]["cactus"](_Ambiguous(), Query("x"))
        assert cand["SMILES"] == "CCO"

    def test_cactus_without_names_keeps_the_structure(self):
        class _Nameless(_Cactus):
            def resolve(self, identifier, representation):
                if representation == "names":
                    raise NotFoundError("no names")
                return "CCO"

        [cand] = LOOKUPS["cas"]["cactus"](_Nameless(), Query("64-17-5"))
        assert cand["SMILES"] == "CCO" and cand["name"] is None
        assert cand["InChIKey"] == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"

    def test_pubchem_adapter_takes_cas_from_synonyms(self):
        cand = candidate_from_pubchem_online(_ASPIRIN_ROW, ["aspirin", "50-78-2"])
        assert cand["source"] == "PubChem (online)"
        assert cand["CAS_candidates"] == ["50-78-2"]
        assert cand["molecular_mass"] == pytest.approx(180.16)

    def test_cactus_adapter_derives_the_inchikey(self):
        cand = candidate_from_cactus("CC(=O)Oc1ccccc1C(O)=O", ["Aspirin", "50-78-2"])
        assert cand["InChIKey"] == _ASPIRIN_IK
        assert cand["CAS_candidates"] == ["50-78-2"]


class TestGetCidsByInchi:
    def test_posts_the_inchi_in_the_body(self, monkeypatch):
        sent = {}

        class _Response:
            def json(self):
                return {"IdentifierList": {"CID": [2244]}}

        def fake_request(self, url, method="GET", data=None, **kwargs):
            sent.update(url=url, method=method, data=data)
            return _Response()

        monkeypatch.setattr(PubChemAPI, "_make_request", fake_request)
        assert PubChemAPI(use_cache=False).get_cids_by_inchi(_ASPIRIN_INCHI) == [2244]
        assert sent["method"] == "POST"
        assert sent["url"].endswith("/compound/inchi/cids/JSON")
        assert sent["data"] == {"inchi": _ASPIRIN_INCHI}


class TestEnrich:
    def test_enrich_carries_the_online_counts(self, online):
        frame = pd.DataFrame({"CAS": ["50-78-2", "50-78-2", "0000-00-0"]})
        out = _search(online_fallback=True).enrich(frame, "CAS")
        assert out.attrs["online_fallbacks"] == 2
        assert out.attrs["online_resolved"] == 1
        assert list(out["provesid_InChIKey"].iloc[:2]) == [_ASPIRIN_IK, _ASPIRIN_IK]
