"""
Offline tests for the CAS Common Chemistry client.

``tests/test_cascommonchem.py`` skips itself entirely without a CAS API key, so
until this file existed the module had no coverage on a developer machine at
all --- and it is the module that gained the most in the move to the shared
transport, having had no pacing and no retry before it.

Nothing here touches the network or needs a real key: the key is only a header,
and ``requests.get`` is stubbed throughout.
"""

import logging

import pytest
import requests

from provesid.cascommonchem import (
    CASCommonChem,
    CASCommonChemError,
    CASCommonChemNotFoundError,
    CASCommonChemTimeoutError,
    _lookup_failed,
)


class _FakeResponse:
    """The slice of ``requests.Response`` the transport and parser touch."""

    def __init__(self, status_code, payload=None, text=""):
        self.status_code = status_code
        self._payload = payload
        self.text = text or str(payload or "")
        self.headers = {}

    def json(self):
        if self._payload is None:
            raise ValueError("no JSON body")
        return self._payload


@pytest.fixture
def cas():
    """A client with a stub key and no cache, so stubs are never bypassed."""
    client = CASCommonChem(api_key="test-key-not-a-real-one", use_cache=False)
    client._http.backoff = 0.0
    client._http.min_interval = 0.0
    return client


WATER = {
    "rn": "7732-18-5",
    "name": "Water",
    "molecularFormula": "H<sub>2</sub>O",
    "molecularMass": "18.02",
    "inchiKey": "InChIKey=XLYOFNOQVPJJNP-UHFFFAOYSA-N",
    "synonyms": ["Water", "Dihydrogen oxide"],
}


# --------------------------------------------------------------------------
# The happy path
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_a_successful_lookup_is_flattened_into_the_result(cas, monkeypatch):
    """Every key CAS sends lands at the top level, beside the status."""
    monkeypatch.setattr(
        requests, "get",
        lambda url, **kwargs: _FakeResponse(200, payload=WATER),
    )
    result = cas.cas_to_detail("7732-18-5")

    assert result["status"] == "Success"
    assert result["found"] is True
    assert result["rn"] == "7732-18-5"
    assert result["synonyms"] == ["Water", "Dihydrogen oxide"]


@pytest.mark.unit
def test_the_api_key_is_sent_on_every_request(cas, monkeypatch):
    """The key lives on the transport now, not in each method."""
    seen = {}

    def record(url, **kwargs):
        seen.update(kwargs.get("headers") or {})
        return _FakeResponse(200, payload=WATER)

    monkeypatch.setattr(requests, "get", record)
    cas.cas_to_detail("7732-18-5")

    assert seen["X-API-KEY"] == "test-key-not-a-real-one"


@pytest.mark.unit
def test_a_name_search_follows_through_to_the_detail_call(cas, monkeypatch):
    """A search gives a CAS number; the detail call gives the substance."""
    urls = []

    def route(url, **kwargs):
        urls.append(url)
        if "/search" in url:
            return _FakeResponse(200, payload={"count": 1, "results": [{"rn": "7732-18-5"}]})
        return _FakeResponse(200, payload=WATER)

    monkeypatch.setattr(requests, "get", route)
    result = cas.name_to_detail("water")

    assert result["found"] is True
    assert result["rn"] == "7732-18-5"
    assert len(urls) == 2


@pytest.mark.unit
def test_a_search_with_no_hits_is_not_found(cas, monkeypatch):
    """An empty result set is absence, reported in the status."""
    monkeypatch.setattr(
        requests, "get",
        lambda url, **kwargs: _FakeResponse(200, payload={"count": 0, "results": []}),
    )
    result = cas.name_to_detail("not-a-substance-12345")

    assert result["status"] == "Not found"
    assert result["found"] is False


# --------------------------------------------------------------------------
# smiles_to_detail: search by InChI, keep the exact match
# --------------------------------------------------------------------------

ETHANOL_INCHI = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
ETHANOL_RECORDS = {
    # The order CAS's search returned them in on 2026-09-23: the dimer first.
    "42845-45-4": {"rn": "42845-45-4", "name": "Ethanol, dimer",
                   "molecularFormula": "(C<sub>2</sub>H<sub>6</sub>O)<sub>2</sub>",
                   "inchi": ETHANOL_INCHI, "synonyms": ["Ethanol, dimer"]},
    "64-17-5": {"rn": "64-17-5", "name": "Ethanol",
                "molecularFormula": "C<sub>2</sub>H<sub>6</sub>O",
                "inchi": ETHANOL_INCHI, "synonyms": ["Ethanol", "Ethyl alcohol"]},
}


def _cas_service(records, searched=None):
    """Stub ``requests.get`` as CAS: ``/search`` lists ``records``, ``/detail`` serves one."""
    def route(url, **kwargs):
        if "/search" in url:
            if searched is not None:
                searched.append(requests.utils.unquote(url.split("?q=", 1)[1]))
            hits = [{"rn": rn, "name": r["name"]} for rn, r in records.items()]
            return _FakeResponse(200, payload={"count": len(hits), "results": hits})
        return _FakeResponse(200, payload=records[url.split("cas_rn=", 1)[1]])
    return route


@pytest.mark.unit
def test_a_smiles_is_searched_as_its_inchi(cas, monkeypatch):
    """CAS matches a SMILES only as its own string; an InChI whatever the spelling."""
    searched = []
    monkeypatch.setattr(requests, "get", _cas_service(ETHANOL_RECORDS, searched))
    for smiles in ["CCO", "OCC", "C(O)C"]:
        cas.smiles_to_detail(smiles)
    assert searched == [ETHANOL_INCHI] * 3


@pytest.mark.unit
def test_a_smiles_skips_the_oligomer_that_shares_its_inchi(cas, monkeypatch):
    monkeypatch.setattr(requests, "get", _cas_service(ETHANOL_RECORDS))
    result = cas.smiles_to_detail("CCO")
    assert result["found"] is True
    assert result["rn"] == "64-17-5"


@pytest.mark.unit
def test_a_smiles_prefers_the_record_with_more_synonyms(cas, monkeypatch, caplog):
    """Sodium chloride and rock salt share InChI and formula."""
    inchi = "InChI=1S/ClH.Na/h1H;/q;+1/p-1"
    records = {
        "14762-51-7": {"rn": "14762-51-7", "name": "Rock salt", "molecularFormula": "ClNa",
                       "inchi": inchi, "synonyms": ["Rock salt", "Halite"]},
        "7647-14-5": {"rn": "7647-14-5", "name": "Sodium chloride", "molecularFormula": "ClNa",
                      "inchi": inchi, "synonyms": ["Sodium chloride", "Salt", "NaCl"]},
    }
    monkeypatch.setattr(requests, "get", _cas_service(records))
    with caplog.at_level(logging.WARNING):
        result = cas.smiles_to_detail("[Na+].[Cl-]")
    assert result["rn"] == "7647-14-5"
    assert "14762-51-7" in caplog.text


@pytest.mark.unit
def test_a_smiles_whose_formula_cas_writes_differently_still_matches(cas, monkeypatch):
    """No record passes the formula test, so the InChI match alone decides."""
    records = {"64-17-5": dict(ETHANOL_RECORDS["64-17-5"], molecularFormula="C2H5OH")}
    monkeypatch.setattr(requests, "get", _cas_service(records))
    assert cas.smiles_to_detail("CCO")["rn"] == "64-17-5"


@pytest.mark.unit
def test_a_smiles_with_no_record_of_its_inchi_is_not_found(cas, monkeypatch):
    """A hit with another InChI (an isotopologue, say) is not the substance."""
    records = {"1516-08-1": {"rn": "1516-08-1", "name": "Ethanol-d6",
                             "molecularFormula": "C<sub>2</sub>D<sub>6</sub>O",
                             "inchi": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3/i1D3,2D2,3D",
                             "synonyms": []}}
    monkeypatch.setattr(requests, "get", _cas_service(records))
    result = cas.smiles_to_detail("CCO")
    assert result["status"] == "Not found"
    assert result["found"] is False


@pytest.mark.unit
@pytest.mark.parametrize("smiles", ["not-a-smiles", "", "   "])
def test_an_unreadable_smiles_is_refused_without_a_request(cas, monkeypatch, smiles):
    def no_network(url, **kwargs):
        raise AssertionError(f"request made: {url}")
    monkeypatch.setattr(requests, "get", no_network)
    result = cas.smiles_to_detail(smiles)
    assert result["status"] == "Invalid SMILES"
    assert result["found"] is False


@pytest.mark.unit
def test_a_failed_smiles_search_reports_the_failure(cas, monkeypatch):
    monkeypatch.setattr(requests, "get", lambda url, **kwargs: _FakeResponse(403))
    assert cas.smiles_to_detail("CCO")["status"] == "Unauthorized - Check API Key"


# --------------------------------------------------------------------------
# Failures are named, and named apart from each other
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_an_unknown_cas_number_is_not_found(cas, monkeypatch):
    """CAS answers an unknown registry number with a 404."""
    monkeypatch.setattr(
        requests, "get",
        lambda url, **kwargs: _FakeResponse(404, text="Not Found"),
    )
    result = cas.cas_to_detail("9999999-99-9")

    assert result["status"] == "Not Found"
    assert result["found"] is False


@pytest.mark.unit
def test_a_rejected_key_is_reported_as_such(cas, monkeypatch, caplog):
    """
    A 401 must not read as "no such substance": the substance may well exist
    and the caller needs to be told to fix their key.
    """
    monkeypatch.setattr(
        requests, "get",
        lambda url, **kwargs: _FakeResponse(401, text="Unauthorized"),
    )
    with caplog.at_level(logging.ERROR):
        result = cas.cas_to_detail("7732-18-5")

    assert result["status"] == "Unauthorized - Check API Key"
    assert result["found"] is False
    assert "authentication failed" in caplog.text


@pytest.mark.unit
def test_a_rejected_key_is_reported_on_a_name_search_too(cas, monkeypatch):
    """Both entry points read the key the same way."""
    monkeypatch.setattr(
        requests, "get",
        lambda url, **kwargs: _FakeResponse(401, text="Unauthorized"),
    )
    assert cas.name_to_detail("water")["status"] == "Unauthorized - Check API Key"


@pytest.mark.unit
def test_a_timeout_is_reported_as_a_timeout(cas, monkeypatch):
    """A timeout is transient, and says so rather than claiming absence."""
    monkeypatch.setattr(
        requests, "get",
        lambda url, **kwargs: (_ for _ in ()).throw(requests.Timeout("gone")),
    )
    result = cas.cas_to_detail("7732-18-5")

    assert result["status"] == "Timeout"
    assert result["found"] is False


@pytest.mark.unit
def test_a_rejected_key_is_never_retried(cas, monkeypatch):
    """Asking again with the same bad key cannot help."""
    attempts = []

    def unauthorized(url, **kwargs):
        attempts.append(url)
        return _FakeResponse(401, text="Unauthorized")

    monkeypatch.setattr(requests, "get", unauthorized)
    cas.cas_to_detail("7732-18-5")

    assert len(attempts) == 1


@pytest.mark.unit
def test_a_busy_service_is_retried(cas, monkeypatch):
    """
    The module had no retry at all before the shared transport; a single 503
    became a permanent "Internal Server Error" in the cache.
    """
    attempts = []

    def busy(url, **kwargs):
        attempts.append(url)
        return _FakeResponse(503, text="Service Unavailable")

    monkeypatch.setattr(requests, "get", busy)
    result = cas.cas_to_detail("7732-18-5")

    assert result["found"] is False
    assert len(attempts) == cas._http.max_retries + 1


@pytest.mark.unit
def test_a_transient_failure_that_clears_is_invisible(cas, monkeypatch):
    """The caller gets the substance, not an apology."""
    attempts = []

    def flaky(url, **kwargs):
        attempts.append(url)
        if len(attempts) == 1:
            return _FakeResponse(503, text="Service Unavailable")
        return _FakeResponse(200, payload=WATER)

    monkeypatch.setattr(requests, "get", flaky)
    result = cas.cas_to_detail("7732-18-5")

    assert result["found"] is True
    assert result["rn"] == "7732-18-5"
    assert len(attempts) == 2


@pytest.mark.unit
def test_no_requests_exception_reaches_the_caller(cas, monkeypatch):
    """The methods here have always promised a dict, never a raise."""
    monkeypatch.setattr(
        requests, "get",
        lambda url, **kwargs: (_ for _ in ()).throw(requests.ConnectionError("down")),
    )
    result = cas.cas_to_detail("7732-18-5")

    assert isinstance(result, dict)
    assert result["found"] is False


# --------------------------------------------------------------------------
# And none of those failures is remembered
# --------------------------------------------------------------------------

@pytest.mark.unit
@pytest.mark.parametrize("result,skipped", [
    ({"found": True, "rn": "7732-18-5"}, False),
    ({"found": False, "status": "Timeout"}, True),
    ({"found": False, "status": "Not Found"}, True),
    ({"status": "Success"}, True),
    (None, True),
])
def test_only_a_found_substance_is_cacheable(result, skipped):
    """Absence and failure look the same from here, so neither is stored."""
    assert _lookup_failed(result) is skipped


@pytest.mark.unit
def test_a_failed_lookup_is_retried_next_time(tmp_path, monkeypatch):
    """
    A throttled request must not become a permanent "no such CAS number".
    """
    import provesid.cache as cache_module
    from provesid.cache import CacheManager

    probe = CacheManager(cache_dir=str(tmp_path / "cas"), service_name="cas")
    monkeypatch.setitem(cache_module._service_caches, "cas", probe)

    client = CASCommonChem(api_key="test-key-not-a-real-one")
    client._http.backoff = 0.0
    client._http.min_interval = 0.0
    client._http.max_retries = 0

    attempts = []

    def flaky(url, **kwargs):
        attempts.append(url)
        if len(attempts) == 1:
            return _FakeResponse(503, text="Service Unavailable")
        return _FakeResponse(200, payload=WATER)

    monkeypatch.setattr(requests, "get", flaky)

    assert client.cas_to_detail("7732-18-5")["found"] is False
    assert client.cas_to_detail("7732-18-5")["found"] is True, \
        "the failure was served from the cache"
    assert len(attempts) == 2


# --------------------------------------------------------------------------
# The exception hierarchy
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_the_exceptions_catch_as_one_service_and_as_all_of_them():
    """A caller can catch CAS, or every web service in the package."""
    from provesid.http import NotFoundError, ServiceError, ServiceTimeoutError

    assert issubclass(CASCommonChemNotFoundError, CASCommonChemError)
    assert issubclass(CASCommonChemTimeoutError, CASCommonChemError)
    assert issubclass(CASCommonChemError, ServiceError)
    assert issubclass(CASCommonChemNotFoundError, NotFoundError)
    assert issubclass(CASCommonChemTimeoutError, ServiceTimeoutError)


# --------------------------------------------------------------------------
# The cache key
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_the_cache_key_is_stable_and_carries_the_schema_version():
    """
    An entry has to be findable in the next process, and retirable when its
    content turns out to be wrong.
    """
    first = CASCommonChem(api_key="test-key-not-a-real-one")
    second = CASCommonChem(api_key="test-key-not-a-real-one")

    assert first.__cache_key__() == second.__cache_key__()
    assert first.__cache_key__()[-1] == CASCommonChem.CACHE_SCHEMA_VERSION
    assert CASCommonChem.CACHE_SCHEMA_VERSION >= 2, \
        "version 1 entries cached failures as answers"


@pytest.mark.unit
def test_the_api_key_is_not_part_of_the_cache_key():
    """
    Two keys reach the same registry and get the same answer, so keying on the
    key would only mean a new key starts from an empty cache.
    """
    one = CASCommonChem(api_key="key-one")
    two = CASCommonChem(api_key="key-two")

    assert one.__cache_key__() == two.__cache_key__()
    assert "key-one" not in str(one.__cache_key__())
