"""
Tests that a failed PubChem fetch is reported, not disguised as "no data".

Offline tests: the transport (``_make_request``) is stubbed, so nothing here
touches the network. They guard the case that motivated the change — PUG-View
answers ``PUGVIEW.ServerBusy`` (HTTP 503) often enough that a client which
returns ``[]`` on failure will, sooner or later, cache "this compound has no
melting point" for a compound that has one.
"""

import pandas as pd
import pytest

from provesid.pubchem import PubChemAPI, PubChemError, PubChemNotFoundError
from provesid.pubchemview import (
    PubChemView,
    PubChemViewError,
    PubChemViewNotFoundError,
)


# --------------------------------------------------------------------------
# PUG-View
# --------------------------------------------------------------------------

@pytest.fixture
def view():
    """A PubChemView with caching disabled so stubs are not served from disk."""
    return PubChemView(use_cache=False)


@pytest.mark.unit
def test_extract_property_data_raises_on_transport_failure(view, monkeypatch):
    """An exhausted retry budget must surface, not become an empty list."""
    def busy(url):
        raise PubChemViewError("Request failed after 4 attempts: 503 ServerBusy")

    monkeypatch.setattr(view, "_make_request", busy)
    with pytest.raises(PubChemViewError):
        view.extract_property_data(2244, "Boiling Point")


@pytest.mark.unit
def test_extract_property_data_returns_empty_on_real_absence(view, monkeypatch):
    """A 404 means the compound genuinely has no such property."""
    def absent(url):
        raise PubChemViewNotFoundError(f"Resource not found: {url}")

    monkeypatch.setattr(view, "_make_request", absent)
    assert view.extract_property_data(2244, "Boiling Point") == []


@pytest.mark.unit
def test_property_table_raises_on_transport_failure(view, monkeypatch):
    """An empty table must never stand in for a failed request."""
    def busy(url):
        raise PubChemViewError("Request failed after 4 attempts: 503 ServerBusy")

    monkeypatch.setattr(view, "_make_request", busy)
    with pytest.raises(PubChemViewError):
        view.get_property_table(2244, "Melting Point")


@pytest.mark.unit
def test_property_table_returns_empty_frame_on_real_absence(view, monkeypatch):
    """A 404 yields an empty frame with the documented columns."""
    def absent(url):
        raise PubChemViewNotFoundError(f"Resource not found: {url}")

    monkeypatch.setattr(view, "_make_request", absent)
    table = view.get_property_table(2244, "Melting Point")
    assert isinstance(table, pd.DataFrame)
    assert table.empty
    assert list(table.columns) == PubChemView.PROPERTY_TABLE_COLUMNS


@pytest.mark.unit
def test_convenience_getter_does_not_cache_a_failure(view, monkeypatch):
    """get_melting_point must raise rather than cache [] for a busy server."""
    def busy(url):
        raise PubChemViewError("Request failed after 4 attempts: 503 ServerBusy")

    monkeypatch.setattr(view, "_make_request", busy)
    with pytest.raises(PubChemViewError):
        view.get_melting_point(2244)


# --------------------------------------------------------------------------
# PUG-REST
# --------------------------------------------------------------------------

@pytest.fixture
def api():
    """A PubChemAPI with caching disabled so stubs are not served from disk."""
    return PubChemAPI(use_cache=False)


@pytest.mark.unit
def test_synonyms_raise_on_transport_failure(api, monkeypatch):
    """A throttled synonym fetch must not look like "no synonyms"."""
    def throttled(url, **kwargs):
        raise PubChemError("HTTP error 429: Too Many Requests")

    monkeypatch.setattr(api, "_make_request", throttled)
    with pytest.raises(PubChemError):
        api.get_compound_synonyms(2244)


@pytest.mark.unit
def test_synonyms_return_empty_on_real_absence(api, monkeypatch):
    """A 404 on the synonym endpoint is a genuine empty result."""
    def absent(url, **kwargs):
        raise PubChemNotFoundError("Resource not found")

    monkeypatch.setattr(api, "_make_request", absent)
    assert api.get_compound_synonyms(2244) == []


@pytest.mark.unit
def test_partial_property_result_records_the_synonym_error(api, monkeypatch):
    """Properties still come back, but the synonym failure is recorded."""
    monkeypatch.setattr(
        api, "_make_request",
        lambda url, **kwargs: _stub_response(
            {"PropertyTable": {"Properties": [{"CID": 2244, "MolecularWeight": "180.16"}]}}
        ),
    )
    monkeypatch.setattr(
        api, "get_compound_synonyms",
        lambda cid, output_format=None: (_ for _ in ()).throw(
            PubChemError("HTTP error 429: Too Many Requests")
        ),
    )

    result = api.get_compound_properties(2244, ["MolecularWeight"], include_synonyms=True)
    assert result["MolecularWeight"] == "180.16"
    assert result["synonyms_error"].startswith("HTTP error 429")


@pytest.mark.unit
def test_partial_property_result_is_not_cached(api, monkeypatch):
    """The skip_if predicate keeps the incomplete result out of the cache."""
    from provesid.pubchem import _synonyms_incomplete

    complete = {"success": True, "CID": 2244, "synonyms": ["aspirin"]}
    partial = {"success": True, "CID": 2244, "synonyms": [],
               "synonyms_error": "HTTP error 429: Too Many Requests"}
    assert _synonyms_incomplete(partial) is True
    assert _synonyms_incomplete(complete) is False


def _stub_response(payload):
    """Build a minimal object that _parse_response can read as JSON."""
    class _Response:
        status_code = 200
        text = ""
        content = b""

        def json(self):
            return payload

    return _Response()


# --------------------------------------------------------------------------
# PUG-View status-code classification
# --------------------------------------------------------------------------

class _FakeResponse:
    """Minimal stand-in for a requests.Response with a chosen status."""

    def __init__(self, status_code, payload=None, text=""):
        self.status_code = status_code
        self._payload = payload
        self.text = text

    def json(self):
        return self._payload

    def raise_for_status(self):
        import requests
        if self.status_code >= 400:
            raise requests.exceptions.HTTPError(f"{self.status_code} Error")


@pytest.mark.unit
@pytest.mark.parametrize("status", [400, 404])
def test_permanent_client_errors_are_absence_and_are_not_retried(view, monkeypatch, status):
    """An unknown compound (404) or heading (400) must not burn retries."""
    import requests
    attempts = []

    def once(url, timeout=None):
        attempts.append(url)
        return _FakeResponse(status, text="PUGVIEW.BadRequest")

    monkeypatch.setattr(requests, "get", once)
    with pytest.raises(PubChemViewNotFoundError):
        view._make_request("https://example.org/heading")
    assert len(attempts) == 1, "a permanent client error was retried"


@pytest.mark.unit
def test_other_client_errors_are_errors_and_are_not_retried(view, monkeypatch):
    """A 403 is permanent too, but it is a failure rather than absence."""
    import requests
    attempts = []

    def once(url, timeout=None):
        attempts.append(url)
        return _FakeResponse(403, text="Forbidden")

    monkeypatch.setattr(requests, "get", once)
    with pytest.raises(PubChemViewError):
        view._make_request("https://example.org/heading")
    assert len(attempts) == 1


@pytest.mark.unit
def test_server_busy_is_retried_then_raises(monkeypatch):
    """503 ServerBusy is transient: retry, then report failure."""
    import requests
    attempts = []
    view = PubChemView(use_cache=False, max_retries=2, backoff_factor=0.0)

    def busy(url, timeout=None):
        attempts.append(url)
        return _FakeResponse(503, text="PUGVIEW.ServerBusy")

    monkeypatch.setattr(requests, "get", busy)
    with pytest.raises(PubChemViewError):
        view._make_request("https://example.org/heading")
    assert len(attempts) == 3, "503 should be retried up to max_retries"


@pytest.mark.unit
def test_throttling_is_retried(monkeypatch):
    """429 is transient and must be retried, not treated as absence."""
    import requests
    attempts = []
    view = PubChemView(use_cache=False, max_retries=1, backoff_factor=0.0)

    def throttled(url, timeout=None):
        attempts.append(url)
        if len(attempts) == 1:
            return _FakeResponse(429, text="Too Many Requests")
        return _FakeResponse(200, payload={"Record": {"Section": []}})

    monkeypatch.setattr(requests, "get", throttled)
    assert view._make_request("https://example.org/heading") == {"Record": {"Section": []}}
    assert len(attempts) == 2


# --------------------------------------------------------------------------
# Fault-code classification: PubChem sheds load behind 4xx statuses too
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_pugview_transient_fault_behind_a_404_is_retried(monkeypatch):
    """A ServerBusy fault is transient whatever status carries it."""
    import requests
    attempts = []
    view = PubChemView(use_cache=False, max_retries=2, backoff_factor=0.0)

    def busy(url, timeout=None):
        attempts.append(url)
        return _FakeResponse(404, payload={"Fault": {"Code": "PUGVIEW.ServerBusy"}})

    monkeypatch.setattr(requests, "get", busy)
    with pytest.raises(PubChemViewError) as excinfo:
        view._make_request("https://example.org/heading")
    assert not isinstance(excinfo.value, PubChemViewNotFoundError)
    assert len(attempts) == 3, "a transient fault must be retried"


@pytest.mark.unit
def test_pugview_notfound_fault_is_absence_and_is_not_retried(monkeypatch):
    """A NotFound fault is the real thing: absence, no retries."""
    import requests
    attempts = []
    view = PubChemView(use_cache=False, max_retries=2, backoff_factor=0.0)

    def absent(url, timeout=None):
        attempts.append(url)
        return _FakeResponse(404, payload={"Fault": {"Code": "PUGVIEW.NotFound"}})

    monkeypatch.setattr(requests, "get", absent)
    with pytest.raises(PubChemViewNotFoundError):
        view._make_request("https://example.org/heading")
    assert len(attempts) == 1


@pytest.mark.unit
def test_pugrest_transient_fault_behind_a_404_is_a_server_error(api, monkeypatch):
    """PUG-REST ServerBusy behind a 404 must not read as an absent record."""
    import requests
    from provesid.pubchem import PubChemServerError

    monkeypatch.setattr(
        requests, "get",
        lambda url, timeout=None, headers=None: _FakeResponse(
            404, payload={"Fault": {"Code": "PUGREST.ServerBusy"}}
        ),
    )
    with pytest.raises(PubChemServerError):
        api._make_request("https://example.org/synonyms")


@pytest.mark.unit
def test_pugrest_notfound_fault_is_absence(api, monkeypatch):
    """A genuine PUGREST.NotFound stays an absence."""
    import requests

    monkeypatch.setattr(
        requests, "get",
        lambda url, timeout=None, headers=None: _FakeResponse(
            404, payload={"Fault": {"Code": "PUGREST.NotFound"}}
        ),
    )
    with pytest.raises(PubChemNotFoundError):
        api._make_request("https://example.org/synonyms")


# --------------------------------------------------------------------------
# Absence is never persisted
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_absent_property_is_not_cached(tmp_path, monkeypatch):
    """An empty extraction must be retried next time, not remembered."""
    import provesid.cache as cache_module
    from provesid.cache import CacheManager

    probe = CacheManager(cache_dir=str(tmp_path / "pv"), service_name="pubchemview")
    monkeypatch.setitem(cache_module._service_caches, "pubchemview", probe)

    view = PubChemView()
    calls = []

    def sometimes_absent(cid, property_name):
        calls.append((cid, property_name))
        if len(calls) == 1:
            raise PubChemViewNotFoundError("No data found")
        return {"Record": {"Section": [{
            "TOCHeading": "Chemical and Physical Properties",
            "Section": [{
                "TOCHeading": "Experimental Properties",
                "Section": [{
                    "TOCHeading": "Melting Point",
                    "Information": [{"Value": {"StringWithMarkup": [{"String": "135 °C"}]}}],
                }],
            }],
        }]}}

    monkeypatch.setattr(view, "get_property", sometimes_absent)

    assert view.extract_property_data(2244, "Melting Point") == []
    second = view.extract_property_data(2244, "Melting Point")
    assert len(second) == 1, "the empty result was served from the cache"
    assert second[0].value == "135 °C"
    assert len(calls) == 2


@pytest.mark.unit
def test_empty_synonym_list_is_not_cached(tmp_path, monkeypatch):
    """A compound reported as having no synonyms is re-checked next time."""
    import provesid.cache as cache_module
    from provesid.cache import CacheManager

    probe = CacheManager(cache_dir=str(tmp_path / "pc"), service_name="pubchem")
    monkeypatch.setitem(cache_module._service_caches, "pubchem", probe)

    api = PubChemAPI()
    calls = []

    def flaky(url, **kwargs):
        calls.append(url)
        if len(calls) == 1:
            raise PubChemNotFoundError("Resource not found (PUGREST.NotFound)")
        return _stub_response(
            {"InformationList": {"Information": [{"Synonym": ["aspirin"]}]}}
        )

    monkeypatch.setattr(api, "_make_request", flaky)

    assert api.get_compound_synonyms(2244) == []
    assert api.get_compound_synonyms(2244) == ["aspirin"]
    assert len(calls) == 2
