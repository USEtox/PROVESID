"""
Tests for the shared HTTP transport.

Fully offline: every test stubs ``requests.get``/``requests.post`` and records
what the client did, because the only interesting questions about a transport
are *how many* requests it made and *how long* it waited between them. Nothing
here touches the network.
"""

import email.utils
import time
from datetime import datetime, timedelta, timezone

import pytest
import requests

from provesid.http import (
    HTTPClient,
    NotFoundError,
    Outcome,
    RateLimitError,
    ServiceError,
    ServiceTimeoutError,
    default_classify,
    retry_after_seconds,
)


class _FakeResponse:
    """A stand-in for ``requests.Response`` carrying just what is classified."""

    def __init__(self, status_code, payload=None, text="", headers=None):
        self.status_code = status_code
        self._payload = payload
        self.text = text if payload is None else str(payload)
        self.headers = headers or {}
        self.content = self.text.encode()

    def json(self):
        if self._payload is None:
            raise ValueError("no JSON body")
        return self._payload


class _Recorder:
    """Records each call and replays a fixed list of responses."""

    def __init__(self, *responses):
        self.responses = list(responses)
        self.calls = []

    def __call__(self, url, **kwargs):
        self.calls.append((url, kwargs))
        result = self.responses[min(len(self.calls) - 1, len(self.responses) - 1)]
        if isinstance(result, Exception):
            raise result
        return result

    @property
    def count(self):
        return len(self.calls)


@pytest.fixture
def sleeps(monkeypatch):
    """Capture every sleep the client asks for instead of serving it."""
    recorded = []

    def fake_sleep(seconds):
        recorded.append(seconds)

    monkeypatch.setattr(time, "sleep", fake_sleep)
    return recorded


# --------------------------------------------------------------------------
# The request itself
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_get_text_returns_the_stripped_body(monkeypatch):
    """A 200 comes back as text with surrounding whitespace removed."""
    monkeypatch.setattr(requests, "get", _Recorder(_FakeResponse(200, text="  CCO\n")))
    assert HTTPClient().get_text("https://example.org/x") == "CCO"


@pytest.mark.unit
def test_get_json_returns_the_parsed_body(monkeypatch):
    """A JSON body is decoded once, by the transport."""
    monkeypatch.setattr(requests, "get", _Recorder(_FakeResponse(200, payload={"ok": True})))
    assert HTTPClient().get_json("https://example.org/x") == {"ok": True}


@pytest.mark.unit
def test_only_url_and_timeout_are_passed_when_nothing_else_is_set(monkeypatch):
    """
    The plain case must read ``requests.get(url, timeout=...)``.

    Several suites stub ``requests.get`` with exactly that signature. A client
    that always passed ``params=None, headers=None`` would break them for no
    gain, so this pins the call shape rather than leaving it to chance.
    """
    calls = []

    def only_url_and_timeout(url, timeout=None):
        calls.append((url, timeout))
        return _FakeResponse(200, text="ok")

    monkeypatch.setattr(requests, "get", only_url_and_timeout)
    assert HTTPClient(timeout=7).get_text("https://example.org/x") == "ok"
    assert calls == [("https://example.org/x", 7)]


@pytest.mark.unit
def test_headers_and_params_are_passed_when_given(monkeypatch):
    """Anything actually supplied reaches ``requests``."""
    recorder = _Recorder(_FakeResponse(200, text="ok"))
    monkeypatch.setattr(requests, "get", recorder)

    client = HTTPClient(headers={"Accept": "text/plain"})
    client.get_text("https://example.org/x", params={"q": "1"}, headers={"X-Key": "k"})

    _, kwargs = recorder.calls[0]
    assert kwargs["params"] == {"q": "1"}
    assert kwargs["headers"] == {"Accept": "text/plain", "X-Key": "k"}


@pytest.mark.unit
def test_post_json_sends_the_body(monkeypatch):
    """POST carries its JSON body and decodes the answer."""
    recorder = _Recorder(_FakeResponse(200, payload={"n": 2}))
    monkeypatch.setattr(requests, "post", recorder)

    assert HTTPClient().post_json("https://example.org/q", json={"n": 1}) == {"n": 2}
    assert recorder.calls[0][1]["json"] == {"n": 1}


@pytest.mark.unit
def test_unsupported_method_is_rejected():
    """Only GET and POST are supported; anything else is a programming error."""
    with pytest.raises(ValueError, match="Unsupported HTTP method"):
        HTTPClient().request("DELETE", "https://example.org/x")


# --------------------------------------------------------------------------
# Pacing
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_rate_limit_paces_successive_calls():
    """Consecutive calls are held at least ``min_interval`` apart."""
    client = HTTPClient(min_interval=0.05)
    start = time.time()
    for _ in range(3):
        client.rate_limit()
    assert time.time() - start >= 0.05 * 2 * 0.9


@pytest.mark.unit
def test_rate_limit_is_a_no_op_when_pacing_is_off(sleeps):
    """``min_interval=0`` never sleeps."""
    client = HTTPClient(min_interval=0)
    for _ in range(5):
        client.rate_limit()
    assert sleeps == []


@pytest.mark.unit
def test_every_attempt_is_paced(monkeypatch):
    """
    A retry is paced like any other request.

    A service already shedding load must not be asked again faster than a
    healthy one, so the pacing sits inside the retry loop rather than in front
    of it.
    """
    recorder = _Recorder(_FakeResponse(503, text="busy"))
    monkeypatch.setattr(requests, "get", recorder)

    client = HTTPClient(min_interval=0.02, max_retries=2, backoff=0.0)
    start = time.time()
    with pytest.raises(ServiceError):
        client.get("https://example.org/x")

    assert recorder.count == 3
    assert time.time() - start >= 0.02 * 2 * 0.9


# --------------------------------------------------------------------------
# Classification
# --------------------------------------------------------------------------

@pytest.mark.unit
@pytest.mark.parametrize("status,expected", [
    (200, Outcome.OK), (204, Outcome.OK),
    (400, Outcome.FATAL), (401, Outcome.FATAL), (403, Outcome.FATAL),
    (404, Outcome.ABSENT),
    (429, Outcome.RETRY), (500, Outcome.RETRY), (502, Outcome.RETRY),
    (503, Outcome.RETRY), (504, Outcome.RETRY),
])
def test_default_classify_reads_the_status(status, expected):
    """The status-only reading, for services that use status codes honestly."""
    assert default_classify(_FakeResponse(status)) is expected


@pytest.mark.unit
def test_404_is_absence_and_is_not_retried(monkeypatch):
    """Absence is a statement about the data; asking again cannot change it."""
    recorder = _Recorder(_FakeResponse(404, text="nope"))
    monkeypatch.setattr(requests, "get", recorder)

    with pytest.raises(NotFoundError):
        HTTPClient(max_retries=3, backoff=0.0).get("https://example.org/x")
    assert recorder.count == 1


@pytest.mark.unit
def test_a_permanent_4xx_is_not_retried(monkeypatch):
    """A rejected key or a malformed request will be rejected again."""
    recorder = _Recorder(_FakeResponse(401, text="bad key"))
    monkeypatch.setattr(requests, "get", recorder)

    with pytest.raises(ServiceError) as excinfo:
        HTTPClient(max_retries=3, backoff=0.0).get("https://example.org/x")
    assert not isinstance(excinfo.value, NotFoundError)
    assert recorder.count == 1
    assert "bad key" in str(excinfo.value)


@pytest.mark.unit
def test_5xx_is_retried_exactly_max_retries_times(monkeypatch):
    """The budget is ``max_retries`` retries after the first attempt."""
    recorder = _Recorder(_FakeResponse(503, text="busy"))
    monkeypatch.setattr(requests, "get", recorder)

    with pytest.raises(ServiceError):
        HTTPClient(max_retries=2, backoff=0.0).get("https://example.org/x")
    assert recorder.count == 3


@pytest.mark.unit
def test_a_transient_failure_that_clears_is_invisible(monkeypatch):
    """A caller never hears about a 503 that the next attempt answered."""
    recorder = _Recorder(_FakeResponse(503), _FakeResponse(200, text="ok"))
    monkeypatch.setattr(requests, "get", recorder)

    assert HTTPClient(max_retries=3, backoff=0.0).get_text("https://example.org/x") == "ok"
    assert recorder.count == 2


@pytest.mark.unit
def test_sustained_throttling_raises_the_rate_limit_class(monkeypatch):
    """429 until the budget runs out is reported as throttling, not as absence."""
    monkeypatch.setattr(requests, "get", _Recorder(_FakeResponse(429, text="slow down")))

    client = HTTPClient(max_retries=1, backoff=0.0, rate_limit_cls=RateLimitError)
    with pytest.raises(RateLimitError):
        client.get("https://example.org/x")


@pytest.mark.unit
def test_a_custom_classifier_decides(monkeypatch):
    """
    The classifier hook is what lets a service override its own status codes.

    Both PubChem and CACTUS answer a permanent condition with a status that
    says otherwise, so this is not a hypothetical.
    """
    recorder = _Recorder(_FakeResponse(500, text="Page not found (404)"))
    monkeypatch.setattr(requests, "get", recorder)

    def body_decides(response):
        if "Page not found" in response.text:
            return Outcome.ABSENT
        return default_classify(response)

    client = HTTPClient(max_retries=3, backoff=0.0, classify=body_decides)
    with pytest.raises(NotFoundError):
        client.get("https://example.org/x")
    assert recorder.count == 1, "a classifier saying ABSENT must stop the retries"


# --------------------------------------------------------------------------
# Back-off and Retry-After
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_backoff_is_exponential_without_a_retry_after(monkeypatch, sleeps):
    """Without guidance from the service, wait 1x, 2x, 4x the base."""
    monkeypatch.setattr(requests, "get", _Recorder(_FakeResponse(503)))

    with pytest.raises(ServiceError):
        HTTPClient(min_interval=0, max_retries=3, backoff=0.5).get("https://example.org/x")
    assert sleeps == [0.5, 1.0, 2.0]


@pytest.mark.unit
def test_retry_after_in_seconds_is_honoured(monkeypatch, sleeps):
    """A service that says when to come back knows better than any curve."""
    monkeypatch.setattr(
        requests, "get",
        _Recorder(_FakeResponse(429, headers={"Retry-After": "3"}), _FakeResponse(200, text="ok")),
    )

    client = HTTPClient(min_interval=0, max_retries=2, backoff=10.0)
    assert client.get_text("https://example.org/x") == "ok"
    assert sleeps == [3.0], "Retry-After must win over the exponential curve"


@pytest.mark.unit
def test_retry_after_as_an_http_date_is_honoured(monkeypatch, sleeps):
    """RFC 9110 allows a date; the wait is the distance to it."""
    when = datetime.now(timezone.utc) + timedelta(seconds=30)
    header = email.utils.format_datetime(when)
    monkeypatch.setattr(
        requests, "get",
        _Recorder(_FakeResponse(503, headers={"Retry-After": header}), _FakeResponse(200, text="ok")),
    )

    client = HTTPClient(min_interval=0, max_retries=2, backoff=1.0, max_backoff=120)
    assert client.get_text("https://example.org/x") == "ok"
    assert len(sleeps) == 1
    assert 25 <= sleeps[0] <= 31


@pytest.mark.unit
def test_retry_after_is_capped_by_max_backoff(monkeypatch, sleeps):
    """A service asking for an hour does not get to stall the caller."""
    monkeypatch.setattr(
        requests, "get",
        _Recorder(_FakeResponse(429, headers={"Retry-After": "3600"}), _FakeResponse(200, text="ok")),
    )

    client = HTTPClient(min_interval=0, max_retries=2, backoff=1.0, max_backoff=5.0)
    assert client.get_text("https://example.org/x") == "ok"
    assert sleeps == [5.0]


@pytest.mark.unit
def test_exponential_backoff_is_capped_too(monkeypatch, sleeps):
    """The curve flattens at ``max_backoff``."""
    monkeypatch.setattr(requests, "get", _Recorder(_FakeResponse(503)))

    with pytest.raises(ServiceError):
        HTTPClient(min_interval=0, max_retries=4, backoff=1.0, max_backoff=2.0).get(
            "https://example.org/x"
        )
    assert sleeps == [1.0, 2.0, 2.0, 2.0]


@pytest.mark.unit
@pytest.mark.parametrize("header,expected", [
    ({"Retry-After": "12"}, 12.0),
    ({"Retry-After": " 4 "}, 4.0),
    ({"Retry-After": "not a number"}, None),
    ({"Retry-After": "Mon, 01 Jan 1990 00:00:00 GMT"}, None),
    ({}, None),
])
def test_retry_after_seconds_parsing(header, expected):
    """Every form the header can take, including the ones to ignore."""
    assert retry_after_seconds(_FakeResponse(200, headers=header)) == expected


# --------------------------------------------------------------------------
# Network failures
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_a_timeout_is_retried_then_raises_the_timeout_class(monkeypatch, sleeps):
    """A timeout is transient until the budget says otherwise."""
    recorder = _Recorder(requests.Timeout("too slow"))
    monkeypatch.setattr(requests, "get", recorder)

    client = HTTPClient(min_interval=0, max_retries=2, backoff=0.0,
                        timeout_cls=ServiceTimeoutError)
    with pytest.raises(ServiceTimeoutError):
        client.get("https://example.org/x")
    assert recorder.count == 3


@pytest.mark.unit
def test_a_connection_error_is_retried_and_can_clear(monkeypatch, sleeps):
    """A connection that fails once and then works is not an error."""
    recorder = _Recorder(requests.ConnectionError("refused"), _FakeResponse(200, text="ok"))
    monkeypatch.setattr(requests, "get", recorder)

    client = HTTPClient(min_interval=0, max_retries=2, backoff=0.0)
    assert client.get_text("https://example.org/x") == "ok"
    assert recorder.count == 2


@pytest.mark.unit
def test_no_raw_requests_exception_reaches_the_caller(monkeypatch):
    """
    Dev-principle 8: a caller sees this package's exceptions, never
    ``requests``'.
    """
    monkeypatch.setattr(requests, "get", _Recorder(requests.TooManyRedirects("loop")))

    with pytest.raises(ServiceError) as excinfo:
        HTTPClient(max_retries=1, backoff=0.0).get("https://example.org/x")
    assert not isinstance(excinfo.value, requests.RequestException)


@pytest.mark.unit
def test_a_non_json_body_is_reported_as_a_service_error(monkeypatch):
    """A 200 carrying HTML has failed as surely as a 500 has."""
    monkeypatch.setattr(requests, "get", _Recorder(_FakeResponse(200, text="<html>oops</html>")))

    with pytest.raises(ServiceError, match="not JSON"):
        HTTPClient().get_json("https://example.org/x")


# --------------------------------------------------------------------------
# What the clients get out of it
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_each_client_raises_its_own_exception_classes(monkeypatch):
    """
    The transport is shared; the exceptions are not.

    A caller keeps catching the service it called, and can now also catch every
    service at once through the shared bases.
    """
    class MyError(ServiceError):
        pass

    class MyNotFound(MyError, NotFoundError):
        pass

    monkeypatch.setattr(requests, "get", _Recorder(_FakeResponse(404)))
    client = HTTPClient(error_cls=MyError, not_found_cls=MyNotFound, backoff=0.0)

    with pytest.raises(MyNotFound):
        client.get("https://example.org/x")
    with pytest.raises(MyError):
        client.get("https://example.org/x")
    with pytest.raises(NotFoundError):
        client.get("https://example.org/x")
    with pytest.raises(ServiceError):
        client.get("https://example.org/x")


@pytest.mark.unit
def test_the_shared_bases_catch_across_modules():
    """Every client's exceptions descend from the shared bases."""
    from provesid.pubchem import PubChemError, PubChemNotFoundError, PubChemTimeoutError
    from provesid.pubchemview import PubChemViewError, PubChemViewNotFoundError
    from provesid.resolver import (
        NCIResolverError,
        NCIResolverNotFoundError,
        NCIResolverTimeoutError,
    )

    for cls in (PubChemError, PubChemViewError, NCIResolverError):
        assert issubclass(cls, ServiceError)
    for cls in (PubChemNotFoundError, PubChemViewNotFoundError, NCIResolverNotFoundError):
        assert issubclass(cls, NotFoundError)
    for cls in (PubChemTimeoutError, NCIResolverTimeoutError):
        assert issubclass(cls, ServiceTimeoutError)

    # And each still catches as its own service, which is what callers wrote.
    assert issubclass(PubChemViewNotFoundError, PubChemViewError)
    assert issubclass(NCIResolverNotFoundError, NCIResolverError)
