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
    RateLimiter,
    RateLimitError,
    ServiceError,
    ServiceTimeoutError,
    default_classify,
    host_limiter,
    release_holds,
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
def test_a_retry_after_longer_than_max_backoff_gives_up_instead_of_asking_early(monkeypatch, sleeps):
    """
    A service asking for an hour does not get to stall the caller --- and does
    not get asked again after five seconds either, which would only earn a
    second refusal. The call gives up after the one request.
    """
    recorder = _Recorder(_FakeResponse(429, headers={"Retry-After": "3600"}), _FakeResponse(200, text="ok"))
    monkeypatch.setattr(requests, "get", recorder)

    client = HTTPClient(min_interval=0, max_retries=2, backoff=1.0, max_backoff=5.0,
                        rate_limit_cls=RateLimitError)
    with pytest.raises(RateLimitError, match="longer than max_backoff=5s"):
        client.get_text("https://example.org/x")
    assert sleeps == []
    assert recorder.count == 1


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


# --------------------------------------------------------------------------
# Pacing is per host, not per object
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_two_clients_on_one_host_share_a_clock():
    """
    PubChem's five per second is a per-IP budget, so objects cannot each keep
    their own clock and still honour it.
    """
    a = HTTPClient(min_interval=0.2, pace_host="https://pubchem.ncbi.nlm.nih.gov/rest/pug")
    b = HTTPClient(min_interval=0.2, pace_host="https://pubchem.ncbi.nlm.nih.gov/rest/pug_view")

    assert a.limiter is b.limiter
    assert a.limiter is host_limiter("pubchem.ncbi.nlm.nih.gov")


@pytest.mark.unit
def test_clients_on_different_hosts_do_not_share_a_clock():
    """One slow service must not pace a fast unrelated one."""
    a = HTTPClient(pace_host="https://pubchem.ncbi.nlm.nih.gov/rest/pug")
    b = HTTPClient(pace_host="https://www.ebi.ac.uk/chebi/backend/api/public")
    assert a.limiter is not b.limiter


@pytest.mark.unit
def test_a_client_with_no_pace_host_paces_alone():
    """A stub, or a service with no shared budget, keeps its own clock."""
    a = HTTPClient(min_interval=0.2)
    b = HTTPClient(min_interval=0.2)
    assert a.limiter is not b.limiter


@pytest.mark.unit
def test_the_shared_clock_actually_delays_the_second_client(monkeypatch, sleeps):
    """
    The point of sharing: a request waits for the last request *anyone* made.
    """
    monkeypatch.setattr(requests, "get", _Recorder(_FakeResponse(200, text="ok")))
    host = "https://shared.example.test/api"
    first = HTTPClient(min_interval=0.5, pace_host=host)
    second = HTTPClient(min_interval=0.5, pace_host=host)

    first.get_text(host)
    sleeps.clear()
    second.get_text(host)

    assert len(sleeps) == 1, "the second client did not wait for the first"
    assert 0 < sleeps[0] <= 0.5


@pytest.mark.unit
def test_last_request_time_stays_per_client():
    """
    A client reports when *it* last asked, even though the clock it measures
    against is shared. A fresh client has never asked.
    """
    host = "https://fresh.example.test/api"
    warm = HTTPClient(pace_host=host)
    warm.rate_limit()

    assert warm.last_request_time > 0
    assert HTTPClient(pace_host=host).last_request_time == 0.0


@pytest.mark.unit
def test_a_bare_limiter_records_even_when_pacing_is_off(sleeps):
    """The request happened, so the next caller has to be able to see it."""
    limiter = RateLimiter()
    assert limiter.last_request_time == 0.0
    limiter.wait(0.0)
    assert limiter.last_request_time > 0
    assert sleeps == []


# --------------------------------------------------------------------------
# The circuit breaker: a Retry-After holds the host, not the request
# --------------------------------------------------------------------------

def _throttled(seconds="30", status=503):
    return _FakeResponse(status, text="busy", headers={"Retry-After": seconds})


@pytest.mark.unit
def test_a_held_host_fails_at_once_without_a_request(monkeypatch, sleeps):
    """
    The point of the breaker. PubChem blocks an IP with ``Retry-After: 30``;
    the first call learns it at the price of one request, and the next call
    is told without asking.
    """
    recorder = _Recorder(_throttled("30"))
    monkeypatch.setattr(requests, "get", recorder)
    client = HTTPClient(min_interval=0, max_retries=3, max_elapsed=10,
                        rate_limit_cls=RateLimitError,
                        retry_exhausted_cls=RateLimitError)

    with pytest.raises(RateLimitError, match="retry budget"):
        client.get("https://example.org/a")
    assert recorder.count == 1

    with pytest.raises(RateLimitError, match="asked for no requests until") as refused:
        client.get("https://example.org/b")
    assert recorder.count == 1, "a held host was asked again"
    assert sleeps == []
    assert refused.value.held_until == client.limiter.not_before


@pytest.mark.unit
def test_only_the_breaker_marks_its_refusals_as_held(monkeypatch, sleeps):
    """
    ``held_until`` is how a batch tells "the host is still held, as already
    reported" from a failure worth a warning of its own.
    """
    monkeypatch.setattr(requests, "get", _Recorder(_throttled("30")))
    client = HTTPClient(min_interval=0, max_retries=0, max_backoff=5)

    with pytest.raises(ServiceError) as refused_by_host:
        client.get("https://example.org/a")
    assert refused_by_host.value.held_until is None

    with pytest.raises(ServiceError) as refused_by_breaker:
        client.get("https://example.org/a")
    assert refused_by_breaker.value.held_until is not None


@pytest.mark.unit
def test_a_plain_exception_class_still_carries_held_until(sleeps):
    """The attribute is set by hand, so any ``rate_limit_cls`` gets it."""
    client = HTTPClient(min_interval=0, max_backoff=1, rate_limit_cls=RuntimeError)
    client.limiter.hold(60)
    with pytest.raises(RuntimeError) as refused:
        client.get("https://example.org/x")
    assert refused.value.held_until == client.limiter.not_before


@pytest.mark.unit
def test_the_hold_binds_every_client_on_the_host(monkeypatch, sleeps):
    """
    ``Retry-After`` is information about the host, so a second client aimed at
    it --- PubChemView beside PubChemAPI --- respects what the first was told.
    """
    recorder = _Recorder(_throttled("30"))
    monkeypatch.setattr(requests, "get", recorder)
    host = "https://breaker.example.test/api"
    first = HTTPClient(min_interval=0, max_retries=0, pace_host=host)
    second = HTTPClient(min_interval=0, max_retries=3, max_backoff=5, pace_host=host,
                        rate_limit_cls=RateLimitError)

    with pytest.raises(ServiceError):
        first.get(host)
    with pytest.raises(RateLimitError, match="breaker.example.test asked for no requests"):
        second.get(host)
    assert recorder.count == 1


@pytest.mark.unit
def test_a_hold_on_one_host_does_not_touch_another(monkeypatch, sleeps):
    """Being throttled by PubChem says nothing about ChEBI."""
    monkeypatch.setattr(requests, "get", _Recorder(_FakeResponse(200, text="ok")))
    host_limiter("https://held.example.test").hold(600)

    other = HTTPClient(min_interval=0, pace_host="https://free.example.test")
    assert other.get_text("https://free.example.test/x") == "ok"


@pytest.mark.unit
def test_a_short_hold_is_waited_out_then_asked(monkeypatch, sleeps):
    """
    A hold the client would have waited anyway is waited, not failed: a
    one-second pause is cheaper than a lost answer.
    """
    recorder = _Recorder(_FakeResponse(200, text="ok"))
    monkeypatch.setattr(requests, "get", recorder)
    client = HTTPClient(min_interval=0, max_backoff=5)
    client.limiter.hold(2)

    assert client.get_text("https://example.org/x") == "ok"
    assert len(sleeps) == 1 and 1.5 < sleeps[0] <= 2
    assert recorder.count == 1


@pytest.mark.unit
def test_the_hold_wait_counts_against_max_elapsed(monkeypatch, sleeps):
    """
    Waiting out a hold is waiting, so it comes out of the same budget a retry
    would spend. Here the hold fits, and the retry it leaves no room for is
    declined.
    """
    recorder = _Recorder(_FakeResponse(503, text="busy"))
    monkeypatch.setattr(requests, "get", recorder)
    client = HTTPClient(min_interval=0, max_retries=3, backoff=4, max_elapsed=10)
    client.limiter.hold(8)

    with pytest.raises(ServiceError, match="retry budget"):
        client.get("https://example.org/x")
    assert recorder.count == 1
    assert len(sleeps) == 1 and 7.5 < sleeps[0] <= 8


@pytest.mark.unit
def test_a_longer_hold_is_never_shortened():
    """The host's latest word does not retract what it told another request."""
    limiter = RateLimiter()
    limiter.hold(60)
    limiter.hold(5)
    assert limiter.held_for() > 55
    limiter.hold(0)
    assert limiter.held_for() > 55


@pytest.mark.unit
def test_a_response_without_retry_after_holds_nothing(monkeypatch, sleeps):
    """A plain 503 is this request's bad luck, not a statement about the host."""
    monkeypatch.setattr(requests, "get", _Recorder(_FakeResponse(503)))
    client = HTTPClient(min_interval=0, max_retries=1, backoff=0)

    with pytest.raises(ServiceError):
        client.get("https://example.org/x")
    assert client.limiter.held_for() == 0.0


@pytest.mark.unit
def test_release_holds_reopens_every_host(monkeypatch, sleeps):
    """For a caller who knows better than the header, and for the test suite."""
    monkeypatch.setattr(requests, "get", _Recorder(_FakeResponse(200, text="ok")))
    host = "https://released.example.test"
    host_limiter(host).hold(600)

    release_holds()
    assert HTTPClient(min_interval=0, pace_host=host).get_text(host) == "ok"
    assert sleeps == []


# --------------------------------------------------------------------------
# The total wait budget
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_max_elapsed_stops_retrying_before_the_attempts_run_out(monkeypatch, sleeps):
    """
    A service that says "come back in 30s" three times must not cost 90s.

    PubChem does exactly this to a throttled IP, which is why the budget
    exists.
    """
    recorder = _Recorder(_FakeResponse(503, text="busy", headers={"Retry-After": "30"}))
    monkeypatch.setattr(requests, "get", recorder)

    client = HTTPClient(max_retries=3, max_elapsed=30, min_interval=0)
    with pytest.raises(ServiceError, match="stopped retrying after 30s of its 30s retry budget"):
        client.get("https://example.org/x")

    assert sleeps == [30.0], "the budget allowed exactly one wait"
    assert recorder.count == 2


@pytest.mark.unit
def test_max_elapsed_does_not_curtail_cheap_retries(monkeypatch, sleeps):
    """A fast back-off curve fits inside the budget and is spent in full."""
    recorder = _Recorder(_FakeResponse(503, text="busy"))
    monkeypatch.setattr(requests, "get", recorder)

    client = HTTPClient(max_retries=3, backoff=0.5, max_elapsed=30, min_interval=0)
    with pytest.raises(ServiceError):
        client.get("https://example.org/x")

    assert sleeps == [0.5, 1.0, 2.0]
    assert recorder.count == 4


@pytest.mark.unit
def test_no_budget_means_max_retries_is_the_only_bound(monkeypatch, sleeps):
    """The default is unchanged: retry until the attempts run out."""
    recorder = _Recorder(_FakeResponse(503, text="busy", headers={"Retry-After": "30"}))
    monkeypatch.setattr(requests, "get", recorder)

    with pytest.raises(ServiceError, match="after 3 attempt"):
        HTTPClient(max_retries=2, min_interval=0).get("https://example.org/x")

    assert sleeps == [30.0, 30.0]
    assert recorder.count == 3


# --------------------------------------------------------------------------
# Telling one failure from another
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_a_busy_service_and_a_bad_request_raise_different_classes(monkeypatch):
    """
    "The service kept failing" and "the request was wrong" are different
    things, and PubChem's callers have always been able to tell them apart.
    """
    class MyError(ServiceError):
        pass

    class MyBusy(MyError):
        pass

    client = HTTPClient(error_cls=MyError, retry_exhausted_cls=MyBusy,
                        max_retries=1, backoff=0.0, min_interval=0)

    monkeypatch.setattr(requests, "get", _Recorder(_FakeResponse(503, text="busy")))
    with pytest.raises(MyBusy):
        client.get("https://example.org/x")

    monkeypatch.setattr(requests, "get", _Recorder(_FakeResponse(403, text="nope")))
    with pytest.raises(MyError) as excinfo:
        client.get("https://example.org/x")
    assert not isinstance(excinfo.value, MyBusy)


@pytest.mark.unit
def test_retry_exhausted_defaults_to_the_error_class(monkeypatch):
    """A client that does not care keeps one exception for both."""
    monkeypatch.setattr(requests, "get", _Recorder(_FakeResponse(503, text="busy")))
    with pytest.raises(ServiceError):
        HTTPClient(max_retries=0, min_interval=0).get("https://example.org/x")


@pytest.mark.unit
@pytest.mark.parametrize("status,payload", [(404, None), (403, None), (503, None)])
def test_the_exception_carries_the_status_and_the_url(monkeypatch, status, payload):
    """
    CAS distinguishes a rejected key from an unknown CAS number by status, so
    the status has to survive the raise.
    """
    monkeypatch.setattr(requests, "get", _Recorder(_FakeResponse(status, text="x")))

    with pytest.raises(ServiceError) as excinfo:
        HTTPClient(max_retries=0, min_interval=0).get("https://example.org/thing")

    assert excinfo.value.status_code == status
    assert excinfo.value.url == "https://example.org/thing"
    assert excinfo.value.response is not None


@pytest.mark.unit
def test_a_timeout_carries_no_status(monkeypatch, sleeps):
    """Nothing arrived, so there is nothing to report a status for."""
    monkeypatch.setattr(requests, "get", _Recorder(requests.Timeout("gone")))

    with pytest.raises(ServiceTimeoutError) as excinfo:
        HTTPClient(max_retries=0, min_interval=0,
                   timeout_cls=ServiceTimeoutError).get("https://example.org/x")

    assert excinfo.value.status_code is None
    assert excinfo.value.response is None


@pytest.mark.unit
def test_a_plain_exception_class_is_still_usable(monkeypatch):
    """
    The detail is keyword-only, so a class that never declared it must not be
    handed it --- that would turn a service failure into a TypeError.
    """
    class Plain(Exception):
        pass

    monkeypatch.setattr(requests, "get", _Recorder(_FakeResponse(403, text="nope")))

    with pytest.raises(Plain):
        HTTPClient(error_cls=Plain, max_retries=0, min_interval=0).get(
            "https://example.org/x"
        )


# --------------------------------------------------------------------------
# Sessions
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_a_session_is_used_when_one_is_given(monkeypatch):
    """
    ChEBI keeps a session for its pooled connection and persistent headers, so
    the transport has to call through it rather than around it.
    """
    recorder = _Recorder(_FakeResponse(200, text="ok"))
    session = requests.Session()
    monkeypatch.setattr(session, "get", recorder)
    monkeypatch.setattr(requests, "get", _Recorder(_FakeResponse(500, text="wrong path")))

    assert HTTPClient(session=session).get_text("https://example.org/x") == "ok"
    assert recorder.count == 1


@pytest.mark.unit
def test_without_a_session_the_module_function_is_called(monkeypatch):
    """
    Patching ``requests.get`` must keep working: several suites rely on it, and
    a ``from requests import get`` or an always-on session would defeat them.
    """
    recorder = _Recorder(_FakeResponse(200, text="ok"))
    monkeypatch.setattr(requests, "get", recorder)

    assert HTTPClient().get_text("https://example.org/x") == "ok"
    assert recorder.count == 1
