# Network behaviour

Every web client in PROVESID — `PubChemAPI`, `PubChemView`,
`NCIChemicalIdentifierResolver`, `ChEBI`, `CASCommonChem` and `OPSIN` —
makes its requests through one shared transport,
[`provesid.http`](../api/http.md). It decides how long to wait between
requests, when a failure is worth retrying, and how long to back off, so the
answer is the same for every service. Successful answers are also cached on
disk; see [Caching](caching.md). `ClassyFireAPI` is the exception: its three
raw calls still use `requests` directly, and only `get_classification` goes
through the transport. The ClassyFire service has classified nothing new
since February 2023.

## Exceptions

Each client raises its own exception classes, and they share bases, so a
caller can catch one service or all of them:

```python
from provesid.http import ServiceError, NotFoundError

try:
    smiles = resolver.resolve("aspirin", "smiles")
except NotFoundError:
    ...          # any service saying "no such record"
except ServiceError:
    ...          # any service failing for any other reason
```

A raw `requests` exception never reaches the caller. Each exception carries the
`status_code`, the `url` and the `response` it came from; `status_code` and
`response` are `None` when nothing arrived, as with a timeout or a connection
failure.

## What a response means

Services disagree about what a response *means*, and several answer a
permanent condition with a status code that says otherwise:

| Service | Says | Means |
|---|---|---|
| PubChem | `404` with `{"Fault": {"Code": "PUGVIEW.ServerBusy"}}` | ask again later |
| PubChem | `400` with `PUGVIEW.BadRequest` | no such heading — permanent |
| NCI CACTUS | `500` with body `<h1>Page not found (404)</h1>` | no such compound — permanent |
| OPSIN | `404` with `{"status": "FAILURE", "message": "…"}` | the answer, and the only place the reason is |

So the retry policy is fixed and the *reading* of a response is per service:
each client passes a `classify` function that maps a response to an
[`Outcome`][provesid.http.Outcome]: `OK`, `ABSENT`, `RETRY` or `FATAL`.
`provesid.pubchem.pugrest_classify` and `pugview_classify` read PubChem's
fault codes, `provesid.resolver.nci_classify` reads CACTUS's body, and
`provesid.opsin.opsin_classify` treats OPSIN's 404 as a success. Services that
use status codes honestly, ChEBI and CAS Common Chemistry, use
`default_classify`.

## Retries

- Retried: a `RETRY` verdict, request timeouts and connection errors.
- Not retried: `ABSENT` (the record does not exist) and `FATAL` (a malformed
  request, a rejected key) — asking again cannot change either.
- A request is made at most `max_retries + 1` times.
- The wait is `backoff * 2 ** attempt`, capped at `max_backoff`, unless the
  service sent a `Retry-After`, in seconds or as an HTTP date, in which case
  that wins. A `Retry-After` longer than `max_backoff` is not cut short: the
  call gives up instead, because asking before the time the service named
  only earns a second refusal.
- Retries are paced like any other request.
- `max_elapsed`, when set, caps the *total* waiting. Both PubChem clients set
  a 10-second budget (`provesid.pubchem.RETRY_WAIT_BUDGET`), which allows
  1 + 2 + 4 seconds for a transient 500 or a timeout but declines PubChem's
  `Retry-After: 30`.

PubChem sends that header when it has throttled or blocked an IP, and such a
block does not lift in thirty seconds: waiting the full thirty and asking again
returns the same 503. For an unattended bulk job that would rather wait a
throttle out, raise the budget:

```python
api._http.max_elapsed = 180
```

## Pacing is per host

A client's `min_interval` is its own limit on how fast it asks. The clock it
is measured against belongs to the **host**, because the limit being respected
does too: PubChem publishes five requests per second *per IP*, not per Python
object. Every client aimed at one host shares one
[`RateLimiter`][provesid.http.RateLimiter]:

```python
>>> from provesid.pubchem import PubChemAPI
>>> from provesid.pubchemview import PubChemView
>>> PubChemAPI()._http.limiter is PubChemView()._http.limiter
True
```

The effective rate for a host is set by its most impatient client: sharing the
clock stops two clients doubling a limit, but not one client configured with
`min_interval=0.01` from exceeding it alone.

## The circuit breaker

A `Retry-After` says when the *host* will next answer, not when one request
may be repeated, so it is kept on the host's shared `RateLimiter` as a "not
before" time. Every client aimed at that host checks it before asking:

- a hold no longer than the client would wait anyway (`max_backoff`, and what
  is left of `max_elapsed`) is waited out, and the request is then sent;
- a longer hold fails the call at once with the client's rate-limit exception
  (`PubChemServerError`, `PubChemViewError`, …), **without sending a
  request**.

So resolving a thousand names against a PubChem that has blocked this IP costs
one refused request, not a thousand; the other 999 fail in microseconds:

```python
>>> api = PubChemAPI()
>>> api.get_compound_synonyms(2244)        # PubChem: 503, Retry-After: 30
PubChemServerError: ... stopped retrying after 0s of its 10s retry budget ...
>>> PubChemView().get_property(2244, "Boiling Point")   # no request sent
PubChemViewError: pubchem.ncbi.nlm.nih.gov asked for no requests until 14:02:31 ...
```

The refusal carries `held_until`, the end of the hold. `Search` logs a held
host at DEBUG; the transport has already warned once, when the hold was
recorded, so a batch does not print a warning per query.

A shorter `Retry-After` never shortens a longer hold already in place, and
only a `Retry-After` sets one: a bare 503 is one request's bad luck, not a
statement about the host. After a network change, say, clear every hold with
[`release_holds`][provesid.http.release_holds], or one host's with
`host_limiter(url).release()`. Holds, like the clock, are per process.

## Sessions

A client can make its calls through a `requests.Session`, for connection
pooling and persistent headers. `ChEBI` does: its `User-Agent` and `Accept`
headers live on the session, and an ontology walk reuses one connection across
many small requests.
