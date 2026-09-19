# HTTP Transport

Every PROVESID web-API client makes its requests through one shared transport.
It is the single place that decides how long to wait between requests, when a
failure is worth retrying, and how long to back off before asking again.

## Why it exists

Each client used to carry its own copy of that logic. The six copies had
drifted apart: none of them honoured `Retry-After`, only one retried at all,
and three modules — ChEBI, CAS Common Chemistry and OPSIN — had no rate
limiting whatsoever. Centralising it is what
[development principle 8](https://github.com/USEtox/PROVESID) asks for, and it
means a fix to the retry policy is a fix everywhere.

## What a caller sees

Nothing changes for a caller. Each client keeps its own exception classes and
passes them to the transport, so `except PubChemViewError` and
`except NCIResolverError` work exactly as before. What is new is that those
classes now also descend from shared bases, so a caller can catch every
service at once:

```python
from provesid.http import ServiceError, NotFoundError

try:
    smiles = resolver.resolve("aspirin", "smiles")
except NotFoundError:
    ...          # any service saying "no such record"
except ServiceError:
    ...          # any service failing for any other reason
```

A raw `requests` exception never reaches the caller.

Each exception also carries the status, the URL and the response it came from,
which is how a client tells one failure from another without re-reading the
wire:

```python
from provesid.cascommonchem import CASCommonChemError

try:
    detail = cas._http.get_json(url)
except CASCommonChemError as exc:
    if exc.status_code == 401:
        ...      # the key was rejected, not the CAS number
```

`status_code` and `response` are `None` when nothing arrived — a timeout or a
connection failure.

## Classifying a response

The one thing services genuinely disagree about is what a response *means*,
and two of the services this package talks to answer a permanent condition
with a status code that says otherwise:

| Service | Says | Means |
|---|---|---|
| PubChem | `404` with `{"Fault": {"Code": "PUGVIEW.ServerBusy"}}` | ask again later |
| PubChem | `400` with `PUGVIEW.BadRequest` | no such heading — permanent |
| NCI CACTUS | `500` with body `<h1>Page not found (404)</h1>` | no such compound — permanent |
| OPSIN | `404` with `{"status": "FAILURE", "message": "…"}` | the answer, and the only place the reason is |

Retrying the CACTUS case costs four requests and several seconds of back-off to
learn an answer that will not change; treating OPSIN's 404 as absence throws
away the explanation of *which part of the name* could not be parsed. So the
policy is fixed and the *reading* is pluggable: a client passes a `classify`
callback mapping a response to an [`Outcome`][provesid.http.Outcome].
`provesid.pubchem.pugrest_classify` and `pugview_classify` read PubChem's fault
codes, `provesid.resolver.nci_classify` reads CACTUS's body, and
`provesid.opsin.opsin_classify` treats OPSIN's 404 as a success. Services that
use status codes honestly — ChEBI, CAS Common Chemistry — get
`default_classify` and need to do nothing.

The two PubChem classifiers differ in one place: a `400` that carries no fault
code. PUG-View takes the heading as a query parameter, so a `400` means "no
such heading" — absence. PUG-REST takes its whole query in the URL path, so a
`400` means the path was wrong, and the caller needs PubChem's own explanation
of which property name it misspelled.

## Retry policy

- Retried: a `RETRY` verdict, request timeouts, and connection errors.
- Not retried: `ABSENT` (the record does not exist) and `FATAL` (a malformed
  request, a rejected key) — asking again cannot change either.
- A request is made at most `max_retries + 1` times.
- The wait is `backoff * 2 ** attempt`, capped at `max_backoff` — unless the
  service sent a `Retry-After`, in seconds or as an HTTP date, in which case
  that wins, also capped.
- Pacing applies to retries too: a service already shedding load is not asked
  again faster than a healthy one.
- `max_elapsed`, when set, caps the *total* waiting: "do not make the caller
  wait longer than this". Both PubChem clients set a 10-second budget
  (`provesid.pubchem.RETRY_WAIT_BUDGET`), which leaves the cheap curve intact —
  1 + 2 + 4 for a transient 500 or a timeout — while declining PubChem's
  `Retry-After: 30`.

  That last part is deliberate. PubChem sends that header when it has throttled
  or blacklisted an IP, and such a block does not lift in thirty seconds:
  measured on 2026-09-19, waiting the full thirty and asking again returned the
  same 503. The wait buys nothing while every call pays it, which for a caller
  resolving a thousand names is hours instead of an immediate "you are blocked".
  A caller who does want to wait a throttle out raises the budget, which is the
  right setting for an unattended bulk job:

  ```python
  api._http.max_elapsed = 180
  ```

## Pacing is per host, not per object

A client's `min_interval` is its own promise about how fast it will ask. The
clock it measures that promise against belongs to the **host**, because the
limit being respected does too: PubChem publishes five requests per second *per
IP*, not per Python object. A `PubChemAPI` and a `PubChemView` in one process —
the ordinary way to use this package — would each have kept their own clock and
together asked twice as fast as PubChem allows.

So a client names the service it shares its pacing with, and
[`host_limiter`][provesid.http.host_limiter] hands every client aimed at that
host the same [`RateLimiter`][provesid.http.RateLimiter]:

```python
>>> from provesid.pubchem import PubChemAPI
>>> from provesid.pubchemview import PubChemView
>>> PubChemAPI()._http.limiter is PubChemView()._http.limiter
True
```

`last_request_time` stays per client — it says when *that* client last asked —
while the wait is measured against the shared clock. A client given no
`pace_host` paces alone, which is what a stub wants.

The effective rate for a host is set by its most impatient client: sharing the
clock stops two clients doubling a limit, but it does not stop one client
configured with `min_interval=0.01` from exceeding it by itself.

## Sessions

`HTTPClient` will make its calls through a `requests.Session` when given one,
for connection pooling and persistent headers. ChEBI uses this: its
`User-Agent` and `Accept` headers live on the session, and an ontology walk
reuses one connection across many small requests. A client with no session
calls `requests.get` and `requests.post` through the module, which is what lets
a test stub them.

::: provesid.http
