# HTTP Transport

Every PROVESID web-API client makes its requests through one shared transport.
It is the single place that decides how long to wait between requests, when a
failure is worth retrying, and how long to back off before asking again.

## Why it exists

Each client used to carry its own copy of that logic. The four copies had
drifted apart: none of them honoured `Retry-After`, only one retried at all,
and three modules had no rate limiting whatsoever. Centralising it is what
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

## Classifying a response

The one thing services genuinely disagree about is what a response *means*,
and two of the services this package talks to answer a permanent condition
with a status code that says otherwise:

| Service | Says | Means |
|---|---|---|
| PubChem | `404` with `{"Fault": {"Code": "PUGVIEW.ServerBusy"}}` | ask again later |
| PubChem | `400` with `PUGVIEW.BadRequest` | no such heading — permanent |
| NCI CACTUS | `500` with body `<h1>Page not found (404)</h1>` | no such compound — permanent |

Retrying the last one costs four requests and several seconds of back-off to
learn an answer that will not change. So the policy is fixed and the *reading*
is pluggable: a client passes a `classify` callback mapping a response to an
[`Outcome`][provesid.http.Outcome]. `provesid.pubchem.pubchem_classify` reads
PubChem's fault codes; `provesid.resolver.nci_classify` reads CACTUS's body.
Services that use status codes honestly get `default_classify` and need to do
nothing.

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

::: provesid.http
