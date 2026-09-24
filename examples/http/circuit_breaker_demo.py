"""
See the circuit breaker stop PubChem requests while PubChem has throttled
this machine, without going online.

PubChem answers an IP it has throttled or blacklisted with ``HTTP 503`` and
``Retry-After: 30``. Before the circuit breaker, every call found this out for
itself: a batch of a thousand names cost a thousand refused requests, and each
one told PubChem the block was still needed.

A ``Retry-After`` is now recorded on the *host's* shared clock
(``provesid.http.host_limiter``). Every client aimed at that host,
``PubChemAPI`` and ``PubChemView`` alike, then:

- waits the hold out, when it is no longer than the client would wait anyway;
- or fails at once **without a request**, when it is longer.

This script stubs ``requests.get`` with a throttled PubChem, so it runs
offline and never touches the real service.

Run with::

    uv run python examples/http/circuit_breaker_demo.py
"""

import logging
import time

import requests

from provesid import PubChemAPI, PubChemView, release_holds
from provesid.http import host_limiter

logging.basicConfig(level=logging.ERROR)

sent = []


class ThrottledResponse:
    """What PubChem sends a throttled IP."""

    status_code = 503
    headers = {"Retry-After": "30"}
    text = '{"Fault": {"Code": "PUGREST.ServerBusy"}}'
    content = text.encode()

    def json(self):
        return {"Fault": {"Code": "PUGREST.ServerBusy"}}


def throttled_pubchem(url, **kwargs):
    sent.append(url)
    return ThrottledResponse()


requests.get = throttled_pubchem

api = PubChemAPI()
view = PubChemView(use_cache=False)

print("1. The first call is refused, and learns the hold:")
start = time.perf_counter()
try:
    api.get_compound_synonyms(2244)
except Exception as exc:
    print(f"   {type(exc).__name__}: {str(exc)[:90]}...")
print(f"   requests sent so far: {len(sent)}, "
      f"{(time.perf_counter() - start) * 1e3:.1f} ms")

limiter = host_limiter("https://pubchem.ncbi.nlm.nih.gov")
print(f"\n   pubchem.ncbi.nlm.nih.gov is held for another "
      f"{limiter.held_for():.0f} s")

print("\n2. A thousand more calls, split over both PubChem clients:")
start = time.perf_counter()
failures = 0
for cid in range(1, 1001):
    try:
        if cid % 2:
            api.get_compound_synonyms(cid)
        else:
            view.get_property(cid, "Boiling Point")
    except Exception:
        failures += 1
print(f"   {failures} failed, requests sent so far: {len(sent)}, "
      f"{(time.perf_counter() - start) * 1e3:.1f} ms in total")

print("\n3. Released by hand, e.g. after a network change, PubChem is asked again:")
release_holds()
try:
    api.get_compound_synonyms(2244)
except Exception:
    pass
print(f"   requests sent so far: {len(sent)}")
