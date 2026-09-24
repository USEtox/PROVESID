#!/usr/bin/env python3
"""
Where the PROVESID cache lives, how to address one service, and how a version
bump retires stale entries.

Three things this demonstrates, all offline --- no network call is made:

1. **The cache persists.** It lives in the per-user cache directory, not the
   system temp directory, so it is still there after a reboot.
2. **One cache per service.** Every cache function takes ``service=``; the
   PubChem cache can be cleared without touching OPSIN's.
3. **A version in every key.** Bumping ``@cached(version=N)``, a client's
   ``CACHE_SCHEMA_VERSION`` or ``CACHE_KEY_VERSION`` makes entries of the
   previous *shape* unreachable rather than deserialising them into something
   the caller cannot read.

The demo redirects the cache into a throwaway directory via
``PROVESID_CACHE_DIR`` so it leaves nothing behind.

Usage:
    python examples/cache/cache_layout_and_versioning_demo.py
"""

import os
import tempfile

# Redirect the cache before provesid is imported, so this demo never touches
# the real one. In your own code you would simply leave the default alone.
_SANDBOX = tempfile.TemporaryDirectory(prefix="provesid-cache-demo-")
os.environ["PROVESID_CACHE_DIR"] = _SANDBOX.name

import provesid  # noqa: E402
from provesid.cache import cached, get_service_cache  # noqa: E402
from provesid.utils import user_cache_path  # noqa: E402


def show_where_the_cache_lives() -> None:
    """Section 1: a real cache directory, overridable by one variable."""
    print("1. Where the cache lives")
    print("-" * 60)
    print(f"   this run (PROVESID_CACHE_DIR): {user_cache_path(ensure_exists=False)}")

    # What it would be without the override.
    override = os.environ.pop("PROVESID_CACHE_DIR")
    print(f"   the default on this platform:  {user_cache_path(ensure_exists=False)}")
    os.environ["PROVESID_CACHE_DIR"] = override

    print("\n   Before, the default was the system temp directory, which most")
    print("   Linux systems clear on boot --- so a 'persistent' cache lasted")
    print("   until the next reboot. Datasets live elsewhere again")
    print("   (PROVESID_DATA_DIR), because those are not disposable.\n")


def show_one_cache_per_service() -> None:
    """Section 2: service= replaces fourteen hand-written functions."""
    print("2. One cache per service")
    print("-" * 60)
    print(f"   services: {', '.join(provesid.CACHE_SERVICES)}\n")

    # Put one entry in two different services.
    get_service_cache("pubchem").set("demo", (2244,), {}, "aspirin")
    get_service_cache("opsin").set("demo", ("ethanol",), {}, "CCO")

    for name, info in provesid.get_all_cache_info().items():
        print(f"   {name:12s} {info['disk_entries']:3d} entries  "
              f"{info['cache_directory']}")

    print("\n   Clearing PubChem leaves OPSIN alone:")
    provesid.clear_cache(service="pubchem")
    print(f"   pubchem: {get_service_cache('pubchem').get('demo', (2244,), {})}")
    print(f"   opsin:   {get_service_cache('opsin').get('demo', ('ethanol',), {})}")

    print("\n   A typo raises rather than quietly using the global cache:")
    try:
        provesid.clear_cache(service="pubchemm")
    except ValueError as exc:
        print(f"   ValueError: {exc}\n")


def show_version_retires_stale_shapes() -> None:
    """Section 3: the cache cannot detect a shape change, so you declare it."""
    print("3. A version in every key")
    print("-" * 60)

    calls = []

    @cached(service="pubchem", version=1)
    def summary(cid):
        calls.append(cid)
        return {"cid": cid}

    print(f"   v1 first call:  {summary(2244)}  (function ran: {len(calls)}x)")
    print(f"   v1 second call: {summary(2244)}  (function ran: {len(calls)}x)")

    # The next release adds a field to what this function returns. Without the
    # bump, the v1 entry above would be served to callers reading ['parsed'].
    @cached(service="pubchem", version=2)
    def summary(cid):  # noqa: F811 - deliberately the next release's version
        calls.append(cid)
        return {"cid": cid, "parsed": True}

    print(f"   v2 first call:  {summary(2244)}  (function ran: {len(calls)}x)")
    print("\n   The v1 entry is still on disk; it is simply unreachable.")
    print("   Three versions, bump the narrowest that covers the change:")
    print("     @cached(version=N)          one function")
    print("     Client.CACHE_SCHEMA_VERSION one client, via __cache_key__")
    print(f"     CACHE_KEY_VERSION           everything (now "
          f"{provesid.CACHE_KEY_VERSION})")

    from provesid import PubChemView
    print(f"\n   PubChemView.__cache_key__() -> {PubChemView().__cache_key__()}")


def main() -> None:
    print("PROVESID cache: layout, services and versioning")
    print("=" * 60 + "\n")
    show_where_the_cache_lives()
    show_one_cache_per_service()
    show_version_retires_stale_shapes()
    print("\nDone. The sandbox cache directory is removed on exit.")


if __name__ == "__main__":
    try:
        main()
    finally:
        _SANDBOX.cleanup()
