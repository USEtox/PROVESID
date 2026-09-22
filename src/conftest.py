"""
Run the docstring examples as doctests::

    pytest --doctest-modules src/provesid

This file lives beside the package rather than inside it, so it is never
installed. ``testpaths`` is ``tests``, so the ordinary suite never collects the
examples; they are checked only when asked for.

Examples that need the network carry ``# doctest: +SKIP``, and this file makes
that a rule rather than a convention: every outbound connection is refused
while the examples run, so an example that forgets the marker fails here
instead of depending on a service being up. Set
``PROVESID_DOCTEST_ALLOW_NETWORK=1`` to lift the block when re-recording an
online example's output.

Examples that read an offline database run against the installed copy, and a
module whose database is not installed is skipped here rather than
downloaded: a doctest run must never start a multi-gigabyte download.
"""

import os
import socket
import tempfile

import pytest

# A stray analysis script that reads a CSV at import time; not package code.
collect_ignore = ["provesid/data"]

# Several examples clear or fill the cache. Keep them out of the user's.
os.environ.setdefault(
    "PROVESID_CACHE_DIR", os.path.join(tempfile.gettempdir(), "provesid_doctest_cache")
)

# The provesid.config examples store and remove an API key. Point them at a
# throwaway directory, unconditionally: an earlier example that ran against
# ~/.config/provesid replaced a real CAS key with its placeholder.
_CONFIG_SANDBOX = tempfile.mkdtemp(prefix="provesid_doctest_config_")
os.environ["XDG_CONFIG_HOME"] = _CONFIG_SANDBOX
os.environ["APPDATA"] = _CONFIG_SANDBOX

# Every example is either offline or marked ``# doctest: +SKIP``. Hold the
# package to that: refuse the connection rather than let an unmarked example
# ask a live service, which passes or fails with the service's mood. Three
# examples reached PubChem and ChEBI unmarked until this went in.
_ALLOW_NETWORK = os.environ.get("PROVESID_DOCTEST_ALLOW_NETWORK") == "1"
_real_connect = socket.socket.connect


def _refuse_connect(self, address, *args, **kwargs):
    """Let a loopback connection through; refuse everything else."""
    host = address[0] if isinstance(address, tuple) else address
    if isinstance(host, str) and (host.startswith("127.") or host in ("::1", "localhost")):
        return _real_connect(self, address, *args, **kwargs)
    raise RuntimeError(
        f"this example asked {host} for something: the doctest run is offline. "
        "Mark it '# doctest: +SKIP', or set PROVESID_DOCTEST_ALLOW_NETWORK=1 "
        "to re-record its output."
    )


if not _ALLOW_NETWORK:
    socket.socket.connect = _refuse_connect


# The offline databases each module's examples read.
_DATASETS_BY_MODULE = {
    "provesid": ["pubchem", "comptox", "chebi", "chembl"],
    "provesid.pubchem_id": ["pubchem"],
    "provesid.comptox": ["comptox"],
    "provesid.chebi_sdf": ["chebi"],
    "provesid.chembl": ["chembl"],
    "provesid.zeropm": ["zeropm"],
    "provesid.search": ["pubchem", "comptox", "chebi", "chembl"],
    "provesid.sqlite_client": ["pubchem"],
    "provesid.sources": ["pubchem", "comptox", "chebi", "chembl"],
}


def pytest_collection_modifyitems(config, items):
    """Skip the examples of a module whose database is not installed."""
    from _pytest.doctest import DoctestItem

    from provesid import datasets

    absent = {}
    for item in items:
        if not isinstance(item, DoctestItem):
            continue
        module = item.dtest.globs.get("__name__", "")
        needed = _DATASETS_BY_MODULE.get(module, [])
        if module not in absent:
            absent[module] = [name for name in needed if not datasets.is_present(name)]
        if absent[module]:
            item.add_marker(pytest.mark.skip(
                reason="dataset not installed: " + ", ".join(absent[module])
            ))
