"""
Run the docstring examples as doctests::

    pytest --doctest-modules src/provesid

This file lives beside the package rather than inside it, so it is never
installed. ``testpaths`` is ``tests``, so the ordinary suite never collects the
examples; they are checked only when asked for.

Examples that need the network carry ``# doctest: +SKIP``. Examples that read
an offline database run against the installed copy, and a module whose
database is not installed is skipped here rather than downloaded: a doctest
run must never start a multi-gigabyte download.
"""

import os
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

# The offline databases each module's examples read.
_DATASETS_BY_MODULE = {
    "provesid": ["pubchem", "comptox", "chebi", "chembl"],
    "provesid.pubchem_id": ["pubchem"],
    "provesid.comptox": ["comptox"],
    "provesid.chebi_sdf": ["chebi"],
    "provesid.chembl": ["chembl"],
    "provesid.zeropm": ["zeropm"],
    "provesid.search": ["pubchem", "comptox", "chebi", "chembl"],
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
