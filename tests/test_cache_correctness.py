"""
Tests for cache-key stability and for never caching a failed lookup.

These are offline unit tests: no test here touches the network. They guard two
defects that made the persistent cache useless and, worse, made transient
network errors permanent:

1. The cache key used to be built with ``json.dumps(..., default=str)`` over the
   call's arguments. For a bound method the first argument is ``self``, whose
   default ``str()`` embeds a memory address, so every instance — and every
   process — produced a different key and no persistent entry was ever reused.
2. ``@cached`` stored whatever a function returned, including the
   ``{'success': False, 'error': ...}`` dicts and empty lists that the clients
   hand back after an HTTP 429/503 or a timeout.
"""

import subprocess
import sys
from pathlib import Path

import pytest

from provesid.cache import (
    CacheManager,
    cached,
    is_failure_result,
    stable_key_part,
)
from provesid import PubChemAPI, PubChemView, NCIChemicalIdentifierResolver


REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def manager(tmp_path):
    """A CacheManager writing into an isolated temporary directory."""
    return CacheManager(cache_dir=str(tmp_path / "cache"), service_name="test")


# --------------------------------------------------------------------------
# 1. Key stability
# --------------------------------------------------------------------------

@pytest.mark.unit
@pytest.mark.parametrize("client_cls", [PubChemAPI, PubChemView, NCIChemicalIdentifierResolver])
def test_key_is_identical_for_two_client_instances(manager, client_cls):
    """Two clients of the same class and endpoint must share a cache key."""
    first, second = client_cls(), client_cls()
    assert first is not second
    assert manager._get_cache_key("f", (first, 2244), {}) == \
           manager._get_cache_key("f", (second, 2244), {})


@pytest.mark.unit
def test_key_does_not_contain_a_memory_address(manager):
    """The normalised argument for a client must not be its default repr."""
    part = stable_key_part(PubChemView())
    assert "0x" not in repr(part)
    assert part[:2] == ["provesid.pubchemview.PubChemView",
                        "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"]


@pytest.mark.unit
def test_client_key_carries_its_cache_schema_version():
    """
    A change to the *shape* of a cached result has to change the key.

    Pickle stores an instance's ``__dict__``, so an entry written before a
    dataclass gained a field restores without that attribute and every caller
    reading it raises. The version in the key is what makes such an entry
    unreachable instead.
    """
    view = PubChemView()
    assert view.CACHE_SCHEMA_VERSION in view.__cache_key__()

    bumped = PubChemView()
    bumped.CACHE_SCHEMA_VERSION = view.CACHE_SCHEMA_VERSION + 1
    assert stable_key_part(view) != stable_key_part(bumped)


@pytest.mark.unit
def test_key_is_stable_across_processes(tmp_path):
    """A key computed in a fresh interpreter must match this one."""
    script = (
        "import sys; sys.path.insert(0, %r)\n"
        "from provesid import PubChemView\n"
        "from provesid.cache import CacheManager\n"
        "cm = CacheManager(cache_dir=%r, service_name='test')\n"
        "print(cm._get_cache_key('f', (PubChemView(), 2244), {}))\n"
        % (str(REPO_ROOT / "src"), str(tmp_path / "cache"))
    )
    out = subprocess.run([sys.executable, "-c", script], capture_output=True,
                         text=True, check=True, cwd=str(REPO_ROOT))
    manager = CacheManager(cache_dir=str(tmp_path / "cache"), service_name="test")
    expected = manager._get_cache_key("f", (PubChemView(), 2244), {})
    assert out.stdout.strip() == expected


@pytest.mark.unit
def test_different_endpoint_gets_a_different_key(manager):
    """__cache_key__ keeps results from two endpoints apart."""
    default = PubChemView()
    mirror = PubChemView(base_url="https://mirror.example.org/pug_view")
    assert manager._get_cache_key("f", (default, 1), {}) != \
           manager._get_cache_key("f", (mirror, 1), {})


@pytest.mark.unit
def test_tuple_and_list_arguments_agree(manager):
    """Normalisation makes an equivalent tuple and list hash the same."""
    assert manager._get_cache_key("f", (("a", "b"),), {}) == \
           manager._get_cache_key("f", (["a", "b"],), {})


@pytest.mark.unit
def test_distinct_arguments_still_get_distinct_keys(manager):
    """Normalisation must not collapse genuinely different calls."""
    keys = {
        manager._get_cache_key("f", (2244,), {}),
        manager._get_cache_key("f", (702,), {}),
        manager._get_cache_key("g", (2244,), {}),
        manager._get_cache_key("f", (2244,), {"fmt": "JSON"}),
        manager._get_cache_key("f", (2244,), {"fmt": "CSV"}),
    }
    assert len(keys) == 5


# --------------------------------------------------------------------------
# 2. A cached None is not a cache miss
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_cached_none_round_trips(manager):
    """A stored None must be reported as a hit, not repeated forever."""
    manager.set("f", (1,), {}, None)
    manager._memory_cache.clear()  # force the disk path
    found, value = manager.get("f", (1,), {})
    assert found is True
    assert value is None


@pytest.mark.unit
def test_absent_entry_is_a_miss(manager):
    """An entry that was never written is a miss."""
    found, value = manager.get("f", ("never-written",), {})
    assert found is False
    assert value is None


# --------------------------------------------------------------------------
# 3. Failures are never stored
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_is_failure_result_recognises_the_in_band_convention():
    """The predicate accepts only the documented failure shape."""
    assert is_failure_result({"success": False, "error": "Server busy"})
    assert not is_failure_result({"success": True, "CID": 2244})
    assert not is_failure_result({"CID": 2244})
    assert not is_failure_result([])
    assert not is_failure_result(None)


@pytest.mark.unit
def test_failure_dict_is_returned_but_not_cached(tmp_path, monkeypatch):
    """A transient failure must not be remembered as the answer."""
    import provesid.cache as cache_module
    probe = CacheManager(cache_dir=str(tmp_path / "probe"), service_name="probe")
    monkeypatch.setitem(cache_module._service_caches, "probe", probe)

    calls = []

    @cached(service="probe")
    def flaky(cid):
        calls.append(cid)
        if len(calls) == 1:
            return {"success": False, "cid": cid, "error": "Server busy"}
        return {"success": True, "cid": cid, "MolecularWeight": "180.16"}

    first = flaky(2244)
    assert first["success"] is False

    second = flaky(2244)
    assert second["success"] is True, "the failure was served from the cache"
    assert calls == [2244, 2244]

    # The success is cached, so a third call does not hit the function.
    assert flaky(2244)["MolecularWeight"] == "180.16"
    assert calls == [2244, 2244]


@pytest.mark.unit
def test_skip_if_prevents_storage(tmp_path, monkeypatch):
    """skip_if covers failures that are not success:False dicts."""
    import provesid.cache as cache_module
    probe = CacheManager(cache_dir=str(tmp_path / "probe2"), service_name="probe2")
    monkeypatch.setitem(cache_module._service_caches, "probe2", probe)

    calls = []

    @cached(service="probe2", skip_if=lambda result: result == [])
    def sometimes_empty(cid):
        calls.append(cid)
        return [] if len(calls) == 1 else ["aspirin"]

    assert sometimes_empty(2244) == []
    assert sometimes_empty(2244) == ["aspirin"]
    assert calls == [2244, 2244]


@pytest.mark.unit
def test_raised_exception_is_not_cached(tmp_path, monkeypatch):
    """An exception propagates and leaves nothing behind."""
    import provesid.cache as cache_module
    probe = CacheManager(cache_dir=str(tmp_path / "probe3"), service_name="probe3")
    monkeypatch.setitem(cache_module._service_caches, "probe3", probe)

    calls = []

    @cached(service="probe3")
    def raises_once(cid):
        calls.append(cid)
        if len(calls) == 1:
            raise RuntimeError("Server busy")
        return {"success": True, "cid": cid}

    with pytest.raises(RuntimeError):
        raises_once(2244)
    assert raises_once(2244)["success"] is True
    assert calls == [2244, 2244]


@pytest.mark.unit
def test_use_cache_false_still_writes_a_success(tmp_path, monkeypatch):
    """use_cache=False means "do not read"; a success is still stored."""
    import provesid.cache as cache_module
    probe = CacheManager(cache_dir=str(tmp_path / "probe4"), service_name="probe4")
    monkeypatch.setitem(cache_module._service_caches, "probe4", probe)

    calls = []

    @cached(service="probe4")
    def counted(cid):
        calls.append(cid)
        return {"success": True, "n": len(calls)}

    counted(2244, use_cache=False)
    assert counted(2244)["n"] == 1
    assert calls == [2244]
