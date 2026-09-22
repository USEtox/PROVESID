"""
Tests for where the cache lives, how it is addressed, and how it is versioned.

Three defects motivate this file:

1. The default cache directory was ``tempfile.gettempdir()/provesid_cache``,
   while ``docs/advanced_caching.md`` promised the cache survived a restart.
   On most Linux systems ``/tmp`` is cleared on boot, so it did not.
2. Fourteen near-identical module-level functions (``clear_pubchem_cache``,
   ``get_cas_cache_info``, ...) said what ``clear_cache(service=...)`` and
   ``get_cache_info(service=...)`` say once. The service list is data.
3. Nothing in the key recorded the *shape* of what was stored, so a release
   that added a field to a cached dataclass served entries that unpickled into
   objects missing that attribute.

Everything here is offline.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

import provesid
from provesid import cache as cache_module
from provesid.cache import (
    CACHE_KEY_VERSION,
    CACHE_SERVICES,
    CacheManager,
    cached,
    clear_cache,
    get_all_cache_info,
    get_cache_info,
    get_service_cache,
)
from provesid.utils import user_cache_path


REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def manager_factory(tmp_path):
    """Build isolated CacheManagers under the test's own directory."""
    def _build(name):
        return CacheManager(cache_dir=str(tmp_path / name), service_name=name)
    return _build


# --------------------------------------------------------------------------
# 1. The cache directory persists
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_default_cache_dir_is_not_the_system_temp_dir(monkeypatch, tmp_path):
    """A cache under /tmp is discarded on reboot; this one is not."""
    import tempfile

    monkeypatch.delenv("PROVESID_CACHE_DIR", raising=False)
    default = user_cache_path(ensure_exists=False)
    system_temp = os.path.realpath(tempfile.gettempdir())

    assert not os.path.realpath(default).startswith(system_temp + os.sep), (
        f"the cache still defaults into the system temp directory: {default}"
    )
    assert "provesid" in default.lower()


@pytest.mark.unit
def test_cache_dir_honours_the_environment_override(monkeypatch, tmp_path):
    """PROVESID_CACHE_DIR relocates the whole cache, service dirs included."""
    monkeypatch.setenv("PROVESID_CACHE_DIR", str(tmp_path / "elsewhere"))

    assert user_cache_path(ensure_exists=False) == str(tmp_path / "elsewhere")
    assert user_cache_path("pubchem", ensure_exists=False) == str(
        tmp_path / "elsewhere" / "pubchem"
    )
    # The override is also what a freshly built manager picks up.
    manager = CacheManager(service_name="pubchem")
    assert manager.cache_dir == tmp_path / "elsewhere" / "pubchem"


@pytest.mark.unit
def test_cache_dir_is_separate_from_the_dataset_dir(monkeypatch, tmp_path):
    """Disposable responses and 2.4 GiB datasets are cleaned independently."""
    from provesid.utils import user_dataset_path

    monkeypatch.delenv("PROVESID_CACHE_DIR", raising=False)
    monkeypatch.delenv("PROVESID_DATA_DIR", raising=False)
    assert user_cache_path(ensure_exists=False) != user_dataset_path(
        ensure_exists=False
    )


@pytest.mark.unit
def test_entries_survive_a_fresh_interpreter(tmp_path):
    """The point of the whole module: write in one process, read in the next."""
    script = (
        "import os, sys\n"
        f"sys.path.insert(0, {str(REPO_ROOT / 'src')!r})\n"
        f"os.environ['PROVESID_CACHE_DIR'] = {str(tmp_path / 'persist')!r}\n"
        "from provesid.cache import cached\n"
        "@cached(service='pubchem')\n"
        "def fetch(cid):\n"
        "    print('CALLED')\n"
        "    return {'success': True, 'cid': cid}\n"
        "fetch(2244)\n"
    )
    first = subprocess.run([sys.executable, "-c", script],
                           capture_output=True, text=True, check=True)
    second = subprocess.run([sys.executable, "-c", script],
                            capture_output=True, text=True, check=True)

    assert "CALLED" in first.stdout
    assert "CALLED" not in second.stdout, (
        "the second process re-ran the function, so nothing persisted"
    )


# --------------------------------------------------------------------------
# 2. The service list is data
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_every_service_gets_its_own_directory():
    """Clearing PubChem must not touch OPSIN."""
    directories = {
        name: get_cache_info(service=name)['cache_directory']
        for name in CACHE_SERVICES
    }
    assert len(set(directories.values())) == len(CACHE_SERVICES)
    assert get_cache_info()['cache_directory'] not in directories.values()


@pytest.mark.unit
def test_service_caches_are_shared_and_built_lazily():
    """Two clients of one service share entries and one manager."""
    assert get_service_cache('opsin') is get_service_cache('opsin')
    assert get_service_cache('opsin') is not get_service_cache('pubchem')
    assert get_service_cache(None) is get_service_cache(None)


@pytest.mark.unit
def test_unknown_service_raises_rather_than_silently_using_the_global_cache():
    """A typo must not scatter entries where clear_cache never looks."""
    with pytest.raises(ValueError, match="Unknown cache service"):
        get_cache_info(service='pubchemm')
    with pytest.raises(ValueError, match="Unknown cache service"):
        clear_cache(service='pubchemm')
    with pytest.raises(ValueError, match="Unknown cache service"):
        @cached(service='not-a-service')
        def f():
            return 1


@pytest.mark.unit
def test_get_all_cache_info_covers_the_global_cache_and_every_service():
    info = get_all_cache_info()
    assert set(info) == {'global', *CACHE_SERVICES}
    assert all('cache_directory' in entry for entry in info.values())


@pytest.mark.unit
def test_clear_cache_addresses_one_service_at_a_time(monkeypatch, tmp_path):
    """clear_cache(service=...) is the whole of the old per-service family."""
    kept = CacheManager(cache_dir=str(tmp_path / "kept"), service_name="opsin")
    dropped = CacheManager(cache_dir=str(tmp_path / "dropped"),
                           service_name="pubchem")
    monkeypatch.setitem(cache_module._service_caches, "opsin", kept)
    monkeypatch.setitem(cache_module._service_caches, "pubchem", dropped)

    kept.set("f", (1,), {}, "keep me")
    dropped.set("f", (1,), {}, "drop me")

    clear_cache(service="pubchem")

    assert kept.get("f", (1,), {}) == (True, "keep me")
    assert dropped.get("f", (1,), {}) == (False, None)


@pytest.mark.unit
def test_clear_cache_all_services_clears_everything(monkeypatch, tmp_path):
    managers = {}
    for name in CACHE_SERVICES:
        manager = CacheManager(cache_dir=str(tmp_path / name), service_name=name)
        manager.set("f", (1,), {}, name)
        monkeypatch.setitem(cache_module._service_caches, name, manager)
        managers[name] = manager

    clear_cache(all_services=True)

    for name, manager in managers.items():
        assert manager.get("f", (1,), {}) == (False, None), name


@pytest.mark.unit
def test_the_per_service_function_family_is_gone():
    """One parameterised function per operation, not one per service."""
    for name in CACHE_SERVICES:
        assert not hasattr(cache_module, f"clear_{name}_cache")
        assert not hasattr(cache_module, f"get_{name}_cache_info")


@pytest.mark.unit
def test_clients_still_expose_their_own_cache_helpers():
    """The per-service functions went; the client methods that used them stay."""
    from provesid import (
        CASCommonChem,
        NCIChemicalIdentifierResolver,
        PubChemAPI,
        PubChemView,
    )
    from provesid.opsin import OPSIN

    for client, service in [
        (PubChemAPI(), 'pubchem'),
        (PubChemView(), 'pubchemview'),
        (NCIChemicalIdentifierResolver(), 'nci'),
        (CASCommonChem(), 'cas'),
        (OPSIN(), 'opsin'),
    ]:
        info = client.get_cache_info()
        assert info['cache_directory'] == get_cache_info(
            service=service)['cache_directory']


# --------------------------------------------------------------------------
# 3. A version in every key
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_the_global_key_version_is_part_of_every_key(manager_factory,
                                                     monkeypatch):
    """Bumping CACHE_KEY_VERSION retires every entry everywhere."""
    manager = manager_factory("global-version")
    before = manager._get_cache_key("f", (1,), {})

    monkeypatch.setattr(cache_module, "CACHE_KEY_VERSION", CACHE_KEY_VERSION + 1)
    after = manager._get_cache_key("f", (1,), {})

    assert before != after


@pytest.mark.unit
def test_a_function_version_retires_only_that_function(manager_factory):
    """@cached(version=N) is the narrowest of the three versions."""
    manager = manager_factory("fn-version")
    v1 = manager._get_cache_key("f", (1,), {}, version=1)
    v2 = manager._get_cache_key("f", (1,), {}, version=2)
    other = manager._get_cache_key("g", (1,), {}, version=1)

    assert v1 != v2
    assert v1 != other
    assert v1 == manager._get_cache_key("f", (1,), {}, version=1)


@pytest.mark.unit
def test_bumping_a_functions_version_makes_the_old_entry_unreachable(
        manager_factory, monkeypatch):
    """The whole point: the stale *shape* is never served after a bump."""
    manager = manager_factory("shape")
    monkeypatch.setitem(cache_module._service_caches, "shape-probe", manager)

    calls = []

    @cached(service="shape-probe", version=1)
    def summary(cid):
        calls.append(cid)
        return {"cid": cid}

    assert summary(2244) == {"cid": 2244}
    assert summary(2244) == {"cid": 2244}
    assert calls == [2244]

    # Same function, new return shape, version bumped alongside it.
    @cached(service="shape-probe", version=2)
    def summary(cid):  # noqa: F811 - deliberately the next release's version
        calls.append(cid)
        return {"cid": cid, "parsed": True}

    assert summary(2244) == {"cid": 2244, "parsed": True}, (
        "the version-1 entry was served to a caller expecting version 2"
    )
    assert calls == [2244, 2244]


@pytest.mark.unit
def test_client_schema_version_sits_between_the_other_two():
    """A client bump retires that client's entries and no others."""
    from provesid import PubChemView
    from provesid.opsin import OPSIN

    assert PubChemView.CACHE_SCHEMA_VERSION in PubChemView().__cache_key__()
    assert OPSIN.CACHE_SCHEMA_VERSION in OPSIN().__cache_key__()


@pytest.mark.unit
def test_a_cached_none_survives_an_export_import_round_trip(tmp_path):
    """`is not _MISS`, not `is not None`: None is a legitimate cached value."""
    source = CacheManager(cache_dir=str(tmp_path / "src"), service_name="rt")
    source.set("f", (1,), {}, None)

    export_file = tmp_path / "export.pkl"
    assert source.export_cache(str(export_file))

    target = CacheManager(cache_dir=str(tmp_path / "dst"), service_name="rt")
    assert target.import_cache(str(export_file))
    assert target.get("f", (1,), {}) == (True, None)


@pytest.mark.unit
def test_entries_written_by_an_earlier_process_are_counted_and_exported(tmp_path):
    """
    metadata.json is flushed every hundred writes, so a process that made two
    never recorded them. Export used to walk the metadata and so returned
    True with nothing in it; the files on disk are the record.
    """
    cache_dir = str(tmp_path / "cache")
    writer = CacheManager(cache_dir=cache_dir, service_name="rt")
    writer.set("f", (1,), {}, "one")
    writer.set("f", (2,), {}, "two")

    later = CacheManager(cache_dir=cache_dir, service_name="rt")
    assert later.get_cache_info()["disk_entries"] == 2

    export_file = tmp_path / "export.pkl"
    assert later.export_cache(str(export_file))
    target = CacheManager(cache_dir=str(tmp_path / "dst"), service_name="rt")
    assert target.import_cache(str(export_file))
    assert target.get("f", (1,), {}) == (True, "one")
    assert target.get("f", (2,), {}) == (True, "two")


@pytest.mark.unit
def test_the_new_cache_api_is_exported():
    for name in ('CACHE_SERVICES', 'CACHE_KEY_VERSION', 'clear_cache',
                 'get_cache_info', 'get_all_cache_info', 'get_cache_size',
                 'get_service_cache', 'export_cache', 'import_cache',
                 'set_cache_warning_threshold', 'enable_cache_warnings'):
        assert hasattr(provesid, name), name

