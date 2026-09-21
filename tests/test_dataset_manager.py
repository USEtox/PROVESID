"""
Tests for the dataset manager --- ``status``/``plan``/``fetch``/``remove`` and
the ``Search(datasets=...)`` policy.

The defect these cover is the one a new user met first: ``Search("cas")`` on a
clean machine constructed four clients that each default to
``auto_download=True``, so one CAS lookup fetched ~32 GB without a word.  The
behaviour that matters is therefore mostly *absence* --- that nothing is
downloaded, that a missing dataset is reported rather than fixed behind the
caller's back, and that the message names the call that would fix it.

Every test passes an explicit ``data_dir``: the manager's whole job is to
report on a directory, and a test that used the real one would report on
whatever the developer happens to have installed.  The files it creates are
empty or near-empty --- the manager reads names and sizes, never contents.
"""

import inspect
import os

import pytest

from provesid import datasets
from provesid.datasets import DATASETS, MissingDatasetError
from provesid.search import Search


def touch(directory, name, size=0):
    """Create a file of ``size`` bytes and return its path."""
    path = os.path.join(str(directory), name)
    with open(path, "wb") as handle:
        handle.write(b"\0" * size)
    return path


# ── The registry itself ───────────────────────────────────────────────────────

class TestRegistry:

    def test_names_match_the_sources_search_queries(self):
        """The two tables are written separately and must not drift apart.

        ``Search`` names its sources in ``_ALL_SOURCE_KEYS`` and the registry
        names the datasets they read.  A dataset added to one and not the other
        would go unreported by ``status`` or unreachable by ``fetch``, with
        nothing to say so.
        """
        assert set(DATASETS) == set(Search._ALL_SOURCE_KEYS)
        assert set(datasets.DEFAULT_DATASETS) == set(Search._DEFAULT_SOURCE_KEYS)

    def test_every_entry_is_described_and_sized(self):
        for name, dataset in DATASETS.items():
            assert dataset.name == name
            assert dataset.title and dataset.role and dataset.source
            assert dataset.download_bytes > 0
            assert dataset.resident_bytes > 0
            assert dataset.peak_bytes >= dataset.resident_bytes

    def test_chembl_peak_exceeds_what_it_installs(self):
        """ChEMBL is built out of files far larger than the one it leaves.

        The archive and the full release both exist before the extract is
        built, so a laptop with 28 GB free still cannot install a dataset that
        ends up costing 2.4 GiB --- exactly the surprise the manager exists to
        prevent.
        """
        chembl = DATASETS["chembl"]
        assert chembl.peak_bytes > 10 * chembl.resident_bytes
        assert all(DATASETS[name].peak_bytes == DATASETS[name].resident_bytes
                   for name in ("comptox", "zeropm"))

    def test_pubchem_is_sized_as_the_ftp_build_the_default_route_runs(self):
        """``fetch`` builds PubChem from FTP; the registry must say what that costs.

        Eight files are downloaded and deleted one at a time, so the transfer
        is several times what is kept, and the peak is the database plus the
        largest file --- not the sum of the files.
        """
        from provesid import PubChemID
        from provesid.pubchem_ftp import SOURCE_FILES

        assert inspect.signature(PubChemID.__init__).parameters["source"].default == "ftp"
        pubchem = DATASETS["pubchem"]
        assert pubchem.download_bytes == sum(source.approx_bytes for source in SOURCE_FILES)
        assert pubchem.download_bytes > 3 * pubchem.resident_bytes
        assert pubchem.peak_bytes == (pubchem.resident_bytes
                                      + max(source.approx_bytes for source in SOURCE_FILES))
        assert pubchem.peak_bytes < pubchem.resident_bytes + pubchem.download_bytes

    def test_chembl_is_sized_as_the_extract_the_default_route_installs(self):
        """The registry must describe the route ``fetch`` actually takes.

        ``fetch`` constructs ``CheMBL`` with its default ``source``, which
        compacts the release and deletes it, so what lands on disk is the
        ~2.4 GiB extract and not the 27.7 GiB download it was built from. If
        that default ever moves, these sizes become a lie told before a 5.8 GB
        transfer.
        """
        from provesid import CheMBL

        assert inspect.signature(CheMBL.__init__).parameters["source"].default \
            == "sqlite"
        chembl = DATASETS["chembl"]
        assert chembl.resident_bytes < chembl.download_bytes
        assert 2 * 1024 ** 3 < chembl.resident_bytes < 4 * 1024 ** 3

    def test_unknown_name_is_rejected_with_the_known_ones(self):
        with pytest.raises(KeyError, match="pubchem"):
            datasets.status("pubchem_id")

    def test_human_bytes(self):
        assert datasets.human_bytes(0) == "0 B"
        assert datasets.human_bytes(512) == "512 B"
        assert datasets.human_bytes(1536) == "1.5 KiB"
        assert datasets.human_bytes(2322595840) == "2.2 GiB"

    def test_fetch_command_is_copy_pasteable(self):
        assert datasets.fetch_command("chembl") == "provesid.datasets.fetch('chembl')"
        assert (datasets.fetch_command(["chebi", "pubchem"])
                == "provesid.datasets.fetch(['pubchem', 'chebi'])")


# ── status ────────────────────────────────────────────────────────────────────

class TestStatus:

    def test_empty_directory(self, tmp_path):
        frame = datasets.status(data_dir=tmp_path)
        assert list(frame["dataset"]) == datasets.dataset_names()
        assert not frame["present"].any()
        assert frame["bytes"].sum() == 0
        assert (frame["path"] == "").all()
        assert frame.attrs["data_dir"] == str(tmp_path)

    def test_counts_the_derived_index_with_the_dataset(self, tmp_path):
        """ChEBI's index is a fifth of what ChEBI costs on disk.

        It is built on first use rather than downloaded, so a status table that
        counted only the downloaded file would understate ChEBI by ~75 MiB.
        """
        touch(tmp_path, "chebi.sdf", 1000)
        touch(tmp_path, "chebi.sdf.index.pkl", 200)
        row = datasets.status("chebi", data_dir=tmp_path).iloc[0]
        assert row["present"] and row["files"] == 2 and row["bytes"] == 1200
        assert row["path"].endswith("chebi.sdf")

    def test_a_partial_download_takes_space_without_being_present(self, tmp_path):
        """An interrupted 2 GB download is reported, and is not the dataset.

        This is the row that explains a full disk after a download that was
        cancelled, and treating it as an installed dataset would be worse than
        not reporting it at all.
        """
        touch(tmp_path, "pubchem_id.db.part", 4096)
        row = datasets.status("pubchem", data_dir=tmp_path).iloc[0]
        assert not row["present"]
        assert row["bytes"] == 4096
        assert row["path"] == ""

    def test_the_chembl_extract_wins_over_the_release_it_came_from(self, tmp_path):
        """Both answer identically, and one costs a tenth of the page cache.

        ``CheMBL`` opens the extract in that case, so the status table has to
        name the same file or it is describing a database nobody will open.
        """
        touch(tmp_path, "chembl_36.db", 100)
        touch(tmp_path, "chembl_36_provesid.db", 10)
        row = datasets.status("chembl", data_dir=tmp_path).iloc[0]
        assert row["path"].endswith("chembl_36_provesid.db")
        assert row["release"] == "36 (extract)"
        assert row["bytes"] == 110          # both are on disk and both count

    def test_a_newer_release_beats_a_smaller_extract(self, tmp_path):
        touch(tmp_path, "chembl_36_provesid.db", 10)
        touch(tmp_path, "chembl_37.db", 100)
        row = datasets.status("chembl", data_dir=tmp_path).iloc[0]
        assert row["path"].endswith("chembl_37.db")
        assert row["release"] == "37"

    def test_zeropm_reports_the_newest_version(self, tmp_path):
        touch(tmp_path, "zeropm-v0-0-3.sqlite", 10)
        touch(tmp_path, "zeropm-v0-0-4.sqlite", 10)
        row = datasets.status("zeropm", data_dir=tmp_path).iloc[0]
        assert row["release"] == "0.0.4"
        assert row["path"].endswith("zeropm-v0-0-4.sqlite")

    def test_release_is_empty_where_the_filename_does_not_say(self, tmp_path):
        touch(tmp_path, "pubchem_id.db", 10)
        assert datasets.status("pubchem", data_dir=tmp_path).iloc[0]["release"] == ""

    def test_missing_lists_only_what_is_absent(self, tmp_path):
        touch(tmp_path, "pubchem_id.db", 10)
        assert datasets.missing(["pubchem", "chembl"], tmp_path) == ["chembl"]
        assert datasets.is_present("pubchem", tmp_path)
        assert not datasets.is_present("chembl", tmp_path)


# ── plan ──────────────────────────────────────────────────────────────────────

class TestPlan:

    def test_totals_cover_only_what_is_missing(self, tmp_path):
        touch(tmp_path, "pubchem_id.db", 10)
        frame = datasets.plan(["pubchem", "chebi"], data_dir=tmp_path)

        actions = dict(zip(frame["dataset"], frame["action"]))
        assert actions == {"pubchem": "present", "chebi": "download"}
        assert frame.attrs["total_download_bytes"] == DATASETS["chebi"].download_bytes
        assert frame.attrs["total_resident_bytes"] == DATASETS["chebi"].resident_bytes

    def test_a_clean_machine_is_told_the_whole_number(self, tmp_path):
        """The 32 GB nobody was asked about, on screen before anything moves."""
        frame = datasets.plan(datasets.DEFAULT_DATASETS, data_dir=tmp_path)
        expected = sum(DATASETS[name].download_bytes
                       for name in datasets.DEFAULT_DATASETS)
        assert frame.attrs["total_download_bytes"] == expected
        assert (frame["action"] == "download").all()
        assert frame.attrs["peak_bytes"] > frame.attrs["total_resident_bytes"]

    def test_force_plans_a_redownload_of_what_is_present(self, tmp_path):
        touch(tmp_path, "pubchem_id.db", 10)
        frame = datasets.plan("pubchem", data_dir=tmp_path, force=True)
        assert frame.iloc[0]["action"] == "download"

    def test_nothing_to_do_is_zero_not_empty(self, tmp_path):
        touch(tmp_path, "pubchem_id.db", 10)
        frame = datasets.plan("pubchem", data_dir=tmp_path)
        assert len(frame) == 1
        assert frame.attrs["total_download_bytes"] == 0


# ── fetch ─────────────────────────────────────────────────────────────────────

class _StubClient:
    """Stands in for a source client: records its kwargs, writes its file."""

    calls = []

    def __init__(self, filename, **kwargs):
        self.filename = filename
        _StubClient.calls.append(kwargs)
        directory = kwargs["data_dir"]
        os.makedirs(directory, exist_ok=True)
        with open(os.path.join(directory, filename), "wb") as handle:
            handle.write(b"\0" * 10)


@pytest.fixture
def stub_clients(monkeypatch):
    """Replace the real clients with stubs that create the expected file."""
    _StubClient.calls = []
    names = {"pubchem": "pubchem_id.db", "comptox": "comptox_chemicals.db",
             "chebi": "chebi.sdf", "chembl": "chembl_37.db",
             "zeropm": "zeropm-v0-0-4.sqlite"}

    def fake_client_class(name):
        filename = names[name]
        return lambda **kwargs: _StubClient(filename, **kwargs)

    monkeypatch.setattr(datasets, "_client_class", fake_client_class)
    return _StubClient


class TestFetch:

    def test_installs_a_missing_dataset_and_returns_its_path(self, tmp_path, stub_clients):
        installed = datasets.fetch("pubchem", data_dir=tmp_path)
        assert installed["pubchem"] == str(tmp_path / "pubchem_id.db")
        assert os.path.exists(installed["pubchem"])

    def test_asks_the_client_to_download(self, tmp_path, stub_clients):
        datasets.fetch("chebi", data_dir=tmp_path)
        assert stub_clients.calls == [
            {"auto_download": True, "data_dir": str(tmp_path), "redownload": False}
        ]

    def test_present_datasets_are_not_touched(self, tmp_path, stub_clients):
        """Calling fetch on a list has to be cheap and repeatable.

        Re-fetching a 27.7 GB dataset because it was named again would make
        ``fetch(["pubchem", "chembl"])`` a trap.
        """
        touch(tmp_path, "pubchem_id.db", 10)
        installed = datasets.fetch(["pubchem", "chebi"], data_dir=tmp_path)

        assert stub_clients.calls == [
            {"auto_download": True, "data_dir": str(tmp_path), "redownload": False}
        ]
        assert set(installed) == {"pubchem", "chebi"}
        assert os.path.getsize(installed["pubchem"]) == 10   # untouched

    def test_force_redownloads_a_present_dataset(self, tmp_path, stub_clients):
        touch(tmp_path, "pubchem_id.db", 999)
        datasets.fetch("pubchem", data_dir=tmp_path, force=True)
        assert stub_clients.calls[0]["redownload"] is True

    def test_announces_the_total_before_transferring(self, tmp_path, stub_clients, caplog):
        with caplog.at_level("INFO", logger="provesid.datasets"):
            datasets.fetch(["pubchem", "chebi"], data_dir=tmp_path)
        announcement = "\n".join(caplog.messages)
        assert "2 dataset(s)" in announcement
        assert "GiB" in announcement          # the size, before the first byte

    def test_unknown_name_downloads_nothing(self, tmp_path, stub_clients):
        with pytest.raises(KeyError):
            datasets.fetch(["pubchem", "nonesuch"], data_dir=tmp_path)
        assert stub_clients.calls == []
        assert os.listdir(tmp_path) == []


# ── remove ────────────────────────────────────────────────────────────────────

class TestRemove:

    def test_dry_run_deletes_nothing(self, tmp_path):
        touch(tmp_path, "chembl_36.db", 100)
        frame = datasets.remove("chembl", data_dir=tmp_path, dry_run=True)

        assert frame.attrs["freed_bytes"] == 100
        assert not frame["removed"].any()
        assert os.path.exists(tmp_path / "chembl_36.db")

    def test_removes_the_dataset_and_everything_derived_from_it(self, tmp_path):
        """A leftover index or archive reclaims a fraction of the space.

        Removing ChEBI and leaving its 75 MiB index behind would also leave
        ``status`` reporting bytes for a dataset that is gone.
        """
        touch(tmp_path, "chebi.sdf", 100)
        touch(tmp_path, "chebi.sdf.index.pkl", 20)
        touch(tmp_path, "chebi.sdf.gz.part", 5)
        touch(tmp_path, "pubchem_id.db", 50)

        frame = datasets.remove("chebi", data_dir=tmp_path)

        assert frame["removed"].all()
        assert frame.attrs["freed_bytes"] == 125
        assert os.listdir(tmp_path) == ["pubchem_id.db"]
        assert not datasets.is_present("chebi", tmp_path)

    def test_removes_every_chembl_release_in_the_directory(self, tmp_path):
        touch(tmp_path, "chembl_36.db", 100)
        touch(tmp_path, "chembl_36_provesid.db", 10)
        frame = datasets.remove("chembl", data_dir=tmp_path)
        assert frame.attrs["freed_bytes"] == 110
        assert os.listdir(tmp_path) == []

    def test_removing_what_is_not_there_is_not_an_error(self, tmp_path):
        frame = datasets.remove("chembl", data_dir=tmp_path)
        assert len(frame) == 0
        assert frame.attrs["freed_bytes"] == 0


# ── require ───────────────────────────────────────────────────────────────────

class TestRequire:

    def test_silent_when_everything_is_present(self, tmp_path):
        touch(tmp_path, "pubchem_id.db", 10)
        datasets.require("pubchem", tmp_path)

    def test_names_each_missing_dataset_its_size_and_the_fetch_call(self, tmp_path):
        touch(tmp_path, "pubchem_id.db", 10)
        with pytest.raises(MissingDatasetError) as raised:
            datasets.require(["pubchem", "chebi", "chembl"], tmp_path)

        message = str(raised.value)
        assert "pubchem" not in message.split("Install them with:")[0]
        assert "ChEBI SDF" in message and "ChEMBL" in message
        assert "provesid.datasets.fetch(['chebi', 'chembl'])" in message
        assert str(tmp_path) in message


# ── Search(datasets=...) ──────────────────────────────────────────────────────

@pytest.fixture
def no_downloads(monkeypatch):
    """Make any bulk download fail the test loudly, wherever it is called from."""
    def refuse(*args, **kwargs):
        raise AssertionError(f"download_file was called: {args[:1]}")

    for module in ("pubchem", "comptox", "zeropm", "chembl", "chebi", "datasets"):
        monkeypatch.setattr(f"provesid.{module}.download_file", refuse,
                            raising=False)


class _Recorder:
    """A source client that records how it was constructed and does nothing."""

    calls = []

    def __init__(self, **kwargs):
        _Recorder.calls.append(kwargs)


@pytest.fixture
def recording_clients(monkeypatch):
    _Recorder.calls = []
    for symbol in ("ChebiSDF", "CompToxID", "PubChemID", "CheMBL", "ZeroPM"):
        monkeypatch.setattr(f"provesid.search.{symbol}", _Recorder)
    return _Recorder


class TestSearchDatasetPolicy:

    def test_present_is_the_default(self):
        assert Search("cas", data_dir="/nonexistent").datasets == "present"

    def test_present_downloads_nothing_and_reports_what_is_missing(
        self, tmp_path, no_downloads, caplog
    ):
        """The defect in one test: a clean machine, one CAS lookup, no 32 GB.

        The search still runs --- ``Search`` already degrades to the sources it
        has --- and the log says which sources it is running without and what
        installing them would cost.
        """
        search = Search("cas", data_dir=tmp_path, show_progress=False)
        with caplog.at_level("WARNING", logger="provesid.search"):
            frame = search.search("50-00-0")

        assert os.listdir(tmp_path) == []
        assert search.sources_available == []
        assert search.sources_unavailable == list(Search._DEFAULT_SOURCE_KEYS)
        assert len(frame) == 1 and frame.iloc[0]["CASRN"] == "50-00-0"

        reported = "\n".join(caplog.messages)
        for name in Search._DEFAULT_SOURCE_KEYS:
            assert f"provesid.datasets.fetch('{name}')" in reported

    def test_present_uses_the_sources_that_are_installed(self, tmp_path, recording_clients):
        """Missing is per dataset, not all-or-nothing."""
        Search("cas", data_dir=tmp_path)._ensure_clients()
        assert len(recording_clients.calls) == len(Search._DEFAULT_SOURCE_KEYS)
        assert all(call["auto_download"] is False for call in recording_clients.calls)

    def test_auto_restores_the_old_behaviour(self, tmp_path, recording_clients):
        Search("cas", datasets="auto", data_dir=tmp_path)._ensure_clients()
        assert recording_clients.calls
        assert all(call["auto_download"] is True for call in recording_clients.calls)

    def test_required_raises_in_the_constructor(self, tmp_path, no_downloads):
        """Before the first query, not after an hour of them."""
        with pytest.raises(MissingDatasetError) as raised:
            Search("cas", datasets="required", data_dir=tmp_path)
        assert "provesid.datasets.fetch(" in str(raised.value)
        assert os.listdir(tmp_path) == []

    def test_required_is_satisfied_by_what_is_present(self, tmp_path, recording_clients):
        for name in ("pubchem_id.db", "comptox_chemicals.db", "chebi.sdf",
                     "chembl_37.db"):
            touch(tmp_path, name, 10)
        Search("cas", datasets="required", data_dir=tmp_path)

    def test_required_ignores_sources_whose_client_was_passed(self, tmp_path):
        """A client handed in has already found its data, wherever that is."""
        given = _Recorder()
        Search("cas", datasets="required", data_dir=tmp_path,
               chebi=given, comptox=given, pubchem=given, chembl=given)

    def test_required_does_not_demand_zeropm_unless_it_is_queried(self, tmp_path):
        for name in ("pubchem_id.db", "comptox_chemicals.db", "chebi.sdf",
                     "chembl_37.db"):
            touch(tmp_path, name, 10)
        Search("cas", datasets="required", data_dir=tmp_path)      # no ZeroPM

        with pytest.raises(MissingDatasetError, match="zeropm"):
            Search("cas", datasets="required", data_dir=tmp_path, use_zeropm=True)

    def test_an_unknown_policy_is_rejected(self):
        with pytest.raises(ValueError, match="datasets must be one of"):
            Search("cas", datasets="yes please")

    def test_redownload_needs_a_policy_that_downloads(self, tmp_path):
        """Silently ignoring it would hand back a stale copy with no sign."""
        with pytest.raises(ValueError, match="redownload=True"):
            Search("cas", redownload=True, data_dir=tmp_path)
        Search("cas", redownload=True, datasets="auto", data_dir=tmp_path)
