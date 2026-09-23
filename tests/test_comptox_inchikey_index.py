"""Tests for the CompTox InChIKey index, built by the first InChIKey lookup.

The downloaded database indexes DTXSID, CASRN and the preferred name only, so
each InChIKey lookup scanned the 1.2 M rows (0.2 s) and each missed skeleton
search took 0.9 s.  Each test builds a small database in the real schema, so
nothing is downloaded and the real file is never touched.
"""

import logging
import os
import shutil
import sqlite3
import stat
import threading

import pytest

import provesid.comptox
from provesid.comptox import INCHIKEY_INDEX, NAME_INDEX_TABLE, CompToxID
from provesid.sources import comptox_skeleton_search

_COLUMNS = (
    "DTXSID", "PREFERRED_NAME", "CASRN", "DTXCID", "INCHIKEY", "IUPAC_NAME", "SMILES",
    "MOLECULAR_FORMULA", "AVERAGE_MASS", "MONOISOTOPIC_MASS", "QSAR_READY_SMILES",
    "MS_READY_SMILES", "IDENTIFIER",
)

ASPIRIN = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
ASPIRIN_ANION = "BSYNRYMUTXBXSQ-UHFFFAOYSA-M"
DEHYDROACETIC_NONSTANDARD = "PGRHXDWITVMQBC-UHFFFAOYNA-N"
DEHYDROACETIC_STANDARD = "PGRHXDWITVMQBC-UHFFFAOYSA-N"

_ROWS = [
    {"DTXSID": "DTXSID5020108", "PREFERRED_NAME": "Aspirin", "CASRN": "50-78-2",
     "INCHIKEY": ASPIRIN, "IDENTIFIER": "Aspirin"},
    {"DTXSID": "DTXSID0000003", "PREFERRED_NAME": "Aspirin anion", "CASRN": "0-00-3",
     "INCHIKEY": ASPIRIN_ANION, "IDENTIFIER": "Aspirin anion"},
    {"DTXSID": "DTXSID6020014", "PREFERRED_NAME": "Dehydroacetic acid", "CASRN": "520-45-6",
     "INCHIKEY": DEHYDROACETIC_NONSTANDARD, "IDENTIFIER": "Dehydroacetic acid"},
    {"DTXSID": "DTXSID0000004", "PREFERRED_NAME": "Mixture", "CASRN": "0-00-4",
     "INCHIKEY": None, "IDENTIFIER": "Mixture"},
]


def _write_database(path):
    """A database in the real schema and with the downloaded file's indexes."""
    connection = sqlite3.connect(path)
    connection.execute(f"CREATE TABLE chemicals ({', '.join(c + ' TEXT' for c in _COLUMNS)})")
    connection.execute("CREATE INDEX idx_dtxsid ON chemicals(DTXSID)")
    connection.execute("CREATE INDEX idx_casrn ON chemicals(CASRN)")
    connection.execute("CREATE INDEX idx_preferred_name ON chemicals(PREFERRED_NAME)")
    connection.executemany(
        f"INSERT INTO chemicals VALUES ({', '.join('?' for _ in _COLUMNS)})",
        [tuple(row.get(c) for c in _COLUMNS) for row in _ROWS],
    )
    connection.commit()
    connection.close()


@pytest.fixture
def db_path(tmp_path):
    path = tmp_path / "comptox_chemicals.db"
    _write_database(path)
    return str(path)


@pytest.fixture
def db(db_path):
    with CompToxID(db_path=db_path, auto_download=False) as client:
        yield client


def _has_inchikey_index(client):
    return client.conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='index' AND name=?", (INCHIKEY_INDEX,)
    ).fetchone() is not None


def _plan(client, sql, params):
    return " ".join(row[3] for row in client.conn.execute(f"EXPLAIN QUERY PLAN {sql}", params))


# ── Building ─────────────────────────────────────────────────────────────────


@pytest.mark.unit
def test_opening_the_database_does_not_build_the_index(db):
    db.get_by_casrn("50-78-2")
    assert not _has_inchikey_index(db)


@pytest.mark.unit
def test_first_inchikey_lookup_builds_the_index_and_says_so_once(db, caplog):
    with caplog.at_level(logging.WARNING, logger="provesid.comptox"):
        assert db.get_by_inchikey(ASPIRIN)["DTXSID"] == "DTXSID5020108"
        assert db.inchikey_to_dtxsid(DEHYDROACETIC_STANDARD) == "DTXSID6020014"
    assert _has_inchikey_index(db)
    building = [r for r in caplog.records if "InChIKey index" in r.message]
    assert len(building) == 1


@pytest.mark.unit
def test_a_database_that_has_the_index_is_left_alone(db_path, caplog):
    with CompToxID(db_path=db_path, auto_download=False) as client:
        client.get_by_inchikey(ASPIRIN)
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="provesid.comptox"):
        with CompToxID(db_path=db_path, auto_download=False) as client:
            assert client.get_by_inchikey(ASPIRIN)["DTXSID"] == "DTXSID5020108"
    assert not [r for r in caplog.records if "InChIKey index" in r.message]


@pytest.mark.unit
def test_concurrent_first_lookups_all_succeed(db):
    results = []
    threads = [
        threading.Thread(target=lambda: results.append(db.inchikey_to_dtxsid(ASPIRIN)))
        for _ in range(4)
    ]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert results == ["DTXSID5020108"] * 4
    assert _has_inchikey_index(db)


@pytest.mark.unit
@pytest.mark.skipif(os.name == "nt" or os.geteuid() == 0,
                    reason="needs POSIX permissions and a non-root user")
def test_read_only_database_scans_instead(db_path, caplog):
    os.chmod(db_path, stat.S_IRUSR)
    try:
        with CompToxID(db_path=db_path, auto_download=False) as client:
            with caplog.at_level(logging.WARNING, logger="provesid.comptox"):
                assert client.inchikey_to_dtxsid(ASPIRIN) == "DTXSID5020108"
                assert client.inchikey_to_dtxsid(DEHYDROACETIC_STANDARD) == "DTXSID6020014"
            assert not _has_inchikey_index(client)
        failures = [r for r in caplog.records if "could not be built" in r.message]
        assert len(failures) == 1
    finally:
        os.chmod(db_path, stat.S_IRUSR | stat.S_IWUSR)


@pytest.mark.unit
def test_download_builds_both_indexes(tmp_path, monkeypatch):
    source = tmp_path / "upstream.db"
    _write_database(source)

    def fake_download(url, destination, verify, **kwargs):
        verify(str(source))
        shutil.copy(source, destination)

    monkeypatch.setattr(provesid.comptox, "download_file", fake_download)
    with CompToxID(db_path=str(tmp_path / "comptox_chemicals.db")) as client:
        assert _has_inchikey_index(client)
        assert client.has_name_index
        assert client.conn.execute(
            f"SELECT COUNT(*) FROM {NAME_INDEX_TABLE}"
        ).fetchone()[0] > 0


# ── Using it ─────────────────────────────────────────────────────────────────


@pytest.mark.unit
def test_the_lookup_searches_the_index(db):
    db.get_by_inchikey(ASPIRIN)
    plan = _plan(db, "SELECT * FROM chemicals WHERE INCHIKEY IN (?, ?)", (ASPIRIN, ASPIRIN))
    assert f"USING INDEX {INCHIKEY_INDEX}" in plan


@pytest.mark.unit
def test_the_skeleton_search_searches_the_index(db):
    db.get_by_inchikey(ASPIRIN)
    plan = _plan(db, "SELECT * FROM chemicals WHERE INCHIKEY GLOB ? LIMIT 20", (f"{ASPIRIN[:14]}*",))
    assert f"USING INDEX {INCHIKEY_INDEX}" in plan
    rows = comptox_skeleton_search(db, ASPIRIN[:14])
    assert sorted(row["PREFERRED_NAME"] for row in rows) == ["Aspirin", "Aspirin anion"]


@pytest.mark.unit
def test_the_index_changes_no_answer(db_path):
    queries = [ASPIRIN, ASPIRIN_ANION, DEHYDROACETIC_STANDARD, DEHYDROACETIC_NONSTANDARD,
               "XXXXXXXXXXXXXX-UHFFFAOYSA-N"]
    with CompToxID(db_path=db_path, auto_download=False) as client:
        client._inchikey_index_checked = True     # scan, as before the index
        scanned = [client.inchikey_to_dtxsid(key) for key in queries]
        assert not _has_inchikey_index(client)
    with CompToxID(db_path=db_path, auto_download=False) as client:
        indexed = [client.inchikey_to_dtxsid(key) for key in queries]
        assert _has_inchikey_index(client)
    assert indexed == scanned == [
        "DTXSID5020108", "DTXSID0000003", "DTXSID6020014", "DTXSID6020014", None,
    ]
