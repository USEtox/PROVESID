"""Tests for the CompTox lookup indexes, each built by the first lookup by its column.

The downloaded database indexes DTXSID, CASRN and the preferred name only, so
each lookup by InChIKey, SMILES, DTXCID or molecular formula scanned the
1.2 M rows (0.14 to 0.2 s) and each missed skeleton search took 0.9 s.  Each
test builds a small database in the real schema, so nothing is downloaded and
the real file is never touched.
"""

import logging
import os
import shutil
import sqlite3
import stat
import threading

import pytest

import provesid.comptox
from provesid.comptox import LOOKUP_INDEXES, NAME_INDEX_TABLE, CompToxID
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
ASPIRIN_SMILES = "CC(=O)OC1=C(C=CC=C1)C(O)=O"
CAFFEIC_SMILES = "OC(=O)C=CC1=CC(O)=C(O)C=C1"

# Two chemicals share a formula, and two share a SMILES (as 893 do in the
# 2025 release), so that the order of the answers is tested too.
_ROWS = [
    {"DTXSID": "DTXSID5020108", "PREFERRED_NAME": "Aspirin", "CASRN": "50-78-2",
     "DTXCID": "DTXCID50108", "INCHIKEY": ASPIRIN, "SMILES": ASPIRIN_SMILES,
     "MOLECULAR_FORMULA": "C9H8O4", "IDENTIFIER": "Aspirin"},
    {"DTXSID": "DTXSID0000003", "PREFERRED_NAME": "Aspirin anion", "CASRN": "0-00-3",
     "DTXCID": "DTXCID3", "INCHIKEY": ASPIRIN_ANION, "SMILES": "CC(=O)OC1=CC=CC=C1C([O-])=O",
     "MOLECULAR_FORMULA": "C9H7O4", "IDENTIFIER": "Aspirin anion"},
    {"DTXSID": "DTXSID6020014", "PREFERRED_NAME": "Dehydroacetic acid", "CASRN": "520-45-6",
     "DTXCID": "DTXCID6014", "INCHIKEY": DEHYDROACETIC_NONSTANDARD,
     "SMILES": "CC(=O)C1C(=O)OC(C)=CC1=O", "MOLECULAR_FORMULA": "C8H8O4",
     "IDENTIFIER": "Dehydroacetic acid"},
    {"DTXSID": "DTXSID0000004", "PREFERRED_NAME": "Mixture", "CASRN": "0-00-4",
     "INCHIKEY": None, "IDENTIFIER": "Mixture"},
    {"DTXSID": "DTXSID0000005", "PREFERRED_NAME": "Caffeic acid", "CASRN": "331-39-5",
     "DTXCID": "DTXCID5", "SMILES": CAFFEIC_SMILES, "MOLECULAR_FORMULA": "C9H8O4",
     "IDENTIFIER": "Caffeic acid"},
    {"DTXSID": "DTXSID0000006", "PREFERRED_NAME": "Caffeic acid, duplicate", "CASRN": "0-00-6",
     "DTXCID": "DTXCID6", "SMILES": CAFFEIC_SMILES, "MOLECULAR_FORMULA": "C9H8O4",
     "IDENTIFIER": "Caffeic acid, duplicate"},
]

# A lookup by each indexed column, and the column's hits and a miss.
LOOKUPS = {
    "INCHIKEY": (
        lambda db, key: db.inchikey_to_dtxsid(key),
        [ASPIRIN, ASPIRIN_ANION, DEHYDROACETIC_STANDARD, DEHYDROACETIC_NONSTANDARD,
         "XXXXXXXXXXXXXX-UHFFFAOYSA-N"],
        ["DTXSID5020108", "DTXSID0000003", "DTXSID6020014", "DTXSID6020014", None],
    ),
    "SMILES": (
        lambda db, smiles: db.smiles_to_dtxsid(smiles),
        [ASPIRIN_SMILES, CAFFEIC_SMILES, "C"],
        ["DTXSID5020108", "DTXSID0000005", None],
    ),
    "DTXCID": (
        lambda db, dtxcid: (db.get_by_dtxcid(dtxcid) or {}).get("DTXSID"),
        ["DTXCID50108", "DTXCID6", "DTXCID0"],
        ["DTXSID5020108", "DTXSID0000006", None],
    ),
    "MOLECULAR_FORMULA": (
        lambda db, formula: [r["DTXSID"] for r in db.search_by_formula(formula)],
        ["C9H8O4", "C8H8O4", "H2O"],
        [["DTXSID5020108", "DTXSID0000005", "DTXSID0000006"], ["DTXSID6020014"], []],
    ),
}


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


def _has_index(client, column):
    return client.conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='index' AND name=?", (LOOKUP_INDEXES[column],)
    ).fetchone() is not None


def _plan(client, sql, params):
    return " ".join(row[3] for row in client.conn.execute(f"EXPLAIN QUERY PLAN {sql}", params))


def _lookup(client, column, query):
    lookup, _, _ = LOOKUPS[column]
    return lookup(client, query)


@pytest.mark.unit
def test_every_indexed_column_has_a_lookup_here():
    assert set(LOOKUPS) == set(LOOKUP_INDEXES)


# ── Building ─────────────────────────────────────────────────────────────────


@pytest.mark.unit
def test_opening_the_database_builds_no_index(db):
    db.get_by_casrn("50-78-2")
    assert not any(_has_index(db, column) for column in LOOKUP_INDEXES)


@pytest.mark.unit
@pytest.mark.parametrize("column", LOOKUP_INDEXES)
def test_first_lookup_builds_its_own_index_only_and_says_so_once(db, column, caplog):
    _, queries, _ = LOOKUPS[column]
    with caplog.at_level(logging.WARNING, logger="provesid.comptox"):
        for query in queries:
            _lookup(db, column, query)
    assert _has_index(db, column)
    assert not any(_has_index(db, other) for other in LOOKUP_INDEXES if other != column)
    building = [r for r in caplog.records if "index" in r.message]
    assert len(building) == 1
    assert f"CompTox {column} index" in building[0].message


@pytest.mark.unit
@pytest.mark.parametrize("column", LOOKUP_INDEXES)
def test_a_database_that_has_the_index_is_left_alone(db_path, column, caplog):
    _, queries, _ = LOOKUPS[column]
    with CompToxID(db_path=db_path, auto_download=False) as client:
        _lookup(client, column, queries[0])
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="provesid.comptox"):
        with CompToxID(db_path=db_path, auto_download=False) as client:
            _lookup(client, column, queries[0])
    assert not [r for r in caplog.records if "index" in r.message]


@pytest.mark.unit
def test_concurrent_first_lookups_all_succeed(db):
    results = []
    lookups = [
        lambda: results.append(db.inchikey_to_dtxsid(ASPIRIN)),
        lambda: results.append(db.smiles_to_dtxsid(ASPIRIN_SMILES)),
        lambda: results.append(db.get_by_dtxcid("DTXCID50108")["DTXSID"]),
        lambda: results.append(db.search_by_formula("C9H8O4")[0]["DTXSID"]),
    ]
    threads = [threading.Thread(target=lookup) for lookup in lookups * 2]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert results == ["DTXSID5020108"] * 8
    assert all(_has_index(db, column) for column in LOOKUP_INDEXES)


@pytest.mark.unit
@pytest.mark.skipif(os.name == "nt" or os.geteuid() == 0,
                    reason="needs POSIX permissions and a non-root user")
def test_read_only_database_scans_instead(db_path, caplog):
    os.chmod(db_path, stat.S_IRUSR)
    try:
        with CompToxID(db_path=db_path, auto_download=False) as client:
            with caplog.at_level(logging.WARNING, logger="provesid.comptox"):
                for column, (lookup, queries, expected) in LOOKUPS.items():
                    assert [lookup(client, q) for q in queries] == expected
                    assert [lookup(client, q) for q in queries] == expected
            assert not any(_has_index(client, column) for column in LOOKUP_INDEXES)
        failures = [r for r in caplog.records if "could not be built" in r.message]
        assert sorted(r.args[0] for r in failures) == sorted(LOOKUP_INDEXES)
    finally:
        os.chmod(db_path, stat.S_IRUSR | stat.S_IWUSR)


@pytest.mark.unit
def test_download_builds_every_index(tmp_path, monkeypatch):
    source = tmp_path / "upstream.db"
    _write_database(source)

    def fake_download(url, destination, verify, **kwargs):
        verify(str(source))
        shutil.copy(source, destination)

    monkeypatch.setattr(provesid.comptox, "download_file", fake_download)
    with CompToxID(db_path=str(tmp_path / "comptox_chemicals.db")) as client:
        assert all(_has_index(client, column) for column in LOOKUP_INDEXES)
        assert client.has_name_index
        assert client.conn.execute(
            f"SELECT COUNT(*) FROM {NAME_INDEX_TABLE}"
        ).fetchone()[0] > 0


# ── Using them ───────────────────────────────────────────────────────────────


@pytest.mark.unit
@pytest.mark.parametrize("column", LOOKUP_INDEXES)
def test_the_lookup_searches_the_index(db, column):
    _, queries, _ = LOOKUPS[column]
    _lookup(db, column, queries[0])
    plan = _plan(db, f"SELECT * FROM chemicals WHERE {column} = ?", (queries[0],))
    assert f"USING INDEX {LOOKUP_INDEXES[column]}" in plan


@pytest.mark.unit
def test_the_inchikey_lookup_of_several_spellings_searches_the_index(db):
    db.get_by_inchikey(ASPIRIN)
    plan = _plan(db, "SELECT * FROM chemicals WHERE INCHIKEY IN (?, ?)", (ASPIRIN, ASPIRIN))
    assert f"USING INDEX {LOOKUP_INDEXES['INCHIKEY']}" in plan


@pytest.mark.unit
def test_the_skeleton_search_searches_the_index(db):
    db.get_by_inchikey(ASPIRIN)
    plan = _plan(db, "SELECT * FROM chemicals WHERE INCHIKEY GLOB ? LIMIT 20", (f"{ASPIRIN[:14]}*",))
    assert f"USING INDEX {LOOKUP_INDEXES['INCHIKEY']}" in plan
    rows = comptox_skeleton_search(db, ASPIRIN[:14])
    assert sorted(row["PREFERRED_NAME"] for row in rows) == ["Aspirin", "Aspirin anion"]


@pytest.mark.unit
@pytest.mark.parametrize("column", LOOKUP_INDEXES)
def test_the_index_changes_no_answer(db_path, column):
    lookup, queries, expected = LOOKUPS[column]
    with CompToxID(db_path=db_path, auto_download=False) as client:
        client._lookup_indexes_checked = frozenset(LOOKUP_INDEXES)   # scan, as before
        scanned = [lookup(client, q) for q in queries]
        assert not _has_index(client, column)
    with CompToxID(db_path=db_path, auto_download=False) as client:
        indexed = [lookup(client, q) for q in queries]
        assert _has_index(client, column)
    assert indexed == scanned == expected
