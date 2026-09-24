"""Tests for the CompTox name index: exact name lookups that find synonyms.

Each test builds a three-chemical database shaped like the real
``comptox_chemicals.db``, so nothing is downloaded and the real file is never
touched.
"""

import logging
import os
import sqlite3
import stat
import threading

import pytest

from provesid.comptox import NAME_INDEX_TABLE, NAME_KINDS, CompToxID, name_key
from provesid.sources import LOOKUPS, Query

_COLUMNS = (
    "DTXSID", "PREFERRED_NAME", "CASRN", "DTXCID", "INCHIKEY", "IUPAC_NAME", "SMILES",
    "MOLECULAR_FORMULA", "AVERAGE_MASS", "MONOISOTOPIC_MASS", "QSAR_READY_SMILES",
    "MS_READY_SMILES", "IDENTIFIER",
)

_ROWS = [
    {
        "DTXSID": "DTXSID5020108", "PREFERRED_NAME": "Aspirin", "CASRN": "50-78-2",
        "INCHIKEY": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N", "IUPAC_NAME": "2-acetyloxybenzoic acid",
        "SMILES": "CC(=O)OC1=CC=CC=C1C(O)=O", "MOLECULAR_FORMULA": "C9H8O4",
        "IDENTIFIER": "50-78-2 | 11126-35-5 | Aspirin | Acetylsalicylic acid | ASA | aspirin",
    },
    {
        "DTXSID": "DTXSID9020112", "PREFERRED_NAME": "Atrazine", "CASRN": "1912-24-9",
        "INCHIKEY": "MXWJVTOOROXGIU-UHFFFAOYSA-N", "IUPAC_NAME": None,
        "SMILES": "CCNC1=NC(Cl)=NC(NC(C)C)=N1", "MOLECULAR_FORMULA": "C8H14ClN5",
        "IDENTIFIER": "1912-24-9 | 39400-72-1 | Atrazin | Gesaprim",
    },
    {
        # Lists "Aspirin" as a synonym: must rank after the chemical called it.
        "DTXSID": "DTXSID0000001", "PREFERRED_NAME": "Aspirin mixture", "CASRN": "0-00-1",
        "INCHIKEY": None, "IUPAC_NAME": None, "SMILES": None, "MOLECULAR_FORMULA": None,
        "IDENTIFIER": "Aspirin | Gesaprim | 99999-99-7",
    },
    {
        # Lists the same stray number as the mixture: ambiguous, so unanswered.
        "DTXSID": "DTXSID0000002", "PREFERRED_NAME": "Unrelated", "CASRN": "0-00-2",
        "INCHIKEY": None, "IUPAC_NAME": None, "SMILES": None, "MOLECULAR_FORMULA": None,
        "IDENTIFIER": "99999-99-7 | 0-00-1",
    },
]


@pytest.fixture
def db_path(tmp_path):
    """A small database in the real schema, with no name index."""
    path = tmp_path / "comptox_chemicals.db"
    connection = sqlite3.connect(path)
    connection.execute(f"CREATE TABLE chemicals ({', '.join(c + ' TEXT' for c in _COLUMNS)})")
    connection.execute("CREATE INDEX idx_preferred_name ON chemicals(PREFERRED_NAME)")
    connection.executemany(
        f"INSERT INTO chemicals VALUES ({', '.join('?' for _ in _COLUMNS)})",
        [tuple(row.get(c) for c in _COLUMNS) for row in _ROWS],
    )
    connection.commit()
    connection.close()
    return str(path)


@pytest.fixture
def db(db_path):
    with CompToxID(db_path=db_path, auto_download=False) as client:
        yield client


def _names(rows):
    return [row["PREFERRED_NAME"] for row in rows]


# ── The key ──────────────────────────────────────────────────────────────────


def test_name_key_strips_and_lowercases():
    assert name_key("  Acetylsalicylic ACID ") == "acetylsalicylic acid"


def test_name_key_folds_beyond_ascii():
    # SQLite's lower() would leave these alone; the index is built in Python.
    assert name_key("ÄTHANOL") == "äthanol"


# ── Building ─────────────────────────────────────────────────────────────────


def test_a_downloaded_database_has_no_index_until_asked(db):
    assert not db.has_name_index


def test_build_name_index_counts_distinct_names_per_chemical(db):
    # Aspirin: aspirin, 2-acetyloxybenzoic acid, 50-78-2, 11126-35-5,
    # acetylsalicylic acid, asa (the second "aspirin" is a duplicate) = 6.
    # Atrazine: atrazine, 1912-24-9, 39400-72-1, atrazin, gesaprim = 5.
    # Mixture: aspirin mixture, aspirin, gesaprim, 99999-99-7 = 4.
    # Unrelated: unrelated, 99999-99-7, 0-00-1 = 3.
    assert db.build_name_index() == 18
    assert db.has_name_index


def test_build_name_index_twice_is_a_no_op(db):
    assert db.build_name_index() == db.build_name_index() == 18


def test_each_name_keeps_the_first_kind_it_appeared_under(db):
    db.build_name_index()
    kinds = dict(db.conn.execute(
        f"SELECT name_key, kind FROM {NAME_INDEX_TABLE} WHERE chemical_rowid = 1"
    ).fetchall())
    assert kinds["aspirin"] == NAME_KINDS["preferred"]
    assert kinds["2-acetyloxybenzoic acid"] == NAME_KINDS["iupac"]
    assert kinds["asa"] == NAME_KINDS["identifier"]


def test_first_exact_lookup_builds_the_index_and_says_so(db, caplog):
    with caplog.at_level(logging.WARNING, logger="provesid.comptox"):
        rows = db.search_by_name("Acetylsalicylic acid", exact=True)
    assert _names(rows) == ["Aspirin"]
    assert db.has_name_index
    assert any("Building the CompTox name index" in r.message for r in caplog.records)


def test_substring_lookup_does_not_build_the_index(db):
    db.search_by_name("spir", exact=False)
    assert not db.has_name_index


def test_concurrent_first_lookups_build_once(db):
    results = []
    threads = [
        threading.Thread(target=lambda: results.append(db.search_by_name("ASA", exact=True)))
        for _ in range(4)
    ]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert [_names(r) for r in results] == [["Aspirin"]] * 4
    assert db.conn.execute(f"SELECT COUNT(*) FROM {NAME_INDEX_TABLE}").fetchone()[0] == 18


@pytest.mark.skipif(os.name == "nt" or os.geteuid() == 0,
                    reason="needs POSIX permissions and a non-root user")
def test_read_only_database_falls_back_to_preferred_names(db_path, caplog):
    os.chmod(db_path, stat.S_IRUSR)
    try:
        with CompToxID(db_path=db_path, auto_download=False) as client:
            with caplog.at_level(logging.WARNING, logger="provesid.comptox"):
                assert _names(client.search_by_name("Aspirin", exact=True)) == ["Aspirin"]
                assert client.search_by_name("ASA", exact=True) == []
                client.search_by_name("ASA", exact=True)
            assert not client.has_name_index
        failures = [r for r in caplog.records if "could not be built" in r.message]
        assert len(failures) == 1
    finally:
        os.chmod(db_path, stat.S_IRUSR | stat.S_IWUSR)


# ── Looking up ───────────────────────────────────────────────────────────────


@pytest.mark.parametrize("query, expected", [
    ("Acetylsalicylic acid", ["Aspirin"]),       # a synonym
    ("ACETYLSALICYLIC ACID", ["Aspirin"]),       # in any case
    ("  asa ", ["Aspirin"]),                     # padded
    ("2-acetyloxybenzoic acid", ["Aspirin"]),    # the IUPAC name
    ("39400-72-1", ["Atrazine"]),                # a former CAS number
    ("Atrazin", ["Atrazine"]),
    ("Acetylsalicylic", []),                     # exact means exact
])
def test_exact_lookup_finds_every_name(db, query, expected):
    assert _names(db.search_by_name(query, exact=True)) == expected


def test_a_chemical_called_the_query_outranks_one_listing_it(db):
    assert _names(db.search_by_name("aspirin", exact=True)) == ["Aspirin", "Aspirin mixture"]


def test_ties_keep_database_order_and_respect_limit(db):
    assert _names(db.search_by_name("Gesaprim", exact=True)) == ["Atrazine", "Aspirin mixture"]
    assert _names(db.search_by_name("Gesaprim", exact=True, limit=1)) == ["Atrazine"]


def test_exact_rows_carry_parsed_identifiers(db):
    (row,) = db.search_by_name("ASA", exact=True)
    assert "Acetylsalicylic acid" in row["identifiers"]


def test_get_by_name_still_means_the_preferred_name(db):
    assert db.get_by_name("Aspirin")["DTXSID"] == "DTXSID5020108"
    assert db.get_by_name("ASA") is None


def test_search_name_cell_reaches_synonyms(db):
    candidates = LOOKUPS["name"]["comptox"](db, Query("Acetylsalicylic acid", k=5, label="name"))
    assert [c["DTXSID"] for c in candidates] == ["DTXSID5020108"]


# ── Retired and alternate CAS numbers ────────────────────────────────────────


def test_a_retired_cas_number_finds_the_chemical_that_lists_it(db):
    assert db.get_by_casrn("39400-72-1") is None
    assert db.get_by_alternate_casrn("39400-72-1")["DTXSID"] == "DTXSID9020112"


def test_a_number_listed_by_two_chemicals_is_not_guessed(db):
    assert db.get_by_alternate_casrn("99999-99-7") is None


def test_an_unknown_number_finds_nothing(db):
    assert db.get_by_alternate_casrn("12345-67-8") is None


def test_only_identifier_tokens_count_as_alternate_numbers(db):
    # The mixture lists "Aspirin" as an identifier token, and it is still no CAS number.
    assert db.get_by_alternate_casrn("Aspirin") is None


def test_the_cas_cell_prefers_a_chemical_s_own_casrn(db):
    # 0-00-1 is the mixture's own CASRN and also listed by "Unrelated".
    (cand,) = LOOKUPS["cas"]["comptox"](db, Query("0-00-1", k=1, label="cas"))
    assert cand["DTXSID"] == "DTXSID0000001"


def test_the_cas_cell_falls_back_to_a_retired_number(db):
    (cand,) = LOOKUPS["cas"]["comptox"](db, Query("39400-72-1", k=1, label="cas"))
    assert cand["DTXSID"] == "DTXSID9020112"


@pytest.mark.skipif(os.name == "nt" or os.geteuid() == 0,
                    reason="needs POSIX permissions and a non-root user")
def test_no_index_means_no_alternate_lookup(db_path):
    os.chmod(db_path, stat.S_IRUSR)
    try:
        with CompToxID(db_path=db_path, auto_download=False) as client:
            assert client.get_by_alternate_casrn("39400-72-1") is None
    finally:
        os.chmod(db_path, stat.S_IRUSR | stat.S_IWUSR)
