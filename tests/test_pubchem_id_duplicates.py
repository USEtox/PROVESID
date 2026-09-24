"""
Repeated rows in ``pubchem_id.db`` must not surface as repeated answers.

The Zenodo copy of the database lists some ``(cid, cas)`` pairs twice ---
aspirin's ``50-78-2`` among them --- and a compound can match a name search
through both its title and a synonym. Both used to come back twice. The
database here is a temporary SQLite file with the real schema.
"""

import sqlite3

import pytest

from provesid.pubchem_id import PubChemID


@pytest.fixture
def db(tmp_path):
    """Aspirin with its CAS number listed twice and a synonym matching its title."""
    path = tmp_path / "pubchem_id.db"
    conn = sqlite3.connect(path)
    conn.executescript("""
        CREATE TABLE compounds (
            cid INTEGER PRIMARY KEY, cmpdname TEXT, mf TEXT, inchi TEXT,
            smiles TEXT, inchikey TEXT, iupacname TEXT, mw REAL,
            exactmass REAL, cidcdate TEXT);
        CREATE TABLE cas_numbers (id INTEGER PRIMARY KEY, cid INTEGER, cas TEXT);
        CREATE TABLE synonyms (id INTEGER PRIMARY KEY, cid INTEGER, synonym TEXT);
        INSERT INTO compounds (cid, cmpdname, mf) VALUES (2244, 'Aspirin', 'C9H8O4');
        INSERT INTO cas_numbers (cid, cas) VALUES
            (2244, '50-78-2'), (2244, '50-78-2'), (2244, '11126-35-5');
        INSERT INTO synonyms (cid, synonym) VALUES
            (2244, 'aspirin'), (2244, 'Aspirin (USP)');
    """)
    conn.commit()
    conn.close()
    database = PubChemID(db_path=str(path), auto_download=False)
    yield database
    database.close()


@pytest.mark.unit
def test_a_repeated_cas_row_is_reported_once_in_database_order(db):
    assert db.cid_to_cas(2244) == ['50-78-2', '11126-35-5']
    assert db.name_to_cas('aspirin') == ['50-78-2', '11126-35-5']


@pytest.mark.unit
def test_a_compound_matching_title_and_synonym_is_returned_once(db):
    assert [r['cid'] for r in db.search_by_name('aspirin')] == [2244]
