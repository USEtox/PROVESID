"""
The InChIKey-skeleton searches in :mod:`provesid.sources`.

They read ``client._conn``, which the SQLite clients stopped having when
their connection handling moved to :class:`~provesid.sqlite_client.SQLiteClient`
(it is ``conn``). The AttributeError was caught and logged, so
``Search(inchikey_skeleton=True)`` silently got nothing from CompTox or
PubChem. These run against real, in-memory SQLite tables.
"""

import sqlite3

import pytest

from provesid.comptox import CompToxID
from provesid.pubchem_id import PubChemID
from provesid.sources import comptox_skeleton_search, pubchem_skeleton_search

ASPIRIN = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
ASPIRIN_ANION = "BSYNRYMUTXBXSQ-UHFFFAOYSA-M"
ETHANOL = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"


def _client(cls, table, key_column, rows):
    """An instance of ``cls`` over one in-memory table, without its __init__."""
    connection = sqlite3.connect(":memory:", check_same_thread=False)
    connection.row_factory = sqlite3.Row
    connection.execute(f"CREATE TABLE {table} (name TEXT, {key_column} TEXT)")
    connection.executemany(f"INSERT INTO {table} VALUES (?, ?)", rows)
    client = object.__new__(cls)
    client._adopt_connection(connection)
    return client


@pytest.mark.unit
def test_comptox_skeleton_search_finds_every_key_on_the_skeleton():
    comptox = _client(CompToxID, "chemicals", "INCHIKEY",
                      [("Aspirin", ASPIRIN), ("Aspirin anion", ASPIRIN_ANION),
                       ("Ethanol", ETHANOL)])
    rows = comptox_skeleton_search(comptox, ASPIRIN[:14])
    assert sorted(row["name"] for row in rows) == ["Aspirin", "Aspirin anion"]


@pytest.mark.unit
def test_pubchem_skeleton_search_finds_every_key_on_the_skeleton():
    pubchem = _client(PubChemID, "compounds", "inchikey",
                      [("Aspirin", ASPIRIN), ("Ethanol", ETHANOL)])
    rows = pubchem_skeleton_search(pubchem, ASPIRIN[:14])
    assert [row["name"] for row in rows] == ["Aspirin"]
