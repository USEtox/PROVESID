"""
Tests for bulk property retrieval and the offline-first property lookup.

Everything here is offline. The PUG-REST transport (``_make_request``) is
stubbed and the local database is a temporary SQLite file carrying the same
schema as ``pubchem_id.db``, so the tests assert the routing decisions — path
versus POST, offline versus online, absent versus unknown — rather than
PubChem's current answers.
"""

import sqlite3

import pandas as pd
import pytest

from provesid.pubchem import (
    PROPERTY_CHUNK_SIZE,
    URL_IDENTIFIER_LIMIT,
    PubChemAPI,
    PubChemNotFoundError,
    PubChemServerError,
)
from provesid.pubchem_id import PubChemID


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

class _FakeResponse:
    """The slice of ``requests.Response`` the parser and error handler touch."""

    def __init__(self, payload):
        self.status_code = 200
        self._payload = payload
        self.text = str(payload)

    def json(self):
        return self._payload


def _property_table(rows):
    """Wrap property rows the way PUG-REST does."""
    return {"PropertyTable": {"Properties": rows}}


@pytest.fixture
def api():
    """A client with caching off, so stubbed answers are never served from disk."""
    return PubChemAPI(use_cache=False)


@pytest.fixture
def recorder(api, monkeypatch):
    """
    Replace the transport with a recorder that answers from a CID table.

    Returns the list of ``(url, method, data)`` calls made, so a test can assert
    both the answer and how many requests it took to get it.
    """
    calls = []
    known = {2244: "C9H8O4", 702: "C2H6O", 5793: "C6H12O6"}

    def transport(url, method="GET", data=None, timeout=30, headers=None):
        calls.append((url, method, data))
        if method == "POST":
            identifiers = data["cid"]
        else:
            # .../compound/cid/<identifiers>/property/<props>/JSON
            identifiers = url.split("/cid/")[1].split("/property/")[0]
        rows = []
        for raw in identifiers.split(","):
            cid = int(raw)
            row = {"CID": cid}
            if cid in known:
                row["MolecularFormula"] = known[cid]
            rows.append(row)
        return _FakeResponse(_property_table(rows))

    monkeypatch.setattr(api, "_make_request", transport)
    monkeypatch.setattr(api, "_rate_limit", lambda: None)
    return calls


@pytest.fixture
def local_db(tmp_path):
    """
    A miniature ``pubchem_id.db`` holding three compounds.

    CID 2244 is complete, CID 233 has a NULL ``xlogp`` — the way PubChem records
    a compound it computes no logP for — and CID 702 is complete. CID 5793 is
    deliberately absent, so it can stand for a compound only the online API
    knows.
    """
    path = tmp_path / "pubchem_id.db"
    conn = sqlite3.connect(path)
    conn.executescript("""
        CREATE TABLE compounds (
            cid INTEGER PRIMARY KEY, cmpdname TEXT, mf TEXT, inchi TEXT,
            smiles TEXT, inchikey TEXT, iupacname TEXT, mw REAL, polararea REAL,
            complexity REAL, xlogp REAL, heavycnt INTEGER, hbonddonor INTEGER,
            hbondacc INTEGER, rotbonds INTEGER, exactmass REAL, charge INTEGER,
            cidcdate TEXT);
        CREATE TABLE cas_numbers (id INTEGER PRIMARY KEY, cid INTEGER, cas TEXT);
        CREATE TABLE synonyms (id INTEGER PRIMARY KEY, cid INTEGER, synonym TEXT);
    """)
    conn.executemany(
        "INSERT INTO compounds (cid, cmpdname, mf, smiles, inchikey, mw, xlogp, "
        "charge, heavycnt) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        [
            (2244, "Aspirin", "C9H8O4", "CC(=O)OC1=CC=CC=C1C(=O)O",
             "BSYNRYMUTXBXSQ-UHFFFAOYSA-N", 180.16, 1.2, 0, 13),
            (233, "Arsenate", "AsHO4-2", "[O-][As](=O)([O-])O", None, 139.93,
             None, -2, 5),
            (702, "Ethanol", "C2H6O", "CCO", "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
             46.07, -0.1, 0, 3),
        ],
    )
    conn.commit()
    conn.close()
    return str(path)


@pytest.fixture
def db(local_db, api, recorder):
    """The offline database wired to the recording online client."""
    return PubChemID(db_path=local_db, auto_download=False, api=api)


# --------------------------------------------------------------------------
# Bulk retrieval: one request for many CIDs
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_many_cids_cost_one_request(api, recorder):
    """The point of the bulk path: N compounds, one round trip."""
    rows = api.get_properties_for_cids([2244, 702, 5793], ["MolecularFormula"])

    assert len(recorder) == 1
    assert {row["CID"]: row["MolecularFormula"] for row in rows} == {
        2244: "C9H8O4", 702: "C2H6O", 5793: "C6H12O6"}


@pytest.mark.unit
def test_short_identifier_list_goes_in_the_url(api, recorder):
    """A list that fits in a URL is fetched with GET, as PubChem prefers."""
    api.get_properties_for_cids([2244, 702], ["MolecularFormula"])

    url, method, data = recorder[0]
    assert method == "GET"
    assert data is None
    assert "/cid/2244,702/property/MolecularFormula/JSON" in url


@pytest.mark.unit
def test_long_identifier_list_switches_to_post(api, recorder):
    """
    Past the URL length ceiling the identifiers move into the request body.

    The URL itself must then carry no identifier segment, or PubChem would see
    the operation where it expects the CIDs.
    """
    cids = list(range(10_000, 10_000 + 400))
    api.get_properties_for_cids(cids, ["MolecularFormula"], chunk_size=1000)

    url, method, data = recorder[0]
    assert method == "POST"
    assert data["cid"].startswith("10000,10001")
    assert url.endswith("/compound/cid/property/MolecularFormula/JSON")
    assert len(",".join(str(c) for c in cids)) > URL_IDENTIFIER_LIMIT


@pytest.mark.unit
def test_chunking_splits_the_request(api, recorder):
    """Chunk size bounds the work lost when one request has to be retried."""
    api.get_properties_for_cids(list(range(500)), ["MolecularFormula"], chunk_size=200)

    assert len(recorder) == 3


@pytest.mark.unit
def test_duplicate_cids_are_asked_once(api, recorder):
    """A repeated CID is one row, and the first-appearance order is kept."""
    rows = api.get_properties_for_cids([702, 2244, 702], ["MolecularFormula"])

    assert [row["CID"] for row in rows] == [702, 2244]
    assert recorder[0][0].count("702") == 1


@pytest.mark.unit
def test_unknown_cid_comes_back_as_a_bare_row(api, recorder):
    """PubChem answers a CID it has no record of with the CID and nothing else."""
    rows = api.get_properties_for_cids([2244, 424242], ["MolecularFormula"])

    by_cid = {row["CID"]: row for row in rows}
    assert by_cid[424242] == {"CID": 424242}


@pytest.mark.unit
def test_empty_inputs(api, recorder):
    """No CIDs is an empty answer; no properties is a programming error."""
    assert api.get_properties_for_cids([], ["MolecularFormula"]) == []
    assert recorder == []

    with pytest.raises(ValueError):
        api.get_properties_for_cids([2244], [])
    with pytest.raises(ValueError):
        api.get_properties_for_cids([2244], ["MolecularFormula"], chunk_size=0)


@pytest.mark.unit
def test_absent_record_is_an_empty_table(api, monkeypatch):
    """A 404 for the whole request is absence, and absence is an empty list."""
    def absent(url, method="GET", data=None, timeout=30, headers=None):
        raise PubChemNotFoundError("Resource not found (PUGREST.NotFound)")

    monkeypatch.setattr(api, "_make_request", absent)
    assert api.get_properties_for_cids([424242], ["MolecularFormula"]) == []


@pytest.mark.unit
def test_transport_failure_is_not_an_empty_table(api, monkeypatch):
    """
    A busy server must not be reported as "these compounds have no properties".

    Partial results are not returned either: a caller seeing a short table would
    have no way to tell it apart from a complete one.
    """
    def busy(url, method="GET", data=None, timeout=30, headers=None):
        raise PubChemServerError("Server busy (PUGREST.ServerBusy)")

    monkeypatch.setattr(api, "_make_request", busy)
    with pytest.raises(PubChemServerError):
        api.get_properties_for_cids([2244], ["MolecularFormula"])


@pytest.mark.unit
def test_batch_reports_a_row_per_requested_cid(api, recorder):
    """
    The legacy-shaped wrapper keeps every CID, answered or not.

    A caller that zips the result against its input must not have rows silently
    dropped, so an unknown CID is reported rather than omitted.
    """
    rows = api.get_compound_properties_batch([2244, 424242, 702], ["MolecularFormula"])

    assert [row["cid"] for row in rows] == [2244, 424242, 702]
    assert rows[0]["success"] is True
    assert rows[0]["MolecularFormula"] == "C9H8O4"
    assert rows[1]["success"] is False
    assert "No such compound" in rows[1]["error"]
    assert len(recorder) == 1


# --------------------------------------------------------------------------
# Offline first
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_local_database_answers_without_network(db, recorder):
    """The whole point of the offline database: no request at all."""
    result = db.properties(2244, ["MolecularFormula", "MolecularWeight"])

    assert result == {"CID": 2244, "Source": "offline",
                      "MolecularFormula": "C9H8O4", "MolecularWeight": 180.16}
    assert recorder == []


@pytest.mark.unit
def test_string_cid_is_accepted(db, recorder):
    """``"2244"`` and ``2244`` name the same compound and share a row."""
    assert db.properties("2244", ["MolecularFormula"])["CID"] == 2244
    assert recorder == []


@pytest.mark.unit
def test_null_column_is_reported_as_absent_not_as_none(db, recorder):
    """
    PubChem omits a property it has no value for, and so does this.

    CID 233 has no InChIKey in this database, so the key is missing rather than
    present-and-None — and missing must not trigger a pointless online request.
    """
    result = db.properties(233, ["MolecularFormula", "InChIKey"])

    assert result["MolecularFormula"] == "AsHO4-2"
    assert "InChIKey" not in result
    assert recorder == []


@pytest.mark.unit
def test_descriptors_go_online_even_when_the_database_has_the_column(db, recorder):
    """
    XLogP is PubChem's model output, not data about the compound.

    This database carries an ``xlogp`` column, as every Zenodo copy does, and
    the lookup still asks PubChem: a descriptor is served from one place, so
    the answer never depends on which CIDs happen to be local.
    """
    def transport(url, method="GET", data=None, timeout=30, headers=None):
        recorder.append((url, method, data))
        return _FakeResponse(_property_table([{"CID": 2244, "XLogP": 1.2}]))

    db.api._make_request = transport
    assert "XLogP" not in db.offline_properties
    assert db.properties(2244, ["XLogP"]) == {"CID": 2244, "Source": "online", "XLogP": 1.2}
    assert len(recorder) == 1


@pytest.mark.unit
def test_cid_absent_offline_falls_back_online(db, recorder):
    """CID 5793 is not in the local database, so it is fetched."""
    result = db.properties(5793, ["MolecularFormula"])

    assert result == {"CID": 5793, "Source": "online",
                      "MolecularFormula": "C6H12O6"}
    assert len(recorder) == 1


@pytest.mark.unit
def test_property_absent_offline_falls_back_online(db, recorder):
    """
    A property with no local column sends the request online.

    This database predates the ``monoisotopicmass`` column, so even a CID held
    locally has to be asked about online.
    """
    monkeypatched = {"CID": 2244, "MonoisotopicMass": "180.04225873"}

    def transport(url, method="GET", data=None, timeout=30, headers=None):
        recorder.append((url, method, data))
        return _FakeResponse(_property_table([monkeypatched]))

    db.api._make_request = transport
    result = db.properties(2244, ["MonoisotopicMass"])

    assert result == {"CID": 2244, "Source": "online",
                      "MonoisotopicMass": 180.04225873}
    assert len(recorder) == 1


@pytest.mark.unit
def test_offline_only_never_reaches_the_network(db, recorder):
    """``use_online_fallback=False`` is a hard guarantee, not a preference."""
    assert db.properties(5793, ["MolecularFormula"], use_online_fallback=False) is None
    assert db.properties(2244, ["MonoisotopicMass"], use_online_fallback=False) is None
    assert recorder == []


@pytest.mark.unit
def test_unknown_to_both_sources_is_none(db, recorder):
    """A CID neither source knows is reported as unknown, not as empty."""
    assert db.properties(424242, ["MolecularFormula"]) is None
    assert len(recorder) == 1


@pytest.mark.unit
def test_bulk_lookup_asks_online_only_for_what_is_missing(db, recorder):
    """
    The local database absorbs most of the work; one request covers the rest.

    Two of these three CIDs are local, so exactly one CID should appear in the
    single request that is made.
    """
    rows = db.properties_for_cids([2244, 5793, 702], ["MolecularFormula"])

    assert [row["CID"] for row in rows] == [2244, 5793, 702]
    assert [row["Source"] for row in rows] == ["offline", "online", "offline"]
    assert len(recorder) == 1
    assert recorder[0][0].count("5793") == 1
    assert "2244" not in recorder[0][0]


@pytest.mark.unit
def test_defaults_to_every_locally_available_property(db, recorder):
    """Naming no properties must not turn a local lookup into a request."""
    result = db.properties(2244)

    assert recorder == []
    assert result["Title"] == "Aspirin"
    assert result["MolecularWeight"] == 180.16
    # This database predates ``monoisotopicmass``, so the default leaves it out
    # rather than sending every default lookup online for it.
    assert "MonoisotopicMass" not in db.offline_properties
    assert set(result) <= set(db.offline_properties) | {"CID", "Source"}


@pytest.mark.unit
def test_values_have_the_same_type_from_either_source(db, recorder):
    """
    A table assembled from both sources has to be usable as one table.

    PUG-REST reports ``MolecularWeight`` as the string ``"46.07"`` while SQLite
    holds a float, so the wrapper normalises both to float.
    """
    def transport(url, method="GET", data=None, timeout=30, headers=None):
        recorder.append((url, method, data))
        return _FakeResponse(_property_table(
            [{"CID": 5793, "MolecularWeight": "180.16", "HeavyAtomCount": "12"}]))

    db.api._make_request = transport
    offline, online = db.properties_for_cids([702, 5793], ["MolecularWeight"])
    assert isinstance(offline["MolecularWeight"], float)
    assert isinstance(online["MolecularWeight"], float)

    counted = db.properties(5793, ["HeavyAtomCount"])
    assert isinstance(counted["HeavyAtomCount"], int)


@pytest.mark.unit
def test_table_has_a_row_for_every_cid_asked_about(db, recorder):
    """Including the ones nothing knows, so the frame can be joined."""
    table = db.properties_table([2244, 5793, 424242], ["MolecularFormula"])

    assert isinstance(table, pd.DataFrame)
    assert list(table.columns) == ["CID", "Source", "MolecularFormula"]
    assert list(table["CID"]) == [2244, 5793, 424242]
    assert list(table["Source"]) == ["offline", "online", "missing"]
    assert pd.isna(table.loc[2, "MolecularFormula"])


@pytest.mark.unit
def test_table_is_strictly_offline_on_request(db, recorder):
    """Every non-local CID reads as missing rather than being fetched."""
    table = db.properties_table([2244, 5793], ["MolecularFormula"],
                                use_online_fallback=False)

    assert list(table["Source"]) == ["offline", "missing"]
    assert recorder == []


@pytest.mark.unit
def test_malformed_cid_is_rejected_before_any_request(db, recorder):
    """
    PubChem answers one bad CID by failing the whole batch, so catch it here.

    A ``PUGREST.BadRequest`` names nothing; a ValueError names the value.
    """
    with pytest.raises(ValueError, match="abc"):
        db.properties("abc")
    with pytest.raises(ValueError):
        db.properties_for_cids([2244, None], ["MolecularFormula"])
    assert recorder == []


@pytest.mark.unit
def test_uncastable_value_is_passed_through(db):
    """A surprising value reaches the caller rather than being dropped."""
    assert PubChemID._cast_property("MolecularWeight", "n/a") == "n/a"
    assert PubChemID._cast_property("MolecularFormula", "C9H8O4") == "C9H8O4"


@pytest.mark.unit
def test_large_offline_lookup_is_one_statement_per_batch(db, recorder):
    """
    More CIDs than SQLite will bind at once still resolve in one call.

    The IN list is filled in batches of ``_SQL_PARAMETER_LIMIT``; asking about
    several times that many CIDs must return the same rows as asking about a
    handful, in the same order.
    """
    assert len(range(1, 2500)) > PubChemID._SQL_PARAMETER_LIMIT * 2
    cids = list(range(1, 2500))
    rows = db.properties_for_cids(cids, ["MolecularFormula"],
                                  use_online_fallback=False)

    assert [row["CID"] for row in rows] == [233, 702, 2244]
    assert recorder == []


@pytest.mark.unit
def test_offline_property_map_names_real_properties():
    """
    Every locally served name must be a property PubChem actually publishes.

    A typo here would send the fallback to PUG-REST with an invalid property
    name, which fails the whole request.
    """
    from provesid.pubchem import CompoundProperties

    published = {value for name, value in vars(CompoundProperties).items()
                 if not name.startswith('_')}
    assert set(PubChemID.OFFLINE_PROPERTIES) <= published
    assert PubChemID.DEFAULT_PROPERTIES == tuple(PubChemID.OFFLINE_PROPERTIES)


@pytest.mark.unit
def test_chunk_size_default_is_shared():
    """The database layer batches its fallback the same way the client does."""
    assert PROPERTY_CHUNK_SIZE > 0
