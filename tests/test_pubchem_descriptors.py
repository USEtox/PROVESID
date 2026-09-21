"""
Tests for computed descriptors on demand: ``PubChemID.descriptors`` and
:func:`provesid.pubchem_id.rdkit_descriptors`.

Everything here is offline. RDKit runs for real; the PUG-REST transport is
stubbed with one that answers whatever properties a request names from a small
table, and records every request, so the tests can assert both the values and
which source produced them --- in particular that ``source="rdkit"`` never
touches the network for a compound the local database holds, and never reads a
stored descriptor column even when a Zenodo copy has one.
"""

import sqlite3

import pytest

from provesid.pubchem import PubChemAPI
from provesid.pubchem_id import (
    PUBCHEM_DESCRIPTORS,
    RDKIT_DESCRIPTORS,
    PubChemID,
    _rdkit_descriptor_functions,
    rdkit_descriptors,
)

ASPIRIN = "CC(=O)OC1=CC=CC=C1C(=O)O"

#: RDKit's values for aspirin, every descriptor.
ASPIRIN_RDKIT = {
    'MolLogP': 1.3101,
    'TPSA': 63.6,
    'HBondDonorCount': 1,
    'HBondAcceptorCount': 3,
    'RotatableBondCount': 2,
    'HeavyAtomCount': 13,
    'Charge': 0,
}

#: What the stub PUG-REST knows. CID 5793 (glucose) is absent from the local
#: database, so it stands for a compound only PubChem has.
ONLINE = {
    2244: {'SMILES': ASPIRIN, 'XLogP': 1.2, 'TPSA': 63.6, 'Complexity': 212,
           'Charge': 0, 'HBondDonorCount': 1, 'HBondAcceptorCount': 4,
           'RotatableBondCount': 3, 'HeavyAtomCount': 13},
    5793: {'SMILES': "C(C1C(C(C(C(O1)O)O)O)O)O", 'XLogP': -2.6, 'TPSA': 110,
           'Complexity': 151, 'Charge': 0, 'HBondDonorCount': 5,
           'HBondAcceptorCount': 6, 'RotatableBondCount': 1,
           'HeavyAtomCount': 12},
}


class _FakeResponse:
    """The slice of ``requests.Response`` the parser touches."""

    def __init__(self, payload):
        self.status_code = 200
        self._payload = payload
        self.text = str(payload)

    def json(self):
        return self._payload


@pytest.fixture
def api():
    """A client with caching off, so stubbed answers are never served from disk."""
    return PubChemAPI(use_cache=False)


@pytest.fixture
def requests_made(api, monkeypatch):
    """
    Answer property requests from :data:`ONLINE`, recording each one.

    Returns the list of ``(cids, properties)`` requested, one entry per
    request.
    """
    calls = []

    def transport(url, method="GET", data=None, timeout=30, headers=None):
        properties = url.split("/property/")[1].split("/")[0].split(",")
        identifiers = (data["cid"] if method == "POST"
                       else url.split("/cid/")[1].split("/property/")[0])
        cids = [int(raw) for raw in identifiers.split(",")]
        calls.append((cids, properties))
        rows = []
        for cid in cids:
            known = ONLINE.get(cid, {})
            rows.append({"CID": cid, **{name: known[name]
                                        for name in properties if name in known}})
        return _FakeResponse({"PropertyTable": {"Properties": rows}})

    monkeypatch.setattr(api, "_make_request", transport)
    monkeypatch.setattr(api, "_rate_limit", lambda: None)
    return calls


def _make_database(path, *, zenodo_columns=False):
    """
    Write a miniature ``pubchem_id.db``.

    The default is the shape an FTP build has: no descriptor columns. With
    ``zenodo_columns=True`` it has the eight descriptor columns of a Zenodo
    copy, filled with PubChem's values for aspirin, so a test can check they
    are not what ``source="rdkit"`` returns.

    CID 2244 is aspirin, 702 ethanol, 999 a compound whose SMILES RDKit
    cannot parse, and 998 one with no structure at all.
    """
    descriptor_columns = (", polararea REAL, complexity REAL, xlogp REAL, "
                          "heavycnt INTEGER, hbonddonor INTEGER, hbondacc INTEGER, "
                          "rotbonds INTEGER, charge INTEGER") if zenodo_columns else ""
    conn = sqlite3.connect(path)
    conn.executescript(f"""
        CREATE TABLE compounds (
            cid INTEGER PRIMARY KEY, cmpdname TEXT, mf TEXT, inchi TEXT,
            smiles TEXT, inchikey TEXT, iupacname TEXT, mw REAL,
            exactmass REAL, monoisotopicmass REAL, cidcdate TEXT
            {descriptor_columns});
        CREATE TABLE cas_numbers (id INTEGER PRIMARY KEY, cid INTEGER, cas TEXT);
        CREATE TABLE synonyms (id INTEGER PRIMARY KEY, cid INTEGER, synonym TEXT);
    """)
    conn.executemany(
        "INSERT INTO compounds (cid, cmpdname, mf, smiles) VALUES (?, ?, ?, ?)",
        [
            (2244, "Aspirin", "C9H8O4", ASPIRIN),
            (702, "Ethanol", "C2H6O", "CCO"),
            (999, "Unreadable", "C", "C1CC"),
            (998, "Structureless", None, None),
        ],
    )
    if zenodo_columns:
        conn.execute("UPDATE compounds SET xlogp = 1.2, polararea = 63.6, "
                     "hbondacc = 4, rotbonds = 3 WHERE cid = 2244")
    conn.commit()
    conn.close()
    return str(path)


@pytest.fixture
def db(tmp_path, api, requests_made):
    """An FTP-shaped local database wired to the recording online client."""
    path = _make_database(tmp_path / "pubchem_id.db")
    with PubChemID(db_path=path, auto_download=False, api=api) as client:
        yield client


# --------------------------------------------------------------------------
# rdkit_descriptors: one structure, no database
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_every_rdkit_descriptor_of_aspirin():
    assert rdkit_descriptors(ASPIRIN) == ASPIRIN_RDKIT


@pytest.mark.unit
def test_descriptors_come_back_in_the_order_asked():
    result = rdkit_descriptors(ASPIRIN, ['HeavyAtomCount', 'TPSA'])
    assert list(result) == ['HeavyAtomCount', 'TPSA']


@pytest.mark.unit
def test_floats_carry_no_floating_point_residue():
    """CalcTPSA sums to 63.60000000000001 for aspirin; the tabulated value is 63.6."""
    assert rdkit_descriptors(ASPIRIN, ['TPSA']) == {'TPSA': 63.6}


@pytest.mark.unit
def test_counts_are_ints():
    result = rdkit_descriptors(ASPIRIN)
    for name in ('HBondDonorCount', 'HBondAcceptorCount', 'RotatableBondCount',
                 'HeavyAtomCount', 'Charge'):
        assert type(result[name]) is int, name


@pytest.mark.unit
def test_tpsa_counts_sulfur_and_phosphorus():
    """PubChem's TPSA includes S and P; Ertl's original, RDKit's default, does not."""
    methanethiol = rdkit_descriptors("CS", ['TPSA'])
    assert methanethiol['TPSA'] > 0


@pytest.mark.unit
def test_charge_is_the_net_formal_charge():
    assert rdkit_descriptors("[O-][As](=O)([O-])O", ['Charge']) == {'Charge': -2}
    assert rdkit_descriptors("C[N+](C)(C)C.[Cl-]", ['Charge']) == {'Charge': 0}


@pytest.mark.unit
@pytest.mark.parametrize("smiles", ["C1CC", "not a molecule", "", None])
def test_unreadable_structure_is_none(smiles):
    assert rdkit_descriptors(smiles) is None


@pytest.mark.unit
def test_rdkit_is_asked_for_xlogp_and_says_mollogp():
    with pytest.raises(ValueError, match="MolLogP"):
        rdkit_descriptors(ASPIRIN, ['XLogP'])


@pytest.mark.unit
def test_rdkit_is_asked_for_complexity_and_says_pubchem():
    with pytest.raises(ValueError, match="source='pubchem'"):
        rdkit_descriptors(ASPIRIN, ['Complexity'])


@pytest.mark.unit
def test_unknown_descriptor_lists_the_choices():
    with pytest.raises(ValueError, match="HeavyAtomCount"):
        rdkit_descriptors(ASPIRIN, ['Polarizability'])


@pytest.mark.unit
def test_empty_descriptor_list_is_an_error():
    with pytest.raises(ValueError, match="at least one"):
        rdkit_descriptors(ASPIRIN, [])


@pytest.mark.unit
def test_every_listed_rdkit_descriptor_has_a_function():
    assert tuple(_rdkit_descriptor_functions()) == RDKIT_DESCRIPTORS


@pytest.mark.unit
def test_every_pubchem_descriptor_is_type_normalised():
    """PUG-REST's values are cast like the rest of the property layer's."""
    assert set(PUBCHEM_DESCRIPTORS) <= set(PubChemID._PROPERTY_CASTS)


# --------------------------------------------------------------------------
# PubChemID.descriptors, source="rdkit"
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_rdkit_answers_from_the_stored_smiles_without_the_network(db, requests_made):
    assert db.descriptors(2244) == {'CID': 2244, 'Source': 'rdkit', **ASPIRIN_RDKIT}
    assert requests_made == []


@pytest.mark.unit
def test_string_cid_is_accepted(db, requests_made):
    assert db.descriptors("702", ['HeavyAtomCount']) == {
        'CID': 702, 'Source': 'rdkit', 'HeavyAtomCount': 3}


@pytest.mark.unit
def test_compound_not_held_locally_has_its_smiles_fetched(db, requests_made):
    """Only the structure goes over the network; the descriptors are still RDKit's."""
    result = db.descriptors(5793, ['HBondDonorCount', 'HeavyAtomCount'])

    assert result == {'CID': 5793, 'Source': 'rdkit',
                      'HBondDonorCount': 5, 'HeavyAtomCount': 12}
    assert requests_made == [([5793], ['SMILES'])]


@pytest.mark.unit
def test_strictly_offline_lookup_of_a_compound_not_held_is_none(db, requests_made):
    assert db.descriptors(5793, use_online_fallback=False) is None
    assert requests_made == []


@pytest.mark.unit
def test_compound_unknown_everywhere_is_none(db, requests_made):
    assert db.descriptors(123456789) is None


@pytest.mark.unit
def test_unreadable_smiles_is_a_record_with_no_values(db, requests_made):
    """The compound exists; RDKit just cannot read it. That is not 'unknown'."""
    assert db.descriptors(999) == {'CID': 999, 'Source': 'rdkit'}
    assert requests_made == []


@pytest.mark.unit
def test_compound_without_structure_is_a_record_with_no_values(db, requests_made):
    assert db.descriptors(998) == {'CID': 998, 'Source': 'rdkit'}
    assert requests_made == []


@pytest.mark.unit
def test_zenodo_descriptor_columns_are_not_what_rdkit_returns(tmp_path, api, requests_made):
    """A Zenodo copy stores PubChem's XLogP and acceptor count; source='rdkit' computes its own."""
    path = _make_database(tmp_path / "zenodo.db", zenodo_columns=True)
    with PubChemID(db_path=path, auto_download=False, api=api) as zenodo:
        result = zenodo.descriptors(2244, ['MolLogP', 'HBondAcceptorCount'])

    assert result == {'CID': 2244, 'Source': 'rdkit',
                      'MolLogP': 1.3101, 'HBondAcceptorCount': 3}
    assert requests_made == []


@pytest.mark.unit
def test_many_compounds_held_locally_cost_no_request(db, requests_made):
    rows = db.descriptors_for_cids([2244, 702, 2244], ['HeavyAtomCount'])

    assert rows == [{'CID': 2244, 'Source': 'rdkit', 'HeavyAtomCount': 13},
                    {'CID': 702, 'Source': 'rdkit', 'HeavyAtomCount': 3}]
    assert requests_made == []


@pytest.mark.unit
def test_bulk_lookup_fetches_only_the_missing_structures(db, requests_made):
    rows = db.descriptors_for_cids([5793, 2244, 702], ['HeavyAtomCount'])

    assert [row['CID'] for row in rows] == [5793, 2244, 702]
    assert {row['Source'] for row in rows} == {'rdkit'}
    assert requests_made == [([5793], ['SMILES'])]


# --------------------------------------------------------------------------
# PubChemID.descriptors, source="pubchem"
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_pubchem_source_returns_pubchems_values_labelled_online(db, requests_made):
    result = db.descriptors(2244, ['XLogP', 'Complexity'], source='pubchem')

    assert result == {'CID': 2244, 'Source': 'online', 'XLogP': 1.2,
                      'Complexity': 212.0}
    assert requests_made == [([2244], ['XLogP', 'Complexity'])]


@pytest.mark.unit
def test_pubchem_source_defaults_to_all_eight(db, requests_made):
    result = db.descriptors(2244, source='pubchem')

    assert set(result) == {'CID', 'Source', *PUBCHEM_DESCRIPTORS}
    assert requests_made == [([2244], list(PUBCHEM_DESCRIPTORS))]


@pytest.mark.unit
def test_pubchem_source_is_asked_for_mollogp_and_says_xlogp(db, requests_made):
    with pytest.raises(ValueError, match="XLogP"):
        db.descriptors(2244, ['MolLogP'], source='pubchem')
    assert requests_made == []


@pytest.mark.unit
def test_pubchem_source_cannot_be_strictly_offline(db, requests_made):
    with pytest.raises(ValueError, match="source='rdkit'"):
        db.descriptors(2244, source='pubchem', use_online_fallback=False)
    assert requests_made == []


@pytest.mark.unit
def test_unknown_source_is_rejected_before_any_request(db, requests_made):
    with pytest.raises(ValueError, match="'rdkit' or 'pubchem'"):
        db.descriptors(2244, source='cactvs')
    assert requests_made == []


@pytest.mark.unit
def test_malformed_cid_is_rejected_before_any_request(db, requests_made):
    with pytest.raises(ValueError, match="integer"):
        db.descriptors_for_cids([2244, "aspirin"])
    assert requests_made == []


# --------------------------------------------------------------------------
# descriptors_table
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_table_has_a_row_for_every_cid_asked_about(db, requests_made):
    table = db.descriptors_table([2244, 123456789, 999], ['TPSA'])

    assert list(table.columns) == ['CID', 'Source', 'TPSA']
    assert table['CID'].tolist() == [2244, 123456789, 999]
    assert table['Source'].tolist() == ['rdkit', 'missing', 'rdkit']
    assert table.loc[0, 'TPSA'] == 63.6
    assert table.loc[1:, 'TPSA'].isna().all()


@pytest.mark.unit
def test_tables_from_either_source_share_their_common_columns(db, requests_made):
    """Switching source keeps the column names wherever the quantity is the same."""
    common = ['TPSA', 'HBondDonorCount', 'HeavyAtomCount']
    rdkit = db.descriptors_table([2244], common)
    pubchem = db.descriptors_table([2244], common, source='pubchem')

    assert list(rdkit.columns) == list(pubchem.columns)
    assert rdkit['Source'].tolist() == ['rdkit']
    assert pubchem['Source'].tolist() == ['online']


@pytest.mark.unit
def test_table_defaults_to_every_descriptor_of_the_source(db, requests_made):
    table = db.descriptors_table([2244])
    assert list(table.columns) == ['CID', 'Source', *RDKIT_DESCRIPTORS]


@pytest.mark.unit
def test_table_rejects_an_unknown_source(db, requests_made):
    with pytest.raises(ValueError, match="'rdkit' or 'pubchem'"):
        db.descriptors_table([2244], source='cactvs')
