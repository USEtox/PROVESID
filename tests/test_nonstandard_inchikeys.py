"""
Non-standard InChIKeys from CompTox and ZeroPM (plan §31.3 item 5).

CompTox stores a non-standard key (flag ``N``, ``...UHFFFAOYNA-N``) for about
11% of its substances, and ZeroPM a non-standard InChI and key for about 5%.
Such a key never equals the standard key PubChem, ChEBI and ChEMBL publish.
Before the fix, ``Search`` returned it as ``InChIKey``, a DTXSID query looked
the other sources up by it and found nothing, and a standard-key query missed
the CompTox row.

The unit tests build a two-row CompTox database in the real schema. The
``Search`` tests use the installed databases and skip without them.
"""

import sqlite3

import pytest

from provesid.comptox import CompToxID
from provesid.search import Search, _cluster_candidates
from provesid.tools import make_candidate, standardize_inchi_and_key
from provesid.utils import inchikey_flag_variants, is_standard_inchikey
from provesid.zeropm import ZeroPM

# Dehydroacetic acid as CompTox stores it, and its standard key.
DEHYDROACETIC_SMILES = "CC(=O)C1C(=O)OC(C)=CC1=O"
DEHYDROACETIC_NONSTANDARD = "PGRHXDWITVMQBC-UHFFFAOYNA-N"
DEHYDROACETIC_STANDARD = "PGRHXDWITVMQBC-UHFFFAOYSA-N"

# DTXSID00891854: the non-standard key differs in the stereo hash too.
STEREO_SMILES = r"COC(=O)[C@@]1(O)[C@@H](Cl)C(=O)C(Cl)=C1\C=C\C"
STEREO_NONSTANDARD = "VFXXEEUGYINLKM-DXABRBBMNA-N"
STEREO_STANDARD = "VFXXEEUGYINLKM-PVXSWLFQSA-N"


@pytest.mark.unit
def test_the_flag_tells_a_standard_key_from_a_non_standard_one():
    assert is_standard_inchikey(DEHYDROACETIC_STANDARD)
    assert not is_standard_inchikey(DEHYDROACETIC_NONSTANDARD)
    assert not is_standard_inchikey(None)
    assert inchikey_flag_variants(DEHYDROACETIC_NONSTANDARD) == [
        DEHYDROACETIC_NONSTANDARD, DEHYDROACETIC_STANDARD,
    ]


@pytest.mark.unit
@pytest.mark.parametrize("smiles,stored,standard", [
    (DEHYDROACETIC_SMILES, DEHYDROACETIC_NONSTANDARD, DEHYDROACETIC_STANDARD),
    (STEREO_SMILES, STEREO_NONSTANDARD, STEREO_STANDARD),
])
def test_a_non_standard_key_is_recomputed_from_the_structure(smiles, stored, standard):
    _, key = standardize_inchi_and_key(smiles, None, stored)
    assert key == standard


@pytest.mark.unit
def test_a_non_standard_inchi_is_recomputed_when_there_is_no_smiles():
    """ZeroPM publishes InChI, not SMILES."""
    inchi, key = standardize_inchi_and_key(
        None, "InChI=1/CH2O/c1-2/h1H2", "WSFSSNUMVMOOMR-UHFFFAOYNA-N"
    )
    assert inchi == "InChI=1S/CH2O/c1-2/h1H2"
    assert key == "WSFSSNUMVMOOMR-UHFFFAOYSA-N"


@pytest.mark.unit
def test_a_non_standard_key_without_a_readable_structure_becomes_none():
    """A non-standard key is never passed on as if it could match."""
    assert standardize_inchi_and_key("not a smiles", None, DEHYDROACETIC_NONSTANDARD) == (None, None)
    assert standardize_inchi_and_key(None, None, DEHYDROACETIC_NONSTANDARD) == (None, None)


@pytest.mark.unit
def test_standard_values_are_not_touched():
    """Nothing is recomputed, even from a structure that disagrees."""
    assert standardize_inchi_and_key("CCO", "InChI=1S/CH2O/c1-2/h1H2", DEHYDROACETIC_STANDARD) == (
        "InChI=1S/CH2O/c1-2/h1H2", DEHYDROACETIC_STANDARD,
    )


@pytest.mark.unit
def test_a_comptox_candidate_clusters_with_a_standard_one_without_skeletons():
    comptox = make_candidate("CompTox", smiles=DEHYDROACETIC_SMILES, inchikey=DEHYDROACETIC_NONSTANDARD)
    pubchem = make_candidate("PubChemID", smiles=DEHYDROACETIC_SMILES, inchikey=DEHYDROACETIC_STANDARD)
    assert comptox["InChIKey"] == DEHYDROACETIC_STANDARD
    clusters = _cluster_candidates([comptox, pubchem], by_skeleton=False)
    assert len(clusters) == 1


# ── CompToxID.get_by_inchikey, on a small database in the real schema ───────

_COLUMNS = (
    "DTXSID", "PREFERRED_NAME", "CASRN", "DTXCID", "INCHIKEY", "IUPAC_NAME", "SMILES",
    "MOLECULAR_FORMULA", "AVERAGE_MASS", "MONOISOTOPIC_MASS", "QSAR_READY_SMILES",
    "MS_READY_SMILES", "IDENTIFIER",
)


@pytest.fixture
def comptox(tmp_path):
    path = tmp_path / "comptox_chemicals.db"
    rows = [
        {"DTXSID": "DTXSID6020014", "PREFERRED_NAME": "Dehydroacetic acid",
         "CASRN": "520-45-6", "INCHIKEY": DEHYDROACETIC_NONSTANDARD, "SMILES": DEHYDROACETIC_SMILES},
        {"DTXSID": "DTXSID0000009", "PREFERRED_NAME": "Stored twice",
         "CASRN": "0-00-9", "INCHIKEY": "VFXXEEUGYINLKM-DXABRBBMSA-N"},
        {"DTXSID": "DTXSID0000010", "PREFERRED_NAME": "Stored twice, non-standard",
         "CASRN": "0-01-0", "INCHIKEY": "VFXXEEUGYINLKM-DXABRBBMNA-N"},
    ]
    connection = sqlite3.connect(path)
    connection.execute(f"CREATE TABLE chemicals ({', '.join(c + ' TEXT' for c in _COLUMNS)})")
    connection.executemany(
        f"INSERT INTO chemicals VALUES ({', '.join('?' for _ in _COLUMNS)})",
        [tuple(row.get(c) for c in _COLUMNS) for row in rows],
    )
    connection.commit()
    connection.close()
    with CompToxID(db_path=str(path), auto_download=False) as client:
        yield client


@pytest.mark.unit
def test_a_standard_key_finds_the_row_stored_under_the_non_standard_one(comptox):
    row = comptox.get_by_inchikey(DEHYDROACETIC_STANDARD)
    assert row["DTXSID"] == "DTXSID6020014"
    assert row["INCHIKEY"] == DEHYDROACETIC_NONSTANDARD  # the record as stored
    assert comptox.inchikey_to_dtxsid(DEHYDROACETIC_STANDARD) == "DTXSID6020014"


@pytest.mark.unit
def test_the_exact_spelling_wins_when_both_are_stored(comptox):
    assert comptox.get_by_inchikey("VFXXEEUGYINLKM-DXABRBBMSA-N")["DTXSID"] == "DTXSID0000009"
    assert comptox.get_by_inchikey("VFXXEEUGYINLKM-DXABRBBMNA-N")["DTXSID"] == "DTXSID0000010"


# ── Search, on the installed databases ─────────────────────────────────────

@pytest.fixture(scope="module")
def installed():
    s = Search("cas", show_progress=False)
    s._ensure_clients()
    missing = [k for k in ("chebi", "comptox", "pubchem") if k in s.sources_unavailable]
    s.close()
    if missing:  # pragma: no cover - environment dependent
        pytest.skip(f"Offline sources unavailable: {', '.join(missing)}")


@pytest.mark.integration
@pytest.mark.parametrize("dtxsid,standard", [
    ("DTXSID6020014", DEHYDROACETIC_STANDARD),         # dehydroacetic acid
    ("DTXSID8020040", "QBYJBZPUGVGKQQ-SJJAEHHWSA-N"),  # aldrin
])
def test_a_dtxsid_query_reaches_the_other_sources(installed, dtxsid, standard):
    """CompTox's key used to be the one the others were asked for."""
    with Search("dtxsid", show_progress=False) as s:
        row = s.search([dtxsid]).iloc[0]
    assert row["InChIKey"] == standard
    assert row["n_source_support"] >= 3


@pytest.mark.integration
def test_a_standard_key_query_finds_the_comptox_row(installed):
    with Search("inchikey", show_progress=False) as s:
        row = s.search([DEHYDROACETIC_STANDARD]).iloc[0]
    assert row["DTXSID"] == "DTXSID6020014"
    assert row["source_details"]["CompTox"]["found"]


@pytest.mark.integration
def test_comptox_is_not_split_off_without_skeleton_clustering(installed):
    with Search("cas", cluster_by_skeleton=False, n_hits="all", show_progress=False) as s:
        df = s.search(["520-45-6"])
    assert len(df) == 1
    assert df.iloc[0]["InChIKey"] == DEHYDROACETIC_STANDARD
    assert df.iloc[0]["n_source_support"] >= 3


# ── ZeroPM's InChI lookups ───────────────────────────────────────────────────

# trans-1,4-Cyclohexanediol (6995-79-5): ZeroPM stores only the non-standard
# InChI, which differs from the standard one by its prefix alone.
DIOL_STANDARD_INCHI = "InChI=1S/C6H12O2/c7-5-1-2-6(8)4-3-5/h5-8H,1-4H2/t5-,6-"
DIOL_STORED_INCHI = "InChI=1/C6H12O2/c7-5-1-2-6(8)4-3-5/h5-8H,1-4H2/t5-,6-"
DIOL_STANDARD_KEY = "VKONPUDBRVKQLM-IZLXSQMJSA-N"
DIOL_STORED_KEY = "VKONPUDBRVKQLM-IZLXSQMJNA-N"
ETHANOL_INCHI = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
ETHANOL_KEY = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"


@pytest.fixture
def zeropm_substances():
    """A ZeroPM client over a ``substances`` table alone, in memory."""
    connection = sqlite3.connect(":memory:")
    connection.execute(
        "CREATE TABLE substances (inchi_id INTEGER PRIMARY KEY, inchi TEXT, inchikey TEXT)"
    )
    connection.executemany("INSERT INTO substances VALUES (?, ?, ?)", [
        (1, DIOL_STORED_INCHI, DIOL_STORED_KEY),
        (2, ETHANOL_INCHI, ETHANOL_KEY),
    ])
    client = object.__new__(ZeroPM)
    client._adopt_connection(connection, db_path=":memory:")
    with client:
        yield client


@pytest.mark.unit
@pytest.mark.parametrize("inchi,inchi_id", [
    (DIOL_STORED_INCHI, 1),                         # as stored
    (DIOL_STANDARD_INCHI, 1),                       # standard, stored non-standard
    ("InChI=1/C2H6O/c1-2-3/h3H,2H2,1H3", 2),        # non-standard, stored standard
    ("InChI=1S/C2H6O/c1-3-2/h1H3", None),           # dimethyl ether: not there
    ("not an InChI", None),
    ("InChI=1S/garbage", None),
])
def test_an_inchi_is_found_by_string_or_by_its_key(zeropm_substances, inchi, inchi_id):
    found = zeropm_substances._find_substance_by_inchi(inchi)
    assert (found[0] if found else None) == inchi_id


@pytest.mark.unit
def test_the_row_found_by_key_is_returned_as_stored(zeropm_substances):
    assert tuple(zeropm_substances._find_substance_by_inchi(DIOL_STANDARD_INCHI)) == (
        1, DIOL_STORED_INCHI, DIOL_STORED_KEY,
    )


@pytest.fixture(scope="module")
def zeropm():
    try:
        client = ZeroPM(auto_download=False)
    except FileNotFoundError:  # pragma: no cover - environment dependent
        pytest.skip("ZeroPM is not installed")
    with client:
        yield client


@pytest.mark.integration
def test_zeropm_finds_a_substance_stored_only_under_its_non_standard_key(zeropm):
    table = zeropm.get_id_table_from_inchikey("ANPMEHLBMHHCML-UHFFFAOYSA-N")
    assert table is not None
    assert "173524-60-2" in table["cas"].tolist()


@pytest.mark.integration
def test_zeropm_finds_a_substance_stored_only_under_its_non_standard_inchi(zeropm):
    table = zeropm.get_id_table_from_inchi(DIOL_STANDARD_INCHI)
    assert table["cas"].dropna().tolist() == ["6995-79-5"]
    assert set(table["inchi"]) == {DIOL_STORED_INCHI}
    assert zeropm.get_cas_from_inchi(DIOL_STANDARD_INCHI) == "6995-79-5"


@pytest.mark.integration
def test_every_zeropm_key_lookup_tries_both_flags(zeropm):
    assert zeropm.get_cas_from_inchikey(DIOL_STANDARD_KEY) == "6995-79-5"
    assert zeropm.get_smiles_from_inchikey(DIOL_STANDARD_KEY) == "O[C@H]1CC[C@H](O)CC1"
    assert zeropm.batch_get_cas_from_inchikey([DIOL_STANDARD_KEY, "XXXXXXXXXXXXXX-UHFFFAOYSA-N"]) == {
        DIOL_STANDARD_KEY: "6995-79-5", "XXXXXXXXXXXXXX-UHFFFAOYSA-N": None,
    }


@pytest.mark.integration
def test_relative_stereo_is_still_out_of_reach(zeropm):
    """Stored as ``/t9-,10+,11-,12-/s2``: no standard InChI hashes the same."""
    standard = ("InChI=1S/C13H20O2/c1-4-15-13(14)12-10-6-5-9(7-10)11(12)8(2)3"
                "/h5-6,8-12H,4,7H2,1-3H3/t9-,10+,11-,12-/m0/s1")
    assert zeropm.get_id_table_from_inchi(standard) is None


@pytest.mark.integration
def test_search_by_inchi_reaches_zeropm(installed, zeropm):
    with Search("inchi", use_zeropm=True, show_progress=False) as s:
        row = s.search([DIOL_STANDARD_INCHI]).iloc[0]
    assert row["source_details"]["ZeroPM"]["found"]
    assert row["InChIKey"] == DIOL_STANDARD_KEY
