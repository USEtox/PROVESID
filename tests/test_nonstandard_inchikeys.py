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


@pytest.mark.integration
def test_zeropm_finds_a_substance_stored_only_under_its_non_standard_key():
    from provesid import ZeroPM
    try:
        zpm = ZeroPM(auto_download=False)
    except FileNotFoundError:  # pragma: no cover - environment dependent
        pytest.skip("ZeroPM is not installed")
    with zpm:
        table = zpm.get_id_table_from_inchikey("ANPMEHLBMHHCML-UHFFFAOYSA-N")
    assert table is not None
    assert "173524-60-2" in table["cas"].tolist()
