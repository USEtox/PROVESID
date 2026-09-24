"""
Tests for building ``pubchem_id.db`` from PubChem's FTP files.

Nothing is stubbed between the builder and the network. A miniature PubChem
release --- the real directory layout, the real file formats, gzipped, with an
``.md5`` beside every file --- is written to a temporary directory and served
over HTTP on localhost, so listing the releases, downloading, checking MD5s,
streaming, filtering and writing the database all run as they would against
``ftp.ncbi.nlm.nih.gov``.

The release is chosen to exercise the decisions the builder makes:

* aspirin (2244) has two CAS rows for the same number, and cross-references of
  every stored type plus one (Wikidata) that is not stored;
* ethanol (702) appears *before* aspirin in some files, so nothing may assume
  the files are sorted;
* chloroform-d (71583) is written ``CHCl3`` by PubChem, so its weight must come
  from the isotope label in its SMILES rather than its formula;
* CID 999001 has only a CAS number that fails the check digit, and CID 999002
  has no CAS number at all --- neither may appear anywhere in the database,
  although both have a line in every file.
"""

import functools
import gzip
import hashlib
import http.server
import os
import sqlite3
import threading

import pytest

from provesid import pubchem_ftp
from provesid.datasets import DownloadError
from provesid.pubchem_id import PubChemID
from provesid.pubchem_ftp import (
    build_pubchem_id_db,
    extras_url,
    list_releases,
    molecular_weight,
    resolve_release,
    selected_files,
)
from provesid.utils import check_CASRN

RELEASE = "2026-09-01"

ASPIRIN_INCHI = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
ETHANOL_INCHI = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
CHLOROFORM_D_INCHI = "InChI=1S/CHCl3/c2-1(3)4/h1H/i1D"

#: One miniature release: file name to lines, as PubChem writes them.
FILES = {
    "CID-Identifiers.tsv.gz": [
        "702\t64-17-5\tCAS",
        "702\t200-578-6\tEuropean Community (EC) Number",
        "702\tDTXSID9020584\tDSSTox Substance ID",
        "2244\t50-78-2\tCAS",
        "2244\t50-78-2\tCAS",
        "2244\tCHEBI:15365\tChEBI ID",
        "2244\tCHEMBL25\tChEMBL ID",
        "2244\tDTXSID5020108\tDSSTox Substance ID",
        "2244\t200-064-1\tEuropean Community (EC) Number",
        "2244\tR16CO5Y76E\tUNII",
        "2244\tQ18216\tWikidata",
        "71583\t865-49-6\tCAS",
        "999001\t50-78-3\tCAS",
        "999001\tCHEBI:1\tChEBI ID",
        "999002\tCHEBI:2\tChEBI ID",
    ],
    "CID-Date.gz": [
        "702\t2004-09-16", "2244\t2004-09-16", "71583\t2005-03-26",
        "999001\t2020-01-01", "999002\t2020-01-01",
    ],
    "CID-Mass.gz": [
        "2244\tC9H8O4\t180.04225873\t180.04225873",
        "702\tC2H6O\t46.041864811\t46.041864811",
        "71583\tCHCl3\t118.9175601\t118.9175601",
        "999001\tC6H6\t78.04695\t78.04695",
        "999002\tC6H6\t78.04695\t78.04695",
    ],
    "CID-SMILES.gz": [
        "702\tCCO",
        "2244\tCC(=O)OC1=CC=CC=C1C(=O)O",
        "71583\t[2H]C(Cl)(Cl)Cl",
        "999001\tC1=CC=CC=C1",
        "999002\tC1=CC=CC=C1",
    ],
    "CID-Title.gz": [
        "2244\tAspirin", "702\tEthanol", "71583\tChloroform-d",
        "999001\tBenzene?", "999002\tBenzene!",
    ],
    "CID-IUPAC.gz": [
        "2244\t2-acetyloxybenzoic acid", "702\tethanol",
        "71583\ttrichloro(deuterio)methane", "999001\tbenzene", "999002\tbenzene",
    ],
    "CID-InChI-Key.gz": [
        f"2244\t{ASPIRIN_INCHI}\tBSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        f"702\t{ETHANOL_INCHI}\tLFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        f"71583\t{CHLOROFORM_D_INCHI}\tHEDRZPFGACZZDS-MICDWDOJSA-N",
        "999001\tInChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H\tUHOVQNZJYSORNB-UHFFFAOYSA-N",
        "999002\tInChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H\tUHOVQNZJYSORNB-UHFFFAOYSA-N",
    ],
    "CID-Synonym-filtered.gz": [
        "2244\tAspirin", "2244\tacetylsalicylic acid", "2244\t50-78-2",
        "702\tEthanol", "702\tethyl alcohol",
        "71583\tChloroform-d",
        "999001\tbenzene", "999002\tbenzene",
    ],
}


def _write_release(directory, files):
    """Write gzipped files, each with the ``.md5`` PubChem publishes beside it."""
    os.makedirs(directory, exist_ok=True)
    for name, lines in files.items():
        path = os.path.join(directory, name)
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write("\n".join(lines) + "\n")
        with open(path, "rb") as handle:
            digest = hashlib.md5(handle.read()).hexdigest()
        with open(path + ".md5", "w") as handle:
            handle.write(f"{digest}  {name}\n")


class _Server:
    """A static file server on localhost that remembers what was asked for."""

    def __init__(self, root):
        self.root = root
        requests_seen = self.requests = []

        class Handler(http.server.SimpleHTTPRequestHandler):
            def do_GET(self):
                requests_seen.append(self.path)
                super().do_GET()

            def log_message(self, *args):
                pass

        self.httpd = http.server.ThreadingHTTPServer(
            ("127.0.0.1", 0), functools.partial(Handler, directory=str(root)))
        self.thread = threading.Thread(target=self.httpd.serve_forever, daemon=True)
        self.thread.start()
        self.url = f"http://127.0.0.1:{self.httpd.server_address[1]}"

    def fetched(self, filename):
        """How many times a file (not its checksum) was downloaded."""
        return sum(1 for path in self.requests if path.endswith("/" + filename))

    def close(self):
        self.httpd.shutdown()
        self.httpd.server_close()


@pytest.fixture
def ftp(tmp_path):
    """
    A miniature PubChem compound root, served over HTTP.

    ``Monthly/`` holds two snapshots, of which only the newer has files;
    ``Extras/`` at the root is the rolling dump. Yields the server; its
    ``url`` goes where :data:`provesid.pubchem_ftp.FTP_ROOT` would.
    """
    root = tmp_path / "pubchem"
    _write_release(root / "Monthly" / RELEASE / "Extras", FILES)
    (root / "Monthly" / RELEASE / "TIMESTAMP").write_text("2026/08/31 18:38:04\n")
    (root / "Monthly" / "2026-08-01").mkdir(parents=True)
    _write_release(root / "Extras", FILES)
    server = _Server(root)
    yield server
    server.close()


@pytest.fixture
def built(ftp, tmp_path):
    """The database built from the miniature release with default options."""
    path = build_pubchem_id_db(str(tmp_path / "data" / "pubchem_id.db"),
                               base_url=ftp.url, progress=False)
    return path


def _query(path, sql, *params):
    connection = sqlite3.connect(path)
    try:
        return connection.execute(sql, params).fetchall()
    finally:
        connection.close()


# ──────────────────────────────────────────────────────────────────────────────
# Releases
# ──────────────────────────────────────────────────────────────────────────────

@pytest.mark.unit
def test_list_releases_is_newest_first_with_current_last(ftp):
    assert list_releases(base_url=ftp.url) == [RELEASE, "2026-08-01", "current"]


@pytest.mark.unit
def test_latest_resolves_to_the_newest_snapshot(ftp):
    assert resolve_release("latest", base_url=ftp.url) == RELEASE


@pytest.mark.unit
def test_a_named_release_costs_no_request(ftp):
    assert resolve_release("2026-08-01", base_url=ftp.url) == "2026-08-01"
    assert resolve_release("current", base_url=ftp.url) == "current"
    assert ftp.requests == []


@pytest.mark.unit
@pytest.mark.parametrize("release", ["2026-9-1", "september", "", "Latest"])
def test_a_bad_release_is_named_as_such(release):
    with pytest.raises(ValueError, match="is not a PubChem release"):
        resolve_release(release)


@pytest.mark.unit
def test_extras_url_for_snapshots_and_the_rolling_dump():
    assert extras_url(RELEASE, base_url="http://x/pc") == f"http://x/pc/Monthly/{RELEASE}/Extras"
    assert extras_url("current", base_url="http://x/pc/") == "http://x/pc/Extras"


@pytest.mark.unit
def test_selected_files_honour_the_options():
    assert [f.key for f in selected_files()] == [
        "identifiers", "date", "mass", "smiles", "title", "iupac", "inchi", "synonyms"]
    assert [f.key for f in selected_files(include_inchi=False, include_synonyms=False)] == [
        "identifiers", "date", "mass", "smiles", "title", "iupac"]


# ──────────────────────────────────────────────────────────────────────────────
# Molecular weight
# ──────────────────────────────────────────────────────────────────────────────

@pytest.mark.unit
@pytest.mark.parametrize("formula, expected", [
    ("C9H8O4", 180.16),        # aspirin; PubChem: 180.16
    ("C9H18NO4+", 204.24),     # PubChem: 204.24 -- IUPAC 2005 H and C, not 2013
    ("C12H22N2O2", 226.32),    # PubChem: 226.32
    ("C6H13O9P", 260.14),      # PubChem: 260.14 -- rounds half up, not to even
    ("H2Se", 80.99),           # PubChem: 80.99 -- Se at 78.971, not 78.96
    ("AsHO4-2", 139.93),       # a charge suffix is not an element
    ("Na+", 22.99),
])
def test_molecular_weight_matches_pubchem(formula, expected):
    assert molecular_weight(formula) == expected


@pytest.mark.unit
@pytest.mark.parametrize("formula", ["", None, "c9h8o4", "C9H8O4)", "Xx2"])
def test_molecular_weight_of_a_bad_formula_is_none(formula):
    assert molecular_weight(formula) is None


# ──────────────────────────────────────────────────────────────────────────────
# The build
# ──────────────────────────────────────────────────────────────────────────────

@pytest.mark.unit
def test_only_compounds_with_a_valid_cas_number_are_in_scope(built):
    assert not check_CASRN("50-78-3")
    cids = [row[0] for row in _query(built, "SELECT cid FROM compounds ORDER BY cid")]
    assert cids == [702, 2244, 71583]
    for table in ("cas_numbers", "synonyms", "xrefs"):
        assert _query(built, f"SELECT COUNT(*) FROM {table} WHERE cid > 999000") == [(0,)]


@pytest.mark.unit
def test_cas_numbers_are_curated_and_deduplicated(built):
    rows = _query(built, "SELECT cid, cas FROM cas_numbers ORDER BY cid")
    assert rows == [(702, "64-17-5"), (2244, "50-78-2"), (71583, "865-49-6")]


@pytest.mark.unit
def test_every_column_is_filled_from_its_file(built):
    connection = sqlite3.connect(built)
    connection.row_factory = sqlite3.Row
    aspirin = dict(connection.execute("SELECT * FROM compounds WHERE cid = 2244").fetchone())
    connection.close()
    assert aspirin == {
        "cid": 2244,
        "cmpdname": "Aspirin",
        "mf": "C9H8O4",
        "inchi": ASPIRIN_INCHI,
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        "iupacname": "2-acetyloxybenzoic acid",
        "mw": 180.16,
        "exactmass": 180.04225873,
        "monoisotopicmass": 180.04225873,
        "cidcdate": "2004-09-16",
    }


@pytest.mark.unit
def test_an_isotope_label_weighs_what_it_is(built):
    """
    PubChem writes chloroform-d as CHCl3; it is one neutron heavier than that.

    PubChem reports 119.37 for chloroform and 120.38 for chloroform-d, two
    weights no single chlorine weight reproduces together; this builder gives
    119.37 and 120.37.
    """
    assert molecular_weight("CHCl3") == 119.37
    assert _query(built, "SELECT mw FROM compounds WHERE cid = 71583") == [(120.37,)]


@pytest.mark.unit
def test_synonyms_keep_pubchems_order(built):
    rows = _query(built, "SELECT synonym FROM synonyms WHERE cid = 2244 ORDER BY id")
    assert [row[0] for row in rows] == ["Aspirin", "acetylsalicylic acid", "50-78-2"]


@pytest.mark.unit
def test_cross_references_keep_only_the_stored_types(built):
    rows = _query(built, "SELECT source, identifier FROM xrefs WHERE cid = 2244 "
                         "ORDER BY source")
    assert rows == [("chebi", "CHEBI:15365"), ("chembl", "CHEMBL25"),
                    ("dtxsid", "DTXSID5020108"), ("ec", "200-064-1"),
                    ("unii", "R16CO5Y76E")]


@pytest.mark.unit
def test_the_indexes_the_lookups_rely_on_exist(built):
    names = {row[0] for row in _query(built, "SELECT name FROM sqlite_master "
                                             "WHERE type = 'index'")}
    assert {"idx_compounds_inchikey", "idx_cas_cas", "idx_synonyms_synonym",
            "idx_xrefs_cid"} <= names


@pytest.mark.unit
def test_provenance_names_the_release_and_every_file(built, ftp):
    entries = dict(_query(built, "SELECT key, value FROM provenance"))
    assert entries["release"] == RELEASE
    assert entries["release_timestamp"] == "2026/08/31 18:38:04"
    assert entries["source_url"] == f"{ftp.url}/Monthly/{RELEASE}/Extras"
    assert entries["inchi_source"] == "pubchem"
    assert entries["compounds"] == "3"
    assert entries["cas_rows_rejected"] == "1"
    assert entries["mw_from_isotopes"] == "1"

    files = _query(built, "SELECT file, url, md5, lines_read FROM provenance_files")
    assert [row[0] for row in files] == [f.filename for f in selected_files()]
    for name, url, md5, lines in files:
        assert url == f"{ftp.url}/Monthly/{RELEASE}/Extras/{name}"
        published = ftp.root / "Monthly" / RELEASE / "Extras" / f"{name}.md5"
        assert md5 == published.read_text().split()[0]
        assert lines == len(FILES[name])


@pytest.mark.unit
def test_downloads_are_deleted_and_nothing_temporary_is_left(built):
    assert os.listdir(os.path.dirname(built)) == ["pubchem_id.db"]


@pytest.mark.unit
def test_kept_downloads_are_reused_rather_than_fetched_again(ftp, tmp_path):
    path = str(tmp_path / "pubchem_id.db")
    build_pubchem_id_db(path, base_url=ftp.url, keep_downloads=True, progress=False)
    assert ftp.fetched("CID-SMILES.gz") == 1
    kept = tmp_path / "pubchem_ftp" / RELEASE / "CID-SMILES.gz"
    assert kept.exists()

    build_pubchem_id_db(path, base_url=ftp.url, keep_downloads=True, force=True,
                        progress=False)
    assert ftp.fetched("CID-SMILES.gz") == 1


@pytest.mark.unit
def test_a_damaged_kept_download_is_fetched_again(ftp, tmp_path):
    path = str(tmp_path / "pubchem_id.db")
    build_pubchem_id_db(path, base_url=ftp.url, keep_downloads=True, progress=False)
    (tmp_path / "pubchem_ftp" / RELEASE / "CID-Title.gz").write_bytes(b"not gzip")

    build_pubchem_id_db(path, base_url=ftp.url, keep_downloads=True, force=True,
                        progress=False)
    assert ftp.fetched("CID-Title.gz") == 2
    assert _query(path, "SELECT cmpdname FROM compounds WHERE cid = 2244") == [("Aspirin",)]


@pytest.mark.unit
def test_an_existing_database_is_not_replaced_without_force(built, ftp):
    ftp.requests.clear()
    with pytest.raises(FileExistsError, match="force=True"):
        build_pubchem_id_db(built, base_url=ftp.url, progress=False)
    assert ftp.requests == []


@pytest.mark.unit
def test_a_failed_build_leaves_the_existing_database_alone(built, ftp):
    """A file that fails its checksum stops the build; the old database stays."""
    before = open(built, "rb").read()
    # Corrupt one file on the server without updating its published MD5.
    with open(ftp.root / "Monthly" / RELEASE / "Extras" / "CID-IUPAC.gz", "ab") as handle:
        handle.write(b"tampered")

    with pytest.raises(DownloadError, match="Checksum mismatch"):
        build_pubchem_id_db(built, base_url=ftp.url, force=True, progress=False)

    assert open(built, "rb").read() == before
    assert not os.path.exists(built + ".tmp")


@pytest.mark.unit
def test_without_inchi_the_file_is_skipped_and_rdkit_fills_the_columns(ftp, tmp_path):
    path = build_pubchem_id_db(str(tmp_path / "pubchem_id.db"), base_url=ftp.url,
                               include_inchi=False, progress=False)
    assert ftp.fetched("CID-InChI-Key.gz") == 0
    assert _query(path, "SELECT inchi, inchikey FROM compounds WHERE cid = 2244") == [
        (ASPIRIN_INCHI, "BSYNRYMUTXBXSQ-UHFFFAOYSA-N")]
    assert dict(_query(path, "SELECT key, value FROM provenance"))["inchi_source"] == "rdkit"


@pytest.mark.unit
def test_without_synonyms_the_file_is_skipped(ftp, tmp_path):
    path = build_pubchem_id_db(str(tmp_path / "pubchem_id.db"), base_url=ftp.url,
                               include_synonyms=False, progress=False)
    assert ftp.fetched("CID-Synonym-filtered.gz") == 0
    assert _query(path, "SELECT COUNT(*) FROM synonyms") == [(0,)]


@pytest.mark.unit
def test_the_rolling_dump_is_read_from_the_root_extras(ftp, tmp_path):
    path = build_pubchem_id_db(str(tmp_path / "pubchem_id.db"), base_url=ftp.url,
                               release="current", progress=False)
    assert all(request.startswith("/Extras/") for request in ftp.requests)
    entries = dict(_query(path, "SELECT key, value FROM provenance"))
    assert entries["release"] == "current"
    assert entries["release_timestamp"] == ""


# ──────────────────────────────────────────────────────────────────────────────
# PubChemID on a database built from FTP
# ──────────────────────────────────────────────────────────────────────────────

@pytest.fixture
def ftp_default(ftp, monkeypatch):
    """Point the builder's default root at the miniature release."""
    monkeypatch.setattr(pubchem_ftp, "FTP_ROOT", ftp.url)
    return ftp


@pytest.mark.unit
def test_pubchemid_builds_from_ftp_by_default(ftp_default, tmp_path):
    with PubChemID(data_dir=str(tmp_path)) as db:
        assert db.source == "ftp"
        assert db.cas_to_cid("50-78-2") == 2244
        assert db.get_by_inchikey("LFQSCWFLJHTTHZ-UHFFFAOYSA-N")["cmpdname"] == "Ethanol"
        assert db.search_by_name("acetylsalicylic acid", exact=True)[0]["cid"] == 2244
    assert ftp_default.fetched("CID-SMILES.gz") == 1


@pytest.mark.unit
def test_pubchemid_serves_monoisotopic_mass_and_not_descriptors_offline(
        ftp_default, tmp_path):
    with PubChemID(data_dir=str(tmp_path)) as db:
        assert "MonoisotopicMass" in db.offline_properties
        assert "XLogP" not in db.offline_properties
        result = db.properties(2244, ["MonoisotopicMass", "MolecularWeight"],
                               use_online_fallback=False)
        assert result == {"CID": 2244, "Source": "offline",
                          "MonoisotopicMass": 180.04225873, "MolecularWeight": 180.16}


@pytest.mark.unit
def test_pubchemid_reports_provenance_and_xrefs(ftp_default, tmp_path):
    with PubChemID(data_dir=str(tmp_path)) as db:
        record = db.provenance()
        assert record["release"] == RELEASE
        assert [f["file"] for f in record["files"]][0] == "CID-Identifiers.tsv.gz"
        assert db.xrefs(2244)["dtxsid"] == ["DTXSID5020108"]
        assert db.xrefs("702") == {"dtxsid": ["DTXSID9020584"], "ec": ["200-578-6"]}
        assert db.xrefs(71583) == {}


@pytest.mark.unit
def test_pubchemid_batch_tables_follow_the_databases_columns(ftp_default, tmp_path):
    with PubChemID(data_dir=str(tmp_path)) as db:
        table = db.get_by_cas_batch(["50-78-2", "0-00-0"])
        assert list(table["cas"]) == ["50-78-2"]
        assert "monoisotopicmass" in table.columns
        assert "xlogp" not in table.columns
        assert list(db.get_by_cas_batch([]).columns) == list(table.columns)
        by_smiles = db.get_by_smiles_batch(["CCO"])
        assert list(by_smiles["cas"]) == ["64-17-5"]


@pytest.mark.unit
def test_pubchemid_zenodo_route_downloads(monkeypatch, tmp_path):
    calls = []

    def fake_download(db_path=None, zenodo_url=None, force=False):
        calls.append((db_path, zenodo_url, force))
        sqlite3.connect(db_path).executescript(
            "CREATE TABLE compounds (cid INTEGER PRIMARY KEY, mf TEXT, xlogp REAL);")

    monkeypatch.setattr(PubChemID, "download_database", staticmethod(fake_download))
    with PubChemID(data_dir=str(tmp_path), source="zenodo") as db:
        assert db.offline_properties == {"MolecularFormula": "mf"}
        assert db.provenance() == {}
        with pytest.raises(RuntimeError, match="build_pubchem_id_db"):
            db.xrefs(2244)
    assert calls == [(str(tmp_path / "pubchem_id.db"), PubChemID.DEFAULT_DB_URL, False)]


@pytest.mark.unit
def test_pubchemid_rejects_an_unknown_source_before_anything_else(tmp_path):
    with pytest.raises(ValueError, match="'ftp', 'zenodo'"):
        PubChemID(data_dir=str(tmp_path), source="sdf")
    assert list(tmp_path.iterdir()) == []


@pytest.mark.unit
def test_pubchemid_opens_whatever_is_on_disk_regardless_of_source(built, monkeypatch):
    """``source`` describes an acquisition; an existing database is just opened."""
    def refuse(*args, **kwargs):
        raise AssertionError("nothing should be fetched")

    monkeypatch.setattr(pubchem_ftp, "build_pubchem_id_db", refuse)
    monkeypatch.setattr(PubChemID, "download_database", staticmethod(refuse))
    for source in PubChemID.SOURCES:
        with PubChemID(db_path=built, source=source) as db:
            assert db.cas_to_cid("64-17-5") == 702


@pytest.mark.unit
def test_the_dataset_manager_counts_and_removes_kept_downloads(ftp, tmp_path):
    """Kept source files are the dataset's; ``datasets.remove`` takes them too."""
    from provesid import datasets

    data = tmp_path / "data"
    build_pubchem_id_db(str(data / "pubchem_id.db"), base_url=ftp.url,
                        keep_downloads=True, progress=False)
    row = datasets.status("pubchem", data_dir=data).iloc[0]
    assert row["present"]
    assert row["files"] == 1 + len(selected_files())

    removed = datasets.remove("pubchem", data_dir=data)
    assert removed["removed"].all()
    assert not any(path.is_file() for path in data.rglob("*"))


@pytest.mark.unit
def test_a_download_dir_elsewhere_loses_only_itself(ftp, tmp_path):
    """Cleaning up an empty download directory must not climb out of it."""
    outside = tmp_path / "scratch" / "empty-parent" / "downloads"
    build_pubchem_id_db(str(tmp_path / "data" / "pubchem_id.db"), base_url=ftp.url,
                        download_dir=str(outside), progress=False)
    assert not outside.exists()
    assert outside.parent.exists()
