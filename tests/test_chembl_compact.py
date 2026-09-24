"""
Tests for ``CheMBL.compact`` — the ChEMBL extract.

A full ChEMBL release is ~30 GB across 74 tables, of which PROVESID reads eight.
``compact`` copies those eight, drops the ``molfile`` column and rebuilds only
the indexes this package's queries need, producing ~2.6 GB that answers
identically.

These tests build a miniature ChEMBL — the eight real tables with a handful of
rows each, plus a bioactivity table standing in for the 66 the package ignores —
so they need neither the 30 GB database nor network access.  The equivalence
check against a real release lives at the bottom, marked ``slow``, and skips
when no full database is present.
"""

import glob
import inspect
import logging
import os
import shutil
import sqlite3
import tarfile
from types import SimpleNamespace

import pytest

from provesid import CheMBL, ChEMBLError


# ── The miniature ChEMBL ──────────────────────────────────────────────────────

def _make_full_chembl(path, n_compounds=25):
    """Create a small but structurally faithful stand-in for a ChEMBL release.

    Carries every table and column ``compact`` reads, a ``molfile`` column with
    enough bulk to be visible in a size comparison, and an ``activities`` table
    representing the bioactivity half of ChEMBL that the package never opens.
    """
    conn = sqlite3.connect(path)
    conn.executescript(
        """
        CREATE TABLE molecule_dictionary (
            molregno INTEGER PRIMARY KEY, pref_name TEXT, chembl_id TEXT,
            max_phase NUMERIC, therapeutic_flag INTEGER, molecule_type TEXT);
        CREATE TABLE compound_structures (
            molregno INTEGER PRIMARY KEY, molfile TEXT, standard_inchi TEXT,
            standard_inchi_key TEXT, canonical_smiles TEXT);
        CREATE TABLE compound_properties (
            molregno INTEGER PRIMARY KEY, mw_freebase REAL, alogp REAL,
            hba INTEGER, hbd INTEGER, psa REAL, full_molformula TEXT);
        CREATE TABLE molecule_synonyms (
            molregno INTEGER, syn_type TEXT, molsyn_id INTEGER, synonyms TEXT);
        CREATE TABLE molecule_hierarchy (
            molregno INTEGER PRIMARY KEY, parent_molregno INTEGER,
            active_molregno INTEGER);
        CREATE TABLE chembl_id_lookup (
            chembl_id TEXT PRIMARY KEY, entity_type TEXT, entity_id INTEGER,
            status TEXT, last_active INTEGER);
        CREATE TABLE pesticide_classification (
            pest_class_id INTEGER PRIMARY KEY, compound_name TEXT, mec_id INTEGER,
            mechanism_comment TEXT, ref_type TEXT, ref_id TEXT, ref_url TEXT);
        CREATE TABLE pesticide_class_mapping (
            mol_pest_id INTEGER PRIMARY KEY, pest_class_id INTEGER,
            molregno INTEGER);
        -- One of the 66 tables PROVESID never reads.
        CREATE TABLE activities (activity_id INTEGER PRIMARY KEY, molregno INTEGER,
            standard_value REAL, standard_units TEXT);
        """
    )
    for i in range(1, n_compounds + 1):
        conn.execute(
            "INSERT INTO molecule_dictionary VALUES (?,?,?,?,?,?)",
            (i, f"COMPOUND {i}", f"CHEMBL{i}", 4, 1, "Small molecule"),
        )
        conn.execute(
            "INSERT INTO compound_structures VALUES (?,?,?,?,?)",
            (i, "MOLBLOCK " * 500, f"InChI=1S/C{i}H{i}", f"KEY{i:016d}-A-N",
             "C" * (i % 7 + 1)),
        )
        conn.execute(
            "INSERT INTO compound_properties VALUES (?,?,?,?,?,?,?)",
            (i, 100.0 + i, 1.5, 2, 1, 40.0, f"C{i}H{i}"),
        )
        conn.execute(
            "INSERT INTO molecule_synonyms VALUES (?,?,?,?)",
            (i, "TRADE_NAME", i, f"synonym-of-{i}"),
        )
        conn.execute("INSERT INTO molecule_hierarchy VALUES (?,?,?)", (i, i, i))
        conn.execute(
            "INSERT INTO chembl_id_lookup VALUES (?,?,?,?,?)",
            (f"CHEMBL{i}", "COMPOUND", i, "ACTIVE", 1),
        )
        # Assay and document ids share the lookup table with compounds; only the
        # COMPOUND rows are ever queried, and only they should be copied.
        conn.execute(
            "INSERT INTO chembl_id_lookup VALUES (?,?,?,?,?)",
            (f"CHEMBL_ASSAY_{i}", "ASSAY", i, "ACTIVE", 1),
        )
        conn.execute("INSERT INTO activities VALUES (?,?,?,?)", (i, i, 1.0, "nM"))
    conn.execute(
        "INSERT INTO pesticide_classification VALUES (1,'ATRAZINE',1,'C1','FRAC','1','u')"
    )
    conn.execute("INSERT INTO pesticide_class_mapping VALUES (1, 1, 1)")
    conn.commit()
    conn.close()
    return str(path)


@pytest.fixture
def full_db(tmp_path):
    """A miniature full ChEMBL release named like the real thing."""
    return _make_full_chembl(tmp_path / "chembl_36.db")


def _tables(path):
    with sqlite3.connect(path) as conn:
        return {
            row[0]
            for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")
        }


# ── What the extract contains ─────────────────────────────────────────────────

class TestCompactContents:
    """The extract holds the eight tables, and nothing else."""

    def test_only_the_eight_used_tables_are_copied(self, full_db):
        extract = CheMBL.compact(full_db)
        tables = _tables(extract)

        assert set(CheMBL.COMPACT_TABLES) <= tables
        assert "activities" not in tables, "bioactivity data must not be copied"

    def test_molfile_column_is_dropped(self, full_db):
        extract = CheMBL.compact(full_db)
        with sqlite3.connect(extract) as conn:
            columns = {
                row[1] for row in conn.execute("PRAGMA table_info(compound_structures)")
            }
        assert "molfile" not in columns
        assert {"molregno", "canonical_smiles", "standard_inchi",
                "standard_inchi_key"} == columns

    def test_only_compound_rows_of_the_lookup_table_survive(self, full_db):
        extract = CheMBL.compact(full_db)
        with sqlite3.connect(extract) as conn:
            kinds = {
                row[0]
                for row in conn.execute("SELECT DISTINCT entity_type FROM chembl_id_lookup")
            }
        assert kinds == {"COMPOUND"}

    def test_the_extract_is_smaller_than_its_source(self, full_db):
        extract = CheMBL.compact(full_db)
        assert os.path.getsize(extract) < os.path.getsize(full_db)

    def test_keep_inchi_false_drops_the_column_and_its_index(self, full_db):
        extract = CheMBL.compact(full_db, keep_inchi=False)
        with sqlite3.connect(extract) as conn:
            columns = {
                row[1] for row in conn.execute("PRAGMA table_info(compound_structures)")
            }
            indexes = {
                row[0]
                for row in conn.execute(
                    "SELECT name FROM sqlite_master WHERE type='index'"
                )
            }
        assert "standard_inchi" not in columns
        assert "ix_cs_inchi" not in indexes
        assert "ix_cs_inchikey" in indexes, "the InChIKey index must survive"

    def test_expected_indexes_are_built(self, full_db):
        extract = CheMBL.compact(full_db)
        with sqlite3.connect(extract) as conn:
            indexes = {
                row[0]
                for row in conn.execute(
                    "SELECT name FROM sqlite_master WHERE type='index' "
                    "AND name IS NOT NULL"
                )
            }
        for name, _ in CheMBL.COMPACT_INDEXES:
            assert name in indexes, f"missing index {name}"


# ── Provenance ────────────────────────────────────────────────────────────────

class TestCompactProvenance:
    """An extract can always say where it came from."""

    def test_provenance_records_the_source_and_the_builder(self, full_db):
        extract = CheMBL.compact(full_db)
        provenance = CheMBL.read_provenance(extract)

        assert provenance["release"] == "36"
        assert provenance["source_database"] == "chembl_36.db"
        assert provenance["source_bytes"] == str(os.path.getsize(full_db))
        assert provenance["schema_version"] == str(CheMBL.COMPACT_SCHEMA_VERSION)
        assert provenance["keep_inchi"] == "1"
        assert provenance["built_at"].endswith("+00:00")

    def test_provenance_records_row_counts_per_table(self, full_db):
        extract = CheMBL.compact(full_db)
        provenance = CheMBL.read_provenance(extract)

        with sqlite3.connect(extract) as conn:
            for table in CheMBL.COMPACT_TABLES:
                actual = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
                assert provenance[f"rows.{table}"] == str(actual)

    def test_a_full_release_has_no_provenance(self, full_db):
        assert CheMBL.read_provenance(full_db) is None
        assert CheMBL.is_compact_database(full_db) is False

    def test_an_extract_is_recognised_as_one(self, full_db):
        assert CheMBL.is_compact_database(CheMBL.compact(full_db)) is True

    def test_missing_file_is_not_an_extract(self, tmp_path):
        assert CheMBL.is_compact_database(str(tmp_path / "nope.db")) is False


# ── Source selection, naming and safety ───────────────────────────────────────

class TestCompactSourceAndSafety:
    """compact picks the right source and never destroys anything by accident."""

    def test_default_destination_sits_beside_the_source(self, full_db):
        assert CheMBL.compact(full_db) == str(
            os.path.join(os.path.dirname(full_db), "chembl_36_provesid.db")
        )

    def test_newest_release_is_compacted_when_no_source_is_given(self, tmp_path):
        _make_full_chembl(tmp_path / "chembl_35.db")
        _make_full_chembl(tmp_path / "chembl_37.db")

        extract = CheMBL.compact(data_dir=str(tmp_path))

        assert os.path.basename(extract) == "chembl_37_provesid.db"

    def test_an_existing_extract_is_never_a_compaction_source(self, tmp_path):
        """With only an extract on disk, there is nothing left to compact."""
        full = _make_full_chembl(tmp_path / "chembl_36.db")
        CheMBL.compact(full)
        os.remove(full)

        with pytest.raises(FileNotFoundError, match="No full ChEMBL database"):
            CheMBL.compact(data_dir=str(tmp_path))

    def test_compacting_an_extract_is_refused(self, full_db):
        extract = CheMBL.compact(full_db)
        with pytest.raises(ChEMBLError, match="already a PROVESID extract"):
            CheMBL.compact(extract)

    def test_existing_extract_is_not_overwritten_without_force(self, full_db):
        CheMBL.compact(full_db)
        with pytest.raises(FileExistsError, match="force=True"):
            CheMBL.compact(full_db)

    def test_force_rebuilds_over_an_existing_extract(self, full_db):
        first = CheMBL.compact(full_db)
        assert CheMBL.compact(full_db, force=True) == first

    def test_the_source_survives_by_default(self, full_db):
        CheMBL.compact(full_db)
        assert os.path.exists(full_db), "compact must not delete 30 GB unasked"

    def test_remove_source_deletes_the_full_database(self, full_db):
        extract = CheMBL.compact(full_db, remove_source=True)
        assert not os.path.exists(full_db)
        assert os.path.exists(extract)

    def test_a_truncated_source_is_rejected_before_anything_is_written(self, tmp_path):
        """A database missing a needed table fails with its name, not a SQL error."""
        path = str(tmp_path / "chembl_36.db")
        conn = sqlite3.connect(path)
        conn.execute("CREATE TABLE molecule_dictionary (molregno INTEGER)")
        conn.commit()
        conn.close()

        with pytest.raises(ChEMBLError, match="compound_structures"):
            CheMBL.compact(path)
        assert not glob.glob(str(tmp_path / "*.tmp"))

    def test_a_missing_column_is_reported_by_name(self, tmp_path):
        path = _make_full_chembl(tmp_path / "chembl_36.db")
        with sqlite3.connect(path) as conn:
            conn.execute("ALTER TABLE compound_structures DROP COLUMN standard_inchi")

        with pytest.raises(ChEMBLError, match="standard_inchi"):
            CheMBL.compact(path)

    def test_nonexistent_source_path_is_reported(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="No ChEMBL database at"):
            CheMBL.compact(str(tmp_path / "absent.db"))

    def test_no_temporary_file_is_left_behind(self, full_db):
        CheMBL.compact(full_db)
        assert not glob.glob(os.path.join(os.path.dirname(full_db), "*.tmp"))


# ── Opening an extract ────────────────────────────────────────────────────────

class TestOpeningAnExtract:
    """CheMBL prefers the extract and reports what it opened."""

    def test_extract_is_preferred_over_the_full_release_beside_it(self, tmp_path):
        full = _make_full_chembl(tmp_path / "chembl_36.db")
        extract = CheMBL.compact(full)

        chembl = CheMBL(data_dir=str(tmp_path), auto_download=False)

        assert chembl.db_path == extract
        assert chembl.is_compact is True
        assert chembl.release == 36

    def test_a_newer_full_release_beats_an_older_extract(self, tmp_path):
        """Release number first, extract-ness only as the tie-break."""
        CheMBL.compact(_make_full_chembl(tmp_path / "chembl_36.db"))
        _make_full_chembl(tmp_path / "chembl_37.db")

        chembl = CheMBL(data_dir=str(tmp_path), auto_download=False)

        assert chembl.release == 37
        assert chembl.is_compact is False

    def test_a_full_release_reports_no_provenance(self, tmp_path):
        _make_full_chembl(tmp_path / "chembl_36.db")
        chembl = CheMBL(data_dir=str(tmp_path), auto_download=False)

        assert chembl.is_compact is False
        assert chembl.provenance is None

    def test_a_stale_extract_warns_with_the_rebuild_instruction(self, tmp_path, caplog):
        extract = CheMBL.compact(_make_full_chembl(tmp_path / "chembl_36.db"))
        with sqlite3.connect(extract) as conn:
            conn.execute(
                f"UPDATE {CheMBL.PROVENANCE_TABLE} SET value = '0' "
                "WHERE key = 'schema_version'"
            )

        with caplog.at_level("WARNING"):
            CheMBL(db_path=extract, auto_download=False)

        assert "compact(force=True)" in caplog.text

    def test_queries_answer_the_same_from_the_extract(self, tmp_path):
        """Every public lookup returns the same content from either database."""
        full_path = _make_full_chembl(tmp_path / "chembl_36.db")
        extract_path = CheMBL.compact(full_path)

        full = CheMBL(db_path=full_path, auto_download=False)
        extract = CheMBL(db_path=extract_path, auto_download=False)

        assert full.chembl_id_to_molregno("CHEMBL7") == extract.chembl_id_to_molregno("CHEMBL7")
        assert full.molregno_to_chembl_id(7) == extract.molregno_to_chembl_id(7)
        assert full.get_compound(7) == extract.get_compound(7)
        assert full.get_properties(7) == extract.get_properties(7)
        assert full.get_molecule_dictionary(7) == extract.get_molecule_dictionary(7)
        assert full.get_molecule_hierarchy(7) == extract.get_molecule_hierarchy(7)
        assert full.search_by_inchikey("KEY0000000000007-A-N") == \
            extract.search_by_inchikey("KEY0000000000007-A-N")
        assert full.search_by_inchi("InChI=1S/C7H7") == extract.search_by_inchi("InChI=1S/C7H7")
        assert full.get_pesticide_classifications(1) == extract.get_pesticide_classifications(1)
        # search_by_name has no ORDER BY, so compare content rather than order.
        assert sorted(c["molregno"] for c in full.search_by_name("COMPOUND 7", exact=True)) == \
            sorted(c["molregno"] for c in extract.search_by_name("COMPOUND 7", exact=True))

    def test_get_compound_no_longer_returns_a_molfile(self, tmp_path):
        """Dropped from the query, so a full release stops returning it too."""
        full_path = _make_full_chembl(tmp_path / "chembl_36.db")
        compound = CheMBL(db_path=full_path, auto_download=False).get_compound(3)

        assert compound is not None
        assert "molfile" not in compound
        assert compound["canonical_smiles"] is not None


# ── Acquisition: source= ──────────────────────────────────────────────────────

def _make_archive(tmp_path, release=36, compounds=25):
    """Package a miniature release the way EBI does, nested two deep.

    ChEMBL's archive holds ``chembl_NN/chembl_NN_sqlite/chembl_NN.db``, and the
    depth has changed between releases, so the extraction walks the tree rather
    than assuming a layout. The tests use the real shape.
    """
    staging = tmp_path / "staging" / f"chembl_{release}" / f"chembl_{release}_sqlite"
    staging.mkdir(parents=True)
    _make_full_chembl(staging / f"chembl_{release}.db", n_compounds=compounds)

    archive = tmp_path / f"chembl_{release}_sqlite.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(staging / f"chembl_{release}.db",
                arcname=f"chembl_{release}/chembl_{release}_sqlite/"
                        f"chembl_{release}.db")
    return str(archive)


@pytest.fixture
def served_release(tmp_path, monkeypatch):
    """Stand in for the 5.8 GB transfer with a local copy of a tiny archive.

    Returns the URL the fake download answers to; ``download_file`` is replaced
    for the duration, so nothing reaches the network and the rest of the
    download path --- extraction, validation, compaction, cleanup --- runs for
    real.
    """
    archive = _make_archive(tmp_path)
    url = "http://example.invalid/chembl_36_sqlite.tar.gz"
    calls = []

    def fake_download_file(requested_url, dest, **kwargs):
        calls.append((requested_url, dest))
        shutil.copyfile(archive, dest)
        return dest

    monkeypatch.setattr("provesid.chembl.download_file", fake_download_file)
    return SimpleNamespace(url=url, calls=calls)


def _install(tmp_path, served_release, **kwargs):
    """Construct a CheMBL that downloads the miniature release."""
    return CheMBL(data_dir=str(tmp_path / "data"), db_url=served_release.url,
                  **kwargs)


class TestSourceValidation:
    """``source`` is checked before anything is fetched."""

    def test_default_is_the_compacting_route(self):
        assert inspect.signature(CheMBL.__init__).parameters["source"].default == "sqlite"

    def test_an_unknown_route_is_refused(self, tmp_path):
        with pytest.raises(ValueError) as exc_info:
            CheMBL(data_dir=str(tmp_path), source="ftp")
        assert "'sqlite'" in str(exc_info.value)

    def test_the_check_happens_before_any_download(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "provesid.chembl.download_file",
            lambda *a, **k: pytest.fail("a bad source must not reach the network"),
        )
        with pytest.raises(ValueError):
            CheMBL(data_dir=str(tmp_path), source="nonsense")


class TestDownloadBuildsTheExtract:
    """``source="sqlite"``: what a first install leaves on disk."""

    def test_only_the_extract_survives(self, tmp_path, served_release):
        chembl = _install(tmp_path, served_release)
        directory = tmp_path / "data"

        assert chembl.db_path == str(directory / "chembl_36_provesid.db")
        assert chembl.is_compact
        assert chembl.release == 36
        assert not os.path.exists(directory / "chembl_36.db"), \
            "the full release must be deleted once the extract verifies"
        assert glob.glob(str(directory / "*.tar.gz")) == []
        assert glob.glob(str(directory / "*.incoming")) == []
        assert glob.glob(str(directory / "*.tmp")) == []

    def test_the_extract_answers_and_knows_where_it_came_from(self, tmp_path,
                                                              served_release):
        chembl = _install(tmp_path, served_release)

        assert chembl.search_by_chembl_id("CHEMBL7")["pref_name"] == "COMPOUND 7"
        assert chembl.provenance["release"] == "36"
        assert chembl.provenance["source_database"] == "chembl_36.db"

    def test_full_keeps_the_release_and_builds_no_extract(self, tmp_path,
                                                          served_release):
        chembl = _install(tmp_path, served_release, source="full")
        directory = tmp_path / "data"

        assert chembl.db_path == str(directory / "chembl_36.db")
        assert not chembl.is_compact
        assert not os.path.exists(directory / "chembl_36_provesid.db")
        assert glob.glob(str(directory / "*.tar.gz")) == []

    def test_both_routes_transfer_the_same_archive(self, tmp_path, served_release):
        """The choice is what is kept, not what is fetched."""
        _install(tmp_path, served_release, source="full")
        assert served_release.calls[0][0] == served_release.url
        assert served_release.calls[0][1].endswith("chembl_36.db.tar.gz")

    def test_a_database_already_on_disk_is_left_alone(self, tmp_path, monkeypatch):
        """``source`` describes a download; it must not compact what is here."""
        directory = tmp_path / "data"
        directory.mkdir()
        _make_full_chembl(directory / "chembl_36.db")
        monkeypatch.setattr(
            "provesid.chembl.download_file",
            lambda *a, **k: pytest.fail("an installed release must not be refetched"),
        )

        chembl = CheMBL(data_dir=str(directory), source="sqlite")

        assert chembl.db_path == str(directory / "chembl_36.db")
        assert not chembl.is_compact
        assert not os.path.exists(directory / "chembl_36_provesid.db")

    def test_redownload_rebuilds_over_an_existing_extract(self, tmp_path,
                                                          served_release):
        """An extract from a previous attempt must not block the next one."""
        chembl = _install(tmp_path, served_release)

        again = _install(tmp_path, served_release, redownload=True)

        assert again.db_path == chembl.db_path
        assert again.is_compact
        assert len(served_release.calls) == 2, "the archive was fetched again"
        # The whole route ran a second time, deletion included.
        assert not os.path.exists(tmp_path / "data" / "chembl_36.db")

    def test_a_failed_compaction_keeps_the_release_and_says_so(
        self, tmp_path, served_release, monkeypatch, caplog
    ):
        """A download is worth too much to throw away over the step after it."""
        def explode(*args, **kwargs):
            raise ChEMBLError("no space left on device")

        monkeypatch.setattr(CheMBL, "compact", staticmethod(explode))

        with caplog.at_level(logging.WARNING, logger="provesid.chembl"):
            chembl = _install(tmp_path, served_release)

        directory = tmp_path / "data"
        assert chembl.db_path == str(directory / "chembl_36.db")
        assert os.path.exists(chembl.db_path)
        assert chembl.search_by_chembl_id("CHEMBL3") is not None
        assert "CheMBL.compact(remove_source=True)" in caplog.text
        assert "no space left on device" in caplog.text


class TestCompactionHint:
    """Opening a large full release mentions the option, once."""

    def test_a_big_full_release_suggests_compacting(self, tmp_path, monkeypatch, caplog):
        _make_full_chembl(tmp_path / "chembl_36.db")
        monkeypatch.setattr(CheMBL, "_COMPACT_HINT_BYTES", 0)

        with caplog.at_level(logging.INFO, logger="provesid.chembl"):
            CheMBL(data_dir=str(tmp_path), auto_download=False)

        assert "CheMBL.compact(remove_source=True)" in caplog.text

    def test_a_small_database_says_nothing(self, tmp_path, caplog):
        """Below the threshold there is nothing to reclaim, so no advice."""
        _make_full_chembl(tmp_path / "chembl_36.db")

        with caplog.at_level(logging.INFO, logger="provesid.chembl"):
            CheMBL(data_dir=str(tmp_path), auto_download=False)

        assert "CheMBL.compact" not in caplog.text

    def test_an_extract_is_never_told_to_compact(self, tmp_path, caplog):
        CheMBL.compact(_make_full_chembl(tmp_path / "chembl_36.db"))
        os.remove(tmp_path / "chembl_36.db")

        with caplog.at_level(logging.INFO, logger="provesid.chembl"):
            CheMBL(data_dir=str(tmp_path), auto_download=False)

        assert "CheMBL.compact" not in caplog.text


# ── Against a real release ────────────────────────────────────────────────────

@pytest.mark.slow
class TestAgainstARealRelease:
    """Equivalence on real data, when a full ChEMBL database happens to be here."""

    @staticmethod
    def _find_full_release():
        directory = os.environ.get("PROVESID_DATA_DIR", "")
        for path in sorted(glob.glob(os.path.join(directory, "chembl_*.db"))):
            if not CheMBL.is_compact_database(path):
                return path
        return None

    def test_extract_answers_identically_on_real_data(self, tmp_path):
        source = self._find_full_release()
        if source is None:
            pytest.skip("no full ChEMBL database available")

        extract_path = CheMBL.compact(
            source, dest_path=str(tmp_path / "chembl_provesid.db")
        )
        full = CheMBL(db_path=source, auto_download=False)
        extract = CheMBL(db_path=extract_path, auto_download=False)

        molregnos = [
            row[0]
            for row in full.cursor.execute(
                "SELECT molregno FROM compound_structures ORDER BY random() LIMIT 200"
            )
        ]
        for molregno in molregnos:
            assert full.get_compound(molregno) == extract.get_compound(molregno)

        assert os.path.getsize(extract_path) < os.path.getsize(source) / 5


class TestVerificationSample:
    """The sample that justifies deleting 30 GB must cover the whole table."""

    def test_every_row_is_compared_when_the_table_is_small(self, full_db):
        with sqlite3.connect(full_db) as conn:
            sample = CheMBL._verification_sample(conn)
        assert sorted(sample) == list(range(1, 26))

    def test_sample_is_capped_and_spread_across_the_range(self, tmp_path, monkeypatch):
        """With more rows than the cap, the draw spans the table, not one end."""
        path = _make_full_chembl(tmp_path / "chembl_36.db", n_compounds=2000)
        monkeypatch.setattr(CheMBL, "_VERIFY_SAMPLE", 200)

        with sqlite3.connect(path) as conn:
            sample = CheMBL._verification_sample(conn)

        assert len(sample) == 200
        assert len(set(sample)) == 200, "no molregno compared twice"
        assert min(sample) < 400 and max(sample) > 1600, "sample hugs one end"

    def test_an_empty_table_yields_an_empty_sample(self, tmp_path):
        path = _make_full_chembl(tmp_path / "chembl_36.db", n_compounds=0)
        with sqlite3.connect(path) as conn:
            assert CheMBL._verification_sample(conn) == []

    def test_verification_catches_a_corrupted_extract(self, full_db, monkeypatch):
        """A silently altered row must stop remove_source from firing."""
        real_build = CheMBL._build_compact.__func__

        def build_then_corrupt(cls, source_path, dest_path, tables, keep_inchi):
            counts = real_build(cls, source_path, dest_path, tables, keep_inchi)
            with sqlite3.connect(dest_path) as conn:
                conn.execute(
                    "UPDATE compound_structures SET canonical_smiles = 'WRONG' "
                    "WHERE molregno = 5"
                )
            return counts

        monkeypatch.setattr(CheMBL, "_build_compact", classmethod(build_then_corrupt))

        with pytest.raises(ChEMBLError, match="disagrees with its source"):
            CheMBL.compact(full_db, remove_source=True)

        assert os.path.exists(full_db), "a failed verification must not delete the source"
        assert not glob.glob(os.path.join(os.path.dirname(full_db), "*provesid*"))

    def test_verification_catches_a_row_count_mismatch(self, full_db, monkeypatch):
        real_build = CheMBL._build_compact.__func__

        def build_then_delete_rows(cls, source_path, dest_path, tables, keep_inchi):
            counts = real_build(cls, source_path, dest_path, tables, keep_inchi)
            with sqlite3.connect(dest_path) as conn:
                conn.execute("DELETE FROM molecule_synonyms WHERE molregno > 10")
            return counts

        monkeypatch.setattr(CheMBL, "_build_compact", classmethod(build_then_delete_rows))

        with pytest.raises(ChEMBLError, match="Row count mismatch in molecule_synonyms"):
            CheMBL.compact(full_db, remove_source=True)
        assert os.path.exists(full_db)
