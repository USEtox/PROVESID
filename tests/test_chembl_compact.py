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
import os
import sqlite3

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
