"""
Tests for ``CheMBL(source="mysql")`` and ``CheMBL.build_from_mysql_dump``.

The MySQL route has to build *the same extract* as the SQLite route (§9.7 of
the 2026-09-20 plan: "build both once, compare row counts and a content hash
per table").  So every test here starts from the miniature ChEMBL of
``test_chembl_compact``, writes it out as ``mysqldump`` would, and builds the
extract both ways: ``compact`` from the SQLite file, ``build_from_mysql_dump``
from the dump.  ``CheMBL.extract_digest`` then has to agree table by table.

A slow test at the bottom does the same with a real release, when both a
``chembl_NN_mysql.tar.gz`` and a ``chembl_NN_provesid.db`` of the same release
are in ``PROVESID_DATA_DIR``.
"""

import glob
import io
import logging
import os
import re
import shutil
import sqlite3
import tarfile
from types import SimpleNamespace

import pytest

from provesid import CheMBL, ChEMBLError
from test_chembl_compact import _make_full_chembl


# ── Writing a miniature ChEMBL out as mysqldump would ─────────────────────────

# The MySQL types ChEMBL's dump uses for the SQLite types of the miniature.
_MYSQL_TYPES = {"INTEGER": "bigint", "TEXT": "varchar(4000)",
                "REAL": "double", "NUMERIC": "decimal(9,2)"}

# mysql_real_escape_string's set: backslash first, so it is not doubled twice.
_ESCAPES = [("\\", "\\\\"), ("'", "\\'"), ('"', '\\"'), ("\0", "\\0"),
            ("\n", "\\n"), ("\r", "\\r"), ("\x1a", "\\Z")]


def _mysql_literal(value, mysql_type):
    if value is None:
        return "NULL"
    if isinstance(value, str):
        for raw, escaped in _ESCAPES:
            value = value.replace(raw, escaped)
        return f"'{value}'"
    if mysql_type.startswith("decimal"):
        return f"{float(value):.2f}"     # decimal(9,2) always prints two places
    return repr(value)


def _write_dump(db_path, dump_path, *, rows_per_insert=7, complete_insert=False,
                skip_tables=(), empty_tables=(), drop_columns=None):
    """Write every table of ``db_path`` in mysqldump's layout."""
    drop_columns = drop_columns or {}
    conn = sqlite3.connect(db_path)
    out = ["-- MySQL dump 10.13  Distrib 8.0.36, for Linux (x86_64)",
           "/*!40101 SET NAMES utf8mb4 */;",
           "/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;"]
    tables = [row[0] for row in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")]
    for table in tables:
        if table in skip_tables:
            continue
        info = [(row[1], _MYSQL_TYPES[row[2]])
                for row in conn.execute(f"PRAGMA table_info({table})")
                if row[1] not in drop_columns.get(table, ())]
        names = [name for name, _ in info]
        out += [f"DROP TABLE IF EXISTS `{table}`;",
                "/*!40101 SET @saved_cs_client     = @@character_set_client */;",
                f"CREATE TABLE `{table}` ("]
        out += [f"  `{name}` {kind} DEFAULT NULL," for name, kind in info]
        out += [f"  PRIMARY KEY (`{names[0]}`),",
                f"  KEY `idx_{table}` (`{names[0]}`)",
                ") ENGINE=InnoDB DEFAULT CHARSET=utf8mb3;",
                f"LOCK TABLES `{table}` WRITE;",
                f"/*!40000 ALTER TABLE `{table}` DISABLE KEYS */;"]
        if table not in empty_tables:
            rows = conn.execute(f"SELECT {', '.join(names)} FROM {table}").fetchall()
            order = list(range(len(names)))
            if complete_insert:
                order.reverse()     # a column list in a different order
            head = f"INSERT INTO `{table}` "
            if complete_insert:
                head += "(" + ", ".join(f"`{names[i]}`" for i in order) + ") "
            for start in range(0, len(rows), rows_per_insert):
                chunk = rows[start:start + rows_per_insert]
                values = ",".join(
                    "(" + ",".join(_mysql_literal(row[i], info[i][1]) for i in order) + ")"
                    for row in chunk
                )
                out.append(f"{head}VALUES {values};")
        out += [f"/*!40000 ALTER TABLE `{table}` ENABLE KEYS */;", "UNLOCK TABLES;"]
    out.append("-- Dump completed on 2026-09-21 12:00:00")
    conn.close()
    with open(dump_path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(out) + "\n")


def _add_awkward_values(db_path):
    """Values that exercise every escape, NULLs, and the NUMERIC 4.0 → 4 rule."""
    conn = sqlite3.connect(db_path)
    awkward = ["it's", 'say "hi"', "back\\slash", "two\nlines", "cr\rlf",
               "nul\0byte", "ctrl\x1aZ", "β-carotène", "trailing\\"]
    for offset, text in enumerate(awkward, start=1):
        conn.execute("INSERT INTO molecule_synonyms VALUES (?,?,?,?)",
                     (offset, "OTHER", 1000 + offset, text))
    conn.execute("INSERT INTO molecule_synonyms VALUES (1, 'OTHER', 2000, NULL)")
    conn.execute("UPDATE molecule_dictionary SET pref_name = NULL WHERE molregno = 2")
    conn.execute("UPDATE molecule_dictionary SET max_phase = 0.5 WHERE molregno = 3")
    conn.execute("UPDATE molecule_dictionary SET max_phase = NULL WHERE molregno = 4")
    conn.commit()
    conn.close()


def _pack(dump_path, archive_path, release=36):
    """Package a dump the way EBI does: chembl_NN/chembl_NN_mysql/*.dmp."""
    with tarfile.open(archive_path, "w:gz") as tar:
        tar.add(dump_path, arcname=f"chembl_{release}/chembl_{release}_mysql/"
                                   f"chembl_{release}_mysql.dmp")
        info = tarfile.TarInfo(f"chembl_{release}/chembl_{release}_mysql/INSTALL_mysql")
        text = b"mysql -u user -p chembl_36 < chembl_36_mysql.dmp\n"
        info.size = len(text)
        tar.addfile(info, io.BytesIO(text))
    return str(archive_path)


@pytest.fixture
def release(tmp_path):
    """One miniature release, as a SQLite file and as a MySQL dump archive."""
    work = tmp_path / "work"
    work.mkdir()
    db = str(work / "chembl_36.db")
    _make_full_chembl(db)
    _add_awkward_values(db)
    _write_dump(db, work / "chembl_36_mysql.dmp")
    archive = _pack(work / "chembl_36_mysql.dmp", work / "chembl_36_mysql.tar.gz")
    return SimpleNamespace(db=db, archive=archive, work=work)


def _variant(tmp_path, release, **dump_options):
    """The same release dumped with different options."""
    dump = tmp_path / "variant.dmp"
    _write_dump(release.db, dump, **dump_options)
    return _pack(dump, tmp_path / "chembl_36_mysql.tar.gz")


def _schema(path):
    conn = sqlite3.connect(path)
    try:
        return {
            name: sorted((row[1], row[2])
                         for row in conn.execute(f"PRAGMA table_info({name})"))
            for name in CheMBL.COMPACT_TABLES
        }
    finally:
        conn.close()


def _index_names(path):
    conn = sqlite3.connect(path)
    try:
        return {row[0] for row in conn.execute(
            "SELECT name FROM sqlite_master WHERE type='index'")}
    finally:
        conn.close()


# ── Equivalence with the SQLite route ─────────────────────────────────────────

class TestSameExtractAsTheSqliteRoute:
    """§9.7: row counts and a content hash per table, both routes."""

    def test_every_table_digests_identically(self, tmp_path, release):
        from_sqlite = CheMBL.compact(release.db, str(tmp_path / "a.db"))
        from_mysql = CheMBL.build_from_mysql_dump(release.archive,
                                                  str(tmp_path / "b.db"))

        expected = CheMBL.extract_digest(from_sqlite)
        assert set(expected) == set(CheMBL.COMPACT_TABLES)
        assert CheMBL.extract_digest(from_mysql) == expected

    def test_keep_inchi_false_matches_too(self, tmp_path, release):
        from_sqlite = CheMBL.compact(release.db, str(tmp_path / "a.db"),
                                     keep_inchi=False)
        from_mysql = CheMBL.build_from_mysql_dump(
            release.archive, str(tmp_path / "b.db"), keep_inchi=False)

        assert CheMBL.extract_digest(from_mysql) == CheMBL.extract_digest(from_sqlite)
        assert "ix_cs_inchi" not in _index_names(from_mysql)

    def test_column_types_and_indexes_match(self, tmp_path, release):
        from_sqlite = CheMBL.compact(release.db, str(tmp_path / "a.db"))
        from_mysql = CheMBL.build_from_mysql_dump(release.archive,
                                                  str(tmp_path / "b.db"))

        assert _schema(from_mysql) == _schema(from_sqlite)
        assert _index_names(from_mysql) == _index_names(from_sqlite)

    def test_numeric_values_are_stored_as_the_release_stores_them(self, tmp_path,
                                                                 release):
        """``decimal`` 4.00 must land as the integer 4, and 0.50 as a real."""
        path = CheMBL.build_from_mysql_dump(release.archive, str(tmp_path / "b.db"))
        conn = sqlite3.connect(path)
        rows = dict(conn.execute(
            "SELECT molregno, typeof(max_phase) FROM molecule_dictionary "
            "WHERE molregno IN (1, 3, 4)"))
        conn.close()
        assert rows == {1: "integer", 3: "real", 4: "null"}

    def test_queries_answer_the_same(self, tmp_path, release):
        from_sqlite = CheMBL(db_path=CheMBL.compact(release.db, str(tmp_path / "a.db")),
                             auto_download=False)
        from_mysql = CheMBL(db_path=CheMBL.build_from_mysql_dump(
            release.archive, str(tmp_path / "b.db")), auto_download=False)

        for molregno in range(1, 26):
            assert from_mysql.get_compound(molregno) == from_sqlite.get_compound(molregno)
            assert (from_mysql.get_molecule_dictionary(molregno)
                    == from_sqlite.get_molecule_dictionary(molregno))
        for name in ["it's", "two\nlines", "β-carotène", "trailing\\", "COMPOUND 7"]:
            assert (from_mysql.search_by_name(name, exact=True)
                    == from_sqlite.search_by_name(name, exact=True))
            assert from_mysql.search_by_name(name, exact=True), name
        assert (from_mysql.search_pesticide_by_name("ATRAZINE")
                == from_sqlite.search_pesticide_by_name("ATRAZINE"))

    def test_a_column_list_in_another_order_is_honoured(self, tmp_path, release):
        archive = _variant(tmp_path, release, complete_insert=True)
        from_sqlite = CheMBL.compact(release.db, str(tmp_path / "a.db"))
        from_mysql = CheMBL.build_from_mysql_dump(archive, str(tmp_path / "b.db"))
        assert CheMBL.extract_digest(from_mysql) == CheMBL.extract_digest(from_sqlite)

    def test_the_unused_tables_are_not_copied(self, tmp_path, release):
        path = CheMBL.build_from_mysql_dump(release.archive, str(tmp_path / "b.db"))
        conn = sqlite3.connect(path)
        tables = {row[0] for row in conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' "
            "AND name NOT LIKE 'sqlite_%'")}
        conn.close()
        assert tables == set(CheMBL.COMPACT_TABLES) | {CheMBL.PROVENANCE_TABLE}


class TestExtractDigest:
    """The comparison has to be able to fail."""

    def test_a_changed_value_changes_the_digest(self, tmp_path, release):
        path = CheMBL.compact(release.db, str(tmp_path / "a.db"))
        before = CheMBL.extract_digest(path)
        conn = sqlite3.connect(path)
        conn.execute("UPDATE compound_structures SET canonical_smiles = 'N' "
                     "WHERE molregno = 5")
        conn.commit()
        conn.close()
        after = CheMBL.extract_digest(path)
        assert after["compound_structures"] != before["compound_structures"]
        assert after["molecule_dictionary"] == before["molecule_dictionary"]

    def test_an_integer_and_an_equal_real_differ(self, tmp_path):
        """The digest compares storage types, not just values: 180 != 180.0."""
        digests = []
        for name, value in (("int.db", 1), ("real.db", 1.0)):
            conn = sqlite3.connect(tmp_path / name)
            # No declared types, so no affinity: each value is stored as given.
            conn.execute("CREATE TABLE pesticide_class_mapping "
                         "(mol_pest_id, pest_class_id, molregno)")
            conn.execute("INSERT INTO pesticide_class_mapping VALUES (1, 1, ?)", (value,))
            conn.commit()
            conn.close()
            digests.append(CheMBL.extract_digest(str(tmp_path / name)))
        assert 1 == 1.0
        assert digests[0] != digests[1]

    def test_row_order_does_not_matter(self, tmp_path, release):
        path = CheMBL.compact(release.db, str(tmp_path / "a.db"))
        before = CheMBL.extract_digest(path)
        conn = sqlite3.connect(path)
        conn.executescript(
            "CREATE TABLE s AS SELECT * FROM molecule_synonyms ORDER BY molsyn_id DESC;"
            "DELETE FROM molecule_synonyms;"
            "INSERT INTO molecule_synonyms SELECT * FROM s; DROP TABLE s;"
        )
        conn.close()
        assert CheMBL.extract_digest(path) == before

    def test_a_full_release_digests_like_its_extract(self, tmp_path, release):
        """molfile and the non-compound lookup rows are outside the digest."""
        path = CheMBL.compact(release.db, str(tmp_path / "a.db"))
        assert CheMBL.extract_digest(release.db) == CheMBL.extract_digest(path)


# ── build_from_mysql_dump: paths, provenance, failure ─────────────────────────

class TestBuild:

    def test_default_destination_sits_beside_the_dump(self, release):
        path = CheMBL.build_from_mysql_dump(release.archive)
        assert path == str(release.work / "chembl_36_provesid.db")

    def test_provenance_says_it_came_from_the_dump(self, release):
        path = CheMBL.build_from_mysql_dump(release.archive)
        provenance = CheMBL.read_provenance(path)
        assert provenance["source_format"] == "mysql"
        assert provenance["source_database"] == "chembl_36_mysql.tar.gz"
        assert provenance["release"] == "36"
        assert provenance["schema_version"] == str(CheMBL.COMPACT_SCHEMA_VERSION)
        assert provenance["rows.chembl_id_lookup"] == "25"   # compounds only

    def test_the_sqlite_route_records_its_format_too(self, tmp_path, release):
        path = CheMBL.compact(release.db, str(tmp_path / "a.db"))
        assert CheMBL.read_provenance(path)["source_format"] == "sqlite"

    def test_the_dump_is_kept_by_default_and_removed_on_request(self, tmp_path,
                                                               release):
        CheMBL.build_from_mysql_dump(release.archive, str(tmp_path / "a.db"))
        assert os.path.exists(release.archive)
        CheMBL.build_from_mysql_dump(release.archive, str(tmp_path / "b.db"),
                                     remove_source=True)
        assert not os.path.exists(release.archive)

    def test_an_existing_extract_needs_force(self, tmp_path, release):
        dest = str(tmp_path / "a.db")
        CheMBL.build_from_mysql_dump(release.archive, dest)
        with pytest.raises(FileExistsError):
            CheMBL.build_from_mysql_dump(release.archive, dest)
        assert CheMBL.build_from_mysql_dump(release.archive, dest, force=True) == dest

    def test_a_missing_dump_is_reported(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            CheMBL.build_from_mysql_dump(str(tmp_path / "nope.tar.gz"))

    @pytest.mark.parametrize("options, message", [
        ({"skip_tables": ("molecule_hierarchy",)}, "no CREATE TABLE for molecule_hierarchy"),
        ({"empty_tables": ("compound_properties",)}, "no rows for compound_properties"),
        ({"drop_columns": {"compound_structures": ("canonical_smiles",)}},
         "missing column(s) canonical_smiles"),
    ])
    def test_an_incomplete_dump_is_refused_by_name(self, tmp_path, release,
                                                   options, message):
        archive = _variant(tmp_path, release, **options)
        dest = tmp_path / "out.db"
        with pytest.raises(ChEMBLError, match=re.escape(message)):
            CheMBL.build_from_mysql_dump(archive, str(dest))
        assert not dest.exists()
        assert glob.glob(str(tmp_path / "*.tmp")) == []

    def test_an_archive_without_a_dump_is_refused(self, tmp_path, release):
        archive = tmp_path / "chembl_36_mysql.tar.gz"
        with tarfile.open(archive, "w:gz") as tar:
            tar.add(release.db, arcname="chembl_36/chembl_36.db")
        with pytest.raises(ChEMBLError, match="no .dmp or .sql"):
            CheMBL.build_from_mysql_dump(str(archive), str(tmp_path / "out.db"))

    def test_a_corrupt_archive_is_refused(self, tmp_path, release):
        data = open(release.archive, "rb").read()
        broken = tmp_path / "chembl_36_mysql.tar.gz"
        broken.write_bytes(data[: len(data) // 2])
        with pytest.raises(ChEMBLError, match="Could not read the MySQL dump"):
            CheMBL.build_from_mysql_dump(str(broken), str(tmp_path / "out.db"))
        assert not (tmp_path / "out.db").exists()

    def test_an_unreadable_value_is_refused(self, tmp_path, release):
        dump = tmp_path / "bad.dmp"
        _write_dump(release.db, dump)
        text = dump.read_text(encoding="utf-8")
        text = text.replace("'COMPOUND 7'", "0x434F4D504F554E442037", 1)
        dump.write_text(text, encoding="utf-8")
        archive = _pack(dump, tmp_path / "chembl_36_mysql.tar.gz")
        with pytest.raises(ChEMBLError, match="molecule_dictionary"):
            CheMBL.build_from_mysql_dump(archive, str(tmp_path / "out.db"))


# ── CheMBL(source="mysql") ────────────────────────────────────────────────────

@pytest.fixture
def served_dump(release, monkeypatch):
    """Stand in for the 2.1 GB transfer; everything after it runs for real."""
    calls = []

    def fake_download_file(requested_url, dest, **kwargs):
        calls.append((requested_url, dest))
        shutil.copyfile(release.archive, dest)
        return dest

    monkeypatch.setattr("provesid.chembl.download_file", fake_download_file)
    return SimpleNamespace(
        url="http://example.invalid/chembl_36_mysql.tar.gz", calls=calls)


class TestMysqlRoute:

    def test_mysql_is_a_route(self):
        assert "mysql" in CheMBL.SOURCES

    def test_only_the_extract_survives(self, tmp_path, served_dump):
        directory = tmp_path / "data"
        chembl = CheMBL(data_dir=str(directory), db_url=served_dump.url,
                        source="mysql")

        assert chembl.db_path == str(directory / "chembl_36_provesid.db")
        assert chembl.is_compact
        assert chembl.release == 36
        assert chembl.provenance["source_format"] == "mysql"
        assert chembl.search_by_chembl_id("CHEMBL7")["pref_name"] == "COMPOUND 7"
        assert sorted(os.listdir(directory)) == ["chembl_36_provesid.db"]

    def test_the_dump_is_what_is_downloaded(self, tmp_path, served_dump):
        CheMBL(data_dir=str(tmp_path / "data"), db_url=served_dump.url,
               source="mysql")
        assert len(served_dump.calls) == 1
        url, dest = served_dump.calls[0]
        assert url == served_dump.url
        assert dest.endswith("chembl_36_mysql.tar.gz")

    def test_the_latest_dump_is_resolved_from_the_listing(self, tmp_path,
                                                          served_dump, monkeypatch):
        listing = SimpleNamespace(
            text='<a href="chembl_38_sqlite.tar.gz">x</a>'
                 '<a href="chembl_38_mysql.tar.gz">x</a>',
            raise_for_status=lambda: None,
        )
        monkeypatch.setattr("provesid.chembl.requests.get", lambda *a, **k: listing)

        CheMBL(data_dir=str(tmp_path / "data"), source="mysql")

        assert served_dump.calls[0][0].endswith("/chembl_38_mysql.tar.gz")

    def test_the_sqlite_route_still_resolves_the_sqlite_archive(self, monkeypatch):
        listing = SimpleNamespace(
            text='<a href="chembl_38_mysql.tar.gz">x</a>'
                 '<a href="chembl_38_sqlite.tar.gz">x</a>',
            raise_for_status=lambda: None,
        )
        monkeypatch.setattr("provesid.chembl.requests.get", lambda *a, **k: listing)
        assert CheMBL.resolve_latest_db_url().endswith("chembl_38_sqlite.tar.gz")
        assert (CheMBL.resolve_latest_db_url(archive="mysql")
                .endswith("chembl_38_mysql.tar.gz"))

    def test_a_database_already_on_disk_is_left_alone(self, tmp_path, monkeypatch):
        directory = tmp_path / "data"
        directory.mkdir()
        _make_full_chembl(directory / "chembl_36.db")
        monkeypatch.setattr(
            "provesid.chembl.download_file",
            lambda *a, **k: pytest.fail("an installed release must not be refetched"),
        )
        chembl = CheMBL(data_dir=str(directory), source="mysql")
        assert chembl.db_path == str(directory / "chembl_36.db")

    def test_redownload_rebuilds_over_an_existing_extract(self, tmp_path, served_dump):
        directory = tmp_path / "data"
        first = CheMBL(data_dir=str(directory), db_url=served_dump.url, source="mysql")
        again = CheMBL(data_dir=str(directory), db_url=served_dump.url,
                       source="mysql", redownload=True)
        assert again.db_path == first.db_path
        assert len(served_dump.calls) == 2
        assert sorted(os.listdir(directory)) == ["chembl_36_provesid.db"]

    def test_a_failed_build_removes_the_dump_and_names_the_other_route(
        self, tmp_path, release, monkeypatch
    ):
        broken = tmp_path / "broken.tar.gz"
        with tarfile.open(broken, "w:gz") as tar:
            tar.add(release.db, arcname="chembl_36/nothing_here.txt")
        monkeypatch.setattr(
            "provesid.chembl.download_file",
            lambda url, dest, **k: shutil.copyfile(broken, dest) or dest,
        )
        directory = tmp_path / "data"

        with pytest.raises(ChEMBLError, match="source='sqlite'"):
            CheMBL(data_dir=str(directory), source="mysql",
                   db_url="http://example.invalid/chembl_36_mysql.tar.gz")

        assert os.listdir(directory) == []

    def test_a_failed_download_is_reported_as_chembl_error(self, tmp_path, monkeypatch):
        from provesid.datasets import DownloadError

        def refuse(*args, **kwargs):
            raise DownloadError("connection refused")

        monkeypatch.setattr("provesid.chembl.download_file", refuse)
        with pytest.raises(ChEMBLError, match="Download failed"):
            CheMBL(data_dir=str(tmp_path / "data"), source="mysql",
                   db_url="http://example.invalid/chembl_36_mysql.tar.gz")


# ── Against a real release ────────────────────────────────────────────────────

@pytest.mark.slow
class TestAgainstARealRelease:
    """The §9.7 comparison on real data, when both inputs are on disk."""

    def test_the_dump_builds_the_same_extract_as_the_release(self, tmp_path, caplog):
        directory = os.environ.get("PROVESID_DATA_DIR", "")
        pairs = []
        for dump in sorted(glob.glob(os.path.join(directory, "chembl_*_mysql.tar.gz"))):
            number = re.search(r"chembl_(\d+)_mysql", dump).group(1)
            extract = os.path.join(directory, f"chembl_{number}_provesid.db")
            if CheMBL.is_compact_database(extract):
                pairs.append((dump, extract))
        if not pairs:
            pytest.skip("needs chembl_NN_mysql.tar.gz and chembl_NN_provesid.db")

        dump, reference = pairs[-1]
        with caplog.at_level(logging.INFO, logger="provesid.chembl"):
            built = CheMBL.build_from_mysql_dump(dump, str(tmp_path / "from_mysql.db"))
        assert CheMBL.extract_digest(built) == CheMBL.extract_digest(reference)
