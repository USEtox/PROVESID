"""
Install ChEMBL as 2.4 GiB instead of 27.7 GiB: ``CheMBL(source=...)``.

A ChEMBL release is published as a 5.8 GB archive that unpacks into a 27.7 GiB
SQLite database of 74 tables.  PROVESID reads eight of them, so
``CheMBL.compact()`` was written to shrink a release already on disk
(``examples/chembl/compact_demo.py``).  This script is about the other half of
that story: a machine that has never had ChEMBL should not have to install
27.7 GiB first and shrink it afterwards by hand.

``CheMBL(source="sqlite")`` -- the default -- finishes the download by building
the extract and deleting the release it came from.  ``CheMBL(source="full")``
keeps the whole release, for anyone who wants the other 66 tables.  Both
transfer exactly the same bytes: the choice is what stays on disk.

The transient cost of those two is real and is worth knowing before you start:
the archive and the unpacked release both exist on the way in, so installing
ChEMBL that way needs about 33.4 GiB free even though it leaves 2.4 GiB behind.
``datasets.plan()`` reports that as ``peak_bytes``.

``CheMBL(source="mysql")`` avoids it.  It downloads ChEMBL's 2.1 GB MySQL dump
instead and reads the extract's eight tables straight out of it, so the 27.7 GiB
release never exists and about 4.5 GiB free is enough.  It builds the same
extract, which this script checks with ``CheMBL.extract_digest``.

This script demonstrates all three routes against a miniature release it builds
and serves itself -- as a SQLite archive and as a MySQL dump -- so it needs no
network and downloads nothing large.

Run with::

    python examples/chembl/download_source_demo.py
"""

import http.server
import logging
import os
import sqlite3
import tarfile
import tempfile
import threading

from provesid import CheMBL, datasets

# The download path reports through the standard logging module, like the rest
# of PROVESID -- turn it on to watch the extract being built and verified.
logging.basicConfig(level=logging.INFO, format="%(levelname)-7s %(message)s")

RELEASE = 36


# ── A miniature ChEMBL, packaged the way EBI packages the real one ───────────

def build_archive(directory, n_compounds=200):
    """Write a small but structurally faithful ``chembl_NN_sqlite.tar.gz``.

    Only the eight tables PROVESID reads, plus an ``activities`` table standing
    in for the 66 it never opens -- which is the whole point of the extract, so
    the demo would be dishonest without one.
    """
    nested = os.path.join(directory, f"chembl_{RELEASE}", f"chembl_{RELEASE}_sqlite")
    os.makedirs(nested, exist_ok=True)
    db_path = os.path.join(nested, f"chembl_{RELEASE}.db")

    conn = sqlite3.connect(db_path)
    conn.executescript(
        """
        CREATE TABLE molecule_dictionary (
            molregno INTEGER PRIMARY KEY, pref_name TEXT, chembl_id TEXT,
            max_phase NUMERIC, therapeutic_flag INTEGER, molecule_type TEXT);
        CREATE TABLE compound_structures (
            molregno INTEGER PRIMARY KEY, molfile TEXT, standard_inchi TEXT,
            standard_inchi_key TEXT, canonical_smiles TEXT);
        CREATE TABLE compound_properties (
            molregno INTEGER PRIMARY KEY, mw_freebase REAL, full_molformula TEXT);
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
        CREATE TABLE activities (
            activity_id INTEGER PRIMARY KEY, molregno INTEGER,
            standard_value REAL, standard_units TEXT);
        """
    )
    for i in range(1, n_compounds + 1):
        conn.execute("INSERT INTO molecule_dictionary VALUES (?,?,?,?,?,?)",
                     (i, f"COMPOUND {i}", f"CHEMBL{i}", 4, 1, "Small molecule"))
        # The molfile column is a quarter of the real database and is read
        # nowhere; give it enough bulk here for the size comparison to mean
        # something.
        conn.execute("INSERT INTO compound_structures VALUES (?,?,?,?,?)",
                     (i, "MOLBLOCK " * 2000, f"InChI=1S/C{i}H{i}",
                      f"KEY{i:016d}-A-N", "C" * (i % 7 + 1)))
        conn.execute("INSERT INTO compound_properties VALUES (?,?,?)",
                     (i, 100.0 + i, f"C{i}H{i}"))
        conn.execute("INSERT INTO molecule_synonyms VALUES (?,?,?,?)",
                     (i, "TRADE_NAME", i, f"synonym-of-{i}"))
        conn.execute("INSERT INTO molecule_hierarchy VALUES (?,?,?)", (i, i, i))
        conn.execute("INSERT INTO chembl_id_lookup VALUES (?,?,?,?,?)",
                     (f"CHEMBL{i}", "COMPOUND", i, "ACTIVE", 1))
        conn.execute("INSERT INTO activities VALUES (?,?,?,?)", (i, i, 1.0, "nM"))
    conn.execute("INSERT INTO pesticide_classification "
                 "VALUES (1,'ATRAZINE',1,'C1','FRAC','1','url')")
    conn.execute("INSERT INTO pesticide_class_mapping VALUES (1, 1, 1)")
    conn.commit()
    conn.close()

    archive = os.path.join(directory, f"chembl_{RELEASE}_sqlite.tar.gz")
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(db_path, arcname=f"chembl_{RELEASE}/chembl_{RELEASE}_sqlite/"
                                 f"chembl_{RELEASE}.db")
    return archive, db_path


def build_mysql_dump(db_path, directory):
    """Write the same release as ``chembl_NN_mysql.tar.gz``, as mysqldump would.

    One ``CREATE TABLE`` per table, then its rows as multi-row ``INSERT``
    statements, one per line, with MySQL's backslash escaping.
    """
    mysql_types = {"INTEGER": "bigint", "TEXT": "varchar(4000)",
                   "REAL": "double", "NUMERIC": "decimal(9,2)"}

    def literal(value, mysql_type):
        if value is None:
            return "NULL"
        if isinstance(value, str):
            for raw, escaped in (("\\", "\\\\"), ("'", "\\'"), ("\n", "\\n")):
                value = value.replace(raw, escaped)
            return f"'{value}'"
        if mysql_type.startswith("decimal"):
            return f"{float(value):.2f}"
        return repr(value)

    dump_path = os.path.join(directory, f"chembl_{RELEASE}_mysql.dmp")
    conn = sqlite3.connect(db_path)
    with open(dump_path, "w", encoding="utf-8") as out:
        out.write("-- MySQL dump 10.13\n/*!40101 SET NAMES utf8mb4 */;\n")
        tables = [row[0] for row in conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table'")]
        for table in tables:
            columns = [(row[1], mysql_types[row[2]])
                       for row in conn.execute(f"PRAGMA table_info({table})")]
            out.write(f"CREATE TABLE `{table}` (\n")
            out.writelines(f"  `{name}` {kind} DEFAULT NULL,\n" for name, kind in columns)
            out.write(f"  PRIMARY KEY (`{columns[0][0]}`)\n) ENGINE=InnoDB;\n")
            rows = conn.execute(f"SELECT * FROM {table}").fetchall()
            for start in range(0, len(rows), 50):
                values = ",".join(
                    "(" + ",".join(literal(value, kind)
                                   for value, (_, kind) in zip(row, columns)) + ")"
                    for row in rows[start:start + 50])
                out.write(f"INSERT INTO `{table}` VALUES {values};\n")
    conn.close()

    archive = os.path.join(directory, f"chembl_{RELEASE}_mysql.tar.gz")
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(dump_path, arcname=f"chembl_{RELEASE}/chembl_{RELEASE}_mysql/"
                                   f"chembl_{RELEASE}_mysql.dmp")
    os.remove(dump_path)
    return archive


def serve(directory):
    """Serve ``directory`` over HTTP on a free port, and return the base URL."""
    handler = type(
        "QuietHandler",
        (http.server.SimpleHTTPRequestHandler,),
        {"log_message": lambda *args: None},
    )
    httpd = http.server.ThreadingHTTPServer(
        ("127.0.0.1", 0),
        lambda *args, **kw: handler(*args, directory=directory, **kw),
    )
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    return httpd, f"http://127.0.0.1:{httpd.server_port}/"


def listing(directory):
    """What is in the data directory, largest first."""
    rows = [(name, os.path.getsize(os.path.join(directory, name)))
            for name in sorted(os.listdir(directory))]
    return "\n".join(f"      {size / 1e6:9.2f} MB  {name}"
                     for name, size in sorted(rows, key=lambda row: -row[1]))


def main():
    work = tempfile.mkdtemp(prefix="provesid-chembl-source-demo-")
    served = os.path.join(work, "server")
    os.makedirs(served)

    print("\nBuilding a miniature ChEMBL release and serving it locally...")
    archive, db_path = build_archive(served)
    dump = build_mysql_dump(db_path, served)
    httpd, base_url = serve(served)
    url = base_url + os.path.basename(archive)
    dump_url = base_url + os.path.basename(dump)
    print(f"    SQLite archive: {os.path.getsize(archive) / 1e6:.2f} MB at {url}")
    print(f"    MySQL dump    : {os.path.getsize(dump) / 1e6:.2f} MB at {dump_url}")

    try:
        # ── 1. The default: download, compact, keep only the extract ─────────
        print("\n1. CheMBL(source='sqlite') -- the default:")
        compact_dir = os.path.join(work, "default")
        chembl = CheMBL(data_dir=compact_dir, db_url=url)

        print(f"\n    open database : {os.path.basename(chembl.db_path)}")
        print(f"    is an extract : {chembl.is_compact}")
        print(f"    release       : {chembl.release}")
        print(f"    built from    : {chembl.provenance['source_database']} "
              f"({int(chembl.provenance['source_bytes']) / 1e6:.2f} MB)")
        print("    what is left on disk:")
        print(listing(compact_dir))

        # The extract is not a different database, only a smaller one.
        compound = chembl.search_by_chembl_id("CHEMBL42")
        print(f"\n    CHEMBL42 -> {compound['pref_name']}, "
              f"{compound['canonical_smiles']}")

        # ── 2. The other route: keep all 74 tables ───────────────────────────
        print("\n2. CheMBL(source='full') -- keep the whole release:")
        full_dir = os.path.join(work, "full")
        full = CheMBL(data_dir=full_dir, db_url=url, source="full")

        print(f"\n    open database : {os.path.basename(full.db_path)}")
        print(f"    is an extract : {full.is_compact}")
        print("    what is left on disk:")
        print(listing(full_dir))

        extract_bytes = os.path.getsize(chembl.db_path)
        full_bytes = os.path.getsize(full.db_path)
        print(f"\n    the extract is {100 * (1 - extract_bytes / full_bytes):.1f}% "
              f"smaller, and answers the same questions:")
        print(f"    CHEMBL42 -> {full.search_by_chembl_id('CHEMBL42')['pref_name']}")

        # ── 3. The MySQL dump: the same extract, no release on the way ───────
        print("\n3. CheMBL(source='mysql') -- build the extract from the dump:")
        mysql_dir = os.path.join(work, "mysql")
        from_dump = CheMBL(data_dir=mysql_dir, db_url=dump_url, source="mysql")

        print(f"\n    open database : {os.path.basename(from_dump.db_path)}")
        print(f"    built from    : {from_dump.provenance['source_database']} "
              f"({from_dump.provenance['source_format']})")
        print("    what is left on disk:")
        print(listing(mysql_dir))

        # Row counts and a content hash per table, typeof included: the two
        # routes have to agree on every value, not just on the answers.
        by_sqlite = CheMBL.extract_digest(chembl.db_path)
        by_mysql = CheMBL.extract_digest(from_dump.db_path)
        differing = sorted(t for t in by_sqlite if by_sqlite[t] != by_mysql.get(t))
        print(f"\n    tables compared: {len(by_sqlite)}, differing: {differing or 'none'}")
        print(f"    CHEMBL42 -> {from_dump.search_by_chembl_id('CHEMBL42')['pref_name']}")

        # ── 4. A route name that does not exist ──────────────────────────────
        print("\n4. A route that is not there says so:")
        try:
            CheMBL(data_dir=full_dir, source="ftp")
        except ValueError as exc:
            print(f"    source='ftp': {exc}")

        # ── 5. What this costs on the real thing ─────────────────────────────
        print("\n5. The real release, before you start it:")
        todo = datasets.plan("chembl", data_dir=full_dir, force=True)
        row = todo.iloc[0]
        print(f"    download  : {row['download']}")
        print(f"    installed : {row['installed']}   (the extract)")
        print(f"    peak      : {datasets.human_bytes(todo.attrs['peak_bytes'])}"
              "   (archive + release, before either is deleted)")
        print(f"    note      : {row['note']}")
        print("""
    from provesid import CheMBL, datasets

    datasets.plan("chembl")      # 5.7 GiB down, 2.4 GiB installed, 33.4 GiB peak
    CheMBL()                     # ...and this installs it, as the extract
    CheMBL(source="mysql")       # the same extract from 2.1 GB, ~4.5 GiB peak
    CheMBL(source="full")        # ...or keeps all 74 tables, 27.7 GiB
""")
    finally:
        httpd.shutdown()

    print(f"(demo files are in {work})")


if __name__ == "__main__":
    main()
