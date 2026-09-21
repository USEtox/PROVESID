"""
Build the PubChem identifier database from PubChem's own FTP files.

``PubChemID`` answers CAS, name, InChIKey and formula lookups from a local
SQLite file, ``pubchem_id.db``. It used to be downloaded from Zenodo, as a copy
someone had built by hand from a CSV export --- with no record of which PubChem
release it came from, and with CAS numbers pulled out of free text by a regular
expression, 1.25% of which fail the CAS check digit.

``PubChemID(source="ftp")`` -- now the default -- builds it instead, from a
dated monthly snapshot of ``Compound/Extras/`` on PubChem's FTP site:

* the compounds are those PubChem itself maps to a CAS number, each checked;
* every source file is checked against the MD5 PubChem publishes beside it;
* the database records its release and every file's MD5 (``provenance()``);
* DSSTox, ChEBI, ChEMBL, EC and UNII cross-references come for free
  (``xrefs()``).

The real build transfers 15.4 GB and takes about 12 minutes on top of the
download, so this script builds a miniature release -- three compounds, the real
file formats and directory layout -- serves it over HTTP on localhost, and
builds from that. The last section asks the real FTP site which releases exist,
which is one small request.

Run with::

    python examples/pubchem/ftp_build_demo.py
"""

import functools
import gzip
import hashlib
import http.server
import logging
import os
import tempfile
import threading

from provesid import PubChemID
from provesid.pubchem_ftp import build_pubchem_id_db, list_releases, molecular_weight

logging.basicConfig(level=logging.INFO, format="%(levelname)-7s %(message)s")

RELEASE = "2026-09-01"

# ── A miniature release, in PubChem's formats ────────────────────────────────

#: The eight files the builder reads, as PubChem writes them: one line per
#: compound (or per synonym, or per identifier), CID first, tab-separated.
#: CID 999001 has only a CAS number that fails the check digit, so it is left
#: out of every table.
FILES = {
    "CID-Identifiers.tsv.gz": [
        "2244\t50-78-2\tCAS", "2244\tCHEBI:15365\tChEBI ID",
        "2244\tCHEMBL25\tChEMBL ID", "2244\tDTXSID5020108\tDSSTox Substance ID",
        "702\t64-17-5\tCAS", "702\tDTXSID9020584\tDSSTox Substance ID",
        "71583\t865-49-6\tCAS",
        "999001\t50-78-3\tCAS",
    ],
    "CID-Date.gz": ["2244\t2004-09-16", "702\t2004-09-16", "71583\t2005-03-26",
                    "999001\t2020-01-01"],
    "CID-Mass.gz": ["2244\tC9H8O4\t180.04225873\t180.04225873",
                    "702\tC2H6O\t46.041864811\t46.041864811",
                    "71583\tCHCl3\t118.9175601\t118.9175601",
                    "999001\tC6H6\t78.04695\t78.04695"],
    "CID-SMILES.gz": ["2244\tCC(=O)OC1=CC=CC=C1C(=O)O", "702\tCCO",
                      "71583\t[2H]C(Cl)(Cl)Cl", "999001\tC1=CC=CC=C1"],
    "CID-Title.gz": ["2244\tAspirin", "702\tEthanol", "71583\tChloroform-d",
                     "999001\tBenzene"],
    "CID-IUPAC.gz": ["2244\t2-acetyloxybenzoic acid", "702\tethanol",
                     "71583\ttrichloro(deuterio)methane", "999001\tbenzene"],
    "CID-InChI-Key.gz": [
        "2244\tInChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
        "\tBSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        "702\tInChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3\tLFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        "71583\tInChI=1S/CHCl3/c2-1(3)4/h1H/i1D\tHEDRZPFGACZZDS-MICDWDOJSA-N",
        "999001\tInChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H\tUHOVQNZJYSORNB-UHFFFAOYSA-N",
    ],
    "CID-Synonym-filtered.gz": ["2244\tAspirin", "2244\tacetylsalicylic acid",
                                "702\tEthanol", "702\tethyl alcohol",
                                "71583\tChloroform-d", "999001\tbenzene"],
}


def write_release(root):
    """Lay the files out as ``Compound/Monthly/<date>/Extras/``, with ``.md5`` files."""
    extras = os.path.join(root, "Monthly", RELEASE, "Extras")
    os.makedirs(extras)
    for name, lines in FILES.items():
        path = os.path.join(extras, name)
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write("\n".join(lines) + "\n")
        with open(path, "rb") as handle:
            digest = hashlib.md5(handle.read()).hexdigest()
        with open(path + ".md5", "w") as handle:
            handle.write(f"{digest}  {name}\n")
    with open(os.path.join(root, "Monthly", RELEASE, "TIMESTAMP"), "w") as handle:
        handle.write("2026/08/31 18:38:04\n")


def serve(root):
    """Serve ``root`` over HTTP on localhost; return the server and its URL."""
    class Quiet(http.server.SimpleHTTPRequestHandler):
        def log_message(self, *args):
            pass

    httpd = http.server.ThreadingHTTPServer(
        ("127.0.0.1", 0), functools.partial(Quiet, directory=root))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    return httpd, f"http://127.0.0.1:{httpd.server_address[1]}"


def main():
    with tempfile.TemporaryDirectory() as workdir:
        ftp_root = os.path.join(workdir, "pubchem")
        write_release(ftp_root)
        httpd, url = serve(ftp_root)

        # ── 1. Build ────────────────────────────────────────────────────────
        print("\n1. Build the database from the newest snapshot the server lists")
        db_path = build_pubchem_id_db(os.path.join(workdir, "data", "pubchem_id.db"),
                                      base_url=url, progress=False)
        leftovers = sorted(os.listdir(os.path.dirname(db_path)))
        print(f"   built {db_path}")
        print(f"   left in the data directory: {leftovers}  (source files deleted)")

        with PubChemID(db_path=db_path) as db:
            # ── 2. Lookups work as they always have ─────────────────────────
            print("\n2. The lookups PubChemID has always offered")
            print(f"   cas_to_cid('50-78-2')        -> {db.cas_to_cid('50-78-2')}")
            print(f"   cas_to_inchikey('64-17-5')   -> {db.cas_to_inchikey('64-17-5')}")
            print(f"   cas_to_cid('50-78-3')        -> {db.cas_to_cid('50-78-3')}"
                  "   (fails the check digit, so never stored)")

            # ── 3. Properties that are data, offline ────────────────────────
            print("\n3. Properties served from disk")
            print(f"   {sorted(db.offline_properties)}")
            result = db.properties(2244, ["MolecularWeight", "MonoisotopicMass"],
                                   use_online_fallback=False)
            print(f"   aspirin: {result}")
            print(f"   chloroform-d is CHCl3 to PubChem, {molecular_weight('CHCl3')} by "
                  f"its formula; the isotope label in its SMILES makes it "
                  f"{db.get_by_cid(71583)['mw']}")

            # ── 4. Where it came from ───────────────────────────────────────
            print("\n4. Provenance, written into the database")
            record = db.provenance()
            for key in ("release", "release_timestamp", "compounds", "cas_numbers",
                        "cas_rows_rejected", "inchi_source"):
                print(f"   {key:18s} {record[key]}")
            for entry in record["files"][:3]:
                print(f"   {entry['file']:26s} md5 {entry['md5']}  "
                      f"{entry['rows_kept']}/{entry['lines_read']} lines kept")
            print("   ...")

            # ── 5. Cross-references ─────────────────────────────────────────
            print("\n5. Identifiers other databases give aspirin, as PubChem links them")
            print(f"   {db.xrefs(2244)}")

        httpd.shutdown()
        httpd.server_close()

    # ── 6. The real thing ──────────────────────────────────────────────────
    print("\n6. Releases on the real PubChem FTP site")
    try:
        print(f"   {list_releases()}")
    except Exception as exc:  # offline machines still get through the demo
        print(f"   could not reach PubChem: {exc}")
    print("   PubChemID() builds from the first of these when pubchem_id.db is missing;")
    print("   PubChemID(source='zenodo') downloads the 2.2 GiB prebuilt copy instead.")


if __name__ == "__main__":
    main()
