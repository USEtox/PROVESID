"""
Computed descriptors on demand: RDKit's from the stored SMILES, or PubChem's.

``pubchem_id.db`` stores identifiers and structures, not descriptors. XLogP,
TPSA and the atom and bond counts are the output of a model run over the
structure, and PubChem's model (Cactvs, XLogP3) is not the only one. So
``PubChemID.descriptors()`` computes them when asked, and says which model
answered:

* ``source="rdkit"`` (default): RDKit, from the SMILES already on disk. No
  network, about half a millisecond per compound. ``Source='rdkit'``.
* ``source="pubchem"``: PubChem's own values over PUG-REST.
  ``Source='online'``. The only way to ``XLogP`` and ``Complexity``.

So that this runs anywhere, the local database here is a miniature one with
three compounds, written to a temporary directory; with the real database,
``PubChemID()`` is all it takes. The last section makes one request to PubChem
and is skipped if PubChem cannot be reached.

Run with::

    python examples/pubchem/descriptors_demo.py
"""

import logging
import os
import sqlite3
import tempfile

import pandas as pd

from provesid import PubChemError, PubChemID
from provesid.pubchem import PUBCHEM_DESCRIPTORS, RDKIT_DESCRIPTORS, rdkit_descriptors

logging.basicConfig(level=logging.INFO, format="%(levelname)-7s %(message)s")
pd.set_option("display.width", 120)

COMPOUNDS = [
    (2244, "Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O"),
    (702, "Ethanol", "CCO"),
    (2519, "Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
]


def miniature_database(directory: str) -> str:
    """Write a ``pubchem_id.db`` with the tables PubChemID reads and three compounds."""
    path = os.path.join(directory, "pubchem_id.db")
    with sqlite3.connect(path) as conn:
        conn.executescript("""
            CREATE TABLE compounds (cid INTEGER PRIMARY KEY, cmpdname TEXT,
                mf TEXT, inchi TEXT, smiles TEXT, inchikey TEXT, iupacname TEXT,
                mw REAL, exactmass REAL, monoisotopicmass REAL, cidcdate TEXT);
            CREATE TABLE cas_numbers (id INTEGER PRIMARY KEY, cid INTEGER, cas TEXT);
            CREATE TABLE synonyms (id INTEGER PRIMARY KEY, cid INTEGER, synonym TEXT);
        """)
        conn.executemany("INSERT INTO compounds (cid, cmpdname, smiles) VALUES (?, ?, ?)",
                         COMPOUNDS)
    return path


with tempfile.TemporaryDirectory() as workdir:
    with PubChemID(db_path=miniature_database(workdir), auto_download=False) as db:

        # ── One compound, every RDKit descriptor ─────────────────────────────
        print("\nRDKit computes:", ", ".join(RDKIT_DESCRIPTORS))
        print("PubChem publishes:", ", ".join(PUBCHEM_DESCRIPTORS))
        print("\ndb.descriptors(2244):")
        print(" ", db.descriptors(2244))

        # ── A table, strictly offline ────────────────────────────────────────
        # use_online_fallback=False guarantees no request: a CID the database
        # does not hold (here 5793, glucose) gets Source='missing' rather than
        # having its SMILES fetched.
        print("\ndb.descriptors_table([2244, 702, 2519, 5793], use_online_fallback=False):")
        print(db.descriptors_table([2244, 702, 2519, 5793],
                                   ["MolLogP", "TPSA", "HBondDonorCount"],
                                   use_online_fallback=False).to_string(index=False))

        # ── Names that belong to the other source ────────────────────────────
        # RDKit's logP is Crippen's model, not XLogP3, so it is not called XLogP.
        try:
            db.descriptors(2244, ["XLogP"])
        except ValueError as error:
            print(f"\nAsking RDKit for XLogP: {error}")

        # ── PubChem's own values, side by side with RDKit's ──────────────────
        common = ["TPSA", "HBondDonorCount", "HBondAcceptorCount", "RotatableBondCount"]
        cids = [cid for cid, _, _ in COMPOUNDS]
        try:
            pubchem = db.descriptors_table(cids, common, source="pubchem")
        except PubChemError as error:
            print(f"\nPubChem could not be reached ({error}); skipping the comparison.")
        else:
            rdkit = db.descriptors_table(cids, common)
            print("\nThe same names, two models (rdkit above, PubChem below):")
            print(pd.concat([rdkit, pubchem]).to_string(index=False))

# ── A structure that is not in PubChem ───────────────────────────────────────
print("\nrdkit_descriptors('OCC(O)CO'):", rdkit_descriptors("OCC(O)CO"))
