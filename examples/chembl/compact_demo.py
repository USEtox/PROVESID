"""
Shrink a ChEMBL release to the part PROVESID actually reads.

A full ChEMBL release is about 30 GB spread over 74 tables.  PROVESID opens
eight of them and reads no bioactivity data at all, so almost the whole file is
dead weight on your disk.  ``CheMBL.compact()`` copies those eight tables, drops
the ``molfile`` column -- a quarter of the entire database on its own, and
consumed nowhere in this package -- and rebuilds only the indexes the queries
here need.

Measured on ChEMBL 36: **29.74 GB to 2.60 GB**, with every public method
returning the same compounds.

Run with::

    python examples/chembl/compact_demo.py

The script never deletes anything.  The last section shows the one call that
does, so you can reclaim the space once you are satisfied the extract answers
your own queries.
"""

import logging
import os

from provesid import CheMBL

# compact() reports its progress through the standard logging module, like the
# rest of PROVESID -- turn it on to watch the tables go by.
logging.basicConfig(level=logging.INFO, format="%(levelname)-7s %(message)s")


def gigabytes(path):
    """Size of a file in GB, for printing."""
    return os.path.getsize(path) / 1e9


def main():
    # ── 1. What is on disk now ────────────────────────────────────────────────
    # Opening CheMBL with no arguments reuses whatever release is already there.
    chembl = CheMBL(auto_download=False)
    print(f"\nOpen database : {chembl.db_path}")
    print(f"Release       : {chembl.release}")
    print(f"Size          : {gigabytes(chembl.db_path):.2f} GB")
    print(f"Is an extract : {chembl.is_compact}")

    if chembl.is_compact:
        print("\nThis is already a PROVESID extract. Its provenance:")
        for key, value in sorted(chembl.provenance.items()):
            if not key.startswith("rows."):
                print(f"  {key:18} {value}")
        return

    # ── 2. A reference answer from the full database ──────────────────────────
    # Keep one result to compare against after compaction.
    aspirin_molregno = chembl.chembl_id_to_molregno("CHEMBL25")
    before = chembl.get_compound(aspirin_molregno)
    print(f"\nAspirin from the full database: {before['pref_name']}, "
          f"{before['canonical_smiles']}")

    # ── 3. Build the extract ──────────────────────────────────────────────────
    # The source is opened read-only and the extract is verified against it
    # before it replaces anything, so this is safe to run at any time.
    print("\nBuilding the extract...")
    extract_path = CheMBL.compact(force=True)

    source_gb = gigabytes(chembl.db_path)
    extract_gb = gigabytes(extract_path)
    print(f"\n  {chembl.db_path}  {source_gb:8.2f} GB")
    print(f"  {extract_path}  {extract_gb:8.2f} GB")
    print(f"  saved {source_gb - extract_gb:.2f} GB "
          f"({100 * (1 - extract_gb / source_gb):.1f}% smaller)")

    # ── 4. The extract answers identically ────────────────────────────────────
    # A plain CheMBL() now prefers the extract over the full release beside it.
    compact = CheMBL(auto_download=False)
    print(f"\nCheMBL() now opens : {compact.db_path}")
    print(f"Is an extract      : {compact.is_compact}")

    after = compact.get_compound(compact.chembl_id_to_molregno("CHEMBL25"))
    print(f"Aspirin from the extract: {after['pref_name']}, "
          f"{after['canonical_smiles']}")
    print(f"Identical result   : {before == after}")

    # get_compound no longer carries a molfile. Build one with RDKit if needed:
    #     from rdkit import Chem
    #     Chem.MolToMolBlock(Chem.MolFromSmiles(after["canonical_smiles"]))
    print(f"molfile in result  : {'molfile' in after}")

    # ── 5. Where it came from ─────────────────────────────────────────────────
    # An extract is a subset, so it records what built it. A later PROVESID that
    # needs a ninth table can then say so instead of failing with 'no such table'.
    print("\nProvenance:")
    for key, value in sorted(compact.provenance.items()):
        if not key.startswith("rows."):
            print(f"  {key:18} {value}")
    print("\nRows copied:")
    for key, value in sorted(compact.provenance.items()):
        if key.startswith("rows."):
            print(f"  {key[5:]:28} {int(value):>12,}")

    # ── 6. Reclaiming the space ───────────────────────────────────────────────
    # Nothing above deleted anything. To free the ~27 GB, either delete the full
    # database yourself or let compact() do it once it has verified the extract:
    #
    #     CheMBL.compact(remove_source=True, force=True)
    #
    # remove_source only ever runs after verification passes.
    print(f"\nThe full database is still at {chembl.db_path} "
          f"({source_gb:.2f} GB).")
    print("Reclaim it with:  CheMBL.compact(remove_source=True, force=True)")


if __name__ == "__main__":
    main()
