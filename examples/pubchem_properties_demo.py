"""
Demo script for PubChem property retrieval: in bulk, and offline first.

Two things are demonstrated:

1. ``PubChemAPI.get_properties_for_cids()`` / ``get_compound_properties_batch()``
   ask PubChem about a whole list of compounds in one request instead of one
   request per compound.
2. ``PubChemID.properties()`` and friends read the local SQLite database first
   and fall back to PUG-REST only for what it cannot answer, which for most
   property lookups means no network traffic at all.

Run it with the local database already downloaded; the first ``PubChemID()``
otherwise fetches ~2.3 GB.
"""

import logging
import time

from provesid import PubChemAPI, PubChemID
from provesid.pubchem import PROPERTY_CHUNK_SIZE

# The offline/online decisions are logged at DEBUG, which is the easiest way to
# see which source answered a given call.
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

print("=" * 80)
print("PubChem Property Retrieval Demo")
print("=" * 80)
print()

# ---------------------------------------------------------------------------
# 1. Bulk retrieval from the online API
# ---------------------------------------------------------------------------
print("1. Bulk retrieval: many CIDs, few requests")
print("-" * 80)

api = PubChemAPI()
cids = [2244, 702, 5793, 1983, 3672]  # aspirin, ethanol, glucose, paracetamol, ibuprofen

start = time.time()
rows = api.get_properties_for_cids(cids, ["MolecularFormula", "MolecularWeight"])
elapsed = time.time() - start

for row in rows:
    print(f"  CID {row['CID']:>7}  {row['MolecularFormula']:<12} {row['MolecularWeight']}")
print(f"  {len(rows)} compounds in one request, {elapsed:.2f}s")
print(f"  (lists longer than {PROPERTY_CHUNK_SIZE} CIDs are split into that many per request)")
print()

# The same request reshaped: one dict per CID, including the ones PubChem has no
# record of, so the result can be zipped against the input.
print("  Reshaped, with an unknown CID included:")
for row in api.get_compound_properties_batch([2244, 999999999], ["MolecularFormula"]):
    if row["success"]:
        print(f"    CID {row['cid']}: {row['MolecularFormula']}")
    else:
        print(f"    CID {row['cid']}: {row['error']}")
print()

# ---------------------------------------------------------------------------
# 2. Offline first
# ---------------------------------------------------------------------------
print("2. Offline first: the local database answers what it can")
print("-" * 80)

db = PubChemID()
print(f"  Database: {db.db_path}")
print()

# Everything asked for here has a column in the local database, so this costs
# no request at all.
start = time.time()
result = db.properties(2244, ["MolecularFormula", "MolecularWeight", "XLogP", "SMILES"])
print(f"  db.properties(2244, [...]) -> Source={result['Source']} "
      f"in {(time.time() - start) * 1000:.1f} ms")
for key, value in result.items():
    if key not in ("CID", "Source"):
        print(f"    {key:<18} {value}")
print()

# MonoisotopicMass has no local column, so this one goes online.
result = db.properties(2244, ["MonoisotopicMass", "MolecularWeight"])
print(f"  db.properties(2244, ['MonoisotopicMass', ...]) -> Source={result['Source']}")
print(f"    MonoisotopicMass   {result['MonoisotopicMass']}")
print()
print("  Locally available properties:")
print(f"    {', '.join(sorted(PubChemID.OFFLINE_PROPERTIES))}")
print()

# A property the compound genuinely has no value for is absent from the result
# rather than None -- the same way PubChem reports it. CID 233 (arsenate) has no
# computed logP.
result = db.properties(233, ["MolecularFormula", "XLogP"])
print(f"  CID 233 has no logP in PubChem: 'XLogP' in result -> {'XLogP' in result}")
print()

# ---------------------------------------------------------------------------
# 3. A bulk lookup that splits itself between the two sources
# ---------------------------------------------------------------------------
print("3. Bulk lookup across both sources")
print("-" * 80)

# A thousand low CIDs: most are in the local database, the rest are fetched in
# one batched request.
many = list(range(1000, 2000))

start = time.time()
table = db.properties_table(many, ["MolecularFormula", "MolecularWeight"])
elapsed = time.time() - start

counts = table["Source"].value_counts().to_dict()
print(f"  {len(table)} CIDs in {elapsed:.2f}s")
for source in ("offline", "online", "missing"):
    if source in counts:
        print(f"    {source:<8} {counts[source]}")
print()
print(table.head().to_string(index=False))
print()

# Strictly local, for a script that must not touch the network.
offline_only = db.properties_table(many, ["MolecularWeight"], use_online_fallback=False)
served = (offline_only["Source"] == "offline").sum()
print(f"  use_online_fallback=False: {served}/{len(offline_only)} answered locally, "
      "no requests made")
print()

print("=" * 80)
print("Demo completed!")
print("=" * 80)
