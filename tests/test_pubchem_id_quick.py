"""Quick test of PubChemID class."""

from provesid import PubChemID

# Initialize
print("Initializing PubChemID...")
db = PubChemID()

# Get stats
print("\nDatabase Statistics:")
stats = db.get_stats()
print(f"  Total compounds: {stats['total_compounds']:,}")
print(f"  Total CAS numbers: {stats['total_cas_numbers']:,}")
print(f"  Compounds with CAS: {stats['compounds_with_cas']:,}")
print(f"  Total synonyms: {stats['total_synonyms']:,}")
print(f"  Compounds with InChIKey: {stats['compounds_with_inchikey']:,}")
print(f"  Database size: {stats['database_size_mb']:.2f} MB")

# Test lookups
print("\n--- Test 1: Lookup by CAS (Aspirin: 50-78-2) ---")
result = db.get_by_cas("50-78-2")
if result:
    print(f"  CID: {result['cid']}")
    print(f"  Name: {result['cmpdname']}")
    print(f"  Formula: {result['mf']}")
    print(f"  InChIKey: {result['inchikey']}")
    print(f"  CAS numbers: {result['cas_numbers']}")
else:
    print("  Not found")

# Test conversion
print("\n--- Test 2: CAS to InChI conversion ---")
inchi = db.cas_to_inchi("50-78-2")
if inchi:
    print(f"  InChI: {inchi[:50]}...")
else:
    print("  Not found")

# Test batch conversion
print("\n--- Test 3: Batch CAS to CID conversion ---")
cas_list = ["50-78-2", "50-00-0", "64-17-5"]  # Aspirin, Formaldehyde, Ethanol
results = db.batch_cas_to_cid(cas_list)
for cas, cid in results.items():
    print(f"  {cas} -> CID {cid}")

# Test identifier table
print("\n--- Test 4: Get identifier table ---")
df = db.get_id_table_from_cas("50-78-2")
if df is not None:
    print(df.to_string())
else:
    print("  Not found")

# Test search by name
print("\n--- Test 5: Search by name (partial match: 'aspirin') ---")
results = db.search_by_name("aspirin", exact=False, limit=3)
for r in results:
    print(f"  CID {r['cid']}: {r['cmpdname']}")

print("\n✓ All tests completed!")
