"""
Demo script for PubChemID CAS conversion methods.

This script demonstrates the new conversion methods added to PubChemID:
- smiles_to_cas()
- name_to_cas()
- formula_to_cas()
- batch_smiles_to_cas()
- batch_name_to_cas()
- batch_formula_to_cas()
"""

from provesid import PubChemID

print("=" * 80)
print("PubChemID CAS Conversion Methods Demo")
print("=" * 80)
print()

# Initialize PubChemID database
print("Initializing PubChemID database...")
db = PubChemID()
print(f"✓ Database loaded: {db.db_path}")
print()

# Get database stats
stats = db.get_stats()
print(f"Database contains {stats['total_compounds']:,} compounds")
print(f"Total CAS numbers: {stats['total_cas_numbers']:,}")
print()

print("=" * 80)
print("1. SMILES to CAS Conversion")
print("=" * 80)
print()

# Test SMILES to CAS
smiles_examples = {
    "C": "Methane",
    "CO": "Methanol",
    "CCO": "Ethanol",
    "CC(=O)OC1=CC=CC=C1C(=O)O": "Aspirin"
}

for smiles, name in smiles_examples.items():
    cas_list = db.smiles_to_cas(smiles)
    if cas_list:
        print(f"{name:15s} ({smiles:30s})")
        print(f"  → CAS: {', '.join(cas_list[:3])}")
        if len(cas_list) > 3:
            print(f"      ... and {len(cas_list) - 3} more")
    else:
        print(f"{name:15s} ({smiles:30s}) → Not found")
    print()

print("=" * 80)
print("2. Chemical Name to CAS Conversion")
print("=" * 80)
print()

# Test name to CAS
names = ["aspirin", "caffeine", "glucose", "ethanol"]

for name in names:
    cas_list = db.name_to_cas(name, exact=False)
    if cas_list:
        print(f"{name.capitalize():15s} → CAS: {', '.join(cas_list[:3])}")
        if len(cas_list) > 3:
            print(f"{'':15s}    ... and {len(cas_list) - 3} more")
    else:
        print(f"{name.capitalize():15s} → Not found")

print()

print("=" * 80)
print("3. Molecular Formula to CAS Conversion")
print("=" * 80)
print()

# Test formula to CAS
formulas = ["H2O", "C9H8O4", "CH4O", "C2H6O"]

for formula in formulas:
    cas_list = db.formula_to_cas(formula, limit=50)
    if cas_list:
        print(f"{formula:10s} → Found {len(cas_list)} compounds")
        print(f"{'':10s}    First 3 CAS: {', '.join(cas_list[:3])}")
        if len(cas_list) > 3:
            print(f"{'':10s}    ... and {len(cas_list) - 3} more")
    else:
        print(f"{formula:10s} → Not found")
    print()

print("=" * 80)
print("4. Batch Conversions")
print("=" * 80)
print()

# Batch SMILES to CAS
print("Batch SMILES to CAS:")
smiles_batch = ["C", "CO", "CCO"]
results = db.batch_smiles_to_cas(smiles_batch)
for smiles, cas_list in results.items():
    if cas_list:
        print(f"  {smiles:6s} → {len(cas_list)} CAS number(s)")
    else:
        print(f"  {smiles:6s} → Not found")
print()

# Batch name to CAS
print("Batch Name to CAS:")
names_batch = ["aspirin", "caffeine", "glucose"]
results = db.batch_name_to_cas(names_batch, exact=False)
for name, cas_list in results.items():
    if cas_list:
        print(f"  {name:10s} → {len(cas_list)} CAS number(s)")
    else:
        print(f"  {name:10s} → Not found")
print()

# Batch formula to CAS
print("Batch Formula to CAS:")
formulas_batch = ["H2O", "CH4", "C9H8O4"]
results = db.batch_formula_to_cas(formulas_batch, limit=50)
for formula, cas_list in results.items():
    if cas_list:
        print(f"  {formula:10s} → {len(cas_list)} CAS number(s)")
    else:
        print(f"  {formula:10s} → Not found")

print()
print("=" * 80)
print("Demo completed!")
print("=" * 80)
