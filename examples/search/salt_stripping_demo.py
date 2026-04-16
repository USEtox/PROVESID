"""Salt stripping demo.

Demonstrates the ``strip_salts=True`` flag, which uses RDKit's SaltRemover
and largest-fragment picker to strip counter-ions before lookup.

The parent SMILES and parent InChIKey are stored in dedicated columns so the
original query is never lost.

Run with::

    uv run python examples/search/salt_stripping_demo.py
"""

from provesid import Search, strip_salts, normalize_structure

# ── Manual strip_salts utility ────────────────────────────────────────────────
salts = [
    "[Na+].[Cl-].CC(=O)O",          # Sodium chloride + acetic acid
    "c1ccccc1.CC(=O)Oc1ccccc1C(=O)O",  # Benzene + aspirin
    "CC(=O)Oc1ccccc1C(=O)O",        # Aspirin alone (no salt)
]

print("=== Manual strip_salts() ===")
for s in salts:
    parent = strip_salts(s)
    norm = normalize_structure(parent) if parent else {}
    print(f"  Input:   {s}")
    print(f"  Parent:  {parent}")
    print(f"  IK:      {norm.get('inchikey')}")
    print()

# ── Integrated into Search ─────────────────────────────────────────────────────
smiles_list = [
    "[Na+].[Cl-].CC(=O)O",
    "CC(=O)Oc1ccccc1C(=O)O",
]

print("=== Search with strip_salts=True ===")
s_search = Search("smiles", strip_salts=True, show_progress=True)
df = s_search.search(smiles_list)
print(
    df[["query", "parent_smiles", "parent_inchikey", "canonical_smiles",
        "InChIKey", "name"]]
    .to_string(index=False)
)
