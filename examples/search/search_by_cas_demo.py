"""Search by CAS Registry Number — basic demo.

Demonstrates the most common use-case: resolving a list of CAS numbers into
a unified DataFrame with cross-database identifiers.

Run with::

    uv run python examples/search/search_by_cas_demo.py
"""

from provesid import Search

cas_list = [
    "50-00-0",   # Formaldehyde
    "64-17-5",   # Ethanol
    "50-78-2",   # Aspirin
    "7732-18-5", # Water
    "58-08-2",   # Caffeine
]

# Create a Search instance for CAS number lookup.
# Clients (ChEBI, CompTox, PubChemID, ZeroPM, ChEMBL) are initialised
# lazily on the first call to search().
s = Search("cas", show_progress=True)

df = s.search(cas_list)

print("Resolved", len(df), "compounds\n")
print(
    df[["query", "CASRN", "name", "canonical_smiles", "InChIKey",
        "confidence", "source", "match_score"]]
    .to_string(index=False)
)
