"""Similarity search via Tanimoto fingerprints.

When ``similarity_threshold > 0``, the Search class uses RDKit Morgan
fingerprints to find structurally similar compounds in the offline databases
if an exact structure match fails.

Run with::

    uv run python examples/search/similarity_search_demo.py
"""

from provesid import Search

# SMILES that may not have an exact database entry — use structural similarity
query_smiles = [
    "CC(=O)Oc1ccccc1C(=O)O",   # Aspirin (exact match expected)
    "CC(=O)Nc1ccc(O)cc1",      # Paracetamol
]

print("=== Similarity search (threshold=0.7) ===")
s = Search(
    "smiles",
    similarity_threshold=0.7,
    show_progress=True,
)
df = s.search(query_smiles)
print(
    df[["query", "name", "CASRN", "canonical_smiles", "InChIKey",
        "confidence", "match_method", "match_score"]]
    .to_string(index=False)
)
