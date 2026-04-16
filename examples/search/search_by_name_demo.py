"""Search by compound name — exact and fuzzy matching demo.

Shows how to resolve chemical names using either exact matching or fuzzy
string matching (RapidFuzz) when the exact spelling is uncertain.

Run with::

    uv run python examples/search/search_by_name_demo.py
"""

from provesid import Search

names_exact = ["aspirin", "caffeine", "ethanol", "paracetamol"]

# ── exact name search ─────────────────────────────────────────────────────────
print("=== Exact name search ===")
s_exact = Search("name", show_progress=True)
df_exact = s_exact.search(names_exact)
print(
    df_exact[["query", "name", "CASRN", "canonical_smiles", "confidence",
              "match_method"]]
    .to_string(index=False)
)

# ── fuzzy name search ─────────────────────────────────────────────────────────
names_fuzzy = ["aspirine", "kafein", "ethaanol"]  # deliberate misspellings

print("\n=== Fuzzy name search (similarity_threshold=0.7) ===")
s_fuzzy = Search("name", fuzzy=True, similarity_threshold=0.7, show_progress=True)
df_fuzzy = s_fuzzy.search(names_fuzzy)
print(
    df_fuzzy[["query", "name", "CASRN", "confidence", "match_method", "match_score"]]
    .to_string(index=False)
)
