"""How corroboration shapes a Search result's confidence.

A hit's confidence is not just "which database answered". It also accounts for
how many *independent* databases carry the same structure — the
``n_source_support`` column — because a lone source trivially agrees with
itself and would otherwise outrank a structure that three databases share.

This demo shows:

1. every candidate structure for a query, with its support and confidence;
2. how ``min_source_support`` refuses to answer without corroboration;
3. that a group record (a SMILES with an attachment point) is never returned.

Run with::

    uv run python examples/search/confidence_and_corroboration_demo.py
"""

import pandas as pd

from provesid import Search

pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 20)

COLUMNS = ["hit_rank", "name", "InChIKey", "source", "n_source_support", "confidence"]

# 108-62-3 is metaldehyde, the molluscicide in slug pellets.  This CAS used to
# come back as "Mefluidide" — a different pesticide entirely.
CAS = "108-62-3"

# ── 1. Every candidate, ranked ────────────────────────────────────────────────
print(f"All candidate structures for CAS {CAS}:\n")
every_hit = Search("cas", n_hits="all", show_progress=False).search(CAS)
print(every_hit[COLUMNS].to_string(index=False))

print(
    "\nThe top hit is the one the most databases agree on.  Confidence is\n"
    "  base × query_term × (0.5 + 0.5 × consensus) × support_factor\n"
    "where support_factor is 0.85 for a single database, 0.95 for two, and\n"
    "1.0 from three upwards.\n"
)

# ── 2. Requiring corroboration ────────────────────────────────────────────────
# Some CAS numbers are carried by only one database — sodium stibogluconate is
# in ChEBI and nowhere else offline.  Whether an uncorroborated answer is good
# enough depends on what you are doing with it, so it is your call to make.
mixed = [
    "108-62-3",     # metaldehyde  — five databases agree
    "16037-91-5",   # sodium stibogluconate — ChEBI only
    "1234567-89-0",  # not a real CAS
]

lenient = Search("cas", show_progress=False).search(mixed)
strict = Search("cas", min_source_support=2, show_progress=False).search(mixed)

print("Default (any hit accepted):\n")
print(lenient[["query", "name", "n_source_support", "confidence"]].to_string(index=False))
print("\nWith min_source_support=2 (two databases must agree):\n")
print(strict[["query", "name", "n_source_support", "confidence"]].to_string(index=False))
print(
    "\nThe single-source hit is dropped rather than returned at face value; it\n"
    "already scored 0.765 rather than the 0.90 base for an exact CAS match.\n"
)

# ── 3. Group records are never substances ─────────────────────────────────────
# ChEBI indexes substituents ("…dienoyl group", SMILES `*C(=O)CCCC=CCC=CCCCCC`)
# alongside real compounds.  A registry number identifies a substance, so those
# records are dropped from identifier lookups; RDKit cannot process them either.
print("\nNo returned structure carries an attachment point:")
resolved = Search("cas", n_hits="all", show_progress=False).search(mixed)
smiles = [s for s in resolved["SMILES"] if isinstance(s, str)]
print(f"  checked {len(smiles)} structures, "
      f"{sum('*' in s for s in smiles)} contain '*'")
