"""Can the structure in a Search row be used as it is?

Everything downstream of ``Search`` --- a QSAR model, a property estimate, a
characterisation factor --- takes the structure in the row at face value.
These columns say whether it should:

1. ``structure_status``: the source's SMILES as given, or rebuilt because it
   did not parse (ChEBI's export drops aromatic ring N-H: valsartan, the
   triptans, allopurinol, 5-fluorouracil), or no usable structure at all;
2. ``identity_consistent``: whether the InChIKey and formula the sources state
   belong to that structure;
3. ``conflict``: whether an independent database names a different compound
   for the same query, and which (``runner_up_*``);
4. ``validated_by``: which checks against your own reference values a hit
   passed, with ``validate=``.

The CAS numbers are from a structure audit of 630 pharmaceuticals, where each
of these went wrong before the checks existed.

Run with::

    uv run python examples/search/structure_qc_demo.py
"""

import pandas as pd

from provesid import Search

pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 20)
pd.set_option("display.max_colwidth", 40)

CAS = {
    "137862-53-4": "valsartan: ChEBI's SMILES lacks the tetrazole N-H",
    "51-21-8": "5-fluorouracil: same, rebuilt from ChEBI's molfile",
    "118443-89-3": "cefquinome sulfate: CompTox and PubChem name different compounds",
    "27314-97-2": "tirapazamine: CompTox's record has the wrong formula",
    "123997-26-2": "eprinomectin: CompTox lists a structure-less mixture",
}

with Search("cas", show_progress=False) as search:
    df = search.search(list(CAS))

# ── 1. Where each structure came from ─────────────────────────────────────────
print("Structure status, and whether the stated identifiers fit it:\n")
print(df[["query", "name", "structure_status", "identity_consistent",
          "molecular_formula", "source_formula"]].to_string(index=False))
print(
    "\n'repaired:molfile' means the source's SMILES did not parse and the\n"
    "structure was rebuilt from its own molfile, accepted only because the\n"
    "result has the InChIKey the source states.  Before, these rows carried an\n"
    "unreadable SMILES, no canonical SMILES, and a normal-looking confidence.\n"
)

# ── 2. Disagreement between databases ─────────────────────────────────────────
print("Conflicts between independent databases:\n")
print(df[["query", "name", "source", "n_source_support", "conflict", "n_clusters",
          "runner_up_sources", "runner_up_confidence", "confidence"]].to_string(index=False))
print(
    "\nA contested answer's confidence is multiplied by its share of the\n"
    "independent databases: cefquinome's one-against-one split halves it.\n"
    "The strict preset returns nothing for a contested query instead:\n"
)
with Search("cas", preset="strict", show_progress=False) as strict:
    refused = strict.search(["118443-89-3", "27314-97-2"])
print(refused[["query", "name", "conflict", "n_clusters"]].to_string(index=False))

# ── 3. A DTXSID query is checked the same way ─────────────────────────────────
print("\nAtrazine, tirapazamine and cefquinome sulfate asked by DTXSID:\n")
with Search("dtxsid", show_progress=False) as by_dtxsid:
    dtx = by_dtxsid.search(["DTXSID9020112", "DTXSID30865361", "DTXSID6048593"])
print(dtx[["query", "name", "molecular_formula", "source", "n_source_support",
           "conflict", "confidence"]].to_string(index=False))
print(
    "\nOnly CompTox holds DTXSIDs, so the other databases are asked for the CAS\n"
    "number CompTox gives, and each answers from its own record: atrazine is\n"
    "confirmed three times, tirapazamine's wrong CompTox formula loses to ChEBI\n"
    "and PubChem, and cefquinome's split is flagged.  Asked for CompTox's own\n"
    "InChIKey instead, they could only agree with it.\n"
)

# ── 4. Checking hits against reference values ─────────────────────────────────
print("\nValidating against the molecular weights a measurement table carries:\n")
table = pd.DataFrame({
    "CAS": ["27314-97-2", "137862-53-4", "118443-89-3"],
    "MW": [178.15, 435.52, 500.0],          # the last matches neither reading
})
with Search("cas", validate={"mw": "MW"}, show_progress=False) as checked:
    validated = checked.search(table, column="CAS")
print(validated[["query", "name", "molecular_mass", "MW", "validated_by"]].to_string(index=False))
print(
    "\nA hit that fails is passed over for the next-ranked compound; when none\n"
    "passes, the query gets a no-match row.\n"
)

# ── 5. Which data produced this frame ─────────────────────────────────────────
print("Dataset provenance recorded in df.attrs['dataset_versions']:\n")
for source, entry in df.attrs["dataset_versions"].items():
    print(f"  {source:9s} {entry}")
