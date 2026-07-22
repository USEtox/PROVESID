"""Attach ChEBI class labels to a set of chemicals and group by class.

This mirrors the downstream use case the chebifier backend is built for: take a
table of chemicals, classify each into ChEBI ontology classes, and attach a
class label so results can be grouped/aggregated (e.g. per-class model metrics).

It also shows:
  * name resolution (ChEBI id -> human-readable name),
  * ``to_labels`` producing the ``{inchikey: label}`` mapping,
  * the on-disk cache (a second run is served from cache).

Install the backend first (Linux/CPU)::

    bash scripts/install_chebifier.sh

Run::

    python examples/chebifier/chebifier_grouping_example.py
"""

import pandas as pd

from provesid.cache import get_chebifier_cache_info
from provesid.taxonomy import ChebifierClassifier, chebifier_available


def main() -> None:
    if not chebifier_available():
        print("chebifier is not installed. Run: bash scripts/install_chebifier.sh")
        return

    # A small "evaluation anchor" table: chemicals we want to group by class.
    chemicals = pd.DataFrame(
        {
            "name": ["benzene", "toluene", "phenol", "glucose", "fructose"],
            "smiles": [
                "c1ccccc1",
                "Cc1ccccc1",
                "Oc1ccccc1",
                "OCC1OC(O)C(O)C(O)C1O",
                "OCC1(O)OCC(O)C(O)C1O",
            ],
        }
    )

    clf = ChebifierClassifier(use_cache=True, resolve_names=True)

    # Classify all structures in one batched call; keep the InChIKey join key.
    taxonomy = clf.classify(chemicals["smiles"].tolist())

    # Use the most specific (first) predicted ChEBI class as a single label for
    # grouping. In a real pipeline you might pick a fixed ontology level instead.
    def first_label(row: pd.Series) -> str:
        ids = (row["chebi_ids"] or "").split("|") if row["chebi_ids"] else []
        names = (row["chebi_names"] or "").split("|") if row["chebi_names"] else []
        if not ids:
            return "unclassified"
        return names[0] if names else ids[0]

    taxonomy["chem_class"] = taxonomy.apply(first_label, axis=1)

    # Attach the class label back onto the input table via SMILES.
    merged = chemicals.merge(
        taxonomy[["smiles", "inchikey", "chem_class"]], on="smiles", how="left"
    )
    print("Chemicals with attached ChEBI class label:")
    print(merged.to_string(index=False))

    # {inchikey: label} mapping — the shape downstream grouping code consumes.
    labels = ChebifierClassifier.to_labels(taxonomy, level="chebi_names")
    print("\n{inchikey: chebi_names} entries:", len(labels))

    # Group the anchor by the assigned class (toy per-class aggregation).
    print("\nCounts per assigned class:")
    print(merged.groupby("chem_class")["name"].count().to_string())

    # Cache visibility — the second classify() is served from disk.
    info = get_chebifier_cache_info()
    print(f"\nchebifier cache: {info.get('disk_entries')} entries, "
          f"{info.get('total_size_mb', 0):.2f} MB")


if __name__ == "__main__":
    main()
