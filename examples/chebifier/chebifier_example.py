"""Classify molecules into ChEBI ontology classes with the chebifier backend.

This is an optional, heavy backend. Install it first (Linux/CPU) with::

    bash scripts/install_chebifier.sh

See ``docs/chebifier.md`` for the full installation story, the model-storage
behaviour, and known issues.

Run::

    python examples/chebifier/chebifier_example.py
"""

from provesid.taxonomy import ChebifierClassifier, chebifier_available


def main() -> None:
    # Feature-detect without importing PyTorch. If the extra is missing, a call
    # to classify() would raise a clear ChebifierMissingError instead.
    if not chebifier_available():
        print("chebifier is not installed. Run: bash scripts/install_chebifier.sh")
        return

    # Construct once and reuse; the (expensive) ensemble loads lazily on the
    # first classify() call and its weights are stored in the shared PROVESID
    # dataset directory (set PROVESID_DATA_DIR to relocate).
    classifier = ChebifierClassifier(
        use_cache=True,        # cache predictions on disk by InChIKey
        resolve_names=False,   # set True to also resolve ChEBI ids -> names
    )

    molecules = {
        "benzene": "c1ccccc1",
        "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
        "glucose": "OCC1OC(O)C(O)C(O)C1O",
    }

    # First call loads the model and downloads weights on first ever run.
    df = classifier.classify(list(molecules.values()))

    print("Tidy taxonomy table:")
    print(df[["smiles", "inchikey", "source"]].to_string(index=False))
    print()

    for name, smiles in molecules.items():
        row = df[df["smiles"] == smiles].iloc[0]
        chebi_ids = (row["chebi_ids"] or "").split("|") if row["chebi_ids"] else []
        print(f"{name} ({smiles}): {len(chebi_ids)} ChEBI classes")
        print("  first few:", chebi_ids[:8])

    # Collapse to the {inchikey: label} mapping downstream code (e.g. per-class
    # model evaluation) typically consumes.
    labels = ChebifierClassifier.to_labels(df, level="chebi_ids")
    print("\n{inchikey: chebi_ids} mapping has", len(labels), "entries")

    # A second call over the same structures is served from the on-disk cache.
    df_again = classifier.classify(list(molecules.values()))
    print("Cache round-trip identical:", df["chebi_ids"].equals(df_again["chebi_ids"]))


if __name__ == "__main__":
    main()
