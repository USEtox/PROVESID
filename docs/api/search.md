# Search API Reference

The `Search` class is the primary entry point for resolving chemical identifiers
across all offline databases (ChEBI, CompTox, PubChemID, ZeroPM, ChEMBL).

::: provesid.search.Search
    options:
      show_source: false
      show_bases: false
      members:
        - __init__
        - search

---

## Module utilities

::: provesid.search.normalize_structure
    options:
      show_source: false

::: provesid.search.strip_salts
    options:
      show_source: false

---

## Output schema

Every row returned by `Search.search()` contains the following columns:

| Column | Type | Description |
|---|---|---|
| `query` | str | Original input value |
| `CASRN` | str \| None | CAS Registry Number |
| `name` | str \| None | Preferred common name |
| `IUPAC_name` | str \| None | IUPAC systematic name |
| `molecular_formula` | str \| None | Molecular formula |
| `SMILES` | str \| None | Original SMILES from source |
| `canonical_smiles` | str \| None | RDKit-canonical SMILES |
| `kekulized_smiles` | str \| None | Kekulized SMILES from RDKit |
| `InChI` | str \| None | InChI string |
| `InChIKey` | str \| None | Full InChIKey (always reported) |
| `DTXSID` | str \| None | CompTox substance identifier |
| `molecular_mass` | float \| None | Molecular weight |
| `Synonyms` | str \| None | Semicolon-separated synonyms |
| `parent_smiles` | str \| None | SMILES after salt stripping (opt-in) |
| `parent_inchikey` | str \| None | InChIKey after salt stripping (opt-in) |
| `foundby` | str | How the match was found |
| `source` | str \| None | Source that provided the primary SMILES |
| `source_details` | dict | Per-source traceability record |
| `confidence` | float | Overall confidence score [0, 1] |
| `match_method` | str | Matching method used |
| `match_score` | float | Cross-source consensus score [0, 1] |
| `consensus_source` | str \| None | Source chosen by consensus algorithm |
| `source_match_scores` | dict | Per-source agreement scores |

---

## Confidence scoring

Confidence is computed as:

$$
\text{confidence} = \text{base} \times (0.5 + 0.5 \times \text{consensus\_score})
$$

| Match method | Base confidence |
|---|---|
| Exact InChIKey | 1.00 |
| Exact canonical SMILES | 0.95 |
| InChI | 0.95 |
| Exact CAS | 0.90 |
| DTXSID | 0.90 |
| Exact name | 0.80 |
| InChIKey skeleton | 0.75 |
| Tanimoto similarity | tanimoto × 0.85 |
| Fuzzy name | rapidfuzz ratio |
| Formula | 0.30 |
