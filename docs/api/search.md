# Search API Reference

The `Search` class is the primary entry point for resolving chemical identifiers
across the offline databases ChEBI, CompTox, PubChemID and ChEMBL.

## Targeted sources

`Search` queries four sources by default: **ChEBI**, **CompTox**, **PubChemID**
and **ChEMBL**.

**ZeroPM is not among them.** It aggregates regulatory inventories rather than
curating compounds, so its name→structure rows are noticeably noisier than the
other four — yet it counted as a full independent vote in the corroboration
ranking that drives `confidence` and `min_source_support`. The
`ZeroPM` class itself is unchanged and remains fully available for direct use
(see the [ZeroPM tutorial](../examples/zeropm/zeropm-example.md)); only the
resolver stopped targeting it.

Pass `use_zeropm=True` to put it back in the pool:

```python
from provesid import Search

# Default: four sources, ZeroPM absent from source_details.
Search("cas").search("50-00-0")

# Opt back in — mainly useful for fuzzy name queries, since ZeroPM is the only
# source that does true fuzzy *retrieval* (e.g. "caffiene" -> caffeine).
Search("name", fuzzy=True, use_zeropm=True).search("caffiene")
```

A `zeropm=` client passed to the constructor is ignored (with a warning) unless
`use_zeropm=True` is also set.

What each source can be asked, and by which identifier, is one table in
`provesid.sources` (see [Source lookup table](#source-lookup-table) below).
A gap in it is a source with no index for that identifier, which `Search`
reaches through one it does have, such as ChEMBL through the SMILES another
source found for a CAS number.

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

## Source lookup table

::: provesid.sources
    options:
      show_source: false
      members:
        - Query
        - LOOKUPS
        - rank_rows_by_completeness

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
| `hit_rank` | int | Rank among the hits for this query (0 = best) |
| `n_source_support` | int | Independent databases carrying this structure |

---

## Confidence scoring

Confidence combines four signals: how strong the match method is, how well the
candidate matches the query itself, how well the databases that answered agree
with each other, and how many of them carried the structure at all:

$$
\text{confidence} = \text{base}
\times (w_q \times \text{query\_score} + (1 - w_q))
\times (0.5 + 0.5 \times \text{consensus\_score})
\times \text{support}(n)
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
| Fuzzy name | rapidfuzz ratio × 0.80 |
| Formula | 0.30 |

For exact-identifier methods `query_score` is 1.0, which collapses the second
term. `w_q` is the `query_weight` argument (default 0.5).

### Corroboration

`support(n)` scales the score by the number of *independent* databases carrying
the structure (`n_source_support` in the output):

| Databases agreeing | Factor |
|---|---|
| 0 (OPSIN-only parse) | 1.00 |
| 1 | 0.85 |
| 2 | 0.95 |
| 3 or more | 1.00 |

This term is what keeps provenance from beating evidence. `consensus_score`
measures how well the sources that answered agree — but a lone source agrees
with itself perfectly, so without `support(n)` a single database hit (0.90)
outranks a structure that three databases all carry (0.8777), and the resolver
returns the compound nothing corroborates.

Use `min_source_support=` (on the constructor or per call) to require
corroboration outright: `Search("cas", min_source_support=2)` returns only
structures that at least two databases agree on.
