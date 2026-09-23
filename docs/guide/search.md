# Resolving identifiers with `Search`

`Search` takes a column of identifiers of one kind — CAS numbers, names,
SMILES, InChIs, InChIKeys, DTXSIDs or formulas — asks every installed offline
database about each, and returns one table: the structure the databases agree
on, every identifier they hold for it, where each came from, and how confident
the answer is. The [Search tutorial](../examples/search/search_tutorial.ipynb)
works through a real dataset with it.

```python
from provesid import Search

with Search("cas") as s:
    df = s.search(["50-78-2", "39400-72-1", "1912-24-9"])

df[["query", "name", "InChIKey", "source", "n_source_support", "confidence"]]
```

```text
        query                  name                     InChIKey   source  n_source_support  confidence
0     50-78-2  acetylsalicylic acid  BSYNRYMUTXBXSQ-UHFFFAOYSA-N    ChEBI                 4       0.892
1  39400-72-1              Atrazine  MXWJVTOOROXGIU-UHFFFAOYSA-N  CompTox                 1       0.765
2   1912-24-9              atrazine  MXWJVTOOROXGIU-UHFFFAOYSA-N    ChEBI                 4       0.900
```

The second query is a retired CAS number for atrazine. Only CompTox still
lists it, so the answer is right but carries one database's support, and its
confidence says so.

`search()` also takes a DataFrame, or a `.csv` or `.parquet` path, with
`column=` naming the identifiers. `enrich(df, column)` adds the result columns
to your own table instead, resolving each distinct value once — the cheap way
to attach identifiers to a measurement table where one compound fills many
rows. The full argument list is on the [API page](../api/search.md).

## Sources

By default `Search` queries four databases: **PubChem** (`PubChemID`),
**CompTox** (`CompToxID`), **ChEBI** (`ChebiSDF`) and **ChEMBL** (`CheMBL`).
It uses the ones that are installed, and logs once which are not; see
[Installing the offline databases](datasets.md).

**ZeroPM is left out** unless `use_zeropm=True`. It aggregates regulatory
inventories rather than curating compounds, so its name→structure rows are
noisier than the others, and as a full vote in the corroboration count they
used to push wrong structures up the ranking. It is also the one source that
retrieves misspelled names, which is why the `"recall"` preset turns it on.

What each source can be asked, and by which identifier, is one table,
`provesid.sources.LOOKUPS`. A gap in it is a source with no index for that
identifier, which `Search` reaches through one it does have — ChEMBL through
the SMILES another source found for a CAS number, for instance.

## Presets

The arguments that decide what counts as a match and what is returned have
three named settings, `Search.PRESETS`:

| preset | what changes from `balanced` | for |
|---|---|---|
| `"balanced"` | nothing: the constructor's defaults | general use |
| `"strict"` | `min_source_support=2` | tables where a wrong structure costs more than a missing one |
| `"recall"` | `fuzzy`, `inchikey_skeleton`, `similarity_threshold=0.7`, `use_zeropm`, `n_hits="all"` | finding candidates to review by hand |

```python
Search.PRESETS["strict"]                          # a plain dict: inspect it
Search("name", preset="strict").search("atrazine")
Search("name", preset="recall", n_hits=3)         # explicit arguments win
df.attrs["preset"], df.attrs["settings"]          # what produced this frame
```

An argument passed explicitly overrides the preset, even when it equals the
`balanced` value: `preset="recall", use_zeropm=False` leaves ZeroPM out.

## Online fallback

`Search` opens no socket by default. With `online_fallback=True`, a query that
no offline source answered, and only such a query, is asked of PubChem
PUG-REST and the NCI/CADD resolver (CACTUS):

```python
df = Search("cas", online_fallback=True).search(["50-78-2", "1912-24-9"])
df["source"]                  # "PubChem (online)" / "CACTUS" on rows the network supplied
df.attrs["online_fallbacks"]  # how many queries went online
df.attrs["online_resolved"]   # how many of those came back with an answer
```

Each service counts as one vote in `n_source_support`, the same as a database.
A service that fails is logged and left out, and formula queries are never
sent. How the requests are paced and retried is described in
[Network behaviour](network.md).

## What each row says

| Column | Description |
|---|---|
| `query` | the input value |
| `CASRN` | CAS Registry Number |
| `name`, `IUPAC_name` | preferred common name, IUPAC name |
| `molecular_formula`, `molecular_mass` | formula and molecular weight |
| `SMILES` | SMILES as the source gave it |
| `canonical_smiles`, `kekulized_smiles` | RDKit's canonical and Kekulé forms |
| `InChI`, `InChIKey` | InChI and full InChIKey |
| `DTXSID` | CompTox substance identifier |
| `Synonyms` | semicolon-separated synonyms |
| `parent_smiles`, `parent_inchikey` | the structure after salt stripping (`strip_salts=True`) |
| `foundby` | how the match was found |
| `source` | the source that provided the primary SMILES |
| `source_details` | per-source record of what each database returned |
| `confidence` | overall confidence in [0, 1], see below |
| `match_method` | the matching method used |
| `match_score` | cross-source consensus score in [0, 1] |
| `consensus_source` | the source the consensus chose |
| `source_match_scores` | per-source agreement scores |
| `hit_rank` | rank among the hits for this query (0 = best) |
| `n_source_support` | independent databases carrying this structure |
| `opsin_smiles` | OPSIN's parse of a name query (`use_opsin=True`) |

With `n_hits` above 1 (or `"all"`) a query can have several rows, one per
candidate structure, ranked by `hit_rank`. The frame's `attrs` record the
preset and settings it was made with, the sources that were and were not
available, and the online-fallback counts, so a saved result says how it was
produced.

## Confidence

Confidence combines four signals: how strong the match method is, how well the
candidate matches the query itself, how well the databases that answered agree
with each other, and how many of them carried the structure at all:

$$
\text{confidence} = \text{base}
\times (w_q \times \text{query\_score} + (1 - w_q))
\times (0.5 + 0.5 \times \text{consensus\_score})
\times \text{support}(n)
$$

| Match method | Base |
|---|---|
| Exact InChIKey | 1.00 |
| OPSIN parse of a name | 0.97 |
| Exact canonical SMILES | 0.95 |
| InChI | 0.95 |
| Exact CAS | 0.90 |
| DTXSID | 0.90 |
| Exact name | 0.80 |
| InChIKey skeleton | 0.75 |
| Tanimoto similarity | Tanimoto × 0.85 |
| Fuzzy name | rapidfuzz ratio × 0.80 |
| Formula | 0.30 |

For exact-identifier methods `query_score` is 1.0, which collapses the second
term. `w_q` is the `query_weight` argument (default 0.5).

`support(n)` scales the score by the number of *independent* databases carrying
the structure, `n_source_support`:

| Databases agreeing | Factor |
|---|---|
| 0 (an OPSIN-only parse) | 1.00 |
| 1 | 0.85 |
| 2 | 0.95 |
| 3 or more | 1.00 |

This term is what keeps provenance from beating evidence. `consensus_score`
measures how well the sources that answered agree — but a lone source agrees
with itself perfectly, so without `support(n)` a single database hit (0.90)
would outrank a structure that three databases all carry (0.8777).

Use `min_source_support=` (on the constructor or per call) to require
corroboration outright: `Search("cas", min_source_support=2)` returns only
structures that at least two databases agree on, which is the `"strict"`
preset.

Because corroboration drives confidence, a run with a database missing scores
lower than a full run on the same input. Scores are comparable only between
runs on the same set of sources, which `df.attrs["sources_available"]`
records; `Search(datasets="required")` refuses to run on fewer.

## Closing

`Search` is a context manager. It closes the database clients it constructed,
and leaves alone any you passed to its constructor — those are yours, and you
may still be using them:

```python
from provesid import Search, PubChemID

db = PubChemID()
with Search("cas", pubchem=db) as s:
    df = s.search("50-00-0")
db.get_by_cas("64-17-5")        # still open
```
