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
0     50-78-2  acetylsalicylic acid  BSYNRYMUTXBXSQ-UHFFFAOYSA-N    ChEBI                 3       0.880
1  39400-72-1              Atrazine  MXWJVTOOROXGIU-UHFFFAOYSA-N  CompTox                 1       0.765
2   1912-24-9              atrazine  MXWJVTOOROXGIU-UHFFFAOYSA-N    ChEBI                 3       0.900
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

**ZeroPM is left out** unless `sources` names it. It aggregates regulatory
inventories rather than curating compounds, so its name→structure rows are
noisier than the others, and as a full vote in the corroboration count they
used to push wrong structures up the ranking. It is also the one source that
retrieves misspelled names, which is why the `"recall"` preset turns it on.

`sources` chooses the databases: a list of keys, one key, or `"all"`.
A database left out is never opened.

```python
Search("cas", sources="all")                      # the four, plus ZeroPM
Search("cas", sources=["pubchem", "comptox"])     # these two and no others
```

What each source can be asked, and by which identifier, is one table,
`provesid.sources.LOOKUPS`. A gap in it is a source with no index for that
identifier, which `Search` reaches through one it does have — ChEMBL through
the InChIKeys the other sources found for a CAS number, for instance. A
record found that way is an echo of the answer that found it, not a second
opinion: it fills fields the others left empty, and is shown in
`source_details`, but it does not vote and does not count in
`n_source_support`.

A DTXSID is CompTox's own identifier, so for a DTXSID query the other
sources are asked for the CAS number CompTox gives the substance. Each that
holds it answers from its own record and votes, as for a CAS query; one that
does not is reached through the InChIKeys found, as an echo. A one-against-one
split goes to the CompTox record the DTXSID names, with `conflict` set.

## Presets

The arguments that decide what counts as a match and what is returned have
three named settings, `Search.PRESETS`:

| preset | what changes from `balanced` | for |
|---|---|---|
| `"balanced"` | nothing: the constructor's defaults | general use |
| `"strict"` | `min_source_support=2`, `on_conflict="refuse"` | tables where a wrong structure costs more than a missing one |
| `"recall"` | `fuzzy`, `inchikey_skeleton`, `similarity_threshold=0.7`, `sources="all"`, `n_hits="all"` | finding candidates to review by hand |

```python
Search.PRESETS["strict"]                          # a plain dict: inspect it
Search("name", preset="strict").search("atrazine")
Search("name", preset="recall", n_hits=3)         # explicit arguments win
df.attrs["preset"], df.attrs["settings"]          # what produced this frame
```

An argument passed explicitly overrides the preset, even when it equals the
`balanced` value: `preset="recall", sources=["chebi", "comptox", "pubchem",
"chembl"]` leaves ZeroPM out.

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
| `CASRN` | the compound's CAS Registry Number, the current one where the sources say which; a search by a retired number reports the current one, and `query` keeps the number searched |
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
| `structure_status` | `valid`, `repaired:molfile` / `repaired:inchi` / `repaired:nh_repair`, `invalid` or `none`; see [Structure quality](#structure-quality) |
| `identity_consistent` | whether the InChIKey and formula the sources state belong to the structure |
| `source_inchikey`, `source_formula` | the InChIKey and formula as the sources state them |
| `n_clusters` | how many distinct compounds the sources returned for the query |
| `conflict` | an independent source's best answer is a different compound; see [Conflicts](#conflicts) |
| `runner_up_inchikey`, `runner_up_sources`, `runner_up_confidence` | the compound this row competed with |

With `n_hits` above 1 (or `"all"`) a query can have several rows, one per
candidate structure, ranked by `hit_rank`. The frame's `attrs` record the
preset and settings it was made with, the sources that were and were not
available, the dataset files they read (`dataset_versions`: path, size,
modification time and release of each, ChEBI's index format, and the provesid
version), and the online-fallback counts, so a saved result says how it was
produced.

## Structure quality

A row is only as useful as its structure, and three things can be wrong with
one that look fine at a glance. `Search` repairs what it can and says what it
could not.

**A source's SMILES does not parse, or is not the compound it names.** 6% of
ChEBI's exported SMILES drop the hydrogen of an aromatic ring nitrogen
(valsartan's tetrazole reads `c2nnnn2`), and RDKit cannot read them; another
10% parse but lost a double bond's configuration, or, like guanine's
`Nc1nc(=O)c2ncnc2n1`, two hydrogens. Each candidate's structure is rebuilt
from the source's own molfile, its InChI, or by putting the ring N–H back,
and the rebuild is kept only when its InChIKey is the one the source states.
`structure_status` says which: `repaired:molfile` for valsartan. A SMILES that
cannot be rebuilt is dropped rather than passed on, so another source's
structure can take its place; when none can, the row is `invalid`, has no
SMILES, scores lower, and is logged at WARNING.

**The identifiers belong to different molecules.** When the row has a
structure, `InChIKey`, `InChI`, `molecular_formula` and `molecular_mass` are
computed from it, so they always describe the `SMILES` of the row. What the
sources stated is kept in `source_inchikey` and `source_formula`, and
`identity_consistent` is False — and a WARNING logged — when they do not match
the structure in skeleton, charge or element counts (stereo is not compared).

**No structure at all.** CompTox lists some substances (eprinomectin,
gentamicin, neomycin) as mixtures with a DTXSID and no structure. Such a record
joins the structured cluster of the same substance when exactly one shares its
DTXSID, name or CAS number, adding its DTXSID and its support. A row that ends
up with no structure (`structure_status` `none` or `invalid`) has its
confidence multiplied by 0.8, so a structured answer wins a tie.

## Conflicts

When independent databases return different compounds for one query, the
answer is a judgement, and the row says so:

```python
df = Search("cas", n_hits=1).search("118443-89-3")   # cefquinome sulfate
df[["name", "source", "conflict", "n_clusters", "runner_up_sources", "confidence"]]
```

```text
                 name   source  conflict  n_clusters runner_up_sources  confidence
0  Cefquinome sulfate  CompTox      True           2       [PubChemID]      0.3825
```

CompTox's record for this CAS number is a 1,2,4-thiadiazole; PubChem's is
cefquinome. A *conflict* is a different compound as the best answer of an
independent database — not a salt, solvate, charge state or stereoisomer of
this one (they compare by the InChIKey skeleton of their largest fragment,
neutralised and in RDKit's canonical tautomer), not
a lower-ranked alternative, and not a fuzzy, Tanimoto or formula match, whose
results are alternatives by design. The confidence of a contested row is
multiplied by its share of the independent databases behind the two readings:
a one-against-one split halves it, two against one keeps two thirds.

`on_conflict="refuse"`, part of the `"strict"` preset, returns a no-match row
with `conflict=True` instead of choosing. That costs answers the majority gets
right — tirapazamine's CAS number, where ChEBI and PubChem agree and CompTox
has the wrong formula — in exchange for never returning a disputed structure.

## Validating hits against your own data

Most measurement tables already carry a molecular weight, a formula or a
structure, and checking each hit against it is the cheapest way to catch a
confident wrong answer. `validate=` names those columns; a hit that fails is
passed over for the next-ranked compound, and the checks each returned hit
passed are in `validated_by`:

```python
table = pd.DataFrame({"CAS": ["27314-97-2"], "MW": [178.15]})
s = Search("cas", validate={"mw": "MW", "mw_tolerance": 0.5})
s.search(table, column="CAS")[["name", "molecular_formula", "validated_by"]]
```

```text
           name molecular_formula validated_by
0  tirapazamine          C7H6N4O2           mw
```

The keys are `"mw"` (a molecular weight), `"smiles"` (a structure, compared
by weight; an identical one is also reported as `"smiles"`), `"formula"`
(element counts must match) and `"name"` (reported, never required). An empty
reference cell skips its check. This is
[`resolve_cascade`](../api/search.md)'s `accept=` check run inside a single
search, so it also works where `Search` is called directly.

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
structures that at least two databases agree on, which is half of the
`"strict"` preset.

Two further factors apply after this one: a contested answer is multiplied by
its share of the independent databases (see [Conflicts](#conflicts)), and an
answer without a structure by 0.8.

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

Passing one client does not stop `Search` from building the others: here
ChEBI, CompTox and ChEMBL are opened as usual on the first search, and
closed at the end of the `with` block. To query only the client you pass,
add `sources="pubchem"`.
