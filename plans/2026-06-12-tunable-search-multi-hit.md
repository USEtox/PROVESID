# Plan: Tunable, multi-hit `Search` with PYOPSIN structure anchoring

**Date:** 2026-06-12
**Author:** Ali A Eftekhari + Claude Code
**Target module:** `src/provesid/search.py` (and helpers in `src/provesid/tools.py`)
**Status:** Proposed

---

## 1. Problem statement

`Search.search()` returns exactly **one** merged result per query, and that result is
not always correct — most often for **name** queries. Two structural causes:

1. **Each source is collapsed to a single candidate before consensus runs.**
   In `_candidates_from_name` (and the other resolvers) we take `rows[0]` (CompTox/
   PubChem/ChEMBL) or the first synonym hit (ChEBI), so the consensus engine
   (`_compute_consensus`) only ever sees **one candidate per source**. If the correct
   compound is hit #2 in a source while a wrong compound is hit #1, the correct
   compound is invisible to scoring. A single wrong #1 hit can dominate the whole
   answer.

2. **The scoring/fuzzy machinery is hard-coded and inflexible.** rapidfuzz cut-off,
   scorer choice, the `0.7` ZeroPM magic constant, the consensus-compatibility
   threshold (`0.35`), and `_PRIORITY_TOLERANCE` (`0.05`) are all baked in. There is
   no way for the caller to tune precision/recall, and no way to use a structure
   oracle (PYOPSIN) to anchor an ambiguous name to a real molecule.

### Design decisions (confirmed with user, 2026-06-12)

- **Multi-hit output = distinct compounds, ranked.** When more than one hit is
  requested, each output row is one *distinct candidate compound* (deduped across
  sources by InChIKey skeleton / canonical structure), ranked by confidence, with a
  new `hit_rank` column. Default remains **one row per query**.
- **PYOPSIN is opt-in, default off** (`use_opsin=False`) because it needs a Java
  runtime. When enabled, IUPAC names are converted to SMILES offline and used as a
  high-confidence structure anchor.
- **Full root-cause refactor**: collect **top-K candidates per source**, pool them,
  cluster by structure, and let scoring choose among *all* candidates rather than one
  pre-chosen winner per source.

---

## 2. Current architecture (as-is)

```
search(queries)
  └─ _resolve_single(q)
       └─ _resolve_name(q)                    # (and _resolve_cas/_smiles/...)
            ├─ _candidates_from_name(q)        # {source: ONE candidate or None}
            │     chebi  -> rows[0]
            │     comptox-> get_by_name or search_by_name[0]
            │     pubchem-> search_by_name[0]
            │     zeropm -> id_table -> one candidate
            │     chembl -> search_by_name[0]
            ├─ (fuzzy fallback) _fuzzy_name_candidates(q)
            └─ _finalise_result(result, candidates, method, details)
                 ├─ _compute_consensus(candidates)   # picks 1 anchor source
                 ├─ apply compatible candidates -> ONE merged result dict
                 ├─ normalize_structure / salts / confidence
                 └─ returns ONE row
```

Key data structure today: `candidates: Dict[str, Optional[dict]]` — **one candidate
per source key**. This is the bottleneck to break.

Relevant source-method capabilities already available (currently under-used):
- `chebi.search_by_name(name, exact)`, `search_by_synonym(name, exact)` → `List[dict]`
- `comptox.search_by_name(name, exact, limit)` → `List[dict]`; `get_by_name` → one
- `pubchem.search_by_name(name, exact, limit)` → `List[dict]`
- `zeropm.query_similar_name(name, number_of_results=5, score_cutoff=80)` (fuzzy!)
- `chembl.search_by_name(name, limit)` → `List[dict]`

---

## 3. Target architecture (to-be)

Introduce a **candidate pool**: a flat list of candidate records (not one-per-source),
each tagged with its originating source and a per-source match score. Consensus and
ranking then operate over the pool and over **structure clusters**.

```
search(queries, n_hits=1, ...)
  └─ _resolve_single(q) -> List[dict]          # now returns a LIST of ranked hits
       └─ _resolve_name(q)
            ├─ pool = _candidate_pool_from_name(q)        # list[candidate], top-K/source
            │     + (opt) PYOPSIN anchor candidate (structure)
            ├─ clusters = _cluster_candidates(pool)       # group by structure identity
            ├─ ranked  = _rank_clusters(clusters, q)      # confidence per cluster
            └─ _finalise_hits(ranked, n_hits)             # 1..N merged rows + hit_rank
  └─ flatten list-of-lists -> DataFrame (n_hits column added)
```

### 3.1 Candidate record (extended)

Reuse the existing `_make_candidate(...)` shape, adding two transient fields used only
during ranking (dropped from final output unless surfaced):

- `query_match_score: float` — how well this candidate matches the *query itself*
  (e.g. name similarity for name queries, 1.0 for exact CAS/InChIKey, Tanimoto for
  structure fallback). This is the missing signal today.
- `_origin_rank: int` — the rank position the candidate had within its source's result
  list (0-based). Used as a tie-breaker and for diagnostics.

### 3.2 Structure clustering / dedup

`_cluster_candidates(pool)` groups candidates that refer to the **same compound**:

1. Primary key: full InChIKey when present.
2. Fallback key: 14-char InChIKey **skeleton** (connectivity) — merges
   stereo/charge/isotope variants when `cluster_by_skeleton=True` (default True;
   matches the spirit of `inchikey_skeleton`).
3. Fallback key: canonical SMILES.
4. Candidates with no structure at all form singleton clusters keyed by normalized
   name (so name-only sources still contribute, but score lower).

Each cluster collects: member candidates, the set of contributing sources, and the
best `query_match_score` among members.

### 3.3 Cluster ranking → confidence

Per cluster, compute a confidence in `[0, 1]` combining:

- **base** from `_BASE_CONFIDENCE[match_method]` (existing table),
- **query agreement** = `query_match_score` (name/Tanimoto/exact signal),
- **cross-source consensus** = number/agreement of distinct sources backing the
  cluster (re-using `_candidate_similarity` pairwise within and across clusters), and
- **PYOPSIN bonus**: when an OPSIN-derived structure anchor exists and a cluster
  matches it (by InChIKey or skeleton), that cluster gets a strong confidence boost —
  this is what fixes wrong name hits, because OPSIN gives a near-ground-truth structure
  for parseable IUPAC names.

Formula (generalizes the current `_compute_confidence`):

```
confidence = clamp(
    base
    * (w_q * query_match_score + (1 - w_q))      # query agreement term
    * (0.5 + 0.5 * consensus_score)              # existing consensus modulation
    + opsin_match_bonus                          # additive, capped
)
```

Weights (`w_q`, `opsin_match_bonus`, etc.) become tunable attributes (see §4).

The **best cluster** (rank 0) is the new "single best match" — equivalent to today's
output when `n_hits=1`, but chosen from the *full* pool rather than from pre-collapsed
per-source winners.

---

## 4. New `Search` attributes (constructor)

All keyword-only, all with defaults preserving current behavior.

| Attribute | Type | Default | Purpose |
|---|---|---|---|
| `n_hits` | `int \| "all"` | `1` | Default number of ranked hits to return per query. Can be overridden per-call in `search()`. |
| `min_confidence` | `float` | `0.0` | Drop hits below this confidence before truncating to `n_hits`. |
| `use_opsin` | `bool` | `False` | Enable PYOPSIN IUPAC→structure anchoring for name queries (needs Java). |
| `opsin_jar_fpath` | `str` | `"default"` | Passed to `PYOPSIN`. |
| `top_k_per_source` | `int` | `5` | How many candidates to pull from each source before pooling. |
| `cluster_by_skeleton` | `bool` | `True` | Merge stereo/charge variants when clustering. |
| `fuzzy_score_cutoff` | `float` | `80.0` | rapidfuzz / ZeroPM score cut-off (0–100). Replaces hard-coded `80`/`0.7`. |
| `fuzzy_scorer` | `str` | `"WRatio"` | rapidfuzz scorer name (`WRatio`,`token_sort_ratio`,`partial_ratio`,…). |
| `consensus_compat_threshold` | `float` | `0.35` | Was hard-coded in `_candidate_compatible_with_consensus`. |
| `query_weight` (`w_q`) | `float` | `0.5` | Weight of query-agreement vs base in confidence. |
| `return_alternatives` | `bool` | `False` | When `n_hits=1`, still attach runner-up summaries in a new `alternatives` column (for inspection without exploding rows). |

`search()` gains per-call overrides:

```python
def search(self, queries, *, column=None, n_hits=None, min_confidence=None): ...
```
(falling back to the instance attributes when `None`).

Validation: `n_hits` must be `>= 1` or the literal `"all"`; raise `ValueError`
otherwise. `fuzzy_scorer` validated against a whitelist mapped to rapidfuzz functions.

---

## 5. Output schema changes

Add three columns to `OUTPUT_COLUMNS`:

- `hit_rank` (int, 0 = best) — position of this row among the hits for its query.
- `n_source_support` (int) — how many distinct databases back this hit's cluster.
- `opsin_smiles` (str, nullable) — SMILES OPSIN produced for the query name (only
  populated when `use_opsin=True` and parse succeeded), for traceability.

Optional column (only when `return_alternatives=True` and `n_hits==1`):
- `alternatives` (list[dict]) — compact `{name, InChIKey, confidence, source}` for the
  next few runner-up clusters.

**Backward compatibility:** with defaults (`n_hits=1`), output is one row per query as
today, plus the three new always-present columns. Existing column order preserved;
new columns appended. Downstream code selecting by name keeps working.

---

## 6. PYOPSIN integration

- Add lazy `self._opsin` (a `PYOPSIN` instance) created on first name query when
  `use_opsin=True`. Wrap construction/calls in try/except; on any failure (no Java,
  py2opsin import error) log a warning **once** and disable OPSIN for the session
  (`self._opsin_available = False`) so we don't spam.
- In `_candidate_pool_from_name`:
  1. `smiles = self._opsin.get_smiles(name)` (offline, fast).
  2. If non-empty, run it through `normalize_structure` to get canonical SMILES +
     InChIKey, build an **anchor candidate** (`source="OPSIN"`, `query_match_score=1.0`,
     `match_method="opsin"`), and add it to the pool.
  3. Additionally use the OPSIN InChIKey to do a direct **structure lookup** across
     sources (`_candidate_pool_from_inchikey`-style), pulling in the *correct* compound
     even when its name spelling differs — this is the main accuracy win.
- Add `"opsin": 0.97` to `_BASE_CONFIDENCE` (just below exact InChIKey).
- OPSIN also usable for `smiles`/`inchi` cross-checks later, but scope here = names.

---

## 7. Implementation steps

Ordered, each independently testable.

1. **Introduce candidate-pool helpers in `tools.py`** (new, additive — don't break
   existing `_compute_consensus` callers in `tools.py`'s `ids_from_*`):
   - `_cluster_candidates(pool, by_skeleton) -> List[Cluster]`
   - `_cluster_confidence(cluster, *, base, weights, opsin_anchor) -> float`
   - `_rank_clusters(clusters, query, ...) -> List[(confidence, cluster)]`
   - Keep these pure/stateless so they're unit-testable without DB clients.

2. **Refactor name resolution to pool-based** in `search.py`:
   - New `_candidate_pool_from_name(name) -> List[candidate]` that pulls `top_k_per_source`
     from each source (using existing `limit`/`exact` params), tags each with
     `query_match_score = _text_similarity(name, cand.name/synonyms)` and `_origin_rank`.
   - Fold the existing fuzzy logic in: when exact pool is weak, widen with
     `exact=False` + `fuzzy_score_cutoff`/`fuzzy_scorer` and ZeroPM
     `query_similar_name(number_of_results=top_k_per_source, score_cutoff=fuzzy_score_cutoff)`.
   - Add OPSIN anchor (§6).

3. **Generalize `_finalise_result` → `_finalise_hits`**:
   - Build merged result dict per cluster (reuse `_apply_candidate_to_result` over the
     cluster's members, anchored on the cluster's best-structure member).
   - Run `normalize_structure`, salt stripping, confidence per cluster.
   - Set `hit_rank`, `n_source_support`, `opsin_smiles`.
   - Return `List[dict]` (ranked, filtered by `min_confidence`, truncated to `n_hits`).
   - Keep a thin `_finalise_result` wrapper returning `[0]` for resolvers not yet
     migrated, so we can land name first.

4. **Update `_resolve_single` + `search()`** to handle list-returning resolvers:
   - `_resolve_single` returns `List[dict]`; `search()` flattens with `tqdm` over
     queries and concatenates. Each query's hits stay contiguous and ordered.
   - Apply per-call `n_hits`/`min_confidence` overrides; validate.
   - DataFrame assembly: extra-column merge (from DataFrame/file input) must now
     **broadcast** original row data across the multiple hit-rows of that query — join
     on a per-query index rather than positional `reset_index`. (Important edge case:
     today's positional concat assumes 1 row per query; multi-hit breaks that.)

5. **Migrate remaining resolvers** to the pool model for consistency (CAS, SMILES,
   InChI, InChIKey, DTXSID, formula). For exact-identifier types the pool is usually
   a single cluster, so behavior is unchanged at `n_hits=1`; but `formula` (inherently
   many hits) and `smiles` Tanimoto fallback benefit immediately from ranked multi-hit.

6. **Constructor + validation** for all new attributes (§4); update class docstring and
   the module-level usage examples.

7. **Confidence table + weights**: add `opsin` base; expose weights; keep
   `_compute_confidence` working for the single-cluster path.

---

## 8. Testing plan

Add to the existing test suite (mirror current test layout; check `tests/`).

- **Unit (no DB, no Java):**
  - `_cluster_candidates`: same-InChIKey merge; skeleton merge on/off; SMILES-only and
    name-only fallbacks; empty pool.
  - `_cluster_confidence` / ranking: ordering monotonic in query score and source
    support; OPSIN anchor pushes the matching cluster to rank 0.
  - `n_hits` validation; `"all"`; `min_confidence` filtering.
- **Integration (DB clients, marked like existing offline tests):**
  - A known ambiguous name where today's code picks the wrong compound (capture a
    real example from the user — see §10) → assert correct compound is rank 0 *or*
    appears within top-`n_hits`.
  - `n_hits="all"` returns >1 row for an ambiguous name; `hit_rank` contiguous 0..N.
  - DataFrame input: extra columns correctly broadcast across multi-hit rows.
- **OPSIN (marked, skip if Java/py2opsin missing — reuse the existing `opsin` pytest
  marker):** IUPAC name → correct InChIKey anchor; graceful no-op when OPSIN
  unavailable (no exception, falls back to name matching).
- **Regression:** default `Search("name").search([...])` still yields exactly one row
  per query; existing columns unchanged; `Search("cas")` outputs identical to before.

---

## 9. Risks & mitigations

| Risk | Mitigation |
|---|---|
| Multi-hit rows break downstream code expecting 1 row/query | Default `n_hits=1` unchanged; new columns appended, not inserted. Document the broadcast behavior for DataFrame inputs. |
| Pulling `top_k_per_source` from 5 sources slows batch runs | `top_k_per_source` default modest (5); clustering is O(pool²) but pool is tiny per query. Profile on a batch; keep `show_progress`. |
| OPSIN/Java absent in CI or user env | Opt-in (`use_opsin=False`); one-time warning + session disable; tests skip via marker. |
| Skeleton clustering over-merges distinct stereoisomers | `cluster_by_skeleton=True` default but configurable; `n_source_support` + full InChIKey preserved in output for the user to disambiguate. |
| Tunable weights make confidence non-comparable across configs | Document that `confidence` is only comparable within a fixed configuration; keep default weights == current behavior at `n_hits=1` as far as possible. |

---

## 10. Open questions for the user (non-blocking)

1. **A concrete failing example.** Please provide 1–3 name queries where the current
   `search` picks the *wrong* compound (and what the right answer is). These become
   regression tests and let us calibrate default weights against reality.
2. **Default `top_k_per_source`** — is 5 reasonable, or do you want it higher for
   recall (at some speed cost)?
3. **`alternatives` column** — useful, or noise? (Off by default either way.)
4. Should the pool refactor also expose a **raw, unranked** debug accessor (e.g.
   `Search.debug_pool(query)` returning every candidate + scores) to help you tune
   weights interactively? Recommended; cheap to add.

---

## 11. Out of scope (future work)

- Vectorized Parquet + fingerprint Tanimoto search (already noted as future in
  `_tanimoto_candidates`).
- Learning/auto-tuning the confidence weights from a labeled set.
- Applying OPSIN to non-name identifier types as a cross-validation oracle.
