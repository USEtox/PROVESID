# Status of the refactor, and what to do next

**Date:** 2026-09-20
**Author:** Ali A Eftekhari + Claude Code
**Scope:** all of `src/provesid/`, `docs/`, `examples/`, `scripts/`, `README.md`
**Companion:** `plans/2026-08-02-package-wide-refactor.md` (§1–23), which this
document assesses rather than replaces.
**Status:** Proposed.

---

## 1. Method

Read every module in `src/provesid/` (20 407 lines), re-measured docstring
coverage with an AST pass, collected the test suite, diffed `examples/` against
`docs/examples/`, and compared the shipped `pubchem_id.db` row-by-row against
PubChem's own CAS↔CID mapping. Probed the PubChem FTP service directly on
2026-09-20 for §8: real file sizes, real row counts, range-request support and
measured throughput from this machine. Every number below was measured, not
estimated, unless it says otherwise.

Decisions confirmed with the user before writing §8 are recorded in §8.2.

---

## 2. Where the refactor stands

§12 of the August plan sequenced seventeen steps. Seven have landed and one is
partly done; all but the three from 2026-08-02 landed in the last two days.

| # | Step | Status |
|---:|---|---|
| 1 | repo hygiene, `pyproject.toml` metadata | **done** (§21) |
| 2 | delete legacy resolvers from `tools.py` | **done** (§21) |
| 3 | `http.py`; migrate `resolver.py`, `pubchemview.py` | **done** (§22) |
| 4 | migrate `pubchem.py`, `chebi.py`, `cascommonchem.py`, `opsin.py` | **done** (§23) |
| 5 | `classyfire.py` raises; rewrite its tests | postponed (§17.0) |
| 6 | extract `pubchemview_parse.py` | **done** (§19, ahead of sequence) |
| 7 | split `pubchem.py` → `+ pubchem_id.py`; `chebi.py` → `+ chebi_sdf.py` | not started |
| 8 | `sources.py`; collapse the `Search` source ladders | not started |
| 9 | fix the fuzzy-name mislabelling | **done** (§16.1) |
| 10 | `cache.py` parameterisation and `_MISS` sentinel | **done** (§17 below) |
| 11 | `zeropm.py` logging, `reach.py` xlsx, `config.py` printing | not started |
| 12 | `enrich`, `resolve_cascade`, `mw_within` | **done** (§16.2) |
| 13 | docstrings with examples, module by module | not started |
| 14 | rebuild `docs/` from docstrings | not started |
| 15 | the 17 notebooks | not started |
| 16 | rewrite `README.md` | not started |
| 17 | `CHANGELOG.md` for the release | not started |

Four unplanned pieces of work also landed and are worth keeping in view because
later steps must not undo them: the cache-correctness fixes (§17), bulk and
offline-first properties (§18), the single typed PUG-View parser (§19), and
`Search` dropping ZeroPM from the default vote (§20).

### 2.1 Measured state, 2026-09-20

| | 2026-08-02 | today |
|---|---:|---:|
| `src/provesid/` lines | ~16 700 | 20 407 |
| Public objects | 415 | 483 |
| …with no docstring | 32 | **34** |
| …with no usage example | 305 | **298** |
| Copies of the HTTP transport | 4 (+3 with none) | **1** |
| Modules handling HTTP 429 | 0 | **all but `classyfire.py`** |
| Copies of the *bulk download* layer | 5 | **5** |
| Tests collected | — | 1035 |
| Notebooks | 0 | **0** |
| Hand-written `docs/api/` lines | ~3 800 | **3 570** |

The package grew by 3 700 lines while the transport shrank to one copy, which is
the right trade: the growth is `http.py` (849), `pubchemview_parse.py` (644) and
the property and enrichment layers, all of which replaced duplicated or absent
behaviour. Docstring coverage is essentially unmoved because step 13 has not
started — the 68 new public objects arrived documented (`http.py` is 19/19 with
examples), which is why the "no docstring" count rose by only two.

---

## 3. What the refactor got right, and should be copied

`http.py` is now the house standard and the rest of this document is largely
"apply it elsewhere". Specifically:

- **One policy object, many contracts.** The transport owns *when to ask again*;
  each client keeps its URLs, its parsing and its exception type. The `classify`
  callback is the one hook, and it exists because PubChem genuinely lies about
  its status codes (§22.1, §23.4).
- **The pacer belongs to the host, not the object.** §23.2 moved the rate limiter
  off the client instance onto a per-host shared clock, because PubChem's limit
  is per IP. Two `PubChemView` instances no longer exceed it together.
- **Evidence in the docstring.** Every non-obvious decision cites the observation
  that forced it. `ServiceError`'s docstring explains why `status_code` is
  keyword-only. This is the standard step 13 should hold every module to.
- **Failures are not answers.** §17.2 and §23.6 stopped the caches storing
  failures; §17.4 gave them a real miss sentinel.

Two more strengths worth stating plainly, because the improvements below must
not damage them: `Search` is genuinely offline — it opens no socket on any code
path — and `PubChemID.properties()` is the one place where dev-principle §9's
offline→online fallback is actually implemented, with a `Source` column naming
which one answered.

---

## 4. Priority improvements

Ordered by what a user hits first. Sizes are S/M/L in the August plan's sense.

### 4.1 First run downloads 32 GB without asking (L)

This is the largest user-friendliness defect in the package, and it is not in the
August plan at all.

```python
from provesid import Search
df = Search("cas").search("50-00-0")
```

On a clean machine `_ensure_clients` (`search.py:651`) constructs `ChebiSDF`,
`CompToxID`, `PubChemID` and `CheMBL`, each with `auto_download=True`, in a loop.
Measured from the copies on this machine:

| Source | Resident | Notes |
|---|---:|---|
| `chebi.sdf` + index | 943 MB | index built on first use |
| `comptox_chemicals.db` | 817 MB | |
| `pubchem_id.db` | 2.2 GB | |
| `chembl_36.db` | **28 GB** | from a 5.8 GB archive, extracted |
| **total** | **~32 GB** | for one CAS lookup |

Nothing announces the total before it starts, nothing asks, and nothing offers to
proceed with the sources already present. ChEMBL alone is 87% of it, and in
`Search` it only *enriches* — it is queried last, by SMILES, after the primary
sources have already produced a structure (`search.py:1361`).

The user constraint is explicit: this package is for researchers on laptops,
where 32 GB is often the whole free disk.

**Proposal — a dataset manager, `src/provesid/datasets.py`.**

```python
from provesid import datasets

datasets.status()          # table: name, present?, bytes on disk, release, path
datasets.plan(["pubchem", "chebi"])   # what would be downloaded, and how big
datasets.fetch("pubchem")  # explicit, resumable, checksummed
datasets.remove("chembl")  # reclaim the space, by name
```

and a `Search(..., datasets="present")` policy argument:

| value | behaviour |
|---|---|
| `"present"` | use whatever is on disk; log which sources are missing. **Proposed default.** |
| `"auto"` | today's behaviour — download whatever is missing |
| `"required"` | raise, naming the missing datasets and the exact `datasets.fetch` call |

`Search` already degrades gracefully when a source is unavailable and already
warns that confidence is not comparable across runs (`search.py:693`), so
`"present"` needs no new machinery in the resolver — only that the clients stop
downloading behind the caller's back. Changing the default is a breaking change,
which dev-principle §1 permits.

`ChEMBL` deserves a second look on its own terms: 28 GB for an enrichment-only
source is a bad bargain. **§9 answers it** — the package uses 8 of ChEMBL's 74
tables, and an extract holding just those is 2.60 GB with byte-identical
results. That changes the arithmetic above from ~32 GB to ~6.5 GB and is the
single highest-value item in this document.

### 4.2 Five copies of the bulk-download layer, none resumable (M)

§23.8 deliberately left these out of the `http.py` migration, correctly: a 100 MB
file with a progress bar wants resumption and checksums, not a 5-per-second
pacer. But the consequence is that the transport that handles the *largest*
transfers in the package is the one with no shared implementation.

| Module | Call site | Retry | Resume | Checksum |
|---|---|:-:|:-:|:-:|
| `pubchem.py` | `download_database`, L2198 | no | no | no |
| `comptox.py` | `download_database`, L168 | no | no | no |
| `zeropm.py` | `download_database`, L142 | no | no | no |
| `chembl.py` | `download_database`, L357 | no | no | no |
| `chebi.py` | `download_sdf`, L1093 | no | no | no |

An interrupted 5.8 GB ChEMBL download starts again from zero — and so does an
interrupted 2.2 GB PubChem one. Corruption is caught unevenly: the two gzipped
downloads (`chembl.py`'s `.tar.gz`, `chebi.py`'s `.sdf.gz`) get truncation
detection for free from gzip's CRC and length trailer, but the three plain
SQLite downloads (`pubchem.py`, `comptox.py`, `zeropm.py`) have none. All three
run a connect-and-query check afterwards, which catches a grossly broken file;
a file truncated in its interior opens cleanly and fails much later, on the
first query that touches a missing page — which the user will read as a data
problem, not a download problem.

**Proposal — `download_file()` in the new `datasets.py`:** streamed, `Range`-resumed
against a `.part` file, checksum-verified when the server publishes one, atomic
rename on success, one progress bar, one logger. PubChem's FTP mirror publishes
an `.md5` beside every file (verified: `CID-SMILES.gz.md5` →
`3659dd5c96fc506bb11b7fd8cce4553d`) and answers `Range` with HTTP 206 (verified),
so both features are real for the largest dataset in §8. Five call sites collapse
to five one-line calls.

### 4.3 Four SQLite clients that cannot be closed (S) — **done, §16**

`PubChemID`, `CompToxID`, `ZeroPM` and `CheMBL` each open a connection in
`__init__` and close it in `__del__`. None has a public `close()`, none supports
`with`, and none passes `check_same_thread=False`.

The consequences are ordinary but real: a notebook that rebinds `db = PubChemID()`
leaks until the collector runs; on Windows the file stays locked, so
`redownload=True` on a live object cannot replace the file it is holding open;
and a connection created in one thread raises `ProgrammingError` if used from
another, which a user who reaches for `ThreadPoolExecutor` over a CAS list will
hit immediately.

**Proposal:** give all four a `close()`, `__enter__`/`__exit__`, and either
`check_same_thread=False` with a documented "read-only, safe to share" contract
or a thread-local connection. Keep `__del__` as a backstop. `Search` closes the
clients it constructed. Size S, and it removes a class of report that is annoying
to diagnose from a bug report.

### 4.4 Offline→online fallback exists in exactly one place (M) — **done, §23**

Dev-principle §9 asks for a two-stage lookup with an opt-out parameter,
documented, logged at DEBUG. `PubChemID.properties()` does this properly
(`pubchem.py:2961`). Nothing else does.

In particular `Search` — the package's flagship entry point, and the thing a user
reaches for first — has no online path at all. When all four offline sources miss
a CAS number, the user gets an empty row, even though `PubChemAPI`,
`NCIChemicalIdentifierResolver` and `OPSIN` are in the same package and would
very often answer. That is the wrong shape for "offline first, online offered
too": offline first is a *performance and traffic* decision, and it should not
cost the answer.

**Proposal:** `Search(..., online_fallback=False)`, off by default so the
offline-only guarantee stays the default and batch runs stay reproducible. When
enabled, queries that produced no candidate from any offline source — and only
those — are retried against PUG-REST and CACTUS, the result is tagged
`source="pubchem-online"` in `source_details`, and `foundby` records it. The
existing `sources_available` machinery already gives the honest reporting this
needs. Log each fallback at DEBUG, and count them in `df.attrs` so a user can see
how much network a batch actually used.

### 4.5 The `Search` source ladders (§5 B.2) (M) — **done, §22**

Unchanged since August, and now the largest remaining piece of duplication in the
package. `search.py` holds seven `_resolve_<type>` methods, each a hand-written
ladder of four or five blocks of exactly this shape:

```python
if self._chebi is not None:
    try:
        rows = self._chebi.search_by_cas(cas)
        if rows:
            candidates["chebi"] = candidate_from_chebi_row(rows[0])
    except Exception as exc:
        log.warning("ChEBI CAS lookup failed for %r: %s", cas, exc)
```

Roughly thirty copies differing only in the client, the method name, the adapter
and the log string. Adding a source means editing seven methods; forgetting one
is silent. §5 B.2's `sources.py` — a table of small functions, one row per
(source, identifier type) — remains the right answer, and §20 already did the
hard part by settling which sources deserve a vote.

This is also the natural home for §4.4: an online rung is one more row in the
table rather than seven more blocks.

### 4.6 `Search.__init__` has 22 keyword arguments (S)

They are all documented, and the docstring is genuinely good — it explains why
`WRatio` is a trap and what `min_source_support` costs. But 22 knobs is not a
user-friendly surface, and the defaults encode a policy the user cannot name.

**Proposal:** keep every argument, add named presets that a user can read and
cite in a methods section:

```python
Search("name", preset="strict")    # min_source_support=2, no fuzzy, no skeleton
Search("name", preset="recall")    # fuzzy, skeleton, similarity, n_hits="all"
Search("cas")                      # preset="balanced", today's defaults
```

A preset is a dict of the same keyword arguments, exposed as
`Search.PRESETS` so it is inspectable, with explicit arguments overriding it.
Size S, and it makes the interesting question ("what settings produced this
table?") answerable in one word.

### 4.7 `cache.py` (S, partly §5 B.7)

Three things, in order of user impact:

1. **The cache does not survive a restart.** The default directory is
   `tempfile.gettempdir()/provesid_cache` (`cache.py:174`), while
   `docs/advanced_caching.md` promises that it does. On most Linux systems
   `/tmp` is cleared on boot. `platformdirs.user_cache_dir` is already an
   effective dependency through `utils.user_dataset_path`; use it.
2. **Fourteen near-identical module-level functions.** `clear_pubchem_cache`,
   `clear_cas_cache`, … and seven `get_<svc>_cache_info` twins
   (`cache.py:613–667`). `clear_cache(service="pubchem")` and
   `get_cache_info(service="pubchem")` say the same thing, and
   `get_all_service_cache_info` already proves the service list is data.
3. **No cache version in the key.** §19.6 found this the hard way when a cached
   return *shape* changed. The fix landed for that one case; the general
   mechanism — a version component in every key, bumped when a return shape
   changes — is still missing.

### 4.8 `classyfire.py` is the last module on raw `requests` (S, step 5)

Unchanged since the August finding that the service has classified nothing since
February 2023. It still has no rate limiting, no retry, no 429 handling, a class
docstring that tells users to `from data_extract import ClassyFireAPI`, and
seventeen tests that pass by asserting on HTTP status codes rather than on a
completed classification. It is still exported from `__init__.py` and still
listed in `docs/index.md` under "Online interfaces" with no warning.

Step 5 as written is right: every method raises `ServiceUnavailableError` with
the February 2023 evidence and a pointer to `ChebifierClassifier`; the URL
construction stays, commented, in case the service returns; the tests assert the
raise. The one thing to add is that `docs/` and `README.md` must stop
advertising it in the same breath as the working services.

### 4.9 Documentation is a second, drifting copy of the code (M, step 14) — **done, §29**

- `docs/api/` is **3 570 hand-written lines**. Six of its eleven pages use
  mkdocstrings at all; the other five — `chebi.md` (360), `cascommonchem.md`
  (675), `classyfire.md` (405), `opsin.md` (214) and `index.md` (322), 1 976
  lines between them — are pure hand-written duplication of docstrings that
  already exist.
- It has already drifted: `docs/api/index.md` documents
  `PubChemPUGViewAPI`, a name the package has not exported for some time.
- **`docs/examples/` is a byte-identical copy of `examples/`** — 19 files, all
  `SAME` on diff today. Two copies of every tutorial, one of which will rot.
- The nav has no page for `CompToxID`, `ZeroPM`, `PubChemID`, `taxonomy`,
  `cache`, `config`, `tools` or `reach` — that is, for most of the offline half
  of the package, which is the half the project is built around.
- The nav ships a **"Modernization"** section pointing at `docs/plans/`, which
  contains internal working documents (`PLAN_search.md` still describes the
  `tools.py` resolvers deleted in §21.4). Internal plans do not belong in user
  documentation.

Step 14 stands as written. Add to it: delete `docs/examples/` and point mkdocs at
`examples/` once, drop `docs/plans/` from the nav, and add the missing offline
pages.

### 4.10 The README sells the wrong package (S, step 16)

The first paragraph describes PROVESID as "Pythonic access to online services",
with offline mentioned last as an aspiration ("also aims to provide an offline
platform when data files are available"). That is the opposite of the design
decision. Concretely:

- **`Search` is not mentioned anywhere in the README.** The package's primary
  entry point is invisible to a new user.
- Three of the tutorial links are dead: `ChEBI_tutorial.ipynb`,
  `zeropm-example.ipynb`, `classyfire_tutorial.ipynb` — the files are `.md`.
- The badge says "Python 3.8+"; `pyproject.toml` says `>=3.12`.
- ClassyFire is listed among the offered services (§4.8).
- The 32 GB of §4.1 appears only as a parenthetical about `uv`.

### 4.11 Notebooks (L, step 15)

Still zero `.ipynb` in the repo. `mkdocs-jupyter` is configured and waiting,
`scripts/convert_notebooks_to_myst.sh` and `generate_notebooks_from_myst.sh`
exist, `docs/quickstart.md` is already a MyST notebook with a
`validate_docs_local.sh --execute` path. The infrastructure is in place; the
content is not. Dev-principle §7 asks for a notebook per feature area.

This is also where §21.5's gap sits: `Search` has no notebook at all, because
B.1 deleted the old `notebooks.py` on the understanding that step 15 would
replace it.

### 4.12 Smaller items, unchanged from the August plan

- ~~**Step 7**, the `pubchem.py` / `chebi.py` splits.~~ **Done, §21.** `pubchem.py` is 3 301 lines
  holding two unrelated things: an online PUG-REST client and an offline SQLite
  database. §8 below makes the split more attractive, since the FTP builder
  belongs beside `PubChemID`, not beside `PubChemAPI`.
- **Step 11**: `zeropm.py` logging, `reach.py` xlsx, `config.py` printing.
- ~~**A circuit breaker** (§23.8).~~ **Done, §27.** A `Retry-After` is information about the *host*,
  so it belongs on the host's clock as a "not before T" that every client
  respects, making a known-throttled host fail at once. The shared `RateLimiter`
  is already the right place.
- ~~**34 public objects still have no docstring**~~ **Done, §28.** Concentrated in `pubchem.py`
  (10), `opsin.py` (7), `search.py` (6) and `cache.py` (2).
- **Stray files in `src/provesid/data/`**: `build_pubchem_id_db.py` (a duplicate
  of `scripts/build_pubchem_id_db.py`), `examining_pubchem.py`,
  `recreate_tables.sql`, `schema_documentation.txt`. A package data directory
  should hold data.
- **`scripts/README.md` is wrong and Windows-specific**: it says a local build is
  "required before using `PubChemID`", which auto-download made untrue, and it
  documents `cd c:\projects\git\PROVESID`.

### 4.13 Passing one source client to `Search` silently disables the others

Found while validating §9, not by reading the code. `Search.__init__` ends with

```python
self._clients_initialized: bool = any(
    c is not None for c in [chebi, comptox, pubchem, zeropm, chembl]
)
```

and `_ensure_clients` opens with `if not self._clients_initialized:`. So handing
`Search` a single pre-built client — the obvious thing to do when you want to
point it at a particular ChEMBL file, which is exactly what §9's validation
needed — marks *all five* as initialised and leaves the other four `None`.

`Search("cas", chembl=CheMBL(db_path=...))` therefore runs on one source. It
says so, at WARNING, through the §20 availability message, and every result
comes back empty because `_resolve_cas` only reaches ChEMBL through a SMILES the
primary sources were supposed to supply. But nothing refuses the call, and a
caller who is not watching the log gets a table of empty rows with
`confidence` 0.0 rather than an error.

The fix is small: track initialisation per client rather than with one flag, so a
caller can supply one source and have the rest built lazily as documented. Size
S. Until then, passing any client means passing all of them.

---

## 5. Classification is deferred

Per the user's instruction, `taxonomy.py` (the chebifier backend) and the
classification algorithm are **out of scope** for this round. They are
experimental and will be revisited at a later stage of the refactor.

One caveat that should not wait, because it is hygiene rather than algorithm:
`examples/chebifier/__pycache__/` is still tracked in git, which §11 asked to
remove and §21 did not reach.

---

## 6. Sequencing

Steps are independently committable and leave the suite green.

| # | Commit | §  | Size |
|---:|---|---|---|
| 1 | ~~**`CheMBL.compact()`: 30 GB → 2.6 GB from a full database already on disk**~~ **done, §11** | **9** | **S** |
| 2 | ~~rewrite `search_by_name` as a `UNION`; add `ORDER BY` — 743 ms → 10 µs~~ **done, §12** | 9.5, 9.3 | S |
| 3 | ~~`datasets.py`: `download_file` with resume + checksum; migrate the 5 download sites~~ **done, §13** | 4.2 | M |
| 4 | ~~`datasets.status/plan/fetch/remove`; `Search(datasets="present")` as the default~~ **done, §14** | 4.1 | M |
| 5 | ~~build the ChEMBL extract during download; `source=` on `CheMBL`~~ **done, §15** | 9.7 | M |
| 6 | ~~`close()`, context managers and thread safety on the four SQLite clients~~ **done, §16** | 4.3 | S |
| 7 | `classyfire.py` raises; rewrite its tests; drop it from the docs' service lists — **postponed, §17.0** | 4.8 | S |
| 8 | ~~`cache.py`: persistent dir, parameterised service functions, key version~~ **done, §17** | 4.7 | S |
| 9 | ~~`pubchem_ftp.py`: the FTP identifier builder~~ **done, §18** | 8 | L |
| 10 | ~~`PubChemID.descriptors()` — RDKit descriptors on demand~~ **done, §19** | 8.6 | M |
| 11 | ~~`CheMBL(source="mysql")` — the streaming route, verified against step 5~~ **done, §20** | 9.7 | M |
| 12 | ~~split `pubchem.py` → `+ pubchem_id.py`; `chebi.py` → `+ chebi_sdf.py`~~ **done, §21** | 4.12 | M |
| 13 | ~~`sources.py`; collapse the `Search` ladders~~ **done, §22** | 4.5 | M |
| 14 | ~~`Search(online_fallback=...)` as a row in the source table~~ **done, §23** | 4.4 | M |
| 15 | ~~`Search.PRESETS`~~ **done, §24** | 4.6 | S |
| 16 | ~~circuit breaker on the shared `RateLimiter`~~ **done, §27** | 4.12 | M |
| 17 | ~~docstrings with examples, module by module~~ **done, §28** | 4.12 | L |
| 18 | ~~rebuild `docs/`; delete `docs/examples/`; drop `docs/plans/` from the nav~~ **done, §29** | 4.9 | M |
| 19 | notebooks, `search/` first | 4.11 | L |
| 20 | rewrite `README.md` around offline-first and `Search` | 4.10 | S |

Step 1 comes first because it is the largest single win in the document, costs
almost nothing, and needs no download — it reclaims 27 GB on a machine that
already has ChEMBL. Step 2 is an unrelated S-sized fix that fell out of measuring
step 1 and pays for itself immediately. Steps 3–4 follow because every later
step's first-run story depends on them, and because §8's builder is written on
top of step 3's downloader. Step 5 before 11, since 11 is verified against 5. Step 13 before 14, and 17 before 18.

## 7. Verification

Unchanged from §13 of the August plan: `pytest` green locally at every commit, no
GitHub Actions test jobs, `mkdocs build --strict` clean, notebooks executed before
their outputs are committed. Two additions:

- §8's builder needs a **small-scale test**: build from a 10 000-CID slice of the
  real FTP files and assert the schema, the row counts and a handful of known
  CAS→CID pairs. The full build is too slow for the suite.
- §4.1's default change needs a test that constructs `Search` with no datasets
  present and asserts that **nothing is downloaded** and the error names the
  `datasets.fetch` call.
- §9's extract needs an **equivalence test**, not just a smoke test: with both a
  full and a compacted ChEMBL present, assert that every public method of
  `CheMBL` returns the same object for a fixed sample of molregnos. The
  prototype ran this over 2 000 random compounds; the committed version can use
  a small fixed list so it runs anywhere, with the 2 000-compound sweep as a
  marked slow test. Route 3 (§9.7) additionally needs a per-table row-count and
  content-hash comparison against route 2 before it is trusted.

---

## 8. The extra section: building the PubChem identifier database from the FTP site

### 8.1 What exists today, and why it should change

`pubchem_id.db` is built by a manual, unreproducible pipeline:

1. A human opens PubChem's classification browser, navigates to
   *Names and Identifiers → Other Identifiers → CAS*, and downloads a CSV.
2. `scripts/build_pubchem_id_db.py` reads that CSV and **regexes CAS numbers,
   InChI and InChIKey out of the `cmpdsynonym` free-text column**
   (`CAS_PATTERN = r'\b\d{2,7}-\d{2}-\d\b'`).
3. The resulting 2.2 GB SQLite file is uploaded to Zenodo by hand, and
   `PubChemID.DEFAULT_DB_URL` is edited to point at it.

Three things are wrong with this, and the third is measurable.

**It cannot be reproduced.** The CSV is a point-and-click artifact. Nothing in
the repo records which PubChem release it came from — the filename says
`202601`, and that is the entire provenance.

**It cannot be refreshed without a person.** Every update is the same manual
sequence, so the database ages between releases.

**It is measurably wrong.** PubChem publishes its own curated CAS↔CID mapping in
`Compound/Extras/CID-Identifiers.tsv.gz` — third-party database identifiers with
an explicit type column. Comparing it against the shipped `pubchem_id.db`
(measured 2026-09-20, against the FTP file regenerated 2026-09-19):

| | shipped DB (regex) | FTP mapping | |
|---|---:|---:|---|
| CAS rows | 1 400 544 | 1 462 322 | |
| distinct (CID, CAS) pairs | 1 389 628 | 1 462 322 | |
| distinct CIDs carrying a CAS | 1 323 167 | **1 431 500** | +108 333 |
| **rows failing the CAS check digit** | **17 379** | **122** | 0.008% vs 1.25% |
| CIDs present only in the other | 16 821 | **125 154** | |

The regex over free-text synonyms produces **17 379 CAS numbers that are not
valid CAS numbers** — 1.25% of the table — because `\d{2,7}-\d{2}-\d` matches
plenty of things that are not registry numbers. `provesid.utils.check_CASRN`
already implements the check-digit test; the build script never calls it.
Meanwhile the authoritative mapping covers **125 154 CIDs the shipped database
has never heard of**.

So the FTP route is not merely more convenient. It is more complete by 8% and
roughly 150× cleaner on the one field the whole database is organised around.

### 8.2 Decisions taken

Confirmed with the user on 2026-09-20:

1. **Drop the eight computed descriptors** (`XLogP`, `TPSA`, `Complexity`,
   `Charge`, `HBondDonorCount`, `HBondAcceptorCount`, `RotatableBondCount`,
   `HeavyAtomCount`) from the stored database. They are *computed from the
   structure*, not data about the substance, and the database's job is
   identifiers and structures. Compute them with RDKit **on demand**, when the
   user asks for them (§8.6).
2. **Scope: CAS-bearing compounds only** — about 1.43 M CIDs, the same shape as
   today. The package targets researchers on laptops; hundreds of gigabytes are
   not affordable, especially beside ChEBI and ChEMBL.
3. **FTP is the default path.** Zenodo stays as an alternative route that the
   user refreshes every few months by hand. Downloads are **on demand, not on a
   schedule** — a weekly refresh would be overkill and unnecessary load.

### 8.3 What the FTP site actually offers

All measurements taken 2026-09-20 from `https://ftp.ncbi.nlm.nih.gov/pubchem/`.
The HTTPS mirror answers `Range` requests with HTTP 206 (verified) and publishes
an `.md5` beside every file (verified), so §4.2's resumable, checksummed
downloader applies directly.

`Compound/Extras/`, regenerated with every PubChem dump (`Last-Modified:
2026-09-19` on every file below):

| File | Compressed | Content |
|---|---:|---|
| `CID-Identifiers.tsv.gz` | **98 MB** | CID, identifier, identifier type — **1 462 322 `CAS` rows** |
| `CID-SMILES.gz` | 1.49 GB | CID → isomeric SMILES |
| `CID-InChI-Key.gz` | 7.37 GB | CID → InChI, InChIKey |
| `CID-Mass.gz` | 1.39 GB | CID → formula, monoisotopic mass, exact mass |
| `CID-Title.gz` | 1.89 GB | CID → compound title (the `cmpdname` column) |
| `CID-IUPAC.gz` | 1.85 GB | CID → computed IUPAC name |
| `CID-Synonym-filtered.gz` | 0.97 GB | CID → synonym, structure-consistent names only |
| **total** | **15.05 GB** | |

Measured throughput from this machine: **12.2 MB/s**, so the full set is
**~21 minutes**. Skipping `CID-InChI-Key.gz` and deriving InChI/InChIKey from the
SMILES with RDKit brings it to 7.7 GB and ~11 minutes, at the cost of using
RDKit's InChI rather than PubChem's — worth offering as a flag, not as the
default.

`CID-Identifiers.tsv.gz` carries far more than CAS, all for the same 98 MB, and
several of these are identifiers `Search` already reconciles across sources:

| type | rows | | type | rows |
|---|---:|---|---|---:|
| Nikkaji Number | 3 584 326 | | European Community (EC) Number | 355 723 |
| ChEMBL ID | 2 929 213 | | NSC Number | 317 185 |
| **CAS** | **1 462 322** | | HMDB ID | 203 807 |
| Wikidata | 1 326 468 | | ChEBI ID | 179 537 |
| DSSTox Substance ID | 1 149 154 | | UNII | 129 279 |

The **1 149 154 DSSTox (DTXSID) rows and 179 537 ChEBI rows** are notable: they
are a free cross-source identifier bridge in a file that is already being
downloaded, and `Search` currently has to *infer* those links through structure
matching.

**Substances.** `Substance/Extras/SID-Map.gz` (3.58 GB) gives SID, source name,
registry identifier and the standardised CID. Sources that register CAS numbers
as their registry identifier yield a substance-level CAS layer, and
`Substance/Extras/Source-Names` (72 KB) maps depositor IDs to display names. This
is the honest route for the "substances with CAS number" half of the current
manual download, but it is a second, larger pipeline and §8.5 keeps it out of the
first version.

**What is not there.** The eight computed descriptors of §8.2 exist only in the
full SDF dump, `Compound/CURRENT-Full/SDF/` — **359 files of ~339 MB each, about
120 GB**, of which 358 contain at least one CAS-bearing CID (they are spread
across the whole CID range: half below CID 14.4 M, but 5% above 138 M). This is
exactly why decision 1 was taken.

**Reproducibility.** `Compound/Monthly/YYYY-MM-01/` holds frozen snapshots — the
same `Extras/` file set plus a `TIMESTAMP` file (`2026-09-01/` →
`2026/08/31 18:38:04`) and `killed-CIDs` / `updated-CIDs` lists. Building against
a dated monthly snapshot rather than the rolling current dump makes a build
reproducible and citable, and matches decision 3's "not weekly" cadence exactly.
**Default to the newest monthly snapshot**, with the rolling `Extras/` available
as an opt-in.

### 8.4 Proposed design: `src/provesid/pubchem_ftp.py`

```python
from provesid.pubchem_ftp import build_pubchem_id_db, list_releases

list_releases()
# ['2026-09-01', '2026-08-01', '2026-07-01', '2026-06-01', 'current']

build_pubchem_id_db(
    release="2026-09-01",   # a monthly snapshot; "current" for the rolling dump
    scope="cas",            # CAS-bearing CIDs only
    include_inchi=True,     # False → derive from SMILES with RDKit, saves 7.4 GB
    include_synonyms=True,
    keep_downloads=False,   # delete each source file once it has been parsed
)
```

and, as the user-facing entry point, `PubChemID` gains a third acquisition route
beside the existing Zenodo download:

```python
PubChemID(source="ftp")      # build from PubChem FTP — proposed default
PubChemID(source="zenodo")   # download the prebuilt artifact
PubChemID(source="local")    # use what is on disk, or raise
```

**The algorithm**, in one pass per file and constant memory:

1. Fetch `CID-Identifiers.tsv.gz` (98 MB). Select the `CAS` rows, validate each
   with `check_CASRN`, and hold the CID set — 1 431 500 integers in a sorted
   `array('l')` is about 11 MB. This set is the scope filter for everything else.
2. For each remaining file, in turn: download with resume and md5 check
   (§4.2), stream-decompress, keep only lines whose CID is in the set, insert in
   batches, then delete the file unless `keep_downloads=True`.

Streaming one file at a time is what makes this fit on a laptop: peak transient
disk is the largest single file (6.9 GiB for `CID-InChI-Key.gz`, or 1.8 GiB if
`include_inchi=False`), not the 15 GB total. The membership test makes no
assumption about the files being CID-sorted — they appear to be, but a build that
silently depends on it would be fragile.

3. Compute `mw` from the formula in `CID-Mass.gz` using standard atomic weights.
   PubChem's own molecular weight is exactly this, so the column stays faithful
   while `CID-Mass.gz` supplies formula, monoisotopic mass and exact mass
   directly.
4. Build the indexes the current schema already defines, then write a
   `provenance` table — release, per-file URL, md5, byte count, row counts,
   build timestamp, builder version — so a database on disk can always say where
   it came from. This is what the current artifact cannot do.

**Schema.** Keep `compounds` / `cas_numbers` / `synonyms` unchanged, minus the
eight descriptor columns, plus the `provenance` table. `PubChemID.OFFLINE_PROPERTIES`
drops to the eight columns that remain data:

| kept (from FTP) | dropped (computed → §8.6) |
|---|---|
| `MolecularFormula`, `MolecularWeight`, `ExactMass` | `XLogP`, `TPSA`, `Complexity` |
| `SMILES`, `InChI`, `InChIKey` | `Charge`, `HeavyAtomCount` |
| `IUPACName`, `Title` | `HBondDonorCount`, `HBondAcceptorCount`, `RotatableBondCount` |

Add `MonoisotopicMass`, which `CID-Mass.gz` supplies for free and the current
database lacks. Optionally add a `xrefs` table from the same
`CID-Identifiers.tsv.gz` pass, carrying the DTXSID, ChEBI, ChEMBL, EC and UNII
rows of §8.3 — 5 M rows, and a direct answer to cross-source questions `Search`
currently resolves by structure matching.

**Estimated cost of a build** (the one set of numbers here that is projected
rather than measured): ~21 minutes of download at the measured 12.2 MB/s, plus a
parse and insert pass dominated by SQLite write throughput — on the order of
20–40 minutes. Resulting database roughly 1.5–2.5 GB, comparable to today's
2.2 GB: the eight dropped columns are small next to the ~10 M synonym rows.
These should be measured on the first real build and written into the
documentation rather than guessed at.

### 8.5 What the first version should not do

- **No substance layer.** `SID-Map.gz` is 3.58 GB and answers a different
  question. Add it once the compound builder is in use, as `scope="cas+sid"`.
- **No SDF harvest.** 120 GB for eight computed columns, per decision 1.
- **No incremental update.** `Compound/Daily/` and `Weekly/` exist, and
  `Monthly/*/killed-CIDs` and `updated-CIDs` make a real incremental update
  possible, but decision 3 sets the cadence at on-demand. A rebuild against a
  newer monthly snapshot is simpler and, at ~40 minutes, cheap enough.
- **No change to Zenodo.** It keeps working, unchanged, as `source="zenodo"`.

### 8.6 The other half of decision 1: descriptors on demand

Dropping eight columns is only acceptable if a user who wants them can still get
them. Two routes, and the user's instruction is to compute:

```python
db = PubChemID()
db.descriptors(2244)                       # RDKit, from the stored SMILES
db.descriptors(2244, source="pubchem")     # PubChem's own values, via PUG-REST
```

`descriptors()` reads the SMILES already in the database and computes TPSA,
`MolLogP`, H-bond donors and acceptors, rotatable bonds, heavy-atom count and
formal charge with RDKit — no network, milliseconds per molecule, and RDKit is
already a hard dependency. The values must be **labelled as RDKit's**, because
they will not always equal PubChem's: `XLogP3` is a different model from Crippen
`MolLogP`, and `Complexity` (Bertz/Hendrickson/Ihlenfeldt) has no RDKit
equivalent at all — `BertzCT` is related but not the same number. So:

- The returned record carries `Source='rdkit'`, exactly as the existing property
  layer carries `Source='offline'` / `'online'` (§18).
- `Complexity` is **not** offered from RDKit. It is PUG-REST or nothing, and the
  docstring says so.
- `descriptors(..., source="pubchem")` routes through the existing offline→online
  property machinery, so a user who needs PubChem's exact numbers keeps them.

This also removes a real inconsistency: today `properties()` returns PubChem's
`XLogP` from the local database for compounds that happen to be in it, and
PubChem's `XLogP` over the network for those that are not, with no way to ask for
either deliberately. After this change the three sources are named and selectable.

### 8.7 What this buys

| | today | after |
|---|---|---|
| How the database is built | manual CSV download + regex | one function call |
| Provenance | a number in a filename | a `provenance` table with md5s per file |
| CIDs with a CAS | 1 323 167 | **1 431 500** |
| Invalid CAS rows | 17 379 (1.25%) | ~122 (0.008%) |
| Refresh | a person, hours | `build_pubchem_id_db(release=...)`, ~40 min |
| Reproducible | no | yes — a dated monthly snapshot |
| Cross-source identifiers | inferred by structure matching | 5 M rows, free |
| Descriptor columns | 8, stored, unlabelled origin | computed on demand, labelled |
| Zenodo dependency | required | optional |

---

## 9. Shrinking ChEMBL from 30 GB to 2.6 GB

§4.1 said ChEMBL deserved a second look on its own terms: 28 GB resident, 87% of
the package's total footprint, for a source that only *enriches* a `Search`
result. This section is that look. It is written from a working prototype, not
from a design sketch — the extract described here was built, measured and
validated against the full database on 2026-09-20.

**Result: 29.74 GB → 2.60 GB, built in 31 seconds, returning the same compounds
from every public method.** Measuring it also turned up a second, unrelated prize
that has nothing to do with disk: `search_by_name` costs 743 ms per call and can
cost 10 µs (§9.5).

### 9.1 The package uses 8 of ChEMBL's 74 tables

Every SQL statement in `chembl.py` touches one of eight tables:

| Table | Rows | Used by |
|---|---:|---|
| `molecule_dictionary` | 2 878 135 | every lookup; `get_molecule_dictionary` returns all 28 columns |
| `compound_structures` | 2 854 815 | `search_by_inchi` / `_inchikey` / `_smiles`, `get_compound` |
| `compound_properties` | 2 858 458 | `get_properties`; `mw_freebase` is what `Search` reads |
| `molecule_hierarchy` | 2 786 844 | `get_molecule_hierarchy` |
| `chembl_id_lookup` | 5 352 300 | `chembl_id_to_molregno`, filtered to `entity_type='COMPOUND'` |
| `molecule_synonyms` | 132 937 | `search_by_name`, `get_compound` |
| `pesticide_classification` | 595 | the four pesticide methods |
| `pesticide_class_mapping` | 593 | the four pesticide methods |

The other 66 tables are never opened. They include everything ChEMBL is actually
famous for: `activities` (24 267 312 rows), `activity_properties` (11 898 213),
`compound_records` (3 774 137), `assays` (1 890 749), and the whole target,
document and binding-site apparatus. PROVESID uses ChEMBL as a structure and
name index, not as a bioactivity database.

### 9.2 Where the 30 GB goes

Measured with `dbstat` against the local `chembl_36.db` (29.74 GB, 4 096-byte
pages, 7 260 702 pages). The eight used tables, with their indexes, come to
**14.87 GB** — so half the file is tables the package never reads, and half is
tables it reads inefficiently:

| Used table (incl. indexes) | Size |
|---|---:|
| `compound_structures` | **13.21 GB** |
| `compound_properties` | 584 MB |
| `chembl_id_lookup` | 477 MB |
| `molecule_dictionary` | 460 MB |
| `molecule_hierarchy` | 126 MB |
| `molecule_synonyms` | 12 MB |
| pesticide tables | 0.2 MB |
| **total** | **14.87 GB** |

`compound_structures` is 89% of that, and one column explains it. Summing column
lengths across all 2.85 M rows:

| Column | Total bytes |
|---|---:|
| **`molfile`** | **7.72 GB** |
| `standard_inchi` | 468 MB |
| `canonical_smiles` | 167 MB |
| `standard_inchi_key` | 77 MB |

**A quarter of the entire 30 GB database is the `molfile` column** — the full MOL
block for every compound. PROVESID reads it in exactly one place: `get_compound`
selects `cs.molfile` (`chembl.py:766`). Nothing consumes it. `Search` never looks
at it; `candidate_from_chembl_row` (`tools.py:751`) reads `pref_name`,
`canonical_smiles`, `standard_inchi`, `standard_inchi_key`, `synonyms` and
`mw_freebase`, and nothing else. No test asserts on it, no example uses it. It is
mentioned in one hand-written doc page (`docs/api/chembl.md:66`) that §4.9 deletes
anyway. And it is redundant: a MOL block is reconstructible from the SMILES with
RDKit, which is already a hard dependency.

The rest is index weight that serves queries the package never issues.
`compound_properties` carries nine indexes (~330 MB) on `alogp`, `psa`, `rtb`,
`hba`, `hbd`, `mw_freebase` and the Ro5 counts — every one of them for range
queries, while `chembl.py` only ever does `WHERE molregno = ?`.

### 9.3 The measured extract

The prototype copies the eight tables, drops `molfile`, keeps only
`entity_type='COMPOUND'` in the lookup table, and builds only the indexes the
package's queries need:

```
molecule_dictionary          2,878,135 rows     2.7s
compound_structures          2,854,815 rows    21.6s
compound_properties          2,858,458 rows     1.9s
molecule_synonyms              132,937 rows     0.1s
molecule_hierarchy           2,786,844 rows     0.8s
chembl_id_lookup             3,038,372 rows    11.4s   (from 5,352,300)
pesticide_classification           595 rows     0.0s
pesticide_class_mapping            593 rows     0.0s
...14 indexes, ANALYZE, VACUUM
chembl_36_lite.db: 2.60 GB in 31s
```

Composition of the result:

| Table (incl. indexes) | Extract | Full | |
|---|---:|---:|---|
| `compound_structures` | 1 679 MB | 13 206 MB | −87% |
| `molecule_dictionary` | 334 MB | 460 MB | −27% |
| `compound_properties` | 290 MB | 584 MB | −50% |
| `chembl_id_lookup` | 199 MB | 477 MB | −58% |
| `molecule_hierarchy` | 87 MB | 126 MB | −31% |
| `molecule_synonyms` | 10 MB | 12 MB | |
| pesticide tables | 0.1 MB | 0.2 MB | |
| **total** | **2.60 GB** | **29.74 GB** | **−91.3%** |

**Validation.** With `cs.molfile` removed from `get_compound`'s SELECT, every
public method of `CheMBL` returns the same *compounds* as the full database:

- all 15 public methods, exercised on aspirin (`CHEMBL25`), an inexact name
  search, an InChI, an InChIKey, a SMILES and a pesticide query — **0
  differences**;
- `get_compound` over **2 000 randomly chosen molregnos** — **0 differences**;
- 500 `search_by_smiles` calls: 0.02 s on both.

One qualification, and it is a finding about the existing code rather than about
the extract. Over 200 random preferred names, `search_by_name(exact=True)`
returned a **different list order** on the extract for 63 of them — and on closer
inspection **every one of those 63 was an ordering difference with identical
content**: re-checked over a fresh sample of 60 names, 25 differed in order and
**0 differed in content**.

The cause is that `search_by_name` issues `SELECT DISTINCT … LIMIT ?` with **no
`ORDER BY`**, so the row order is whatever the query plan happens to produce.
That is unstable on the *current* database too — a `VACUUM`, an index change or a
SQLite upgrade can reorder it — and when a query matches more than `limit` rows
it silently changes *which* rows come back, not just their order. The extract did
not introduce this; it exposed it.

Two consequences. `search_by_name` should get a deterministic `ORDER BY
md.molregno` (an S-sized fix, and one worth making regardless of §9). And the
equivalence test of §7 must compare **sets of molregnos**, not lists, or it will
fail for the wrong reason.

### 9.4 The one code change it requires

```diff
--- a/src/provesid/chembl.py
@@ get_compound
                 cs.canonical_smiles,
                 cs.standard_inchi,
-                cs.standard_inchi_key,
-                cs.molfile
+                cs.standard_inchi_key
             FROM molecule_dictionary md
```

plus the two docstring lines that mention it (`chembl.py:12`, `chembl.py:739`).
That is the whole migration. Under dev-principle §1 no shim is needed; callers
who want a MOL block get it from RDKit, and the `get_compound` docstring should
say so.

### 9.5 A second prize: `search_by_name` is 743 ms per call, and needn't be

Measuring the extract turned up something the disk-space question would never
have surfaced. `search_by_name(exact=True)` costs **743 ms per lookup** on the
full database — 200 random preferred names took 148.6 s. On the extract it is
602.8 ms, which is only the reward for scanning a smaller table.

`EXPLAIN QUERY PLAN` says why, and it is the same on both:

```
SCAN md USING INDEX sqlite_autoindex_molecule_dictionary_1
SEARCH ms USING COVERING INDEX ... (molregno=?) LEFT-JOIN
```

A full scan of 2.88 M rows, every call. The query is

```sql
SELECT DISTINCT md.molregno FROM molecule_dictionary md
LEFT JOIN molecule_synonyms ms ON md.molregno = ms.molregno
WHERE LOWER(md.pref_name) = LOWER(?) OR LOWER(ms.synonyms) = LOWER(?)
```

and it cannot use an index for two compounding reasons: `LOWER(col) = ?` defeats
an index on `col`, and an `OR` whose arms live in *different tables* across a
`LEFT JOIN` defeats indexing altogether. Adding expression indexes on
`lower(pref_name)` and `lower(synonyms)` changes nothing on its own — I added
them to the prototype and the plan still shows a full scan. **The index is not
the fix; the query shape is.**

Rewritten as a union of two independently indexable lookups:

```sql
SELECT molregno FROM molecule_dictionary WHERE lower(pref_name) = lower(?)
UNION
SELECT molregno FROM molecule_synonyms  WHERE lower(synonyms)  = lower(?)
```

the plan becomes two index searches, and the measurement is:

| | per exact name lookup |
|---|---:|
| full database, current query | 743.2 ms |
| extract, current query | 602.8 ms |
| **extract, `UNION` + expression indexes** | **10 µs** |

Measured over 20 000 lookups (0.19 s total). That is a **~77 000× speedup**, and
the rewrite returns identical molregnos — 0 disagreements over 200 names. The two
expression indexes cost a few tens of MB and are already inside the 2.60 GB
above.

This matters beyond `CheMBL` itself: `Search` reaches ChEMBL by SMILES rather
than by name, so it does not pay this today — but any name-driven enrichment over
a few thousand rows currently costs 743 ms each, which is forty minutes for a
3 000-row dataset. The rewrite is an S-sized change to one method and is worth
making whether or not the extract lands.

Note that the `UNION` also fixes the ordering instability of §9.3 only partly —
`LIMIT` without `ORDER BY` is still nondeterministic when a name matches more
rows than the limit, so the `ORDER BY md.molregno` recommended there is still
needed.

### 9.6 A leaner variant, and why not to take it

Dropping `standard_inchi` as well — deriving it from the SMILES with RDKit when a
caller wants it, and routing `search_by_inchi` through an InChI→InChIKey
conversion — gives **1.55 GB**, another 1.05 GB saved (the column is 468 MB and
its two indexes are ~1.17 GB).

I do not recommend it as the default. `candidate_from_chembl_row` passes ChEMBL's
`standard_inchi` straight into the candidate record, so dropping it either
degrades `Search` output or puts an RDKit InChI round-trip on the hot path for
every ChEMBL hit. The 1.05 GB is not worth that. Keep it available as
`compact(keep_inchi=False)` for users genuinely short of disk, and document what
it costs.

### 9.7 How the user gets the extract

This is the part that needs a decision, because the extract has to be built from
*something*, and ChEMBL does not publish the eight tables as flat files. Probed
on 2026-09-20 (`https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/`,
now at release 37; `Range` answered with HTTP 206):

| File | Size | Usable? |
|---|---:|---|
| `chembl_37_sqlite.tar.gz` | 5.76 GB | yes — today's path |
| `chembl_37_mysql.tar.gz` | 2.10 GB | yes — plain-text SQL, one multi-row `INSERT INTO` per table |
| `chembl_37_postgresql.tar.gz` | 2.09 GB | **no** — `PGDMP` custom format, needs `pg_restore` |
| `chembl_37_chemreps.txt.gz` | 293 MB | partial — structures only, no names or properties |
| `chembl_37.sdf.gz` | 931 MB | partial — its only SD tag is `chembl_id` |

`chemreps` and the SDF are tempting and insufficient: neither carries
`pref_name`, `mw_freebase`, synonyms, the hierarchy or the pesticide tables.
Those exist only inside the SQL dumps. So three routes, in the order I would
build them:

**Route 1 — `CheMBL.compact()`, from a full database already on disk.**
No download at all. 31 seconds, then delete the 30 GB. This is the immediate win
and should land first: it turns an existing installation from 30 GB into 2.6 GB
today, and it is the same code path the other two routes reuse.

**Route 2 — build during download, delete the full file.** The `sqlite.tar.gz`
route, unchanged except that `download_database` finishes by building the extract
and removing the intermediate. Steady state 2.6 GB instead of 30 GB. The cost is
a transient peak of ~31 GB while the full file exists, which is exactly the
laptop problem §4.1 is about — so the peak must be *stated before the download
starts*, and `datasets.plan()` (§4.1) is where it belongs.

**Route 3 — stream the MySQL dump.** 2.10 GB downloaded, decompressed on the
fly, with a parser that recognises `INSERT INTO` for the eight tables and ignores
everything else. Peak transient disk ≈ 0 beyond the 2.6 GB result, and the
download is 3.7 GB smaller than route 2. The cost is a MySQL-dump parser —
mysqldump's escaping is well defined (`\'`, `\\`, `\n`, `\r`, `\0`, `\Z`, bare
`NULL`), but it is still a parser where routes 1 and 2 have none.

Route 3 is the right end state for the constraint the user actually stated —
laptops, direct from source, no hosting — and its correctness is checkable
against route 2: build both once, compare row counts and a content hash per
table. I would not write it until routes 1 and 2 are in use, and I would not
ship it without that comparison as a test.

**Not recommended: hosting a prebuilt 2.6 GB extract.** It is the easiest thing
for users, but it re-creates exactly the hosting dependency §8.2 decided to move
away from, and it makes PROVESID a redistributor of ChEMBL rather than a reader
of it. ChEMBL's licence (CC BY-SA) permits redistribution with attribution, so
this stays available as a fallback if routes 2 and 3 both prove impractical —
but it should be the last choice, not the first.

### 9.8 Proposed API

```python
from provesid import CheMBL

# Route 1 — shrink what is already on disk.
CheMBL.compact()                       # find the full db, build the extract
CheMBL.compact(remove_source=True)     # …and delete the 30 GB when it verifies
CheMBL.compact(keep_inchi=False)       # the 1.55 GB variant (§9.6)

# Routes 2 and 3 — acquisition.
CheMBL(source="mysql")    # stream the 2.10 GB dump           (route 3)
CheMBL(source="sqlite")   # download 5.76 GB, extract, compact (route 2)
CheMBL(source="full")     # today's behaviour: keep the 30 GB database
```

`compact` writes a `provenance` table in the same shape §8.4 proposes for
`pubchem_id.db`: release, source URL, source checksum, per-table row counts,
build timestamp, and the PROVESID version that built it. `CheMBL.__init__` reads
it and logs the release, so a compacted database can always say what it came
from — and so a later PROVESID that needs a ninth table can detect that the
extract predates it and say so, instead of failing with `no such table`.

That last point is the real risk of this change and deserves a guard: an extract
is a *subset*, and subsets go stale differently from copies. A
`schema_version` in `provenance`, bumped whenever the eight-table list changes,
turns "no such table: compound_records" into "this ChEMBL extract was built for
PROVESID 0.7 and is missing tables 0.9 needs; run `CheMBL.compact(rebuild=True)`".

### 9.9 Effect

| | today | after |
|---|---:|---:|
| ChEMBL on disk | 29.74 GB | **2.60 GB** |
| …as a share of PROVESID's ~32 GB | 87% | 40% of a new ~6.5 GB total |
| Tables present | 74 | 8 |
| Download to acquire | 5.76 GB | 2.10 GB (route 3) |
| Peak transient disk | ~36 GB | ~2.6 GB (route 3), ~31 GB (route 2) |
| Build time from a local full db | — | 31 s |
| `search_by_name(exact=True)` | 743 ms | **10 µs** (§9.5) |
| Public methods affected | — | none (§9.3 on result *order*) |
| Code changes | — | one SELECT, two docstrings, one query rewrite |

Together with §8, this is what makes the package usable on a laptop: PubChem
~2 GB, ChEMBL 2.6 GB, ChEBI 0.9 GB, CompTox 0.8 GB — **about 6.5 GB for all four
`Search` sources**, down from ~32 GB.

---

## 10. Explicitly not doing

- The classification algorithm and `taxonomy.py`, per §5.
- Splitting `zeropm.py`. Unchanged from the August plan: it is one coherent
  SQLite interface.
- The confidence/consensus algorithm in `Search` beyond the §4.5 restructuring,
  which must be behaviour-preserving and is pinned by
  `tests/test_search_precision_regression.py`.
- Reviving ClassyFire.
- Backward-compatibility shims, per dev-principle §1.

---

## 11. Landed on 2026-09-20 — step 1, `CheMBL.compact()` (§9)

§9 was written from a throwaway prototype. This is the same thing as shipped
code, with the differences that mattered.

### 11.1 What landed

- `CheMBL.compact()`, a classmethod, with `source_path`, `dest_path`,
  `data_dir`, `keep_inchi`, `remove_source` and `force`.
- `CheMBL.is_compact_database()`, `read_provenance()` and `compact_path_for()`
  as the public surface around it.
- `COMPACT_TABLES`, `COMPACT_INDEXES`, `COMPACT_SCHEMA_VERSION`,
  `COMPACT_SUFFIX` and `PROVENANCE_TABLE` as documented class attributes, so
  *which eight tables* is data rather than something buried in a method.
- `CheMBL.is_compact` and `CheMBL.provenance` on every instance.
- `_find_local_database` now prefers an extract over a full release **of the
  same release number**, and a newer full release over an older extract.
- `get_compound` no longer selects `cs.molfile` (§9.4), so a full release stops
  returning it too. That is the only behaviour change for existing callers.
- 35 tests in `tests/test_chembl_compact.py`, an example in
  `examples/chembl/compact_demo.py`, and a `CHANGELOG.md` entry.

Measured on the real ChEMBL 36, through the shipped code rather than the
prototype: **29.74 GB → 2.60 GB, 91.3% smaller, in 61 s** (the prototype's 31 s
plus verification and `ANALYZE`).

### 11.2 Verification is the feature, not the copy

The prototype proved a 2.6 GB file could answer the same questions. Shipping it
means `remove_source=True` destroys 30 GB on the strength of that claim, so the
build now proves itself before anything is deleted:

1. SQLite's `quick_check` on the result.
2. A row count per table against the source **under the same filter** — the
   `chembl_id_lookup` count has to match the source's `WHERE entity_type =
   'COMPOUND'` count, not its total.
3. 500 compounds re-read from both databases and compared column by column.

The extract is written to `<dest>.tmp` and only `os.replace`d into position once
all three pass; a failure removes the partial file and leaves the source alone.
Two tests corrupt a verified extract on purpose — one altering a SMILES, one
deleting rows — and assert that `remove_source=True` raises and the 30 GB
survives.

### 11.3 Sampling had to be fixed to be worth anything

The first implementation drew its 500 compounds from `ORDER BY molregno LIMIT
10000`, which is fast but samples only the lowest 0.35% of the table — a
verification that would miss any damage to the other 99.65%. `ORDER BY
random()` is unbiased but scans 2.9 M rows.

`_verification_sample` now reads `MIN`, `MAX` and `COUNT`, draws random values
from that range and keeps the ones that exist. A few hundred index seeks, spread
across the whole table, with a fixed seed so a failure reproduces. A table
smaller than the sample size is compared in full, which is what makes the
miniature ChEMBL in the tests exercise the same code path as the real one.

### 11.4 Three things §9 did not anticipate

**`_find_local_database` had to change, or the extract would never be opened.**
Its regex was `chembl_(\d+)\.db$`, which `chembl_36_provesid.db` does not match
— it would have scored −1 and lost to the 30 GB file every time. The ordering is
now (release, is-extract), so the extract wins a tie and a newer release still
wins outright.

**The provenance table needed a schema version, and §9.8 was right to insist.**
An extract is a subset, so a later PROVESID needing a ninth table would meet
`no such table` with no clue the database was merely old.
`_warn_if_extract_is_stale` compares `COMPACT_SCHEMA_VERSION` on open and names
`CheMBL.compact(force=True)`. It warns rather than raises: an old extract still
answers everything the package asked of it when it was built.

**`compact()` must refuse to compact an extract.** Running it twice in a
directory would otherwise produce `chembl_36_provesid_provesid.db` from a source
that has already lost `molfile`. `_resolve_compact_source` skips extracts when
searching, and an explicit extract path raises.

### 11.4a Two hazards found while validating, both now fixed

**A doctest in `download_database` really downloads 5.8 GB.** Running
`pytest --doctest-modules src/provesid/chembl.py` executes

```
>>> chembl = CheMBL(auto_download=False)
>>> chembl.download_database(force=True)
```

which is not an illustration but an instruction. It was caught here only
because it started writing `chembl_36_provesid.db.tar.gz.tmp` into the dataset
directory — that is, it had resolved `db_path` to the *extract* and was on
course to overwrite it. The example predates this work and was always
destructive; what changed is which file it would have destroyed. Both lines are
now `# doctest: +SKIP`, with a comment saying why.

`pyproject.toml` sets `testpaths = ["tests"]`, so the ordinary suite never
collects these. `pubchem.py` has 26 docstring examples with the same shape
(§23.7a of the August plan), and they should be audited before anything starts
running doctests across `src/`.

**`with sqlite3.connect(...)` does not close the connection.** It manages a
transaction, not the handle. `is_compact_database` used that form and is called
once per candidate file by `_find_local_database`, so it leaked a connection per
call — harmless on Linux, a held file lock on Windows, where it would have
blocked `redownload=True` from replacing the very file it had just probed. The
three helpers that read a database now close explicitly in a `finally`.

### 11.5 Validation

- `pytest tests/` — **1008 passed, 34 skipped, 0 failed** (5m51s). The three
  PubChem failures §23.7a recorded on 2026-09-20 now pass; that IP is evidently
  no longer throttled.
- `tests/test_chembl_compact.py` — 34 offline tests plus one `slow` test that
  compacts the real release and compares 200 compounds, skipped when no full
  database is present.
- All 15 public `CheMBL` methods compared full vs extract on ChEMBL 36, plus
  `get_compound` over **3 000 random compounds** — 0 differences.
- `Search("cas")` over eight CAS numbers with all four sources, once against the
  full ChEMBL and once against the extract: identical `name`,
  `canonical_smiles`, `InChI`, `InChIKey`, `molecular_mass`, `confidence`,
  `n_source_support`, `source` and `match_method`, with ChEMBL voting
  (`n_source_support` up to 4) in both.

### 11.6 Still open

- The 30 GB original is **not** deleted unless asked. On this machine both files
  now exist; `CheMBL.compact(remove_source=True, force=True)` reclaims the
  27.14 GB whenever you want it.
- ~~Step 2 (§9.5, the `search_by_name` `UNION` rewrite)~~ — **done, §12**.
- ~~§9.3's ordering instability~~ — **done, §12**. The tests added in step 1
  still compare sets rather than lists, which remains correct.
- §4.13 — passing one client to `Search` disables the others — was found doing
  this work and is not fixed.
- `keep_inchi=False` is implemented and tested but unmeasured on real data at
  this commit; §9.6's 1.55 GB comes from the prototype.
- Routes 2 and 3 of §9.7 (build during download, stream the MySQL dump) are not
  started. `compact()` is deliberately the whole of step 1. ~~Route 2~~ —
  **done, §15**.

---

## 12. Landed on 2026-09-20 — step 2, the `search_by_name` rewrite (§9.5, §9.3)

§9.5 predicted 743 ms → 10 µs from a `UNION` rewrite, and §9.3 asked for an
`ORDER BY`. Both landed together, because they are the same method and the
second is only meaningful once the first has settled what "the matching rows"
means.

### 12.1 What landed

- `CheMBL.search_by_name` now issues a `UNION` of two single-table lookups —
  one on `molecule_dictionary(lower(pref_name))`, one on
  `molecule_synonyms(lower(synonyms))` — instead of an `OR` across a
  `LEFT JOIN`, and closes with `ORDER BY molregno LIMIT ?`.
- The `COMPACT_INDEXES` comment no longer says the two `lower(...)` expression
  indexes are useless to this method. They are now the whole point of it.
- 32 tests in `tests/test_chembl_name_search.py`, an example in
  `examples/chembl/name_search_demo.py`, and a `CHANGELOG.md` entry.

Both modes changed shape; only `exact=True` changed speed by orders of
magnitude, because a leading-wildcard `LIKE` is a scan whatever the query
shape.

### 12.2 Measured, on ChEMBL 36

200 real search terms — 150 preferred names and 50 synonyms, drawn at random
with a fixed seed — against both databases on this machine:

| | exact lookup, per call |
|---|---:|
| full release, old query | 764.4 ms |
| full release, new query | 220.9 ms |
| extract, old query | 559.7 ms |
| **extract, new query** | **10.1 µs** |

| | substring lookup, per call |
|---|---:|
| extract, old query | 665.2 ms |
| extract, new query | 274.8 ms |

So §9.5's headline holds: **764 ms → 10 µs, ~76 000×**, against the full
release that is today's default. What §9.5 did not say is that the rewrite
pays on a full release too — 764 ms → 221 ms, 3.5×, with no indexes at all —
because removing the join removes the automatic covering index SQLite was
building on `molecule_synonyms` for every call. A user who has not run
`compact()` still gets most of an order of magnitude.

The plan on the extract is now what §9.5 wanted:

```
MERGE (UNION)
LEFT   SEARCH molecule_dictionary USING INDEX ix_md_pref_lower (<expr>=?)
RIGHT  SEARCH molecule_synonyms  USING INDEX ix_ms_syn_lower  (<expr>=?)
```

against the old

```
SCAN md USING INDEX ix_md_molregno
BLOOM FILTER ON ms (molregno=?)
SEARCH ms USING AUTOMATIC COVERING INDEX (molregno=?) LEFT-JOIN
```

### 12.3 Equivalence

Compared as *sets* of molregnos with a limit high enough that neither query
truncates — because the old query's truncation is exactly the thing that is not
reproducible:

- 200 exact terms — **0 differences**;
- 25 substring fragments — **0 differences**;
- 50 exact terms checked for sortedness — **0 unsorted results**.

### 12.4 One semantic difference, and why it does not matter

The `LEFT JOIN` drove everything from `molecule_dictionary`, so a synonym row
whose `molregno` is not in `molecule_dictionary` could never match. The synonym
arm of the `UNION` has no such anchor and would return that orphan molregno.

`get_compound` selects `FROM molecule_dictionary`, so an orphan yields `None`
and is dropped from the result either way. ChEMBL has no orphans in practice —
the 200-term comparison above would have found them — but the outcome is
identical by construction rather than by luck, which is why no extra filter was
added.

### 12.5 The tests

`tests/test_chembl_name_search.py` builds a miniature ChEMBL of 200 compounds,
laid out for the cases that matter rather than for realism: 40 compounds share
one preferred name so that `limit` has to *choose*, rows are inserted in
descending `molregno` so a sorted result cannot be an accident of physical row
order, and one name is one compound's `pref_name` and a different compound's
synonym so the two arms have to merge.

Three things are pinned that were not pinned before:

- the **legacy query is kept in the test file** and the rewrite is asserted to
  find the same compounds, per term, as the shape it replaced;
- a truncated result is asserted **stable across a `VACUUM`**, which is the
  concrete form of §9.3's instability;
- the **query plan itself** is asserted, via SQLite's trace hook, so a future
  edit that quietly reintroduces a join or a scan fails a test rather than a
  benchmark nobody runs. The hook reports statements with parameters already
  substituted, which is what makes the recovered SQL plannable as it stands.

### 12.6 Validation

- `pytest tests/` — **1068 passed, 34 skipped, 0 failed** (8m01s).
- `examples/chembl/name_search_demo.py` run against the real extract.

One process note, since it will recur at every step: `uv run pytest` re-resolves
`uv.lock`, pulling in the whole chebifier/`chebai` optional extra — 1 883 added
lines, CUDA wheels included — which has nothing to do with the change under
test. The lockfile was reverted here. It is worth deciding, before step 3, that
either the lock is committed once deliberately or the test command stops
writing it.

### 12.7 Still open

- §4.13 — passing one client to `Search` disables the others — remains unfixed.
- ~~Step 3~~ — **done, §13**. Step 4 (the dataset manager) is next.

---

## 13. Landed on 2026-09-20 — step 3, `datasets.py` (§4.2)

§4.2 asked for one `download_file`: streamed, `Range`-resumed against a `.part`
file, checksum-verified when the server publishes one, atomic rename on
success, one progress bar, one logger — and five call sites collapsing to five
one-line calls. That is what landed, with one addition §4.2 did not ask for and
one it did not anticipate.

### 13.1 What landed

- `src/provesid/datasets.py`: `download_file`, `read_checksum`, `md5_of_file`,
  `DownloadError`, and the `CHUNK_SIZE` / `PART_SUFFIX` constants.
- All five download sites migrated: `pubchem.py`, `comptox.py`, `zeropm.py`,
  `chembl.py`, `chebi.py`.
- 26 tests in `tests/test_datasets.py`, an example in
  `examples/datasets/resumable_download_demo.py`, and two `CHANGELOG.md`
  entries.
- `DownloadError`, `download_file`, `md5_of_file` and `read_checksum` exported
  from `provesid`.

`DownloadError` subclasses `ServiceError`, so the whole family stays catchable
through one base. The separate class earns its place on recovery rather than on
taxonomy: an API call that fails is retried or abandoned, while a download that
fails has usually left a resumable `.part` file on disk.

### 13.2 The addition: `verify`, and checking before the rename

§4.2 listed retry, resume and checksum. Reading the five sites turned up a
fourth defect it did not mention, and it is the one with teeth.

Only `pubchem.py` validated the downloaded file *before* moving it into place.
`comptox.py` renamed first and then looked for the `chemicals` table, so a
failed check left the broken file exactly where the working database had been.
`zeropm.py` had the same shape and then deleted the database, leaving nothing.

So `download_file` takes a `verify` callback, handed the finished `.part` file,
and the rename happens only after it returns. All five sites now behave the way
the best of them did. A download can no longer destroy a working database, and
the five inconsistent cleanup paths became one.

### 13.3 What the tests found: a resumed download restarts at a chunk boundary

The resumption tests were written asserting that a transfer cut at 40 000 bytes
resumes at `bytes=40000-`. They failed, resuming at zero, and the reason is
worth recording because it is a property of the module rather than a bug in it.

`iter_content` yields whole chunks. A connection that breaks mid-chunk raises
before that chunk is handed over, so those bytes never reach the disk and the
next request asks from the last *complete* chunk. With the 1 MB default that
costs at most a megabyte per interruption — nothing against 5.8 GB — but it
means the resumed offset is a multiple of the chunk size, and a test that
assumed otherwise was testing urllib3's buffering rather than this module.

The tests now pin the real property (`resumed_from % chunk_size == 0`, and
progress is monotonic across repeated interruptions) and say why in the
docstring.

A second finding from the same tests: a file rejected by `verify` or by its
checksum must be **deleted**, not kept. It downloaded completely, so resuming
it would finish instantly and fail the same check again — an endless loop over
a 5.8 GB file. A transfer that merely stopped is kept, which is the whole
point.

### 13.3a A partial file has to say where it came from

Resumption has a hole that §4.2 does not mention: a `.part` left over from a
different release — or from a different dataset that happens to share a
destination — would be resumed, splicing two files together into something
that looks plausible and is not. A checksum catches it, but only PubChem's FTP
mirror publishes one; the four Zenodo and EBI downloads have no backstop at
all.

So `download_file` writes the URL to a `<dest>.part.source` marker and resumes
only a partial that came from the same URL. Anything else is discarded, with a
line in the log saying so. The marker is removed when the file lands.

### 13.4 The five sites

| Module | Before | After |
|---|---|---|
| `pubchem.py` | 60 lines, verified before rename | one call + a `verify` |
| `comptox.py` | 45 lines, verified *after* rename | one call + a `verify` |
| `zeropm.py` | 45 lines, verified after rename, deleted on failure | one call + a `verify` |
| `chembl.py` | 40 lines, then extraction inline | one call; extraction moved to `_extract_database` |
| `chebi.py` | 50 lines, then gunzip inline | one call; gunzip kept, gzip's CRC still the truncation check |

`chembl.py` and `chebi.py` keep their archives beside the destination rather
than in a temporary directory, so an interrupted download's `.part` file is
found and resumed next time. `chembl.py`'s inline extraction moved into
`_extract_database`, which left `download_database` short enough to read.

ChEMBL needed one more change to keep the same contract as the other four.
`verify` covers a file that `download_file` puts in place, but ChEMBL's
database is written by the *extraction*, which renamed straight onto
`db_path` — so a `force=True` re-download that failed to extract destroyed the
release already on disk. The archive now extracts to `<db_path>.incoming`,
answers a query there, and is moved into position only then.

`requests` and `tqdm` are no longer imported by `comptox.py` or `zeropm.py` at
all.

### 13.5 Verification

- `tests/test_datasets.py` runs a real HTTP server on localhost rather than
  stubbing `requests`, because the behaviour that matters is in the `Range`
  request and the response status — the parts a stub would have to fake, and so
  the parts a stub would let us get wrong. The server can be told to answer a
  status, hang up after N bytes, ignore `Range`, or serve a corrupt body.
- What is pinned: resumption from a leftover `.part`, resumption after a
  dropped connection, convergence over three successive drops, a partial from
  another URL or with no marker being discarded, a server that ignores
  `Range`, HTTP 416, retry on 503/502, a finite budget, a fatal 404 that is not
  retried, checksum match and mismatch, `read_checksum` against a
  coreutils-format sidecar, and that no failure path ever writes to the
  destination. 29 tests.
- One existing test had to be rewritten:
  `test_download_database_force_parameter` patched
  `provesid.zeropm.requests.get`, which no longer exists. It was worth
  rewriting on its own terms — it wrapped the call in `except Exception: pass`,
  so it asserted nothing at all. It now patches `download_file`, checks the URL
  and destination it is handed, and a companion test exercises the `verify`
  callback against a file that is not SQLite.

### 13.6 Validation

- `pytest tests/` — **1098 passed, 34 skipped, 0 failed** (7m39s).
- `examples/datasets/resumable_download_demo.py` run end to end: a dropped
  connection resumed at a chunk boundary, a leftover `.part` resumed, a
  `.part` from another URL discarded, a published checksum fetched and matched,
  and a rejected file leaving the previous good one in place.

### 13.7 Still open

- The real datasets were not re-downloaded to test this end to end; that costs
  ~9 GB and the local server exercises the same code paths. The first real
  download will be the proof.
- §8's `pubchem_ftp.py` builder (step 9) is written on top of this, and is
  where `checksum_url` finally pays for itself: PubChem publishes an `.md5`
  beside every FTP file, so every one of the seven source files is verified.
- Step 4 — `datasets.status/plan/fetch/remove` and `Search(datasets="present")`
  — is next, and is what stops a first run downloading 32 GB without asking.
  **Landed; see §14.**

---

## 14. Landed on 2026-09-20 — step 4, the dataset manager (§4.1)

§4.1 called the 32 GB first run "the largest user-friendliness defect in the
package". It asked for four functions and one keyword argument. Both landed,
with two additions the section did not anticipate and one deliberate omission.

### 14.1 What landed

- `src/provesid/datasets.py` grew a second half: the `Dataset` dataclass, the
  `DATASETS` registry, `DEFAULT_DATASETS`, `MissingDatasetError`, and
  `status`, `plan`, `fetch`, `remove`, `require`, `missing`, `is_present`,
  `dataset_files`, `dataset_names`, `data_directory`, `fetch_command` and
  `human_bytes`.
- `Search(datasets=...)` with the three policies of §4.1's table, defaulting
  to `"present"`.
- 40 tests in `tests/test_dataset_manager.py`, an example in
  `examples/datasets/dataset_manager_demo.py`, a §8 in `docs/quickstart.md`,
  and two `CHANGELOG.md` entries.
- `datasets`, `DATASETS`, `Dataset`, `MissingDatasetError` and `human_bytes`
  exported from `provesid`.

The module is now "the bulk datasets: what they are, and how they get here" —
`download_file` moves the bytes, the registry says which bytes and whether
they are already on disk. Keeping them together is what lets `status` count a
leftover `.part` against the dataset it belongs to.

### 14.2 The registry cannot import the clients, and should not

The obvious home for "how big is ChEMBL, and where does it live" is
`CheMBL` — and it is the wrong one, because constructing a client is exactly
the act the manager exists to avoid, and because all five client modules
import `datasets` for `download_file`. So the facts live in a table of frozen
dataclasses, and `fetch` imports the client it needs through
`importlib.import_module` at the moment it needs it.

The cost is that two tables now describe the same five sources: `DATASETS`
here and `Search._ALL_SOURCE_KEYS` there. A test pins them equal, because a
dataset added to one and not the other would be unreportable by `status` or
unreachable by `fetch`, with nothing to say so.

### 14.3 Measured sizes, and the peak nobody counts

The sizes in the registry were measured from the copies on this machine rather
than taken from §4.1's table, which rounded:

| Dataset | Download | Installed |
|---|---:|---:|
| PubChem | 2.2 GiB | 2.2 GiB |
| CompTox | 816.6 MiB | 816.6 MiB |
| ChEBI | 250 MiB (gz) | 954.2 MiB (SDF + 74.5 MiB index) |
| ChEMBL 36 | 5.7 GiB | 27.7 GiB |
| ZeroPM | 438.7 MiB | 438.7 MiB |
| **four default sources** | **8.9 GiB** | **31.6 GiB** |

Two of those rows are worth the trouble they cost.

**ChEBI's index is a fifth of ChEBI.** It is built on first use rather than
downloaded, so a registry that listed only the downloaded file would understate
ChEBI by 74 MiB and `remove` would leave it behind. Derived files are therefore
a separate `extras` tuple: counted in `status`, deleted by `remove`, and never
proof that the dataset is installed.

**A laptop with 28 GB free still cannot install ChEMBL.** The archive sits
beside the database it extracts into, so the worst moment needs 33.4 GiB, not
27.7 GiB. `plan()` reports that as `attrs["peak_bytes"]` — everything else
already installed, plus the largest single transient overhead, since the
datasets install one after another rather than at once. For the four default
sources it is 37.3 GiB against 31.6 GiB installed. §4.1 does not mention this,
and it is precisely the failure that wastes an afternoon: 5.8 GB transferred,
then no space to unpack it.

### 14.4 Two things §4.1 did not ask for

**A leftover `.part` is reported without being a dataset.** An interrupted
2.2 GB download leaves 2.2 GB on disk that no source can read. `status` counts
it (it explains a full disk after a cancelled download) and `present` stays
False (it is not the dataset). Both halves matter: treating it as installed
would be worse than not reporting it at all.

**`status` names the file a client would actually open.** With `chembl_36.db`
and `chembl_36_provesid.db` side by side, `CheMBL` opens the extract, so a
table naming the 27.7 GB original would be describing a database nobody will
open. `_preferred_file` mirrors `CheMBL._find_local_database`'s rule — newest
release, extract at equal releases — but reads it off the filename rather than
opening each candidate, so a status table stays a directory listing. The
divergence, documented where it is implemented: an extract renamed to hide its
`_provesid` suffix is misreported here and opened anyway by `CheMBL`, which
checks the file itself.

### 14.5 `redownload=True` had to become an error

`Search(redownload=True)` means "fetch a fresh copy". Under the new default it
cannot, and the two plausible readings — ignore it, or treat it as an implicit
`datasets="auto"` — are both worse than refusing. Ignoring it hands back a
stale database with no indication; promoting it would let a keyword that used
to mean "re-download the copy I already have" silently authorise the 32 GB
first download that this step exists to prevent. So it raises `ValueError`
naming the two ways forward. Dev-principle §1 permits the break.

### 14.6 `"required"` raises in the constructor

Clients are constructed lazily on the first `search()`, which is where the
natural place to enforce a policy would be. `"required"` is checked in
`__init__` instead, so a batch script fails on the line that configured it
rather than an hour into a run. The check is a `glob` over one directory: no
client is constructed and nothing is downloaded, which is the whole point.

It skips any source whose client the caller passed in — that client has
already found its data, wherever it put it — and it does not demand ZeroPM
unless `use_zeropm=True` put it among the queried sources.

### 14.7 What this does not do

`ChebiSDF`, `CompToxID`, `PubChemID`, `ZeroPM` and `CheMBL` still default to
`auto_download=True` when constructed directly. `Search` no longer relies on
that default, which is where the 32 GB came from, but `CheMBL()` on a clean
machine still fetches 5.8 GB unasked. Changing the five client defaults is a
wider break than §4.1 asked for and belongs with step 6, which is already
opening those constructors for `close()` and context-manager support.

`fetch` installs by constructing the client, which for ChEBI also builds the
index — a few minutes with no progress bar. That is the right definition of
"installed" (the dataset is usable when `fetch` returns), but it means `fetch`
is not purely a download.

### 14.8 Verification

- `tests/test_dataset_manager.py`, 40 tests, every one against an explicit
  `data_dir` in `tmp_path`: the manager's job is to report on a directory, and
  a test that used the real one would report on whatever the developer happens
  to have installed. The files it creates are empty — the manager reads names
  and sizes, never contents.
- §7 asked for a test that constructs `Search` with no datasets present and
  asserts that nothing is downloaded and the error names the `datasets.fetch`
  call. Both exist: `test_present_downloads_nothing_and_reports_what_is_missing`
  runs a real CAS search in an empty directory with every module's
  `download_file` monkeypatched to fail the test if called, then asserts the
  directory is still empty and that each of the four warnings carries its own
  `provesid.datasets.fetch('...')`; `test_required_raises_in_the_constructor`
  covers the raising policy.
- What else is pinned: the registry against `Search._ALL_SOURCE_KEYS`, ChEMBL's
  peak exceeding its installed size, a `.part` counted but not present, the
  extract preferred over the release it came from, a newer release preferred
  over an older extract, ZeroPM's newest version chosen out of two, `plan`
  totals over a mixed present/missing set, `fetch` skipping what is present and
  passing `auto_download=True` for what is missing, `fetch` announcing the
  total before transferring, `remove` taking derived files and every ChEMBL
  release with it, `dry_run` deleting nothing, and each of the three policies
  reaching the clients as the right `auto_download` value.

### 14.9 Validation

- `pytest tests/` — **1135 passed, 34 skipped, 3 failed** (9m30s). The three
  failures are live PubChem 503s in `test_pubchem.py` and `test_pubchemview.py`;
  they fail identically on the unmodified tree (verified with `git stash`).
- `mkdocs build --strict` — clean.
- `examples/datasets/dataset_manager_demo.py` run end to end: the status table
  for this machine, the 8.9 GiB / 31.6 GiB / 37.3 GiB plan for a clean one, a
  real CAS search in an empty directory that downloaded nothing and reported
  four missing sources, and `"required"`'s message.
- `Search("cas").search("50-00-0")` against the real datasets still returns
  formaldehyde at confidence 0.9 from all four sources.

### 14.10 Still open

- The five client defaults (§14.7), which step 6 should take.
- `fetch` was exercised against stub clients, not against the real Zenodo and
  EBI downloads — that costs ~9 GB, and the transport underneath it is what
  `tests/test_datasets.py` already covers against a real HTTP server. The first
  real install will be the proof.
- No `docs/api/datasets.md` yet; the module has none from step 3 either, and
  step 18 rebuilds the API reference.
- Step 5 — building the ChEMBL extract during download, and `source=` on
  `CheMBL` (§9.7) — is next. It is what turns the 31.6 GiB in §14.3 into
  ~6.5 GB for a user who never needed the full release in the first place.
  **Landed; see §15**, at 6.3 GiB.

---

## 15. Landed on 2026-09-20 — step 5, `CheMBL(source=...)` (§9.7 route 2)

§9.7 asked for "the `sqlite.tar.gz` route, unchanged except that
`download_database` finishes by building the extract and removing the
intermediate". That is what landed. Step 1 (§11) could only shrink a release
that was *already* on disk, which left the worst case untouched: a machine that
had never had ChEMBL still installed 27.7 GiB and then had to be told to shrink
it. A first install now costs **2.42 GiB**.

### 15.1 What landed

- `CheMBL(source=...)`, with `SOURCES = ("sqlite", "full")` as a documented
  class attribute and `"sqlite"` as the default. `"sqlite"` compacts the
  release as the last step of the download and deletes it; `"full"` is the
  behaviour of every release until now.
- `_compact_after_download`, called from `download_database`, and
  `_validate_source`, called from `__init__` before anything is fetched.
- `_UNIMPLEMENTED_SOURCES`, so that `source="mysql"` — route 3, step 11 — says
  it is not implemented rather than reading as a typo.
- `_suggest_compacting_a_full_release`: one INFO line when a full release over
  1 GiB is opened, naming `compact`. The constructor does not compact what it
  finds (§15.4), so saying the option exists is the only thing left to do.
- `datasets.DATASETS["chembl"]` resized to what an install now leaves, plus
  `chembl_*.db.tmp` among its `extras`.
- 15 tests (14 in `tests/test_chembl_compact.py`, 1 in
  `tests/test_dataset_manager.py`), `examples/chembl/download_source_demo.py`,
  a `CHANGELOG.md` entry, and the size claims in `search.py`, `datasets.py`,
  `docs/quickstart.md`, `docs/api/chembl.md` and `examples/chembl/README.md`
  brought in line.

### 15.2 The archive is deleted before the extract is built, and that is the whole peak question

§9.7 called route 2's transient peak "exactly the laptop problem §4.1 is
about". The peak is decided by one line's position. Three files can exist
during an install — the 5.7 GiB archive, the 27.7 GiB release, the 2.42 GiB
extract — and the order they are removed in decides which two overlap:

| moment | archive | release | extract | total |
|---|---:|---:|---:|---:|
| extracting | 5.7 | 27.7 | — | **33.4 GiB** |
| compacting, archive already gone | — | 27.7 | 2.4 | 30.1 GiB |
| compacting, archive kept | 5.7 | 27.7 | 2.4 | 35.8 GiB |

So `os.remove(archive_path)` moved to *before* the compaction rather than
after it, and the result is that **route 2 needs no more free disk than the
full route always did**: the peak stays the extraction moment, 33.4 GiB, and
what changes is only what is left at the end. Had the archive been kept until
the download function returned — the obvious ordering, and the one that reads
more naturally — the new default would have raised the requirement by 2.4 GiB
for every user, including those with the disk to spare.

### 15.3 A compaction that fails does not throw away the download

`_compact_after_download` logs a warning and returns; it does not raise. By the
time it runs, the release is in place and answers every query, so a failure
costs disk rather than function — and the alternative is discarding a 5.8 GB
transfer over a step that `CheMBL.compact(remove_source=True)` repeats in one
call. The warning names that call and the path the release is at, and
`db_path` stays on the release, so the object the constructor returns is
usable either way.

Two things make that safe rather than merely convenient. `compact` already
verifies the extract against its source before deleting anything (§11.2), so
"the compaction failed" never means "the release was deleted anyway". And it
builds into `<dest>.tmp` and `os.replace`s, so the `force=True` this path
passes cannot destroy a good extract from an earlier release of the same
number when the rebuild fails.

### 15.4 `source` describes a download, not a directory

`CheMBL(source="sqlite")` on a machine that already has `chembl_36.db`
downloads nothing and compacts nothing. That is deliberate and is pinned by a
test: deleting 27 GB is something to ask for (`compact(remove_source=True)`),
not something a constructor should do because a keyword argument has a default.
The same test forbids any network access, since "already installed" must also
mean "not re-fetched".

What is left is the hint of §15.1 — and it is silent below 1 GiB, which keeps
it out of the miniature databases the tests build and, more to the point, stops
it firing where there is nothing to reclaim.

### 15.5 What the registry had to say

§14.3 measured ChEMBL at 27.7 GiB installed and called its 33.4 GiB peak "the
peak nobody counts". With the extract as the default, the two swap roles: the
peak is now more than ten times the installed size, and it is the only number
that decides whether the install succeeds.

| | before | after |
|---|---:|---:|
| ChEMBL download | 5.7 GiB | 5.7 GiB |
| ChEMBL installed | 27.7 GiB | **2.42 GiB** |
| ChEMBL peak | 33.4 GiB | 33.4 GiB |
| four default sources, installed | 31.6 GiB | **6.3 GiB** |
| four default sources, peak | 37.3 GiB | 37.3 GiB |

`resident_bytes` is the measured size of `chembl_36_provesid.db` on this
machine (2 599 391 232 bytes), not §9.3's rounded 2.60 GB. A test pins the
registry's ChEMBL row against `CheMBL.__init__`'s `source` default, because
these sizes are shown to a user *before* a 5.8 GB transfer: if the default ever
moves back, the number that made them agree to it would be a lie.

### 15.6 Verification

- `tests/test_chembl_compact.py` grew two fixtures and four classes. The tests
  do not stub the download path: a miniature release is packaged as a real
  `chembl_36/chembl_36_sqlite/chembl_36.db` tar.gz, `download_file` is replaced
  by a copy, and everything after it — extraction, validation, compaction,
  verification, deletion — runs for real.
- What is pinned: that the default leaves only the extract and no `.tar.gz`,
  `.incoming` or `.tmp` behind; that the extract answers and carries its
  provenance; that `"full"` keeps the release and builds no extract; that both
  routes request the same archive; that a release already on disk is neither
  re-fetched nor compacted; that `redownload=True` rebuilds over an existing
  extract; that a failed compaction keeps the release, keeps it queryable and
  names `compact` in the warning; that a bad `source` raises before the network
  is touched; and that `"mysql"` says it is unimplemented.
- The compaction hint has its own three tests, including that an extract is
  never told to compact itself.

### 15.7 Validation

- `pytest tests/` — **1150 passed, 34 skipped, 3 failed** (10m33s). The three
  failures are the live PubChem 503s of §14.9, unchanged and unrelated:
  `test_error_handling_invalid_cid`, `test_malformed_property_names` and
  `TestPubChemView::test_error_handling`, each of which asserts on an error
  message and gets `HTTP 503` from PubChem instead.
- `mkdocs build --strict` — clean.
- `examples/chembl/download_source_demo.py` run end to end: it builds a
  miniature release, serves it over a real HTTP server on localhost, installs
  it both ways, prints what each route leaves on disk (0.18 MB against
  3.76 MB, 95.3% smaller), shows the two route errors, and ends with
  `datasets.plan("chembl")` for the real thing.
- `datasets.plan()` on this machine now reports 8.9 GiB / 6.3 GiB / 37.3 GiB
  for the four default sources, against 8.9 / 31.6 / 37.3 before.
- `Search("cas").search("50-00-0")` against the real datasets still returns
  formaldehyde, `WSFSSNUMVMOOMR-UHFFFAOYSA-N`, at confidence 0.9. Nothing in
  this step touches a query path: the changes are the constructor's argument
  check, the tail of `download_database`, one log line, and sizes in a table.

### 15.8 Still open

- ~~**Route 3 (step 11), `source="mysql"`.**~~ **done, §20.** 2.1 GB transferred instead of
  5.7 GiB and no 27.7 GiB intermediate at all, which is the version of this
  that a laptop with 20 GB free can actually run. §9.7 wants it verified
  against route 2 table by table; route 2 is now the thing to verify against.
- The real 5.8 GB download has still not been run end to end (§13.7, §14.10).
  Every step of it is exercised, but against a miniature release.
- `keep_inchi` is not exposed on the constructor. §9.6's 1.55 GB variant is
  reachable only by compacting by hand afterwards, which is the right default
  surface for something that degrades `Search` output.
- The five client defaults (§14.7) are unchanged: `CheMBL()` on a clean machine
  still downloads without being asked — 2.42 GiB now rather than 27.7 GiB.
  Step 6.

---

## 16. Landed on 2026-09-20 — step 6, `close()`, `with` and threads (§4.3)

§4.3 asked for "a `close()`, `__enter__`/`__exit__`, and either
`check_same_thread=False` with a documented 'read-only, safe to share'
contract or a thread-local connection", with `__del__` kept as a backstop and
`Search` closing the clients it constructed. All of that landed, as one mixin
rather than four copies. Of the two options §4.3 offered, only the second was
available — §16.3.

### 16.1 What landed

- `src/provesid/sqlite_client.py` (406 lines): `SQLiteClient`, the mixin, and
  `DatabaseClosedError`. `_open_database` is the one call a client makes;
  `conn`, `cursor`, `closed` and `db_file` are properties; `close`,
  `__enter__`, `__exit__` and `__del__` are the lifetime. Both are exported
  from `provesid`.
- `PubChemID`, `CompToxID`, `ZeroPM` and `CheMBL` inherit it. Each lost its
  `self.conn = sqlite3.connect(...)` and its `__del__`; none of their 186
  `self.conn` / `self.cursor` call sites changed.
- `Search.close`, `__enter__`, `__exit__`, and `_owned_clients` — the source
  keys whose client this instance built, which is what it may close (§16.4).
  A closed `Search` raises from `_ensure_clients` rather than searching on.
- `_adopt_connection`, for an instance built without `__init__` (§16.6).
- 50 tests in `tests/test_sqlite_lifecycle.py`,
  `examples/sqlite_clients/` (a demo and a README), `docs/api/sqlite_clients.md`
  in the nav and the API index, a `docs/quickstart.md` section, and a
  `CHANGELOG.md` entry.

### 16.2 The hazard is `self.cursor`, not `self.conn`

§4.3 described the threading problem as a connection created in one thread
being used in another. That is the *symptom*; it is also the benign half. The
four clients hold 186 references to `self.conn` and `self.cursor` between
them, and 166 of those are `self.cursor` — 135 in `zeropm.py`, 31 in
`chembl.py` — written as two statements:

```python
self.cursor.execute("SELECT ... WHERE cas = ?", (cas,))
result = self.cursor.fetchone()
```

That is why `check_same_thread=False`, the first of §4.3's two options, is the
wrong one. It removes the check and leaves every thread sharing one cursor, so
two threads interleaving those two statements read each other's rows. The
exception it would silence is the thing protecting the caller from a wrong
answer, and a wrong answer in a chemical identifier resolver is a far worse
outcome than a traceback. So the unit of thread affinity has to be the cursor,
not just the connection, and `SQLiteClient` hands out one of each per thread.

`check_same_thread=False` is still passed — but for a different reason, and
only because per-thread connections make it safe: `close()` runs on whichever
thread owns the object, and it has to be able to close the connections the
workers opened.

### 16.3 "Read-only, safe to share" was not available

§4.3's other option assumed these are read-only databases. `ZeroPM` is not:
`create_indexes` and `create_view` both execute DDL and commit
(`zeropm.py:1663`, `zeropm.py:2001`). Opening with `mode=ro` would have broken
them, and a contract that says "read-only" while one client writes is worse
than no contract. Per-thread connections need no such promise — SQLite
serialises the writes itself — which is the second reason the choice was made
the way it was.

### 16.4 `Search` closes what it built, and nothing else

`_ensure_clients` records a source key in `_owned_clients` only when it
successfully constructed that client. A client passed to the constructor is
never in the list, so `Search.close()` leaves it open: it belongs to the
caller, who may still be using it, and closing someone else's database on the
way out of a `with` block is not a service. Both halves are pinned by tests.

A closed `Search` raises `DatabaseClosedError` from `_ensure_clients` rather
than rebuilding its clients. Rebuilding would be the surprising choice — a
`with` block that quietly re-opens 6.3 GiB of databases after it ends is not
what the block said.

### 16.5 Threads make the pool work; they do not make it fast

Worth recording because the obvious reading of this step is "PROVESID is now
faster on many cores", and it is not. Measured on this machine against the
real `pubchem_id.db`:

| workload | serial | 4 threads | 8 threads |
|---|---:|---:|---:|
| 5 000 `get_by_cas` | 0.29 s | 4.80 s | 14.65 s |
| 400 `get_by_cas`, each + 20 ms wait | 8.52 s | — | 1.17 s |

A warm local lookup takes tens of microseconds, which is less than the GIL
handoff around the `sqlite3` call costs, so a pool over nothing but lookups
loses badly. This is not something the mixin introduced: a plain
`sqlite3.connect` per thread, with no PROVESID in the picture, measures the
same (0.06 s serial against 2.94 s on eight threads), and chunking the work
does not help. A pool pays when each item also waits on something, which is
the case a user actually reaches for it in — a lookup feeding a request, a
file read, an RDKit call — and there it is 7.3×.

So the value of this step is that the pool no longer *fails*, plus the two
lifetime features. The docs, the example and the CHANGELOG all say so rather
than leaving the reader to infer a speedup.

### 16.6 Two things §4.3 did not anticipate

**`conn` became a property, which broke two test helpers.**
`tests/test_search_new_methods.py` builds a client with `object.__new__` and
then assigns `obj.conn = mock_conn`, which a read-only property refuses. A
setter would have been the easy fix and the wrong one: assigning to `conn`
leaves the replaced connection open and off the registry, so nothing would
ever close it. `_adopt_connection` is that assignment with the bookkeeping
kept — it registers the connection, so `close()` still reaches it — and it
also makes an in-memory database a supported way to build a client.

**`close()` has to work on a half-built object.** `__del__` calls it, and
`__del__` runs on an instance whose constructor raised before
`_open_database` — a missing database file, for instance, which is an
ordinary path with `auto_download=False`. Every attribute access in `close`
and `closed` therefore goes through `getattr(self, ..., default)`, and
`closed` reports `True` for an object that owns nothing.

### 16.7 Verification

- `tests/test_sqlite_lifecycle.py`, 50 tests in two halves. The first builds a
  three-row database in `tmp_path` and exercises the mixin through a
  ten-line `TinyClient`, so closing, re-closing, `with`, the per-thread
  handles, the `__del__` backstop, the row factory and an adopted connection
  are all tested without an installed dataset. The second runs the same
  contract against all four real clients, parametrised and skipped per
  dataset, plus five `Search` tests.
- Two of the threading tests needed a `threading.Barrier` to be tests at all:
  without one, `ThreadPoolExecutor` finishes each task before starting the
  next worker, so four submissions run on one thread and "each thread gets its
  own connection" passes trivially against a single connection.
- What is pinned: that a closed client names itself and its file in the error;
  that the error is a `RuntimeError` subclass, so existing handlers still
  catch it; that `close` is idempotent and reaches connections opened on other
  threads; that a worker thread sees `DatabaseClosedError` rather than
  sqlite3's own; that `__exit__` closes while letting the exception through;
  that `ZeroPM`'s tuple rows survive; that a write on one thread is visible on
  another; and that `conn` has no setter.

### 16.8 Validation

- `pytest tests/` — **1200 passed, 34 skipped, 3 failed** of 1 237 collected
  (8m08s). The three failures are the live PubChem 503s of §14.9 and §15.7 —
  `test_error_handling_invalid_cid`, `test_malformed_property_names` and
  `TestPubChemView::test_error_handling` — unchanged and unrelated. An
  intermediate run also failed the six `test_search_new_methods.py` mock
  tests of §16.6; those are fixed and pass.
- `mkdocs build --strict` — clean, after `__exit__`'s three parameters were
  annotated: griffe warns on an unannotated parameter, and `--strict` turns
  that into a failure.
- `examples/sqlite_clients/connection_lifetime_demo.py` run end to end: it
  builds a ten-row database, shows `with`, `close`, the error message, 2 000
  lookups on eight threads over nine connections, and `close()` reaching all
  five connections a barriered pool left open, then queries the real
  PubChemID, CompToxID and CheMBL.
- `Search("cas").search("50-00-0")` against the real datasets still returns
  formaldehyde, `WSFSSNUMVMOOMR-UHFFFAOYSA-N`, at confidence 0.9.

### 16.9 Still open

- `ChebiSDF` is not an `SQLiteClient` — it holds a pickled index in memory
  rather than a file handle — so `Search.close()` skips it. It is the one
  source whose memory a `with` block does not hand back.
- Nothing calls `close()` on the client an instance-level
  `download_database(force=True)` is about to replace. On Windows that
  `os.replace` still fails against a live object; the remedy is now available
  (`close()` first) but not automatic. Wiring it would mean a `reopen()` on
  the mixin, which is a step of its own.
- The five client defaults (§14.7, §15.8) are still unchanged: `CheMBL()` on a
  clean machine downloads without being asked. That is `Search`'s default
  today, not the clients'.

---

## 17. Landed on 2026-09-20 — step 8, `cache.py` (§4.7)

### 17.0 Step 7 is postponed, not dropped

`classyfire.py` stays as it is for now. The decision is the user's and the
reason is throughput: every method of that client blocks on a service whose
API is as slow as it gets, so the module is the least rewarding thing in the
queue to touch and the last that anything else depends on. §4.8's finding is
unchanged and the work is unchanged — every method raises
`ServiceUnavailableError` citing February 2023, the URL construction stays
commented, the seventeen status-code tests assert the raise, and `docs/` and
`README.md` stop listing it beside the working services. Nothing in steps 8–20
blocks on it; §4.8 sequenced it early only because it is small. It keeps its
entry in `CACHE_SERVICES` in the meantime, because the cache does not care
whether a service answers.

### 17.1 What landed

§4.7 asked for three things in order of user impact, and all three landed in
one commit because they are all one file.

- `provesid.utils.user_cache_path`, the cache twin of `user_dataset_path`:
  `platformdirs.user_cache_dir`, `PROVESID_CACHE_DIR` as the override.
- `cache.py` 691 → 884 lines, and **29 module-level functions → 14**. The
  eighteen removed are the seven `clear_<svc>_cache`, the seven
  `get_<svc>_cache_info`, `export_service_cache`, `import_service_cache`,
  `get_all_service_cache_info` and `clear_all_service_caches`. The three added
  are `get_service_cache`, `get_all_cache_info` and `_require_known_service`.
- `CACHE_KEY_VERSION`, `CACHE_SERVICES`, `get_service_cache` and
  `get_all_cache_info` exported from `provesid`.
- 18 tests in `tests/test_cache_layout_and_versioning.py`;
  `examples/cache/` (the new demo, the old `cache_demo.py` moved in, a README);
  `docs/advanced_caching.md` rewritten in four places; a `CHANGELOG.md` entry.

### 17.2 The directory was the whole of defect 1, and it was one line

`tempfile.gettempdir()/provesid_cache` → `user_cache_path()`. What makes it
worth a section is what it changes about everything *else* in the file: with
the cache under `/tmp`, a stale entry, a wrong key and a test that called
`provesid.clear_cache()` all cost nothing, because the next reboot fixed them.
None of that is true now.

Three consequences had to be handled in the same commit:

- **`tests/conftest.py` sets `PROVESID_CACHE_DIR`.** Several tests call
  `provesid.clear_cache()` outright. Before, they deleted a directory that was
  about to be deleted anyway; now they would delete responses the developer
  paid network time for.
- **`examples/cache/cache_demo.py` sandboxes itself** for the same reason — it
  clears the cache and drops the size-warning threshold to 1 MB.
- **The versioning of defect 3 stops being optional.** A wrong-shaped entry is
  now permanent until someone clears the cache by hand.

There is no migration of entries from the old location. They were never
reliably there to migrate, which is the defect.

### 17.3 Laziness fell out of the parameterisation

`_service_caches` was a module-level dict literal of seven `CacheManager()`
constructions, evaluated at import, and `CacheManager.__init__` calls
`mkdir(parents=True, exist_ok=True)`. So `import provesid` created eight
directories — seven services and the global cache — on a machine that never
made a single call. Under `/tmp` nobody noticed. Under `~/.cache` it is litter.

`get_service_cache(service)` builds each manager on first use and memoises it.
Two callers of the same service still share one manager and therefore share
the in-memory tier, which is the property `@cached` relies on.

### 17.4 An unknown service has to raise

`_service_caches.get(service, _global_cache)` was the old lookup: a typo — or a
service name that used to exist — silently wrote to the global cache, where
`clear_cache(service=...)` for that name would never look, and
`get_cache_info` would report zero entries for a cache that was filling up
somewhere else. `get_service_cache` raises `ValueError` naming the known
services instead, and `@cached(service=...)` raises at *decoration* time rather
than on the first call, so a bad name is an import error rather than a runtime
surprise.

One exception, and it is deliberate: a name already present in
`_service_caches` counts as known even if it is not in `CACHE_SERVICES`. Four
existing tests and two new ones inject a throwaway `CacheManager` under a
synthetic name (`"probe"`, `"probe2"`, …) to test the decorator without
touching a real service directory. Rejecting those would have meant
monkeypatching `CACHE_SERVICES` in every one of them, which tests the patch
rather than the code.

### 17.5 Three versions, and why not one

§4.7 asked for "a version component in every key, bumped when a return shape
changes". One number would have done that, but every bump would discard every
entry of every service — for a change to one function's return dict. So the
key carries three, and the rule is to bump the narrowest that covers the
change:

| Bump | Retires | Lives in |
|---|---|---|
| `@cached(version=N)` | one function's entries | the decorator |
| `Client.CACHE_SCHEMA_VERSION` | one client's entries | `__cache_key__` |
| `CACHE_KEY_VERSION` | every entry, everywhere | `cache.py` |

The middle one already existed — §19.6 of the August plan added it to
`PubChemView`, `OPSIN` and `CASCommonChem` after `PropertyData` gained a field
— but it only reaches methods of a client that declares it. A module-level
`@cached` function has no `self`, so nothing versioned its key at all; that is
the loophole `version=` closes. `CACHE_KEY_VERSION` is the one to bump if the
key derivation itself changes, since that is not a shape change any narrower
number describes.

Retired entries are not deleted. Finding them would mean reading every pickle
in the directory, which is the cost the version exists to avoid; `clear_cache`
reclaims the space when the user wants it back.

### 17.6 One bug found while reading the file

`export_cache` skipped any entry whose value was `None`:

```python
value = self._load_from_disk(cache_key)
if value is not None:          # meant: is not _MISS
    export_data['cache_data'][cache_key] = value
```

`_load_from_disk` returns the `_MISS` sentinel for a missing entry precisely so
that a cached `None` is distinguishable from no entry (§17.4 of the August
plan). This line predates the sentinel and was never updated, so a function
that legitimately returns `None` had its entries dropped from every export.
Fixed, with a test that round-trips a cached `None`.

### 17.7 Verification

- `tests/test_cache_layout_and_versioning.py`, 18 tests in three groups. What
  is pinned: that the default directory is not under the system temp directory
  and is not the dataset directory; that `PROVESID_CACHE_DIR` moves both the
  root and the service subdirectories; that an entry written by one
  interpreter is read by the next **without the function running again**
  (a real subprocess pair, which is the only honest test of persistence);
  that the seven services get seven distinct directories, none of them the
  global one; that `clear_cache(service=...)` leaves the other services
  untouched and `all_services=True` does not; that a bad service name raises
  from `get_cache_info`, from `clear_cache` and from `@cached`; that the
  eighteen removed functions are gone; that the client `clear_cache` /
  `get_cache_info` methods still point at the right service; that each of the
  three versions changes the key and that bumping a function's version stops
  the previous shape being served; and the cached-`None` round trip of §17.6.
- Doctests: six examples in `cache.py` are marked `# doctest: +SKIP` because
  they clear caches, write files or mutate thresholds. Under `/tmp` that was
  harmless; it is not now.
- `examples/cache/cache_layout_and_versioning_demo.py` run end to end. It
  prints the sandbox directory beside the real platform default, fills two
  service caches and clears one, shows the `ValueError` for `'pubchemm'`, and
  runs the same function at `version=1` and `version=2` to show the second
  call re-running rather than reading the first's entry.

### 17.8 Validation

- `pytest tests/` — **1201 passed, 35 skipped, 19 failed** of 1 255 collected
  (10m02s). PubChem was down for the duration: all nineteen failures are live
  calls to `pubchem.ncbi.nlm.nih.gov` returning HTTP 503, four in
  `test_pubchem.py` and fifteen in `test_pubchemview.py`, and a bare `curl` of
  both the PUG and PUG-View endpoints returned 503 at the same moment. No
  offline test failed, and no test in either file failed for any other reason.
  This is a wider outage than the three flaky 503s of §14.9, §15.7 and §16.8,
  not a different one. Re-running with those two files excluded: **1189
  passed, 33 skipped, 0 failed** (11m05s).
- `mkdocs build --strict` — clean.
- `examples/cache/cache_demo.py` and the new demo both run; neither leaves
  anything in `~/.cache/provesid`.

### 17.9 Still open

- **`get_cache_size` is O(n) in the number of entries and is called every 100
  writes.** `docs/chebifier.md` already warns about this: one pickle per entry
  in one flat directory, re-globbed and re-`stat`ed on every hundredth write,
  which is O(n²) over a long run and serialises processes that share the
  directory. Parameterising the functions did not touch it. The fix is a size
  counter maintained in the metadata rather than a directory walk, or a
  two-level directory fan-out by key prefix.
- **Metadata is written every 100 operations, so `disk_entries` on disk is
  usually stale.** It is accurate within a process, because the dict is in
  memory; a process that exits between saves loses the tail of its metadata,
  though not the entries themselves.
- **No entry ever expires.** Chemical data rarely changes, which is the
  justification, but there is no TTL and no way to ask for one.
- **The `_memory_cache` is unbounded and never evicts.** A long run over
  100k structures holds every result in RAM as well as on disk.
- Nothing in this step touched `is_empty_result` / `skip_if`, which remain the
  per-function way to say "absence is not an answer".

---

## 18. Landed on 2026-09-21 — step 9, `pubchem_ftp.py` (§8)

§8.4 asked for one function that builds `pubchem_id.db` from a dated monthly
snapshot of PubChem's FTP site, and for `PubChemID(source="ftp")` as the
default route to it. That is what landed, with four departures, each argued
below: the molecular weight is *not* PubChem's (§18.3), the scope is a
trade and not a strict gain (§18.4), the descriptors stop being served from
disk even by an old Zenodo copy (§18.5), and there is no `scope=` or
`source="local"` (§18.6).

### 18.1 What landed

- `src/provesid/pubchem_ftp.py`: `build_pubchem_id_db`, `list_releases`,
  `resolve_release`, `extras_url`, `selected_files`, `molecular_weight`, and
  the data they run on --- `SOURCE_FILES`, `XREF_TYPES`, `ATOMIC_WEIGHTS`.
  `release="latest"` (default) resolves to the newest `Monthly/YYYY-MM-01/`
  by reading the directory listing; `"current"` is the rolling `Extras/`; a
  date names a snapshot and costs no request.
- `PubChemID(source=...)`, `SOURCES = ("ftp", "zenodo")`, `"ftp"` the default,
  checked before anything is fetched. Like `CheMBL(source=...)` (§15.4) it
  describes an acquisition, not a file: whatever is on disk is opened.
- `PubChemID.offline_properties`, per instance, from the columns the open
  database actually has; `provenance()`; `xrefs(cid)`.
- `get_by_cas_batch` / `get_by_smiles_batch` return the database's own
  columns instead of a hard-coded list naming the eight descriptors.
- `datasets.DATASETS["pubchem"]` sized for the FTP route, with the download
  directory and the `.tmp` among its `extras`, so `status` counts and `remove`
  deletes a kept or interrupted build.
- `scripts/build_pubchem_id_db.py` rewritten as an `argparse` wrapper for the
  Zenodo refresh; the CSV-and-regex original and its stale duplicate in
  `src/provesid/data/` deleted.
- 46 tests in `tests/test_pubchem_ftp.py`, `examples/pubchem/ftp_build_demo.py`,
  `docs/api/pubchem.md` (a new "The Local Database" section and the module's
  API), `CHANGELOG.md`, `scripts/README.md`, `src/provesid/data/README_PUBCHEM.md`.

**The schema** is §8.4's: `compounds` / `cas_numbers` / `synonyms` as before
minus the eight descriptor columns, plus `monoisotopicmass`, plus `xrefs`,
`provenance` and `provenance_files`. Two things §8.4 did not list:
`cidcdate` is kept, from `CID-Date.gz` (315 MB), because the column existed
and the file is small; it is ISO `YYYY-MM-DD` where the CSV had `YYYYMMDD`.
And `synonyms` now holds PubChem's *filtered* list verbatim, CAS-shaped
strings included --- the old builder moved anything matching the regex into
`cas_numbers`, which is how the 17 379 invalid CAS numbers got there.

### 18.2 How the build is shaped

One pass per file, one file at a time, as §8.4 said. Three details decided
the shape:

- **Scope is a Python `set`, not §8.4's sorted `array('l')`.** 1.43 M ints are
  about 75 MB as a set against 11 MB as an array, and a set lookup is one hash
  where bisection is twenty comparisons, on each of ~150 M lines per large
  file. Readability and speed both favour the set, and 75 MB is not a laptop
  problem.
- **Rows go straight to SQLite and are deduplicated there.** Collecting the
  1.46 M `(cid, cas)` pairs in a dict first cost ~250 MB; `INSERT` as read and
  one `DELETE ... WHERE rowid NOT IN (SELECT MIN(rowid) ... GROUP BY ...)`
  costs nothing. The identifier file is read twice --- scope first, then the
  cross-references filtered by it --- because nothing guarantees it is sorted,
  and 98 MB twice is seconds.
- **Columns are filled by `UPDATE ... WHERE cid = ?`**, one file at a time, on
  the `INTEGER PRIMARY KEY`. Lines are handled as bytes and decoded only when
  they survive the scope filter: ~1% of a large file does.

A failed build deletes its `.tmp` and leaves any existing database untouched;
a successful one renames it into place and removes the emptied download
directory. `keep_downloads=True` keeps the files, and a later build reuses a
kept file after checking it against the published MD5 --- a damaged one is
fetched again (tested). Two bugs the tests found before any real data did: the
database's own directory was never created (it only existed because the
default download directory is inside it), and the empty-directory cleanup
would have climbed out of a `download_dir` placed outside the data directory.

### 18.3 The molecular weight cannot be PubChem's

§8.4 step 3 said "PubChem's own molecular weight is exactly this" --- the sum
of standard atomic weights over the formula. Tested against the 1.58 M
weights in the shipped database, it is not, for two separate reasons.

**The atomic weights.** RDKit's table (IUPAC 2013-ish, `C 12.011`,
`H 1.008`) reproduces PubChem at any precision for 74% of compounds. IUPAC
2005 (`C 12.0107`, `H 1.00794`) does better for C, H, N, O and P. Per
element, choosing whichever of the two reproduces more PubChem weights, and
trying the obvious alternatives for the rest, settles on IUPAC 2005 with twelve
elements moved: S 32.067, Se 78.971, Ge 72.630, Mo 95.95, B 10.812,
Yb 173.045, Hg 200.592, Hf 178.486, Cd 112.414, Cl 35.45, Si 28.085,
I 126.904. That is `ATOMIC_WEIGHTS`, and its comment says it was fitted.

**The rounding.** PubChem reports weights to between zero and three decimals
--- aspirin 180.16, C₂₀H₃₂O₄ 336.5, tetramethyllead 267 --- and no rule tried
(summed uncertainties, root-sum-square, largest single term) separates the
cases. 7% of weights match *no* rounding of the computed value under any one
table. PubChem itself is not consistent: chloroform is 119.37, which needs
Cl 35.45, and chloroform-*d* is 120.38, which needs Cl 35.453.

So `mw` is computed with the fitted table and rounded half-up to two decimals.
Against the shipped weights: on the 1 401 407 compounds the two
databases share, 62.0% are identical, 74.3% within 0.01, 97.1% within 0.05,
and 337 differ by more than 0.5 --- almost all compounds of lead, which PubChem
rounds to a whole number (tetrabutyllead: 435 against 435.66), and
technetium, which PubChem weighs as ⁹⁷Tc and RDKit's fallback as 98. Rounding
by magnitude instead (three decimals below 50, five significant figures) was
tried and moves the exact-match rate by less than half a point either way, so
the simple rule stayed. Where the two differ, PubChem
has usually rounded harder, and this is the more precise number. The
docstrings, the docs and `provenance["mw_source"]` say it is computed.

Isotopically labelled compounds need the SMILES: PubChem writes chloroform-*d*
as `CHCl3`, so the formula weighs it 119.37. After all columns are in, rows
whose SMILES carries an isotope label are recomputed atom by atom, the label's
isotope mass for labelled atoms and `ATOMIC_WEIGHTS` for the rest
(8 534 compounds).

### 18.4 The scope is a trade, not the 8% gain §8.1 described

§8.1 compared CIDs *carrying a CAS number* and found the FTP mapping ahead by
108 333. Measured against the whole shipped database, 2026-09-01 snapshot:

| | shipped (regex) | FTP (curated) |
|---|---:|---:|
| compounds | 1 589 910 | 1 430 379 |
| compounds with a CAS number | 1 323 167 | 1 430 379 |
| CAS rows failing the check digit | 17 379 | 0 (123 rejected) |
| distinct valid CAS numbers | 1 352 757 | 1 371 578 |

and the difference, both ways:

- **188 503 compounds leave.** 172 408 never had a CAS row: they were in the
  CSV because PubChem's classification browser files them under "CAS", most
  likely through a *substance* that carries one --- which is §8.5's deferred
  substance layer. The other 16 095 had only regex-found CAS numbers PubChem's
  mapping does not give them.
- **28 972 compounds arrive** that the shipped database never had, and
  123 307 more compounds carry a CAS number than did.
- **34 611 valid CAS numbers stop resolving.** They came from synonyms and are
  not in PubChem's mapping; 33 094 of them are in PubChem's own *filtered*
synonym list (names PubChem keeps because they are consistent with the
structure), 25 734 on a compound already in scope. **53 432 resolve that did not.**
- **10 039 CAS numbers resolve to a different first CID** than the shipped
  database's `get_by_cas` returns.

That is fewer compounds and more CAS numbers, every one of them valid, and
every one of them a link PubChem vouches for. For `Search`, which also asks
CompTox and ChEBI for a CAS number, the lost 2.6% is partly covered. It is the
trade decision 2 of §8.2 implied, but it should be a decision rather than a
side effect; §18.9 lists the option.

### 18.5 Descriptors are online-only, whatever is on disk

§8.4 dropped the eight descriptors from `OFFLINE_PROPERTIES`. The question that
raised is what a Zenodo copy --- which still has the columns --- should do. It
could serve them: the columns are there. It does not, for §8.6's reason: before
this change `properties(cid, ["XLogP"])` returned a months-old snapshot for
compounds in the file and live PubChem for the rest, with no way to ask for
either. Now a descriptor has one source, PUG-REST, labelled `online`, until
step 10 adds RDKit's as a named third.

`offline_properties` is per instance because the two kinds of database differ
in the other direction too: `MonoisotopicMass` is offline only in an FTP
build, and making it part of a Zenodo copy's *default* property list would
send every default lookup online for it.

### 18.6 What was left out of §8.4's sketch

- **`scope="cas"`.** One value, and §8.5 keeps the second (`"cas+sid"`) out of
  this version; a parameter that accepts one string is noise. When the
  substance layer lands it can arrive with its parameter.
- **`PubChemID(source="local")`.** `auto_download=False` already means "use what
  is on disk, or raise", and has since before this plan.
- **A release argument on `PubChemID`.** The constructor builds from the newest
  snapshot; a user who wants another calls `build_pubchem_id_db(release=...)`,
  which is one line and says what it does.

### 18.7 The real build

Built once, end to end, against the real 2026-09-01 snapshot
(`TIMESTAMP` 2026/08/31 18:38:04), with `keep_downloads=True` so it could be
repeated without the network.

| file | compressed | lines read | kept | download | read + write |
|---|---:|---:|---:|---:|---:|
| `CID-Identifiers.tsv.gz` | 98 MB | 12 745 309 | 3 238 798 | 2m05s | 22s |
| `CID-Date.gz` | 330 MB | 124 599 998 | 1 430 378 | 6m48s | 67s |
| `CID-Mass.gz` | 1.39 GB | 124 599 998 | 1 430 378 | 32m57s | 83s |
| `CID-SMILES.gz` | 1.49 GB | 124 599 998 | 1 430 378 | 42m16s | 81s |
| `CID-Title.gz` | 1.89 GB | 124 599 721 | 1 430 378 | 54m49s | 87s |
| `CID-IUPAC.gz` | 1.85 GB | 123 699 553 | 1 416 370 | 30m59s | 87s |
| `CID-InChI-Key.gz` | 7.36 GB | 124 599 995 | 1 430 378 | 1h51m30s | 166s |
| `CID-Synonym-filtered.gz` | 0.97 GB | 117 897 265 | 15 770 496 | 17m31s | 97s |
| indexes | | | | | 21s |
| **total** | **15.37 GB** | | | **4h59m** | **11.9 min** |

PubChem served 0.85 MB/s all afternoon --- a fourteenth of the 12.2 MB/s §8.3
measured the day before, and a bare `curl` alongside got the same --- so the
build took **310.8 minutes, 96% of it waiting for bytes**. At §8.3's rate the
download is 21 minutes and the whole build about half an hour, which is what
§8.4 projected; the documentation therefore quotes the processing time and the
transfer size, not a wall-clock time that the network decides. A second build
from the kept files took **11.7 minutes** including re-verifying all 15 GB of
MD5s, and produced an identical `compounds` table. Peak RSS **299 MB**.

The database is **2.48 GB** (2 484 281 344 bytes): 1 430 379 compounds,
1 461 053 CAS rows, 15 770 496 synonyms, 1 777 745 cross-references
(DTXSID 1 129 293, EC 335 218, ChEMBL 145 712, UNII 120 113, ChEBI 47 409).
Synonyms and their text index are half of it (1.27 GB); `compounds` is
0.56 GB, `xrefs` with its indexes 0.12 GB. §8.4's estimate of 1.5–2.5 GB
held. §8.3's "5 M rows" of cross-references was for every CID; restricted to
the scope it is 1.78 M. One CID with a CAS number is in no structure file
(`compounds_without_smiles = 1`); 14 009 have no computed IUPAC name, as in
PubChem.

`datasets.DATASETS["pubchem"]` now carries these numbers: 15 374 591 318 bytes
to download, 2 484 281 344 resident, and a peak of the database plus the
7 361 682 757-byte InChI file. For the four default sources `datasets.plan()`
reports 21.0 GiB to download, 6.5 GiB installed, 37.4 GiB at peak --- the peak
still ChEMBL's.

### 18.8 Verification

- `tests/test_pubchem_ftp.py`, 46 tests. Nothing between the builder and the
  network is stubbed: a miniature release --- the real directory layout, the
  real line formats, gzipped, an `.md5` beside every file, a `TIMESTAMP` --- is
  served by `http.server` on localhost, and the builder lists, resolves,
  downloads, verifies, streams, filters and writes against it. The release is
  built to exercise the decisions: a duplicated CAS row, a Wikidata row that
  must not reach `xrefs`, files in different CID orders, chloroform-*d*, a CID
  whose only CAS fails the check digit and one with no CAS at all, both
  present in every file and required to appear in no table.
- What is pinned: the scope; deduplication; every `compounds` column of
  aspirin exactly; the isotope correction; PubChem's synonym order; the stored
  cross-reference types; the indexes; provenance keys, and that every
  recorded MD5 is the published one; that nothing is left beside the database;
  reuse of kept downloads, and re-fetching of a damaged one; `FileExistsError`
  before any request; that a checksum failure leaves the existing database
  byte-identical and no `.tmp`; that `include_inchi=False` and
  `include_synonyms=False` never request their file; the rolling dump's paths;
  seven PubChem weights computed from their formulas, charge suffixes and bad
  formulas; `PubChemID` building by default, serving `MonoisotopicMass`
  offline and not `XLogP`, `provenance()` and `xrefs()`, the batch tables'
  columns, the Zenodo route, a bad `source` raising before anything is written,
  and an existing database opened unchanged by either route; and that
  `datasets.status` counts and `datasets.remove` deletes kept downloads.
- `tests/test_pubchem_properties.py`: three tests rewritten for descriptors
  going online, one added (`XLogP` goes online although the column exists).
  `tests/test_dataset_manager.py`: PubChem's registry row pinned to the FTP
  route and to `PubChemID.__init__`'s default, as §15.5 did for ChEMBL.
- `examples/pubchem/ftp_build_demo.py` run end to end.

### 18.9 Validation

- `pytest tests/` — **1249 passed, 35 skipped, 19 failed** (8m20s). The
  nineteen are exactly §17.8's: live calls to `pubchem.ncbi.nlm.nih.gov` that
  returned HTTP 503, four in `test_pubchem.py` and fifteen in
  `test_pubchemview.py`, and a bare `curl` of PUG-REST returned 503 at the same
  moment. Re-running with those two files excluded: **1237 passed,
  33 skipped, 0 failed** (6m56s).
- `mkdocs build --strict` — clean, after one docstring example was reworded
  because autorefs read its list output as a cross-reference.
- `pytest --doctest-modules src/provesid/pubchem_ftp.py` — 5 passed, 1 skipped.
- Against the shipped database, on the 1 401 407 CIDs both have: formula and
  InChIKey identical for all, SMILES for all but 55, exact mass within 0.001
  for all, title for 94.8% (PubChem retitles compounds between releases).
  Molecular weight as in §18.3.
- `PubChemID(db_path=<new>)`: `get_by_cas` for formaldehyde, ethanol,
  aspirin, water and benzene returns the same CIDs and names as the shipped
  database; `search_by_name("aspirin", exact=True)` returns 2244; `cas_to_cid`
  runs in 57 µs against 59 µs; `properties(2244)` offline carries
  `MonoisotopicMass`; `xrefs(712)` returns formaldehyde's ChEBI, ChEMBL,
  DTXSID, EC and UNII. Formaldehyde shows the rounding difference: 30.03
  against PubChem's 30.026.

### 18.10 Still open

- **The lost CAS numbers of §18.4.** A `cas_from_synonyms=True` option --- CAS-shaped
filtered synonyms that pass the check digit, added to `cas_numbers` for
compounds already in scope --- would win back 25 734 of the 34 611 at no extra
download, since the synonym file is read anyway. It is a precision/recall
choice PubChem has already made one way, so it is the user's call rather than
a default; if taken, the rows should carry their origin so a lookup can tell a
curated link from a synonym. The 172 408 compounds with
  no compound-level CAS are §8.5's substance layer (`SID-Map.gz`, 3.58 GB).
- **The Zenodo copy is still the CSV build.** It should be replaced by an FTP
  build (`scripts/build_pubchem_id_db.py`), which also gives Zenodo users
  provenance and cross-references; `DEFAULT_DB_URL` then moves.
- **Step 10**, `PubChemID.descriptors()`, is what makes §18.5 whole: until it
  lands, a user wanting XLogP without the network has no route at all.
- **`Search` does not use `xrefs` yet.** DTXSID and ChEBI links from PubChem
  are exactly what `Search` infers by structure matching today; wiring them in
  belongs with step 13's source table.
- **A failed build starts over.** The database is rebuilt from nothing on
  rerun, though downloads resume and kept files are reused. Resuming at the
  file level would need the `.tmp` to record which files it has absorbed.
- **Download and parse are sequential.** Fetching the next file while parsing
  the current one would roughly halve the wall-clock time, at the price of a
  second file on disk at the peak.

---

## 19. Landed on 2026-09-21 — step 10, `PubChemID.descriptors()` (§8.6)

§8.6 asked for the other half of decision 1: the eight descriptors dropped
from the database, computed with RDKit when asked, labelled as RDKit's, with
PubChem's own values still available by name. That is what landed, with the
one naming decision §8.6 left implicit (§19.2), and with the RDKit definitions
chosen by measurement rather than by default (§19.3).

### 19.1 What landed

- `PubChemID.descriptors(cid, descriptors=None, source="rdkit",
  use_online_fallback=True)`, with `descriptors_for_cids` and
  `descriptors_table`, shaped exactly like `properties` /
  `properties_for_cids` / `properties_table`. The row-per-CID reshaping they
  share moved into `PubChemID._table`.
- `source="rdkit"` reads the SMILES through `properties_for_cids(cids,
  ["SMILES"])`, so it gets offline-first-then-online for free. A compound the
  database holds costs no request, and one it does not hold has only its
  SMILES fetched. The result says `Source='rdkit'` either way.
- `source="pubchem"` is `properties_for_cids(cids, names)`, which goes straight
  online for these names and labels the result `Source='online'`, as
  `properties()` already did.
- Module level in `pubchem.py`: `rdkit_descriptors(smiles, descriptors=None)`
  for any structure, and `RDKIT_DESCRIPTORS` / `PUBCHEM_DESCRIPTORS`.
- 36 tests in `tests/test_pubchem_descriptors.py`,
  `examples/pubchem/descriptors_demo.py`, a "Descriptors on Demand" section in
  `docs/api/pubchem.md`, `CHANGELOG.md`.

### 19.2 `MolLogP`, not `XLogP`

§8.6 said the values must be "labelled as RDKit's" and that `Source='rdkit'`
does the labelling. For six of the seven that is enough: a TPSA is a TPSA and a
donor count is a donor count, even when two programs count them differently.
Keeping PubChem's names means a table can switch source without renaming its
columns (tested).

The logP is different. `XLogP` is not the name of a quantity but of a model,
XLogP3, and RDKit's number comes from another model, Crippen's. Returning
Crippen's value under the key `XLogP` would put a wrong claim in the column
name, and `Source='rdkit'` would not undo it once the column is out of the
table. So RDKit's logP is `MolLogP`, RDKit's own name. Asking RDKit for
`XLogP` raises an error pointing to `MolLogP` and to `source="pubchem"`. Asking
RDKit for `Complexity` gets the same kind of redirect, and so does asking PubChem
for `MolLogP`.

### 19.3 The definitions were measured against PubChem

RDKit offers more than one definition of several of these descriptors.
The shipped Zenodo copy still stores PubChem's values, so each candidate was
compared against them on 20 000 random CAS-bearing compounds (4 of the SMILES
did not parse):

| descriptor | RDKit function | equal to PubChem |
|---|---|---:|
| `HeavyAtomCount` | `GetNumHeavyAtoms` | 100.00% |
| `Charge` | `GetFormalCharge` | 99.31% → 100% (§19.4) |
| `HBondDonorCount` | `CalcNumHBD` (= `Lipinski.NumHDonors`) | 93.52% |
| | `CalcNumLipinskiHBD` | 81.57% |
| `RotatableBondCount` | `CalcNumRotatableBonds`, default | 73.90% |
| | `StrictLinkages` | 76.13% |
| | `NonStrict` | 70.84% |
| `TPSA` | `CalcTPSA(includeSandP=True)` | **69.80%** |
| | `CalcTPSA()` default | 60.27% |
| `HBondAcceptorCount` | `CalcNumHBA` (= `Lipinski.NumHAcceptors`) | 63.47% |
| | `CalcNumLipinskiHBA` | 48.41% |
| `MolLogP` vs `XLogP` | `Crippen.MolLogP` | 61.6% within 0.5; median gap 0.37, p90 1.21 |

TPSA is the one where the choice matters: PubChem counts sulfur and phosphorus,
Ertl's original definition and RDKit's default leave them out, and counting
them gains ten points. Rotatable bonds stayed on RDKit's default. The
`StrictLinkages` variant is two points closer, but it is not the usual
definition and not the one a user would expect. The remaining gaps are Cactvs
counting differently from RDKit, not a wrong choice. PubChem counts, for
instance, the bond to a CF3 group as rotatable, and fluorine and halide
counter-ions as acceptors. So the table is in `rdkit_descriptors`' docstring and in
the docs, and "use one source throughout an analysis" is the advice.

`MolLogP` and `TPSA` are rounded to four decimals. Both are sums of tabulated
per-atom contributions with at most four decimals, so this drops the
floating-point residue (aspirin's TPSA came out 63.60000000000001) and nothing
else.

### 19.4 One bug found in the Zenodo copy

Every `Charge` "disagreement" in §19.3 was the Zenodo copy's fault, not
RDKit's. The CSV build stored negative charges as unsigned 32-bit integers:
**8 617 compounds carry `4294967295`-style charges, and not one row has a
negative charge**. Until step 9, `properties(cid, ["Charge"])` served that
column from disk, so every anion in the database got a charge of about four
billion. Step 9 stopped serving descriptors from disk (§18.5), so nothing
returns these values now. It is one more reason to replace the Zenodo copy with
an FTP build (§18.10).

### 19.5 Cost

On the 1.59 M-compound Zenodo copy: **0.47 ms** for one compound, and
**8.0 s for 10 000** random compounds (0.8 ms each), with no requests. Of the
~0.5 ms RDKit spends per molecule, parsing is ~145 µs and Crippen's logP
~245 µs; the other six descriptors together are ~115 µs. Asking only for the
descriptors needed is the one real saving. The wrapper's overhead (validating
names, blocking RDKit's log) is under 5 µs.

### 19.6 Verification

- `tests/test_pubchem_descriptors.py`, 36 tests. RDKit runs for real, and the
  PUG-REST transport is replaced by a recorder that answers any requested
  property from a small table. Pinned: every RDKit descriptor of aspirin
  exactly; order, types, the S-and-P TPSA, net charge; unparsable and empty
  SMILES; the three redirecting errors and an unknown name, all before any
  request; that `source="rdkit"` makes no request for held compounds and one
  `SMILES`-only request for the rest; strictly offline returning None; a
  compound with an unreadable or no structure coming back as a record with no
  values, not as unknown; that a Zenodo-shaped database's stored `xlogp` and
  `hbondacc` are *not* what `source="rdkit"` returns; `source="pubchem"`'s
  values, label and request; its refusal of `use_online_fallback=False`; the
  table's `missing` rows and shared columns across sources.
- `examples/pubchem/descriptors_demo.py` run end to end. Its PubChem section
  took the skip path: PUG-REST answered 503, as in §18.9.

### 19.7 Validation

- `pytest tests/test_pubchem_descriptors.py tests/test_pubchem_properties.py
  tests/test_pubchem_ftp.py`: **110 passed**.
- `pytest tests/`: **1284 passed, 36 skipped, 19 failed** (12m20s). The
  nineteen are exactly §18.9's: live calls to `pubchem.ncbi.nlm.nih.gov` in
  `test_pubchem.py` (4) and `test_pubchemview.py` (15), all HTTP 503. PubChem
  is still refusing this machine a day later.
- `mkdocs build --strict`: clean.
- `pytest --doctest-modules src/provesid/pubchem.py`: the new examples pass. The
  25 failures there are the same 25 as before this change: older examples that
  assume a default database or the network.

### 19.8 Still open

- **`Search` does not offer descriptors.** Its output could carry them as
  optional columns, computed from the consensus SMILES, which would also cover
  compounds PubChem does not hold. That belongs with step 15's presets.
- **`properties()` still accepts `XLogP` and friends** and sends them online,
  which is now the same as `descriptors(source="pubchem")`. Two routes to one
  answer is tolerable while `properties()` is PubChem's property vocabulary. If
  it ever narrows to what is data, those names should redirect to
  `descriptors()`.
- **The Zenodo copy's wrapped charges** (§19.4) are one more reason for the FTP
  rebuild already listed in §18.10.

---

## 20. Landed on 2026-09-21 — step 11, `CheMBL(source="mysql")` (§9.7 route 3)

§9.7 described route 3 as "2.10 GB downloaded, decompressed on the fly, with a
parser that recognises `INSERT INTO` for the eight tables and ignores everything
else", and said it should not ship without being compared with route 2 "table
by table". The parser landed, and so did the comparison. It was run on real
data, but not on EBI's own dump, because `ftp.ebi.ac.uk` refused connections
from this machine all day (§20.6). There is one departure from §9.7: the dump
is downloaded to disk before it is read (§20.2).

### 20.1 What landed

- `src/provesid/mysqldump.py`: `read_statements`, `parse_values`, `unescape`,
  `sqlite_affinity`, and `DumpFormatError`. It has no ChEMBL knowledge: it reads
  `CREATE TABLE` and `INSERT` statements for the tables it is asked about and
  skips the lines of every other table.
- `CheMBL.build_from_mysql_dump(dump_path, dest_path=None, *, keep_inchi,
  remove_source, force)`, with the same signature and contract as `compact`:
  it builds into `<dest>.tmp`, checks the result, then `os.replace`s it into place.
- `CheMBL(source="mysql")`: `SOURCES` is now `("sqlite", "mysql", "full")`,
  and `_UNIMPLEMENTED_SOURCES` is gone. `resolve_latest_db_url(archive=...)`
  finds either archive in the `latest/` listing.
- `CheMBL.extract_digest(db_path)`: a row count and a content hash per table,
  independent of row order and sensitive to each value's SQLite type.
- `_index_and_seal`, split out of `_build_compact`, so the two routes build the
  same indexes and run the same `ANALYZE` and `VACUUM`. The provenance record
  gains `source_format`.
- 54 tests in `tests/test_mysqldump.py`, 33 in `tests/test_chembl_mysql.py`
  (one `slow`), and one in `tests/test_dataset_manager.py`. The third route is
  in `examples/chembl/download_source_demo.py`. Also updated:
  `docs/api/chembl.md` (new section "Installing from the MySQL dump"),
  `examples/chembl/README.md` and `CHANGELOG.md`.

### 20.2 The dump is downloaded, then read; not read off the socket

§9.7 promised a "peak transient disk ≈ 0 beyond the 2.6 GB result". That
would mean parsing the HTTP response directly, which loses the most useful
property of `download_file`: resumption (§13). A dropped connection at 90% of
2.1 GB would then cost the whole transfer again. So the archive is saved with
`download_file`, and `tarfile` then streams it (`r|gz`) without extracting
it. The price is 2.1 GB of transient disk.

The archive is deleted as soon as it has been read, **before** indexing and
`VACUUM` (`remove_source=True`, which is what the constructor passes). `VACUUM`
writes a second copy of the database, so deleting the archive afterwards would
put three large files on disk at once. The peak is therefore about the archive
plus the unindexed extract, or two copies of the extract during `VACUUM`,
whichever is larger: **~4.5 GiB against route 2's 33.4 GiB**. This figure is
reasoned from the file sizes, not measured on EBI's archive.

### 20.3 Types come from the dump, through SQLite's own rule

Route 2's extract gets its column types from `CREATE TABLE … AS SELECT`, which
writes the *affinity* of each source column (`INT`, `TEXT`, `NUM`, `REAL`), and
that affinity decides how each value is stored. A `NUMERIC(2,1)` `max_phase`
of 4.0 is stored as the integer 4. 9 992 `mw_freebase` values are integers,
and the rest are reals. To match that, the mysql route must use the same
affinities. `sqlite_affinity` applies SQLite's documented rule (§3.1 of its
datatype page) to the MySQL type name. `bigint` gives `INT`, `varchar(20)`
gives `TEXT`, `decimal(9,2)` gives `NUM`. The dump's `4.00` is parsed as a
`float` and stored as 4 by the column, which is what happens in route 2. A
test pins `typeof(max_phase)` for an integer, a half and a NULL, and
`extract_digest` hashes `typeof` alongside each value, so an integer-vs-real
mismatch anywhere would fail the comparison.

Projection and filtering happen in SQLite as well. Each table gets an
`INSERT … SELECT <kept> FROM (SELECT ? AS c1, …) WHERE <where>`, built from
the same `COMPACT_TABLES` entry `compact` uses. `molfile` is bound and
dropped, and only `COMPOUND` rows of `chembl_id_lookup` are kept, with the
same `where` text in both routes.

### 20.4 Refusing, not guessing

The parser handles exactly what `mysqldump` emits: quoted strings with
MySQL's escapes, `NULL`, and unquoted numbers, all on one line per `INSERT`.
Anything else inside a *wanted* table raises `DumpFormatError` with the table
name and offset. That includes a hex literal (`--hex-blob`), a `_binary`
introducer, a missing `;`, and an unterminated string. Lines of unwanted
tables are rejected on their first bytes and never decoded, so a blob in
`activities` costs nothing. `build_from_mysql_dump` also refuses a dump with no
`.dmp`/`.sql` member, a missing table or column, or a wanted table with **zero
rows**. An empty table is the likeliest symptom of `INSERT` lines the reader
failed to recognise, so it is treated as an error rather than accepted. In the
constructor, any failure deletes the archive and names `source="sqlite"` as
the route that does not depend on the dump's format.

The string pattern uses possessive quantifiers (`[^'\\]++`, `*+`, available
since Python 3.11). The ordinary `(?:[^'\\]+|\\.)*` backtracks
exponentially on an unterminated string, which is exactly what a truncated line
contains. A test feeds it 150 000 characters with no closing quote and requires
the refusal to arrive within a second.

### 20.5 One adjacent bug: the registry could not see a ChEMBL `.part`

`datasets.DATASETS["chembl"].extras` listed `chembl_*_sqlite.tar.gz`, but
route 2 saves its archive as `db_path + ".tar.gz"`, which is
`chembl_36.db.tar.gz`. An interrupted 5.8 GB download was therefore invisible
to `status()` and survived `remove("chembl")`. The globs now name the files
both routes actually write. A test creates both routes' `.part` files; against
the old registry it freed 0 of their 50 bytes.

### 20.6 Validation

- **The §9.7 comparison, on real data.** EBI's FTP server refused every
  connection (HTTP and HTTPS, `Connection refused`) from this machine for the
  whole session, while other hosts answered. So the 2.1 GB dump could not be
  fetched. As a stand-in, the eight tables of the local ChEMBL 36 release
  (`chembl_36.db`, molfile included) were written out in mysqldump's layout and
  escaping. Rows were in primary-key order, as `mysqldump` writes them, and an
  unwanted table with a hex literal was added at the front. The result was a
  9.60 GB dump, 1.84 GB gzipped. `build_from_mysql_dump` read it in **359 s at
  185 MB peak RSS**, 5.66 MB/s of compressed input, and wrote 2.60 GB. Against
  the `chembl_36_provesid.db` that route 1 built on 2026-09-20,
  **`extract_digest` agrees on all eight tables: 14 550 749 rows with
  identical counts and content hashes.**
- What this does *not* prove: that EBI's dump has the layout this parser
  expects (one statement per line, backticked names, `ENGINE=` closing lines,
  no `--hex-blob`), or how long the 66 unwanted tables take to decompress and
  skip. Those are properties of EBI's file and need EBI's file.
- `pytest tests/`: **1372 passed, 36 skipped, 18 failed** (13m03s). The 18 are
  the live `pubchem.ncbi.nlm.nih.gov` HTTP 503s of §19.7 in `test_pubchem.py`
  and `test_pubchemview.py`, one fewer than yesterday.
- `mkdocs build --strict`: clean. `pytest --doctest-modules
  src/provesid/mysqldump.py`: 4 passed.
- `examples/chembl/download_source_demo.py` ran end to end over a local HTTP
  server with all three routes. The `sqlite` and `mysql` extracts
  digest identically.

### 20.7 Still open

- **Run the comparison on EBI's own dump.** When `ftp.ebi.ac.uk` answers,
  download `releases/chembl_36/chembl_36_mysql.tar.gz`, put it in
  `PROVESID_DATA_DIR` beside `chembl_36_provesid.db`, and run
  `pytest -m slow tests/test_chembl_mysql.py`. Only then should
  `CheMBL.__init__`'s default move from `"sqlite"` to `"mysql"`, together with
  the registry's `download_bytes` and `peak_bytes` and the test that pins them
  to the default (§15.5).
- **No checksum.** EBI publishes `checksums.txt` (SHA-256) per release, while
  `download_file` only verifies MD5. The byte-count check, and the tar stream failing on a short read,
  catch truncation, but not substitution.
- **The peak is reasoned, not measured** (§20.2).


---

## 21. Landed on 2026-09-21 — step 12, the `pubchem.py` and `chebi.py` splits (§4.12)

The August plan's B.3 and B.4, done as specified. `pubchem.py` was 3 730
lines holding an online PUG-REST client and an offline SQLite database. It is
now 2 029 lines of client, beside `pubchem_id.py` at 1 793. `chebi.py` was
1 614 lines and is now 947 lines of REST client, beside `chebi_sdf.py` at 687.
`from provesid import PubChemID, ChebiSDF` is unchanged. Imports that named a
submodule are not, and per dev-principle §1 there is no shim.

### 21.1 What landed

- `src/provesid/pubchem_id.py`: `PubChemID`, plus the descriptor block that
  §19 put in `pubchem.py`, which is `RDKIT_DESCRIPTORS`, `PUBCHEM_DESCRIPTORS`,
  `rdkit_descriptors`, `_rdkit_descriptor_functions` and
  `_check_descriptor_names`. It imports `PubChemAPI` and `PROPERTY_CHUNK_SIZE`
  from `pubchem.py`. Nothing imports the other way, so there is no cycle.
- `src/provesid/chebi_sdf.py`: `ChebiSDF`. `get_chebi_entity` and
  `search_chebi` stay with the online client, as B.4 said.
- A module docstring for each of the four files. `pubchem.py` had only a
  two-line comment before.
- Docstrings for the nine PUG-REST constant classes that B.3 asked for, from
  `Domain` to `OutputFormat`. Four of them carry a doctest.
- Callers updated: `__init__.py`, `search.py`, `datasets._CLIENTS`, eight test
  files, `examples/pubchem/descriptors_demo.py`,
  `examples/ChEBI/chebi_sdf_tutorial.md`, `docs/api/pubchem.md` (a new
  `::: provesid.pubchem_id`) and `docs/api/chebi.md` (a new section on
  `ChebiSDF` with `::: provesid.chebi_sdf`). `CHANGELOG.md` has an entry under
  *Changed*, marked breaking.

### 21.2 Verbatim, and checked to be

B.3 said "verbatim", and a move is only safe if it really is. Every top-level
function, class and assignment in HEAD's two files was compared, by AST span
and source text, with its counterpart in the four new files. The comparison
found no definition missing, none added, and none with a changed body. There
is one exception, made afterwards and on purpose: `ChebiSDF.download_sdf`
gained `-> str` (§21.4). The only other edits to the old files removed the 16
imports that only the moved classes used. One of those was `from rdkit import
Chem` in `chebi.py`, so `import provesid.chebi` no longer loads RDKit.

### 21.3 A test guard that would have gone quiet

`tests/test_dataset_manager.py::no_downloads` makes any bulk download fail the
test. It does this by patching `download_file` in each module that downloads,
using `raising=False`. After the split, `provesid.pubchem` and `provesid.chebi`
no longer have a `download_file`. Leaving the list alone would therefore not
have failed a single test. It would have patched two attributes that nothing
reads and left the real call sites in `pubchem_id` and `chebi_sdf` unguarded.
The list now names the new modules. This is the case where a move breaks
nothing visibly and still removes a safety check, so it is worth recording.

### 21.4 `ChebiSDF` is rendered for the first time

`docs/api/chebi.md` had never had a `:::` directive, so no `ChebiSDF`
docstring had ever been through griffe. Adding one made `mkdocs build
--strict` fail on `download_sdf`, whose docstring promised a `str` without an
annotation. The annotation is correct because the method returns
`self.sdf_path`, and it is the only change to moved code.

### 21.5 Validation

- `pytest` on the affected files before the docstring edits, namely
  `test_pubchem_descriptors`, `test_pubchem_properties`, `test_pubchem_id`,
  `test_pubchem_id_quick`, `test_pubchem_ftp`, `test_chebi_sdf`,
  `test_chebi`, `test_dataset_manager`, `test_search_new_methods` and
  `test_sqlite_lifecycle`: **352 passed, 4 skipped**.
- `pytest tests/`: **1372 passed, 36 skipped, 18 failed** (12m30s). This is
  §20.6's result to the test. The 18 are the same live
  `pubchem.ncbi.nlm.nih.gov` HTTP 503s in `test_pubchem.py` and
  `test_pubchemview.py`. The run began before the constant-class docstrings
  were added, so the offline PubChem, ChEBI SDF and HTTP files were run again
  on the final code: 139 passed.
- `mkdocs build --strict`: clean.
- `pytest --doctest-modules` on the two new modules: 23 failures, and HEAD
  has **the same 23** on the same docstrings in their old location. They are
  step 17's problem. The examples call methods on a database the doctest never
  opened, or ask PubChem, which answered 503 again today (§19.7). The four new
  doctests on the constant classes pass.

### 21.6 Still open

- `docs/api/pubchem.md` and `chebi.md` are hand-written pages with a `:::`
  block dropped into them. Step 18 rebuilds them from the docstrings.
- The historical table in `datasets.py`'s module docstring still lists
  `pubchem.py` and `chebi.py` as download sites. It records what was true when
  §13 found the five copies, so it was left as history rather than rewritten.

---

## 22. Landed on 2026-09-21 — step 13, `sources.py` and the `Search` ladders (§4.5)

The August plan's B.2, done in the shape it proposed. `search.py` was 3 392
lines and is now 2 721. `sources.py` is 378, so the package is 293 lines
shorter. `search.py` had 55 hand-written `if self._x is not None: try: ...
except Exception: log.warning(...)` blocks. Of those, 47 in nine live methods
became one table and one driver, six more in two dead methods were deleted,
and two remain in `_tanimoto_candidates` (§22.3).

### 22.1 What landed

- `src/provesid/sources.py`. `LOOKUPS[kind][source](client, query)` returns
  that source's candidates, best first, already adapted by the `tools`
  `candidate_from_*` helpers, and an empty list on a miss. There are nine
  kinds: `cas`, `inchikey`, `inchikey_skeleton`, `inchi`, `smiles`, `dtxsid`,
  `name`, `fuzzy_name` and `formula`. The module docstring prints the table as
  a grid, so the gaps can be read at a glance. Lookups do not catch
  exceptions. The three skeleton searches that read `_conn` directly and
  `rank_rows_by_completeness` moved here from `search.py`.
- `Query(value, label, k, fuzzy_cutoff)`, a frozen dataclass, instead of the
  bare `(client, query)` lambdas B.2 sketched. B.2's sketch had no way to
  express two things the ladders really did: name and formula lookups take
  `top_k_per_source` rows, and a ZeroPM candidate is named after the user's
  query, which for DTXSID is not the InChIKey ZeroPM was asked by (§22.2).
- `Search._collect(kind, value, *, label, k, sources)` is the only place a
  source is queried. A failing source is logged as `"<Source> <kind> lookup
  failed for <value>: <error>"` and left out, and the others still vote.
  `Search._pool(hits, method, score)` flattens the result in source order and
  then rank order. `score` is either a number or a function of the candidate,
  which covers exact identifiers (1.0), names (`_name_score`) and formulas
  (`_completeness_score`) without three copies of an `add()` closure.
- The five `self._chebi` … `self._chembl` attributes became one
  `self._clients` dict, and `getattr(self, f"_{key}")` is gone. The factory
  map stays local to `_ensure_clients`, and not only out of habit (§22.5).
- The resolvers now read as what they are. The CAS resolver, for example,
  does this:

  ```python
  hits = self._collect("cas", cas)
  smiles = _first_smiles_from_candidates(hits)
  if not is_missing(smiles):
      hits.update(self._collect("smiles", str(smiles), sources=["chembl"]))
  return result, self._pool(hits, "exact_cas"), None
  ```

  `_resolve_cas` went from 66 lines to 24, and most of what remains is its
  docstring. `_skeleton_candidates` is gone, because it was exactly one row
  of the table.
- Deleted as dead, as B.2 asked: `_candidates_from_name`,
  `_fuzzy_name_candidates`, `_most_complete_row` and its three tests, the
  `rapidfuzz.process` import, and the six unused `tools` imports B.2 listed.
  Four adapters that only the ladders used are now imported by `sources.py`
  instead.
- `tests/test_sources.py` has 21 tests: the table's shape, one lookup per
  adapter shape, `Query`, the driver's failure isolation, and pool ordering.
  `docs/api/search.md` gained a "Source lookup table" section rendered from
  the module. `CHANGELOG.md` has an entry under *Changed*.

### 22.2 Cross-source routes stay in the resolver

Not every cell of the grid is filled. ChEMBL has no CAS numbers, ChEBI is
asked about a SMILES by its InChIKey, CompTox is asked about an InChI by its
InChIKey, and the DTXSID resolver asks everyone except CompTox by the
InChIKey that CompTox returned. Those routes depend on what another source
answered first, so they are not table rows. They are one `_collect(...,
sources=[...])` line each in the resolver, which is where the dependency is
visible. Putting them in the table would have needed a small language for
"ask after", and there are only five such routes.

The DTXSID route is the reason `Query.label` exists. HEAD named ZeroPM's
candidate after the DTXSID even though ZeroPM was asked by InChIKey. The
name is meaningless either way, but a refactor that changes output is not a
refactor. So the label is carried through, and a test pins it.

### 22.3 Deliberately left alone

- `_tanimoto_candidates`. It scores each row against a fingerprint, applies
  the threshold per source and returns the best score. That is a similarity
  search with two source-specific probes, not a lookup, and turning it into a
  table row would hide the scoring. It now reads `self._clients` and returns
  the new hits shape, and is otherwise unchanged.
- **§4.13 is still open.** `_clients_initialized` is still one flag, so
  passing one client still disables lazy construction of the others. Fixing
  it needs a way to tell "not passed" apart from "passed `None`", because
  three test helpers pass `chebi=None, ...` precisely to get a one-source
  `Search`. That is a behaviour change with its own design question, so it
  was kept out of a commit whose claim is "no behaviour change".
- §4.4's online fallback is now one more row per kind plus a flag, as §4.5
  predicted. That is step 14.

### 22.4 Equivalence on the real databases

The stub tests show that the code paths run. They cannot show that 47 blocks
were transcribed correctly, and in particular they cannot show which
`exact=`, `limit=` or `or`-fallback each source used. So HEAD was checked
out into a separate worktree, and the same script ran on both trees against
the installed ChEBI, CompTox, PubChemID, ChEMBL 36 and ZeroPM. It covered
7 configurations × 7 identifier types with `n_hits="all"`, so every
candidate that survives clustering is compared and not just the winner. The
configurations were default, `use_zeropm`, `fuzzy`, `fuzzy`+`use_zeropm`,
`inchikey_skeleton`, `similarity_threshold=0.5`, and
`min_source_support=2`+`strip_salts`. The 47 queries included misses, typos,
salts, an unparsable SMILES, InChIKeys with made-up second blocks (for
the skeleton path) and a made-up DTXSID.

**All 437 result rows are identical in every column**, with one exception
that is not a difference. With ZeroPM on, the `Synonyms` string of two SMILES
queries lists the same names in a different order. Running HEAD twice gives
two different orders as well, because `candidate_from_zeropm_smiles` collects
synonyms through a set and Python randomises string hashing per process. That
is an existing nondeterminism in `tools.py`. It is worth fixing, but not
here.

### 22.5 Validation

- Baseline before any edit: the eight `Search`-touching test files,
  **409 passed**.
- `tests/test_sources.py`, `test_search.py` and `test_search_multihit.py`
  after the change: 154 passed.
- `pytest tests/`: **1388 passed, 35 skipped, 21 failed** (85 min). 18 of
  the failures are §20.6's live `pubchem.ncbi.nlm.nih.gov` failures in
  `test_pubchem.py` and `test_pubchemview.py`, the same set §21.5 recorded.
  The other three were this change's, and are fixed. One read
  `getattr(search, f"_{key}")` in `test_sqlite_lifecycle.py`. The other two
  are worth recording. The first draft moved the client factory map onto the
  class as `_CLIENT_FACTORIES`. That captured `ChebiSDF`, `PubChemID` and the
  rest at import time. `test_dataset_manager.py`'s `recording_clients`
  fixture patches `provesid.search.PubChemID` and friends, so the two
  dataset-policy tests silently built the **real** clients against a
  `tmp_path`, instead of recording the calls. That is why the run took 85
  minutes rather than §21.5's 12. The map is back inside `_ensure_clients`,
  where the lookup happens at call time, with a comment saying why. After the
  fix, `test_dataset_manager`, `test_sqlite_lifecycle`, `test_sources`,
  `test_search` and `test_search_multihit` gave 247 passed in 3.7 s.
- `pytest --doctest-modules src/provesid/sources.py`: 2 passed.
- `ruff check --select F,E9` on both modules: clean.
- `mkdocs build --strict`: clean.

### 22.6 Still open

- §4.13, above.
- The ZeroPM synonym order (§22.4) should be sorted or insertion-ordered in
  `tools.candidate_from_zeropm_smiles`.
- `comptox_skeleton_search` and `pubchem_skeleton_search` still reach into
  `_conn`. A public `search_by_inchikey_prefix` on each client would let the
  table call a method, as every other row does.

---

## 23. Landed on 2026-09-22 — step 14, `Search(online_fallback=...)` (§4.4)

§4.4 as proposed, in the shape §22.3 predicted. The fallback is two new
columns in the source table plus a flag. No resolver changed.

### 23.1 What landed

- `Search(..., online_fallback=False)`. It is off by default, so a run still
  opens no socket. When it is on, a query whose resolver returned an empty
  pool is asked again of PubChem PUG-REST and CACTUS. The check sits in
  `_resolve_single`, after the resolver returns, so none of the seven
  resolvers knows the fallback exists. An empty pool is the precise meaning
  of "no candidate from any offline source". A Tanimoto hit, a skeleton
  match or a fuzzy hit is a candidate, so it prevents the fallback.
- `sources.py` gained `ONLINE_SOURCE_KEYS = ["pubchem_online", "cactus"]`,
  their display names (`"PubChem (online)"`, `"CACTUS"`) and one cell per
  kind they can answer. The online question is the query asked as its own
  type, since the identifier types and the lookup kinds share their names.
  Cells by kind:

  | kind | PubChem (online) | CACTUS |
  |---|---|---|
  | `cas`, `name` | `get_cids_by_name(name_type="complete")` | yes |
  | `smiles`, `inchikey` | `get_cids_by_<kind>` | yes |
  | `inchi` | `get_cids_by_inchi` (new) | yes |
  | `dtxsid` | name index, via PubChem's DSSTox synonyms | -- |
  | `formula`, `fuzzy_name`, `inchikey_skeleton` | -- | -- |

  A formula names thousands of PubChem compounds, and neither service has a
  fuzzy or prefix search worth a round trip. So a formula query is never
  sent online and is not counted as a fallback.
- The online sources are asked only through `sources=`. `_collect`'s default
  is still the offline `_SOURCE_KEYS`, so no cross-source route in a resolver
  can reach the network by accident. `_pool`, `_build_source_details` and
  the merge loop in `_build_result_for_cluster` iterate `_SOURCE_KEYS +
  _ONLINE_KEYS`. `_ONLINE_KEYS` is empty when the flag is off, so an offline
  run's `source_details` has exactly the keys it had before.
- The provenance §4.4 asked for goes in `source` and `source_details`, per
  row, and in two new `df.attrs` counters, `online_fallbacks` (queries sent
  online) and `online_resolved` (those answered). Both are reset on each
  `search()` call, and `enrich` carries them. §4.4 also said `foundby`
  should record it. That was not done: `foundby` names the identifier type
  queried, and every other column that describes a row's origin already
  says "PubChem (online)".
- Clients are built by `_ensure_online_clients` on the first query that
  needs them, not in `_ensure_clients`. A run answered entirely offline
  therefore never constructs them. The factories are looked up at call
  time, as §22.5 found necessary, so the tests patch
  `provesid.search.PubChemAPI` and `NCIChemicalIdentifierResolver`.
- New helpers: `PubChemAPI.get_cids_by_inchi`, which POSTs because an InChI's
  slashes cannot travel in a URL path, `tools.candidate_from_pubchem_online`
  and `tools.candidate_from_cactus`. The PubChem candidate is built from one
  property table for all CIDs plus each CID's synonyms, which is where
  PubChem keeps CAS numbers. The CACTUS candidate is built from `smiles`
  plus `names`, with the InChIKey derived by RDKit.

### 23.2 Two votes, not one answer

Each service is an independent vote in `n_source_support`, as a database
is. Two services that agree therefore score like two databases that agree,
and `min_source_support=2` requires both of them. The alternative was to
treat the network as a single "online" source, or to discount it. That
would have written a trust judgement into the scoring that nothing
measured. CACTUS draws on NCI's own collection, not PubChem's, so the two
really are separate witnesses. The per-row `source_details` shows the
difference either way.

### 23.3 "Not found" is a miss, and everything else is a failure

`sources.py` says its lookups do not catch exceptions. The online lookups
make one exception: `NotFoundError` from `http.py` becomes `[]`. Without
it, every offline miss that is also unknown online would log
"lookup failed" at WARNING, which describes a miss as a fault. Any other
error, a 503 or a timeout, still reaches `_collect`, is logged at WARNING and
costs that service's vote. PubChem also answers some unknown identifiers
with CID 0 instead of a 404, and that is filtered out as a miss.

### 23.4 Validation

- The six `Search` test files before the change: 258 passed, and the same
  after it, with `test_sources.py`'s two table-shape tests widened to the
  online keys.
- `tests/test_search_online_fallback.py`: 27 tests, and no socket is opened.
  They cover off-by-default (no client built, `source_details` unchanged),
  offline hits not retried, only the misses going online, all six identifier
  types asked as themselves, name scoring and `top_k_per_source` CIDs,
  formula never sent, the DEBUG lines, not-found versus failure, one failing
  service leaving the other's vote, CID 0, CIDs with no properties, CACTUS
  ambiguity and a missing name list, the POST body of `get_cids_by_inchi`,
  and `enrich` carrying the counters.
- `pytest tests/`: 1427 passed, 34 skipped, 11 failed (9 min 17 s). §20.6's
  live PubChem failures did not recur this time. All 11 failures were
  `test_search_enrichment.py`, which stubs `search()` with a frame that has
  no `attrs`. The first draft of `enrich` copied the provenance from the
  result's `attrs`. It now reads it from the instance, as before, and those
  11 tests and the eight `Search` files give 326 passed. Suite total:
  **1438 passed, 34 skipped, 0 failed**.
- Doctests for `sources.py` and `tools.py`: 13 passed. `ruff check
  --select F,E9`: clean. `mkdocs build --strict`: clean.
- **Live, 2026-09-22.** PubChem answered every request from this machine
  with `503 PUGREST.ServerBusy`, the condition behind §20.6's failures, so
  its column could not be exercised against the real service. The 503 was
  logged and cost PubChem's vote. CACTUS answered, so the fallback delivered
  the answer PubChem could not. On the installed ChEBI, CompTox, PubChemID
  and ChEMBL 36, 18 current CAS numbers chosen to be obscure all resolved
  offline, so the fallback is rarely needed for current numbers. **All 8
  retired CAS numbers tried missed offline** (for example `39400-72-1` and
  `11121-31-6` for atrazine, `11126-35-5` for aspirin), and CACTUS resolved
  all 8 to the correct InChIKey. Retired numbers are what old datasets carry, and they are the
  realistic case for this feature. `examples/search/online_fallback_demo.py`
  uses one.

### 23.5 Found while validating: `CASRN` is the smallest CAS string, not the best one

With no CAS in the query, the live runs reported atrazine as `11121-31-6`
and aspirin as `11126-35-5`. Both are real but retired numbers. CACTUS
lists the current number first. `extract_cas_values` sorts,
`make_candidate` sorts again, and `first_cas` takes element 0, so the
lexicographically smallest number wins. Nothing about that is
online-specific. CompTox's own `CASRN` column and PubChem's synonym order
are discarded the same way offline. Fixing it means keeping each source's
order in `CAS_candidates`, which changes offline output. That is its own
commit with its own equivalence check, so it was kept out of this one.

### 23.6 Still open

- §23.5, CAS order.
- PubChem's column has not been run against the live service (§23.4).
  DTXSID-by-synonym rests on PubChem's documented DSSTox deposit and has not
  been observed from here.
- A batch sent online while a service is down logs one WARNING per query.
  Step 16's circuit breaker is the fix: once a host is known to be
  throttled, `RateLimiter` should fail its requests at once, with one
  message.
- §4.13 and §22.6 are unchanged.

---

## 24. Landed on 2026-09-22 — step 15, `Search.PRESETS` (§4.6)

§4.6 as proposed. `balanced`, `strict` and `recall` are dicts of constructor
arguments. Explicit arguments override them, and each result frame records
which one it ran under.

### 24.1 What landed

- `Search.PRESETS`, a class attribute with three entries over the same
  thirteen keys: `fuzzy`, `fuzzy_score_cutoff`, `fuzzy_scorer`,
  `inchikey_skeleton`, `similarity_threshold`, `use_zeropm`,
  `top_k_per_source`, `cluster_by_skeleton`, `consensus_compat_threshold`,
  `query_weight`, `n_hits`, `min_confidence`, `min_source_support`. These are
  the arguments that decide what counts as a match and what is returned.

  | preset | differs from `balanced` in |
  |---|---|
  | `balanced` | nothing; it holds the constructor defaults |
  | `strict` | `min_source_support=2` |
  | `recall` | `fuzzy`, `inchikey_skeleton`, `similarity_threshold=0.7`, `use_zeropm`, `n_hits="all"` |

- `Search(..., preset="balanced")`. The thirteen arguments now default to
  `None`, meaning "take the preset's value", which is the convention
  `search()` already used for its per-call overrides. With the old literal
  defaults, `preset="recall", fuzzy=False` would have been
  indistinguishable from not passing `fuzzy`. The thirteen defaults are now
  written only in `PRESETS["balanced"]`, and each docstring entry says
  "Balanced: X".
- `Search.settings`, a property giving the values in force, keyed as
  `PRESETS`, and `Search.preset`.
- `df.attrs["preset"]` and `df.attrs["settings"]` on every `search()` and
  `enrich()` frame. `settings` includes that call's `n_hits`,
  `min_confidence` and `min_source_support`, so it describes the call and
  not just the instance. The four existing provenance entries now come from
  one `_provenance()` helper. `enrich` still reads them from the instance,
  not from `results.attrs` (§23.4).
- Left out of the presets on purpose: `online_fallback` (a preset should not
  open a socket), `use_opsin` (needs Java), `strip_salts` and
  `return_alternatives` (they change output columns, not what matches),
  `datasets` and the client arguments.
- `examples/search/presets_demo.py`; a "Presets" section and the two new
  members in `docs/api/search.md`.

### 24.2 Measured on the installed databases

200 queries per sample, drawn by `rowid` stride from CompTox (whose
structure is the truth), and scored by InChIKey skeleton as in
`test_search_precision_regression.py`. `n_hits=1` for every preset,
recall included, so each gets one answer to be judged. Synonyms come from
the test's `_pick_synonym`. Typos are the preferred name with one interior
character deleted (seed 0).

| sample | preset | correct | wrong | not found | time |
|---|---|---:|---:|---:|---:|
| CAS | balanced | 200 | 0 | 0 | 0.5 s |
| | strict | 199 | 0 | 1 | 0.5 s |
| | recall | 200 | 0 | 0 | 7.7 s |
| synonym | balanced | 28 | 1 | 171 | 36 s |
| | strict | 2 | 0 | 198 | 35 s |
| | recall | 197 | 3 | 0 | 628 s |
| typo | balanced | 0 | 3 | 197 | 33 s |
| | strict | 0 | 1 | 199 | 33 s |
| | recall | 39 | 60 | 101 | 730 s |

What this says about the presets:

- On CAS numbers all three are right, and strict costs one answer in 200:
  a structure only one database holds.
- `strict` does what it says. Its answers are never wrong in these samples,
  and on names it answers almost nothing. Most of these synonyms are known
  to one database only.
- `recall`'s top hit on a typo is wrong 60 times in 99. That is why
  `recall` returns `n_hits="all"` and its docstring says to read
  `confidence` and `n_source_support`: it finds candidates for review, not
  answers. It is also 20× slower than balanced on names.

### 24.3 Found while measuring: CompTox's exact name lookup ignores its synonyms

Ablating the synonym row shows the gain is fuzzy widening's alone:
`balanced + use_zeropm` scores 31 correct, while `balanced + fuzzy` scores
197 correct, 2 wrong and 1 not found in 607 s. The queries are not typos.
They are names like `Acetaldoxime` and `2,3-dimethylvaleraldehyde` that
CompTox holds verbatim. `CompToxID.search_by_name(exact=True)` and
`get_by_name` compare against `PREFERRED_NAME` only
(`comptox.py:429`, `:372`). The synonyms sit in the pipe-separated
`IDENTIFIER` column, and only the `LIKE '%…%'` path reads it. That path is
a full-table scan at about 3 s a query. So balanced and strict miss about
85 % of CompTox synonyms, and fuzzy finds them slowly by substring.

The fix is a synonym table built with the CompTox database, one row per
`(name, DTXSID)` and indexed, as §12 did for ChEMBL. The exact `name` cell
in `sources.py` would then consult it. That changes offline output and
needs its own equivalence check, so it is not part of this step. The
sample also flatters the fix, because every synonym was drawn from
CompTox.

### 24.4 Validation

- `tests/test_search_presets.py`: 18 tests, no database opened. They cover
  the table's shape; every key being a constructor argument that defaults
  to `None`; balanced equalling a bare `Search()`; strict differing from
  balanced in one key; each preset reproduced by `settings`; explicit
  overrides, including one equal to the balanced value; unknown presets and
  scorers refused; `PRESETS` not aliased by an instance; strict dropping an
  uncorroborated hit that balanced keeps; and `attrs` recording per-call
  settings in `search` and `enrich`.
- The eight existing `Search` files plus `test_sources.py`: 394 passed,
  unchanged by the new `None` defaults.
- `pytest tests/`: **1456 passed, 34 skipped, 0 failed** (7 min 17 s), which
  is §23.4's 1438 plus the 18 new tests.
- The `settings` doctest passes. `ruff check --select F,E9` is clean, and
  `mkdocs build --strict` is clean with `PRESETS` and `settings` rendered.
- `examples/search/presets_demo.py` ran end to end on the installed data.
  Recall resolves the typo `atrazin` to atrazine (4 sources, 0.77); balanced
  and strict leave it unresolved.

### 24.5 Still open

- ~~§24.3, the CompTox synonym index.~~ **Done, §25.**
- §19.8's descriptor columns were marked "with step 15's presets". They
  are not a preset matter: they add output columns and change no matching
  setting. They belong with `strip_salts` and `return_alternatives`, as an
  output option.
- Presets name policies, but the three were chosen by reasoning, and §24.2
  measured them afterwards. No sample has yet tested whether, for example,
  `cluster_by_skeleton=False` belongs in `strict`.
- §23.5, §4.13 and §22.6 are unchanged.

---

## 25. Landed on 2026-09-22 — the CompTox name index (§24.3)

This was done ahead of step 16, because §24.3 showed it was the largest
recall gap in the default preset. It uses the approach §12 took for ChEMBL:
a lookup that had to scan gets an index to search instead.

### 25.1 What landed

- `CompToxID.build_name_index()` adds a table, `chemical_names`, to
  `comptox_chemicals.db`. It is keyed `(name_key, kind, chemical_rowid)`,
  `WITHOUT ROWID`, with one row per distinct name of each chemical. The
  names are its `PREFERRED_NAME`, its `IUPAC_NAME` and every `|`-separated
  token of `IDENTIFIER`. `kind` records where each name came from (0
  preferred, 1 IUPAC, 2 identifier; `NAME_KINDS`). `name_key` is
  `str.strip().lower()`, applied in Python on both the build and the lookup
  side, so it folds non-ASCII case where SQLite's `lower()` would not.
- `search_by_name(name, exact=True)` reads that table. It is now
  case-insensitive, matches synonyms, IUPAC names and retired CAS numbers,
  and orders results by `kind` and then by database row. A chemical *called*
  the query therefore comes before one that only lists it. `exact=False` is
  unchanged: a substring scan of about 3 s. `get_by_name` still means the
  preferred name only.
- **When it is built.** `download_database` builds it straight after the
  download. On a database downloaded before this change, the first exact
  lookup builds it and logs one WARNING (~20 s). `has_name_index` reports
  whether it exists. On a read-only file the build fails once, is logged
  once, and exact lookups fall back to the old case-sensitive
  `PREFERRED_NAME = ?`. `BEGIN IMMEDIATE` takes the write lock before the
  existence check, so two processes cannot both build it. A lock on the
  instance does the same for two threads.
- `sources.py`'s CompTox `name` cell no longer has its
  `or get_by_name(...)` fallback. That ran the same `PREFERRED_NAME = ?`
  query as the old exact path, so it could never add anything.
- The `comptox` registry entry reports 1.1 GiB resident, with a note, and
  `Search`'s docstrings say ~6.7 GiB installed instead of ~6.5.
- `tests/test_comptox_name_index.py` (22 tests, on a three-chemical
  database), `examples/comptox/name_index_demo.py`, and a section in
  `docs/api/sqlite_clients.md`.

### 25.2 Layout, measured

On the installed release (1 246 399 chemicals), the index has 5 128 983
names:

| layout | added | build |
|---|---:|---:|
| ordinary table plus index on `name_key` | 540 MiB | 9 s |
| **`WITHOUT ROWID`, key-first primary key** | **290 MiB** | **22 s** |
| same, with a 64-bit hash for the key | 100 MiB | 24 s |

The hash layout saves 190 MiB, about 3 % of a full install. In exchange,
every hit would need a second check against the row's names, and the
table could not be read with plain SQL. Per the readability principle, the
text key was kept. An exact lookup takes 30–110 µs.

### 25.3 Before and after, on the four default databases

The samples are §24.2's, plus 200 preferred names. "Before" is this code
with the index forced off, which is exactly the old lookup. "After" is an
indexed copy of the same file.

| sample | preset | before (right / wrong / none) | after | answers changed |
|---|---|---|---|---|
| CAS | balanced | 200 / 0 / 0 | 200 / 0 / 0 | 0 |
| CAS | strict | 199 / 0 / 1 | 199 / 0 / 1 | 0 |
| preferred name | balanced | 200 / 0 / 0 | 200 / 0 / 0 | 0 |
| preferred name | strict | 189 / 0 / 11 | 189 / 0 / 11 | 0 |
| synonym | balanced | 28 / 1 / 171 | **200 / 0 / 0** | 172, all to right |
| synonym | strict | 2 / 0 / 198 | 28 / 0 / 172 | 26, all to right |
| typo | balanced | 0 / 3 / 197 | 0 / 3 / 197 | 0 |
| typo | strict | 0 / 1 / 199 | 0 / 1 / 199 | 0 |

No answer went from right to wrong, and nothing outside synonym queries
moved. The one synonym that had been answered wrongly is now answered
correctly: exact CompTox candidates now outrank whatever the other sources
had matched loosely. Timing is unchanged, because the other three sources
dominate a name query.

The sample is still CompTox's own synonyms, so 200/200 flatters the
change. A synonym only another source holds is unaffected. What the table
shows reliably is that nothing regressed.

### 25.4 Validation

- `tests/test_comptox_name_index.py`: 22 passed. They cover the key, the
  row count and kinds, idempotence, the lazy build and its WARNING,
  substring lookups not triggering a build, four threads building once, the
  read-only fallback warning once, synonym, case, IUPAC and retired-CAS
  matches, ranking, `limit`, parsed identifiers, `get_by_name`'s meaning,
  and the `sources.py` cell.
- `tests/test_dataset_manager.py`: 43 passed with the new registry size.
- `examples/comptox/name_index_demo.py` was run against the indexed copy,
  not the default path, so that running the demo would not build the index
  into the user's database as a side effect of validation.
- `pytest tests/`: 1476 passed, 34 skipped, 2 failed (8 min 48 s). The run
  was also the lazy path's first real use. `conftest.py` points the suite
  at the git-ignored copy in `src/provesid/data/`, and the first exact name
  lookup built the index there (856 MB → 1.16 GB).
  - Both failures were `test_search_precision_regression.py` treating
    `asprin` as a typo that no database lists. CompTox lists it among
    aspirin's synonyms, so an exact lookup now finds aspirin, correctly and
    correctly labelled `exact_name`. The two tests now use `aspirn`, which no
    source lists and which fuzzy matching still resolves to aspirin. A new
    test, `test_a_listed_misspelling_is_an_exact_match`, pins the `asprin`
    behaviour. The file: 11 passed.
  - Suite total: **1479 passed, 34 skipped, 0 failed**.
- `ruff check --select F,E9` is clean on the new and changed files, apart
  from an unused `Union` import in `comptox.py` that was already there. The
  new `name_key` doctest passes. The module docstring's example, which reads
  `result['preferred_name']`, was already failing and is left for step 17.
  `mkdocs build --strict` is clean.

### 25.5 Still open

- ~~**Retired CAS numbers, offline.**~~ **Done, §26.** The index already maps `39400-72-1` to
  atrazine. §23.4 found that all 8 retired numbers it tried missed offline,
  and that only CACTUS resolved them. Consulting the index from the `cas`
  cell when `get_by_casrn` misses would answer them without the network.
  That changes CAS output, so it is its own step, with a check that no
  current number moves.
- **A synonym shared by several substances.** A generic synonym like `ASA`
  now reaches every chemical that lists it, ranked by `kind` and row order.
  None of the 200 was ambiguous in a way that hurt, but a sample of
  deliberately generic synonyms has not been run.
- The installed `~/.local/share/provesid/comptox_chemicals.db` on this
  machine does not have the index yet. It will be built on the first exact
  name lookup, or when `CompToxID().build_name_index()` is called. The
  suite's copy in `src/provesid/data/` has it.
- Found while choosing a replacement typo: `Search("name", fuzzy=True,
  use_zeropm=True)` resolves `asprine` to **Fenyramidol**, which is
  PHENYRAMIDOL, the "Evasprin" compound of the regression this test file
  was written for. This happened before §25 as well. The fuzzy-name path
  still lets a substring-sharing synonym outrank the intended compound.

---

## 26. Landed on 2026-09-22 — retired CAS numbers, offline (§25.5)

§25's name index already held every CAS number CompTox lists, current or
not. This change lets a CAS query use it. It is the offline counterpart of
the job §23.4 found only CACTUS could do.

### 26.1 What landed

- `CompToxID.get_by_alternate_casrn(casrn)` returns the chemical that lists
  `casrn` among its `IDENTIFIER` tokens (name-index `kind` 2), but only
  when **exactly one** chemical does. It returns `None` when the input is
  not shaped like a CAS number, when no chemical lists it, when two do, or
  when there is no index. `get_by_casrn` is unchanged.
- `sources.py`'s CompTox `cas` cell is now `get_by_casrn(v) or
  get_by_alternate_casrn(v)`. A number that is some chemical's own `CASRN`
  therefore belongs to that chemical, whatever else lists it. The new path
  runs only on a miss, so no current number can change answer.
- `examples/search/online_fallback_demo.py` had used `39400-72-1` as its
  example of what only CACTUS knows. That number is now answered offline,
  and the comments say so. `examples/comptox/name_index_demo.py` and
  `docs/api/sqlite_clients.md` show the new method.

### 26.2 How ambiguous the listed numbers are

On the installed release, 1 303 984 CAS-shaped identifier tokens are
listed, and 1 220 051 of them are also some chemical's `CASRN`. That leaves
**83 933** that appear only in `IDENTIFIER`, and so only the new path can
answer them. 83 929 are listed by a single chemical. Four are listed by two
unrelated ones, for example `5990-67-0` by both tetrandrine and a
piperidinyl pyridinecarboxylate. They look like data errors, and they go
unanswered rather than guessed. One of the 83 933 fails the CAS checksum.
It is still answered, because a typo'd number CompTox recorded is a
number some dataset may carry.

A first draft took any single match. The test suite caught the flaw: an
identifier token is not necessarily a CAS number. "Aspirin" is listed by
exactly one chemical as a *synonym*, so the method returned that chemical as
the owner of the "CAS number" `Aspirin`. The input is now required to be
CAS-shaped.

### 26.3 Before and after

300 numbers drawn with seed 0 from the 83 929 unambiguous ones, which is
the only population whose answers can change, plus the three retired
numbers §23.4 named, whose structures CACTUS had confirmed. "Before" is
this code with `get_by_alternate_casrn` returning `None`, which is exactly
the old lookup. Truth is the structure of the chemical CompTox lists the
number under. That is circular for CompTox's own answers, so the
informative rows are agreement with the other sources and the three known
numbers.

| preset | before (right / wrong / no structure / none) | after |
|---|---|---|
| balanced | 6 / 1 / 1 / 292 | 244 / 0 / 53 / 3 |
| strict | 1 / 0 / 1 / 298 | 12 / 0 / 1 / 287 |

"No structure" covers two cases. 45 are substances CompTox gives a
DTXSID but no InChIKey, such as mixtures and UVCBs. 8 are structures from
another source for chemicals CompTox holds no InChIKey for.

- **Agreement.** Before the change, another source answered 8 of the 300.
  Afterwards CompTox agrees with it on 7, and 6 of those gain a
  corroborating vote. That vote is what moves 11 answers into `strict`.
- **The one disagreement** is `68459-97-2`. PubChem draws it as a 1:2:1
  zinc chloride / diazo-amine compound. CompTox lists it as an alternate
  number of the bis-diazonium tetrachlorozincate (`15280-31-6`). That is
  the same zinc chloride diazonium salt with a different stoichiometry.
  Both answers have one source, so the tie now goes to CompTox. Neither is
  clearly wrong.
- **The three known numbers** (`39400-72-1` and `11121-31-6` for atrazine,
  `11126-35-5` for aspirin) went from none to right in balanced. They stay
  none in strict, correctly, because CompTox is their only witness.
- Timing: 300 queries took 0.4 s before and 0.6 s after.

### 26.4 Validation

- `tests/test_comptox_name_index.py`: 29 passed (22 + 7). The new tests
  cover a retired number found, an ambiguous number unanswered, an unknown
  number, a listed synonym not treated as a CAS number, the `cas` cell
  preferring a chemical's own `CASRN` over another chemical's listing, the
  cell's fallback, and no index meaning no answer.
- `_CompToxStub` in `test_search.py` gained `get_by_alternate_casrn`.
  Without it, every CAS miss against the stub would have been logged as a
  failed lookup.
- `pytest tests/`: **1486 passed, 34 skipped, 0 failed** (7 min 53 s).
- Both demos were run. Search on the four default sources with
  `online_fallback=True` answers `39400-72-1` from CompTox and sends only
  the nonsense number online.

### 26.5 Still open

- **`CASRN` reports the number queried, not the current one.**
  `_resolve_cas` sets `result["CASRN"] = cas`, so atrazine found by
  `39400-72-1` is reported as `CASRN = 39400-72-1`. The current `1912-24-9`
  is available in the candidate. The online fallback does the same, so
  this is §23.5's question, and the two should be settled together: should
  a row say what was asked, or what the compound is? A `CASRN_current`
  column, or a note in `source_details`, would give both.
- **A CAS example for the online fallback.** Of 3 000 ZeroPM CAS numbers,
  only 4 miss all four default databases, and CACTUS answered none of
  them. Only one of the four was even sent online. The other three produced
  an offline candidate with neither a structure nor a DTXSID, and under
  §23's rule that blocks the fallback. Whether a candidate with no
  identifiers at all should count as "found" deserves a look.
- §24.5's other items, §23.5, §4.13 and §22.6 are unchanged.

---

## 27. Landed on 2026-09-22 — step 16, the circuit breaker (§4.12)

§23.8 of the August plan and §23.6 here asked for the same thing. A
`Retry-After` is information about the *host*, so it belongs on the host's
shared `RateLimiter` as a "not before T" that every client respects. A host
known to be throttled should then fail at once, with one message.

### 27.1 What landed

- `RateLimiter` gained `not_before`, `hold(seconds)`, `held_for()` and
  `release()`, and a `host` for messages. `host_limiter` names its limiters.
  `provesid.release_holds()` clears every host's hold in the process.
- `HTTPClient.request` records every `Retry-After` it is sent on the
  limiter, before it decides whether to retry. A hold only moves later.
  A shorter `Retry-After` never shortens a longer one already in place.
- **On entry**, a call checks the hold. If the hold is no longer than the
  client would wait anyway, which is `max_backoff` capped by `max_elapsed`,
  the call waits it out, and the wait counts against `max_elapsed`. A longer
  hold raises the client's `rate_limit_cls` **without a request**, with a
  message naming the host and the time. The exception carries
  `held_until`.
- **Within a call**, a `Retry-After` longer than `max_backoff` now ends the
  retries. Before, it was cut down to `max_backoff`, and the client asked
  again before the time the host had named.
- `Search._collect` logs a lookup refused by the breaker at DEBUG. The
  transport already warned once, when the hold was recorded.
- `tests/conftest.py` releases all holds after every test. Holds live on
  process-wide limiters, and a stub answering `Retry-After: 30` under
  PubChem's real host would otherwise fail every later PubChem test.
- `docs/api/http.md` gained "A throttle belongs to the host too".
  `examples/http/circuit_breaker_demo.py` shows the breaker against a stubbed,
  throttled PubChem, offline.

### 27.2 Checked only on entry, not before each retry

The retry loop already waits a `Retry-After` it can afford, so checking
the hold again before each retry would only re-read a time the loop has
just slept through. It would also make the retry depend on the wall clock.
The suite replaces `time.sleep` and not `time.time`, so a stubbed sleep would
leave the hold in place, and every retry test would refuse its own second
attempt. What the breaker has to stop is the *next* call, and that call
checks on entry.

### 27.3 The one behaviour change a caller can see

`test_retry_after_is_capped_by_max_backoff` asserted that a service asking
for an hour got asked again after `max_backoff` (5 s), and that the second
answer was used. Under the breaker, that is the request the host asked not
to receive. The test is now
`test_a_retry_after_longer_than_max_backoff_gives_up_instead_of_asking_early`:
one request, no sleep, `rate_limit_cls`. For the clients as configured this
matters little. `max_backoff` is 60 s for all of them, and both PubChem
clients already declined `Retry-After: 30` through their 10 s `max_elapsed`.

### 27.4 Effect

Against a stubbed PubChem answering `503` with `Retry-After: 30`
(`examples/http/circuit_breaker_demo.py`), 1 000 calls split over `PubChemAPI`
and `PubChemView` after the first refusal send **no requests** and take
18 ms in total. Before, each of them sent one request. At PubChem's 0.2 s
pacing that is about 200 s of refused requests, each one telling PubChem
that the block is still needed. That figure is arithmetic, not a
measurement. The breaker has not been seen tripping against the live
PubChem.

### 27.5 Validation

- `tests/test_http.py`: 70 passed (60 + 10). The new tests cover: a held host failing
  with no request; the hold binding a second client on the same host; one
  host's hold not touching another; a short hold waited out; the hold's wait
  counting against `max_elapsed`; a longer hold never shortened; no hold
  without `Retry-After`; `release_holds`; `held_until` set only by the
  breaker; `held_until` on a plain exception class. The capped-`Retry-After`
  test was rewritten, as §27.3 describes.
- `tests/test_pubchem_failure_reporting.py`: a throttle seen by PUG-REST
  stops a PUG-View call without a request.
- `tests/test_search_online_fallback.py`: a held host is logged at DEBUG,
  not WARNING.
- Doctests in `http.py`: 20 passed, 5 skipped (the network examples).
- `pytest tests/`: **1498 passed, 34 skipped, 0 failed** (7 min 23 s).
- The demo was run: 1 000 refused calls, one request sent, 18 ms.

### 27.6 Still open

- **The hold is per process**, like the pacing clock. Two processes on one
  machine each learn a throttle for themselves.
- **Bulk downloads do not see it.** `datasets.download_file` has its own
  retry and reads `Retry-After` itself. Its hosts (the FTP and Zenodo
  mirrors) are not the API hosts, so there is nothing to share yet.
- **Only `Retry-After` trips the breaker.** A service that fails every
  request with a bare 503 still costs one retry curve per call. Counting
  consecutive failures per host would catch that. It would also have to
  decide when a host is healthy again, which `Retry-After` decides for us.
  Nothing here needs it yet.
- §26.5, §24.5's other items, §23.5, §4.13 and §22.6 are unchanged.

## 28. Landed on 2026-09-22 — step 17, docstrings with examples (§4.12)

August's workstream C asked for every public object to carry a Google-style
docstring with a runnable example, checked as a doctest, with offline
examples using real values and network examples showing real but skipped
output. That is done for all of `src/provesid/`.

### 28.1 What landed

- **Coverage.** All 521 public objects have a docstring with an example. The
  same AST pass at the start of the step found 13 with no docstring and 305
  with no example (it then counted property setters too). Every module has a
  module docstring;
  `__init__`, `opsin`, `cascommonchem`, `classyfire`, `pubchemview`,
  `resolver`, `utils` and `zeropm` had none.
- **A way to run them.** `pytest --doctest-modules src/provesid`, set up by
  `src/conftest.py`, which is not installed. It skips a module's examples
  when the database they read is absent, so a doctest run never downloads,
  it sandboxes the cache and the config directory, and it refuses every
  outbound connection so that an example needing the network has to say so
  (§28.6). `doctest_optionflags` gains `NORMALIZE_WHITESPACE` and `ELLIPSIS`.
  `TESTING.md` says how to run it, and how to lift the block
  (`PROVESID_DOCTEST_ALLOW_NETWORK=1`) to re-record an online output.
- **Result.** 477 passed, 87 skipped, 0 failed, with the network refused,
  where the baseline was 133 passed, 59 failed, 35 skipped. With
  `PROVESID_DATA_DIR` pointing at an empty directory: 281 passed, 283
  skipped, nothing downloaded --- the directory gains one empty `chebifier/`
  from a classifier that is constructed but never run. A skipped item is a
  docstring every one of whose examples is `+SKIP`, not a docstring without
  examples.
- **Outputs are real.** Offline outputs come from the installed databases.
  Online ones were recorded against the live services on 2026-09-22 and
  marked `+SKIP`. CAS Common Chemistry is the exception, see §28.4.
- **`Examples:`, not `Example:`.** griffe renders a Google `Examples:`
  section as a code block but an `Example:` section as a Markdown
  admonition, where `rows[0]["cas"]` is read as a reference link. August's
  house style said `Example:`; every header is now `Examples:`, and
  `mkdocs build --strict` is clean. Trailing whitespace was stripped from
  every module touched, in the same commits (`git diff -w` shows the rest).
- **Validation.** `pytest tests/`: 1511 passed, 33 skipped, 0 failed (9 min).
  `pytest --doctest-modules src/provesid`: 476 passed, 88 skipped.
  `mkdocs build --strict`: clean.

### 28.2 Defects found by writing the examples

Each is fixed with a test that fails on the old code.

| Module | Defect |
|---|---|
| `pubchem_id` | The Zenodo copy repeats (cid, cas) pairs: `cid_to_cas(2244)` was `['50-78-2', '50-78-2']`. `search_by_name` could return a compound twice. |
| `zeropm` | `query_name_regex(case_sensitive=True)` used `LIKE`, which ignores case; now `GLOB`. `get_id_table_from_zeropm_id` repeated rows (no `DISTINCT`). |
| `cache` | `metadata.json` is flushed every 100 writes and never at exit, so `export_cache` returned `True` having exported nothing from a short-lived process, and `disk_entries` read 0 beside 37 files. Both now read the files. |
| `chebi` | All seven text-body calls (`calculate_*`, `depict_structure`) got HTTP 406 because the session sends `Accept: application/json`; they returned None. The mocks never checked headers. They also take a molfile, not SMILES. |
| `resolver` | `quote(safe='')` sent `/` as `%2F`, which CACTUS's front end answers with 404: every InChI and every stereo SMILES was "not found", including through `Search`'s online fallback. The live test skipped that failure as "not supported for this format". |
| `tools` | `candidate_from_zeropm_smiles("CCO")` reported 13C-labelled ethanol's InChIKey: it took the first row of whichever CAS sorted first. |
| `opsin` | `PYOPSIN.get_id_from_list` gave the statuses `'S'`, `'U'`, `'C'`: it indexed one status string per name. |
| `cascommonchem` | CAS answers a bad key with 403; only 401 mapped to "Unauthorized - Check API Key". |
| `sources` | `comptox_skeleton_search` and `pubchem_skeleton_search` read `client._conn`, gone since `SQLiteClient`; the AttributeError was logged, so `inchikey_skeleton=True` drew on ChEBI alone. |

### 28.3 Docstrings that were wrong

Among others: `PubChemID.properties(use_online_fallback=False)` returns None,
not a partial record, when the request needs the network; nine ChEMBL
examples used molregno 15 for aspirin (it is 1280; 15 is CHEMBL6214);
`extract_identifiers_from_synonyms`' `ec_number` holds Enzyme Commission
numbers, not EC inventory numbers; `get_compound_properties` reports a bad
property in its dict rather than raising; `get_cas_by_substructure` searches
only the first 10 000 ZeroPM structures; `Search`'s module example claimed
`fuzzy=True` rescues "caffiene", which it does only with ZeroPM in the pool
(`preset="recall"`); `mw_within` and `resolve_cascade` printed counts from
no dataset. ClassyFire's docstring now says it has classified nothing new
since February 2023 (step 7 is still postponed).

### 28.4 Side effects of this step on the working machine

- **The CAS key.** The baseline doctest run executed the old
  `set_cas_api_key("your-cas-api-key-here")` example against the real
  `~/.config/provesid/config.json` and replaced the stored key. The earlier
  key cannot be recovered from the file. `src/conftest.py` now sandboxes the
  config directory. The stored placeholder outranks `CAS_API_KEY` in the
  environment, so it has to be replaced or removed before CAS works again.
  For the same reason no live CAS output was recorded.
- **The CompTox name index** was built in the installed database (~290 MiB)
  by an exact name lookup made while collecting outputs, as the library does
  on first use.

### 28.5 Still open

- §23.5 shows in the examples: `first_cas` picks the retired `11126-35-5`
  for aspirin from CompTox, and a fuzzy "asprin" returns it as `CASRN`.
- `ChEBI.get_complete_entity` and `batch_get_entities` are compatibility
  aliases, which dev-principle 1 says not to keep.
- `NCIChemicalIdentifierResolver.get_molecular_data`, `resolve_multiple` and
  `batch_resolve` are cached even when some representations failed.
- `ClassyFireAPI.query_status` is cached, so polling with the cache on sees
  the first answer forever (documented, not changed).
- The online examples are skipped by design, so they will drift as the
  services change; re-record them before a release
  (`PROVESID_DOCTEST_ALLOW_NETWORK=1`, and read the diff: a service that has
  changed its answer is the point of the exercise).

### 28.6 Three examples were asking live services

Re-running the step's own verification with outbound connections refused
turned up three examples that reached the network without
`# doctest: +SKIP`: `PubChemAPI.get_properties_for_cids`,
`PubChemAPI.get_compound_properties_batch` and `ChEBI.get_compound`. They
passed when they were written because PubChem and ChEBI were up; PubChem
answered 503 on the re-run and two of them failed, which is how they were
found. `ChEBI.get_compound`'s example also asserted nothing --- it bound
`water` and stopped --- so it could only ever fail by raising.

All three are marked now, with output recorded live on 2026-09-22
(`get_properties_for_cids` and the ChEBI record were re-checked against the
services the same day), and `ChEBI.get_compound` asserts the name and
accession it gets back. The block in `src/conftest.py` is what keeps the
fourth one from happening: an unmarked example now fails with a message
naming the host it asked, rather than passing whenever the service is up.

The block also settles what the remaining 87 skips are for: downloads, index
builds and other mutations; the provenance and cross-reference tables, which
the Zenodo copy of `pubchem_id.db` does not carry; CompTox's exact name
searches, which build a ~20 s index on a machine that has never run one; and
the online services. Three markers did not belong there and are gone:
`pubchem_id`'s module example (`cas_to_cid("50-78-2")`), `sqlite_client`'s
module examples --- the context manager and the eight-thread pool the module
is about, both of which now run --- and `CompToxID`'s class example, which
asserted nothing and now checks the DTXSID it gets back. `sqlite_client`
joins the dataset map in `src/conftest.py` for the same reason the others
are in it: without the PubChem database its examples would ask for a
download.

---

## 29. Landed on 2026-09-22 — step 18, the documentation rebuilt (§4.9)

§4.9 asked for four things: stop hand-writing what the docstrings say, delete
`docs/examples/`, take `docs/plans/` out of the nav, and give the offline half
of the package pages. All four are done. Step 17 is what made the first
possible: every public object now has a docstring with a checked example, so
the API pages can be the docstrings.

### 29.1 What landed

- **API reference: 21 pages, each a lead and a `:::` directive.** `docs/api/`
  goes from 4 045 hand-written lines to 207. New pages: `PubChemID` (with
  `pubchem_ftp`), `CompToxID`, `ChebiSDF`, `ZeroPM`, `REACHDossierID`,
  `datasets`, `taxonomy`, `cache`, `config`, and `tools` with `utils`;
  `CheMBL` gains `mysqldump`, `PubChemView` keeps `pubchemview_parse`. The nav
  groups them as Search, Datasets, Offline databases, Online services,
  Taxonomy, Infrastructure.
- **Eight guides in `docs/guide/`**, holding what the docstrings do not:
  installing the offline databases (new — the registry, `status/plan/fetch/remove`,
  where files go, `Search(datasets=...)`, PubChem's two routes and ChEMBL's
  three); resolving identifiers with `Search` (from the old API page: sources,
  presets, online fallback, the output columns, confidence); using the local
  databases directly (closing, threads, PubChem properties and descriptors,
  CompTox names and retired CAS numbers); experimental properties from PubChem
  (`PropertyData`, `ParsedValue`, property tables); network behaviour (from the
  old `http.md`); caching, API keys and Chebifier (moved and corrected). The
  history in the old pages ("before this module…", "why it exists") was left
  out: a user needs the behaviour, and the history is in this file and the
  changelog.
- **Home and quick start lead offline.** `index.md` opens with a `Search`
  call and its real output and says nothing is downloaded until asked;
  `quickstart.md` goes datasets → `Search` → a client directly → the online
  services. `data_methods.md` is gone; its one idea is the home page's.
- **`docs/examples/` is deleted.** It was nine symlinks into `examples/`, not
  a copy (§4.9 said "byte-identical copy", which is what the symlinks looked
  like to `diff`), so nothing had rotted yet, but it would have taken a tenth
  folder to a symlink and back to add a tutorial. A MkDocs hook,
  `scripts/mkdocs_hooks.py`, adds every `.md` and `.ipynb` under `examples/`
  except `README.md` to the site at the same `examples/...` paths, ahead of
  mkdocs-jupyter so it still renders them as notebooks. The `.py` demos,
  which the symlinks had also published as notebook pages outside the nav, are
  not pages any more.
- **`docs/plans/` moved to `plans/`**, and the "Modernization" nav section is
  gone. `exclude_docs` had nothing left to exclude and went too.
- **Validation.** `mkdocs.yml` gains `validation:` at `warn` for omitted
  files, absolute and unrecognised links and missing anchors, so `--strict`
  fails on a dead link or a dead `#heading` (checked by breaking one).
- Links to the moved pages updated in `README.md`, `scripts/README.md`,
  `scripts/install_chebifier.sh`, `scripts/validate_docs_local.sh` (the
  tutorial list now reads `examples/`), `src/provesid/taxonomy.py`,
  `src/provesid/data/README_PUBCHEM.md`, `pyproject.toml`, two example READMEs
  and a demo. `TESTING.md` says how to build the docs. The deploy workflow's
  path filter gains `scripts/mkdocs_hooks.py`; no test job was touched.

### 29.2 Rendering the whole package found two parser problems

Only six modules had ever been rendered. Rendering all of them failed
`--strict`, and the warnings were two real defects in how the docs had been
parsed all along, not noise:

- **Four modules are NumPy-style.** `chembl`, `comptox`, `reach` and `zeropm`
  use `Parameters\n----------`; `mkdocs.yml` said `docstring_style: google`.
  `CheMBL` had been rendered that way on its page all along: its sections
  came out as Markdown setext headings, and ZeroPM's example outputs
  (`['Alcohol', 'ETHANOL', ...]`) as reference links, which is what
  `--strict` caught. A per-page `docstring_style: numpy` does not work:
  mkdocstrings-python applies it to the directive's own object and not to a
  class's members, because the package is loaded, and parsed, once. The fix
  is `docstring_style: auto` with Google as the default, which detects the
  style per docstring.
- **Every wrapped `Returns:` rendered as several returned values.** griffe's
  Google parser starts a new item at each line of the section's base
  indentation, so a description wrapped over four lines was four rows in the
  Returns table — `Search.search` had 17, `cas_to_detail` 17. About 200 docstrings,
  private ones included, were affected. `returns_multiple_items: false` makes the section one value.
  `returns_named_value` stays at its default: turned off, it splits the first
  line on its first colon, which would cut 32 prose descriptions (`"Mapping
  of canonical ``CHEBI:<id>`` string…"`) in half.

`auto` with per-style options triggers a mkdocstrings-python 2.0 quirk: it
fills in defaults for the Sphinx style and then warns, once per module, that
griffe's Sphinx parser does not take one of them (`warn_missing_types`). The
package has no Sphinx-style docstrings. The hook filters that one message and
says why; remove the filter when mkdocstrings-python stops emitting it.
`warn_missing_types: false` is set for Google and NumPy: most online clients'
signatures carry no annotations, and the docstrings give the types in prose.

### 29.3 What the old pages said that was wrong

- `api-keys.md` said the config file's permissions "are set to be readable
  only by your user account". `config.py` sets no permissions. The guide now
  says so, and says the stored key outranks `CAS_API_KEY` (§28.4 is how that
  bit this machine).
- `chembl.md`'s manual download ran `cd src/provesid/data`; datasets have not
  lived there since the per-user data directory. The section is gone;
  `datasets.fetch` and `CheMBL(source=...)` are the documented routes.
- `http.md` listed ClassyFire among the clients on the shared transport; it
  is still on raw `requests` (step 7 is postponed).
- `advanced_caching.md` implied every online client is cached. The online
  `ChEBI` client is not, and `use_cache=False` is a constructor argument on
  four clients and a per-call one on `ClassyFireAPI`.
- `search.md`'s confidence table left out OPSIN's base (0.97) and its column
  table `opsin_smiles`.
- `api/index.md` still imported `PubChemPUGViewAPI`.

### 29.4 Validation

- `mkdocs build --strict`: clean, 0 warnings.
- Every tutorial in `scripts/validate_docs_local.sh` round-trips through
  jupytext at its new path.
- The offline cells of the new quick start run against the installed
  databases; the example outputs in `index.md`, `guide/search.md`,
  `guide/datasets.md` (`plan()` on an empty `PROVESID_DATA_DIR`, and the
  `MissingDatasetError` text) and `guide/local-databases.md` were produced on
  this machine on 2026-09-22.
- `pytest --doctest-modules src/provesid/taxonomy.py` (the one module whose
  text changed): 10 passed, 5 skipped. No code changed.

### 29.5 Still open

- **Sphinx roles render literally.** Docstrings use `:class:`, `:func:`,
  `:data:` and `:attr:`; mkdocstrings shows them as `:data: OUTPUT_COLUMNS`.
  They should become Markdown cross-references (`` [`Search`][provesid.Search] ``)
  module by module, which `--strict` would then check.
- A typed Google return (`dict: Keyed by …`) renders the type in the *Name*
  column, since `returns_named_value` stays on (§29.2).
- Bulleted lists inside a `Returns:` description (`cas_to_detail`,
  `get_molecular_data`) run together into one paragraph: they need a blank
  line before the list.
- The tutorials were moved, not reviewed. They are step 19's.
- The quick start's online cells were not executed (the block from §28.6 does
  not apply to MyST pages, but no live run was made); `--execute` in
  `validate_docs_local.sh` does it.
