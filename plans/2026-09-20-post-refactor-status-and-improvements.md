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
| 5 | `classyfire.py` raises; rewrite its tests | not started |
| 6 | extract `pubchemview_parse.py` | **done** (§19, ahead of sequence) |
| 7 | split `pubchem.py` → `+ pubchem_id.py`; `chebi.py` → `+ chebi_sdf.py` | not started |
| 8 | `sources.py`; collapse the `Search` source ladders | not started |
| 9 | fix the fuzzy-name mislabelling | **done** (§16.1) |
| 10 | `cache.py` parameterisation and `_MISS` sentinel | partly (§17.4) |
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

### 4.3 Four SQLite clients that cannot be closed (S)

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

### 4.4 Offline→online fallback exists in exactly one place (M)

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

### 4.5 The `Search` source ladders (§5 B.2, still open) (M)

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

### 4.9 Documentation is a second, drifting copy of the code (M, step 14)

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

- **Step 7**, the `pubchem.py` / `chebi.py` splits. `pubchem.py` is 3 301 lines
  holding two unrelated things: an online PUG-REST client and an offline SQLite
  database. §8 below makes the split more attractive, since the FTP builder
  belongs beside `PubChemID`, not beside `PubChemAPI`.
- **Step 11**: `zeropm.py` logging, `reach.py` xlsx, `config.py` printing.
- **A circuit breaker** (§23.8). A `Retry-After` is information about the *host*,
  so it belongs on the host's clock as a "not before T" that every client
  respects, making a known-throttled host fail at once. The shared `RateLimiter`
  is already the right place.
- **34 public objects still have no docstring**, concentrated in `pubchem.py`
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
| 4 | `datasets.status/plan/fetch/remove`; `Search(datasets="present")` as the default | 4.1 | M |
| 5 | build the ChEMBL extract during download; `source=` on `CheMBL` | 9.7 | M |
| 6 | `close()`, context managers and thread safety on the four SQLite clients | 4.3 | S |
| 7 | `classyfire.py` raises; rewrite its tests; drop it from the docs' service lists | 4.8 | S |
| 8 | `cache.py`: persistent dir, parameterised service functions, key version | 4.7 | S |
| 9 | **`pubchem_ftp.py`: the FTP identifier builder** | **8** | **L** |
| 10 | `PubChemID.descriptors()` — RDKit descriptors on demand | 8.6 | M |
| 11 | `CheMBL(source="mysql")` — the streaming route, verified against step 5 | 9.7 | M |
| 12 | split `pubchem.py` → `+ pubchem_id.py`; `chebi.py` → `+ chebi_sdf.py` | 4.12 | M |
| 13 | `sources.py`; collapse the `Search` ladders | 4.5 | M |
| 14 | `Search(online_fallback=...)` as a row in the source table | 4.4 | M |
| 15 | `Search.PRESETS` | 4.6 | S |
| 16 | circuit breaker on the shared `RateLimiter` | 4.12 | M |
| 17 | docstrings with examples, module by module | 4.12 | L |
| 18 | rebuild `docs/`; delete `docs/examples/`; drop `docs/plans/` from the nav | 4.9 | M |
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
  started. `compact()` is deliberately the whole of step 1.

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
