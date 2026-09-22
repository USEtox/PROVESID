# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **`Search(online_fallback=True)`: ask PubChem and CACTUS what no offline
  database holds.** Off by default, so a run still opens no socket. When on,
  a query that produced no candidate from any offline source, and only such
  a query, is asked of PubChem's PUG-REST service and the NCI/CADD
  resolver (CACTUS). Their answers are pooled, clustered and scored like
  offline ones, and each service counts as one vote in `n_source_support`:

  ```python
  df = Search("cas", online_fallback=True).search(["50-78-2", "1912-24-9"])
  df["source"]                  # "PubChemID", ..., "PubChem (online)" or "CACTUS"
  df.attrs["online_fallbacks"]  # queries that went online
  df.attrs["online_resolved"]   # ...and those the services answered
  ```

  The two services are two more columns in `provesid.sources.LOOKUPS`,
  `pubchem_online` and `cactus` (`ONLINE_SOURCE_KEYS`), with rows for CAS,
  name, SMILES, InChI, InChIKey and, on PubChem only, DTXSID. Formula queries are never
  sent online. A service's "not found" is a miss. Any other failure is logged
  at WARNING and costs only that service's vote. Each fallback is logged at
  DEBUG. New helpers: `PubChemAPI.get_cids_by_inchi`, which POSTs because an
  InChI cannot travel in a URL path, and `tools.candidate_from_pubchem_online`
  and `tools.candidate_from_cactus`. See
  `examples/search/online_fallback_demo.py`.

- **ChEMBL from its MySQL dump: `CheMBL(source="mysql")`.** Builds the same
  2.4 GiB extract as the default route from ChEMBL's 2.1 GB MySQL dump, read
  as a stream, so the 27.7 GiB release never exists on disk:

  ```python
  CheMBL(source="mysql")                                  # 2.1 GB down, ~4.5 GiB free
  CheMBL.build_from_mysql_dump("chembl_37_mysql.tar.gz")  # from a dump you already have
  CheMBL.extract_digest(path)                             # {table: (rows, hash)}
  ```

  The new `provesid.mysqldump` module reads exactly what `mysqldump` writes
  and raises `DumpFormatError` on anything else rather than guess a value.
  Column types come from the dump's `CREATE TABLE` and map to the affinity
  SQLite gives them, so values are stored the way the `sqlite` route stores
  them: `decimal` 4.00 becomes the integer 4 in both. `extract_digest`
  fingerprints each table with a row count and an order-independent content
  hash that includes each value's type, and is how the routes were compared.
  On the real ChEMBL 36, written out as a mysqldump, all eight tables
  (14.6 M rows) digest identically to the `sqlite` route's extract. The build
  took 6 min and 185 MB of memory. Provenance now records `source_format`.
  The default stays `source="sqlite"` until the same comparison has been run on
  EBI's own dump. `examples/chembl/download_source_demo.py` shows all three
  routes.

- **Computed descriptors on demand: `PubChemID.descriptors()`.** The database
  no longer stores XLogP, TPSA or the atom and bond counts; this computes them
  when asked, and names the model that answered:

  ```python
  db = PubChemID()
  db.descriptors(2244)                                  # RDKit, from the stored SMILES
  # {'CID': 2244, 'Source': 'rdkit', 'MolLogP': 1.3101, 'TPSA': 63.6, ...}
  db.descriptors(2244, ["XLogP"], source="pubchem")     # PubChem's, over PUG-REST
  db.descriptors_table(cids)                            # a row per CID
  rdkit_descriptors("CCO")                              # any structure, no database
  ```

  `source="rdkit"` (default) needs no network for a compound the database holds
  and costs about 0.5 ms per compound; for one it does not hold, only the SMILES
  is fetched. Names are PubChem's where the quantity is the same, except the
  logP: RDKit's is Crippen's, not XLogP3, so it is `MolLogP`, and asking RDKit
  for `XLogP` or `Complexity` raises an error saying where to get it. The
  values are RDKit's, not PubChem's Cactvs values. Against 20 000 compounds,
  heavy atoms and charge agree for all, donors for 94%, rotatable bonds 74%,
  TPSA 70% (counted with S and P, as PubChem does), acceptors 63%.
  `examples/pubchem/descriptors_demo.py`.

- **`pubchem_id.db` is built from PubChem's FTP site, and says where it came
  from.** The new `provesid.pubchem_ftp` module builds the database
  `PubChemID` reads from a dated monthly snapshot of `Compound/Extras/`, and
  `PubChemID(source=...)` chooses between that and the Zenodo download:

  ```python
  from provesid import PubChemID
  from provesid.pubchem_ftp import build_pubchem_id_db, list_releases

  db = PubChemID()                    # source="ftp" (default): build if missing
  db = PubChemID(source="zenodo")     # download the 2.2 GiB prebuilt copy
  list_releases()                     # ['2026-09-01', ..., 'current']
  build_pubchem_id_db(release="2026-09-01", keep_downloads=True)

  db.provenance()["release"]          # '2026-09-01', and every file's MD5
  db.xrefs(2244)                      # {'chebi': [...], 'dtxsid': [...], ...}
  ```

  *Which compounds.* Those PubChem's own identifier file maps to a CAS number,
  each checked against the check digit: 1 430 379 in the 2026-09-01 snapshot,
  with 123 CAS rows rejected. The Zenodo database found its CAS numbers by
  running `\d{2,7}-\d{2}-\d` over free-text synonyms, and 17 379 of them
  (1.25%) are not valid CAS numbers.

  *How.* Eight files, each downloaded resumably, checked against the MD5
  PubChem publishes beside it, streamed once to keep the compounds in scope and
  deleted, so the free disk needed is the database plus the largest file
  (7.4 GB), not the 15.4 GB total. The database is built beside its destination
  and moved into place only when complete. `include_inchi=False` skips the
  InChI file and computes InChI with RDKit; `include_synonyms=False` leaves out
  the synonyms; `keep_downloads=True` keeps the files, and a later build reuses
  any whose MD5 still matches.

  *What it records.* A `provenance` table (release, snapshot timestamp, build
  time, row counts, builder version) and a `provenance_files` table (URL, MD5,
  bytes, lines read and kept, per file). `PubChemID.provenance()` returns both.

  *What it adds.* `monoisotopicmass`, and an `xrefs` table of DSSTox, ChEBI,
  ChEMBL, EC and UNII identifiers from the same file, which `PubChemID.xrefs()`
  reads.

  *Molecular weight* is in none of the FTP files, so it is computed from the
  formula (`provesid.pubchem_ftp.molecular_weight`) with
  `provesid.pubchem_ftp.ATOMIC_WEIGHTS` --- IUPAC 2005 with twelve elements
  fitted to PubChem's own weights --- and from the SMILES for isotopically
  labelled compounds, whose formula PubChem writes unlabelled. PubChem's
  weights cannot be reproduced exactly: it rounds them to between zero and
  three decimals by a rule no function of the formula fits. This agrees with
  PubChem to the second decimal for most compounds and is the more precise
  number where PubChem rounds harder.

  `scripts/build_pubchem_id_db.py` is now a command-line wrapper over the
  builder, for refreshing the Zenodo copy; the CSV-and-regex version it
  replaces, and its duplicate in `src/provesid/data/`, are **removed**.
  `examples/pubchem/ftp_build_demo.py` builds a miniature release served on
  localhost.

- **The cache is in a cache directory now, one per service, with a version in
  every key.** Three changes to `provesid.cache`:

  ```python
  import provesid

  provesid.get_cache_info(service="pubchem")["cache_directory"]
  # '/home/you/.cache/provesid/pubchem'   (was /tmp/provesid_cache/pubchem)

  provesid.clear_cache(service="pubchem")     # one service
  provesid.clear_cache(all_services=True)     # all of them
  provesid.get_all_cache_info()               # {'global': ..., 'pubchem': ...}
  ```

  *It persists.* The default was `tempfile.gettempdir()/provesid_cache` while
  `docs/advanced_caching.md` promised the cache survived a restart. Most Linux
  systems clear `/tmp` on boot, so it did not. It now uses
  `platformdirs.user_cache_dir` through the new `provesid.utils.user_cache_path`
  — `~/.cache/provesid/` on Linux, `~/Library/Caches/provesid/` on macOS,
  `%LOCALAPPDATA%\USEtox\provesid\Cache\` on Windows — overridable with
  `PROVESID_CACHE_DIR`, and separate from `PROVESID_DATA_DIR` because cached
  responses are disposable and 2.4 GiB datasets are not. Entries written under
  the old default are not migrated; they were never reliably there to migrate.

  *The service list is data.* Fourteen near-identical module-level functions
  (`clear_pubchem_cache`, `clear_cas_cache`, … and seven `get_<svc>_cache_info`
  twins) plus `export_service_cache` and `import_service_cache` are **removed**.
  Every cache function — `clear_cache`, `get_cache_info`, `get_cache_size`,
  `export_cache`, `import_cache`, `set_cache_warning_threshold`,
  `enable_cache_warnings` — now takes `service=`, with the names in
  `provesid.CACHE_SERVICES`; `get_all_cache_info()` and
  `clear_cache(all_services=True)` cover the lot. An unknown service raises
  `ValueError` instead of silently writing to the global cache, where the
  matching `clear_cache` would never look. `get_service_cache()` returns the
  underlying `CacheManager`, and builds it on first use, so importing
  `provesid` no longer creates eight cache directories.

  *Keys carry a version.* Pickle stores an instance's `__dict__`, so an entry
  written before a cached dataclass gained a field comes back as an object
  missing that attribute, and every caller that reads it raises on a machine
  where only the package version changed. There are now three versions folded
  into the key, and you bump the narrowest one that covers the change:
  `@cached(version=N)` for one function, a client's `CACHE_SCHEMA_VERSION` for
  one client, `provesid.cache.CACHE_KEY_VERSION` for everything. Retired
  entries are unreachable, not deleted; `clear_cache` reclaims the space.

  Also fixed: `export_cache` dropped any entry whose cached value was `None`
  (`is not None` where it meant `is not _MISS`), so a legitimately cached
  `None` did not survive an export/import round trip.

- **The four SQLite clients can be closed, used in a `with` block, and queried
  from several threads.** `PubChemID`, `CompToxID`, `ZeroPM` and `CheMBL` each
  opened a connection in `__init__` and closed it in `__del__`, and nowhere
  else. All four now inherit `provesid.SQLiteClient`:

  ```python
  with PubChemID() as db:                 # released at the end of the block
      inchi = db.cas_to_inchi("50-78-2")

  db = CompToxID()
  db.close()                              # or by hand; idempotent
  db.closed                               # True

  with PubChemID() as db:                 # a connection per thread
      with ThreadPoolExecutor(8) as pool:
          rows = list(pool.map(db.get_by_cas, cas_numbers))
  ```

  The threading case is the one users met first: a single connection is bound
  to the thread that created it, so the obvious way to resolve ten thousand
  CAS numbers against a local 2.2 GB database died on the first worker with
  `sqlite3.ProgrammingError: SQLite objects created in a thread can only be
  used in that same thread`. Each thread now gets its own connection and
  cursor, opened on its first query; `close()` closes all of them, from
  whichever thread calls it.

  Per-thread connections rather than `check_same_thread=False`, which removes
  the *check* and not the problem: every thread would still share one
  `self.cursor`, and these clients `execute` in one statement and `fetchone`
  in the next, so two threads interleaving those would read each other's rows.
  A wrong answer is a worse failure than the exception it replaces.

  Threads are now *possible*, which is not the same as faster: a tight loop of
  nothing but local lookups is slower on a pool than in a plain loop (5 000
  `get_by_cas` calls: 0.29 s serially, 14.65 s on eight threads), because each
  query takes tens of microseconds and the GIL handoff around it costs more.
  That is `sqlite3` under CPython, not this change — plain `sqlite3.connect`
  per thread measures the same. A pool pays when each item also waits on
  something: 400 lookups each followed by 20 ms of waiting took 8.52 s
  serially and 1.17 s on eight threads.

  `Search` is a context manager too, and closes the source clients it
  constructed — but not one passed to its constructor, which belongs to the
  caller. Querying a closed client or a closed `Search` raises
  `DatabaseClosedError`, a `RuntimeError` subclass that names the client and
  the file, rather than sqlite3's `Cannot operate on a closed database`.
  `__del__` remains as a backstop, so code that never closes anything behaves
  as it did.

- **Installing ChEMBL now costs 2.4 GiB, not 27.7 GiB: `CheMBL(source=...)`.**
  `CheMBL.compact()` could already shrink a release *already on disk*, which
  left the worst case untouched — a machine that had never had ChEMBL still
  installed 27.7 GiB and had to be told to shrink it afterwards.

  The download now finishes by building the extract and deleting the release
  it came from:

  ```python
  CheMBL()                  # 5.8 GB transferred, 2.4 GiB installed (default)
  CheMBL(source="full")     # ...or keep all 74 tables, 27.7 GiB
  ```

  Both routes transfer the same archive: the choice is what stays on disk. The
  extract is verified against the release — `quick_check`, a row count per
  table, and 500 compounds compared column by column — *before* the release is
  deleted, so the file that survives is the one that was checked.

  `source` only describes a download. A release already on disk is opened as
  it is: shrinking 27 GB is `compact()`'s job, not something a constructor
  should do unasked. Opening a full release now logs the one line that says
  the option exists.

  A compaction that fails after a successful download is logged, not raised.
  The release is in place by then and answers every query, so the failure
  costs disk rather than function, and throwing away a 5.8 GB download over a
  step that can be repeated with one call would be the worse trade.

- **`datasets.plan()` now tells the truth about ChEMBL.** Its registry entry
  said 27.7 GiB installed; it is 2.4 GiB, from a 5.8 GB download, with a peak
  of 33.4 GiB while the archive and the full release both exist. That peak is
  the number a laptop actually fails on, and it is now the one that is
  reported: for the four default sources, 8.9 GiB downloaded, **6.3 GiB
  installed** (down from 31.6 GiB), 37.3 GiB needed at the worst moment.

- **A dataset manager, and `Search` no longer downloads 32 GB behind your
  back.** On a clean machine `Search("cas").search("50-00-0")` constructed
  four source clients that each default to `auto_download=True`, so one CAS
  lookup fetched ChEBI, CompTox, PubChem and ChEMBL — about 32 GB, of which
  ChEMBL is 87% and only *enriches* a structure the other sources already
  found. Nothing announced the total, nothing asked, and nothing offered to
  proceed with the sources already present. This package is for researchers on
  laptops, where 32 GB is often the whole free disk.

  `provesid.datasets` now describes the five datasets without opening them, so
  the decision can be made before the first byte moves:

  ```python
  from provesid import datasets

  datasets.status()                    # present? size on disk? release? path?
  datasets.plan(["pubchem", "chebi"])  # what a download would transfer
  datasets.fetch("pubchem")            # install by name; skips what is present
  datasets.remove("chembl")            # reclaim the space, by name
  ```

  `status()` and `plan()` return DataFrames and read only filenames and
  `stat`, so they are instant with 30 GB of ChEMBL in the directory. `plan()`
  also reports the *peak* requirement, which for ChEMBL exceeds the installed
  size by the 5.8 GB archive sitting beside the database it extracts into — a
  laptop with 28 GB free and 27.7 GB of ChEMBL to install still fails.

- **`Search(datasets=...)`, defaulting to `"present"`.**

  | value | behaviour |
  |---|---|
  | `"present"` | use whatever is on disk; name each missing dataset, its size and the `fetch` call that installs it. **The new default.** |
  | `"auto"` | the behaviour before this release — download whatever is missing |
  | `"required"` | raise `MissingDatasetError` in the constructor, before any query, naming the missing datasets and the exact `fetch` call |

  `Search` already degraded gracefully when a source was unavailable and
  already warned that confidence is not comparable across runs, so `"present"`
  needed no new machinery in the resolver — only that the clients stop
  downloading behind the caller's back.

- **One resumable, checksummed downloader for every bulk dataset:
  `provesid.datasets.download_file`.** Five modules had each grown their own
  copy of "stream the response into a temporary file with a progress bar", and
  all five shared the same three defects — no retry, no resumption, no
  checksum. An interrupted 5.8 GB ChEMBL download started again from zero, and
  so did an interrupted 2.2 GB PubChem one.

  Corruption was caught unevenly, too. The two gzipped downloads got truncation
  detection for free from gzip's CRC and length trailer; the three plain SQLite
  downloads had none. A file truncated in its interior opens cleanly and fails
  much later, on the first query that touches a missing page — which a user
  reads as a data problem, not a download problem.

  `download_file` streams into a `.part` file beside the destination, resumes
  it with an HTTP `Range` request, and applies four checks in order, each of
  which leaves the destination untouched if it fails:

  1. the byte count against the size the server declared;
  2. an MD5, from `expected_md5` or fetched from `checksum_url` — PubChem's FTP
     mirror publishes an `.md5` beside every file;
  3. the caller's own `verify` callback, handed the finished file;
  4. an atomic rename.

  ```python
  from provesid.datasets import download_file

  download_file(url, "/data/CID-SMILES.gz", checksum_url=url + ".md5")
  ```

  A file rejected by its checksum or by `verify` is deleted rather than kept:
  it is already complete, so resuming it would fail the same check again. A
  transfer that merely stopped is kept, and the next call continues from it —
  but only when it came from the same URL, recorded in a `.part.source` marker
  beside it. A partial left by a different download is discarded instead of
  being spliced onto this one, which a checksum would catch but only PubChem's
  FTP mirror publishes one.

  `PubChemID.download_database`, `CompToxID.download_database`,
  `ZeroPM.download_database`, `CheMBL.download_database` and
  `ChebiSDF.download_sdf` are now one call each into this.

- **`CheMBL.compact()` — a ChEMBL release without the 27 GB you never read.**
  A full ChEMBL release is ~30 GB across 74 tables. PROVESID opens eight of
  them (`molecule_dictionary`, `compound_structures`, `compound_properties`,
  `molecule_synonyms`, `molecule_hierarchy`, `chembl_id_lookup` and the two
  pesticide tables) and reads no bioactivity data at all: `activities` alone is
  24.3 M rows that nothing in this package has ever queried.

  Worse, a quarter of the whole database is a single column. `molfile` — the
  MOL block for every compound — totals 7.72 GB, and `get_compound` was the
  only thing that selected it, for callers that never used it.

  `CheMBL.compact()` copies the eight tables, drops `molfile`, keeps only the
  `COMPOUND` rows of `chembl_id_lookup`, and rebuilds just the indexes this
  package's queries need. Measured on ChEMBL 36: **29.74 GB → 2.60 GB, 91.3%
  smaller**, in about a minute, with every public method returning the same
  compounds — verified over 3 000 randomly chosen compounds.

  ```python
  from provesid import CheMBL

  CheMBL.compact()                      # build chembl_37_provesid.db beside the original
  CheMBL.compact(remove_source=True)    # …and delete the 30 GB once it verifies
  CheMBL.compact(keep_inchi=False)      # ~1.55 GB; see the caveat below
  ```

  The source is opened read-only, the extract is written to a temporary file
  and checked against its source — SQLite's `quick_check`, a row count per
  table, and 500 compounds re-read from both databases — before it replaces
  anything. `remove_source` deletes the original only after that passes.

  A later `CheMBL()` opens the extract in preference to a full release of the
  same number, so nothing else has to change. `CheMBL.is_compact` and
  `CheMBL.provenance` report what was opened; an extract records the release it
  came from, the source file and its size, the build time, the PROVESID version
  and its row counts. It also records a schema version, so a future PROVESID
  that needs a ninth table warns and names the rebuild command instead of
  failing with `no such table`.

  `keep_inchi=False` saves a further ~1.0 GB by dropping `standard_inchi` and
  its indexes, but `search_by_inchi` then has no column to match and `Search`
  loses the InChI it reads straight from ChEMBL. Leave it alone unless disk is
  genuinely short.

- **One shared HTTP transport, and `Retry-After` honoured for the first time.**
  Every web-API client used to carry its own copy of "pause, request, decide
  what the status code meant, maybe give up". The six copies had drifted
  apart: only `pubchemview` retried at all, none honoured `Retry-After`, and
  `chebi`, `cascommonchem` and `opsin` had no rate limiting whatsoever.

  `provesid.http.HTTPClient` is now the single place that decides when to ask
  again. It paces requests, retries HTTP 429, 5xx, timeouts and connection
  errors with exponential back-off capped at `max_backoff`, and prefers a
  `Retry-After` the service sent — in seconds or as an HTTP date — over its own
  curve. The pacing applies to retries too: a service already shedding load is
  not asked again faster than a healthy one. No raw `requests` exception
  reaches a caller.

  Every web-API client is migrated: `resolver`, `pubchemview`, `pubchem`,
  `chebi`, `cascommonchem` and `opsin`. Only `classyfire` still holds
  `requests` calls of its own, and only because its service has been down
  since February 2023. The bulk database downloads
  (`ChEBISDF.download_sdf`, `chembl`, `comptox`, `zeropm`, `pubchem`'s Zenodo
  dataset) are deliberately left alone: a hundred-megabyte streamed file wants
  resumption and checksums, not a 5-per-second pacer.

  Nothing changes for a caller. Each client passes its own exception classes to
  the transport, so `except NCIResolverError` and `except PubChemViewError`
  work exactly as before; what is new is that those classes also descend from
  `provesid.http.ServiceError` / `NotFoundError` / `ServiceTimeoutError`, so a
  caller can now catch every service at once.

- **Pacing is now shared per host, not per client object.** PubChem publishes
  five requests per second **per IP**, and a `PubChemAPI` plus a `PubChemView`
  in one process — which is the ordinary way to use this package, and what
  `Search` does — each kept their own clock. Each believed it was pacing
  correctly while together they asked twice as fast as PubChem allows.

  `min_interval` stays the client's own promise about how fast it will ask; the
  clock it is measured against now belongs to the host.
  `provesid.http.RateLimiter` holds one host's clock behind a lock, and
  `provesid.http.host_limiter(url_or_host)` hands every client aimed at that
  host the same one. ChEBI and OPSIN are both served from `www.ebi.ac.uk`, so
  they share a clock too — which is correct, the limit being the host's.

  `last_request_time` still reports when *that* client last asked. The
  effective rate for a host is set by its most impatient client, so this stops
  two clients from doubling a limit but not one client configured at
  `min_interval=0.01` from exceeding it alone.

- **`max_elapsed`, a ceiling on the total time spent waiting between retries** —
  "do not make the caller wait longer than this". PubChem answers a throttled or
  blacklisted IP with `Retry-After: 30`, which the transport honours, so three
  retries would be a ninety-second call;
  `tests/test_pubchemview.py::TestPubChemView::test_error_handling` measured
  92.5 s. Dropping the retry would be worse, because a caller resolving ten
  thousand compounds loses one to every transient 503.

  Both PubChem clients set `provesid.pubchem.RETRY_WAIT_BUDGET = 10.0`, which
  leaves the cheap curve intact — 1 + 2 + 4 for a transient 500 or a timeout —
  while declining the 30-second throttle. That block does not lift in thirty
  seconds: measured on 2026-09-19, waiting the full thirty and asking again
  returned the same 503, so the wait buys nothing while every call pays it. A
  throttled `get_compound_by_cid` now fails in 1.1 s. Raise the budget with
  `api._http.max_elapsed = 180` for an unattended bulk job that would rather
  wait.

  `max_elapsed` defaults to `None`, so every other client is unaffected.

- **`HTTPClient(session=...)`**, so a client that keeps a `requests.Session`
  for connection pooling and persistent headers makes its calls through it.
  `chebi` uses this: its `User-Agent` and `Accept` live on the session, and an
  ontology walk reuses one connection. Without a session the transport calls
  `requests.get`/`requests.post` through the module, which is what lets a test
  stub them.

- **Exceptions carry the response detail.** `ServiceError` now has
  `status_code`, `url` and `response`, keyword-only and defaulting to None, so
  `raise PubChemError("...")` by hand is unchanged. CAS Common Chemistry needs
  this: it reports a rejected key as 401 and an unknown CAS number as 404, and
  both have to become different strings in the dict it returns.

- **`retry_exhausted_cls`**, raised when a transient condition outlives the
  retry budget, as distinct from a permanent error. PubChem passes
  `PubChemServerError`, because its callers catch that one to skip and
  `PubChemError` to fail. Defaults to `error_cls`, so no other client notices.

- **New exception classes**, all exported from `provesid` and all descending
  from the shared bases: `ChEBINotFoundError`, `ChEBITimeoutError`,
  `CASCommonChemError`, `CASCommonChemNotFoundError`,
  `CASCommonChemTimeoutError`, `OPSINError`, `OPSINNotFoundError`,
  `OPSINTimeoutError`. `PubChemServerError` and `PubChemTimeoutError` are now
  exported too. `ChEBIError` was a bare `Exception` and is now a
  `ServiceError`, so `except ServiceError` catches ChEBI as well.

- **`tests/test_cascommonchem_offline.py`**, 19 offline tests.
  `tests/test_cascommonchem.py` skips itself entirely without a CAS API key, so
  the module that gained the most in this change had no coverage on a developer
  machine at all. The key is only a header, so a stub key and a stubbed
  `requests.get` cover the lot.

- **`provesid.cache.is_empty_result`**, the `skip_if` predicate that keeps an
  empty answer out of the cache. It was written twice, identically, in
  `pubchem` and `pubchemview`; it now lives beside `is_failure_result`, which
  is where `skip_if` belongs.

- **One parser for PUG-View's free-text values, with typed output.** PubChem
  reports every experimental property as prose written by whoever deposited it —
  `"138-140 °C"`, `"8.5X10-5 mm Hg at 25 °C"`, `"greater than or equal to 100
  mg/mL"`, `"Vapor pressure, kPa at 20°C: 24"` — and PROVESID had two parsers
  for it that disagreed with each other: a 25-line one behind `PropertyData` and
  a 390-line property-specific regex cascade behind `get_property_table`, both
  returning strings.

  The new `provesid.pubchemview_parse` module holds the only one.
  `parse_value(text, heading)` returns a `ParsedValue` carrying:

  - `value`, or `value_min`/`value_max` for a range (a single value fills both
    bounds too, so a numeric filter needs no special case);
  - `unit`, spelling normalised — `torr` and `mm Hg` both report as `mmHg`;
  - `value_si`/`value_min_si`/`value_max_si`/`unit_si`, the same quantity in
    `K`, `Pa`, `kg/m³`, `mol/m³`, `Pa·s`, `m²/s` or `N/m`;
  - `temperature_c`, the temperature the measurement was made *at*, which is a
    condition rather than the value and was the largest source of wrong answers;
  - `operator` (`>`, `<`, `>=`, `<=`, `~`) when the entry bounds the value;
  - `qualitative` (`insoluble`, `miscible`, `negligible`) when a word replaced
    the number;
  - `text`, always — nothing the parser cannot read is lost.

  The `heading` argument settles what a string alone cannot: a bare `138` under
  "Melting Point" is 138 °C, while a bare `1.19` under "LogP" is dimensionless
  and complete.

  Measured on 366 real value strings across ten compounds and ten properties:
  90% yield a number, 8% a qualitative term, and the remaining 2% are entries
  where PubChem states no value at all. Aspirin's melting point entries, which
  are deposited variously as `135 °C` and `275 °F`, now all read 408.15 K.

  Nothing is invented: an unrecognised unit is reported as written with no SI
  conversion, and `%`, `ppm` and `ppb` are never converted, since a composition
  needs a density to become a concentration. `M` for molar is not recognised
  either — it is indistinguishable from metres and from a stray capital.

- **`get_property_table` now returns numbers.** New columns `ValueMin`,
  `ValueMax`, `ValueSI`, `ValueMinSI`, `ValueMaxSI`, `UnitSI`, `Operator`,
  `Qualitative` and `Heading` join the existing ones, and the full list is
  published as `PubChemView.PROPERTY_TABLE_COLUMNS`.

- **Every PUG-View heading works, not only the experimental ones.** Both
  response parsers hard-coded the path *Chemical and Physical Properties →
  Experimental Properties → \*​*, but PUG-View nests a requested heading
  wherever it sits in the compound's table of contents: "GHS Classification"
  under *Safety and Hazards*, "Drug Indication" under *Drug and Medication
  Information*. Those returned a full response and an empty result — data
  reported as absent while it was right there. `_iter_information` walks the
  record instead, and `get_property_table(2244, "GHS Classification")` now
  returns its 18 rows.

- `get_property_summary` gained `numeric_values`, `numeric_values_si` and
  `units_si`; `export_properties_to_dict` gained the parsed numbers as flat,
  JSON-serialisable keys. `PropertyData` gained `heading` and `parsed`.

- `PubChemView.CACHE_SCHEMA_VERSION` is part of the client's cache key, so an
  entry written before the shape of a result changed becomes unreachable rather
  than being deserialised into the wrong structure. This mattered immediately: a
  `PropertyData` pickled by the previous version restores without a `parsed`
  attribute, and everything reading it would raise. Bump the constant whenever a
  cached return shape changes.

- **Bulk property retrieval.** `PubChemAPI.get_properties_for_cids(cids,
  properties)` asks PubChem about a whole list of compounds in one request
  rather than one request per compound: 450 compounds now cost 3 requests and
  about 2 seconds, where the per-CID loop cost 450 requests and a minute and a
  half. Lists are split into `PROPERTY_CHUNK_SIZE` (200) CIDs per request, so a
  retry redoes a chunk rather than everything.

- **POST for long identifier lists.** PUG-REST caps a URL at roughly 2000
  characters, which a few hundred CIDs exceed. `_build_post_url` builds the
  same endpoint without the identifier segment, and the bulk property path
  switches to POST — identifiers in the body — once the list outgrows
  `URL_IDENTIFIER_LIMIT`. Below that it keeps using GET, which PubChem prefers.

- **Offline-first property lookup.** `PubChemID.properties(cid, properties)`
  reads the local SQLite database first and consults PUG-REST only for what it
  cannot answer, per the two-stage lookup the development principles ask for.
  For most property work this means no network traffic at all: 473 of 1000 low
  CIDs were served from disk in the same call, and the remaining 527 in three
  batched requests.

  `properties_for_cids()` does the same for a list, and `properties_table()`
  returns a DataFrame with a row for every CID asked about. Each record carries
  a `Source` key reading `offline`, `online` or (in the table) `missing`.
  `use_online_fallback=False` keeps a lookup strictly local.

  `PubChemID.OFFLINE_PROPERTIES` maps the 16 property names the database can
  answer — formula, weight, exact mass, isomeric SMILES, InChI, InChIKey, IUPAC
  name, title, XLogP, TPSA, complexity, charge and the H-bond, rotatable-bond
  and heavy-atom counts — to their columns. Anything else
  (`MonoisotopicMass`, `ConnectivitySMILES`, the 3D descriptors, the patent and
  literature counts) is online-only, and asking for one sends the whole request
  online rather than assembling a row from two different PubChem snapshots.

  A property with no value is absent from the result rather than `None`, which
  is how PubChem reports it — `'XLogP' not in result` means PubChem computes no
  logP for that compound, not that the lookup fell short. Values are normalised
  to one type across both sources, since PUG-REST returns `MolecularWeight` as
  a string where the database holds a float.

- `PubChemID(api=...)` accepts the `PubChemAPI` used for fallback; one is
  created on first use otherwise, so a strictly offline session never builds
  one.

### Fixed
- **`datasets.status()` and `remove("chembl")` could not see an interrupted
  ChEMBL download.** The registry looked for `chembl_NN_sqlite.tar.gz`, but the
  download is saved as `chembl_NN.db.tar.gz`. A half-finished 5.8 GB `.part`
  was therefore neither counted nor removed. Both routes' real archive names
  are now listed.

- **`CheMBL.search_by_name` scanned all 2.9 M compounds on every call, and its
  truncated results were not reproducible.** The method asked one question of
  two tables:

  ```sql
  SELECT DISTINCT md.molregno FROM molecule_dictionary md
  LEFT JOIN molecule_synonyms ms ON md.molregno = ms.molregno
  WHERE LOWER(md.pref_name) = LOWER(?) OR LOWER(ms.synonyms) = LOWER(?)
  LIMIT ?
  ```

  Neither arm can use an index in that shape — `LOWER(col) = ?` defeats an index
  on `col`, and an `OR` whose arms live in different tables across a `LEFT JOIN`
  defeats indexing altogether — so SQLite scanned `molecule_dictionary` in full
  for every lookup. Adding expression indexes changes nothing on its own; the
  query shape is the problem.

  It is now a `UNION` of two independently indexable lookups, one per table.
  Measured on ChEMBL 36, exact lookups:

  | | per call |
  |---|---:|
  | full release, old query | 764 ms |
  | full release, new query | 221 ms |
  | `compact()` extract, old query | 560 ms |
  | **`compact()` extract, new query** | **10 µs** |

  The 10 µs comes from the `lower(pref_name)` and `lower(synonyms)` expression
  indexes that `CheMBL.compact()` already builds into every extract; a full
  release has no such indexes and still scans, but it scans two small queries
  instead of a join. `exact=False` cannot be indexed either way — a
  leading-wildcard `LIKE` is a scan by construction — and improves from 665 ms
  to 275 ms on the extract for the same reason the join is gone.

  The rewrite returns the same compounds: 0 differences over 200 real names and
  25 substring fragments.

- **`CheMBL.search_by_name` returned a different set of compounds run to run.**
  `LIMIT` with no `ORDER BY` returns whatever the query plan happens to produce,
  so a name matching more compounds than `limit` silently changed *which* ones
  came back — after a `VACUUM`, an index change or a SQLite upgrade. Comparing
  a `compact()` extract against its source found 25 of 60 names differing in
  order, and none in content, which is what put the defect in view.

  Results are now ordered by `molregno`, so a truncated result is the `limit`
  lowest `molregno` values and is the same on every call and every copy of the
  database.

- **ChEMBL's archive extracted straight onto the database it was replacing.**
  `download_database(force=True)` renamed the extracted file onto `db_path`
  and validated afterwards, so an extraction or validation failure destroyed
  the release already on disk — and then deleted what was left. The archive now
  extracts to `<db_path>.incoming`, answers a query there, and is moved into
  position only once it has.

- **Three downloads renamed the file into place before checking it.**
  `CompToxID.download_database` moved the downloaded file onto the destination
  and only then looked for the `chemicals` table; a failed check therefore left
  the broken file where the working database had been. `ZeroPM` had the same
  shape and deleted the database on failure, leaving nothing at all.

  Every download now validates the `.part` file *before* the rename, which is
  what `PubChemID` alone already did. A failed download cannot replace or
  remove a working database.

- **OPSIN threw away the reason for every failure it reported.** OPSIN answers a
  name it cannot parse with HTTP 404 and a complete JSON body —
  `{"status": "FAILURE", "message": "notachemical12345 was uninterpretable due
  to the following section of the name: ..."}` — and the body is the only place
  that explanation exists. `OPSIN.get_id` mapped the status code through a
  lookup table and returned without reading it, so `message` was empty for
  every failure the module ever reported, including in the WARNING
  `get_id_from_list` logs, which therefore always read
  `Failed to get ID for x: ` with nothing after the colon.

  `provesid.opsin.opsin_classify` treats a 404 as a success for exactly this
  reason, and `get_id` reads `status` and `message` out of the body. The
  status-code table, `OPSIN.responses`, is gone — there is nothing left for it
  to do.

- **CAS Common Chemistry cached its own failures as answers.** `cas_to_detail`
  and `name_to_detail` report a failure in-band, by returning
  `{"status": "Timeout", "found": False, ...}`, and both were `@cached` with no
  `skip_if`. One timed-out request became a permanent "no such CAS number" on
  disk. This is the same defect fixed elsewhere in this release:
  `is_failure_result` looks for `success: False`, and CAS says `found: False`.
  Both methods now skip the cache unless `found` is True — absence included, for
  the reason `is_empty_result` gives.

- **Old OPSIN and CAS cache entries are retired.** Both clients cached failures
  before this release, so a name or CAS number looked up during one momentary
  outage sits on disk as a permanent failure for a record that exists — and
  every cached OPSIN failure carries an empty `message`, because the old code
  never read the body that explains it. Neither is fixable in place, so both
  clients gained a `CACHE_SCHEMA_VERSION = 2` inside a new `__cache_key__`,
  which makes version 1 entries unreachable. The CAS key deliberately excludes
  the API key: two keys reach the same registry and get the same answer.

- **A momentary OPSIN outage raised a bare `KeyError`.** `get_id` looked the
  status code up in a three-entry table, so a 503 — not one of the three —
  failed inside the lookup rather than being reported. It is retried now, and
  reported as a `"FAILURE"` whose `message` says what happened.

- `docs/api/pubchemview.md` documented `PubChemView(pause_time=...)`, a
  parameter that has never existed — both snippets raised `TypeError`.
  `docs/api/nci_resolver.md` described the default pacing as "3 requests per
  second"; it is 0.1 s between requests, so at most 10.

- **The persistent cache never hit across processes.** `CacheManager._get_cache_key`
  built its key with `json.dumps(..., default=str)` over the call's arguments.
  For a bound method the first argument is `self`, and the default `str()` of a
  client object embeds its memory address
  (`<PubChemView object at 0x7f...>`), so every instance — and every
  interpreter — produced a different key. No entry written by one run was ever
  reachable from the next: every "cached" method in `pubchem`, `pubchemview`,
  `resolver`, `cascommonchem`, `classyfire` and `opsin` re-fetched from the
  network on every call while the cache directory filled with unreachable
  duplicates.

  Arguments are now normalised by `cache.stable_key_part` into a form that is
  stable across processes: objects reduce to the value they declare through
  `__cache_key__()`, or to their fully qualified class name. `PubChemAPI`,
  `PubChemView` and `NCIChemicalIdentifierResolver` declare
  `__cache_key__` as `(class path, base_url)`, so two clients pointing at the
  same endpoint share cache entries while two different endpoints stay apart.

  **Existing cache entries are orphaned** by the new key scheme. They are inert,
  not incorrect; run `provesid.cache.clear_cache()` (or
  `clear_all_service_caches()`) to reclaim the disk space.

- **A cached `None` was indistinguishable from a cache miss.** `_load_from_disk`
  returned `None` both for "no entry" and for "the stored value is `None`", so
  any function whose result is `None` re-ran on every call. It now returns a
  `_MISS` sentinel.

- **Failed lookups were cached as answers.** `@cached` stored whatever a
  function returned, including the `{'success': False, 'error': ...}` dicts and
  the empty lists that the clients hand back after an HTTP 429, a PUG-View
  `ServerBusy` (503) or a timeout. A single transient error therefore became a
  permanent "this compound has no data". `@cached` now skips storage for any
  result that `cache.is_failure_result` recognises, and takes an optional
  `skip_if` predicate for clients that signal failure some other way.

- **`PubChemView` reported a failed fetch as an absent property.**
  `extract_property_data` caught `PubChemViewError` — the base class, which
  covers an exhausted retry budget — logged "not found", and returned `[]`.
  Combined with the caching bug above, a 503 during
  `get_melting_point(2244)` was stored as "aspirin has no melting point".
  Transport failures now propagate; only a genuine 404 yields `[]`.
  `get_property_table` makes the same distinction: an empty frame always means
  "no such data", never "the fetch failed".

- **`PubChemAPI.get_compound_synonyms` swallowed every error** into `[]`, with
  the same consequence. A 404 still returns `[]`; anything else raises.
  `get_compound_properties` keeps returning the properties it retrieved when the
  follow-up synonym request fails, but now records the failure under
  `synonyms_error` and, via `skip_if`, is not cached while incomplete.

- **PubChem sheds load behind 4xx statuses, so absence is now decided by the
  fault code, not the HTTP status.** Both services describe every error in the
  body — `{"Fault": {"Code": "PUGVIEW.NotFound", ...}}` for a compound that
  genuinely has no such data, `PUGVIEW.BadRequest` for an unknown heading,
  `...ServerBusy` when the service is merely busy. `pubchem.fault_code` and
  `pubchemview.fault_code` read it: a transient code is retried whatever status
  carries it, `NotFound`/`BadRequest` is absence and is not retried, and any
  other 4xx is a non-retryable error. Previously an unknown heading (HTTP 400)
  cost four requests to learn a permanent answer.

- **Absence is no longer persisted at all.** Even with the classification above,
  a wrongly-reported empty result would be permanent once cached, while
  re-fetching an empty one costs a single cheap request. Every cached extraction
  method in `pubchemview` and `PubChemAPI.get_compound_synonyms` now carries a
  `skip_if` predicate, so only positive results are stored. The raw response
  fetchers (`get_property`, `get_experimental_properties`) are unaffected: they
  raise on absence and so never had an empty to store.

- Removed an unreachable `return cache_info` left behind in
  `PubChemAPI.get_cache_info`.

### Changed
- **The documentation is rebuilt around the docstrings.** Every API page is
  now a short lead and a `:::` directive: `docs/api/` goes from 4 045
  hand-written lines, which had drifted from the code, to 207, and every
  public module has a page, including the
  offline half that had none (`PubChemID`, `CompToxID`, `ZeroPM`,
  `REACHDossierID`, `datasets`, `taxonomy`, `cache`, `config`, `tools`).
  What the docstrings do not say moved into eight guides under
  `docs/guide/`: installing the offline databases (new), `Search`, using the
  local databases directly, experimental properties from PubChem, network
  behaviour, caching, API keys and Chebifier. The home page and the quick
  start lead with `Search` and `datasets` rather than the online clients.
  - `docs/examples/` (symlinks to `examples/`) is deleted; a MkDocs hook,
    `scripts/mkdocs_hooks.py`, serves the tutorials straight from
    `examples/` at the same URLs. The `.py` demos are no longer rendered as
    notebook pages.
  - `docs/plans/` moved to `plans/`, and the "Modernization" section left
    the nav.
  - Docstrings are parsed with `docstring_style: auto`: `chembl`, `comptox`,
    `reach` and `zeropm` are NumPy-style and had been rendered as Google,
    with their example outputs read as Markdown links. A Google `Returns:`
    section is one value, where each wrapped line used to render as a
    separate returned value.
  - `mkdocs build --strict` now also fails on a link to a missing page or
    heading.
  - Docstrings cross-reference with Markdown, `` [`Search`][provesid.search.Search] ``,
    in place of the Sphinx roles (`:class:`, `:meth:` …) that the site showed
    literally; 591 of them are links now, and `--strict` fails on one whose
    target is gone. Constants documented with `#:` comments, which
    mkdocstrings does not read, carry attribute docstrings instead, so
    `DATASETS`, `Search.PRESETS`, `LOOKUPS`, `OUTPUT_COLUMNS` and 74 others
    now appear in the API reference. A typed return reads `(dict): …`, which
    renders the type in the Type column, not as a name; bulleted lists inside
    a section render as lists.
- **`Search` asks its sources through one table, `provesid.sources.LOOKUPS`.**
  The 47 hand-written `if client: try: ... except: log` blocks in nine
  methods of `search.py` became one row per (lookup kind, source) and a single
  driver, `Search._collect`. Adding a source is now one column in
  `sources.py` rather than an edit to every resolver. Results are unchanged:
  437 result rows over seven identifier types and seven configurations
  (default, ZeroPM, fuzzy, fuzzy with ZeroPM, InChIKey skeleton, Tanimoto,
  strict with salt stripping) were compared with the previous code on the real
  databases, and all of them are identical. The dead `_candidates_from_name`,
  `_fuzzy_name_candidates` and `_most_complete_row` were deleted, and
  `search.py` shrank from 3 392 to 2 721 lines. Only the wording of the
  warning logged when one source fails changed, to
  `"<Source> <kind> lookup failed for <query>: <error>"`.
- **Breaking: `PubChemID` and `ChebiSDF` have modules of their own.** The
  offline database clients were split off from the online clients they shared
  a file with: `PubChemID`, `rdkit_descriptors`, `RDKIT_DESCRIPTORS` and
  `PUBCHEM_DESCRIPTORS` moved from `provesid.pubchem` to `provesid.pubchem_id`,
  and `ChebiSDF` from `provesid.chebi` to `provesid.chebi_sdf`. The code
  moved unchanged. `from provesid import PubChemID, ChebiSDF` works as before,
  while `from provesid.pubchem import PubChemID` must become
  `from provesid.pubchem_id import PubChemID`. A side effect is that
  `import provesid.chebi` no longer imports RDKit.
- **Breaking: `PubChemID()` builds a missing database from PubChem's FTP site
  instead of downloading it from Zenodo.** That is 15.4 GB transferred and
  about 12 minutes of processing, against 2.2 GiB; pass `source="zenodo"` for the
  download. A database already on disk is opened as before, whichever way it
  was made, and `datasets.plan("pubchem")` now reports the FTP route's cost.
- **Breaking: `PubChemID.properties()` no longer serves the eight computed
  descriptors from disk.** `XLogP`, `TPSA`, `Complexity`, `Charge`,
  `HBondDonorCount`, `HBondAcceptorCount`, `RotatableBondCount` and
  `HeavyAtomCount` are PubChem's model outputs, not data about a compound; they
  now always come from PUG-REST, labelled `Source='online'`, even from a Zenodo
  copy that still stores them. Before, the same property came from a months-old
  snapshot for local compounds and from live PubChem for the rest, with no way
  to tell which. `MonoisotopicMass` moves the other way, and is served from
  disk by a database built from FTP. Which properties an open database can
  answer is `db.offline_properties`, and is also the default property list.
- **`PubChemID.get_by_cas_batch()` and `get_by_smiles_batch()` return the
  database's own columns** rather than a fixed list that named the descriptor
  columns, so a database built from FTP yields `monoisotopicmass` and no
  columns of `None`.
- **Breaking: `Search` no longer downloads missing datasets by default.** Pass
  `datasets="auto"` for the old behaviour. `redownload=True` now requires
  `datasets="auto"` and raises otherwise, rather than being silently ignored
  and handing back a stale copy.
- **`CheMBL.get_compound()` no longer returns `molfile`.** It is a quarter of a
  full ChEMBL database, nothing in PROVESID consumed it, and it is absent from
  the extract above. Build a MOL block from `canonical_smiles` with RDKit when
  you need one — RDKit is already a hard dependency.
- **`CheMBL` prefers a compacted extract.** With `chembl_37.db` and
  `chembl_37_provesid.db` side by side, `CheMBL()` opens the extract. Release
  number still wins first: an older extract does not beat a newer full release.

- **`OPSIN.base_url` is `https://www.ebi.ac.uk/opsin/ws/`.** The Cambridge
  address it used to name, `opsin.ch.cam.ac.uk`, answers every request with a
  301 to it (verified live, 2026-09-19), so the old URL worked and cost a
  redirect per name. `OPSIN.responses` is removed with it.

- **PUG-REST and PUG-View read a bare 400 differently, so they have separate
  classifiers.** `pubchem_classify` is now `pugrest_classify` and
  `pugview_classify`. PUG-View takes the heading as a query parameter, so a 400
  means "no such heading" — absence. PUG-REST takes its whole query in the URL
  path, so a 400 means the path was wrong, and the caller needs PubChem's own
  explanation of which property name it misspelled, which absence would discard.
  The fault code still decides first for both, so the split only governs a 400
  carrying no fault at all.

- **`PubChemAPI.pause_time` and `last_request_time` are properties** over the
  shared transport rather than plain attributes. `pause_time` is still settable
  mid-batch and still takes effect on the next request, retries included;
  `last_request_time` is read-only.

- **ChEBI paces itself at 10 requests per second and retries twice** with a
  half-second base back-off, where it previously did neither. Two retries rather
  than the transport's three: an ontology walk makes many small requests to a
  fast service, where a 1-2-4 second curve costs more than the request it is
  protecting. EBI publishes no per-IP figure for the ChEBI 2.0 API, so the
  pacing is politeness rather than a quoted limit.

  A ChEBI timeout message changed with the move: it now names the URL and the
  number of attempts.

- **CAS Common Chemistry paces itself at 5 requests per second and retries**,
  where it previously did neither. A 401 is never retried, because asking again
  with the same key cannot help.

- **An unresolvable identifier now costs the NCI resolver one request instead
  of four.** CACTUS answers an identifier it cannot resolve with **HTTP 500**
  carrying the body `<h1>Page not found (404)</h1>` — verified live on
  2026-09-19 against `this_is_definitely_not_a_chemical_12345`, while
  `α-glucose` answers 200. Its status code is not a reliable guide, so
  `provesid.resolver.nci_classify` reads the body: a 5xx carrying that page is
  absence and is raised at once, every other 5xx keeps its usual retryable
  reading. Two live tests that took 10.1 s and 9.6 s now take 0.57 s each.

  The exception type changes with it: an unresolvable identifier now raises
  `NCIResolverNotFoundError`, which is what it always meant, rather than
  `NCIResolverError("Internal server error")`. `NCIResolverNotFoundError` is a
  subclass, so `except NCIResolverError` is unaffected.

  Conversely, a *genuine* transient failure from CACTUS — a 429, a real 5xx, a
  timeout — is now retried, where before it was raised on the first attempt.

- **PubChem's fault-code classification has one home.** `fault_code` and
  `TRANSIENT_FAULT_CODES` were duplicated verbatim in `pubchem` and
  `pubchemview`, with `pubchemview` carrying the larger set. Both now live in
  `provesid.pubchem`, joined by `ABSENCE_FAULT_CODES` and by
  `pubchem_classify`, the classifier both services share.
  `provesid.pubchemview.fault_code` is gone; import it from `provesid.pubchem`.

- `NCIChemicalIdentifierResolver.pause_time` and
  `PubChemView.min_request_interval` are now properties over the transport's
  interval. Reading them is unchanged; setting one still takes effect on the
  next request, retries included.

- `NCIChemicalIdentifierResolver` logs to a module logger instead of the root
  logger, and `download_image` gains the retry and pacing the rest of the
  client already had.

- **`get_property_table`'s `ExperimentalValue` column is now a float**, not a
  string. An entry that reports a range leaves it NaN and fills `ValueMin` and
  `ValueMax` instead, rather than putting `"138-140"` in a column callers were
  expected to plot. `Unit` is normalised, and `Temperature` is a float in °C
  instead of a string like `"25°C"`.

- **`PubChemView._extract_experimental_value_and_unit` and `_parse_value_string`
  are gone**, replaced by `parse_value`. Between them they were 415 lines of
  property-specific regex cascade, and the two disagreed: the same string parsed
  one way through `get_melting_point` and another through `get_property_table`.

- Entries that state no value — PubChem's pointers to its own external tables,
  such as `{"Value": {"ExternalTableName": "iupacpka"}}` — are left out of
  property results instead of appearing as blank rows.

- The four module-level convenience functions in `pubchemview`
  (`get_experimental_property`, `get_all_experimental_properties`,
  `get_property_values_only`, `get_property_table`) are no longer cached
  themselves. Each delegates to a cached `PubChemView` method, so the second
  cache only kept a duplicate copy of the same payload on disk — under a key
  that did not carry `CACHE_SCHEMA_VERSION`, which is exactly how an upgrade
  would have served a stale shape.

- **`get_compound_properties_batch` now issues one request per 200 CIDs**
  instead of one per CID, having been rebuilt on `get_properties_for_cids`. The
  return shape is unchanged — one dict per CID with the properties plus
  `success`, `cid` and `error` — except that a CID PubChem has no record of is
  now reported with `success=False` and `"No such compound"` rather than
  whatever the single-CID path happened to return. The `output_format`
  parameter is gone: only JSON can be reshaped into per-CID dicts, so the
  parameter promised something it never delivered. Synonyms are not included,
  as they need a request per compound; the method never returned them.

- **`Search` no longer targets the ZeroPM database.** ZeroPM harvests regulatory
  inventories rather than curating compounds, so its name→structure rows are
  noisier than ChEBI/CompTox/PubChem/ChEMBL — while counting as a full
  independent vote in the corroboration ranking that 0.6.0 made drive
  `confidence` and `min_source_support`. The resolver now queries four sources
  by default, and ZeroPM's database is not even opened.

  The `ZeroPM` class is untouched and stays fully available for direct use; only
  `Search` stopped consulting it. Pass `use_zeropm=True` to restore the previous
  five-source behaviour. A `zeropm=` client handed to the constructor is ignored
  (with a warning) unless `use_zeropm=True` is set too, so "disabled" does not
  depend on how the caller happened to build the instance.

  Consequences: `source_details` no longer carries a `"ZeroPM"` entry,
  `sources_available` lists four keys, and confidence values shift slightly
  wherever ZeroPM used to vote. Fuzzy name queries lose recall — ZeroPM was the
  only source doing true fuzzy *retrieval*, so a typo sharing no substring with
  the real name (e.g. `"caffiene"`) now returns no match instead of a guess.
  Precision is unaffected: a misspelling still never resolves to a *different*
  compound, which `tests/test_search_precision_regression.py` now asserts
  explicitly for the default source set.

- `Search._SOURCE_KEYS` is now a per-instance attribute reflecting the sources
  that instance targets; the full catalogue lives in `Search._ALL_SOURCE_KEYS`
  and the default set in `Search._DEFAULT_SOURCE_KEYS`.

## [0.7.0] - 2026-08-17

### Fixed
- **`ChebifierClassifier.classify()` wrote model data into the caller's working
  directory.** The ensemble is *built* under `data_dir` because chemlog_extra and
  the smoother resolve paths relative to cwd, but `classify()` ran
  `predict_smiles_list` without that `chdir`, so predicting also created
  `data/chebi_v244/` and `data/chebi_v200/` wherever the process happened to be
  running. Now under `data_dir`, like every other path.
- **ChEMBL was unreachable, and `Search` lost a source without saying so.**
  `CheMBL.DEFAULT_DB_URL` pinned a release number *inside* EBI's moving
  `latest/` path (`latest/chembl_36_sqlite.tar.gz`), so the download started
  404ing the day ChEMBL 37 shipped — and would break again at 38. The release is
  now resolved from the `latest/` directory listing
  (`CheMBL.resolve_latest_db_url()`), with `DEFAULT_DB_URL` kept as a pinned
  fallback (`CheMBL.FALLBACK_RELEASE`) for when the listing cannot be read.
  `db_url=`, `db_name=` and `db_path=` still override everything.

  Because `Search` catches source-initialisation failures and continues, and
  0.6.0 made corroboration drive confidence, every run since ChEMBL 37 appeared
  scored lower than a full-source run and satisfied `min_source_support` less
  often, with nothing in the result to show why.

### Changed
- **`CheMBL()` no longer pins a database filename.** `db_name` defaults to
  `None`: the newest `chembl_*.db` already in the data directory is reused (so a
  new ChEMBL release does not trigger a multi-gigabyte re-download — pass
  `redownload=True` for that), and otherwise the name is derived from the
  resolved archive (`chembl_37_sqlite.tar.gz` → `chembl_37.db`). New
  `CheMBL.release` attribute reports the release number in use.
- **Simplified the chebifier install to the two commands upstream now supports.**
  chebifier 1.2.2 ships a `models` extra that pins the whole model stack
  (`chebai` 1.2.0, `chebai-graph` 1.0.0, `chemlog-extra` 1.0.1, `c3p` 0.5.0), all
  on PyPI, so `scripts/install_chebifier.sh` is now a thin wrapper around:

  ```bash
  uv pip install "chebifier[models]"
  uv pip install torch==2.12.0 torch_scatter torch_geometric \
      -f https://data.pyg.org/whl/torch-2.12.0+cpu.html
  ```

  Consequences, all verified by running the full ensemble on CPU (benzene and
  aspirin → sensible ChEBI classes, every model incl. the GNNs loaded):
  - The `provesid[chebifier]` extra is now `chebifier[models]==1.2.2` (was bare
    `chebifier==1.2.1`), so `pip install 'provesid[chebifier]'` alone gets the
    transformer and rule-based models working.
  - **No more git installs** for `chemlog-extra` and `c3p`, and no explicit
    `chebi-utils` — `chebai-graph` 1.0.0 does not need it.
  - **torch is no longer capped at 2.11.** Only `torch_scatter` is required (not
    `torch_sparse`/`torch_cluster`/`pyg_lib`); `torch_cluster`, which has no wheel
    past torch 2.11, was the sole reason for the old pin. torch **2.12.0** now.
  - **No index patching needed.** `chebai-graph` 1.0.0 predates the property-index
    drift that broke the `v244` GNN checkpoints, so the installer no longer
    rewrites files inside site-packages. `ensure_v244_indices()` stays as a
    runtime safety net (reports `ok` on a clean install) for a hand-upgraded
    `chebai-graph`.
  - The installer takes torch from the PyTorch CPU index by default
    (`TORCH_INDEX_URL=""` to opt out): **1.6 GB** of site-packages instead of
    **5.4 GB**, since plain PyPI torch adds 2.7 GB of CUDA wheels plus triton. It
    also verifies with the interpreter it installed into, rather than whatever
    `python3` resolves to (this failed when `VIRTUAL_ENV` was set but the
    environment was not activated).
- `CHEBIFIER_PINNED_VERSION` is `"1.2.2"`.

### Added
- **`ChebifierClassifier(with_scores=True)` and `predict_with_scores()`** —
  per-label confidence from the chebifier ensemble. `predict_smiles_list`
  computes a smoothed net score, thresholds it at 0 to pick the surviving
  classes, then discards it; `predict_with_scores` runs the same four steps
  (`gather_predictions` → `consolidate_predictions` → smoother → `> 0`) and keeps
  it, so `classify()`'s `confidence` column is populated for no extra cost — the
  models still run exactly once. Verified to reproduce `predict_smiles_list`'s
  label sets exactly (120 molecules, 0 mismatches); scores land in `(0, 1]`.
- **`ChebifierClassifier(exclude_models=[...])`** — build the ensemble without
  named models. The motivating case: chebai's tokenizer returns `None` for SMILES
  outside its vocabulary and the collator then dies on `len(None)`, taking the
  **whole batch** down. Only `electra` does this, and a reduced ensemble
  classifies those structures fine — which turned 414 hard failures into 0 over a
  100k-compound run. Both new options participate in the cache key, so a
  score-less or reduced-ensemble prediction is never served to a caller who asked
  for something else.
- **`chebi_class_names()`** — ChEBI id → name for all ~204k terms, read from the
  ontology snapshot the ensemble already loads. Offline, and the practical way to
  label thousands of predicted classes; the alternative was one HTTP call per id.
  Ids are bare, matching what the ensemble predicts.
- **`Search.sources_available` / `Search.sources_unavailable`**, mirrored on the
  result frame as `df.attrs["sources_available"]` / `["sources_unavailable"]`,
  plus a single warning naming the missing sources. A degraded run is now
  identifiable after the fact instead of looking like a full one.
- Tests: ChEMBL release resolution against a mocked listing (including the
  fallback and the "existing database is reused without network access" case),
  a live check that the resolved archive URL is downloadable, and a smoke test
  asserting all five offline `Search` sources initialise.

## [0.6.0] - 2026-08-04

### Fixed
- **ChEBI record lookups returned a *neighbouring* compound.** `ChebiSDF` built
  its index with file offsets counted in text mode, where universal-newline
  translation collapses `\r\n` to `\n`. The ChEBI SDF mixes line endings
  (~59 000 CRLF lines in the 2026 release), so every one of them under-counted a
  byte and the drift grew to ~59 kB by the end of the file. `get_compound_by_id`
  then seeked into an adjacent record and returned it — silently, since a
  neighbouring record parses perfectly well. Offsets are now computed and
  consumed in binary. The ChEBI *data* was never wrong; only the offsets were.

  This is what made `Search("cas")` return the wrong structure for 18 of 65
  pesticide CAS numbers (ChEBI answered "Mefluidide" for metaldehyde's
  108-62-3). All 18 now resolve correctly.
- **A cached ChEBI index is now validated against the SDF** (format version +
  file size) and rebuilt when it does not match, instead of being trusted
  blindly. Without this, an index written by an affected release keeps returning
  wrong compounds after the code is fixed.
- **Corroboration now counts in `Search`'s confidence.** The score was driven by
  *which* database answered rather than *how many* agreed: a lone ChEBI hit
  scored a flat 0.90 while a structure CompTox, PubChem and ZeroPM all carried
  scored 0.8777, so an uncorroborated hit outranked a three-source consensus.
  `consensus_score` cannot express this on its own — a single source agrees with
  itself perfectly — so confidence is now multiplied by a corroboration factor
  (1 source ×0.85, 2 ×0.95, 3+ ×1.0).
- **Group records are no longer returned as compounds.** A SMILES with an
  attachment point (`*C(=O)CCCC=CCC=CCCCCC`, a ChEBI *group*) is a substituent,
  never the substance a CAS or name denotes; such candidates are dropped unless
  the query is itself a group SMILES. RDKit could not process them either
  (`Unsupported in this mode element '*'`).

### Added
- **`Search(min_source_support=...)`** (also a per-call `search()` override) —
  require a structure to be carried by at least this many independent databases
  before it is returned, trading recall for precision.
- `tests/test_search_scoring_truth.py` — scoring-system regression suite: unit
  tests for the confidence rules, ChEBI record round-trip tests that would have
  caught the offset drift, and CompTox-truth samples for every identifier type
  (CAS, name, InChIKey, DTXSID, SMILES, InChI) plus `enrich()` and
  `resolve_cascade()`. Measured on 1500 sampled CompTox CASRNs: 1498 correct,
  2 disagreements (both genuine cross-database structure differences, not
  ranking defects).
- `examples/search/confidence_and_corroboration_demo.py`.

## [0.5.0] - 2026-08-02

### Added
- **`Search.enrich(df, column)`** — attach resolved identifier columns to a
  DataFrame, searching each *distinct* value once and broadcasting the result to
  every row that carries it. Row order and index are preserved; added columns are
  namespaced with a configurable `prefix` (default `provesid_`).
- **`resolve_cascade(df, stages, ...)`** — resolve rows through an ordered list of
  `Search` stages, passing each stage only the rows still unresolved, so every row
  is resolved by the most reliable identifier it actually has. Records
  `resolved_by` and `validated_by` per row, with an optional RDKit fallback that
  derives identifiers from the row's own structure.
- **`mw_within(tolerance, reference_column=...)`** — validator factory for
  `resolve_cascade`'s `accept` argument: accepts a hit only when its molecular
  weight agrees with the structure the dataset already carried, and reports any
  additional SMILES/name agreement. This is what stops a confident-but-wrong
  identifier match from ending a cascade.
- `CheMBL.search_by_name(..., exact=True)` for exact (case-insensitive) name and
  synonym matching.
- `ZeroPM.match_similar_name()` and `ZeroPM.get_id_table_from_similar_name()` —
  fuzzy name matching that reports *what* matched and *how well*, instead of
  discarding it.
- `ZeroPM.zeropm_id_to_inchi_id()` — the reverse of `get_zeropm_id`.
- `taxonomy.ensure_element_class_mappings()`, `default_ensemble_available()` and
  `missing_ensemble_modules()`. The last two report on the *whole* default
  ensemble (transformer, graph, rule-based and c3p models each live in a separate
  package), so a partial install can be detected up front instead of failing with
  a bare `ModuleNotFoundError` at predict time.

### Fixed
- **`Search` resolved misspelled names to unrelated compounds and labelled them
  exact matches.** `Search("name", fuzzy=True).search("asprin")` returned
  PHENYRAMIDOL with `match_method="exact_name"`. Three defects combined:
  - `CheMBL.search_by_name` had no exact mode, so `Search`'s exact pass received a
    substring match — PHENYRAMIDOL carries the synonym `"Evasprin"`, which
    contains `"asprin"` — and tagged it `exact_name`.
  - The "did the exact pass find a strong match?" test used a fuzzy score, and
    `WRatio("asprin", "Evasprin")` is 85.7, clearing the cut-off. That suppressed
    the fuzzy widening which would have found aspirin. The test now requires an
    actual (case- and whitespace-insensitive) name equality.
  - A fuzzy match's confidence base was the raw similarity (up to 1.0) while an
    exact name match was pinned at 0.80, so **a typo could score higher than the
    correct spelling**. The fuzzy base is now scaled by the exact-name base, so an
    approximate match can never outrank an exact one.
- `Search`'s ZeroPM fuzzy branch was dead code: it acted only on a `DataFrame`,
  but `query_similar_name` returns a list of ids, so ZeroPM — the only source
  doing true fuzzy *retrieval* — never contributed to fuzzy name search.
- ChEMBL was queried in `Search`'s exact name pass but omitted from the fuzzy
  widening pass; it now participates in both.
- **`ZeroPM.get_pm_probabilities`, `batch_get_pm_probabilities` and
  `get_all_zeropm_chemicals` could never return P/M probability data.** All three
  keyed `pm_probabilities` on `zeropm_id`, but that table is keyed on `inchi_id`,
  so every call raised `sqlite3.OperationalError: no such column: zeropm_id`.
  `get_pm_probabilities` now translates a `zeropm_id` via the new
  `zeropm_id_to_inchi_id`, and the two joins use `inchi_id`.
- **`ChebifierClassifier` could not build the default ensemble at all.**
  `chemlog_extra` reads its element-class mapping files from a *working-directory
  relative* path and rebuilds them from the ChEBI graph when absent — but that
  rebuild crashes, because 288 of the graph's 205k nodes carry `name: None` and
  the builder does `" molecular entity" in properties["name"]`. Every
  `classify()` call failed with `TypeError: argument of type 'NoneType' is not
  iterable`. The new `ensure_element_class_mappings` writes both files into the
  PROVESID chebifier data directory using upstream's own derivation rules
  (skipping unnamed nodes), and the ensemble is now constructed with that
  directory as the working directory.

### Fixed — tests
- `test_pubchem_id.py::test_init_nonexistent_path` omitted `auto_download=False`,
  so instead of asserting `FileNotFoundError` it **downloaded the ~2.3 GB
  database into the repository root** on every run. It now passes the flag and
  uses `tmp_path`.
- `test_zeropm.py::test_get_cas_from_name_integration` asserted a
  name → CAS → same-CAS round trip that the data model does not support (names
  are many-to-many with CAS numbers), and its fallback clause compared a CAS
  against a list of *names*, so it could never hold. It passed only because its
  `SELECT ... LIMIT 1` had no `ORDER BY`: creating indexes — which another test in
  the file does to the shared database — changed which row came back and broke
  it. Replaced with a deterministic test of the guarantee that does hold.
- The three P/M probability tests repeated the same wrong join key in their own
  SQL, and their `else: pytest.skip(...)` branches turned the resulting error
  into a silent skip — which is how the production bug survived. They now assert
  that the fixture query found data, and check the returned values against the
  database so a wrong join cannot pass again.
- `test_zeropm.py` asserted `dtype == object` for string columns; pandas 3 gives
  `StringDtype`. Now uses `pandas.api.types.is_string_dtype`.
- `test_search.py::test_exact_inchikey_with_low_consensus` expected 0.5 for a
  zero consensus score. Zero consensus only occurs when no source matched at all
  (one source scores 1.0; two fully disagreeing sources score 0.5), so 0.0 is
  correct and the test encoded the un-special-cased formula. Replaced with tests
  for both the zero and the partial-agreement cases.

### Removed — tests
- The end-to-end chebifier classification tests (`TestLiveClassification`).
  chebifier stays an optional extra, and its full model stack (transformer,
  graph/GNN, rule-based and c3p models — separate packages, some git-only) is
  awkward enough to install that the test suite should not depend on it. The
  remaining `test_taxonomy.py` tests all pass with the extra absent, verified by
  running them with the stack made unimportable. The now-unused `chebifier`
  pytest marker was dropped from `pyproject.toml`.

### Changed
- **`Search`'s default `fuzzy_scorer` is now `"ratio"`, was `"WRatio"`.** `WRatio`
  includes a partial-ratio term that scores a short name highly whenever it
  appears anywhere inside the query, which makes `fuzzy_score_cutoff` stop
  discriminating: `WRatio("caffiene", "ne")` and
  `WRatio("zzzznotachemical", "Mica")` are both 90, while `ratio` puts both at 40
  and still scores the genuine typo `ratio("caffiene", "caffeine")` at 87.5. Pass
  `fuzzy_scorer="WRatio"` to restore the old behaviour.
  Confidence values for `fuzzy_name` matches change as a result.
- `ZeroPM` "not found in database" messages moved from `WARNING` to `DEBUG` and
  onto the instance logger. They fired on the root logger during successful
  `Search` runs.

### Removed
- GitHub Actions workflows for running tests (`test.yml`,
  `test-with-api-keys.yml`) and for releasing (`release.yml`). Tests are run
  locally; releases are made manually with `twine`. Only the documentation
  deploy workflow remains.

## [0.3.0] - 2026-04-16

### Added
- **`Search` class** — unified cross-database chemical identifier resolver (`provesid.Search`)
  - Accepts CAS, name, SMILES, InChI, InChIKey, DTXSID, or molecular formula as input
  - Queries all five offline databases (ChEBI, CompTox, PubChemID, ZeroPM, ChEMBL) in a single call
  - Returns a `pandas.DataFrame` with 23 standardised columns including `confidence`, `match_method`, `source_details`, and `source_match_scores`
  - Optional **salt stripping** (`strip_salts=True`): uses RDKit `SaltRemover` + largest-fragment picker; parent SMILES and InChIKey stored in dedicated columns
  - Optional **fuzzy name matching** (`fuzzy=True`): RapidFuzz-based candidate search with configurable similarity threshold
  - Optional **Tanimoto similarity** search (`similarity_threshold>0`): Morgan fingerprint fallback when no exact match is found
  - Optional **InChIKey skeleton** search (`inchikey_skeleton=True`): matches stereoisomers via 14-character skeleton prefix
  - Accepts `str`, `list[str]`, `pd.DataFrame` (with `column=` kwarg), or a CSV/Parquet file path as input
  - Confidence scoring model: base per match method × (0.5 + 0.5 × cross-source consensus score)
  - Full per-source traceability in `source_details` column
- **`normalize_structure(smiles)`** — RDKit helper returning canonical SMILES, Kekulized SMILES, InChI, InChIKey, and molecular weight
- **`strip_salts(smiles)`** — standalone salt-stripping utility exported from `provesid`
- **`OUTPUT_COLUMNS`** — list of all 23 column names in the `Search` result schema, exported from `provesid`
- Example scripts: `examples/search/search_by_cas_demo.py`, `search_by_name_demo.py`, `salt_stripping_demo.py`, `similarity_search_demo.py`
- API documentation page `docs/api/search.md`

### Fixed
- `strip_salts`: when `SaltRemover` strips all fragments, fall back to the largest fragment of the original molecule instead of returning an empty string

## [0.2.0] - 2025-09-29

### Added
- **🚀 Unlimited Caching System**: Complete overhaul of caching infrastructure
  - Unlimited cache by default (no size limits)
  - Persistent cache storage across sessions
  - 5GB warning threshold with configurable monitoring
  - Memory + disk hybrid caching for optimal performance
  - SHA256 cache keys for security and uniqueness
  - Import/Export functionality for team collaboration (pickle and JSON formats)
  - Global cache management functions: `clear_cache()`, `get_cache_info()`, `export_cache()`, `import_cache()`

- **📊 Comprehensive API Caching**: All major APIs now support unlimited caching
  - **PubChemAPI**: 19 cached methods including `get_compounds()`, `get_properties()`, `get_synonyms()`, etc.
  - **CASCommonChem**: 2 cached methods (`cas_to_detail()`, `name_to_detail()`)
  - **NCIChemicalIdentifierResolver**: 15 cached methods including all convenience functions
  - **PubChemView**: 15+ cached methods for experimental property extraction
  - **ClassyFireAPI**: 3 cached methods (`submit_query()`, `query_status()`, `get_query()`)
  - **OPSIN**: 2 cached methods (`get_id()`, `get_id_from_list()`)

- **🔧 Cache Management Methods**: Each API class now includes:
  - `clear_cache()`: Clear cached data for that specific API
  - `get_cache_info()`: Get detailed cache statistics and information

- **📈 Performance Improvements**:
  - Significant speed improvements for repeated API calls
  - Reduced API rate limiting issues
  - Offline capability when APIs are unavailable
  - Cross-session data persistence

### Changed
- **Breaking**: Removed `cache_size` parameter from PubChemAPI constructor (now unlimited by default)
- **Breaking**: Replaced all `@lru_cache(maxsize=X)` decorators with unlimited `@cached` decorator
- Cache behavior is now consistent across all APIs
- Cache storage moved from memory-only to persistent disk storage

### Enhanced
- **Test Coverage**: Added comprehensive caching tests
  - 168 tests passing with new caching system
  - Cache persistence, import/export, and size monitoring tests
- **Documentation**: Updated API documentation to reflect caching capabilities
- **Error Handling**: Improved cache error handling and recovery

### Technical Details
- New `cache.py` module with `CacheManager` class
- Automatic cache directory creation in system temp folder
- Thread-safe cache operations
- Configurable warning thresholds and cache policies
- Full backward compatibility (existing code continues to work)

### Performance Metrics
- Cache hit rates: Near 100% for repeated identical requests
- Memory usage: Efficient hybrid memory/disk storage
- Disk usage: Automatic monitoring with configurable warnings
- Speed improvement: 10-100x faster for cached requests

## [0.1.0] - Initial Release

### Added
- Initial implementation of PROVESID package
- PubChemAPI for PubChem REST API access
- CASCommonChem for CAS Common Chemistry API
- NCIChemicalIdentifierResolver for NCI resolver
- PubChemView for experimental properties
- ClassyFireAPI for chemical classification
- OPSIN for IUPAC name to structure conversion
- ChEBI API integration
- Basic caching with lru_cache (limited size)
- Comprehensive test suite
- Documentation and examples

[0.2.0]: https://github.com/USEtox/PROVESID/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/USEtox/PROVESID/releases/tag/v0.1.0