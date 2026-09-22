"""
The bulk datasets PROVESID reads: what they are, and how they get here.

Two halves. [`download_file`][provesid.datasets.download_file] is the transport
--- one resumable, checksummed downloader shared by every dataset in the
package. The registry below is the manager:
[`status`][provesid.datasets.status] says what is on disk,
[`plan`][provesid.datasets.plan] what a download would cost,
[`fetch`][provesid.datasets.fetch] installs a dataset by name and
[`remove`][provesid.datasets.remove] reclaims its space. The second half exists
because the first one worked too well: a clean machine running one CAS lookup
through [`Search`][provesid.search.Search] used to spend ~32 GB without asking
anyone.

[`provesid.http`][provesid.http] is the transport for web *APIs* --- small
requests, paced at five a second, retried a few times.  The bulk datasets are a
different problem entirely: ChEMBL's archive is 5.8 GB, PubChem's identifier
database 2.2 GB, and what those transfers want is not a pacer but *resumption*,
a checksum and an atomic rename.  So they were deliberately left out of the
``http.py`` migration --- which left the largest transfers in the package as
the only ones with no shared implementation.

Five modules had grown their own copy of "stream the response into a temporary
file with a progress bar", and all five shared the same defects:

======================  =====  ======  ========
Call site               Retry  Resume  Checksum
======================  =====  ======  ========
``pubchem.py``          no     no      no
``comptox.py``          no     no      no
``zeropm.py``           no     no      no
``chembl.py``           no     no      no
``chebi.py``            no     no      no
======================  =====  ======  ========

An interrupted 5.8 GB ChEMBL download started again from zero.  Corruption was
caught unevenly: the two gzipped downloads got truncation detection for free
from gzip's CRC trailer, while the three plain SQLite files had none at all ---
a file truncated in its interior opens cleanly and fails much later, on the
first query that touches a missing page, which a user reads as a data problem
rather than as a download problem.

[`download_file`][provesid.datasets.download_file] is the one implementation.
It resumes against a ``.part`` file with an HTTP ``Range`` request, verifies
an MD5 when the server publishes one, checks the byte count it actually
received against the one the server declared, hands the finished file to the
caller's own validity check, and only then moves it into place.  Nothing is
ever renamed onto the destination until every one of those has passed, so a
failed download can never replace a good database with a broken one.

Examples:
    >>> from provesid.datasets import download_file
    >>> download_file(                                    # doctest: +SKIP
    ...     "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz",
    ...     "/data/CID-SMILES.gz",
    ...     checksum_url="https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz.md5",
    ... )
    '/data/CID-SMILES.gz'

    >>> from provesid import datasets
    >>> datasets.status()                                 # doctest: +SKIP
    >>> datasets.plan(["pubchem", "chebi"])               # doctest: +SKIP
    >>> datasets.fetch("pubchem")                         # doctest: +SKIP
    >>> datasets.remove("chembl")                         # doctest: +SKIP
"""

import glob
import hashlib
import importlib
import logging
import os
import re
import time
from dataclasses import dataclass
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple, Union
from urllib.parse import urlsplit

import pandas as pd
import requests
from tqdm import tqdm

from .http import RETRYABLE_STATUS, ServiceError, retry_after_seconds
from .utils import user_dataset_path

logger = logging.getLogger(__name__)

CHUNK_SIZE = 1024 * 1024
"""Bytes read from the socket, and from disk while checksumming, at a time.
A megabyte is large enough that the per-chunk work disappears against the
transfer and small enough that the progress bar still moves.
"""

PART_SUFFIX = ".part"
"""Suffix of the partial file a download accumulates into.  It sits beside the
destination rather than in a temporary directory so that a resumed download
finds it, and so that the final move is a rename within one filesystem.
"""

SOURCE_SUFFIX = PART_SUFFIX + ".source"
"""Suffix of the marker recording which URL a partial file came from.  Without
it, a ``.part`` left over from a different release --- or a different
dataset that happens to share a destination --- would be resumed, splicing
two files together.  A checksum would catch that, but only the PubChem FTP
mirror publishes one; the Zenodo downloads have no such backstop.
"""


class DownloadError(ServiceError):
    """
    A bulk dataset could not be downloaded, or arrived damaged.

    A [`ServiceError`][provesid.http.ServiceError] like any other failure this package
    reports from a remote service, so a caller can catch the whole family; the
    separate class exists because the recovery differs. An API call that fails
    is retried or abandoned, while a download that fails has usually left a
    resumable ``.part`` file on disk and will pick up where it stopped on the
    next call.

    The exception is a checksum mismatch, which deletes the partial file: the
    bytes on disk are known to be wrong, so resuming from them would only
    produce the same wrong file again.

    Examples:
        >>> try:                                                 # doctest: +SKIP
        ...     download_file("https://zenodo.org/records/0/files/missing.db", "/tmp/x.db")
        ... except DownloadError as exc:
        ...     print("will resume from", exc)
        >>> issubclass(DownloadError, ServiceError)
        True
    """


def md5_of_file(path: str, chunk_size: int = CHUNK_SIZE) -> str:
    """
    MD5 digest of a file, read in chunks.

    Args:
        path: File to digest.
        chunk_size: Bytes to read at a time. The default keeps a multi-gigabyte
            file's digest off the heap.

    Returns:
        The digest as lower-case hexadecimal.

    Examples:
        >>> import tempfile, os
        >>> handle, path = tempfile.mkstemp()
        >>> _ = os.write(handle, b"provesid"); os.close(handle)
        >>> md5_of_file(path)
        '302c7456c0426cb916dfd4920a9c25bc'
        >>> os.remove(path)
    """
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_checksum(url: str, *, session: Optional[requests.Session] = None,
                  timeout: float = 30) -> str:
    """
    Fetch a checksum published beside a file.

    PubChem's FTP mirror publishes an ``.md5`` next to every file, in the
    format ``coreutils`` writes: the digest, whitespace, then the filename.
    Some servers publish the digest alone. Both are accepted.

    Args:
        url: URL of the checksum file.
        session: ``requests.Session`` to fetch through, for connection reuse.
        timeout: Seconds to wait for the response.

    Returns:
        The digest as lower-case hexadecimal.

    Raises:
        DownloadError: If the checksum file cannot be fetched or does not look
            like a digest.

    Examples:
        >>> read_checksum(                                          # doctest: +SKIP
        ...     "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz.md5")
        '3659dd5c96fc506bb11b7fd8cce4553d'
    """
    getter = session.get if session is not None else requests.get
    try:
        response = getter(url, timeout=timeout)
        response.raise_for_status()
    except requests.RequestException as exc:
        raise DownloadError(f"Could not fetch checksum from {url}: {exc}",
                            url=url) from exc

    digest = response.text.strip().split()[0].lower() if response.text.strip() else ""
    if len(digest) != 32 or not all(c in "0123456789abcdef" for c in digest):
        raise DownloadError(
            f"{url} does not contain an MD5 digest: {response.text[:80]!r}",
            url=url,
        )
    return digest


def download_file(
    url: str,
    dest: str,
    *,
    expected_md5: Optional[str] = None,
    checksum_url: Optional[str] = None,
    verify: Optional[Callable[[str], None]] = None,
    resume: bool = True,
    max_retries: int = 4,
    backoff: float = 2.0,
    max_backoff: float = 60.0,
    timeout: float = 60.0,
    chunk_size: int = CHUNK_SIZE,
    progress: bool = True,
    description: Optional[str] = None,
    session: Optional[requests.Session] = None,
    log: Optional[logging.Logger] = None,
) -> str:
    """
    Download a large file, resuming and verifying it, and move it into place.

    The download accumulates into ``dest + '.part'``. If that file is already
    there from an interrupted attempt, the request carries a ``Range`` header
    and the transfer continues from where it stopped; a server that ignores the
    header and answers 200 simply starts the file again, which is correct
    though slower. Nothing is written to ``dest`` until the transfer has
    completed and every check below has passed, so an interrupted or damaged
    download can never replace a good file.

    A partial file is resumed only when it came from the same ``url``. The URL
    is recorded in a ``.part.source`` marker beside it, and a partial left by
    some other download is discarded instead of being spliced onto this one.

    Four checks, in order, each of which leaves ``dest`` untouched if it fails:

    1. **Byte count.** When the server declares a size, a short file is a
       truncated transfer and is retried rather than accepted. This is the
       check the three SQLite downloads never had --- a file truncated in its
       interior opens cleanly and fails much later, on the first query to touch
       a missing page.
    2. **Checksum**, when one is available from ``expected_md5`` or
       ``checksum_url``.
    3. **The caller's own** ``verify``, which is handed the finished file and
       raises if it is not usable --- typically opening it and querying a table
       it must contain.
    4. **Atomic rename** onto ``dest``.

    Args:
        url: What to download.
        dest: Where it ends up. Parent directories are created. Overwritten
            only on success.
        expected_md5: Digest the finished file must have. Takes precedence over
            ``checksum_url``.
        checksum_url: URL of a published checksum to fetch and use, such as the
            ``.md5`` PubChem writes beside every FTP file. Ignored when
            ``expected_md5`` is given.
        verify: Called with the path of the finished file before it is moved
            into place; raise from it to reject the download. Any exception
            propagates unchanged, so a caller keeps its own error type. A
            rejected file is deleted rather than kept for resumption: it is
            already complete, so resuming it would fail the same check again.
        resume: Whether to continue an existing ``.part`` file. False discards
            it and starts over. A partial file is resumed only when it came
            from this same ``url``, which is recorded in a marker beside it;
            one left by a different download is discarded rather than spliced
            onto the new one.
        max_retries: Retries *after* the first attempt, so the file is fetched
            at most ``max_retries + 1`` times. Each retry resumes from what is
            already on disk, so a flaky connection makes progress rather than
            restarting.
        backoff: Base for the exponential wait, ``backoff * 2 ** attempt``
            seconds. A ``Retry-After`` header overrides it when larger.
        max_backoff: Ceiling on any single wait, including one ``Retry-After``
            asks for.
        timeout: Seconds to wait for the response to begin. Not a ceiling on
            the transfer, which may legitimately run for an hour.
        chunk_size: Bytes read from the socket at a time.
        progress: Whether to show a progress bar. It starts at the resumed
            offset, so a resumed download reports its true position.
        description: Label for the progress bar. Defaults to the filename.
        session: ``requests.Session`` to download through.
        log: Logger for the progress and retry messages. Defaults to this
            module's.

    Returns:
        ``dest``.

    Raises:
        DownloadError: If the transfer could not be completed within the retry
            budget, if the checksum does not match, or if the server answered
            with a status that is not worth retrying.

    Examples:
        >>> import sqlite3
        >>> def must_be_a_database(path):
        ...     sqlite3.connect(path).execute("SELECT 1 FROM compounds LIMIT 1")
        >>> download_file(                                       # doctest: +SKIP
        ...     "https://zenodo.org/records/1234/files/pubchem_id.db",
        ...     "/data/pubchem_id.db",
        ...     verify=must_be_a_database,
        ... )
        '/data/pubchem_id.db'
    """
    log = log or logger
    dest = os.path.abspath(os.path.expanduser(dest))
    part = dest + PART_SUFFIX
    parent = os.path.dirname(dest)
    if parent:
        os.makedirs(parent, exist_ok=True)

    if expected_md5 is None and checksum_url is not None:
        expected_md5 = read_checksum(checksum_url, session=session, timeout=timeout)
        log.debug("Expecting MD5 %s from %s", expected_md5, checksum_url)

    label = description or os.path.basename(dest)
    marker = dest + SOURCE_SUFFIX
    if not resume:
        _discard_partial(part, marker)
    elif os.path.exists(part) and _recorded_source(marker) != url:
        log.info("Discarding a partial %s left by a different download", label)
        _discard_partial(part, marker)

    log.info("Downloading %s from %s", label, url)

    with open(marker, "w") as handle:
        handle.write(url)

    for attempt in range(max_retries + 1):
        already = os.path.getsize(part) if os.path.exists(part) else 0
        if already and attempt == 0:
            log.info("Resuming %s from %.1f MB already on disk", label, already / 1e6)
        try:
            _stream_into(url, part, already, label=label, timeout=timeout,
                         chunk_size=chunk_size, progress=progress,
                         session=session, log=log)
            break
        except _Retryable as exc:
            if exc.restart and os.path.exists(part):
                # The server cannot serve the range we asked for, so what is on
                # disk is not a prefix of the file it is offering.
                os.remove(part)
            if attempt == max_retries:
                raise DownloadError(
                    f"Download of {url} failed after {max_retries + 1} attempts: "
                    f"{exc}. {_resume_hint(part)}",
                    url=url,
                    status_code=exc.status_code,
                ) from exc
            wait = min(max(backoff * 2 ** attempt, exc.retry_after or 0), max_backoff)
            log.warning("Download of %s interrupted (%s); retrying in %.1fs "
                        "(attempt %d of %d)", label, exc, wait, attempt + 2,
                        max_retries + 1)
            time.sleep(wait)

    if expected_md5 is not None:
        log.info("Verifying checksum of %s", label)
        actual = md5_of_file(part, chunk_size)
        if actual != expected_md5:
            _discard_partial(part, marker)
            raise DownloadError(
                f"Checksum mismatch for {url}: expected MD5 {expected_md5}, "
                f"got {actual}. The partial file has been deleted; run the "
                f"download again to start fresh.",
                url=url,
            )

    if verify is not None:
        log.info("Checking that %s is usable", label)
        try:
            verify(part)
        except Exception:
            # The transfer completed; the bytes are simply not what they should
            # be. Resuming from them would complete instantly and fail the same
            # check again, so the partial file goes -- as it does for a
            # checksum mismatch, and for the same reason.
            _discard_partial(part, marker)
            raise

    os.replace(part, dest)
    if os.path.exists(marker):
        os.remove(marker)
    log.info("%s ready at %s (%.2f GB)", label, dest, os.path.getsize(dest) / 1e9)
    return dest


def _recorded_source(marker: str) -> Optional[str]:
    """The URL a partial file came from, or None if it was not recorded."""
    try:
        with open(marker) as handle:
            return handle.read().strip()
    except OSError:
        return None


def _discard_partial(part: str, marker: str) -> None:
    """Delete a partial file and the marker naming where it came from."""
    for path in (part, marker):
        if os.path.exists(path):
            os.remove(path)


class _Retryable(Exception):
    """
    A download attempt failed in a way that another attempt might survive.

    Carries what the retry loop needs to decide *how* to try again:
    ``retry_after`` when the server named a delay, and ``restart`` when what is
    on disk must be thrown away first rather than resumed.
    """

    def __init__(self, message: str, *, status_code: Optional[int] = None,
                 retry_after: Optional[float] = None, restart: bool = False):
        super().__init__(message)
        self.status_code = status_code
        self.retry_after = retry_after
        self.restart = restart


def _stream_into(url: str, part: str, already: int, *, label: str, timeout: float,
                 chunk_size: int, progress: bool, session: Optional[requests.Session],
                 log: logging.Logger) -> None:
    """
    Append one attempt's worth of ``url`` to ``part``, starting at ``already``.

    Raises:
        _Retryable: On a connection failure, a timeout, a retryable status, or
            a transfer that ended short of the size the server declared.
        DownloadError: On a status that another attempt would not change.
    """
    getter = session.get if session is not None else requests.get
    headers = {"Range": f"bytes={already}-"} if already else {}

    try:
        response = getter(url, stream=True, timeout=timeout, headers=headers)
    except requests.Timeout as exc:
        raise _Retryable(f"timed out after {timeout}s") from exc
    except requests.RequestException as exc:
        raise _Retryable(f"connection failed: {exc}") from exc

    with response:
        status = response.status_code
        if status == 416:
            # Requested Range Not Satisfiable: the .part file is at least as
            # long as the resource, so it belongs to a different file or a
            # previous release.
            raise _Retryable("the partial file does not match the server's copy "
                             "(HTTP 416); starting over", status_code=416,
                             restart=True)
        if status in RETRYABLE_STATUS:
            raise _Retryable(f"HTTP {status}", status_code=status,
                             retry_after=retry_after_seconds(response))
        if status >= 400:
            raise DownloadError(f"Download of {url} failed with HTTP {status}",
                                url=url, status_code=status, response=response)

        resuming = status == 206
        if already and not resuming:
            # The server ignored the Range header, so its bytes start at zero
            # and appending them would corrupt the file.
            log.info("%s does not support resumption; downloading %s in full",
                     _host_of(url), label)
            already = 0

        declared = response.headers.get("content-length")
        remaining = int(declared) if declared is not None and declared.isdigit() else None
        total = already + remaining if remaining is not None else None

        mode = "ab" if already else "wb"
        bar = tqdm(total=total, initial=already, unit="B", unit_scale=True,
                   desc=f"Downloading {label}", disable=not progress)
        written = already
        try:
            with open(part, mode) as handle:
                for chunk in response.iter_content(chunk_size=chunk_size):
                    if not chunk:
                        continue
                    handle.write(chunk)
                    written += len(chunk)
                    bar.update(len(chunk))
        except requests.RequestException as exc:
            raise _Retryable(f"transfer interrupted: {exc}") from exc
        finally:
            bar.close()

    if total is not None and written < total:
        # A connection dropped mid-transfer without raising; the bytes on disk
        # are a valid prefix, so the next attempt resumes from here.
        raise _Retryable(
            f"transfer ended at {written} bytes of {total} declared by the server"
        )


def _resume_hint(part: str) -> str:
    """Tell the user whether a retry would pick up where this one stopped."""
    if os.path.exists(part):
        return (f"{os.path.getsize(part) / 1e6:.1f} MB are on disk at {part}; "
                f"running the download again resumes from there.")
    return "Nothing usable was left on disk; running the download again starts fresh."


def _host_of(url: str) -> str:
    """Host of a URL, for a log line that names who misbehaved."""
    return urlsplit(url).netloc or url


# ──────────────────────────────────────────────────────────────────────────────
# The dataset registry
#
# Everything above moves bytes; everything below answers "which bytes, and do I
# already have them?".  On a clean machine one CAS lookup through
# ``Search`` used to fetch ~32 GB without asking, because every
# source client defaults to ``auto_download=True`` and Search constructs all of
# them.  The package is for researchers on laptops, where 32 GB is often the
# whole free disk, so the download has to become something the user asks for.
#
# That needs a place where the five datasets are described *without* opening
# them: their filenames, their sizes, and what each one is for.  The client
# classes cannot serve as that place --- constructing one is exactly the act
# we are trying to avoid --- so the facts live here, in a table, and the
# clients are imported lazily inside ``fetch``.
# ──────────────────────────────────────────────────────────────────────────────

_MB = 1024 ** 2


@dataclass(frozen=True)
class Dataset:
    """
    One bulk dataset PROVESID can read offline, described without opening it.

    Attributes:
        name: Key used everywhere in this module, and the same string
            [`Search`][provesid.search.Search] uses for the source it feeds.
        title: Human-readable name for logs and tables.
        role: One line on what the dataset contributes to a search.
        patterns: Filenames, relative to the data directory, whose presence
            means the dataset is installed. Globs rather than plain names,
            because ChEMBL names its file after the release and ZeroPM after
            the version.
        extras: Globs for files that belong to the dataset but do not prove it
            is installed --- a derived index, a leftover ``.part``. They count
            towards the space it occupies and are removed with it.
        download_bytes: Size of the transfer, measured.
        resident_bytes: Size on disk once installed, measured, including
            anything built on first use.
        peak_bytes: Most disk needed at any one moment during installation.
            Larger than ``resident_bytes`` only for ChEMBL, which downloads a
            5.8 GB archive, extracts a 27.7 GiB database beside it and only
            then compacts that into the 2.4 GiB it keeps.
        source: Where the file comes from, for messages that have to tell a
            user what is about to be fetched.
        note: Anything a user deciding whether to fetch this should know.

    Examples:
        >>> chebi = DATASETS["chebi"]
        >>> chebi.title, chebi.patterns
        ('ChEBI SDF', ('chebi.sdf',))
        >>> human_bytes(DATASETS["chembl"].peak_bytes)
        '33.4 GiB'
    """

    name: str
    title: str
    role: str
    patterns: Tuple[str, ...]
    download_bytes: int
    resident_bytes: int
    source: str
    extras: Tuple[str, ...] = ()
    peak_bytes: int = 0
    note: str = ""

    def __post_init__(self) -> None:
        if not self.peak_bytes:
            object.__setattr__(self, "peak_bytes",
                               max(self.download_bytes, self.resident_bytes))


PUBCHEM_FTP_DOWNLOAD = 15_374_591_318
"""The PubChem identifier database as
[`provesid.pubchem_ftp`][provesid.pubchem_ftp] builds it, measured on the first
real build (2026-09-01 snapshot): the eight source files, the finished
database, and the largest file, ``CID-InChI-Key.gz``, which is the only one on
disk beside the database at the worst moment.
"""
PUBCHEM_FTP_RESIDENT = 2_484_281_344
PUBCHEM_FTP_LARGEST_FILE = 7_361_682_757

DATASETS: Dict[str, Dataset] = {
    "pubchem": Dataset(
        name="pubchem",
        title="PubChem identifiers",
        role="CAS, name, InChIKey and formula lookups; the broadest source",
        patterns=("pubchem_id.db",),
        # ``pubchem_ftp/<release>/`` holds the source files of a build: while
        # it runs, after an interrupted one, or for good with
        # ``keep_downloads=True``. The ``.part`` pair is the Zenodo route's.
        extras=("pubchem_id.db.tmp", "pubchem_ftp/*/*",
                "pubchem_id.db.part", "pubchem_id.db.part.source"),
        # Sized for the route ``fetch`` takes, ``PubChemID(source="ftp")``:
        # eight files from the newest monthly snapshot, read one at a time and
        # deleted, so the peak is the finished database plus the largest file.
        download_bytes=PUBCHEM_FTP_DOWNLOAD,
        resident_bytes=PUBCHEM_FTP_RESIDENT,
        peak_bytes=PUBCHEM_FTP_RESIDENT + PUBCHEM_FTP_LARGEST_FILE,
        source="PubChem FTP (Compound/Monthly/<newest>/Extras)",
        note="built from eight PubChem FTP files, one at a time: about 12 "
             "minutes of processing on top of a 15.4 GB download; "
             "PubChemID(source='zenodo') downloads a 2.2 GiB prebuilt copy "
             "instead",
    ),
    "comptox": Dataset(
        name="comptox",
        title="EPA CompTox chemicals",
        role="DTXSID lookups, and curated CAS-name pairs",
        patterns=("comptox_chemicals.db",),
        extras=("comptox_chemicals.db.part", "comptox_chemicals.db.part.source"),
        download_bytes=856293376,
        # The download plus the name index CompToxID builds into it (290 MiB),
        # measured 2026-09-22.
        resident_bytes=1160740864,
        source="Zenodo record 18833587",
        note="a name index (~20 s to build, 290 MiB) is added to the file "
             "after the download, so exact name lookups find synonyms",
    ),
    "chebi": Dataset(
        name="chebi",
        title="ChEBI SDF",
        role="curated structures, synonyms and ChEBI IDs",
        patterns=("chebi.sdf",),
        # The index is built on first use and is a fifth of the total; a status
        # table that ignored it would understate ChEBI by 74 MiB.
        extras=("chebi.sdf.index.pkl", "chebi.sdf.gz",
                "chebi.sdf.gz.part", "chebi.sdf.gz.part.source", "chebi.sdf.tmp"),
        download_bytes=250 * _MB,
        resident_bytes=922455587 + 78104601,
        source="EBI FTP (chebi.sdf.gz)",
        note="downloaded gzipped (~250 MB) and expanded; the search index is "
             "built on first use and takes a few minutes",
    ),
    "chembl": Dataset(
        name="chembl",
        title="ChEMBL",
        role="enrichment only --- adds ChEMBL IDs to structures already found",
        patterns=("chembl_*.db",),
        # The ``sqlite`` route saves its archive as ``chembl_NN.db.tar.gz``;
        # the ``mysql`` route keeps the dump's own name. Either can be left
        # half-downloaded as a ``.part`` pair.
        extras=("chembl_*.db.tar.gz", "chembl_*.db.tar.gz.part",
                "chembl_*.db.tar.gz.part.source",
                "chembl_*_mysql.tar.gz", "chembl_*_mysql.tar.gz.part",
                "chembl_*_mysql.tar.gz.part.source",
                "chembl_*.db.incoming", "chembl_*.db.tmp"),
        download_bytes=5800 * _MB,
        # What is left when the install finishes: the extract, not the release.
        # ``CheMBL(source="sqlite")`` -- the default -- compacts the 27.7 GiB
        # database into 2.42 GiB and deletes it, so what this dataset costs to
        # install and what it costs to keep are an order of magnitude apart,
        # and ``peak_bytes`` below is the number that decides whether it fits.
        resident_bytes=2599391232,
        peak_bytes=29739835392 + 5800 * _MB,
        source="EBI FTP (chembl_NN_sqlite.tar.gz)",
        note="installs as a 2.4 GiB extract, but the 5.8 GB archive and the "
             "27.7 GiB release it is built from both exist on the way in --- "
             "33.4 GiB has to be free. CheMBL(source='mysql') builds the "
             "same extract from the 2.1 GB MySQL dump and needs ~4.5 GiB; "
             "CheMBL(source='full') keeps the release",
    ),
    "zeropm": Dataset(
        name="zeropm",
        title="ZeroPM inventory",
        role="regulatory inventories; off unless Search(use_zeropm=True)",
        patterns=("zeropm-*.sqlite",),
        extras=("zeropm-*.sqlite.part", "zeropm-*.sqlite.part.source"),
        download_bytes=459968512,
        resident_bytes=459968512,
        source="ZeroPM-H2020 GitHub repository",
    ),
}
"""The five datasets, in the order [`Search`][provesid.search.Search] benefits
from them: the three primary sources first, then ChEMBL, which only enriches a
structure the others already found, then ZeroPM, which is off by default.

Sizes were measured on 2026-09-20 from the copies on a machine that had all
five: ChEBI SDF 879.7 MiB plus a 74.5 MiB index, CompTox 816.6 MiB plus a 290.3
MiB name index, PubChem 2.31 GiB as the FTP build leaves it (see
``PUBCHEM_FTP_RESIDENT``), ChEMBL 36 2.42 GiB as the extract an install now
keeps (27.7 GiB as the full release it is built from, out of a 5.8 GB archive),
ZeroPM 438.7 MiB. They are advisory --- a later release is a little larger ---
and are used to tell the user what a download will cost before it starts, not
to check anything.
"""

DEFAULT_DATASETS: Tuple[str, ...] = ("pubchem", "comptox", "chebi", "chembl")
"""Datasets [`Search`][provesid.search.Search] queries unless ``use_zeropm=True``."""


class MissingDatasetError(ServiceError):
    """
    A dataset was needed and is not on disk.

    Raised by [`require`][provesid.datasets.require] --- and so by
    ``Search(datasets="required")`` --- instead of downloading tens of
    gigabytes on the user's behalf. The message names every missing dataset,
    what it costs, and the exact [`fetch`][provesid.datasets.fetch] call that
    would install it.

    A [`ServiceError`][provesid.http.ServiceError] so that the whole family stays
    catchable through one base, as
    [`DownloadError`][provesid.datasets.DownloadError] is.

    Examples:
        >>> import tempfile
        >>> require("chembl", tempfile.mkdtemp())
        Traceback (most recent call last):
        ...
        provesid.datasets.MissingDatasetError: 1 dataset(s) missing from ...
    """


def human_bytes(count: float) -> str:
    """
    Format a byte count the way a user reading a size wants it.

    Args:
        count: Number of bytes.

    Returns:
        The size in the largest unit that leaves a number above 1, binary
        units, one decimal place.

    Examples:
        >>> human_bytes(2322595840)
        '2.2 GiB'
        >>> human_bytes(0)
        '0 B'
    """
    if count < 1024:
        return f"{int(count)} B"
    for unit in ("KiB", "MiB", "GiB", "TiB"):
        count /= 1024
        if count < 1024 or unit == "TiB":
            return f"{count:.1f} {unit}"
    return f"{count:.1f} TiB"  # pragma: no cover - unreachable, kept explicit


def dataset_names() -> List[str]:
    """
    Names of every dataset in the registry, in registry order.

    Returns:
        The keys of [`DATASETS`][provesid.datasets.DATASETS].

    Examples:
        >>> dataset_names()
        ['pubchem', 'comptox', 'chebi', 'chembl', 'zeropm']
    """
    return list(DATASETS)


def _resolve_names(names: Optional[Union[str, Iterable[str]]]) -> List[str]:
    """
    Normalise a name, an iterable of names, or None into a list of valid names.

    Args:
        names: One dataset name, several, or None for every dataset in the
            registry.

    Returns:
        Dataset names in registry order, without duplicates.

    Raises:
        KeyError: If a name is not in the registry. The message lists the
            names that are, since a typo here is otherwise reported much later
            as an empty table.
    """
    if names is None:
        return dataset_names()
    if isinstance(names, str):
        names = [names]
    requested = list(dict.fromkeys(names))
    unknown = [name for name in requested if name not in DATASETS]
    if unknown:
        raise KeyError(
            f"Unknown dataset(s): {', '.join(sorted(unknown))}. "
            f"Known datasets: {', '.join(dataset_names())}."
        )
    return [name for name in dataset_names() if name in requested]


def data_directory(data_dir: Optional[str] = None) -> str:
    """
    The directory the datasets live in.

    Args:
        data_dir: An explicit directory, which is returned as given (expanded
            and made absolute). None for the per-user default, which honours
            ``PROVESID_DATA_DIR``.

    Returns:
        Absolute path to the dataset directory.

    Examples:
        >>> data_directory("/data/provesid")
        '/data/provesid'
        >>> data_directory() == user_dataset_path()
        True
    """
    if data_dir is not None:
        return os.path.abspath(os.path.expanduser(str(data_dir)))
    return user_dataset_path()


def dataset_files(name: str, data_dir: Optional[str] = None,
                  *, include_extras: bool = True) -> List[str]:
    """
    Files belonging to one dataset that are actually on disk.

    Args:
        name: Dataset name.
        data_dir: Directory to look in; None for the default.
        include_extras: Include derived and leftover files --- ChEBI's index, a
            ``.part`` from an interrupted download. They occupy real space, so
            [`status`][provesid.datasets.status] and
            [`remove`][provesid.datasets.remove] want them;
            [`is_present`][provesid.datasets.is_present] does not.

    Returns:
        Absolute paths, sorted, of the files that exist.

    Raises:
        KeyError: If ``name`` is not a known dataset.

    Examples:
        >>> import tempfile
        >>> directory = tempfile.mkdtemp()
        >>> dataset_files("zeropm", directory)
        []
        >>> open(os.path.join(directory, "zeropm-v0-0-4.sqlite"), "w").close()
        >>> [os.path.basename(path) for path in dataset_files("zeropm", directory)]
        ['zeropm-v0-0-4.sqlite']
    """
    dataset = DATASETS[name]
    directory = data_directory(data_dir)
    globs = dataset.patterns + (dataset.extras if include_extras else ())
    found: List[str] = []
    for pattern in globs:
        found.extend(glob.glob(os.path.join(directory, pattern)))
    return sorted(dict.fromkeys(os.path.abspath(path) for path in found))


def is_present(name: str, data_dir: Optional[str] = None) -> bool:
    """
    Whether a dataset is installed and usable.

    A leftover ``.part`` does not count, and neither does ChEBI's index on its
    own: both are files the dataset leaves behind rather than the dataset.

    Args:
        name: Dataset name.
        data_dir: Directory to look in; None for the default.

    Returns:
        True when at least one file matching the dataset's own patterns exists.
        Only the name is checked; an empty or damaged file counts.

    Raises:
        KeyError: If ``name`` is not a known dataset.

    Examples:
        >>> import tempfile
        >>> is_present("chembl", tempfile.mkdtemp())
        False
    """
    return bool(dataset_files(name, data_dir, include_extras=False))


def missing(names: Optional[Union[str, Iterable[str]]] = None,
            data_dir: Optional[str] = None) -> List[str]:
    """
    Which of the named datasets are not installed.

    Args:
        names: Dataset name, names, or None for all of them.
        data_dir: Directory to look in; None for the default.

    Returns:
        Names of the datasets that are absent, in registry order.

    Raises:
        KeyError: If a name is not in the registry.

    Examples:
        >>> missing(["pubchem", "chembl"])            # doctest: +SKIP
        ['chembl']
    """
    return [name for name in _resolve_names(names)
            if not is_present(name, data_dir)]


def fetch_command(names: Union[str, Iterable[str]]) -> str:
    """
    The exact call that installs the given datasets, as a string.

    Error messages that tell a user what went wrong should also tell them what
    to type; this builds that line so every message spells it the same way.

    Args:
        names: Dataset name or names.

    Returns:
        A copy-pasteable Python call.

    Examples:
        >>> fetch_command("chembl")
        "provesid.datasets.fetch('chembl')"
        >>> fetch_command(["pubchem", "chebi"])
        "provesid.datasets.fetch(['pubchem', 'chebi'])"
    """
    resolved = _resolve_names(names)
    if len(resolved) == 1:
        return f"provesid.datasets.fetch({resolved[0]!r})"
    return f"provesid.datasets.fetch({resolved!r})"


def _release_of(name: str, path: str) -> str:
    """
    The release a file on disk belongs to, read from its name.

    Only two of the five datasets carry a version in the filename, and both
    matter to a user reading a status table: ChEMBL because the release number
    changes what is in it, ZeroPM because the inventory is versioned. The
    others are single-file snapshots whose version lives in the Zenodo record
    they came from, and inventing a version for them would be worse than
    leaving the column empty.

    Args:
        name: Dataset name.
        path: Path to one of its files.

    Returns:
        A short release string, or ``""`` when the filename does not say.
    """
    if not path:
        return ""
    basename = os.path.basename(path)
    if name == "chembl":
        match = re.search(r"chembl_(\d+)", basename)
        if not match:
            return ""
        return (f"{match.group(1)} (extract)" if "_provesid" in basename
                else match.group(1))
    if name == "zeropm":
        match = re.search(r"v(\d+)-(\d+)-(\d+)", basename)
        return ".".join(match.groups()) if match else ""
    return ""


def _preferred_file(name: str, paths: List[str]) -> str:
    """
    The one file a client would open, out of several a dataset may have.

    Only ChEMBL and ZeroPM can have more than one: ChEMBL names its file after
    the release and leaves older ones in place, and
    [`CheMBL.compact`][provesid.chembl.CheMBL.compact] writes an extract beside
    the full release it was built from. The rule here is
    `CheMBL._find_local_database`'s --- newest release wins, and at equal
    releases the extract, since the two answer identically and one costs a
    tenth of the page cache.

    It reads the rule off the filename rather than opening each candidate, so
    a status table stays a directory listing. An extract renamed to hide its
    ``_provesid`` suffix would be misreported here and opened anyway by
    [`CheMBL`][provesid.chembl.CheMBL], which checks the file itself.

    Args:
        name: Dataset name.
        paths: Candidate paths, as
            [`dataset_files`][provesid.datasets.dataset_files] returns them.

    Returns:
        The preferred path, or ``""`` when there are no candidates.
    """
    if not paths:
        return ""
    if name == "chembl":
        def key(path: str) -> Tuple[int, int]:
            match = re.search(r"chembl_(\d+)", os.path.basename(path))
            return (int(match.group(1)) if match else -1,
                    1 if "_provesid" in os.path.basename(path) else 0)
        return max(paths, key=key)
    if name == "zeropm":
        return max(paths, key=lambda path: _release_of(name, path))
    return paths[0]


def status(names: Optional[Union[str, Iterable[str]]] = None,
           data_dir: Optional[str] = None) -> pd.DataFrame:
    """
    What is on disk, dataset by dataset.

    The first thing to run on a machine whose disk is filling up, and the
    answer to "will this search use all four sources?". Nothing is downloaded,
    nothing is opened --- the table is built from filenames and ``stat`` calls,
    so it is instant even with 30 GB of ChEMBL in the directory.

    Args:
        names: Dataset name, names, or None for every dataset.
        data_dir: Directory to look in; None for the per-user default.

    Returns:
        A DataFrame with one row per dataset and the columns:

        ``dataset``
            Registry name, the same string ``Search`` uses for the source.
        ``title``
            Human-readable name.
        ``present``
            Whether the dataset itself is installed. A leftover ``.part`` or a
            stale index does not make this True, though both are counted in
            ``bytes``.
        ``files``
            Number of files found, including derived and partial ones.
        ``bytes`` / ``size``
            Space occupied, as an integer and as a readable string.
        ``release``
            Release read from the filename, where the filename says.
        ``path``
            The dataset's main file, or ``""`` when it is absent.

        ``df.attrs`` carries ``data_dir`` and ``total_bytes``.

    Raises:
        KeyError: If a name is not in the registry.

    Examples:
        >>> from provesid import datasets
        >>> datasets.status()[["dataset", "present", "size"]]   # doctest: +SKIP
          dataset  present      size
        0 pubchem     True   2.2 GiB
        1 comptox     True     1.1 GiB
        2   chebi     True   954.2 MiB
        3  chembl    False         0 B
        4  zeropm     True   438.7 MiB
    """
    directory = data_directory(data_dir)
    rows: List[Dict[str, Any]] = []
    for name in _resolve_names(names):
        dataset = DATASETS[name]
        files = dataset_files(name, directory)
        primary = dataset_files(name, directory, include_extras=False)
        total = sum(os.path.getsize(path) for path in files
                    if os.path.exists(path))
        rows.append({
            "dataset": name,
            "title": dataset.title,
            "present": bool(primary),
            "files": len(files),
            "bytes": total,
            "size": human_bytes(total),
            "release": _release_of(name, _preferred_file(name, primary)),
            "path": _preferred_file(name, primary),
        })

    frame = pd.DataFrame(rows, columns=["dataset", "title", "present", "files",
                                        "bytes", "size", "release", "path"])
    frame.attrs["data_dir"] = directory
    frame.attrs["total_bytes"] = int(frame["bytes"].sum()) if rows else 0
    return frame


def plan(names: Optional[Union[str, Iterable[str]]] = None,
         data_dir: Optional[str] = None, *, force: bool = False) -> pd.DataFrame:
    """
    What [`fetch`][provesid.datasets.fetch] would download, and how much disk
    it would take.

    The question §4.1 of the refactor plan says nobody was asked: a clean
    machine used to spend 32 GB on one CAS lookup without a word. Run this
    first and the number is on screen before anything is transferred.

    The sizes are the measured ones in
    [`DATASETS`][provesid.datasets.DATASETS], so they are advisory: a newer
    ChEMBL release is a little larger than the one they were taken from. They
    are the right order of magnitude, which is what the decision turns on.

    Args:
        names: Dataset name, names, or None for every dataset.
        data_dir: Directory the datasets would go in; None for the default.
        force: Plan a re-download of datasets that are already installed, as
            ``fetch(force=True)`` would.

    Returns:
        A DataFrame with one row per dataset and the columns ``dataset``,
        ``action`` (``"download"`` or ``"present"``), ``download`` /
        ``installed`` (readable sizes), ``download_bytes`` /
        ``resident_bytes`` / ``peak_bytes``, ``role`` and ``note``.

        ``df.attrs`` carries ``data_dir``, ``total_download_bytes``,
        ``total_resident_bytes`` and ``peak_bytes`` --- the last being the most
        disk needed at any one moment, which for ChEMBL is more than ten times
        the installed size, because the 2.4 GiB extract is built from a
        27.7 GiB release that is downloaded, unpacked and then deleted.

    Raises:
        KeyError: If a name is not in the registry.

    Examples:
        >>> from provesid import datasets
        >>> todo = datasets.plan(["pubchem", "chebi"])        # doctest: +SKIP
        >>> datasets.human_bytes(                             # doctest: +SKIP
        ...     todo.attrs["total_download_bytes"])
        '2.4 GiB'
    """
    directory = data_directory(data_dir)
    rows: List[Dict[str, Any]] = []
    for name in _resolve_names(names):
        dataset = DATASETS[name]
        would_download = force or not is_present(name, directory)
        rows.append({
            "dataset": name,
            "action": "download" if would_download else "present",
            "download": human_bytes(dataset.download_bytes) if would_download else "-",
            "installed": human_bytes(dataset.resident_bytes) if would_download else "-",
            "download_bytes": dataset.download_bytes if would_download else 0,
            "resident_bytes": dataset.resident_bytes if would_download else 0,
            "peak_bytes": dataset.peak_bytes if would_download else 0,
            "role": dataset.role,
            "note": dataset.note,
        })

    frame = pd.DataFrame(rows, columns=["dataset", "action", "download", "installed",
                                        "download_bytes", "resident_bytes",
                                        "peak_bytes", "role", "note"])
    frame.attrs["data_dir"] = directory
    frame.attrs["total_download_bytes"] = int(frame["download_bytes"].sum()) if rows else 0
    frame.attrs["total_resident_bytes"] = int(frame["resident_bytes"].sum()) if rows else 0
    # Peak disk is not the sum of the peaks: the datasets are installed one
    # after another, so the worst moment is everything else already on disk
    # plus the largest transient overhead of a single install -- ChEMBL's, whose
    # archive and full release are both deleted once the extract is built.
    overhead = (frame["peak_bytes"] - frame["resident_bytes"]).max() if rows else 0
    frame.attrs["peak_bytes"] = int(frame["resident_bytes"].sum() + overhead) if rows else 0
    return frame


_CLIENTS: Dict[str, Tuple[str, str]] = {
    "pubchem": (".pubchem_id", "PubChemID"),
    "comptox": (".comptox", "CompToxID"),
    "chebi": (".chebi_sdf", "ChebiSDF"),
    "chembl": (".chembl", "CheMBL"),
    "zeropm": (".zeropm", "ZeroPM"),
}
"""Import paths of the client classes, resolved only when
[`fetch`][provesid.datasets.fetch] runs. The client modules import *this* one
for [`download_file`][provesid.datasets.download_file], so the registry cannot
import them back at module level.
"""


def _client_class(name: str) -> type:
    """
    Import and return the client class that installs a dataset.

    Args:
        name: Dataset name.

    Returns:
        The class, e.g. [`PubChemID`][provesid.pubchem_id.PubChemID] for ``"pubchem"``.
    """
    module_name, class_name = _CLIENTS[name]
    module = importlib.import_module(module_name, package=__package__)
    return getattr(module, class_name)


def fetch(names: Union[str, Iterable[str]], data_dir: Optional[str] = None,
          *, force: bool = False, progress: bool = True) -> Dict[str, str]:
    """
    Download and install datasets, by name.

    Each dataset is installed by constructing its client with
    ``auto_download=True``, which is the one code path that knows how to
    finish the job: ChEBI's SDF has to be expanded from gzip and indexed,
    ChEMBL's archive extracted and checked, and all five are verified before
    anything is moved into place. The client is closed again --- the point here
    is the files, not the connection.

    Datasets already present are skipped unless ``force=True``, so calling this
    on a list is cheap and repeatable.

    Args:
        names: Dataset name or names. There is no "all" default: fetching
            everything transfers ~21 GiB and needs ~37 GiB free while ChEMBL is
            unpacked, which is a decision that has to be spelled out.
        data_dir: Directory to install into; None for the per-user default.
        force: Re-download datasets that are already installed. For ChEMBL this
            fetches the current release, which may be a newer one.
        progress: Show the per-file progress bar.

    Returns:
        Mapping of dataset name to the path of its main file.

    Raises:
        KeyError: If a name is not in the registry.
        DownloadError: If a transfer could not be completed.

    Examples:
        >>> from provesid import datasets
        >>> datasets.fetch(["pubchem", "chebi"])            # doctest: +SKIP
        {'pubchem': '/home/me/.local/share/provesid/pubchem_id.db',
         'chebi': '/home/me/.local/share/provesid/chebi.sdf'}
    """
    directory = data_directory(data_dir)
    os.makedirs(directory, exist_ok=True)
    requested = _resolve_names(names)
    todo = [name for name in requested
            if force or not is_present(name, directory)]

    if todo:
        # The total, before the first byte moves. This is the announcement
        # whose absence made a first run cost 32 GB unasked.
        upcoming = plan(todo, directory, force=force)
        logger.info(
            "Fetching %d dataset(s) into %s: %s. Download %s, %s on disk when "
            "done, %s needed at peak.",
            len(todo), directory, ", ".join(todo),
            human_bytes(upcoming.attrs["total_download_bytes"]),
            human_bytes(upcoming.attrs["total_resident_bytes"]),
            human_bytes(upcoming.attrs["peak_bytes"]),
        )

    installed: Dict[str, str] = {}
    for name in requested:
        if name not in todo:
            found = _preferred_file(name, dataset_files(name, directory,
                                                        include_extras=False))
            logger.info("%s already present at %s", DATASETS[name].title, found)
            installed[name] = found
            continue

        logger.info("Fetching %s (%s, %s)", DATASETS[name].title,
                    human_bytes(DATASETS[name].download_bytes),
                    DATASETS[name].source)
        client = _client_class(name)(
            auto_download=True, data_dir=directory, redownload=force,
        )
        try:
            found = _preferred_file(name, dataset_files(name, directory,
                                                        include_extras=False))
            if not found:  # pragma: no cover - the client would have raised
                raise DownloadError(
                    f"{DATASETS[name].title} reported success but left no file "
                    f"matching {DATASETS[name].patterns} in {directory}"
                )
            installed[name] = found
        finally:
            connection = getattr(client, "conn", None)
            if connection is not None:
                connection.close()
            del client

    return installed


def remove(names: Union[str, Iterable[str]], data_dir: Optional[str] = None,
           *, dry_run: bool = False) -> pd.DataFrame:
    """
    Delete datasets from disk, by name, and report the space reclaimed.

    Deletes the dataset's own files and everything derived from them --- an
    index, a half-finished ``.part``, ChEMBL's extracted archive --- because
    leaving those behind reclaims a fraction of the space and confuses the next
    [`status`][provesid.datasets.status]. For ChEMBL that means *every* release
    and extract in the directory, not only the one a client would open;
    ``dry_run=True`` lists them first.

    This is destructive and there is no undo beyond fetching again, so pass
    ``dry_run=True`` first to see the list. Naming the datasets explicitly is
    deliberate: there is no "all".

    Args:
        names: Dataset name or names.
        data_dir: Directory to delete from; None for the per-user default.
        dry_run: List what would go without deleting anything.

    Returns:
        A DataFrame with one row per file and the columns ``dataset``,
        ``path``, ``bytes``, ``size`` and ``removed``. ``df.attrs`` carries
        ``freed_bytes`` --- what was reclaimed, or what would be.

    Raises:
        KeyError: If a name is not in the registry.

    Examples:
        >>> from provesid import datasets
        >>> datasets.remove("chembl", dry_run=True)          # doctest: +SKIP
        >>> datasets.remove("chembl")                        # doctest: +SKIP
    """
    directory = data_directory(data_dir)
    rows: List[Dict[str, Any]] = []
    for name in _resolve_names(names):
        for path in dataset_files(name, directory):
            size = os.path.getsize(path) if os.path.exists(path) else 0
            removed = False
            if not dry_run:
                try:
                    os.remove(path)
                    removed = True
                except OSError as exc:
                    # A Windows client still holding the database open is the
                    # common case; say which file and carry on with the rest.
                    logger.warning("Could not remove %s: %s", path, exc)
            rows.append({
                "dataset": name,
                "path": path,
                "bytes": size,
                "size": human_bytes(size),
                "removed": removed,
            })

    frame = pd.DataFrame(rows, columns=["dataset", "path", "bytes", "size", "removed"])
    freed = int(frame.loc[frame["removed"] | dry_run, "bytes"].sum()) if rows else 0
    frame.attrs["data_dir"] = directory
    frame.attrs["freed_bytes"] = freed
    logger.info("%s %s across %d file(s)",
                "Would free" if dry_run else "Freed", human_bytes(freed), len(rows))
    return frame


def require(names: Union[str, Iterable[str]], data_dir: Optional[str] = None) -> None:
    """
    Raise unless every named dataset is installed.

    What ``Search(datasets="required")`` calls, and what any code that must not
    silently run on fewer sources should call. The message is the useful part:
    it names each missing dataset with its size and ends with the exact
    [`fetch`][provesid.datasets.fetch] call, so the user never has to look one up.

    Args:
        names: Dataset name or names.
        data_dir: Directory to look in; None for the per-user default.

    Raises:
        MissingDatasetError: If any named dataset is absent. Nothing is
            downloaded.
        KeyError: If a name is not in the registry.

    Examples:
        >>> from provesid import datasets
        >>> datasets.require(["pubchem", "chembl"])          # doctest: +SKIP
        Traceback (most recent call last):
        provesid.datasets.MissingDatasetError: 1 dataset is missing from ...
    """
    absent = missing(names, data_dir)
    if not absent:
        return

    directory = data_directory(data_dir)
    lines = [
        f"  {DATASETS[name].title} ({name}): {human_bytes(DATASETS[name].download_bytes)}"
        f" to download, {human_bytes(DATASETS[name].resident_bytes)} on disk"
        f" --- {DATASETS[name].role}"
        for name in absent
    ]
    raise MissingDatasetError(
        f"{len(absent)} dataset(s) missing from {directory}:\n"
        + "\n".join(lines)
        + f"\n\nInstall them with:\n  {fetch_command(absent)}\n"
        "Or pass datasets='present' to run on whatever is already on disk, "
        "or datasets='auto' to download automatically."
    )
