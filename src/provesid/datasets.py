"""
One resumable, checksummed downloader for the bulk datasets PROVESID reads.

:mod:`provesid.http` is the transport for web *APIs* --- small requests, paced
at five a second, retried a few times.  The bulk datasets are a different
problem entirely: ChEMBL's archive is 5.8 GB, PubChem's identifier database
2.2 GB, and what those transfers want is not a pacer but *resumption*, a
checksum and an atomic rename.  So they were deliberately left out of the
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

:func:`download_file` is the one implementation.  It resumes against a
``.part`` file with an HTTP ``Range`` request, verifies an MD5 when the server
publishes one, checks the byte count it actually received against the one the
server declared, hands the finished file to the caller's own validity check,
and only then moves it into place.  Nothing is ever renamed onto the
destination until every one of those has passed, so a failed download can
never replace a good database with a broken one.

Example:
    >>> from provesid.datasets import download_file
    >>> download_file(                                    # doctest: +SKIP
    ...     "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz",
    ...     "/data/CID-SMILES.gz",
    ...     checksum_url="https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz.md5",
    ... )
    '/data/CID-SMILES.gz'
"""

import hashlib
import logging
import os
import time
from typing import Callable, Optional
from urllib.parse import urlsplit

import requests
from tqdm import tqdm

from .http import RETRYABLE_STATUS, ServiceError, retry_after_seconds

logger = logging.getLogger(__name__)

#: Bytes read from the socket, and from disk while checksumming, at a time.
#: A megabyte is large enough that the per-chunk work disappears against the
#: transfer and small enough that the progress bar still moves.
CHUNK_SIZE = 1024 * 1024

#: Suffix of the partial file a download accumulates into.  It sits beside the
#: destination rather than in a temporary directory so that a resumed download
#: finds it, and so that the final move is a rename within one filesystem.
PART_SUFFIX = ".part"

#: Suffix of the marker recording which URL a partial file came from.  Without
#: it, a ``.part`` left over from a different release --- or a different
#: dataset that happens to share a destination --- would be resumed, splicing
#: two files together.  A checksum would catch that, but only the PubChem FTP
#: mirror publishes one; the Zenodo downloads have no such backstop.
SOURCE_SUFFIX = PART_SUFFIX + ".source"


class DownloadError(ServiceError):
    """
    A bulk dataset could not be downloaded, or arrived damaged.

    A :class:`~provesid.http.ServiceError` like any other failure this package
    reports from a remote service, so a caller can catch the whole family; the
    separate class exists because the recovery differs. An API call that fails
    is retried or abandoned, while a download that fails has usually left a
    resumable ``.part`` file on disk and will pick up where it stopped on the
    next call.

    The exception is a checksum mismatch, which deletes the partial file: the
    bytes on disk are known to be wrong, so resuming from them would only
    produce the same wrong file again.
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

    Example:
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

    Example:
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

    Example:
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
