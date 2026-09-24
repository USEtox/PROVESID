"""
Tests for ``provesid.datasets.download_file`` — the bulk downloader.

The five modules that download a dataset each had their own copy of "stream the
response into a temporary file with a progress bar", none of which retried,
resumed or checksummed.  ``download_file`` is the one implementation, and what
has to be pinned is exactly the behaviour those copies lacked: that an
interrupted transfer continues rather than restarting, that a damaged file is
caught before it replaces a good one, and that nothing lands at the destination
unless every check passed.

These tests run a real HTTP server on localhost rather than stubbing
``requests``, because the interesting behaviour is in the ``Range`` request and
the response status — the parts a stub would have to fake, and so the parts a
stub would let us get wrong.
"""

import hashlib
import http.server
import os
import threading

import pytest

from provesid.datasets import (
    DownloadError,
    PART_SUFFIX,
    SOURCE_SUFFIX,
    download_file,
    md5_of_file,
    read_checksum,
)


PAYLOAD = bytes(range(256)) * 400          # 102 400 bytes, and not compressible
PAYLOAD_MD5 = hashlib.md5(PAYLOAD).hexdigest()


# ── A server that can be made to misbehave ───────────────────────────────────

class _Handler(http.server.BaseHTTPRequestHandler):
    """Serve ``PAYLOAD``, with the faults a download has to survive.

    The server's ``plan`` decides what each request gets:

    ``statuses``
        A list of statuses to answer with before serving anything, one per
        request — ``[503, 503]`` makes the first two attempts fail.
    ``cut_after``
        Bytes to send before hanging up, once per entry, simulating a dropped
        connection with the full length already declared.
    ``ignore_range``
        Answer 200 with the whole file even when a ``Range`` was asked for,
        which is what a server without resumption support does.
    ``corrupt``
        Serve a byte-for-byte wrong file of the right length.
    """

    def log_message(self, *args):
        pass  # keep the test output readable

    @property
    def plan(self):
        return self.server.plan

    def do_GET(self):
        if self.path.endswith(".md5"):
            body = f"{PAYLOAD_MD5}  payload.bin\n".encode()
            self.send_response(200)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return

        self.plan.setdefault("requests", []).append(self.headers.get("Range"))

        if self.plan.get("statuses"):
            status = self.plan["statuses"].pop(0)
            self.send_response(status)
            self.send_header("Content-Length", "0")
            self.end_headers()
            return

        payload = bytes(len(PAYLOAD)) if self.plan.get("corrupt") else PAYLOAD

        start = 0
        requested = self.headers.get("Range")
        if requested and not self.plan.get("ignore_range"):
            start = int(requested.split("=")[1].split("-")[0])
            if start >= len(payload):
                self.send_response(416)
                self.send_header("Content-Length", "0")
                self.end_headers()
                return
            self.send_response(206)
            self.send_header("Content-Range",
                             f"bytes {start}-{len(payload) - 1}/{len(payload)}")
        else:
            self.send_response(200)

        body = payload[start:]
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()

        if self.plan.get("cut_after"):
            cut = self.plan["cut_after"].pop(0)
            self.wfile.write(body[:cut])
            self.wfile.flush()
            self.close_connection = True
            return
        self.wfile.write(body)


@pytest.fixture
def server():
    """A localhost HTTP server whose ``plan`` the test can rewrite."""
    httpd = http.server.HTTPServer(("127.0.0.1", 0), _Handler)
    httpd.plan = {}
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    httpd.url = f"http://127.0.0.1:{httpd.server_port}/payload.bin"
    yield httpd
    httpd.shutdown()
    httpd.server_close()


@pytest.fixture
def dest(tmp_path):
    return str(tmp_path / "payload.bin")


# ── The ordinary case ────────────────────────────────────────────────────────

class TestAPlainDownload:

    def test_the_file_arrives_intact(self, server, dest):
        assert download_file(server.url, dest, progress=False) == dest
        assert open(dest, "rb").read() == PAYLOAD

    def test_the_partial_file_is_gone_afterwards(self, server, dest):
        download_file(server.url, dest, progress=False)
        assert not os.path.exists(dest + PART_SUFFIX)

    def test_missing_parent_directories_are_created(self, server, tmp_path):
        nested = str(tmp_path / "a" / "b" / "payload.bin")
        download_file(server.url, nested, progress=False)
        assert os.path.exists(nested)

    def test_an_existing_file_is_replaced(self, server, dest):
        open(dest, "wb").write(b"stale")
        download_file(server.url, dest, progress=False)
        assert open(dest, "rb").read() == PAYLOAD


# ── Resumption, which is the point of the module ─────────────────────────────

class TestResumption:

    def test_an_interrupted_transfer_continues_from_disk(self, server, dest):
        """The first attempt is cut short; the second asks for the rest."""
        server.plan["cut_after"] = [40_000]
        download_file(server.url, dest, progress=False, backoff=0, chunk_size=8192)
        assert open(dest, "rb").read() == PAYLOAD
        first, second = server.plan["requests"]
        assert first is None
        assert second.startswith("bytes=")
        assert int(second[len("bytes="):-1]) > 0

    def test_resumption_picks_up_from_what_was_written_not_what_arrived(
        self, server, dest
    ):
        """A broken connection loses the chunk it was in the middle of.

        ``iter_content`` yields whole chunks, so bytes received after the last
        complete chunk are never handed over and never reach the disk. The
        resumed request therefore asks from the last chunk boundary, not from
        the last byte the server managed to send. With a megabyte chunk that
        costs at most a megabyte per interruption, which is nothing against a
        5.8 GB archive -- but it means the offset is a multiple of the chunk
        size, and a test that assumed otherwise would be testing urllib3's
        buffering rather than this module.
        """
        server.plan["cut_after"] = [40_000]
        download_file(server.url, dest, progress=False, backoff=0, chunk_size=8192)
        resumed_from = int(server.plan["requests"][1][len("bytes="):-1])
        assert resumed_from % 8192 == 0
        assert 0 < resumed_from <= 40_000
        assert open(dest, "rb").read() == PAYLOAD

    def test_a_leftover_part_file_is_resumed_not_restarted(self, server, dest):
        open(dest + PART_SUFFIX, "wb").write(PAYLOAD[:30_000])
        open(dest + SOURCE_SUFFIX, "w").write(server.url)
        download_file(server.url, dest, progress=False)
        assert open(dest, "rb").read() == PAYLOAD
        assert server.plan["requests"] == ["bytes=30000-"]

    def test_resume_false_discards_what_is_on_disk(self, server, dest):
        open(dest + PART_SUFFIX, "wb").write(PAYLOAD[:30_000])
        open(dest + SOURCE_SUFFIX, "w").write(server.url)
        download_file(server.url, dest, progress=False, resume=False)
        assert open(dest, "rb").read() == PAYLOAD
        assert server.plan["requests"] == [None]

    def test_a_server_that_ignores_range_still_produces_the_right_file(
        self, server, dest
    ):
        """Answering 200 to a Range request means the bytes start at zero."""
        open(dest + PART_SUFFIX, "wb").write(PAYLOAD[:30_000])
        open(dest + SOURCE_SUFFIX, "w").write(server.url)
        server.plan["ignore_range"] = True
        download_file(server.url, dest, progress=False)
        assert open(dest, "rb").read() == PAYLOAD

    def test_a_part_file_from_another_url_is_discarded(self, server, dest):
        """Splicing two files together would produce a plausible-looking ruin.

        Only PubChem's FTP mirror publishes a checksum; the Zenodo downloads
        have no backstop, so the URL a partial came from is recorded beside it.
        """
        open(dest + PART_SUFFIX, "wb").write(b"bytes from some other dataset")
        open(dest + SOURCE_SUFFIX, "w").write("http://example.invalid/other.bin")
        download_file(server.url, dest, progress=False)
        assert open(dest, "rb").read() == PAYLOAD
        assert server.plan["requests"] == [None]

    def test_a_part_file_with_no_marker_is_discarded(self, server, dest):
        """A partial from before this was recorded cannot be vouched for."""
        open(dest + PART_SUFFIX, "wb").write(PAYLOAD[:30_000])
        download_file(server.url, dest, progress=False)
        assert open(dest, "rb").read() == PAYLOAD
        assert server.plan["requests"] == [None]

    def test_the_marker_is_removed_once_the_file_is_in_place(self, server, dest):
        download_file(server.url, dest, progress=False)
        assert not os.path.exists(dest + SOURCE_SUFFIX)

    def test_a_part_file_longer_than_the_resource_starts_over(self, server, dest):
        """HTTP 416: what is on disk cannot be a prefix of what is offered."""
        open(dest + PART_SUFFIX, "wb").write(PAYLOAD + b"extra")
        open(dest + SOURCE_SUFFIX, "w").write(server.url)
        download_file(server.url, dest, progress=False, backoff=0)
        assert open(dest, "rb").read() == PAYLOAD

    def test_several_interruptions_still_converge(self, server, dest):
        """Each attempt makes progress, so a flaky link finishes eventually."""
        server.plan["cut_after"] = [20_000, 20_000, 20_000]
        download_file(server.url, dest, progress=False, backoff=0, chunk_size=8192)
        assert open(dest, "rb").read() == PAYLOAD

        offsets = [
            0 if r is None else int(r[len("bytes="):-1])
            for r in server.plan["requests"]
        ]
        assert len(offsets) == 4
        assert offsets == sorted(offsets)
        assert offsets[0] == 0 and offsets[-1] > 0


# ── Retry, and giving up ─────────────────────────────────────────────────────

class TestRetry:

    def test_a_transient_status_is_retried(self, server, dest):
        server.plan["statuses"] = [503, 502]
        download_file(server.url, dest, progress=False, backoff=0)
        assert open(dest, "rb").read() == PAYLOAD

    def test_the_budget_is_finite_and_the_error_says_what_is_on_disk(
        self, server, dest
    ):
        server.plan["statuses"] = [503] * 10
        with pytest.raises(DownloadError, match="after 3 attempts"):
            download_file(server.url, dest, progress=False, backoff=0, max_retries=2)

    def test_a_fatal_status_is_not_retried(self, server, dest):
        server.plan["statuses"] = [404]
        with pytest.raises(DownloadError, match="HTTP 404"):
            download_file(server.url, dest, progress=False, backoff=0)
        assert server.plan["requests"] == [None]

    def test_nothing_lands_at_the_destination_when_the_download_fails(
        self, server, dest
    ):
        open(dest, "wb").write(b"the good file")
        server.plan["statuses"] = [404]
        with pytest.raises(DownloadError):
            download_file(server.url, dest, progress=False, backoff=0)
        assert open(dest, "rb").read() == b"the good file"


# ── Checksums ────────────────────────────────────────────────────────────────

class TestChecksums:

    def test_a_matching_checksum_passes(self, server, dest):
        download_file(server.url, dest, expected_md5=PAYLOAD_MD5, progress=False)
        assert os.path.exists(dest)

    def test_a_published_checksum_is_fetched_and_used(self, server, dest):
        download_file(server.url, dest, checksum_url=server.url + ".md5",
                      progress=False)
        assert open(dest, "rb").read() == PAYLOAD

    def test_read_checksum_takes_the_digest_from_a_coreutils_line(self, server):
        assert read_checksum(server.url + ".md5") == PAYLOAD_MD5

    def test_a_mismatch_is_refused(self, server, dest):
        server.plan["corrupt"] = True
        with pytest.raises(DownloadError, match="Checksum mismatch"):
            download_file(server.url, dest, expected_md5=PAYLOAD_MD5, progress=False)
        assert not os.path.exists(dest)

    def test_a_mismatch_deletes_the_partial_file(self, server, dest):
        """Resuming from bytes known to be wrong would only repeat the failure."""
        server.plan["corrupt"] = True
        with pytest.raises(DownloadError):
            download_file(server.url, dest, expected_md5=PAYLOAD_MD5, progress=False)
        assert not os.path.exists(dest + PART_SUFFIX)

    def test_an_explicit_digest_wins_over_a_published_one(self, server, dest):
        with pytest.raises(DownloadError, match="Checksum mismatch"):
            download_file(server.url, dest, expected_md5="0" * 32,
                          checksum_url=server.url + ".md5", progress=False)

    def test_md5_of_file_reads_in_chunks(self, tmp_path):
        path = tmp_path / "f"
        path.write_bytes(PAYLOAD)
        assert md5_of_file(str(path), chunk_size=17) == PAYLOAD_MD5


# ── The caller's own check ───────────────────────────────────────────────────

class TestVerifyCallback:

    def test_verify_sees_the_finished_file_before_it_moves(self, server, dest):
        seen = {}

        def check(path):
            seen["path"] = path
            seen["bytes"] = os.path.getsize(path)

        download_file(server.url, dest, verify=check, progress=False)
        assert seen["path"] == dest + PART_SUFFIX
        assert seen["bytes"] == len(PAYLOAD)

    def test_a_rejected_file_never_reaches_the_destination(self, server, dest):
        open(dest, "wb").write(b"the good file")

        def reject(path):
            raise ValueError("not a database")

        with pytest.raises(ValueError, match="not a database"):
            download_file(server.url, dest, verify=reject, progress=False)
        assert open(dest, "rb").read() == b"the good file"

    def test_a_rejected_file_is_deleted_rather_than_resumed(self, server, dest):
        """A complete file that fails its check would fail it again on resume."""

        def reject(path):
            raise ValueError("not a database")

        with pytest.raises(ValueError):
            download_file(server.url, dest, verify=reject, progress=False)
        assert not os.path.exists(dest + PART_SUFFIX)

    def test_the_callers_exception_type_is_preserved(self, server, dest):
        class MyError(Exception):
            pass

        def reject(path):
            raise MyError("mine")

        with pytest.raises(MyError):
            download_file(server.url, dest, verify=reject, progress=False)
