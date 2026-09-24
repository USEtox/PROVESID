"""
Download a large dataset file so that an interruption costs nothing.

Every bulk dataset PROVESID reads -- ChEMBL's 5.8 GB archive, PubChem's 2.2 GB
identifier database, CompTox, ZeroPM, the ChEBI SDF -- now goes through
``provesid.datasets.download_file``.  It streams into a ``.part`` file beside
the destination, resumes that file with an HTTP ``Range`` request if the
transfer is interrupted, verifies an MD5 when the server publishes one, checks
the byte count against what the server declared, hands the finished file to
your own validity check, and only then moves it into place.

The practical consequence: losing a connection at 4 GB of a 5.8 GB download
costs you the seconds since the last chunk, not the four gigabytes.

This script demonstrates it against a local server it starts itself, so it
needs no network and downloads nothing large.  The last section shows the calls
you would make against the real thing.

Run with::

    python examples/datasets/resumable_download_demo.py
"""

import hashlib
import http.server
import logging
import os
import sqlite3
import tempfile
import threading

from provesid.datasets import DownloadError, download_file, md5_of_file

# download_file reports through the standard logging module, like the rest of
# PROVESID -- turn it on to watch the resume happen.
logging.basicConfig(level=logging.INFO, format="%(levelname)-7s %(message)s")

PAYLOAD = bytes(range(256)) * 4000          # 1 024 000 bytes
PAYLOAD_MD5 = hashlib.md5(PAYLOAD).hexdigest()


# ── A server that drops the connection the first time ────────────────────────

class FlakyHandler(http.server.BaseHTTPRequestHandler):
    """Serve PAYLOAD, hanging up mid-response on the first request only."""

    already_failed = False

    def log_message(self, *args):
        pass

    def do_GET(self):
        if self.path.endswith(".md5"):
            body = f"{PAYLOAD_MD5}  payload.bin\n".encode()
            self.send_response(200)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return

        # A Range header is how download_file says "continue from here".
        requested = self.headers.get("Range")
        start = int(requested.split("=")[1].split("-")[0]) if requested else 0
        print(f"    server: request for bytes {start}- "
              f"({'resumed' if start else 'from the beginning'})")

        body = PAYLOAD[start:]
        self.send_response(206 if start else 200)
        if start:
            self.send_header("Content-Range",
                             f"bytes {start}-{len(PAYLOAD) - 1}/{len(PAYLOAD)}")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()

        if not FlakyHandler.already_failed:
            # Send a third of it, then hang up -- a dropped connection.
            FlakyHandler.already_failed = True
            self.wfile.write(body[: len(body) // 3])
            self.wfile.flush()
            self.close_connection = True
            return
        self.wfile.write(body)


def main():
    httpd = http.server.HTTPServer(("127.0.0.1", 0), FlakyHandler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    url = f"http://127.0.0.1:{httpd.server_port}/payload.bin"

    work = tempfile.mkdtemp(prefix="provesid-download-demo-")
    dest = os.path.join(work, "payload.bin")

    # ── 1. A download that survives a dropped connection ─────────────────────
    # chunk_size is small here only so the demo is visible: iter_content
    # yields whole chunks, so with the 1 MB default the whole payload would
    # arrive as one chunk and the interrupted attempt would have written
    # nothing to resume from.  Against a 5.8 GB archive the default is right.
    print("\n1. Downloading from a server that hangs up once:")
    download_file(url, dest, progress=False, backoff=0, chunk_size=8192)
    print(f"    got {os.path.getsize(dest)} bytes, intact: "
          f"{open(dest, 'rb').read() == PAYLOAD}")
    print(f"    no .part file left behind: {not os.path.exists(dest + '.part')}")

    # ── 2. Resuming a partial file from a previous run ───────────────────────
    # This is what an interrupted session leaves on disk: the bytes, and a
    # marker naming the URL they came from.  download_file resumes a partial
    # only when that URL matches, so a .part left by some other download is
    # discarded rather than spliced onto this one.
    print("\n2. Resuming a .part file left by an earlier run:")
    os.remove(dest)
    with open(dest + ".part", "wb") as handle:
        handle.write(PAYLOAD[:700_000])
    with open(dest + ".part.source", "w") as handle:
        handle.write(url)
    download_file(url, dest, progress=False)
    print(f"    intact: {open(dest, 'rb').read() == PAYLOAD}")

    # The same partial, but claiming a different origin, is thrown away.
    print("\n2b. A .part file left by a *different* download:")
    os.remove(dest)
    with open(dest + ".part", "wb") as handle:
        handle.write(b"bytes from some other dataset")
    with open(dest + ".part.source", "w") as handle:
        handle.write("http://example.invalid/something-else.bin")
    download_file(url, dest, progress=False)
    print(f"    intact: {open(dest, 'rb').read() == PAYLOAD}")

    # ── 3. A published checksum ──────────────────────────────────────────────
    # PubChem's FTP mirror publishes an .md5 beside every file.  Point
    # checksum_url at it and the digest is fetched and enforced.
    print("\n3. Verifying against the checksum the server publishes:")
    os.remove(dest)
    download_file(url, dest, checksum_url=url + ".md5", progress=False)
    print(f"    md5 matches: {md5_of_file(dest) == PAYLOAD_MD5}")

    # ── 4. A damaged file never replaces a good one ──────────────────────────
    # verify= is handed the finished file *before* it is moved into place, so
    # rejecting it leaves whatever was already at the destination alone.
    print("\n4. Refusing a file that fails your own check:")

    def must_be_a_database(path):
        sqlite3.connect(path).execute("SELECT 1 FROM compounds LIMIT 1")

    try:
        download_file(url, dest, verify=must_be_a_database, progress=False)
    except sqlite3.DatabaseError as exc:
        print(f"    rejected: {exc}")
    print(f"    the previous good file is untouched: "
          f"{open(dest, 'rb').read() == PAYLOAD}")

    # A rejected file is deleted rather than kept: it downloaded completely,
    # so resuming it would finish instantly and fail the same check again.
    print(f"    the rejected download was not kept: "
          f"{not os.path.exists(dest + '.part')}")

    # A wrong checksum is refused the same way, for the same reason.
    try:
        download_file(url, dest, expected_md5="0" * 32, progress=False)
    except DownloadError as exc:
        print(f"    checksum mismatch reported: {str(exc)[:60]}...")

    httpd.shutdown()

    # ── 5. Against the real datasets ─────────────────────────────────────────
    print("\n5. The real thing -- these are the calls the clients now make:")
    print("""
    from provesid import CheMBL, PubChemID

    CheMBL()                    # 5.8 GB archive, resumable
    PubChemID.download_database()   # 2.2 GB, resumable, checked before it lands

    # Or directly, for a file PubChem publishes a checksum beside:
    from provesid.datasets import download_file
    base = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/"
    download_file(base + "CID-SMILES.gz", "CID-SMILES.gz",
                  checksum_url=base + "CID-SMILES.gz.md5")
    """)
    print(f"(demo files are in {work})")


if __name__ == "__main__":
    main()
