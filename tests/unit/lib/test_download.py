"""
Tests for `lib.download.fetch_to`, `write_through_part` and `download_file` against a loopback server: the
User-Agent, POST bodies, short reads, a server that stalls and a response failing its content check; and the
content checks themselves.
"""

import concurrent.futures
import gzip
import os
import socket
import threading
import time
from urllib.error import ContentTooShortError, URLError

import pytest

from pocketmapper.lib import download
from pocketmapper.lib.download import (
    USER_AGENT,
    DownloadJob,
    UnexpectedContentError,
    check_gzip,
    check_json,
    check_mmcif,
    download_file,
    fetch_to,
    is_transient_error,
    write_through_part,
)

BODY = b"0123456789"


def response(body, content_length):
    """An HTTP/1.1 200 response carrying `body` and declaring `content_length` bytes."""
    head = f"HTTP/1.1 200 OK\r\nContent-Length: {content_length}\r\nConnection: close\r\n\r\n"
    return head.encode() + body


class OneShotServer:
    """
    Loopback server that accepts one connection, records the request and sends `reply`, or nothing if None.
    """

    def __init__(self, reply):
        """Listen on an ephemeral port and serve one connection in a background thread."""
        self.reply = reply
        self.request = b""
        self.release = threading.Event()
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.bind(("127.0.0.1", 0))
        self.sock.listen(1)
        self.url = f"http://127.0.0.1:{self.sock.getsockname()[1]}/file"
        self.thread = threading.Thread(target=self.serve, daemon=True)
        self.thread.start()

    def serve(self):
        """Read the request headers and body, then send the reply, or hold the connection open until released."""
        conn, _ = self.sock.accept()
        with conn:
            conn.settimeout(5)
            while b"\r\n\r\n" not in self.request:
                self.request += conn.recv(65536)
            head, _, body = self.request.partition(b"\r\n\r\n")
            length = 0
            for line in head.split(b"\r\n"):
                if line.lower().startswith(b"content-length:"):
                    length = int(line.split(b":", 1)[1])
            while len(body) < length:
                body += conn.recv(65536)
            self.request = head + b"\r\n\r\n" + body
            if self.reply is None:
                self.release.wait(5)
            else:
                conn.sendall(self.reply)

    def close(self):
        """Release a held connection and stop listening."""
        self.release.set()
        self.thread.join(5)
        self.sock.close()


@pytest.fixture
def serve():
    """A function starting a `OneShotServer` with a given reply; every server is closed after the test."""
    servers = []

    def start(reply):
        """Start and return a server sending `reply`."""
        server = OneShotServer(reply)
        servers.append(server)
        return server

    yield start
    for server in servers:
        server.close()


@pytest.mark.parametrize(
    ("data", "method"),
    [(None, b"GET"), (b"1abc,2xyz", b"POST")],
)
def test_fetch_to_writes_body(serve, tmp_path, data, method):
    """A complete response lands in the file, and the request carries the User-Agent, method and body."""
    server = serve(response(BODY, len(BODY)))
    out = os.path.join(str(tmp_path), "out")
    fetch_to(server.url, out, data=data)
    server.thread.join(5)
    with open(out, "rb") as f:
        assert f.read() == BODY
    assert server.request.startswith(method + b" /file ")
    assert f"User-Agent: {USER_AGENT}".encode() in server.request
    if data:
        assert server.request.endswith(data)


@pytest.mark.parametrize(
    ("reply", "expected"),
    [
        (response(BODY, 100), ContentTooShortError),
        (None, (TimeoutError, URLError)),
    ],
    ids=["short_read", "stall"],
)
def test_fetch_to_failures_are_transient(serve, tmp_path, reply, expected):
    """A short body or a stalled server raises an error that `is_transient_error` retries."""
    server = serve(reply)
    with pytest.raises(expected) as info:
        fetch_to(server.url, os.path.join(str(tmp_path), "out"), timeout=0.2)
    assert is_transient_error(info.value)


def copy_reversed(src_fpath, dst_fpath):
    """A transform writing the bytes of `src_fpath` reversed to `dst_fpath`."""
    with open(src_fpath, "rb") as f_in, open(dst_fpath, "wb") as f_out:
        f_out.write(f_in.read()[::-1])


def failing_transform(src_fpath, dst_fpath):
    """A transform that writes part of its output and then raises."""
    with open(dst_fpath, "wb") as f_out:
        f_out.write(b"half")
    raise ValueError("transform failed")


@pytest.mark.parametrize(
    ("reply", "transform", "error", "expected"),
    [
        pytest.param(response(BODY, len(BODY)), None, None, BODY, id="success"),
        pytest.param(response(BODY, len(BODY)), copy_reversed, None, BODY[::-1], id="success_transform"),
        pytest.param(response(BODY, 100), None, ContentTooShortError, None, id="short_read"),
        pytest.param(response(BODY, 100), copy_reversed, ContentTooShortError, None, id="short_read_transform"),
        pytest.param(None, None, (TimeoutError, URLError), None, id="stall"),
        pytest.param(response(BODY, len(BODY)), failing_transform, ValueError, None, id="transform_fails"),
    ],
)
def test_write_through_part_leaves_no_scratch(serve, tmp_path, monkeypatch, reply, transform, error, expected):
    """Success leaves only the final file, failure leaves nothing: no scratch file survives either way."""

    def quick_fetch(url, raw_fpath, data=None):
        """`fetch_to` with a short timeout."""
        fetch_to(url, raw_fpath, data, timeout=0.2)

    monkeypatch.setattr(download, "fetch_to", quick_fetch)
    server = serve(reply)
    out = os.path.join(str(tmp_path), "sub", "out.json")
    if error is None:
        write_through_part(server.url, out, transform, None)
        assert os.listdir(os.path.dirname(out)) == ["out.json"]
        with open(out, "rb") as f:
            assert f.read() == expected
    else:
        with pytest.raises(error):
            write_through_part(server.url, out, transform, None)
        assert os.listdir(os.path.dirname(out)) == []


def test_download_many_interrupt_cancels_queued(tmp_path, monkeypatch):
    """A KeyboardInterrupt mid-run cancels the queued jobs instead of waiting for all of them."""
    lock = threading.Lock()
    calls = []

    def fetch(url, out_fpath, transform=None, check=None, log_extra=None):
        """Count the call and take a moment, as a download would."""
        with lock:
            calls.append(url)
        time.sleep(0.05)
        return True

    def interrupted(futures):
        """Yield the first future to finish, then raise as Ctrl-C would."""
        yield next(concurrent.futures.as_completed(futures))
        raise KeyboardInterrupt

    monkeypatch.setattr(download, "as_completed", interrupted)
    jobs = {i: DownloadJob(f"http://x/{i}", os.path.join(str(tmp_path), f"{i}.txt")) for i in range(50)}
    with pytest.raises(KeyboardInterrupt):
        download.download_many(jobs, fetch=fetch, max_workers=2)
    assert len(calls) < 10


HTML_BODY = b"<html><body>Down for maintenance</body></html>"
GOOD_CONTENT = {
    check_gzip: gzip.compress(b"data_x\n"),
    check_json: b'{"a": 1}',
    check_mmcif: b"# a comment\n\n#\ndata_x\n_cell.length_a 10\n",
}


@pytest.mark.parametrize(
    ("check", "content", "accepted"),
    [
        pytest.param(check, content, accepted, id=f"{check.__name__}_{name}")
        for check in GOOD_CONTENT
        for name, content, accepted in [
            ("good", GOOD_CONTENT[check], True),
            ("html", HTML_BODY, False),
            ("empty", b"", False),
        ]
    ],
)
def test_content_checks(tmp_path, check, content, accepted):
    """Each check accepts its kind of file and rejects an HTML page or an empty file, without a retry."""
    fpath = os.path.join(str(tmp_path), "response")
    with open(fpath, "wb") as f:
        f.write(content)
    if accepted:
        assert check(fpath) is None
    else:
        with pytest.raises(UnexpectedContentError) as info:
            check(fpath)
        assert fpath in str(info.value)
        assert not is_transient_error(info.value)


def test_download_file_unexpected_content(serve, tmp_path, monkeypatch):
    """A 200 HTML page failing `check_gzip` costs one request and leaves nothing in the output directory."""
    calls = []
    real_fetch_to = download.fetch_to

    def counting_fetch(url, raw_fpath, data=None):
        """Count the request, then run the real `fetch_to` with a short timeout."""
        calls.append(url)
        real_fetch_to(url, raw_fpath, data, timeout=0.2)

    monkeypatch.setattr(download, "fetch_to", counting_fetch)
    server = serve(response(HTML_BODY, len(HTML_BODY)))
    out_dir = os.path.join(str(tmp_path), "pdb")
    out = os.path.join(out_dir, "4Q5J.cif.gz")
    assert download_file(server.url, out, base_delay=0, check=check_gzip) is False
    assert calls == [server.url]
    assert os.listdir(out_dir) == []
