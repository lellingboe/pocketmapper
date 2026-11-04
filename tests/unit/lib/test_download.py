"""
Tests for `lib.download.fetch_to` and `write_through_part` against a loopback server: the User-Agent, POST
bodies, short reads and a server that stalls.
"""

import os
import socket
import threading
from urllib.error import ContentTooShortError, URLError

import pytest

from pocketmapper.lib import download
from pocketmapper.lib.download import USER_AGENT, fetch_to, is_transient_error, write_through_part

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


def test_write_through_part_stall_leaves_nothing(serve, tmp_path, monkeypatch):
    """A stalled fetch leaves neither the final file nor its `.part` behind."""

    def quick_fetch(url, raw_fpath, data=None):
        """`fetch_to` with a short timeout."""
        fetch_to(url, raw_fpath, data, timeout=0.2)

    monkeypatch.setattr(download, "fetch_to", quick_fetch)
    server = serve(None)
    out = os.path.join(str(tmp_path), "sub", "out.json")
    with pytest.raises((TimeoutError, URLError)):
        write_through_part(server.url, out, None, None)
    assert not os.path.exists(out)
    assert not os.path.exists(f"{out}.part")
