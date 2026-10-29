"""
Shared fixtures for the pockets unit tests: the seeded PISA cache, a PDB cache holding 4Q5J, and guards
against the network and against downloads.
"""

import logging
import os
import shutil
import socket

import pytest

from pocketmapper.lib.constants import PACKAGE_LOGGER
from pocketmapper.pockets import pisa_downloader

TESTS_DIR = os.path.join(os.path.dirname(__file__), "..", "..")
SEED_POCKET_DIR = os.path.join(TESTS_DIR, "data", "seed_cache", "pockets")
STRUCTURE_4Q5J = os.path.join(TESTS_DIR, "integration", "fixtures", "4Q5J.cif.gz")


class NetworkCalled(RuntimeError):
    """Raised in place of a network connection."""


@pytest.fixture
def network_attempts(monkeypatch):
    """Block and record every connection and name lookup, and fail the test if any was attempted."""
    attempts = []

    def blocked(*args, **kwargs):
        """Record the attempt and raise in place of the network call."""
        attempts.append(args)
        raise NetworkCalled(f"network call attempted: {args}")

    monkeypatch.setattr(socket.socket, "connect", blocked)
    monkeypatch.setattr(socket, "getaddrinfo", blocked)
    yield attempts
    assert attempts == []


@pytest.fixture
def no_downloads(monkeypatch, network_attempts):
    """Make every download function `pisa_downloader` reaches record the call and raise, and fail the test on any."""
    calls = []

    def refuse(*args, **kwargs):
        """Record the call and raise in place of the download."""
        calls.append(args)
        raise AssertionError(f"download attempted: {args}")

    for name in ("download_api", "download_file", "download_many"):
        monkeypatch.setattr(pisa_downloader, name, refuse)
    # The source table holds its own references to the fetch functions
    sources = {source: (url_for, refuse) for source, (url_for, _) in pisa_downloader.ASSEMBLY_SOURCES.items()}
    monkeypatch.setattr(pisa_downloader, "ASSEMBLY_SOURCES", sources)
    yield calls
    assert calls == []


@pytest.fixture
def pocket_dir(tmp_path):
    """A copy of the seeded pocket cache, with its PISA files, as the test's own pocket directory."""
    pocket_dir = os.path.join(str(tmp_path), "pockets")
    shutil.copytree(SEED_POCKET_DIR, pocket_dir)
    return pocket_dir


@pytest.fixture
def pdb_dir(work_dir):
    """`work_dir/pdb`, where the `processor` fixture's PDB records point, holding the 4Q5J structure."""
    pdb_dir = os.path.join(work_dir, "pdb")
    os.makedirs(pdb_dir)
    shutil.copy(STRUCTURE_4Q5J, os.path.join(pdb_dir, "4Q5J.cif.gz"))
    return pdb_dir


@pytest.fixture
def warnings(caplog):
    """A function returning the messages of the package's warnings captured so far."""
    caplog.set_level(logging.WARNING, logger=PACKAGE_LOGGER)

    def messages():
        """The messages of the captured WARNING records logged under the package logger."""
        return [
            record.getMessage()
            for record in caplog.records
            if record.levelno == logging.WARNING and record.name.startswith(PACKAGE_LOGGER)
        ]

    return messages
