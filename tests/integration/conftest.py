"""
Shared fixtures for the integration tests: resource gating, the seeded cache, the network guard, `run_cli`
and the assertions on a run's outputs.

The assertion helpers are fixtures that return functions, as in `tests/unit/conftest.py`, so a test module
needs no import from this directory, whatever pytest's import mode.
"""

import csv
import os
import shutil
import socket

import pytest

from pocketmapper import cli
from pocketmapper.lib.entries import read_entries

TESTS_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIXTURES_DIR = os.path.join(TESTS_DIR, "integration", "fixtures")
SEED_DIR = os.path.join(TESTS_DIR, "data", "seed_cache")
NETWORK_CACHE_DIR = os.path.join(TESTS_DIR, "integration", "network_cache")

# The structure the seed shares with fixtures/, copied into the cache rather than committed twice
SEED_STRUCTURE = "4Q5J.cif.gz"

# Output files of a run, by their default names under results_dir
COMPARISON_NAME = "pocket_comparison.tsv"
ALIGNMENT_NAME = "alignment.tsv"


class NetworkBlocked(RuntimeError):
    """
    Raised in place of a network call. `lib.download.is_transient_error` returns False for it, so no retry waits.
    """


def pytest_runtest_setup(item):
    """Skip a test whose marker names a resource this machine lacks, saying which."""
    # Most specific resource first, so a test needing several is skipped for the rarest one
    if item.get_closest_marker("huge") and os.environ.get("POCKETMAPPER_TEST_HUGE") != "1":
        pytest.skip("downloads the full-PDB Foldseek database; set POCKETMAPPER_TEST_HUGE=1")
    if item.get_closest_marker("pdb_fsdb"):
        fsdb_path = os.environ.get("POCKETMAPPER_PDB_FSDB")
        if not fsdb_path or not os.path.exists(fsdb_path):
            pytest.skip("needs a PDB Foldseek database; set POCKETMAPPER_PDB_FSDB to an existing path")
    if item.get_closest_marker("foldseek") and shutil.which("foldseek") is None:
        pytest.skip("foldseek is not on PATH")
    if item.get_closest_marker("network") and os.environ.get("POCKETMAPPER_TEST_NETWORK") != "1":
        pytest.skip("needs the network; set POCKETMAPPER_TEST_NETWORK=1")


@pytest.fixture
def seed_cache(tmp_path):
    """A fresh copy of the committed seed cache, with 4Q5J.cif.gz in `pdb_structures/`, as a string path."""
    cache_dir = os.path.join(str(tmp_path), "cache")
    shutil.copytree(SEED_DIR, cache_dir)
    pdb_dir = os.path.join(cache_dir, "pdb_structures")
    os.makedirs(pdb_dir, exist_ok=True)
    shutil.copy2(os.path.join(FIXTURES_DIR, SEED_STRUCTURE), os.path.join(pdb_dir, SEED_STRUCTURE))
    return cache_dir


@pytest.fixture
def network_cache():
    """The cache shared by `network` tests: `POCKETMAPPER_TEST_CACHE`, else `tests/integration/network_cache/`."""
    cache_dir = os.path.abspath(os.environ.get("POCKETMAPPER_TEST_CACHE") or NETWORK_CACHE_DIR)
    os.makedirs(cache_dir, exist_ok=True)
    return cache_dir


@pytest.fixture(autouse=True)
def no_network(request, monkeypatch):
    """Block and record every connection and name lookup of a test not marked `network`; fail if any was tried."""
    if request.node.get_closest_marker("network"):
        yield
        return

    attempts = []

    def blocked(target):
        """Record an attempt on `target` and raise NetworkBlocked."""
        attempts.append(str(target))
        raise NetworkBlocked(f"network access blocked in an offline test: {target}")

    def connect(self, address):
        """Stand-in for `socket.socket.connect`."""
        blocked(address)

    def connect_ex(self, address):
        """Stand-in for `socket.socket.connect_ex`."""
        blocked(address)

    def getaddrinfo(host, port, *args, **kwargs):
        """Stand-in for `socket.getaddrinfo`, where a lookup starts before any connection."""
        blocked(f"{host}:{port}")

    def create_connection(address, *args, **kwargs):
        """Stand-in for `socket.create_connection`."""
        blocked(f"{address[0]}:{address[1]}")

    monkeypatch.setattr(socket.socket, "connect", connect)
    monkeypatch.setattr(socket.socket, "connect_ex", connect_ex)
    monkeypatch.setattr(socket, "getaddrinfo", getaddrinfo)
    monkeypatch.setattr(socket, "create_connection", create_connection)
    yield
    # Download code may catch NetworkBlocked and carry on, so the record, not the raise, is the check
    assert not attempts, (
        f"{len(attempts)} network attempt(s) in an offline test, so the seed cache lacks a file the run needs: "
        f"{', '.join(sorted(set(attempts)))}"
    )


@pytest.fixture
def run_cli(monkeypatch):
    """A function running `cli.cli` on its arguments from a working directory and returning the exit code."""

    def run(*argv, cwd=None):
        """Run `cli.cli(list(argv))` from `cwd` (default `fixtures/`); return 0 or the `SystemExit` code."""
        monkeypatch.chdir(cwd or FIXTURES_DIR)
        try:
            cli.cli(list(argv))
        except SystemExit as e:
            return 0 if e.code is None else e.code
        return 0

    return run


@pytest.fixture
def comparison_rows():
    """A function returning the data rows of a run's pocket comparison table."""

    def rows(results_dir, name=COMPARISON_NAME):
        """The rows of `name` under `results_dir`, each a dict keyed by column."""
        path = os.path.join(results_dir, name)
        assert os.path.isfile(path), f"no {path}"
        with open(path, newline="") as f:
            return list(csv.DictReader(f, delimiter="\t"))

    return rows


@pytest.fixture
def distinct_queries(comparison_rows):
    """A function returning the distinct values of the `query` column of a run's pocket comparison table."""

    def queries(results_dir, name=COMPARISON_NAME):
        """The set of `query` values in `name` under `results_dir`."""
        return {row["query"] for row in comparison_rows(results_dir, name)}

    return queries


@pytest.fixture
def failed():
    """A function telling whether an entries table marks an entry as failed for a given reason."""

    def marked(entries_path, reason, pocket_id):
        """True if a row of `entries_path` has `pocket_id` and `failed_reason` equal to `reason`."""
        return any(
            row["failed_reason"] == reason and row["pocket_id"] == pocket_id for row in read_entries(entries_path)
        )

    return marked


@pytest.fixture
def same_outputs():
    """A function telling whether two runs wrote the same comparison and alignment tables, ignoring row order."""

    def same(dir_a, dir_b):
        """True if `pocket_comparison.tsv` and `alignment.tsv` under `dir_a` and `dir_b` hold the same lines."""
        for name in (COMPARISON_NAME, ALIGNMENT_NAME):
            path_a = os.path.join(dir_a, name)
            path_b = os.path.join(dir_b, name)
            assert os.path.isfile(path_a), f"no {path_a}"
            assert os.path.isfile(path_b), f"no {path_b}"
            with open(path_a) as f_a, open(path_b) as f_b:
                if sorted(f_a) != sorted(f_b):
                    return False
        return True

    return same
