"""
Unit tests for `pocketmapper.pockets.pisa`: PISA pockets by bond and by buried surface, on the seeded 4Q5J
interfaces.
"""

import os
import shutil
import socket

import pytest

from pocketmapper.pockets.pisa import pisa_pockets_from_records

SEED_POCKET_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "data", "seed_cache", "pockets")

# The threshold the e2e case used; it removes residues from both seeded pockets
BSA_THRESHOLD = 25

# (pocket_id, residues by bond, buried at 0, buried past BSA_THRESHOLD), counted by hand from the seed
POCKET_SIZES = [
    pytest.param("4Q5J:A_E", 11, 19, 14, id="4Q5J:A_E"),
    pytest.param("4Q5J:B_F", 7, 18, 14, id="4Q5J:B_F"),
]


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
def interface_dir(tmp_path):
    """A copy of the seeded PISA interface files, in the test's own pocket directory."""
    pocket_dir = os.path.join(str(tmp_path), "pockets")
    shutil.copytree(SEED_POCKET_DIR, pocket_dir)
    return os.path.join(pocket_dir, "pisa", "interface_pairs")


@pytest.fixture
def records(processor):
    """The 4Q5J query records, built by the parser."""
    records = []
    for entry in ("4Q5J:A_E", "4Q5J:B_F"):
        ((_, record, reason),) = processor.process_qt_cmdline_input(entry, "query")
        assert reason is None
        records.append(record)
    return records


def residues(records, interface_dir, bsa_threshold=None):
    """pocket_id -> the set of residue ids `pisa_pockets_from_records` gives under `bsa_threshold`."""
    pockets = pisa_pockets_from_records(records, interface_dir, bsa_threshold=bsa_threshold)
    return {pocket_id: set(pocket.res_auth_ids) for pocket_id, pocket in pockets.items()}


@pytest.mark.parametrize("pocket_id, bonded, buried, past_threshold", POCKET_SIZES)
def test_bsa_threshold_narrows_pocket(
    network_attempts, records, interface_dir, pocket_id, bonded, buried, past_threshold
):
    """A threshold keeps a strict subset of the buried residues, and with 0 every bonded residue is buried."""
    by_bond = residues(records, interface_dir)[pocket_id]
    at_zero = residues(records, interface_dir, 0)[pocket_id]
    narrowed = residues(records, interface_dir, BSA_THRESHOLD)[pocket_id]
    assert (len(by_bond), len(at_zero), len(narrowed)) == (bonded, buried, past_threshold)
    assert narrowed < at_zero
    assert by_bond <= at_zero
