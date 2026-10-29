"""
Unit tests for `pocketmapper.pockets.pisa`: PISA pockets by bond and by buried surface, on the seeded 4Q5J
interfaces, and the PISA downloader with downloads off.
"""

import os

import pytest

from pocketmapper.pockets.pisa import pisa_pockets, pisa_pockets_from_records
from pocketmapper.pockets.pisa_downloader import PisaDownloader

# The start of the one warning `download_missing_interfaces` logs with downloads off
SKIPPED_WARNING = "Downloads are off"

# (pdb_list, cached files removed first, count of entries skipped). 4Q5J is the seeded entry; 1ABC and
# 2XYZ are not cached. With only its interface file removed, 4Q5J still has its summary and assemblies.
OFFLINE_CASES = [
    pytest.param(["4Q5J"], [], 0, id="all_cached"),
    pytest.param(["4Q5J", "1abc"], [], 1, id="mixed"),
    pytest.param(["1ABC", "2XYZ", "1abc"], [], 2, id="none_cached"),
    pytest.param(["4Q5J"], [os.path.join("interface_pairs", "4q5j.json")], 1, id="interfaces_missing"),
]

# The threshold the e2e case used; it removes residues from both seeded pockets
BSA_THRESHOLD = 25

# (pocket_id, residues by bond, buried at 0, buried past BSA_THRESHOLD), counted by hand from the seed
POCKET_SIZES = [
    pytest.param("4Q5J:A_E", 11, 19, 14, id="4Q5J:A_E"),
    pytest.param("4Q5J:B_F", 7, 18, 14, id="4Q5J:B_F"),
]


@pytest.fixture
def interface_dir(pocket_dir):
    """The seeded PISA interface files, in the test's own pocket directory."""
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


def tree(root):
    """Every file under `root`, relative to it."""
    return {os.path.relpath(os.path.join(d, f), root) for d, _, files in os.walk(root) for f in files}


@pytest.mark.parametrize("pdb_list, removed, skipped", OFFLINE_CASES)
def test_downloader_offline(no_downloads, warnings, pocket_dir, pdb_list, removed, skipped):
    """With downloads off, entries not in `interface_pairs/` are skipped with one warning, and nothing is written."""
    pisa_dir = os.path.join(pocket_dir, "pisa")
    for path in removed:
        os.remove(os.path.join(pisa_dir, path))
    before = tree(pisa_dir)

    failures = PisaDownloader(pisa_dir, download=False).download_missing_interfaces(pdb_list)

    assert failures == {}
    assert tree(pisa_dir) == before
    assert not os.path.exists(os.path.join(pisa_dir, "errors.json"))
    skipped_warnings = [message for message in warnings() if message.startswith(SKIPPED_WARNING)]
    if skipped:
        assert len(skipped_warnings) == 1
        assert f" {skipped} " in skipped_warnings[0]
    else:
        assert skipped_warnings == []


def test_pisa_pockets_offline(no_downloads, processor, pdb_dir, pocket_dir):
    """With downloads off, a cached entry's pocket is built and an uncached entry's is left out."""
    records = []
    for entry in ("4Q5J:A_E", "1ABC:A_B"):
        ((_, record, reason),) = processor.process_qt_cmdline_input(entry, "query")
        assert reason is None
        records.append(record)

    pockets = pisa_pockets(records, PisaDownloader(os.path.join(pocket_dir, "pisa"), download=False))

    assert list(pockets) == ["4Q5J:A_E"]
    assert pockets["4Q5J:A_E"].pocket_exists
    assert not os.path.exists(os.path.join(pocket_dir, "pisa", "errors.json"))
