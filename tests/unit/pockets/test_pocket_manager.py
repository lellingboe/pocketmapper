"""
Unit tests for `pocketmapper.pockets.pocket_manager`: building pockets, expanding Foldseek PDB-database hits
into pisa records, and the pocket files. Offline, on the seeded 4Q5J PISA data and the `pdb_dir` fixture.
"""

import os

import pandas as pd
import pytest

from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.pockets.pocket import Pocket, PocketResidue
from pocketmapper.pockets.pocket_manager import (
    POCKETS_FILE_VERSION,
    PocketManager,
    read_pockets_file,
    write_pockets_file,
    write_pockets_table,
)

# (entry, pocket method the parser gives it): one entry per pocket method, all on 4Q5J
METHOD_CASES = [
    pytest.param("4Q5J:A_E", "pisa", id="pisa"),
    pytest.param("4Q5J:A:1101,1104", "passthrough", id="passthrough"),
    pytest.param("4Q5J:B_F:1", "vdw", id="vdw"),
    pytest.param("4Q5J:A", "whole_chain", id="whole_chain"),
]

# (bsa_threshold, residues in the 4Q5J:A_E pocket), as test_pisa's POCKET_SIZES counts them from the seed
BSA_CASES = [
    pytest.param(None, 11, id="by_bond"),
    pytest.param(0, 19, id="buried_at_0"),
    pytest.param(25, 14, id="buried_past_25"),
]

# The seed has 4Q5J's interfaces A,E and B,F. 1ABC is PDB-style but not cached; the last is not PDB-style.
HIT_NAMES = ["4q5j-assembly1_A", "4q5j-assembly1_B", "4q5j-assembly2_B", "1abc-assembly1_A", "AF-P06493-F1-model_v4_A"]

# (pocket_id, preprocess_name) of every record `expand_pdb_hits` gives for HIT_NAMES
EXPECTED_HIT_RECORDS = [
    ("4Q5J:A_E", "4q5j-assembly1_A"),
    ("4Q5J:B_F", "4q5j-assembly1_B"),
    ("4Q5J:B_F", "4q5j-assembly2_B"),
]

# hit name -> pocket_ids, as `expand_pdb_hits` maps HIT_NAMES
EXPECTED_HIT_CHAINS = {
    "4q5j-assembly1_A": ["4Q5J:A_E"],
    "4q5j-assembly1_B": ["4Q5J:B_F"],
    "4q5j-assembly2_B": ["4Q5J:B_F"],
    "1abc-assembly1_A": [],
    "AF-P06493-F1-model_v4_A": [],
}

# The directory left out of the manager; `expand_pdb_hits` needs all three
MISSING_DIR_CASES = [pytest.param(name, id=name) for name in ("pdb_dir", "alphafold_dir", "fsdb_dir")]

# (file content, id): pockets files `read_pockets_file` must reject
BAD_POCKETS_FILES = [
    pytest.param(f'{{"version": {POCKETS_FILE_VERSION - 1}, "pockets": {{}}, "chains": {{}}}}', id="wrong_version"),
    pytest.param("not json at all", id="not_json"),
    pytest.param(f'{{"version": {POCKETS_FILE_VERSION}, "chains": {{}}}}', id="missing_pockets"),
    pytest.param(f'{{"version": {POCKETS_FILE_VERSION}, "pockets": {{}}}}', id="missing_chains"),
]


def parse(processor, entry):
    """The QTRecord dict the parser gives for `entry`, which must parse."""
    ((_, record, reason),) = processor.process_qt_cmdline_input(entry, "query")
    assert reason is None
    return record


@pytest.fixture
def make_manager(pocket_dir, cache_dirs):
    """A factory of PocketManagers over the test's pocket directory and cache directories, offline by default."""

    def make(**kwargs):
        """A PocketManager with `kwargs` overriding the test's directories and `download=False`."""
        return PocketManager(pocket_dir, **(cache_dirs | {"download": False} | kwargs))

    return make


@pytest.fixture
def hit_records(no_downloads, pdb_dir, make_manager):
    """The records `expand_pdb_hits` gives for HIT_NAMES, with downloads off."""
    records, _, _ = make_manager().expand_pdb_hits(HIT_NAMES)
    return records


@pytest.mark.parametrize("entry, pocket_method", METHOD_CASES)
def test_build_each_method(no_downloads, processor, pdb_dir, pocket_dir, entry, pocket_method):
    """Each method's pocket is built with coordinates, and `pocket_dir` holds nothing but `pisa/` afterwards."""
    record = parse(processor, entry)
    assert record["pocket_method"] == pocket_method

    pockets = PocketManager(pocket_dir).build([record])

    assert list(pockets) == [entry]
    assert pockets[entry].has_coords
    assert os.listdir(pocket_dir) == ["pisa"]


@pytest.mark.parametrize("bsa_threshold, size", BSA_CASES)
def test_build_bsa_threshold(no_downloads, processor, pdb_dir, pocket_dir, bsa_threshold, size):
    """The manager's `bsa_threshold` reaches the pisa builder."""
    pockets = PocketManager(pocket_dir, bsa_threshold=bsa_threshold).build([parse(processor, "4Q5J:A_E")])

    assert len(pockets["4Q5J:A_E"].res_auth_ids) == size


def test_build_offline_uncached(no_downloads, network_attempts, processor, pdb_dir, pocket_dir):
    """With downloads off, a pisa entry whose PISA data is not cached is left out, and nothing is requested."""
    records = [parse(processor, entry) for entry in ("4Q5J:A_E", "1ABC:A_B")]

    pockets = PocketManager(pocket_dir, download=False).build(records)

    assert list(pockets) == ["4Q5J:A_E"]
    assert network_attempts == []


def test_expand_pdb_hits(no_downloads, network_attempts, pdb_dir, make_manager):
    """Each cached hit gives one record per interface, keyed by its hit name, and every hit is a `chains` key."""
    records, chains, failures = make_manager().expand_pdb_hits(HIT_NAMES)

    assert sorted((r["pocket_id"], r["preprocess_name"]) for r in records) == EXPECTED_HIT_RECORDS
    assert all(r["pocket_method"] == "pisa" for r in records)
    assert chains == EXPECTED_HIT_CHAINS
    assert failures == []
    assert network_attempts == []


def test_expand_pdb_hits_no_structures(no_downloads, local_root, make_manager):
    """With no hit structure available, `expand_pdb_hits` raises."""
    empty_pdb_dir = os.path.join(local_root, "empty_pdb")
    os.makedirs(empty_pdb_dir)

    with pytest.raises(PocketMapperError):
        make_manager(pdb_dir=empty_pdb_dir).expand_pdb_hits(HIT_NAMES)


@pytest.mark.parametrize("missing", MISSING_DIR_CASES)
def test_expand_pdb_hits_missing_dir(no_downloads, processor, pdb_dir, pocket_dir, make_manager, missing):
    """A manager without a directory builds pockets, but `expand_pdb_hits` raises before downloading anything."""
    pisa_files = sorted(os.walk(pocket_dir))
    manager = make_manager(download=True, **{missing: None})

    assert list(manager.build([parse(processor, "4Q5J:A_E")])) == ["4Q5J:A_E"]
    with pytest.raises(PocketMapperError, match=missing):
        manager.expand_pdb_hits(HIT_NAMES)
    assert sorted(os.walk(pocket_dir)) == pisa_files


def test_pockets_file_round_trip(local_root):
    """A written pockets file reads back equal, keeping `res_auth_ids` order, `seq_pos`, None fields and None pockets."""
    residues = {
        "25": PocketResidue(res_code="GLY", res_code_single="G", seq_pos=20, ca_coords=[1.0, 2.0, 3.0]),
        "3": PocketResidue(res_code="ALA", res_code_single="A", seq_pos=-1, uniprot_pos=7),
        "10": PocketResidue(res_code_single="k"),
    }
    pockets = {
        "4Q5J:A_E": Pocket(
            res_auth_ids=[3, 10, 25], residues=residues, ca_sequence="AKG", pocket_exists=True, has_coords=True
        ),
        "1ABC:A_B": None,
    }
    chains = {"4Q5J_A_hash": ["4Q5J:A_E"], "1ABC_A_hash": ["1ABC:A_B"], "empty_chain": []}
    path = os.path.join(local_root, "out", "pockets.json")

    write_pockets_file(pockets, chains, path)
    read_pockets, read_chains = read_pockets_file(path)

    assert read_pockets == pockets
    assert read_pockets["4Q5J:A_E"].res_auth_ids == [3, 10, 25]
    assert list(read_pockets["4Q5J:A_E"].residues) == ["25", "3", "10"]
    assert read_chains == chains


@pytest.mark.parametrize("content", BAD_POCKETS_FILES)
def test_read_pockets_file_rejects(write_file, content):
    """A pockets file of the wrong version, not JSON or missing a key raises PocketMapperError."""
    path = write_file("pockets.json", content)

    with pytest.raises(PocketMapperError):
        read_pockets_file(path)


@pytest.mark.parametrize("built", [pytest.param(True, id="built"), pytest.param(False, id="not_built")])
def test_write_pockets_table_once_per_pocket(hit_records, pocket_dir, local_root, built):
    """A pocket_id given by two hit records is written once, with its residues if built and empty if not."""
    records = [record for record in hit_records if record["pocket_id"] == "4Q5J:B_F"]
    assert len(records) == 2
    pockets = PocketManager(pocket_dir, download=False).build(records) if built else {"4Q5J:B_F": None}
    path = os.path.join(local_root, "out", "pockets.tsv")

    write_pockets_table(records, pockets, path)
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

    res_auth_ids = ",".join(map(str, pockets["4Q5J:B_F"].res_auth_ids)) if built else ""
    assert table.to_dict("records") == [
        {"pocket_id": "4Q5J:B_F", "chain": "B", "res_auth_ids": res_auth_ids, "method": "pisa"}
    ]
