"""
Unit tests for unusable pocket residues: one the chain lacks or with no CA atom is dropped from its pocket with one
warning, and a pocket left with none of its named residues is not built. On 4Q5J chain A, where 9999 is not in the
chain and 1201 is its ANP ligand, which has no CA atom in the polymer.
"""

import json
import os

import pytest

from pocketmapper.lib.structure_manager import StructureManager
from pocketmapper.pockets.pisa import pisa_pockets
from pocketmapper.pockets.pisa_downloader import PisaDownloader
from pocketmapper.pockets.structure import parse_pocket_from_struct, passthrough_pockets

# The start of the one warning `parse_pocket_from_struct` logs per pocket for its dropped residues
DROPPED_WARNING = "Dropped "

# (entry, residue ids kept or None when the entry is left out, text the dropped-residue warning holds)
PASSTHROUGH_CASES = [
    pytest.param("4Q5J:A:1101,9999,1104", ["1101", "1104"], "not in the chain: 9999", id="not_in_chain"),
    pytest.param("4Q5J:A:1101,1201,1104", ["1101", "1104"], "no CA atom in the polymer: 1201", id="no_ca"),
    pytest.param("4Q5J:A:9999", None, "not in the chain: 9999", id="only_not_in_chain"),
    pytest.param("4Q5J:A:1201,9999", None, "not in the chain: 9999; no CA atom in the polymer: 1201", id="none_left"),
]

# A residue of chain A that the seeded A,E interface lists both as bonded and as buried at threshold 0
RENAMED_RESIDUE = "1139"


def parsed_record(processor, entry):
    """The record the parser builds for `entry`, which must parse."""
    ((_, record, reason),) = processor.process_qt_cmdline_input(entry, "query")
    assert reason is None
    return record


def dropped_warnings(warnings):
    """The captured dropped-residue warnings."""
    return [message for message in warnings() if message.startswith(DROPPED_WARNING)]


def assert_all_positioned(pocket):
    """Every listed residue of `pocket` is in `residues` with a usable `seq_pos`."""
    assert set(pocket.res_auth_ids) == set(pocket.residues)
    assert all(isinstance(residue.seq_pos, int) and residue.seq_pos >= 0 for residue in pocket.residues.values())


@pytest.mark.parametrize("entry, kept, reason", PASSTHROUGH_CASES)
def test_passthrough_drops_unusable(network_attempts, processor, pdb_dir, warnings, entry, kept, reason):
    """Unusable residues are dropped with one warning; an entry left with none is not built."""
    pockets = passthrough_pockets([parsed_record(processor, entry)], StructureManager(pdb_dir, download=False))

    dropped = dropped_warnings(warnings)
    assert len(dropped) == 1
    assert reason in dropped[0]
    if kept is None:
        assert pockets == {}
        return
    assert pockets[entry].res_auth_ids == kept
    assert_all_positioned(pockets[entry])


def test_empty_residue_list_kept(pdb_dir, warnings):
    """An empty residue list gives a pocket with `pocket_exists` False, not None, and drops nothing."""
    pocket = parse_pocket_from_struct(os.path.join(pdb_dir, "4Q5J.cif.gz"), "A", [])

    assert pocket is not None
    assert not pocket.pocket_exists
    assert pocket.res_auth_ids == []
    assert dropped_warnings(warnings) == []


@pytest.fixture
def renamed_pisa_dir(pocket_dir):
    """The seeded PISA cache with chain A's RENAMED_RESIDUE renumbered to 9999 in the A,E interface."""
    pisa_dir = os.path.join(pocket_dir, "pisa")
    path = os.path.join(pisa_dir, "interface_pairs", "4q5j.json")
    with open(path) as f:
        entry = json.load(f)
    interface = entry["interfaces"]["A,E"]
    molecule = next(mol for mol in interface["molecules"] if mol["chain_id"] == "A")
    seq_id_lists = [molecule["residue_seq_ids"]]
    for value in interface.values():
        if isinstance(value, dict):
            seq_id_lists.append(value[f"atom_site_{molecule['molecule_id']}_seq_nums"])
    for seq_ids in seq_id_lists:
        seq_ids[:] = ["9999" if seq_id == RENAMED_RESIDUE else seq_id for seq_id in seq_ids]
    with open(path, "w") as f:
        json.dump(entry, f)
    return pisa_dir


@pytest.mark.parametrize("bsa_threshold", [pytest.param(None, id="bonded"), pytest.param(0, id="buried")])
def test_pisa_drops_residue_not_in_chain(no_downloads, processor, pdb_dir, renamed_pisa_dir, warnings, bsa_threshold):
    """An interface residue the chain lacks is dropped with one warning, and the rest keep their seq_pos."""
    record = parsed_record(processor, "4Q5J:A_E")
    pockets = pisa_pockets(
        [record],
        StructureManager(pdb_dir, download=False),
        PisaDownloader(renamed_pisa_dir, download=False),
        bsa_threshold=bsa_threshold,
    )

    pocket = pockets["4Q5J:A_E"]
    assert pocket.pocket_exists
    assert "9999" not in pocket.res_auth_ids
    assert len(pocket.res_auth_ids) > 0
    assert_all_positioned(pocket)
    dropped = dropped_warnings(warnings)
    assert len(dropped) == 1
    assert "not in the chain: 9999" in dropped[0]


def test_pisa_empty_under_rule_kept(no_downloads, processor, pdb_dir, pocket_dir, warnings):
    """A pisa pocket with no residue past the threshold is kept with `pocket_exists` False."""
    record = parsed_record(processor, "4Q5J:A_E")
    pockets = pisa_pockets(
        [record],
        StructureManager(pdb_dir, download=False),
        PisaDownloader(os.path.join(pocket_dir, "pisa"), download=False),
        bsa_threshold=100,
    )

    assert list(pockets) == ["4Q5J:A_E"]
    assert not pockets["4Q5J:A_E"].pocket_exists
    assert pockets["4Q5J:A_E"].res_auth_ids == []
    assert dropped_warnings(warnings) == []
