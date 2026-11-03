"""
Unit tests for residue keys: a residue is its number and insertion code, so 100, 100A and 100B are three residues.

The synthetic cases renumber 4Q5J residues in memory to give them insertion codes, or build a small PISA interface
dict. The 4C59 cases run the builders on a real entry: its chain B numbers 111A-111E with no plain 111
(`tests/data/insertion_codes/`).
"""

import os
import shutil

import gemmi
import pytest

from pocketmapper.lib.residues import ca_residues, residue_key, residue_sort_key
from pocketmapper.lib.structure_manager import StructureManager, load_structure
from pocketmapper.pockets.pisa import bonded_residues, buried_residues, pisa_pockets
from pocketmapper.pockets.pisa_downloader import PisaDownloader
from pocketmapper.pockets.structure import parse_pocket_from_struct, passthrough_pockets, whole_chain_pockets
from pocketmapper.pockets.vdw import vdw_contact_residues

TESTS_DIR = os.path.join(os.path.dirname(__file__), "..", "..")
STRUCTURE_4Q5J = os.path.join(TESTS_DIR, "integration", "fixtures", "4Q5J.cif.gz")
INSERTION_CODES_DIR = os.path.join(TESTS_DIR, "data", "insertion_codes")
STRUCTURE_4C59 = os.path.join(INSERTION_CODES_DIR, "4C59.cif.gz")
INTERFACES_4C59 = os.path.join(INSERTION_CODES_DIR, "4c59.json")

# (chain, old number) -> new SeqId. 855 and 856 follow 854 in chain A; 1105 follows 1104, and both contact chain E.
RENUMBERED = {("A", 855): (854, "A"), ("A", 856): (854, "B"), ("A", 1105): (1104, "A")}

BOND_TYPES = ["hydrogen_bonds", "salt_bridges", "disulfide_bonds", "covalent_bonds", "other_bonds"]


@pytest.fixture
def renumbered():
    """4Q5J with chain A's 855, 856 and 1105 renumbered to 854A, 854B and 1104A."""
    st = load_structure(STRUCTURE_4Q5J, {})
    for chain in st[0]:
        for res in chain:
            new = RENUMBERED.get((chain.name, res.seqid.num))
            if new is not None:
                res.seqid = gemmi.SeqId(*new)
    return st


def test_residue_key_forms():
    """gemmi's space and PISA's None both mean no insertion code; a code keeps its case."""
    assert residue_key(100, " ") == residue_key(100, None) == residue_key(100) == "100"
    assert residue_key(100, "A") == "100A"
    assert residue_key(100, "a") == "100a"
    assert residue_key(-3, None) == "-3"


def test_residue_sort_key_order():
    """Keys sort by number, then insertion code, negative numbers included."""
    keys = ["101", "100B", "-3", "100", "2", "100A", "-10"]
    assert sorted(keys, key=residue_sort_key) == ["-10", "-3", "2", "100", "100A", "100B", "101"]


def test_whole_chain_keys_distinct(renumbered):
    """A whole-chain pocket lists every residue once, each at its own seq_pos."""
    pocket = parse_pocket_from_struct(renumbered, "A", None)
    ids = pocket.res_auth_ids
    assert len(set(ids)) == len(ids) == len(pocket.residues) == len(pocket.ca_sequence)
    assert {"854", "854A", "854B"} <= set(ids)
    seq_positions = [residue.seq_pos for residue in pocket.residues.values()]
    assert len(set(seq_positions)) == len(seq_positions)


def test_passthrough_keys_seq_pos(renumbered):
    """Listed residues sharing a number get the seq_pos of their own residue in `ca_residues`."""
    expected = {
        residue_key(res.seqid.num, res.seqid.icode): seq_pos
        for seq_pos, res in enumerate(ca_residues(renumbered[0]["A"]))
    }
    pocket = parse_pocket_from_struct(renumbered, "A", sorted(["854B", "854", "854A"], key=residue_sort_key))
    assert pocket.res_auth_ids == ["854", "854A", "854B"]
    assert {key: pocket.residues[key].seq_pos for key in pocket.res_auth_ids} == {
        key: expected[key] for key in ("854", "854A", "854B")
    }


def test_vdw_contact_keys(renumbered):
    """Contact residues come back keyed with their insertion code."""
    model = renumbered[0]
    contacts = vdw_contact_residues(ca_residues(model["A"]), model["E"].get_polymer().first_conformer())
    assert "1104" in contacts
    assert "1104A" in contacts
    assert "1105" not in contacts


def interface_dict(seq_nums, ins_codes):
    """A one-molecule interface with one hydrogen bond per residue and no other bonds."""
    interface = {
        bond_type: {
            "atom_site_1_seq_nums": [],
            "atom_site_1_inscodes": [],
            "atom_site_1_residues": [],
            "atom_site_1_unp_nums": [],
        }
        for bond_type in BOND_TYPES
    }
    interface["hydrogen_bonds"] = {
        "atom_site_1_seq_nums": seq_nums,
        "atom_site_1_inscodes": ins_codes,
        "atom_site_1_residues": ["ALA"] * len(seq_nums),
        "atom_site_1_unp_nums": [None] * len(seq_nums),
    }
    return interface


def molecule_dict(seq_ids, ins_codes):
    """An interface molecule whose residues are all half buried."""
    return {
        "residue_seq_ids": seq_ids,
        "residue_ins_codes": ins_codes,
        "residue_label_comp_ids": ["ALA"] * len(seq_ids),
        "accessible_surface_areas": [10.0] * len(seq_ids),
        "buried_surface_areas": [5.0] * len(seq_ids),
    }


@pytest.mark.parametrize(
    "residues",
    [
        pytest.param(lambda nums, codes: bonded_residues(interface_dict(nums, codes), "1"), id="bonded"),
        pytest.param(lambda nums, codes: buried_residues(molecule_dict(nums, codes), 0), id="buried"),
    ],
)
def test_pisa_keys_keep_insertion_codes(residues):
    """A residue with an insertion code is kept beside the one without it, and sorts after it."""
    keys = residues([100, 101, 100], ["A", None, None])
    assert set(keys) == {"100", "100A", "101"}
    assert sorted(keys, key=residue_sort_key) == ["100", "100A", "101"]


# 4C59, a real entry with insertion codes


@pytest.fixture
def pdb_dir_4c59(local_root):
    """`local_root/pdb`, where the `processor` fixture's PDB records point, holding the 4C59 structure."""
    pdb_dir = os.path.join(local_root, "pdb")
    os.makedirs(pdb_dir)
    shutil.copy(STRUCTURE_4C59, os.path.join(pdb_dir, "4C59.cif.gz"))
    return pdb_dir


@pytest.fixture
def pisa_dir_4c59(tmp_path):
    """A PISA cache holding only 4C59's interface file, which is all a downloader with downloads off reads."""
    pisa_dir = os.path.join(str(tmp_path), "pisa")
    os.makedirs(os.path.join(pisa_dir, "interface_pairs"))
    shutil.copy(INTERFACES_4C59, os.path.join(pisa_dir, "interface_pairs", "4c59.json"))
    return pisa_dir


@pytest.fixture
def positions_4c59():
    """Residue key -> index in `ca_residues` for 4C59 chain B: the seq_pos every builder must give."""
    st = load_structure(STRUCTURE_4C59, {})
    return {residue_key(res.seqid.num, res.seqid.icode): i for i, res in enumerate(ca_residues(st[0]["B"]))}


def parsed_record(processor, entry):
    """The record the parser builds for `entry`, which must parse."""
    ((_, record, reason),) = processor.process_qt_cmdline_input(entry, "query")
    assert reason is None
    return record


def assert_keyed_by_position(pocket, positions):
    """Every residue of `pocket` is listed once, in display order, at the seq_pos of its own key."""
    ids = pocket.res_auth_ids
    assert len(set(ids)) == len(ids) == len(pocket.residues)
    assert {key: pocket.residues[key].seq_pos for key in ids} == {key: positions[key] for key in ids}


def test_4c59_whole_chain(processor, pdb_dir_4c59, positions_4c59):
    """The whole chain lists 111A-111E as five residues, in chain order."""
    record = parsed_record(processor, "4C59:B")
    pockets = whole_chain_pockets([record], StructureManager(pdb_dir_4c59, download=False))
    pocket = pockets["4C59:B"]
    assert pocket.res_auth_ids == list(positions_4c59)
    assert {"111A", "111B", "111C", "111D", "111E"} <= set(pocket.res_auth_ids)
    assert "111" not in pocket.residues
    assert_keyed_by_position(pocket, positions_4c59)


def test_4c59_passthrough(processor, pdb_dir_4c59, positions_4c59):
    """Listed ids with insertion codes are built in display order, each at its own seq_pos."""
    record = parsed_record(processor, "4C59:B:111C,110,111A,111B,112")
    pockets = passthrough_pockets([record], StructureManager(pdb_dir_4c59, download=False))
    pocket = pockets["4C59:B:111C,110,111A,111B,112"]
    assert pocket.res_auth_ids == ["110", "111A", "111B", "111C", "112"]
    assert_keyed_by_position(pocket, positions_4c59)


def test_4c59_passthrough_bare_number(processor, pdb_dir_4c59, warnings):
    """111 alone names no residue of 4C59 chain B, which has only 111A-111E, so the entry is left out."""
    record = parsed_record(processor, "4C59:B:110,111")
    pockets = passthrough_pockets([record], StructureManager(pdb_dir_4c59, download=False))
    assert pockets == {}
    assert any("Residue(s) 111 of 4C59:B:110,111" in message for message in warnings())


# (bsa_threshold, residue count, insertion-code residues), read from 4C59's interface file. At 0, 111C is the one
# coded residue its partner does not bury.
PISA_4C59_CASES = [
    pytest.param(None, 13, ["111A"], id="bonded"),
    pytest.param(0, 20, ["111A", "111B", "111D", "111E"], id="buried"),
]


@pytest.mark.parametrize("bsa_threshold, count, coded", PISA_4C59_CASES)
def test_4c59_pisa(no_downloads, processor, pdb_dir_4c59, pisa_dir_4c59, positions_4c59, bsa_threshold, count, coded):
    """A PISA pocket keeps each insertion-code residue apart, sorted after its number, at its own seq_pos."""
    record = parsed_record(processor, "4C59:B_A")
    pockets = pisa_pockets(
        [record],
        StructureManager(pdb_dir_4c59, download=False),
        PisaDownloader(pisa_dir_4c59, download=False),
        bsa_threshold=bsa_threshold,
    )
    pocket = pockets["4C59:B_A"]
    ids = pocket.res_auth_ids
    assert len(ids) == count
    assert [key for key in ids if not key.isdigit()] == coded
    assert ids == sorted(ids, key=residue_sort_key)
    assert_keyed_by_position(pocket, positions_4c59)
