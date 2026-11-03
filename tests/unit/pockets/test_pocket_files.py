"""
Unit tests for `pocketmapper.pockets.pocket_files`: the pockets file round trip and its rejection of bad files.
"""

import json
import os
from dataclasses import replace

import pytest

from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.pockets.pocket import Pocket, PocketResidue
from pocketmapper.pockets.pocket_files import (
    POCKETS_FILE_VERSION,
    check_pocket,
    read_pockets_file,
    serialise_pockets,
    write_pockets_file,
)

# (file content, id): pockets files `read_pockets_file` must reject
BAD_POCKETS_FILES = [
    pytest.param(f'{{"version": {POCKETS_FILE_VERSION - 1}, "pockets": {{}}, "chains": {{}}}}', id="wrong_version"),
    pytest.param("not json at all", id="not_json"),
    pytest.param(f'{{"version": {POCKETS_FILE_VERSION}, "chains": {{}}}}', id="missing_pockets"),
    pytest.param(f'{{"version": {POCKETS_FILE_VERSION}, "pockets": {{}}}}', id="missing_chains"),
]


def make_pocket():
    """A sound pocket on the chain "AKG": residues "3" (A), "10" (lowercase k) and "25" (G), in file order."""
    residues = {
        "25": PocketResidue(res_code="GLY", res_code_single="G", seq_pos=2, ca_coords=[1.0, 2.0, 3.0]),
        "3": PocketResidue(res_code="ALA", res_code_single="A", seq_pos=0, uniprot_pos=7),
        "10": PocketResidue(res_code_single="k", seq_pos=1),
    }
    return Pocket(
        res_auth_ids=["3", "10", "25"], residues=residues, ca_sequence="AKG", pocket_exists=True, has_coords=True
    )


def break_residue(pocket, res_id, **fields):
    """Replace fields of one residue of `pocket` in place."""
    pocket.residues[res_id] = replace(pocket.residues[res_id], **fields)


def duplicate_id(pocket):
    """Repeat an id in `res_auth_ids`."""
    pocket.res_auth_ids.append("3")


def drop_residue(pocket):
    """Remove residue "10" from `residues`, leaving its id listed."""
    del pocket.residues["10"]


# (breaking function, id): a sound pocket made invalid, one per check
BREAKERS = [
    pytest.param(duplicate_id, id="duplicate_id"),
    pytest.param(drop_residue, id="id_missing_from_residues"),
    pytest.param(lambda p: break_residue(p, "3", seq_pos=None), id="seq_pos_none"),
    pytest.param(lambda p: break_residue(p, "3", seq_pos=-1), id="seq_pos_negative"),
    pytest.param(lambda p: break_residue(p, "3", seq_pos=3), id="seq_pos_too_large"),
    pytest.param(lambda p: break_residue(p, "3", seq_pos="0"), id="seq_pos_not_int"),
    pytest.param(lambda p: break_residue(p, "3", seq_pos=1), id="code_mismatch"),
    pytest.param(lambda p: break_residue(p, "3", res_code_single=None), id="code_none"),
]


def test_pockets_file_round_trip(local_root):
    """A written pockets file reads back equal, keeping res_auth_ids order, seq_pos, None fields and None pockets."""
    pocket = make_pocket()
    pockets = {
        "4Q5J:A_E": pocket,
        "1ABC:A_B": None,
    }
    chains = {"4Q5J_A_hash": ["4Q5J:A_E"], "1ABC_A_hash": ["1ABC:A_B"], "empty_chain": []}
    path = os.path.join(local_root, "out", "pockets.json")

    write_pockets_file(pockets, chains, path)
    read_pockets, read_chains = read_pockets_file(path)

    assert read_pockets == pockets
    assert read_pockets["4Q5J:A_E"].res_auth_ids == ["3", "10", "25"]
    assert list(read_pockets["4Q5J:A_E"].residues) == ["25", "3", "10"]
    assert read_chains == chains


@pytest.mark.parametrize("content", BAD_POCKETS_FILES)
def test_read_pockets_file_rejects(write_file, content):
    """A pockets file of the wrong version, not JSON or missing a key raises PocketMapperError."""
    path = write_file("pockets.json", content)

    with pytest.raises(PocketMapperError):
        read_pockets_file(path)


@pytest.mark.parametrize("breaker", BREAKERS)
def test_write_pockets_file_rejects_invalid_pocket(local_root, breaker):
    """A pocket that fails a `check_pocket` check raises on write, names the pocket, and leaves no file."""
    pocket = make_pocket()
    breaker(pocket)
    path = os.path.join(local_root, "pockets.json")

    with pytest.raises(PocketMapperError, match="4Q5J:A_E"):
        write_pockets_file({"1ABC:A_B": None, "4Q5J:A_E": pocket}, {}, path)

    assert not os.path.exists(path)


@pytest.mark.parametrize("breaker", BREAKERS)
def test_read_pockets_file_rejects_invalid_pocket(write_file, breaker):
    """A hand-written file holding a pocket that fails a check raises, naming the pocket, the file and a rerun."""
    pocket = make_pocket()
    breaker(pocket)
    content = {"version": POCKETS_FILE_VERSION, "pockets": serialise_pockets({"4Q5J:A_E": pocket}), "chains": {}}
    path = write_file("pockets.json", json.dumps(content))

    with pytest.raises(PocketMapperError, match=r"4Q5J:A_E.*pockets\.json.*rerun pockets"):
        read_pockets_file(path)


@pytest.mark.parametrize(
    "pocket",
    [
        pytest.param(make_pocket(), id="sound"),
        pytest.param(Pocket(res_auth_ids=[], residues={}, ca_sequence="AKG", pocket_exists=False), id="empty"),
        pytest.param(Pocket(), id="default"),
    ],
)
def test_check_pocket_accepts(pocket):
    """A sound pocket, and one built but empty under its rule, pass."""
    assert check_pocket("4Q5J:A_E", pocket) is None


@pytest.mark.parametrize("code", ["k", "K"])
def test_check_pocket_ignores_case(code):
    """A residue code matches `ca_sequence` whichever case either is in."""
    pocket = make_pocket()
    break_residue(pocket, "10", res_code_single=code)

    assert check_pocket("4Q5J:A_E", pocket) is None
    pocket.ca_sequence = "AkG"
    assert check_pocket("4Q5J:A_E", pocket) is None


def test_check_pocket_ignores_extra_residue_keys():
    """A key of `residues` outside `res_auth_ids` is not an error: nothing reads it."""
    pocket = make_pocket()
    pocket.res_auth_ids.remove("25")

    assert check_pocket("4Q5J:A_E", pocket) is None
