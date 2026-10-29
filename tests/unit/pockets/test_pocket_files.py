"""
Unit tests for `pocketmapper.pockets.pocket_files`: the pockets file round trip and its rejection of bad files.
"""

import os

import pytest

from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.pockets.pocket import Pocket, PocketResidue
from pocketmapper.pockets.pocket_files import POCKETS_FILE_VERSION, read_pockets_file, write_pockets_file

# (file content, id): pockets files `read_pockets_file` must reject
BAD_POCKETS_FILES = [
    pytest.param(f'{{"version": {POCKETS_FILE_VERSION - 1}, "pockets": {{}}, "chains": {{}}}}', id="wrong_version"),
    pytest.param("not json at all", id="not_json"),
    pytest.param(f'{{"version": {POCKETS_FILE_VERSION}, "chains": {{}}}}', id="missing_pockets"),
    pytest.param(f'{{"version": {POCKETS_FILE_VERSION}, "pockets": {{}}}}', id="missing_chains"),
]


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
