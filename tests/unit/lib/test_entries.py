"""
Unit tests for `lib.entries`: setting aside the records whose structure is not on disk.
"""

import os

import pytest

from pocketmapper.lib.entries import split_missing_structures
from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.structure_manager import StructureManager


@pytest.mark.parametrize(
    "struct_type, struct_info, written, kept",
    [
        pytest.param("pdb", "4q5j", os.path.join("pdb", "4Q5J.cif.gz"), True, id="pdb_present"),
        pytest.param("pdb", "4Q5J", None, False, id="pdb_missing"),
        pytest.param("alphafold", "p12345", os.path.join("af", "P12345.cif.gz"), True, id="alphafold_present"),
        pytest.param("local_file", "mine.cif", "mine.cif", True, id="local_present"),
        pytest.param("local_file", "mine.cif", None, False, id="local_missing"),
        pytest.param(
            "foldseek_db",
            "human_domains",
            os.path.join("fsdb", "pdb"),
            True,
            id="fsdb_present",
        ),
        pytest.param("foldseek_db", "pdb", None, False, id="fsdb_missing"),
    ],
)
def test_split_missing_structures(write_file, local_root, struct_type, struct_info, written, kept):
    """A record is kept only if the file `record_path` gives exists."""
    if written is not None:
        write_file(written, "")
    manager = StructureManager(
        os.path.join(local_root, "pdb"),
        os.path.join(local_root, "af"),
        os.path.join(local_root, "fsdb"),
        local_root,
    )
    record = {"pocket_id": "q1", "struct_type": struct_type, "struct_info": struct_info}

    result = split_missing_structures([record], manager, "query", "align", None, {})
    assert result == ([record] if kept else [])


def test_split_missing_structures_needs_dir():
    """A record whose cache directory the manager lacks raises."""
    record = {"struct_type": "pdb", "struct_info": "4Q5J"}
    with pytest.raises(PocketMapperError, match="needs pdb_dir"):
        split_missing_structures([record], StructureManager(), "query", "align", None, {})
