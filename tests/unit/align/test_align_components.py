"""
Unit tests for the align components: they locate structures through the manager, not the record's `struct_path`.
"""

import os
import shutil

import pytest

from pocketmapper.align.sequence_aligner import SequenceAligner
from pocketmapper.align.structure_preprocessor import StructurePreprocessor
from pocketmapper.lib.structure_manager import StructureManager

FIXTURE = os.path.join(os.path.dirname(__file__), "..", "..", "integration", "fixtures", "4Q5J.cif.gz")


@pytest.fixture
def manager(local_root):
    """A StructureManager whose local files resolve against `local_root`."""
    return StructureManager(local_root=local_root)


@pytest.fixture
def record(local_root):
    """A local-file record naming a copy of 4Q5J under `local_root`, with a `struct_path` that does not exist."""
    shutil.copyfile(FIXTURE, os.path.join(local_root, "copy.cif.gz"))
    return {
        "struct_type": "local_file",
        "struct_info": "copy.cif.gz",
        "struct_path": os.path.join(local_root, "missing.cif.gz"),
        "chain_info": "A_E",
        "pocket_id": "copy:A_E",
        "preprocess_name": "copy_A_x",
    }


@pytest.mark.parametrize("chain, expected", [("A_E", True), ("Z_E", False)])
def test_preprocess_records_locates_with_manager(manager, record, local_root, chain, expected):
    """The reference structure comes from `record_path`; a chain the structure lacks fails."""
    record["chain_info"] = chain
    cache, search = os.path.join(local_root, "cache"), os.path.join(local_root, "search")
    os.makedirs(search)
    status = StructurePreprocessor(manager).preprocess_records([record], cache, search)
    assert status == {record["pocket_id"]: expected}
    assert os.path.isfile(os.path.join(search, "copy_A_x.cif.gz")) is expected


def test_align_records_locates_with_manager(manager, record):
    """The sequence is read from `record_path`, so a record aligns perfectly with itself."""
    table = SequenceAligner(manager).align_records([record], [record])
    assert len(table) == 1
    assert table.iloc[0]["fident"] == 1
