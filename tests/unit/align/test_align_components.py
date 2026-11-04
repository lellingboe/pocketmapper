"""
Unit tests for the align components: structure location through the manager and `SequenceAligner.align_seqs`.
"""

import os
import shutil

import pytest

from pocketmapper.align.sequence_aligner import SequenceAligner, build_aligner
from pocketmapper.align.structure_preprocessor import StructurePreprocessor
from pocketmapper.lib.structure_manager import StructureManager

FIXTURE = os.path.join(os.path.dirname(__file__), "..", "..", "integration", "fixtures", "4Q5J.cif.gz")


@pytest.fixture
def manager(local_root):
    """A StructureManager whose local files resolve against `local_root`."""
    return StructureManager(local_root=local_root)


@pytest.fixture
def record(local_root):
    """A local-file record naming a copy of 4Q5J under `local_root`."""
    shutil.copyfile(FIXTURE, os.path.join(local_root, "copy.cif.gz"))
    return {
        "struct_type": "local_file",
        "struct_info": "copy.cif.gz",
        "chain_info": "A_E",
        "pocket_id": "copy:A_E",
        "preprocess_name": "copy_A_x",
    }


@pytest.mark.parametrize("chain, expected", [("A_E", True), ("Z_E", False)])
def test_preprocess_records_locates_with_manager(manager, record, local_root, chain, expected):
    """The reference structure comes from `record_path`; a chain the structure lacks fails; no scratch is left."""
    record["chain_info"] = chain
    cache, search = os.path.join(local_root, "cache"), os.path.join(local_root, "search")
    os.makedirs(search)
    status = StructurePreprocessor(manager).preprocess_records([record], cache, search)
    assert status == {record["pocket_id"]: expected}
    assert os.path.isfile(os.path.join(search, "copy_A_x.cif.gz")) is expected
    if expected:
        assert os.listdir(cache) == ["copy_A_x.cif.gz"]


def test_align_records_locates_with_manager(manager, record):
    """The sequence is read from `record_path`, so a record aligns perfectly with itself."""
    table = SequenceAligner(manager).align_records([record], [record])
    assert len(table) == 1
    assert table.iloc[0]["fident"] == 1
    assert table.iloc[0]["mismatch"] == 0


# (query, target, qaln, taln, qstart, qend, tstart, tend, mismatch, gapopen)
ALIGN_SEQS_CASES = [
    pytest.param(
        "ACDEFRHIKLM", "ACDEFWHIKLM", "ACDEFRHIKLM", "ACDEFWHIKLM", 1, 11, 1, 11, 1, 0, id="r-w-mismatch-not-gapped"
    ),
    pytest.param("PPPPACDEFHIK", "ACDEF", "ACDEF", "ACDEF", 5, 9, 1, 5, 0, 0, id="flanks-on-query-only"),
    pytest.param(
        "GGGGACDEFRHIKLMGGG",
        "PPACDEFWHIKLMPPPP",
        "ACDEFRHIKLM",
        "ACDEFWHIKLM",
        5,
        15,
        3,
        13,
        1,
        0,
        id="flanks-on-both-sequences",
    ),
    pytest.param(
        "ACDEFGHIKLMNPQRSTVWY",
        "ACDEFGHIKLMSTVWY",
        "ACDEFGHIKLMNPQRSTVWY",
        "ACDEFGHIKLM----STVWY",
        1,
        20,
        1,
        16,
        0,
        1,
        id="internal-gap",
    ),
]


@pytest.mark.parametrize("query, target, qaln, taln, qstart, qend, tstart, tend, mismatch, gapopen", ALIGN_SEQS_CASES)
def test_align_seqs_local_columns(query, target, qaln, taln, qstart, qend, tstart, tend, mismatch, gapopen):
    """Only the aligned region is kept, start..end locate it 1-based inclusive, and mismatch is a count."""
    columns = SequenceAligner(None).align_seqs(query, target, build_aligner())
    assert (columns["qaln"], columns["taln"]) == (qaln, taln)
    assert (columns["qstart"], columns["qend"], columns["tstart"], columns["tend"]) == (qstart, qend, tstart, tend)
    # The gapped strings hold exactly the residues start..end, no more and no fewer
    assert columns["qaln"].replace("-", "") == query[qstart - 1 : qend]
    assert columns["taln"].replace("-", "") == target[tstart - 1 : tend]
    assert columns["alnlen"] == len(qaln)
    assert columns["fident"] == sum(a == b for a, b in zip(qaln, taln)) / len(qaln)
    assert type(columns["mismatch"]) is int
    assert (columns["mismatch"], columns["gapopen"]) == (mismatch, gapopen)


def test_align_seqs_no_local_alignment():
    """Sequences with no positive-scoring residue pair have no local alignment, rather than an error."""
    assert SequenceAligner(None).align_seqs("W", "G", build_aligner()) is None
