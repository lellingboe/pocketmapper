"""
Unit tests for `superpose.structure_aligner`: locating structures through a `StructureManager`, and the records of
Foldseek-database hits.
"""

import os
import shutil

import numpy as np
import pandas as pd
import pytest

from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.files import safe_filename
from pocketmapper.lib.structure_manager import StructureManager
from pocketmapper.superpose import structure_aligner
from pocketmapper.superpose.structure_aligner import StructureAligner

TESTS_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIXTURE_STRUCTURE = os.path.join(TESTS_DIR, "integration", "fixtures", "4Q5J.cif.gz")
IDENTITY_U = ",".join(str(x) for x in np.eye(3).flatten())
ZERO_T = "0,0,0"


@pytest.fixture
def manager(local_root, cache_dirs):
    """A StructureManager with every cache directory under `local_root`, which it resolves local files against."""
    return StructureManager(**cache_dirs, local_root=local_root, download=False)


@pytest.fixture
def fake_extract(monkeypatch):
    """Replace `extract_fsdb_structures` with a copy of the fixture per entry; return the calls it got."""
    calls = []

    def extract(db_path, entry_names, out_dir, threads, log_extra=None):
        """Copy the fixture structure in place of each entry and return entry name -> its path."""
        calls.append({"db_path": db_path, "entry_names": entry_names, "out_dir": out_dir})
        struct_dir = os.path.join(out_dir, "fsdb_structures")
        os.makedirs(struct_dir, exist_ok=True)
        paths = {name: os.path.join(struct_dir, f"{name}.cif.gz") for name in entry_names}
        for path in paths.values():
            shutil.copy2(FIXTURE_STRUCTURE, path)
        return paths

    monkeypatch.setattr(structure_aligner, "extract_fsdb_structures", extract)
    return calls


@pytest.mark.parametrize(
    "struct_type, struct_info, stored_at",
    [
        ("pdb", "4q5j", os.path.join("pdb", "4Q5J.cif.gz")),
        ("alphafold", "P12345", os.path.join("af", "P12345.cif.gz")),
        ("local_file", os.path.join("structs", "mine.cif.gz"), os.path.join("structs", "mine.cif.gz")),
        ("local_file", os.path.join("{root}", "elsewhere", "mine.cif.gz"), os.path.join("elsewhere", "mine.cif.gz")),
    ],
)
def test_transform_reads_record_path(local_root, manager, struct_type, struct_info, stored_at):
    """`transform` reads each record's structure where the manager locates it, an absolute local file included."""
    struct_path = os.path.join(local_root, stored_at)
    os.makedirs(os.path.dirname(struct_path), exist_ok=True)
    shutil.copy2(FIXTURE_STRUCTURE, struct_path)
    record = {"struct_type": struct_type, "struct_info": struct_info.format(root=local_root), "chain_info": "A_F"}
    records = [dict(record, pocket_id="query"), dict(record, pocket_id="target")]
    out_path = os.path.join(local_root, "out", "aligned.pdb")

    written = StructureAligner(manager).transform(records, [None, (np.eye(3), np.zeros(3))], out_path)

    assert written == out_path
    with open(out_path) as f:
        text = f.read()
    assert "MOLECULE: query;" in text
    assert "MOLECULE: target;" in text


def test_transform_missing_dir_raises(local_root):
    """A record whose type needs a directory the manager lacks raises, rather than being dropped."""
    record = {"pocket_id": "4Q5J:A", "struct_type": "pdb", "struct_info": "4Q5J", "chain_info": "A"}
    aligner = StructureAligner(StructureManager(local_root=local_root))

    with pytest.raises(PocketMapperError, match="needs pdb_dir"):
        aligner.transform([record], [None], os.path.join(local_root, "aligned.pdb"))


@pytest.mark.parametrize(
    "target_records, target_ids, entry_names",
    [
        # Records given: ids are pocket_ids, read through preprocess_name; an id no record covers is dropped
        (
            [{"pocket_id": "4Q5J:B_F", "preprocess_name": "4q5j-assembly1_B"}],
            ["4Q5J:B_F", "1ABC:A"],
            {"4Q5J:B_F": "4q5j-assembly1_B"},
        ),
        # No records: ids are the database's entry names
        ([], ["AF-P1-F1", "AF-P2-F1"], {"AF-P1-F1": "AF-P1-F1", "AF-P2-F1": "AF-P2-F1"}),
    ],
)
def test_fsdb_target_records(local_root, monkeypatch, fake_extract, target_records, target_ids, entry_names):
    """Each hit record names its rebuilt structure as an absolute `local_file`, whatever the out_dir or local_root."""
    monkeypatch.chdir(local_root)
    out_dir = os.path.abspath("aligned")
    elsewhere = StructureManager(local_root=os.path.join(local_root, "elsewhere"))

    records = StructureAligner(elsewhere).fsdb_target_records("db", target_records, target_ids, "aligned", 1)

    assert fake_extract == [{"db_path": "db", "entry_names": list(entry_names.values()), "out_dir": out_dir}]
    assert records == {
        target_id: {
            "pocket_id": target_id,
            "preprocess_name": entry_name,
            "chain_info": None,
            "struct_type": "local_file",
            "struct_info": os.path.join(out_dir, "fsdb_structures", f"{entry_name}.cif.gz"),
        }
        for target_id, entry_name in entry_names.items()
    }
    for record in records.values():
        assert elsewhere.record_path(record) == record["struct_info"]


def test_align_structs_fsdb(local_root, manager, fake_extract):
    """Against a Foldseek database, the rebuilt hits are superposed onto a query read through the manager."""
    shutil.copy2(FIXTURE_STRUCTURE, os.path.join(local_root, "query.cif.gz"))
    query = {
        "pocket_id": "Q:A",
        "preprocess_name": "query_A",
        "struct_type": "local_file",
        "struct_info": "query.cif.gz",
        "chain_info": "A",
    }
    comparison = pd.DataFrame(
        {
            "query": ["Q:A", "Q:A"],
            "target": ["AF-P1-F1", "AF-P2-F1"],
            "overlap_count": [5, 0],
            "jaccard_index": [0.5, 0.0],
            "min_overlap_similarity": [0.5, 0.0],
        }
    )
    alignment = pd.DataFrame({"query": ["query_A"], "target": ["AF-P1-F1"], "u": [IDENTITY_U], "t": [ZERO_T]})
    out_dir = os.path.join(local_root, "aligned")

    out_paths = StructureAligner(manager).align_structs(
        [query], [], comparison, out_dir, "foldseek", alignment=alignment, fsdb_path="db"
    )

    out_path = os.path.join(out_dir, f"{safe_filename('Q:A')}.pdb")
    assert out_paths == {"Q:A": out_path}
    assert fake_extract[0]["entry_names"] == ["AF-P1-F1"]
    with open(out_path) as f:
        text = f.read()
    assert "MOLECULE: Q:A;" in text
    assert "MOLECULE: AF-P1-F1;" in text
