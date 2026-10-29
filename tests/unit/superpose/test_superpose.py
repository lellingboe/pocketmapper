"""
Unit tests for `superpose.superpose`: the Foldseek database path handed to the aligner.
"""

import os

import pytest

from pocketmapper.lib.foldseek import BUNDLED_HUMAN_DOMAINS_DB
from pocketmapper.lib.structure_manager import StructureManager
from pocketmapper.superpose import superpose
from pocketmapper.superpose.structure_aligner import StructureAligner


@pytest.fixture
def aligner_calls(monkeypatch):
    """Replace `StructureAligner.align_structs` with a recorder; return the calls, each its manager and arguments."""
    calls = []

    def align_structs(self, **kwargs):
        """Record the aligner's manager and the arguments, and write nothing."""
        calls.append(dict(kwargs, manager=self.manager))
        return {}

    monkeypatch.setattr(StructureAligner, "align_structs", align_structs)
    return calls


@pytest.mark.parametrize(
    "struct_info, fsdb_path",
    [
        ("pdb", os.path.join("{fsdb_dir}", "pdb")),
        ("human_domains", BUNDLED_HUMAN_DOMAINS_DB),
        (os.path.join("dbs", "mydb"), os.path.join("{local_root}", "dbs", "mydb")),
    ],
)
def test_fsdb_path_from_manager(write_file, local_root, cache_dirs, aligner_calls, struct_info, fsdb_path):
    """The aligner gets the database's path from the manager: a bundled name, or a user path under `local_root`."""
    write_file("query.cif.gz", "")
    alignment = write_file("alignment.tsv", "query\ttarget\nquery_A\tAF-P1-F1\n")
    query = {"pocket_id": "Q:A", "struct_type": "local_file", "struct_info": "query.cif.gz"}
    database = {"pocket_id": "db", "struct_type": "foldseek_db", "struct_info": struct_info}
    manager = StructureManager(**cache_dirs, local_root=local_root, download=False)

    superpose.superpose_top_targets(
        {"query": [query], "target": [database]},
        manager,
        None,
        "comparison.tsv",
        alignment,
        None,
        os.path.join(local_root, "aligned"),
        "foldseek",
        10,
        1,
    )

    expected = fsdb_path.format(fsdb_dir=cache_dirs["fsdb_dir"], local_root=local_root)
    [call] = aligner_calls
    assert os.path.normpath(call["fsdb_path"]) == os.path.normpath(expected)
    assert call["manager"] is manager
    assert call["query_records"] == [query]
    assert call["target_records"] == []
