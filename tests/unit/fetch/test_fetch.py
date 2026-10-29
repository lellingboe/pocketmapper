"""
Unit tests for `fetch.fetch`: a missing Foldseek database with downloads turned off.
"""

import os
import re

import pytest

from pocketmapper.fetch import fetch
from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.structure_manager import StructureManager


@pytest.fixture
def no_foldseek(monkeypatch):
    """Make any foldseek probe or run in `fetch` fail the test."""

    def forbidden(*args, **kwargs):
        """Fail the test: nothing may be downloaded."""
        raise AssertionError("foldseek was called")

    monkeypatch.setattr(fetch, "require_foldseek", forbidden)
    monkeypatch.setattr(fetch, "run_foldseek", forbidden)


def test_fetch_missing_fsdb_no_download(local_root, no_foldseek):
    """A missing database with `download=False` raises, without probing or running foldseek."""
    fsdb_path = os.path.join(local_root, "fsdb", "pdb")
    record = {"struct_type": "foldseek_db", "struct_info": "pdb"}
    manager = StructureManager(fsdb_dir=os.path.join(local_root, "fsdb"), download=False)

    message = f"The Foldseek database 'PDB' is not at {fsdb_path}, and downloads are turned off"
    with pytest.raises(PocketMapperError, match=f"^{re.escape(message)}"):
        fetch.fetch_missing_fsdb(record, manager, 1, os.path.join(local_root, "tmp"))
    assert not os.path.exists(os.path.dirname(fsdb_path))


def test_fetch_missing_fsdb_cached_no_download(write_file, local_root, no_foldseek):
    """A database already on disk passes with `download=False`."""
    write_file(os.path.join("fsdb", "pdb"), "")
    record = {"struct_type": "foldseek_db", "struct_info": "pdb"}
    manager = StructureManager(fsdb_dir=os.path.join(local_root, "fsdb"), download=False)

    fetch.fetch_missing_fsdb(record, manager, 1, os.path.join(local_root, "tmp"))
