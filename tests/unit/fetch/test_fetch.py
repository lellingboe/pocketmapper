"""
Unit tests for `fetch.fetch`: a missing Foldseek database with downloads turned off.
"""

import os
import re

import pytest

from pocketmapper.fetch import fetch
from pocketmapper.lib.exceptions import PocketMapperError


@pytest.fixture
def no_foldseek(monkeypatch):
    """Make any foldseek probe or run in `fetch` fail the test."""

    def forbidden(*args, **kwargs):
        """Fail the test: nothing may be downloaded."""
        raise AssertionError("foldseek was called")

    monkeypatch.setattr(fetch, "require_foldseek", forbidden)
    monkeypatch.setattr(fetch, "run_foldseek", forbidden)


def test_fetch_missing_fsdb_no_download(work_dir, no_foldseek):
    """A missing database with `download=False` raises, without probing or running foldseek."""
    fsdb_path = os.path.join(work_dir, "fsdb", "human_domains", "human_domains")
    record = {"struct_type": "foldseek_db", "struct_info": "human_domains", "struct_path": fsdb_path}

    message = f"The Foldseek database 'HUMAN_DOMAINS' is not at {fsdb_path}, and downloads are turned off"
    with pytest.raises(PocketMapperError, match=f"^{re.escape(message)}"):
        fetch.fetch_missing_fsdb(record, 1, os.path.join(work_dir, "tmp"), download=False)
    assert not os.path.exists(os.path.dirname(fsdb_path))


def test_fetch_missing_fsdb_cached_no_download(write_file, work_dir, no_foldseek):
    """A database already on disk passes with `download=False`."""
    fsdb_path = write_file(os.path.join("fsdb", "human_domains", "human_domains"), "")
    record = {"struct_type": "foldseek_db", "struct_info": "human_domains", "struct_path": fsdb_path}

    fetch.fetch_missing_fsdb(record, 1, os.path.join(work_dir, "tmp"), download=False)
