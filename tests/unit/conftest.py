"""
Shared fixtures for the unit tests: a scratch working directory and the objects built against it.
"""

import os

import pytest

from pocketmapper.parse.qt_processor import QTProcessor


@pytest.fixture
def work_dir(tmp_path):
    """The test's scratch directory, as a string."""
    return str(tmp_path)


@pytest.fixture
def processor(work_dir):
    """A QTProcessor with relative cache directories, resolved against `work_dir`."""
    return QTProcessor(pdb_dir="pdb", alphafold_dir="af", fsdb_dir="fsdb", work_dir=work_dir)


@pytest.fixture
def write_file(work_dir):
    """A factory that writes `content` to `name` under `work_dir` and returns the file's absolute path."""

    def write(name, content):
        """Write `content` (str or bytes) to `name` under `work_dir` and return the absolute path."""
        path = os.path.join(work_dir, name)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        mode = "wb" if isinstance(content, bytes) else "w"
        with open(path, mode) as f:
            f.write(content)
        return path

    return write


@pytest.fixture
def cache_dirs(work_dir):
    """The cache directory arguments of the step functions, each under `work_dir`."""
    return {
        "pdb_dir": os.path.join(work_dir, "pdb"),
        "alphafold_dir": os.path.join(work_dir, "af"),
        "fsdb_dir": os.path.join(work_dir, "fsdb"),
    }
