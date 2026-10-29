"""
Shared fixtures for the unit tests: a scratch working directory and the objects built against it.
"""

import os

import pytest

from pocketmapper.parse.qt_processor import QTProcessor


@pytest.fixture
def local_root(tmp_path):
    """The test's scratch directory, as a string."""
    return str(tmp_path)


@pytest.fixture
def processor(local_root):
    """A QTProcessor resolving local files against `local_root`."""
    return QTProcessor(local_root=local_root)


@pytest.fixture
def write_file(local_root):
    """A factory that writes `content` to `name` under `local_root` and returns the file's absolute path."""

    def write(name, content):
        """Write `content` (str or bytes) to `name` under `local_root` and return the absolute path."""
        path = os.path.join(local_root, name)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        mode = "wb" if isinstance(content, bytes) else "w"
        with open(path, mode) as f:
            f.write(content)
        return path

    return write


@pytest.fixture
def cache_dirs(local_root):
    """The cache directory arguments of the step functions, each under `local_root`."""
    return {
        "pdb_dir": os.path.join(local_root, "pdb"),
        "alphafold_dir": os.path.join(local_root, "af"),
        "fsdb_dir": os.path.join(local_root, "fsdb"),
    }
