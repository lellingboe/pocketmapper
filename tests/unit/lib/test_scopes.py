"""
Unit tests for `lib.scopes.delete_temp_dir`: what `delete_tmp` keeps and what it removes.
"""

import os

import pytest

from pocketmapper.lib.constants import TEMP_SUBDIRS
from pocketmapper.lib.scopes import delete_temp_dir

# What each TEMP_SUBDIRS subdirectory holds, so a deletion has to recurse
SUBDIR_FILE = os.path.join("nested", "structure.cif")


def make_temp_dir(temp_dir, extras):
    """Create `temp_dir` with every TEMP_SUBDIRS subdirectory, each holding a file, plus the `extras` paths."""
    for name in TEMP_SUBDIRS.values():
        write(os.path.join(temp_dir, name, SUBDIR_FILE))
    for extra in extras:
        write(os.path.join(temp_dir, extra))


def write(path):
    """Write a small file at `path`, making its directory."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("x\n")


def listing(temp_dir):
    """Every file under `temp_dir`, as sorted relative paths, or None if `temp_dir` does not exist."""
    if not os.path.isdir(temp_dir):
        return None
    paths = []
    for root, _, names in os.walk(temp_dir):
        paths.extend(os.path.relpath(os.path.join(root, name), temp_dir) for name in names)
    return sorted(paths)


@pytest.fixture
def temp_dir(local_root):
    """The path of a scratch directory under `local_root`, not yet created."""
    return os.path.join(local_root, "tmp")


@pytest.mark.parametrize("extras", [[], ["user_file.txt", os.path.join("user_dir", "notes.txt")]])
def test_delete_tmp_0_keeps_everything(temp_dir, extras):
    """With delete_tmp 0, every subdirectory and every other file is kept."""
    make_temp_dir(temp_dir, extras)
    before = listing(temp_dir)
    delete_temp_dir(temp_dir, 0)
    assert listing(temp_dir) == before


def test_delete_tmp_1_removes_dir_left_empty(temp_dir):
    """With delete_tmp 1, the TEMP_SUBDIRS go, and then `temp_dir` itself, as nothing else is in it."""
    make_temp_dir(temp_dir, [])
    delete_temp_dir(temp_dir, 1)
    assert not os.path.exists(temp_dir)


@pytest.mark.parametrize(
    "extras",
    [
        ["user_file.txt"],
        [os.path.join("user_dir", "notes.txt")],
        ["user_file.txt", os.path.join("user_dir", "notes.txt")],
    ],
)
def test_delete_tmp_1_keeps_other_files(temp_dir, extras):
    """With delete_tmp 1, only the TEMP_SUBDIRS go; other files keep `temp_dir` in place."""
    make_temp_dir(temp_dir, extras)
    delete_temp_dir(temp_dir, 1)
    assert listing(temp_dir) == sorted(extras)
    for name in TEMP_SUBDIRS.values():
        assert not os.path.exists(os.path.join(temp_dir, name))


def test_delete_tmp_1_partial_subdirs(temp_dir):
    """With delete_tmp 1, a `temp_dir` holding only some of the TEMP_SUBDIRS is still removed."""
    write(os.path.join(temp_dir, TEMP_SUBDIRS["foldseek"], SUBDIR_FILE))
    delete_temp_dir(temp_dir, 1)
    assert not os.path.exists(temp_dir)


def test_delete_tmp_1_empty_dir(temp_dir):
    """With delete_tmp 1, an empty `temp_dir` is removed."""
    os.makedirs(temp_dir)
    delete_temp_dir(temp_dir, 1)
    assert not os.path.exists(temp_dir)


@pytest.mark.parametrize("delete_tmp", [0, 1])
def test_delete_temp_dir_never_made(temp_dir, delete_tmp):
    """A `temp_dir` that was never made is left absent, without error."""
    delete_temp_dir(temp_dir, delete_tmp)
    assert not os.path.exists(temp_dir)
