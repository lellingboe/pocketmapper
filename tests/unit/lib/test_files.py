"""
Unit tests for `lib.files.part_path`.
"""

import os

from pocketmapper.lib.files import part_path


def test_part_path_unique_beside_destination(local_root):
    """Each call gives a new `.part` name in the destination's directory, starting with its name, creating nothing."""
    dest = os.path.join(local_root, "sub", "4Q5J.cif.gz")
    names = [part_path(dest) for _ in range(3)]
    assert len(set(names)) == len(names)
    for name in names:
        assert os.path.dirname(name) == os.path.dirname(dest)
        assert name.startswith(f"{dest}.")
        assert name.endswith(".part")
        assert not os.path.exists(name)
