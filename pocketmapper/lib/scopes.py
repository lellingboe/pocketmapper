"""
Run-state scopes: the scratch directory a block owns, and the command a run names.

The only module state in the package's helpers lives here: `HELD_TEMP_DIRS`, the scratch directories an
open `temp_dir_scope` holds, and `OUTER_COMMAND`, the command an open `run_scope` names.
"""

import logging
import os
import shutil
from contextlib import contextmanager

from pocketmapper.lib.constants import TEMP_SUBDIRS

logger = logging.getLogger(__name__)

# Absolute paths of the scratch directories an open temp_dir_scope holds
HELD_TEMP_DIRS = set()

# The command of the outermost open run_scope; empty when none is open
OUTER_COMMAND = []


def empty_temp_dir(temp_dir):
    """
    Delete the scratch subdirectories so that whatever uses them next starts with them empty.

    Only the subdirectories in TEMP_SUBDIRS are deleted, never `temp_dir` itself or anything else in
    it, so a `temp_dir` pointed at a directory in use loses nothing.

    Args:
        temp_dir (str): The scratch directory. Not created here.

    Returns:
        None
    """
    for name in TEMP_SUBDIRS.values():
        shutil.rmtree(os.path.join(temp_dir, name), ignore_errors=True)


def delete_temp_dir(temp_dir, delete_tmp):
    """
    Delete a scratch directory once it is no longer needed.

    Nothing is deleted if `delete_tmp` is 0 or the directory was never created. Only the subdirectories
    in TEMP_SUBDIRS are deleted, then `temp_dir` if that leaves it empty.

    Args:
        temp_dir (str): The scratch directory.
        delete_tmp (int): 1 to delete it, 0 to keep it.

    Returns:
        None
    """
    log_extra = {"stage": "Cleaning Up"}

    # Only a step that used scratch space created it
    if not os.path.isdir(temp_dir):
        logger.debug(f"No temp_dir was created at {temp_dir}; nothing to delete", extra=log_extra)
        return

    if delete_tmp == 0:
        logger.info(f"delete_tmp is 0; keeping {temp_dir}", extra=log_extra)
        return

    empty_temp_dir(temp_dir)
    if not os.listdir(temp_dir):
        os.rmdir(temp_dir)
    else:
        logger.debug(f"Leaving {temp_dir} in place: it holds files PocketMapper did not make", extra=log_extra)


@contextmanager
def temp_dir_scope(temp_dir, delete_tmp):
    """
    Own a scratch directory for the length of a `with` block: empty it on entry, delete it on exit.

    The emptying and deletion are `empty_temp_dir`'s and `delete_temp_dir`'s. A block that raises
    deletes nothing, leaving the scratch for inspection. A scope for a directory an enclosing scope
    already holds does nothing, on entry or exit: the outer scope owns it.

    Args:
        temp_dir (str): The scratch directory.
        delete_tmp (int): 1 to delete it on exit, 0 to keep it.

    Yields:
        None
    """
    key = os.path.abspath(temp_dir)
    if key in HELD_TEMP_DIRS:
        yield
        return
    HELD_TEMP_DIRS.add(key)
    # Released however the block ends, or a failed run would leave the next one in this process
    # neither emptying nor deleting it
    try:
        empty_temp_dir(temp_dir)
        yield
        delete_temp_dir(temp_dir, delete_tmp)
    finally:
        HELD_TEMP_DIRS.discard(key)


@contextmanager
def run_scope(command):
    """
    Mark a command as running for the length of a `with` block, unless an enclosing command is.

    A scope opened inside another does nothing: the outer command owns the run.

    Args:
        command (str): The command, e.g. "search" or "align".

    Yields:
        bool: True for the outermost scope, False for one nested in another.
    """
    if OUTER_COMMAND:
        yield False
        return
    OUTER_COMMAND.append(command)
    # Released however the block ends, or a failed run would leave the next one in this process nested
    try:
        yield True
    finally:
        OUTER_COMMAND.clear()


def outer_command():
    """
    The command the outermost open `run_scope` names.

    Returns:
        str: The command, or None when no scope is open.
    """
    return OUTER_COMMAND[0] if OUTER_COMMAND else None
