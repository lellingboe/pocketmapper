"""
Run-state scopes: the scratch directory a block owns, and the command a run names.

The only module state in the package's helpers lives here: `HELD_TEMP_DIRS`, the scratch directories an
open `temp_dir_scope` holds, and `OUTER_COMMAND`, the command an open `run_scope` names.
"""

import logging
import os
import shutil
from contextlib import contextmanager

logger = logging.getLogger(__name__)

# Absolute paths of the scratch directories an open temp_dir_scope holds
HELD_TEMP_DIRS = set()

# The command of the outermost open run_scope; empty when none is open
OUTER_COMMAND = []


def is_within(path, roots):
    """
    Report whether a path resolves to somewhere inside one of `roots`.

    Args:
        path (str): The path to test.
        roots (list): Candidate containing directories; only one has to match.

    Returns:
        bool: True if `path` is inside (or equal to) any of `roots`, both sides resolved with
            `os.path.realpath` so `..` segments and symlinks cannot walk out of a root.
    """
    real_path = os.path.realpath(path)
    for root in roots:
        real_root = os.path.realpath(root)
        # commonpath raises on a mix of absolute and relative paths; realpath has made both
        # absolute, so the only remaining raiser is a different drive, which cannot be a match.
        try:
            if os.path.commonpath([real_path, real_root]) == real_root:
                return True
        except ValueError:
            continue
    return False


def empty_temp_dir(temp_dir, roots):
    """
    Delete a scratch directory so that whatever uses it next starts with it empty.

    Only a `temp_dir` inside one of `roots` is deleted, so a mistyped path costs a stray directory,
    not its contents. One outside them is reused as it is, with a warning if it exists.

    Args:
        temp_dir (str): The scratch directory. Not created here.
        roots (list): The directories it may be deleted under.

    Returns:
        None
    """
    log_extra = {"stage": "Configuring Workflow"}

    if is_within(temp_dir, roots):
        shutil.rmtree(temp_dir, ignore_errors=True)
    elif os.path.exists(temp_dir):
        logger.warning(
            f"Reusing temp_dir {temp_dir} without emptying it: it is outside {' and '.join(roots)}. "
            "Empty it yourself if a previous run left anything there.",
            extra=log_extra,
        )


def delete_temp_dir(temp_dir, delete_tmp, roots):
    """
    Delete a scratch directory once it is no longer needed.

    Nothing is deleted if `delete_tmp` is 0, the directory was never created, or it resolves
    outside every one of `roots` (warned about).

    Args:
        temp_dir (str): The scratch directory.
        delete_tmp (int): 1 to delete it, 0 to keep it.
        roots (list): The directories it may be deleted under.

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

    if not is_within(temp_dir, roots):
        logger.warning(
            f"Not deleting temp_dir {temp_dir}: it is outside {' and '.join(roots)}. "
            "Remove it yourself if that was intended.",
            extra=log_extra,
        )
        return
    shutil.rmtree(temp_dir)


@contextmanager
def temp_dir_scope(temp_dir, delete_tmp, roots):
    """
    Own a scratch directory for the length of a `with` block: empty it on entry, delete it on exit.

    The emptying and deletion are `empty_temp_dir`'s and `delete_temp_dir`'s, with their guards. A
    block that raises deletes nothing, leaving the scratch for inspection. A scope for a directory an
    enclosing scope already holds does nothing, on entry or exit: the outer scope owns it.

    Args:
        temp_dir (str): The scratch directory.
        delete_tmp (int): 1 to delete it on exit, 0 to keep it.
        roots (list): The directories it may be emptied or deleted under.

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
        empty_temp_dir(temp_dir, roots)
        yield
        delete_temp_dir(temp_dir, delete_tmp, roots)
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
