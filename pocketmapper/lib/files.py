"""
Stateless file and filename helpers.
"""

import gzip
import hashlib
import logging
import os
import re
import shutil

from pocketmapper.exceptions import PocketMapperError

logger = logging.getLogger(__name__)

UNSAFE_FILENAME_CHARS = re.compile(r"[^A-Za-z0-9._-]")


def make_dir(path, log_extra):
    """
    Create a directory and any missing parents, if it does not already exist.

    Args:
        path (str): The directory. An empty string, the directory of a bare filename, is the
            working directory and is left alone.
        log_extra (dict): Logging `extra` for the failure message.

    Returns:
        None

    Raises:
        PocketMapperError: If the directory cannot be created.
    """
    if not path:
        return
    try:
        os.makedirs(path, exist_ok=True)
    except OSError as e:
        logger.critical(f"Error creating directory {path}", extra=log_extra)
        raise PocketMapperError(f"Error creating directory {path}") from e


def gzip_file(src_fpath, dst_fpath):
    """
    Write a gzip-compressed copy of a file.

    Args:
        src_fpath (str): File to read.
        dst_fpath (str): Path to write the compressed copy to; overwritten if it already exists.

    Returns:
        None
    """
    with open(src_fpath, "rb") as f_in:
        with gzip.open(dst_fpath, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def safe_filename(name, max_len=80):
    """
    Build a filesystem-safe filename stem from a pocket_id.

    Args:
        name (str): The pocket_id, or any raw input string, which may be a path.
        max_len (int): Ceiling on the result's length. Defaults to 80.

    Returns:
        str: The basename of `name` with every character outside [A-Za-z0-9._-] replaced by "_",
            truncated, and suffixed with an md5 of the full `name` so distinct names never collide.
            At most max_len characters (or 32, if max_len leaves no room).
    """
    # Basename only: a stem with a directory in it would resolve outside the directory it is joined onto
    name_hash = hashlib.md5(name.encode()).hexdigest()
    stem = UNSAFE_FILENAME_CHARS.sub("_", os.path.basename(name.rstrip("/")))
    stem = stem[: max(max_len - len(name_hash) - 1, 0)].rstrip("_")
    return f"{stem}_{name_hash}" if stem else name_hash


def jsonify_dict(item):
    """
    Recursively turn sets into lists so a dict becomes JSON-serialisable.

    Args:
        item: Any value; dicts are walked and their keys coerced to strings.

    Returns:
        The same structure with every set replaced by a list.
    """
    if isinstance(item, set):
        return list(item)
    elif isinstance(item, dict):
        return {str(k): jsonify_dict(v) for k, v in item.items()}
    else:
        return item
