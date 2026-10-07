"""
The files steps hand each other: records, the cache manifest and the failed-entries log.

A records file is a JSON list of `QTRecord` dicts. The cache manifest, `cache_dirs.json` in a
results directory, names the cache directories a chain of steps shares. `failed_entries.json`
collects every record a step dropped, with the reason.

Also derives what a target records file says about the run: whether the target is a Foldseek
database and, if so, how its pockets are resolved.
"""

import json
import logging
import os

from pocketmapper.exceptions import PocketMapperError
from pocketmapper.lib import make_dir

logger = logging.getLogger(__name__)

# The cache manifest's fixed name within a results directory, and the directories it names
CACHE_MANIFEST_NAME = "cache_dirs.json"
CACHE_MANIFEST_KEYS = (
    "cache_dir",
    "pdb_dir",
    "alphafold_dir",
    "pocket_dir",
    "foldseek_preprocessed_structure_dir",
    "fsdb_dir",
)


def require_file(path, what):
    """
    Check that a file a step reads exists.

    Args:
        path (str): The file.
        what (str): What the file holds, for the error message, e.g. "alignment".

    Returns:
        None

    Raises:
        PocketMapperError: If there is no file at `path`.
    """
    if not os.path.isfile(path):
        msg = f"No {what} file at {path}"
        logger.critical(msg)
        raise PocketMapperError(msg)


def read_json(path, what):
    """
    Read a JSON file a step needs as input.

    Args:
        path (str): The file.
        what (str): What the file holds, for the error message, e.g. "query records".

    Returns:
        The parsed JSON.

    Raises:
        PocketMapperError: If the file is missing or is not JSON.
    """
    require_file(path, what)
    try:
        with open(path) as f:
            return json.load(f)
    except ValueError as e:
        msg = f"Could not read the {what} file {path}: {e}"
        logger.critical(msg)
        raise PocketMapperError(msg) from e


def write_json(value, path):
    """
    Write a value as JSON, creating the file's directory.

    Args:
        value: Anything `json.dump` accepts.
        path (str): The file to write.

    Returns:
        None

    Raises:
        PocketMapperError: If the directory cannot be created.
    """
    make_dir(os.path.dirname(path), {"stage": "Writing Files"})
    with open(path, "w") as f:
        json.dump(value, f, indent=4)


def read_records(path):
    """
    Read a records file.

    Args:
        path (str): The records file.

    Returns:
        list: QTRecord dicts, in file order.

    Raises:
        PocketMapperError: If the file is missing or is not JSON.
    """
    return read_json(path, "records")


def write_records(records, path):
    """
    Write a records file, creating its directory.

    Args:
        records (list): QTRecord dicts.
        path (str): The records file.

    Returns:
        None
    """
    write_json(records, path)
    logger.debug(f"Wrote {len(records)} records to {path}")


def write_cache_manifest(results_dir, cache_dirs):
    """
    Write `cache_dirs.json` into a results directory.

    Args:
        results_dir (str): The results directory.
        cache_dirs (dict): Each of CACHE_MANIFEST_KEYS -> its directory. Stored absolute.

    Returns:
        None
    """
    write_json({key: os.path.abspath(cache_dirs[key]) for key in CACHE_MANIFEST_KEYS}, cache_manifest_path(results_dir))


def read_cache_manifest(results_dir):
    """
    Read `cache_dirs.json` from a results directory.

    Args:
        results_dir (str): The results directory.

    Returns:
        dict: Each of CACHE_MANIFEST_KEYS -> its absolute directory.

    Raises:
        PocketMapperError: If the manifest is missing, unreadable or lacks a directory.
    """
    path = cache_manifest_path(results_dir)
    manifest = read_json(path, "cache manifest")
    missing = [key for key in CACHE_MANIFEST_KEYS if key not in manifest]
    if missing:
        msg = f"Cache manifest {path} does not name {', '.join(missing)}; rerun parse"
        logger.critical(msg)
        raise PocketMapperError(msg)
    return manifest


def cache_manifest_path(results_dir):
    """
    Where a results directory's cache manifest lives.

    Args:
        results_dir (str): The results directory.

    Returns:
        str: The path of `cache_dirs.json` in it.
    """
    return os.path.join(results_dir, CACHE_MANIFEST_NAME)


def failed_entry(pocket_id, step, reason, source, record=None, detail=None):
    """
    Build one `failed_entries.json` entry.

    Args:
        pocket_id (str): The entry as typed.
        step (str): The step that dropped it, e.g. "fetch".
        reason (str): Why, e.g. "structure_not_found".
        source (str): The records file the record came from, or the query/target input it was
            parsed from.
        record (dict, optional): The dropped record, whose fields are appended. Defaults to None,
            for an entry that never became a record.
        detail (str, optional): A human-readable explanation. Defaults to None, which adds none.

    Returns:
        dict: The entry.
    """
    entry = {"pocket_id": pocket_id, "step": step, "reason": reason, "source": source}
    if detail is not None:
        entry["detail"] = detail
    if record is not None:
        entry |= {key: value for key, value in record.items() if key not in entry}
    return entry


def start_failed_entries(path):
    """
    Start a fresh `failed_entries.json`, holding an empty list.

    Args:
        path (str): The failed-entries file.

    Returns:
        None
    """
    write_json([], path)


def append_failed_entries(path, entries):
    """
    Add entries to `failed_entries.json`, creating it if missing.

    Args:
        path (str): The failed-entries file.
        entries (list): Entries from `failed_entry`. Nothing is written when empty.

    Returns:
        None
    """
    if not entries:
        return
    existing = read_json(path, "failed entries") if os.path.isfile(path) else []
    write_json(existing + list(entries), path)
    logger.info(f"{len(entries)} entries dropped; see {path}")


def unique_by(records, *fields):
    """
    Keep the first record for each distinct combination of some fields.

    Args:
        records (list): QTRecord dicts.
        *fields (str): The fields to tell records apart by.

    Returns:
        list: The first record carrying each combination, in record order.
    """
    unique = {}
    for record in records:
        unique.setdefault(tuple(record[field] for field in fields), record)
    return list(unique.values())


def fsdb_record(target_records):
    """
    The Foldseek-database record among the target records, if the target is a database.

    Args:
        target_records (list): Target QTRecord dicts.

    Returns:
        dict: The `foldseek_db` record, or None when the targets are structures.
    """
    return next((record for record in target_records if record["struct_type"] == "foldseek_db"), None)


def synthesise_target_pockets(target_records):
    """
    Whether target pockets must be synthesised from the alignment rather than looked up.

    Args:
        target_records (list): Target QTRecord dicts.

    Returns:
        bool: True for a Foldseek-database target that was not expanded into PISA records.
    """
    record = fsdb_record(target_records)
    return record is not None and record["fsdb_pockets"] != "pisa"


def preproc_to_ids(records):
    """
    Map each aligned chain to the pockets on it.

    Args:
        records (list): QTRecord dicts, query and target together.

    Returns:
        dict: preprocess_name -> the pocket_ids on that chain, in record order, without repeats.
    """
    # Keyed by preprocess_name: one chain can carry several pockets (several pocket_ids), so
    # membership must be tested on the key, not on the pocket_id, or each chain keeps only its last.
    mapping = {}
    for record in records:
        pocket_ids = mapping.setdefault(record["preprocess_name"], [])
        if record["pocket_id"] not in pocket_ids:
            pocket_ids.append(record["pocket_id"])
    return mapping
