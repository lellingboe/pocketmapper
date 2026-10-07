"""
What the steps share about parsed entries and the ones they skip: the failed-entries log and helpers
over QTRecord dicts.

No step hands another its records: each re-derives them from the entries
(`steps.parse.parse_entries`). `failed_entries.json` collects every entry a step left out, with the
reason; nothing reads it back.

Also derives what the target records say about the run: whether the target is a Foldseek database.
"""

import json
import logging
import os

from pocketmapper.exceptions import PocketMapperError
from pocketmapper.lib import make_dir
from pocketmapper.lib import outer_command

logger = logging.getLogger(__name__)


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


def failed_entry(pocket_id, step, reason, source, record=None, detail=None):
    """
    Build one `failed_entries.json` entry.

    Args:
        pocket_id (str): The entry as typed.
        step (str): The step that dropped it, e.g. "fetch_structures".
        reason (str): Why, e.g. "structure_not_found".
        source (str): The query/target input the entry was parsed from, or the file it was read
            from.
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
    Add entries to `failed_entries.json`, creating it if missing, skipping any already listed.

    An entry is already listed when the file, or an earlier entry in `entries`, holds its
    `(pocket_id, reason)`: a chain of steps that each re-parse the entries meets each failure once per
    step.

    Args:
        path (str): The failed-entries file.
        entries (list): Entries from `failed_entry`. Nothing is written when none is new.

    Returns:
        list: The entries that were new, in order.
    """
    existing = read_json(path, "failed entries") if os.path.isfile(path) else []
    new = unlisted(entries, existing)
    if new:
        write_json(existing + new, path)
        logger.info(f"{len(new)} entries skipped; see {path}")
    return new


def unlisted(entries, listed):
    """
    The failure entries not already listed.

    Args:
        entries (list): Entries from `failed_entry`.
        listed (list): Entries already recorded.

    Returns:
        list: Each entry of `entries` whose `(pocket_id, reason)` is neither in `listed` nor on an
            earlier entry of `entries`, in order.
    """
    seen = {(entry["pocket_id"], entry["reason"]) for entry in listed}
    new = []
    for entry in entries:
        key = (entry["pocket_id"], entry["reason"])
        if key not in seen:
            seen.add(key)
            new.append(entry)
    return new


def report_failures(path, entries, log_extra, message=None):
    """
    Add entries to the failed-entries file and log them: new ones as warnings, ones already listed at DEBUG.

    Args:
        path (str or None): The failed-entries file, or None to write none, every entry then being
            new unless it repeats an earlier one.
        entries (list): Entries from `failed_entry`.
        log_extra (dict): Logging `extra` for the messages.
        message (str, optional): What the entries have in common, logged once ahead of their
            pocket_ids. Defaults to None, which logs each entry's own `detail`.

    Returns:
        list: The entries that were new, in order.
    """
    new = append_failed_entries(path, entries) if path is not None else unlisted(entries, [])
    new_ids = {id(entry) for entry in new}
    repeats = [entry for entry in entries if id(entry) not in new_ids]
    for level, group, suffix in (
        (logging.WARNING, new, ""),
        (logging.DEBUG, repeats, " (already listed as failed)"),
    ):
        if not group:
            continue
        if message is None:
            for entry in group:
                logger.log(level, f"{entry['detail']}; skipping this entry{suffix}", extra=log_extra)
        else:
            pocket_ids = ", ".join(dict.fromkeys(entry["pocket_id"] for entry in group))
            logger.log(level, f"{message}{suffix}: {pocket_ids}", extra=log_extra)
    return new


def split_missing_structures(records, step, source, failed_entries_path, log_extra):
    """
    Set aside the records whose structure, or Foldseek database, is not on disk.

    Reports them through `report_failures` as `structure_not_found`, with a hint to run
    fetch_structures unless running inside search, which already has.

    Args:
        records (list): One side's QTRecord dicts.
        step (str): The step skipping them, for the failure entries.
        source (str): The input they were parsed from, for the failure entries.
        failed_entries_path (str): The failed-entries file, appended to.
        log_extra (dict): Logging `extra` for the messages.

    Returns:
        list: The records whose structure is on disk.
    """
    kept = []
    failures = []
    for record in records:
        if os.path.exists(record["struct_path"]):
            kept.append(record)
        else:
            failures.append(failed_entry(record["pocket_id"], step, "structure_not_found", source, record))
    hint = "" if outer_command() == "search" else "; run fetch_structures first"
    report_failures(failed_entries_path, failures, log_extra, f"No structure on disk, skipping{hint}")
    return kept


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
