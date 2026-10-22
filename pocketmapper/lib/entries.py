"""
What the steps share about parsed entries: the entries table and helpers over QTRecord dicts.

`parse` writes the entries table, `entries.tsv`: a row per query and target entry, holding its QTRecord
fields, the input it came from and, once a step has dropped it, which step, why and how. Every later step
reads its records from the table (`open_entries`) and marks the ones it drops (`report_failures`), so no
step parses the entries again or meets a failure an earlier step already reported.

A step clears its own and every later step's marks when it opens the table, so rerunning a step gives the
table it gave the first time.
"""

import logging
import os

import pandas as pd

from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.files import make_dir
from pocketmapper.lib.scopes import outer_command

logger = logging.getLogger(__name__)

# The fields of `parse.qt_processor.QTRecord`, in its order; a record dict holds exactly these
RECORD_FIELDS = (
    "pocket_id",
    "struct_info",
    "chain_info",
    "residue_info",
    "struct_type",
    "struct_path",
    "preprocess_name",
    "pocket_method",
)

# The entries table's columns. The failed_* columns are empty until a step drops the entry.
ENTRY_COLUMNS = ("side",) + RECORD_FIELDS + ("source", "failed_step", "failed_reason", "failed_detail")

# The steps in the order they run; a step clears the marks of itself and every step after it
STEP_ORDER = ("parse", "fetch_structures", "align", "pockets", "compare", "superpose")


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


def entry_row(side, record, source, failed_step=None, failed_reason=None, failed_detail=None):
    """
    Build one row of the entries table.

    Args:
        side (str): "query" or "target".
        record (dict): A QTRecord dict, or for an entry that never became a record, one holding only
            its `pocket_id`.
        source (str): The input or entries file the entry was read from.
        failed_step (str, optional): The step that dropped it. Defaults to None, for not dropped.
        failed_reason (str, optional): Why, e.g. "invalid_entry". Defaults to None.
        failed_detail (str, optional): A human-readable explanation. Defaults to None.

    Returns:
        dict: Column -> value, for every one of ENTRY_COLUMNS.
    """
    row = {"side": side}
    row |= {key: record.get(key) for key in RECORD_FIELDS}
    row |= {
        "source": source,
        "failed_step": failed_step,
        "failed_reason": failed_reason,
        "failed_detail": failed_detail,
    }
    return row


def write_entries(rows, path):
    """
    Write the entries table, creating its directory.

    Args:
        rows (list): Rows from `entry_row`.
        path (str or None): The table, tab-separated, or None to write none.

    Returns:
        None

    Raises:
        PocketMapperError: If the directory cannot be created.
    """
    if path is None:
        return
    make_dir(os.path.dirname(path), {"stage": "Writing Files"})
    pd.DataFrame(rows, columns=list(ENTRY_COLUMNS)).to_csv(path, sep="\t", index=False)


def read_entries(path):
    """
    Read the entries table.

    Args:
        path (str): The table.

    Returns:
        list: A row dict per entry, in table order, None for each empty cell.

    Raises:
        PocketMapperError: If the file is missing, unreadable or lacks a column.
    """
    require_file(path, "entries")
    try:
        table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except ValueError as e:
        msg = f"Could not read the entries file {path}: {e}"
        logger.critical(msg)
        raise PocketMapperError(msg) from e
    missing = [column for column in ENTRY_COLUMNS if column not in table.columns]
    if missing:
        msg = f"The entries file {path} lacks the column(s) {', '.join(missing)}; rerun parse to rewrite it"
        logger.critical(msg)
        raise PocketMapperError(msg)
    return [
        {column: (row[column] if row[column] != "" else None) for column in ENTRY_COLUMNS}
        for row in table.to_dict("records")
    ]


def entry_sides(rows):
    """
    The records of the entries no step has dropped, by side.

    Args:
        rows (list): Rows of the entries table.

    Returns:
        dict: "query" and "target" -> that side's QTRecord dicts, in table order.
    """
    sides = {"query": [], "target": []}
    for row in rows:
        if row["failed_step"] is None:
            sides[row["side"]].append({key: row[key] for key in RECORD_FIELDS})
    return sides


def open_entries(path, step):
    """
    Read the records a step works on from the entries table.

    Clears the marks of `step` and of every step after it, rewriting the table if any was set, so a
    rerun starts from what the earlier steps left.

    Args:
        path (str): The entries table.
        step (str): The step opening it, one of STEP_ORDER.

    Returns:
        dict: "query" and "target" -> the QTRecord dicts of the entries no earlier step dropped.

    Raises:
        PocketMapperError: If the table cannot be read, or a side has no entry left.
    """
    log_extra = {"stage": "Reading Entries"}

    rows = read_entries(path)
    cleared = STEP_ORDER[STEP_ORDER.index(step) :]
    reset = [row for row in rows if row["failed_step"] in cleared]
    for row in reset:
        row["failed_step"] = row["failed_reason"] = row["failed_detail"] = None
    if reset:
        write_entries(rows, path)

    sides = entry_sides(rows)
    errors = []
    for name, records in sides.items():
        if not records:
            logger.critical(f"No {name} entries left in {path}", extra=log_extra)
            errors.append(f"no {name} entries left")
    if errors:
        raise PocketMapperError("; ".join(errors))
    return sides


def failure(side, record, reason, detail=None):
    """
    Describe one entry a step drops.

    Args:
        side (str or None): The record's side, or None for a record that is no entry (a Foldseek hit),
            which is logged but not marked.
        record (dict): The QTRecord dict.
        reason (str): Why, e.g. "structure_not_found".
        detail (str, optional): A human-readable explanation. Defaults to None.

    Returns:
        dict: "side", "pocket_id", "reason" and "detail".
    """
    return {"side": side, "pocket_id": record["pocket_id"], "reason": reason, "detail": detail}


def report_failures(path, step, failures, log_extra, message=None):
    """
    Log the entries a step drops as warnings, and mark them in the entries table.

    Each entry's rows are marked only if no step has marked them yet, so the first reason stands.

    Args:
        path (str or None): The entries table, or None to mark nothing.
        step (str): The step dropping them.
        failures (list): Entries from `failure`. Nothing is written when empty.
        log_extra (dict): Logging `extra` for the messages.
        message (str, optional): What the entries have in common, logged once ahead of their
            pocket_ids. Defaults to None, which logs each entry's own `detail`.

    Returns:
        None

    Raises:
        PocketMapperError: If the table cannot be read.
    """
    if not failures:
        return
    if message is None:
        for entry in failures:
            logger.warning(f"{entry['detail']}; skipping this entry", extra=log_extra)
    else:
        pocket_ids = ", ".join(dict.fromkeys(entry["pocket_id"] for entry in failures))
        logger.warning(f"{message}: {pocket_ids}", extra=log_extra)

    marks = {(entry["side"], entry["pocket_id"]): entry for entry in failures if entry["side"] is not None}
    if path is None or not marks:
        return
    rows = read_entries(path)
    for row in rows:
        entry = marks.get((row["side"], row["pocket_id"]))
        if entry is not None and row["failed_step"] is None:
            row["failed_step"] = step
            row["failed_reason"] = entry["reason"]
            row["failed_detail"] = entry["detail"]
    write_entries(rows, path)
    logger.info(f"{len(marks)} entries skipped; see {path}", extra=log_extra)


def split_missing_structures(records, side, step, entries_path, log_extra):
    """
    Set aside the records whose structure, or Foldseek database, is not on disk.

    Reports them through `report_failures` as `structure_not_found`, with a hint to run
    fetch_structures (or pass --fetch_missing 1) unless running inside search, which already has.

    Args:
        records (list): One side's QTRecord dicts.
        side (str): Which side they are, "query" or "target".
        step (str): The step skipping them.
        entries_path (str): The entries table, marked.
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
            failures.append(failure(side, record, "structure_not_found"))
    hint = "" if outer_command() == "search" else "; run fetch_structures first, or pass --fetch_missing 1"
    report_failures(entries_path, step, failures, log_extra, f"No structure on disk, skipping{hint}")
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
