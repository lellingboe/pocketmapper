"""
Step 1: parse the query and target inputs, checking each entry against its pocket method.

A dry run: no network. Writes the entries table every later step reads its records from (see
`lib.entries`), and its settings dump is the job file the rest of a chain takes.
"""

import logging

import pandas as pd

from pocketmapper.lib.entries import ENTRY_COLUMNS, entry_row, entry_sides, fsdb_record, write_entries
from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.logs import log_to_file
from pocketmapper.lib.scopes import run_scope
from pocketmapper.lib.settings import dump_settings, layer_settings, require_setting, resolve_paths
from pocketmapper.parse.qt_processor import QTProcessor

logger = logging.getLogger(__name__)


def parse(
    query=None,
    target=None,
    job_file=None,
    results_dir=None,
    work_dir=None,
    verbosity=None,
    log_path=None,
    job_settings_path=None,
    query_pocket_method=None,
    target_pocket_method=None,
    cache_dir=None,
    pdb_dir=None,
    alphafold_dir=None,
    pocket_dir=None,
    foldseek_preprocessed_structure_dir=None,
    fsdb_dir=None,
    entries_path=None,
):
    """
    Parse the query and target inputs, checking each entry against its pocket method. No network.

    Writes the entries table afresh, the entries that could not be parsed marked as dropped, even when
    the entries are then rejected. Without a `results_dir`, writes no log, settings or entries table
    unless that file's own path is given. Logs nothing about the table it returns: printing it is the
    caller's business.

    Args:
        query (str, optional): Query entry, or a file of one entry per line. Required here or in
            `job_file`, not both.
        target (str, optional): Target entry, a file of them, or a Foldseek database. As `query`.
        job_file (str or dict, optional): JSON job file of job key -> value, or the same
            already loaded. Any argument given overrides it.
        results_dir (str, optional): Where the log, the settings and the entries table go by default.
            Defaults to None, for none.
        work_dir (str, optional): Directory that entries and relative paths resolve against.
            Defaults to the working directory.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Defaults to <results_dir>/info.log.
        job_settings_path (str, optional): Where these settings are written, as a job file for later
            steps. Defaults to <results_dir>/parse_settings.json. Not written when run inside search.
        query_pocket_method (str, optional): Pocket method to force on every query entry, or "auto"
            (the default) to infer it per entry.
        target_pocket_method (str, optional): As `query_pocket_method`, for the target side.
        cache_dir (str, optional): Defaults to DEFAULT_CACHE_DIR.
        pdb_dir (str, optional): Defaults to <cache_dir>/pdb_structures.
        alphafold_dir (str, optional): Defaults to <cache_dir>/alphafold_structures.
        pocket_dir (str, optional): Defaults to <cache_dir>/pockets.
        foldseek_preprocessed_structure_dir (str, optional): Defaults to
            <cache_dir>/foldseek_preprocessed_structures.
        fsdb_dir (str, optional): Defaults to <cache_dir>/fsdb.
        entries_path (str, optional): Where the entries table is written. Defaults to
            <results_dir>/entries.tsv.

    Returns:
        pandas.DataFrame: The entries table, a row per entry with the columns of
            `lib.entries.ENTRY_COLUMNS`.

    Raises:
        PocketMapperError: If the job file cannot be read, query or target is given both ways or
            neither way, a pocket method is unknown, a side has no valid entries, or a
            Foldseek-database target is not the only target entry.
    """
    values = layer_settings(
        job_file,
        {
            "query": query,
            "target": target,
            "results_dir": results_dir,
            "work_dir": work_dir,
            "verbosity": verbosity,
            "log_path": log_path,
            "job_settings_path": job_settings_path,
            "query_pocket_method": query_pocket_method,
            "target_pocket_method": target_pocket_method,
            "cache_dir": cache_dir,
            "pdb_dir": pdb_dir,
            "alphafold_dir": alphafold_dir,
            "pocket_dir": pocket_dir,
            "foldseek_preprocessed_structure_dir": foldseek_preprocessed_structure_dir,
            "fsdb_dir": fsdb_dir,
            "entries_path": entries_path,
        },
    )
    for key in ("query", "target"):
        require_setting(values, key)
    values = resolve_paths(values, "parse")
    with run_scope("parse") as outermost, log_to_file(values["log_path"], values["verbosity"]):
        if outermost:
            dump_settings(values)
        _, rows = parse_entries(
            values["query"],
            values["target"],
            values["query_pocket_method"],
            values["target_pocket_method"],
            values,
            values["work_dir"],
            values["entries_path"],
        )
        return pd.DataFrame(rows, columns=list(ENTRY_COLUMNS))


def parse_entries(
    query,
    target,
    query_pocket_method,
    target_pocket_method,
    cache_dirs,
    work_dir,
    entries_path=None,
):
    """
    Parse both sides' entries into records, and check that the two sides make a search.

    Every entry that could not be parsed, and every Foldseek-database query entry, is marked as dropped,
    `invalid_entry`, and logged as a warning. The entries table is written before the sides are
    checked, so it lists them even when the check fails.

    Args:
        query (str or list): Query entry, a file of one entry per line, or a list of either.
        target (str or list): Target entry, a file of them, a Foldseek database, or a list of these.
        query_pocket_method (str): Pocket method to force on every query entry, or "auto".
        target_pocket_method (str): As `query_pocket_method`, for the target side.
        cache_dirs (dict): Holds "pdb_dir", "alphafold_dir" and "fsdb_dir", which give the
            records' structure paths.
        work_dir (str): Directory that entries files, local structure files, a user Foldseek
            database and relative cache directories resolve against.
        entries_path (str, optional): Where the entries table is written. Defaults to None, for
            nowhere.

    Returns:
        tuple: (sides, rows). `sides` is "query" and "target" -> the QTRecord dicts of that side's
            entries kept, in input order; `rows` the entries table, as `lib.entries.entry_row` builds
            it, every entry included.

    Raises:
        PocketMapperError: If a pocket method is unknown, an entries file cannot be read, either side
            has no valid entries, or a Foldseek-database target is not the only target entry.
    """
    log_extra = {"stage": "Determine Query/Target Types"}

    rows = []
    for name, entries, pocket_method in (
        ("query", query, query_pocket_method),
        ("target", target, target_pocket_method),
    ):
        rows += parse_side(entries, name, pocket_method, cache_dirs, work_dir)
    for row in rows:
        # A database can only be searched
        if row["side"] == "query" and row["struct_type"] == "foldseek_db" and row["failed_step"] is None:
            row["failed_step"] = "parse"
            row["failed_reason"] = "invalid_entry"
            row["failed_detail"] = f"A Foldseek database cannot be a query entry: {row['pocket_id']}"
    for row in rows:
        if row["failed_step"] is not None:
            logger.warning(f"{row['failed_detail']}; skipping this entry", extra=log_extra)
    write_entries(rows, entries_path)

    sides = entry_sides(rows)
    logger.info(f"Parsed {len(sides['query'])} query and {len(sides['target'])} target entries", extra=log_extra)
    errors = []
    for name, records in sides.items():
        if not records:
            logger.critical(f"No valid {name} entries after processing", extra=log_extra)
            errors.append(f"no valid {name} entries")
    if errors:
        raise PocketMapperError("; ".join(errors))

    if fsdb_record(sides["target"]) is not None and len(sides["target"]) > 1:
        msg = (
            "A Foldseek database target must be the only target entry, but "
            f"{len(sides['target'])} target entries were given"
        )
        logger.critical(msg, extra=log_extra)
        raise PocketMapperError(msg)
    return sides, rows


def parse_side(entries, name, pocket_method, cache_dirs, work_dir):
    """
    Parse one side's entries into rows of the entries table.

    Args:
        entries (str or list): An entry, a file of one entry per line, or a list of either.
        name (str): Which side this is, "query" or "target".
        pocket_method (str): Pocket method to force on every entry, or "auto".
        cache_dirs (dict): Holds "pdb_dir", "alphafold_dir" and "fsdb_dir", which give the
            records' structure paths.
        work_dir (str): Directory that entries files, local structure files, a user Foldseek
            database and relative cache directories resolve against.

    Returns:
        list: A row per entry, in input order; one that could not be parsed is marked as dropped by
            parse, `invalid_entry`, its record fields empty.

    Raises:
        PocketMapperError: If `pocket_method` is unknown or an entries file cannot be read.
    """
    qtprocessor = QTProcessor(
        pdb_dir=cache_dirs["pdb_dir"],
        alphafold_dir=cache_dirs["alphafold_dir"],
        fsdb_dir=cache_dirs["fsdb_dir"],
        work_dir=work_dir,
    )
    rows = []
    for source in [entries] if isinstance(entries, str) else entries:
        results = qtprocessor.process_qt_cmdline_input(qt_input=source, name=name, pocket_method=pocket_method)
        for entry, record, reason in results:
            if record is None:
                rows.append(entry_row(name, {"pocket_id": entry}, source, "parse", "invalid_entry", reason))
            else:
                rows.append(entry_row(name, record, source))
    return rows
