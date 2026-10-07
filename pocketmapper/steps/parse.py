"""
Step 1: parse the query and target inputs, checking each entry against its pocket method.

No network. Every later step parses the entries again, from its own settings, through
`parse_entries`; parse is the cheap check before anything is fetched, and its settings dump is the
job file the rest of a chain takes.
"""

import logging

from pocketmapper.entries import failed_entry
from pocketmapper.entries import fsdb_record
from pocketmapper.entries import report_failures
from pocketmapper.entries import start_failed_entries
from pocketmapper.exceptions import PocketMapperError
from pocketmapper.lib import log_to_file
from pocketmapper.lib import run_scope
from pocketmapper.qt_processor import QTProcessor
from pocketmapper.settings import dump_settings
from pocketmapper.settings import layer_settings
from pocketmapper.settings import require_setting
from pocketmapper.settings import resolve_paths

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
    failed_entries_path=None,
    query_pocket_method=None,
    target_pocket_method=None,
    cache_dir=None,
    pdb_dir=None,
    alphafold_dir=None,
    pocket_dir=None,
    foldseek_preprocessed_structure_dir=None,
    fsdb_dir=None,
):
    """
    Parse the query and target inputs, checking each entry against its pocket method. No network.

    Starts `failed_entries_path` afresh with the entries that could not be parsed. Without a
    `results_dir`, writes no log, settings or failed entries unless that file's own path is given.

    Args:
        query (str, optional): Query entry, or a file of one entry per line. Required here or in
            `job_file`, not both.
        target (str, optional): Target entry, a file of them, or a Foldseek database. As `query`.
        job_file (str or dict, optional): JSON job file of job key -> value, or the same
            already loaded. Any argument given overrides it.
        results_dir (str, optional): Where the log, the settings and the failed entries go by
            default. Defaults to None, for none.
        work_dir (str, optional): Directory that entries and relative paths resolve against.
            Defaults to the working directory.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Defaults to <results_dir>/info.log.
        job_settings_path (str, optional): Where these settings are written, as a job file for later
            steps. Defaults to <results_dir>/parse_settings.json. Not written when run inside search.
        failed_entries_path (str, optional): Defaults to <results_dir>/failed_entries.json.
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

    Returns:
        None

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
            "failed_entries_path": failed_entries_path,
            "query_pocket_method": query_pocket_method,
            "target_pocket_method": target_pocket_method,
            "cache_dir": cache_dir,
            "pdb_dir": pdb_dir,
            "alphafold_dir": alphafold_dir,
            "pocket_dir": pocket_dir,
            "foldseek_preprocessed_structure_dir": foldseek_preprocessed_structure_dir,
            "fsdb_dir": fsdb_dir,
        },
    )
    for key in ("query", "target"):
        require_setting(values, key)
    values = resolve_paths(values, "parse")
    with run_scope("parse") as outermost, log_to_file(values["log_path"], values["verbosity"]):
        if outermost:
            dump_settings(values)
        parse_inputs(
            values["query"],
            values["target"],
            values["query_pocket_method"],
            values["target_pocket_method"],
            values,
            values["work_dir"],
            values["failed_entries_path"],
        )


def parse_inputs(
    query,
    target,
    query_pocket_method,
    target_pocket_method,
    cache_dirs,
    work_dir,
    failed_entries_path,
):
    """
    Parse both sides' input, starting `failed_entries_path` afresh with the entries that fail.

    Every entry that could not be parsed and every Foldseek-database query entry is listed as
    `invalid_entry`.

    Args:
        query (str): Query entry, or a file of one entry per line.
        target (str): Target entry, a file of them, or a Foldseek database.
        query_pocket_method (str): Pocket method to force on every query entry, or "auto".
        target_pocket_method (str): As `query_pocket_method`, for the target side.
        cache_dirs (dict): Holds "pdb_dir", "alphafold_dir" and "fsdb_dir".
        work_dir (str): Directory that entries files, local structure files, a user Foldseek
            database and relative cache directories resolve against.
        failed_entries_path (str or None): The failed-entries file, or None to write none.

    Returns:
        dict: "query" and "target" -> that side's QTRecord dicts, as `parse_entries` returns them.

    Raises:
        PocketMapperError: If a pocket method is unknown, either side has no valid entries, or a
            Foldseek-database target is not the only target entry.
    """
    log_extra = {"stage": "Determine Query/Target Types"}

    if failed_entries_path is not None:
        start_failed_entries(failed_entries_path)
    sides, _ = parse_entries(
        query,
        target,
        query_pocket_method,
        target_pocket_method,
        cache_dirs,
        work_dir,
        failed_entries_path=failed_entries_path,
    )
    logger.info(f"Parsed {len(sides['query'])} query and {len(sides['target'])} target entries", extra=log_extra)
    return sides


def parse_job_entries(values, step):
    """
    Parse both sides' entries as a step's settings name them, adding failures to its failed-entries file.

    Args:
        values (dict): Job key -> value, from `settings.resolve_paths`. Reads query, target, both
            pocket methods, the cache directories, `work_dir` and `failed_entries_path`.
        step (str): The step parsing them, named in the failure entries.

    Returns:
        dict: "query" and "target" -> that side's QTRecord dicts.

    Raises:
        PocketMapperError: As `parse_entries`.
    """
    sides, _ = parse_entries(
        values["query"],
        values["target"],
        values["query_pocket_method"],
        values["target_pocket_method"],
        values,
        values["work_dir"],
        step=step,
        failed_entries_path=values["failed_entries_path"],
    )
    return sides


def parse_entries(
    query,
    target,
    query_pocket_method,
    target_pocket_method,
    cache_dirs,
    work_dir,
    step="parse",
    failed_entries_path=None,
):
    """
    Parse both sides' entries into records, and check that the two sides make a search.

    Every entry that could not be parsed, and every Foldseek-database query entry, is a failure,
    `invalid_entry`, reported through `entries.report_failures` (a warning, or DEBUG for one already
    in `failed_entries_path`) before the sides are checked, so it is listed even when the check
    fails.

    Args:
        query (str or list): Query entry, a file of one entry per line, or a list of either.
        target (str or list): Target entry, a file of them, a Foldseek database, or a list of these.
        query_pocket_method (str): Pocket method to force on every query entry, or "auto".
        target_pocket_method (str): As `query_pocket_method`, for the target side.
        cache_dirs (dict): Holds "pdb_dir", "alphafold_dir" and "fsdb_dir", which give the
            records' structure paths.
        work_dir (str): Directory that entries files, local structure files, a user Foldseek
            database and relative cache directories resolve against.
        step (str, optional): The step parsing them, named in the failure entries. Defaults to
            "parse".
        failed_entries_path (str, optional): The failed-entries file, appended to. Defaults to None,
            which writes none and warns about every failure.

    Returns:
        tuple: (sides, failures). `sides` is "query" and "target" -> that side's QTRecord dicts, in
            input order; `failures` the failed-entries entries for the entries left out.

    Raises:
        PocketMapperError: If a pocket method is unknown, an entries file cannot be read, either side
            has no valid entries, or a Foldseek-database target is not the only target entry.
    """
    log_extra = {"stage": "Determine Query/Target Types"}

    sides = {}
    failures = []
    for name, entries, pocket_method in (
        ("query", query, query_pocket_method),
        ("target", target, target_pocket_method),
    ):
        sides[name], rejected = parse_side(entries, name, pocket_method, cache_dirs, work_dir)
        failures += [
            failed_entry(entry, step, "invalid_entry", source, detail=reason) for entry, reason, source in rejected
        ]

    # A database can only be searched, not searched with
    for record in [record for record in sides["query"] if record["struct_type"] == "foldseek_db"]:
        reason = f"A Foldseek database cannot be a query entry: {record['pocket_id']}"
        failures.append(failed_entry(record["pocket_id"], step, "invalid_entry", query, record, reason))
        sides["query"].remove(record)
    report_failures(failed_entries_path, failures, log_extra)

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
    return sides, failures


def parse_structure_side(entries, name, struct_type, cache_dirs, work_dir):
    """
    Resolve the structure of each of one side's entries, ignoring their chain and residue parts.

    Args:
        entries (str or list): An entry or bare structure id, a file of them, or a list of either.
        name (str): Which side this is, e.g. "query", for the log and error messages.
        struct_type (str): Structure type to force on every entry, or "auto".
        cache_dirs (dict): Holds "pdb_dir", "alphafold_dir" and "fsdb_dir", which give the
            structure paths.
        work_dir (str): Directory that entries files, local structure files, a user Foldseek
            database and relative cache directories resolve against.

    Returns:
        tuple: (structures, rejected). `structures` holds a dict per entry resolved, as
            `QTProcessor.parse_structure` returns it, in input order; `rejected` an
            (entry, reason, source) triple for each entry that could not be.

    Raises:
        PocketMapperError: If `struct_type` is unknown or an entries file cannot be read.
    """
    qtprocessor = QTProcessor(
        pdb_dir=cache_dirs["pdb_dir"],
        alphafold_dir=cache_dirs["alphafold_dir"],
        fsdb_dir=cache_dirs["fsdb_dir"],
        work_dir=work_dir,
    )
    structures = []
    rejected = []
    for source in [entries] if isinstance(entries, str) else entries:
        resolved, failed = qtprocessor.process_structure_input(source, name, struct_type)
        structures += resolved
        rejected += [(entry, reason, source) for entry, reason in failed]
    return structures, rejected


def parse_side(entries, name, pocket_method, cache_dirs, work_dir):
    """
    Parse one side's entries into records.

    Args:
        entries (str or list): An entry, a file of one entry per line, or a list of either.
        name (str): Which side this is, e.g. "query", for the log and error messages.
        pocket_method (str): Pocket method to force on every entry, or "auto".
        cache_dirs (dict): Holds "pdb_dir", "alphafold_dir" and "fsdb_dir", which give the
            records' structure paths.
        work_dir (str): Directory that entries files, local structure files, a user Foldseek
            database and relative cache directories resolve against.

    Returns:
        tuple: (records, rejected). `records` holds the QTRecord dicts parsed, in input order;
            `rejected` an (entry, reason, source) triple for each entry that could not be, `source`
            being the entry or file it came from.

    Raises:
        PocketMapperError: If `pocket_method` is unknown or an entries file cannot be read.
    """
    qtprocessor = QTProcessor(
        pdb_dir=cache_dirs["pdb_dir"],
        alphafold_dir=cache_dirs["alphafold_dir"],
        fsdb_dir=cache_dirs["fsdb_dir"],
        work_dir=work_dir,
    )
    records = []
    rejected = []
    for source in [entries] if isinstance(entries, str) else entries:
        parsed, failed = qtprocessor.process_qt_cmdline_input(qt_input=source, name=name, pocket_method=pocket_method)
        records += parsed
        rejected += [(entry, reason, source) for entry, reason in failed]
    return records, rejected
