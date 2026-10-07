"""
Step 1: parse the query and target inputs into records files.

No network. Writes the cache manifest, so later steps resolve every cache path the way this one did.
"""

import logging
import os

from pocketmapper.exceptions import PocketMapperError
from pocketmapper.lib import log_to_file
from pocketmapper.qt_processor import QTProcessor
from pocketmapper.records import CACHE_MANIFEST_KEYS
from pocketmapper.records import append_failed_entries
from pocketmapper.records import failed_entry
from pocketmapper.records import fsdb_record
from pocketmapper.records import start_failed_entries
from pocketmapper.records import write_cache_manifest
from pocketmapper.records import write_records
from pocketmapper.settings import layer_settings
from pocketmapper.settings import require_setting
from pocketmapper.settings import resolve_paths

logger = logging.getLogger(__name__)


def parse(
    query=None,
    target=None,
    job_file=None,
    results_dir=None,
    verbosity=None,
    log_path=None,
    failed_entries_path=None,
    query_pocket_method=None,
    target_pocket_method=None,
    cache_dir=None,
    pdb_dir=None,
    alphafold_dir=None,
    pocket_dir=None,
    foldseek_preprocessed_structure_dir=None,
    fsdb_dir=None,
    query_records_path=None,
    target_records_path=None,
):
    """
    Parse the query and target inputs into records files. No network.

    Writes the records files and `cache_dirs.json`, naming every cache directory absolute, and
    starts `failed_entries_path` afresh with the entries that could not be parsed.

    Args:
        query (str, optional): Query entry, or a file of one entry per line. Required here or in
            `job_file`, not both.
        target (str, optional): Target entry, a file of them, or a Foldseek database. As `query`.
        job_file (str or dict, optional): JSON job file of job key -> value, or the same
            already loaded. Any argument given overrides it.
        results_dir (str, optional): Defaults to pocketmapper_results_<YYMMDD_HHMMSS>.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Defaults to <results_dir>/info.log.
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
        query_records_path (str, optional): Defaults to <results_dir>/query_records.json.
        target_records_path (str, optional): Defaults to <results_dir>/target_records.json.

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
            "verbosity": verbosity,
            "log_path": log_path,
            "failed_entries_path": failed_entries_path,
            "query_pocket_method": query_pocket_method,
            "target_pocket_method": target_pocket_method,
            "cache_dir": cache_dir,
            "pdb_dir": pdb_dir,
            "alphafold_dir": alphafold_dir,
            "pocket_dir": pocket_dir,
            "foldseek_preprocessed_structure_dir": foldseek_preprocessed_structure_dir,
            "fsdb_dir": fsdb_dir,
            "query_records_path": query_records_path,
            "target_records_path": target_records_path,
        },
    )
    for key in ("query", "target"):
        require_setting(values, key)
    values = resolve_paths(values)
    with log_to_file(values["log_path"], values["verbosity"]):
        parse_inputs(
            values["query"],
            values["target"],
            values["query_pocket_method"],
            values["target_pocket_method"],
            {key: values[key] for key in CACHE_MANIFEST_KEYS},
            values["results_dir"],
            values["query_records_path"],
            values["target_records_path"],
            values["failed_entries_path"],
        )


def parse_inputs(
    query,
    target,
    query_pocket_method,
    target_pocket_method,
    cache_dirs,
    results_dir,
    query_records_path,
    target_records_path,
    failed_entries_path,
):
    """
    Parse both sides' input into records, and write them with the cache manifest.

    Starts `failed_entries_path` afresh, then adds every entry that could not be parsed and every
    Foldseek-database query entry, as `invalid_entry`.

    Args:
        query (str): Query entry, or a file of one entry per line.
        target (str): Target entry, a file of them, or a Foldseek database.
        query_pocket_method (str): Pocket method to force on every query entry, or "auto".
        target_pocket_method (str): As `query_pocket_method`, for the target side.
        cache_dirs (dict): Each of `records.CACHE_MANIFEST_KEYS` -> its directory. Made absolute for
            the record paths and the manifest.
        results_dir (str): Where the manifest is written.
        query_records_path (str): Where the query records are written.
        target_records_path (str): Where the target records are written.
        failed_entries_path (str): The failed-entries file.

    Returns:
        None

    Raises:
        PocketMapperError: If a pocket method is unknown, either side has no valid entries, or a
            Foldseek-database target is not the only target entry.
    """
    log_extra = {"stage": "Determine Query/Target Types"}

    start_failed_entries(failed_entries_path)
    cache_dirs = {key: os.path.abspath(path) for key, path in cache_dirs.items()}
    qtprocessor = QTProcessor(
        pdb_dir=cache_dirs["pdb_dir"],
        alphafold_dir=cache_dirs["alphafold_dir"],
        fsdb_dir=cache_dirs["fsdb_dir"],
    )

    sides = {}
    failures = []
    for name, qt_input, pocket_method in (
        ("query", query, query_pocket_method),
        ("target", target, target_pocket_method),
    ):
        records, rejected = qtprocessor.process_qt_cmdline_input(
            qt_input=qt_input, name=name, pocket_method=pocket_method
        )
        failures += [
            failed_entry(entry, "parse", "invalid_entry", qt_input, detail=reason) for entry, reason in rejected
        ]
        sides[name] = records

    # A database can only be searched, not searched with
    for record in [record for record in sides["query"] if record["struct_type"] == "foldseek_db"]:
        reason = f"A Foldseek database cannot be a query entry: {record['pocket_id']}"
        logger.warning(f"{reason}; skipping this entry", extra=log_extra)
        failures.append(failed_entry(record["pocket_id"], "parse", "invalid_entry", query, record, reason))
        sides["query"].remove(record)
    append_failed_entries(failed_entries_path, failures)

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

    write_records(sides["query"], query_records_path)
    write_records(sides["target"], target_records_path)
    write_cache_manifest(results_dir, cache_dirs)
    logger.info(f"Parsed {len(sides['query'])} query and {len(sides['target'])} target entries", extra=log_extra)
