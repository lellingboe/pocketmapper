"""
Step 2: download the structures and bundled Foldseek database the records need.

A Foldseek PDB database's hits are not known yet; the align step fetches theirs. PISA interfaces are
fetched by the pockets step.
"""

import logging
import os

from pocketmapper.downloads.structure_downloader import StructureDownloader
from pocketmapper.exceptions import PocketMapperError
from pocketmapper.foldseek import run_foldseek
from pocketmapper.lib import log_to_file
from pocketmapper.lib import make_dir
from pocketmapper.lib import temp_dir_scope
from pocketmapper.records import append_failed_entries
from pocketmapper.records import failed_entry
from pocketmapper.records import read_cache_manifest
from pocketmapper.records import read_records
from pocketmapper.records import unique_by
from pocketmapper.records import write_records
from pocketmapper.settings import input_path
from pocketmapper.settings import layer_settings
from pocketmapper.settings import require_foldseek
from pocketmapper.settings import require_setting
from pocketmapper.settings import resolve_delete_tmp
from pocketmapper.settings import resolve_paths
from pocketmapper.settings import resolve_threads

logger = logging.getLogger(__name__)


def fetch_structures(
    job_file=None,
    results_dir=None,
    work_dir=None,
    verbosity=None,
    log_path=None,
    failed_entries_path=None,
    query_records=None,
    target_records=None,
    query_records_path=None,
    target_records_path=None,
    threads=None,
    temp_dir=None,
    delete_tmp=None,
):
    """
    Download the structures and Foldseek database the records need.

    Reads the cache directories from `results_dir`'s cache manifest. Drops the records whose
    structure cannot be fetched, adding them to `failed_entries_path`. Empties `temp_dir` on the way
    in and, unless `delete_tmp` is 0, deletes it on the way out, unless an enclosing call holds it.

    Args:
        job_file (str or dict, optional): JSON job file of job key -> value, or the same
            already loaded. Any argument given overrides it.
        results_dir (str, optional): The results directory `parse` wrote to. Required here or in
            `job_file`.
        work_dir (str, optional): Directory that entries and relative paths resolve against.
            Defaults to the working directory.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Defaults to <results_dir>/info.log.
        failed_entries_path (str, optional): Defaults to <results_dir>/failed_entries.json.
        query_records (str, optional): Defaults to the job file's query_records_path, else
            <results_dir>/query_records.json.
        target_records (str, optional): As `query_records`, for the target side.
        query_records_path (str, optional): Where the query records left are written. Defaults to
            `query_records`.
        target_records_path (str, optional): As `query_records_path`. Defaults to `target_records`.
        threads (int, optional): Defaults to one per available core.
        temp_dir (str, optional): Defaults to <results_dir>/tmp.
        delete_tmp (int, optional): 1 deletes `temp_dir` at the end; 0 keeps it. Defaults to
            DEFAULT_DELETE_TMP.

    Returns:
        None

    Raises:
        PocketMapperError: If the job file cannot be read, no results_dir is given, the manifest or a
            records file is missing, a setting is invalid, no structure for a side could be fetched,
            or a Foldseek database cannot be downloaded.
    """
    # The records paths rewrite the inputs, so they are not layered: an explicit one would become
    # its input's default too
    values = layer_settings(
        job_file,
        {
            "results_dir": results_dir,
            "work_dir": work_dir,
            "verbosity": verbosity,
            "log_path": log_path,
            "failed_entries_path": failed_entries_path,
            "threads": threads,
            "temp_dir": temp_dir,
            "delete_tmp": delete_tmp,
        },
    )
    require_setting(values, "results_dir")
    values = resolve_paths(values)
    query_records = input_path(values, query_records, "query_records_path")
    target_records = input_path(values, target_records, "target_records_path")
    with log_to_file(values["log_path"], values["verbosity"]):
        cache_dirs = read_cache_manifest(values["results_dir"])
        threads = resolve_threads(values["threads"])
        delete_tmp = resolve_delete_tmp(values["delete_tmp"])

        roots = [cache_dirs["cache_dir"], values["results_dir"]]
        with temp_dir_scope(values["temp_dir"], delete_tmp, roots):
            fetch_inputs(
                query_records,
                target_records,
                query_records_path if query_records_path is not None else query_records,
                target_records_path if target_records_path is not None else target_records,
                values["failed_entries_path"],
                threads,
                values["temp_dir"],
            )


def fetch_inputs(
    query_records,
    target_records,
    query_records_path,
    target_records_path,
    failed_entries_path,
    threads,
    temp_dir,
):
    """
    Download the structures and Foldseek database both sides' records need.

    Writes each structure to its record's `struct_path`. Records whose structure cannot be fetched
    are dropped and added to `failed_entries_path` as `structure_not_found`.

    Args:
        query_records (str): The query records file.
        target_records (str): The target records file.
        query_records_path (str): Where the query records left are written. May be `query_records`.
        target_records_path (str): As `query_records_path`, for the target side.
        failed_entries_path (str): The failed-entries file, appended to.
        threads (int): Thread count for a Foldseek database download.
        temp_dir (str): Scratch directory; a Foldseek database download works under it.

    Returns:
        None

    Raises:
        PocketMapperError: If a records file cannot be read, no structure for a side could be fetched,
            or a Foldseek database cannot be downloaded.
    """
    log_extra = {"stage": "Downloading Structures"}

    kept = {}
    for name, in_path in (("query", query_records), ("target", target_records)):
        records = read_records(in_path)
        for record in records:
            if record["struct_type"] == "foldseek_db":
                fetch_missing_fsdb(record, threads, os.path.join(temp_dir, "foldseek_tmp"))
        structure_records = [record for record in records if record["struct_type"] != "foldseek_db"]
        found = fetch_missing_structures(name, structure_records) if structure_records else {}

        failures = [
            failed_entry(record["pocket_id"], "fetch_structures", "structure_not_found", in_path, record)
            for record in structure_records
            if not found[record["struct_info"]]
        ]
        append_failed_entries(failed_entries_path, failures)
        kept[name] = [
            record for record in records if record["struct_type"] == "foldseek_db" or found[record["struct_info"]]
        ]
        if not kept[name]:
            logger.critical(f"Insufficient {name} structures after fetching", extra=log_extra)
            raise PocketMapperError(f"Insufficient {name} structures after fetching. No valid {name} entries remain.")

    write_records(kept["query"], query_records_path)
    write_records(kept["target"], target_records_path)


def fetch_missing_structures(name, records):
    """
    Download the reference structures a list of records needs.

    Writes each structure to its record's `struct_path`, unless a file is already there.

    Args:
        name (str): Which records these are, e.g. "query" or "foldseek hit". Used in logging.
        records (list): QTRecord dicts.

    Returns:
        dict: struct_info -> whether its structure is available.
    """
    log_extra = {"stage": "Downloading Structures"}

    # One download per structure, however many pockets sit on it
    found = StructureDownloader().download_missing_structures(unique_by(records, "struct_info"))
    logger.debug(f"Structure fetcher results: {found}", extra=log_extra)

    logger.info(f"{sum(found.values())}/{len(found)} {name} required structures available", extra=log_extra)
    missing = list(dict.fromkeys(record["pocket_id"] for record in records if not found[record["struct_info"]]))
    if missing:
        logger.warning(f"Missing structures for {name}(s): {', '.join(missing)}", extra=log_extra)
    return found


def fetch_missing_fsdb(record, threads, tmp_dir):
    """
    Download a bundled Foldseek database if it is not already on disk.

    Probes the foldseek binary before downloading.

    Args:
        record (dict): The `foldseek_db` record; the database is named by its `struct_info` and
            written to its `struct_path`.
        threads (int): Thread count for the download.
        tmp_dir (str): Scratch directory for `foldseek databases`, created if needed.

    Returns:
        None

    Raises:
        PocketMapperError: If foldseek cannot run, the destination or `tmp_dir` cannot be created, or
            -- from `run_foldseek` -- if the download fails.
    """
    log_extra = {"stage": "Fetching Missing Foldseek Database"}
    fsdb_name = record["struct_info"].upper()
    fsdb_path = record["struct_path"]
    if os.path.exists(fsdb_path):
        return

    require_foldseek(f"Downloading the Foldseek database '{fsdb_name}' needs foldseek")
    logger.info(f"Fetching bundled Foldseek database '{fsdb_name}' to {fsdb_path}", extra=log_extra)
    make_dir(os.path.dirname(fsdb_path), log_extra)
    make_dir(tmp_dir, log_extra)
    run_foldseek(["databases", fsdb_name, fsdb_path, tmp_dir, "--threads", str(threads)], log_extra)
    logger.info(f"Successfully fetched Foldseek database '{fsdb_name}'", extra=log_extra)
