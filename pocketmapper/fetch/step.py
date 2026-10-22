"""
Step 2: download the structures and bundled Foldseek database the entries table names.

A Foldseek PDB database's hits are not known yet; the pockets step fetches theirs. PISA interfaces are
fetched by the pockets step.
"""

import logging
import os

from pocketmapper.fetch.structure_downloader import download_structures
from pocketmapper.lib.constants import TEMP_SUBDIRS
from pocketmapper.lib.entries import failure, open_entries, report_failures, unique_by
from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.files import make_dir
from pocketmapper.lib.foldseek import run_foldseek
from pocketmapper.lib.logs import log_to_file
from pocketmapper.lib.scopes import run_scope, temp_dir_scope
from pocketmapper.lib.settings import (
    layer_settings,
    require_foldseek,
    resolve_delete_tmp,
    resolve_step_paths,
    resolve_threads,
)

logger = logging.getLogger(__name__)


def fetch_structures(
    job_file=None,
    entries=None,
    verbosity=None,
    log_path=None,
    threads=None,
    temp_dir=None,
    delete_tmp=None,
):
    """
    Download the structures and Foldseek databases the entries table names.

    A local file is checked, not fetched. Marks the entries whose structure cannot be fetched
    (`structure_not_found`) in the entries table. Empties the scratch subdirectories of `temp_dir` on the
    way in and, unless `delete_tmp` is 0, deletes them on the way out, unless an enclosing call holds it.

    Args:
        job_file (str or dict, optional): JSON job file of job key -> value, or the same already
            loaded, e.g. parse's settings. Any argument given overrides it.
        entries (str, optional): The entries table. Defaults to the job file's entries_path.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Where the log is appended. Defaults to the job file's, else none.
        threads (int, optional): Defaults to one per available core.
        temp_dir (str, optional): Scratch space for a Foldseek database download. Defaults to the job
            file's.
        delete_tmp (int, optional): 1 deletes `temp_dir` at the end; 0 keeps it. Defaults to
            DEFAULT_DELETE_TMP.

    Returns:
        None

    Raises:
        PocketMapperError: If the job file cannot be read, a path is not given, the entries table cannot
            be read or leaves a side empty, a setting is invalid, no structure for a side could be
            fetched, or a Foldseek database cannot be downloaded.
    """
    values = layer_settings(
        job_file,
        {
            "entries_path": entries,
            "verbosity": verbosity,
            "log_path": log_path,
            "threads": threads,
            "temp_dir": temp_dir,
            "delete_tmp": delete_tmp,
        },
    )
    values = resolve_step_paths(values, ("entries_path", "temp_dir"), ("log_path",))
    with run_scope("fetch_structures"), log_to_file(values["log_path"], values["verbosity"]):
        threads = resolve_threads(values["threads"])
        delete_tmp = resolve_delete_tmp(values["delete_tmp"])
        sides = open_entries(values["entries_path"], "fetch_structures")
        with temp_dir_scope(values["temp_dir"], delete_tmp):
            fetch_inputs(sides, values["entries_path"], threads, values["temp_dir"])


def fetch_inputs(sides, entries_path, threads, temp_dir):
    """
    Download the structures and Foldseek databases the given records name.

    Writes each structure to its `struct_path`. Marks those whose structure cannot be fetched in the
    entries table as `structure_not_found`.

    Args:
        sides (dict): "query" and "target" -> that side's QTRecord dicts.
        entries_path (str or None): The entries table, marked, or None for none.
        threads (int): Thread count for a Foldseek database download.
        temp_dir (str): Scratch directory; a Foldseek database download works under it.

    Returns:
        None

    Raises:
        PocketMapperError: If no structure for a side could be fetched, or a Foldseek database cannot
            be downloaded.
    """
    log_extra = {"stage": "Downloading Structures"}

    for name, records in sides.items():
        for record in records:
            if record["struct_type"] == "foldseek_db":
                fetch_missing_fsdb(record, threads, os.path.join(temp_dir, TEMP_SUBDIRS["foldseek"]))
        structure_records = [record for record in records if record["struct_type"] != "foldseek_db"]
        found = fetch_missing_structures(name, structure_records) if structure_records else {}

        failures = [
            failure(name, record, "structure_not_found")
            for record in structure_records
            if not found[record["struct_info"]]
        ]
        report_failures(entries_path, "fetch_structures", failures, log_extra, f"Missing structures for {name}(s)")
        if len(failures) == len(records):
            msg = f"Insufficient {name} structures after fetching. No valid {name} entries remain."
            logger.critical(msg, extra=log_extra)
            raise PocketMapperError(msg)


def fetch_missing_entries(sides):
    """
    Download the structure of every record whose structure is not yet on disk.

    Nothing is reported for a structure that cannot be fetched: the caller skips what is still
    missing. Foldseek databases are left alone.

    Args:
        sides (dict): Side name, e.g. "query" -> that side's QTRecord dicts.

    Returns:
        None
    """
    for name, records in sides.items():
        structure_records = [record for record in records if record["struct_type"] != "foldseek_db"]
        if structure_records:
            fetch_missing_structures(name, structure_records)


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
    found = download_structures(unique_by(records, "struct_info"))
    logger.debug(f"Structure fetcher results: {found}", extra=log_extra)

    logger.info(f"{sum(found.values())}/{len(found)} {name} required structures available", extra=log_extra)
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
