"""
Step 2: download the structures and bundled Foldseek database the query and target entries need.

A Foldseek PDB database's hits are not known yet; the align step fetches theirs. PISA interfaces are
fetched by the pockets step.
"""

import logging
import os

from pocketmapper.downloads.structure_downloader import StructureDownloader
from pocketmapper.entries import failed_entry
from pocketmapper.entries import report_failures
from pocketmapper.entries import unique_by
from pocketmapper.exceptions import PocketMapperError
from pocketmapper.foldseek import run_foldseek
from pocketmapper.lib import log_to_file
from pocketmapper.lib import make_dir
from pocketmapper.lib import run_scope
from pocketmapper.lib import temp_dir_scope
from pocketmapper.settings import dump_settings
from pocketmapper.settings import layer_settings
from pocketmapper.settings import require_foldseek
from pocketmapper.settings import require_setting
from pocketmapper.settings import resolve_delete_tmp
from pocketmapper.settings import resolve_paths
from pocketmapper.settings import resolve_threads
from pocketmapper.steps.parse import parse_job_entries

logger = logging.getLogger(__name__)


def fetch_structures(
    job_file=None,
    results_dir=None,
    work_dir=None,
    verbosity=None,
    log_path=None,
    job_settings_path=None,
    failed_entries_path=None,
    cache_dir=None,
    pdb_dir=None,
    alphafold_dir=None,
    pocket_dir=None,
    foldseek_preprocessed_structure_dir=None,
    fsdb_dir=None,
    threads=None,
    temp_dir=None,
    delete_tmp=None,
):
    """
    Download the structures and Foldseek database the query and target entries need.

    Parses the entries the job file names. Adds the entries whose structure cannot be fetched to
    `failed_entries_path`. Empties `temp_dir` on the way in and, unless `delete_tmp` is 0, deletes it
    on the way out, unless an enclosing call holds it.

    Args:
        job_file (str or dict, optional): JSON job file of job key -> value, or the same already
            loaded, e.g. parse's settings. Any argument given overrides it. Must set query and target.
        results_dir (str, optional): Required here or in `job_file`.
        work_dir (str, optional): Directory that entries and relative paths resolve against.
            Defaults to the working directory.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Defaults to <results_dir>/info.log.
        job_settings_path (str, optional): Where these settings are written, as a job file for later
            steps. Defaults to <results_dir>/fetch_structures_settings.json. Not written when run
            inside search.
        failed_entries_path (str, optional): Defaults to <results_dir>/failed_entries.json.
        cache_dir (str, optional): Defaults to DEFAULT_CACHE_DIR.
        pdb_dir (str, optional): Defaults to <cache_dir>/pdb_structures.
        alphafold_dir (str, optional): Defaults to <cache_dir>/alphafold_structures.
        pocket_dir (str, optional): Defaults to <cache_dir>/pockets.
        foldseek_preprocessed_structure_dir (str, optional): Defaults to
            <cache_dir>/foldseek_preprocessed_structures.
        fsdb_dir (str, optional): Defaults to <cache_dir>/fsdb.
        threads (int, optional): Defaults to one per available core.
        temp_dir (str, optional): Defaults to <results_dir>/tmp.
        delete_tmp (int, optional): 1 deletes `temp_dir` at the end; 0 keeps it. Defaults to
            DEFAULT_DELETE_TMP.

    Returns:
        None

    Raises:
        PocketMapperError: If the job file cannot be read, query, target or results_dir is not given,
            the entries are rejected as `parse` rejects them, a setting is invalid, no structure for a
            side could be fetched, or a Foldseek database cannot be downloaded.
    """
    values = layer_settings(
        job_file,
        {
            "results_dir": results_dir,
            "work_dir": work_dir,
            "verbosity": verbosity,
            "log_path": log_path,
            "job_settings_path": job_settings_path,
            "failed_entries_path": failed_entries_path,
            "cache_dir": cache_dir,
            "pdb_dir": pdb_dir,
            "alphafold_dir": alphafold_dir,
            "pocket_dir": pocket_dir,
            "foldseek_preprocessed_structure_dir": foldseek_preprocessed_structure_dir,
            "fsdb_dir": fsdb_dir,
            "threads": threads,
            "temp_dir": temp_dir,
            "delete_tmp": delete_tmp,
        },
    )
    for key in ("query", "target", "results_dir"):
        require_setting(values, key)
    values = resolve_paths(values, "fetch_structures")
    with run_scope("fetch_structures") as outermost, log_to_file(values["log_path"], values["verbosity"]):
        if outermost:
            dump_settings(values)
        threads = resolve_threads(values["threads"])
        delete_tmp = resolve_delete_tmp(values["delete_tmp"])
        sides = parse_job_entries(values, "fetch_structures")

        roots = [values["cache_dir"], values["results_dir"]]
        with temp_dir_scope(values["temp_dir"], delete_tmp, roots):
            fetch_inputs(
                sides,
                {"query": values["query"], "target": values["target"]},
                values["failed_entries_path"],
                threads,
                values["temp_dir"],
            )


def fetch_inputs(sides, sources, failed_entries_path, threads, temp_dir):
    """
    Download the structures and Foldseek database both sides' records need.

    Writes each structure to its record's `struct_path`. Records whose structure cannot be fetched
    are added to `failed_entries_path` as `structure_not_found`.

    Args:
        sides (dict): "query" and "target" -> that side's QTRecord dicts.
        sources (dict): "query" and "target" -> the input the side was parsed from, for the failure
            entries.
        failed_entries_path (str): The failed-entries file, appended to.
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
                fetch_missing_fsdb(record, threads, os.path.join(temp_dir, "foldseek_tmp"))
        structure_records = [record for record in records if record["struct_type"] != "foldseek_db"]
        found = fetch_missing_structures(name, structure_records) if structure_records else {}

        failures = [
            failed_entry(record["pocket_id"], "fetch_structures", "structure_not_found", sources[name], record)
            for record in structure_records
            if not found[record["struct_info"]]
        ]
        report_failures(failed_entries_path, failures, log_extra, f"Missing structures for {name}(s)")
        if len(failures) == len(records):
            logger.critical(f"Insufficient {name} structures after fetching", extra=log_extra)
            raise PocketMapperError(f"Insufficient {name} structures after fetching. No valid {name} entries remain.")


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
    found = StructureDownloader().download_missing_structures(unique_by(records, "struct_info"))
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
