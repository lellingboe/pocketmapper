"""
Step 2: download the structures and bundled Foldseek database the entries name: some given, or the
query and target.

A Foldseek PDB database's hits are not known yet; the pockets step fetches theirs. PISA interfaces are
fetched by the pockets step.
"""

import logging
import os

import pandas as pd

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
from pocketmapper.steps.parse import parse_structure_side

logger = logging.getLogger(__name__)


def fetch_structures(
    entries=None,
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
    struct_type=None,
    out_dir=None,
    structures_tsv_path=None,
    threads=None,
    temp_dir=None,
    delete_tmp=None,
):
    """
    Download the structures and Foldseek databases some entries name.

    Only each entry's structure is resolved; its chain and residue parts are ignored. A local file is
    listed, not fetched. Adds the entries whose structure cannot be resolved (`invalid_entry`) or
    fetched (`structure_not_found`) to `failed_entries_path`. Empties `temp_dir` on the way in and,
    unless `delete_tmp` is 0, deletes it on the way out, unless an enclosing call holds it. Without a
    `results_dir`, writes no log, settings or failed entries unless that file's own path is given.

    Args:
        entries (list, optional): Entries or bare structure ids (`4Q5J`, `P12345`, `pdb`), or files
            of them. Defaults, when None or empty, to the job file's entries, else its query and
            target, where a side whose pocket method is foldseek_db is fetched as a Foldseek database
            and a Foldseek-database query entry is rejected, as `parse` rejects it.
        job_file (str or dict, optional): JSON job file of job key -> value, or the same already
            loaded, e.g. parse's settings. Any argument given overrides it.
        results_dir (str, optional): Where the log, the settings, the failed entries and `temp_dir` go
            by default. Defaults to None, for none.
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
        struct_type (str, optional): Structure type to force on every entry -- "pdb", "alphafold" or
            "foldseek_db" -- or "auto" to infer it. Defaults to DEFAULT_STRUCT_TYPE.
        out_dir (str, optional): One directory to write every structure and database to, as
            <out_dir>/<ID>.cif.gz, in place of the cache directories. Defaults to None, the cache.
        structures_tsv_path (str, optional): Where a table of every structure resolved is written:
            id, struct_type, path and whether it is on disk (ok). Defaults to None, for none.
        threads (int, optional): Defaults to one per available core.
        temp_dir (str, optional): Defaults to <results_dir>/tmp, or <cache_dir>/tmp without a
            `results_dir`.
        delete_tmp (int, optional): 1 deletes `temp_dir` at the end; 0 keeps it. Defaults to
            DEFAULT_DELETE_TMP.

    Returns:
        None

    Raises:
        PocketMapperError: If the job file cannot be read, neither entries nor query and target are
            given, a setting is invalid, no structure for a side could be fetched, or a Foldseek
            database cannot be downloaded.
    """
    values = layer_settings(
        job_file,
        {
            # The command line gives [] for none
            "entries": entries or None,
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
            "struct_type": struct_type,
            "out_dir": out_dir,
            "structures_tsv_path": structures_tsv_path,
            "threads": threads,
            "temp_dir": temp_dir,
            "delete_tmp": delete_tmp,
        },
    )
    if not values["entries"]:
        for key in ("query", "target"):
            require_setting(values, key)
    values = resolve_paths(values, "fetch_structures")
    with run_scope("fetch_structures") as outermost, log_to_file(values["log_path"], values["verbosity"]):
        if outermost:
            dump_settings(values)
        threads = resolve_threads(values["threads"])
        delete_tmp = resolve_delete_tmp(values["delete_tmp"])
        sides, sources = resolve_job_structures(values)

        roots = [root for root in (values["cache_dir"], values["results_dir"]) if root is not None]
        with temp_dir_scope(values["temp_dir"], delete_tmp, roots):
            fetch_inputs(
                sides,
                sources,
                values["failed_entries_path"],
                threads,
                values["temp_dir"],
                values["structures_tsv_path"],
            )


def resolve_job_structures(values):
    """
    Resolve the structures a fetch_structures run's settings name, reporting the entries rejected.

    Args:
        values (dict): Job key -> value, from `settings.resolve_paths`.

    Returns:
        tuple: (sides, sources). `sides` maps "entries", or "query" and "target", to that side's
            structure dicts; `sources` maps the same names to the input each was parsed from.

    Raises:
        PocketMapperError: If a structure type is unknown or an entries file cannot be read.
    """
    log_extra = {"stage": "Processing Inputs"}

    cache_dirs = values
    if values["out_dir"] is not None:
        cache_dirs = dict.fromkeys(("pdb_dir", "alphafold_dir", "fsdb_dir"), values["out_dir"])
    if values["entries"]:
        inputs = {"entries": (values["entries"], values["struct_type"])}
    else:
        # A database path is a file too, and would otherwise be read as an entries file
        inputs = {
            name: (
                values[name],
                "foldseek_db" if values[f"{name}_pocket_method"] == "foldseek_db" else values["struct_type"],
            )
            for name in ("query", "target")
        }

    sides = {}
    sources = {}
    failures = []
    for name, (entries, struct_type) in inputs.items():
        structures, rejected = parse_structure_side(entries, name, struct_type, cache_dirs, values["work_dir"])
        failures += [
            failed_entry(entry, "fetch_structures", "invalid_entry", source, detail=reason)
            for entry, reason, source in rejected
        ]
        sides[name] = structures
        sources[name] = entries if isinstance(entries, str) else " ".join(entries)

    # A database can only be searched, not searched with
    for structure in [structure for structure in sides.get("query", []) if structure["struct_type"] == "foldseek_db"]:
        reason = f"A Foldseek database cannot be a query entry: {structure['pocket_id']}"
        failures.append(
            failed_entry(
                structure["pocket_id"], "fetch_structures", "invalid_entry", sources["query"], structure, reason
            )
        )
        sides["query"].remove(structure)
    report_failures(values["failed_entries_path"], failures, log_extra)
    return sides, sources


def fetch_inputs(sides, sources, failed_entries_path, threads, temp_dir, structures_tsv_path=None):
    """
    Download the structures and Foldseek databases the given records or structures name.

    Writes each structure to its `struct_path`. Those whose structure cannot be fetched are added to
    `failed_entries_path` as `structure_not_found`.

    Args:
        sides (dict): Side name, e.g. "query" -> that side's QTRecord or structure dicts; each needs
            `pocket_id`, `struct_info`, `struct_type` and `struct_path`.
        sources (dict): Side name -> the input the side was parsed from, for the failure entries.
        failed_entries_path (str or None): The failed-entries file, appended to, or None for none.
        threads (int): Thread count for a Foldseek database download.
        temp_dir (str): Scratch directory; a Foldseek database download works under it.
        structures_tsv_path (str, optional): Where a table of every structure is written: id,
            struct_type, path and ok, whether it is on disk. Defaults to None, for none.

    Returns:
        None

    Raises:
        PocketMapperError: If no structure for a side could be fetched, a Foldseek database cannot
            be downloaded, or the table's directory cannot be created.
    """
    log_extra = {"stage": "Downloading Structures"}

    # The table is written however the fetching ends, so it shows what did arrive
    try:
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
                msg = f"Insufficient {name} structures after fetching. No valid {name} entries remain."
                logger.critical(msg, extra=log_extra)
                raise PocketMapperError(msg)
    finally:
        if structures_tsv_path is not None:
            write_structures_table(sides, structures_tsv_path)


def write_structures_table(sides, path):
    """
    Write a table of every structure: id, struct_type, path and ok, whether it is on disk.

    Args:
        sides (dict): Side name -> structure or record dicts. A structure several entries name is
            listed once.
        path (str): The table to write, tab-separated.

    Returns:
        None

    Raises:
        PocketMapperError: If the table's directory cannot be created.
    """
    rows = {}
    for record in [record for records in sides.values() for record in records]:
        rows.setdefault(
            record["struct_path"],
            {
                "id": record["struct_info"],
                "struct_type": record["struct_type"],
                "path": record["struct_path"],
                "ok": int(os.path.exists(record["struct_path"])),
            },
        )
    make_dir(os.path.dirname(path), {"stage": "Downloading Structures"})
    pd.DataFrame(list(rows.values()), columns=["id", "struct_type", "path", "ok"]).to_csv(path, sep="\t", index=False)
    logger.info(f"Structures listed in {path}", extra={"stage": "Downloading Structures"})


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
