"""
Step 2: download the structures and bundled Foldseek database the entries table names.

A Foldseek PDB database's hits are not known yet; the pockets step fetches theirs. PISA interfaces are
fetched by the pockets step. With `offline` 1 nothing is downloaded: the cache is only checked.
"""

import logging
import os

from pocketmapper.lib.constants import TEMP_SUBDIRS
from pocketmapper.lib.entries import failure, open_entries, report_failures
from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.logs import log_to_file
from pocketmapper.lib.scopes import run_scope, temp_dir_scope
from pocketmapper.lib.settings import (
    check_offline_refresh,
    layer_settings,
    resolve_delete_tmp,
    resolve_offline,
    resolve_refresh_structures,
    resolve_step_paths,
    resolve_threads,
)
from pocketmapper.lib.structure_manager import StructureManager, fetch_missing_fsdb, fetch_missing_structures

logger = logging.getLogger(__name__)


def fetch_structures(
    job_file=None,
    entries=None,
    verbosity=None,
    log_path=None,
    threads=None,
    temp_dir=None,
    delete_tmp=None,
    offline=None,
    refresh_structures=None,
    pdb_dir=None,
    alphafold_dir=None,
    fsdb_dir=None,
    local_root=None,
):
    """
    Download the structures and Foldseek databases the entries table names.

    A local file is checked, not fetched, and so is everything with `offline` 1. Marks the entries whose
    structure cannot be fetched (`structure_not_found`) in the entries table. Empties the scratch
    subdirectories of `temp_dir` on the way in and, unless `delete_tmp` is 0, deletes them on the way out,
    unless an enclosing call holds it.

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
        offline (int, optional): 1 downloads nothing: an uncached structure is marked, an uncached
            Foldseek database is an error. Defaults to DEFAULT_OFFLINE.
        refresh_structures (int, optional): 1 re-downloads PDB and AlphaFold structures even when they
            are cached. Defaults to DEFAULT_REFRESH_STRUCTURES.
        pdb_dir (str, optional): Cache of PDB structures. Defaults to the job file's.
        alphafold_dir (str, optional): Cache of AlphaFold structures. Defaults to the job file's.
        fsdb_dir (str, optional): Cache of Foldseek databases. Defaults to the job file's.
        local_root (str, optional): Directory that relative local structure files and Foldseek databases
            resolve against. Defaults to the job file's, else the current directory.

    Returns:
        None

    Raises:
        PocketMapperError: If the job file cannot be read, a path is not given, the entries table cannot
            be read or leaves a side empty, a setting is invalid, `offline` and `refresh_structures` are
            both 1, no structure for a side could be fetched, or a Foldseek database cannot be
            downloaded (or, offline, is not cached).
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
            "offline": offline,
            "refresh_structures": refresh_structures,
            "pdb_dir": pdb_dir,
            "alphafold_dir": alphafold_dir,
            "fsdb_dir": fsdb_dir,
            "local_root": local_root,
        },
    )
    values = resolve_step_paths(
        values, ("entries_path", "temp_dir", "pdb_dir", "alphafold_dir", "fsdb_dir"), ("log_path", "local_root")
    )
    with run_scope("fetch_structures"), log_to_file(values["log_path"], values["verbosity"]):
        threads = resolve_threads(values["threads"])
        delete_tmp = resolve_delete_tmp(values["delete_tmp"])
        offline = resolve_offline(values["offline"])
        refresh_structures = resolve_refresh_structures(values["refresh_structures"])
        check_offline_refresh(offline, refresh_structures)
        sides = open_entries(values["entries_path"], "fetch_structures")
        with temp_dir_scope(values["temp_dir"], delete_tmp):
            manager = StructureManager(
                values["pdb_dir"],
                values["alphafold_dir"],
                values["fsdb_dir"],
                values["local_root"],
                download=not offline,
                refresh=bool(refresh_structures),
            )
            fetch_inputs(sides, manager, values["entries_path"], threads, values["temp_dir"])


def fetch_inputs(sides, manager, entries_path, threads, temp_dir):
    """
    Download the structures and Foldseek databases the given records name.

    Writes each structure and database to the path `manager` gives. Marks those whose structure cannot be
    fetched in the entries table as `structure_not_found`.

    Args:
        sides (dict): "query" and "target" -> that side's QTRecord dicts.
        manager (StructureManager): Locates and fetches the structures and databases, under its download
            and refresh policy.
        entries_path (str or None): The entries table, marked, or None for none.
        threads (int): Thread count for a Foldseek database download.
        temp_dir (str): Scratch directory; a Foldseek database download works under it.

    Returns:
        None

    Raises:
        PocketMapperError: If no structure for a side could be fetched, a Foldseek database cannot be
            downloaded or, with `manager` not downloading, is not cached, or `manager` lacks a directory a
            record needs.
    """
    log_extra = {"stage": "Downloading Structures"}

    for name, records in sides.items():
        for record in records:
            if record["struct_type"] == "foldseek_db":
                fetch_missing_fsdb(record, manager, threads, os.path.join(temp_dir, TEMP_SUBDIRS["foldseek"]))
        structure_records = [record for record in records if record["struct_type"] != "foldseek_db"]
        found = fetch_missing_structures(name, structure_records, manager) if structure_records else {}

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
