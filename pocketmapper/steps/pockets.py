"""
Step 4: build the pocket of every record into one pockets file.

A record whose pocket cannot be built is dropped from its records file, so every record left has a
pocket for the comparison to read. Against a PDB-named Foldseek database, each hit's PISA interfaces
are expanded into pockets of their own here, keyed in the pockets file by hit name.
"""

import logging
import os
import re
from dataclasses import asdict

import pandas as pd

from pocketmapper.exceptions import PocketMapperError
from pocketmapper.lib import fsdb_pocket_mode
from pocketmapper.lib import log_to_file
from pocketmapper.lib import make_dir
from pocketmapper.lib import parse_foldseek_pdb_entry_name
from pocketmapper.lib import run_scope
from pocketmapper.pockets.pisa import PisaParser
from pocketmapper.pockets.pisa import download_pisa_interfaces
from pocketmapper.pockets.pocket_fetcher import PocketFetcher
from pocketmapper.pockets.pocket_fetcher import write_pockets_file
from pocketmapper.qt_processor import QTProcessor
from pocketmapper.records import append_failed_entries
from pocketmapper.records import failed_entry
from pocketmapper.records import fsdb_record
from pocketmapper.records import preproc_to_ids
from pocketmapper.records import read_cache_manifest
from pocketmapper.records import read_records
from pocketmapper.records import require_file
from pocketmapper.records import write_records
from pocketmapper.settings import dump_settings
from pocketmapper.settings import input_path
from pocketmapper.settings import layer_settings
from pocketmapper.settings import require_setting
from pocketmapper.settings import resolve_paths
from pocketmapper.settings import resolve_pisa_source
from pocketmapper.steps.fetch_structures import fetch_missing_structures

logger = logging.getLogger(__name__)


def pockets(
    records=None,
    job_file=None,
    results_dir=None,
    work_dir=None,
    verbosity=None,
    log_path=None,
    job_settings_path=None,
    failed_entries_path=None,
    alignment=None,
    pockets_path=None,
    pisa_source=None,
):
    """
    Build the pocket of every record in some records files.

    Reads the cache directories from `results_dir`'s cache manifest. Overwrites `pockets_path`
    with the pockets of every file named, and rewrites each file in place without the records whose
    pocket could not be built, adding those to `failed_entries_path`. Downloads any PISA data not
    already cached and, against a PDB-named Foldseek database, the structures of its hits.

    Args:
        records (list, optional): The records files. The pockets file holds only the pockets of the
            files named. Defaults, when None or empty, to the query and target records files: the
            job file's query_records_path and target_records_path, else the standard files under
            <results_dir>.
        job_file (str or dict, optional): JSON job file of job key -> value, or the same
            already loaded. Any argument given overrides it.
        results_dir (str, optional): The results directory `parse` wrote to. Required here or in
            `job_file`.
        work_dir (str, optional): Directory that entries and relative paths resolve against.
            Defaults to the working directory.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Defaults to <results_dir>/info.log.
        job_settings_path (str, optional): Where these settings are written, as a job file for later
            steps. Defaults to <results_dir>/pockets_settings.json. Not written when run inside search.
        failed_entries_path (str, optional): Defaults to <results_dir>/failed_entries.json.
        alignment (str, optional): The alignment table, read only for a Foldseek-database target.
            Defaults to the job file's alignment_path, else <results_dir>/alignment.tsv.
        pockets_path (str, optional): Defaults to <results_dir>/pockets.json.
        pisa_source (str, optional): "ftp" or "api". Defaults to DEFAULT_PISA_SOURCE.

    Returns:
        None

    Raises:
        PocketMapperError: If the job file cannot be read, no results_dir is given, the manifest or a
            records file is missing, `pisa_source` is invalid, or a Foldseek-database target has no
            alignment.
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
            "pockets_path": pockets_path,
            "pisa_source": pisa_source,
        },
    )
    require_setting(values, "results_dir")
    values = resolve_paths(values, "pockets")
    # The command line gives [] for none
    if not records:
        records = [values["query_records_path"], values["target_records_path"]]
    with run_scope("pockets") as outermost, log_to_file(values["log_path"], values["verbosity"]):
        if outermost:
            dump_settings(values)
        cache_dirs = read_cache_manifest(values["results_dir"])
        pisa_source = resolve_pisa_source(values["pisa_source"])
        build_pockets(
            records,
            input_path(values, alignment, "alignment_path"),
            values["pockets_path"],
            values["failed_entries_path"],
            cache_dirs,
            pisa_source,
        )


def build_pockets(records_paths, alignment, pockets_path, failed_entries_path, cache_dirs, pisa_source):
    """
    Build a Pocket for every record across several records files, and write them to `pockets_path`.

    Rewrites each records file in place without the records whose pocket could not be built, and
    adds those to `failed_entries_path` as `pocket_not_built`. Foldseek-database records have no
    pocket and are kept. Against a PDB-named Foldseek database, also builds a pisa pocket per PISA
    interface of each hit (see `expand_fsdb_pdb_targets`). The pockets file lists every record and
    hit pocket, `null` where it was not built, and maps each `preprocess_name` and hit name to its
    pocket_ids. Downloads any PISA data not already cached under `pocket_dir/pisa/`, and overwrites
    the per-method pocket files under `pocket_dir`.

    Args:
        records_paths (list): The records files, e.g. query then target. Records are built in this
            order, then the hits, which decides which of several identical pisa records names the
            shared pocket.
        alignment (str): The alignment table. Read only when a record is a Foldseek database.
        pockets_path (str): Where the pockets are written, overwriting any there.
        failed_entries_path (str): The failed-entries file, appended to.
        cache_dirs (dict): Each of `records.CACHE_MANIFEST_KEYS` -> its absolute directory.
        pisa_source (str): Where PISA interfaces are fetched from: "ftp" or "api".

    Returns:
        None

    Raises:
        PocketMapperError: If a records file cannot be read, a Foldseek-database target has no
            alignment, or a PDB-named database's hits have interfaces but none of their structures
            could be fetched.
    """
    log_extra = {"stage": "Getting Pockets"}

    files = {path: read_records(path) for path in records_paths}
    buildable = [record for records in files.values() for record in records if record["struct_type"] != "foldseek_db"]
    chains = preproc_to_ids(buildable)

    hits = []
    if fsdb_record([record for records in files.values() for record in records]) is not None:
        require_file(alignment, "alignment")
        hit_names = pd.read_csv(alignment, sep="\t", usecols=["target"], dtype=str)["target"].unique().tolist()
        if fsdb_pocket_mode(hit_names) == "pisa":
            hits, hit_chains, failures = expand_fsdb_pdb_targets(hit_names, cache_dirs, pisa_source, alignment)
            append_failed_entries(failed_entries_path, failures)
            chains |= hit_chains

    built = PocketFetcher().fetch_pockets(
        buildable + hits,
        cache_dirs["pocket_dir"],
        builder_options={"pisa": {"pisa_source": pisa_source}},
    )

    # Checked on the merged pockets: a later method can overwrite an earlier one's pocket with None
    failures = []
    for path, records in files.items():
        kept = []
        for record in records:
            if record["struct_type"] == "foldseek_db" or built.get(record["pocket_id"]) is not None:
                kept.append(record)
            else:
                failures.append(failed_entry(record["pocket_id"], "pockets", "pocket_not_built", path, record))
        files[path] = kept
    failures += [
        failed_entry(record["pocket_id"], "pockets", "pocket_not_built", alignment, record)
        for record in hits
        if built.get(record["pocket_id"]) is None
    ]
    if failures:
        logger.warning(
            f"No pocket for {', '.join(dict.fromkeys(entry['pocket_id'] for entry in failures))}; skipping them",
            extra=log_extra,
        )
    append_failed_entries(failed_entries_path, failures)

    # Every pocket given or expanded is listed, None where it was not built
    pocket_ids = [record["pocket_id"] for record in buildable] + [pid for ids in chains.values() for pid in ids]
    pockets = {pocket_id: built.get(pocket_id) for pocket_id in pocket_ids}
    make_dir(os.path.dirname(pockets_path), log_extra)
    write_pockets_file(pockets, chains, pockets_path)
    logger.info(f"Pockets written to {pockets_path}", extra=log_extra)
    for path, records in files.items():
        write_records(records, path)


def expand_fsdb_pdb_targets(hit_names, cache_dirs, pisa_source, source):
    """
    Turn the hits of a PDB Foldseek-database search into pisa records, one per PISA interface.

    Downloads the PISA data and the structures those records need. Each record's `preprocess_name` is
    the hit name rather than the one `QTProcessor` derives.

    Args:
        hit_names (list): The database entry names the search hit.
        cache_dirs (dict): Each of `records.CACHE_MANIFEST_KEYS` -> its absolute directory.
        pisa_source (str): Where PISA interfaces are fetched from: "ftp" or "api".
        source (str): The alignment table, for the failure entries.

    Returns:
        tuple: (records, chains, failure entries). `records` holds those whose structure is
            available. `chains` maps every hit name to the pocket_ids of its records, including those
            whose structure could not be fetched (each a `structure_not_found` failure); it is empty
            for a hit that is not PDB-style or has no interface the pisa pocket method can spell.

    Raises:
        PocketMapperError: If interfaces were found but none of their structures could be fetched.
    """
    log_extra = {"stage": "Expanding Foldseek DB Targets"}

    chains = {hit_name: [] for hit_name in hit_names}
    hits = {}  # foldseek entry name -> (pdb_id, chain_id)
    for hit_name in hit_names:
        resolved = parse_foldseek_pdb_entry_name(hit_name)
        if resolved is not None:
            hits[hit_name] = resolved

    pdb_list = sorted({pdb_id for pdb_id, _ in hits.values()})
    logger.info(
        f"Retrieving PISA interfaces for {len(pdb_list)} PDB entries behind {len(hits)} Foldseek hits",
        extra=log_extra,
    )
    interface_dir = download_pisa_interfaces(pdb_list, cache_dirs["pocket_dir"], pisa_source)

    # Building one record per interface the hit chain takes part in
    parser = PisaParser()
    qtprocessor = QTProcessor(
        pdb_dir=cache_dirs["pdb_dir"],
        alphafold_dir=cache_dirs["alphafold_dir"],
        fsdb_dir=cache_dirs["fsdb_dir"],
    )
    # PISA stores chain pairs the input grammar cannot spell (multi-character chain ids); those are
    # counted and skipped here rather than each rejected with a warning.
    pisa_pattern = qtprocessor.pocket_methods["pisa"][0]
    unspellable = 0
    records = []
    for hit_name, (pdb_id, chain_id) in hits.items():
        for partner in parser.get_interface_partners(pdb_id, chain_id, interface_dir):
            if not re.match(pisa_pattern, f"{chain_id}_{partner}"):
                unspellable += 1
                continue
            record, _ = qtprocessor.parse_individual_qt(f"{pdb_id}:{chain_id}_{partner}", pocket_method="pisa")
            if record is None:
                continue  # the reason is logged
            # The alignment is keyed by the Foldseek entry name
            record.preprocess_name = hit_name
            records.append(asdict(record))
            chains[hit_name].append(record.pocket_id)
    if unspellable:
        logger.info(
            f"Skipped {unspellable} PISA interfaces whose chain ids the pisa pocket method cannot express",
            extra=log_extra,
        )
    if not records:
        logger.warning("No PISA interfaces found for any Foldseek hit", extra=log_extra)
        return [], chains, []

    # Fetched last, so only entries with an interface are downloaded. Failed fetches are dropped.
    found = fetch_missing_structures("foldseek hit", records)
    if not any(found.values()):
        logger.critical("Insufficient foldseek hit structures after fetching", extra=log_extra)
        raise PocketMapperError(
            "Insufficient foldseek hit structures after fetching. No valid foldseek hit entries remain."
        )
    failures = [
        failed_entry(record["pocket_id"], "pockets", "structure_not_found", source, record)
        for record in records
        if not found[record["struct_info"]]
    ]
    records = [record for record in records if found[record["struct_info"]]]

    logger.info(
        f"Added {len(records)} PISA target pockets from {len({r['preprocess_name'] for r in records})} Foldseek hits",
        extra=log_extra,
    )
    return records, chains, failures
