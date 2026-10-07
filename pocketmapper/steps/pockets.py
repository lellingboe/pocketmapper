"""
Step 4: build the pocket of every query and target entry into one pockets file.

Every entry given is listed, `null` where its pocket could not be built. Against a PDB-named Foldseek
database, each hit's PISA interfaces are expanded into pockets of their own here, keyed in the
pockets file by hit name.
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
from pocketmapper.records import require_file
from pocketmapper.records import split_missing_structures
from pocketmapper.settings import dump_settings
from pocketmapper.settings import input_path
from pocketmapper.settings import layer_settings
from pocketmapper.settings import require_setting
from pocketmapper.settings import resolve_paths
from pocketmapper.settings import resolve_pisa_source
from pocketmapper.steps.fetch_structures import fetch_missing_structures
from pocketmapper.steps.parse import parse_job_entries

logger = logging.getLogger(__name__)


def pockets(
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
    alignment=None,
    pockets_path=None,
    pisa_source=None,
):
    """
    Build the pocket of every query and target entry into one pockets file.

    Parses the entries the job file names. Overwrites `pockets_path`, and adds the entries whose
    structure is missing or whose pocket could not be built to `failed_entries_path`. Downloads any
    PISA data not already cached and, against a PDB-named Foldseek database, the structures of its
    hits.

    Args:
        job_file (str or dict, optional): JSON job file of job key -> value, or the same already
            loaded, e.g. parse's settings. Any argument given overrides it. Must set query and target.
        results_dir (str, optional): Required here or in `job_file`.
        work_dir (str, optional): Directory that entries and relative paths resolve against.
            Defaults to the working directory.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Defaults to <results_dir>/info.log.
        job_settings_path (str, optional): Where these settings are written, as a job file for later
            steps. Defaults to <results_dir>/pockets_settings.json. Not written when run inside search.
        failed_entries_path (str, optional): Defaults to <results_dir>/failed_entries.json.
        cache_dir (str, optional): Defaults to DEFAULT_CACHE_DIR.
        pdb_dir (str, optional): Defaults to <cache_dir>/pdb_structures.
        alphafold_dir (str, optional): Defaults to <cache_dir>/alphafold_structures.
        pocket_dir (str, optional): Defaults to <cache_dir>/pockets.
        foldseek_preprocessed_structure_dir (str, optional): Defaults to
            <cache_dir>/foldseek_preprocessed_structures.
        fsdb_dir (str, optional): Defaults to <cache_dir>/fsdb.
        alignment (str, optional): The alignment table, read only for a Foldseek-database target.
            Defaults to the job file's alignment_path, else <results_dir>/alignment.tsv.
        pockets_path (str, optional): Defaults to <results_dir>/pockets.json.
        pisa_source (str, optional): "ftp" or "api". Defaults to DEFAULT_PISA_SOURCE.

    Returns:
        None

    Raises:
        PocketMapperError: If the job file cannot be read, query, target or results_dir is not given,
            the entries are rejected as `parse` rejects them, `pisa_source` is invalid, or a
            Foldseek-database target has no alignment.
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
            "pockets_path": pockets_path,
            "pisa_source": pisa_source,
        },
    )
    for key in ("query", "target", "results_dir"):
        require_setting(values, key)
    values = resolve_paths(values, "pockets")
    with run_scope("pockets") as outermost, log_to_file(values["log_path"], values["verbosity"]):
        if outermost:
            dump_settings(values)
        pisa_source = resolve_pisa_source(values["pisa_source"])
        sides = parse_job_entries(values, "pockets")
        build_pockets(
            sides,
            {"query": values["query"], "target": values["target"]},
            input_path(values, alignment, "alignment_path"),
            values["pockets_path"],
            values["failed_entries_path"],
            values,
            pisa_source,
        )


def build_pockets(sides, sources, alignment, pockets_path, failed_entries_path, cache_dirs, pisa_source):
    """
    Build a Pocket for every record of both sides, and write them to `pockets_path`.

    Adds the records whose structure is missing (`structure_not_found`) or whose pocket could not be
    built (`pocket_not_built`) to `failed_entries_path`. Foldseek-database records have no pocket.
    Against a PDB-named Foldseek database, also builds a pisa pocket per PISA interface of each hit
    (see `expand_fsdb_pdb_targets`). The pockets file lists every record and hit pocket, `null`
    where it was not built, and maps each `preprocess_name` and hit name to its pocket_ids.
    Downloads any PISA data not already cached under `pocket_dir/pisa/`, and overwrites the
    per-method pocket files under `pocket_dir`.

    Args:
        sides (dict): "query" and "target" -> that side's QTRecord dicts. Built query first, then
            target, then the hits, which decides which of several identical pisa records names the
            shared pocket.
        sources (dict): "query" and "target" -> the input the side was parsed from, for the failure
            entries.
        alignment (str): The alignment table. Read only when the target is a Foldseek database.
        pockets_path (str): Where the pockets are written, overwriting any there.
        failed_entries_path (str): The failed-entries file, appended to.
        cache_dirs (dict): Holds "pdb_dir", "alphafold_dir", "fsdb_dir" and "pocket_dir".
        pisa_source (str): Where PISA interfaces are fetched from: "ftp" or "api".

    Returns:
        None

    Raises:
        PocketMapperError: If a Foldseek-database target has no alignment, or a PDB-named database's
            hits have interfaces but none of their structures could be fetched.
    """
    log_extra = {"stage": "Getting Pockets"}

    # (record, the input it came from), in build order
    given = [
        (record, sources[name])
        for name, records in sides.items()
        for record in records
        if record["struct_type"] != "foldseek_db"
    ]
    chains = preproc_to_ids([record for record, _ in given])
    buildable = []
    for name, records in sides.items():
        kept, missing = split_missing_structures(
            [record for record in records if record["struct_type"] != "foldseek_db"],
            "pockets",
            sources[name],
            log_extra,
        )
        buildable += [(record, sources[name]) for record in kept]
        append_failed_entries(failed_entries_path, missing)

    if fsdb_record(sides["target"]) is not None:
        require_file(alignment, "alignment")
        hit_names = pd.read_csv(alignment, sep="\t", usecols=["target"], dtype=str)["target"].unique().tolist()
        if fsdb_pocket_mode(hit_names) == "pisa":
            hits, hit_chains, failures = expand_fsdb_pdb_targets(hit_names, cache_dirs, pisa_source, alignment)
            append_failed_entries(failed_entries_path, failures)
            buildable += [(record, alignment) for record in hits]
            chains |= hit_chains

    built = PocketFetcher().fetch_pockets(
        [record for record, _ in buildable],
        cache_dirs["pocket_dir"],
        builder_options={"pisa": {"pisa_source": pisa_source}},
    )

    # Checked on the merged pockets: a later method can overwrite an earlier one's pocket with None
    failures = [
        failed_entry(record["pocket_id"], "pockets", "pocket_not_built", source, record)
        for record, source in buildable
        if built.get(record["pocket_id"]) is None
    ]
    if failures:
        logger.warning(
            f"No pocket for {', '.join(dict.fromkeys(entry['pocket_id'] for entry in failures))}; skipping them",
            extra=log_extra,
        )
    append_failed_entries(failed_entries_path, failures)

    # Every pocket given or expanded is listed, None where it was not built
    pocket_ids = [record["pocket_id"] for record, _ in given] + [pid for ids in chains.values() for pid in ids]
    pockets = {pocket_id: built.get(pocket_id) for pocket_id in pocket_ids}
    make_dir(os.path.dirname(pockets_path), log_extra)
    write_pockets_file(pockets, chains, pockets_path)
    logger.info(f"Pockets written to {pockets_path}", extra=log_extra)


def expand_fsdb_pdb_targets(hit_names, cache_dirs, pisa_source, source):
    """
    Turn the hits of a PDB Foldseek-database search into pisa records, one per PISA interface.

    Downloads the PISA data and the structures those records need. Each record's `preprocess_name` is
    the hit name rather than the one `QTProcessor` derives.

    Args:
        hit_names (list): The database entry names the search hit.
        cache_dirs (dict): Holds "pdb_dir", "alphafold_dir", "fsdb_dir" and "pocket_dir".
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
