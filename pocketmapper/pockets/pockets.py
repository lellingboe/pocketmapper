"""
Step 4: build the pocket of every query and target entry into one pockets file.

Every entry no earlier step dropped is listed, `null` where its pocket could not be built. Against a
PDB-named Foldseek database, each hit's PISA interfaces are expanded into pockets of their own here,
keyed in the pockets file by hit name.
"""

import logging
import os
import re
from dataclasses import asdict

import pandas as pd

from pocketmapper.fetch.fetch import fetch_missing_entries, fetch_missing_structures
from pocketmapper.lib.entries import (
    failure,
    fsdb_record,
    open_entries,
    preproc_to_ids,
    report_failures,
    require_file,
    split_missing_structures,
)
from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.files import make_dir
from pocketmapper.lib.foldseek import fsdb_pocket_mode, parse_foldseek_pdb_entry_name
from pocketmapper.lib.logs import log_to_file
from pocketmapper.lib.residues import split_chain_info
from pocketmapper.lib.scopes import run_scope
from pocketmapper.lib.settings import (
    layer_settings,
    resolve_fetch_missing,
    resolve_offline,
    resolve_pisa_bsa_threshold,
    resolve_pisa_source,
    resolve_step_paths,
)
from pocketmapper.lib.structure_manager import StructureManager
from pocketmapper.parse.qt_processor import QTProcessor
from pocketmapper.pockets.pisa import download_pisa_interfaces, interface_partners, load_interfaces
from pocketmapper.pockets.pocket_fetcher import fetch_pockets, write_pockets_file

logger = logging.getLogger(__name__)


def pockets(
    job_file=None,
    entries=None,
    alignment=None,
    verbosity=None,
    log_path=None,
    pdb_dir=None,
    alphafold_dir=None,
    pocket_dir=None,
    fsdb_dir=None,
    pockets_path=None,
    pockets_tsv_path=None,
    pisa_source=None,
    pisa_bsa_threshold=None,
    fetch_missing=None,
    offline=None,
):
    """
    Build the pocket of every query and target entry into one pockets file.

    Overwrites `pockets_path`, and marks the entries whose structure is missing or whose pocket could
    not be built in the entries table. Unless `offline` is 1, downloads any PISA data not already cached
    and, against a PDB-named Foldseek database target, the structures of its hits.

    Args:
        job_file (str or dict, optional): JSON job file of job key -> value, or the same already
            loaded, e.g. parse's settings. Any argument given overrides it.
        entries (str, optional): The entries table. Defaults to the job file's entries_path.
        alignment (str, optional): The alignment table, read only for a Foldseek-database target.
            Defaults to the job file's alignment_path.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Where the log is appended. Defaults to the job file's, else none.
        pdb_dir (str, optional): Where a Foldseek hit's PDB structure is fetched to. Defaults to the
            job file's.
        alphafold_dir (str, optional): Defaults to the job file's.
        pocket_dir (str, optional): Cache of parsed pockets and PISA data. Defaults to the job file's.
        fsdb_dir (str, optional): Defaults to the job file's.
        pockets_path (str, optional): Defaults to the job file's.
        pockets_tsv_path (str, optional): Where a table of every pocket is written: pocket_id, chain,
            res_auth_ids and method. Defaults to None, for none.
        pisa_source (str, optional): "ftp" or "api". Defaults to DEFAULT_PISA_SOURCE.
        pisa_bsa_threshold (float, optional): Build pisa pockets from buried surface: a residue joins
            when the partner chain buries more than this percent of its accessible surface, at least 0
            and below 100. Defaults to the job file's, else None, which builds them from interface bonds.
        fetch_missing (int, optional): 1 downloads an entry structure missing from the cache, 0 skips
            the entry as `structure_not_found`. Defaults to DEFAULT_FETCH_MISSING.
        offline (int, optional): 1 downloads nothing, whatever `fetch_missing` says: an uncached
            structure is skipped as with `fetch_missing` 0, and a pisa pocket whose PISA data is not
            cached is not built. Defaults to DEFAULT_OFFLINE.

    Returns:
        None

    Raises:
        PocketMapperError: If the job file cannot be read, a path is not given, the entries table cannot
            be read or leaves a side empty, a setting is invalid, or a Foldseek-database target has no
            alignment.
    """
    values = layer_settings(
        job_file,
        {
            "entries_path": entries,
            "alignment_path": alignment,
            "verbosity": verbosity,
            "log_path": log_path,
            "pdb_dir": pdb_dir,
            "alphafold_dir": alphafold_dir,
            "pocket_dir": pocket_dir,
            "fsdb_dir": fsdb_dir,
            "pockets_path": pockets_path,
            "pockets_tsv_path": pockets_tsv_path,
            "pisa_source": pisa_source,
            "pisa_bsa_threshold": pisa_bsa_threshold,
            "fetch_missing": fetch_missing,
            "offline": offline,
        },
    )
    values = resolve_step_paths(
        values,
        ("entries_path", "alignment_path", "pdb_dir", "alphafold_dir", "pocket_dir", "fsdb_dir", "pockets_path"),
        ("log_path", "pockets_tsv_path"),
    )
    with run_scope("pockets"), log_to_file(values["log_path"], values["verbosity"]):
        pisa_source = resolve_pisa_source(values["pisa_source"])
        pisa_bsa_threshold = resolve_pisa_bsa_threshold(values["pisa_bsa_threshold"])
        fetch_missing = resolve_fetch_missing(values["fetch_missing"])
        offline = resolve_offline(values["offline"])
        sides = open_entries(values["entries_path"], "pockets")
        if fetch_missing and not offline:
            fetch_missing_entries(sides, StructureManager())
        build_pockets(
            sides,
            values["entries_path"],
            values["alignment_path"],
            values["pockets_path"],
            values,
            pisa_source,
            values["pockets_tsv_path"],
            pisa_bsa_threshold=pisa_bsa_threshold,
            download=not offline,
        )


def build_pockets(
    sides,
    entries_path,
    alignment,
    pockets_path,
    cache_dirs,
    pisa_source,
    pockets_tsv_path=None,
    pisa_bsa_threshold=None,
    download=True,
):
    """
    Build a Pocket for every record of every side, and write them to `pockets_path`.

    Marks the records whose structure is missing (`structure_not_found`) or whose pocket could not be
    built (`pocket_not_built`) in the entries table. Foldseek-database records have no pocket. Against a
    PDB-named Foldseek database, also builds a pisa pocket per PISA interface of each hit (see
    `expand_fsdb_pdb_targets`); a hit is no entry, so one that fails is only logged. The pockets file
    lists every record and hit pocket, `null` where it was not built, and maps each `preprocess_name`
    and hit name to its pocket_ids. Downloads any PISA data not already cached under `pocket_dir/pisa/`
    unless `download` is False, and overwrites the per-method pocket files under `pocket_dir`.

    Args:
        sides (dict): "query" and "target" -> that side's QTRecord dicts. Built in order, then the hits,
            which decides which of several identical pisa records names the shared pocket.
        entries_path (str or None): The entries table, marked, or None for none.
        alignment (str): The alignment table. Read only when the target is a Foldseek database.
        pockets_path (str): Where the pockets are written, overwriting any there.
        cache_dirs (dict): Holds "pdb_dir", "alphafold_dir", "fsdb_dir" and "pocket_dir".
        pisa_source (str): Where PISA interfaces are fetched from: "ftp" or "api".
        pockets_tsv_path (str, optional): Where a table of every pocket is written: pocket_id, chain,
            res_auth_ids and method. Defaults to None, for none.
        pisa_bsa_threshold (float, optional): Percent of a residue's accessible surface the partner
            chain must bury, beyond which it joins a pisa pocket; applies to every pisa pocket, the
            Foldseek hits' included. Defaults to None, which builds pisa pockets from interface bonds.
        download (bool, optional): False downloads nothing: neither PISA data nor the structures of a
            PDB-named database's hits, so only what is cached is used. Defaults to True.

    Returns:
        None

    Raises:
        PocketMapperError: If a Foldseek-database target has no alignment file, or a PDB-named database's
            hits have interfaces but none of their structures could be fetched.
    """
    log_extra = {"stage": "Getting Pockets"}

    entry_records = {
        name: [record for record in records if record["struct_type"] != "foldseek_db"]
        for name, records in sides.items()
    }
    given = [(name, record) for name, records in entry_records.items() for record in records]
    chains = preproc_to_ids([record for _, record in given])
    # (side, record) in build order; the side is None for a Foldseek hit
    buildable = []
    for name, records in entry_records.items():
        kept = split_missing_structures(records, name, "pockets", entries_path, log_extra)
        buildable += [(name, record) for record in kept]

    hits = []
    if fsdb_record(sides["target"]) is not None:
        require_file(alignment, "alignment")
        hit_names = pd.read_csv(alignment, sep="\t", usecols=["target"], dtype=str)["target"].unique().tolist()
        if fsdb_pocket_mode(hit_names) == "pisa":
            hits, hit_chains, failures = expand_fsdb_pdb_targets(hit_names, cache_dirs, pisa_source, download=download)
            report_failures(None, "pockets", failures, log_extra, "Missing structures for foldseek hit(s)")
            buildable += [(None, record) for record in hits]
            chains |= hit_chains

    built = fetch_pockets(
        [record for _, record in buildable],
        cache_dirs["pocket_dir"],
        builder_options={
            "pisa": {"pisa_source": pisa_source, "bsa_threshold": pisa_bsa_threshold, "download": download}
        },
    )

    # Looked up by pocket_id in the merged pockets, where a later method's pocket replaces an earlier
    # one's of the same id. Two methods can share an id only across sides or with a Foldseek hit, since
    # each side chooses its method independently.
    failures = [
        failure(name, record, "pocket_not_built")
        for name, record in buildable
        if built.get(record["pocket_id"]) is None
    ]
    report_failures(entries_path, "pockets", failures, log_extra, "No pocket built, skipping")

    # Every pocket given or expanded is listed, None where it was not built
    pocket_ids = [record["pocket_id"] for _, record in given] + [pid for ids in chains.values() for pid in ids]
    pockets = {pocket_id: built.get(pocket_id) for pocket_id in pocket_ids}
    make_dir(os.path.dirname(pockets_path), log_extra)
    write_pockets_file(pockets, chains, pockets_path)
    logger.info(f"Pockets written to {pockets_path}", extra=log_extra)
    if pockets_tsv_path is not None:
        write_pockets_table([record for _, record in given] + hits, pockets, pockets_tsv_path)


def write_pockets_table(records, pockets, path):
    """
    Write a table of every pocket: pocket_id, chain, res_auth_ids and method.

    Args:
        records (list): The QTRecord dicts the pockets were built for; a pocket_id listed twice (on both
            sides, or as a given record and a Foldseek hit) is written once, from its first record.
        pockets (dict): pocket_id -> Pocket, or None for one not built, whose res_auth_ids are left
            empty.
        path (str): The table to write, tab-separated.

    Returns:
        None

    Raises:
        PocketMapperError: If the table's directory cannot be created.
    """
    log_extra = {"stage": "Getting Pockets"}
    rows = {}
    for record in records:
        pocket = pockets.get(record["pocket_id"])
        rows.setdefault(
            record["pocket_id"],
            {
                "pocket_id": record["pocket_id"],
                "chain": split_chain_info(record["chain_info"])[0],
                "res_auth_ids": ",".join(map(str, pocket.res_auth_ids)) if pocket is not None else "",
                "method": record["pocket_method"],
            },
        )
    make_dir(os.path.dirname(path), log_extra)
    pd.DataFrame(list(rows.values()), columns=["pocket_id", "chain", "res_auth_ids", "method"]).to_csv(
        path, sep="\t", index=False
    )
    logger.info(f"Pockets listed in {path}", extra=log_extra)


def expand_fsdb_pdb_targets(hit_names, cache_dirs, pisa_source, download=True):
    """
    Turn the hits of a PDB Foldseek-database search into pisa records, one per PISA interface.

    Downloads the PISA data and the structures those records need, unless `download` is False. Each
    record's `preprocess_name` is the hit name rather than the one `QTProcessor` derives.

    Args:
        hit_names (list): The database entry names the search hit.
        cache_dirs (dict): Holds "pdb_dir", "alphafold_dir", "fsdb_dir" and "pocket_dir".
        pisa_source (str): Where PISA interfaces are fetched from: "ftp" or "api".
        download (bool, optional): False downloads nothing: a hit entry whose PISA data is not cached
            gives no records, and a record whose structure is not cached is a failure. Defaults to True.

    Returns:
        tuple: (records, chains, failures). `records` holds those whose structure is available.
            `chains` maps every hit name to the pocket_ids of its records, including those whose
            structure could not be fetched; it is empty for a hit that is not PDB-style or has no
            interface the pisa pocket method can spell. `failures` has one `structure_not_found`
            failure per such record, as `lib.entries.failure` builds it, with no side.

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
    interface_dir = download_pisa_interfaces(pdb_list, cache_dirs["pocket_dir"], pisa_source, download=download)

    # Each entry's interfaces are read once, however many hits share it
    hit_chains = {}  # pdb_id -> the chain ids hit on it, as dict keys in hit order
    for pdb_id, chain_id in hits.values():
        hit_chains.setdefault(pdb_id, {})[chain_id] = None
    partners = {}  # (pdb_id, chain_id) -> partner chain ids
    for pdb_id, chain_ids in hit_chains.items():
        interfaces = load_interfaces(pdb_id, interface_dir)
        if interfaces is None:
            logger.debug(f"Could not load PISA data for {pdb_id}", extra=log_extra)
            continue
        for chain_id in chain_ids:
            partners[(pdb_id, chain_id)] = interface_partners(interfaces, chain_id)
            if not partners[(pdb_id, chain_id)]:
                logger.debug(f"No PISA interface involving chain {chain_id} of {pdb_id}", extra=log_extra)

    # Building one record per interface the hit chain takes part in
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
        for partner in partners.get((pdb_id, chain_id), []):
            if not re.match(pisa_pattern, f"{chain_id}_{partner}"):
                unspellable += 1
                continue
            record, reason = qtprocessor.parse_individual_qt(f"{pdb_id}:{chain_id}_{partner}", pocket_method="pisa")
            if record is None:
                logger.warning(f"{reason}; skipping this interface", extra=log_extra)
                continue
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
    found = fetch_missing_structures("foldseek hit", records, StructureManager(download=download))
    if not any(found.values()):
        logger.critical("Insufficient foldseek hit structures after fetching", extra=log_extra)
        raise PocketMapperError(
            "Insufficient foldseek hit structures after fetching. No valid foldseek hit entries remain."
        )
    failures = [failure(None, record, "structure_not_found") for record in records if not found[record["struct_info"]]]
    records = [record for record in records if found[record["struct_info"]]]

    logger.info(
        f"Added {len(records)} PISA target pockets from {len({r['preprocess_name'] for r in records})} Foldseek hits",
        extra=log_extra,
    )
    return records, chains, failures
