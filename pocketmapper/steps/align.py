"""
Step 3: align the query chains against the target chains into an alignment table.

With the foldseek aligner each chain is first cut out of its structure and cached. Against a
Foldseek database, the hits then decide how target pockets are resolved: a PDB-named database's
hits become pisa target records, any other database's pockets are synthesised from the alignment.
"""

import logging
import os
import re
from dataclasses import asdict

import pandas as pd

from pocketmapper.constants import FOLDSEEK_FORMAT_OUTPUT
from pocketmapper.exceptions import PocketMapperError
from pocketmapper.foldseek import run_foldseek
from pocketmapper.lib import make_dir
from pocketmapper.lib import parse_foldseek_pdb_entry_name
from pocketmapper.pockets.pisa import PisaParser
from pocketmapper.pockets.pisa import download_pisa_interfaces
from pocketmapper.qt_processor import QTProcessor
from pocketmapper.records import append_failed_entries
from pocketmapper.records import failed_entry
from pocketmapper.records import fsdb_record
from pocketmapper.records import read_records
from pocketmapper.records import unique_by
from pocketmapper.records import write_records
from pocketmapper.sequence_aligner import SequenceAligner
from pocketmapper.steps.fetch import fetch_missing_structures
from pocketmapper.structure_preprocessor import StructurePreprocessor

logger = logging.getLogger(__name__)


def align_chains(
    query_records,
    target_records,
    query_records_path,
    target_records_path,
    alignment_path,
    failed_entries_path,
    aligner,
    threads,
    verbosity,
    temp_dir,
    cache_dirs,
    pisa_source,
):
    """
    Align every query chain against every target chain and write the table to `alignment_path`.

    Drops, and adds to `failed_entries_path`, each record whose structure is missing
    (`structure_not_found`) or cannot be cut down to its chain (`structure_preprocessing_failed`).
    Against a Foldseek database, sets its record's `fsdb_pockets` and, for a PDB-named database,
    appends a pisa record per PISA interface of each hit, downloading their PISA data and structures.
    Records a previous expansion appended are replaced, so the step can rerun on its own output.

    Args:
        query_records (str): The query records file.
        target_records (str): The target records file.
        query_records_path (str): Where the query records left are written. May be `query_records`.
        target_records_path (str): As `query_records_path`, for the target side.
        alignment_path (str): Where the alignment table is written.
        failed_entries_path (str): The failed-entries file, appended to.
        aligner (str): "foldseek" or "seq". A Foldseek-database target needs "foldseek".
        threads (int): Thread count for Foldseek.
        verbosity (int): The run's verbosity; Foldseek's own is capped at 3.
        temp_dir (str): Scratch directory; Foldseek's inputs and databases are built under it.
        cache_dirs (dict): Each of `records.CACHE_MANIFEST_KEYS` -> its absolute directory. Expanded
            records' paths are resolved against them, and their PISA data cached in `pocket_dir`.
        pisa_source (str): Where expanded records' PISA interfaces are fetched from: "ftp" or "api".

    Returns:
        None

    Raises:
        PocketMapperError: If a records file cannot be read, a side has no usable records, a
            Foldseek invocation fails, or every expanded hit's structure fails to download.
    """
    log_extra = {"stage": "Alignment"}

    sides = {"query": read_records(query_records), "target": read_records(target_records)}
    sources = {"query": query_records, "target": target_records}

    # Only the database record survives a previous expansion
    database = fsdb_record(sides["target"])
    if database is not None:
        database = dict(database, fsdb_pockets=None)
        sides["target"] = [database]

    # A structure fetch never produced, or one removed since
    failures = []
    for name, records in sides.items():
        sides[name] = []
        for record in records:
            if record["struct_type"] == "foldseek_db" or os.path.isfile(record["struct_path"]):
                sides[name].append(record)
                continue
            logger.warning(
                f"No structure at {record['struct_path']} for {name} {record['pocket_id']}; skipping it",
                extra=log_extra,
            )
            failures.append(failed_entry(record["pocket_id"], "align", "structure_not_found", sources[name], record))

    tmp_dirs = {
        "query": os.path.join(temp_dir, "query_structures"),
        "target": os.path.join(temp_dir, "target_structures"),
    }
    if aligner == "foldseek":
        logger.info("Preprocessing structures for Foldseek...", extra=log_extra)
        for name in ("query",) if database is not None else ("query", "target"):
            sides[name], preprocess_failures = foldseek_preprocessing(
                sides[name], cache_dirs["foldseek_preprocessed_structure_dir"], tmp_dirs[name], sources[name]
            )
            failures += preprocess_failures
        logger.info("Finished preprocessing structures", extra={"stage": "Preprocessing Structures"})
    append_failed_entries(failed_entries_path, failures)

    errors = []
    for name, records in sides.items():
        if not records:
            logger.critical(f"No usable {name} records to align", extra=log_extra)
            errors.append(f"no usable {name} records")
    if errors:
        raise PocketMapperError("; ".join(errors))

    if aligner == "foldseek":
        logger.info("Running Foldseek easy-search...", extra=log_extra)
        target_db_path = database["struct_path"] if database is not None else None
        foldseek_alignment(
            tmp_dirs, target_db_path, alignment_path, os.path.join(temp_dir, "foldseek_tmp"), threads, verbosity
        )
    else:
        logger.info("Running local pairwise aligner...", extra=log_extra)
        local_alignment(sides["query"], sides["target"], alignment_path)

    if database is not None:
        expanded, failures = expand_fsdb_pdb_targets(alignment_path, cache_dirs, pisa_source, target_records)
        append_failed_entries(failed_entries_path, failures)
        if expanded is None:
            database["fsdb_pockets"] = "whole_chain"
        else:
            database["fsdb_pockets"] = "pisa"
            sides["target"] += expanded

    write_records(sides["query"], query_records_path)
    write_records(sides["target"], target_records_path)


def foldseek_preprocessing(records, cache_dir, search_dir, source):
    """
    Write each record's alignment chain as a single-chain structure for Foldseek to index.

    Caches each copy as `<cache_dir>/<preprocess_name>.cif.gz` and copies it into `search_dir`,
    creating it. Records sharing a chain are preprocessed once.

    Args:
        records (list): One side's QTRecord dicts. Foldseek-database records are passed through.
        cache_dir (str): Directory the single-chain copies are cached in.
        search_dir (str): The side's scratch directory for Foldseek's input.
        source (str): The records file they came from, for the failure entries.

    Returns:
        tuple: (records left, failure entries). Every record on a chain that could not be
            preprocessed is dropped as `structure_preprocessing_failed`.

    Raises:
        PocketMapperError: If `search_dir` cannot be created.
    """
    log_extra = {"stage": "Preprocessing Structures"}

    unique_records = unique_by(records, "preprocess_name", "chain_info")
    logger.debug(f"Records to preprocess: {unique_records}", extra=log_extra)
    make_dir(search_dir, log_extra)

    results = StructurePreprocessor().preprocess_records(
        records=unique_records, cache_dir=cache_dir, search_dir=search_dir
    )
    logger.debug(f"Preprocessing results: {results}", extra=log_extra)

    # The preprocessor reports only the record it processed for each chain
    failed_chains = {
        (record["preprocess_name"], record["chain_info"])
        for record in unique_records
        if not results[record["pocket_id"]]
    }
    kept = []
    failures = []
    for record in records:
        if (record["preprocess_name"], record["chain_info"]) in failed_chains:
            failures.append(
                failed_entry(record["pocket_id"], "align", "structure_preprocessing_failed", source, record)
            )
        else:
            kept.append(record)
    return kept, failures


def foldseek_alignment(tmp_dirs, target_db_path, alignment_path, foldseek_tmp_dir, threads, verbosity):
    """
    Build the query (and, unless the target is a database, target) Foldseek DB, then search them.

    Writes the databases into the side's scratch directories.

    Args:
        tmp_dirs (dict): "query" and "target" -> the side's scratch directory of preprocessed chains.
        target_db_path (str): The Foldseek database to search, or None to build one from the
            target chains.
        alignment_path (str): Where the alignment table is written.
        foldseek_tmp_dir (str): Scratch directory for `easy-search`, created if needed.
        threads (int): Thread count for Foldseek.
        verbosity (int): The run's verbosity; Foldseek's own is capped at 3.

    Returns:
        None

    Raises:
        PocketMapperError: If a directory cannot be created or a Foldseek invocation fails.
    """
    log_extra = {"stage": "Foldseek Alignment"}
    logger.info("Running Foldseek alignment...", extra=log_extra)

    query_db_path = os.path.join(tmp_dirs["query"], "query_db")
    run_foldseek(["createdb", tmp_dirs["query"], query_db_path, "--threads", str(threads)], log_extra)

    if target_db_path is not None:
        logger.debug(f"Targeting Foldseek DB at {target_db_path}", extra=log_extra)
    else:
        target_db_path = os.path.join(tmp_dirs["target"], "target_db")
        run_foldseek(["createdb", tmp_dirs["target"], target_db_path, "--threads", str(threads)], log_extra)

    make_dir(os.path.dirname(alignment_path), log_extra)
    make_dir(foldseek_tmp_dir, log_extra)
    query_target_align_cmd = [
        "easy-search",
        query_db_path,
        target_db_path,
        alignment_path,
        foldseek_tmp_dir,
        "--format-output",
        FOLDSEEK_FORMAT_OUTPUT,
        "--format-mode",
        "4",
        "-e",
        "0.001",
        "--file-include",
        r".*\.cif\.gz",
        "--max-seqs",
        "5000",
        "--threads",
        str(threads),
        "-v",
        str(min(3, verbosity)),  # capped at info: Foldseek's debug output is very verbose
    ]
    run_foldseek(query_target_align_cmd, log_extra)
    logger.debug("Foldseek alignment completed successfully", extra=log_extra)


def local_alignment(query_records, target_records, alignment_path):
    """
    Align every query chain against every target chain by sequence.

    Args:
        query_records (list): Query QTRecord dicts. Records sharing a chain are aligned once.
        target_records (list): Target QTRecord dicts, as `query_records`.
        alignment_path (str): Where the alignment table is written.

    Returns:
        None

    Raises:
        PocketMapperError: If the directory of `alignment_path` cannot be created.
    """
    log_extra = {"stage": "Local Alignment"}
    logger.info("Running local sequence alignments...", extra=log_extra)

    alignment = SequenceAligner().align_records(
        unique_by(query_records, "preprocess_name"), unique_by(target_records, "preprocess_name")
    )
    make_dir(os.path.dirname(alignment_path), log_extra)
    alignment.to_csv(alignment_path, index=False, sep="\t")


def expand_fsdb_pdb_targets(alignment_path, cache_dirs, pisa_source, source):
    """
    Turn the hits of a PDB Foldseek-database search into pisa target records.

    Builds one record per PISA interface of each hit chain, downloading the PISA data and the
    structures those records need. Each record's `preprocess_name` is the Foldseek entry name rather
    than the one `QTProcessor` derives, and its preprocess paths are None.

    Args:
        alignment_path (str): The alignment table, whose `target` column names the hits.
        cache_dirs (dict): Each of `records.CACHE_MANIFEST_KEYS` -> its absolute directory.
        pisa_source (str): Where PISA interfaces are fetched from: "ftp" or "api".
        source (str): The target records file, for the failure entries.

    Returns:
        tuple: (records, failure entries). `records` is None when no hit is PDB-named, so the database
            is not a PDB database, and may be empty when no hit has a usable interface. Records whose
            structure cannot be fetched are left out, as `structure_not_found` failures.

    Raises:
        PocketMapperError: If interfaces were found but none of their structures could be fetched.
    """
    log_extra = {"stage": "Expanding Foldseek DB Targets"}

    alignment_df = pd.read_csv(alignment_path, sep="\t", engine="c")
    hits = {}  # foldseek entry name -> (pdb_id, chain_id)
    for hit_name in alignment_df["target"].unique().tolist():
        resolved = parse_foldseek_pdb_entry_name(hit_name)
        if resolved is not None:
            hits[hit_name] = resolved
    if not hits:
        logger.info(
            "Foldseek database target is not a PDB database; keeping whole-chain target pockets",
            extra=log_extra,
        )
        return None, []

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
    if unspellable:
        logger.info(
            f"Skipped {unspellable} PISA interfaces whose chain ids the pisa pocket method cannot express",
            extra=log_extra,
        )
    if not records:
        logger.warning("No PISA interfaces found for any Foldseek hit", extra=log_extra)
        return [], []

    # Fetched last, so only entries with an interface are downloaded. Failed fetches are dropped.
    found = fetch_missing_structures("foldseek hit", records)
    if not any(found.values()):
        logger.critical("Insufficient foldseek hit structures after fetching", extra=log_extra)
        raise PocketMapperError(
            "Insufficient foldseek hit structures after fetching. No valid foldseek hit entries remain."
        )
    failures = [
        failed_entry(record["pocket_id"], "align", "structure_not_found", source, record)
        for record in records
        if not found[record["struct_info"]]
    ]
    records = [record for record in records if found[record["struct_info"]]]

    logger.info(
        f"Added {len(records)} PISA target pockets from {len({r['preprocess_name'] for r in records})} Foldseek hits",
        extra=log_extra,
    )
    return records, failures
