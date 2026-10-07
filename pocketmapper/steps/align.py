"""
Step 3: align the query chains against the target chains into an alignment table.

With the foldseek aligner each chain is first cut out of its structure and cached. A Foldseek
database target is searched as it is; its hits' pockets are the pockets step's business.
"""

import logging
import os

from pocketmapper.constants import FOLDSEEK_FORMAT_OUTPUT
from pocketmapper.exceptions import PocketMapperError
from pocketmapper.foldseek import run_foldseek
from pocketmapper.lib import log_to_file
from pocketmapper.lib import make_dir
from pocketmapper.lib import temp_dir_scope
from pocketmapper.records import append_failed_entries
from pocketmapper.records import failed_entry
from pocketmapper.records import fsdb_record
from pocketmapper.records import read_cache_manifest
from pocketmapper.records import read_records
from pocketmapper.records import unique_by
from pocketmapper.records import write_records
from pocketmapper.sequence_aligner import SequenceAligner
from pocketmapper.settings import check_fsdb_aligner
from pocketmapper.settings import layer_settings
from pocketmapper.settings import require_foldseek
from pocketmapper.settings import require_setting
from pocketmapper.settings import resolve_aligner
from pocketmapper.settings import resolve_delete_tmp
from pocketmapper.settings import resolve_paths
from pocketmapper.settings import resolve_threads
from pocketmapper.structure_preprocessor import StructurePreprocessor

logger = logging.getLogger(__name__)


def align(
    job_file=None,
    results_dir=None,
    verbosity=None,
    log_path=None,
    failed_entries_path=None,
    query_records=None,
    target_records=None,
    query_records_path=None,
    target_records_path=None,
    alignment_path=None,
    aligner=None,
    threads=None,
    temp_dir=None,
    delete_tmp=None,
):
    """
    Align the query chains against the target chains into an alignment table.

    Reads the cache directories from `results_dir`'s cache manifest. Drops the records whose
    structure is missing or cannot be preprocessed, adding them to `failed_entries_path`. Empties
    `temp_dir` on the way in and, unless `delete_tmp` is 0, deletes it on the way out, unless an
    enclosing call holds it.

    Args:
        job_file (str or dict, optional): JSON job file of Settings field name -> value, or the same
            already loaded. Any argument given overrides it.
        results_dir (str, optional): The results directory `parse` wrote to. Required here or in
            `job_file`.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Defaults to <results_dir>/info.log.
        failed_entries_path (str, optional): Defaults to <results_dir>/failed_entries.json.
        query_records (str, optional): Defaults to the job file's query_records_path, else
            <results_dir>/query_records.json.
        target_records (str, optional): As `query_records`, for the target side.
        query_records_path (str, optional): Where the query records left are written. Defaults to
            `query_records`.
        target_records_path (str, optional): As `query_records_path`. Defaults to `target_records`.
        alignment_path (str, optional): Defaults to <results_dir>/alignment.tsv.
        aligner (str, optional): "foldseek" or "seq". Defaults to DEFAULT_ALIGNER.
        threads (int, optional): Defaults to one per available core.
        temp_dir (str, optional): Defaults to <results_dir>/tmp.
        delete_tmp (int, optional): 1 deletes `temp_dir` at the end; 0 keeps it. Defaults to
            DEFAULT_DELETE_TMP.

    Returns:
        None

    Raises:
        PocketMapperError: If the job file cannot be read, no results_dir is given, the manifest or a
            records file is missing, a setting is invalid, foldseek is needed but cannot run, a
            Foldseek-database target is aligned with "seq", a side has no usable records, or a
            Foldseek invocation fails.
    """
    # The records paths rewrite the inputs, so they are not layered: an explicit one would become
    # its input's default too
    values = layer_settings(
        job_file,
        {
            "results_dir": results_dir,
            "verbosity": verbosity,
            "log_path": log_path,
            "failed_entries_path": failed_entries_path,
            "alignment_path": alignment_path,
            "aligner": aligner,
            "threads": threads,
            "temp_dir": temp_dir,
            "delete_tmp": delete_tmp,
        },
    )
    require_setting(values, "results_dir")
    values = resolve_paths(values)
    query_records = query_records if query_records is not None else values["query_records_path"]
    target_records = target_records if target_records is not None else values["target_records_path"]
    with log_to_file(values["log_path"], values["verbosity"]):
        cache_dirs = read_cache_manifest(values["results_dir"])
        aligner = resolve_aligner(values["aligner"])
        if fsdb_record(read_records(target_records)) is not None:
            check_fsdb_aligner(aligner)
        if aligner == "foldseek":
            require_foldseek(
                "The foldseek aligner was selected",
                "Or pass --aligner seq to use the built-in BLOSUM62 sequence aligner.",
            )
        threads = resolve_threads(values["threads"])
        delete_tmp = resolve_delete_tmp(values["delete_tmp"])

        roots = [cache_dirs["cache_dir"], values["results_dir"]]
        with temp_dir_scope(values["temp_dir"], delete_tmp, roots):
            align_chains(
                query_records,
                target_records,
                query_records_path if query_records_path is not None else query_records,
                target_records_path if target_records_path is not None else target_records,
                values["alignment_path"],
                values["failed_entries_path"],
                aligner,
                threads,
                values["verbosity"],
                values["temp_dir"],
                cache_dirs["foldseek_preprocessed_structure_dir"],
            )


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
    preprocessed_dir,
):
    """
    Align every query chain against every target chain and write the table to `alignment_path`.

    Drops, and adds to `failed_entries_path`, each record whose structure is missing
    (`structure_not_found`) or cannot be cut down to its chain (`structure_preprocessing_failed`).

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
        preprocessed_dir (str): Cache of the single-chain structures Foldseek is given.

    Returns:
        None

    Raises:
        PocketMapperError: If a records file cannot be read, a side has no usable records, or a
            Foldseek invocation fails.
    """
    log_extra = {"stage": "Alignment"}

    sides = {"query": read_records(query_records), "target": read_records(target_records)}
    sources = {"query": query_records, "target": target_records}

    database = fsdb_record(sides["target"])

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
                sides[name], preprocessed_dir, tmp_dirs[name], sources[name]
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
