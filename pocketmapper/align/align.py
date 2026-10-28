"""
Step 3: align the query chains against the target chains into an alignment table.

With the foldseek aligner each chain is first cut out of its structure and cached. A Foldseek
database target is searched as it is; its hits' pockets are the pockets step's business.
"""

import logging
import os

from pocketmapper.align.sequence_aligner import SequenceAligner
from pocketmapper.align.structure_preprocessor import StructurePreprocessor
from pocketmapper.fetch.fetch import fetch_missing_entries, fetch_missing_fsdb
from pocketmapper.lib.constants import FOLDSEEK_FORMAT_OUTPUT, TEMP_SUBDIRS
from pocketmapper.lib.entries import (
    failure,
    fsdb_record,
    open_entries,
    report_failures,
    split_missing_structures,
    unique_by,
)
from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.files import make_dir
from pocketmapper.lib.foldseek import run_foldseek
from pocketmapper.lib.logs import log_to_file
from pocketmapper.lib.scopes import run_scope, temp_dir_scope
from pocketmapper.lib.settings import (
    check_fsdb_aligner,
    layer_settings,
    require_foldseek,
    resolve_aligner,
    resolve_delete_tmp,
    resolve_fetch_missing,
    resolve_step_paths,
    resolve_threads,
)

logger = logging.getLogger(__name__)


def align(
    job_file=None,
    entries=None,
    verbosity=None,
    log_path=None,
    foldseek_preprocessed_structure_dir=None,
    alignment_path=None,
    aligner=None,
    threads=None,
    fetch_missing=None,
    temp_dir=None,
    delete_tmp=None,
):
    """
    Align the query chains against the target chains into an alignment table.

    Marks the entries whose structure is missing or cannot be preprocessed in the entries table.
    Empties the scratch subdirectories of `temp_dir` on the way in and, unless `delete_tmp` is 0,
    deletes them on the way out, unless an enclosing call holds it.

    Args:
        job_file (str or dict, optional): JSON job file of job key -> value, or the same already
            loaded, e.g. parse's settings. Any argument given overrides it.
        entries (str, optional): The entries table. Defaults to the job file's entries_path.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Where the log is appended. Defaults to the job file's, else none.
        foldseek_preprocessed_structure_dir (str, optional): Cache of the single-chain structures
            Foldseek is given. Defaults to the job file's.
        alignment_path (str, optional): Defaults to the job file's.
        aligner (str, optional): "foldseek" or "seq". Defaults to DEFAULT_ALIGNER.
        threads (int, optional): Defaults to one per available core.
        fetch_missing (int, optional): 1 downloads an entry structure or bundled Foldseek database
            missing from the cache, 0 skips the entry as `structure_not_found`. Defaults to
            DEFAULT_FETCH_MISSING.
        temp_dir (str, optional): Defaults to the job file's.
        delete_tmp (int, optional): 1 deletes `temp_dir` at the end; 0 keeps it. Defaults to
            DEFAULT_DELETE_TMP.

    Returns:
        None

    Raises:
        PocketMapperError: If the job file cannot be read, a path is not given, the entries table cannot
            be read or leaves a side empty, a setting is invalid, foldseek is needed but cannot run, a
            Foldseek-database target is aligned with "seq", a side has no usable records, or a Foldseek
            invocation fails.
    """
    values = layer_settings(
        job_file,
        {
            "entries_path": entries,
            "verbosity": verbosity,
            "log_path": log_path,
            "foldseek_preprocessed_structure_dir": foldseek_preprocessed_structure_dir,
            "alignment_path": alignment_path,
            "aligner": aligner,
            "threads": threads,
            "fetch_missing": fetch_missing,
            "temp_dir": temp_dir,
            "delete_tmp": delete_tmp,
        },
    )
    values = resolve_step_paths(
        values,
        ("entries_path", "foldseek_preprocessed_structure_dir", "alignment_path", "temp_dir"),
        ("log_path",),
    )
    with run_scope("align"), log_to_file(values["log_path"], values["verbosity"]):
        aligner = resolve_aligner(values["aligner"])
        sides = open_entries(values["entries_path"], "align")
        if fsdb_record(sides["target"]) is not None:
            check_fsdb_aligner(aligner)
        if aligner == "foldseek":
            require_foldseek(
                "The foldseek aligner was selected",
                "Or pass --aligner seq to use the built-in BLOSUM62 sequence aligner.",
            )
        threads = resolve_threads(values["threads"])
        delete_tmp = resolve_delete_tmp(values["delete_tmp"])
        fetch_missing = resolve_fetch_missing(values["fetch_missing"])

        with temp_dir_scope(values["temp_dir"], delete_tmp):
            database = fsdb_record(sides["target"])
            if fetch_missing:
                fetch_missing_entries(sides)
                if database is not None:
                    fetch_missing_fsdb(database, threads, os.path.join(values["temp_dir"], TEMP_SUBDIRS["foldseek"]))
            align_chains(
                sides,
                values["entries_path"],
                values["alignment_path"],
                aligner,
                threads,
                values["verbosity"],
                values["temp_dir"],
                values["foldseek_preprocessed_structure_dir"],
            )


def align_chains(
    sides,
    entries_path,
    alignment_path,
    aligner,
    threads,
    verbosity,
    temp_dir,
    preprocessed_dir,
):
    """
    Align every query chain against every target chain and write the table to `alignment_path`.

    Leaves out, marking it in the entries table, each record whose structure or Foldseek database is
    missing (`structure_not_found`) or whose structure cannot be cut down to its chain
    (`structure_preprocessing_failed`).

    Args:
        sides (dict): "query" and "target" -> that side's QTRecord dicts.
        entries_path (str or None): The entries table, marked, or None for none.
        alignment_path (str): Where the alignment table is written.
        aligner (str): "foldseek" or "seq". A Foldseek-database target needs "foldseek".
        threads (int): Thread count for Foldseek.
        verbosity (int): The run's verbosity; Foldseek's own is capped at 3.
        temp_dir (str): Scratch directory; Foldseek's inputs and databases are built under it.
        preprocessed_dir (str): Cache of the single-chain structures Foldseek is given.

    Returns:
        None

    Raises:
        PocketMapperError: If a side has no usable records, or a Foldseek invocation fails.
    """
    log_extra = {"stage": "Alignment"}

    database = fsdb_record(sides["target"])

    # A structure fetch never produced, or one removed since
    usable = {
        name: split_missing_structures(records, name, "align", entries_path, log_extra)
        for name, records in sides.items()
    }

    tmp_dirs = {name: os.path.join(temp_dir, TEMP_SUBDIRS[name]) for name in ("query", "target")}
    if aligner == "foldseek":
        logger.info("Preprocessing structures for Foldseek...", extra=log_extra)
        for name in ("query",) if database is not None else ("query", "target"):
            usable[name], failures = foldseek_preprocessing(usable[name], name, preprocessed_dir, tmp_dirs[name])
            report_failures(entries_path, "align", failures, log_extra, "Structure preprocessing failed, skipping")
        logger.info("Finished preprocessing structures", extra={"stage": "Preprocessing Structures"})

    errors = []
    for name, records in usable.items():
        if not records:
            logger.critical(f"No usable {name} records to align", extra=log_extra)
            errors.append(f"no usable {name} records")
    if errors:
        raise PocketMapperError("; ".join(errors))

    if aligner == "foldseek":
        logger.info("Running Foldseek easy-search...", extra=log_extra)
        target_db_path = database["struct_path"] if database is not None else None
        foldseek_alignment(
            tmp_dirs,
            target_db_path,
            alignment_path,
            os.path.join(temp_dir, TEMP_SUBDIRS["foldseek"]),
            threads,
            verbosity,
        )
    else:
        logger.info("Running local pairwise aligner...", extra=log_extra)
        local_alignment(usable["query"], usable["target"], alignment_path)


def foldseek_preprocessing(records, side, cache_dir, search_dir):
    """
    Write each record's alignment chain as a single-chain structure for Foldseek to index.

    Caches each copy as `<cache_dir>/<preprocess_name>.cif.gz` and copies it into `search_dir`,
    creating it. Records sharing a chain are preprocessed once.

    Args:
        records (list): One side's QTRecord dicts. Foldseek-database records are passed through.
        side (str): Which side they are, for the failures.
        cache_dir (str): Directory the single-chain copies are cached in.
        search_dir (str): The side's scratch directory for Foldseek's input.

    Returns:
        tuple: (records left, failures, as `lib.entries.failure` builds them). Every record on a chain
            that could not be preprocessed is dropped as `structure_preprocessing_failed`.

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
            failures.append(failure(side, record, "structure_preprocessing_failed"))
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
