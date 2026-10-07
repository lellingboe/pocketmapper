"""
Step 3: align the query chains against the target chains into an alignment table.

With the foldseek aligner each chain is first cut out of its structure and cached. A Foldseek
database target is searched as it is; its hits' pockets are the pockets step's business.
"""

import logging
import os

from pocketmapper.constants import FOLDSEEK_FORMAT_OUTPUT
from pocketmapper.entries import failed_entry
from pocketmapper.entries import fsdb_record
from pocketmapper.entries import report_failures
from pocketmapper.entries import split_missing_structures
from pocketmapper.entries import unique_by
from pocketmapper.exceptions import PocketMapperError
from pocketmapper.foldseek import run_foldseek
from pocketmapper.lib import log_to_file
from pocketmapper.lib import make_dir
from pocketmapper.lib import run_scope
from pocketmapper.lib import temp_dir_scope
from pocketmapper.sequence_aligner import SequenceAligner
from pocketmapper.settings import check_fsdb_aligner
from pocketmapper.settings import dump_settings
from pocketmapper.settings import layer_settings
from pocketmapper.settings import require_foldseek
from pocketmapper.settings import require_setting
from pocketmapper.settings import resolve_aligner
from pocketmapper.settings import resolve_delete_tmp
from pocketmapper.settings import resolve_paths
from pocketmapper.settings import resolve_threads
from pocketmapper.steps.parse import parse_job_entries
from pocketmapper.structure_preprocessor import StructurePreprocessor

logger = logging.getLogger(__name__)


def align(
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
    alignment_path=None,
    aligner=None,
    threads=None,
    temp_dir=None,
    delete_tmp=None,
):
    """
    Align the query chains against the target chains into an alignment table.

    Parses the entries the job file names. Adds the entries whose structure is missing or cannot be
    preprocessed to `failed_entries_path`. Empties `temp_dir` on the way in and, unless `delete_tmp`
    is 0, deletes it on the way out, unless an enclosing call holds it.

    Args:
        job_file (str or dict, optional): JSON job file of job key -> value, or the same already
            loaded, e.g. parse's settings. Any argument given overrides it. Must set query and target.
        results_dir (str, optional): Required here or in `job_file`.
        work_dir (str, optional): Directory that entries and relative paths resolve against.
            Defaults to the working directory.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Defaults to <results_dir>/info.log.
        job_settings_path (str, optional): Where these settings are written, as a job file for later
            steps. Defaults to <results_dir>/align_settings.json. Not written when run inside search.
        failed_entries_path (str, optional): Defaults to <results_dir>/failed_entries.json.
        cache_dir (str, optional): Defaults to DEFAULT_CACHE_DIR.
        pdb_dir (str, optional): Defaults to <cache_dir>/pdb_structures.
        alphafold_dir (str, optional): Defaults to <cache_dir>/alphafold_structures.
        pocket_dir (str, optional): Defaults to <cache_dir>/pockets.
        foldseek_preprocessed_structure_dir (str, optional): Defaults to
            <cache_dir>/foldseek_preprocessed_structures.
        fsdb_dir (str, optional): Defaults to <cache_dir>/fsdb.
        alignment_path (str, optional): Defaults to <results_dir>/alignment.tsv.
        aligner (str, optional): "foldseek" or "seq". Defaults to DEFAULT_ALIGNER.
        threads (int, optional): Defaults to one per available core.
        temp_dir (str, optional): Defaults to <results_dir>/tmp.
        delete_tmp (int, optional): 1 deletes `temp_dir` at the end; 0 keeps it. Defaults to
            DEFAULT_DELETE_TMP.

    Returns:
        None

    Raises:
        PocketMapperError: If the job file cannot be read, query, target or results_dir is not given,
            the entries are rejected as `parse` rejects them, a setting is invalid, foldseek is
            needed but cannot run, a Foldseek-database target is aligned with "seq", a side has no
            usable records, or a Foldseek invocation fails.
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
            "alignment_path": alignment_path,
            "aligner": aligner,
            "threads": threads,
            "temp_dir": temp_dir,
            "delete_tmp": delete_tmp,
        },
    )
    for key in ("query", "target", "results_dir"):
        require_setting(values, key)
    values = resolve_paths(values, "align")
    with run_scope("align") as outermost, log_to_file(values["log_path"], values["verbosity"]):
        if outermost:
            dump_settings(values)
        aligner = resolve_aligner(values["aligner"])
        sides = parse_job_entries(values, "align")
        if fsdb_record(sides["target"]) is not None:
            check_fsdb_aligner(aligner)
        if aligner == "foldseek":
            require_foldseek(
                "The foldseek aligner was selected",
                "Or pass --aligner seq to use the built-in BLOSUM62 sequence aligner.",
            )
        threads = resolve_threads(values["threads"])
        delete_tmp = resolve_delete_tmp(values["delete_tmp"])

        roots = [values["cache_dir"], values["results_dir"]]
        with temp_dir_scope(values["temp_dir"], delete_tmp, roots):
            align_chains(
                sides,
                {"query": values["query"], "target": values["target"]},
                values["alignment_path"],
                values["failed_entries_path"],
                aligner,
                threads,
                values["verbosity"],
                values["temp_dir"],
                values["foldseek_preprocessed_structure_dir"],
            )


def align_chains(
    sides,
    sources,
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

    Leaves out, adding to `failed_entries_path`, each record whose structure or Foldseek database is
    missing (`structure_not_found`) or whose structure cannot be cut down to its chain
    (`structure_preprocessing_failed`).

    Args:
        sides (dict): "query" and "target" -> that side's QTRecord dicts.
        sources (dict): "query" and "target" -> the input the side was parsed from, for the failure
            entries.
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
        PocketMapperError: If a side has no usable records, or a Foldseek invocation fails.
    """
    log_extra = {"stage": "Alignment"}

    database = fsdb_record(sides["target"])

    # A structure fetch never produced, or one removed since
    usable = {
        name: split_missing_structures(records, "align", sources[name], failed_entries_path, log_extra)
        for name, records in sides.items()
    }

    tmp_dirs = {
        "query": os.path.join(temp_dir, "query_structures"),
        "target": os.path.join(temp_dir, "target_structures"),
    }
    if aligner == "foldseek":
        logger.info("Preprocessing structures for Foldseek...", extra=log_extra)
        for name in ("query",) if database is not None else ("query", "target"):
            usable[name], failures = foldseek_preprocessing(
                usable[name], preprocessed_dir, tmp_dirs[name], sources[name]
            )
            report_failures(failed_entries_path, failures, log_extra, "Structure preprocessing failed, skipping")
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
            tmp_dirs, target_db_path, alignment_path, os.path.join(temp_dir, "foldseek_tmp"), threads, verbosity
        )
    else:
        logger.info("Running local pairwise aligner...", extra=log_extra)
        local_alignment(usable["query"], usable["target"], alignment_path)


def foldseek_preprocessing(records, cache_dir, search_dir, source):
    """
    Write each record's alignment chain as a single-chain structure for Foldseek to index.

    Caches each copy as `<cache_dir>/<preprocess_name>.cif.gz` and copies it into `search_dir`,
    creating it. Records sharing a chain are preprocessed once.

    Args:
        records (list): One side's QTRecord dicts. Foldseek-database records are passed through.
        cache_dir (str): Directory the single-chain copies are cached in.
        search_dir (str): The side's scratch directory for Foldseek's input.
        source (str): The input they were parsed from, for the failure entries.

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
