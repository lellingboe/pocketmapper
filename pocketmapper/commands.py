"""
One function per pipeline step, each runnable on its own: `parse`, `fetch`, `align`, `pockets`,
`compare` and `superpose`. Run in that order in one `results_dir`, they do what
`PocketMapper.search` does.

Each command resolves its path defaults under `results_dir`, logs to `log_path` for the length of
the call, validates its settings, manages its scratch directory and calls its step in
`pocketmapper.steps`. A command reads only files: the previous step's outputs and, for `fetch`,
`align` and `pockets`, the cache manifest `parse` wrote, so a chain cannot split its cache.

Defaults: an input defaults to the standard file under `results_dir`, and an output that rewrites
an input (`fetch` and `align` records) defaults to that input's path.
"""

import logging
import os
from contextlib import contextmanager

import pandas as pd

from pocketmapper.constants import DEFAULT_ALIGN_COUNT
from pocketmapper.constants import DEFAULT_ALIGN_STRUCT_METHOD
from pocketmapper.constants import DEFAULT_ALIGNER
from pocketmapper.constants import DEFAULT_CACHE_DIR
from pocketmapper.constants import DEFAULT_DELETE_TMP
from pocketmapper.constants import DEFAULT_PISA_SOURCE
from pocketmapper.constants import DEFAULT_POCKET_METHOD
from pocketmapper.constants import DEFAULT_VERBOSITY
from pocketmapper.lib import delete_temp_dir
from pocketmapper.lib import empty_temp_dir
from pocketmapper.lib import log_to_file
from pocketmapper.lib import make_dir
from pocketmapper.records import fsdb_record
from pocketmapper.records import read_cache_manifest
from pocketmapper.records import read_records
from pocketmapper.records import require_file
from pocketmapper.settings import cache_path
from pocketmapper.settings import check_fsdb_align_struct_method
from pocketmapper.settings import check_fsdb_aligner
from pocketmapper.settings import default_results_dir
from pocketmapper.settings import require_foldseek
from pocketmapper.settings import resolve_align_struct_method
from pocketmapper.settings import resolve_aligner
from pocketmapper.settings import resolve_delete_tmp
from pocketmapper.settings import resolve_pisa_source
from pocketmapper.settings import resolve_threads
from pocketmapper.settings import results_path
from pocketmapper.steps.align import align_chains
from pocketmapper.steps.compare import compare_aligned_pockets
from pocketmapper.steps.fetch import fetch_inputs
from pocketmapper.steps.parse import parse_inputs
from pocketmapper.steps.pockets import build_pockets
from pocketmapper.steps.superpose import superpose_top_targets

logger = logging.getLogger(__name__)


@contextmanager
def command_log(log_path, verbosity):
    """
    Log the package to `log_path` for the length of a `with` block, creating its directory.

    The file is appended to, so a chain of commands in one results directory writes one log.

    Args:
        log_path (str): The log file.
        verbosity (int): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR.

    Yields:
        None

    Raises:
        PocketMapperError: If the log's directory cannot be created.
    """
    make_dir(os.path.dirname(log_path), {"stage": "Configuring Settings"})
    with log_to_file(log_path, verbosity):
        yield


def parse(
    query,
    target,
    results_dir=None,
    verbosity=DEFAULT_VERBOSITY,
    log_path=None,
    failed_entries_path=None,
    query_pocket_method=DEFAULT_POCKET_METHOD,
    target_pocket_method=DEFAULT_POCKET_METHOD,
    cache_dir=DEFAULT_CACHE_DIR,
    pdb_dir=None,
    alphafold_dir=None,
    pocket_dir=None,
    foldseek_preprocessed_structure_dir=None,
    fsdb_dir=None,
    query_records_path=None,
    target_records_path=None,
):
    """
    Parse the query and target inputs into records files. No network.

    Writes the records files and `cache_dirs.json`, naming every cache directory absolute, and
    starts `failed_entries_path` afresh with the entries that could not be parsed.

    Args:
        query (str): Query entry, or a file of one entry per line.
        target (str): Target entry, a file of them, or a Foldseek database.
        results_dir (str, optional): Defaults to pocketmapper_results_<YYMMDD_HHMMSS>.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Defaults to <results_dir>/info.log.
        failed_entries_path (str, optional): Defaults to <results_dir>/failed_entries.json.
        query_pocket_method (str, optional): Pocket method to force on every query entry, or "auto"
            (the default) to infer it per entry.
        target_pocket_method (str, optional): As `query_pocket_method`, for the target side.
        cache_dir (str, optional): Defaults to DEFAULT_CACHE_DIR.
        pdb_dir (str, optional): Defaults to <cache_dir>/pdb_structures.
        alphafold_dir (str, optional): Defaults to <cache_dir>/alphafold_structures.
        pocket_dir (str, optional): Defaults to <cache_dir>/pockets.
        foldseek_preprocessed_structure_dir (str, optional): Defaults to
            <cache_dir>/foldseek_preprocessed_structures.
        fsdb_dir (str, optional): Defaults to <cache_dir>/fsdb.
        query_records_path (str, optional): Defaults to <results_dir>/query_records.json.
        target_records_path (str, optional): Defaults to <results_dir>/target_records.json.

    Returns:
        None

    Raises:
        PocketMapperError: If a pocket method is unknown, a side has no valid entries, or a
            Foldseek-database target is not the only target entry.
    """
    if results_dir is None:
        results_dir = default_results_dir()
    cache_dirs = {
        "cache_dir": cache_dir,
        "pdb_dir": cache_path(cache_dir, "pdb_dir", pdb_dir),
        "alphafold_dir": cache_path(cache_dir, "alphafold_dir", alphafold_dir),
        "pocket_dir": cache_path(cache_dir, "pocket_dir", pocket_dir),
        "foldseek_preprocessed_structure_dir": cache_path(
            cache_dir, "foldseek_preprocessed_structure_dir", foldseek_preprocessed_structure_dir
        ),
        "fsdb_dir": cache_path(cache_dir, "fsdb_dir", fsdb_dir),
    }
    with command_log(results_path(results_dir, "log_path", log_path), verbosity):
        parse_inputs(
            query,
            target,
            query_pocket_method,
            target_pocket_method,
            cache_dirs,
            results_dir,
            results_path(results_dir, "query_records_path", query_records_path),
            results_path(results_dir, "target_records_path", target_records_path),
            results_path(results_dir, "failed_entries_path", failed_entries_path),
        )


def fetch(
    results_dir,
    verbosity=DEFAULT_VERBOSITY,
    log_path=None,
    failed_entries_path=None,
    query_records=None,
    target_records=None,
    query_records_path=None,
    target_records_path=None,
    pisa_source=DEFAULT_PISA_SOURCE,
    threads=None,
    temp_dir=None,
    delete_tmp=DEFAULT_DELETE_TMP,
):
    """
    Download the structures, Foldseek database and PISA interfaces the records need.

    Reads the cache directories from `results_dir`'s cache manifest. Drops the records whose
    structure cannot be fetched, adding them to `failed_entries_path`. Empties `temp_dir` on the way
    in and, unless `delete_tmp` is 0, deletes it on the way out.

    Args:
        results_dir (str): The results directory `parse` wrote to.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Defaults to <results_dir>/info.log.
        failed_entries_path (str, optional): Defaults to <results_dir>/failed_entries.json.
        query_records (str, optional): Defaults to <results_dir>/query_records.json.
        target_records (str, optional): Defaults to <results_dir>/target_records.json.
        query_records_path (str, optional): Where the query records left are written. Defaults to
            `query_records`.
        target_records_path (str, optional): As `query_records_path`. Defaults to `target_records`.
        pisa_source (str, optional): "ftp" or "api". Defaults to DEFAULT_PISA_SOURCE.
        threads (int, optional): Defaults to one per available core.
        temp_dir (str, optional): Defaults to <results_dir>/tmp.
        delete_tmp (int, optional): 1 deletes `temp_dir` at the end; 0 keeps it. Defaults to
            DEFAULT_DELETE_TMP.

    Returns:
        None

    Raises:
        PocketMapperError: If the manifest or a records file is missing, a setting is invalid, no
            structure for a side could be fetched, or a Foldseek database cannot be downloaded.
    """
    query_records = results_path(results_dir, "query_records_path", query_records)
    target_records = results_path(results_dir, "target_records_path", target_records)
    temp_dir = results_path(results_dir, "temp_dir", temp_dir)
    with command_log(results_path(results_dir, "log_path", log_path), verbosity):
        cache_dirs = read_cache_manifest(results_dir)
        threads = resolve_threads(threads)
        delete_tmp = resolve_delete_tmp(delete_tmp)
        pisa_source = resolve_pisa_source(pisa_source)

        roots = [cache_dirs["cache_dir"], results_dir]
        empty_temp_dir(temp_dir, roots)
        fetch_inputs(
            query_records,
            target_records,
            query_records_path or query_records,
            target_records_path or target_records,
            results_path(results_dir, "failed_entries_path", failed_entries_path),
            cache_dirs["pocket_dir"],
            pisa_source,
            threads,
            temp_dir,
        )
        delete_temp_dir(temp_dir, delete_tmp, roots)


def align(
    results_dir,
    verbosity=DEFAULT_VERBOSITY,
    log_path=None,
    failed_entries_path=None,
    query_records=None,
    target_records=None,
    query_records_path=None,
    target_records_path=None,
    alignment_path=None,
    aligner=DEFAULT_ALIGNER,
    threads=None,
    pisa_source=DEFAULT_PISA_SOURCE,
    temp_dir=None,
    delete_tmp=DEFAULT_DELETE_TMP,
):
    """
    Align the query chains against the target chains into an alignment table.

    Reads the cache directories from `results_dir`'s cache manifest. Drops the records whose
    structure is missing or cannot be preprocessed, adding them to `failed_entries_path`. Against a
    PDB Foldseek database, appends a pisa target record per interface of each hit, downloading what
    they need. Empties `temp_dir` on the way in and, unless `delete_tmp` is 0, deletes it on the way
    out.

    Args:
        results_dir (str): The results directory `parse` wrote to.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Defaults to <results_dir>/info.log.
        failed_entries_path (str, optional): Defaults to <results_dir>/failed_entries.json.
        query_records (str, optional): Defaults to <results_dir>/query_records.json.
        target_records (str, optional): Defaults to <results_dir>/target_records.json.
        query_records_path (str, optional): Where the query records left are written. Defaults to
            `query_records`.
        target_records_path (str, optional): As `query_records_path`. Defaults to `target_records`.
        alignment_path (str, optional): Defaults to <results_dir>/alignment.tsv.
        aligner (str, optional): "foldseek" or "seq". Defaults to DEFAULT_ALIGNER.
        threads (int, optional): Defaults to one per available core.
        pisa_source (str, optional): "ftp" or "api", for a PDB Foldseek database's hits. Defaults to
            DEFAULT_PISA_SOURCE.
        temp_dir (str, optional): Defaults to <results_dir>/tmp.
        delete_tmp (int, optional): 1 deletes `temp_dir` at the end; 0 keeps it. Defaults to
            DEFAULT_DELETE_TMP.

    Returns:
        None

    Raises:
        PocketMapperError: If the manifest or a records file is missing, a setting is invalid,
            foldseek is needed but cannot run, a Foldseek-database target is aligned with "seq", a
            side has no usable records, or a Foldseek invocation fails.
    """
    query_records = results_path(results_dir, "query_records_path", query_records)
    target_records = results_path(results_dir, "target_records_path", target_records)
    temp_dir = results_path(results_dir, "temp_dir", temp_dir)
    with command_log(results_path(results_dir, "log_path", log_path), verbosity):
        cache_dirs = read_cache_manifest(results_dir)
        aligner = resolve_aligner(aligner)
        if fsdb_record(read_records(target_records)) is not None:
            check_fsdb_aligner(aligner)
        if aligner == "foldseek":
            require_foldseek(
                "The foldseek aligner was selected",
                "Or pass --aligner seq to use the built-in BLOSUM62 sequence aligner.",
            )
        threads = resolve_threads(threads)
        pisa_source = resolve_pisa_source(pisa_source)
        delete_tmp = resolve_delete_tmp(delete_tmp)

        roots = [cache_dirs["cache_dir"], results_dir]
        empty_temp_dir(temp_dir, roots)
        align_chains(
            query_records,
            target_records,
            query_records_path or query_records,
            target_records_path or target_records,
            results_path(results_dir, "alignment_path", alignment_path),
            results_path(results_dir, "failed_entries_path", failed_entries_path),
            aligner,
            threads,
            verbosity,
            temp_dir,
            cache_dirs,
            pisa_source,
        )
        delete_temp_dir(temp_dir, delete_tmp, roots)


def pockets(
    records,
    results_dir,
    verbosity=DEFAULT_VERBOSITY,
    log_path=None,
    failed_entries_path=None,
    pockets_path=None,
    pisa_source=DEFAULT_PISA_SOURCE,
):
    """
    Build the pocket of every record in some records files.

    Reads the pocket cache directory from `results_dir`'s cache manifest. Overwrites `pockets_path`
    with the pockets of every file named, and rewrites each file in place without the records whose
    pocket could not be built, adding those to `failed_entries_path`. Downloads any PISA data not
    already cached.

    Args:
        records (list): The records files. Name both sides' files in one call: the pockets file
            holds only the pockets of the files named.
        results_dir (str): The results directory `parse` wrote to.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Defaults to <results_dir>/info.log.
        failed_entries_path (str, optional): Defaults to <results_dir>/failed_entries.json.
        pockets_path (str, optional): Defaults to <results_dir>/pockets.json.
        pisa_source (str, optional): "ftp" or "api". Defaults to DEFAULT_PISA_SOURCE.

    Returns:
        None

    Raises:
        PocketMapperError: If the manifest or a records file is missing, or `pisa_source` is invalid.
    """
    with command_log(results_path(results_dir, "log_path", log_path), verbosity):
        cache_dirs = read_cache_manifest(results_dir)
        pisa_source = resolve_pisa_source(pisa_source)
        build_pockets(
            records,
            results_path(results_dir, "pockets_path", pockets_path),
            results_path(results_dir, "failed_entries_path", failed_entries_path),
            cache_dirs["pocket_dir"],
            pisa_source,
        )


def compare(
    results_dir,
    verbosity=DEFAULT_VERBOSITY,
    log_path=None,
    query_records=None,
    target_records=None,
    alignment=None,
    pockets=None,
    pocket_comparison_path=None,
):
    """
    Compare the pockets of every aligned query/target pair into a pocket comparison table.

    Writes unknown_ids.json and incorrect_mapping.json to `results_dir` when either has anything to
    report, deleting any left there by an earlier run.

    Args:
        results_dir (str): The results directory the inputs default to.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Defaults to <results_dir>/info.log.
        query_records (str, optional): Defaults to <results_dir>/query_records.json.
        target_records (str, optional): Defaults to <results_dir>/target_records.json.
        alignment (str, optional): Defaults to <results_dir>/alignment.tsv.
        pockets (str, optional): Defaults to <results_dir>/pockets.json.
        pocket_comparison_path (str, optional): Defaults to <results_dir>/pocket_comparison.tsv.

    Returns:
        None

    Raises:
        PocketMapperError: If an input is missing or unreadable, or a record has no pocket.
    """
    with command_log(results_path(results_dir, "log_path", log_path), verbosity):
        compare_aligned_pockets(
            results_path(results_dir, "query_records_path", query_records),
            results_path(results_dir, "target_records_path", target_records),
            results_path(results_dir, "alignment_path", alignment),
            results_path(results_dir, "pockets_path", pockets),
            results_path(results_dir, "pocket_comparison_path", pocket_comparison_path),
            results_dir,
        )


def superpose(
    results_dir,
    verbosity=DEFAULT_VERBOSITY,
    log_path=None,
    query_records=None,
    target_records=None,
    pocket_comparison=None,
    alignment=None,
    aligned_structure_dir=None,
    align_struct_method=DEFAULT_ALIGN_STRUCT_METHOD,
    align_count=DEFAULT_ALIGN_COUNT,
    threads=None,
):
    """
    Superpose the top targets of each query onto it, one PDB per query.

    The aligner is read off the alignment table: the seq aligner leaves its `u` column "-". Against
    a Foldseek database, the target structures are rebuilt out of it, which needs foldseek.

    Args:
        results_dir (str): The results directory the inputs default to.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Defaults to <results_dir>/info.log.
        query_records (str, optional): Defaults to <results_dir>/query_records.json.
        target_records (str, optional): Defaults to <results_dir>/target_records.json.
        pocket_comparison (str, optional): Defaults to <results_dir>/pocket_comparison.tsv.
        alignment (str, optional): Defaults to <results_dir>/alignment.tsv.
        aligned_structure_dir (str, optional): Defaults to <results_dir>/aligned_structures.
        align_struct_method (str, optional): "pocket", "foldseek" or "auto", which picks "foldseek"
            for a Foldseek alignment and "pocket" for a seq one. Defaults to
            DEFAULT_ALIGN_STRUCT_METHOD.
        align_count (int, optional): Most targets to superpose onto each query; 0 disables.
            Defaults to DEFAULT_ALIGN_COUNT.
        threads (int, optional): Defaults to one per available core.

    Returns:
        None

    Raises:
        PocketMapperError: If an input is missing, a setting is invalid or does not suit the
            alignment or target, or foldseek is needed but cannot run.
    """
    log_extra = {"stage": "Structural Alignment"}

    alignment = results_path(results_dir, "alignment_path", alignment)
    pocket_comparison = results_path(results_dir, "pocket_comparison_path", pocket_comparison)
    target_records = results_path(results_dir, "target_records_path", target_records)
    with command_log(results_path(results_dir, "log_path", log_path), verbosity):
        require_file(alignment, "alignment")
        require_file(pocket_comparison, "pocket comparison")
        transforms = pd.read_csv(alignment, sep="\t", usecols=["u"], dtype=str)["u"]
        if transforms.empty:
            logger.warning(f"{alignment} holds no alignments; nothing to superpose", extra=log_extra)
            return
        aligner = "seq" if transforms.iloc[0] == "-" else "foldseek"
        align_struct_method = resolve_align_struct_method(align_struct_method, aligner)
        if fsdb_record(read_records(target_records)) is not None:
            check_fsdb_align_struct_method(align_struct_method)
            require_foldseek("Rebuilding the target structures out of the Foldseek database needs foldseek")
        threads = resolve_threads(threads)

        superpose_top_targets(
            results_path(results_dir, "query_records_path", query_records),
            target_records,
            pocket_comparison,
            alignment,
            results_path(results_dir, "aligned_structure_dir", aligned_structure_dir),
            align_struct_method,
            align_count,
            threads,
        )
