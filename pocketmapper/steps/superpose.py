"""
Step 6: superpose the top targets of each query onto it, one PDB per query.
"""

import logging

import pandas as pd

from pocketmapper.lib import fsdb_pocket_mode
from pocketmapper.lib import log_to_file
from pocketmapper.lib import run_scope
from pocketmapper.pockets.pocket_fetcher import read_pockets_file
from pocketmapper.records import fsdb_record
from pocketmapper.records import require_file
from pocketmapper.records import split_missing_structures
from pocketmapper.settings import check_fsdb_align_struct_method
from pocketmapper.settings import dump_settings
from pocketmapper.settings import input_path
from pocketmapper.settings import layer_settings
from pocketmapper.settings import require_foldseek
from pocketmapper.settings import require_setting
from pocketmapper.settings import resolve_align_struct_method
from pocketmapper.settings import resolve_paths
from pocketmapper.settings import resolve_threads
from pocketmapper.steps.parse import parse_job_entries
from pocketmapper.structure_aligner import StructureAligner

logger = logging.getLogger(__name__)


def superpose(
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
    pocket_comparison=None,
    alignment=None,
    pockets=None,
    aligned_structure_dir=None,
    align_struct_method=None,
    align_count=None,
    threads=None,
):
    """
    Superpose the top targets of each query onto it, one PDB per query.

    Parses the entries the job file names, and adds those whose structure is missing to
    `failed_entries_path`. The aligner is read off the alignment table: the seq aligner leaves its
    `u` column "-". Against a Foldseek database, the target structures are rebuilt out of it, which
    needs foldseek.

    Args:
        job_file (str or dict, optional): JSON job file of job key -> value, or the same already
            loaded, e.g. parse's settings. Any argument given overrides it. Must set query and target.
        results_dir (str, optional): The results directory the inputs default to. Required here or
            in `job_file`.
        work_dir (str, optional): Directory that entries and relative paths resolve against.
            Defaults to the working directory.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Defaults to <results_dir>/info.log.
        job_settings_path (str, optional): Where these settings are written, as a job file for later
            steps. Defaults to <results_dir>/superpose_settings.json. Not written when run inside search.
        failed_entries_path (str, optional): Defaults to <results_dir>/failed_entries.json.
        cache_dir (str, optional): Defaults to DEFAULT_CACHE_DIR.
        pdb_dir (str, optional): Defaults to <cache_dir>/pdb_structures.
        alphafold_dir (str, optional): Defaults to <cache_dir>/alphafold_structures.
        pocket_dir (str, optional): Defaults to <cache_dir>/pockets.
        foldseek_preprocessed_structure_dir (str, optional): Defaults to
            <cache_dir>/foldseek_preprocessed_structures.
        fsdb_dir (str, optional): Defaults to <cache_dir>/fsdb.
        pocket_comparison (str, optional): Defaults to the job file's pocket_comparison_path, else
            <results_dir>/pocket_comparison.tsv.
        alignment (str, optional): As `pocket_comparison`, from alignment_path.
        pockets (str, optional): As `pocket_comparison`, from pockets_path. Read only against a
            PDB-named Foldseek database, for its hits' pockets.
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
        PocketMapperError: If the job file cannot be read, query, target or results_dir is not given,
            the entries are rejected as `parse` rejects them, an input is missing, a setting is invalid or does not suit the alignment or target, or foldseek is
            needed but cannot run.
    """
    log_extra = {"stage": "Structural Alignment"}

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
            "aligned_structure_dir": aligned_structure_dir,
            "align_struct_method": align_struct_method,
            "align_count": align_count,
            "threads": threads,
        },
    )
    for key in ("query", "target", "results_dir"):
        require_setting(values, key)
    values = resolve_paths(values, "superpose")
    pocket_comparison = input_path(values, pocket_comparison, "pocket_comparison_path")
    alignment = input_path(values, alignment, "alignment_path")
    pockets = input_path(values, pockets, "pockets_path")
    with run_scope("superpose") as outermost, log_to_file(values["log_path"], values["verbosity"]):
        if outermost:
            dump_settings(values)
        require_file(alignment, "alignment")
        require_file(pocket_comparison, "pocket comparison")
        transforms = pd.read_csv(alignment, sep="\t", usecols=["u"], dtype=str)["u"]
        if transforms.empty:
            logger.warning(f"{alignment} holds no alignments; nothing to superpose", extra=log_extra)
            return
        aligner = "seq" if transforms.iloc[0] == "-" else "foldseek"
        align_struct_method = resolve_align_struct_method(values["align_struct_method"], aligner)
        sides = parse_job_entries(values, "superpose")
        if fsdb_record(sides["target"]) is not None:
            check_fsdb_align_struct_method(align_struct_method)
            require_foldseek("Rebuilding the target structures out of the Foldseek database needs foldseek")
        threads = resolve_threads(values["threads"])

        superpose_top_targets(
            sides,
            {"query": values["query"], "target": values["target"]},
            values["failed_entries_path"],
            pocket_comparison,
            alignment,
            pockets,
            values["aligned_structure_dir"],
            align_struct_method,
            values["align_count"],
            threads,
        )


def superpose_top_targets(
    sides,
    sources,
    failed_entries_path,
    pocket_comparison,
    alignment,
    pockets,
    aligned_structure_dir,
    align_struct_method,
    align_count,
    threads,
):
    """
    Superpose the top `align_count` targets of each query onto it, into `aligned_structure_dir`.

    Leaves out, adding to `failed_entries_path`, each record whose structure is missing
    (`structure_not_found`). Against a Foldseek database, the target structures are rebuilt out of it
    under `aligned_structure_dir`.

    Args:
        sides (dict): "query" and "target" -> that side's QTRecord dicts.
        sources (dict): "query" and "target" -> the input the side was parsed from, for the failure
            entries.
        failed_entries_path (str): The failed-entries file, appended to.
        pocket_comparison (str): The pocket comparison table.
        alignment (str): The alignment table; its transforms are read with the "foldseek" method.
        pockets (str): The pockets file; read only against a PDB-named Foldseek database, whose
            `chains` name the pockets on each hit.
        aligned_structure_dir (str): Where the superposed structures are written.
        align_struct_method (str): "pocket" or "foldseek", already resolved.
        align_count (int): Most targets to superpose onto each query; 0 writes nothing.
        threads (int): Thread count for rebuilding structures out of a Foldseek database.

    Returns:
        None

    Raises:
        PocketMapperError: If the pockets file cannot be read, or Foldseek fails rebuilding
            structures.
    """
    log_extra = {"stage": "Structural Alignment"}

    # A structure fetch never produced, or one removed since
    query = split_missing_structures(sides["query"], "superpose", sources["query"], failed_entries_path, log_extra)
    target = sides["target"]

    # With a Foldseek database, a hit's transform fits the database's own structure, so that is
    # what is superposed. Given no target records, the aligner reads target ids as entry names: only
    # a PDB database's hits have pockets of their own, each found on its hit by `chains`.
    database = fsdb_record(target)
    fsdb_path = database["struct_path"] if database is not None else None
    if database is None:
        target = split_missing_structures(target, "superpose", sources["target"], failed_entries_path, log_extra)
    else:
        hit_names = pd.read_csv(alignment, sep="\t", usecols=["target"], dtype=str)["target"].unique()
        target = []
        if fsdb_pocket_mode(hit_names) == "pisa":
            require_file(pockets, "pockets")
            _, chains = read_pockets_file(pockets)
            hit_names = set(hit_names)
            target = [
                {"pocket_id": pocket_id, "preprocess_name": name}
                for name, pocket_ids in chains.items()
                if name in hit_names
                for pocket_id in pocket_ids
            ]

    StructureAligner().align_structs(
        query_records=query,
        target_records=target,
        pocket_comparison=pocket_comparison,
        out_dir=aligned_structure_dir,
        method=align_struct_method,
        align_count=align_count,
        alignment=alignment,
        threads=threads,
        fsdb_path=fsdb_path,
    )
