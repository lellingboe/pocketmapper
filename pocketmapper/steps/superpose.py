"""
Step 6: superpose the top targets of each query onto it, one PDB per query.
"""

import logging

import pandas as pd

from pocketmapper.lib import fsdb_pocket_mode
from pocketmapper.lib import log_to_file
from pocketmapper.pockets.pocket_fetcher import read_pockets_file
from pocketmapper.records import fsdb_record
from pocketmapper.records import read_records
from pocketmapper.records import require_file
from pocketmapper.settings import check_fsdb_align_struct_method
from pocketmapper.settings import layer_settings
from pocketmapper.settings import require_foldseek
from pocketmapper.settings import require_setting
from pocketmapper.settings import resolve_align_struct_method
from pocketmapper.settings import resolve_paths
from pocketmapper.settings import resolve_threads
from pocketmapper.structure_aligner import StructureAligner

logger = logging.getLogger(__name__)


def superpose(
    job_file=None,
    results_dir=None,
    verbosity=None,
    log_path=None,
    query_records=None,
    target_records=None,
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

    The aligner is read off the alignment table: the seq aligner leaves its `u` column "-". Against
    a Foldseek database, the target structures are rebuilt out of it, which needs foldseek.

    Args:
        job_file (str or dict, optional): JSON job file of Settings field name -> value, or the same
            already loaded. Any argument given overrides it.
        results_dir (str, optional): The results directory the inputs default to. Required here or
            in `job_file`.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Defaults to <results_dir>/info.log.
        query_records (str, optional): Defaults to the job file's query_records_path, else
            <results_dir>/query_records.json.
        target_records (str, optional): As `query_records`, from target_records_path.
        pocket_comparison (str, optional): As `query_records`, from pocket_comparison_path.
        alignment (str, optional): As `query_records`, from alignment_path.
        pockets (str, optional): As `query_records`, from pockets_path. Read only against a PDB-named
            Foldseek database, for its hits' pockets.
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
        PocketMapperError: If the job file cannot be read, no results_dir is given, an input is
            missing, a setting is invalid or does not suit the alignment or target, or foldseek is
            needed but cannot run.
    """
    log_extra = {"stage": "Structural Alignment"}

    values = layer_settings(
        job_file,
        {
            "results_dir": results_dir,
            "verbosity": verbosity,
            "log_path": log_path,
            "aligned_structure_dir": aligned_structure_dir,
            "align_struct_method": align_struct_method,
            "align_count": align_count,
            "threads": threads,
        },
    )
    require_setting(values, "results_dir")
    values = resolve_paths(values)
    query_records = query_records if query_records is not None else values["query_records_path"]
    target_records = target_records if target_records is not None else values["target_records_path"]
    pocket_comparison = pocket_comparison if pocket_comparison is not None else values["pocket_comparison_path"]
    alignment = alignment if alignment is not None else values["alignment_path"]
    pockets = pockets if pockets is not None else values["pockets_path"]
    with log_to_file(values["log_path"], values["verbosity"]):
        require_file(alignment, "alignment")
        require_file(pocket_comparison, "pocket comparison")
        transforms = pd.read_csv(alignment, sep="\t", usecols=["u"], dtype=str)["u"]
        if transforms.empty:
            logger.warning(f"{alignment} holds no alignments; nothing to superpose", extra=log_extra)
            return
        aligner = "seq" if transforms.iloc[0] == "-" else "foldseek"
        align_struct_method = resolve_align_struct_method(values["align_struct_method"], aligner)
        if fsdb_record(read_records(target_records)) is not None:
            check_fsdb_align_struct_method(align_struct_method)
            require_foldseek("Rebuilding the target structures out of the Foldseek database needs foldseek")
        threads = resolve_threads(values["threads"])

        superpose_top_targets(
            query_records,
            target_records,
            pocket_comparison,
            alignment,
            pockets,
            values["aligned_structure_dir"],
            align_struct_method,
            values["align_count"],
            threads,
        )


def superpose_top_targets(
    query_records,
    target_records,
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

    Against a Foldseek database, the target structures are rebuilt out of it under
    `aligned_structure_dir`.

    Args:
        query_records (str): The query records file.
        target_records (str): The target records file.
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
        PocketMapperError: If a records file or the pockets file cannot be read, or Foldseek fails
            rebuilding structures.
    """
    query = read_records(query_records)
    target = read_records(target_records)

    # With a Foldseek database, a hit's transform fits the database's own structure, so that is
    # what is superposed. Given no target records, the aligner reads target ids as entry names: only
    # a PDB database's hits have pockets of their own, each found on its hit by `chains`.
    database = fsdb_record(target)
    fsdb_path = database["struct_path"] if database is not None else None
    if database is not None:
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
