"""
Step 6: superpose the top targets of each query onto it, one PDB per query.
"""

import logging

import pandas as pd

from pocketmapper.fetch.fetch import fetch_missing_entries
from pocketmapper.lib.entries import fsdb_record, open_entries, require_file, split_missing_structures
from pocketmapper.lib.foldseek import fsdb_pocket_mode
from pocketmapper.lib.logs import log_to_file
from pocketmapper.lib.scopes import run_scope
from pocketmapper.lib.settings import (
    check_fsdb_align_struct_method,
    layer_settings,
    require_foldseek,
    resolve_align_struct_method,
    resolve_fetch_missing,
    resolve_offline,
    resolve_step_paths,
    resolve_threads,
)
from pocketmapper.lib.structure_manager import StructureManager
from pocketmapper.pockets.pocket_manager import read_pockets_file
from pocketmapper.superpose.structure_aligner import StructureAligner

logger = logging.getLogger(__name__)


def superpose(
    job_file=None,
    entries=None,
    pocket_comparison=None,
    alignment=None,
    pockets=None,
    verbosity=None,
    log_path=None,
    aligned_structure_dir=None,
    align_struct_method=None,
    align_count=None,
    threads=None,
    fetch_missing=None,
    offline=None,
    pdb_dir=None,
    alphafold_dir=None,
    fsdb_dir=None,
    local_root=None,
):
    """
    Superpose the top targets of each query onto it, one PDB per query.

    Marks the entries whose structure is missing in the entries table. The aligner is read off the
    alignment table: the seq aligner leaves its `u` column "-". Against a Foldseek database, the target
    structures are rebuilt out of it, which needs foldseek.

    Args:
        job_file (str or dict, optional): JSON job file of job key -> value, or the same already
            loaded, e.g. parse's settings. Any argument given overrides it.
        entries (str, optional): The entries table. Defaults to the job file's entries_path.
        pocket_comparison (str, optional): Defaults to the job file's pocket_comparison_path.
        alignment (str, optional): Defaults to the job file's alignment_path.
        pockets (str, optional): Defaults to the job file's pockets_path. Read only against a
            PDB-named Foldseek database, for its hits' pockets.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Where the log is appended. Defaults to the job file's, else none.
        aligned_structure_dir (str, optional): Defaults to the job file's.
        align_struct_method (str, optional): "pocket", "foldseek" or "auto", which picks "foldseek"
            for a Foldseek alignment and "pocket" for a seq one. Defaults to
            DEFAULT_ALIGN_STRUCT_METHOD.
        align_count (int, optional): Most targets to superpose onto each query; 0 disables.
            Defaults to DEFAULT_ALIGN_COUNT.
        threads (int, optional): Defaults to one per available core.
        fetch_missing (int, optional): 1 downloads an entry structure missing from the cache, 0 skips
            the entry as `structure_not_found`. Defaults to DEFAULT_FETCH_MISSING.
        offline (int, optional): 1 downloads nothing, whatever `fetch_missing` says, so an uncached
            structure is skipped as with `fetch_missing` 0. Defaults to DEFAULT_OFFLINE.
        pdb_dir (str, optional): Cache of PDB structures. Defaults to the job file's.
        alphafold_dir (str, optional): Cache of AlphaFold structures. Defaults to the job file's.
        fsdb_dir (str, optional): Cache of Foldseek databases. Defaults to the job file's.
        local_root (str, optional): Directory that relative local structure files and Foldseek databases
            resolve against. Defaults to the job file's, else the current directory.

    Returns:
        None

    Raises:
        PocketMapperError: If the job file cannot be read, a path is not given, the entries table cannot
            be read or leaves a side empty, an input is missing, a setting is invalid or does not suit
            the alignment or target, or foldseek is needed but cannot run.
    """
    log_extra = {"stage": "Structural Alignment"}

    values = layer_settings(
        job_file,
        {
            "entries_path": entries,
            "pocket_comparison_path": pocket_comparison,
            "alignment_path": alignment,
            "pockets_path": pockets,
            "verbosity": verbosity,
            "log_path": log_path,
            "aligned_structure_dir": aligned_structure_dir,
            "align_struct_method": align_struct_method,
            "align_count": align_count,
            "threads": threads,
            "fetch_missing": fetch_missing,
            "offline": offline,
            "pdb_dir": pdb_dir,
            "alphafold_dir": alphafold_dir,
            "fsdb_dir": fsdb_dir,
            "local_root": local_root,
        },
    )
    values = resolve_step_paths(
        values,
        (
            "entries_path",
            "pocket_comparison_path",
            "alignment_path",
            "pockets_path",
            "aligned_structure_dir",
            "pdb_dir",
            "alphafold_dir",
            "fsdb_dir",
        ),
        ("log_path", "local_root"),
    )
    alignment = values["alignment_path"]
    with run_scope("superpose"), log_to_file(values["log_path"], values["verbosity"]):
        require_file(alignment, "alignment")
        require_file(values["pocket_comparison_path"], "pocket comparison")
        transforms = pd.read_csv(alignment, sep="\t", usecols=["u"], dtype=str)["u"]
        if transforms.empty:
            logger.warning(f"{alignment} holds no alignments; nothing to superpose", extra=log_extra)
            return
        aligner = "seq" if transforms.iloc[0] == "-" else "foldseek"
        align_struct_method = resolve_align_struct_method(values["align_struct_method"], aligner)
        fetch_missing = resolve_fetch_missing(values["fetch_missing"])
        offline = resolve_offline(values["offline"])
        sides = open_entries(values["entries_path"], "superpose")
        if fetch_missing and not offline:
            fetch_missing_entries(sides, StructureManager())
        if fsdb_record(sides["target"]) is not None:
            check_fsdb_align_struct_method(align_struct_method)
            require_foldseek("Rebuilding the target structures out of the Foldseek database needs foldseek")
        threads = resolve_threads(values["threads"])

        superpose_top_targets(
            sides,
            values["entries_path"],
            values["pocket_comparison_path"],
            alignment,
            values["pockets_path"],
            values["aligned_structure_dir"],
            align_struct_method,
            values["align_count"],
            threads,
        )


def superpose_top_targets(
    sides,
    entries_path,
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

    Leaves out, marking it in the entries table, each record whose structure is missing
    (`structure_not_found`). Against a Foldseek database, the target structures are rebuilt out of it
    under `aligned_structure_dir`.

    Args:
        sides (dict): "query" and "target" -> that side's QTRecord dicts.
        entries_path (str or None): The entries table, marked, or None for none.
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
    query = split_missing_structures(sides["query"], "query", "superpose", entries_path, log_extra)
    target = sides["target"]

    # With a Foldseek database, a hit's transform fits the database's own structure, so that is
    # what is superposed. Given no target records, the aligner reads target ids as entry names: only
    # a PDB database's hits have pockets of their own, each found on its hit by `chains`.
    database = fsdb_record(target)
    fsdb_path = database["struct_path"] if database is not None else None
    if database is None:
        target = split_missing_structures(target, "target", "superpose", entries_path, log_extra)
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
