"""
Step 4: build the pocket of every query and target entry into one pockets file.

Every entry no earlier step dropped is listed, `null` where its pocket could not be built. Against a
PDB-named Foldseek database, each hit's PISA interfaces are expanded into pockets of their own here,
keyed in the pockets file by hit name.
"""

import logging
import os

import pandas as pd

from pocketmapper.fetch.fetch import fetch_missing_entries
from pocketmapper.lib.entries import (
    failure,
    fsdb_record,
    open_entries,
    preproc_to_ids,
    report_failures,
    require_file,
    split_missing_structures,
)
from pocketmapper.lib.files import make_dir
from pocketmapper.lib.foldseek import fsdb_pocket_mode
from pocketmapper.lib.logs import log_to_file
from pocketmapper.lib.scopes import run_scope
from pocketmapper.lib.settings import (
    layer_settings,
    resolve_fetch_missing,
    resolve_offline,
    resolve_pisa_bsa_threshold,
    resolve_pisa_source,
    resolve_step_paths,
)
from pocketmapper.pockets.pocket_manager import PocketManager, write_pockets_file, write_pockets_table

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
    local_root=None,
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
        pocket_dir (str, optional): Parent of the PISA cache, `pocket_dir/pisa/`. Defaults to the job file's.
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
        local_root (str, optional): Directory that relative local structure files and Foldseek databases
            resolve against. Defaults to the job file's, else the current directory.

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
            "local_root": local_root,
        },
    )
    values = resolve_step_paths(
        values,
        ("entries_path", "alignment_path", "pdb_dir", "alphafold_dir", "pocket_dir", "fsdb_dir", "pockets_path"),
        ("log_path", "pockets_tsv_path", "local_root"),
    )
    with run_scope("pockets"), log_to_file(values["log_path"], values["verbosity"]):
        pisa_source = resolve_pisa_source(values["pisa_source"])
        pisa_bsa_threshold = resolve_pisa_bsa_threshold(values["pisa_bsa_threshold"])
        fetch_missing = resolve_fetch_missing(values["fetch_missing"])
        offline = resolve_offline(values["offline"])
        sides = open_entries(values["entries_path"], "pockets")
        manager = PocketManager(
            values["pocket_dir"],
            pdb_dir=values["pdb_dir"],
            alphafold_dir=values["alphafold_dir"],
            fsdb_dir=values["fsdb_dir"],
            pisa_source=pisa_source,
            bsa_threshold=pisa_bsa_threshold,
            download=not offline,
        )
        if fetch_missing and not offline:
            fetch_missing_entries(sides, manager.structures)
        build_pockets(
            sides,
            values["entries_path"],
            values["alignment_path"],
            values["pockets_path"],
            manager,
            values["pockets_tsv_path"],
        )


def build_pockets(sides, entries_path, alignment, pockets_path, manager, pockets_tsv_path=None):
    """
    Build a Pocket for every record of every side, and write them to `pockets_path`.

    Marks the records whose structure is missing (`structure_not_found`) or whose pocket could not be
    built (`pocket_not_built`) in the entries table. Foldseek-database records have no pocket. Against a
    PDB-named Foldseek database, also builds a pisa pocket per PISA interface of each hit (see
    `PocketManager.expand_pdb_hits`); a hit is no entry, so one that fails is only logged. The pockets
    file lists every record and hit pocket, `null` where it was not built, and maps each
    `preprocess_name` and hit name to its pocket_ids. Downloads any PISA data, and the hits' structures,
    not already cached, unless `manager`'s downloads are off.

    Args:
        sides (dict): "query" and "target" -> that side's QTRecord dicts. Built in order, then the hits,
            which decides which of several identical pisa records names the shared pocket.
        entries_path (str or None): The entries table, marked, or None for none.
        alignment (str): The alignment table. Read only when the target is a Foldseek database.
        pockets_path (str): Where the pockets are written, overwriting any there.
        manager (PocketManager): Builds the pockets and expands the hits, with its PISA source,
            `bsa_threshold` (every pisa pocket, the Foldseek hits' included) and download policy.
        pockets_tsv_path (str, optional): Where a table of every pocket is written: pocket_id, chain,
            res_auth_ids and method. Defaults to None, for none.

    Returns:
        None

    Raises:
        PocketMapperError: If a Foldseek-database target has no alignment file, a PDB-named database's
            hits have interfaces but none of their structures could be fetched, `manager` lacks a
            directory the hits need, or an output directory cannot be created.
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
            hits, hit_chains, failures = manager.expand_pdb_hits(hit_names)
            report_failures(None, "pockets", failures, log_extra, "Missing structures for foldseek hit(s)")
            buildable += [(None, record) for record in hits]
            chains |= hit_chains

    built = manager.build([record for _, record in buildable])

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
        make_dir(os.path.dirname(pockets_tsv_path), log_extra)
        write_pockets_table([record for _, record in given] + hits, pockets, pockets_tsv_path)
