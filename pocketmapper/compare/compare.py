"""
Step 5: compare the pockets of every aligned query/target pair into a pocket comparison table.
"""

import json
import logging
import os
from importlib.resources import files

import pandas as pd

from pocketmapper.compare.pocket_comparison import compare_pockets
from pocketmapper.lib.entries import fsdb_record, open_entries, require_file
from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.files import jsonify_dict, make_dir
from pocketmapper.lib.foldseek import bundled_offset_table, fsdb_pocket_mode
from pocketmapper.lib.logs import log_to_file
from pocketmapper.lib.scopes import run_scope
from pocketmapper.lib.settings import layer_settings, resolve_step_paths
from pocketmapper.pockets.pocket_fetcher import read_pockets_file

logger = logging.getLogger(__name__)

# The packaged BLAST-format similarity matrix the comparison scores with. Anchored on the parent package:
# `data` has no __init__.py, and importlib.resources only resolves a namespace package from 3.12
BLOSUM_PATH = str(files("pocketmapper").joinpath("data", "blosum62.bla"))


def compare(
    job_file=None,
    entries=None,
    alignment=None,
    pockets=None,
    verbosity=None,
    log_path=None,
    pocket_comparison_path=None,
):
    """
    Compare the pockets of every aligned query/target pair into a pocket comparison table.

    Reads the entries table for whether the target is a Foldseek database. Writes unknown_ids.json and
    incorrect_mapping.json beside `pocket_comparison_path` when either has anything to report, deleting
    any left there by an earlier run.

    Args:
        job_file (str or dict, optional): JSON job file of job key -> value, or the same already
            loaded, e.g. parse's settings. Any argument given overrides it.
        entries (str, optional): The entries table. Defaults to the job file's entries_path.
        alignment (str, optional): Defaults to the job file's alignment_path.
        pockets (str, optional): Defaults to the job file's pockets_path.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Where the log is appended. Defaults to the job file's, else none.
        pocket_comparison_path (str, optional): Defaults to the job file's.

    Returns:
        None

    Raises:
        PocketMapperError: If the job file cannot be read, a path is not given, the entries table cannot
            be read or leaves a side empty, an input is missing or unreadable, or the alignment names a
            chain the pockets file does not.
    """
    values = layer_settings(
        job_file,
        {
            "entries_path": entries,
            "alignment_path": alignment,
            "pockets_path": pockets,
            "verbosity": verbosity,
            "log_path": log_path,
            "pocket_comparison_path": pocket_comparison_path,
        },
    )
    values = resolve_step_paths(
        values, ("entries_path", "alignment_path", "pockets_path", "pocket_comparison_path"), ("log_path",)
    )
    with run_scope("compare"), log_to_file(values["log_path"], values["verbosity"]):
        sides = open_entries(values["entries_path"], "compare")
        compare_aligned_pockets(
            sides["target"],
            values["alignment_path"],
            values["pockets_path"],
            values["pocket_comparison_path"],
        )


def compare_aligned_pockets(target_records, alignment, pockets, pocket_comparison_path):
    """
    Compare the pockets of every aligned query/target pair and write `pocket_comparison_path`.

    Writes unknown_ids.json and incorrect_mapping.json beside `pocket_comparison_path` when either has
    anything to report, and deletes any left there by an earlier run either way.

    Args:
        target_records (list): The target QTRecord dicts; read for whether the target is a Foldseek
            database.
        alignment (str): The alignment table.
        pockets (str): The pockets file. Its `chains` must name every chain in the alignment, but
            the hits of a Foldseek database whose pockets are synthesised.
        pocket_comparison_path (str): Where the comparison table is written.

    Returns:
        None

    Raises:
        PocketMapperError: If an input is missing or unreadable, the alignment names a chain the
            pockets file does not, the bundled database's offset table is missing from the
            installation, or an output directory cannot be created.
    """
    log_extra = {"stage": "Comparing Pockets Based on Alignment"}

    # Both are written only when non-empty, so an old copy would otherwise outlive this run
    report_dir = os.path.dirname(pocket_comparison_path)
    report_paths = {name: os.path.join(report_dir, f"{name}.json") for name in ("unknown_ids", "incorrect_mapping")}
    for path in report_paths.values():
        if os.path.isfile(path):
            os.remove(path)

    require_file(pockets, "pockets")
    pocket_dict, chain_pockets = read_pockets_file(pockets)
    require_file(alignment, "alignment")

    logger.info("Reading alignment results...", extra=log_extra)
    alignment_df = pd.read_csv(alignment, sep="\t", engine="c")
    logger.info(f"{len(alignment_df)} alignment pairs to compare", extra=log_extra)
    logger.debug(f"Alignment pairs: \n{alignment_df.head()}", extra=log_extra)

    logger.debug(f"Preprocessed name to pocket ID mapping: {chain_pockets}", extra=log_extra)

    # A PDB Foldseek database's hits have PISA pockets; any other database's hits have none, so
    # their pockets must be synthesised
    database = fsdb_record(target_records)
    synthesise = database is not None and fsdb_pocket_mode(alignment_df["target"].unique()) == "whole_chain"
    offset_table_path = None
    if synthesise:
        offset_table_path = bundled_offset_table(database["struct_path"])
        if offset_table_path is None:
            logger.info(
                "Foldseek database ships no offset table; target residue ids will be positions "
                "within each database entry rather than UniProt coordinates",
                extra=log_extra,
            )
        elif not os.path.isfile(offset_table_path):
            # Packaged with the database, so absence means a broken install. No fallback: the
            # bundled database's residue ids are documented as UniProt coordinates.
            msg = f"Bundled human-domains offset table is missing from the installation: {offset_table_path}"
            logger.critical(msg, extra=log_extra)
            raise PocketMapperError(msg)
        else:
            logger.info(f"Mapping target residue ids to UniProt coordinates using {offset_table_path}", extra=log_extra)

    check_coverage(alignment_df, chain_pockets, synthesise, pockets)

    pockets_df, unknown_alias, incorrect_mapping = compare_pockets(
        alignment_df,
        pocket_dict,
        preproc_to_ids=chain_pockets,
        blosum_path=BLOSUM_PATH,
        synthesise_target_pockets=synthesise,
        offset_table_path=offset_table_path,
    )

    # Logging cases where a residue was given a single char name unfamiliar to pocketmapper
    if len(unknown_alias) > 0:
        make_dir(report_dir, log_extra)
        logger.warning(f"Unknown Foldseek Alias, see {report_paths['unknown_ids']}", extra=log_extra)
        with open(report_paths["unknown_ids"], "w") as f:
            json.dump(jsonify_dict(dict(unknown_alias)), f)

    # logging cases where foldseek mapping had low sequence identity to the parsed structure
    if len(incorrect_mapping) > 0:
        make_dir(report_dir, log_extra)
        logger.warning(
            f"Foldseek mapping with low sequence identity to parsed structure, see {report_paths['incorrect_mapping']}",
            extra=log_extra,
        )
        with open(report_paths["incorrect_mapping"], "w") as f:
            json.dump(jsonify_dict(dict(incorrect_mapping)), f)

    make_dir(os.path.dirname(pocket_comparison_path), log_extra)
    pockets_df.to_csv(pocket_comparison_path, index=False, sep="\t")
    logger.info(f"Pocket comparison results saved to {pocket_comparison_path}", extra=log_extra)


def check_coverage(alignment_df, chain_pockets, synthesise, pockets):
    """
    Check that the pockets file has seen every chain the alignment names.

    A chain the pockets file lists with no built pocket is fine: its rows are skipped. One it does not
    list at all was never given to pockets, and would otherwise silently give no rows.

    Args:
        alignment_df (pandas.DataFrame): The alignment table.
        chain_pockets (dict): The pockets file's `chains`.
        synthesise (bool): Whether target pockets are synthesised, which exempts the target names.
        pockets (str): The pockets file, for the message.

    Returns:
        None

    Raises:
        PocketMapperError: If an alignment name is not a key of `chain_pockets`.
    """
    names = alignment_df["query"].tolist() + ([] if synthesise else alignment_df["target"].tolist())
    unseen = [name for name in dict.fromkeys(names) if name not in chain_pockets]
    if unseen:
        shown = ", ".join(unseen[:5]) + (f" and {len(unseen) - 5} more" if len(unseen) > 5 else "")
        msg = (
            f"{pockets} does not cover the alignment: no entry for {shown}. Run pockets on the entries "
            "that were aligned, unchanged since align"
        )
        logger.critical(msg, extra={"stage": "Comparing Pockets Based on Alignment"})
        raise PocketMapperError(msg)
