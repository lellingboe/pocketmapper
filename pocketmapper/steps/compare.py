"""
Step 5: compare the pockets of every aligned query/target pair into a pocket comparison table.
"""

import json
import logging
import os

import pandas as pd

from pocketmapper.exceptions import PocketMapperError
from pocketmapper.foldseek import bundled_offset_table
from pocketmapper.lib import fsdb_pocket_mode
from pocketmapper.lib import jsonify_dict
from pocketmapper.lib import log_to_file
from pocketmapper.lib import make_dir
from pocketmapper.lib import run_scope
from pocketmapper.pocket_comparison import compare_pockets
from pocketmapper.pockets.pocket_fetcher import read_pockets_file
from pocketmapper.records import fsdb_record
from pocketmapper.records import require_file
from pocketmapper.settings import dump_settings
from pocketmapper.settings import input_path
from pocketmapper.settings import layer_settings
from pocketmapper.settings import require_setting
from pocketmapper.settings import resolve_paths
from pocketmapper.steps.parse import parse_job_entries

logger = logging.getLogger(__name__)

# The packaged BLAST-format similarity matrix the comparison scores with
BLOSUM_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "blosum62.bla")


def compare(
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
    alignment=None,
    pockets=None,
    pocket_comparison_path=None,
):
    """
    Compare the pockets of every aligned query/target pair into a pocket comparison table.

    Parses the entries the job file names, for whether the target is a Foldseek database. Writes
    unknown_ids.json and incorrect_mapping.json beside `pocket_comparison_path` when either has
    anything to report, deleting any left there by an earlier run.

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
            steps. Defaults to <results_dir>/compare_settings.json. Not written when run inside search.
        failed_entries_path (str, optional): Where entries rejected on parsing are listed. Defaults
            to <results_dir>/failed_entries.json.
        cache_dir (str, optional): Defaults to DEFAULT_CACHE_DIR.
        pdb_dir (str, optional): Defaults to <cache_dir>/pdb_structures.
        alphafold_dir (str, optional): Defaults to <cache_dir>/alphafold_structures.
        pocket_dir (str, optional): Defaults to <cache_dir>/pockets.
        foldseek_preprocessed_structure_dir (str, optional): Defaults to
            <cache_dir>/foldseek_preprocessed_structures.
        fsdb_dir (str, optional): Defaults to <cache_dir>/fsdb.
        alignment (str, optional): Defaults to the job file's alignment_path, else
            <results_dir>/alignment.tsv.
        pockets (str, optional): As `alignment`, from pockets_path.
        pocket_comparison_path (str, optional): Defaults to <results_dir>/pocket_comparison.tsv.

    Returns:
        None

    Raises:
        PocketMapperError: If the job file cannot be read, query, target or results_dir is not given,
            the entries are rejected as `parse` rejects them, an input is missing or unreadable, or
            the alignment names a chain the pockets file does not.
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
            "pocket_comparison_path": pocket_comparison_path,
        },
    )
    for key in ("query", "target", "results_dir"):
        require_setting(values, key)
    values = resolve_paths(values, "compare")
    with run_scope("compare") as outermost, log_to_file(values["log_path"], values["verbosity"]):
        if outermost:
            dump_settings(values)
        sides = parse_job_entries(values, "compare")
        compare_aligned_pockets(
            sides["target"],
            input_path(values, alignment, "alignment_path"),
            input_path(values, pockets, "pockets_path"),
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
