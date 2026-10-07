"""
Step 5: compare the pockets of every aligned query/target pair into a pocket comparison table.
"""

import json
import logging
import os

import pandas as pd

from pocketmapper.exceptions import PocketMapperError
from pocketmapper.foldseek import bundled_offset_table
from pocketmapper.lib import jsonify_dict
from pocketmapper.lib import log_to_file
from pocketmapper.lib import make_dir
from pocketmapper.pocket_comparison import compare_pockets
from pocketmapper.pockets.pocket_fetcher import load_pockets
from pocketmapper.records import fsdb_record
from pocketmapper.records import preproc_to_ids
from pocketmapper.records import read_records
from pocketmapper.records import require_file
from pocketmapper.records import synthesise_target_pockets
from pocketmapper.settings import layer_settings
from pocketmapper.settings import require_setting
from pocketmapper.settings import resolve_paths

logger = logging.getLogger(__name__)

# The packaged BLAST-format similarity matrix the comparison scores with
BLOSUM_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "blosum62.bla")


def compare(
    job_file=None,
    results_dir=None,
    verbosity=None,
    log_path=None,
    query_records=None,
    target_records=None,
    alignment=None,
    pockets=None,
    pocket_comparison_path=None,
):
    """
    Compare the pockets of every aligned query/target pair into a pocket comparison table.

    Writes unknown_ids.json and incorrect_mapping.json beside `pocket_comparison_path` when either has
    anything to report, deleting any left there by an earlier run.

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
        alignment (str, optional): As `query_records`, from alignment_path.
        pockets (str, optional): As `query_records`, from pockets_path.
        pocket_comparison_path (str, optional): Defaults to <results_dir>/pocket_comparison.tsv.

    Returns:
        None

    Raises:
        PocketMapperError: If the job file cannot be read, no results_dir is given, an input is
            missing or unreadable, a record has no pocket, or a Foldseek-database target has not been
            through align.
    """
    values = layer_settings(
        job_file,
        {
            "results_dir": results_dir,
            "verbosity": verbosity,
            "log_path": log_path,
            "pocket_comparison_path": pocket_comparison_path,
        },
    )
    require_setting(values, "results_dir")
    values = resolve_paths(values)
    with log_to_file(values["log_path"], values["verbosity"]):
        compare_aligned_pockets(
            query_records if query_records is not None else values["query_records_path"],
            target_records if target_records is not None else values["target_records_path"],
            alignment if alignment is not None else values["alignment_path"],
            pockets if pockets is not None else values["pockets_path"],
            values["pocket_comparison_path"],
        )


def compare_aligned_pockets(query_records, target_records, alignment, pockets, pocket_comparison_path):
    """
    Compare the pockets of every aligned query/target pair and write `pocket_comparison_path`.

    Writes unknown_ids.json and incorrect_mapping.json beside `pocket_comparison_path` when either has
    anything to report, and deletes any left there by an earlier run either way.

    Args:
        query_records (str): The query records file.
        target_records (str): The target records file.
        alignment (str): The alignment table.
        pockets (str): The pockets file. Must hold a pocket for every record but a Foldseek database.
        pocket_comparison_path (str): Where the comparison table is written.

    Returns:
        None

    Raises:
        PocketMapperError: If an input is missing or unreadable, a record has no pocket, a
            Foldseek-database target has not been through align, the bundled database's offset table
            is missing from the installation, or an output directory cannot be created.
    """
    log_extra = {"stage": "Comparing Pockets Based on Alignment"}

    # Both are written only when non-empty, so an old copy would otherwise outlive this run
    report_dir = os.path.dirname(pocket_comparison_path)
    report_paths = {name: os.path.join(report_dir, f"{name}.json") for name in ("unknown_ids", "incorrect_mapping")}
    for path in report_paths.values():
        if os.path.isfile(path):
            os.remove(path)

    query = read_records(query_records)
    target = read_records(target_records)
    require_file(pockets, "pockets")
    try:
        pocket_dict = load_pockets(pockets)
    except (ValueError, TypeError, AttributeError, KeyError) as e:
        # Not JSON, or JSON that is not a pockets file
        msg = f"Could not read the pockets file {pockets}: {e}"
        logger.critical(msg, extra=log_extra)
        raise PocketMapperError(msg) from e
    require_file(alignment, "alignment")

    # A record with no pocket would silently give no rows
    unbuilt = [
        record["pocket_id"]
        for record in query + target
        if record["struct_type"] != "foldseek_db" and pocket_dict.get(record["pocket_id"]) is None
    ]
    if unbuilt:
        msg = (
            f"{pockets} holds no pocket for {', '.join(dict.fromkeys(unbuilt))}; build the pockets of both "
            "records files first"
        )
        logger.critical(msg, extra=log_extra)
        raise PocketMapperError(msg)

    logger.info("Reading alignment results...", extra=log_extra)
    alignment_df = pd.read_csv(alignment, sep="\t", engine="c")
    logger.info(f"{len(alignment_df)} alignment pairs to compare", extra=log_extra)
    logger.debug(f"Alignment pairs: \n{alignment_df.head()}", extra=log_extra)

    chain_pockets = preproc_to_ids(query + target)
    logger.debug(f"Preprocessed name to pocket ID mapping: {chain_pockets}", extra=log_extra)

    # A PDB Foldseek database's hits have PISA pockets; any other database has no target records,
    # so its pockets must be synthesised
    synthesise = synthesise_target_pockets(target)
    offset_table_path = None
    if synthesise:
        offset_table_path = bundled_offset_table(fsdb_record(target)["struct_path"])
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
