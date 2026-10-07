"""
Step 1: parse the query and target inputs into records files.

No network. Writes the cache manifest, so later steps resolve every cache path the way this one did.
"""

import logging
import os

from pocketmapper.exceptions import PocketMapperError
from pocketmapper.qt_processor import QTProcessor
from pocketmapper.records import append_failed_entries
from pocketmapper.records import failed_entry
from pocketmapper.records import fsdb_record
from pocketmapper.records import start_failed_entries
from pocketmapper.records import write_cache_manifest
from pocketmapper.records import write_records

logger = logging.getLogger(__name__)


def parse_inputs(
    query,
    target,
    query_pocket_method,
    target_pocket_method,
    cache_dirs,
    results_dir,
    query_records_path,
    target_records_path,
    failed_entries_path,
):
    """
    Parse both sides' input into records, and write them with the cache manifest.

    Starts `failed_entries_path` afresh, then adds every entry that could not be parsed and every
    Foldseek-database query entry, as `invalid_entry`.

    Args:
        query (str): Query entry, or a file of one entry per line.
        target (str): Target entry, a file of them, or a Foldseek database.
        query_pocket_method (str): Pocket method to force on every query entry, or "auto".
        target_pocket_method (str): As `query_pocket_method`, for the target side.
        cache_dirs (dict): Each of `records.CACHE_MANIFEST_KEYS` -> its directory. Made absolute for
            the record paths and the manifest.
        results_dir (str): Where the manifest is written.
        query_records_path (str): Where the query records are written.
        target_records_path (str): Where the target records are written.
        failed_entries_path (str): The failed-entries file.

    Returns:
        None

    Raises:
        PocketMapperError: If a pocket method is unknown, either side has no valid entries, or a
            Foldseek-database target is not the only target entry.
    """
    log_extra = {"stage": "Determine Query/Target Types"}

    start_failed_entries(failed_entries_path)
    cache_dirs = {key: os.path.abspath(path) for key, path in cache_dirs.items()}
    qtprocessor = QTProcessor(
        pdb_dir=cache_dirs["pdb_dir"],
        alphafold_dir=cache_dirs["alphafold_dir"],
        fsdb_dir=cache_dirs["fsdb_dir"],
    )

    sides = {}
    failures = []
    for name, qt_input, pocket_method in (
        ("query", query, query_pocket_method),
        ("target", target, target_pocket_method),
    ):
        records, rejected = qtprocessor.process_qt_cmdline_input(
            qt_input=qt_input, name=name, pocket_method=pocket_method
        )
        failures += [
            failed_entry(entry, "parse", "invalid_entry", qt_input, detail=reason) for entry, reason in rejected
        ]
        sides[name] = records

    # A database can only be searched, not searched with
    for record in [record for record in sides["query"] if record["struct_type"] == "foldseek_db"]:
        reason = f"A Foldseek database cannot be a query entry: {record['pocket_id']}"
        logger.warning(f"{reason}; skipping this entry", extra=log_extra)
        failures.append(failed_entry(record["pocket_id"], "parse", "invalid_entry", query, record, reason))
        sides["query"].remove(record)
    append_failed_entries(failed_entries_path, failures)

    errors = []
    for name, records in sides.items():
        if not records:
            logger.critical(f"No valid {name} entries after processing", extra=log_extra)
            errors.append(f"no valid {name} entries")
    if errors:
        raise PocketMapperError("; ".join(errors))

    if fsdb_record(sides["target"]) is not None and len(sides["target"]) > 1:
        msg = (
            "A Foldseek database target must be the only target entry, but "
            f"{len(sides['target'])} target entries were given"
        )
        logger.critical(msg, extra=log_extra)
        raise PocketMapperError(msg)

    write_records(sides["query"], query_records_path)
    write_records(sides["target"], target_records_path)
    write_cache_manifest(results_dir, cache_dirs)
    logger.info(f"Parsed {len(sides['query'])} query and {len(sides['target'])} target entries", extra=log_extra)
