"""
Step 4: build the pocket of every record into one pockets file.

A record whose pocket cannot be built is dropped from its records file, so every record left has a
pocket for the comparison to read.
"""

import logging
import os

from pocketmapper.lib import log_to_file
from pocketmapper.lib import make_dir
from pocketmapper.pockets.pocket_fetcher import PocketFetcher
from pocketmapper.pockets.pocket_fetcher import dump_pockets
from pocketmapper.records import append_failed_entries
from pocketmapper.records import failed_entry
from pocketmapper.records import read_cache_manifest
from pocketmapper.records import read_records
from pocketmapper.records import write_records
from pocketmapper.settings import layer_settings
from pocketmapper.settings import require_setting
from pocketmapper.settings import resolve_paths
from pocketmapper.settings import resolve_pisa_source

logger = logging.getLogger(__name__)


def pockets(
    records=None,
    job_file=None,
    results_dir=None,
    verbosity=None,
    log_path=None,
    failed_entries_path=None,
    pockets_path=None,
    pisa_source=None,
):
    """
    Build the pocket of every record in some records files.

    Reads the pocket cache directory from `results_dir`'s cache manifest. Overwrites `pockets_path`
    with the pockets of every file named, and rewrites each file in place without the records whose
    pocket could not be built, adding those to `failed_entries_path`. Downloads any PISA data not
    already cached.

    Args:
        records (list, optional): The records files. The pockets file holds only the pockets of the
            files named. Defaults, when None or empty, to the query and target records files: the
            job file's query_records_path and target_records_path, else the standard files under
            <results_dir>.
        job_file (str or dict, optional): JSON job file of Settings field name -> value, or the same
            already loaded. Any argument given overrides it.
        results_dir (str, optional): The results directory `parse` wrote to. Required here or in
            `job_file`.
        verbosity (int, optional): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. Defaults to DEFAULT_VERBOSITY.
        log_path (str, optional): Defaults to <results_dir>/info.log.
        failed_entries_path (str, optional): Defaults to <results_dir>/failed_entries.json.
        pockets_path (str, optional): Defaults to <results_dir>/pockets.json.
        pisa_source (str, optional): "ftp" or "api". Defaults to DEFAULT_PISA_SOURCE.

    Returns:
        None

    Raises:
        PocketMapperError: If the job file cannot be read, no results_dir is given, the manifest or a
            records file is missing, or `pisa_source` is invalid.
    """
    values = layer_settings(
        job_file,
        {
            "results_dir": results_dir,
            "verbosity": verbosity,
            "log_path": log_path,
            "failed_entries_path": failed_entries_path,
            "pockets_path": pockets_path,
            "pisa_source": pisa_source,
        },
    )
    require_setting(values, "results_dir")
    values = resolve_paths(values)
    # The command line gives [] for none
    if not records:
        records = [values["query_records_path"], values["target_records_path"]]
    with log_to_file(values["log_path"], values["verbosity"]):
        cache_dirs = read_cache_manifest(values["results_dir"])
        pisa_source = resolve_pisa_source(values["pisa_source"])
        build_pockets(
            records,
            values["pockets_path"],
            values["failed_entries_path"],
            cache_dirs["pocket_dir"],
            pisa_source,
        )


def build_pockets(records_paths, pockets_path, failed_entries_path, pocket_dir, pisa_source):
    """
    Build a Pocket for every record across several records files, and write them to `pockets_path`.

    Rewrites each records file in place without the records whose pocket could not be built, and
    adds those to `failed_entries_path` as `pocket_not_built`. Foldseek-database records have no
    pocket and are kept. Downloads any PISA data not already cached under `pocket_dir/pisa/`, and
    overwrites the per-method pocket files under `pocket_dir`.

    Args:
        records_paths (list): The records files, e.g. query then target. Records are built in this
            order, which decides which of several identical pisa records names the shared pocket.
        pockets_path (str): Where the pockets are written, overwriting any there.
        failed_entries_path (str): The failed-entries file, appended to.
        pocket_dir (str): Pocket cache directory.
        pisa_source (str): Where PISA interfaces are fetched from: "ftp" or "api".

    Returns:
        None

    Raises:
        PocketMapperError: If a records file cannot be read.
    """
    log_extra = {"stage": "Getting Pockets"}

    files = {path: read_records(path) for path in records_paths}
    buildable = [record for records in files.values() for record in records if record["struct_type"] != "foldseek_db"]
    pockets = PocketFetcher().fetch_pockets(
        buildable,
        pocket_dir,
        builder_options={"pisa": {"pisa_source": pisa_source}},
    )

    # Checked on the merged pockets: a later method can overwrite an earlier one's pocket with None
    failures = []
    for path, records in files.items():
        kept = []
        for record in records:
            if record["struct_type"] == "foldseek_db" or pockets.get(record["pocket_id"]) is not None:
                kept.append(record)
            else:
                failures.append(failed_entry(record["pocket_id"], "pockets", "pocket_not_built", path, record))
        files[path] = kept
    if failures:
        logger.warning(
            f"No pocket for {', '.join(dict.fromkeys(entry['pocket_id'] for entry in failures))}; skipping them",
            extra=log_extra,
        )
    append_failed_entries(failed_entries_path, failures)

    make_dir(os.path.dirname(pockets_path), log_extra)
    dump_pockets({pocket_id: pocket for pocket_id, pocket in pockets.items() if pocket is not None}, pockets_path)
    logger.info(f"Pockets written to {pockets_path}", extra=log_extra)
    for path, records in files.items():
        write_records(records, path)
