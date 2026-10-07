"""
Step 4: build the pocket of every record into one pockets file.

A record whose pocket cannot be built is dropped from its records file, so every record left has a
pocket for the comparison to read.
"""

import logging
import os

from pocketmapper.lib import make_dir
from pocketmapper.pockets.pocket_fetcher import PocketFetcher
from pocketmapper.pockets.pocket_fetcher import dump_pockets
from pocketmapper.records import append_failed_entries
from pocketmapper.records import failed_entry
from pocketmapper.records import read_records
from pocketmapper.records import write_records

logger = logging.getLogger(__name__)


def build_pockets(records_paths, pockets_path, failed_entries_path, pocket_dir, pisa_source, download_pisa=True):
    """
    Build a Pocket for every record across several records files, and write them to `pockets_path`.

    Rewrites each records file in place without the records whose pocket could not be built, and
    adds those to `failed_entries_path` as `pocket_not_built`. Foldseek-database records have no
    pocket and are kept. Also overwrites the per-method pocket files under `pocket_dir`.

    Args:
        records_paths (list): The records files, e.g. query then target. Records are built in this
            order, which decides which of several identical pisa records names the shared pocket.
        pockets_path (str): Where the pockets are written, overwriting any there.
        failed_entries_path (str): The failed-entries file, appended to.
        pocket_dir (str): Pocket cache directory.
        pisa_source (str): Where PISA interfaces are fetched from: "ftp" or "api".
        download_pisa (bool, optional): False reads only the PISA data already cached. Defaults to True.

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
        builder_options={"pisa": {"pisa_source": pisa_source, "download": download_pisa}},
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
