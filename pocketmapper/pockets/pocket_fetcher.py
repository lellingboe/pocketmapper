"""
Dispatch of records to their pocket method's builder.

`POCKET_BUILDERS` is the whole table of pocket methods. Every builder takes `(records, pocket_dir)`
and returns pocket_id -> Pocket, so the fetcher needs to know nothing about any one method.
"""

import json
import logging
import os
from dataclasses import asdict

from pocketmapper.pockets.pisa import pisa_pockets
from pocketmapper.pockets.structure import passthrough_pockets
from pocketmapper.pockets.structure import whole_chain_pockets
from pocketmapper.pockets.vdw import vdw_pockets

logger = logging.getLogger(__name__)

# pocket_method -> (log label, dedup keys, builder), in merge order.
# Only pisa dedups: its pocket is fully determined by structure and chain. The others must not --
# two passthrough pockets on one chain differ only in `residue_info`.
POCKET_BUILDERS = {
    "pisa": ("PISA", ("struct_info", "chain_info"), pisa_pockets),
    "passthrough": ("passthrough", None, passthrough_pockets),
    "vdw": ("VDW", None, vdw_pockets),
    "whole_chain": ("whole chain", None, whole_chain_pockets),
}


def dump_pockets(pockets, path):
    """
    Write a pocket collection to `path` as compact JSON.

    Args:
        pockets (dict): pocket_id -> Pocket. A None value is written as `null`.
        path (str): File to write.

    Returns:
        None: Writes a file.
    """
    serialisable = {pid: asdict(pocket) if pocket is not None else None for pid, pocket in pockets.items()}
    with open(path, "w") as f:
        json.dump(serialisable, f)


class PocketFetcher:
    """
    Builds the pockets a set of records names, dispatching each record to its pocket method's builder.
    """

    def fetch_pockets(self, records, pocket_dir):
        """
        Build a Pocket for each record, by the builder of its `pocket_method`.

        Downloads any PISA data not already cached under `pocket_dir`, and overwrites
        `pocket_dir/<method>_pockets.json` for every method present. PISA records sharing a structure
        and chain are built once, under the first one's pocket_id.

        Args:
            records (list): QTRecord dicts; reads `pocket_id`, `pocket_method`, `struct_info`,
                `struct_path`, `chain_info` and `residue_info`. Query and target records may be mixed.
                Used as given: filter out any whose structure is unavailable before calling.
            pocket_dir (str): Pocket cache directory.

        Returns:
            dict: pocket_id -> Pocket. A record whose pocket cannot be built is left out, or for the
                pisa and vdw methods may map to None.
        """
        log_extra = {"stage": "Getting Pockets"}
        logger.info("Starting pocket retrieval...", extra=log_extra)

        pockets = {}
        for pocket_method, (label, dedup_keys, builder) in POCKET_BUILDERS.items():
            method_log_extra = {"stage": f"Retrieving {label} Pockets"}
            logger.info(f"Checking for {label} pockets...", extra=method_log_extra)

            method_records = [record for record in records if record["pocket_method"] == pocket_method]
            if dedup_keys is not None:
                seen = set()
                unique_records = []
                for record in method_records:
                    key = tuple(record[k] for k in dedup_keys)
                    if key not in seen:
                        seen.add(key)
                        unique_records.append(record)
                method_records = unique_records

            if not method_records:
                logger.info(f"No {label} pockets to retrieve", extra=method_log_extra)
                continue
            logger.info(f"{len(method_records)} {label} pockets to retrieve", extra=method_log_extra)

            method_pockets = builder(method_records, pocket_dir)
            dump_pockets(method_pockets, os.path.join(pocket_dir, f"{pocket_method}_pockets.json"))
            logger.debug(f"Extracted {label} pockets: {method_pockets}", extra=log_extra)
            pockets |= method_pockets

        logger.debug(f"Combined pockets: {pockets}", extra=log_extra)
        return pockets
