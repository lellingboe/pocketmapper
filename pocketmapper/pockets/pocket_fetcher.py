"""
Dispatch of records to their pocket method's builder.

`POCKET_BUILDERS` is the whole table of pocket methods. Every builder takes `(records, pocket_dir)`,
plus any keyword options its caller passes for that method, and returns pocket_id -> Pocket, so the
fetcher needs to know nothing about any one method.

Also the pocket file formats: `dump_pockets`/`load_pockets` for a plain pocket_id -> Pocket mapping
(the per-method cache files), `write_pockets_file`/`read_pockets_file` for a versioned pockets file
that also maps each aligned chain to the pockets on it.
"""

import json
import logging
import os
from dataclasses import asdict

from pocketmapper.exceptions import PocketMapperError
from pocketmapper.pockets.pisa import pisa_pockets
from pocketmapper.pockets.pocket import Pocket
from pocketmapper.pockets.pocket import PocketResidue
from pocketmapper.pockets.structure import passthrough_pockets
from pocketmapper.pockets.structure import whole_chain_pockets
from pocketmapper.pockets.vdw import vdw_pockets

logger = logging.getLogger(__name__)

# The pockets file format `write_pockets_file` writes and `read_pockets_file` accepts
POCKETS_FILE_VERSION = 2

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
    with open(path, "w") as f:
        json.dump(serialise_pockets(pockets), f)


def load_pockets(path):
    """
    Read a pocket collection written by `dump_pockets`.

    Args:
        path (str): File to read.

    Returns:
        dict: pocket_id -> Pocket, or None where the file holds `null`. Each Pocket's residues are
            rebuilt as PocketResidues, keeping `res_auth_ids` in file order and every None field.
    """
    with open(path) as f:
        return deserialise_pockets(json.load(f))


def write_pockets_file(pockets, chains, path):
    """
    Write a pockets file: the pockets plus the map from aligned chain to the pockets on it.

    Args:
        pockets (dict): pocket_id -> Pocket, or None for an entry given and not built.
        chains (dict): preprocess_name -> the pocket_ids on that chain; empty for a chain with none.
        path (str): File to write.

    Returns:
        None: Writes a file.
    """
    with open(path, "w") as f:
        json.dump({"version": POCKETS_FILE_VERSION, "pockets": serialise_pockets(pockets), "chains": chains}, f)


def read_pockets_file(path):
    """
    Read a pockets file written by `write_pockets_file`.

    Args:
        path (str): File to read.

    Returns:
        tuple: (pockets, chains), as `write_pockets_file` takes them; pockets rebuilt as by
            `load_pockets`.

    Raises:
        PocketMapperError: If the file is not JSON, is not a pockets file, or is of another version.
    """
    try:
        with open(path) as f:
            serialised = json.load(f)
    except ValueError as e:
        msg = f"Could not read the pockets file {path}: {e}"
        logger.critical(msg)
        raise PocketMapperError(msg) from e
    version = serialised.get("version") if isinstance(serialised, dict) else None
    if version != POCKETS_FILE_VERSION:
        msg = f"{path} is not a version {POCKETS_FILE_VERSION} pockets file (version {version}); rerun pockets"
        logger.critical(msg)
        raise PocketMapperError(msg)
    try:
        return deserialise_pockets(serialised["pockets"]), dict(serialised["chains"])
    except (TypeError, AttributeError, KeyError, ValueError) as e:
        msg = f"Could not read the pockets file {path}: {e}"
        logger.critical(msg)
        raise PocketMapperError(msg) from e


def serialise_pockets(pockets):
    """
    Turn a pocket collection into JSON-ready dicts.

    Args:
        pockets (dict): pocket_id -> Pocket, or None.

    Returns:
        dict: pocket_id -> the Pocket's fields, or None.
    """
    return {pid: asdict(pocket) if pocket is not None else None for pid, pocket in pockets.items()}


def deserialise_pockets(serialised):
    """
    Rebuild a pocket collection from `serialise_pockets` output.

    Args:
        serialised (dict): pocket_id -> a Pocket's fields, or None.

    Returns:
        dict: pocket_id -> Pocket, or None. Residues are rebuilt as PocketResidues, keeping
            `res_auth_ids` in order and every None field.
    """
    pockets = {}
    for pid, fields in serialised.items():
        if fields is None:
            pockets[pid] = None
            continue
        residues = {res_id: PocketResidue(**residue) for res_id, residue in fields["residues"].items()}
        pockets[pid] = Pocket(**(fields | {"residues": residues}))
    return pockets


class PocketFetcher:
    """
    Builds the pockets a set of records names, dispatching each record to its pocket method's builder.
    """

    def fetch_pockets(self, records, pocket_dir, builder_options=None):
        """
        Build a Pocket for each record, by the builder of its `pocket_method`.

        Downloads any PISA data not already cached under `pocket_dir`, and overwrites
        `pocket_dir/<method>_pockets.json` for every method present, creating `pocket_dir` if needed. PISA records sharing a structure
        and chain are built once, under the first one's pocket_id.

        Args:
            records (list): QTRecord dicts; reads `pocket_id`, `pocket_method`, `struct_info`,
                `struct_path`, `chain_info` and `residue_info`. Query and target records may be mixed.
                Used as given: filter out any whose structure is unavailable before calling.
            pocket_dir (str): Pocket cache directory.
            builder_options (dict, optional): pocket_method -> keyword arguments for that method's
                builder, e.g. `{"pisa": {"pisa_source": "api"}}`. Defaults to None, which passes none.

        Returns:
            dict: pocket_id -> Pocket. A record whose pocket cannot be built is left out, or for the
                pisa and vdw methods may map to None.
        """
        log_extra = {"stage": "Getting Pockets"}
        logger.info("Starting pocket retrieval...", extra=log_extra)
        builder_options = builder_options or {}

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

            method_pockets = builder(method_records, pocket_dir, **builder_options.get(pocket_method, {}))
            os.makedirs(pocket_dir, exist_ok=True)
            dump_pockets(method_pockets, os.path.join(pocket_dir, f"{pocket_method}_pockets.json"))
            logger.debug(f"Extracted {label} pockets: {method_pockets}", extra=log_extra)
            pockets |= method_pockets

        logger.debug(f"Combined pockets: {pockets}", extra=log_extra)
        return pockets
