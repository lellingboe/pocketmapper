"""
The pocket files: the versioned pockets file and a table of every pocket.

`write_pockets_file`/`read_pockets_file` handle the versioned pockets file, which also maps each aligned
chain to the pockets on it, and `write_pockets_table` writes a table of every pocket. This module imports
no builder or manager, so readers of `pockets.json` stay light.
"""

import json
import logging
import os
from dataclasses import asdict

import pandas as pd

from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.residues import split_chain_info
from pocketmapper.pockets.pocket import Pocket, PocketResidue

logger = logging.getLogger(__name__)

# The pockets file format `write_pockets_file` writes and `read_pockets_file` accepts
POCKETS_FILE_VERSION = 2


def write_pockets_file(pockets, chains, path):
    """
    Write a pockets file: the pockets plus the map from aligned chain to the pockets on it.

    Creates the file's directory if needed.

    Args:
        pockets (dict): pocket_id -> Pocket, or None for an entry given and not built.
        chains (dict): preprocess_name -> the pocket_ids on that chain; empty for a chain with none.
        path (str): File to write.

    Returns:
        None
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        json.dump({"version": POCKETS_FILE_VERSION, "pockets": serialise_pockets(pockets), "chains": chains}, f)


def read_pockets_file(path):
    """
    Read a pockets file written by `write_pockets_file`.

    Only the shape and version are checked. A hand-written file with a wrong `seq_pos` reads fine and
    gives silently wrong overlap.

    Args:
        path (str): File to read.

    Returns:
        tuple: (pockets, chains), as `write_pockets_file` takes them; pockets rebuilt as by
            `deserialise_pockets`.

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


def write_pockets_table(records, pockets, path):
    """
    Write a table of every pocket: pocket_id, chain, res_auth_ids and method.

    Creates the table's directory if needed.

    Args:
        records (list): The QTRecord dicts the pockets were built for; a pocket_id listed twice (on both
            sides, or as a given record and a Foldseek hit) is written once, from its first record.
        pockets (dict): pocket_id -> Pocket, or None for one not built, whose res_auth_ids are left
            empty.
        path (str): The table to write, tab-separated.

    Returns:
        None
    """
    log_extra = {"stage": "Getting Pockets"}
    rows = {}
    for record in records:
        pocket = pockets.get(record["pocket_id"])
        rows.setdefault(
            record["pocket_id"],
            {
                "pocket_id": record["pocket_id"],
                "chain": split_chain_info(record["chain_info"])[0],
                "res_auth_ids": ",".join(map(str, pocket.res_auth_ids)) if pocket is not None else "",
                "method": record["pocket_method"],
            },
        )
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    pd.DataFrame(list(rows.values()), columns=["pocket_id", "chain", "res_auth_ids", "method"]).to_csv(
        path, sep="\t", index=False
    )
    logger.info(f"Pockets listed in {path}", extra=log_extra)
