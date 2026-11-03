"""
The pocket files: the versioned pockets file and a table of every pocket.

`write_pockets_file`/`read_pockets_file` handle the versioned pockets file, which also maps each aligned
chain to the pockets on it, and `write_pockets_table` writes a table of every pocket. This module imports
no builder or manager, so readers of `pockets.json` stay light. `write_pockets_file` and `read_pockets_file`
pass every pocket through `check_pocket` (`write_pockets_table` does not), so a pocket whose residues do not
line up with its chain's sequence is refused on the way in and out.
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

    Creates the file's directory if needed. Every pocket is checked by `check_pocket` before anything is
    written, so a refused pocket leaves no file behind.

    Args:
        pockets (dict): pocket_id -> Pocket, or None for an entry given and not built.
        chains (dict): preprocess_name -> the pocket_ids on that chain; empty for a chain with none.
        path (str): File to write.

    Returns:
        None

    Raises:
        PocketMapperError: If a pocket fails `check_pocket`.
    """
    for pocket_id, pocket in pockets.items():
        if pocket is not None:
            check_pocket(pocket_id, pocket)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        json.dump({"version": POCKETS_FILE_VERSION, "pockets": serialise_pockets(pockets), "chains": chains}, f)


def read_pockets_file(path):
    """
    Read a pockets file written by `write_pockets_file`.

    The shape and version are checked, and every pocket is passed through `check_pocket`, so a corrupt or
    stale file (a wrong `seq_pos`, say) is refused rather than giving silently wrong overlap.

    Args:
        path (str): File to read.

    Returns:
        tuple: (pockets, chains), as `write_pockets_file` takes them; pockets rebuilt as by
            `deserialise_pockets`.

    Raises:
        PocketMapperError: If the file is not JSON, is not a pockets file, is of another version, or holds a
            pocket that fails `check_pocket`.
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
        pockets, chains = deserialise_pockets(serialised["pockets"]), dict(serialised["chains"])
    except (TypeError, AttributeError, KeyError, ValueError) as e:
        msg = f"Could not read the pockets file {path}: {e}"
        logger.critical(msg)
        raise PocketMapperError(msg) from e
    for pocket_id, pocket in pockets.items():
        if pocket is not None:
            check_pocket(pocket_id, pocket, path)
    return pockets, chains


def check_pocket(pocket_id, pocket, path=None):
    """
    Raise if a built pocket breaks an invariant `compare` depends on.

    A built pocket (see `structure.parse_pocket_from_struct`) has unique `res_auth_ids`, each a key of
    `residues`, and each residue's `seq_pos` indexes `ca_sequence` at its own one-letter code (compared
    case-insensitively, as `compare_pockets` does). The id order is not checked: it only affects display.
    A key of `residues` outside `res_auth_ids` is not flagged either: nothing reads it, so it cannot change
    a result. A pocket with `pocket_exists` False and no ids passes. Linear in the pocket's residues.

    Args:
        pocket_id (str): The pocket's id, named in the error.
        pocket (Pocket): The pocket to check.
        path (str or None, optional): The pockets file the pocket was read from; None when it is about to
            be written.

    Returns:
        None

    Raises:
        PocketMapperError: If an invariant is broken; the message names the pocket and residue, and for a
            pocket read from `path` says the file is corrupt or stale and pockets should be rerun.
    """
    problem = pocket_problem(pocket)
    if problem is None:
        return
    if path is None:
        msg = f"Refusing to write pocket {pocket_id}: {problem}"
    else:
        msg = f"Pocket {pocket_id} in {path} is invalid ({problem}); the file is corrupt or stale, rerun pockets"
    logger.critical(msg)
    raise PocketMapperError(msg)


def pocket_problem(pocket):
    """
    Describe the first invariant a pocket breaks.

    Args:
        pocket (Pocket): The pocket to check.

    Returns:
        str or None: The problem, naming the residue where there is one; None when the pocket is sound.
    """
    ids = pocket.res_auth_ids
    seen = set()
    for res_id in ids:
        if res_id in seen:
            return f"residue {res_id} is repeated in res_auth_ids"
        seen.add(res_id)
    n = len(pocket.ca_sequence)
    for res_id in ids:
        residue = pocket.residues.get(res_id)
        if residue is None:
            return f"residue {res_id} is in res_auth_ids but not in residues"
        seq_pos = residue.seq_pos
        if isinstance(seq_pos, bool) or not isinstance(seq_pos, int) or not 0 <= seq_pos < n:
            return f"residue {res_id} has seq_pos {seq_pos}, not an index into the {n}-residue ca_sequence"
        # A built pocket takes this code and ca_sequence from the same one_letter_code call, so None is a fault
        code = residue.res_code_single
        if not isinstance(code, str) or pocket.ca_sequence[seq_pos].upper() != code.upper():
            return f"residue {res_id} has code {code} but ca_sequence[{seq_pos}] is {pocket.ca_sequence[seq_pos]}"
    return None


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
