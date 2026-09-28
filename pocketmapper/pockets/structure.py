"""
Construction of a Pocket from a structure's chain.

`parse_pocket_from_struct` both produces and extends one, so a pocket built by any of the pisa,
passthrough, vdw or whole_chain methods is interchangeable downstream. The shape itself -- which
fields exist, which are optional and why -- is declared in `pocket.py`.

Also holds the builders of the two methods that need nothing beyond it: `passthrough_pockets` and
`whole_chain_pockets`.
"""

import logging
import os

import gemmi
from tqdm import tqdm

from pocketmapper.lib import one_letter_code
from pocketmapper.lib import split_chain_info
from pocketmapper.pockets.pocket import Pocket
from pocketmapper.pockets.pocket import PocketResidue

logger = logging.getLogger(__name__)


def parse_pocket_from_struct(struct, chain_id, pocket_residues, pocket=None):
    """
    Build or extend a Pocket from a structure's chain.

    A `pocket` passed in is modified in place: its residues keep the fields they already carry and
    gain `seq_pos` and coordinates, and its `whole_chain` and `ca_sequence` are overwritten.

    Args:
        struct (gemmi.Structure | str): A parsed structure, or a path gemmi will read.
        chain_id (str): The chain to extract.
        pocket_residues (list | None): Author seqids in the pocket, or None to take every CA-bearing
            residue of the chain. Recorded on the returned Pocket as `whole_chain`.
        pocket (Pocket | None): An existing Pocket to extend; a new one is built when None.

    Returns:
        Pocket: The pocket, or None if the file is missing or the chain is not in the structure.
    """
    log_extra = {"stage": "Parsing Pocket from Structure"}

    # Accept a parsed structure or a path
    if isinstance(struct, gemmi.Structure):
        st = struct
    else:
        if not os.path.exists(struct):
            logger.warning(f"Structure file {struct} does not exist.", extra=log_extra)
            return None
        st = gemmi.read_structure(struct)

    # Verify the specified chain exists and get it
    chain = st[0].find_chain(chain_id)  # first model only
    if not isinstance(chain, gemmi.Chain):
        logger.critical(f"Chain {chain_id} not found in structure {struct}.", extra=log_extra)
        return None

    # seq_pos is the residue's index among the chain's CA-bearing residues -- the alignment's coordinate
    # system. Starts at -1 so the first CA-bearing residue is 0.
    seq_pos = -1
    # With no residue list, res_auth_ids is filled with every CA-bearing residue as the chain is walked
    whole_chain = pocket_residues is None
    if whole_chain:
        pocket_residues = []
    if pocket is None:
        pocket = Pocket(res_auth_ids=[] if whole_chain else [str(x) for x in pocket_residues])
    pocket.whole_chain = whole_chain
    ca_sequence = []
    # first_conformer: a microheterogeneous position is one residue, as Foldseek reads it
    for res in chain.first_conformer():
        res_id = res.seqid.num
        ca_atom = res.get_ca()
        if ca_atom is None:  # only CA-bearing residues are indexed
            if res_id in pocket_residues:
                logger.debug(
                    f"{st.name}:{chain_id}:{res_id} ({res.name}) does not have CA coords and cannot be compared",
                    extra=log_extra,
                )
                if str(res_id) in pocket.residues:
                    # -1 marks a pocket residue with no CA, which has no alignment position
                    pocket.residues[str(res_id)].seq_pos = -1
            continue
        seq_pos += 1
        res_single_code = one_letter_code(res.name)
        ca_sequence.append(res_single_code)
        if whole_chain:
            pocket.res_auth_ids.append(str(res_id))
        elif res_id not in pocket_residues:  # Only recording residue info for pocket residues
            continue

        # setdefault: a residue already on `pocket` keeps the fields it carries
        residue = pocket.residues.setdefault(str(res_id), PocketResidue())
        residue.res_code = res.name
        residue.res_code_single = res_single_code
        residue.seq_pos = seq_pos
        residue.ca_coords = list(ca_atom.pos)

        pocket.pocket_exists = True
        pocket.has_coords = True
    pocket.ca_sequence = "".join(ca_sequence)
    return pocket


def passthrough_pockets(records, pocket_dir):
    """
    Build a Pocket per record from the residue ids listed in its `residue_info`.

    Args:
        records (list): QTRecord dicts with `pocket_method == "passthrough"`.
        pocket_dir (str): Unused.

    Returns:
        dict: pocket_id -> Pocket, with `res_auth_ids` in ascending order whatever order they were
            listed in. A record whose structure or chain cannot be read, or that names a residue the
            chain lacks or that has no CA atom, is skipped with a warning.
    """
    log_extra = {"stage": "Retrieving passthrough Pockets"}

    pockets = {}
    for record in tqdm(records):
        domain_chain, _ = split_chain_info(record["chain_info"])
        pocket = parse_pocket_from_struct(
            struct=record["struct_path"],
            chain_id=domain_chain,
            # Ascending, whatever order the ids were typed in
            pocket_residues=sorted(int(x) for x in record["residue_info"].split(",")),
        )
        # A missing structure or chain gives None
        if pocket is None:
            logger.warning(
                f"Could not parse chain {domain_chain} of {record['struct_info']} for {record['pocket_id']}, "
                "skipping this entry",
                extra=log_extra,
            )
            continue
        # residues holds only the requested ids found with a CA atom
        unusable = [res_id for res_id in pocket.res_auth_ids if res_id not in pocket.residues]
        if unusable:
            logger.warning(
                f"Residue(s) {','.join(unusable)} of {record['pocket_id']} are not in chain {domain_chain} "
                "or have no CA atom, skipping this entry",
                extra=log_extra,
            )
            continue
        pockets[record["pocket_id"]] = pocket
    return pockets


def whole_chain_pockets(records, pocket_dir):
    """
    Build a Pocket per record from every CA-bearing residue of its chain.

    Args:
        records (list): QTRecord dicts with `pocket_method == "whole_chain"`.
        pocket_dir (str): Unused.

    Returns:
        dict: pocket_id -> Pocket. A record whose structure or chain cannot be read is skipped with a
            warning.
    """
    log_extra = {"stage": "Retrieving whole chain Pockets"}

    pockets = {}
    for record in tqdm(records):
        domain_chain, _ = split_chain_info(record["chain_info"])
        pocket = parse_pocket_from_struct(
            struct=record["struct_path"],
            chain_id=domain_chain,
            pocket_residues=None,  # None means the whole chain
        )
        # A missing structure or chain gives None
        if pocket is None:
            logger.warning(
                f"Could not parse chain {record['chain_info']} of {record['struct_info']} for {record['pocket_id']}, "
                "skipping this entry",
                extra=log_extra,
            )
            continue
        pockets[record["pocket_id"]] = pocket
    return pockets
