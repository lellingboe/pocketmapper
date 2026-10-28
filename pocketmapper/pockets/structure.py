"""
Construction of a Pocket from a structure's chain.

`parse_pocket_from_struct` both produces and extends one, so a pocket built by any of the pisa,
passthrough, vdw or whole_chain methods is interchangeable downstream. It is the one place `seq_pos` is
set, from `lib.residues.ca_residues`. The shape itself -- which fields exist, which are optional and
why -- is declared in `pocket.py`.

Also holds the structure reading every builder shares (`records_by_structure`, through
`lib.structure_manager.load_structure`) and the builders of the two methods that need nothing beyond it:
`passthrough_pockets` and `whole_chain_pockets`. Both run `residue_list_pockets` and differ only in the
residue list they pass it.
"""

import logging

from tqdm import tqdm

from pocketmapper.lib.residues import ca_residues, one_letter_code, split_chain_info
from pocketmapper.lib.structure_manager import load_structure
from pocketmapper.pockets.pocket import Pocket, PocketResidue

logger = logging.getLogger(__name__)


def records_by_structure(records, log_extra):
    """
    Group records by `struct_path` and read each group's structure once.

    Groups come in order of each path's first appearance, and keep the records' order within a group.
    Each structure is read only when its group is reached, so one is held in memory at a time.

    Args:
        records (list): QTRecord dicts; `struct_path` is read.
        log_extra (dict): `extra` for `load_structure`'s warning about a missing or unreadable file.

    Yields:
        tuple: (structure, group): the `load_structure` result for one path, None if the file is
            missing or unusable, and the list of records with that path.
    """
    groups = {}
    for record in records:
        groups.setdefault(record["struct_path"], []).append(record)
    for path, group in tqdm(groups.items()):
        yield load_structure(path, log_extra), group


def parse_pocket_from_struct(struct, chain_id, pocket_residues, pocket=None):
    """
    Build or extend a Pocket from a structure's chain.

    A `pocket` passed in is modified in place: its residues keep the fields they already carry and
    gain `seq_pos` and coordinates, and its `whole_chain` and `ca_sequence` are overwritten. A structure
    passed in has `setup_entities()` called on it (see `lib.structure_manager.load_structure`).

    Args:
        struct (gemmi.Structure or str): A parsed structure, or a path gemmi will read.
        chain_id (str): The chain to extract.
        pocket_residues (list or None): Author seqids (int) in the pocket, or None to take every residue
            of `lib.residues.ca_residues`. Recorded on the returned Pocket as `whole_chain`.
        pocket (Pocket or None, optional): An existing Pocket to extend; a new one is built when None.

    Returns:
        Pocket: The pocket, or None if the file is missing or unusable, or the chain is not in the structure.
    """
    log_extra = {"stage": "Parsing Pocket from Structure"}

    st = load_structure(struct, log_extra)
    if st is None:
        return None
    chain = st[0].find_chain(chain_id)  # first model only
    if chain is None:
        logger.warning(f"Chain {chain_id} not found in structure {st.name}.", extra=log_extra)
        return None

    # With no residue list, res_auth_ids is filled with every CA residue as the chain is walked
    whole_chain = pocket_residues is None
    wanted = set() if whole_chain else set(pocket_residues)
    if pocket is None:
        pocket = Pocket(res_auth_ids=[] if whole_chain else [str(x) for x in pocket_residues])
    pocket.whole_chain = whole_chain

    # seq_pos is the residue's index among the CA residues -- the alignment's coordinate system
    reached = set()
    ca_sequence = []
    for seq_pos, res in enumerate(ca_residues(chain)):
        res_id = res.seqid.num
        res_single_code = one_letter_code(res.name)
        ca_sequence.append(res_single_code)
        if whole_chain:
            pocket.res_auth_ids.append(str(res_id))
        elif res_id not in wanted:  # Only recording residue info for pocket residues
            continue
        reached.add(res_id)

        # setdefault: a residue already on `pocket` keeps the fields it carries
        residue = pocket.residues.setdefault(str(res_id), PocketResidue())
        residue.res_code = res.name
        residue.res_code_single = res_single_code
        residue.seq_pos = seq_pos
        residue.ca_coords = list(res.get_ca().pos)

        pocket.pocket_exists = True
        pocket.has_coords = True
    pocket.ca_sequence = "".join(ca_sequence)

    # -1 marks a requested residue the chain has but the walk did not reach (no CA atom, or not polymer),
    # which has no alignment position. One the chain lacks keeps seq_pos None.
    unreached = wanted - reached
    for res in chain.first_conformer():
        res_id = res.seqid.num
        if res_id not in unreached:
            continue
        logger.debug(
            f"{st.name}:{chain_id}:{res_id} ({res.name}) has no CA in the polymer and cannot be compared",
            extra=log_extra,
        )
        if str(res_id) in pocket.residues:
            pocket.residues[str(res_id)].seq_pos = -1
    return pocket


def residue_list_pockets(records, residue_ids, log_extra):
    """
    Build a Pocket per record from one structure chain and a residue list, reading each structure once.

    Args:
        records (list): QTRecord dicts.
        residue_ids (callable): Takes a record and returns the author seqids (int) of its pocket in
            ascending order, or None for the whole chain.
        log_extra (dict): `extra` for the warnings, naming the builder's stage.

    Returns:
        dict: pocket_id -> Pocket, for built pockets only. A record whose structure or chain cannot be
            read, or that names a residue the chain lacks or that has no CA atom, is left out with a
            warning.
    """
    pockets = {}
    for st, group in records_by_structure(records, log_extra):
        for record in group:
            domain_chain, _ = split_chain_info(record["chain_info"])
            pocket = None if st is None else parse_pocket_from_struct(st, domain_chain, residue_ids(record))
            if pocket is None:
                logger.warning(
                    f"Could not read chain {domain_chain} of {record['struct_info']} for {record['pocket_id']}, "
                    "skipping this entry",
                    extra=log_extra,
                )
                continue
            # residues holds only the requested ids found with a CA atom, so a whole chain never fails here
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


def passthrough_pockets(records, pocket_dir):
    """
    Build a Pocket per record from the residue ids listed in its `residue_info`.

    Args:
        records (list): QTRecord dicts with `pocket_method == "passthrough"`.
        pocket_dir (str): Unused.

    Returns:
        dict: pocket_id -> Pocket, with `res_auth_ids` in ascending order whatever order they were
            listed in. Records that cannot be built are left out (see `residue_list_pockets`).
    """
    return residue_list_pockets(
        records,
        # Ascending, whatever order the ids were typed in
        lambda record: sorted(int(x) for x in record["residue_info"].split(",")),
        {"stage": "Retrieving Passthrough Pockets"},
    )


def whole_chain_pockets(records, pocket_dir):
    """
    Build a Pocket per record from every CA residue of its chain.

    Args:
        records (list): QTRecord dicts with `pocket_method == "whole_chain"`.
        pocket_dir (str): Unused.

    Returns:
        dict: pocket_id -> Pocket. Records that cannot be built are left out (see
            `residue_list_pockets`).
    """
    return residue_list_pockets(records, lambda record: None, {"stage": "Retrieving Whole Chain Pockets"})
