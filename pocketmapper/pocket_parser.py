"""
Construction of the pocket dict, the one shape every pocket method returns.

`parse_pocket_from_struct` both produces and extends that dict, so a pocket built by any of the
pisa, passthrough, vdw or whole_chain methods is interchangeable downstream. The dict holds the
top-level keys `res_auth_ids`, `ca_sequence`, `pocket_exists`, `has_coords` and `whole_chain`,
plus one entry per residue keyed by the author seqid as a **string**.
"""

import gemmi
import logging
import os

from pocketmapper.constants import SINGLE_AA_CODE


def parse_pocket_from_struct(struct, chain_id, pocket_residues, pocket=None):
    """
    Build or extend a pocket dict from a structure's chain.

    Walks the chain once, recording `seq_pos` -- the residue's index among the CA-bearing residues of
    that chain -- for every pocket residue. That index is what maps a pocket into the alignment, so a
    caller computing it any other way gets zero overlap with no error.

    Args:
        struct (gemmi.Structure | str): A parsed structure, or a path gemmi will read.
        chain_id (str): The chain to extract.
        pocket_residues (list | None): Author seqids in the pocket, or None for an open search, where
            every CA-bearing residue of the chain becomes the pocket. Which of the two it was is
            recorded on the returned dict as `whole_chain`.
        pocket (dict | None): An existing pocket dict to extend; a new one is built when None.

    Returns:
        dict: The pocket dict, or None if the file is missing or the chain is not in the structure.
    """
    stage = {"stage": "Parsing Pocket from Structure"}

    # Ensure st is a gemmi.Structure object
    if isinstance(struct, gemmi.Structure):
        st = struct
    else:
        if not os.path.exists(struct):
            logging.warning(f"Structure file {struct} does not exist.", extra=stage)
            return None
        st = gemmi.read_structure(struct)

    # Verify the specified chain exists and get it
    chain = st[0].find_chain(chain_id)  # Assuming we are interested in the first model
    if not isinstance(chain, gemmi.Chain):
        logging.critical(f"Chain {chain_id} not found in structure {struct}.", extra=stage)
        return None

    seq_pos = (
        -1
    )  # To keep track of the position in the sequence, starting at -1 so that the first residue is 0 after incrementing
    # An open search has no residue list to seed res_auth_ids from -- it is filled in as the chain is
    # walked, so it holds exactly the CA-bearing residues, in chain order.
    whole_chain = pocket_residues is None
    if whole_chain:
        pocket_residues = []
    if pocket is None:
        pocket = {
            "res_auth_ids": [] if whole_chain else [str(x) for x in pocket_residues],
            "pocket_exists": False,
            "has_coords": False,
        }
    pocket["whole_chain"] = whole_chain
    ca_sequence = []
    for res in chain:
        res_id = res.seqid.num
        ca_atom = res.get_ca()
        if ca_atom is None:  # Foldseek only uses residues with CA atom coords
            if res_id in pocket_residues:
                logging.warning(
                    f"{st.name}:{chain_id}:{res_id} ({res.name}) does not have CA coords and will be excluded from the comparison",
                    extra=stage,
                )
                if str(res_id) in pocket:
                    pocket[str(res_id)][
                        "seq_pos"
                    ] = -1  # will be ignore in the later comparison since foldseek ignores it
            continue
        seq_pos += 1
        res_single_code = SINGLE_AA_CODE.get(res.name, "X")
        ca_sequence.append(res_single_code)
        if whole_chain:
            pocket["res_auth_ids"].append(str(res_id))
        elif res_id not in pocket_residues:  # Only recording residue info for pocket residues
            continue

        if pocket.get(str(res_id)) is None:
            pocket[str(res_id)] = {}
        pocket[str(res_id)].update(
            {
                "res_code": res.name,
                "res_code_single": res_single_code,
                "seq_pos": seq_pos,  # We don't have the info to map this to a seq pos, but we can still use it in the comparison based on res_id
                "ca_coords": list(ca_atom.pos),
            }
        )

        pocket["pocket_exists"] = (
            True  # If at least one pocket residue has CA coords, we can include this pocket in the comparison
        )
        pocket["has_coords"] = True
    pocket["ca_sequence"] = "".join(ca_sequence)
    return pocket
