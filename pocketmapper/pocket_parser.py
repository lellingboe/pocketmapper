"""
Construction of a Pocket from a structure's chain.

`parse_pocket_from_struct` both produces and extends one, so a pocket built by any of the pisa,
passthrough, vdw or whole_chain methods is interchangeable downstream. The shape itself -- which
fields exist, which are optional and why -- is declared in `pocket.py`.
"""

import gemmi
import logging
import os

from pocketmapper.constants import SINGLE_AA_CODE
from pocketmapper.pocket import Pocket, PocketResidue


def parse_pocket_from_struct(struct, chain_id, pocket_residues, pocket=None):
    """
    Build or extend a Pocket from a structure's chain.

    Walks the chain once, recording `seq_pos` -- the residue's index among the CA-bearing residues of
    that chain -- for every pocket residue. That index is what maps a pocket into the alignment, so a
    caller computing it any other way gets zero overlap with no error.

    Args:
        struct (gemmi.Structure | str): A parsed structure, or a path gemmi will read.
        chain_id (str): The chain to extract.
        pocket_residues (list | None): Author seqids in the pocket, or None for an open search, where
            every CA-bearing residue of the chain becomes the pocket. Which of the two it was is
            recorded on the returned Pocket as `whole_chain`.
        pocket (Pocket | None): An existing Pocket to extend; a new one is built when None. This is
            how a PISA pocket, which starts with residue codes but no coordinates, acquires them.

    Returns:
        Pocket: The pocket, or None if the file is missing or the chain is not in the structure.
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
        pocket = Pocket(res_auth_ids=[] if whole_chain else [str(x) for x in pocket_residues])
    pocket.whole_chain = whole_chain
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
                if str(res_id) in pocket.residues:
                    # -1 falls outside every aligned region, so the residue is ignored in the later
                    # comparison -- which matches Foldseek, which never saw it either.
                    pocket.residues[str(res_id)].seq_pos = -1
            continue
        seq_pos += 1
        res_single_code = SINGLE_AA_CODE.get(res.name, "X")
        ca_sequence.append(res_single_code)
        if whole_chain:
            pocket.res_auth_ids.append(str(res_id))
        elif res_id not in pocket_residues:  # Only recording residue info for pocket residues
            continue

        # setdefault rather than assignment: on the pisa path the residue already exists and carries
        # res_code/uniprot_pos from the interface, which must survive being given coordinates.
        residue = pocket.residues.setdefault(str(res_id), PocketResidue())
        residue.res_code = res.name
        residue.res_code_single = res_single_code
        residue.seq_pos = seq_pos
        residue.ca_coords = list(ca_atom.pos)

        pocket.pocket_exists = (
            True  # If at least one pocket residue has CA coords, we can include this pocket in the comparison
        )
        pocket.has_coords = True
    pocket.ca_sequence = "".join(ca_sequence)
    return pocket
