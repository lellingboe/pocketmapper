import gemmi
import logging
import os

from pocketmapper.constants import SINGLE_AA_CODE


def parse_pocket_from_struct(struct, chain_id, pocket_residues, pocket=None):
    """
    Reads a structure file and return a dictionary with ca_sequence and ca_coords.

    Struct: gemmi.Structure object or path to file which will be read by gemmi.read_structure
    chain_id: the chain id to extract from the structure
    pocket_residues: list of residue ids that are in the pocket

    """
    stage = {"stage": "Parsing Pocket from Structure"}

    # Ensure st is a gemmi.Structure object
    if isinstance(struct, gemmi.Structure):
        st = struct
    else:
        if not os.path.exists(struct):
            logging.critical(f"Structure file {struct} does not exist.", extra=stage)
            exit(1)
        st = gemmi.read_structure(struct)

    # Verify the specified chain exists and get it
    chain = st[0].find_chain(chain_id)  # Assuming we are interested in the first model
    if not isinstance(chain, gemmi.Chain):
        logging.critical(f"Chain {chain_id} not found in structure {struct}.", extra=stage)
        exit(1)

    seq_pos = (
        -1
    )  # To keep track of the position in the sequence, starting at -1 so that the first residue is 0 after incrementing
    if pocket is None:
        pocket = {
            "res_auth_ids": [str(x) for x in pocket_residues],
            "pocket_exists": False,
            "has_coords": False,
        }
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
        if res_id not in pocket_residues:  # Only recording residue info for pocket residues
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
