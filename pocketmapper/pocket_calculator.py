"""
Docstring for pocketmapper.pocket_calculator
"""

import os
import logging
import gemmi
from itertools import product
from numpy.linalg import norm

from pocketmapper.constants import SINGLE_AA_CODE

# from pocketmapper.lib import jsonify_dict


class PocketCalculator:
    def __init__(self):
        self.single_aa_code = {
            "CYS": "C",
            "ASP": "D",
            "SER": "S",
            "GLN": "Q",
            "LYS": "K",
            "ILE": "I",
            "PRO": "P",
            "THR": "T",
            "PHE": "F",
            "ASN": "N",
            "GLY": "G",
            "HIS": "H",
            "LEU": "L",
            "ARG": "R",
            "TRP": "W",
            "ALA": "A",
            "VAL": "V",
            "GLU": "E",
            "TYR": "Y",
            "MET": "M",
            "SEP": "S",  # phospho
            "TPO": "T",  # phospho
        }

    def sets_to_lists(self, item):
        """
        Recursively looks for sets in a dictionary and turns then into lists
        This allows dicts with sets to become JSON serializeable
        """
        if isinstance(item, set):
            return list(item)
        elif isinstance(item, dict):
            return {k: self.sets_to_lists(v) for k, v in item.items()}
        else:
            return item

    def pocket_overlap(self, structure, domain_chain, motif_chain):
        """
        structure1: Gemmi structure or path to mmcif file
        chain1, chain2 : Strings -> Chain IDs
        """
        stage = {"stage": "VdW pocket calculation"}

        # Ensure st is a gemmi.Structure object
        if isinstance(structure, gemmi.Structure):
            pass
        else:
            if not os.path.exists(structure):
                logging.warning(f"Structure file {structure} does not exist.", extra=stage)
                return None
            structure = gemmi.read_structure(structure)
        structure.setup_entities()

        domain_residues = structure[0][domain_chain].get_polymer()
        motif_residues = structure[0][motif_chain].get_polymer()

        pocket_data = {}
        ca_num = 0
        ca_sequence = []
        for res1 in domain_residues:
            if "CA" not in res1:  #
                continue
            res_single_code = SINGLE_AA_CODE.get(res1.name, "X")
            ca_sequence.append(res_single_code)
            for res2 in motif_residues:
                # atom ordering per residue: ['N', 'CA', 'C', 'O', 'CB', R1, R1, ...]
                # if res1.name == "GLY":
                #    backbone1 = [0, 2, 3]
                # else:
                #    backbone1 = [0, 1, 2, 3]
                # if res2.name == "GLY":
                #    backbone2 = [0, 2, 3]
                # else:
                #    backbone2 = [0, 1, 2, 3]
                # Count contacts between residue and ATP
                contacts = 0
                for atom1, atom2 in product(res1, res2):
                    distance = norm(list(atom1.pos - atom2.pos))
                    if distance > 20.0:
                        break  # skip distant atoms to save time
                    vdw_range = atom1.element.vdw_r + atom2.element.vdw_r
                    overlap = vdw_range - distance
                    # if distance < 4.0:
                    if overlap > -0.4:
                        contacts += 1
                        continue

                # If contacts, add to pocket data
                if contacts > 0:
                    res_data = {
                        "res_code": res1.name,
                        "res_code_single": self.single_aa_code.get(res1.name, "X"),
                        "uniprot_pos": -1,
                        "seq_pos": ca_num,
                        "ca_coords": list(res1.get_ca().pos),
                    }
                    pocket_data[str(res1.seqid.num)] = res_data

            # Counts CA-bearing domain residues, so it must advance once per res1 (in step with
            # ca_sequence above), not once per res1/res2 pair. Downstream, seq_pos is the residue's
            # index within this chain's CA sequence, which is what maps it into the alignment.
            ca_num += 1

        pocket_data["res_auth_ids"] = [str(x) for x in pocket_data.keys()]
        pocket_data["id_pos_codes_match"] = True
        pocket_data["pocket_exists"] = True
        pocket_data["has_coords"] = True
        pocket_data["ca_sequence"] = "".join(ca_sequence)
        return pocket_data

    def atp_pocket_overlap(self, struct_path, atp_chain_id, name):
        structure = gemmi.read_structure(struct_path, format=gemmi.CoorFormat.Mmcif)
        structure.setup_entities()
        chain = structure[0][atp_chain_id]

        # Identify ATP residue
        atp_residue = None
        for residue in chain.whole():
            if residue.name == "ATP":
                atp_residue = residue
        if atp_residue is None:
            raise ValueError("No ATP residue found in the specified chain.")

        ca_num = 0
        pocket_data = {}
        for residue in chain.get_polymer():
            # FOldseek skips residues with no CA atom
            if "CA" not in residue:
                continue

            # Count contacts between residue and ATP
            contacts = 0
            for atom1, atom2 in product(atp_residue, residue):
                distance = norm(list(atom1.pos - atom2.pos))
                vdw_range = atom1.element.vdw_r + atom2.element.vdw_r
                overlap = vdw_range - distance
                # if distance < 4.0:
                if overlap > -0.4:
                    contacts += 1
                    continue

            # If contacts, add to pocket data
            if contacts > 0:
                res_data = {
                    "res_code": residue.name,
                    "res_code_single": self.single_aa_code.get(residue.name, "X"),
                    "uniprot_pos": -1,
                    "seq_pos": ca_num,
                    "ca_coords": list(residue.get_ca().pos),
                }
                pocket_data[str(residue.seqid.num)] = res_data
            ca_num += 1

        pocket_data["res_auth_ids"] = [str(x) for x in pocket_data.keys()]
        pocket_data["id_pos_codes_match"] = True
        pocket_data["pocket_exists"] = True
        pocket_data["has_coords"] = True
        return {name: pocket_data}


if __name__ == "__main__":
    pc = PocketCalculator()
    struct_path = (
        r"/Users/lellingboe/Work/data/kinase_edit/atp_pocket/pocketmapper_cache/divided_structs/4WB5_A_B.cif.gz"
    )
    atp_chain_id = "A"
    name = "4WB5_A_B"
    pocket = pc.atp_pocket_overlap(struct_path, atp_chain_id, name)
