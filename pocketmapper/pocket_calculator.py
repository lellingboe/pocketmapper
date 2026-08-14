"""
Docstring for pocketmapper.pocket_calculator
"""

import os
import logging
import gemmi
from itertools import product
from numpy.linalg import norm

from pocketmapper.constants import SINGLE_AA_CODE


class PocketCalculator:
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
                # Count contacts between the two residues
                contacts = 0
                for atom1, atom2 in product(res1, res2):
                    distance = norm(list(atom1.pos - atom2.pos))
                    if distance > 20.0:
                        break  # skip distant atoms to save time
                    vdw_range = atom1.element.vdw_r + atom2.element.vdw_r
                    overlap = vdw_range - distance
                    if overlap > -0.4:
                        contacts += 1
                        continue

                # If contacts, add to pocket data
                if contacts > 0:
                    res_data = {
                        "res_code": res1.name,
                        "res_code_single": res_single_code,
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
        """
        Pocket residues of a chain that contact its own bound ATP ligand.

        Unlike pocket_overlap, which takes contacts between two polymer chains, this walks the
        polymer against the ATP HETATM residue in the same chain.

        NOTE: retained deliberately -- not currently called by search() or any pocket method, and
        kept for planned ATP-pocket work. Do not remove as dead code.
        """
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
            # Foldseek skips residues with no CA atom
            if "CA" not in residue:
                continue

            # Count contacts between residue and ATP
            contacts = 0
            for atom1, atom2 in product(atp_residue, residue):
                distance = norm(list(atom1.pos - atom2.pos))
                vdw_range = atom1.element.vdw_r + atom2.element.vdw_r
                overlap = vdw_range - distance
                if overlap > -0.4:
                    contacts += 1
                    continue

            # If contacts, add to pocket data
            if contacts > 0:
                res_data = {
                    "res_code": residue.name,
                    "res_code_single": SINGLE_AA_CODE.get(residue.name, "X"),
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
