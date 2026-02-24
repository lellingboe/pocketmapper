"""
Docstring for pocketmapper.pocket_calculator
"""

import gemmi
from itertools import product
from numpy.linalg import norm

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
        self.vdw_radii = {"C": 1.88, "N": 1.64, "O": 1.46, "S": 1.77, "P": 1.87, "H": 1.0}
        # https://www.cgl.ucsf.edu/chimerax/docs/user/commands/clashes.html

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

    # TODO rewrite with get polymer
    def pocket_overlap(self, structure, domain_chain, motif_chain):
        """
        structure1, structure2: Biopython models
        chain1, chain2 : Strings -> Chain IDs
        """

        model = structure[0]

        pocket_res_ids = dict()
        motif_res_ids = dict()
        full_interaction = dict()

        problem_atoms = set()
        problem_residues = set()

        # Filter out hetatoms
        domain_residues = [x for x in model[domain_chain].get_residues() if x.id[0] == " "]
        motif_residues = [x for x in model[motif_chain].get_residues()]

        for res1, res2 in product(domain_residues, motif_residues):
            # atom ordering per residue: ['N', 'CA', 'C', 'O', 'CB', R1, R1, ...]
            if res1.get_resname() == "GLY":
                backbone1 = [0, 2, 3]
            else:
                backbone1 = [0, 1, 2, 3]
            if res2.get_resname() == "GLY":
                backbone2 = [0, 2, 3]
            else:
                backbone2 = [0, 1, 2, 3]

            for (pos1, atom1), (pos2, atom2) in product(enumerate(res1.get_atoms()), enumerate(res2.get_atoms())):
                distance = atom1 - atom2
                if distance > 5:
                    continue

                # Skipping pocket residues not in the standard 20
                if atom1.parent.resname not in self.single_aa_code:
                    problem_residues.add(res1.resname)
                    continue

                # VDW Radii
                try:
                    vdw1 = self.vdw_radii[atom1.id[0]]
                except KeyError:
                    problem_atoms.add(atom1.id)
                    continue
                try:
                    vdw2 = self.vdw_radii[atom2.id[0]]
                except KeyError:
                    problem_atoms.add(atom2.id)
                    continue

                vdw_range = vdw1 + vdw2
                overlap = vdw_range - distance
                if overlap > -0.4:
                    (full_interaction.setdefault(res1.id[1], dict()).setdefault(res2.id[1], set())).add(
                        (pos1 not in backbone1, pos2 not in backbone2)
                    )

                    pocket_res_ids.setdefault(res1.id[1], False)
                    if pos1 not in backbone1:
                        pocket_res_ids[res1.id[1]] = True
                    motif_res_ids.setdefault(res2.id[1], False)
                    if pos2 not in backbone2:
                        motif_res_ids[res2.id[1]] = True

        # Dict for mapping residue id to sequence position
        res_id_to_pos = {}
        res_pos_coords = {}
        for i, res in enumerate(domain_residues):
            # if res.id[1] in pocket_res_ids:
            res_id_to_pos[res.id[1]] = i
            res_pos_coords[i] = list(res.get_atoms())[1].coord.tolist()

        # mapping pocket ids to sequence position for foldseek
        pocket_res_pos = {res_id_to_pos[k]: v for k, v in pocket_res_ids.items()}

        full_interaction = self.sets_to_lists(full_interaction)  # sets are not JSON serializable

        return {
            "pocket_exists": len(pocket_res_ids) > 0,
            "pocket_res_ids": pocket_res_ids,
            "pocket_res_pos": pocket_res_pos,
            "res_id_to_pos": res_id_to_pos,
            "pocket_to_motif_sidechain_overlap": full_interaction,
            "res_pos_coords": res_pos_coords,
        }

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

    def get_pocket_data(self, struct_path, chain_id, name, pocket_ids):
        """
        Docstring for get_pocket_data

        :param self: Description
        :param struct_path: Description
        :param chain_id: Description
        :param name: Description
        :param pocket_ids: Description
        """

        structure = gemmi.read_structure(struct_path, format=gemmi.CoorFormat.Mmcif)
        structure.setup_entities()

        pass


if __name__ == "__main__":
    pc = PocketCalculator()
    struct_path = (
        r"/Users/lellingboe/Work/data/kinase_edit/atp_pocket/pocketmapper_cache/divided_structs/4WB5_A_B.cif.gz"
    )
    atp_chain_id = "A"
    name = "4WB5_A_B"
    pocket = pc.atp_pocket_overlap(struct_path, atp_chain_id, name)
