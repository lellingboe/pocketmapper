"""
Van der Waals contact pockets, computed directly from coordinates.

This is the `vdw` pocket method: rather than reading a precomputed interface, it walks two chains
atom by atom and keeps the residues whose van der Waals radii approach within 0.4 A. That makes it
the only method available for structures with no PISA data -- AlphaFold models and local files.
"""

import os
import logging
import gemmi
from itertools import product
from numpy.linalg import norm

from pocketmapper.constants import SINGLE_AA_CODE


class PocketCalculator:
    """
    Computes pockets from van der Waals contacts between chains.
    """

    def pocket_overlap(self, structure, domain_chain, motif_chain):
        """
        Residues of the domain chain that make van der Waals contact with the motif chain.

        `seq_pos` is counted over CA-bearing residues of the domain chain only, which is the index the
        alignment is keyed on -- see `pocket_parser.parse_pocket_from_struct`.

        Atom pairs more than 20 A apart end the scan for that residue pair, which is a large speedup and
        safe because no van der Waals radii reach that far.

        Args:
            structure (gemmi.Structure | str): A parsed structure, or a path gemmi will read.
            domain_chain (str): Chain the pocket belongs to.
            motif_chain (str): Chain it is in contact with.

        Returns:
            dict: A pocket dict with `res_auth_ids`, `ca_sequence`, `pocket_exists`, `has_coords` and one
                entry per contacting residue. None if the structure file does not exist.
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

        Unlike `pocket_overlap`, which takes contacts between two polymer chains, this walks the polymer
        against the ATP HETATM residue in the same chain.

        NOTE: retained deliberately -- not currently called by `search()` or any pocket method, and kept
        for planned ATP-pocket work. Do not remove as dead code.

        Args:
            struct_path (str): Path to an mmCIF structure.
            atp_chain_id (str): Chain holding both the polymer and the ATP residue.
            name (str): Key the pocket is returned under.

        Returns:
            dict: {name: pocket dict}.

        Raises:
            ValueError: If the chain contains no ATP residue.
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
