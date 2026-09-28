"""
Van der Waals contact pockets, computed directly from coordinates.

This is the `vdw` pocket method: rather than reading a precomputed interface, it walks two chains
atom by atom and keeps the residues whose van der Waals radii approach within 0.4 A. That makes it
the only interface method available for a local file, which has no PISA data. It needs two chains,
so it is not offered for an AlphaFold model, which is always a single chain. `vdw_pockets` is the
method's builder.
"""

import logging
import os
from itertools import product

import gemmi
from numpy.linalg import norm
from tqdm import tqdm

from pocketmapper.lib import one_letter_code
from pocketmapper.lib import split_chain_info
from pocketmapper.pockets.pocket import Pocket
from pocketmapper.pockets.pocket import PocketResidue

logger = logging.getLogger(__name__)


class VdWCalculator:
    """
    Computes pockets from van der Waals contacts between chains.
    """

    def pocket_overlap(self, structure, domain_chain, motif_chain):
        """
        Residues of the domain chain that make van der Waals contact with the motif chain.

        A `gemmi.Structure` passed in has `setup_entities()` called on it.

        Args:
            structure (gemmi.Structure | str): A parsed structure, or a path gemmi will read.
            domain_chain (str): Chain the pocket belongs to.
            motif_chain (str): Chain it is in contact with.

        Returns:
            Pocket: One residue per contacting residue of the domain chain. None if the structure file
                does not exist.
        """
        log_extra = {"stage": "VdW Pocket Calculation"}

        # Accept a parsed structure or a path
        if isinstance(structure, gemmi.Structure):
            pass
        else:
            if not os.path.exists(structure):
                logger.warning(f"Structure file {structure} does not exist.", extra=log_extra)
                return None
            structure = gemmi.read_structure(structure)
        structure.setup_entities()

        domain_residues = structure[0][domain_chain].get_polymer()
        motif_residues = structure[0][motif_chain].get_polymer()

        pocket = Pocket()
        ca_num = 0
        ca_sequence = []
        for res1 in domain_residues:
            if "CA" not in res1:  # only CA-bearing residues are indexed
                continue
            res_single_code = one_letter_code(res1.name)
            ca_sequence.append(res_single_code)
            for res2 in motif_residues:
                # Count contacts between the two residues
                contacts = 0
                for atom1, atom2 in product(res1, res2):
                    distance = norm(list(atom1.pos - atom2.pos))
                    if distance > 20.0:
                        break  # no van der Waals radii reach this far, so skip the rest of the pair
                    vdw_range = atom1.element.vdw_r + atom2.element.vdw_r
                    overlap = vdw_range - distance
                    if overlap > -0.4:
                        contacts += 1
                        continue

                # If contacts, add to pocket data
                if contacts > 0:
                    pocket.residues[str(res1.seqid.num)] = PocketResidue(
                        res_code=res1.name,
                        res_code_single=res_single_code,
                        uniprot_pos=-1,
                        seq_pos=ca_num,
                        ca_coords=list(res1.get_ca().pos),
                    )

            # Once per CA-bearing domain residue, in step with ca_sequence, not once per res1/res2 pair
            ca_num += 1

        # Residues are added in chain order, so this is ascending
        pocket.res_auth_ids = list(pocket.residues)
        pocket.pocket_exists = True
        pocket.has_coords = True
        pocket.ca_sequence = "".join(ca_sequence)
        return pocket

    # Not called by any pocket method; kept deliberately for planned ATP-pocket work. Do not remove
    # as dead code.
    def atp_pocket_overlap(self, struct_path, atp_chain_id, name):
        """
        Residues of a chain's polymer that contact the ATP residue in the same chain.

        Args:
            struct_path (str): Path to an mmCIF structure.
            atp_chain_id (str): Chain holding both the polymer and the ATP residue.
            name (str): Key the pocket is returned under.

        Returns:
            dict: {name: Pocket}.

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
        pocket = Pocket()
        ca_sequence = []
        for residue in chain.get_polymer():
            # Only CA-bearing residues are indexed
            if "CA" not in residue:
                continue
            ca_sequence.append(one_letter_code(residue.name))

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
                pocket.residues[str(residue.seqid.num)] = PocketResidue(
                    res_code=residue.name,
                    res_code_single=one_letter_code(residue.name),
                    uniprot_pos=-1,
                    seq_pos=ca_num,
                    ca_coords=list(residue.get_ca().pos),
                )
            ca_num += 1

        pocket.res_auth_ids = list(pocket.residues)
        pocket.pocket_exists = True
        pocket.has_coords = True
        pocket.ca_sequence = "".join(ca_sequence)
        return {name: pocket}


def vdw_pockets(records, pocket_dir):
    """
    Build a Pocket per record from the residues of its first chain in van der Waals contact with its second.

    Args:
        records (list): QTRecord dicts with `pocket_method == "vdw"`.
        pocket_dir (str): Unused.

    Returns:
        dict: pocket_id -> Pocket. A record whose structure file is missing maps to None.
    """
    pockets = {}
    pc = VdWCalculator()
    for record in tqdm(records):
        domain_chain, motif_chain = split_chain_info(record["chain_info"])
        pockets[record["pocket_id"]] = pc.pocket_overlap(
            structure=record["struct_path"],
            domain_chain=domain_chain,
            motif_chain=motif_chain,
        )
    return pockets
