"""
Van der Waals contact pockets, computed directly from coordinates.

This is the `vdw` pocket method: rather than reading a precomputed interface, it compares the atoms of two
chains and keeps the residues whose van der Waals radii approach within 0.4 A. That makes it the only
interface method available for a local file, which has no PISA data. It needs two chains, so it is not
offered for an AlphaFold model, which is always a single chain. `vdw_pockets` is the method's builder.
It finds the contact residues and leaves the pocket itself to `structure.parse_pocket_from_struct`, so
`seq_pos` and coordinates are set as for every other method.

Every structure gets `setup_entities()` (through `lib.structure_manager.load_structure`), because
`lib.residues.ca_residues` and `get_polymer()` read the entity types it sets.
"""

import logging

import numpy as np

from pocketmapper.lib.residues import ca_residues, split_chain_info
from pocketmapper.lib.structure_manager import load_structure
from pocketmapper.pockets.structure import parse_pocket_from_struct, records_by_structure

logger = logging.getLogger(__name__)

# Two atoms are in contact when their van der Waals radii reach within this distance (A) of each other
CONTACT_TOLERANCE = 0.4


def atom_arrays(residues):
    """
    Coordinates and van der Waals radii of every atom of some residues, altlocs included.

    Args:
        residues (iterable): gemmi.Residue objects.

    Returns:
        tuple: (coords, radii): an (n, 3) array of positions (A) and an (n,) array of radii (A).
    """
    atoms = [atom for res in residues for atom in res]
    coords = np.array([atom.pos.tolist() for atom in atoms], dtype=float).reshape(-1, 3)
    radii = np.array([atom.element.vdw_r for atom in atoms], dtype=float)
    return coords, radii


def vdw_contact_residues(residues, partner_residues):
    """
    The residues with an atom in van der Waals contact with an atom of the partner residues.

    Two atoms are in contact when `r1 + r2 - d > -CONTACT_TOLERANCE`. Every atom of a residue counts,
    altlocs included.

    Args:
        residues (iterable): gemmi.Residue objects to test, in chain order.
        partner_residues (iterable): gemmi.Residue objects they may contact.

    Returns:
        list: Author seqids (int) of the contacting residues, in the order of `residues`, each once.
    """
    partner_coords, partner_radii = atom_arrays(partner_residues)
    if len(partner_radii) == 0:
        return []
    contacts = []
    for res in residues:
        coords, radii = atom_arrays([res])
        distances = np.linalg.norm(coords[:, None, :] - partner_coords[None, :, :], axis=-1)
        if (radii[:, None] + partner_radii[None, :] - distances > -CONTACT_TOLERANCE).any():
            contacts.append(res.seqid.num)
    return list(dict.fromkeys(contacts))


def vdw_pockets(records, structures):
    """
    Build a Pocket per record from the residues of its first chain in van der Waals contact with its second.

    The domain side is the chain's CA residues (`lib.residues.ca_residues`); the motif side is every atom
    of the second chain's first-conformer polymer residues.

    Args:
        records (list): QTRecord dicts with `pocket_method == "vdw"`.
        structures (StructureManager): Locates the records' structures.

    Returns:
        dict: pocket_id -> Pocket, for built pockets only. A record whose structure or either chain cannot
            be read is left out with a warning. One with no contacts is kept, with `pocket_exists=False`.
    """
    log_extra = {"stage": "Retrieving VdW Pockets"}
    pockets = {}
    for st, group in records_by_structure(records, structures, log_extra):
        for record in group:
            names = split_chain_info(record["chain_info"])
            chains = [None if st is None else st[0].find_chain(name) for name in names]  # first model only
            missing = [name for name, chain in zip(names, chains) if chain is None]
            if missing:
                logger.warning(
                    f"Could not read chain {','.join(missing)} of {record['struct_info']} for "
                    f"{record['pocket_id']}, skipping this entry",
                    extra=log_extra,
                )
                continue
            domain, motif = chains
            # first_conformer: a microheterogeneous position is one residue, as Foldseek reads it
            contact_ids = vdw_contact_residues(ca_residues(domain), motif.get_polymer().first_conformer())
            # Contacts come in chain order, so res_auth_ids is ascending
            pockets[record["pocket_id"]] = parse_pocket_from_struct(st, domain.name, contact_ids)
    return pockets


# Not called by any pocket method; kept deliberately for planned ATP-pocket work. Do not remove
# as dead code.
def atp_pocket_overlap(struct_path, atp_chain_id, name):
    """
    Residues of a chain's polymer that contact the ATP residue in the same chain.

    Args:
        struct_path (str): Path to a structure file, its format inferred from the extension.
        atp_chain_id (str): Chain holding both the polymer and the ATP residue.
        name (str): Key the pocket is returned under.

    Returns:
        dict: {name: Pocket}.

    Raises:
        ValueError: If the structure file is missing or unreadable, or the chain contains no ATP residue.
    """
    structure = load_structure(struct_path, {"stage": "Retrieving ATP Pockets"})
    if structure is None:
        raise ValueError(f"Structure file {struct_path} is missing or unreadable.")
    chain = structure[0][atp_chain_id]

    atp_residue = None
    for residue in chain.whole():
        if residue.name == "ATP":
            atp_residue = residue
    if atp_residue is None:
        raise ValueError("No ATP residue found in the specified chain.")

    contact_ids = vdw_contact_residues(ca_residues(chain), [atp_residue])
    return {name: parse_pocket_from_struct(structure, atp_chain_id, contact_ids)}
