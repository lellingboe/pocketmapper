"""
Residue naming and selection: one-letter codes, a chain's CA residues and the chain_info field of an
input entry.
"""

from pocketmapper.lib.constants import FOLDSEEK_AA_CODES


def ca_residues(chain):
    """
    The residues of a chain that have an alignment position, in chain order.

    These are the first-conformer polymer residues with a CA atom of element C (`res.get_ca()` is not
    None, so `res.get_ca().pos` is safe to read). A residue's index in this list is its `seq_pos`.
    `get_polymer()` reads the entity types that `gemmi.Structure.setup_entities()` sets: on a structure
    without them, such as a gemmi-written mmCIF with no entity data, this returns no residue.

    Args:
        chain (gemmi.Chain): The chain to walk.

    Returns:
        list: gemmi.Residue objects, valid while the structure they belong to is alive.
    """
    # first_conformer: a microheterogeneous position is one residue, as Foldseek reads it
    return [res for res in chain.get_polymer().first_conformer() if res.get_ca() is not None]


def one_letter_code(res_name):
    """
    Map a three-letter residue name to its one-letter code, as Foldseek does.

    Modified residues map to their parent amino acid (CME -> C). Names Foldseek does not list,
    including ligands, water and nucleotides, map to X.

    Args:
        res_name (str): Residue name as it appears in the structure, e.g. "ALA".

    Returns:
        str: A single uppercase letter.
    """
    return FOLDSEEK_AA_CODES.get(res_name, "X")


def split_chain_info(chain_info):
    """
    Split a chain_info field into its domain chain and its motif chain.

    Args:
        chain_info (str): The middle field of an input entry: a lone chain "A", or "A_B" where A carries
            the pocket and B is its binding partner. Chains of any length split correctly. Not None.

    Returns:
        tuple: (domain_chain, motif_chain); motif_chain is None when the entry named no partner.
    """
    chains = chain_info.split("_")
    domain_chain = chains[0]
    motif_chain = chains[1] if len(chains) > 1 else None
    return domain_chain, motif_chain
