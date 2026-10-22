"""
Residue naming: one-letter codes and the chain_info field of an input entry.
"""

from pocketmapper.constants import FOLDSEEK_AA_CODES


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
