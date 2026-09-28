"""
The declared shape of a pocket: the one structure every pocket method returns.

This module imports nothing, for the same reason `constants.py` gives: a shape that five producers
and one consumer must agree on belongs beside none of them. `pocket_parser` builds a Pocket from a
structure, `pisa_parser` and `pocket_calculator` from interface data and coordinates respectively,
and `pocket_comparison` synthesises one for a Foldseek-database hit -- all four then hand the same
thing to `pocket_comparison.compare_pockets`.

Every field carries a default. A producer that legitimately cannot fill one -- a PISA pocket before
`pocket_parser.parse_pocket_from_struct` enriches it, a synthesised database pocket with no
coordinates -- leaves it at its default rather than omitting it, so no consumer has to guess whether
to read a field with `.get` or straight indexing.
"""

from dataclasses import dataclass
from dataclasses import field


@dataclass
class PocketResidue:
    """
    One residue of a pocket, keyed in `Pocket.residues` by its author seqid as a string.

    Fields a producer cannot fill are left at None.
    """

    res_code: str | None = None
    res_code_single: str | None = None
    # Index among the chain's CA-bearing residues, which is the residue's position in the alignment.
    # -1 for a pocket residue with no CA atom; None if the residue was never found in the chain.
    seq_pos: int | None = None
    ca_coords: list | None = None
    # UniProt position, where the producer has one; otherwise -1 or None.
    uniprot_pos: str | int | None = None


@dataclass
class Pocket:
    """
    A pocket on a single chain. The chain itself is implicit in the pocket_id this is stored under.

    Read-only once built: consumers must never write to one, since the same instance is read
    without copying wherever the pocket is used.
    """

    # Pocket residue ids in ascending order. Not `list(residues)`: its order can differ, and it may
    # name ids that `residues` lacks.
    res_auth_ids: list = field(default_factory=list)
    residues: dict = field(default_factory=dict)
    # CA sequence of the WHOLE chain, not of the pocket.
    ca_sequence: str = ""
    # False until at least one pocket residue is found.
    pocket_exists: bool = False
    has_coords: bool = False
    # True when the whole chain stands in for a pocket (an open search).
    whole_chain: bool = False
