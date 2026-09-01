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

from dataclasses import dataclass, field


@dataclass
class PocketResidue:
    """
    One residue of a pocket, keyed in `Pocket.residues` by its author seqid as a string.

    Which fields are populated depends on the producer, and the gaps are load-bearing rather than
    accidental: a pocket synthesised for a Foldseek-database hit carries `seq_pos` alone, which is
    what suppresses the descriptive pocket_2_* columns and the RMSD block downstream.
    """

    res_code: str | None = None
    res_code_single: str | None = None
    # The residue's index among the CA-bearing residues of its chain -- the value that maps a pocket
    # into the alignment, so a producer computing it any other way yields zero overlap with no error.
    # -1 marks a pocket residue with no CA atom, which Foldseek ignores and which therefore falls
    # outside every aligned region. None means no producer ever reached this residue: PISA lists
    # interface residues that need not exist in the parsed chain, and `_map_pocket_into_alignment`
    # will raise on those rather than skip them. Pre-existing, and left as a hard failure.
    seq_pos: int | None = None
    ca_coords: list | None = None
    # Written by the PISA and vdw producers, read by nothing. Kept because PISA supplies real values.
    uniprot_pos: str | int | None = None


@dataclass
class Pocket:
    """
    A pocket on a single chain. The chain itself is implicit in the pocket_id this is stored under.

    `pocket_comparison` must never write to one of these. Each side's projection onto an alignment is
    returned as a `_MappedPocket` instead, which is what lets the same Pocket be read straight out of
    the pocket collection on every alignment row rather than deep-copied.
    """

    # The ordered pocket residue list, and NOT the same thing as `list(residues)`. The two diverge on
    # the PISA path: this is seeded from the PISA interface while `residues` is filled in chain order
    # and may miss interface residues absent from the parsed chain. The order is load-bearing --
    # `_overlap_ids` returns ids in this order and the two sides' overlap lists have to correspond
    # position for position.
    res_auth_ids: list = field(default_factory=list)
    residues: dict = field(default_factory=dict)
    # The CA sequence of the WHOLE chain, not of the pocket. `_seq_identity` compares it against the
    # sequence the aligner reported for that chain, so a pocket-only sequence would never match.
    ca_sequence: str = ""
    # False until at least one pocket residue is found with CA coordinates; a pocket that never gets
    # there is skipped by the comparison rather than compared as empty.
    pocket_exists: bool = False
    has_coords: bool = False
    # True for an open search, where the whole chain stands in for a pocket. Flagged on the pocket
    # itself rather than on the run so one search can mix open and pocketed targets.
    whole_chain: bool = False
