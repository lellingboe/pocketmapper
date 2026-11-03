"""
BLOSUM62 scoring of aligned pocket sequences.
"""

import math

from Bio.Align import substitution_matrices

BLOSUM62 = substitution_matrices.load("BLOSUM62")

# Letters scored as themselves. Any other letter is scored as X, including "*", whose -4 against every
# residue would take full_similarity below -1.
SCORED_LETTERS = frozenset(BLOSUM62.alphabet) - {"*"}


def normalise(seq):
    """
    Uppercase a sequence and turn every letter outside SCORED_LETTERS into X.

    Args:
        seq (str): An aligned sequence.

    Returns:
        str: The sequence, ready to index BLOSUM62.
    """
    return "".join(res if res in SCORED_LETTERS else "X" for res in seq.upper())


def scored_pairs(seq_a, seq_b):
    """
    The aligned residue pairs both scores use: every position whose seq_a residue has a positive self-score.

    This drops the positions where seq_a is X, the only residue BLOSUM62 scores below zero against itself.

    Args:
        seq_a (str): Sequence whose residues are kept or dropped; must be the same length as seq_b.
        seq_b (str): Second sequence.

    Returns:
        list: Normalised (seq_a residue, seq_b residue) tuples, in order.
    """
    return [(a, b) for a, b in zip(normalise(seq_a), normalise(seq_b)) if BLOSUM62[a, a] > 0]


def binary_similarity(seq_a, seq_b):
    """
    Fraction of scored positions where two aligned sequences score above zero.

    Args:
        seq_a (str): First sequence; must be the same length as seq_b. Positions where it is X are skipped.
        seq_b (str): Second sequence.

    Returns:
        float: Score in 0..1: how much of the sequence is conservatively substituted, not how strongly.
            NaN when every position is skipped.
    """
    pairs = scored_pairs(seq_a, seq_b)
    if not pairs:
        return math.nan
    return sum(BLOSUM62[a, b] > 0 for a, b in pairs) / len(pairs)


def full_similarity(seq_a, seq_b):
    """
    Substitution score of two aligned sequences, normalised per position by seq_a's self-score.

    Args:
        seq_a (str): First sequence; must be the same length as seq_b. Positions where it is X are skipped.
        seq_b (str): Second sequence.

    Returns:
        float: Mean over scored positions of the score divided by seq_a's residue scored against itself.
            In -1..1, and 1 for an identical pair however strongly its residues are conserved. NaN when
            every position is skipped.
    """
    pairs = scored_pairs(seq_a, seq_b)
    if not pairs:
        return math.nan
    return sum(BLOSUM62[a, b] / BLOSUM62[a, a] for a, b in pairs) / len(pairs)
