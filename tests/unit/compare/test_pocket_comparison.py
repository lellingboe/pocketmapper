"""
Unit tests for `pocketmapper.compare.pocket_comparison`: residues are paired by alignment position, not by label.
"""

import math
from dataclasses import replace

import pandas as pd
import pytest
from Bio.SVDSuperimposer import SVDSuperimposer

from pocketmapper.compare.pocket_comparison import (
    AlignmentRow,
    Context,
    aln_positions,
    compare_pocket_pair,
    compare_pockets,
    first_sequence_mismatch,
    map_pocket_into_alignment,
    score_overlap,
)
from pocketmapper.lib.constants import ALIGNMENT_COLUMNS
from pocketmapper.pockets.pocket import Pocket, PocketResidue

CHAIN_SEQ = "ACDEFGHIKLMN"

# Chain positions (0-based seq_pos) that make up the pocket on both sides
POCKET_POSITIONS = range(2, 10)


def helix_ca(k):
    """CA coordinates of chain position `k` on an ideal helix, so no three points are collinear."""
    angle = math.radians(100 * k)
    return [2.3 * math.cos(angle), 2.3 * math.sin(angle), 1.5 * k]


def make_pocket(labels):
    """A pocket on the CHAIN_SEQ chain, labelling chain position `k` as `labels[k]`, ids sorted numerically."""
    residues = {
        labels[k]: PocketResidue(res_code_single=CHAIN_SEQ[k], seq_pos=k, ca_coords=helix_ca(k))
        for k in POCKET_POSITIONS
    }
    return Pocket(
        res_auth_ids=sorted(residues, key=int),
        residues=residues,
        ca_sequence=CHAIN_SEQ,
        pocket_exists=True,
        has_coords=True,
    )


def identity_alignment():
    """An ungapped alignment row of CHAIN_SEQ against itself."""
    values = dict.fromkeys(ALIGNMENT_COLUMNS, "-")
    end = len(CHAIN_SEQ)
    values.update(query="q", target="t", qstart=1, qend=end, tstart=1, tend=end, qaln=CHAIN_SEQ, taln=CHAIN_SEQ)
    return AlignmentRow(**values)


ASCENDING = [str(k + 1) for k in range(len(CHAIN_SEQ))]

# Target labels by chain position. The query is always ASCENDING.
LABEL_TABLE = [
    pytest.param(ASCENDING, id="ascending"),
    # A middle block renumbered high, as in a T4L fusion
    pytest.param([str(5010 + k) if 4 <= k <= 6 else str(k + 1) for k in range(len(CHAIN_SEQ))], id="block-renumbered"),
    # Labels rotated along the chain, as in a circular permutant
    pytest.param([str((k + 6) % len(CHAIN_SEQ) + 1) for k in range(len(CHAIN_SEQ))], id="circular-permutant"),
]


@pytest.mark.parametrize("target_labels", LABEL_TABLE)
def test_pairs_by_alignment_position(target_labels):
    """Identical chains superpose exactly, and each overlap id pair shares a chain position, whatever the labels."""
    p1 = make_pocket(ASCENDING)
    p2 = make_pocket(target_labels)
    aln = identity_alignment()
    ctx = Context(superimposer=SVDSuperimposer(), descriptions={}, mismatches={}, offsets={})
    q_mapped = map_pocket_into_alignment(p1, aln.qaln, aln_positions(aln.qaln), aln.qstart, aln.qend)
    t_mapped = map_pocket_into_alignment(p2, aln.taln, aln_positions(aln.taln), aln.tstart, aln.tend)

    output = compare_pocket_pair(aln, "q", p1, q_mapped, "t", p2, t_mapped, ctx)

    assert output["overlap_count"] == len(POCKET_POSITIONS)
    assert output["rmsd"] == pytest.approx(0, abs=1e-6)
    assert all(float(dist) == pytest.approx(0, abs=1e-3) for dist in output["ca_dists"].split(","))
    query_ids = output["query_overlap_ids"].split(",")
    target_ids = output["target_overlap_ids"].split(",")
    assert [p1.residues[res].seq_pos for res in query_ids] == [p2.residues[res].seq_pos for res in target_ids]
    assert output["query_seq_overlap"] == output["target_seq_overlap"]


def test_score_overlap_one_side_nan():
    """An all-X query overlap leaves `_1_2` NaN, and min and max both take `_2_1`."""
    values = dict.fromkeys(ALIGNMENT_COLUMNS, "-")
    values.update(qaln="XXX", taln="WAC")
    output = score_overlap(AlignmentRow(**values), [0, 1, 2])

    assert math.isnan(output["overlap_similarity_1_2"])
    assert not math.isnan(output["overlap_similarity_2_1"])
    assert output["min_overlap_similarity"] == output["overlap_similarity_2_1"]
    assert output["max_overlap_similarity"] == output["overlap_similarity_2_1"]


LONG_SEQ = "ACDEFGHIKLMNPQRSTVWY" * 5
# Position 90 of 100 is near the chain end, where an identity threshold would let a shift through
LATE_SUBSTITUTION = LONG_SEQ[:90] + ("A" if LONG_SEQ[90] != "A" else "C") + LONG_SEQ[91:]

SEQUENCE_TABLE = [
    pytest.param(LONG_SEQ, LONG_SEQ, None, id="identical"),
    pytest.param("ACDEF", "acdef", None, id="aligner-lowercase"),
    pytest.param("acdef", "ACDEF", None, id="pocket-lowercase"),
    pytest.param("ACXEF", "ACDEF", None, id="x-in-pocket"),
    pytest.param("ACDEF", "ACXEF", None, id="x-in-aligner"),
    pytest.param("ACdEF", "ACxEF", None, id="lowercase-x"),
    pytest.param("ACDEF", "ACDGF", 3, id="substitution"),
    pytest.param(LONG_SEQ, LATE_SUBSTITUTION, 90, id="late-substitution"),
    pytest.param("ACDEF", "ACD", 3, id="aligner-shorter"),
    pytest.param("ACD", "ACDEF", 3, id="pocket-shorter"),
    pytest.param("ACDEF", "", 0, id="aligner-empty"),
    pytest.param("AXDEF", "ACDGF", 3, id="x-does-not-excuse-later-difference"),
]


@pytest.mark.parametrize(("pocket_seq", "aln_seq", "expected"), SEQUENCE_TABLE)
def test_first_sequence_mismatch(pocket_seq, aln_seq, expected):
    """The first differing index is returned, None for a match, and the result is memoised."""
    pocket = Pocket(ca_sequence=pocket_seq)
    cache = {}

    assert first_sequence_mismatch("p", pocket, "dom", aln_seq, cache) == expected
    assert cache == {("p", "dom"): expected}
    # Served from the memo, whatever sequence is passed
    assert first_sequence_mismatch("p", pocket, "dom", "ZZZ", cache) == expected


COMPARE_MISMATCH_TABLE = [
    pytest.param(True, False, {"q": "query"}, id="query"),
    pytest.param(False, True, {"t": "target"}, id="target"),
    pytest.param(True, True, {"q": "query", "t": "target"}, id="both"),
]


@pytest.mark.parametrize(("query_bad", "target_bad", "recorded"), COMPARE_MISMATCH_TABLE)
def test_compare_pockets_records_sequence_mismatch(query_bad, target_bad, recorded):
    """A pocket whose sequence differs from the aligner's is recorded, every failing side too, and gives no row."""
    shifted = CHAIN_SEQ[:10] + "WW"
    p1 = make_pocket(ASCENDING)
    p2 = make_pocket(ASCENDING)
    if query_bad:
        p1 = replace(p1, ca_sequence=shifted)
    if target_bad:
        p2 = replace(p2, ca_sequence=shifted)
    aln = identity_alignment()._replace(qseq=CHAIN_SEQ, tseq=CHAIN_SEQ)

    pockets_df, _, incorrect_mapping = compare_pockets(
        pd.DataFrame([list(aln)], columns=ALIGNMENT_COLUMNS), {"q": p1, "t": p2}, {"q": ["q"], "t": ["t"]}
    )

    assert pockets_df.empty
    assert incorrect_mapping == {
        pocket_id: {
            "side": side,
            "reason": "sequence_mismatch",
            "first_mismatch": 10,
            "pocket_seq": shifted,
            "aligner_seq": CHAIN_SEQ,
        }
        for pocket_id, side in recorded.items()
    }


def test_compare_pockets_matching_sequences_give_a_row():
    """Pockets whose sequences match the aligner's, case aside, are compared with nothing recorded."""
    aln = identity_alignment()._replace(qseq=CHAIN_SEQ, tseq=CHAIN_SEQ.lower())

    pockets_df, _, incorrect_mapping = compare_pockets(
        pd.DataFrame([list(aln)], columns=ALIGNMENT_COLUMNS),
        {"q": make_pocket(ASCENDING), "t": make_pocket(ASCENDING)},
        {"q": ["q"], "t": ["t"]},
    )

    assert len(pockets_df) == 1
    assert pockets_df.loc[0, "overlap_count"] == len(POCKET_POSITIONS)
    assert not incorrect_mapping
