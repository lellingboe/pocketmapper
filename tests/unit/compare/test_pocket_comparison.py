"""
Unit tests for `pocketmapper.compare.pocket_comparison`: residues are paired by alignment position, not by label.
"""

import math

import pytest
from Bio.SVDSuperimposer import SVDSuperimposer

from pocketmapper.compare.pocket_comparison import (
    AlignmentRow,
    Context,
    aln_positions,
    compare_pocket_pair,
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
    ctx = Context(superimposer=SVDSuperimposer(), descriptions={}, identities={}, offsets={})
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
