"""
Pocket comparison: map two pockets onto a shared alignment and score their overlap.

This is step 6 of the pipeline -- it consumes the alignment table written by the Foldseek or
BLOSUM62 aligner together with the pocket dicts produced by the pocket methods, and returns the
rows that become pocket_comparison.tsv.

The alignment table's column order is a positional contract, declared once as
constants.ALIGNMENT_COLUMNS: _foldseek_alignment passes constants.FOLDSEEK_FORMAT_OUTPUT to
Foldseek's --format-output, SequenceAligner.align_records pins its DataFrame to the same list, and
each row is unpacked here into an AlignmentRow. Reorder the constant and all three move together;
edit any producer in isolation and the comparison breaks silently.

Nothing here mutates the pocket dicts it is given. Each side's projection onto the alignment is
returned as a _MappedPocket instead of being written back into the pocket, which is what lets the
same pocket be read straight out of pocket_dict on every alignment row rather than deep-copied.
"""

import json
import logging
from collections import defaultdict, namedtuple
from itertools import product
from typing import NamedTuple

import pandas as pd
from Bio.SVDSuperimposer import SVDSuperimposer
from numpy import array
from numpy import linalg as LA
from tqdm import tqdm

from pocketmapper.constants import ALIGNMENT_COLUMNS
from pocketmapper.lib import binary_similarity, full_similarity, read_blast_similarity_matrix

# One alignment row, unpacked positionally. Built with AlignmentRow(*values), so it depends on the
# column order exactly as the old row[12]-style indexing did -- it just says which column it means.
AlignmentRow = namedtuple("AlignmentRow", ALIGNMENT_COLUMNS)

# Below this sequence identity between a pocket's own CA sequence and the sequence the aligner
# reported for that chain, the two are assumed to be numbering different things (typically an
# assembly vs the asymmetric unit) and every comparison involving that pocket is dropped.
MIN_SEQ_IDENTITY = 0.8

# Every column compare_pockets can produce, in output order.
#
# A comparison stops early whenever there is nothing further to compute -- no overlapping residues, no CA
# coordinates, fewer than three overlapping residues to superpose, or an open (whole-chain) target with no
# pocket 2 to describe. Those rows leave the remaining fields empty rather than dropping the columns, so the
# written table always has this exact schema and consumers can rely on a column existing regardless of what
# any individual comparison found.
POCKET_COMPARISON_COLUMNS = [
    "pocket_1",
    "pocket_2",
    "evalue",
    "lddt",
    "pocket_1_res_ids",
    "pocket_1_len",
    "pocket_1_seq",
    "pocket_1_pct_aln",
    "pocket_2_res_ids",
    "pocket_2_len",
    "pocket_2_seq",
    "pocket_2_pct_aln",
    "overlap_count",
    "pocket_1_overlap_ids",
    "pocket_2_overlap_ids",
    "jaccard_index",
    "pocket_1_seq_overlap",
    "pocket_2_seq_overlap",
    "overlap_identity",
    "overlap_similarity_binary",
    "overlap_similarity_1_2",
    "overlap_similarity_2_1",
    "min_overlap_similarity",
    "max_overlap_similarity",
    "p2_to_p1_u",
    "p2_to_p1_t",
    "p1_to_p2_u",
    "p1_to_p2_t",
    "rmsd",
    "ca_dists",
]


class _MappedPocket(NamedTuple):
    """
    A pocket's residues projected onto one alignment row.

    positions are indices into the gapped alignment string, in res_auth_ids order; pos_by_res is the
    same thing keyed by author seqid; in_aln_count is how many of the pocket's residues landed inside
    the aligned region at all (the numerator of pocket_N_pct_aln). code_mismatches holds the residues
    whose code disagreed with the aligner's, as (aligner code, tri-code, author seqid).

    A projection depends only on the pocket and the alignment row, so it is reused across every
    pairing on that row -- but unknown_ids records one entry per *pairing*, which is why the
    mismatches are carried here rather than written straight out.
    """

    positions: list
    pos_by_res: dict
    in_aln_count: int
    code_mismatches: list


def _aln_positions(aln_seq):
    """
    Map each non-gap position of an aligned sequence to its index within the gapped string.

    The keys would be the contiguous range 0..n-1, so this is a list rather than a dict: callers
    index it with a value already range-checked against the aligned region.
    """
    return [i for i, res in enumerate(aln_seq) if res != "-"]


def _map_pocket_into_alignment(pocket, aln_seq, aln_positions, start, end):
    """
    Project a pocket's residues onto one side of an alignment row.

    seq_pos is the residue's index among its chain's CA-bearing residues; the aligner reports the
    region it actually aligned as 1-based start..end, so shifting by 1 - start puts a residue on the
    aligner's coordinates and residues outside 0..(end - start) fall outside the aligned region.

    Residues whose single-letter code disagrees with the one the aligner reported at the same
    position are collected for unknown_ids. A synthesised Foldseek-DB pocket carries no residue
    codes, so there is nothing to check against there.

    Reads `pocket` only -- the projection is returned rather than written back into it.
    """
    adj = 1 - start
    aligned_len = end - start + 1

    positions = []
    pos_by_res = []
    code_mismatches = []
    for res in pocket["res_auth_ids"]:
        entry = pocket[res]
        adj_pos = int(entry["seq_pos"]) + adj
        if not -1 < adj_pos < aligned_len:
            continue
        aln_pos = aln_positions[adj_pos]
        positions.append(aln_pos)
        pos_by_res.append((res, aln_pos))

        aln_res_code = aln_seq[aln_pos]
        if "res_code_single" in entry and aln_res_code != entry["res_code_single"]:
            code_mismatches.append((aln_res_code, entry["res_code"], res))

    return _MappedPocket(
        positions=positions,
        pos_by_res=dict(pos_by_res),
        in_aln_count=len(positions),
        code_mismatches=code_mismatches,
    )


def _record_code_mismatches(mapped, self_id, other_id, unknown_ids):
    """Fold one projection's code disagreements into unknown_ids, keyed by this particular pairing."""
    for aln_res_code, res_code, res in mapped.code_mismatches:
        unknown_ids[aln_res_code][res_code].add(f"{other_id},{self_id},{res}")


def _synthesise_target_pocket(aln):
    """
    A whole-chain stand-in for a Foldseek-database hit that has no pocket record of its own.

    Only for a database whose entries are not PDB chains (human_domains); a PDB database gets real
    PISA pockets instead, via _expand_fsdb_pdb_targets. The residues are alignment indices rather
    than author seqids and carry no codes or coordinates, so the pocket_2_* columns and the RMSD
    block are both suppressed downstream.
    """
    pocket = {
        "res_auth_ids": [str(k) for k in range(aln.tend)],
        "id_pos_codes_match": True,
        "pocket_exists": True,
        "has_coords": False,
        "whole_chain": True,
        "ca_sequence": aln.tseq,
    }
    pocket.update({str(k): {"seq_pos": k} for k in range(aln.tend)})
    return pocket


def _describe_pocket(pocket_id, pocket, cache):
    """
    The residue list, length and sequence of a pocket, memoised.

    None of the three depend on the alignment row, so on a search with thousands of hits against one
    query this is computed once instead of once per row.
    """
    described = cache.get(pocket_id)
    if described is None:
        res_ids = pocket["res_auth_ids"]
        described = (
            ",".join(res_ids),
            len(res_ids),
            "".join([pocket[res]["res_code_single"] for res in res_ids]),
        )
        cache[pocket_id] = described
    return described


def _seq_identity(pocket_id, pocket, domain, aln_seq, cache):
    """
    Identity between a pocket's own CA sequence and the sequence the aligner reported for its chain.

    Memoised on (pocket_id, domain): the aligner's qseq/tseq is a property of the aligned chain, so
    it is the same on every row that names that chain.
    """
    key = (pocket_id, domain)
    identity = cache.get(key)
    if identity is None:
        ca_sequence = pocket["ca_sequence"]
        identity = sum(map(str.__eq__, aln_seq, ca_sequence)) / len(ca_sequence)
        cache[key] = identity
    return identity


def _overlap_ids(pocket, mapped, overlap_positions):
    """The pocket's author seqids that landed on an overlapping alignment position, in pocket order."""
    pos_by_res = mapped.pos_by_res
    return [res for res in pocket["res_auth_ids"] if pos_by_res.get(res, -1) in overlap_positions]


def _superpose(p1, p2, p1_overlap_ids, p2_overlap_ids, overlap_count, sup):
    """
    Superpose the two pockets on their overlapping residues.

    Returns the transforms both ways, the RMSD and the per-residue CA distances, or nothing when
    either side has no coordinates or there are too few points to fit a rotation.
    """
    if not p1["has_coords"] or not p2["has_coords"] or overlap_count < 3:
        return {}

    x = array([p1[res]["ca_coords"] for res in p1_overlap_ids])
    y = array([p2[res]["ca_coords"] for res in p2_overlap_ids])

    sup.set(x, y)
    sup.run()
    u, t = sup.get_rotran()
    fields = {"p2_to_p1_u": u.flatten().tolist(), "p2_to_p1_t": t.tolist()}

    # TODO do this with matrix algebra instead of doing it twice
    sup.set(y, x)
    sup.run()
    u, t = sup.get_rotran()
    fields["p1_to_p2_u"] = u.flatten().tolist()
    fields["p1_to_p2_t"] = t.tolist()
    fields["rmsd"] = sup.get_rms()

    ca_dists = LA.norm(sup.get_transformed() - y, axis=1)
    fields["ca_dists"] = ",".join([str(round(dist, 3)) for dist in ca_dists])
    return fields


def parse_pocket_transform(u_cell, t_cell):
    """
    Turn a p2_to_p1_u / p2_to_p1_t cell of pocket_comparison.tsv into a gemmi-convention (u, t).

    Two conversions, both easy to get silently wrong:

    _superpose stores Biopython's SVDSuperimposer.get_rotran() output verbatim, and that rotation is
    RIGHT-multiplying (dot(coords, u) + t). gemmi.Transform LEFT-multiplies (u @ v + t), which is the
    convention Foldseek's own u already uses, so the matrix must be TRANSPOSED here. Checked against
    the two fits of the same pair in tests/e2e/e2e_results/test_core_1: |u.T - foldseek_u|max = 0.049,
    |u - foldseek_u|max = 1.08. Never hand a raw cell to gemmi.

    The cells are Python list reprs rather than the comma-joined strings alignment.tsv uses for the
    same quantities, because _superpose writes lists and to_csv reprs them.

    Returns None when the pair has no transform -- fewer than three overlapping residues, or a pocket
    with no coordinates, both of which leave _superpose returning nothing and the cells empty.
    """
    if not isinstance(u_cell, str) or not isinstance(t_cell, str):
        return None
    u = array(json.loads(u_cell)).reshape((3, 3)).T
    t = array(json.loads(t_cell))
    return u, t


def _score_overlap(aln, overlap_positions, similarity_matrix):
    """Sequence identity and the three BLOSUM62 similarity scores over the overlapping residues."""
    p1_aln_seq = "".join([aln.qaln[pos] for pos in overlap_positions])
    p2_aln_seq = "".join([aln.taln[pos] for pos in overlap_positions])

    similarity_1_2 = full_similarity(p1_aln_seq, p2_aln_seq, similarity_matrix)
    similarity_2_1 = full_similarity(p2_aln_seq, p1_aln_seq, similarity_matrix)

    return {
        "pocket_1_seq_overlap": p1_aln_seq,
        "pocket_2_seq_overlap": p2_aln_seq,
        "overlap_identity": sum(map(str.__eq__, p1_aln_seq, p2_aln_seq)) / len(overlap_positions),
        "overlap_similarity_binary": binary_similarity(p1_aln_seq, p2_aln_seq, similarity_matrix),
        "overlap_similarity_1_2": similarity_1_2,
        "overlap_similarity_2_1": similarity_2_1,
        "min_overlap_similarity": min(similarity_1_2, similarity_2_1),
        "max_overlap_similarity": max(similarity_1_2, similarity_2_1),
    }


def _compare_pocket_pair(aln, pocket_id_1, p1, p1_mapped, pocket_id_2, p2, p2_mapped, ctx):
    """
    Score one pocket against one other on a single alignment row.

    Both pockets are assumed to exist -- the caller drops empty ones before projecting them. A pair
    that simply does not overlap still returns a row: the descriptor columns are worth having, and
    the missing fields are filled in from POCKET_COMPARISON_COLUMNS at the end.
    """
    output = {
        "pocket_1": pocket_id_1,
        "pocket_2": pocket_id_2,
        "evalue": aln.evalue,
        "lddt": aln.lddt,
    }

    (
        output["pocket_1_res_ids"],
        output["pocket_1_len"],
        output["pocket_1_seq"],
    ) = _describe_pocket(pocket_id_1, p1, ctx.descriptions)
    output["pocket_1_pct_aln"] = p1_mapped.in_aln_count / output["pocket_1_len"]

    # An open search -- a whole chain rather than a pocket on it -- has no pocket 2 to describe, and
    # the chain's length would swamp both these columns and the union the Jaccard index normalises
    # by. That is flagged on the pocket itself, so one run can mix open and pocketed targets.
    p2_is_whole_chain = p2.get("whole_chain")
    if not p2_is_whole_chain:
        (
            output["pocket_2_res_ids"],
            output["pocket_2_len"],
            output["pocket_2_seq"],
        ) = _describe_pocket(pocket_id_2, p2, ctx.descriptions)
        output["pocket_2_pct_aln"] = p2_mapped.in_aln_count / output["pocket_2_len"]

    # Kept in pocket-1 order: the overlap sequences are built by indexing the alignment strings with it.
    p2_positions = set(p2_mapped.positions)
    overlap_positions = [pos for pos in p1_mapped.positions if pos in p2_positions]
    output["overlap_count"] = len(overlap_positions)
    if not overlap_positions:
        return output

    overlap_set = set(overlap_positions)
    p1_overlap_ids = _overlap_ids(p1, p1_mapped, overlap_set)
    p2_overlap_ids = _overlap_ids(p2, p2_mapped, overlap_set)
    output["pocket_1_overlap_ids"] = ",".join(p1_overlap_ids)
    output["pocket_2_overlap_ids"] = ",".join(p2_overlap_ids)

    if not p2_is_whole_chain:
        union_size = len(p1["res_auth_ids"]) + len(p2["res_auth_ids"]) - len(overlap_positions)
        output["jaccard_index"] = len(overlap_positions) / union_size

    output.update(_score_overlap(aln, overlap_positions, ctx.similarity_matrix))
    output.update(_superpose(p1, p2, p1_overlap_ids, p2_overlap_ids, len(overlap_positions), ctx.superimposer))
    return output


class _Context(NamedTuple):
    """The scoring state shared by every comparison in one call."""

    similarity_matrix: dict
    superimposer: SVDSuperimposer
    descriptions: dict
    identities: dict


def _resolve_pockets(domain, pocket_dict, preproc_to_ids):
    """The pockets sitting on one aligned chain. One chain can carry several pockets."""
    return {
        pocket_id: pocket_dict[pocket_id] for pocket_id in preproc_to_ids.get(domain) or [] if pocket_id in pocket_dict
    }


def compare_pockets(
    alignment_df,
    pocket_dict,
    preproc_to_ids,
    blosum_path,
    synthesise_target_pockets=False,
):
    """
    Compare two pockets based on the alignment that bridges them.

    blosum_path: path to a BLAST-format similarity matrix. The packaged one is
    os.path.join(os.path.dirname(pocketmapper.__file__), "blosum62.bla").

    synthesise_target_pockets: build a whole-chain pseudo-pocket per alignment row instead of looking
    the target up in pocket_dict. Only for a Foldseek database that has no target records of its own
    (see _expand_fsdb_pdb_targets); it is a property of the job.

    Whether the pocket_2_* descriptor columns and the jaccard_index are written, by contrast, is a
    property of each pocket -- see _compare_pocket_pair.

    Returns the comparison table, the residue codes the aligner and pocketmapper disagreed on, and
    the pockets whose sequence did not match the aligner's well enough to be trusted.
    """
    ctx = _Context(
        similarity_matrix=read_blast_similarity_matrix(blosum_path),
        superimposer=SVDSuperimposer(),
        descriptions={},
        identities={},
    )

    unknown_ids = defaultdict(lambda: defaultdict(set))  # for saving tri-code ids which are unknown
    incorrect_mapping = defaultdict(dict)  # for saving cases where foldseek mapping doesn't match pocketmapper sequence

    existing_calcs = set()
    output_rows = []

    for values in tqdm(alignment_df.itertuples(index=False, name=None)):
        aln = AlignmentRow(*values)
        try:
            pockets_1 = _resolve_pockets(aln.query, pocket_dict, preproc_to_ids)
            if not pockets_1:
                continue

            if synthesise_target_pockets:
                pockets_2 = {aln.target: _synthesise_target_pocket(aln)}
            else:
                pockets_2 = _resolve_pockets(aln.target, pocket_dict, preproc_to_ids)
                if not pockets_2:
                    continue

            q_positions = _aln_positions(aln.qaln)
            t_positions = _aln_positions(aln.taln)

            # Each side is projected onto this row once, then reused across every pairing.
            mapped_1 = {}
            mapped_2 = {}

            for pocket_id_1, pocket_id_2 in product(pockets_1, pockets_2):
                # Checking for A-B comparison if B-A has already been calculated
                if (pocket_id_1, pocket_id_2) in existing_calcs:
                    continue
                existing_calcs.add((pocket_id_1, pocket_id_2))

                if pocket_id_1 in incorrect_mapping or pocket_id_2 in incorrect_mapping:
                    continue

                p1 = pockets_1[pocket_id_1]
                p2 = pockets_2[pocket_id_2]
                if not p1["pocket_exists"] or not p2["pocket_exists"]:
                    continue

                p1_identity = _seq_identity(pocket_id_1, p1, aln.query, aln.qseq, ctx.identities)
                if p1_identity < MIN_SEQ_IDENTITY:
                    incorrect_mapping[pocket_id_1] = {
                        "p1_seq_identity": p1_identity,
                        "p1_seq": p1["ca_sequence"],
                        "fs_seq": aln.qseq,
                    }

                p2_identity = _seq_identity(pocket_id_2, p2, aln.target, aln.tseq, ctx.identities)
                if p2_identity < MIN_SEQ_IDENTITY:
                    incorrect_mapping[pocket_id_2] = {
                        "p2_seq_identity": p2_identity,
                        "p2_seq": p2["ca_sequence"],
                        "fs_seq": aln.tseq,
                    }

                if pocket_id_1 not in mapped_1:
                    mapped_1[pocket_id_1] = _map_pocket_into_alignment(p1, aln.qaln, q_positions, aln.qstart, aln.qend)
                if pocket_id_2 not in mapped_2:
                    mapped_2[pocket_id_2] = _map_pocket_into_alignment(p2, aln.taln, t_positions, aln.tstart, aln.tend)
                p1_mapped = mapped_1[pocket_id_1]
                p2_mapped = mapped_2[pocket_id_2]

                # The projections are shared across pairings, but unknown_ids names both pockets, so
                # each pairing contributes its own entries.
                _record_code_mismatches(p1_mapped, pocket_id_1, pocket_id_2, unknown_ids)
                _record_code_mismatches(p2_mapped, pocket_id_2, pocket_id_1, unknown_ids)

                output_rows.append(
                    _compare_pocket_pair(aln, pocket_id_1, p1, p1_mapped, pocket_id_2, p2, p2_mapped, ctx)
                )

        except Exception:
            logging.exception(
                f"Uncontrolled error calculating {aln.query} and {aln.target}",
                extra={"stage": "Pocket Comparison"},
            )
            raise

    # Reindex onto the declared schema so every row carries every column. Any column produced but not
    # declared is kept and flagged rather than silently dropped, so the list above can't drift out of date.
    pockets_df = pd.DataFrame.from_dict(output_rows)
    undeclared = [column for column in pockets_df.columns if column not in POCKET_COMPARISON_COLUMNS]
    if undeclared:
        logging.warning(
            f"Pocket comparison produced undeclared column(s) {undeclared}; add them to POCKET_COMPARISON_COLUMNS",
            extra={"stage": "Pocket Comparison"},
        )
    pockets_df = pockets_df.reindex(columns=POCKET_COMPARISON_COLUMNS + undeclared)

    return pockets_df, unknown_ids, incorrect_mapping
