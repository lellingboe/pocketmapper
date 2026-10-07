"""
Pocket comparison: map two pockets onto a shared alignment and score their overlap.

Takes an alignment table and the Pockets on its chains, and returns one comparison row per pocket
pair. Each alignment row is unpacked by position into an AlignmentRow, so the column order declared
as constants.ALIGNMENT_COLUMNS is a contract.

Nothing here mutates the Pockets it is given. Each side's projection onto the alignment is returned
as a MappedPocket instead, so the same Pocket is read on every alignment row without a copy.
"""

import logging
from collections import defaultdict
from collections import namedtuple
from itertools import product
from typing import NamedTuple

import pandas as pd
from Bio.SVDSuperimposer import SVDSuperimposer
from numpy import array
from numpy import linalg as LA
from tqdm import tqdm

from pocketmapper.constants import ALIGNMENT_COLUMNS
from pocketmapper.exceptions import PocketMapperError
from pocketmapper.lib import binary_similarity
from pocketmapper.lib import full_similarity
from pocketmapper.lib import read_blast_similarity_matrix
from pocketmapper.lib import read_offset_table
from pocketmapper.lib import seq_to_uniprot_map
from pocketmapper.pockets.pocket import Pocket
from pocketmapper.pockets.pocket import PocketResidue

logger = logging.getLogger(__name__)

# One alignment row. Built with AlignmentRow(*values), so it depends on the column order.
AlignmentRow = namedtuple("AlignmentRow", ALIGNMENT_COLUMNS)

# Below this sequence identity between a pocket's own CA sequence and the sequence the aligner
# reported for that chain, the two are assumed to be numbering different things (typically an
# assembly vs the asymmetric unit) and every comparison involving that pocket is dropped.
MIN_SEQ_IDENTITY = 0.8

# Every column compare_pockets can produce, in output order. A comparison that stops early leaves the
# remaining fields empty, so every row carries every column.
POCKET_COMPARISON_COLUMNS = [
    "query",
    "target",
    "evalue",
    "lddt",
    "overlap_count",
    "jaccard_index",
    "overlap_identity",
    "max_overlap_similarity",
    "rmsd",
    "ca_dists",
    "query_pct_aln",
    "target_pct_aln",
    "query_len",
    "target_len",
    "query_seq",
    "target_seq",
    "query_seq_overlap",
    "target_seq_overlap",
    "query_res_ids",
    "target_res_ids",
    "query_overlap_ids",
    "target_overlap_ids",
    "overlap_similarity_binary",
    "overlap_similarity_1_2",
    "overlap_similarity_2_1",
    "min_overlap_similarity",
    "target_to_query_u",
    "target_to_query_t",
    "query_to_target_u",
    "query_to_target_t",
]


class MappedPocket(NamedTuple):
    """
    A pocket's residues projected onto one alignment row.

    positions are indices into the gapped alignment string, in res_auth_ids order; pos_by_res is the
    same thing keyed by author seqid; in_aln_count is how many of the pocket's residues landed inside
    the aligned region at all (the numerator of pocket_N_pct_aln). code_mismatches holds the residues
    whose code disagreed with the aligner's, as (aligner code, tri-code, author seqid).

    Mismatches are carried rather than recorded: a projection is reused across every pairing on its
    row, but unknown_ids records one entry per pairing.
    """

    positions: list
    pos_by_res: dict
    in_aln_count: int
    code_mismatches: list


def aln_positions(aln_seq):
    """
    Map each non-gap position of an aligned sequence to its index within the gapped string.

    Args:
        aln_seq (str): One side of an alignment row, with "-" for gaps.

    Returns:
        list: Index into the gapped string for each ungapped position, in order.
    """
    return [i for i, res in enumerate(aln_seq) if res != "-"]


def map_pocket_into_alignment(pocket, aln_seq, aln_positions, start, end):
    """
    Project a pocket's residues onto one side of an alignment row.

    Args:
        pocket (Pocket): The pocket to project.
        aln_seq (str): This side's gapped alignment string.
        aln_positions (list): As returned by `aln_positions` for `aln_seq`.
        start (int): 1-based first aligned residue on this side.
        end (int): 1-based last aligned residue on this side.

    Returns:
        MappedPocket: The projection. Residues outside the aligned region are left out, and a residue
            with no single-letter code is never counted as a code mismatch.
    """
    # seq_pos indexes the chain's CA-bearing residues from 0; the aligned region is 1-based start..end
    adj = 1 - start
    aligned_len = end - start + 1

    positions = []
    pos_by_res = []
    code_mismatches = []
    for res in pocket.res_auth_ids:
        entry = pocket.residues[res]
        adj_pos = int(entry.seq_pos) + adj
        if not -1 < adj_pos < aligned_len:
            continue
        aln_pos = aln_positions[adj_pos]
        positions.append(aln_pos)
        pos_by_res.append((res, aln_pos))

        aln_res_code = aln_seq[aln_pos]
        if entry.res_code_single is not None and aln_res_code != entry.res_code_single:
            code_mismatches.append((aln_res_code, entry.res_code, res))

    return MappedPocket(
        positions=positions,
        pos_by_res=dict(pos_by_res),
        in_aln_count=len(positions),
        code_mismatches=code_mismatches,
    )


def record_code_mismatches(mapped, self_id, other_id, unknown_ids):
    """
    Fold one projection's code disagreements into unknown_ids, keyed by this particular pairing.

    Args:
        mapped (MappedPocket): The projection whose mismatches are being recorded.
        self_id (str): pocket_id of the pocket that was projected.
        other_id (str): pocket_id of the pocket it is being compared against.
        unknown_ids (dict): Nested accumulator, mutated in place.

    Returns:
        None: Mutates `unknown_ids`.
    """
    for aln_res_code, res_code, res in mapped.code_mismatches:
        unknown_ids[aln_res_code][res_code].add(f"{other_id},{self_id},{res}")


def synthesise_target_pocket(aln, ctx):
    """
    A whole-chain stand-in for a Foldseek-database hit that has no pocket record of its own.

    Args:
        aln (AlignmentRow): The row to build the pseudo-pocket from; reads `target`, `tend` and `tseq`.
        ctx (Context): Scoring state, read for `offsets`.

    Returns:
        Pocket: `whole_chain` set, `has_coords` false, and residues carrying `seq_pos` alone. Residue ids
            are 1-indexed UniProt positions when `ctx.offsets` is set, else 0-indexed positions within
            the entry.

    Raises:
        PocketMapperError: The database ships an offset table but it does not describe this hit --
            the table and the database have drifted apart.
    """
    # An entry is a domain carved out of a UniProt sequence, so its own numbering means nothing outside
    # PocketMapper; UniProt positions do
    res_ids = [str(k) for k in range(aln.tend)]
    if ctx.offsets:
        res_ids = uniprot_res_ids(aln, ctx.offsets)

    return Pocket(
        res_auth_ids=res_ids,
        # The ids are labels and seq_pos the 0-indexed position, which is why the ids can be renumbered.
        # A repeated id would collapse two residues into one entry and give the survivor the wrong
        # seq_pos. No codes or coordinates, which is what suppresses the code-mismatch check and the
        # RMSD block. uniprot_pos stays unset: it would duplicate the key.
        residues={res_id: PocketResidue(seq_pos=k) for k, res_id in enumerate(res_ids)},
        pocket_exists=True,
        has_coords=False,
        whole_chain=True,
        ca_sequence=aln.tseq,
    )


def uniprot_res_ids(aln, offsets):
    """
    The UniProt residue numbers of a database entry's first `tend` positions.

    Args:
        aln (AlignmentRow): The row, read for `target` and `tend`.
        offsets (dict): entry name -> region spec, from `read_offset_table`.

    Returns:
        list: Residue ids as strings, one per position 0..tend-1.

    Raises:
        PocketMapperError: The entry is absent from the table, its spec is malformed, or it is shorter
            than the alignment reaches. All three mean the table and the database have drifted apart.
    """
    # Not memoised by entry: across the five human_domains e2e cases an entry appears on about one row
    # per query (704 rows / 704 entries, 791/791, 830/830, 796/796, 2,219/809), so a memo would miss
    log_extra = {"stage": "Pocket Comparison"}
    domain = offsets.get(aln.target)
    if domain is None:
        msg = (
            f"Foldseek database entry {aln.target} is missing from the offset table shipped with the "
            "database; the two have drifted apart. Refresh the offset table alongside the database."
        )
        logger.critical(msg, extra=log_extra)
        raise PocketMapperError(msg)

    try:
        uniprot_map = seq_to_uniprot_map(domain)
    except ValueError as error:
        msg = f"Malformed offset table entry for {aln.target}: {domain!r} ({error})"
        logger.critical(msg, extra=log_extra)
        raise PocketMapperError(msg)

    # Only a spec that is too short is detectable: the row carries no tlen. Without this check the
    # failure would be a bare KeyError naming neither the entry nor the table.
    if len(uniprot_map) < aln.tend:
        msg = (
            f"Offset table entry for {aln.target} spans {len(uniprot_map)} residues but the alignment "
            f"reaches position {aln.tend}; the table and the database have drifted apart."
        )
        logger.critical(msg, extra=log_extra)
        raise PocketMapperError(msg)

    return [str(uniprot_map[k]) for k in range(aln.tend)]


def describe_pocket(pocket_id, pocket, cache):
    """
    The residue list, length and sequence of a pocket, memoised.

    Args:
        pocket_id (str): Cache key.
        pocket (Pocket): The pocket to describe.
        cache (dict): Memo, updated in place.

    Returns:
        tuple: (comma-joined res_auth_ids, residue count, single-letter sequence).
    """
    described = cache.get(pocket_id)
    if described is None:
        res_ids = pocket.res_auth_ids
        described = (
            ",".join(res_ids),
            len(res_ids),
            "".join([pocket.residues[res].res_code_single for res in res_ids]),
        )
        cache[pocket_id] = described
    return described


def seq_identity(pocket_id, pocket, domain, aln_seq, cache):
    """
    Identity between a pocket's own CA sequence and the sequence the aligner reported for its chain.

    Args:
        pocket_id (str): Half the cache key.
        pocket (Pocket): The pocket, read for `ca_sequence`.
        domain (str): The chain's `preprocess_name`; the other half of the cache key.
        aln_seq (str): The ungapped sequence the aligner reported for that chain, the same on every
            row that names it.
        cache (dict): Memo keyed by (pocket_id, domain), updated in place.

    Returns:
        float: Fraction of positions that agree.
    """
    key = (pocket_id, domain)
    identity = cache.get(key)
    if identity is None:
        ca_sequence = pocket.ca_sequence
        identity = sum(map(str.__eq__, aln_seq, ca_sequence)) / len(ca_sequence)
        cache[key] = identity
    return identity


def overlap_ids(pocket, mapped, overlap_positions):
    """
    The pocket's author seqids that landed on an overlapping alignment position, in pocket order.

    Args:
        pocket (Pocket): The pocket, read for `res_auth_ids`.
        mapped (MappedPocket): Its projection onto this row.
        overlap_positions (set): Alignment positions shared by both pockets.

    Returns:
        list: Author seqids, ordered to match the other side's list position for position.
    """
    pos_by_res = mapped.pos_by_res
    return [res for res in pocket.res_auth_ids if pos_by_res.get(res, -1) in overlap_positions]


def format_vector(values):
    """
    Comma-join values at three decimals, the format Foldseek writes its `u` and `t` in.

    Args:
        values (iterable): Floats to format.

    Returns:
        str: e.g. `0.846,-0.524,0.910`.
    """
    return ",".join(f"{value:.3f}" for value in values)


def superpose(p1, p2, p1_overlap_ids, p2_overlap_ids, overlap_count, sup):
    """
    Superpose the two pockets on their overlapping residues.

    Args:
        p1 (Pocket): Query pocket, read for `has_coords` and CA coordinates.
        p2 (Pocket): Target pocket.
        p1_overlap_ids (list): Query author seqids in overlap order.
        p2_overlap_ids (list): Target author seqids in the same order.
        overlap_count (int): How many residues overlap.
        sup (Bio.SVDSuperimposer.SVDSuperimposer): Reused across calls.

    Returns:
        dict: The transforms both ways, the RMSD and the per-residue CA distances. Empty when either
            side has no coordinates or there are fewer than three points to fit a rotation.
    """
    if not p1.has_coords or not p2.has_coords or overlap_count < 3:
        return {}

    x = array([p1.residues[res].ca_coords for res in p1_overlap_ids])
    y = array([p2.residues[res].ca_coords for res in p2_overlap_ids])

    sup.set(x, y)
    sup.run()
    u, t = sup.get_rotran()
    fields = {"target_to_query_u": format_vector(u.flatten()), "target_to_query_t": format_vector(t)}

    # TODO do this with matrix algebra instead of doing it twice
    sup.set(y, x)
    sup.run()
    u, t = sup.get_rotran()
    fields["query_to_target_u"] = format_vector(u.flatten())
    fields["query_to_target_t"] = format_vector(t)
    fields["rmsd"] = sup.get_rms()

    ca_dists = LA.norm(sup.get_transformed() - y, axis=1)
    fields["ca_dists"] = ",".join([str(round(dist, 3)) for dist in ca_dists])
    return fields


def parse_pocket_transform(u_cell, t_cell):
    """
    Turn a target_to_query_u / target_to_query_t cell of pocket_comparison.tsv into a gemmi-convention (u, t).

    Args:
        u_cell (str): The `target_to_query_u` cell, nine comma-joined floats.
        t_cell (str): The `target_to_query_t` cell, three comma-joined floats.

    Returns:
        tuple: (u, t) ready for gemmi, or None when the cells are empty because the pair has no
            transform.
    """
    # The cells hold SVDSuperimposer.get_rotran() verbatim, a RIGHT-multiplying rotation
    # (dot(coords, u) + t). gemmi.Transform LEFT-multiplies (u @ v + t), as Foldseek's u does, so u is
    # transposed. Checked against both fits of one pair in e2e test_core_1: |u.T - foldseek_u|max =
    # 0.049, |u - foldseek_u|max = 1.08.
    if not isinstance(u_cell, str) or not isinstance(t_cell, str):
        return None
    u = array([float(x) for x in u_cell.split(",")]).reshape((3, 3)).T
    t = array([float(x) for x in t_cell.split(",")])
    return u, t


def score_overlap(aln, overlap_positions, similarity_matrix):
    """
    Sequence identity and the three BLOSUM62 similarity scores over the overlapping residues.

    Args:
        aln (AlignmentRow): The row, read for `qaln` and `taln`.
        overlap_positions (list): Alignment positions shared by both pockets.
        similarity_matrix (dict): BLOSUM62, as returned by `lib.read_blast_similarity_matrix`.

    Returns:
        dict: The two overlap sequences plus the identity and similarity columns.
    """
    p1_aln_seq = "".join([aln.qaln[pos] for pos in overlap_positions])
    p2_aln_seq = "".join([aln.taln[pos] for pos in overlap_positions])

    similarity_1_2 = full_similarity(p1_aln_seq, p2_aln_seq, similarity_matrix)
    similarity_2_1 = full_similarity(p2_aln_seq, p1_aln_seq, similarity_matrix)

    return {
        "query_seq_overlap": p1_aln_seq,
        "target_seq_overlap": p2_aln_seq,
        "overlap_identity": sum(map(str.__eq__, p1_aln_seq, p2_aln_seq)) / len(overlap_positions),
        "overlap_similarity_binary": binary_similarity(p1_aln_seq, p2_aln_seq, similarity_matrix),
        "overlap_similarity_1_2": similarity_1_2,
        "overlap_similarity_2_1": similarity_2_1,
        "min_overlap_similarity": min(similarity_1_2, similarity_2_1),
        "max_overlap_similarity": max(similarity_1_2, similarity_2_1),
    }


def compare_pocket_pair(aln, pocket_id_1, p1, p1_mapped, pocket_id_2, p2, p2_mapped, ctx):
    """
    Score one pocket against one other on a single alignment row.

    Args:
        aln (AlignmentRow): The row bridging the two pockets.
        pocket_id_1 (str): Query pocket id.
        p1 (Pocket): Query pocket, with `pocket_exists` set, as must be `p2`.
        p1_mapped (MappedPocket): Its projection onto this row.
        pocket_id_2 (str): Target pocket id.
        p2 (Pocket): Target pocket.
        p2_mapped (MappedPocket): Its projection onto this row.
        ctx (Context): Scoring state shared across the call.

    Returns:
        dict: One comparison row, holding only the descriptor columns when the pockets do not overlap.
            The target_* descriptor columns and `jaccard_index` are omitted when p2 is a whole chain
            rather than a pocket on one.
    """
    output = {
        "query": pocket_id_1,
        "target": pocket_id_2,
        "evalue": aln.evalue,
        "lddt": aln.lddt,
    }

    (
        output["query_res_ids"],
        output["query_len"],
        output["query_seq"],
    ) = describe_pocket(pocket_id_1, p1, ctx.descriptions)
    output["query_pct_aln"] = p1_mapped.in_aln_count / output["query_len"]

    # A whole chain has no pocket to describe, and its length would swamp these columns and the
    # Jaccard union. Flagged per pocket, so one run can mix open and pocketed targets.
    p2_is_whole_chain = p2.whole_chain
    if not p2_is_whole_chain:
        (
            output["target_res_ids"],
            output["target_len"],
            output["target_seq"],
        ) = describe_pocket(pocket_id_2, p2, ctx.descriptions)
        output["target_pct_aln"] = p2_mapped.in_aln_count / output["target_len"]

    # Kept in pocket-1 order: the overlap sequences are built by indexing the alignment strings with it.
    p2_positions = set(p2_mapped.positions)
    overlap_positions = [pos for pos in p1_mapped.positions if pos in p2_positions]
    output["overlap_count"] = len(overlap_positions)
    if not overlap_positions:
        return output

    overlap_set = set(overlap_positions)
    p1_overlap_ids = overlap_ids(p1, p1_mapped, overlap_set)
    p2_overlap_ids = overlap_ids(p2, p2_mapped, overlap_set)
    output["query_overlap_ids"] = ",".join(p1_overlap_ids)
    output["target_overlap_ids"] = ",".join(p2_overlap_ids)

    if not p2_is_whole_chain:
        union_size = len(p1.res_auth_ids) + len(p2.res_auth_ids) - len(overlap_positions)
        output["jaccard_index"] = len(overlap_positions) / union_size

    output.update(score_overlap(aln, overlap_positions, ctx.similarity_matrix))
    output.update(superpose(p1, p2, p1_overlap_ids, p2_overlap_ids, len(overlap_positions), ctx.superimposer))
    return output


class Context(NamedTuple):
    """
    The scoring state shared by every comparison in one call.

    `offsets` is the database's offset table (entry name -> region spec), empty when the target is not
    a Foldseek database or ships no table.
    """

    similarity_matrix: dict
    superimposer: SVDSuperimposer
    descriptions: dict
    identities: dict
    offsets: dict


def resolve_pockets(domain, pocket_dict, preproc_to_ids):
    """
    The pockets sitting on one aligned chain. One chain can carry several pockets.

    Args:
        domain (str): The chain's `preprocess_name`.
        pocket_dict (dict): pocket_id -> Pocket, or None for one not built, for every pocket in the run.
        preproc_to_ids (dict): preprocess_name -> the pocket_ids on that chain.

    Returns:
        dict: pocket_id -> Pocket, for this chain's built pockets only.
    """
    return {
        pocket_id: pocket_dict[pocket_id]
        for pocket_id in preproc_to_ids.get(domain) or []
        if pocket_dict.get(pocket_id) is not None
    }


def compare_pockets(
    alignment_df,
    pocket_dict,
    preproc_to_ids,
    blosum_path,
    synthesise_target_pockets=False,
    offset_table_path=None,
):
    """
    Compare every pair of pockets on the chains each alignment row bridges.

    Args:
        alignment_df (pandas.DataFrame): The alignment table, columns in `ALIGNMENT_COLUMNS` order.
        pocket_dict (dict): pocket_id -> Pocket, or None for one not built.
        preproc_to_ids (dict): preprocess_name -> the pocket_ids sitting on that chain; bridges the
            alignment's keys to the pockets'.
        blosum_path (str): Path to a BLAST-format similarity matrix. The packaged one is
            os.path.join(os.path.dirname(pocketmapper.__file__), "blosum62.bla").
        synthesise_target_pockets (bool): Build a whole-chain pseudo-pocket per alignment row instead
            of looking the target up in `pocket_dict`, for a Foldseek database with no target records.
        offset_table_path (str | None): Path to the offset table shipped with that database, which
            renumbers the synthesised targets' residues into UniProt coordinates. Consulted only
            alongside `synthesise_target_pockets`; None leaves them as positions within the entry.

    Returns:
        tuple: (comparison table, the residue codes the aligner and pocketmapper disagreed on, the
            pockets whose sequence did not match the aligner's well enough to be trusted).

    Raises:
        PocketMapperError: If the offset table does not describe a hit. Any other error comparing a
            row is logged and re-raised.
    """
    ctx = Context(
        similarity_matrix=read_blast_similarity_matrix(blosum_path),
        superimposer=SVDSuperimposer(),
        descriptions={},
        identities={},
        offsets=read_offset_table(offset_table_path) if offset_table_path else {},
    )

    log_extra = {"stage": "Pocket Comparison"}
    unknown_ids = defaultdict(lambda: defaultdict(set))  # for saving tri-code ids which are unknown
    incorrect_mapping = defaultdict(dict)  # for saving cases where foldseek mapping doesn't match pocketmapper sequence

    existing_calcs = set()
    output_rows = []

    for values in tqdm(alignment_df.itertuples(index=False, name=None)):
        aln = AlignmentRow(*values)
        # Foldseek lowercases residues it masked for seeding (CA B-factor below --mask-bfactor-threshold)
        aln = aln._replace(qaln=aln.qaln.upper(), taln=aln.taln.upper(), qseq=aln.qseq.upper(), tseq=aln.tseq.upper())
        try:
            pockets_1 = resolve_pockets(aln.query, pocket_dict, preproc_to_ids)
            if not pockets_1:
                continue

            if synthesise_target_pockets:
                pockets_2 = {aln.target: synthesise_target_pocket(aln, ctx)}
            else:
                pockets_2 = resolve_pockets(aln.target, pocket_dict, preproc_to_ids)
                if not pockets_2:
                    continue

            q_positions = aln_positions(aln.qaln)
            t_positions = aln_positions(aln.taln)

            # Each side is projected onto this row once, then reused across every pairing.
            mapped_1 = {}
            mapped_2 = {}

            for pocket_id_1, pocket_id_2 in product(pockets_1, pockets_2):
                # A pair already scored on an earlier row is skipped
                if (pocket_id_1, pocket_id_2) in existing_calcs:
                    continue
                existing_calcs.add((pocket_id_1, pocket_id_2))

                if pocket_id_1 in incorrect_mapping or pocket_id_2 in incorrect_mapping:
                    continue

                p1 = pockets_1[pocket_id_1]
                p2 = pockets_2[pocket_id_2]
                if not p1.pocket_exists or not p2.pocket_exists:
                    continue

                p1_identity = seq_identity(pocket_id_1, p1, aln.query, aln.qseq, ctx.identities)
                if p1_identity < MIN_SEQ_IDENTITY:
                    incorrect_mapping[pocket_id_1] = {
                        "p1_seq_identity": p1_identity,
                        "p1_seq": p1.ca_sequence,
                        "fs_seq": aln.qseq,
                    }

                p2_identity = seq_identity(pocket_id_2, p2, aln.target, aln.tseq, ctx.identities)
                if p2_identity < MIN_SEQ_IDENTITY:
                    incorrect_mapping[pocket_id_2] = {
                        "p2_seq_identity": p2_identity,
                        "p2_seq": p2.ca_sequence,
                        "fs_seq": aln.tseq,
                    }

                # Both sides are checked first so each failing pocket is recorded, then the pair is dropped
                if pocket_id_1 in incorrect_mapping or pocket_id_2 in incorrect_mapping:
                    continue

                if pocket_id_1 not in mapped_1:
                    mapped_1[pocket_id_1] = map_pocket_into_alignment(p1, aln.qaln, q_positions, aln.qstart, aln.qend)
                if pocket_id_2 not in mapped_2:
                    mapped_2[pocket_id_2] = map_pocket_into_alignment(p2, aln.taln, t_positions, aln.tstart, aln.tend)
                p1_mapped = mapped_1[pocket_id_1]
                p2_mapped = mapped_2[pocket_id_2]

                # The projections are shared across pairings, but unknown_ids names both pockets, so
                # each pairing contributes its own entries.
                record_code_mismatches(p1_mapped, pocket_id_1, pocket_id_2, unknown_ids)
                record_code_mismatches(p2_mapped, pocket_id_2, pocket_id_1, unknown_ids)

                output_rows.append(
                    compare_pocket_pair(aln, pocket_id_1, p1, p1_mapped, pocket_id_2, p2, p2_mapped, ctx)
                )

        except Exception:
            logger.exception(
                f"Uncontrolled error calculating {aln.query} and {aln.target}",
                extra=log_extra,
            )
            raise

    # Reindex onto the declared schema so every row carries every column. An undeclared column is kept
    # and flagged rather than silently dropped.
    pockets_df = pd.DataFrame.from_dict(output_rows)
    undeclared = [column for column in pockets_df.columns if column not in POCKET_COMPARISON_COLUMNS]
    if undeclared:
        logger.warning(
            f"Pocket comparison produced undeclared column(s) {undeclared}; add them to POCKET_COMPARISON_COLUMNS",
            extra=log_extra,
        )
    pockets_df = pockets_df.reindex(columns=POCKET_COMPARISON_COLUMNS + undeclared)

    return pockets_df, unknown_ids, incorrect_mapping
