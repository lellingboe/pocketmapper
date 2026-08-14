import hashlib
import logging
from copy import deepcopy
from Bio.SVDSuperimposer import SVDSuperimposer
from itertools import product
from tqdm import tqdm
from collections import defaultdict
import pandas as pd
from numpy import array
from numpy import linalg as LA

# Every column compare_pockets can produce, in output order.
#
# A comparison stops early whenever there is nothing further to compute -- no overlapping residues, no CA
# coordinates, fewer than three overlapping residues to superpose, or a foldseek-db/alphafold target with no
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
    "pocket_1_pct_overlap",
    "pocket_2_pct_overlap",
    "min_pct_overlap",
    "max_pct_overlap",
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


def jsonify_dict(item):
    """
    Recursively looks for sets in a dictionary and turns then into lists
    This allows dicts with sets to become JSON serializeable
    """
    if isinstance(item, set):
        return list(item)
    elif isinstance(item, dict):
        return {str(k): jsonify_dict(v) for k, v in item.items()}
    else:
        return item


def safe_filename(name, max_len=80):
    """
    Build a filesystem-safe filename stem from a pocket_id.

    pocket_ids can embed long comma-separated residue lists (e.g. for
    passthrough/VDW queries), which can exceed OS filename length limits.
    Names longer than max_len are truncated and given an md5 suffix to keep
    them unique.
    """
    safe_name = name.replace(":", "_").replace(",", "_")
    if len(safe_name) <= max_len:
        return safe_name
    name_hash = hashlib.md5(name.encode()).hexdigest()
    return f"{safe_name[:max_len]}_{name_hash}"


"""
Foldseek output format used for comparison:
0 query
1 target
2 fident
3 alnlen
4 mismatch
5 gapopen
6 qstart
7 qend
8 tstart
9 tend
10 evalue
11 lddt
12 qaln
13 taln
14 u
15 t
16 qseq
17 tseq
"""


def compare_pockets(
    alignment_df,
    pocket_dict,
    preproc_to_ids,
    blosum_path=r"/home/data/motif_aligner/blosum62.bla",
    alphafold=False,
):
    """
    Compare two pockets based on foldseek alignment
    """

    blosum_similarity_matrix = read_blast_similarity_matrix(blosum_path)

    unknown_ids = defaultdict(lambda: defaultdict(set))  # for saving tri-code ids which are unknown
    incorrect_mapping = defaultdict(dict)  # for saving cases where foldseek mapping doesn't match pocketmapper sequence

    # Setting up vars for use later
    existing_calcs = set()
    output_rows = []
    sup = SVDSuperimposer()

    # TODO divide this into common things for each pocket and a cross-comparison
    for row in tqdm(alignment_df.itertuples(index=False)):
        domain_1 = row[0]
        domain_2 = row[1]
        try:

            # MAPPING SEQUENCE POSITION TO ALIGNMENT POSITION #####
            p1_seq_to_aln = {}
            i = 0
            for j, res in enumerate(row[12]):
                if res != "-":
                    p1_seq_to_aln[i] = j
                    i += 1

            p2_seq_to_aln = {}
            i = 0
            for j, res in enumerate(row[13]):
                if res != "-":
                    p2_seq_to_aln[i] = j
                    i += 1

            # GETTING POCKETS WHICH CORRESPOND TO THE FOLDSEEK NAME
            pockets_1 = {}
            for pocket_id in preproc_to_ids.get(domain_1):
                if pocket_id in pocket_dict:
                    pockets_1[pocket_id] = pocket_dict[pocket_id]
            if len(pockets_1) == 0:
                continue

            if alphafold:
                pockets_2 = {
                    domain_2: {
                        "res_auth_ids": [str(k) for k in range(row[9])],
                        "id_pos_codes_match": True,
                        "pocket_exists": True,
                        "has_coords": False,
                        "ca_sequence": row[17],
                    }
                }
                pockets_2[domain_2].update({str(k): {"seq_pos": k} for k in range(row[9])})
            else:
                pockets_2 = {}
                for pocket_id in preproc_to_ids.get(domain_2):
                    if pocket_id in pocket_dict:
                        pockets_2[pocket_id] = pocket_dict[pocket_id]
                if len(pockets_2) == 0:
                    continue

            # Iterating through aligned pairs
            for pocket_id_1, pocket_id_2 in product(pockets_1.keys(), pockets_2.keys()):

                # Checking for A-B comparison if B-A has already been calculated
                if (pocket_id_1, pocket_id_2) in existing_calcs:
                    continue
                existing_calcs.add((pocket_id_1, pocket_id_2))

                if pocket_id_1 in incorrect_mapping or pocket_id_2 in incorrect_mapping:
                    continue

                # Starting the output list
                output = {
                    "pocket_1": pocket_id_1,
                    "pocket_2": pocket_id_2,
                    "evalue": row[10],
                    "lddt": row[11],
                }

                # Getting the pockets
                p1 = deepcopy(pockets_1.get(pocket_id_1))
                p2 = deepcopy(pockets_2.get(pocket_id_2))
                if not p1["pocket_exists"] or not p2["pocket_exists"]:
                    continue

                p1_seq_identity = sum(map(str.__eq__, row[16], p1["ca_sequence"])) / len(p1["ca_sequence"])
                if p1_seq_identity < 0.8:
                    incorrect_mapping[pocket_id_1] = {
                        "p1_seq_identity": p1_seq_identity,
                        "p1_seq": p1["ca_sequence"],
                        "fs_seq": row[16],
                    }

                p2_seq_identity = sum(map(str.__eq__, row[17], p2["ca_sequence"])) / len(p2["ca_sequence"])
                if p2_seq_identity < 0.8:
                    incorrect_mapping[pocket_id_2] = {
                        "p2_seq_identity": p2_seq_identity,
                        "p2_seq": p2["ca_sequence"],
                        "fs_seq": row[17],
                    }

                # Calculated Metrics
                # pocket 1
                p1_adj = 1 - row[6]  # row[6] is qstart
                p1_fs_len = row[7] - row[6] + 1  # row[7] is qend
                p1_in_aln_region_count = 0
                p1["foldseek_pos"] = []
                for res in p1["res_auth_ids"]:
                    fs_start_adj_pos = int(p1[res]["seq_pos"]) + p1_adj
                    in_aln_region = -1 < fs_start_adj_pos < p1_fs_len
                    p1[res]["in_fs_aln_region"] = in_aln_region
                    if in_aln_region:
                        p1_in_aln_region_count += 1
                        fs_res_pos = p1_seq_to_aln[fs_start_adj_pos]
                        p1["foldseek_pos"].append(fs_res_pos)
                        p1[res]["fs_pos"] = fs_res_pos
                        p1[res]["fs_res_code"] = row[12][fs_res_pos]

                        # Checking fs single res code and pocektmapper single res codes match
                        if p1[res]["fs_res_code"] != p1[res]["res_code_single"]:
                            unknown_ids[p1[res]["fs_res_code"]][p1[res]["res_code"]].add(
                                ",".join([pocket_id_2, pocket_id_1, res])
                            )

                output["pocket_1_res_ids"] = ",".join(p1["res_auth_ids"])
                output["pocket_1_len"] = len(p1["res_auth_ids"])
                output["pocket_1_seq"] = "".join([p1[x]["res_code_single"] for x in p1["res_auth_ids"]])
                output["pocket_1_pct_aln"] = p1_in_aln_region_count / len(p1["res_auth_ids"])

                # pocket 2
                p2_adj = 1 - row[8]
                p2_fs_len = row[9] - row[8] + 1
                p2_in_aln_region_count = 0
                p2["foldseek_pos"] = []
                for res in p2["res_auth_ids"]:
                    fs_start_adj_pos = int(p2[res]["seq_pos"]) + p2_adj
                    in_aln_region = -1 < fs_start_adj_pos < p2_fs_len
                    p2[res]["in_fs_aln_region"] = in_aln_region
                    if in_aln_region:
                        p2_in_aln_region_count += 1
                        fs_res_pos = p2_seq_to_aln[fs_start_adj_pos]
                        p2["foldseek_pos"].append(fs_res_pos)
                        p2[res]["fs_pos"] = fs_res_pos
                        p2[res]["fs_res_code"] = row[13][fs_res_pos]

                        # Checking fs single res code and pocektmapper single res codes match

                        if not alphafold and (p2[res]["fs_res_code"] != p2[res]["res_code_single"]):
                            debug_id = ",".join([pocket_id_1, pocket_id_2, res])
                            unknown_ids[p2[res]["fs_res_code"]][p2[res]["res_code"]].add(debug_id)

                if not alphafold:
                    output["pocket_2_res_ids"] = ",".join(p2["res_auth_ids"])
                    output["pocket_2_len"] = len(p2["res_auth_ids"])
                    output["pocket_2_seq"] = "".join([p2[x]["res_code_single"] for x in p2["res_auth_ids"]])
                    output["pocket_2_pct_aln"] = p2_in_aln_region_count / len(p2["res_auth_ids"])

                # OVERLAP
                overlapping_residues = [x for x in p1["foldseek_pos"] if x in p2["foldseek_pos"]]
                output["overlap_count"] = len(overlapping_residues)
                if len(overlapping_residues) == 0:  # If no overlapping resides, end here
                    output_rows.append(output)
                    continue

                # overlap ids
                p1_overlap_ids = []
                for res in p1["res_auth_ids"]:
                    if p1[res].get("fs_pos", -1) in overlapping_residues:
                        p1_overlap_ids.append(res)
                output["pocket_1_overlap_ids"] = ",".join(p1_overlap_ids)
                p2_overlap_ids = []
                for res in p2["res_auth_ids"]:
                    if p2[res].get("fs_pos", -1) in overlapping_residues:
                        p2_overlap_ids.append(res)
                output["pocket_2_overlap_ids"] = ",".join(p2_overlap_ids)

                # percent overlap
                p1_pct_overlap = len(overlapping_residues) / len(p1["res_auth_ids"])
                output["pocket_1_pct_overlap"] = p1_pct_overlap

                if not alphafold:
                    p2_pct_overlap = len(overlapping_residues) / len(p2["res_auth_ids"])
                    output["pocket_2_pct_overlap"] = p2_pct_overlap
                    output["min_pct_overlap"] = min(p1_pct_overlap, p2_pct_overlap)
                    output["max_pct_overlap"] = max(p1_pct_overlap, p2_pct_overlap)

                #####
                # Aligned residues as a sequence
                p1_aln_seq = "".join([row[12][x] for x in overlapping_residues])
                p2_aln_seq = "".join([row[13][x] for x in overlapping_residues])
                output["pocket_1_seq_overlap"] = p1_aln_seq
                output["pocket_2_seq_overlap"] = p2_aln_seq

                # Identity
                overlap_identity = sum(map(str.__eq__, p1_aln_seq, p2_aln_seq)) / len(overlapping_residues)
                output["overlap_identity"] = overlap_identity

                # BLOSUM62 similarity
                similarity = binary_similarity(p1_aln_seq, p2_aln_seq, blosum_similarity_matrix)
                output["overlap_similarity_binary"] = similarity

                similarity_1_2 = full_similarity(p1_aln_seq, p2_aln_seq, blosum_similarity_matrix)
                similarity_2_1 = full_similarity(p2_aln_seq, p1_aln_seq, blosum_similarity_matrix)
                output["overlap_similarity_1_2"] = similarity_1_2
                output["overlap_similarity_2_1"] = similarity_2_1
                output["min_overlap_similarity"] = min(similarity_1_2, similarity_2_1)
                output["max_overlap_similarity"] = max(similarity_1_2, similarity_2_1)

                # RMSD dings
                if p1["has_coords"] and p2["has_coords"]:
                    x = array([p1[str(x)]["ca_coords"] for x in p1_overlap_ids])
                    y = array([p2[str(x)]["ca_coords"] for x in p2_overlap_ids])

                    if len(overlapping_residues) > 2:
                        sup.set(x, y)
                        sup.run()
                        u, t = sup.get_rotran()
                        output["p2_to_p1_u"] = u.flatten().tolist()
                        output["p2_to_p1_t"] = t.tolist()

                        # TODO do this with matrix algebra isntead of doing it twice
                        sup.set(y, x)
                        sup.run()
                        u, t = sup.get_rotran()
                        output["p1_to_p2_u"] = u.flatten().tolist()
                        output["p1_to_p2_t"] = t.tolist()
                        output["rmsd"] = sup.get_rms()

                        x_on_y = sup.get_transformed()
                        ca_dists = LA.norm(x_on_y - y, axis=1)
                        output["ca_dists"] = ",".join([str(round(x, 3)) for x in ca_dists])

                output_rows.append(output)

        except KeyError:
            logging.warning(
                f"Uncontrolled KeyError calculating {domain_1} and {domain_2}",
                extra={"stage": "Pocket Comparison"},
            )
            raise
        except Exception:
            logging.exception(
                f"Uncontrolled error calculating {domain_1} and {domain_2}",
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


def binary_similarity(seqA, seqB, similarity_matrix):
    """
    Similarity of A->B and B->A
    score = 1 if full similarity score > 0 else score = 0
    normalized to the length of the seqence
    """
    seqA = seqA.replace("U", "X").upper()
    seqB = seqB.replace("U", "X").upper()

    similarity = [(similarity_matrix[x][y] > 0) for x, y in zip(seqA, seqB)]  # True or False
    similarity_score = sum(similarity) / len(similarity)
    return similarity_score


def full_similarity(seqA, seqB, similarity_matrix):
    """
    Similarity of A->B
    blosum scores normalized by
    normalized to the length of the seqence
    """
    seqA = seqA.replace("U", "X").upper()
    seqB = seqB.replace("U", "X").upper()

    similarity = [similarity_matrix[x][y] for x, y in zip(seqA, seqB)]
    similarity_max = [similarity_matrix[x][x] for x in seqA]
    similarity_normalized = [x / y for x, y in zip(similarity, similarity_max)]
    return sum(similarity_normalized) / len(similarity_normalized)


def read_blast_similarity_matrix(similarity_matrix_path, delimiter=" "):
    similarity_matrix = {}
    file_content = open(similarity_matrix_path).read().strip().split("\n")
    header = None
    idx_to_aa = None
    max_score = 0
    row_counter = 0
    for line in file_content:
        # Skip comment lines
        if line[0] == "#":
            continue

        # parsing the first lines into a dict
        if header is None:
            if delimiter == " ":
                header = line.strip().split()  # splits on whitespace and discards empty results
            else:
                header = line.strip().split(delimiter)
            idx_to_aa = dict(list(zip(list(range(0, len(header))), header)))
        else:
            if delimiter == " ":
                fields = line.strip().split()
            else:
                fields = line.strip().split(delimiter)
            from_aa = idx_to_aa[row_counter]
            row_counter += 1
            similarity_matrix[from_aa] = {}
            for idx, score in enumerate(fields):
                to_aa = idx_to_aa[idx]
                similarity_matrix[from_aa][to_aa] = float(score)
                if to_aa not in similarity_matrix:
                    similarity_matrix[to_aa] = {}
                similarity_matrix[to_aa][from_aa] = float(score)
                if float(score) > max_score:
                    max_score = float(score)
            similarity_matrix[from_aa]["-"] = -4

    # Add and construct the "-" entry
    similarity_matrix["-"] = {}
    for aa in similarity_matrix:
        if aa == "-":
            similarity_matrix["-"][aa] = -1
        else:
            similarity_matrix["-"][aa] = -4
    return similarity_matrix
