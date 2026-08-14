"""
Generic, stateless helpers shared across PocketMapper.

Nothing here knows about the pipeline, Settings, or the pocket dict shape -- each function takes
plain values and returns plain values. Workflow logic belongs in the component modules
(pocket_comparison, sequence_aligner, ...) rather than here.
"""

import hashlib


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
