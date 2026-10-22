"""
Substitution-matrix scoring of aligned pocket sequences.
"""


def read_blast_similarity_matrix(similarity_matrix_path, delimiter=" "):
    """
    Read a BLAST-format substitution matrix into a nested dict.

    Args:
        similarity_matrix_path (str): Path to the matrix file. Comment lines are skipped and the first
            remaining line is the residue header, which also fixes the row order.
        delimiter (str): Column separator. The default " " splits on any whitespace run.

    Returns:
        dict: Nested residue -> residue -> float score, stored symmetrically. Adds a "-" gap row and
            column the file does not carry: -4 against any residue, -1 against itself.
    """
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


def binary_similarity(seqA, seqB, similarity_matrix):
    """
    Fraction of positions where two aligned sequences score above zero.

    Args:
        seqA (str): First sequence; must be the same length as seqB.
        seqB (str): Second sequence.
        similarity_matrix (dict): Nested residue -> residue -> score.

    Returns:
        float: Score in 0..1: how much of the sequence is conservatively substituted, not how strongly.
    """
    seqA = seqA.replace("U", "X").upper()
    seqB = seqB.replace("U", "X").upper()

    similarity = [(similarity_matrix[x][y] > 0) for x, y in zip(seqA, seqB)]  # True or False
    similarity_score = sum(similarity) / len(similarity)
    return similarity_score


def full_similarity(seqA, seqB, similarity_matrix):
    """
    Substitution score of two aligned sequences, normalised per position.

    Args:
        seqA (str): First sequence; must be the same length as seqB.
        seqB (str): Second sequence.
        similarity_matrix (dict): Nested residue -> residue -> score.

    Returns:
        float: Mean over positions of the score divided by seqA's residue scored against itself, so a
            perfect match scores 1 however strongly that residue is conserved.
    """
    seqA = seqA.replace("U", "X").upper()
    seqB = seqB.replace("U", "X").upper()

    similarity = [similarity_matrix[x][y] for x, y in zip(seqA, seqB)]
    similarity_max = [similarity_matrix[x][x] for x in seqA]
    similarity_normalized = [x / y for x, y in zip(similarity, similarity_max)]
    return sum(similarity_normalized) / len(similarity_normalized)
