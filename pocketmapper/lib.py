"""
Generic, stateless helpers shared across PocketMapper.

Nothing here knows about the pipeline, Settings, or the Pocket shape -- each function takes
plain values and returns plain values. Workflow logic belongs in the component modules
(pocket_comparison, sequence_aligner, ...) rather than here.
"""

import hashlib
import os
import re

_UNSAFE_FILENAME_CHARS = re.compile(r"[^A-Za-z0-9._-]")


def jsonify_dict(item):
    """
    Recursively turn sets into lists so a dict becomes JSON-serialisable.

    Args:
        item: Any value; dicts are walked and their keys coerced to strings.

    Returns:
        The same structure with every set replaced by a list.
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

    A pocket_id is a raw input string, so it may be a path ("/data/foo.cif.gz:B_F") and may embed a
    long comma-separated residue list (passthrough/VDW queries). Only the basename is kept -- a stem
    with a directory component in it resolves outside the directory the caller joins it onto -- and
    every remaining character outside [A-Za-z0-9._-] becomes "_".

    Both of those steps are lossy, and truncation to max_len is lossy again, so an md5 of the full
    original name is always appended: without it two distinct pocket_ids can reduce to one filename
    and silently overwrite each other. A consequence is that the resulting files are not greppable
    for the input string -- match on the `MOLECULE` records inside them instead.

    Args:
        name (str): The pocket_id, or any raw input string.
        max_len (int): Ceiling on the result's length. Defaults to 80.

    Returns:
        str: A safe stem, at most max_len characters (or 32, if max_len leaves no room).
    """
    name_hash = hashlib.md5(name.encode()).hexdigest()
    stem = _UNSAFE_FILENAME_CHARS.sub("_", os.path.basename(name.rstrip("/")))
    stem = stem[: max(max_len - len(name_hash) - 1, 0)].rstrip("_")
    return f"{stem}_{name_hash}" if stem else name_hash


def binary_similarity(seqA, seqB, similarity_matrix):
    """
    Fraction of positions where two aligned sequences score above zero.

    Each position scores 1 if its substitution score is positive and 0 otherwise, so this measures
    how much of the sequence is conservatively substituted rather than how strongly.

    Args:
        seqA (str): First sequence; must be the same length as seqB.
        seqB (str): Second sequence.
        similarity_matrix (dict): Nested residue -> residue -> score.

    Returns:
        float: Score in 0..1, normalised by sequence length.
    """
    seqA = seqA.replace("U", "X").upper()
    seqB = seqB.replace("U", "X").upper()

    similarity = [(similarity_matrix[x][y] > 0) for x, y in zip(seqA, seqB)]  # True or False
    similarity_score = sum(similarity) / len(similarity)
    return similarity_score


def full_similarity(seqA, seqB, similarity_matrix):
    """
    Substitution score of two aligned sequences, normalised per position.

    Each position's score is divided by the score of seqA's residue against itself, so a perfect
    match scores 1 regardless of how strongly that residue is conserved.

    Args:
        seqA (str): First sequence; must be the same length as seqB.
        seqB (str): Second sequence.
        similarity_matrix (dict): Nested residue -> residue -> score.

    Returns:
        float: Mean normalised score over the sequence.
    """
    seqA = seqA.replace("U", "X").upper()
    seqB = seqB.replace("U", "X").upper()

    similarity = [similarity_matrix[x][y] for x, y in zip(seqA, seqB)]
    similarity_max = [similarity_matrix[x][x] for x in seqA]
    similarity_normalized = [x / y for x, y in zip(similarity, similarity_max)]
    return sum(similarity_normalized) / len(similarity_normalized)


def read_blast_similarity_matrix(similarity_matrix_path, delimiter=" "):
    """
    Read a BLAST-format substitution matrix into a nested dict.

    Comment lines are skipped and the first remaining line is taken as the residue header, which also
    fixes the row order. Scores are stored symmetrically, so either lookup order works.

    A gap row and column are added that the file does not carry: "-" against any residue scores -4
    and against itself -1, so a caller can score an alignment containing gaps without special-casing
    them.

    Args:
        similarity_matrix_path (str): Path to the matrix file.
        delimiter (str): Column separator. The default " " splits on any whitespace run.

    Returns:
        dict: Nested residue -> residue -> float score, including the "-" gap entries.
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


# Foldseek's PDB database names its entries "<pdbid>-assembly<N>_<chain>", with a "-<copy>" suffix on
# chains duplicated within an assembly (e.g. "5ian-assembly1_B-2").
_FOLDSEEK_PDB_ENTRY = re.compile(r"^(?P<pdb>[0-9A-Za-z]{4})-assembly(?P<assembly>\d+)_(?P<chain>.+)$")


def parse_foldseek_pdb_entry_name(name):
    """
    Resolve a Foldseek PDB-database entry name into the PDB ID and chain it came from.

    The assembly number and any chain-copy suffix are discarded, so "5ian-assembly1_B-2" and
    "5ian-assembly2_B" both resolve to ("5IAN", "B"). The PDB ID is upper-cased to match the form
    QTProcessor derives from user input, so a structure already fetched for a query is reused rather
    than downloaded a second time.

    Args:
        name (str): A Foldseek database entry name.

    Returns:
        tuple: (pdb_id, chain), or None for anything that is not a PDB-style entry name -- which is
            also how a caller tells a PDB Foldseek database apart from one built on something else
            (e.g. human_domains).
    """
    match = _FOLDSEEK_PDB_ENTRY.match(name)
    if match is None:
        return None
    return match.group("pdb").upper(), match.group("chain").split("-")[0]
