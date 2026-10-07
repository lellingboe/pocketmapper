"""
Generic, stateless helpers shared across PocketMapper.

Nothing here knows about the pipeline, Settings, or the Pocket shape -- each function takes
plain values and returns plain values, or does one thing to the filesystem or the package logger.
Workflow logic belongs in the component modules rather than here.

The one piece of state is `HELD_TEMP_DIRS`, the scratch directories an open `temp_dir_scope` holds.
"""

import gzip
import hashlib
import logging
import os
import re
import shutil
from contextlib import contextmanager

from pocketmapper.constants import FOLDSEEK_AA_CODES
from pocketmapper.constants import LOG_FORMAT
from pocketmapper.constants import PACKAGE_LOGGER
from pocketmapper.exceptions import PocketMapperError

logger = logging.getLogger(__name__)

UNSAFE_FILENAME_CHARS = re.compile(r"[^A-Za-z0-9._-]")

# Absolute paths of the scratch directories an open temp_dir_scope holds
HELD_TEMP_DIRS = set()


class StageFilter(logging.Filter):
    """
    Supply a missing `stage` attribute from the name of the function that logged the record.

    `LOG_FORMAT` interpolates `%(stage)s`, which is not a stock LogRecord attribute, so without this
    filter a record logged with no `extra={"stage": ...}` fails to format.
    """

    def filter(self, record):
        """
        Default `record.stage` to `record.funcName` and keep the record.

        Args:
            record (logging.LogRecord): The record about to be emitted; modified in place.

        Returns:
            bool: Always True -- nothing is filtered out.
        """
        if not hasattr(record, "stage"):
            record.stage = record.funcName
        return True


def format_handler(handler):
    """
    Give a handler the package log format and the `StageFilter` that format needs.

    Args:
        handler (logging.Handler): The handler to configure; modified in place.

    Returns:
        logging.Handler: The same handler.
    """
    handler.setFormatter(logging.Formatter(LOG_FORMAT))
    handler.addFilter(StageFilter())
    return handler


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

    Args:
        name (str): The pocket_id, or any raw input string, which may be a path.
        max_len (int): Ceiling on the result's length. Defaults to 80.

    Returns:
        str: The basename of `name` with every character outside [A-Za-z0-9._-] replaced by "_",
            truncated, and suffixed with an md5 of the full `name` so distinct names never collide.
            At most max_len characters (or 32, if max_len leaves no room).
    """
    # Basename only: a stem with a directory in it would resolve outside the directory it is joined onto
    name_hash = hashlib.md5(name.encode()).hexdigest()
    stem = UNSAFE_FILENAME_CHARS.sub("_", os.path.basename(name.rstrip("/")))
    stem = stem[: max(max_len - len(name_hash) - 1, 0)].rstrip("_")
    return f"{stem}_{name_hash}" if stem else name_hash


def gzip_file(src_fpath, dst_fpath):
    """
    Write a gzip-compressed copy of a file.

    Args:
        src_fpath (str): File to read.
        dst_fpath (str): Path to write the compressed copy to; overwritten if it already exists.

    Returns:
        None
    """
    with open(src_fpath, "rb") as f_in:
        with gzip.open(dst_fpath, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def is_within(path, roots):
    """
    Report whether a path resolves to somewhere inside one of `roots`.

    Args:
        path (str): The path to test.
        roots (list): Candidate containing directories; only one has to match.

    Returns:
        bool: True if `path` is inside (or equal to) any of `roots`, both sides resolved with
            `os.path.realpath` so `..` segments and symlinks cannot walk out of a root.
    """
    real_path = os.path.realpath(path)
    for root in roots:
        real_root = os.path.realpath(root)
        # commonpath raises on a mix of absolute and relative paths; realpath has made both
        # absolute, so the only remaining raiser is a different drive, which cannot be a match.
        try:
            if os.path.commonpath([real_path, real_root]) == real_root:
                return True
        except ValueError:
            continue
    return False


def make_dir(path, log_extra):
    """
    Create a directory and any missing parents, if it does not already exist.

    Args:
        path (str): The directory. An empty string, the directory of a bare filename, is the
            working directory and is left alone.
        log_extra (dict): Logging `extra` for the failure message.

    Returns:
        None

    Raises:
        PocketMapperError: If the directory cannot be created.
    """
    if not path:
        return
    try:
        os.makedirs(path, exist_ok=True)
    except OSError as e:
        logger.critical(f"Error creating directory {path}", extra=log_extra)
        raise PocketMapperError(f"Error creating directory {path}") from e


def empty_temp_dir(temp_dir, roots):
    """
    Delete a scratch directory so that whatever uses it next starts with it empty.

    Only a `temp_dir` inside one of `roots` is deleted, so a mistyped path costs a stray directory,
    not its contents. One outside them is reused as it is, with a warning if it exists.

    Args:
        temp_dir (str): The scratch directory. Not created here.
        roots (list): The directories it may be deleted under.

    Returns:
        None
    """
    log_extra = {"stage": "Configuring Workflow"}

    if is_within(temp_dir, roots):
        shutil.rmtree(temp_dir, ignore_errors=True)
    elif os.path.exists(temp_dir):
        logger.warning(
            f"Reusing temp_dir {temp_dir} without emptying it: it is outside {' and '.join(roots)}. "
            "Empty it yourself if a previous run left anything there.",
            extra=log_extra,
        )


def delete_temp_dir(temp_dir, delete_tmp, roots):
    """
    Delete a scratch directory once it is no longer needed.

    Nothing is deleted if `delete_tmp` is 0, the directory was never created, or it resolves
    outside every one of `roots` (warned about).

    Args:
        temp_dir (str): The scratch directory.
        delete_tmp (int): 1 to delete it, 0 to keep it.
        roots (list): The directories it may be deleted under.

    Returns:
        None
    """
    log_extra = {"stage": "Cleaning Up"}

    # Only a step that used scratch space created it
    if not os.path.isdir(temp_dir):
        logger.debug(f"No temp_dir was created at {temp_dir}; nothing to delete", extra=log_extra)
        return

    if delete_tmp == 0:
        logger.info(f"delete_tmp is 0; keeping {temp_dir}", extra=log_extra)
        return

    if not is_within(temp_dir, roots):
        logger.warning(
            f"Not deleting temp_dir {temp_dir}: it is outside {' and '.join(roots)}. "
            "Remove it yourself if that was intended.",
            extra=log_extra,
        )
        return
    shutil.rmtree(temp_dir)


@contextmanager
def temp_dir_scope(temp_dir, delete_tmp, roots):
    """
    Own a scratch directory for the length of a `with` block: empty it on entry, delete it on exit.

    The emptying and deletion are `empty_temp_dir`'s and `delete_temp_dir`'s, with their guards. A
    block that raises deletes nothing, leaving the scratch for inspection. A scope for a directory an
    enclosing scope already holds does nothing, on entry or exit: the outer scope owns it.

    Args:
        temp_dir (str): The scratch directory.
        delete_tmp (int): 1 to delete it on exit, 0 to keep it.
        roots (list): The directories it may be emptied or deleted under.

    Yields:
        None
    """
    key = os.path.abspath(temp_dir)
    if key in HELD_TEMP_DIRS:
        yield
        return
    HELD_TEMP_DIRS.add(key)
    # Released however the block ends, or a failed run would leave the next one in this process
    # neither emptying nor deleting it
    try:
        empty_temp_dir(temp_dir, roots)
        yield
        delete_temp_dir(temp_dir, delete_tmp, roots)
    finally:
        HELD_TEMP_DIRS.discard(key)


@contextmanager
def log_to_file(log_path, verbosity):
    """
    Log the package to a file, at a level set by `verbosity`, for the length of a `with` block.

    Creates the file's directory, sets the `pocketmapper` logger's level and adds a file handler
    appending to `log_path`. Both are undone on exit, however the block ends. If the logger already
    writes to `log_path`, as inside an enclosing call, only the level is set.

    Args:
        log_path (str): The log file.
        verbosity (int): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR.

    Yields:
        None

    Raises:
        PocketMapperError: If the log's directory cannot be created.
    """
    make_dir(os.path.dirname(log_path), {"stage": "Configuring Settings"})
    if verbosity == 4:
        log_level = logging.DEBUG
    elif verbosity == 3:
        log_level = logging.INFO
    elif verbosity == 2:
        log_level = logging.WARNING
    else:
        log_level = logging.ERROR

    package_logger = logging.getLogger(PACKAGE_LOGGER)
    previous_log_level = package_logger.level
    package_logger.setLevel(log_level)
    # A second handler on the same file would write every line twice
    handler = None
    log_path = os.path.abspath(log_path)
    if not any(
        isinstance(existing, logging.FileHandler) and existing.baseFilename == log_path
        for existing in package_logger.handlers
    ):
        handler = format_handler(logging.FileHandler(log_path))
        package_logger.addHandler(handler)
    try:
        yield
    finally:
        if handler is not None:
            package_logger.removeHandler(handler)
            handler.close()
        package_logger.setLevel(previous_log_level)


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


# Foldseek's PDB database names its entries "<pdbid>-assembly<N>_<chain>", with a "-<copy>" suffix on
# chains duplicated within an assembly (e.g. "5ian-assembly1_B-2").
FOLDSEEK_PDB_ENTRY = re.compile(r"^(?P<pdb>[0-9A-Za-z]{4})-assembly(?P<assembly>\d+)_(?P<chain>.+)$")


def parse_foldseek_pdb_entry_name(name):
    """
    Resolve a Foldseek PDB-database entry name into the PDB ID and chain it came from.

    Args:
        name (str): A Foldseek database entry name.

    Returns:
        tuple: (pdb_id, chain), with pdb_id upper-cased and the assembly number and any chain-copy
            suffix discarded -- "5ian-assembly1_B-2" and "5ian-assembly2_B" both give ("5IAN", "B").
            None for a name that is not PDB-style.
    """
    match = FOLDSEEK_PDB_ENTRY.match(name)
    if match is None:
        return None
    return match.group("pdb").upper(), match.group("chain").split("-")[0]


def split_chain_info(chain_info):
    """
    Split a chain_info field into its domain chain and its motif chain.

    Args:
        chain_info (str): The middle field of an input entry: a lone chain "A", or "A_B" where A carries
            the pocket and B is its binding partner. Chains of any length split correctly. Not None.

    Returns:
        tuple: (domain_chain, motif_chain); motif_chain is None when the entry named no partner.
    """
    chains = chain_info.split("_")
    domain_chain = chains[0]
    motif_chain = chains[1] if len(chains) > 1 else None
    return domain_chain, motif_chain


def one_letter_code(res_name):
    """
    Map a three-letter residue name to its one-letter code, as Foldseek does.

    Modified residues map to their parent amino acid (CME -> C). Names Foldseek does not list,
    including ligands, water and nucleotides, map to X.

    Args:
        res_name (str): Residue name as it appears in the structure, e.g. "ALA".

    Returns:
        str: A single uppercase letter.
    """
    return FOLDSEEK_AA_CODES.get(res_name, "X")


def seq_to_uniprot_map(domain):
    """
    Map each 0-indexed position of a domain onto its 1-indexed UniProt position.

    Args:
        domain (str): start1-stop1_start2-stop2_... in UniProt coords, 1-indexed and inclusive. A
            domain need not be contiguous, hence a list of regions rather than one offset.

    Returns:
        dict: 0-indexed domain position -> 1-indexed UniProt position. Increasing and injective as long
            as the regions are increasing and non-overlapping.
    """
    regions = domain.split("_")
    seq_pos_to_uniprot_pos = {}
    pos = 0
    for region in regions:
        start, stop = region.split("-")
        for i in range(int(start), int(stop) + 1):
            seq_pos_to_uniprot_pos[pos] = i
            pos += 1
    return seq_pos_to_uniprot_pos


def read_offset_table(path):
    """
    Read an offset table: the UniProt coordinates of every entry in a Foldseek database.

    Args:
        path (str): Path to the tab-separated table. Its header row is discarded.

    Returns:
        dict: database entry name -> region spec, as `seq_to_uniprot_map` parses.
    """
    offsets = {}
    with open(path) as f:
        next(f, None)  # header: unique_id, uniprot_domain
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            entry, domain = line.split("\t")
            offsets[entry] = domain
    return offsets
