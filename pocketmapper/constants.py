"""
Shared constants and declared table schemas.

This module imports nothing, which is why the tables several modules must agree on live here
rather than beside any one of their users: `ALIGNMENT_COLUMNS` and `FOLDSEEK_FORMAT_OUTPUT` are a
positional contract between the two aligners and the comparison, and `SINGLE_AA_CODE` is read by
both the pocket parser and the local aligner. Each constant carries its own rationale above it.
"""

# TODO keep phospho information
SINGLE_AA_CODE = {
    "CYS": "C",
    "ASP": "D",
    "SER": "S",
    "GLN": "Q",
    "LYS": "K",
    "ILE": "I",
    "PRO": "P",
    "THR": "T",
    "PHE": "F",
    "ASN": "N",
    "GLY": "G",
    "HIS": "H",
    "LEU": "L",
    "ARG": "R",
    "TRP": "W",
    "ALA": "A",
    "VAL": "V",
    "GLU": "E",
    "TYR": "Y",
    "MET": "M",
    "SEP": "S",  # phosphoserine
    "TPO": "T",  # phosphothreonine
    "PTR": "Y",  # phosphotyrosine
    "MSE": "M",  # selenomethionine
}

# Appended to every error/warning about a missing foldseek binary, so the install line is
# written once. Foldseek is an optional external dependency and is never bundled.
FOLDSEEK_INSTALL_HINT = (
    "Install it with: conda install -c conda-forge -c bioconda foldseek "
    "(precompiled binaries: https://dev.mmseqs.com/foldseek/)."
)

# The alignment table's columns, in order. This is a positional contract shared by three modules:
# _foldseek_alignment passes FOLDSEEK_FORMAT_OUTPUT to Foldseek's --format-output, the local
# SequenceAligner builds the same columns in the same order, and pocket_comparison unpacks each row
# positionally into an AlignmentRow. Reordering this list moves all three together; editing any one
# of them in isolation breaks the comparison silently, which is why the list lives here.
ALIGNMENT_COLUMNS = [
    "query",
    "target",
    "fident",
    "alnlen",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "tstart",
    "tend",
    "evalue",
    "lddt",
    "qaln",
    "taln",
    "u",
    "t",
    "qseq",
    "tseq",
]

FOLDSEEK_FORMAT_OUTPUT = ",".join(ALIGNMENT_COLUMNS)

# The structural-alignment methods step 7 accepts. "auto" is resolved to one of the other two by
# _resolve_align_struct_method before anything downstream reads it.
ALIGN_STRUCT_METHODS = ("auto", "pocket", "foldseek")

# The chain used when an entry names a structure but no chain at all ("4Q5J"). AlphaFold models are
# always a single chain A, and it is the first chain of most PDB entries.
DEFAULT_CHAIN = "A"


# The one part of `search --help` that argparse cannot generate: the examples. Every option,
# including the twelve paths, is now built from the parser in cli.py, so nothing about them is
# duplicated here. It hangs off the `search` subparser only -- the bare `pocketmapper --help` lists
# subcommands and nothing else. Kept to 80 columns. Anything longer -- the input grammar, the
# databases, the output columns, the Foldseek fallback -- lives in the README, which the footer
# points at.
CLI_SEARCH_EPILOG = """
Examples:
  # One pair, using Foldseek when the binary is installed and the built-in
  # BLOSUM62 aligner when it is not.
  pocketmapper search 4Q5J:B_F 4Q5J:A_E --results_dir ./out

  # Search a pocket against the bundled Foldseek DB of human domains.
  pocketmapper search 4Q5J:B_F human_domains

  # Batch mode: one entry per line in each file.
  pocketmapper search queries.txt targets.txt --settings config.json

Input grammar, databases, output columns and the Foldseek fallback are
documented in the README:
    https://github.com/lellingboe/pocketmapper
"""
