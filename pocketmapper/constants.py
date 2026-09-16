"""
Shared constants and declared table schemas.

This module imports nothing, which is why the tables several modules must agree on live here
rather than beside any one of their users: `ALIGNMENT_COLUMNS` and `FOLDSEEK_FORMAT_OUTPUT` are a
positional contract between the two aligners and the comparison, and `FOLDSEEK_AA_CODES` is the
table behind `lib.one_letter_code`. Each constant carries its own rationale above it.
"""

# Three-letter residue names to one-letter codes, copied verbatim from Foldseek's threeToOneAA
# (src/strucclustutils/GemmiWrapper.cpp) so local sequences match the ones Foldseek builds. Names
# not listed map to X there, and lib.one_letter_code does the same. gemmi's tabulated codes are not
# a substitute: they differ on 14 of these entries. Refresh this whenever Foldseek changes its table.
FOLDSEEK_AA_CODES = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ABA": "A",
    "ASP": "D",
    "ASX": "B",
    "CYS": "C",
    "CSH": "S",
    "GLN": "Q",
    "GLU": "E",
    "GLX": "Z",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "MSE": "M",
    "ORN": "A",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRY": "W",
    "TRP": "W",
    "TYR": "Y",
    "UNK": "X",
    "VAL": "V",
    "SEC": "C",
    "PYL": "O",
    "SEP": "S",
    "TPO": "T",
    "PCA": "E",
    "CSO": "C",
    "PTR": "Y",
    "KCX": "K",
    "CSD": "C",
    "LLP": "K",
    "CME": "C",
    "MLY": "K",
    "DAL": "A",
    "TYS": "Y",
    "OCS": "C",
    "M3L": "K",
    "FME": "M",
    "ALY": "K",
    "HYP": "P",
    "CAS": "C",
    "CRO": "T",
    "CSX": "C",
    "DPR": "P",
    "DGL": "E",
    "DVA": "V",
    "CSS": "C",
    "DPN": "F",
    "DSN": "S",
    "DLE": "L",
    "HIC": "H",
    "NLE": "L",
    "MVA": "V",
    "MLZ": "K",
    "CR2": "G",
    "SAR": "G",
    "DAR": "R",
    "DLY": "K",
    "YCM": "C",
    "NRQ": "M",
    "CGU": "E",
    "0TD": "D",
    "MLE": "L",
    "DAS": "D",
    "DTR": "W",
    "CXM": "M",
    "TPQ": "Y",
    "DCY": "C",
    "DSG": "N",
    "DTY": "Y",
    "DHI": "H",
    "MEN": "N",
    "DTH": "T",
    "SAC": "S",
    "DGN": "Q",
    "AIB": "A",
    "SMC": "C",
    "IAS": "D",
    "CIR": "R",
    "BMT": "T",
    "DIL": "I",
    "FGA": "E",
    "PHI": "F",
    "CRQ": "Q",
    "SME": "M",
    "GHP": "G",
    "MHO": "M",
    "NEP": "H",
    "TRQ": "W",
    "TOX": "W",
    "ALC": "A",
    "SCH": "C",
    "MDO": "A",
    "MAA": "A",
    "GYS": "S",
    "MK8": "L",
    "CR8": "H",
    "KPI": "K",
    "SCY": "C",
    "DHA": "S",
    "OMY": "Y",
    "CAF": "C",
    "0AF": "W",
    "SNN": "N",
    "MHS": "H",
    "SNC": "C",
    "PHD": "D",
    "B3E": "E",
    "MEA": "F",
    "MED": "M",
    "OAS": "S",
    "GL3": "G",
    "FVA": "V",
    "PHL": "F",
    "CRF": "T",
    "BFD": "D",
    "MEQ": "Q",
    "DAB": "A",
    "AGM": "R",
    "4BF": "Y",
    "B3A": "A",
    "B3D": "D",
    "B3K": "K",
    "B3Y": "Y",
    "BAL": "A",
    "DBZ": "A",
    "GPL": "K",
    "HSK": "H",
    "HY3": "P",
    "HZP": "P",
    "KYN": "W",
    "MGN": "Q",
}

# The root log format, shared by the CRITICAL-only handler PocketMapper installs at construction
# and by the dictConfig configure_logging replaces it with. `stage` is not a stock LogRecord
# attribute: lib.StageFilter supplies it from the emitting function's name for any record that does
# not carry one, and both handlers must run that filter or an outside record fails to format.
LOG_FORMAT = "%(levelname)s: %(stage)s - %(msg)s"


# Appended to every error/warning about a missing foldseek binary, so the install line is
# written once. Foldseek is an optional external dependency and is never bundled.
FOLDSEEK_INSTALL_HINT = (
    "Install it with: conda install -c conda-forge -c bioconda foldseek "
    "(precompiled binaries: https://dev.mmseqs.com/foldseek/)."
)

# The alignment table's columns, in order. This is a positional contract shared by three modules:
# foldseek_alignment passes FOLDSEEK_FORMAT_OUTPUT to Foldseek's --format-output, the local
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

# The transform sources step 7 can actually use, and the methods the setting accepts. "auto" is
# resolved to one of the other two by resolve_align_struct_method before anything downstream reads it,
# so StructureAligner.align_structs validates against the resolved pair rather than the whole set.
RESOLVED_ALIGN_STRUCT_METHODS = ("pocket", "foldseek")
ALIGN_STRUCT_METHODS = ("auto",) + RESOLVED_ALIGN_STRUCT_METHODS

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
