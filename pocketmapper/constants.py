"""
Shared constants and declared table schemas.

Imports nothing, so every module can import it. Tables that several modules must agree on live here
rather than beside any one of their users.
"""

# Three-letter residue names to one-letter codes, copied verbatim from Foldseek's threeToOneAA
# (src/strucclustutils/GemmiWrapper.cpp); names not listed map to X there. gemmi's tabulated codes
# are not a substitute: they differ on 14 of these entries. Refresh whenever Foldseek changes its table.
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

# The logger every module's logger sits under
PACKAGE_LOGGER = "pocketmapper"

# The format of every handler the package builds. `stage` is not a stock LogRecord attribute, so a
# handler using this format needs lib.StageFilter.
LOG_FORMAT = "%(levelname)s: %(stage)s - %(msg)s"


# Install instructions for the foldseek binary, an optional dependency that is never bundled
FOLDSEEK_INSTALL_HINT = (
    "Install it with: conda install -c conda-forge -c bioconda foldseek "
    "(precompiled binaries: https://dev.mmseqs.com/foldseek/)."
)

# The alignment table's columns, in order. The order is a positional contract, as binding as the
# names: add, remove or reorder a column here, never in one producer or reader alone.
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

# The transforms a structure can be superposed by, and the values the setting accepts. "auto" is
# resolved to one of the other two at run time.
RESOLVED_ALIGN_STRUCT_METHODS = ("pocket", "foldseek")
ALIGN_STRUCT_METHODS = ("auto",) + RESOLVED_ALIGN_STRUCT_METHODS

# The chain aligners the `aligner` setting accepts: Foldseek, or the local BLOSUM62 sequence aligner.
ALIGNERS = ("foldseek", "seq")

# The values the `delete_tmp` setting accepts: 1 deletes temp_dir at the end of the run, 0 keeps it.
DELETE_TMP_VALUES = (0, 1)

# Where PISA interfaces come from: static files on the EBI FTP server, or one paced PDBe API call per
# assembly. Both serve identical JSON.
PISA_SOURCES = ("ftp", "api")

# The structure types a structure-only parse can be forced to, and "auto", which infers one per entry.
# A local file is inferred only: forcing it would add nothing to the path check that infers it.
STRUCT_TYPES = ("auto", "pdb", "alphafold", "foldseek_db")

# Defaults for the settings that have a static value. Settings whose default depends on the run
# (results_dir, threads and the derived paths) default to None and are resolved at run time.
DEFAULT_CACHE_DIR = "pocketmapper_cache"
DEFAULT_VERBOSITY = 3
DEFAULT_ALIGN_COUNT = 10
DEFAULT_ALIGN_STRUCT_METHOD = "auto"
# "auto" infers the pocket method from each entry.
DEFAULT_POCKET_METHOD = "auto"
DEFAULT_ALIGNER = "foldseek"
DEFAULT_DELETE_TMP = 1
DEFAULT_PISA_SOURCE = "ftp"

# The chain used when an entry names a structure but no chain at all ("4Q5J"). AlphaFold models are
# always a single chain A, and it is the first chain of most PDB entries.
DEFAULT_CHAIN = "A"


# The examples for `search --help`, the one part argparse cannot generate. Kept to 80 columns;
# anything longer belongs in the README, which the footer points at.
CLI_SEARCH_EPILOG = """
Examples:
  # One pair, aligned with Foldseek (the binary must be installed).
  pocketmapper search 4Q5J:B_F 4Q5J:A_E --results_dir ./out

  # The same pair with the built-in BLOSUM62 sequence aligner instead.
  pocketmapper search 4Q5J:B_F 4Q5J:A_E --aligner seq --results_dir ./out

  # Search a pocket against the bundled Foldseek DB of human domains.
  pocketmapper search 4Q5J:B_F human_domains

  # Batch mode: one entry per line in each file.
  pocketmapper search queries.txt targets.txt --job_file job.json

  # Everything, query and target included, from a job file.
  pocketmapper search --job_file job.json

Input grammar, databases, output columns and the choice of aligner are
documented in the README:
    https://github.com/lellingboe/pocketmapper
"""

# The examples for every other command's --help, kept to 80 columns like CLI_SEARCH_EPILOG
CLI_STEP_EPILOG = """
Examples:
  # The search workflow one step at a time. parse writes the settings the
  # rest of the chain takes as its job file.
  pocketmapper parse 4Q5J:B_F 4Q5J:A_E --results_dir ./out
  pocketmapper fetch_structures --job_file ./out/parse_settings.json
  pocketmapper align --job_file ./out/parse_settings.json
  pocketmapper pockets --job_file ./out/parse_settings.json
  pocketmapper compare --job_file ./out/parse_settings.json
  pocketmapper superpose --job_file ./out/parse_settings.json

  # After a search, superpose again with other settings.
  pocketmapper superpose --job_file ./out/job_settings.json --align_count 3 \\
      --aligned_structure_dir ./out/top3

Each step parses the query and target entries its job file names, and
reads the earlier steps' files from --results_dir. The README documents each
step's inputs and outputs:
    https://github.com/lellingboe/pocketmapper
"""

# Command -> the examples for its --help
CLI_COMMAND_EPILOGS = {
    "parse": CLI_STEP_EPILOG,
    "fetch_structures": CLI_STEP_EPILOG,
    "align": CLI_STEP_EPILOG,
    "pockets": CLI_STEP_EPILOG,
    "compare": CLI_STEP_EPILOG,
    "superpose": CLI_STEP_EPILOG,
    "search": CLI_SEARCH_EPILOG,
}
