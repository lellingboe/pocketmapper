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


# The --help text. Kept to 80 columns and to one sentence per option: anything longer -- the input
# grammar, the databases, the output columns, the Foldseek fallback -- lives in the README, which
# the footer points at. Every option here must match the Settings dataclass and the README's Options
# table; nothing generates one from the other.
HELP_MESSAGE = """
PocketMapper - compare the binding surfaces of protein-protein interactions.

Usage:
    pocketmapper search [OPTIONS]

Options:
  --query STR         Query entry, or a file with one entry per line.
                      STRUCT[:CHAIN[:RESIDUES]], e.g. 4Q5J:B_F. (required)
  --target STR        Target entry, a file with one entry per line, or a
                      Foldseek DB name: human_domains, pdb. (required)
  --settings PATH     JSON file of {"option": value}; CLI args override it.
                      (default: none)
  --cache_dir DIR     Where structures, pockets and PISA responses are cached.
                      (default: pocketmapper_cache)
  --results_dir DIR   Where results are written.
                      (default: pocketmapper_results_<YYMMDD_HHMMSS>)
  --verbosity INT     Log level: 4=DEBUG, 3=INFO, 2=WARNING, else ERROR.
                      (default: 3)
  --foldseek BOOL     Require the Foldseek aligner (True) or forbid it (False);
                      unset auto-detects the binary. (default: unset)
  --align_count INT   How many top-scoring targets to superpose onto each
                      query; 0 disables. (default: 10)
  --align_struct_method STR
                      Which transform superposes a target onto its query:
                      auto, pocket or foldseek. (default: auto)
  --query_pocket_method STR
                      Force the query pocket method rather than inferring it
                      from the entry: pisa, passthrough, vdw or whole_chain.
                      (default: unset)
  --target_pocket_method STR
                      As --query_pocket_method, for targets; also accepts
                      foldseek_db. (default: unset)
  --help              Show this message and exit.

Advanced options, settable only in the settings JSON. All are paths; the
defaults below write <cache> for cache_dir and <results> for results_dir:
  structure_dir                        <cache>/ref_structures
  pocket_dir                           <cache>/pockets
  foldseek_tmp_dir                     <cache>/foldseek_tmp
  foldseek_preprocessed_structure_dir  <cache>/foldseek_preprocessed_structures
  fsdb_dir                             <cache>/fsdb
  query_dir                            <results>/query_structures
  target_dir                           <results>/target_structures
  aligned_structure_dir                <results>/aligned_structures
  alignment_path                       <results>/alignment.tsv
  pocket_comparison_path               <results>/pocket_comparison.tsv
  job_settings_path                    <results>/job_settings.json
  log_path                             <results>/info.log

Examples:
  # One pair, using Foldseek when the binary is installed and the built-in
  # BLOSUM62 aligner when it is not.
  pocketmapper search --query 4Q5J:B_F --target 4Q5J:A_E --results_dir ./out

  # Search a pocket against the bundled Foldseek DB of human domains.
  pocketmapper search --query 4Q5J:B_F --target human_domains

  # Batch mode: one entry per line in each file.
  pocketmapper search --query queries.txt --target targets.txt \\
      --settings config.json

Input grammar, databases, output columns and the Foldseek fallback are
documented in the README:
    https://github.com/lellingboe/pocketmapper
"""
