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

HELP_MESSAGE = """
    PocketMapper - A tool for mapping and analyzing protein pockets.

    Usage:
        pocketmapper search [OPTIONS]

    Primary options:
        --query QUERY            Query identifier or path. Accepts:
                    - 'STRUCT:CHAINS[:RESIDUES]' (e.g., 1ABC:A_B, 1ABC:A:10,11,12, P12345:A:10,11)
                    - path to a file listing such entries, one per line
        --target TARGET          Target identifier or path. Accepts:
                    - 'STRUCT:CHAINS[:RESIDUES]' (e.g., 2XYZ:C_D)
                    - path to a file listing such entries, one per line
                    - a bundled Foldseek DB name ('human_domains', 'pdb'), which requires the foldseek binary
        --settings FILE          Path to a JSON file of {"ARG": "VALUE", ...} (overridden by explicit CLI args)
        --cache_dir DIR          Directory for caching files
        --results_dir DIR        Directory for writing results
        --verbosity LEVEL        Set verbosity level (4=DEBUG, 3=INFO, 2=WARNING, else ERROR)
        --foldseek BOOL          Whether to use foldseek for structure alignment instead of local sequence
                                 alignment. Foldseek is an external binary that must be on PATH; if it is not
                                 installed it can be added with:
                                     conda install -c conda-forge -c bioconda foldseek
                                 Left unset, foldseek is used when the binary is available and the local
                                 BLOSUM62 aligner is used with a warning when it is not.
                                   True  -- require foldseek; a missing binary is an error.
                                   False -- always use the local aligner. Structural superposition is
                                            unavailable there, so aligned_structures/*.pdb hold only
                                            the query.
        --align_count N          Number of top targets to superpose onto each query (0 disables, default 10)
        --query_pocket_method M  Force the query pocket method ('pisa', 'passthrough' or 'vdw')
        --target_pocket_method M Force the target pocket method ('pisa', 'passthrough' or 'vdw')
        --help                   Show this help message and exit

    Advanced Options (set via settings JSON):
        structure_dir                          Directory to store downloaded/available structures
        pocket_dir                             Directory to store calculated pockets
        foldseek_preprocessed_structure_dir    Directory for preprocessed/divided structures
        foldseek_tmp_dir                       Scratch directory for Foldseek
        fsdb_dir                               Directory for downloaded Foldseek databases
        query_dir                              Temporary directory for query divided structures
        target_dir                             Temporary directory for target divided structures
        aligned_structure_dir                  Directory to write superposed structures to
        alignment_path                         Path to write alignment TSV
        pocket_comparison_path                 Path to write pocket comparison TSV
        job_settings_path                      Path to write the resolved settings to
        log_path                               Path to write the run log to

    Description:
        Orchestrates fetching/preprocessing of structures, runs local or Foldseek alignments,
        derives pockets (PISA, explicit residue lists, or VdW contacts), extracts atom coordinates
        from mmCIF files, compares pockets using alignments and scoring, and writes results to the
        results directory.

    Examples:
        # Single pair with default settings (uses Foldseek when the binary is installed)
        pocketmapper search --query 1ABC:A_B --target 2XYZ:C_D --results_dir ./out

        # Batch mode using files with one entry per line
        pocketmapper search --query queries.txt --target targets.txt --settings config.json

        # Search the built-in human_domains Foldseek DB (requires the foldseek binary)
        pocketmapper search --query 1ABC:A_B --target human_domains --results_dir ./out_fs

        # Force the local BLOSUM62 aligner even when foldseek is installed
        pocketmapper search --query 1ABC:A_B --target 2XYZ:C_D --foldseek False --results_dir ./out_local

        # Override cache and set verbosity to debug
        pocketmapper search --query 1ABC:A_B --target 2XYZ:C_D --cache_dir /tmp/cache --verbosity 4

    Notes:
        - Query/target inputs are interpreted either as single 'STRUCT:CHAINS[:RESIDUES]' strings or as file paths.
        - Boolean settings can be provided on the command line (e.g., --foldseek False).
        - Use a settings JSON to persist complex configurations; CLI options override settings file values.

    For more information, see the project README or the github repository.
"""
