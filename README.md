PocketMapper
    ============

    PocketMapper is a command-line tool to fetch protein structures, compute PISA-derived pockets,
    extract atomic coordinates from mmCIF files, perform local or Foldseek alignments, compare
    pockets across structures and write results to disk. It is intended for comparative analysis
    of binding pockets between query and target protein chains.

    Features
    - Download and cache mmCIF files
    - Preprocess/mmCIF splitting using gemmi
    - Retrieve PISA interface/pocket information and store pocket data
    - Extract CA coordinates from divided structures
    - Perform local alignments or Foldseek-based alignments
    - Compare pockets using alignment and substitution scoring (BLOSUM62)
    - Save tabular results and auxiliary JSON files to a results directory

    Requirements
    - Python 3.8+
    - pandas
    - gemmi
    - pisa (project-specific downloader wrapper used by this package)
    - foldseek (optional, required only when using Foldseek alignment)
    - Additional dependencies: lib (project helper module), LocalAligner class (local_aligner.py)
    - Command-line wrapper: fire

    Installation
    1. Clone the repository (or copy the project into your workspace).
    2. Create a virtualenv and install dependencies:
        python -m venv .venv
        source .venv/bin/activate
        pip install -r requirements.txt

    3. Ensure external tools are available:
        - foldseek (if using Foldseek): install and available on PATH

    Quick start / Usage
    - Basic local alignment run for a single pair:
      pocketmapper search --query 1ABC_A_B --target 2XYZ_C_D --results_dir ./out

    - Batch mode using files with one PDB_CHAIN_CHAIN per line:
      pocketmapper search --query queries.txt --target targets.txt --settings config.json

    - Use bundled Foldseek DB:
      pocketmapper search --query 1ABC_A_B --target ted --foldseek True --results_dir ./out_fs

    Options
    - --query: Query identifier or path. Accepts a single PDB_CHAIN_CHAIN (e.g. 1ABC_A_B) or a file with one per line.
    - --target: Target identifier or path. Accepts single PDB_CHAIN_CHAIN, file, or 'ted' for bundled Foldseek DB.
    - --settings: Path to JSON settings file. CLI args override settings file.
    - --cache_dir: Directory for caching downloaded or intermediate files.
    - --results_dir: Directory to write results and temporary divided structures.
    - --verbose / --debug: Increase log verbosity.
    - --foldseek: If true, run Foldseek alignments (requires foldseek binary and appropriate DB).
    - --pisa_pockets: Whether to retrieve PISA pockets (default: true).

    Configuration (settings JSON)
    The settings JSON may include keys such as:
    - cache_dir, structure_dir, pocket_dir, pisa_dir, divided_struct_dir
    - results_dir, query_dir, target_dir, alignment_path, pocket_comparison_path
    - foldseek (bool), pisa_pockets (bool), structure (bool)

    Outputs
    - alignment.tsv: Alignment report (Foldseek or local aligner)
    - pocket_comparison.tsv: Final pocket comparison table
    - pisa_pockets and intermediate JSON snapshots under pisa_dir
    - unknown_ids.json (if unknown Foldseek aliases are encountered)
    - Divided mmCIF files and temporary directories under results_dir

    Design / Workflow
    1. Parse CLI args and settings
    2. Determine types of query/target (single PDB_CHAIN_CHAIN, file, or foldseek DB)
    3. Fetch mmCIF structures to cache (structure_dir)
    4. Preprocess and divide structures (gemmi) into per-domain files (divided_struct_dir)
    5. Retrieve PISA interface data and compute pockets (pisa)
    6. Extract CA coordinates from divided mmCIFs
    7. Perform alignment (local or foldseek)
    8. Compare pockets using alignments and BLOSUM scoring
    9. Save results and clean up temporary directories

    Extending or Debugging
    - Increase verbosity with --verbose or --debug to get more details in logs.
    - The library 'lib' contains helper functions for fetching mmCIFs, preprocessing and comparing pockets.
    - Local alignment logic is in local_aligner.py (LocalAligner).
    - The default logging writes to test.log in the current working directory.

    License
    - Add your preferred license information here.

    Contributing
    - Report issues or open pull requests against the repository.
    - Add tests for new functionality and keep changes small and focused.

    Contact / Authors
    - See project repository for maintainer and contributor information.