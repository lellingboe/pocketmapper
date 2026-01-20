# PocketMapper

PocketMapper is a command-line tool to compare the binding surfaces of protein-protein interactions.
PocketMapper fetches protein structures from the PDB, fetches contact residues from PDBe (PISA), aligns structures based on sequence (pairwise BLOSUM62) or structure (Foldseek) alignment, and returns information about overlapping surfaces including the conservation of surface residues and their RMSD.
It is intended for comparative analysis of binding pockets between query and target protein chains.

## Installation
PocketMapper has been tested with Python 3.12 and is available on [PyPI](https://pypi.org/project/pocketmapper/). Optionally, [Foldseek](https://github.com/steineggerlab/foldseek) can be installed to enable effecient structural alignment.
```
# Setup conda environment for pocketmapper
conda create --name=pocketmapper python=3.12

# Pip installation of PocketMapepr
pip install pocketmapper

# Optional - Conda installation of Foldseek
conda install -c conda-forge -c bioconda foldseek
```
Foldseek also has precompiled binaries available at https://dev.mmseqs.com/foldseek/

## Usage
### Example commands
Basic sequence alignment run for a single pair:
```
pocketmapper search --query 4Q5J_B_F --target 4Q5J_A_E --results_dir ./out
```
Batch mode using files with one PDB_CHAIN_CHAIN per line:
```
pocketmapper search --query queries.txt --target targets.txt --settings config.json
```
Using the bundled Foldseek DB:
```
pocketmapper search --query 4Q5J_B_F --target human_domains --foldseek True --results_dir ./out_fs
```

### Options
--query: Query identifier or path. Accepts a single PDB_CHAIN_CHAIN (e.g. 4Q5J_B_F) or a file with one per line.\
--target: Target identifier or path. Accepts single PDB_CHAIN_CHAIN, file, or 'human_domains' for bundled Foldseek DB.\
--settings: Path to JSON settings file. CLI args override settings file. Same available arguments as command line.\
--cache_dir: Directory for caching downloaded or intermediate files.\
--results_dir: Directory to write results and temporary divided structures.\
--verbose / --debug: Increase log verbosity.\
--foldseek: If true, run Foldseek alignments (requires foldseek binary).\
--help: Display the help message

### Features
- Download and cache mmCIF files
- Preprocess/mmCIF splitting using gemmi
- Retrieve PISA interface/pocket information and store pocket data
- Extract CA coordinates from divided structures
- Perform local alignments or Foldseek-based alignments
- Compare pockets using alignment and substitution scoring (BLOSUM62)
- Save tabular results and auxiliary JSON files to a results directory

### Outputs
- alignment.tsv: Alignment report (Foldseek or local aligner)
- pocket_comparison.tsv: Final pocket comparison table
- pisa_pockets and intermediate JSON snapshots under pisa_dir
- unknown_ids.json (if unknown Foldseek aliases are encountered e.g., MSE -> M)
- Divided mmCIF files and temporary directories under results_dir

## Contact / Authors
See project repository for maintainer and contributor information.