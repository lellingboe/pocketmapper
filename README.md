# PocketMapper

PocketMapper is a command-line tool to compare the binding surfaces of protein-protein interactions.
PocketMapper fetches protein structures from the PDB, fetches contact residues from PDBe (PISA), aligns structures based on sequence (pairwise BLOSUM62) or structure (Foldseek) alignment, and returns information about overlapping surfaces including the conservation of surface residues and their RMSD.
It is intended for comparative analysis of binding pockets between query and target protein chains.

## Installation
PocketMapper has been tested with Python 3.12 and is available on [PyPI](https://pypi.org/project/pocketmapper/). Optionally, [Foldseek](https://github.com/steineggerlab/foldseek) can be installed to enable efficient structural alignment.
```
# Setup conda environment for pocketmapper
conda create --name=pocketmapper python=3.12
conda activate pocketmapper

# Pip installation of PocketMapper
pip install pocketmapper

# Optional - Conda installation of Foldseek
conda install -c conda-forge -c bioconda foldseek
```
Foldseek also has precompiled binaries available at https://dev.mmseqs.com/foldseek/

## Usage
### Input format
Query and target entries are colon-separated:
```
STRUCTURE:CHAIN[:RESIDUES]
```
- **STRUCTURE** — a 4-character PDB ID (`4Q5J`), a UniProt accession (`P04637`, fetched from AlphaFold), or a path to a local structure file.
- **CHAIN** — the chain the pocket belongs to. For the interface-based methods this is two chains joined by an underscore (`B_F`), where the first is the chain carrying the pocket and the second is its binding partner.
- **RESIDUES** — optional comma-separated author residue numbers (`10,11,12`).

How the pocket residues are derived is inferred from the shape of the entry:

| Example | Method | Meaning |
| --- | --- | --- |
| `4Q5J:B_F` | pisa | Interface residues of chain B against chain F, taken from PDBe PISA |
| `4Q5J:A:10,11,12` | passthrough | Exactly the residues listed, on chain A |
| `4Q5J:A_B:10,11,12` | vdw | Van der Waals contact residues on chain A against chain B |
| `P04637:A:10,11,12` | passthrough | AlphaFold model for the accession, residues listed |
| `human_domains` | — | Bundled Foldseek database (valid as a target only) |

PISA is only available for PDB entries, and AlphaFold/local-file entries only support `passthrough` and
`vdw`. A passthrough entry needs an explicit residue list. The inferred method can be overridden with
`--query_pocket_method` / `--target_pocket_method`.

For batch runs, pass a path to a file containing one such entry per line instead of a single entry.

### Example commands
Structural (Foldseek) alignment for a single pair:
```
pocketmapper search --query 4Q5J:B_F --target 4Q5J:A_E --foldseek True --results_dir ./out_fs
```
Searching a chain against the bundled Foldseek DB of human domains:
```
pocketmapper search --query 4Q5J:B_F --target human_domains --foldseek True --results_dir ./out_hd
```
Batch mode using files with one entry per line:
```
pocketmapper search --query queries.txt --target targets.txt --settings config.json
```

### Options
--query: Query entry, or path to a file with one entry per line (see Input format).\
--target: Target entry, path to a file with one entry per line, or the name of a bundled Foldseek DB ('human_domains', 'pdb').\
--settings: Path to JSON settings file. CLI args override settings file. Same available arguments as command line.\
--cache_dir: Directory for caching downloaded or intermediate files.\
--results_dir: Directory to write results and temporary divided structures.\
--verbosity: Log level, 4=DEBUG, 3=INFO (default), 2=WARNING, anything else=ERROR.\
--foldseek: If true, run Foldseek alignments (requires foldseek binary). Required when the target is a Foldseek DB.\
--query_pocket_method / --target_pocket_method: Force the pocket method ('pisa', 'passthrough', 'vdw') instead of inferring it from the entry.\
--align_count: Number of top-scoring targets to superpose onto each query (default 10, 0 disables).\
--help: Display the help message

### Features
- Download and cache mmCIF files
- Retrieve PISA interface definition
- Extract CA coordinates from PDB structures
- Perform local alignments or Foldseek-based alignments
- Compare pockets using alignment and substitution scoring (BLOSUM62)
- Save tabular results and auxiliary JSON files to a results directory

### Outputs
- alignment.tsv: Alignment report (Foldseek or local aligner)
- pocket_comparison.tsv: Final pocket comparison table
- pisa_pockets.json / passthrough_pockets.json / vdw_pockets.json and cached PISA API responses under the cache pocket directory
- unknown_ids.json (if unknown Foldseek aliases are encountered e.g., MSE -> M)
- Divided mmCIF files and temporary directories under results_dir

## Contact / Authors
See project repository for maintainer and contributor information.

## TBD
- Align function to generate aligned structures
- Information on results table
- Some sort of visual output
- Information of querying alphafold domains with human_domains