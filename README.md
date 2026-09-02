# PocketMapper

PocketMapper is a command-line tool to compare the binding surfaces of protein-protein interactions.
PocketMapper fetches protein structures from the PDB, fetches contact residues from PDBe (PISA), aligns structures based on sequence (pairwise BLOSUM62) or structure (Foldseek) alignment, and returns information about overlapping surfaces including the conservation of surface residues and their RMSD.
It is intended for comparative analysis of binding pockets between query and target protein chains.

## Installation
PocketMapper supports Python 3.10 to 3.14 and is available on [PyPI](https://pypi.org/project/pocketmapper/). [Foldseek](https://github.com/steineggerlab/foldseek) is an optional external binary, but installing it is recommended: PocketMapper uses it by default when it is on `PATH`, and falls back to the built-in BLOSUM62 sequence aligner (with a warning) when it is not. Foldseek is required to search a bundled Foldseek database. Structural superposition works either way: with Foldseek it can use the whole-chain fit, and the local aligner superposes on the pocket instead (see `--align_struct_method`).
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
STRUCTURE[:CHAIN[:RESIDUES]]
```
- **STRUCTURE** — a 4-character PDB ID (`4Q5J`), a UniProt accession (`P04637`, fetched from AlphaFold), or a path to a local structure file.
- **CHAIN** — the chain the pocket belongs to. For the interface-based methods this is two chains joined by an underscore (`B_F`), where the first is the chain carrying the pocket and the second is its binding partner. Optional; omitted, it defaults to chain `A`.
- **RESIDUES** — optional comma-separated author residue numbers (`10,11,12`). Omitting the whole pocket part makes the entry an open search — see below.

How the pocket residues are derived is inferred from the shape of the entry:

| Example | Method | Meaning |
| --- | --- | --- |
| `4Q5J:B_F` | pisa | Interface residues of chain B against chain F, taken from PDBe PISA |
| `4Q5J:A:10,11,12` | passthrough | Exactly the residues listed, on chain A |
| `4Q5J:A_B:10,11,12` | vdw | Van der Waals contact residues on chain A against chain B |
| `P04637:A:10,11,12` | passthrough | AlphaFold model for the accession, residues listed |
| `4Q5J:B` | whole_chain | Open search — the whole of chain B is treated as the pocket |
| `4Q5J` | whole_chain | Open search over the default chain, `A` |
| `human_domains` | — | Bundled Foldseek database (valid as a target only) |

PISA is only available for PDB entries, and AlphaFold/local-file entries only support `passthrough`,
`vdw` and `whole_chain`. A passthrough entry needs an explicit residue list — without one the entry is
an open search instead. The inferred method can be overridden with `--query_pocket_method` /
`--target_pocket_method`.

### Open searches
An entry that names a structure but no pocket asks a different question: *does the query pocket resemble
anything on this chain at all?* The whole chain becomes the pocket, so the answer is carried by
`overlap_count` and `pocket_1_overlap_ids` — how many of the query pocket's residues the target chain
covers, and which ones.

Because there is no pocket on the target to describe, the descriptive and length-normalised columns
(`pocket_2_res_ids`, `pocket_2_len`, `pocket_2_seq`, `pocket_2_pct_aln`, `jaccard_index`) are left empty —
the same shape a Foldseek-database row has. `jaccard_index` in particular needs a second pocket to size
against, and a whole chain's length would swamp the union it normalises by. The overlap itself is still
fully reported: `pocket_2_overlap_ids` gives the author residue numbers the query pocket maps onto, and
because a named structure has real coordinates the superposition columns (`rmsd`, `ca_dists`, the
transforms) are populated too, which a Foldseek-database row cannot offer.

A Foldseek-database row differs in one further way: `pocket_2_overlap_ids` reports UniProt residue
numbers when the database ships an `offset_table.tsv` beside it, as the bundled `human_domains` DB
does. Its entries are domains carved out of UniProt sequences, so the position within an entry means
nothing outside PocketMapper. A database with no such table — any you supply yourself — reports
0-indexed positions within the entry instead, and says so in the log. The query side
(`pocket_1_overlap_ids`) is always author residue numbers either way.

Open and pocketed targets can be mixed freely in one batch file.

For batch runs, pass a path to a file containing one such entry per line instead of a single entry.

### Example commands
A single pair, using Foldseek if it is installed and the local aligner otherwise:
```
pocketmapper search --query 4Q5J:B_F --target 4Q5J:A_E --results_dir ./out_fs
```
Forcing the local BLOSUM62 aligner even when Foldseek is available:
```
pocketmapper search --query 4Q5J:B_F --target 4Q5J:A_E --foldseek False --results_dir ./out_local
```
Searching a chain against the bundled Foldseek DB of human domains (requires the foldseek binary):
```
pocketmapper search --query 4Q5J:B_F --target human_domains --results_dir ./out_hd
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
--foldseek: Whether to run Foldseek alignments instead of the local BLOSUM62 aligner. Left unset, Foldseek is used when its binary is on `PATH` and the local aligner is used with a warning when it is not. `True` requires Foldseek and errors if the binary is missing; `False` always uses the local aligner. A Foldseek DB target always needs the binary.\
--query_pocket_method / --target_pocket_method: Force the pocket method ('pisa', 'passthrough', 'vdw', 'whole_chain', or 'foldseek_db' for a target) instead of inferring it from the entry.\
--align_count: Number of top-scoring targets to superpose onto each query (default 10, 0 disables).\
--align_struct_method: Which transform superposes a target onto its query. `foldseek` uses Foldseek's whole-chain fit and needs the binary; `pocket` uses the fit of the two pockets on their overlapping residues, taken from `pocket_comparison.tsv`, so the pockets rather than the chains end up on top of each other; `auto` (the default) means `foldseek` when Foldseek is in use and `pocket` with the local aligner, which produces no chain transform. `pocket` is unavailable against a Foldseek DB target.\
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
- pisa_pockets.json / passthrough_pockets.json / vdw_pockets.json / whole_chain_pockets.json and cached PISA API responses under the cache pocket directory. The pocket files are written for inspection only and never read back; each pocket is an object of metadata fields plus a `residues` map keyed by author residue number
- unknown_ids.json (if unknown Foldseek aliases are encountered e.g., MSE -> M)
- Divided mmCIF files and temporary directories under results_dir

## Contact / Authors
See project repository for maintainer and contributor information.

## TBD
- Align function to generate aligned structures
- Information on results table
- Some sort of visual output
- Information of querying alphafold domains with human_domains