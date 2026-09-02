# PocketMapper

PocketMapper is a command-line tool to compare the binding surfaces of protein-protein interactions.
It fetches protein structures from the PDB, fetches contact residues from PDBe (PISA), aligns structures
by sequence (pairwise BLOSUM62) or structure (Foldseek), and reports the overlap between two binding
surfaces — which residues are shared, how well they are conserved, and how closely they superpose.
It is intended for comparative analysis of binding pockets between query and target protein chains.

## Installation

### Dependencies
PocketMapper supports **Python 3.10 to 3.14** and is published on
[PyPI](https://pypi.org/project/pocketmapper/). Its Python dependencies (biopython, fire, numpy, pandas,
tqdm, gemmi) are installed by pip.

[Foldseek](https://github.com/steineggerlab/foldseek) is an optional external binary, but installing it is
recommended:

- PocketMapper uses Foldseek by default whenever it is on `PATH`, and falls back to the built-in BLOSUM62
  sequence aligner — with a warning — when it is not.
- Foldseek is **required** to search a Foldseek database target (`human_domains`, `pdb`).
- Structural superposition works either way: with Foldseek it can use the whole-chain fit, and the local
  aligner superposes on the pocket instead (see [Advanced options](#advanced-options)).

### Install with pip
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

```
pocketmapper search --query 4Q5J:B_F --target 4Q5J:A_E --results_dir ./out
```

This compares the PISA interface of chain B against chain F of 4Q5J with the interface of chain A against
chain E of the same entry. Under the hood, one `search` run:

1. Resolves the query and target entries and works out how each pocket should be derived.
2. Downloads and caches the mmCIF files it needs (or an AlphaFold model, for a UniProt accession).
3. Retrieves the PISA interface definitions, or computes van der Waals contacts, or takes the residues
   you listed.
4. Aligns the query and target chains with Foldseek or the local BLOSUM62 aligner.
5. Projects both pockets onto that alignment, scores their overlap, and superposes them on their shared
   residues.
6. Writes `alignment.tsv`, `pocket_comparison.tsv` and the auxiliary JSON files into `--results_dir`, and
   superposed structures for the top hits into `aligned_structures/`.

Downloads are cached in `--cache_dir`, so a second run over the same structures is fast.

### Options

| Option | Type | Default | Summary |
| --- | --- | --- | --- |
| `--query` | str | *required* | Query entry, or a file with one entry per line (see [Input format](#input-format)). |
| `--target` | str | *required* | Target entry, a file with one entry per line, or a Foldseek DB name (`human_domains`, `pdb`). |
| `--settings` | path | none | JSON file of `{"option": value}`; explicit CLI arguments override it. |
| `--cache_dir` | str | `pocketmapper_cache` | Where structures, pockets and PISA responses are cached. |
| `--results_dir` | str | `pocketmapper_results_<YYMMDD_HHMMSS>` | Where results are written. |
| `--verbosity` | int | `3` | Log level: 4=DEBUG, 3=INFO, 2=WARNING, anything else=ERROR. |
| `--foldseek` | bool | unset (auto) | Require the Foldseek aligner (`True`) or forbid it (`False`); unset auto-detects the binary. |
| `--align_count` | int | `10` | How many top-scoring targets to superpose onto each query; `0` disables. |
| `--align_struct_method` | str | `auto` | Which transform superposes a target onto its query: `auto`, `pocket` or `foldseek`. |
| `--query_pocket_method` | str | unset | Force the query pocket method instead of inferring it: `pisa`, `passthrough`, `vdw`, `whole_chain`. |
| `--target_pocket_method` | str | unset | As `--query_pocket_method`, for targets; also accepts `foldseek_db`. |
| `--help` | flag | — | Show the help message and exit. |

The path settings not listed here are covered under [Advanced options](#advanced-options).

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
`vdw` and `whole_chain`. A local file with an interface-style chain spec (`4Q5J.cif.gz:B_F`) therefore
resolves to `vdw`, not `pisa`. A passthrough entry needs an explicit residue list — without one the entry
is an open search instead. The inferred method can be overridden with `--query_pocket_method` /
`--target_pocket_method`.

For batch runs, pass a path to a file containing one such entry per line instead of a single entry.
Query and target files are read independently, and every query is compared against every target.

### Closed vs Open Search
The two shapes of entry ask different questions.

A **closed search** names a pocket on both sides (`4Q5J:B_F` against `4Q5J:A_E`) and asks *how do these two
pockets compare?* Every column of `pocket_comparison.tsv` is populated: both pockets are described, the
overlap is scored, and `jaccard_index` sizes the shared residues against their union.

An **open search** names a structure but no pocket (`4Q5J:B`, or just `4Q5J`) and asks *does the query
pocket resemble anything on this chain at all?* The whole chain becomes the pocket, so the answer is
carried by `overlap_count` and `pocket_1_overlap_ids` — how many of the query pocket's residues the target
chain covers, and which ones.

Because there is no pocket on the target to describe, the descriptive and length-normalised columns
(`pocket_2_res_ids`, `pocket_2_len`, `pocket_2_seq`, `pocket_2_pct_aln`, `jaccard_index`) are left empty —
the same shape a Foldseek-database row has. `jaccard_index` in particular needs a second pocket to size
against, and a whole chain's length would swamp the union it normalises by. The overlap itself is still
fully reported: `pocket_2_overlap_ids` gives the author residue numbers the query pocket maps onto, and
because a named structure has real coordinates the superposition columns (`rmsd`, `ca_dists`, the
transforms) are populated too, which a Foldseek-database row cannot offer.

Open and closed targets can be mixed freely in one batch file.

### Databases
A target can name a Foldseek database instead of a structure, which searches the query pocket against
every entry in it. Both bundled databases require the foldseek binary, and
`--align_struct_method pocket` is rejected against either (there is no second pocket to fit to).

**`human_domains`** ships inside the package — no download. Its entries are structural domains carved out
of human UniProt sequences. There is no interface to compute on a domain, so the target side behaves like
an open search: no `pocket_2_*` descriptors, no `jaccard_index`, no superposition columns. The database
ships an offset table, so `pocket_2_overlap_ids` is reported in **UniProt residue numbers** rather than
positions within the domain. The query side (`pocket_1_overlap_ids`) is always author residue numbers.

**`pdb`** is downloaded on first use (`foldseek databases PDB`) into `<cache_dir>/fsdb/pdb`. Hits are
real PDB chains, so PocketMapper fetches PISA data for each one and compares against a real interface
pocket; hits with no usable PISA data are dropped rather than compared against a stand-in. **This is
slow the first time**: `4Q5J:B_F` returns roughly 4,970 hits across 3,620 entries, and PISA is fetched
per entry behind a rate-limiting sleep, so the first run takes hours. Both counts are logged before the
fetching starts, and reruns are fast from the interface cache.

A Foldseek database you build yourself also works as a target (pass its path, with
`--target_pocket_method foldseek_db`), but without an offset table beside it the target residue ids are
0-indexed positions within the entry. That is logged at INFO, because the same column then means
different things on different runs.

### Outputs
Everything is written under `--results_dir`:

| File | Contents |
| --- | --- |
| `pocket_comparison.tsv` | The main result: one row per pocket pair. Columns below. |
| `alignment.tsv` | The chain alignments the comparison is built on, from Foldseek or the local aligner. |
| `aligned_structures/*.pdb` | The top `--align_count` targets superposed onto each query. Named after the query, so identify a file by its `MOLECULE` records rather than its filename. |
| `job_settings.json` | The fully resolved settings for the run. |
| `info.log` | The run log. |
| `unknown_ids.json` | Residue codes the aligner and the structure disagreed on (e.g. MSE -> M). Only written if any were seen. |
| `incorrect_mapping.json` | Pockets dropped because their own sequence disagreed with the aligner's for that chain — typically assembly vs asymmetric unit numbering. Only written if any were dropped. |

Under `--cache_dir` and surviving between runs: the downloaded mmCIF files (`ref_structures/`), the
raw PISA API responses (`pockets/pisa_responses/`), and the derived pockets themselves
(`pockets/pisa_pockets.json`, `passthrough_pockets.json`, `vdw_pockets.json`,
`whole_chain_pockets.json`). The pocket files are written for inspection only and never read back;
each pocket is an object of metadata fields plus a `residues` map keyed by author residue number.

#### `pocket_comparison.tsv` columns
Pocket 1 is always the query, pocket 2 the target. The table always has all of these columns: a
comparison that stops early — no overlapping residues, no coordinates, fewer than three residues to
superpose, or an open target — leaves the remaining fields **empty** rather than dropping them.

*Identity*

| Column | Meaning |
| --- | --- |
| `pocket_1`, `pocket_2` | The two pocket ids, in the input-entry form (`4Q5J:B_F`). |
| `evalue` | Alignment E-value for the underlying chain pair. `-` with the local aligner. |
| `lddt` | Alignment lDDT reported by Foldseek. `-` with the local aligner. |

*Pocket description*

| Column | Meaning |
| --- | --- |
| `pocket_1_res_ids`, `pocket_2_res_ids` | The pocket's author residue numbers, comma-separated. |
| `pocket_1_len`, `pocket_2_len` | Number of residues in the pocket. |
| `pocket_1_seq`, `pocket_2_seq` | The pocket's residues as single-letter codes, in `res_ids` order. |
| `pocket_1_pct_aln`, `pocket_2_pct_aln` | Fraction of the pocket's residues that fall inside the aligned region at all. A low value means the pocket sits largely outside the alignment. |

*Overlap*

| Column | Meaning |
| --- | --- |
| `overlap_count` | How many alignment positions both pockets occupy. The headline number. |
| `pocket_1_overlap_ids` | The query residues in the overlap, as author residue numbers. |
| `pocket_2_overlap_ids` | The target residues they map onto, in the same order. UniProt numbers for a `human_domains` target; see [Databases](#databases). |
| `jaccard_index` | `overlap_count` divided by the union of the two pockets. Empty for an open or database target. |

*Overlap scoring (BLOSUM62)*

| Column | Meaning |
| --- | --- |
| `pocket_1_seq_overlap`, `pocket_2_seq_overlap` | The two overlap sequences, aligned position for position. |
| `overlap_identity` | Fraction of overlap positions where the two residues are identical. |
| `overlap_similarity_binary` | Fraction of overlap positions with a positive BLOSUM62 score — how much of the overlap is conservatively substituted, ignoring how strongly. |
| `overlap_similarity_1_2` | Mean BLOSUM62 score over the overlap, normalised per position by the query residue's self-score, so an identical pair scores 1.0. |
| `overlap_similarity_2_1` | The same, normalised against the target residue instead. The two differ because self-scores differ between residues. |
| `min_overlap_similarity`, `max_overlap_similarity` | The smaller and larger of the two directional scores. |

*Superposition* — populated only when both pockets have coordinates and at least three residues overlap.

| Column | Meaning |
| --- | --- |
| `p2_to_p1_u`, `p2_to_p1_t` | Rotation (nine values) and translation (three) that put pocket 2 onto pocket 1, fitted on the overlapping CA atoms. |
| `p1_to_p2_u`, `p1_to_p2_t` | The same fit in the other direction. |
| `rmsd` | RMSD of the overlapping CA atoms after superposition, in Å. |
| `ca_dists` | Per-residue CA distance after superposition, comma-separated, in overlap order. |

The rotation matrices are stored in Biopython's right-multiplying convention
(`coords · u + t`). If you feed one to gemmi, which left-multiplies, transpose it first — or use
`pocketmapper.pocket_comparison.parse_pocket_transform`, which does the conversion for you.

### Examples
A single pair, using Foldseek if it is installed and the local aligner otherwise:
```
pocketmapper search --query 4Q5J:B_F --target 4Q5J:A_E --results_dir ./out_fs
```
Forcing the local BLOSUM62 aligner even when Foldseek is available:
```
pocketmapper search --query 4Q5J:B_F --target 4Q5J:A_E --foldseek False --results_dir ./out_local
```
An open search — is this pocket anywhere on chain A of 4Q5J at all?
```
pocketmapper search --query 4Q5J:B_F --target 4Q5J:A --results_dir ./out_open
```
Searching a chain against the bundled Foldseek DB of human domains (requires the foldseek binary):
```
pocketmapper search --query 4Q5J:B_F --target human_domains --results_dir ./out_hd
```
Batch mode using files with one entry per line:
```
pocketmapper search --query queries.txt --target targets.txt --settings config.json
```

### Advanced options

**The settings JSON.** `--settings config.json` takes a flat `{"option": value}` object using the same
names as the CLI options. Settings are layered lowest to highest: built-in defaults, then the JSON file,
then any explicit CLI argument. It is the only way to set the path options below, all of which default to
a location under `--cache_dir` or `--results_dir`:

| Setting | Default |
| --- | --- |
| `structure_dir` | `<cache_dir>/ref_structures` |
| `pocket_dir` | `<cache_dir>/pockets` |
| `foldseek_tmp_dir` | `<cache_dir>/foldseek_tmp` |
| `foldseek_preprocessed_structure_dir` | `<cache_dir>/foldseek_preprocessed_structures` |
| `fsdb_dir` | `<cache_dir>/fsdb` |
| `query_dir` | `<results_dir>/query_structures` |
| `target_dir` | `<results_dir>/target_structures` |
| `aligned_structure_dir` | `<results_dir>/aligned_structures` |
| `alignment_path` | `<results_dir>/alignment.tsv` |
| `pocket_comparison_path` | `<results_dir>/pocket_comparison.tsv` |
| `job_settings_path` | `<results_dir>/job_settings.json` |
| `log_path` | `<results_dir>/info.log` |

`query_dir` and `target_dir` are deleted at the end of a run, as is `foldseek_tmp_dir` when Foldseek was used.

**`--foldseek` is three-valued.** Left unset it means *auto*: Foldseek is used when its binary is on
`PATH`, and the local BLOSUM62 aligner is used with a warning when it is not. `True` makes Foldseek a hard
requirement, so a missing binary is an error rather than a silent change of method — worth setting in a
pipeline where the two aligners are not interchangeable. `False` always uses the local aligner. A Foldseek
database target needs the binary whatever this is set to.

**`--align_struct_method` picks the transform used to write `aligned_structures/`.**

- `foldseek` — Foldseek's whole-chain fit. Puts the two *chains* on top of each other. Needs the binary.
- `pocket` — the fit of the two pockets on their overlapping residues, taken from
  `pocket_comparison.tsv`. Puts the two *pockets* on top of each other, which is usually what you want
  when the chains are otherwise unrelated. Unavailable against a Foldseek database target, where the
  target pocket has no coordinates to fit.
- `auto` (the default) — `foldseek` when Foldseek is in use, `pocket` with the local aligner, which
  produces no whole-chain transform at all.

**Using PocketMapper as a library.** `PocketMapper().search(...)` runs the same workflow as the CLI and
writes the same files; results come back through `results_dir`, not as a return value. The individual
components (`qt_processor`, `structure_fetcher`, `pisa_downloader`, `sequence_aligner`,
`structure_aligner`, `pocket_calculator`, ...) are each usable on their own. Note that `search()`
reconfigures the root logger and deletes its temporary directories on the way out.

## Contact / Authors
PocketMapper is developed by Lachlan Ellingboe (Lachlan.Ellingboe@icr.ac.uk).
Source, issues and feature requests: https://github.com/lellingboe/pocketmapper
