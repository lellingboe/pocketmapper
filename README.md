# PocketMapper

PocketMapper is a command-line tool to compare the binding surfaces of protein-protein interactions.
It fetches protein structures from the PDB, fetches contact residues from PDBe (PISA), aligns structures
by sequence (pairwise BLOSUM62) or structure (Foldseek), and reports the overlap between two binding
surfaces — which residues are shared, how well they are conserved, and how closely they superpose.
It is intended for comparative analysis of binding pockets between query and target protein chains.

## Installation

### Dependencies

PocketMapper supports **Python 3.10 to 3.14** and is published on
[PyPI](https://pypi.org/project/pocketmapper/). Its Python dependencies (biopython, numpy, pandas, tqdm,
gemmi) are installed by pip.

[Foldseek](https://github.com/steineggerlab/foldseek) is an external binary. It is not installed by pip,
but it is the default aligner:

- PocketMapper uses Foldseek unless you pass `--aligner seq`, which selects the built-in BLOSUM62
  sequence aligner. Without `--aligner seq`, a missing or unrunnable binary is an error.
- Foldseek is **required** to search a Foldseek database target (`human_domains`, `pdb`).
- Structural superposition works either way: with Foldseek it can use the whole-chain fit, and the local
  aligner superposes on the pocket instead (see [Advanced options](#advanced-options)).

### Install with pip

```bash
# Pip installation of PocketMapper
pip install pocketmapper

# Conda installation of Foldseek (not needed with --aligner seq)
conda install -c conda-forge -c bioconda foldseek

# Or Homebrew installation of foldseek
brew tap brewsci/bio
brew trust brewsci/bio
brew install foldseek
```

Foldseek also has precompiled binaries available at <https://dev.mmseqs.com/foldseek/>.

## Usage

```bash
pocketmapper search 4Q5J:B_F 4Q5J:A_E --results_dir ./out
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

Downloads are cached in `--cache_dir`, so a second run over the same structures is fast. Each of these
steps is also a command of its own; see [Running one step at a time](#running-one-step-at-a-time).

### Arguments

Both are required, and are given positionally in this order — there is no `--query`/`--target`
spelling. Either may instead be set in a [job file](#advanced-options), in which case it must be
left off the command line.

| Argument | Summary |
| --- | --- |
| `QUERY` | Query entry, or a file with one entry per line (see [Input format](#input-format)). |
| `TARGET` | Target entry, a file with one entry per line, or a Foldseek DB name (`human_domains`, `pdb`). |

### Options

| Option | Type | Default | Summary |
| --- | --- | --- | --- |
| `--job_file` | path | none | JSON file of `{"option": value}`, query and target included. Arguments override it. |
| `--verbosity` | int | `3` | Log level: 4=DEBUG, 3=INFO, 2=WARNING, anything else=ERROR. |
| `--aligner` | str | `foldseek` | Chain aligner: `foldseek` (needs the binary) or `seq` (built-in BLOSUM62 sequence aligner). |
| `--query_pocket_method` | str | `auto` | Query pocket method: `auto` infers it from each entry; `pisa`, `passthrough`, `vdw`, `whole_chain` force it. Each entry is still checked against the method — see [Input format](#input-format). |
| `--target_pocket_method` | str | `auto` | As `--query_pocket_method`, for targets; also accepts `foldseek_db`. |
| `--threads` | int | one per core | Cap on the cores Foldseek uses. |
| `--pisa_source` | str | `ftp` | Where PISA interfaces are fetched from: `ftp` (EBI FTP server, concurrent) or `api` (PDBe API, one paced request per assembly). Both fill the same cache. |
| `--help` | flag | — | Show the help message and exit. |

#### In options

| Option | Default | Summary |
| --- | --- | --- |
| `--work_dir` | the current directory | Directory that entries, entries files and every relative path resolve against. Recorded absolute in the settings dump, so a step run later from another directory resolves them the same way. |

#### Aligned structure options

| Option | Type | Default | Summary |
| --- | --- | --- | --- |
| `--align_count` | int | `10` | How many top-scoring targets to superpose onto each query; `0` disables. |
| `--align_struct_method` | str | `auto` | Which transform superposes a target onto its query: `auto`, `pocket` or `foldseek`. |

#### Cache options

What survives a run and makes the next one fast. Each defaults to a location under `--cache_dir`.

| Option | Default | Summary |
| --- | --- | --- |
| `--cache_dir` | `pocketmapper_cache` | Where structures, pockets and PISA responses are cached. |
| `--pdb_dir` | `<cache_dir>/pdb_structures` | Cache of fetched PDB structures. |
| `--alphafold_dir` | `<cache_dir>/alphafold_structures` | Cache of fetched AlphaFold structures. |
| `--pocket_dir` | `<cache_dir>/pockets` | Cache of parsed pockets. |
| `--foldseek_preprocessed_structure_dir` | `<cache_dir>/foldseek_preprocessed_structures` | Cache of the single-chain structures Foldseek is given. |
| `--fsdb_dir` | `<cache_dir>/fsdb` | Cache of bundled Foldseek databases. |

#### Out options

What the run produces. Each defaults to a location under `--results_dir`, so you only need these to
split a run's outputs across directories.

| Option | Default | Summary |
| --- | --- | --- |
| `--results_dir` | `pocketmapper_results_<YYMMDD_HHMMSS>` | Where results are written. |
| `--aligned_structure_dir` | `<results_dir>/aligned_structures` | Where superposed structures for the top hits are written. |
| `--alignment_path` | `<results_dir>/alignment.tsv` | Where the alignment table is written. |
| `--pocket_comparison_path` | `<results_dir>/pocket_comparison.tsv` | Where the pocket comparison table is written. |
| `--pockets_path` | `<results_dir>/pockets.json` | Where the pockets are written. |
| `--failed_entries_path` | `<results_dir>/failed_entries.json` | Where the entries dropped along the way are listed, with the reason. |
| `--job_settings_path` | `<results_dir>/job_settings.json` | Where this run's resolved settings are dumped. |
| `--log_path` | `<results_dir>/info.log` | Where the run log is written. It is appended to. |

#### Temp options

Scratch space for a single run: emptied before the run uses it and deleted when it ends — set
`--delete_tmp 0` to keep it for inspection. The per-run query and target structures and Foldseek's
own scratch all live inside it, so it is only created with `--aligner foldseek`, or to download a
bundled Foldseek database.

| Option | Default | Summary |
| --- | --- | --- |
| `--temp_dir` | `<results_dir>/tmp` | Per-run scratch, emptied before use and deleted at the end of the run. |
| `--delete_tmp` | `1` | `1` deletes `--temp_dir` at the end of the run; `0` keeps it. |

### Running one step at a time

`search` runs six steps, and each is also a command of its own. Chained, they produce what `search`
produces. So you can stop after any step, look at its output and carry on, or rerun one step with
other settings.

```bash
pocketmapper parse 4Q5J:B_F 4Q5J:A_E --results_dir ./out
pocketmapper fetch_structures --job_file ./out/parse_settings.json
pocketmapper align --job_file ./out/parse_settings.json
pocketmapper pockets --job_file ./out/parse_settings.json
pocketmapper compare --job_file ./out/parse_settings.json
pocketmapper superpose --job_file ./out/parse_settings.json
```

| Command | Does | Reads | Writes |
| --- | --- | --- | --- |
| `parse` | Checks how each query and target entry parses, against its pocket method. No network. | `QUERY`, `TARGET` | starts `failed_entries.json` afresh |
| `fetch_structures` | Downloads the structures and a bundled Foldseek database the entries need. Also works on its own: `fetch_structures 4Q5J P24941 --out_dir ./structures`. | `ENTRY ...`, else the entries | the cache, or `--out_dir`; optionally `--structures_tsv_path` |
| `align` | Aligns the query chains against the target chains. | the entries, the cache | `alignment.tsv` |
| `pockets` | Downloads the PISA interfaces the entries need, then builds the pocket of every entry. Against the `pdb` database, also downloads each hit's PISA interfaces and structure and builds a pocket per interface. | the entries, the cache, `alignment.tsv` (database target only) | `pockets.json` |
| `compare` | Compares the pockets of every aligned pair. | the target entry, `alignment.tsv`, `pockets.json` | `pocket_comparison.tsv`, `unknown_ids.json`, `incorrect_mapping.json` |
| `superpose` | Superposes the top targets onto each query. | the entries, the cache, `pocket_comparison.tsv`, `alignment.tsv`, `pockets.json` (`pdb` database target only) | `aligned_structures/` |

- **Every step parses the entries again**, from the `query`, `target` and pocket methods in its job
  file, and rejects what `parse` rejects. `fetch_structures` resolves only each entry's structure,
  ignoring its chain and residue parts, and takes entries of its own too. No step reads another's records. So every step after
  `parse` needs a job file naming the entries: `parse` writes one, `parse_settings.json`, and so does
  every step run on its own (`<command>_settings.json`), each a valid job file for any later step.
- **Every step but `parse` and `fetch_structures` requires `--results_dir`**, here or in a job file,
  since its inputs default to files in it. Without one, `parse` and `fetch_structures` write no log,
  settings or `failed_entries.json`, only a file whose own path is given (and `fetch_structures` its
  scratch under `<cache_dir>/tmp`). An option without `_path` (`--alignment`, `--pockets`) names a file the step reads,
  one ending in `_path` a file it writes; both default to the names above under `--results_dir`.
- **The cache is just a setting.** Every step takes the [cache options](#cache-options); a chain
  shares a cache by taking them from the same job file. Entries, entries files and relative paths
  resolve against `--work_dir`, which the settings record absolute, so the commands after `parse`
  can run from any working directory.
- **Only `fetch_structures` downloads the entries' structures**, unless asked. `align`, `pockets` and
  `superpose` skip an entry whose structure is not in the cache, listing it as `structure_not_found`;
  with `--fetch_missing 1` they download it instead (and `align` a missing bundled database). `pockets`
  downloads PISA interfaces and, against the `pdb` database, its hits' structures, which only it needs.
- **`compare` stops if `alignment.tsv` names a chain `pockets.json` does not list**, rather than
  silently returning no rows for it — say an entry's local file changed between `align` and
  `pockets`. An entry `pockets` listed but could not build is skipped.
- **`superpose` reads the aligner off `alignment.tsv`**, so `--align_struct_method auto` resolves as it
  would have in `search` without restating `--aligner`.
- **Each step skips what it cannot use** and lists it in `failed_entries.json` with the step and the
  reason; `parse` starts that file afresh.
- **A `search` writes the same result files** and its settings, `job_settings.json`, so any step can
  be rerun on a finished search:

```bash
pocketmapper superpose --job_file ./out/job_settings.json --align_count 3 --aligned_structure_dir ./out/top3
```

Every step takes `--job_file`, `--verbosity`, `--log_path` (appended to, so a chain writes one log),
`--work_dir`, `--results_dir`, `--job_settings_path` and `--failed_entries_path`; every step after
`parse` also takes the [cache options](#cache-options). Beyond those, each takes the options below;
they mean what they do for `search`.

| Command | Options |
| --- | --- |
| `parse` | `QUERY`, `TARGET` (or from the job file), `--query_pocket_method`, `--target_pocket_method`, the [cache options](#cache-options) |
| `fetch_structures` | `ENTRY ...` (entries or bare ids such as `4Q5J`, `P12345`, `pdb`, or files of them; default the job file's query and target), `--struct_type` (`auto`, `pdb`, `alphafold` or `foldseek_db`), `--out_dir` (every file into one directory instead of the cache), `--structures_tsv_path` (a table of `id`, `struct_type`, `path`, `ok`), `--threads`, the [temp options](#temp-options) |
| `align` | `--threads`, the [temp options](#temp-options), `--aligner`, `--alignment_path` and `--fetch_missing` |
| `pockets` | `--alignment`, `--pockets_path`, `--pisa_source`, `--fetch_missing` |
| `compare` | `--alignment`, `--pockets`, `--pocket_comparison_path` |
| `superpose` | `--pocket_comparison`, `--alignment`, `--pockets`, the [aligned structure options](#aligned-structure-options), `--aligned_structure_dir`, `--threads`, `--fetch_missing` |

Each step writes the settings it ran with to `<results_dir>/<command>_settings.json` (move it with
`--job_settings_path`), a job file any later step takes; `search` writes only `job_settings.json`. See
[Advanced options](#advanced-options).

### Input format

Query and target entries are colon-separated:

```text
STRUCTURE[:CHAIN[:RESIDUES]]
```

- **STRUCTURE** — a 4-character PDB ID (`4Q5J`), a UniProt accession (`P04637`, fetched from
  AlphaFold), or a path to a local structure file.
- **CHAIN** — the chain the pocket belongs to. For the interface-based methods this is two chains joined
  by an underscore (`B_F`), where the first is the chain carrying the pocket and the second is its
  binding partner. Optional; omitted, it defaults to chain `A`.
- **RESIDUES** — optional comma-separated author residue numbers (`10,11,12`). Omitting the whole pocket
  part makes the entry an open search — see below.

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

Which methods an entry can reach depends on its structure type:

| Structure type | Methods available |
| --- | --- |
| PDB entry | `pisa`, `passthrough`, `vdw`, `whole_chain` |
| Local file | `passthrough`, `vdw`, `whole_chain` |
| AlphaFold accession | `passthrough`, `whole_chain` |

PISA is only available for PDB entries, so a local file with an interface-style chain spec
(`4Q5J.cif.gz:B_F`) resolves to `vdw`, not `pisa`. An AlphaFold model is a single chain, so the
methods that need a binding partner are not offered for one at all. A passthrough entry needs an
explicit residue list — without one the entry is an open search instead.

The inferred method (`auto`, the default) can be overridden with `--query_pocket_method` /
`--target_pocket_method`, which does not relax any of the above: the entry must still spell what the
method reads, and the method must be available for that structure type. An entry that cannot satisfy
both is skipped with a warning naming what was missing, before anything is downloaded, and the run
fails only if it leaves a side with no usable entries — so one bad line in a batch file does not stop
the rest. An unrecognised method name is rejected outright.

A passthrough entry's residue ids are checked against the structure: an entry naming a residue the chain
does not have, or has without a CA atom, is skipped with a warning rather than compared, and the rest of
the run continues. A residue id listed more than once is collapsed to one, also with a warning.

For batch runs, pass a path to a file containing one such entry per line instead of a single entry.
Query and target files are read independently, and every query is compared against every target.

### Closed vs Open Search

The two shapes of entry ask different questions.

A **closed search** names a pocket on both sides (`4Q5J:B_F` against `4Q5J:A_E`) and asks *how do these two
pockets compare?* Every column of `pocket_comparison.tsv` is populated: both pockets are described, the
overlap is scored, and `jaccard_index` sizes the shared residues against their union.

An **open search** names a structure but no pocket (`4Q5J:B`, or just `4Q5J`) and asks *does the query
pocket resemble anything on this chain at all?* The whole chain becomes the pocket, so the answer is
carried by `overlap_count` and `query_overlap_ids` — how many of the query pocket's residues the target
chain covers, and which ones.

Because there is no pocket on the target to describe, the descriptive and length-normalised columns
(`target_res_ids`, `target_len`, `target_seq`, `target_pct_aln`, `jaccard_index`) are left empty —
the same shape a Foldseek-database row has. `jaccard_index` in particular needs a second pocket to size
against, and a whole chain's length would swamp the union it normalises by. The overlap itself is still
fully reported: `target_overlap_ids` gives the author residue numbers the query pocket maps onto, and
because a named structure has real coordinates the superposition columns (`rmsd`, `ca_dists`, the
transforms) are populated too, which a Foldseek-database row cannot offer.

Open and closed targets can be mixed freely in one batch file.

### Databases

A target can name a Foldseek database instead of a structure, which searches the query pocket against
every entry in it. Both bundled databases require the foldseek binary, and
`--align_struct_method pocket` is rejected against either (there is no second pocket to fit to).

**`human_domains`** ships inside the package — no download. Its entries are structural domains carved out
of human UniProt sequences. There is no interface to compute on a domain, so the target side behaves like
an open search: no `target_*` descriptors, no `jaccard_index`, no superposition columns. The database
ships an offset table, so `target_overlap_ids` is reported in **UniProt residue numbers** rather than
positions within the domain. The query side (`query_overlap_ids`) is always author residue numbers.

**`pdb`** is downloaded on first use (`foldseek databases PDB`) into `<cache_dir>/fsdb/pdb`. Hits are
real PDB chains, so PocketMapper fetches PISA data for each one and compares against a real interface
pocket; hits with no usable PISA data are dropped rather than compared against a stand-in. **This is
slow the first time**: `4Q5J:B_F` returns roughly 4,970 hits across 3,620 entries. The default
`--pisa_source ftp` downloads their PISA files concurrently; with `--pisa_source api` each assembly is a
separate rate-limited request, so the first run takes hours. Both counts are logged before the fetching
starts, and reruns are fast from the interface cache.

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
| `job_settings.json` | The fully resolved settings for the run. A step run on its own writes `<command>_settings.json` instead. |
| `info.log` | The run log. |
| `pockets.json` | Every entry's pocket, keyed by pocket id (`null` where it could not be built), and `chains`: each aligned chain's name -> the pocket ids on it. Versioned; a file from an older version is rejected, so rerun `pockets`. |
| `failed_entries.json` | Every entry skipped along the way: the entry, the step that skipped it, the reason (`invalid_entry`, `structure_not_found`, `structure_preprocessing_failed`, `pocket_not_built`) and the parsed entry. Empty when nothing was skipped. |
| `unknown_ids.json` | Residue codes the aligner and the structure disagreed on (e.g. MSE -> M). Only written if any were seen, beside `pocket_comparison.tsv`. |
| `incorrect_mapping.json` | Pockets dropped because their own sequence disagreed with the aligner's for that chain — typically assembly vs asymmetric unit numbering. Only written if any were dropped, beside `pocket_comparison.tsv`. |

Under `--cache_dir` and surviving between runs: the downloaded mmCIF files (`pdb_structures/` and
`alphafold_structures/`), the raw PISA responses (`pockets/pisa/`), and the derived
pockets themselves (`pockets/pisa_pockets.json`, `passthrough_pockets.json`, `vdw_pockets.json`,
`whole_chain_pockets.json`). The pocket files are written for inspection only and never read back;
each pocket is an object of metadata fields plus a `residues` map keyed by author residue number.

#### `pocket_comparison.tsv` columns

The `query_*` columns describe the query pocket, the `target_*` columns the target. The table always
has all of these columns: a comparison that stops early — no overlapping residues, no coordinates,
fewer than three residues to superpose, or an open target — leaves the remaining fields **empty**
rather than dropping them.

##### Identity

| Column | Meaning |
| --- | --- |
| `query`, `target` | The two pocket ids, in the input-entry form (`4Q5J:B_F`). |
| `evalue` | Alignment E-value for the underlying chain pair. `-` with the local aligner. |
| `lddt` | Alignment lDDT reported by Foldseek. `-` with the local aligner. |

##### Pocket description

| Column | Meaning |
| --- | --- |
| `query_res_ids`, `target_res_ids` | The pocket's author residue numbers, comma-separated. |
| `query_len`, `target_len` | Number of residues in the pocket. |
| `query_seq`, `target_seq` | The pocket's residues as single-letter codes, in `res_ids` order. |
| `query_pct_aln`, `target_pct_aln` | Fraction of the pocket's residues that fall inside the aligned region at all. A low value means the pocket sits largely outside the alignment. |

##### Overlap

| Column | Meaning |
| --- | --- |
| `overlap_count` | How many alignment positions both pockets occupy. The headline number. |
| `query_overlap_ids` | The query residues in the overlap, as author residue numbers. |
| `target_overlap_ids` | The target residues they map onto, in the same order. UniProt numbers for a `human_domains` target; see [Databases](#databases). |
| `jaccard_index` | `overlap_count` divided by the union of the two pockets. Empty for an open or database target. |

##### Overlap scoring (BLOSUM62)

| Column | Meaning |
| --- | --- |
| `query_seq_overlap`, `target_seq_overlap` | The two overlap sequences, aligned position for position. |
| `overlap_identity` | Fraction of overlap positions where the two residues are identical. |
| `overlap_similarity_binary` | Fraction of overlap positions with a positive BLOSUM62 score — how much of the overlap is conservatively substituted, ignoring how strongly. |
| `overlap_similarity_1_2` | Mean BLOSUM62 score over the overlap, normalised per position by the query residue's self-score, so an identical pair scores 1.0. |
| `overlap_similarity_2_1` | The same, normalised against the target residue instead. The two differ because self-scores differ between residues. |
| `min_overlap_similarity`, `max_overlap_similarity` | The smaller and larger of the two directional scores. |

##### Superposition

Populated only when both pockets have coordinates and at least three residues overlap.

| Column | Meaning |
| --- | --- |
| `target_to_query_u`, `target_to_query_t` | Rotation (nine values) and translation (three) that put the target pocket onto the query pocket, fitted on the overlapping CA atoms. Comma-separated at three decimals, like `alignment.tsv`'s `u` and `t`. |
| `query_to_target_u`, `query_to_target_t` | The same fit in the other direction. |
| `rmsd` | RMSD of the overlapping CA atoms after superposition, in Å. |
| `ca_dists` | Per-residue CA distance after superposition, comma-separated, in overlap order. |

The rotation matrices are stored in Biopython's right-multiplying convention
(`coords · u + t`). If you feed one to gemmi, which left-multiplies, transpose it first — or use
`pocketmapper.pocket_comparison.parse_pocket_transform`, which does the conversion for you.

### Examples

A single pair, aligned with Foldseek:

```bash
pocketmapper search 4Q5J:B_F 4Q5J:A_E --results_dir ./out_fs
```

The same pair with the local BLOSUM62 aligner, which needs no Foldseek binary:

```bash
pocketmapper search 4Q5J:B_F 4Q5J:A_E --aligner seq --results_dir ./out_local
```

An open search — is this pocket anywhere on chain A of 4Q5J at all?

```bash
pocketmapper search 4Q5J:B_F 4Q5J:A --results_dir ./out_open
```

Searching a chain against the bundled Foldseek DB of human domains (requires the foldseek binary):

```bash
pocketmapper search 4Q5J:B_F human_domains --results_dir ./out_hd
```

Batch mode using files with one entry per line:

```bash
pocketmapper search queries.txt targets.txt --job_file job.json
```

### Advanced options

**The job file.** `--job_file job.json` takes a flat `{"option": value}` object using the same names
as `search`'s options, minus the leading `--`, plus `query` and `target`. Every command takes one. An
option given on the command line wins; anything it leaves out comes from the job file, then from the
defaults, and unset paths are derived last. Every option can be given either way — the job file is for
keeping a long invocation reproducible, never the only route to a setting. An unrecognised key is an
error rather than being ignored; a key the command does not use is ignored. `query` and `target` must
each come from exactly one place: giving one both positionally and in the job file is an error, so
`pocketmapper search --job_file job.json` is a complete invocation when the file sets both.

```json
{"query": "queries.txt", "target": "human_domains", "align_count": 3}
```

A finished run's `job_settings.json` is itself a valid job file, and the way to rerun any step with
that run's settings: `pocketmapper superpose --job_file out/job_settings.json --align_count 3`. So is
every step's own `<command>_settings.json`. None of them sets `job_settings_path`, and a job file that
does has it ignored, so a step run from another command's settings never writes over them. A dump
sets *every* option, paths included, so a command-line option changes only the setting it names:
`--results_dir` beside it moves no path, and `--alignment_path` moves only the alignment. One thing does
follow `--results_dir` all the same: where `--temp_dir` may be emptied. So to reuse a run's files, give
its `results_dir` or none. With `search`, a new `--results_dir` beside it writes over that run's outputs.

**`--temp_dir` is emptied on the way in and deleted on the way out.** It holds `query_structures/`,
`target_structures/` and `foldseek_tmp/`, each created only when used, so a rerun into the same
`--results_dir` never hands Foldseek what the previous run left behind. Both the emptying and the
deletion are skipped with a warning when `temp_dir` does not resolve to somewhere under `--cache_dir`
or `--results_dir` — a mistyped `--temp_dir` costs you a stray directory, not its contents.
`--delete_tmp 0` keeps it: it holds the per-run inputs the aligner was actually given, which is what
you want when a run returns nothing and you need to see why.

**`--aligner` picks how chains are aligned.** `foldseek`, the default, needs the Foldseek binary on
`PATH` and runnable (checked by running `foldseek -h`); if it is not, the run stops before anything is
downloaded rather than silently switching method. `seq` uses the built-in BLOSUM62 sequence aligner,
which needs no binary but produces no whole-chain transform and cannot search a Foldseek database.

**`--align_struct_method` picks the transform used to write `aligned_structures/`.**

- `foldseek` — Foldseek's whole-chain fit. Puts the two *chains* on top of each other. Needs the binary.
- `pocket` — the fit of the two pockets on their overlapping residues, taken from
  `pocket_comparison.tsv`. Puts the two *pockets* on top of each other, which is usually what you want
  when the chains are otherwise unrelated. Unavailable against a Foldseek database target, where the
  target pocket has no coordinates to fit.
- `auto` (the default) — `foldseek` with `--aligner foldseek`, `pocket` with `--aligner seq`, which
  produces no whole-chain transform at all.

**Using PocketMapper as a library.** `PocketMapper().search(...)` runs the same workflow as the CLI and
writes the same files; results come back through `results_dir`, not as a return value. Each step
command is a function, `pocketmapper.steps.<step>.<step>`, taking its CLI options as keyword arguments;
its `job_file` also takes a dict, such as the settings `search()` returns. The
individual components (`qt_processor`, `downloads.structure_downloader`, `downloads.pisa_downloader`,
`sequence_aligner`, `structure_aligner`, `pockets.pocket_fetcher`, ...) are each usable on their own.
Note that `search()` deletes its temporary directories on the way out.

Logging goes through the `pocketmapper` logger and never touches the root logger, so your own logging
setup decides what is shown. `search()` still writes `info.log` and sets the `pocketmapper` logger's level
from `verbosity` for the length of the call. For the CLI's console format, use
`pocketmapper.lib.format_handler(logging.StreamHandler())` and add the handler to
`logging.getLogger("pocketmapper")`.

The last step is the one you can defer: run with `align_count=0` and no aligned structures are written,
then run `pocketmapper.steps.superpose.superpose(job_file=settings)`. For finer control, hand
`StructureAligner.align_structs` the run's entries, parsed by `pocketmapper.steps.parse.parse_entries`,
and its result files: it takes `query_ids` and
`target_ids` to superpose only part of the run, and `overwrite=False` to leave alone what it has already
written. (Against a Foldseek database target it also needs `fsdb_path`; `steps.superpose` handles that.)

```python
from pocketmapper import PocketMapper
from pocketmapper.steps.parse import parse_entries
from pocketmapper.structure_aligner import StructureAligner

settings = PocketMapper().search("4Q5J:A_E", "4Q5J:B_F", align_count=0, results_dir="results")
sides, _ = parse_entries(
    settings["query"],
    settings["target"],
    settings["query_pocket_method"],
    settings["target_pocket_method"],
    settings,  # the cache directories
    settings["work_dir"],
)

StructureAligner().align_structs(
    query_records=sides["query"],
    target_records=sides["target"],
    pocket_comparison=settings["pocket_comparison_path"],
    alignment=settings["alignment_path"],
    out_dir=settings["aligned_structure_dir"],
    method=settings["align_struct_method"],
    align_count=10,
    query_ids=["4Q5J:A_E"],
)
```

## Contact / Authors

PocketMapper is developed by Lachlan Ellingboe (<Lachlan.Ellingboe@icr.ac.uk>).

Source, issues and feature requests: <https://github.com/lellingboe/pocketmapper>
