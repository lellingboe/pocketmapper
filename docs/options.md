# Command-line reference

Every option of `pocketmapper search` and of the six step commands. For what an entry looks like, see
[Entries and targets](input.md); for what a run writes, see [Outputs](outputs.md).

## The search command

```bash
pocketmapper search QUERY TARGET [options]
```

### Arguments

Both are required, and are given positionally in this order — there is no `--query`/`--target`
spelling. Either may instead be set in a [job file](#job-files-and-aligners), in which case it must be
left off the command line.

| Argument | Summary |
| --- | --- |
| `QUERY` | Query entry, or a file with one entry per line (see [Input format](input.md#input-format)). |
| `TARGET` | Target entry, a file with one entry per line, or a Foldseek DB name (`human_domains`, `pdb`). |

### Options

| Option | Type | Default | Summary |
| --- | --- | --- | --- |
| `-j`, `--job_file` | path | none | JSON file of `{"option": value}`, query and target included. Arguments override it. |
| `-v`, `--verbosity` | int | `3` | Log level: 4=DEBUG, 3=INFO, 2=WARNING, anything else=ERROR. |
| `-a`, `--aligner` | str | `foldseek` | Chain aligner: `foldseek` (needs the binary) or `seq` (built-in BLOSUM62 sequence aligner). |
| `-q`, `--query_pocket_method` | str | `auto` | Query pocket method: `auto` infers it from each entry; `pisa`, `passthrough`, `vdw`, `whole_chain` force it. Each entry is still checked against the method — see [Input format](input.md#input-format). |
| `-t`, `--target_pocket_method` | str | `auto` | As `--query_pocket_method`, for targets; also accepts `foldseek_db`. |
| `-T`, `--threads` | int | one per core | Cap on the cores Foldseek uses. |
| `--pisa_source` | str | `ftp` | Where PISA interfaces are fetched from: `ftp` (EBI FTP server, concurrent) or `api` (PDBe API, one paced request per assembly). Both fill the same cache. |
| `--pisa_bsa_threshold` | float | unset | Build PISA pockets from buried surface: a residue joins when the partner chain buries more than this percent of its accessible surface (at least 0, below 100). Unset, they come from interface bonds. |
| `--offline` | int | `0` | `1` downloads nothing: no structures, PISA interfaces or Foldseek databases. An entry whose structure is not cached is skipped and marked, as a failed download is; a Foldseek database that is not cached is an error. |
| `--refresh_structures` | int | `0` | `1` re-downloads the PDB and AlphaFold structures in `fetch_structures` even when they are cached. PISA interfaces and Foldseek databases are not refreshed. Not with `--offline 1`. |
| `-h`, `--help` | flag | — | Show the help message and exit. |

#### In options

| Option | Default | Summary |
| --- | --- | --- |
| `--local_root` | the current directory | Directory that relative local structure files and user Foldseek database paths resolve against. Every other path, entries files included, resolves against the current directory. Recorded absolute in the settings dump, so a step run later from another directory resolves them the same way. |

#### Aligned structure options

| Option | Type | Default | Summary |
| --- | --- | --- | --- |
| `--align_count` | int | `10` | How many top-scoring targets to superpose onto each query; `0` disables. |
| `--align_struct_method` | str | `auto` | Which transform superposes a target onto its query: `auto`, `pocket` or `foldseek`. |

#### Cache options

What survives a run and makes the next one fast. Each defaults to a location under `--cache_dir`.

| Option | Default | Summary |
| --- | --- | --- |
| `--cache_dir` | `pocketmapper_cache` | Where structures and PISA responses are cached. |
| `--pdb_dir` | `<cache_dir>/pdb_structures` | Cache of fetched PDB structures. |
| `--alphafold_dir` | `<cache_dir>/alphafold_structures` | Cache of fetched AlphaFold structures. |
| `--pocket_dir` | `<cache_dir>/pockets` | Cache of PISA data. |
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
| `--entries_path` | `<results_dir>/entries.tsv` | Where the entries table is written: how each entry parsed, and which step dropped it and why. |
| `--job_settings_path` | `<results_dir>/job_settings.json` | Where this run's resolved settings are dumped. |
| `--log_path` | `<results_dir>/info.log` | Where the run log is written. It is appended to. |

#### Temp options

Scratch space for a single run: emptied before the run uses it and deleted when it ends successfully (a
failed run leaves it in place) — set `--delete_tmp 0` to keep it for inspection. The per-run query and
target structures and Foldseek's own scratch all live inside it, so it is only created with
`--aligner foldseek`, or to download a bundled Foldseek database.

| Option | Default | Summary |
| --- | --- | --- |
| `--temp_dir` | `<results_dir>/tmp` | Per-run scratch, emptied before use and deleted at the end of a successful run. |
| `--delete_tmp` | `1` | `1` deletes `--temp_dir` at the end of the run; `0` keeps it. |

## Step commands

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
| `parse` | Checks how each query and target entry parses, against its pocket method, and prints the entries table. No network. | `QUERY`, `TARGET` | `entries.tsv` (also printed, after the log), `parse_settings.json` |
| `fetch_structures` | Downloads the structures and a bundled Foldseek database the entries need. | `entries.tsv` | the cache |
| `align` | Aligns the query chains against the target chains. | `entries.tsv`, the cache | `alignment.tsv` |
| `pockets` | Downloads the PISA interfaces the entries need, then builds the pocket of every entry. Against the `pdb` database, also downloads each hit's PISA interfaces and structure and builds a pocket per interface. | `entries.tsv`, the cache, `alignment.tsv` (database target only) | `pockets.json` |
| `compare` | Compares the pockets of every aligned pair. | `entries.tsv`, `alignment.tsv`, `pockets.json` | `pocket_comparison.tsv`, `unknown_ids.json`, `incorrect_mapping.json` |
| `superpose` | Superposes the top targets onto each query. | `entries.tsv`, the cache, `pocket_comparison.tsv`, `alignment.tsv`, `pockets.json` (`pdb` database target only) | `aligned_structures/` |

- **Only `parse` reads the query and target inputs.** It writes `entries.tsv`, a row per query and target entry
  (see [Outputs](outputs.md#outputs)), and every later step reads its records from there. A step skips the
  entries an earlier step dropped, and marks the ones it drops itself with the step and the reason.
  When a step starts, it clears its own marks and those of every later step, so a rerun gives the same
  table.
- **Only `parse` and `search` take `--results_dir` and `--cache_dir`.** They resolve every path from
  them and write it into their settings, `parse_settings.json` or `job_settings.json`. Every other step
  takes each path it uses from that job file, and an argument moves only the path it names. An option
  without `_path` (`--entries`, `--alignment`, `--pockets`) names a file the step reads, and one ending
  in `_path` names a file it writes. A relative path given to a step resolves against the current
  directory, except that a relative local structure file or user Foldseek database named in an entry
  resolves against `--local_root`. The settings hold absolute paths, so the commands after `parse` can
  run from any working directory.
- **Without `--results_dir`, `parse` writes nothing** but a file whose own path is given. The steps after
  it need its `entries.tsv`.
- **Only `fetch_structures` downloads the entries' structures**, unless asked. `align`, `pockets` and
  `superpose` skip an entry whose structure is not in the cache and mark it `structure_not_found`.
  With `--fetch_missing 1` they download it instead (and `align` a missing bundled database). `pockets`
  downloads PISA interfaces and, against the `pdb` database, its hits' structures, which only it needs.
- **With `--offline 1` no step downloads anything**, whatever `--fetch_missing` says. A structure that is
  not cached is skipped and marked `structure_not_found`, and a PISA pocket whose interfaces are not
  cached is marked `pocket_not_built`, as after a failed download.
- **`compare` stops if `alignment.tsv` names a chain `pockets.json` does not list**, rather than
  silently returning no rows for it — say `pockets` ran on a different entries table than `align`.
  An entry `pockets` listed but could not build is skipped.
- **`superpose` reads the aligner off `alignment.tsv`**, so `--align_struct_method auto` resolves as it
  would have in `search` without restating `--aligner`.
- **A `search` writes the same result files** and its settings, `job_settings.json`, so any step can
  be rerun on a finished search:

```bash
pocketmapper superpose --job_file ./out/job_settings.json --align_count 3 --aligned_structure_dir ./out/top3
```

Every step takes `--job_file`, `--verbosity` and `--log_path` (appended to, so a chain writes one log).
Every step after `parse` also takes `--entries`. Beyond those, each takes the options below; they mean
what they do for `search`.

| Command | Options |
| --- | --- |
| `parse` | `QUERY`, `TARGET` (or from the job file), `--query_pocket_method`, `--target_pocket_method`, `--local_root`, the [cache options](#cache-options), `--results_dir`, `--entries_path`, `--job_settings_path` |
| `fetch_structures` | `--local_root`, `--pdb_dir`, `--alphafold_dir`, `--fsdb_dir`, `--threads`, the [temp options](#temp-options), `--offline`, `--refresh_structures` |
| `align` | `--local_root`, `--pdb_dir`, `--alphafold_dir`, `--fsdb_dir`, `--aligner`, `--threads`, `--fetch_missing`, `--foldseek_preprocessed_structure_dir`, `--alignment_path`, the [temp options](#temp-options), `--offline` |
| `pockets` | `--alignment`, `--local_root`, `--pdb_dir`, `--alphafold_dir`, `--pocket_dir`, `--fsdb_dir`, `--pockets_path`, `--pockets_tsv_path` (a table of `pocket_id`, `chain`, `res_auth_ids`, `method`), `--pisa_source`, `--pisa_bsa_threshold`, `--fetch_missing`, `--offline` |
| `compare` | `--local_root`, `--fsdb_dir`, `--alignment`, `--pockets`, `--pocket_comparison_path` |
| `superpose` | `--local_root`, `--pdb_dir`, `--alphafold_dir`, `--fsdb_dir`, `--pocket_comparison`, `--alignment`, `--pockets`, the [aligned structure options](#aligned-structure-options), `--aligned_structure_dir`, `--threads`, `--fetch_missing`, `--offline` |

Only `parse` and `search` write their settings, as a job file for the steps after them. See
[Job files and aligners](#job-files-and-aligners).

## Job files and aligners

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
`parse_settings.json`. Neither sets `job_settings_path`, and a job file that does has it ignored, so a
command run from another's settings never writes over them. A dump sets *every* option, paths
included, so a command-line option changes only the setting it names: `--results_dir` beside it moves
no path, and `--alignment_path` moves only the alignment. With `search`, a new `--results_dir` beside it
writes over that run's outputs.

**`--temp_dir` is emptied on the way in and deleted on the way out of a successful run.** It holds
`query_structures/`, `target_structures/` and `foldseek_tmp/`, each created only when used, so a rerun
into the same `--results_dir` never hands Foldseek what the previous run left behind. Only those three
are ever emptied or deleted, and `--temp_dir` itself only once that leaves it empty, so a mistyped
`--temp_dir` costs you a stray directory, not its contents. `--delete_tmp 0` keeps it: it holds the
per-run inputs the aligner was actually given, which is what you want when a run returns nothing and you
need to see why.

**`--aligner` picks how chains are aligned.** `foldseek`, the default, needs the Foldseek binary on
`PATH` and runnable (checked by running `foldseek -h`); if it is not, the run stops before anything is
downloaded rather than silently switching method. `seq` uses the built-in BLOSUM62 sequence aligner (local
alignment, gap open -11, extend -1), which needs no binary but produces no whole-chain transform and cannot
search a Foldseek database. Like Foldseek, it reports only the aligned region of each chain.

**`--align_struct_method` picks the transform used to write `aligned_structures/`.**

- `foldseek` — Foldseek's whole-chain fit. Puts the two *chains* on top of each other. Needs the binary.
- `pocket` — the fit of the two pockets on their overlapping residues, taken from
  `pocket_comparison.tsv`. Puts the two *pockets* on top of each other, which is usually what you want
  when the chains are otherwise unrelated. Unavailable against a Foldseek database target, where the
  target pocket has no coordinates to fit.
- `auto` (the default) — `foldseek` with `--aligner foldseek`, `pocket` with `--aligner seq`, which
  produces no whole-chain transform at all.
