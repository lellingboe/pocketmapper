# Stateless steps

Design draft. Goal: steps stop handing each other record JSONs, and `fetch_structures` and `pockets`
become useful on their own, outside a stage-by-stage pipeline run. Nothing here is implemented.

## Today

`parse` turns the query and target inputs into `query_records.json` / `target_records.json` (lists of
`QTRecord` dicts) and writes `cache_dirs.json`. Every later step reads the records files, and
`fetch_structures`, `align` and `pockets` rewrite them in place, dropping the records they fail. A
standalone step therefore needs a results directory that `parse` already populated.

## Principle: the entry string is the record

Every `QTRecord` field is a pure function of these inputs (`QTProcessor.parse_individual_qt`):

- the entry string (`struct_info[:chain_info[:residue_info]]`), kept verbatim as `pocket_id`;
- the side's pocket method setting (`auto` or forced);
- the cache dirs (`pdb_dir`, `alphafold_dir`, `fsdb_dir`), which give `struct_path`;
- a local file's bytes, hashed into `preprocess_name`;
- the directory a local file, an entries file, a user FSDB and a relative cache dir resolve against,
  and the filesystem itself (`determine_struct_type` calls `os.path.isfile`).

Today parse makes every path absolute once and stores it in the records and the manifest. Without
those files, a step run from another directory must still resolve paths as the first step did; the
`work_dir` setting pins that directory (see Settings).

Only two pieces of record state are not derivable from the entry, and both come from align against
a Foldseek database:

- `fsdb_pockets` (`"pisa"` or `"whole_chain"`): decided by whether any hit name in the alignment
  parses as a PDB entry (`lib.parse_foldseek_pdb_entry_name`). Derivable from `alignment.tsv`.
- The expanded pisa target records (`steps.align.expand_fsdb_pdb_targets`), whose
  `preprocess_name` is the Foldseek hit name. Derivable from `alignment.tsv` plus the PISA cache.

So the records files and `cache_dirs.json` carry nothing a step cannot recompute. Proposal: every
step re-derives its records from the entries, in memory, and hands on only computed artefacts.

## Files

| File | Fate |
|---|---|
| `query_records.json`, `target_records.json` | Gone. Steps re-derive records from the entries. |
| `cache_dirs.json` | Gone. Cache dirs are ordinary settings on every command that touches the cache. |
| `alignment.tsv` | Unchanged. |
| `pocket_comparison.tsv` (+ `unknown_ids.json`, `incorrect_mapping.json`) | Unchanged. |
| `pockets.json` | Gains a version key and the alignment join map (below). |
| `failed_entries.json` | Kept as a log. Started fresh by parse and search; standalone steps append. Nothing reads it back. |
| `job_settings.json` | Kept for search. Every other command writes `<command>_settings.json`; a chain passes `parse_settings.json` on as `--job_file`. `job_settings_path` is never read from a job file (see Settings). |
| Cached structures, preprocessed chains, PISA data | Unchanged. The structure cache is the real hand-off from fetch to align/pockets. |

Also gone: `QTRecord.fsdb_pockets`, `records.synthesise_target_pockets`, and every in-place records
rewrite.

### `pockets.json` v2

Today: `{pocket_id: Pocket}`. `compare` builds `preprocess_name -> [pocket_id]`
(`records.preproc_to_ids`) from the records files, because the alignment is keyed by
`preprocess_name` and the pockets by `pocket_id`. Without records files, `pockets.json` carries the
join itself:

```json
{
  "version": 2,
  "pockets": {"4Q5J:B_F": {...Pocket...}, "4Q5J:A:9999": null},
  "chains":  {"4Q5J_B_<md5>": ["4Q5J:B_F"], "4q5j-assembly1_B": ["4Q5J:B_F"], "1abc-assembly1_A": []}
}
```

- `load_pockets` rejects a missing or different `version` with "rerun pockets", so a v1 file gets a
  clear error rather than a generic read failure.
- `null` = pockets was given the entry and did not build it, for any reason: the build failed, or
  its structure is missing. Absent = pockets was never given it. Every entry pockets is given, built
  or not, has its `preprocess_name` in `chains`. `dump_pockets` already writes None, but
  `build_pockets` filters Nones out before dumping (today's files hold none); it stops filtering.
- `chains` covers structure entries and expanded Foldseek PDB hits alike, so `compare` needs no
  records at all.
- **Every hit of a PDB-named Foldseek database is a `chains` key.** `[]` for a hit whose name does
  not parse as a PDB entry, has no PISA interface, or has only unspellable chain pairs; a hit whose
  structure download failed lists its pocket_ids, with `null` pockets. So absent means only
  "pockets never saw this name". Today these hits are dropped and `resolve_pockets` skips them
  silently; the coverage check in compare would otherwise fail every PDB-FSDB run. A database whose
  hits are not PDB-named (synthesis mode) adds no hit to `chains`.
- A wrapper rather than a field on `Pocket`: `Pocket` stays free of alignment concerns, and a field
  would need to be a list, since one pocket can sit behind two hit names (`4q5j-assembly1_B`,
  `-assembly2_B`).
- `res_auth_ids` order and `seq_pos` round-trip as today.

## Commands

Every command keeps taking a job file, with arguments winning.

**Network.** Only `fetch_structures` downloads the structures that entries name (and bundled
Foldseek databases). align, pockets and superpose skip an entry whose structure is missing, as
`structure_not_found`, with a hint to run `fetch_structures` (omitted under search, which already
ran it) — unless `--fetch_missing` is given, which makes them download it through the same helper.
pockets always downloads its own inputs, PISA data and Foldseek PDB-hit structures, since no
earlier step can fetch them.

### `fetch_structures`

```
fetch_structures ENTRY... [--struct_type auto|pdb|alphafold|foldseek_db] [--out_dir DIR] [--structures_tsv_path PATH]
```

- `ENTRY` is a bare id (`4Q5J`, `P12345`, `pdb`), a full entry (chain and residue parts ignored), or
  a file of them. With no `ENTRY`, uses the job file's query and target entries together. A side
  whose pocket method is `foldseek_db` is fetched as struct type `foldseek_db`: a user database path
  otherwise passes `isfile` and would parse as a local file.
- Entries go through a structure-only parse (`determine_struct_type`, `determine_ref_struct_path`,
  the bundled-database lookup), not `build_record`, which rejects an entry whose pocket part fits no
  method (e.g. `P12345:A_B`).
- `--struct_type` forces the type, as `--*_pocket_method` forces the method today; `auto` infers.
  Needs a forced-type parameter in `QTProcessor` beside the forced pocket method.
- `--out_dir` writes every file to one directory (`<out_dir>/<ID>.cif.gz`, or the database).
  Default: the cache dirs, so a standalone fetch warms the pipeline's cache. Implemented as a
  `QTProcessor` with `pdb_dir`, `alphafold_dir` and `fsdb_dir` all set to `out_dir`; the downloader
  already writes to each record's `struct_path`.
- `--structures_tsv_path`: optional TSV of `id, struct_type, path, ok`, for use outside the
  pipeline. A local file is listed, not fetched.
- `results_dir` is optional. Without it, no log file, settings dump or `failed_entries.json` is
  written unless its path is given. Needs `resolve_paths` to leave results paths unset when there is
  no `results_dir` (today it fills in a timestamped one) and `log_to_file` to do nothing without a
  log path.
- Failures are logged and appended to `failed_entries.json`; nothing is rewritten.

### `pockets`

```
pockets ENTRY... [--pocket_method M] [--alignment ALN --target FSDB] [--fetch_missing 0|1]
        [--pockets_path P] [--pockets_tsv_path T]
```

- Builds the pocket of each entry into `pockets.json` v2. With no `ENTRY`, uses the job file's
  query and target entries, each side with its own pocket method. `ENTRY` and `--pocket_method` are
  job keys (`entries`, `pocket_method`), not `Settings` fields, so the dump replays a positional run.
- A missing entry structure is skipped unless `--fetch_missing` (see Network); the entry is
  recorded as `null`.
- **Foldseek PDB-hit expansion moves here from align.** Given an FSDB target and an alignment, it
  reads the hits, downloads their PISA data and structures, and builds one pisa pocket per
  interface, keyed in `chains` by hit name. align then never touches PISA or hit structures. An FSDB
  target without an alignment is an error.
- `--pockets_tsv_path`: optional human-readable table (`pocket_id, chain, res_auth_ids, method`).

### `align`

Inputs: query and target entries, cache options, aligner, `--fetch_missing`. Output:
`alignment.tsv`. Re-parses entries, skips (or with `--fetch_missing` fetches) missing structures and
a missing bundled database, preprocesses chains, aligns. No PISA, no records written; loses
`--pisa_source`. An FSDB target's path comes from its entry.

### `compare`

Inputs: `alignment.tsv`, `pockets.json`, the target entry and the cache options. Output:
`pocket_comparison.tsv`.

- The target entry says whether the target is an FSDB (lookup vs synthesis) and locates its offset
  table; resolving a bundled name needs `fsdb_dir`.
- `preproc_to_ids` comes from `pockets.json` `chains`.
- Synthesis (whole-chain target pockets) = FSDB target and no hit is PDB-named. One shared helper,
  e.g. `fsdb_pocket_mode(alignment_df)`, replaces the `fsdb_pockets` field; pockets, compare and
  superpose all call it.
- **Coverage check moves to the alignment.** An alignment name with no `chains` entry is an error
  ("run pockets on X"); a `null` pocket or an empty list is a skip. Since pockets lists every entry
  it was given, the error fires only for a name pockets never saw. Synthesis mode is exempt: its
  target names have no `chains` entry. Today `preprocess_name` is fixed at parse and carried in the
  records, so align and compare always agree. Re-hashing local files in every step breaks that: a
  file edited between steps gets a new name, and `resolve_pockets` would silently give zero rows.
  This check turns that into an error.

### `superpose`

Inputs: query and target entries, `alignment.tsv`, `pocket_comparison.tsv`, `pockets.json` (FSDB PDB
hits only), cache options (structure paths and the FSDB path), `--fetch_missing`. Output: aligned
structures. Records re-derived from entries; an expanded hit's `pocket_id -> hit name` comes from
`chains`.

### `parse`

Becomes a dry run: reports how each entry resolves (struct type, pocket method, path, or why it is
rejected) as a TSV. Nothing downstream reads it.

- `parse` returns the table and `cli()` prints it to stdout after the run. `cli()` sends the log to
  stdout too, so printing during the run would interleave the two. `--entries_path` also writes it
  to a file.
- `results_dir` is optional, with the same `resolve_paths` / `log_to_file` note as
  `fetch_structures`. Given one, as the head of a chain, parse starts `failed_entries.json` fresh and
  writes `parse_settings.json`, the `--job_file` the rest of the chain takes.
- `search` still runs it first, so a bad entry or an empty side fails before any download. Under
  `search` the table is logged, not printed.

### `search`

Same order, same job dict passed to each step. Loses the records and manifest writes. It runs
`fetch_structures` itself, so it leaves `fetch_missing` off. It logs parse's table. The FSDB guards
(`check_fsdb_aligner`, `check_fsdb_align_struct_method`) read the target entry instead of a records
file.

## Settings and job files

- **Job keys are separate from `Settings`.** `Settings` stays search's resolved configuration. A job
  key set, the union of every command's options, is what a job file may hold and what
  `layer_settings` seeds from (today: `fields(Settings)`). `read_job_file` checks keys against it
  instead of `fields(Settings)`; unknown keys stay rejected. A dump from before phase 2 (it names
  `query_records_path`) is then rejected, with a message naming the removed keys. Step-only options (`entries`, `pocket_method`, `struct_type`, `out_dir`,
  `structures_tsv_path`, `pockets_tsv_path`, `entries_path`, `fetch_missing`) are job keys only, so
  search takes none of them. `Settings` drops `query_records_path` / `target_records_path` and gains
  `work_dir`. The project-overview checklist for a new option changes with it: a step-only option
  goes in the job key set, not `Settings`.
- **A chained step gets its entries and cache dirs only from `--job_file` or arguments.** Nothing in
  `results_dir` supplies them. A chain:
  `parse Q T --results_dir R ; fetch_structures --job_file R/parse_settings.json ; align --job_file R/parse_settings.json ; …`.
- **A dump is the layered input.** Every command writes what it received, job file then arguments
  then defaults, resolved, to `job_settings_path`, leaving out `job_settings_path` itself. Any dump
  in a chain is therefore a valid job file for any later step, and each step is rerunnable from its
  own dump. Positional entries carry over: `pockets --job_file R/fetch_structures_settings.json`
  builds only the entries fetch was given.
- **`job_settings_path` is never read from a job file.** Each command's default is
  `<results_dir>/<command>_settings.json` (`job_settings.json` for search), unless passed as an
  argument. Otherwise a step chained off a dump would overwrite it: `compare --job_file
  search/job_settings.json` (`test_steps_15`) would replace search's dump. The per-command default
  lives outside `RESULTS_PATH_DEFAULTS`, which holds one default per key.
- **`work_dir`**: default the working directory, dumped absolute. Entries, the paths inside an
  entries file, a user FSDB and every relative path setting resolve against it, not against the
  step's cwd; `pocket_id` stays as typed. A step run elsewhere from a dump then resolves as the first
  step did, including a local-file entry typed on the command line (`4Q5J.cif.gz:B_F`), which no
  `abspath` of the `query` string could fix without changing its `pocket_id`.
- `search` writes only `job_settings.json`; the steps it runs write no dump. The nesting works like
  `lib.temp_dir_scope`: a scope already held by an enclosing call does nothing.
- Every command that touches the cache takes the cache options. With `fetch_missing` 0, a cache
  split across a chain skips entries as `structure_not_found`, or fails when a side ends up empty;
  an FSDB fetched into one `fsdb_dir` and searched from another fails loudly at align.
- `fetch_missing`: 0 or 1, default 0 (`constants.DEFAULT_FETCH_MISSING`), checked by a
  `settings.resolve_fetch_missing` modelled on `resolve_delete_tmp`.
- `results_dir` stays required by align, pockets, compare and superpose; it is optional for
  `fetch_structures` and `parse`. It remains the default root for outputs and `temp_dir`.

## Failure model

Today: "a records file holds only usable records". Proposed: "each step skips what it cannot use and
logs it".

- `failed_entries.json` is started fresh by `parse` and `search`; standalone steps append. A chain
  rerun in the same `results_dir` therefore starts clean, as search does.
- Every step re-parses the entries, so a rejected entry would be logged as `invalid_entry` once per
  step. `append_failed_entries` skips an entry whose `(pocket_id, reason)` is already in the file,
  and the step logs it at DEBUG rather than WARNING. Each rejection and missing structure therefore
  warns once per run, in a chain as under search. `parse_individual_qt` warns itself today; the
  warning moves to the caller, which knows whether the failure is new.
- A missing structure is a skip, not a fetch, outside `fetch_structures` (see Network), so a step
  without `--fetch_missing` never reaches the network for an entry's structure.
- No step rewrites its inputs, so every step is idempotent on rerun. (Today align must keep only the
  `foldseek_db` record on entry to rerun on its own output; that special case goes.)

## Library surface

- Core functions take record lists in memory, not records paths.
- `parse_entries(query, target, query_pocket_method, target_pocket_method, cache_dirs, work_dir)`
  wraps `QTProcessor` for library callers and every step. It owns the side checks now in
  `steps.parse.parse_inputs` (empty side, FSDB not the only target, FSDB query entry), which need
  both sides, so every step rejects what parse rejects. It is built on a one-side helper taking a
  list of entries or an entries file, which positional `ENTRY...` uses; `QTProcessor` takes one
  string or one file today.
- `records.py` holds no files once the records files go; rename it (e.g. `entries.py`).
- `pockets` and `fetch_structures` cores are then directly usable: entries in, files out.

## Invariants touched

- `preprocess_name` is recomputed per step (local files hashed each time). Consistency now rests on
  the inputs not changing between steps; the new coverage check in compare catches a mismatch.
- Entry paths resolve against `work_dir`, dumped absolute, not against each step's cwd.
- `pockets.json` round-trip of `res_auth_ids` order and `seq_pos`: unchanged.
- Open-search filter in `StructureAligner.align_structs` (query-only `pocket_id` in `target`):
  unchanged; `chains` mixes sides exactly as `preproc_to_ids` does today.
- One `pocket_id` behind two hit names (`4q5j-assembly1_B`, `-assembly2_B`): `chains` lists both,
  same as today.

## Suggested phasing

Each phase leaves the e2e suite green.

1. Move FSDB expansion from align to pockets; `pockets.json` gains `version` and `chains` (every
   PDB-database hit listed, empty included; every entry given recorded, `null` if not built).
   compare takes `preproc_to_ids` from `chains` from this phase on, since the expanded records exist
   nowhere else. compare and superpose derive the FSDB mode from the alignment. Records files still
   exist.
2. Steps take entries and re-derive records; drop the records files and `cache_dirs.json`. The job
   key set, `work_dir`, per-command dumps of the layered input, the `job_settings_path` rule and the
   nesting scope; chains pass `--job_file`. `parse_entries` owns the side checks. align and superpose
   skip missing structures as `structure_not_found` themselves (today they rely on fetch having
   dropped those records); failed-entry de-duplication and DEBUG on repeats.
3. Generic `fetch_structures` / `pockets` interfaces: positional entries, `--struct_type`,
   `--out_dir`, `--structures_tsv_path`, `--pockets_tsv_path`, `--fetch_missing`; `parse` as a dry
   run with its table printed by `cli()`.

## Also changes

- `cli.OPTIONS` / `COMMANDS`, README options, step table and chain examples. `cli()` prints parse's
  returned table.
- e2e chained cases (`test_steps_*`): every step after parse gains
  `--job_file @OUT@/parse_settings.json`; cases asserting records files or `same=` on them change.
  `test_steps_8` (`cd @OUT@` mid-chain) passes through `work_dir`: `testfile.txt` line 2 is the
  relative local file `4Q5J.cif.gz:B_F`, which would otherwise not resolve after the `cd`.
  `test_steps_15` writes `compare_settings.json`, leaving search's dump untouched. `test_steps_9`
  (align rerun on its own output) and `_13` (target file align never rewrote) are obsolete. `_10`
  fails through the coverage check instead. `_11` / `_12` should hold for align too.
- project-overview: hand-off files, "target side's shape", Settings (job key set, `work_dir`,
  per-option checklist).

## Decisions

| Question | Decision | Why |
|---|---|---|
| Missing structure in align / pockets / superpose | Skip as `structure_not_found`; opt-in `--fetch_missing` fetches | Network stays in one step by default; the flag keeps standalone use convenient |
| Scope of `--fetch_missing` | Entry structures and bundled databases only; pockets always fetches PISA and hit structures | Those are pockets' own inputs; no earlier step can fetch them |
| `parse` | Kept as a dry run | Cheap validation before any download |
| parse table | Returned to `cli()`, printed after the run; `--entries_path` also writes it | `cli()` logs to stdout, so printing during the run would interleave |
| parse table under search | Logged, not printed | search's stdout stays the run log |
| Join map | `chains` wrapper in `pockets.json`, every PDB-database hit listed | `Pocket` stays alignment-free; one pocket can have two hit names; absent means unseen |
| `null` pocket | Given to pockets and not built, for any reason | compare errors only on names pockets never saw |
| pockets / fetch_structures inputs | Positional `ENTRY...`, job file's query + target when none | Reads naturally standalone; the pipeline still supplies both sides |
| Where options live | Job key set separate from `Settings`; unknown keys rejected | `Settings` stays search's; search takes no step-only options |
| Settings dumps | `<command>_settings.json` per step, holding the layered input minus `job_settings_path`; `search` writes only `job_settings.json` | Any dump in a chain feeds any later step; each step rerunnable from its own dump |
| `job_settings_path` | Never read from a job file | A step chained off a dump would overwrite it |
| Chain inputs | Explicit `--job_file <dump>` on every chained step | Stateless; no hand-off file in `results_dir` |
| Relative paths | `work_dir`, default cwd, dumped absolute; entries and path settings resolve against it | A cwd change mid-chain resolves as the first step did; `pocket_id` stays as typed |
| Replaying positional runs | `entries` and `pocket_method` are job keys | Every run replays from its dump |
| Repeat failures | DEBUG when already in `failed_entries.json` | One warning per failure per run, chain or search |
| `results_dir` | Required except by `fetch_structures` and `parse` | Those two are useful with no results directory |
| `failed_entries.json` start | `parse` and `search` | A chain rerun starts clean, as search does |
| Forced struct type | `--struct_type ... foldseek_db` | Same spelling as the record value and pocket method |
| Structures TSV | `--structures_tsv_path` | Parallels `--pockets_tsv_path`; "manifest" retires with `cache_dirs.json` |
| `pockets.json` format | `"version": 2`, checked by `load_pockets` | An old file gets a clear error |
