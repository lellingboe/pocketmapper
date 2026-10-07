# Subcommands: option A

Spec for splitting `search` into one command per step: `parse`, `fetch`, `align`, `pockets`,
`compare`, `superpose`, plus `search` itself. A spec, not an implementation. Code names were checked
against the tree at `6c615ff`.

## Goal and constraints

- Run the pipeline one step at a time, with a useful result after each. A step needs only the
  previous step's output files. A user can write `alignment.tsv` or `pockets.json` by hand; records
  files carry an md5 `preprocess_name`, so they come from `parse`.
- No `self.` hand-off between steps: each command reads and writes named files.
- `search` stays row-identical: same `pocket_comparison.tsv` / `alignment.tsv` rows (compare sorted;
  row order is nondeterministic). Deliberate changes are listed under Behaviour changes.
- `cli.py` stays the only module that knows about argv or exit codes.
- Only `search` takes a job file. Only `parse` and `search` take query/target strings; every other
  command takes records files.

## Current state

`PocketMapper.search()` (`pocketmapper.py`) runs eight steps that share state on `self`. Each piece
of state and what replaces it:

| `self.` state | Set by | Read by | Under option A |
| --- | --- | --- | --- |
| `settings` | `configure_workflow` | every step | `search`: unchanged. Other commands: their own arguments plus the cache manifest. |
| `query_df`, `target_df` | `configure_query_target`; `success`/`failure_reason` mutated by `fetch_missing_structures`, `foldseek_preprocessing`; `target_df` appended by `expand_fsdb_pdb_targets` | align, pockets, compare, superpose | `query_records.json`, `target_records.json`; failed rows go to `failed_entries.json` |
| `fsdb_target` | `configure_query_target` | fetch, preprocessing, `foldseek_alignment`, expand, compare, `align_structs` | Derived: a `foldseek_db` target record is present (Target-side shapes) |
| `fsdb_pdb_target` | `expand_fsdb_pdb_targets` | compare, `align_structs` | Derived: that record's `fsdb_pockets` is `"pisa"` |
| `query_tmp_dir`, `target_tmp_dir`, `foldseek_tmp_dir` | `configure_temp_dir` | `fetch_missing_fsdb`, preprocessing, `foldseek_alignment` | Per command (`fetch`, `align`) |
| `query_db_path`, `target_db_path` | `foldseek_alignment` | only itself | Locals |
| `log_handler`, `previous_log_level` | `configure_logging` | `reset_logging` | A logging context manager per command |
| `pockets` (memory only) | `get_pockets` | `compare_pockets_based_on_alignment` | `pockets.json` |

Already on disk today and unchanged: `alignment.tsv`, `pocket_comparison.tsv`, `unknown_ids.json`,
`incorrect_mapping.json`, `aligned_structures/`, `job_settings.json`, `info.log`.
`pocket_dir/<method>_pockets.json` stays a per-method cache, overwritten per run and never read back.
It is not the hand-off file.

## Layout

- `pocketmapper/settings.py`: `Settings`, `resolve_paths`, and the `resolve_*` validators as module
  functions. `resolve_aligner` is split into value validation and the foldseek probe, since `fetch`
  and `superpose` use the probe alone.
- `pocketmapper/records.py`: read/write records JSON and the cache manifest, append to
  `failed_entries.json`, the target-side derivations, `preproc_to_ids`.
- `pocketmapper/steps/`: `parse.py`, `fetch.py`, `align.py`, `pockets.py`, `compare.py`,
  `superpose.py`. One function per step, explicit arguments, no `Settings`. A step reads its input
  files, does the work and writes its output files.
- `pocketmapper/commands.py`: one function per command. Each resolves path defaults, validates,
  opens logging, manages its temp dir, and calls its step.
- `pocketmapper/pocketmapper.py`: `PocketMapper.search`. It builds `Settings`, opens logging once,
  and calls the steps in order through the same files. `chain == search` therefore holds by
  construction.
- `pocketmapper/cli.py`: one subparser per command, shared option groups via `parents=`, one
  dispatch table (command → function, called with the parsed args minus `command`). This replaces
  `cli()`'s kwarg block.

## Conventions

- **Option names:** a bare name is a file the command reads (`--alignment`, `--pockets`,
  `--query_records`). A `*_path` name is a file it writes (`--alignment_path`, `--pockets_path`),
  matching the existing `*_path` options.
- **Defaults:** input and output defaults resolve under `--results_dir`. An output that rewrites an
  input (`fetch`, `align` records) defaults to that input's path, so a chain in one `results_dir`
  updates the records in place. Pass a `*_path` to keep the old copy. `pockets` always rewrites its
  input records files in place and has no `*_path` for them.
- **`--results_dir`:** `parse` and `search` default to `pocketmapper_results_<YYMMDD_HHMMSS>`. Every
  other command requires it, since its default inputs live under it.
- **Cache manifest** (`<results_dir>/cache_dirs.json`): written by `parse` and `search`. Holds the
  absolute `cache_dir`, `pdb_dir`, `alphafold_dir`, `pocket_dir`,
  `foldseek_preprocessed_structure_dir` and `fsdb_dir`. `fetch`, `align` and `pockets` read their
  cache directories from it and take no cache options, so a chain cannot split its cache. A missing
  manifest is `logger.critical` + `PocketMapperError`. Fixed name, no option, like
  `unknown_ids.json`. `--pisa_source` is not in it: it says how to download, not where to cache.
- **Shared on every command:** `--results_dir`, `-v/--verbosity`, `--log_path`
  (`<results_dir>/info.log`; `FileHandler` appends, so a chain writes one log).
- **`--failed_entries_path`** (`<results_dir>/failed_entries.json`): on every command that can
  drop a record (`parse`, `fetch`, `align`, `pockets`, `search`). `parse` and `search` start a
  chain, so they truncate it on entry; the others append.
- **Settings dump:** only `search` writes one (`job_settings.json`). Other commands write none; the
  cache manifest holds the cache directories a chain shares.
- **Temp:** `--temp_dir`/`--delete_tmp` only on `fetch`, `align`, `search`. Each command empties
  its scratch on entry and deletes it on exit (`lib.is_within`-guarded, as today; `fetch` and
  `align` take `cache_dir` from the manifest). `superpose`'s FSDB rebuild writes under
  `aligned_structure_dir`, as today.
- **Directories:** each is made by its first writer (unchanged rule).
- **Absolute paths:** commands pass `os.path.abspath` of every cache directory to `QTProcessor`, so
  `struct_path` is stored absolute and a chain can run
  from another cwd. Two `struct_path`s are taken from the input as typed and need their own
  `os.path.abspath`: a local file's (`determine_ref_struct_path`) and a non-bundled Foldseek
  database's (the `foldseek_db` branch of `parse_individual_qt`). `pocket_id` and `struct_info` stay
  as typed. This is safe for `preprocess_name`, which hashes the basename and contents, not the path
  (`parse_individual_qt`).

## Hand-off files

| File | Written by | Read by | Format |
| --- | --- | --- | --- |
| `query_records.json`, `target_records.json` | parse; rewritten by fetch, align, pockets | fetch, align, pockets, compare, superpose | JSON list of `QTRecord` dicts |
| `cache_dirs.json` | parse, search | fetch, align, pockets | JSON object, cache directory name → absolute path |
| `failed_entries.json` | parse, search (truncate); fetch, align, pockets (append) | the user | JSON list of failure entries |
| `alignment.tsv` | align | compare, superpose | Unchanged, `ALIGNMENT_COLUMNS` |
| `pockets.json` | pockets | compare | `dump_pockets` format |
| `pocket_comparison.tsv` | compare | superpose | Unchanged, `POCKET_COMPARISON_COLUMNS` |

Records are JSON, not TSV, to keep `None` and bools. A missing input file is `logger.critical` +
`PocketMapperError`.

### `QTRecord`

- Drops `success` and `failure_reason`. **A records file holds only usable records:** `fetch`,
  `align` and `pockets` each drop the records they fail.
- Adds `fsdb_pockets` (`str | None`): how `align` resolved a Foldseek-database target. `"pisa"` for
  a PDB-named database, `"whole_chain"` for any other. Set only on the `foldseek_db` record, and only
  by `align`; None everywhere else. It is what tells an empty expansion (zero rows) from a non-PDB
  database (synthesised rows).

Dropping failed rows does not change output. Today they reach `preproc_to_ids` and
`align_structs`, but they have no pocket, so `resolve_pockets` skips them and they never reach
`pocket_comparison.tsv`. A record `pockets` drops has no entry in `pockets.json`, which
`resolve_pockets` skips the same way. The "chain == search" e2e case confirms it.

### Target-side shapes

Every former `self.` flag is derived from which of these the target records file holds:

| Shape | Target records | Produced by |
| --- | --- | --- |
| Structures | No `foldseek_db` record | parse |
| FSDB, unaligned | Exactly one record, `foldseek_db`, `fsdb_pockets` None | parse |
| FSDB, whole-chain | That record, `fsdb_pockets` `"whole_chain"` | align: a non-PDB database |
| FSDB-PDB, expanded | That record, `fsdb_pockets` `"pisa"`, plus pisa records: `struct_type` `pdb`, `preprocess_name` = Foldseek entry name, preprocess paths None. There may be no pisa records. | align: a PDB-named database |

- `synthesise_target_pockets` ⇔ a `foldseek_db` record is present and its `fsdb_pockets` is not
  `"pisa"`.
- `StructureAligner.align_structs`' `fsdb_path` is the `foldseek_db` record's `struct_path` in both
  FSDB shapes, and None otherwise. With `foldseek`, an expanded hit's transform fits the DB's
  assembly, not the AU in its `struct_path`, so superpose rebuilds the structures from the DB.
  `target_records`:
  - whole-chain: `[]` (ids are entry names);
  - expanded: every target record, as today. `fsdb_target_records` already leaves out the
    `foldseek_db` record, whose `preprocess_name` is None.
- **Expansion appends to the `foldseek_db` record; it does not replace it.** Given a target file
  that already holds expanded records, `align` keeps only the `foldseek_db` record and expands
  again, so a rerun on its own output gives the same records. An empty expansion (`"pisa"`, no pisa
  records) gives zero comparison rows, as today (warning, no error).
- "Target row 0 is the FSDB record" and `parse_foldseek_pdb_entry_name` leave the hand-off
  contract. The latter is used only inside `align`, to decide whether the DB is PDB-named.
- `preproc_to_ids` is built from both records files, exactly as in
  `compare_pockets_based_on_alignment`.

### `pockets.json`

- `dump_pockets` format, written to `--pockets_path` (`<results_dir>/pockets.json`).
- New `load_pockets` beside `dump_pockets` in `pockets/pocket_fetcher.py`. It rebuilds
  `Pocket`/`PocketResidue` and keeps their None fields (pisa, vdw). JSON lists keep `res_auth_ids`
  order, and `residues` keys are already strings.
- A record is `pocket_not_built` (failed_entries.json) when the merged dict (after
  `pockets |= method_pockets`) has no entry for its `pocket_id` or maps it to None. Checked after
  the merge, since a later method can overwrite an earlier one's pocket with None. Its record is
  dropped, and `pockets.json` holds no `null`.

### `failed_entries.json`

Record failures only. `incorrect_mapping.json` and `unknown_ids.json` stay as they are;
`pockets/pisa/errors.json` stays separate (see Behaviour changes).

- An entry is `{"pocket_id", "step", "reason", "source"}` plus the record's fields when there is a
  record. `source` is the records file the record came from, or the query/target input for `parse`.
- `parse` and `search` truncate it on entry, so a rerun into the same `results_dir` does not
  duplicate entries. `fetch`, `align` and `pockets` append.

| Step | `reason` | When |
| --- | --- | --- |
| parse | `invalid_entry` | `parse_individual_qt` rejects the entry. It must return the reason it logs today instead of a bare None; `determine_struct_type` and `parse_residue_info` pass theirs up, and `expand_fsdb_pdb_targets` handles the new return. Also a `foldseek_db` query entry. |
| fetch | `structure_not_found` | `download_missing_structures` reports False |
| align | `structure_preprocessing_failed` | `preprocess_records` reports False (foldseek aligner). It reports only the record it processed for a `(preprocess_name, chain_info)` group; every record in that group is dropped. |
| align | `structure_not_found` | a record's `struct_path` is missing (`fetch` skipped, or its file removed); or an expanded hit whose structure cannot be fetched |
| pockets | `pocket_not_built` | its `pocket_id` is missing or None in the merged pockets. Also removed from the input records file. |

An unrecognised pocket method *name* still raises up front, as today.

## Commands

Options are listed with their defaults. "cache group" etc. are the shared argparse parents.

### `parse QUERY TARGET`

Parse the input grammar into records. No network. "Check my input and the inferred pocket methods."

- **Reads:** query/target strings or files. Both are required.
- **Writes:** `--query_records_path`, `--target_records_path`; `cache_dirs.json`;
  `failed_entries.json` (truncated first).
- **Options:**
  - `-q/--query_pocket_method`, `-t/--target_pocket_method` (`auto`).
  - Cache group: `--cache_dir`, `--pdb_dir`, `--alphafold_dir`, `--pocket_dir`,
    `--foldseek_preprocessed_structure_dir`, `--fsdb_dir`. Baked into the record paths and the
    manifest, so later commands do not take them.
- **Validates:** pocket method names; a side left empty is fatal; a `foldseek_db` target entry must
  be the only target entry; no query entry may be a `foldseek_db` (today accepted, and it does
  nothing).
- **Reuses:** `QTProcessor.process_qt_cmdline_input`, which now returns `(records, rejected)`:
  `QTRecord` dicts and the rejected entries with their reasons.
- **Notes:**
  - The two FSDB checks of `configure_query_target` move out of parse (Validator placement).

### `fetch`

Download everything the records need that is knowable before alignment: structures, the bundled
FSDB, PISA interfaces. "Pre-warm a cache on a networked node."

- **Reads:** `--query_records`, `--target_records` (standard paths); `cache_dirs.json`.
- **Writes:** structures into each record's `struct_path`; the FSDB into its `struct_path`; PISA
  under `pocket_dir/pisa/`; `--query_records_path`, `--target_records_path` (default: the input
  path); `failed_entries.json`.
- **Options:** `--pisa_source`, `--threads`, temp group.
- **Validates:** manifest present; probes foldseek only when a `foldseek_db` record must be
  downloaded. A side left empty after fetching is fatal, as today.
- **Reuses:** `StructureDownloader.download_missing_structures`, `fetch_missing_fsdb` (as a step
  function), `pockets.pisa.download_pisa_interfaces`.
- **Notes:**
  - Structures first, then PISA only for the pisa records that survived. Today PISA is fetched
    only for records that also passed preprocessing; `fetch` runs before preprocessing, so a pisa
    record whose chain is missing also has its PISA data fetched. Output is unchanged.
  - Later commands are offline afterwards, with two exceptions. FSDB-PDB hits are unknown until
    `align`, which fetches their PISA data and structures. And PISA caches only successes, so a
    standalone `pockets` retries any entry `fetch` failed.
  - A pisa record whose PISA data cannot be fetched is kept; `pockets` reports it as
    `pocket_not_built`. Dropping it here would remove its chain from `alignment.tsv`.
  - `fetch` takes no `--aligner`, so it downloads an FSDB that `align` will then reject with
    `--aligner seq`. Only `search` checks before the download.
  - Rerunning `fetch` on its own output cannot retry a dropped structure: the record is gone. Rerun
    `parse`, or keep the input by passing `--query_records_path` / `--target_records_path`.

### `align`

Align query chains against target chains. Owns preprocessing, its scratch, and FSDB-PDB expansion.

- **Reads:** `--query_records`, `--target_records`; `cache_dirs.json`.
- **Writes:** `--alignment_path`; `--query_records_path`, `--target_records_path` (default: the
  input path); `failed_entries.json`.
- **Options:** `-a/--aligner` (`foldseek`), `--threads`, `--pisa_source` (for expansion), temp
  group.
- **Validates:** manifest present; aligner value and foldseek probe; an FSDB target requires
  `--aligner foldseek`. Every non-`foldseek_db` record's `struct_path` exists; a missing one is
  dropped as `structure_not_found`, so a skipped `fetch` fails cleanly rather than letting a gemmi
  error from `preprocess_records` or `SequenceAligner.align_records` escape `cli()`. A side left
  empty is fatal. Expansion whose hits all fail to fetch is fatal, as today.
- **Reuses:** `foldseek_preprocessing`, `foldseek_alignment`, `local_alignment`,
  `expand_fsdb_pdb_targets`, as step functions. Expansion builds its `QTProcessor` from the
  manifest.
- **Notes:**
  - Sets the `foldseek_db` record's `fsdb_pockets`: `"pisa"` and expansion for a PDB-named
    database, `"whole_chain"` (INFO log) for any other.
  - Preprocessing cache: a record is preprocessed unless
    `<foldseek_preprocessed_structure_dir>/<preprocess_name>.cif.gz` exists. This
    replaces `StructurePreprocessor`'s `set_output_directory` → `update_cache` sequence and its
    unenforced call order.
  - Rerunnable on its own output. Given a target file already holding expanded records, `align`
    keeps only the `foldseek_db` record and expands again (Target-side shapes).

### `pockets RECORDS [RECORDS ...]`

Build pockets. "Give me the PISA/vdw pocket of `4Q5J:B_F`."

- **Reads:** one or more records files, named positionally. No default. `cache_dirs.json`.
- **Writes:** `--pockets_path`, overwritten with the union of every input's pockets; each input
  records file, rewritten in place without its `pocket_not_built` records (no option; copy a file
  first to keep it); `pocket_dir/<method>_pockets.json` (cache, as today); `failed_entries.json`.
- **Options:** `--pisa_source` (downloads any PISA data `fetch` did not, and retries its failures).
  `--results_dir` is required, as on every command but `parse` and `search`.
- **Validates:** each named file exists; manifest present.
- **Reuses:** `PocketFetcher.fetch_pockets`, `dump_pockets`; new `load_pockets`.
- **Notes:**
  - `foldseek_db` records are skipped and kept: no builder, not a failure.
  - No `struct_path` check: the builders already map a missing structure to None, which becomes
    `pocket_not_built`.
  - Name both records files in one call: `pockets.json` holds only the files named, and `compare`
    rejects a record with no pocket.
  - Pre-existing, unchanged: one `pocket_id` built by two methods (`4Q5J:B_F` forced vdw on one
    side, pisa on the other) keeps one pocket, as `pockets |= method_pockets` does today. If the
    later method gives None, every record with that `pocket_id` is `pocket_not_built`.

### `compare`

Compare the pockets of every aligned pair.

- **Reads:** `--query_records`, `--target_records`, `--alignment`, `--pockets`.
- **Writes:** `--pocket_comparison_path`; `unknown_ids.json`, `incorrect_mapping.json` under
  `results_dir`.
- **Options:** none beyond the inputs. Without edited or external inputs, a rerun gives the same
  output.
- **Validates:**
  - Inputs exist.
  - Every non-`foldseek_db` record in both files has a non-null pocket in `pockets.json`. A
    missing one is fatal: it means `pockets` was not run on that file, and would otherwise give
    zero rows silently.
  - The offset table is present when the DB ships one (as today).
- **Reuses:** `pocket_comparison.compare_pockets`.
- **Notes:**
  - No `--fsdb_dir`: a small helper in `foldseek.py` finds the offset table by resolved DB path.
    Only `human_domains` ships one, and its path does not depend on `fsdb_dir`.
  - Deletes stale `unknown_ids.json` / `incorrect_mapping.json` on entry. Both are written only
    when non-empty, so a rerun would otherwise keep the old copies.

### `superpose`

Superpose the top targets onto each query. "Rerun with another `--align_count` or
`--align_struct_method`."

- **Reads:** `--query_records`, `--target_records`, `--pocket_comparison`, `--alignment`.
- **Writes:** `--aligned_structure_dir`.
- **Options:** `--align_struct_method` (`auto`), `--align_count`, `--threads`.
- **Validates:**
  - The aligner is inferred from `alignment.tsv`: `u` is `"-"` (the seq aligner's placeholder)
    → `seq`, else `foldseek`. `resolve_align_struct_method(method, inferred)` resolves `auto` and
    rejects `foldseek` with `seq`. An empty `alignment.tsv` has nothing to superpose.
  - For an FSDB target (a `foldseek_db` record is present): rejects `pocket`. Its alignment always
    comes from foldseek.
  - Probes foldseek only for an FSDB rebuild.
- **Reuses:** `StructureAligner.align_structs`, unchanged.

### `search QUERY TARGET`

Unchanged workflow and options, plus:

- `--query_records_path`, `--target_records_path`, `--pockets_path`, `--failed_entries_path`.
  Each is a `Settings` field, so the places rule applies (Docs and tests).
- Writes the records files, `cache_dirs.json`, `pockets.json` and `failed_entries.json` (truncated
  first), so standalone commands can rerun a step in its `results_dir`.
- The only command with `-j/--job_file` and `job_settings.json`.
- Keeps the FSDB checks (`foldseek` aligner required, `pocket` rejected) before the FSDB download,
  so a bad combination fails before the download starts.
- Its pockets step passes `download=False` to the pisa builder (through `builder_options`): its
  fetch step already tried every entry, and PISA caches only successes, so each entry is attempted
  once, as today.

## Option matrix

Columns: P = parse, F = fetch, A = align, Pk = pockets, C = compare, Sp = superpose, S = search.
"in"/"out" mark an option that names a file read or written.

| Option | P | F | A | Pk | C | Sp | S |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `QUERY`, `TARGET` | ✓ | | | | | | ✓ |
| `RECORDS ...` | | | | in, out | | | |
| `-j/--job_file` | | | | | | | ✓ |
| `-v/--verbosity`, `--log_path`, `--results_dir` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `--failed_entries_path` | ✓ | ✓ | ✓ | ✓ | | | ✓ |
| `-q`, `-t` pocket methods | ✓ | | | | | | ✓ |
| `--cache_dir`, `--pdb_dir`, `--alphafold_dir`, `--pocket_dir`, `--foldseek_preprocessed_structure_dir`, `--fsdb_dir` | ✓ | | | | | | ✓ |
| `--pisa_source` | | ✓ | ✓ | ✓ | | | ✓ |
| `--threads` | | ✓ | ✓ | | | ✓ | ✓ |
| `--temp_dir`, `--delete_tmp` | | ✓ | ✓ | | | | ✓ |
| `-a/--aligner` | | | ✓ | | | | ✓ |
| `--align_count`, `--align_struct_method`, `--aligned_structure_dir` | | | | | | ✓ | ✓ |
| `--query_records`, `--target_records` | | in | in | | in | in | |
| `--query_records_path`, `--target_records_path` | out | out | out | | | | out |
| `--alignment` / `--alignment_path` | | | out | | in | in | out |
| `--pockets` / `--pockets_path` | | | | out | in | | out |
| `--pocket_comparison` / `--pocket_comparison_path` | | | | | out | in | out |
| `--job_settings_path` | | | | | | | ✓ |

Every current `Settings` field appears in at least one command besides `search`, except
`job_settings_path` (search only, like `-j/--job_file`, by decision).

## Validator placement

| Check (today's site) | Commands |
| --- | --- |
| Job file, query/target given once (`configure_workflow`) | search |
| `resolve_paths` | search; others resolve their own defaults |
| Cache manifest present (new) | fetch, align, pockets |
| Pocket method names (`process_qt_cmdline_input`) | parse, search |
| Side empty after parsing (`configure_query_target`) | parse, search |
| FSDB entry is the only target entry (new; today only row 0 decides) | parse, search |
| No `foldseek_db` query entry (new) | parse, search |
| Side empty after fetching (`fetch_missing_structures`) | fetch, search |
| Every record's `struct_path` exists (new) | align |
| Aligner value (`resolve_aligner`) | align, search |
| Foldseek probe (`resolve_aligner`) | align, search; fetch only for an FSDB download; superpose only for an FSDB rebuild |
| `resolve_align_struct_method` | superpose (aligner inferred from `alignment.tsv`), search |
| `resolve_threads` | fetch, align, superpose, search |
| `resolve_delete_tmp` | fetch, align, search |
| `resolve_pisa_source` | fetch, align, pockets, search |
| FSDB requires `--aligner foldseek` (`configure_query_target`) | align, search |
| FSDB rejects `pocket` (`configure_query_target`) | superpose, search (before the download) |
| Every record has a pocket (new) | compare |
| `results_dir` required | fetch, align, pockets, compare, superpose |

## Invariants at risk

Each is documented in project-overview, Invariants; this lists how option A touches them.

- `res_auth_ids` order must survive `pockets.json`, and so must any hand-written one: overlap pairing
  is positional and nothing checks the order.
- `preprocess_name` stays the alignment join key. Expanded records keep the Foldseek entry name.
- A local file's `preprocess_name` hashes its contents at `parse` time. Editing the file before
  `align` leaves a stale name.
- `seq_pos` must come from `parse_pocket_from_struct`. A hand-written `pockets.json` that computes
  it any other way gives zero overlap and no error.
- `compare_pockets` must not write to a `Pocket`.
- Step 8 stays in `StructureAligner.align_structs`.
- Per-step `success` filtering is replaced by "records files hold only usable records". Every step
  that filtered (`local_alignment`, `get_pockets`, `preprocess_records`) now takes every record.

## Behaviour changes vs today

- `QTRecord`: `success` and `failure_reason` removed, `fsdb_pockets` added.
- New files in `results_dir`: `query_records.json`, `target_records.json`, `cache_dirs.json`,
  `pockets.json`, `failed_entries.json`. `parse` and `search` truncate `failed_entries.json` on
  entry.
- Record paths are absolute, including local files and user Foldseek databases. A relative path to
  the bundled `human_domains` DB therefore now matches its offset table, so its target residue ids
  become UniProt coordinates; today only the exact bundled path string matches.
- A `foldseek_db` target given beside other target entries is rejected. Today only row 0 decides.
- A `foldseek_db` query entry is rejected.
- `compare` deletes stale `unknown_ids.json` / `incorrect_mapping.json`.
- `fetch` (and so `search`) downloads PISA data before alignment instead of while building pockets.
  Same files, earlier, and also for records that later fail preprocessing; `search` still attempts
  each entry once.
- `pockets/pisa/errors.json` (written only on failure, overwritten): against an FSDB-PDB target its
  last writer is now `align`'s expansion, not the pockets step. When both query and hit entries
  fail, it holds the hit failures; today it holds the query failures.
- **Library breaks:**
  - `pm.query_df` / `pm.target_df` are gone. Deferring step 8 (`search(align_count=0)`, then
    `StructureAligner.align_structs`) reads the records files instead.
  - `QTRecord.success` / `failure_reason` removed.
  - `StructurePreprocessor` ignores `success` and drops its `set_output_directory` →
    `update_cache` sequence.
  - The `resolve_*` methods move off `PocketMapper` into `settings.py`.
  - `pisa_pockets` gains a `download` option (default True).
  - `QTProcessor.parse_individual_qt` returns the rejection reason instead of None.
  - `QTProcessor.process_qt_cmdline_input` returns `(records, rejected)`: a list of `QTRecord`
    dicts, not a DataFrame, and the rejected entries with their reasons.

Everything else must be row-identical.

## Docs and tests (when implemented)

- README: Usage, Options, Outputs (the new files).
- project-overview: Pipeline, Settings, As a library. The five-places rule becomes four places per
  command: `Settings` (search only), the function signature, the parser and README Options. The
  dispatch table replaces `cli()`'s kwarg block.
- `constants.CLI_SEARCH_EPILOG`, plus an epilog per command.
- Check: per command, its function's `inspect.signature` and its subparser's `_actions` dests must
  match.
- e2e:
  - One case per command.
  - A "chained commands == `search`" case (sort rows before diffing).
  - `failed_entries.json` contents for the deliberately wrong fixture lines
    (`invalid_residues.txt` line 2, `forced_pisa_mixed.txt` line 2).
  - `align` rerun on its own output against an FSDB-PDB target gives the same rows.
  - `pockets` with a single records file, then `compare` rejecting the other side's records.
  - A `foldseek_db` target beside a structure target, rejected by `parse`.
  - A `foldseek_db` query, rejected by `parse`.
  - A standalone command run from another cwd, reading the manifest.
  - `superpose` after `align -a seq`, with no aligner restated.
