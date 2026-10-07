# Project overview
Project implementation specifics. Cross-module and derived facts only. Anything a single docstring or comment already states appears here as a pointer to that site, never as a second copy.

## Pipeline

`cli.py` is **the only module that knows about argv or exit codes**. Subcommands: `search` plus one per
step (`parse`, `fetch_structures`, `align`, `pockets`, `compare`, `superpose`). The steps of `search()` are in the
`pocketmapper.py` module docstring. Parser details (`OPTIONS` table, per-command `COMMANDS` layout,
optional query/target positionals for `search` and `parse`, every default None) are documented in
`cli.py`. Every command takes a job file; parsed arguments are layered over it, so an unset option must
parse to None or it would hide the job file's value.

### Steps and hand-off files

- **Two layers.** `steps/<step>.py` holds an entry function named after the step (job file + arguments,
  path defaults under `results_dir`, logging, validation, temp dir) and a core function (explicit
  paths and values, no `Settings`, no state); conventions in the `steps` package docstring.
  `PocketMapper.search` builds `Settings`, opens the log and temp scope once, and calls every entry
  function with `asdict(settings)` as the job file, so chained commands == `search` by construction.
- **Hand-off files** (`records.py` reads/writes them): `query_records.json` / `target_records.json`
  (JSON list of `QTRecord` dicts, JSON to keep None/bools), `cache_dirs.json` (absolute cache dirs;
  written by parse and search, read by fetch_structures, align, pockets, which take no cache options),
  `failed_entries.json` (truncated by parse/search, appended by fetch_structures/align/pockets),
  `alignment.tsv`,
  `pockets.json` (`write_pockets_file`/`read_pockets_file`: `version` 2, `pockets` with `null` for a record
  given and not built, `chains` = `preprocess_name` -> pocket_ids; per-method cache files keep the plain
  `dump_pockets` mapping), `pocket_comparison.tsv` (compare's
  `unknown_ids.json`/`incorrect_mapping.json` go beside it).
- **A records file holds only usable records.** fetch_structures, align and pockets drop the records they fail
  into `failed_entries.json` (reasons in the `records.failed_entry` callers). No step filters on a
  `success` flag; `QTRecord` has none.
- **The target side's shape replaces the old `self.` flags.** A `foldseek_db` record, always the only
  target entry from parse, means an FSDB target. Its pocket mode comes from the alignment, not from
  any record: `lib.fsdb_pocket_mode` gives `"pisa"` when any hit name is PDB-style, else
  `"whole_chain"` (pockets, compare and superpose each call it). No step appends hit records to a
  records file; the hits' pockets and names live only in `pockets.json`.
- **Record paths are absolute**: parse passes absolute cache dirs to `QTProcessor`, which resolves an
  entries file, a local file and a user FSDB against `work_dir`. `pocket_id`/`struct_info` stay as
  typed.
- `--pisa_source` is not in the manifest; only pockets takes it. fetch_structures and align fetch no
  PISA: the pockets step downloads it, so a pisa record whose PISA fails is aligned, then dropped as
  `pocket_not_built`.

### Input grammar

`struct_info[:chain_info[:residue_info]]`, or a file of such lines. Forms: README "Input format".

- A local file like `4Q5J.cif.gz:B_F` resolves to `vdw`, not `pisa` (PISA is PDB-only). This is how the
  mixed-input e2e fixtures reach vdw.
- **Two tables in `QTProcessor.__init__` are the whole grammar**: `pocket_methods` (method → pattern,
  warning phrase) and `struct_type_pocket_methods` (struct type → methods, in try order).
  `determine_pocket_method` takes the first match; `validate_pocket_method` checks every record (forced or
  inferred) against its method's pattern, so forced and inferred cannot drift.
- **vdw's partner chain is required.** When optional, a forced `vdw` on `4Q5J:A` reached gemmi as chain `None`.
  Tightening is inference-neutral (`whole_chain` claims bare chains first) — verified by replaying the old
  per-struct_type ladder against the loop over 400 generated entries.
- **A rejected entry is skipped, not fatal**: `parse_individual_qt` returns `(None, reason)` and parse
  lists it in `failed_entries.json` as `invalid_entry`. `steps.parse.parse_entries` (both sides, on the
  one-side `parse_side`, which takes an entry, a file or a list) raises only when a
  side ends up empty or a `foldseek_db` target is not the only target entry; a `foldseek_db` query entry
  is itself an `invalid_entry`. An unrecognised method *name* raises up front in
  `process_qt_cmdline_input`.

### Pocket shape

`pockets.pocket.Pocket` declares the fields; `parse_pocket_from_struct` derives `seq_pos`/`whole_chain`.
`residues` is keyed by author seqid string. `res_auth_ids` is not `list(residues)`: it is the ordered list
the comparison walks (PISA seeds it from the interface; `residues` is chain order).

**Every producer must emit `res_auth_ids` ascending; nothing checks.** `overlap_ids` + `superpose` pair the
two sides position by position, so misordering gives wrong `rmsd`, `ca_dists` and transforms while
`overlap_count` and identities stay correct — no warning. vdw/whole_chain walk the chain; pisa/passthrough sort.

### Open searches

A query and target on one chain share a `preprocess_name`, so some rows carry a query-only `pocket_id` in
`target`. `StructureAligner.align_structs` filters them before target lookup (else a bare pandas
`KeyError`). Output shape: README "Open searches".

### Foldseek-DB targets

With a `foldseek_db` target record, no target structures are fetched; `foldseek.extract_fsdb_structures`
rebuilds the superpose selection via `createsubdb` + `convert2pdb`.

- **`createsubdb` rejects `--threads`** (exits non-zero); it is the only one of the five subcommands that
  does. Hence `run_foldseek` adds no flags; each caller builds its own list.
- `align_structs` interprets target ids once per call: given target records → `pocket_id`s mapped through
  `preprocess_name` (PDB DB; `steps.superpose` builds them from `pockets.json` `chains`, hit names in
  alignment order); none → entry names (other DBs). So `fsdb_structures/` never mixes the two.
- Target pocket: `steps.pockets.expand_fsdb_pdb_targets` (PDB DB; needs `--alignment`, an FSDB target
  without one is an error) or `synthesise_target_pocket` (other). **Every hit of a PDB DB is a
  `chains` key**: `[]` for a non-PDB name, no interface or only unspellable pairs; a hit whose
  structure download failed lists its pocket_ids with `null` pockets.
- **One `pocket_id` can sit behind two `preprocess_name`s** (`4q5j-assembly1_B`, `-assembly2_B` →
  `4Q5J:B_F`). `existing_calcs` scores only the first; the transform is whichever assembly Foldseek listed first.
- **Pockets come from the AU, Foldseek's `tseq` from the assembly.** Verified to agree normally (4Q5J
  self-comparison: `overlap_count == pocket_len`, identity 1.0, RMSD ~1e-14); a populated
  `incorrect_mapping.json` signals divergent numbering.
- `--align_struct_method pocket` is rejected for any FSDB target (`settings.check_fsdb_align_struct_method`,
  in search before the download and in superpose); `--aligner seq` likewise (`check_fsdb_aligner`, search
  and align).

**UniProt renumbering** (`offset_table.tsv`; resolved in `steps.compare`, applied in
`synthesise_target_pocket`):
- Only bundled `human_domains`, looked up by resolved DB path in `foldseek.bundled_offset_table`. User
  DBs keep 0-indexed positions, logged at INFO.
- Only `target_overlap_ids` changes — verified against a same-environment baseline.
- Missing entry or short spec aborts the run (no per-row fallback: it would mix coordinate systems in one
  column). **Refresh the table whenever `BUNDLED_HUMAN_DOMAINS_DB` moves.**

**Bundled DB ships without `.source`**; strip it from any refresh. It duplicates `.lookup` and nothing
reads it (verified: `easy-search`, `createsubdb`, `convert2pdb` all work without it); saves 1.4 MB.
`.lookup` must stay (`extract_fsdb_structures` reads it).

The DB is otherwise at its floor: `_ca` is 53 of 77 MB, 8.4M residues at 6.33 B each
(`--coord-store-mode 2`, smallest mode). zstd -19 saves only 30% and foldseek cannot read a compressed DB.

**No cap on enriched hits**, by choice. `4Q5J:B_F` vs bundled `pdb`: ~4,970 hits / ~3,620 entries, hours
on first run with `pisa_source` `api` (per-assembly PISA behind a sleep); `ftp` is concurrent but untimed
at this scale; reruns hit the cache. Add a cap here if needed.

## Pockets

`PocketFetcher.fetch_pockets` is the entry point; `steps.pockets.build_pockets` calls it on every
non-`foldseek_db` record and drops, as `pocket_not_built`, any record whose `pocket_id` is missing or
None in the merged dict (checked after the merge: a later method can overwrite an earlier one with None).

- **`POCKET_BUILDERS` in `pocket_fetcher` is the whole method table.** Builders live beside their
  primitive and all take `(records, pocket_dir)`, plus keyword options the caller passes per method
  through `fetch_pockets(builder_options=...)` (pisa's `pisa_source` today). A new method =
  one row + one builder.
- Records are dicts, not DataFrames. The fetcher builds every record it is given.
- **`expand_fsdb_pdb_targets` lives in `steps/pockets.py`, not in `pockets/`** (it builds records and
  downloads). It shares `pocket_dir/pisa/` with `pisa_pockets`; both go through
  `pockets.pisa.pisa_cache_layout`, the one place the cache layout is spelled out.

## Downloads

See `pocketmapper/downloads/CLAUDE.md` before changing fetching, pacing, retries or caching.

## Invariants

Breaking these gives silently wrong output. Each is documented at its code site; this is the map plus
checks that live nowhere else.

- **`seq_pos`** — declared on `PocketResidue`, set in `parse_pocket_from_struct`, used in
  `map_pocket_into_alignment`. Computed any other way → zero overlap, no error. Check: pocket vs itself
  gives `overlap_count == pocket_len`. It is not the reported id; `synthesise_target_pocket` keeps
  `seq_pos` 0-indexed while keying by UniProt position, which is what makes renumbering safe.
- **Residue letters mirror Foldseek, not gemmi.** `constants.FOLDSEEK_AA_CODES` is Foldseek's
  `threeToOneAA` verbatim (139 names, rest `X`). gemmi 0.7.5's `one_letter_code` disagrees on 14 (`SEC`→`U`;
  `BAL`, `KYN`, `HZP`→`X`).
- **Foldseek letters compare case-insensitively.** `createdb` lowercases residues with CA B-factor below
  threshold (default 0); 4ER4, 2ER6, 2ER9 hit this wholesale. `compare_pockets` uppercases the four
  sequence columns; `alignment.tsv` keeps Foldseek's casing.
- **Every chain walk reads `first_conformer()`** (`parse_pocket_from_struct`, `vdw_pockets`,
  `SequenceAligner`). Microheterogeneity (4Z0Y:A 252: `HS8`+`HIS`) otherwise shifts later `seq_pos` by
  one. `test_core_10` covers it.
- **A pocket failing `MIN_SEQ_IDENTITY` yields no rows** (goes to `incorrect_mapping.json`).
- **`preprocess_name` is the alignment join key** (`QTProcessor.parse_individual_qt`);
  `pockets.json` `chains` (built by `records.preproc_to_ids` plus the FSDB hits) bridges it to
  `pocket_id`. compare's `check_coverage` errors on an alignment name `chains` lacks (target names
  exempt when synthesising); a `null` pocket or `[]` is a skip. A local file is hashed at parse time; editing it
  before align leaves a stale name. **Local files are hashed with their contents**: by
  basename, `/dir1/structure.pdb:A` and `/dir2/structure.pdb:A` collide and one is silently mis-aligned
  (`test_core_9`, `test_local_7`). Hence identical copies share a name, in-place edits get a new one, and
  local `4Q5J.cif.gz:B` ≠ PDB `4Q5J:B`.
- **`chain_info` is split only by `lib.split_chain_info`.** Indexing `chain_info[0]` breaks on
  multi-char chains; patterns reject `4Q5J:AA_BB`, but library callers bypass them.
- **Passthrough `res_auth_ids` must all be keys of `residues`** — `passthrough_pockets` skips the entry
  otherwise (else a `KeyError` aborts the run). Syntax and repeats are caught in
  `QTProcessor.parse_residue_info`; a repeat would silently misalign overlap lists (wrong `pocket_len`,
  `jaccard_index`, `rmsd`, `ca_dists`).
- **Declared schemas**: `constants.ALIGNMENT_COLUMNS`, `POCKET_COMPARISON_COLUMNS`. New columns go in the
  constant.
- **`compare_pockets` must not write to a `Pocket`.**
- **Aligned structures are named `lib.safe_filename(query_id)`**, not by `pocket_id`; grep the
  `MOLECULE` records instead.
- **Two transform sources by `align_struct_method`** (see `StructureAligner` docstring). Only
  `parse_pocket_transform` may read the pocket transform; never pass a raw `target_to_query_*` cell to gemmi.
- **Superposition lives in `StructureAligner.align_structs`**; `steps.superpose` only derives the target
  shape.
- **`pockets.json` round-trips `res_auth_ids` order and `seq_pos`**; a hand-written one that gets
  either wrong gives silently wrong overlap.

**Changing the comparison without changing behaviour**: capture `compare_pockets`' arguments from a real
run and diff old vs new output. Nothing else covers it.

## Logging and errors

**The package never touches the root logger.** Modules use `logging.getLogger(__name__)` under
`constants.PACKAGE_LOGGER`, which has a `NullHandler`. Never call `logging.info(...)` etc. directly.
Handlers are added in two places, each removed in a `finally`:
- **`cli()` adds stdout** before dispatching, because `layer_settings` can `critical` on a bad job
  file before the log is open (inherited root WARNING lets it through).
- **`search()` and each step add `info.log`** and set the level from `verbosity`, through the
  `lib.log_to_file` context manager, which undoes both on exit, so a second call doesn't write into the
  first log. The handler appends, so a chain of commands writes one log. Nested on the same file (a
  step inside `search`), it adds no second handler, else every line is written twice. Handlers have no
  level of their own.

A library caller gets no console output unless it prints `pocketmapper.*`; `lib.format_handler` applies
the CLI format plus `lib.StageFilter`, which fills `%(stage)s` from `funcName` when absent. Declared
stages are Title Case step names; a function name in the stage column means none was declared. Omitting
`extra` is fine where the function name is the best label.

`extra` is always named **`log_extra`**. Single-stage classes (`StructureDownloader`, `PisaDownloader`,
`StructureAligner`, `StructurePreprocessor`) set `self.log_extra` once and never mutate it; others build
it per function. `QTProcessor`'s `.update()` in `process_qt_cmdline_input` is the one deliberate
exception (call-scoped side name). `PocketMapper` holds no state at all.

Errors: `logger.critical(...)` then `raise PocketMapperError(...)`; `cli()` catches and exits 1. Foldseek
goes through `run_foldseek` so failures follow this too. **No `exit()`/`sys.exit()` outside `cli()`.**

## Settings

`Settings` (`settings.py`) is built only by search. **Every `Settings` field is reachable from the CLI**;
the job file sets nothing the CLI can't. Layering is `settings.layer_settings` for every command
(arguments not None over job file over `SETTING_DEFAULTS`, then `resolve_paths`); search's resolution is
in `resolve_settings` (`resolve_*`). `Settings` has no defaults and is built once. Pocket methods keep
`"auto"` (inferred per entry); `align_struct_method` is resolved per run. The steps take a job dict (or
file) keyed by job keys, never a `Settings`, and call the same `settings.resolve_*` validators.

**Job keys are separate from `Settings`.** `settings.JOB_KEYS` = `fields(Settings)` + `STEP_ONLY_KEYS`
(options only some steps take) is what a job file may hold and what `layer_settings` seeds. Unknown
keys are rejected; one in `REMOVED_JOB_KEYS` gets its replacement named. Search drops keys `Settings`
lacks before building one.

**`work_dir`** (default cwd) is stored absolute by `resolve_paths`, which joins every path setting onto
it; step input arguments (`--alignment`, ...) resolve against it through `settings.input_path`. So
every resolved path is absolute, and a step run from a dump elsewhere resolves as the first did.

**Settings dumps.** `settings.dump_settings` writes every command's values, minus `job_settings_path`,
to `job_settings_path`: search's resolved `Settings` to `job_settings.json`, a step's layered values
(paths resolved, before its validators) to `<results_dir>/<command>_settings.json` (the default comes
from `resolve_paths(values, command)`, not `RESULTS_PATH_DEFAULTS`). `layer_settings` drops
`job_settings_path` from a job file, so a step chained off another's dump never overwrites it. A step
nested in search writes none: `lib.run_scope` records the outermost command and nested scopes yield
False (`lib.outer_command()` reads it).

**A reused `job_settings.json` names every path**, so an argument moves only the path it names.
`results_dir` is not inert beside it all the same: it locates `cache_dirs.json` (fetch_structures, align,
pockets) and is a `temp_dir` emptying root. In fetch_structures and align, explicit
`query_records_path`/`target_records_path` stay out of `layer_settings` and default to their input;
layered in, an explicit output would also become the input's default (breaks `test_steps_9`/`_13`).

A new option goes in **four hand-maintained places per command**: the function signature and its
`layer_settings` arguments dict (for search also `Settings`; for a step-only option
`settings.STEP_ONLY_KEYS` instead), `OPTIONS` plus the command's `COMMANDS` row
in `cli.py`, and README Options (for a step, its row in the step table). A fifth for an option with a
static default: `settings.SETTING_DEFAULTS` (from `constants.DEFAULT_*`, which the help also quotes).
Signature and parser defaults are all None. Path defaults are `settings.CACHE_PATH_DEFAULTS` /
`RESULTS_PATH_DEFAULTS` (`cache_path`, `results_path`). `cli()` passes parsed args straight through its
dispatch table, so there is no kwarg block. Miss one → argparse fails or the option is silently
ignored. Check, per command: `inspect.signature(steps.<x>.<x>)` == subparser `_actions` dests; for
search also `dataclasses.fields(Settings)` == signature == `arguments` keys (modulo `job_file`);
`SETTING_DEFAULTS` keys ⊆ `JOB_KEYS`.

Options are grouped by lifetime in `cli.COMMANDS` and README alike: `in`, `aligned structure`, `cache`,
`out`, `temp`, `advanced`. A file a command reads has a bare name (`--alignment`); one it writes ends in
`_path`. A new path setting picks a group in both.

**Every directory is made by its first writer**; there is no bulk creation step. Only the log's
directory is made up front, by `lib.log_to_file` (the file handler opens at once). Pipeline-side writers use `lib.make_dir`
(critical + `PocketMapperError`); components use bare `os.makedirs`. A new writer into a configurable
path must make its directory, or a path option pointed elsewhere fails.

**`temp_dir` is one option, three dirs** (`query_structures/`, `target_structures/`, `foldseek_tmp/`),
named in `steps/align.py` (all three) and `steps/fetch_structures.py` (`foldseek_tmp/`, for a DB download). The two
structure dirs are made in `foldseek_preprocessing`, `foldseek_tmp/` before each Foldseek call that uses
it (MMseqs2 makes only one level of tmp dir), so a `seq` run makes none. Ownership is
`lib.temp_dir_scope` (search, fetch_structures, align): emptied on entry so reruns can't feed `createdb`
stale structures, deleted on a clean exit unless `delete_tmp` is 0. Scopes nest: search's owns temp
for the whole run and the steps' are no-ops inside it, so `search --delete_tmp 0` keeps every step's
scratch. The held dir is released in a `finally`, else a failed library run leaves the next one in that
process neither emptying nor deleting it. Emptying is `lib.is_within`-guarded against `cache_dir` and
`results_dir`, creation is not. It runs after the log is open so the skip warning reaches `info.log`.
`lib.delete_temp_dir` does nothing if it was never made.

`search --help` is generated from `help=` strings; `CLI_SEARCH_EPILOG` holds only examples. Resolution
order, the `aligner` check and tri-state `align_struct_method`: `Settings` docstring and `# 4b.`/`# 4c.` in
`resolve_settings`. `aligner` has no auto mode; a broken `foldseek` fails at 4b before fetching.

## Python versions

**3.10 – 3.14**, verified by full e2e at both ends. Four places must agree: `requires-python`, the
classifiers, `[tool.black] target-version`, README Installation; CI `compat` matrix is a fifth.

- **Floor is 3.10**: `match` statements, PEP 604 dataclass annotations (no `from __future__ import
  annotations`), and biopython itself.
- **`compat` guards the floor, not `lint`** (flake8 parses with its own interpreter). `compileall` at 3.10
  catches syntax; the import step catches annotation/`importlib.resources` failures. It uses
  `pkgutil.walk_packages` with a prefix (`iter_modules` skips subpackage modules) and asserts a module count.
- **Only 3.10 resolves pandas 2.x**; 3.11+ get 3.x. Hence e2e runs 3.10 and 3.14. Both give identical
  row counts on every non-`huge` case.

## Repo layout

- No unit tests; `tests/e2e/` is the whole suite (`pocketmapper-e2e` skill).
- **Deliberately wrong fixture lines — don't fix them:** `invalid_residues.txt` line 2 (`4Q5J:A:9999`,
  `test_invalid_1`, asserted `pocket_not_built`); `forced_pisa_mixed.txt` line 2 (`4Q5J:A`, no partner;
  `test_invalid_8`, asserted `invalid_entry`); `fsdb_mixed_target.txt` (a DB beside a structure;
  `test_steps_11`).
- **`fixtures/job_file.json` must never set a path** (the runner's `--results_dir` would not move it, so
  the case would write outside its dir). Its `align_count: 3` vs `test_settings_2`'s `5` shows priority
  only in `job_settings.json` (5 wins) — check by hand.
- `test_settings_5`/`_6` share `job_file_qt.json`: file-only query/target vs repeated positionally (rejected).
- `test_settings_1` only catches rejected/crashing flags; the runner asserts only on
  `pocket_comparison.tsv`, so silently dropped path options need the per-command signature check.
- **Delete `build/` before building or testing a wheel.** Stale, gitignored, but setuptools reuses
  `build/lib/`, so `pip install .` ships dead modules (`align.py`, `local_aligner.py`, a 3.12-only
  `pisa.py`, top-level `structure_fetcher.py`/`pisa_downloader.py` shadowing `downloads/`). Never edit it.
- Structure parsing is gemmi; biopython only for pairwise alignment and SVD superposition.
- `StructurePreprocessor` caches on `<foldseek_preprocessed_structure_dir>/<preprocess_name>.cif.gz`
  existing and writes its own `.part`;
  `StructureDownloader` uses each record's `struct_path` and `lib_download`'s `.part`.
- **`StructureDownloader` makes a destination's parent directory** just before downloading into it;
  one it cannot make surfaces as `structure_not_found`.
- **e2e cases can chain commands** (`parse ... ; fetch_structures ; cd DIR ; ...`) and assert `files=`, `failed=`,
  `same=`; format in the `run_e2e.sh` header. `test_steps_9` needs `POCKETMAPPER_PDB_FSDB`: any
  Foldseek DB of `<pdb>-assembly<N>.cif.gz` files (`createdb`) has PDB-style entry names, so a few
  cached mmCIFs make a small local one.

## As a library

`PocketMapper().search(...)`, a step's entry or core function in `steps`, or any component directly
(`qt_processor`, `downloads.*`, `structure_preprocessor`, `pockets.pocket_fetcher` or one builder,
`sequence_aligner`, `structure_aligner`, `foldseek`). No component or step takes a `Settings`; the
entry functions take its dict as `job_file`.

- **Superposition can be deferred**: `search(align_count=0)`, then `steps.superpose.superpose(results_dir=...)`,
  or `StructureAligner.align_structs` with `records.read_records` of the records paths in the returned
  settings (FSDB targets also need `fsdb_path`; `steps.superpose` derives it). Records point at
  `pdb_dir`/`alphafold_dir`, which temp deletion never touches. `query_ids`/`target_ids`/`overwrite`
  exist for batching this.
- `__init__` exports only `PocketMapper` and `__version__`; always import submodules explicitly.
- `search(job_file=...)` needs no query/target if the file sets them; both ways is rejected.
- **`search()` side effects**: logger level + `info.log` handler for the call; empties `temp_dir` on
  entry and `rmtree`s it at the end unless `delete_tmp=0`. Both guarded by `lib.is_within` (under
  `cache_dir` or `results_dir`) — a safety net, not a licence. `fetch_structures` and `align` do the same with
  their own `temp_dir` when called outside `search`.
- Returns only resolved `Settings` as a dict; results are in `pocket_comparison.tsv` / `alignment.tsv`.
