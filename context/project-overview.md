# Project overview
Project implementation specifics. Cross-module and derived facts only. Anything a single docstring or comment already states appears here as a pointer to that site, never as a second copy.

## Pipeline

`cli.py` is **the only module that knows about argv or exit codes**; `search` is its one subcommand. The
eight steps of `search()` are in the `pocketmapper.py` module docstring. Parser details (optional
query/target positionals, defaults from `constants`) are documented at the parser. The job file is layered
on top of parsed arguments, so nothing distinguishes a default from a typed value.

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
- **A rejected entry is skipped, not fatal**; `configure_query_target` raises only when a side ends up
  empty. An unrecognised method *name* raises up front in `process_qt_cmdline_input`.

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

With `self.fsdb_target` set, no target structures are fetched; `foldseek.extract_fsdb_structures` rebuilds
the step-8 selection via `createsubdb` + `convert2pdb`.

- **`createsubdb` rejects `--threads`** (exits non-zero); it is the only one of the five subcommands that
  does. Hence `run_foldseek` adds no flags; each caller builds its own list.
- `align_structs` interprets target ids once per call: given target records → `pocket_id`s mapped through
  `preprocess_name` (PDB DB); none → entry names (other DBs). So `fsdb_structures/` never mixes the two.
- Target pocket: `expand_fsdb_pdb_targets` (PDB DB) or `synthesise_target_pocket` (other). **PDB hits
  with no usable PISA data are dropped.**
- **One `pocket_id` can sit behind two `preprocess_name`s** (`4q5j-assembly1_B`, `-assembly2_B` →
  `4Q5J:B_F`). `existing_calcs` scores only the first; the transform is whichever assembly Foldseek listed first.
- **Pockets come from the AU, Foldseek's `tseq` from the assembly.** Verified to agree normally (4Q5J
  self-comparison: `overlap_count == pocket_len`, identity 1.0, RMSD ~1e-14); a populated
  `incorrect_mapping.json` signals divergent numbering.
- `--align_struct_method pocket` is rejected for any FSDB target in `configure_query_target`.

**UniProt renumbering** (`offset_table.tsv`; resolved in `compare_pockets_based_on_alignment`, applied in
`synthesise_target_pocket`):
- Only bundled `human_domains` (the one `bundled_foldseek_dbs` entry with `offset_path`), looked up by
  resolved DB path. User DBs keep 0-indexed positions, logged at INFO.
- Only `target_overlap_ids` changes — verified against a same-environment baseline.
- Missing entry or short spec aborts the run (no per-row fallback: it would mix coordinate systems in one
  column). **Refresh the table whenever `BUNDLED_HUMAN_DOMAINS_DB` moves.**

**Bundled DB ships without `.source`**; strip it from any refresh. It duplicates `.lookup` and nothing
reads it (verified: `easy-search`, `createsubdb`, `convert2pdb` all work without it); saves 1.7 MB.
`.lookup` must stay (`extract_fsdb_structures` reads it).

The DB is otherwise at its floor: `_ca` is 70 of 98 MB, 11.2M residues at 6.33 B each
(`--coord-store-mode 2`, smallest mode). zstd -19 saves only 16% and foldseek cannot read a compressed DB.

**No cap on enriched hits**, by choice. `4Q5J:B_F` vs bundled `pdb`: ~4,970 hits / ~3,620 entries, hours
on first run with `pisa_source` `api` (per-assembly PISA behind a sleep); `ftp` is concurrent but untimed
at this scale; reruns hit the cache. Add a cap here if needed.

## Pockets

`PocketFetcher.fetch_pockets` is the entry point; `PocketMapper.get_pockets` filters to `success` rows.

- **`POCKET_BUILDERS` in `pocket_fetcher` is the whole method table.** Builders live beside their
  primitive and all take `(records, pocket_dir)`, plus keyword options the caller passes per method
  through `fetch_pockets(builder_options=...)` (only `pisa_source` today). A new method = one row + one builder.
- Records are dicts, not DataFrames. The fetcher ignores `success`; filtering is the caller's job.
- **`expand_fsdb_pdb_targets` lives on `PocketMapper`, not in `pockets/`** (it builds records and
  downloads). It shares `pocket_dir/pisa/` with `pisa_pockets`; both go through
  `pockets.pisa.download_pisa_interfaces`, the one place the cache layout is spelled out.

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
  `preproc_to_ids` bridges it to `pocket_id`. **Local files are hashed with their contents**: by
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
- **Step 8 lives in `StructureAligner.align_structs`**; `PocketMapper.align_structs` only unpacks settings.

**Changing step 7 without changing behaviour**: capture `compare_pockets`' arguments from a real run and
diff old vs new output. Nothing else covers it.

## Logging and errors

**The package never touches the root logger.** Modules use `logging.getLogger(__name__)` under
`constants.PACKAGE_LOGGER`, which has a `NullHandler`. Never call `logging.info(...)` etc. directly.
Handlers are added in two places, each removed in a `finally`:
- **`cli()` adds stdout** before `search()`, because `configure_workflow` can `critical` on a bad job file
  before `configure_logging` (inherited root WARNING lets it through).
- **`search()` adds `info.log`** and sets the level from `verbosity`; `reset_logging` undoes both, so a
  second `search()` doesn't write into the first log. Handlers have no level of their own.

A library caller gets no console output unless it prints `pocketmapper.*`; `lib.format_handler` applies
the CLI format plus `lib.StageFilter`, which fills `%(stage)s` from `funcName` when absent. Declared
stages are Title Case step names; a function name in the stage column means none was declared. Omitting
`extra` is fine where the function name is the best label.

`extra` is always named **`log_extra`**. Single-stage classes (`StructureDownloader`, `PisaDownloader`,
`StructureAligner`, `StructurePreprocessor`) set `self.log_extra` once and never mutate it; others build
it per function. `QTProcessor`'s `.update()` in `process_qt_cmdline_input` is the one deliberate
exception (call-scoped side name). `PocketMapper` holds no stage state.

Errors: `logger.critical(...)` then `raise PocketMapperError(...)`; `cli()` catches and exits 1. Foldseek
goes through `run_foldseek` so failures follow this too. **No `exit()`/`sys.exit()` outside `cli()`.**

## Settings

**Every `Settings` field is reachable from the CLI**; the job file sets nothing the CLI can't. Layering
is in `configure_workflow` (job file over args, then `resolve_*`). `Settings` has no defaults and is built
once. Pocket methods keep `"auto"` (inferred per entry); `align_struct_method` is resolved per run. A reused
`job_settings.json` names every field, so it pins everything.

A new option goes in **five hand-maintained places**: `Settings`, the `search()` signature + its
`arguments` dict, the parser in `cli.py`, `cli()`'s kwarg block (easiest to forget), and README Options.
Static defaults are `constants.DEFAULT_*` shared by parser and signature; run-dependent defaults are `None`
and resolved in `configure_workflow`. Miss one → option silently ignored. Check: `dataclasses.fields(Settings)`,
`inspect.signature(PocketMapper.search)`, `arguments` keys, subparser `_actions` dests and `cli()`'s
`x=args.x` lines must match (modulo `job_file`).

Options are grouped by lifetime in `build_parser` and README alike: `aligned structure`, `cache`,
`out`, `temp`. A new path setting picks a group in both.

**`temp_dir` is one option, three dirs** (`query_structures/`, `target_structures/`, `foldseek_tmp/`),
computed in `configure_temp_dir` and held on `PocketMapper`, not `Settings`. Only the first two are
created (Foldseek makes its own). `temp_dir` is emptied on entry so reruns can't feed `createdb` stale
structures; emptying is `lib.is_within`-guarded, creation is not. It runs after `configure_logging` so the
skip warning reaches `info.log`.

`results_dir` must stay in `dirs_to_create` in its own right, or `--aligned_structure_dir` elsewhere
leaves `info.log` with no directory.

`search --help` is generated from `help=` strings; `CLI_SEARCH_EPILOG` holds only examples. Resolution
order, the `aligner` check and tri-state `align_struct_method`: `Settings` docstring and `# 4b.`/`# 4c.` in
`configure_workflow`. `aligner` has no auto mode; a broken `foldseek` fails at 4b before fetching.

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
  `test_invalid_1`); `forced_pisa_mixed.txt` line 2 (`4Q5J:A`, no partner; `test_invalid_8`).
- **`fixtures/job_file.json` must never set a path** (it would beat the runner's `--results_dir`). Its
  `align_count: 3` vs `test_settings_2`'s `5` shows priority only in `job_settings.json` — check by hand.
- `test_settings_5`/`_6` share `job_file_qt.json`: file-only query/target vs repeated positionally (rejected).
- `test_settings_1` only catches rejected/crashing flags; the runner asserts only on
  `pocket_comparison.tsv`, so silently dropped path options need the five-way check.
- **Delete `build/` before building or testing a wheel.** Stale, gitignored, but setuptools reuses
  `build/lib/`, so `pip install .` ships dead modules (`align.py`, `local_aligner.py`, a 3.12-only
  `pisa.py`, top-level `structure_fetcher.py`/`pisa_downloader.py` shadowing `downloads/`). Never edit it.
- Structure parsing is gemmi; biopython only for pairwise alignment and SVD superposition.
- `StructurePreprocessor` requires `set_output_directory()` → `update_cache()` → `preprocess_records()`
  (unenforced) and writes its own `.part`; `StructureDownloader` uses each record's `struct_path` and
  `lib_download`'s `.part`.
- **Nothing creates a structure's parent directory**; a missing one surfaces as `structure_not_found`.
  `configure_workflow` makes `pdb_dir`/`alphafold_dir`; library callers must.

## As a library

`PocketMapper().search(...)` or any component directly (`qt_processor`, `downloads.*`,
`structure_preprocessor`, `pockets.pocket_fetcher` or one builder, `sequence_aligner`, `structure_aligner`,
`foldseek`). No component takes a `Settings`; `pocketmapper.py` unpacks it per call site.

- **Step 8 can be deferred**: `search(align_count=0)`, then `StructureAligner.align_structs` with
  `pm.query_df`/`pm.target_df` records and the result paths from `pm.settings`. Verified byte-identical
  PDBs on both structure and FSDB paths (records point at `pdb_dir`/`alphafold_dir`, which `delete_tmp`
  never touches). `query_ids`/`target_ids`/`overwrite` exist for batching this.
- `__init__` exports only `PocketMapper` and `__version__`; always import submodules explicitly.
- `search(job_file=...)` needs no query/target if the file sets them; both ways is rejected.
- **`search()` side effects**: logger level + `info.log` handler for the call; empties `temp_dir` on
  entry and `rmtree`s it at the end unless `delete_tmp=0`. Both guarded by `lib.is_within` (under
  `cache_dir` or `results_dir`) — a safety net, not a licence.
- Returns only resolved `Settings` as a dict; results are in `pocket_comparison.tsv` / `alignment.tsv`.
