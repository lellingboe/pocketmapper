# Project overview
Cross-module and derived facts only. Anything one docstring or comment already states is left to that site
and not copied here.

## Project Status

- This is an early development project, features may change or break between patches.
- Names/location of outputs may change, but typically the underlying data should remain fairly consistent.
- The project is currently going through refactoring and reorganization, some things may be inconsistent
  while this process is happening

## Where to start

- `pocketmapper.py` module docstring: the steps of `search()` and what each writes under `results_dir`.
- `steps/__init__.py`: each step's entry/core layers, settings layering, and how every step re-derives its
  records from the entries (no records are handed on).
- `cli.py`: the only module that knows about argv or exit codes; why every option defaults to None.
- `entries.py`: `failed_entries.json`, and why each failure warns once per run.
- `qt_processor.py`: the input grammar; the two tables in `QTProcessor.__init__` are the whole of it.
- `pockets/CLAUDE.md`, `downloads/CLAUDE.md`: read before changing those packages.

## Invariants

Breaking these gives silently wrong output. Each is explained at its code site; this is the map, plus
checks that live nowhere else.

- **`seq_pos`**: declared on `PocketResidue`, set in `parse_pocket_from_struct`, read in
  `map_pocket_into_alignment`. Computed any other way, it gives zero overlap and no error. Check: a pocket
  compared with itself gives `overlap_count == pocket_len`.
- **`res_auth_ids` ascending** in every producer (`Pocket`); nothing checks it. `pockets.json` round-trips
  this order and `seq_pos` (`read_pockets_file`).
- **Residue letters mirror Foldseek, not gemmi** (`constants.FOLDSEEK_AA_CODES`). They compare
  case-insensitively (`compare_pockets`), because Foldseek lowercases some residues; 4ER4, 2ER6 and 2ER9
  are lowercased wholesale.
- **Every chain walk reads `first_conformer()`** (`parse_pocket_from_struct`, `vdw_pockets`,
  `SequenceAligner`), and a new walk must too. `test_core_10` covers it.
- **`preprocess_name` is the alignment join key** (`QTProcessor.build_record`), bridged to `pocket_id` by
  `pockets.json` `chains`. Every step recomputes it (re-hashing local files), so consistency depends on the
  inputs not changing between steps. `steps.compare.check_coverage` catches a mismatch. `test_core_9` and
  `test_local_7` cover local-file names.
- **`chain_info` is split only by `lib.residues.split_chain_info`.** `chain_info[0]` breaks on multi-character
  chains; the patterns reject them, but library callers bypass the patterns.

**Changing the comparison without changing behaviour**: capture `compare_pockets`' arguments from a real
run and diff old vs new output. Nothing else covers it.

## Logging and errors

- The package never touches the root logger: `__init__` gives `constants.PACKAGE_LOGGER` a `NullHandler`.
- Handlers are added only in `cli()` (stdout) and `lib.logs.log_to_file` (`info.log`, for `search` and every
  step), each removed in a `finally`.
- Errors go through `PocketMapperError`, which `cli()` turns into exit 1. Foldseek runs through
  `run_foldseek` so its failures raise it too.

## Python versions

**3.10 – 3.14**, verified by full e2e at both ends. Four places must agree: `requires-python`, the
classifiers, `[tool.ruff] target-version`, README Installation; CI `compat` matrix is a fifth.

- **Floor is 3.10**: `match` statements, PEP 604 dataclass annotations (no `from __future__ import
  annotations`), and biopython itself.
- **`compat` guards the floor, not `lint`** (ruff flags some post-3.10 syntax from `target-version`, but is
  not the 3.10 parser). `compileall` at 3.10 catches syntax; the import step catches
  annotation/`importlib.resources` failures. It uses `pkgutil.walk_packages` with a prefix (`iter_modules`
  skips subpackage modules) and asserts a module count.
- **Only 3.10 resolves pandas 2.x**; 3.11+ get 3.x. Hence e2e runs 3.10 and 3.14. Both give identical
  row counts on every non-`huge` case.

## Repo layout

- No unit tests; `tests/e2e/` is the whole suite.
- **Deliberately wrong fixture lines — don't fix them:** `invalid_residues.txt` line 2 (`4Q5J:A:9999`,
  `test_invalid_1`, asserted `pocket_not_built`); `forced_pisa_mixed.txt` line 2 (`4Q5J:A`, no partner;
  `test_invalid_8`, asserted `invalid_entry`); `fsdb_mixed_target.txt` (a DB beside a structure, also via
  `job_file_fsdb_mixed.json` in `test_steps_22`; `test_steps_11`).
- **`fixtures/job_file.json` must never set a path** (the runner's `--results_dir` would not move it, so
  the case would write outside its dir). Its `align_count: 3` vs `test_settings_2`'s `5` shows priority
  only in `job_settings.json` (5 wins) — check by hand.
- **Delete `build/` before building or testing a wheel.** Stale, gitignored, but setuptools reuses
  `build/lib/`, so `pip install .` ships dead modules (`align.py`, `local_aligner.py`, a 3.12-only
  `pisa.py`, top-level `structure_fetcher.py`/`pisa_downloader.py` shadowing `downloads/`). Never edit it.

## As a library

`PocketMapper().search(...)`, a step's entry or core function in `steps`, or any component directly
(`qt_processor`, `downloads.*`, `structure_preprocessor`, `pockets.pocket_fetcher` or one builder,
`sequence_aligner`, `structure_aligner`, `foldseek`). No component or step takes a `Settings`; the
entry functions take its dict as `job_file`.

- **Records come from `steps.parse.parse_entries`**: the side checks parse makes, with no files written.
- **Superposition can be deferred**: `search(align_count=0)`, then
  `steps.superpose.superpose(job_file=settings)`, or `StructureAligner.align_structs` with
  `parse_entries` of the returned settings (FSDB targets also need `fsdb_path`; `steps.superpose`
  derives it). Records point at `pdb_dir`/`alphafold_dir`, which temp deletion never touches.
  `query_ids`/`target_ids`/`overwrite` exist for batching this.
