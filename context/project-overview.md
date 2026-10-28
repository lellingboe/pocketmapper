# Project overview
Cross-module and derived facts only. Anything one docstring or comment already states is left to that site
and not copied here.

## Project Status

- This is an early development project, features may change or break between patches.
- Names/location of outputs may change, but typically the underlying data should remain fairly consistent.
- The project is currently going through refactoring and reorganization, some things may be inconsistent
  while this process is happening

## Where to start

- `search.py` module docstring: the steps of `search()` and what each writes under `results_dir`;
  each step's entry/core layers, settings layering, and how every later step reads its records from
  `entries.tsv` (no records are handed on in memory).
- `cli.py`: the only module that knows about argv or exit codes; why every option defaults to None.
- `lib/entries.py`: `entries.tsv`, why a dropped entry warns only once, and why a rerun step gives the same
  table.
- `parse/qt_processor.py`: the input grammar; the two tables in `QTProcessor.__init__` are the whole of it.
- `parse/CLAUDE.md`, `pockets/CLAUDE.md`, `lib/CLAUDE.md`: read before changing those packages.

## Invariants

Breaking these gives silently wrong output. Each is explained at its code site; this is the map, plus
checks that live nowhere else.

- **`seq_pos`**: declared on `PocketResidue`, set only in `parse_pocket_from_struct` (every builder reaches it,
  vdw included) as the index in `lib.residues.ca_residues`, and read in `map_pocket_into_alignment`. -1 marks
  a requested residue the chain has but `ca_residues` leaves out (no CA atom, or not polymer). Computed any
  other way, it gives zero overlap and no error. Check: a pocket compared with itself gives
  `overlap_count == pocket_len` (`test_self_comparison`).
- **`res_auth_ids` ascending** in every producer (`Pocket`); nothing checks it. `pockets.json` round-trips
  this order and `seq_pos` (`read_pockets_file`).
- **Residue letters mirror Foldseek, not gemmi** (`lib.constants.FOLDSEEK_AA_CODES`). They compare
  case-insensitively (`compare_pockets`), because Foldseek lowercases some residues; 4ER4, 2ER6 and 2ER9
  are lowercased wholesale.
- **Every chain walk reads `first_conformer()`**, and a new walk must too. The CA walk is
  `lib.residues.ca_residues` (`parse_pocket_from_struct`, `vdw_pockets`' domain side, `SequenceAligner`); it
  needs the entity types `setup_entities()` sets (`lib.structure_manager.load_structure` calls it). vdw's motif
  side reads every atom of `get_polymer().first_conformer()`. `test_core_10` (`network`) covers it.
- **`preprocess_name` is the alignment join key** (`QTProcessor.build_record`), bridged to `pocket_id` by
  `pockets.json` `chains`. Only parse computes it (hashing local files); later steps read it from
  `entries.tsv`, so `align` and `pockets` must run on the same table, and a caller that re-runs
  `parse_entries` needs unchanged inputs. `compare.compare.check_coverage` catches a mismatch.
  `test_preprocess_name_local_contents` and `test_local_7` cover local-file names.
- **`chain_info` is split only by `lib.residues.split_chain_info`.** `chain_info[0]` breaks on multi-character
  chains; the patterns reject them, but library callers bypass the patterns.

**Changing the comparison without changing behaviour**: capture `compare_pockets`' arguments from a real
run and diff old vs new output. The tests assert invariants (self-comparison overlap, rows present, steps
matching search), not column values, so nothing else covers it.

## Logging and errors

- The package never touches the root logger: `__init__` gives `lib.constants.PACKAGE_LOGGER` a `NullHandler`.
- Handlers are added only in `cli()` (stdout) and `lib.logs.log_to_file` (`info.log`, for `search` and every
  step), each removed in a `finally`.
- Errors go through `PocketMapperError`, which `cli()` turns into exit 1. Foldseek runs through
  `run_foldseek` so its failures raise it too.

## Python versions

**3.10 – 3.14**; CI's `e2e` job runs pytest with foldseek, `network` tests included, at both ends. Four
places must agree: `requires-python`, the classifiers, `[tool.ruff] target-version`, README Installation; CI
`compat` matrix is a fifth.

- **Floor is 3.10**: `match` statements, PEP 604 dataclass annotations (no `from __future__ import
  annotations`), and biopython itself.
- **`compat` guards the floor, not `lint`** (ruff flags some post-3.10 syntax from `target-version`, but is
  not the 3.10 parser). `compileall` at 3.10 catches syntax; the import step catches
  annotation/`importlib.resources` failures. It uses `pkgutil.walk_packages` with a prefix (`iter_modules`
  skips subpackage modules) and asserts a module count.
- **Only 3.10 resolves pandas 2.x**; 3.11+ get 3.x. Hence CI's `e2e` job runs 3.10 and 3.14. Both gave
  identical row counts on every non-`huge` case of the former bash suite.

## Repo layout

- `tests/unit/` and `tests/integration/` are the suite (`python -m pytest`). Integration tests run the CLI
  in-process, offline from a copy of `tests/data/seed_cache/` unless marked `network`; the markers and
  their env vars are in `pyproject.toml` and `tests/integration/conftest.py`.
- **Deliberately wrong fixture line — don't fix it:** `tests/integration/fixtures/invalid_residues.txt`
  line 2 (`4Q5J:A:9999`, `test_steps_23`, asserted `pocket_not_built`).
- **`tests/integration/fixtures/4Q5J.cif.gz` is also the seed's PDB structure** (`seed_cache` copies it
  into `pdb_structures/`), and `testfile.txt` names it by a relative path; don't move or rename it.
- **Delete `build/` before building or testing a wheel.** Stale, gitignored, but setuptools reuses
  `build/lib/`, so `pip install .` ships dead modules (`align.py`, `local_aligner.py`, a 3.12-only
  `pisa.py`, top-level `structure_fetcher.py`/`pisa_downloader.py`). Never edit it.

## As a library

`PocketMapper().search(...)`, a step's entry or core function in `<step>.<step>`, or any component directly
(`parse.qt_processor`, `lib.structure_manager` (`StructureManager.get` fetches and loads one structure),
`pockets.pisa_downloader`, `align.structure_preprocessor`, `pockets.pocket_fetcher.fetch_pockets` or one
builder, `pockets.pisa` (`load_interfaces`, `interface_partners`, `pisa_pockets_from_records`),
`pockets.vdw.vdw_contact_residues`, `align.sequence_aligner`, `superpose.structure_aligner`, `lib.foldseek`).
No component or step takes a `Settings`; the entry functions take its dict as `job_file`.

- **Records come from `parse.parse.parse_entries`**: the side checks parse makes, with no files written.
  A finished run's own records, its dropped entries left out, are
  `lib.entries.entry_sides(lib.entries.read_entries(entries_path))`.
- **Superposition can be deferred**: `search(align_count=0)`, then
  `superpose.superpose.superpose(job_file=settings)`, or `StructureAligner.align_structs` with the records
  of the run's `entries.tsv`, as above (FSDB targets also need `fsdb_path`; `superpose.superpose` derives it).
  Records point at `pdb_dir`/`alphafold_dir`, which temp deletion never touches.
  `query_ids`/`target_ids`/`overwrite` exist for batching this.
