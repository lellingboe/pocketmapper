# Project overview
Cross-module facts and pointers only (`coding-standards.md` "Context docs").

**Status:** early development, mid-refactor; output names and locations may change between patches. Not
ready for production use.

## Where to start

- `search.py` module docstring: the steps, what each writes, and how every later step reads its records
  from `entries.tsv`.
- `cli.py` module docstring: argv, exit codes, why every option defaults to None.
- `lib/options.py`: the option registry. `context/settings.md` says what a new option needs.
- `lib/entries.py` module docstring: `entries.tsv`.
- `parse/qt_processor.py`: the input grammar, mostly the two tables in `QTProcessor.__init__`.

## Invariants

Breaking these gives silently wrong output. Each is explained at its code site; this is the map, plus
checks that live nowhere else.

- **`seq_pos`** (`PocketResidue.seq_pos`): every built pocket's `seq_pos` is set in
  `parse_pocket_from_struct`, as the index in `lib.residues.ca_residues`, and read in
  `map_pocket_into_alignment`. The one other setter is `compare.pocket_comparison.synthesise_target_pocket`,
  for a pocket that is never built. `pockets.pocket_files.check_pocket` enforces the rule on every pocket
  written to or read from `pockets.json`. Check: a pocket compared with itself gives
  `overlap_count == pocket_len` (`test_self_comparison`).
- **Residue letters mirror Foldseek, not gemmi** (`lib.constants.FOLDSEEK_AA_CODES`) and compare
  case-insensitively; 4ER4, 2ER6 and 2ER9 are lowercased wholesale.
- **Every chain walk reads `first_conformer()`**, and a new walk must too. The walks: `lib.residues.ca_residues`
  (used by `parse_pocket_from_struct`, `vdw_pockets`' domain side and `SequenceAligner`), which needs the
  entity types `lib.structure_manager.load_structure` sets; and vdw's motif side. `test_core_10` (`network`)
  covers it.
- **`preprocess_name` is the alignment join key**, bridged to `pocket_id` by `pockets.json` `chains`. Parse
  computes it (`QTProcessor.build_record`, hashing local files) and `PocketManager.expand_pdb_hits` replaces
  it with the Foldseek hit name; every other step reads it from `entries.tsv`. So `align` and `pockets` must
  run on the same table, and a caller that re-runs `parse_entries` needs unchanged inputs.
  `compare.compare.check_coverage` catches a mismatch. `test_preprocess_name_local_contents` and
  `test_local_7` cover local-file names. More in `parse/CLAUDE.md`.
- **A pocket's `ca_sequence` must equal the aligner's sequence for its chain**
  (`compare.pocket_comparison.first_sequence_mismatch`); a mismatch drops the pocket into
  `incorrect_mapping.json`. Exactness is the point: a threshold lets a shift near the chain end through.
- **`chain_info` is split only by `lib.residues.split_chain_info`.** `chain_info[0]` breaks on multi-character
  chains; the patterns reject them, but library callers bypass the patterns (`pockets/CLAUDE.md`).

**Changing the comparison without changing behaviour**: capture `compare_pockets`' arguments from a real
run and diff old vs new output. The tests mostly assert invariants (self-comparison overlap, rows present,
steps matching search), not column values, so nothing else covers it.

## Logging and errors

- The package never touches the root logger (`pocketmapper/__init__.py`).
- Handlers are added only in `cli()` and `lib.logs.log_to_file`, each removed in a `finally`.
- Errors go through `PocketMapperError`, which `cli()` turns into exit 1; `lib.foldseek.run_foldseek` raises
  it for Foldseek.

## Python versions

**3.10 – 3.14.** Five places must agree: `requires-python`, the classifiers, `[tool.ruff] target-version`,
README Installation, and the CI `compat` matrix; the e2e matrix runs both ends.

- **Floor is 3.10**: PEP 604 dataclass annotations (no `from __future__ import annotations`).
- **`compat` guards the floor, not `lint`**: ruff is not the 3.10 parser; `compileall` at 3.10 catches syntax.
- **Only 3.10 resolves pandas 2.x**; 3.11+ get 3.x. Hence the e2e jobs run 3.10 and 3.14.

## Repo layout

- Tests: `tests/unit/` and `tests/integration/`; markers and env vars in `pyproject.toml` and
  `tests/integration/conftest.py`.
- **Deliberately wrong fixture line — don't fix it:** `tests/integration/fixtures/invalid_residues.txt`
  line 2 (`4Q5J:A:9999`, `test_steps_23`).
- **`tests/integration/fixtures/4Q5J.cif.gz` is also the seed's PDB structure**, and `testfile.txt` names it
  by a relative path; don't move or rename it.
- **Delete `build/` before building or testing a wheel**: setuptools reuses a stale `build/lib/`.
- **Imports run one way.** No `lib` module imports a step package, `search` or `cli`, and no step imports
  another step's entry module, only its components (`tests/unit/test_imports.py`).

## As a library

`PocketMapper().search(...)`, a step's entry or core function in `<step>.<step>`, or a component:
`parse.qt_processor`, `lib.structure_manager`, `pockets.pisa_downloader`, `align.structure_preprocessor`,
`pockets.pocket_manager.PocketManager` or one builder, `pockets.pocket_files`, `pockets.pisa`,
`pockets.vdw.vdw_contact_residues`, `align.sequence_aligner`, `superpose.structure_aligner`, `lib.foldseek`.
No component or step takes a settings object; the entry functions take a job dict as `job_file`.

- **Records**: `parse.parse.parse_entries` makes them without writing files. A finished run's own records,
  its dropped entries left out, are `lib.entries.entry_sides(lib.entries.read_entries(entries_path))`.
- **Superposition can be deferred**: `search(align_count=0)`, then `superpose.superpose.superpose(job_file=...)`,
  or `StructureAligner.align_structs` with the records above (FSDB targets also need `fsdb_path`). The
  structures it reads sit in the cache dirs, which temp deletion never touches.
