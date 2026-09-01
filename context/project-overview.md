# Project overview
Project implementation specifics. Cross-module and derived facts only. Anything a single docstring or comment already states appears here as a pointer to that site, never as a second copy.

## Pipeline

`main()` → `fire.Fire(PocketMapper())`, so **every public method on `PocketMapper` is a CLI subcommand** —
hence the leading underscores on all internals, to keep them out of fire's help. The seven steps of
`search()` are listed in the `pocketmapper.py` module docstring.

### Input grammar

`struct_info[:chain_info[:residue_info]]`; either side may instead be a file with one such string per line.
README's "Input format" table documents the forms; the `qt_processor.py` module docstring points at the two
methods that implement them.

One consequence neither states: a local-file entry like `4Q5J.cif.gz:B_F` resolves to `vdw`, not `pisa` —
`B_F` matches the vdw regex and PISA is PDB-only. That is how the mixed-input e2e fixtures reach the vdw
code.

### Pocket shape

Every pocket method returns a `pocket.Pocket` — the dataclass declares which fields exist, which are
optional and why, and `pocket_parser.parse_pocket_from_struct` shows how `seq_pos` and `whole_chain`
are derived. Residues live under `residues`, keyed by author seqid as a string.

One thing the class states that no producer would: `res_auth_ids` is not `list(residues)`. It is the
ordered residue list the comparison walks, and on the PISA path it is seeded from the interface while
`residues` is filled in chain order.

### Open searches

README's "Open searches" covers the output shape; `_retrieve_whole_chain_pockets` and
`_compare_pocket_pair` cover the per-pocket suppression of the `pocket_2_*` columns.

A `pocket_2` value is not guaranteed to be a target. A query and target sharing a chain share a
`preprocess_name`, so `compare_pockets` pairs every pocket on that chain with every other and some rows
carry a query-only `pocket_id` in `pocket_2`. `_align_structs` filters those out before its `.loc` lookup;
without that it raises a bare pandas `KeyError`.

### Foldseek-DB targets

When the target is a bundled Foldseek DB, `self.fsdb_target` is set: no target structures are fetched or
preprocessed, and `_align_structs` reconstructs target PDBs via `foldseek createsubdb` + `convert2pdb`.
What the target "pocket" is depends on the DB — `_expand_fsdb_pdb_targets` for a PDB DB, and
`pocket_comparison._synthesise_target_pocket` for any other. On the PDB path, **hits with no usable PISA
data are dropped**, not compared against a stand-in.

Three consequences of the PDB path that are not visible from any single file:

- **One `pocket_id` can sit behind two `preprocess_name`s** (`4q5j-assembly1_B` and `4q5j-assembly2_B` both
  resolve to `4Q5J:B_F`). The pocket is computed once and `compare_pockets`' `existing_calcs` scores only
  the first assembly's alignment row, so the transform used is whichever assembly Foldseek reported first.
- **Pockets come from the wwPDB asymmetric unit while Foldseek's `tseq` comes from the assembly.** These
  agree in the ordinary case (verified: a 4Q5J self-comparison through a PDB-named DB gives
  `overlap_count == pocket_len`, identity 1.0, RMSD ~1e-14), and the `MIN_SEQ_IDENTITY` guard catches them
  when they don't — a populated `incorrect_mapping.json` signals that an entry's assembly and AU numbering
  diverged.
- **`--align_struct_method pocket` is rejected for any Foldseek-DB target**, in `_configure_query_target`
  before anything is fetched. The rejection site gives the reason for both kinds of DB.

**No cap on how many hits get enriched**, by choice. `4Q5J:B_F` against the bundled `pdb` DB returns ~4,970
hits across ~3,620 entries, and PISA is fetched per entry behind a sleep, so the first run takes hours.
Reruns are cheap from the interface cache, and `_expand_fsdb_pdb_targets` logs both counts before starting
so the wait is legible. Add a cap here if that becomes untenable.

## Invariants

Breaking one of these generally produces silently wrong output rather than an error. Each is documented at
its code site; what follows is the map of where, plus the checks that live nowhere else.

- **`seq_pos` is the value everything hinges on** — declared on `pocket.PocketResidue`, set in
  `pocket_parser.parse_pocket_from_struct`, used in `pocket_comparison._map_pocket_into_alignment`. A new
  pocket method computing it any other way yields zero overlap with no error. Check it by comparing a
  pocket against itself: `overlap_count == pocket_len`.
- **`preprocess_name` is the alignment join key** — computed in `QTProcessor.parse_individual_qt`.
  Alignments are keyed by it, pockets by `pocket_id`, and `_compare_pockets_based_on_alignment` builds
  `preproc_to_ids` to bridge them.
- **Two tables have declared schemas** — `constants.ALIGNMENT_COLUMNS` and
  `pocket_comparison.POCKET_COMPARISON_COLUMNS`. A new column goes into the constant, never into one
  producer alone; see the note above `ALIGNMENT_COLUMNS`.
- **`compare_pockets` must not write to a `Pocket`** — stated on the `Pocket` class itself.
- **Aligned structures are named by `lib.safe_filename(query_id)`, not by `pocket_id`** — so
  `aligned_structures/*.pdb` filenames aren't greppable for an input string. Match on the `MOLECULE`
  records inside instead.
- **Two transform sources, chosen by `align_struct_method`** — `StructureAligner`'s class docstring names
  them; `pocket_comparison.parse_pocket_transform` is the only legitimate reader of the pocket transform
  and carries the measured evidence. Never hand a raw `p2_to_p1_*` cell to gemmi.

**Changing step 6 without changing behaviour**: capture `compare_pockets`' arguments from a real run and
diff old output against new. Nothing else covers that path.

## Logging and errors

Every log call must pass `extra={"stage": "..."}` or the record fails to format against the root formatter
(`PocketMapper.__init__`). Every call currently does; nothing enforces it.

Errors are `logging.critical(...)` then `raise PocketMapperError(...)`; `main()` catches and exits 1.

**No `exit()`/`sys.exit()` inside modules** — deliberately removed, which no code comment can show. The one
survivor is the CLI affordance `_check_help_search`, so `PocketMapper().search(help=True)` will kill a host
process — library callers must not pass it.

## Settings

A new option goes in **five hand-maintained places**: the `Settings` dataclass, the `cli_overrides` dict in
`_configure_workflow`, the `search()` signature, `HELP_MESSAGE` (`constants.py`) and the README's Options
list. Neither of the last two is generated from the dataclass. Miss one and the option is silently ignored.

Resolution order and the tri-state `foldseek` / `align_struct_method` settings are documented where they are
resolved — the `Settings` docstring and the `# 4b.` / `# 4c.` comments in `_configure_workflow`, which give
the reasons those call sites are load-bearing. Keep new resolution logic there.

## Python versions

Supported: **3.10 – 3.14**, verified by running the full e2e suite on each end. Four hand-maintained
places have to agree: `requires-python` in `pyproject.toml`, the `Programming Language :: Python` classifiers
beside it, `[tool.black] target-version`, and the README's Installation line. The `compat` CI job pins the
range in one more place, as a matrix.

- **The floor is 3.10 and going lower buys nothing.** Three `match` statements (`qt_processor.py` x2,
  `structure_fetcher.py`) and the PEP 604 `str | None` field annotations on `Pocket`, `PocketResidue`,
  `QTRecord` and `Settings` all require it. No module carries `from __future__ import annotations`, so those
  annotations are evaluated at import rather than deferred. Rewriting all of that for 3.9 would still fail:
  biopython requires >=3.10.
- **`compat` is what guards the floor, not `lint`.** flake8 parses with whatever interpreter runs it, so lint
  at 3.12 cannot see a 3.12-only construct. `compileall` at 3.10 is what catches syntax; the import step is
  what catches the annotation and `importlib.resources` failures that compileall cannot.
- **3.10 is the only version pip resolves to pandas 2.x** — 3.11 and up get pandas 3.x. That is why the e2e
  matrix covers 3.10 and 3.14 rather than the middle. Both produce identical comparison row counts across
  every non-`huge` case.
- The bundled Foldseek DB is resolved through `files("pocketmapper")`, not through the data directory, for
  the reason given at that call site in `qt_processor`.

## Repo layout

Each module's own docstring states its remit. Not stated anywhere in the code:

- There are no unit tests. `tests/e2e/` is the whole suite; the `pocketmapper-e2e` skill covers running it.
- `build/` and `dist/` are stale artifacts of an older version. Both are gitignored and untracked, so a
  fresh clone and CI never see them — but setuptools reuses `build/lib/` in place rather than clearing it,
  so on a machine that has one, `pip install .` silently ships whatever dead modules it still holds
  (`align.py`, `local_aligner.py`, `pisa.py`) on top of the current sources. `pisa.py` still carries the
  3.12-only f-string that `pisa_downloader.py` no longer does, so an import-everything check passes in CI
  and fails locally. Delete `build/` before building or testing a wheel; never edit `build/lib/pocketmapper/`.
- **Structure parsing is gemmi throughout** (`.cif.gz` on disk). Biopython is used only for pairwise
  alignment (`sequence_aligner.py`) and SVD superposition (`pocket_comparison.py`).
- `StructureFetcher` and `StructurePreprocessor` share a required call order that nothing enforces; both
  classes' docstrings say so. Both cache on bare filenames and write through a `.part` file, for reasons
  their `update_cache` docstrings give.

## As a library

fire wraps `PocketMapper` only at the `main()` boundary, so nothing here needs a terminal. Two levels of
entry: `PocketMapper().search(...)` does the same work as the CLI, or drive a component directly —
`qt_processor`, `structure_fetcher`, `structure_preprocessor`, `pisa_downloader`, `pisa_parser`,
`sequence_aligner`, `structure_aligner`, `pocket_calculator` are each separately usable.

- **A component reaching into a `Settings` can't be used without building one, and hides which fields it
  depends on** — so no component takes one. The `Settings` is unpacked at each call site in
  `pocketmapper.py` into the values that component needs.
- **`pocketmapper/__init__.py` only exports `main` and `__version__`.** Submodules are reachable as
  `pocketmapper.lib` etc. only as a side effect of `pocketmapper.pocketmapper` importing them — always use
  explicit `from pocketmapper.<module> import <name>`.
- **Always call `Settings(...).resolve_paths()`** if you build one yourself — the failure mode is in that
  method's docstring. `search()` does this for you.
- **`search()` has global side effects**: `logging.config.dictConfig` reconfigures the *root* logger and
  stomps on a host app's logging setup, and `_delete_tmp` `shutil.rmtree`s
  `query_dir`/`target_dir`/`foldseek_tmp_dir` at the end (marked `# TODO this is unsafe`).
- Results come back through files — `search()` returns `None`, so read `pocket_comparison.tsv` /
  `alignment.tsv` from `results_dir` (paths available on `Settings`).
