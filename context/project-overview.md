# Project overview
Project implementation specifics. Cross-module and derived facts only. Anything a single docstring or comment already states appears here as a pointer to that site, never as a second copy.

## Pipeline

`cli.py` holds the argparse parser and the console-script `main()`; it is **the only module that knows
about argv or exit codes**, and `search()` is its one subcommand. The seven steps of `search()` are listed
in the `pocketmapper.py` module docstring.

Two parsing details are load-bearing and documented at the parser: query and target are required
positionals with no `--query`/`--target` spelling — so a settings file's `query`/`target` can never
win on the CLI path, and those two fields exist for library callers alone — and every option defaults
to `None` rather than to a `Settings` default, which is what leaves the JSON settings file
overridable.

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

README's "Open searches" covers the output shape; `retrieve_whole_chain_pockets` and
`compare_pocket_pair` cover the per-pocket suppression of the `pocket_2_*` columns.

A `pocket_2` value is not guaranteed to be a target. A query and target sharing a chain share a
`preprocess_name`, so `compare_pockets` pairs every pocket on that chain with every other and some rows
carry a query-only `pocket_id` in `pocket_2`. `align_structs` filters those out before its `.loc` lookup;
without that it raises a bare pandas `KeyError`.

### Foldseek-DB targets

When the target is a bundled Foldseek DB, `self.fsdb_target` is set: no target structures are fetched or
preprocessed, and `align_structs` reconstructs target PDBs via `foldseek createsubdb` + `convert2pdb`.
What the target "pocket" is depends on the DB — `expand_fsdb_pdb_targets` for a PDB DB, and
`pocket_comparison.synthesise_target_pocket` for any other. On the PDB path, **hits with no usable PISA
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
- **`--align_struct_method pocket` is rejected for any Foldseek-DB target**, in `configure_query_target`
  before anything is fetched. The rejection site gives the reason for both kinds of DB.

**Target residue ids are UniProt coordinates when the DB ships an offset table.** A non-PDB DB's
entries are domains carved out of UniProt sequences, so `synthesise_target_pocket`'s residue *labels*
are renumbered through `offset_table.tsv` — resolved in `compare_pockets_based_on_alignment`, applied
in `synthesise_target_pocket`, both of which carry the reasoning. Three things follow that no single
file states:

- **Only the bundled `human_domains` DB is renumbered.** `foldseek.bundled_human_domains_offset_table`
  matches on the resolved DB path, so a DB you supply yourself keeps 0-indexed positions within the
  entry — logged at INFO, because the same column then means different things on different runs.
- **`pocket_2_overlap_ids` is the only column affected**, since a whole-chain pocket already suppresses
  the other `pocket_2_*` columns and has no coordinates to superpose. Verified: against a
  same-environment baseline, a `human_domains` search changes that column and nothing else.
- **The table and the DB are now coupled.** A hit whose entry is missing from the table, or whose spec
  is shorter than the alignment reaches, aborts the run rather than falling back — a per-row fallback
  would mix two coordinate systems inside one column with nothing in the row to tell them apart.
  Refresh the table whenever `BUNDLED_HUMAN_DOMAINS_DB` moves.

**The bundled DB ships without its `.source` file**, and a refreshed one must be stripped the same way.
Foldseek's `createdb` writes `.source` alongside `.lookup`, but it duplicates the same key-to-name mapping
and nothing reads it: verified by running `easy-search`, `createsubdb` and `convert2pdb` against a copy with
it removed. `.lookup` is the one that must survive — `align_structs` reads it to turn entry names into
database keys. Dropping `.source` saves 1.7 MB in the repo and in both distributions.

The DB is otherwise at its floor. `_ca` is 70 of its 98 MB, holding 11.2M residues at 6.33 bytes each, which
is `foldseek createdb --coord-store-mode 2` (uint16 deltas), the default and the smallest of the three modes.
Being delta-encoded it is already near entropy — zstd -19 takes only 16% off it — and no foldseek module
writes or reads a compressed structure DB, so there is nothing to gain by compressing what ships.

**No cap on how many hits get enriched**, by choice. `4Q5J:B_F` against the bundled `pdb` DB returns ~4,970
hits across ~3,620 entries, and PISA is fetched per entry behind a sleep, so the first run takes hours.
Reruns are cheap from the interface cache, and `expand_fsdb_pdb_targets` logs both counts before starting
so the wait is legible. Add a cap here if that becomes untenable.

## Invariants

Breaking one of these generally produces silently wrong output rather than an error. Each is documented at
its code site; what follows is the map of where, plus the checks that live nowhere else.

- **`seq_pos` is the value everything hinges on** — declared on `pocket.PocketResidue`, set in
  `pocket_parser.parse_pocket_from_struct`, used in `pocket_comparison.map_pocket_into_alignment`. A new
  pocket method computing it any other way yields zero overlap with no error. Check it by comparing a
  pocket against itself: `overlap_count == pocket_len`. It is also **not** the reported residue id:
  `synthesise_target_pocket` keys its residues by UniProt position while leaving `seq_pos` the
  0-indexed alignment coordinate, and that separation is the only reason renumbering is safe.
- **`preprocess_name` is the alignment join key** — computed in `QTProcessor.parse_individual_qt`.
  Alignments are keyed by it, pockets by `pocket_id`, and `compare_pockets_based_on_alignment` builds
  `preproc_to_ids` to bridge them.
- **`chain_info` is split in exactly one place** — `lib.split_chain_info`, which nine call sites across
  seven modules now share. Four of them used to index the string (`chain_info[0]`), which is the domain
  chain only while a chain id is one character. `QTProcessor`'s regexes guarantee that, but a forced
  `--query_pocket_method` / `--target_pocket_method` skips them entirely, so `4Q5J:AA_BB` silently became
  chain `A`. Never re-derive a domain or motif chain inline.
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

That convention now holds on the Foldseek path too, which is most of what `foldseek.run_foldseek` buys:
five of the six invocations used to be a bare `subprocess.run(..., check=True)`, so a failing Foldseek
surfaced as a `CalledProcessError` traceback rather than a message. Exit code was 1 either way.

**No `exit()`/`sys.exit()` inside modules** — deliberately removed, which no code comment can show. There
are now none: the last survivor was `_check_help_search`, deleted along with `search()`'s `help` parameter
when argparse took over `--help`. The only `sys.exit` in the package is in `cli.py`'s `main()`, which is
the boundary and is meant to have one.

## Settings

**Every `Settings` field is reachable from the command line**, and the settings JSON sets nothing the
CLI cannot. The file is a convenience for keeping a long invocation reproducible, never the only route
to a setting; the layering that makes it one is in `configure_workflow`.

A new option goes in **five hand-maintained places**: the `Settings` dataclass, the `search()` signature
together with the `cli_overrides` dict directly beneath it (one site — the dict mirrors the signature and
sits next to it precisely so the two cannot drift), the parser in `cli.py`, `main()`'s kwarg block in the
same file, and the README's Options tables. None is generated from the dataclass. Miss one and the option
is silently ignored — `main()` is the one that reads like boilerplate and is easiest to forget.

Nothing enforces the agreement, but it is checkable in a few lines: `dataclasses.fields(Settings)`,
`inspect.signature(PocketMapper.search)`, the `cli_overrides` keys, the subparser's `_actions` dests and
`main()`'s `x=args.x` lines must all name the same fields (modulo `settings`, which is a parser-side
spelling rather than a field). The `query` and `target` positionals carry those dests, so they line up
with the rest.

Options are grouped by lifetime in both places a human reads them — argparse's argument groups in
`build_parser`, and the README's matching subsections. The path fields alone roughly double the
option count, so leaving them ungrouped would bury `--foldseek` and `--query_pocket_method` among
them. The four groups are `aligned structure options`, `cache options` (what survives a run),
`out options` (what the run produces) and `temp options` (what `delete_tmp` removes at the end);
`--cache_dir` and `--results_dir` head the group whose defaults derive from them. Adding a path
setting means picking one of those groups in both places.

Lifetime wins over where the default comes from, and `--foldseek_tmp_dir` is the one place the two
disagree: it defaults under `cache_dir` but a Foldseek run deletes it, so it sits in `temp options`
rather than with the caches. Moving it back to `cache options` on the strength of its default would
put a directory that does not survive the run under a heading that promises it does.

`search --help` *is* generated, from the parser's `help=` strings; `constants.CLI_SEARCH_EPILOG`
carries only the examples, which is all argparse cannot produce. It hangs off the `search` subparser
alone; the bare `pocketmapper --help` is the subcommand list and nothing more.

Resolution order and the tri-state `foldseek` / `align_struct_method` settings are documented where they are
resolved — the `Settings` docstring and the `# 4b.` / `# 4c.` comments in `configure_workflow`, which give
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
  the reason given at that call site in `foldseek`.

## Repo layout

Each module's own docstring states its remit. Not stated anywhere in the code:

- There are no unit tests. `tests/e2e/` is the whole suite; the `pocketmapper-e2e` skill covers running it.
- **`fixtures/settings_paths.json` is deliberately wrong, and JSON cannot say so.** It sets a `results_dir`
  that must never be used: `test_settings_2` relies on the runner appending its own `--results_dir` after
  the case args, so if CLI-over-file layering ever broke, `pocket_comparison.tsv` would land at the file's
  path and the existing assertion would fail. Its `align_count: 3` is the other half — a value nothing on
  the command line sets, so seeing it in `job_settings.json` proves the file was read at all. Change either
  value and the case stops testing anything.
- **What `test_settings_1` cannot catch.** The runner only ever asserts on `$case_out/pocket_comparison.tsv`,
  so a path option that argparse accepts and something downstream silently drops still passes. The case
  catches a rejected or crashing flag and nothing subtler; the five-way agreement check under "Settings" is
  what covers the rest, by hand.
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

The CLI is confined to `cli.py`, so nothing else here needs a terminal. Two levels of
entry: `PocketMapper().search(...)` does the same work as the CLI, or drive a component directly —
`qt_processor`, `structure_fetcher`, `structure_preprocessor`, `pisa_downloader`, `pisa_parser`,
`sequence_aligner`, `structure_aligner`, `pocket_calculator`, `foldseek` are each separately usable.

- **A component reaching into a `Settings` can't be used without building one, and hides which fields it
  depends on** — so no component takes one. The `Settings` is unpacked at each call site in
  `pocketmapper.py` into the values that component needs.
- **`pocketmapper/__init__.py` only exports `main` and `__version__`**, now from `pocketmapper.cli`.
  Submodules are reachable as `pocketmapper.lib` etc. only as a side effect of that import chain — always
  use explicit `from pocketmapper.<module> import <name>`.
- **Always call `Settings(...).resolve_paths()`** if you build one yourself — the failure mode is in that
  method's docstring. `search()` does this for you.
- **`search()` has global side effects**: `logging.config.dictConfig` reconfigures the *root* logger and
  stomps on a host app's logging setup, and `delete_tmp` `shutil.rmtree`s
  `query_dir`/`target_dir`/`foldseek_tmp_dir` at the end unless `delete_tmp=False`, which keeps all
  three. The rmtree is guarded rather than unconditional: all three are settable, so `lib.is_within`
  skips (with a warning) any that does not resolve under `cache_dir` or `results_dir`. The guard
  bounds the damage from a mistyped path; it is not a reason to point those settings at a directory
  you care about.
- Results come back through files — `search()` returns `None`, so read `pocket_comparison.tsv` /
  `alignment.tsv` from `results_dir` (paths available on `Settings`).
