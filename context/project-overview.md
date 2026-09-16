# Project overview
Project implementation specifics. Cross-module and derived facts only. Anything a single docstring or comment already states appears here as a pointer to that site, never as a second copy.

## Pipeline

`cli.py` holds the argparse parser and the console-script `cli()`; it is **the only module that knows
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

**Two tables in `QTProcessor.__init__` are the whole grammar**, and both directions read them:
`pocket_methods` maps a method to its pocket-info pattern and to the phrase a warning uses, and
`struct_type_pocket_methods` maps a structure type to the methods it supports, in the order they are tried.
`determine_pocket_method` takes the first method whose pattern matches; `validate_pocket_method` checks an
entry against the pattern of the method it was *given*. That second direction is what a forced
`--query_pocket_method` / `--target_pocket_method` goes through, so a forced method is now constrained
exactly as an inferred one is and the two cannot drift. It runs for inferred methods too, where it is a
tautology, so the invariant holds for every record rather than for the forced ones alone.

Two things that follow:

- **The vdw pattern's partner chain is required, and was not always.** It used to be optional, which
  inference could never exercise — `whole_chain` is tried first and claims any bare chain — but which made
  a forced `vdw` on `4Q5J:A` validate and then reach gemmi as chain `None`. Tightening it is inference-neutral;
  verified by replaying the old per-struct_type ladder against the new loop over 400 generated entries.
- **A rejected entry is skipped, not fatal.** `validate_pocket_method` warns and returns False, the record
  is dropped, and `configure_query_target` raises only when a side ends up empty — so one bad line in a
  batch file does not stop the rest. The exception is an unrecognised method *name*, which is a whole-run
  setting rather than one entry: `process_qt_cmdline_input` raises on it before parsing anything.

### Pocket shape

Every pocket method returns a `pocket.Pocket` — the dataclass declares which fields exist, which are
optional and why, and `pocket_parser.parse_pocket_from_struct` shows how `seq_pos` and `whole_chain`
are derived. Residues live under `residues`, keyed by author seqid as a string.

One thing the class states that no producer would: `res_auth_ids` is not `list(residues)`. It is the
ordered residue list the comparison walks, and on the PISA path it is seeded from the interface while
`residues` is filled in chain order.

**Every producer must emit it in ascending residue order**, and nothing checks that. `overlap_ids`
returns each side's ids in its own `res_auth_ids` order and `superpose` pairs the two lists position
for position, so a pocket ordered any other way is superposed against the wrong residues -- wrong
`rmsd`, `ca_dists` and transforms, with `overlap_count` and every identity column still correct,
and no warning. Only the passthrough method takes its order from user input; the rest walk the chain
(`pocket_calculator`, the whole-chain path) or sort (`pisa_parser`, `retrieve_passthrough_pockets`).

### Open searches

README's "Open searches" covers the output shape; `retrieve_whole_chain_pockets` and
`compare_pocket_pair` cover the per-pocket suppression of the `pocket_2_*` columns.

A `pocket_2` value is not guaranteed to be a target. A query and target sharing a chain share a
`preprocess_name`, so `compare_pockets` pairs every pocket on that chain with every other and some rows
carry a query-only `pocket_id` in `pocket_2`. `StructureAligner.align_structs` filters those out before
looking a target record up; without that it raises a bare pandas `KeyError`.

### Foldseek-DB targets

When the target is a bundled Foldseek DB, `self.fsdb_target` is set: no target structures are fetched or
preprocessed, and `foldseek.extract_fsdb_structures` reconstructs target PDBs via `createsubdb` +
`convert2pdb` for whichever entries step 7 selected.
Of the five Foldseek subcommands the package runs, **`createsubdb` is the only one that takes no
`--threads`** — it accepts just `--subdb-mode`, `--id-mode` and `-v`, and passing the flag makes it exit
non-zero. `run_foldseek` therefore adds no flags of its own; every caller builds its own argument list.

`StructureAligner.align_structs` reaches that helper through one `fsdb_path` argument, and reads a
selected target id two ways depending on whether it was given any target records at all — records mean
the ids are `pocket_id`s to map through `preprocess_name` (the PDB DB), no records mean the ids are
entry names themselves (any other DB). The choice is made once for the whole call rather than per id,
so one run's `fsdb_structures/` can never mix entries resolved both ways.

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

- **Only the bundled `human_domains` DB is renumbered.** It is the one entry in
  `foldseek.bundled_foldseek_dbs` with a non-None `offset_path`, and `compare_pockets_based_on_alignment`
  looks that up by the target's resolved DB path rather than by name — so a DB you supply yourself keeps
  0-indexed positions within the entry, logged at INFO because the same column then means different
  things on different runs.
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
it removed. `.lookup` is the one that must survive — `extract_fsdb_structures` reads it to turn entry
names into database keys. Dropping `.source` saves 1.7 MB in the repo and in both distributions.

The DB is otherwise at its floor. `_ca` is 70 of its 98 MB, holding 11.2M residues at 6.33 bytes each, which
is `foldseek createdb --coord-store-mode 2` (uint16 deltas), the default and the smallest of the three modes.
Being delta-encoded it is already near entropy — zstd -19 takes only 16% off it — and no foldseek module
writes or reads a compressed structure DB, so there is nothing to gain by compressing what ships.

**No cap on how many hits get enriched**, by choice. `4Q5J:B_F` against the bundled `pdb` DB returns ~4,970
hits across ~3,620 entries, and PISA is fetched per entry behind a sleep, so the first run takes hours.
Reruns are cheap from the interface cache, and `expand_fsdb_pdb_targets` logs both counts before starting
so the wait is legible. Add a cap here if that becomes untenable.

## Downloads

`downloads/` holds what fetches bytes over HTTP, and nothing else. `foldseek.py` and
`PocketMapper.fetch_missing_fsdb` stay outside it: they download a database by shelling out to the
`foldseek` binary, so they share none of the machinery here, and moving them would split `foldseek.py`
away from `check_foldseek` and `bundled_foldseek_dbs`, which `qt_processor` imports for reasons
unrelated to downloading.

`lib_download` offers **two entry points, not one**, because the two things the package fetches differ
in every operational respect. `download_file` is for bulk structure files: unpaced, stateless, and
called from `download_missing_structures`' pool. `download_api` is for REST endpoints: it
paces its requests and, on a transient failure, doubles that host's pacing and never lowers it
again. Both write through a `.part` file and share one retry test.

Four consequences no single file states:

- **The download pool ignores `--threads`, on purpose.** `fetch_missing_structures` builds
  `StructureDownloader` with its default width, `DOWNLOAD_WORKERS` (8), so `--threads` governs Foldseek
  alone. A worker here waits on a socket rather than on a core, so a thread count says nothing useful
  about how wide the pool should be. A library caller can still pass its own `max_workers`.
- **The pacing delay and the backoff delay are the same number.** That is what "carry the backoff
  forward" means here: the escalation one retry needed becomes the pace of every later request to that
  host. It is also why `download_missing_summaries` and `download_missing_assemblies` no longer sleep
  themselves — the helper owns pacing, and a caller-side sleep would double it.
- **The delay registry is module-level and never decays**, so it outlives any one `PisaDownloader` —
  which matters, because `download_pisa_interfaces` builds a fresh one on each of its two calls per run.
  It equally outlives a whole `search()`, so a library caller running several in one process carries an
  elevated delay across all of them; `reset_host_delays` is the escape hatch.
- **Only 5xx and 408/425/429 are retried.** A whitelist among 4xx rather than a blacklist of 404, so an
  unrecognised 4xx costs one request instead of the whole budget. The old PISA code caught bare
  `Exception`, so every entry PISA lacked cost 5 requests and ~3.75s of sleeping — on the full-PDB path
  that is thousands of entries.

A leftover `.part` is inert in every cache directory: `download_missing_interfaces` globs `*.json`,
which cannot match `x.json.part`, the other PISA stages check an exact path, and
`StructureDownloader` tests for an exact `.cif.gz` destination.

**The PISA failure report is the caller's file, not the downloader's.** `download_missing_interfaces`
returns what each stage could not handle and writes `error_path` only when there is something to write;
`download_pisa_interfaces` hands both of its two calls per run the same
`pocket_dir/pisa_responses/errors.json`. So a second call with failures replaces the first call's
report, and a clean second call leaves the first's file in place. Date it by its mtime, not by its
existence.

Two things that file does not say about itself:

- **Most of what it lists under `assembly_parsing` is routine, not broken.** An interface is skipped
  for a multi-character chain id, and those dominate: over a 13,063-assembly cache, 25,733 of 50,333
  interfaces have one, while none had a molecule count other than two. The entries are
  distinguishable only by shape — `<pdb_code>_<assembly_id>_<interface_id>` for a skip,
  `<pdb_code>_<assembly_id>` or a `_parse_error` suffix for a genuine failure.
- **An entry whose interfaces are all skipped gets no `<pdb_code>.json`.** It therefore never enters
  the interface cache and is reparsed on every later run — from cached assemblies, so at no request
  cost, but it is also why such an entry reappears in every report.

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
  chain only while a chain id is one character. `QTProcessor`'s patterns guarantee that on both the
  inferred and the forced path now (see "Input grammar"), so `4Q5J:AA_BB` is rejected rather than silently
  becoming chain `A` — but the split stays in one place regardless, because a library caller can build a
  record itself and reach the same call sites. Never re-derive a domain or motif chain inline.
- **A passthrough pocket's `res_auth_ids` are all keys of its `residues`** — enforced in
  `retrieve_passthrough_pockets`, which skips the whole entry when they are not.
  `map_pocket_into_alignment` and `describe_pocket` both index `residues` by every `res_auth_ids` id,
  so an id the chain cannot supply used to surface as a `KeyError` out of `compare_pockets`' re-raise
  -- one typo aborting the run. Syntax and repeats are caught earlier, in
  `QTProcessor.parse_residue_info`, before anything is fetched. The repeat is the dangerous one: it
  reached `res_auth_ids` twice and paired the two sides' overlap lists off by one, so `pocket_len`,
  `jaccard_index`, `rmsd` and `ca_dists` all came out wrong with nothing in the row to show it.
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
- **Step 7 is `StructureAligner.align_structs`, not a pipeline method** — `PocketMapper.align_structs`
  only unpacks the `Settings` and the two `fsdb_*` flags into it. Selection, transform lookup and
  writing all live in the component, so a change to any of them belongs there and is reachable without
  running a search.

**Changing step 6 without changing behaviour**: capture `compare_pockets`' arguments from a real run and
diff old output against new. Nothing else covers that path.

## Logging and errors

The root format interpolates `%(stage)s` (`constants.LOG_FORMAT`), which is not a stock LogRecord
attribute. `lib.StageFilter` supplies it from `record.funcName` for any record that arrives without one,
so a missing `stage` degrades to the emitting function's name instead of failing to format. It is
attached **to the handlers**, in both places a handler is built — `PocketMapper.__init__` and the
`configure_logging` dictConfig. Handler-level rather than logger-level because a handler filter also
sees records propagating up from third-party loggers, which never pass an `extra`; a root-logger filter
would let those reach the formatter unprotected.

Three consequences, none visible from one file:

- **The `__init__` handler is load-bearing, not a placeholder.** `configure_workflow` can `logging.critical`
  on a bad settings file *before* it calls `configure_logging`, so the root logger needs a formatting
  handler from construction. It is built by hand rather than via `basicConfig(format=...)` precisely so
  the filter can be attached to it.
- **A declared stage and a defaulted one look different on purpose.** Declared stages are Title Case
  phrases naming a pipeline step ("Foldseek Alignment"); a defaulted one is a function name
  (`write_through_part`). The difference is the signal that nothing declared a stage there.
- **Passing no `extra` is now a legitimate choice**, taken where the function name is already the best
  label — the `"Initialized"`/`"Started"` debug lines in three constructors, and the optional
  `log_extra` parameters of `foldseek.run_foldseek` and the two `downloads.lib_download` entry points.

How the `extra` is built follows one rule, and there is exactly one spelling for it: **`log_extra`**.
A class with a single coherent stage sets `self.log_extra` once in `__init__` and never mutates it
(`StructureDownloader`, `PisaDownloader`, `StructureAligner`, `StructurePreprocessor`). Anything spanning
several stages builds a local `log_extra` per function, or passes the dict inline when the function has
only one call (`pocketmapper.py`, `pisa_parser`, `pocket_parser`, `pocket_comparison`, ...). `QTProcessor`
is the one deliberate `.update()`: `process_qt_cmdline_input` names the side being processed and the
`determine_*` helpers it drives all log under that name, which is call-scoped context rather than drift.
`PocketMapper` itself holds no logging state — it used to, and the stage a step logged under then depended
on which earlier step had last updated it.

All logging is module-level `logging.debug(...)`. No module instantiates a logger; `getLogger` appears
nowhere in the package.

Errors are `logging.critical(...)` then `raise PocketMapperError(...)`; `cli()` catches and exits 1.

That convention now holds on the Foldseek path too, which is most of what `foldseek.run_foldseek` buys:
five of the six invocations used to be a bare `subprocess.run(..., check=True)`, so a failing Foldseek
surfaced as a `CalledProcessError` traceback rather than a message. Exit code was 1 either way.

**No `exit()`/`sys.exit()` inside modules** — deliberately removed, which no code comment can show. There
are now none: the last survivor was `_check_help_search`, deleted along with `search()`'s `help` parameter
when argparse took over `--help`. The only `sys.exit` in the package is in `cli.py`'s `cli()`, which is
the boundary and is meant to have one.

## Settings

**Every `Settings` field is reachable from the command line**, and the settings JSON sets nothing the
CLI cannot. The file is a convenience for keeping a long invocation reproducible, never the only route
to a setting; the layering that makes it one is in `configure_workflow`.

A new option goes in **five hand-maintained places**: the `Settings` dataclass, the `search()` signature
together with the `cli_overrides` dict directly beneath it (one site — the dict mirrors the signature and
sits next to it precisely so the two cannot drift), the parser in `cli.py`, `cli()`'s kwarg block in the
same file, and the README's Options tables. None is generated from the dataclass. Miss one and the option
is silently ignored — `cli()` is the one that reads like boilerplate and is easiest to forget.

Nothing enforces the agreement, but it is checkable in a few lines: `dataclasses.fields(Settings)`,
`inspect.signature(PocketMapper.search)`, the `cli_overrides` keys, the subparser's `_actions` dests and
`cli()`'s `x=args.x` lines must all name the same fields (modulo `settings`, which is a parser-side
spelling rather than a field). The `query` and `target` positionals carry those dests, so they line up
with the rest.

Options are grouped by lifetime in both places a human reads them — argparse's argument groups in
`build_parser`, and the README's matching subsections. The path fields alone roughly double the
option count, so leaving them ungrouped would bury `--foldseek` and `--query_pocket_method` among
them. The four groups are `aligned structure options`, `cache options` (what survives a run),
`out options` (what the run produces) and `temp options` (what `delete_tmp` removes at the end);
`--cache_dir`, `--results_dir` and `--temp_dir` head the group whose defaults derive from them. Adding
a path setting means picking one of those groups in both places.

**`temp_dir` has exactly one option and three directories.** `query_structures/`, `target_structures/`
and `foldseek_tmp/` are computed in `configure_temp_dir` and held on the `PocketMapper` instance, not
declared on `Settings` — a `Settings` field has to be reachable from the command line (above), and
placing these individually is what `--temp_dir` replaced. Only the first two are created: Foldseek
makes its own. `temp_dir` is emptied there rather than merely created, so a rerun into the same
`results_dir` cannot hand `createdb` the previous run's structures; the emptying is guarded by
`lib.is_within` and the creation deliberately is not, since a run pointed outside both roots still
needs somewhere to work. That step sits after `configure_logging` for the same reason 4b/4c do —
the skip warning would otherwise be swallowed by the CRITICAL-only root logger.

`results_dir` is in `configure_workflow`'s `dirs_to_create` in its own right, and has to stay there.
Every other path in that list is settable away from `results_dir`, so without it `configure_logging`'s
file handler is one `--aligned_structure_dir` away from opening a log in a directory nothing made.

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
  `downloads/structure_downloader.py`) and the PEP 604 `str | None` field annotations on `Pocket`, `PocketResidue`,
  `QTRecord` and `Settings` all require it. No module carries `from __future__ import annotations`, so those
  annotations are evaluated at import rather than deferred. Rewriting all of that for 3.9 would still fail:
  biopython requires >=3.10.
- **`compat` is what guards the floor, not `lint`.** flake8 parses with whatever interpreter runs it, so lint
  at 3.12 cannot see a 3.12-only construct. `compileall` at 3.10 is what catches syntax; the import step is
  what catches the annotation and `importlib.resources` failures that compileall cannot. That step walks with
  `pkgutil.walk_packages` and a prefix, not `iter_modules`: `iter_modules` stops at the top level, so once
  `downloads` became a subpackage it would import that package's `__init__` and silently skip the three
  modules under it. It asserts a module count for the same reason.
- **3.10 is the only version pip resolves to pandas 2.x** — 3.11 and up get pandas 3.x. That is why the e2e
  matrix covers 3.10 and 3.14 rather than the middle. Both produce identical comparison row counts across
  every non-`huge` case.
- The bundled Foldseek DB is resolved through `files("pocketmapper")`, not through the data directory, for
  the reason given at that call site in `foldseek`.

## Repo layout

Each module's own docstring states its remit. Not stated anywhere in the code:

- There are no unit tests. `tests/e2e/` is the whole suite; the `pocketmapper-e2e` skill covers running it.
- **`fixtures/invalid_residues.txt`'s second line is deliberately wrong.** `4Q5J:A:9999` names a
  residue chain A does not have, and `test_invalid_1` expects the run to succeed anyway on the first
  line -- so "correcting" the 9999 silently removes the only case covering the skip.
- **`fixtures/forced_pisa_mixed.txt`'s second line is deliberately wrong too.** `4Q5J:A` names no partner
  chain, so it is the entry `--query_pocket_method pisa` must reject; the other two lines are what
  `test_invalid_8` asserts still produce rows. Give the middle line a partner and the case stops proving
  that a forced method skips entries rather than aborting the run.
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
  3.12-only f-string that `downloads/pisa_downloader.py` no longer does, so an import-everything check
  passes in CI and fails locally. A stale copy now also holds `structure_fetcher.py` and
  `pisa_downloader.py` at their old top-level paths, which shadow the `downloads` package versions and
  hide a missed import update. Delete `build/` before building or testing a wheel; never edit
  `build/lib/pocketmapper/`.
- **Structure parsing is gemmi throughout** (`.cif.gz` on disk). Biopython is used only for pairwise
  alignment (`sequence_aligner.py`) and SVD superposition (`pocket_comparison.py`).
- **The two cached-output classes no longer have the same shape.** `StructurePreprocessor` still has
  the required `set_output_directory()` -> `update_cache()` -> `preprocess_records()` order that nothing
  enforces, and caches on bare filenames; its class and `update_cache` docstrings say so.
  `StructureDownloader` has none of that — every record carries its own `struct_path`, so it tests that
  exact path and takes no output directory at all. Both still write through a `.part` file, but
  `StructureDownloader` gets that from `downloads.lib_download` while `StructurePreprocessor` keeps its
  own, since it writes a file it computed rather than one it fetched.
- **Nothing creates a structure's parent directory.** `StructureDownloader` dropped the `makedirs` that
  `set_output_directory` used to do, and `write_through_part` does not add one, so a missing directory
  is a non-transient failure that surfaces as `structure_not_found` rather than an error. On the
  pipeline path `configure_workflow` has already created `structure_dir`; a library caller driving the
  component directly has to create it.

## As a library

The CLI is confined to `cli.py`, so nothing else here needs a terminal. Two levels of
entry: `PocketMapper().search(...)` does the same work as the CLI, or drive a component directly —
`qt_processor`, `downloads.structure_downloader`, `structure_preprocessor`, `downloads.pisa_downloader`,
`pisa_parser`, `sequence_aligner`, `structure_aligner`, `pocket_calculator`, `foldseek` are each
separately usable.

- **Step 7 is the one step that can be deferred.** `search(align_count=0)` writes everything but the
  aligned structures, and `StructureAligner.align_structs` then produces them from the run's own outputs
  — `pm.query_df` / `pm.target_df` as records, plus the two result paths off `pm.settings`. Verified: the
  PDBs come out byte-identical to those of a normal run, on both the structure and the Foldseek-DB path.
  It works because records point at `structure_dir`, which `delete_tmp` never touches, and it is why
  `align_structs` takes `query_ids`, `target_ids` and `overwrite` — a deferred caller superposes a few
  queries at a time rather than all of them.
- **A component reaching into a `Settings` can't be used without building one, and hides which fields it
  depends on** — so no component takes one. The `Settings` is unpacked at each call site in
  `pocketmapper.py` into the values that component needs.
- **`pocketmapper/__init__.py` only exports `PocketMapper` and `__version__`.** The console script
  loads `cli` from `pocketmapper.cli` directly, so the package does not re-export it. Submodules are
  reachable as `pocketmapper.lib` etc. only as a side effect of its importing
  `pocketmapper.pocketmapper` — for anything else, always use explicit
  `from pocketmapper.<module> import <name>`.
- **Always call `Settings(...).resolve_paths()`** if you build one yourself — the failure mode is in that
  method's docstring. `search()` does this for you.
- **`search()` has global side effects**: `logging.config.dictConfig` reconfigures the *root* logger and
  stomps on a host app's logging setup, and `delete_tmp` `shutil.rmtree`s `temp_dir` at the end unless
  `delete_tmp=False`, which keeps it. `configure_workflow` also empties `temp_dir` on the way in. Both
  rmtrees are guarded rather than unconditional: `temp_dir` is settable, so `lib.is_within` skips (with
  a warning) a path that does not resolve under `cache_dir` or `results_dir`. The guard bounds the
  damage from a mistyped path; it is not a reason to point the setting at a directory you care about.
- Results come back through files — `search()` returns `None`, so read `pocket_comparison.tsv` /
  `alignment.tsv` from `results_dir` (paths available on `Settings`).
