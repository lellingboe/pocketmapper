# Settings

How settings, job keys and options are layered, resolved and dumped, and where a new option goes. Moved
out of `context/project-overview.md` so it loads only when working on settings or options.

`Settings` (`lib/settings.py`) is built only by search. **Every `Settings` field is reachable from the CLI**;
the job file sets nothing the CLI can't. Layering is `lib.settings.layer_settings` for every command
(arguments not None over job file over `SETTING_DEFAULTS`), then `resolve_paths` for search and parse or
`resolve_step_paths` for every other step; search's resolution is
in `resolve_settings` (`resolve_*`). `Settings` has no defaults and is built once. Pocket methods keep
`"auto"` (inferred per entry); `align_struct_method` is resolved per run. The steps take a job dict (or
file) keyed by job keys, never a `Settings`, and call the same `lib.settings.resolve_*` validators.

**Job keys are separate from `Settings`.** `lib.settings.JOB_KEYS` = `fields(Settings)` + `STEP_ONLY_KEYS`
(options only some steps take; `STEP_ONLY_PATH_KEYS` among them resolve against `work_dir` in
`resolve_paths`) is what a job file may hold and what `layer_settings` seeds. Unknown
keys are rejected; one in `REMOVED_JOB_KEYS` gets its replacement named. Search drops keys `Settings`
lacks before building one.

**Only search and parse take `cache_dir`, `results_dir` and `work_dir`.** `resolve_paths` derives every
path from them and stores it absolute (`work_dir`, default cwd, is what relative ones join onto). Every
other step takes each path it uses from its job file or arguments (`resolve_step_paths`: required ones
checked, a relative one made absolute against the cwd). A step's input argument (`--entries`,
`--alignment`, ...) is layered over the matching `*_path` key, so `entries` sets `entries_path`. A step
run from a dump elsewhere resolves as the first did, since the dump holds absolute paths.

**Settings dumps.** Only search and parse write one. `lib.settings.dump_settings` writes the values,
minus `job_settings_path`, to `job_settings_path`: search's resolved `Settings` to `job_settings.json`,
parse's layered values (paths resolved) to `<results_dir>/parse_settings.json` (the default comes from
`resolve_paths(values, command)`, not `RESULTS_PATH_DEFAULTS`). `layer_settings` drops
`job_settings_path` from a job file, so a command chained off another's dump never overwrites it. parse
nested in search writes none: `lib.scopes.run_scope` records the outermost command and nested scopes yield
False (`lib.scopes.outer_command()` reads it).

**`results_dir` is optional for parse**: `resolve_paths` then leaves every results path not given None
(`temp_dir` → `<cache_dir>/tmp`), and `log_to_file(None)`, `dump_settings` and `lib.entries.write_entries`
do nothing for a None path. Only search defaults it to a timestamped name.

**A reused dump names every path**, so an argument moves only the path it names.

A new option goes in **four hand-maintained places per command**: the function signature and its
`layer_settings` arguments dict (for search also `Settings`; for a step-only option
`lib.settings.STEP_ONLY_KEYS` instead; for a step's path, its `resolve_step_paths` call), `OPTIONS` plus the command's `COMMANDS` row
in `cli.py`, and README Options (for a step, its row in the step table). A fifth for an option with a
static default: `lib.settings.SETTING_DEFAULTS` (from `lib.constants.DEFAULT_*`, which the help also quotes).
Signature and parser defaults are all None. Path defaults are `lib.settings.CACHE_PATH_DEFAULTS` /
`RESULTS_PATH_DEFAULTS` (`cache_path`, `results_path`). `cli()` passes parsed args straight through its
dispatch table, so there is no kwarg block. Miss one → argparse fails or the option is silently
ignored. Check, per command: `inspect.signature(<step>.step.<command>)` == subparser `_actions` dests; for
search also `dataclasses.fields(Settings)` == signature == `arguments` keys (modulo `job_file`);
`SETTING_DEFAULTS` keys ⊆ `JOB_KEYS`.

Options are grouped by lifetime in `cli.COMMANDS` and README alike: `in`, `aligned structure`, `cache`,
`out`, `temp`, `advanced` (README lists advanced options in the main table: an `Advanced options`
heading there would take the `#advanced-options` anchor). A file a command reads has a bare name (`--alignment`); one it writes ends in
`_path`. A new path setting picks a group in both.

**Every directory is made by its first writer**; there is no bulk creation step. Only the log's
directory is made up front, by `lib.logs.log_to_file` (the file handler opens at once). A writer that does
not make its directory fails once a path option points elsewhere.

**`temp_dir` is one option, three dirs**, `lib.constants.TEMP_SUBDIRS` (`query_structures/`,
`target_structures/`, `foldseek_tmp/`), used in `align/step.py` (all three) and `fetch/step.py`
(`foldseek_tmp/`, for a DB download). The two
structure dirs are made in `foldseek_preprocessing`, `foldseek_tmp/` before each Foldseek call that uses
it (MMseqs2 makes only one level of tmp dir), so a `seq` run makes none. Ownership is
`lib.scopes.temp_dir_scope` (search, fetch_structures, align): emptied on entry so reruns can't feed `createdb`
stale structures, deleted on a clean exit unless `delete_tmp` is 0. Scopes nest: search's owns temp
for the whole run and the steps' are no-ops inside it, so `search --delete_tmp 0` keeps every step's
scratch. The held dir is released in a `finally`, else a failed library run leaves the next one in that
process neither emptying nor deleting it. Only the `TEMP_SUBDIRS` are emptied or deleted, and `temp_dir`
itself only once that leaves it empty, so no root guard is needed and a step needs no `cache_dir` or
`results_dir`. `lib.scopes.delete_temp_dir` does nothing if it was never made.

`search --help` is generated from `help=` strings; `CLI_SEARCH_EPILOG` holds only examples. Resolution
order, the `aligner` check and tri-state `align_struct_method`: `Settings` docstring and `# 4b.`/`# 4c.` in
`resolve_settings`. `aligner` has no auto mode; a broken `foldseek` fails at 4b before fetching.
