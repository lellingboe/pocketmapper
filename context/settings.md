# Settings

How settings, job keys and options are layered, resolved and dumped, and where a new option goes. Moved
out of `context/project-overview.md` so it loads only when working on settings or options.

`Settings` (`settings.py`) is built only by search. **Every `Settings` field is reachable from the CLI**;
the job file sets nothing the CLI can't. Layering is `settings.layer_settings` for every command
(arguments not None over job file over `SETTING_DEFAULTS`, then `resolve_paths`); search's resolution is
in `resolve_settings` (`resolve_*`). `Settings` has no defaults and is built once. Pocket methods keep
`"auto"` (inferred per entry); `align_struct_method` is resolved per run. The steps take a job dict (or
file) keyed by job keys, never a `Settings`, and call the same `settings.resolve_*` validators.

**Job keys are separate from `Settings`.** `settings.JOB_KEYS` = `fields(Settings)` + `STEP_ONLY_KEYS`
(options only some steps take; `STEP_ONLY_PATH_KEYS` among them resolve against `work_dir` in
`resolve_paths`) is what a job file may hold and what `layer_settings` seeds. Unknown
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

**`results_dir` is optional for parse and fetch_structures**: `resolve_paths` then leaves every
results path not given None (`temp_dir` → `<cache_dir>/tmp`), and `log_to_file(None)`,
`dump_settings` and the failed-entries writers do nothing for a None path. Only search defaults it
to a timestamped name; the other steps `require_setting` it.

**A reused dump names every path**, so an argument moves only the path it names. `results_dir` is not
inert beside it all the same: it is a `temp_dir` emptying root.

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
`out`, `temp`, `advanced` (README lists advanced options in the main table: an `Advanced options`
heading there would take the `#advanced-options` anchor). A file a command reads has a bare name (`--alignment`); one it writes ends in
`_path`. A new path setting picks a group in both.

**Every directory is made by its first writer**; there is no bulk creation step. Only the log's
directory is made up front, by `lib.log_to_file` (the file handler opens at once). A writer that does
not make its directory fails once a path option points elsewhere.

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
