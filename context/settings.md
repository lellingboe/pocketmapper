# Settings

How settings, job keys and options are layered, resolved and dumped, and where a new option goes. Moved
out of `context/project-overview.md` so it loads only when working on settings or options.

Settings travel as a dict, job key -> value; there is no settings record. **Every job key is reachable from
the CLI**; the job file sets nothing the CLI can't. Every entry function, search's included, begins with
`values = lib.settings.layer_command("<cmd>", job_file, locals())`: it raises unless the parameters are
exactly `command_options(cmd)` (`job_file` aside), sets each under its job key and layers them
(`lib.settings.layer_settings`: arguments not None over job file over `SETTING_DEFAULTS`). Then
`resolve_paths` for search and parse or `resolve_step_paths` for every other step. The next step is
`lib.settings.validate_options(values, cmd)` (each option's `validate`, in registry order), inside the
command's log scope so a failure reaches the log, before the entries table is opened, anything is
downloaded or any early return; for search it is in `PocketMapper.resolve_settings`. Then the cross checks
(`check_offline_refresh`, `check_fsdb_*`). Pocket methods keep `"auto"` (inferred per entry);
`align_struct_method` is resolved per run (`resolve_auto_align_struct_method`, after validation). The
single-option validators are `lib.options.OPTIONS[name].validate`.

**Job keys.** `lib.options.JOB_KEYS`, derived from the registry `lib.options.OPTIONS` in its order, is
search's keys (`command_keys("search")`) plus `STEP_ONLY_KEYS` (options only some steps take;
`STEP_ONLY_PATH_KEYS` among them resolve against the cwd in `resolve_paths`). It is what a job file may
hold and what `layer_settings` seeds. Unknown keys are rejected; one in `REMOVED_JOB_KEYS` gets its
replacement named. Search keeps only `command_keys("search")` once its settings are checked.

**Only search and parse take `cache_dir` and `results_dir`.** `resolve_paths` derives every
path from them and stores it absolute against the cwd. `local_root` (default cwd, absolute) is the base
only for relative local structure files and user Foldseek database paths; every other path, `local_root`
included, resolves against the cwd. Every other step takes each path it uses from its job file or arguments
(`resolve_step_paths`: required ones checked, a relative one made absolute against the cwd, an unset
optional `local_root` the cwd). A step's input argument (`--entries`,
`--alignment`, ...) is layered over the matching `*_path` key, so `entries` sets `entries_path`. A step
run from a dump elsewhere resolves as the first did, since the dump holds absolute paths.

**Settings dumps.** Only search and parse write one. `lib.settings.dump_settings` writes the values,
minus `job_settings_path`, to `job_settings_path`: search's resolved `command_keys("search")`, in registry
order, to `job_settings.json`,
parse's layered values (paths resolved) to `<results_dir>/parse_settings.json` (the default comes from
`resolve_paths(values, command)`, not `RESULTS_PATH_DEFAULTS`). `layer_settings` drops
`job_settings_path` from a job file, so a command chained off another's dump never overwrites it. parse
nested in search writes none: `lib.scopes.run_scope` records the outermost command and nested scopes yield
False (`lib.scopes.outer_command()` reads it).

**`results_dir` is optional for parse**: `resolve_paths` then leaves every results path not given None
(`temp_dir` → `<cache_dir>/tmp`); `dump_settings` and `lib.entries.write_entries` do nothing for a None
path, and `log_to_file(None)` adds no file handler but still sets the level. Only search defaults it to
a timestamped name.

**A reused dump names every path**, so an argument moves only the path it names.

A new option goes in **five hand-maintained places per command**: the function signature and its
`layer_settings` arguments dict (search has none: `layer_command` checks its signature; for a step's path,
its `resolve_step_paths` call), `OPTIONS` plus the command's `COMMANDS` row in `cli.py`, README Options (for a step, its row in the
step table), and its `Option` in `lib.options.OPTIONS`, from which `SETTING_DEFAULTS` (its `default`, from
`lib.constants.DEFAULT_*`, which the help also quotes), `JOB_KEYS`, `STEP_ONLY_KEYS` and the path defaults
are derived. Signature and parser defaults are all None. Path defaults are `lib.options.CACHE_PATH_DEFAULTS`
/ `RESULTS_PATH_DEFAULTS` (an `Option`'s `under`; `cache_path`, `results_path`). `cli()` passes parsed args
straight through its dispatch table, so there is no kwarg block. Miss one → argparse fails or the option is
silently ignored. `tests/unit/test_cli.py` checks, per command, that the options are exactly the parameters
of its `cli.DISPATCH` function and all default to None; for search also signature ==
`command_options("search")` (modulo `job_file`); and `SETTING_DEFAULTS` keys ⊆ `JOB_KEYS`. No test
checks every step's `arguments` dict, or README.

Options are grouped by lifetime in `cli.COMMANDS` and README alike: `in`, `aligned structure`, `cache`,
`out`, `temp`, `advanced` (README lists advanced options in the main table: an `Advanced options`
heading there would take the `#advanced-options` anchor). A file a command reads has a bare name (`--alignment`); one it writes ends in
`_path`. A new path setting picks a group in both.

**Every directory is made by its first writer**; there is no bulk creation step. Only the log's
directory is made up front, by `lib.logs.log_to_file` (the file handler opens at once). A writer that does
not make its directory fails once a path option points elsewhere.

**`temp_dir` is one option, three dirs**, `lib.constants.TEMP_SUBDIRS` (`query_structures/`,
`target_structures/`, `foldseek_tmp/`), used in `align/align.py` (all three) and `fetch/fetch.py`
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

**`offline`, `fetch_missing` and `refresh_structures`.** `offline` 1 stops every download in each step
that takes it: structures, PISA interfaces and Foldseek databases. It wins over `fetch_missing`, so a
step downloads a missing entry only with `fetch_missing and not offline` (align, pockets and superpose
skip the fetch call; fetch_structures builds `StructureManager(download=False)`). An uncached structure
is marked `structure_not_found`, as a failed download is; an uncached Foldseek database is an error
(fetch_structures: `lib.structure_manager.fetch_missing_fsdb` raises with a `StructureManager(download=False)`;
align: no usable target).
pockets passes `download=not offline` to `PocketManager`, for PISA and the Foldseek hits' structures.
`refresh_structures` 1 re-downloads only PDB and AlphaFold structures, and only in the fetch step
(`StructureManager(refresh=True)`); PISA interfaces and Foldseek databases are never refreshed. Both 1 is
an error, checked before anything is fetched or dumped by one helper, `lib.settings.check_offline_refresh`,
from search's `resolve_settings` and `fetch_structures`; `StructureManager.__init__` rejects the same
pair for a library caller.

`search --help` is generated from `help=` strings; `CLI_SEARCH_EPILOG` holds only examples. Resolution
order, the `aligner` check and tri-state `align_struct_method` (`lib.settings.resolve_auto_align_struct_method`):
the comments in search's `resolve_settings`. `aligner` has no auto mode; a broken `foldseek` fails there, after
every other setting is checked and before fetching.
