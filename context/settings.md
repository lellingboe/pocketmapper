# Settings

Cross-module rules for options, job keys, path settings and their layering.

## Where it is stated

- Layering, path resolution and validation: `lib.settings` module docstring, `layer_command`,
  `resolve_paths`, `resolve_step_paths`, `validate_options`.
- Job keys and what is derived from the registry: `lib.options` module docstring, `JOB_KEYS`.
- Dumps: `lib.settings.dump_settings`; parse's in `parse.parse`; nested parse writes none
  (`lib.scopes.run_scope`).
- Temp dirs: `lib.scopes` (`temp_dir_scope`, `delete_temp_dir`), `lib.constants.TEMP_SUBDIRS`.
- Search's resolution order, the `aligner` check and `align_struct_method` "auto": the comments in
  `PocketMapper.resolve_settings`. `--help` comes from `help=` strings; `lib.constants.CLI_COMMAND_EPILOGS`
  holds only examples.

## Rules

- **Every job key is reachable from the CLI**; a job file sets nothing the CLI can't.
- **Validate first.** Each entry function validates right after layering and resolving its paths, inside its
  log scope, before it opens the entries table, downloads anything or returns early
  (`test_step_validates_before_reading_entries`). The cross checks (`check_offline_refresh`,
  `check_fsdb_*`) follow.
- **Every directory is made by its first writer**; there is no bulk creation step. Only the log's directory
  is made up front (`lib.logs.log_to_file`). A writer that does not make its directory fails once a path
  option points elsewhere.
- **Temp.** `lib.scopes.temp_dir_scope` is owned by search, fetch_structures and align; nested scopes are
  no-ops, so `search --delete_tmp 0` keeps every step's scratch. align uses all three `TEMP_SUBDIRS`,
  fetch only `foldseek_tmp/` (for a DB download). The structure dirs are made in
  `foldseek_preprocessing`, `foldseek_tmp/` before each Foldseek call that uses it, so a `seq` run makes
  none.
- **`offline`, `fetch_missing`, `refresh_structures`.** align, pockets and superpose each fetch a missing
  entry only with `fetch_missing and not offline`; fetch_structures builds `StructureManager(download=False)`
  and pockets passes `download=not offline` to `PocketManager`. Offline, an uncached Foldseek database is
  an error. Both `offline` and `refresh_structures` 1 is rejected by one helper,
  `lib.settings.check_offline_refresh` (search's `resolve_settings` and `fetch_structures`), and by
  `StructureManager.__init__` for a library caller.

## A new option needs

1. One `Option` in `lib.options.OPTIONS`, with its default from `lib.constants.DEFAULT_*` (which the help
   also quotes). Everything else is derived from it; parser defaults are all None.
2. The parameter and its docstring `Args:` entry in each entry function that takes it. `layer_command`
   raises unless the parameters are exactly `command_options(cmd)` (`tests/unit/test_cli.py` checks this per
   command); nothing checks the `Args:` entry.
3. The code that reads it, and for a step's path its `resolve_step_paths` call. Nothing checks this: an
   option that is parsed and validated but never read is silently ignored.
4. The README: its table row (for a step, its entry in the step table's Options column, or a group link),
   in the same group as in `lib.options.GROUPS`. `tests/unit/test_readme.py` compares only flags with
   `command_options`, plus `SHARED_STEP_FLAGS`, because the README wording is richer than `--help`.

README lists advanced options in the main table: an `Advanced options` heading there would take the
`#advanced-options` anchor. A file a command reads has a bare name (`--alignment`); one it writes ends in
`_path`.
