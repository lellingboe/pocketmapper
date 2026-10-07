"""
The pipeline's steps, one module each, in the order `search` runs them: `parse`, `fetch_structures`,
`align`, `pockets`, `compare` and `superpose`.

Each module has two layers. Its entry function, named after the step, is what the CLI and `search`
call: it layers its settings, resolves its path defaults under `results_dir`, logs to `log_path` for
the length of the call, validates its settings, manages its scratch directory and calls the core
function. The core function takes explicit paths and values, never a `Settings`, and keeps no state
between calls. Core functions take records in memory, never a records file.

Settings: an argument given beats the job file (a JSON path, or the same already loaded as a dict),
which beats `settings.SETTING_DEFAULTS`. A job file keyed by job keys, such as a run's
job_settings.json, works as-is; it names every path, so an argument moves only the path it names.
`results_dir` is required, as an argument or in the job file, by every step but `parse`. Run on its
own, each step writes its settings to `<results_dir>/<command>_settings.json`, a job file for the
next; inside `search` it writes none.

Entries: every step parses the query and target entries its settings name, through
`steps.parse.parse_job_entries`, and rejects what `parse` rejects. No step reads another's records;
the cache and the files under `results_dir` are all the steps hand on. Every step takes the cache
options, so a chain shares its cache by taking them from the same job file.

Defaults: an input defaults to the matching `*_path` setting (the standard file under `results_dir`
unless the job file names another).
"""
