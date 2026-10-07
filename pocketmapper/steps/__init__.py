"""
The pipeline's steps, one module each, in the order `search` runs them: `parse`, `fetch_structures`,
`align`, `pockets`, `compare` and `superpose`.

Each module has two layers. Its entry function, named after the step, is what the CLI and `search`
call: it layers its settings, resolves its path defaults under `results_dir`, logs to `log_path` for
the length of the call, validates its settings, manages its scratch directory and calls the core
function. The core function takes explicit paths and values, never a `Settings`, and keeps no state
between calls. The files the steps hand each other are described in `pocketmapper.records`.

Settings: an argument given beats the job file (a JSON path, or the same already loaded as a dict),
which beats `settings.SETTING_DEFAULTS`. A job file keyed by `Settings` field names, such as a run's
job_settings.json, works as-is; it names every path, so an argument moves only the path it names.
`results_dir` is required, as an argument or in the job file, by every step but `parse`.

Defaults: an input defaults to the matching `*_path` setting (the standard file under `results_dir`
unless the job file names another), and an output that rewrites an input (`fetch_structures` and
`align` records) defaults to that input's path.

The cache is chosen once, by `parse`, which writes the cache manifest into `results_dir`.
`fetch_structures`, `align` and `pockets` read it from there and ignore any cache setting in a job
file, so a chain cannot split its cache.
"""
