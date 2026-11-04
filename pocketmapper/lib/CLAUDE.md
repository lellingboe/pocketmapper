# lib

Cross-module facts about `lib/`.

## Downloads

The `lib.download` module docstring covers the two entry points (`download_file` unpaced, `download_api`
paced), `.part` files, content checks, retries, timeouts and `download_many`. What follows from them:

- **`download_many` has two callers**: `StructureManager.fetch` and `PisaDownloader.download_missing_assemblies`.
  PISA summaries call `download_api` directly.
- **A host that stalls every time costs `max_retries` x 60 s per file**, and under `download_api` those files
  run one after another. Ctrl-C waits for the jobs in flight, each through its full retry loop.
- **The download pools ignore `--threads`, on purpose.** No step passes it into `max_workers`, so every pool is
  `DOWNLOAD_WORKERS` wide and `--threads` governs Foldseek alone. A library caller can pass `max_workers`.
- **Callers never sleep between requests.** `download_api` owns pacing, and its pacing delay is the backoff
  delay, so a caller-side sleep would double it.
- **The delay registry never decays and outlives a run.** Each pockets run builds its own `PocketManager` and
  `PisaDownloader`, and nothing in `search()` resets it, so a library caller running several searches in one
  process carries an elevated delay across all of them; `reset_host_delays` is the escape hatch.

**A leftover scratch file (`<name>.<hex>.part`) is inert in every cache directory.**
`download_missing_interfaces` offline globs `*.json`, which cannot match `x.json.<hex>.part`, and everything
else checks an exact path: `download_many` its job's `out_fpath`, `StructureManager.on_disk`, the online
`is_complete` and the other PISA stages their own files, the align step's preprocessing its
`<preprocess_name>.cif.gz`. The same holds for a leftover `summaries/_batch.json.<hex>.part`, which nothing
reads. Being unique, they are never overwritten either: a run killed mid-write leaves one per file it was
writing, until deleted.
