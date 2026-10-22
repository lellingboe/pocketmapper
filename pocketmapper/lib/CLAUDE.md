# lib

Cross-module facts about `lib/`, moved out of `context/project-overview.md` so they load only when
working here.

## Downloads

`lib.download` offers **two entry points, not one**, because the two kinds of thing it fetches differ
in every operational respect. `download_file` is for bulk static files: unpaced and stateless.
`download_api` is for REST endpoints: it paces its requests and, on a transient failure, doubles that
host's pacing and never lowers it again. Both write through a `.part` file, create the destination's
directory, and share one retry test.

**`download_many` is the one download pool**, used by `fetch.structure_downloader.download_structures`
and `PisaDownloader.download_missing_assemblies`. It skips a `DownloadJob` whose file already exists, so
the file is the cache. **With `fetch=download_api` it runs one job at a time, whatever `max_workers`
says**: `download_api` sleeps the host delay outside the lock, so concurrent calls would each sleep the
same delay and then fire together. The rule sits beside the pacer, so no caller repeats it.

Four consequences no single file states:

- **The download pools ignore `--threads`, on purpose.** `fetch.step.fetch_missing_structures` calls
  `download_structures`, and `pockets.pisa.download_pisa_interfaces` builds `PisaDownloader`, each with
  the default width, `DOWNLOAD_WORKERS` (8), so `--threads` governs Foldseek alone. A worker here waits
  on a socket rather than on a core, so a thread count says nothing useful about how wide a pool should
  be. A library caller can still pass its own `max_workers`.
- **The pacing delay and the backoff delay are the same number.** That is what "carry the backoff
  forward" means here: the escalation one retry needed becomes the pace of every later request to that
  host. It is also why `download_missing_summaries` and `download_missing_assemblies` no longer sleep
  themselves — the helper owns pacing, and a caller-side sleep would double it.
- **The delay registry is module-level and never decays**, so it outlives any one `PisaDownloader` —
  which matters, because `pisa_pockets` and `expand_fsdb_pdb_targets` each build a fresh one per run.
  It equally outlives a whole `search()`, so a library caller running several in one process carries an
  elevated delay across all of them; `reset_host_delays` is the escape hatch.
- **Only 5xx and 408/425/429 are retried.** A whitelist among 4xx rather than a blacklist of 404, so an
  unrecognised 4xx costs one request instead of the whole budget. The old PISA code caught bare
  `Exception`, so every entry PISA lacked cost 5 requests and ~3.75s of sleeping — on the full-PDB path
  that is thousands of entries.

A leftover `.part` or `.raw.part` is inert in every cache directory: `download_missing_interfaces`
globs `*.json`, which cannot match `x.json.part`, and everything else checks an exact path —
`download_many` its job's `out_fpath` (a structure's `.cif.gz`, a PISA assembly's `.json`), the other
PISA stages their own files. The same holds for a leftover `summaries/_batch.json`, since no PDB code
starts with `_`.
