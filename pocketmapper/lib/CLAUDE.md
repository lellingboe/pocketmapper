# lib

Cross-module facts about `lib/`, moved out of `context/project-overview.md` so they load only when
working here.

## Downloads

`lib.download` offers **two entry points, not one**, because the two kinds of thing it fetches differ
in every operational respect. `download_file` is for bulk static files: unpaced and stateless.
`download_api` is for REST endpoints: it paces its requests and, on a transient failure, doubles that
host's pacing and never lowers it again. Both write through a `.part` file named uniquely per attempt by
`lib.files.part_path`, so two runs sharing a cache never write the same scratch file. Both also create the
destination's directory, share one retry test, and fetch through `fetch_to` (a User-Agent, and a
`DOWNLOAD_TIMEOUT` of 60 s per socket operation). A host that stalls every time therefore costs
`max_retries` x 60 s per file, and for `download_api` those files run one after another.

**`download_many` is the one download pool**, used by `lib.structure_manager.StructureManager.fetch`
and `PisaDownloader.download_missing_assemblies`. It skips a `DownloadJob` whose file already exists, so
the file is the cache; `overwrite=True` fetches it anyway (a structure refresh), and a failed fetch then
leaves the cached file in place. **The file existing is still the only cache test, but only checked content
reaches the final name**: a job's content check (`check_gzip` for a PDB structure, `check_mmcif` for an
AlphaFold model, `check_json` for a PISA file) runs on the response before it is placed. A 200 response
that fails it (a proxy or maintenance page) raises `UnexpectedContentError`, which is not transient, so it
costs one request, is not cached, and is asked for again on the next run. A file cached before the checks
existed is not re-checked. **With `fetch=download_api` it runs one job at a time, whatever `max_workers`
says**: `download_api` sleeps the host delay outside the lock, so concurrent calls would each sleep the
same delay and then fire together. The rule sits beside the pacer, so no caller repeats it.
**Ctrl-C cancels the queued jobs, not the running ones**: `download_many` waits for the at most
`max_workers` jobs in flight, each of which still finishes its own retries, then re-raises.

Four consequences no single file states:

- **The download pools ignore `--threads`, on purpose.** `lib.structure_manager.fetch_missing_structures`
  calls `StructureManager.fetch`, and `pockets.pocket_manager.PocketManager` builds one `PisaDownloader` and
  one `StructureManager`, each with the default width, `DOWNLOAD_WORKERS` (8), so `--threads` governs
  Foldseek alone. A worker here waits on a socket rather than on a core, so a thread count says nothing
  useful about how wide a pool should be. A library caller can still pass its own `max_workers`;
  `PocketManager(max_workers=...)` sets the width of both its pools.
- **The pacing delay and the backoff delay are the same number.** That is what "carry the backoff
  forward" means here: the escalation one retry needed becomes the pace of every later request to that
  host. It is also why `download_missing_summaries` and `download_missing_assemblies` no longer sleep
  themselves — the helper owns pacing, and a caller-side sleep would double it. A `Retry-After` (capped at
  `max_delay`) lengthens a single wait without changing the host's pace; `download_file` also jitters each
  sleep by 0.8 to 1.2 so concurrent failures do not retry in step.
- **The delay registry is module-level and never decays**, so it outlives any one `PisaDownloader` —
  which matters, because each `PocketManager` builds its own, one per pockets run.
  It equally outlives a whole `search()`, so a library caller running several in one process carries an
  elevated delay across all of them; `reset_host_delays` is the escape hatch.
- **Only 5xx and 408/425/429 are retried.** A whitelist among 4xx rather than a blacklist of 404, so an
  unrecognised 4xx costs one request instead of the whole budget. The old PISA code caught bare
  `Exception`, so every entry PISA lacked cost 5 requests and ~3.75s of sleeping — on the full-PDB path
  that is thousands of entries.

A leftover scratch file (`<name>.<hex>.part`) is inert in every cache directory:
`download_missing_interfaces` offline globs `*.json`, which cannot match `x.json.<hex>.part`, and everything
else checks an exact path — `download_many` its job's `out_fpath` (a structure's `.cif.gz`, a PISA assembly's
`.json`), the online `is_complete` and the other PISA stages their own files, the align step's preprocessing its
`<preprocess_name>.cif.gz`.
The same holds for a leftover `summaries/_batch.json.<hex>.part`, which nothing reads. Being unique, they
are never overwritten either: a run killed mid-write leaves one per file it was writing, until deleted.
