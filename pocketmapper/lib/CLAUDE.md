# lib

Cross-module facts about `lib/`, moved out of `context/project-overview.md` so they load only when
working here.

## Downloads

`lib.download` offers **two entry points, not one**, because the two kinds of thing it fetches differ
in every operational respect. `download_file` is for bulk static files: unpaced, stateless, and
called from the pools of `download_missing_structures` and `download_missing_ftp_assemblies`.
`download_api` is for REST endpoints: it paces its requests and, on a transient failure, doubles that
host's pacing and never lowers it again. Both write through a `.part` file and share one retry test.

Four consequences no single file states:

- **The download pools ignore `--threads`, on purpose.** `fetch.step.fetch_missing_structures` and
  `pockets.pisa.download_pisa_interfaces` build `StructureDownloader` and `PisaDownloader` with their
  default width, `DOWNLOAD_WORKERS` (8), so `--threads` governs Foldseek alone. A worker here waits on
  a socket rather than on a core, so a thread count says nothing useful about how wide a pool should be.
  A library caller can still pass its own `max_workers`.
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

A leftover `.part` is inert in every cache directory: `download_missing_interfaces` globs `*.json`,
which cannot match `x.json.part`, the other PISA stages check an exact path, and
`StructureDownloader` tests for an exact `.cif.gz` destination. The same holds for a leftover
`summaries/_batch.json`, since no PDB code starts with `_`.
