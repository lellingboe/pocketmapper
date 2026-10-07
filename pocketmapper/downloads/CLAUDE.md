# Downloads

Cross-module facts about `downloads/`, moved out of `context/project-overview.md` so they load only
when working here.

`downloads/` holds what fetches bytes over HTTP, and nothing else. `foldseek.py` and
`PocketMapper.fetch_missing_fsdb` stay outside it: they download a database by shelling out to the
`foldseek` binary, so they share none of the machinery here, and moving them would split `foldseek.py`
away from `check_foldseek` and `bundled_foldseek_dbs`, which `qt_processor` imports for reasons
unrelated to downloading.

`lib_download` offers **two entry points, not one**, because the two things the package fetches differ
in every operational respect. `download_file` is for bulk structure files: unpaced, stateless, and
called from `download_missing_structures`' pool. `download_api` is for REST endpoints: it
paces its requests and, on a transient failure, doubles that host's pacing and never lowers it
again. Both write through a `.part` file and share one retry test.

Four consequences no single file states:

- **The download pool ignores `--threads`, on purpose.** `steps.fetch.fetch_missing_structures` builds
  `StructureDownloader` with its default width, `DOWNLOAD_WORKERS` (8), so `--threads` governs Foldseek
  alone. A worker here waits on a socket rather than on a core, so a thread count says nothing useful
  about how wide the pool should be. A library caller can still pass its own `max_workers`.
- **The pacing delay and the backoff delay are the same number.** That is what "carry the backoff
  forward" means here: the escalation one retry needed becomes the pace of every later request to that
  host. It is also why `download_missing_summaries` and `download_missing_assemblies` no longer sleep
  themselves — the helper owns pacing, and a caller-side sleep would double it.
- **The delay registry is module-level and never decays**, so it outlives any one `PisaDownloader` —
  which matters, because the fetch step, `pisa_pockets` and `expand_fsdb_pdb_targets` each build a fresh
  one per run.
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

**PISA assembly interfaces have two sources, one cache.** `source="ftp"` (the default) fetches
`ftp.ebi.ac.uk/.../pdb-assemblies-analysis/split/<code[1:3]>/<code>_assembly<asm>_interfaces.json`
through `download_file` in a `DOWNLOAD_WORKERS` pool, since those are static files; `source="api"` keeps
the paced `download_api` loop. Summaries, which give the assembly ids, come from the API
either way. Both write `assemblies/<code>_<asm>.json`, because the two serve the same JSON (checked on
2026-09-29: 6rkw assembly 1 byte-for-byte, 4q5j and 12 4dx9 assemblies equal as parsed JSON). So a
cached file is source-agnostic, and with a warm cache the toggle fetches nothing. **If the formats ever
diverge, the caches must split.** An assembly the FTP lacks is a 404, not retried, and is reported
under `assembly_downloading` like an API failure; there is no fallback to the API.

**Entry summaries are the one batched request.** `download_missing_summaries` POSTs
`summary_batch_size` ids at a time to `/pdb/entry/summary/` and splits the response into the same
per-entry files a single GET would give, so ~3,620 requests become ~73. The PISA API interface endpoint
takes one assembly per call and has no batched form. Three behaviours of the API, measured on
2026-09-28, shape the code:

- **An unknown id is silently omitted**, not reported, so a code missing from its batch's response is
  counted as a failure. A batch with *no* known id returns 404, which fails the whole batch.
- **Somewhere between 800 and 900 ids the API returns 500.** That is transient to `is_transient_error`,
  so an oversized batch would be retried five times and permanently double the host's pacing. The
  default of 50 stays well clear.
- The summary is the same whether fetched singly or in a batch (verified over 450 cached entries),
  so switching left `parse_summaries` and the existing cache untouched.

**The PISA failure report is the caller's file, not the downloader's.** `download_missing_interfaces`
returns what each stage could not handle and writes `error_path` only when there is something to write;
the fetch step, `pisa_pockets` (unless given `download=False`) and `expand_fsdb_pdb_targets` all reach
it through `download_pisa_interfaces`, so all use the same `pocket_dir/pisa/errors.json`. In `search`,
fetch writes it, then `expand_fsdb_pdb_targets` against a PDB database. A second call with failures replaces the first call's
report, and a clean second call leaves the first's file in place. Date it by its mtime, not by its
existence.

Two things that file does not say about itself:

- **Its `assembly_parsing` entries are distinguishable only by shape** —
  `<pdb_code>_<assembly_id>_<interface_id>` for an interface skipped for not having exactly two
  molecules, `<pdb_code>_<assembly_id>` or a `_parse_error` suffix for a genuine failure. Skips used to
  dominate, for multi-character chain ids; interfaces are now keyed by comma-joined chain ids
  (`"A,B-2"`), so none is skipped for that, and over a 13,063-assembly cache none had a molecule count
  other than two.
- **An entry whose interfaces are all skipped gets no `<pdb_code>.json`.** It therefore never enters
  the interface cache and is reparsed on every later run — from cached assemblies, so at no request
  cost, but it is also why such an entry reappears in every report.

**The cache lives under `pocket_dir/pisa/`; an existing `pocket_dir/pisa_responses/` is dead** and can
be deleted. Nothing migrates it, so the first run after the move refetches every entry from PISA. Its
`interfaces/` held files keyed by concatenated chain ids (`"BF"`), which the parser no longer reads, and
the flattened cache keeps the name `interface_pairs/` so it can never be mistaken for one.

**The input grammar still allows only single-character chains**, so a multi-character interface is
reachable from the cache but not from an input entry. On the Foldseek PDB-DB path,
`expand_fsdb_pdb_targets` checks each chain pair against the `pisa` pattern and skips the ones it
cannot spell, logging one count rather than a warning per record.
