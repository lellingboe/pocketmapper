# Pockets

Facts about `pockets/`, moved out of `context/project-overview.md` so they load only when working here.

- **`POCKET_BUILDERS` in `pocket_fetcher` is the whole method table.** Builders live beside their
  primitive and all take `(records, pocket_dir)`, plus keyword options the caller passes per method
  through the module function `fetch_pockets(builder_options=...)` (pisa's `pisa_source`,
  `bsa_threshold` and `download` today). A new method = one row + one builder.
- **One builder contract.** A builder returns pocket_id -> Pocket for the pockets it built only. A record
  it cannot build (missing or unreadable structure, missing chain, unusable residue, unresolved interface)
  is left out with a warning, never mapped to None; the step then marks it `pocket_not_built` in
  `entries.tsv`. A pocket built but empty (pisa with no residue under the rule, vdw with no contacts) is
  kept, with `pocket_exists=False`. Each builder reads every structure once per call
  (`records_by_structure`, through `lib.structure_manager.load_structure`, which gives None for a file it
  cannot use).
- Records are dicts, not DataFrames. The fetcher passes every record it is given to its method's
  builder, so a caller filters out records whose structure is unavailable first.
- **`expand_fsdb_pdb_targets` lives in `pockets.py`, not beside the builders** (it builds records and
  downloads). It shares `pocket_dir/pisa/` with `pisa_pockets`; both go through
  `pockets.pisa.download_pisa_interfaces`, which picks `pisa/`, and `PisaDownloader` names its contents.

## PISA

`pisa_downloader` fetches through `lib.download`; its pacing, retry and `.part` behaviour is in
`lib/CLAUDE.md`. The files it writes itself (the split summaries, `interface_pairs/<code>.json` and
`errors.json`) go through `write_json`, which also writes a `.part` file and moves it into place, so an
interrupted run cannot leave a truncated interface file that a later run takes as cached.

**PISA assembly interfaces have two sources, one cache.** `source="ftp"` (the default) fetches
`ftp.ebi.ac.uk/.../pdb-assemblies-analysis/split/<code[1:3]>/<code>_assembly<asm>_interfaces.json`
with `download_file`, since those are static files; `source="api"` uses the paced `download_api`. Both
are one `download_many` call over `ASSEMBLY_SOURCES`, pooled `max_workers` wide for FTP and serial for
the API. Summaries, which give the assembly ids, come from the API either way. Both write
`assemblies/<code>_<asm>.json`, because the two serve the same JSON (checked on 2026-09-29: 6rkw
assembly 1 byte-for-byte, 4q5j and 12 4dx9 assemblies equal as parsed JSON). So a cached file is
source-agnostic, and with a warm cache the toggle fetches nothing. **If the formats ever
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
`pisa_pockets` and `expand_fsdb_pdb_targets` both reach it through `download_pisa_interfaces`, so both
use the same `pocket_dir/pisa/errors.json`. In `search`, `expand_fsdb_pdb_targets` writes it against a
PDB database, then the pockets step's `pisa_pockets`. A second call with failures replaces the first call's
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

**Offline, "cached" means `interface_pairs/<code>.json` exists, and nothing else is read.**
`PisaDownloader(download=False)` (through `download_pisa_interfaces` and the `pisa_pockets` builder option
`download`) runs no stage for an entry without that file: one warning gives their count, no request is
sent, and nothing is written, `errors.json` included. `pisa_pockets_from_records` then finds no interface
file for those records and leaves them out, so the step marks them `pocket_not_built`. Accepted
limitation: an entry with cached summaries and assemblies but no `interface_pairs/<code>.json` (an
interrupted run, or every interface skipped, as above) is not reparsed offline, although that would cost
no request. The download stages called directly ignore the flag.

**The cache lives under `pocket_dir/pisa/`; an existing `pocket_dir/pisa_responses/` is dead** and can
be deleted. Nothing migrates it, so the first run after the move refetches every entry from PISA. Its
`interfaces/` held files keyed by concatenated chain ids (`"BF"`), which `load_interfaces` does not read, and
the flattened cache keeps the name `interface_pairs/` so it can never be mistaken for one.

**A pisa pocket is the domain chain's bonded residues, or with `bsa_threshold` its buried ones.** The
buried rule reads ASA and BSA from the interface molecule already in `interface_pairs/`, so it needs no
download and leaves the cache unchanged. A residue joins when `100 * bsa / asa` is strictly greater than
the threshold, and one with ASA 0 never joins. So 0 gives a superset of the bond pocket: on the warm e2e
cache every bonded residue has ASA > 0 and BSA > 0. When a seq id appears twice in a molecule
(microheterogeneity, e.g. 4Z0Y 252 as HS8 then HIS), its first row decides.

**The input grammar still allows only single-character chains**, so a multi-character interface is
reachable from the cache but not from an input entry. On the Foldseek PDB-DB path,
`expand_fsdb_pdb_targets` checks each chain pair against the `pisa` pattern and skips the ones it
cannot spell, logging one count rather than a warning per record. The patterns also reject a symmetry
copy (`A-2`), which a library caller can still name: as the domain chain it is left out of `pisa_pockets`
with a warning, because the structure holds only the asymmetric unit.
