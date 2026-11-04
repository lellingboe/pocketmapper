# Pockets

Cross-module facts about `pockets/`.

## Where it is stated

- The method table, the builder contract and builder options: the `pocket_manager` module docstring.
  Residue dropping, empty vs all-dropped pockets and per-structure reading: the `structure` module docstring
  and `parse_pocket_from_struct`. A new method: `parse/CLAUDE.md` "A new pocket method needs".
- `pockets.json` checks: the `pocket_files` module docstring.
- The PISA cache layout, `interface_pairs/<code>.json`, `is_complete`, offline behaviour and the
  `assembly_parsing` entry shapes: the `pisa_downloader` module docstring and its stages. Pacing, retries and
  `.part` files: `lib/CLAUDE.md`.
- The pocket rules (bonded, or buried with `bsa_threshold`): the `pisa` module docstring.

## PISA

- **Two sources, one cache.** `ftp` and `api` write the same `assemblies/<code>_<asm>.json`, because they serve
  the same JSON (checked 2026-09-29). **If the formats ever diverge, the caches must split.** An assembly the
  FTP lacks is a 404, reported like an API failure; there is no fallback to the API.
- **Summaries are the one batched request**; the interface endpoint has no batched form. A batch with no
  known id returns 404 and fails as a whole. Somewhere between 800 and 900 ids the API returns 500, which
  `is_transient_error` retries five times, permanently doubling the host's pacing; hence the default
  `summary_batch_size` of 50.
- **`pocket_dir/pisa/errors.json` is written by each `download_missing_interfaces` call with failures.** The
  pockets step against a PDB DB calls it twice, in `expand_pdb_hits` and then in `pisa_pockets` (for the
  same hit codes), so the second call usually replaces the first's report. A clean second call leaves the
  first's file in place: date it by its mtime.
- **Skipped interfaces are reported once, read errors every run.** An entry whose interfaces were all skipped
  is still written and complete, so it is not parsed again; an assembly that could not be read is never
  listed as used, so its entry is parsed, and reported, on every online run.
- **The cost of the completeness rule**: an assembly PISA permanently lacks (an FTP 404) is requested again
  on every online run, and twice per pockets run against a PDB DB (12 of ~5,800 assemblies in the
  network-test cache, 2026-11-04).
- **Offline, an uncached entry gives no pocket.** With `download=False` no stage runs for it, so
  `pisa_pockets_from_records` finds no interfaces and the step marks the record `pocket_not_built`. In
  `expand_pdb_hits` an uncached hit entry gives no records, so its `chains` key is `[]`; an uncached hit
  structure is `structure_not_found`, and if no hit structure is available the step raises.
- **`bsa_threshold` reads data already in `interface_pairs/`**, so changing it needs no download and leaves
  the cache unchanged. Observed on the warm e2e cache: every bonded residue has ASA > 0 and BSA > 0, so a
  threshold of 0 gives a superset of the bond pocket.

## Single-character chains

The input grammar allows only single-character chains, so a multi-character interface is reachable from the
cache but not from an input entry. On the Foldseek PDB-DB path, `PocketManager.expand_pdb_hits` checks each
chain pair against the `pisa` pattern and skips the ones it cannot spell, logging one count. The patterns
also reject a symmetry copy (`A-2`), which a library caller can still name: as the domain chain it is left
out of `pisa_pockets` with a warning, because the structure holds only the asymmetric unit. To allow longer
chain ids, change the patterns in `parse.qt_processor` and keep `lib.residues.split_chain_info` the only
splitter.
