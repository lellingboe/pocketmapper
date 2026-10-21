# Pockets

Facts about `pockets/`, moved out of `context/project-overview.md` so they load only when working here.

- **`POCKET_BUILDERS` in `pocket_fetcher` is the whole method table.** Builders live beside their
  primitive and all take `(records, pocket_dir)`, plus keyword options the caller passes per method
  through `fetch_pockets(builder_options=...)` (pisa's `pisa_source` today). A new method =
  one row + one builder.
- Records are dicts, not DataFrames. The fetcher builds every record it is given.
- **`expand_fsdb_pdb_targets` lives in `steps/pockets.py`, not in `pockets/`** (it builds records and
  downloads). It shares `pocket_dir/pisa/` with `pisa_pockets`; both go through
  `pockets.pisa.pisa_cache_layout`, the one place the cache layout is spelled out.
