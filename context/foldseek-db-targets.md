# Foldseek-DB targets

Cross-module rules for `foldseek_db` target records and the bundled DBs.

## Where it is stated

- Superpose's target structures (`createsubdb` + `convert2pdb`): `lib.foldseek.extract_fsdb_structures`.
- How target ids are read: `StructureAligner.align_structs`.
- Rejected combinations (`--align_struct_method pocket`, `--aligner seq`): `lib.settings.check_fsdb_*`.
- Target pockets: `PocketManager.expand_pdb_hits` (PDB DB) and
  `compare.pocket_comparison.synthesise_target_pocket` (other DBs).
- UniProt renumbering: `lib.foldseek.bundled_offset_table`, `compare.pocket_comparison.uniprot_res_ids`.
- The `.dbtype` marker and staged download: `lib.structure_manager.fetch_missing_fsdb`.

## Rules

- **fetch and align fetch no target structures**; against a PDB DB the pockets step fetches the hits'
  asymmetric-unit structures.
- **`createsubdb` rejects `--threads`**; it is the only one of the five Foldseek subcommands that does, which is
  why `run_foldseek` adds no flags.
- **Every hit of a PDB DB is a `chains` key**, `[]` when it has no pocket, because
  `compare.compare.check_coverage` raises for an alignment target missing from `chains`.
- **One `pocket_id` can sit behind two `preprocess_name`s** (`4q5j-assembly1_B`, `-assembly2_B` →
  `4Q5J:B_F`). `existing_calcs` scores only the first; the transform is whichever assembly Foldseek listed
  first.
- **Pockets come from the asymmetric unit, Foldseek's `tseq` from the assembly.** They agree normally (4Q5J
  self-comparison: `overlap_count == pocket_len`, RMSD ~1e-14); a populated `incorrect_mapping.json` signals
  divergent numbering.
- **UniProt renumbering applies only to the bundled `human_domains`** and changes only `target_overlap_ids`.
  A missing entry aborts the run: a per-row fallback would mix coordinate systems in one column. **Refresh
  `offset_table.tsv` whenever `BUNDLED_HUMAN_DOMAINS_DB` moves.**
- **The `pdb` DB counts as cached when `<path>.dbtype` exists.** Known limit: a partial DB written in place
  by older code, with its `.dbtype`, is still trusted; delete it to refetch.
- **The bundled DB ships without `.source`**; strip it from any refresh. Nothing reads it (`easy-search`,
  `createsubdb` and `convert2pdb` work without it). `.lookup` must stay. The DB is otherwise at its size
  floor: `--coord-store-mode 2` is the smallest mode, and Foldseek cannot read a compressed DB.
- **No cap on enriched hits**, by choice, beyond Foldseek's own `--max-seqs 5000` in `align.py`. `4Q5J:B_F`
  against the bundled `pdb` gives ~4,970 hits, which takes hours on a first run with `pisa_source` `api`;
  reruns hit the cache. Add a cap here if needed.
