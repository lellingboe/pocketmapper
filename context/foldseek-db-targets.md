# Foldseek-DB targets

How a `foldseek_db` target record is aligned, pocketed and renumbered, and how the bundled DB is kept.
Moved out of `context/project-overview.md` so it loads only when working on Foldseek-DB targets.

With a `foldseek_db` target record, no target structures are fetched; `lib.foldseek.extract_fsdb_structures`
rebuilds the superpose selection via `createsubdb` + `convert2pdb`.

- **`createsubdb` rejects `--threads`** (exits non-zero); it is the only one of the five subcommands that
  does. Hence `run_foldseek` adds no flags; each caller builds its own list.
- `align_structs` interprets target ids once per call: given target records → `pocket_id`s mapped through
  `preprocess_name` (PDB DB; `superpose.step` builds them from `pockets.json` `chains`, hit names in
  alignment order); none → entry names (other DBs). So `fsdb_structures/` never mixes the two.
- Target pocket: `pockets.step.expand_fsdb_pdb_targets` (PDB DB; needs `--alignment`, an FSDB target
  without one is an error) or `synthesise_target_pocket` (other). **Every hit of a PDB DB is a
  `chains` key**: `[]` for a non-PDB name, no interface or only unspellable pairs; a hit whose
  structure download failed lists its pocket_ids with `null` pockets.
- **One `pocket_id` can sit behind two `preprocess_name`s** (`4q5j-assembly1_B`, `-assembly2_B` →
  `4Q5J:B_F`). `existing_calcs` scores only the first; the transform is whichever assembly Foldseek listed first.
- **Pockets come from the AU, Foldseek's `tseq` from the assembly.** Verified to agree normally (4Q5J
  self-comparison: `overlap_count == pocket_len`, identity 1.0, RMSD ~1e-14); a populated
  `incorrect_mapping.json` signals divergent numbering.
- `--align_struct_method pocket` is rejected for any FSDB target (`lib.settings.check_fsdb_align_struct_method`,
  in search before the download and in superpose); `--aligner seq` likewise (`check_fsdb_aligner`, search
  and align).

**UniProt renumbering** (`offset_table.tsv`; resolved in `compare.step`, applied in
`synthesise_target_pocket`):
- Only bundled `human_domains`, looked up by resolved DB path in `lib.foldseek.bundled_offset_table`. User
  DBs keep 0-indexed positions, logged at INFO.
- Only `target_overlap_ids` changes — verified against a same-environment baseline.
- Missing entry or short spec aborts the run (no per-row fallback: it would mix coordinate systems in one
  column). **Refresh the table whenever `BUNDLED_HUMAN_DOMAINS_DB` moves.**

**Bundled DB ships without `.source`**; strip it from any refresh. It duplicates `.lookup` and nothing
reads it (verified: `easy-search`, `createsubdb`, `convert2pdb` all work without it); saves 1.4 MB.
`.lookup` must stay (`extract_fsdb_structures` reads it).

The DB is otherwise at its floor: `_ca` is 53 of 77 MB, 8.4M residues at 6.33 B each
(`--coord-store-mode 2`, smallest mode). zstd -19 saves only 30% and foldseek cannot read a compressed DB.

**No cap on enriched hits**, by choice. `4Q5J:B_F` vs bundled `pdb`: ~4,970 hits / ~3,620 entries, hours
on first run with `pisa_source` `api` (per-assembly PISA behind a sleep); `ftp` is concurrent but untimed
at this scale; reruns hit the cache. Add a cap here if needed.
