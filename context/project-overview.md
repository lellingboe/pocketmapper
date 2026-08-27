# Project overview

A map, not a manual. If a comment at the code site already explains something, it does not belong here;
what is left is what spans modules, was measured rather than written down, or was left out on purpose.
Every claim names the file or symbol that proves it.

## What this is

PocketMapper compares binding surfaces (pockets) of protein chains. CLI (`pocketmapper`) and importable
library (see "As a library" in `coding-standards.md`). It fetches structures (PDB/AlphaFold), derives
pocket residues (PISA interfaces, explicit residue lists, VdW contacts, or a whole chain for an open
search), aligns query to target chains (BLOSUM62 or Foldseek), maps pocket residues through the alignment,
writes a comparison table.

## Pipeline

`main()` → `fire.Fire(PocketMapper())`, so **every public method on `PocketMapper` is a CLI subcommand** —
hence the leading underscores on all internals, to keep them out of fire's help. `search()` in
`pocketmapper/pocketmapper.py` is the whole workflow, top to bottom:

1. `_configure_workflow` → `Settings`, directories, `job_settings.json`, logging.
2. `_configure_query_target` → `QTProcessor` → one DataFrame of `QTRecord`s per side.
3. `_fetch_missing_structures` (or `_fetch_missing_fsdb`) → mmCIF into `structure_dir`.
4. `_alignment` → foldseek or the local aligner → `alignment.tsv`.
5. `_get_pockets` → `_retrieve_{pisa,passthrough,vdw,whole_chain}_pockets`, merged into one dict.
6. `_compare_pockets_based_on_alignment` → `pocket_comparison.compare_pockets` → `pocket_comparison.tsv`.
7. `_align_structs` → superposes the top `align_count` targets per query into `aligned_structures/`.

### Input grammar

`struct_info[:chain_info[:residue_info]]`, colon-separated; either side may instead be a file with one such
string per line. The README's "Input format" table documents the forms and
`QTProcessor.determine_struct_type` / `determine_pocket_method` implement them (regexes at `qt_processor.py:68-79`). Two things neither states:

- The `whole_chain` regex is tested **first**, because a bare chain also matches the passthrough and vdw
  patterns.
- A local-file entry like `4Q5J.cif.gz:B_F` resolves to `vdw`, not `pisa` — `B_F` matches the vdw regex and
  PISA is PDB-only. That is how the mixed-input e2e fixtures reach the vdw code.

The original input string is kept as `pocket_id`, the identifier used throughout the results.

### Pocket dict shape

Every pocket method returns the same nested dict, produced and extended by
`lib_struct.parse_pocket_from_struct`: top-level `res_auth_ids`, `ca_sequence`, `pocket_exists`,
`has_coords`, `whole_chain`, plus one entry per residue keyed by the **string** author seqid.

Residues without a CA get `seq_pos = -1` and are excluded from the comparison, because Foldseek only sees
CA-bearing residues. `whole_chain` records whether the method passed a residue list or `None`, and
`compare_pockets` branches on it — so it is a property of each pocket, not of the run (see "Open searches").

### Open searches

An entry naming a structure but no pocket (`4Q5J:B`, or `4Q5J` for the default chain) is an *open search*:
`_retrieve_whole_chain_pockets` treats every CA-bearing residue of the chain as the pocket. The result is an
ordinary pocket, so nothing downstream special-cases it. Only the output differs — with no target pocket to
describe (and its length would dilute every ratio), `_compare_pocket_pair` leaves the `pocket_2_*` descriptor
columns and `jaccard_index` empty.

That suppression is **per pocket**: one run can mix an open and a pocketed target and emit both row shapes.
The per-run flag is `synthesise_target_pockets`, which means only "the target side has no records at all,
build a pseudo-pocket from the alignment row" (`pocket_comparison.py:384-389`).

A `pocket_2` value is not guaranteed to be a target. A query and target sharing a chain share a
`preprocess_name`, so `compare_pockets` pairs every pocket on that chain with every other and some rows
carry a query-only `pocket_id` in `pocket_2`. `_align_structs` filters those out before its `.loc` lookup;
without that it raises a bare pandas `KeyError`.

### Foldseek-DB targets

When the target is a bundled Foldseek DB, `self.fsdb_target` is set: no target structures are fetched or
preprocessed, and `_align_structs` reconstructs target PDBs via `foldseek createsubdb` + `convert2pdb`.
What the target "pocket" is depends on the DB:

- **A PDB DB** (`pdb`, or a local prebuilt one) yields real PDB chains, so real PISA pockets.
  `_expand_fsdb_pdb_targets` appends one ordinary `pisa` record per interface to `_target_df`, so no
  separate pocket code exists and every `pocket_2_*` column is populated. **Hits with no usable PISA data
  are dropped**, not compared against a stand-in.
- **Any other DB** (`human_domains`, from AlphaFold models) keeps `synthesise_target_pockets=True`:
  one whole-chain pseudo-pocket per hit, `pocket_2_*` empty.

Three consequences of the PDB path that are not visible from any single file:

- **One `pocket_id` can sit behind two `preprocess_name`s** (`4q5j-assembly1_B` and `4q5j-assembly2_B` both
  resolve to `4Q5J:B_F`). The pocket is computed once and `compare_pockets`' `existing_calcs` scores only
  the first assembly's alignment row, so the transform used is whichever assembly Foldseek reported first.
- **Pockets come from the wwPDB asymmetric unit while Foldseek's `tseq` comes from the assembly.** These
  agree in the ordinary case (verified: a 4Q5J self-comparison through a PDB-named DB gives
  `overlap_count == pocket_len`, identity 1.0, RMSD ~1e-14), and the `MIN_SEQ_IDENTITY` guard catches them
  when they don't — a populated `incorrect_mapping.json` signals that an entry's assembly and AU numbering
  diverged.
- **`--align_struct_method pocket` is rejected for any Foldseek-DB target**, in `_configure_query_target`
  before anything is fetched: a `human_domains` hit has no coordinates, and on the PDB path the pocket is
  fitted on AU coordinates while the structure superposed is the assembly `convert2pdb` extracts. Which kind
  of DB it is isn't known until `_expand_fsdb_pdb_targets` has read the hit names, hence one early rejection
  covering both.

**No cap on how many hits get enriched**, by choice. `4Q5J:B_F` against the bundled `pdb` DB returns ~4,970
hits across ~3,620 entries, and PISA is fetched per entry behind a sleep, so the first run takes hours.
Reruns are cheap from the interface cache, and `_expand_fsdb_pdb_targets` logs both counts before starting
so the wait is legible. Add a cap here if that becomes untenable.

## Invariants

Breaking one of these generally produces silently wrong output rather than an error. See "Extending the
pipeline" in `coding-standards.md` for what to do about them when adding code.

- **`seq_pos` is the value everything hinges on** — the residue's index among the CA-bearing residues *of
  its own chain*, and what maps a pocket residue into the alignment (`lib_struct.parse_pocket_from_struct`).
- **`preprocess_name` is the alignment join key** — `<basename>_<chain><md5>`, computed once in
  `QTProcessor.parse_individual_qt`. Alignments are keyed by it, pockets by `pocket_id`, and
  `_compare_pockets_based_on_alignment` builds `preproc_to_ids` to bridge them. One `preprocess_name` can
  map to several `pocket_id`s (same chain, different pockets).
- **Two tables have declared schemas**: `constants.ALIGNMENT_COLUMNS` (a positional contract shared by
  `_foldseek_alignment`, `SequenceAligner.align_records` and `AlignmentRow`) and
  `pocket_comparison.POCKET_COMPARISON_COLUMNS` (every column a comparison can produce, always present).
- **Aligned structures are named by `lib.safe_filename(query_id)`, not by `pocket_id`** — the name is
  sanitised down to a basename plus an md5, so `aligned_structures/*.pdb` filenames aren't greppable for an
  input string. Match on the `MOLECULE` records inside instead.
- **Two transform sources, chosen by `align_struct_method`**: `foldseek` reads whole-chain `u`/`t` from
  `alignment.tsv`; `pocket` reads `p2_to_p1_u`/`p2_to_p1_t` from `pocket_comparison.tsv`. `p1` is always the
  query, so `p2_to_p1_*` needs no inversion, and the two cells use different conventions and serialisations
  — `pocket_comparison.parse_pocket_transform` is the only legitimate reader of the latter.

## Repo layout notes

- `pocketmapper.py` — `Settings`, `search()` and every pipeline step; the only module that knows `Settings`.
- `pocket_comparison.py` — step 6 in full, plus `parse_pocket_transform`, which step 7 uses to read
  `_superpose`'s output back.
- `lib.py` — generic stateless helpers only. `constants.py` — imports nothing, which is why the shared
  tables and `SINGLE_AA_CODE` live there. `lib_struct.py` — the pocket dict.
- Components (`qt_processor`, `structure_fetcher`, `structure_preprocessor`, `pisa_downloader`,
  `pisa_parser`, `sequence_aligner`, `structure_aligner`, `pocket_calculator`) are separately usable and
  take explicit values, never a `Settings`.
- `PocketCalculator.atp_pocket_overlap` is uncalled but retained for planned ATP-pocket work — don't prune.
- There are no unit tests. `tests/e2e/` is the whole suite; the `pocketmapper-e2e` skill covers running it.
- `build/` and `dist/` are stale checked-in artifacts of an older version. Ignore them; never edit
  `build/lib/pocketmapper/`.
