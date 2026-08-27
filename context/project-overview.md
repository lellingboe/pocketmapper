# Project overview

## What this is

PocketMapper compares binding surfaces (pockets) of protein chains. CLI (`pocketmapper`) and importable
library (see "As a library"). It fetches structures (PDB/AlphaFold), derives pocket residues (PISA
interfaces, explicit residue lists, VdW contacts, or a whole chain for an open search), aligns query to
target chains (BLOSUM62 or Foldseek), maps pocket residues through the alignment, writes a comparison table.
## Pipeline

`main()` → `fire.Fire(PocketMapper())`, so **every public method on `PocketMapper` is a CLI subcommand** —
hence the leading underscores on all internals, to keep them out of fire's help. `search()` in
`pocketmapper/pocketmapper.py` is the whole workflow, top to bottom:

1. `_configure_workflow` → `Settings`, creates dirs, dumps `job_settings.json`, configures logging.
2. `_configure_query_target` → `QTProcessor` parses `--query`/`--target` into two DataFrames of `QTRecord`s
   (one `process_qt_cmdline_input` call per side).
3. `_fetch_missing_structures` (or `_fetch_missing_fsdb`) → mmCIF from wwPDB/AlphaFold into `structure_dir`.
4. `_alignment` → `_foldseek_preprocessing` + `_foldseek_alignment`, or `_local_alignment`. Both write `alignment.tsv`.
5. `_get_pockets` → union of `_retrieve_pisa_pockets` | `_retrieve_passthrough_pockets` | `_retrieve_vdw_pockets`
   | `_retrieve_whole_chain_pockets`.
6. `_compare_pockets_based_on_alignment` → `pocket_comparison.compare_pockets` → `pocket_comparison.tsv`.
7. `_align_structs` → superposes the top `align_count` targets onto each query into `aligned_structures/`.

### Input grammar

`struct_info[:chain_info[:residue_info]]`, colon-separated (the README's `4Q5J_B_F` form is stale). Either
side may instead be a path to a file with one such string per line. `QTProcessor` infers:

- `struct_type` — `pdb` (4-char ID), `alphafold` (UniProt accession), `local_file` (existing path),
  `foldseek_db` (bundled DB name: `human_domains`, `pdb`).
- `pocket_method` — `whole_chain` (`A`, or nothing), `pisa` (`A_B`), `passthrough` (`A:1,2,3`), `vdw`
  (`A_B:1,2,3`); chosen by regex, constrained by `struct_type`, overridable with `--query_pocket_method` /
  `--target_pocket_method`.

Both trailing parts are optional; dropping them means "no pocket" — an **open search** (below). `chain_info`
then defaults to `constants.DEFAULT_CHAIN` (`"A"`), so `4Q5J` is `4Q5J:A`. The `whole_chain` regex is checked
**first** because a bare chain also matches passthrough and vdw; `passthrough_regex` requires ≥1 residue, so
an empty residue list can't reach `_retrieve_passthrough_pockets` and die on `None.split(",")`.

The original input string is kept as `pocket_id`, the identifier used throughout the results.

A local-file entry like `4Q5J.cif.gz:B_F` resolves to `vdw`, not `pisa` — `B_F` matches the vdw regex and
PISA is PDB-only. That's how the mixed-input fixtures reach the vdw code.

### Pocket dict shape

Every pocket method returns the same nested dict, produced/extended by `lib_struct.parse_pocket_from_struct`:
top-level `res_auth_ids` (author seqids as strings), `ca_sequence`, `pocket_exists`, `has_coords`, plus one
entry per residue keyed by the **string** author seqid holding `res_code`, `res_code_single`, `seq_pos`
(0-based index among CA-bearing residues — what maps into the alignment) and `ca_coords`. Residues without a
CA get `seq_pos = -1` and are excluded, because Foldseek only sees CA-bearing residues.

Top level also carries `whole_chain`, set by `parse_pocket_from_struct` from whether it got a residue list or
`None`. `compare_pockets` branches on it to suppress the `pocket_2_*` columns, so it is a property of each
pocket, not of the run (see "Open searches").

**`compare_pockets` never writes to a pocket dict.** It used to `deepcopy` both sides per alignment row to
hang `fs_pos`/`fs_res_code`/`foldseek_pos` off them; that projection is now a returned `_MappedPocket`, and
pockets are read straight out of `pocket_dict`. Keep it that way — the copy was step 6's largest cost, and a
pocket method relying on mutation would now silently see nothing. The projection is computed once per pocket
per alignment row and shared across pairings, but `unknown_ids` names *both* pockets per entry, so its
records are replayed per pairing (`_record_code_mismatches`).

### Open searches

An entry naming a structure but no pocket (`4Q5J:B`, or `4Q5J` for the default chain) is an *open search*:
`_retrieve_whole_chain_pockets` calls `parse_pocket_from_struct(..., pocket_residues=None)`, treating every
CA-bearing residue of the chain as the pocket. Otherwise an ordinary pocket — keyed by `pocket_id`, joined
through `preproc_to_ids`, carrying residue codes and CA coords — so nothing downstream special-cases it and
`_align_structs` superposes these targets normally.

Only the output differs: with no target pocket to describe (and its length would dilute every ratio),
`compare_pockets` leaves `pocket_2_res_ids`, `pocket_2_len`, `pocket_2_seq`, `pocket_2_pct_aln` and
`jaccard_index` empty — the same shape as a `human_domains` row. Everything else is written, including
`pocket_2_overlap_ids` (real author seqids, only the overlapping ones) and the RMSD/transform columns, which
the synthesised Foldseek-DB pocket cannot produce for lack of coordinates.

**Suppression is per pocket, not per run.** It branches on `p2.get("whole_chain")`, so one run can mix an
open and a pocketed target and emit both row shapes. The remaining global flag, `synthesise_target_pockets`
(formerly `alphafold`), now means only "the target side has no records at all, build a pseudo-pocket from the
alignment row" — the non-PDB Foldseek-DB case, nothing else. That synthesised pocket sets `whole_chain` on
itself, which is how it keeps its column suppression.

The residue-code sanity check against Foldseek's alignment is gated on `"res_code_single" in p2[res]`, not on
`whole_chain`: the synthesised pocket has no residue codes and must skip it; a real whole-chain pocket has
them and the check is worth running.

A `pocket_2` value is not guaranteed to be a target. A query and target sharing a chain share a
`preprocess_name`, so `compare_pockets` pairs every pocket on that chain with every other and some rows carry
a query-only `pocket_id` in `pocket_2`. `_align_structs` filters those out before its `.loc` lookup; without
that it raises a bare pandas `KeyError`.

### Foldseek-DB targets

When the target is a bundled Foldseek DB, `self.fsdb_target` is set: no target structures are fetched or
preprocessed, and `_align_structs` reconstructs target PDBs via `foldseek createsubdb` + `convert2pdb`. A
`foldseek_db` target with foldseek off is a hard error whose message splits on `self._foldseek_available` —
binary missing vs. user-disabled have different fixes.

What the target "pocket" is depends on the DB:

- **A PDB DB** (`pdb`, or a local prebuilt one) yields real PDB chains, so real PISA pockets.
  `_expand_fsdb_pdb_targets` runs first in `_get_pockets`: reads hit names from `alignment.tsv`, resolves each
  via `lib.parse_foldseek_pdb_entry_name`, asks PISA which chains each hit chain touches, and appends one
  ordinary `pisa` record per interface to `_target_df` — `_retrieve_pisa_pockets` then handles them like any
  other, so no separate pocket code exists. `self._fsdb_pdb_target` is set and `compare_pockets` runs with
  `synthesise_target_pockets=False`, populating every `pocket_2_*` column. **Hits with no usable PISA data are
  dropped**, not compared against a stand-in.
- **Any other DB** (`human_domains`, from AlphaFold models) keeps the old behaviour:
  `compare_pockets(synthesise_target_pockets=True)` synthesises a whole-chain pocket per hit; `pocket_2_*`
  stays empty.

Three load-bearing points on the PDB path:

**Generated records carry the Foldseek entry name as `preprocess_name`**, not the one `QTProcessor` derives.
That field is the alignment join key — it links these pockets to their alignment rows (via `preproc_to_ids`)
and to their Foldseek transforms (`foldseek_transform` looks up `u`/`t` by it). The records are otherwise
built by `QTProcessor.parse_individual_qt` on a synthesised `"<PDB>:<chain>_<partner>"` string, so that method
is now part of the pipeline, not just CLI parsing — its output shape is depended on in two places.

**The assembly id is deliberately discarded.** `4q5j-assembly1_B` and `4q5j-assembly2_B` both resolve to
`4Q5J:B_F`, so one `pocket_id` can sit behind two `preprocess_name`s. The pocket is computed once and
`compare_pockets`'s `existing_calcs` set scores only the first assembly's alignment row — so the transform
used is whichever assembly Foldseek reported first. `_align_structs` de-duplicates on `pocket_id` when mapping
back to entry names for this reason; without it the `.loc` lookup returns extra rows and superposes twice.

**Pockets come from the wwPDB asymmetric unit while Foldseek's `tseq` comes from the assembly.** These agree
in the ordinary case (verified: a 4Q5J self-comparison through a PDB-named DB gives `overlap_count ==
pocket_len`, identity 1.0, RMSD ~1e-14), and `compare_pockets`'s 0.8 sequence-identity guard catches them when
they don't — a populated `incorrect_mapping.json` signals that an entry's assembly and AU numbering diverged.

`--align_struct_method pocket` is rejected for any Foldseek-DB target (leaving `foldseek`), in `_configure_query_target` before
anything is fetched: a `human_domains` hit has no coordinates at all, and on the PDB path the pocket is fitted
on asymmetric-unit coordinates while the structure superposed is the assembly `convert2pdb` extracts, so the
transform would be applied in the wrong frame. Which kind of DB it is isn't known until `_expand_fsdb_pdb_targets`
has read the hit names, hence the single early rejection covering both.

**No cap on how many hits get enriched**, by choice. `4Q5J:B_F` against the bundled `pdb` DB returns ~4,970
hits across ~3,620 entries, and PISA is fetched per entry behind a sleep, so the first run takes hours. Reruns
are cheap from the interface cache, and `_expand_fsdb_pdb_targets` logs both counts before starting so the
wait is legible. Add a cap here if that becomes untenable.

## Invariants

Breaking one of these generally produces silently wrong output rather than an error. The Foldseek-DB path
carries three more of its own (above).

### Identifiers and join keys

**`seq_pos` is the value everything hinges on** — the residue's index among the CA-bearing residues *of its
own chain*, and what maps a pocket residue into the alignment. A pocket method computing it differently from
`lib_struct.parse_pocket_from_struct` silently produces zero overlap. When adding one, check a pocket against
itself: self-comparison must yield `overlap_count == pocket_len`.

**`preprocess_name` is the join key** — `<basename>_<chain><md5-of-that>` (e.g. `4Q5J_B_<hash>`), computed once
in `QTProcessor.parse_individual_qt`. Alignments are keyed by it, pockets by `pocket_id`;
`_compare_pockets_based_on_alignment` builds `preproc_to_ids` to bridge them. One `preprocess_name` can map to
several `pocket_id`s (same chain, different pockets).

**Aligned structures are named by `lib.safe_filename(query_id)`, not by `pocket_id`.** A `pocket_id` is raw
user input — possibly a path or a long residue list — and `safe_filename` keeps only the basename, sanitises
it and appends an md5 of the original. So `aligned_structures/*.pdb` filenames aren't greppable for an input
string; match on the `MOLECULE` records inside instead.

### Table schemas

**The alignment table's column order is a positional contract**, declared once as `constants.ALIGNMENT_COLUMNS`.
All three parties derive from it: `_foldseek_alignment` passes `constants.FOLDSEEK_FORMAT_OUTPUT` (the same
list, comma-joined) to Foldseek's `--format-output`; `SequenceAligner.align_records` pins its DataFrame to
`columns=ALIGNMENT_COLUMNS`; `pocket_comparison` unpacks each row positionally into an `AlignmentRow`
namedtuple built from the same list. Reorder the constant and all three move together; reads stay positional
(`AlignmentRow(*values)`), they just say which column they mean. Adding a column to one producer alone still
breaks silently — add it to `ALIGNMENT_COLUMNS`.

**The comparison table has a fixed schema.** `pocket_comparison.POCKET_COMPARISON_COLUMNS` declares every
column `compare_pockets` can produce, and the result is reindexed onto it, so all 30 exist on every run — rows
that stop early (no overlap, no coordinates, an open/foldseek-db target with no pocket 2) leave later fields
empty rather than dropping them. Add new output fields to that list as well as to the row dict; a column
produced but not declared is kept and logged as a warning rather than silently dropped.

### Structural alignment (step 7)

**`_align_structs` only superposes targets that actually overlap the query.** It filters on `overlap_count > 0`
before ranking: a target sharing no pocket residues has no common residue set to superpose on and empty
overlap metrics that would sort arbitrarily. It ranks on `jaccard_index`, then `min_overlap_similarity` — and
since a whole-chain target has no `jaccard_index`, an open or Foldseek-DB search sorts every candidate into the
NaN block and is ordered by similarity alone. Queries with no overlapping target are skipped with a log line;
a run where nothing overlaps returns early.

**Two transform sources, chosen by `align_struct_method`.** `foldseek` is the whole-chain `u`/`t` from
`alignment.tsv` (`StructureAligner.foldseek_transform`); `pocket` is the fit of the two pockets on their
overlapping residues, read out of `pocket_comparison.tsv`'s `p2_to_p1_u`/`p2_to_p1_t` and handed to
`StructureAligner.transform`. `SequenceAligner` writes `"-"` for `u`/`t`, so `foldseek` is impossible on the
local-aligner path — rejected in `_resolve_align_struct_method` at settings time, not per record at write
time. `p1` is always the query, so `p2_to_p1_*` needs no inversion.

**`p2_to_p1_u` is Biopython's *right*-multiplying rotation; `gemmi.Transform` *left*-multiplies, as Foldseek's
`u` already does.** `pocket_comparison.parse_pocket_transform` transposes it, and is the only place that may —
never hand a raw cell to gemmi. It also parses the list-repr serialisation those cells use (`_superpose` writes
lists, `to_csv` reprs them), where `alignment.tsv` comma-joins the same quantities.

**Under `pocket` the candidate filter is `overlap_count > 0` *and* a present `p2_to_p1_u`**, applied before
`head(align_count)`. `_superpose` fits nothing below three overlapping residues, so without it those targets
would consume slots and the run would quietly write fewer structures than asked for.

**`StructureAligner.transform` takes transforms positionally, not keyed by `pocket_id`.** A query compared
against itself gives the reference and a target the same `pocket_id`. A record whose transform is `None` is
dropped, and the COMPND header is built from the records that survived — it used to be built from all of them,
naming models the file did not contain.

## Repo layout notes

- `lib.py` holds only generic stateless helpers — `jsonify_dict`, `safe_filename`, the BLOSUM62 matrix reader,
  the similarity scorers. Nothing in it knows about `Settings`, the pipeline or the pocket dict shape; keep it
  that way and put workflow logic in a component module. It used to be a grab-bag: superseded copies of the
  preprocessing, pocket-calculation and PISA-download logic were deleted in favour of the class-based modules,
  and `compare_pockets` moved to `pocket_comparison.py`.
- `PocketCalculator.atp_pocket_overlap` is uncalled but retained for planned ATP-pocket work — don't prune it.
- `pocket_comparison.py` owns step 6 — `POCKET_COMPARISON_COLUMNS`, `compare_pockets` and its helpers
  (`_map_pocket_into_alignment`, `_compare_pocket_pair`, `_superpose`, …), plus `parse_pocket_transform`, which
  step 7 uses to read `_superpose`'s output back. Its column-order contract lives in
  `constants.ALIGNMENT_COLUMNS` because three modules share it. Several caches are local to one
  `compare_pockets` call (`_describe_pocket`, `_seq_identity`): they hold values depending only on a pocket, not
  the alignment row, which matters when thousands of rows name the same query. With no unit tests,
  behaviour-preserving changes here are best checked by capturing `compare_pockets`' arguments from a real run
  and diffing old against new output.
- `constants.py` holds `SINGLE_AA_CODE` (the one three-to-one letter table — duplicates elsewhere were deleted;
  it maps `SEP`/`TPO`/`PTR`/`MSE` and every caller defaults unknowns to `"X"`), `HELP_MESSAGE` (the `--help`
  text, which also documents the settings-file-only "Advanced Options"), and
  `ALIGNMENT_COLUMNS`/`FOLDSEEK_FORMAT_OUTPUT`. It's the right home precisely because it imports nothing.
- `tests/e2e/fixtures/` holds the batch input files and a local `4Q5J.cif.gz`. Cases run with that directory as
  their working directory because `testfile.txt` refers to `4Q5J.cif.gz` by relative path — which is also what
  makes it a local-file-input test. Keep that relative reference if you edit the fixtures.
- `build/` and `dist/` are stale checked-in artifacts of an older version. Ignore them; never edit
  `build/lib/pocketmapper/`.
