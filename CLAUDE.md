# CLAUDE.md

Guidance for Claude Code (claude.ai/code) in this repo.

**Keep this file terse.** Facts, not prose. When editing it, compress rather than expand: one idea per
sentence, no throat-clearing, no restating what the code already says. Add a line only if breaking the
rule would cause silently wrong output rather than an error.

## What this is

PocketMapper compares binding surfaces (pockets) of protein chains. CLI (`pocketmapper`) and importable
library (see "As a library"). It fetches structures (PDB/AlphaFold), derives pocket residues (PISA
interfaces, explicit residue lists, VdW contacts, or a whole chain for an open search), aligns query to
target chains (BLOSUM62 or Foldseek), maps pocket residues through the alignment, writes a comparison table.

## Commands

```bash
pip install -e ".[dev]"            # dev install
black ./ -l 120                    # format (CI: black ./ --check -l 120)
flake8                             # lint (max-line-length 120, E501 ignored)
pre-commit install                 # black + flake8 on commit
bump2version rc|patch|minor|major  # bumps pyproject.toml + pocketmapper/__init__.py, commits and tags
```

**No unit tests**, and the `test` job in `.github/workflows/test_and_deploy.yml` only installs deps — lint is CI's only real gate. Testing is the
end-to-end suite in `tests/e2e/`:

```bash
tests/e2e/run_e2e.sh --list          # the 19 cases and their tags
tests/e2e/run_e2e.sh -t core         # ~35s, no human_domains searches
tests/e2e/run_e2e.sh -o /tmp/pm_e2e  # everything, results under a chosen directory
tests/e2e/run_e2e.sh test_7          # one case by name
```

Cases are grouped and numbered in that order: 1–9 structure-vs-structure pairs (all `core`), 10–14
human_domains DB targets, 15–16 larger Foldseek DB targets, 17–19 the local aligner. Blank lines separate
groups in the `CASES` heredoc and are skipped by the runner (a `#` comment there would *not* be). Adding a
case means renumbering what follows, so describe cases by behaviour, not number.

Each case shells out to the real CLI against live wwPDB / AlphaFold / PDBe PISA — no mocks, network
required. It asserts exit status plus presence (and, where hits are expected, non-emptiness) of
`pocket_comparison.tsv`.

**Always run against the existing cache at `tests/e2e/e2e_results/pocketmapper_cache`** — via
`POCKETMAPPER_E2E_CACHE="$PWD/tests/e2e/e2e_results/pocketmapper_cache"`, or `--cache-dir` for a direct
`pocketmapper search`. It holds several GB of structures, PISA responses and Foldseek DBs; a cold cache is
dramatically slower (PISA is fetched per entry behind a rate-limiting sleep — a full-PDB run is thousands of
calls at ~3/s). The runner's `-o` default is `$PWD/e2e_results`, *not* relative to the script: running from
the repo root silently creates a second empty cache and re-downloads everything. To test cold-cache
behaviour, point at a throwaway dir *inside* `tests/e2e/e2e_results/` so downloads stay reusable.

Missing resources SKIP rather than fail: `needs-pdb-fsdb` wants `POCKETMAPPER_PDB_FSDB` pointed at a
prebuilt Foldseek PDB database; `needs-pdb-download` downloads the full PDB Foldseek DB (2GB download, 7GB
unzipped) so it runs only when named explicitly.

Foldseek is an optional external binary (`conda install -c conda-forge -c bioconda foldseek`), invoked via
`subprocess.run`, required for any `foldseek_db` target, and the **default aligner** (see "The `foldseek`
setting"). So a case is assumed to need it and is skipped when absent; a case opts out with an explicit
`--foldseek False` in its `args`. Cases tagged `local` carry that flag. Keep it that way — when every case
ran Foldseek the suite couldn't see the local branch, which is how it shipped broken. The skip gate matches
`--foldseek False` *before* the catch-all, keeping the two forms distinct.

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

**It does not work at all on the local-aligner path.** `SequenceAligner` writes `"-"` for the `u`/`t`
transforms, which `foldseek_transform` catches and logs per record, so a non-Foldseek run completes but its
aligned PDB contains only the query.

## As a library

Nothing depends on being launched from a terminal — fire wraps `PocketMapper` only at the `main()` boundary.
Two levels of entry:

```python
# 1. Whole pipeline. Same work as the CLI; writes the same files into results_dir.
from pocketmapper.pocketmapper import PocketMapper
PocketMapper().search(query="4Q5J:B_F", target="4Q5J:A_E", results_dir="./out")

# 2. Individual components. None take a Settings -- they are driven by explicit values.
from pocketmapper.qt_processor import QTProcessor
from pocketmapper.structure_fetcher import StructureFetcher
from pocketmapper.lib_struct import parse_pocket_from_struct
qtprocessor = QTProcessor(
    structure_dir="./structs", foldseek_preprocessed_structure_dir="./preproc", fsdb_dir="./fsdb"
)
# One call per side; `name` labels it in logs and errors.
query_df = qtprocessor.process_qt_cmdline_input("4Q5J:B_F", name="query")
target_df = qtprocessor.process_qt_cmdline_input("4Q5J:A_E", name="target")
fetcher = StructureFetcher()
fetcher.set_output_directory("./structs"); fetcher.update_cache()
fetcher.fetch_structures([{"struct_type": "pdb", "struct_info": "4Q5J"}])
pocket = parse_pocket_from_struct("./structs/4Q5J.cif.gz", "B", [100, 101, 102])
```

Level 2 covers `QTProcessor`, `StructureFetcher`, `StructurePreprocessor`, `PisaDownloader`, `PisaParser`,
`SequenceAligner`, `StructureAligner`, `PocketCalculator`, and the `lib`/`lib_struct`/`pocket_comparison`
functions. **`Settings` stays inside `pocketmapper.py`** — resolved there and unpacked at each call site into
the values a component needs. Keep it that way: a component reaching into a `Settings` can't be used without
building one, and hides which fields it depends on.

- **`pocketmapper/__init__.py` only exports `main` and `__version__`.** Submodules are reachable as
  `pocketmapper.lib` etc. only as a side effect of `pocketmapper.pocketmapper` importing them — always use
  explicit `from pocketmapper.<module> import <name>`.
- **Always call `Settings(...).resolve_paths()`** if you build a `Settings` yourself. Unresolved derived paths
  are `None`, giving an opaque `TypeError: expected str, bytes or os.PathLike object, not NoneType` from deep
  inside `os.path.join`. `search()` does this for you.
- **`search()` has global side effects**: `logging.config.dictConfig` reconfigures the *root* logger and will
  stomp on a host app's logging setup. It creates directories and, at the end, `shutil.rmtree`s
  `query_dir`/`target_dir`/`foldseek_tmp_dir` (`_delete_tmp`, marked `# TODO this is unsafe`) — don't point
  those at anything you want to keep.
- **Failures raise `PocketMapperError` rather than exiting**, which is what makes embedding viable; catch it
  instead of the `sys.exit(1)` `main()` turns it into. Preserve this when adding error paths.
- Results come back through files, not return values — `search()` returns `None`, so read
  `pocket_comparison.tsv` / `alignment.tsv` from `results_dir` (paths available on `Settings`).

## Conventions

- **Logging.** Root formatter is `"%(levelname)s: %(stage)s - %(msg)s"`, so *every* log call must pass
  `extra={"stage": "..."}` (usually a local `stage` dict or `self._log_extra`) or the record fails to format.
  Verbosity is numeric: 4=DEBUG, 3=INFO (default), 2=WARNING, else ERROR.
- **Errors.** `logging.critical(...)` then `raise PocketMapperError(...)` (from `pocketmapper.exceptions`);
  `main()` catches it and exits 1. No `exit()`/`sys.exit()` inside modules — that was deliberately removed. The
  one remaining `exit()` is `_check_help_search`, which prints `HELP_MESSAGE` and quits on `--help`; it's a CLI
  affordance, so `PocketMapper().search(help=True)` will kill a host process — library callers must not pass it.
- **Settings.** A frozen-by-convention dataclass updated via `dataclasses.replace`. Resolution order: defaults
  → `--settings` JSON file → explicit CLI args → `resolve_paths()` (which only fills derived paths still `None`,
  so a settings file can pin any individual path). New options go on the dataclass *and* in the `cli_overrides`
  dict in `_configure_workflow` *and* in the `search()` signature — and, if user-facing, in `HELP_MESSAGE`
  (`constants.py`) and the README's Options list, neither generated from the dataclass.
- **The `foldseek` setting is tri-state.** `Settings.foldseek` defaults to `None` ("auto"); `_resolve_foldseek`
  collapses it to a `bool` — `None` → foldseek if `shutil.which("foldseek")` finds it, else a warning and the
  local aligner; `True` → hard requirement, `PocketMapperError` if absent; `False` → local, binary never probed.
  Everything downstream (`_configure_query_target`, `_alignment`, `_delete_tmp`) sees only `True`/`False`, so
  branch on `self._settings.foldseek` rather than re-testing for `None`.

  Its call site is load-bearing: in `_configure_workflow` **after** `_configure_logging` (the root logger is at
  `CRITICAL` until then, so the fallback warning would vanish) and **before** settings are logged and dumped, so
  `job_settings.json` records the resolved value not `null`. It also runs ahead of `_configure_query_target` and
  all fetching, so an unmet `--foldseek True` fails immediately rather than as a raw `FileNotFoundError` from the
  first subprocess after everything has downloaded. Keep new foldseek-availability logic there, not at the
  subprocess call sites.
- **Fetcher/preprocessor API.** `StructureFetcher` and `StructurePreprocessor` follow `set_output_directory()`
  → `update_cache()` → `fetch_*`/`preprocess_records()`. Caching is a plain `os.listdir` snapshot, so
  `update_cache()` must be called after the output dir is set and before work begins.
- **Structure parsing** is gemmi throughout (`.cif.gz` on disk); Biopython only for pairwise alignment and SVD
  superposition.

## Repo layout notes

- `lib.py` holds only generic stateless helpers — `jsonify_dict`, `safe_filename`, the BLOSUM62 matrix reader,
  the similarity scorers. Nothing in it knows about `Settings`, the pipeline or the pocket dict shape; keep it
  that way and put workflow logic in a component module. It used to be a grab-bag: superseded copies of the
  preprocessing, pocket-calculation and PISA-download logic were deleted in favour of the class-based modules,
  and `compare_pockets` moved to `pocket_comparison.py`.
- `PocketCalculator.atp_pocket_overlap` is uncalled but retained for planned ATP-pocket work — don't prune it.
- `pocket_comparison.py` owns step 6 — `POCKET_COMPARISON_COLUMNS`, `compare_pockets` and its helpers
  (`_map_pocket_into_alignment`, `_compare_pocket_pair`, `_superpose`, …). Its column-order contract lives in
  `constants.ALIGNMENT_COLUMNS` because three modules share it. Several caches are local to one
  `compare_pockets` call (`_describe_pocket`, `_seq_identity`): they hold values depending only on a pocket, not
  the alignment row, which matters when thousands of rows name the same query. With no unit tests,
  behaviour-preserving changes here are best checked by capturing `compare_pockets`' arguments from a real run
  and diffing old against new output.
- `constants.py` holds `SINGLE_AA_CODE` (the one three-to-one letter table — duplicates elsewhere were deleted;
  it maps `SEP`/`TPO`/`PTR`/`MSE` and every caller defaults unknowns to `"X"`), `HELP_MESSAGE` (the `--help`
  text, which also documents the settings-file-only "Advanced Options"), and
  `ALIGNMENT_COLUMNS`/`FOLDSEEK_FORMAT_OUTPUT`. It's the right home precisely because it imports nothing.
- `blosum62.bla` and `human_domains/` (a bundled Foldseek DB) ship as package data — see
  `[tool.setuptools.package-data]` in `pyproject.toml`.
- `tests/e2e/fixtures/` holds the batch input files and a local `4Q5J.cif.gz`. Cases run with that directory as
  their working directory because `testfile.txt` refers to `4Q5J.cif.gz` by relative path — which is also what
  makes it a local-file-input test. Keep that relative reference if you edit the fixtures.
- `build/` and `dist/` are stale checked-in artifacts of an older version. Ignore them; never edit
  `build/lib/pocketmapper/`.
