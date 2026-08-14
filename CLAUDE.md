# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

PocketMapper compares the binding surfaces (pockets) of protein chains. It ships as a CLI (`pocketmapper`)
and is equally usable as an importable Python library — see "Using it as a library" below. It fetches
structures (PDB / AlphaFold), derives pocket residues (PISA interfaces, explicit residue lists, or VdW
contact calculation), aligns query chains to target chains (BLOSUM62 sequence alignment or Foldseek),
maps pocket residues through the alignment, and writes a comparison table.

## Commands

```bash
pip install -e ".[dev]"          # dev install
black ./ -l 120                  # format (CI runs `black ./ --check -l 120`)
flake8                           # lint (max-line-length 120, E501 ignored)
pre-commit install               # black + flake8 on commit
bump2version rc|patch|minor|major  # bumps pyproject.toml + pocketmapper/__init__.py, commits and tags
```

There are **no unit tests**, and `.github/workflows/test_and_deploy.yml` has a `test` job that only installs
dependencies — CI's only real gate is lint. What exists is an end-to-end suite in `tests/e2e/`:

```bash
tests/e2e/run_e2e.sh --list          # the 13 cases and their tags
tests/e2e/run_e2e.sh -t core         # ~35s, no human_domains searches
tests/e2e/run_e2e.sh -o /tmp/pm_e2e  # everything, results under a chosen directory
tests/e2e/run_e2e.sh test_4          # one case by name
```

Each case shells out to the real CLI and hits live wwPDB / AlphaFold / PDBe PISA — there are no mocks, so the
suite needs network access. It asserts exit status plus the presence (and, where a pair is known to produce
hits, the non-emptiness) of `pocket_comparison.tsv`. Reuse one `--cache-dir` across runs; a warm cache makes
reruns dramatically faster. Cases needing resources that aren't present are reported as SKIP rather than
failing: `test_10` wants `POCKETMAPPER_PDB_FSDB` pointed at a prebuilt Foldseek PDB database, and `test_13`
downloads the full PDB Foldseek DB (tens of GB) so it only runs when named explicitly.

Whether a case uses Foldseek is decided by `--foldseek` in its own `args` field, not by the runner — cases
tagged `local` omit it to exercise the BLOSUM62 aligner, and only the Foldseek ones are skipped when the
binary is missing. Keep it that way: when every case forced `--foldseek`, the suite could not see the local
branch at all, which is how it shipped broken.

### Invariants worth knowing

**The comparison table has a fixed schema.** `lib.POCKET_COMPARISON_COLUMNS` declares every column
`compare_pockets` can produce, and the result is reindexed onto it before returning, so all 33 columns exist
on every run — rows that stop early (no overlap, no coordinates, an alphafold/foldseek-db target with no
pocket 2) leave the later fields empty rather than dropping them. Add a new output field to that list as well
as to the row dict; a column produced but not declared is kept and logged as a warning rather than silently
dropped, so the list cannot quietly drift.

**`_align_structs` only superposes targets that actually overlap the query.** It filters on
`overlap_count > 0` before ranking, because a target sharing no pocket residues has no common residue set to
superpose on and empty overlap metrics that would sort arbitrarily. Queries with no overlapping target are
skipped with a log line, and a run where nothing overlaps returns early instead of proceeding.

Note that `seq_pos` is the single value everything hinges on: it is the residue's index among the CA-bearing
residues *of its own chain*, and it is what maps a pocket residue into the alignment. Any pocket method that
computes it differently from `lib_struct.parse_pocket_from_struct` will silently produce zero overlap rather
than an error — that was the `vdw` bug below. When adding a pocket method, check a pocket against itself:
self-comparison must yield `overlap_count == pocket_len`.

### Recently fixed

**A zero-overlap result crashed the final stage.** `compare_pockets` built each row as a dict and let pandas
infer the columns, so a run in which *every* row had zero overlap never created `pocket_1_pct_overlap` /
`min_overlap_similarity` at all, and `_align_structs` sorting by them raised
`KeyError: 'pocket_1_pct_overlap'`. Mixed runs were unaffected, which is why it only showed up on small
single-pair comparisons. Fixed by the two invariants above: the fixed schema means the columns always exist,
and the `overlap_count > 0` filter means the sort only ever sees rows with real metrics. Beyond the crash,
the run used to abort before `_delete_tmp()`, leaking `query_structures/`/`target_structures/` into the
results directory. `test_7` and `test_8` cover this.

**`vdw` pockets always scored zero overlap.** `ca_num += 1` in `PocketCalculator.pocket_overlap` sat inside
the inner `for res2 in motif_residues` loop, so it counted residue *pairs* and inflated `seq_pos` by the
partner chain's length (3681 where the reference gives 230). Every vdw residue therefore fell outside the
alignment region — `pct_aln` 0.0, `overlap_count` 0 on every row, no warning, and otherwise plausible output.
Dedenting the increment to the outer loop fixes it; `seq_pos` and `ca_sequence` now match
`parse_pocket_from_struct` exactly. Note `atp_pocket_overlap` never had the bug (it has only one residue
loop), which is likely how it arose: `pocket_overlap` looks derived from it by wrapping a `for res2` loop
around the body without dedenting the counter.

This was also the main trigger for the `KeyError` above: a run over vdw pockets scored zero everywhere and so
always crashed in `_align_structs`. Local-file entries like `4Q5J.cif.gz:B_F` resolve to vdw (`B_F` matches
the vdw regex, and PISA is PDB-only), which is how the mixed-input fixtures reach it.

`_local_alignment` used to read each record's `preprocess_path_gz`, which only exists after
`_foldseek_preprocessing()` — a step that runs solely on the Foldseek branch — so every run without
`--foldseek` died with `FileNotFoundError`. `SequenceAligner.align_records` now reads `struct_path` (the full
reference structure) and selects the chain itself, which is what it was already doing to the pre-split copy.
Verified to produce results identical to the Foldseek path on the same pair. `test_14`/`test_15` are the
regression tests; both fail against the previous code. Note that `_align_structs` still can't superpose on the
local path — `SequenceAligner` writes `"-"` for the `u`/`t` transforms, which `foldseek_transform` catches and
logs per record, so the run completes with an aligned PDB containing only the query.

Foldseek is an optional external binary (`conda install -c conda-forge -c bioconda foldseek`); it is
invoked via `subprocess.run` and is required for `--foldseek True` and for any `foldseek_db` target.

## Pipeline

`main()` → `fire.Fire(PocketMapper())`, so **every public method on `PocketMapper` is a CLI subcommand**
(hence the leading underscores on all internals — that is deliberate, to keep them out of fire's help).
`search()` in `pocketmapper/pocketmapper.py` is the whole workflow and reads top to bottom:

1. `_configure_workflow` → `Settings` (see below), creates directories, dumps `job_settings.json`, configures logging.
2. `_configure_query_target` → `QTProcessor` parses `--query`/`--target` into two DataFrames of `QTRecord`s.
3. `_fetch_missing_structures` (or `_fetch_missing_fsdb`) → downloads mmCIF from wwPDB / AlphaFold into `structure_dir`.
4. `_alignment` → either `_foldseek_preprocessing` + `_foldseek_alignment`, or `_local_alignment`. Both write `alignment.tsv`.
5. `_get_pockets` → union of `_retrieve_pisa_pockets` | `_retrieve_passthrough_pockets` | `_retrieve_vdw_pockets`.
6. `_compare_pockets_based_on_alignment` → `lib.compare_pockets` → `pocket_comparison.tsv`.
7. `_align_structs` → superposes the top `align_count` targets onto each query into `aligned_structures/`.

### Input grammar

Query/target strings are `struct_info:chain_info:residue_info` (colon-separated; the README's
`4Q5J_B_F` form is stale). Either side may instead be a path to a file with one such string per line.
`QTProcessor` infers two things from the string:

- `struct_type` — `pdb` (4-char ID), `alphafold` (UniProt accession), `local_file` (existing path), or
  `foldseek_db` (a bundled DB name: `human_domains`, `pdb`).
- `pocket_method` — `pisa` (`A_B`), `passthrough` (`A:1,2,3`), `vdw` (`A_B:1,2,3`), chosen by regex and
  constrained by `struct_type`. Overridable with `--query_pocket_method` / `--target_pocket_method`.

The original input string is kept as `pocket_id` and is the identifier used throughout the results.

### Two invariants that hold the pipeline together

**The alignment table's column order is a positional contract.** `SequenceAligner.align_records` builds the
exact same 18 columns that the Foldseek `--format-output` flag requests
(`query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,lddt,qaln,taln,u,t,qseq,tseq`),
and `lib.compare_pockets` reads them **by index** (`row[0]`…`row[17]`). Changing or reordering columns in one
producer without the other, or without updating the indices in `compare_pockets`, breaks silently.

**`preprocess_name` is the join key.** It is `<basename>_<chain><md5-of-that>` (e.g. `4Q5J_B_<hash>`), computed
once in `QTProcessor.parse_individual_qt`. Alignments are keyed by it; pockets are keyed by `pocket_id`;
`_compare_pockets_based_on_alignment` builds `preproc_to_ids` to bridge the two. One `preprocess_name` can map
to several `pocket_id`s (same chain, different pockets).

### Pocket dict shape

Every pocket method returns the same nested dict, produced or extended by `lib_struct.parse_pocket_from_struct`:
top-level `res_auth_ids` (list of author seqids as strings), `ca_sequence`, `pocket_exists`, `has_coords`, plus
one entry per residue keyed by the **string** author seqid holding `res_code`, `res_code_single`, `seq_pos`
(0-based index among CA-bearing residues — this is what maps into the alignment) and `ca_coords`.
Residues without a CA atom get `seq_pos = -1` and are excluded, because Foldseek only sees CA-bearing residues.

### Foldseek-DB targets

When the target is a bundled Foldseek DB, `self.fsdb_target` is set and several branches change: no target
structures are fetched or preprocessed, `compare_pockets(alphafold=True)` synthesises a whole-chain "pocket" for
each target hit rather than looking one up, and `_align_structs` reconstructs target PDBs from the DB via
`foldseek createsubdb` + `convert2pdb`. A `foldseek_db` target without `--foldseek True` is a hard error.

## Using it as a library

Nothing in the package depends on being launched from a terminal — fire only wraps `PocketMapper` at the
`main()` boundary. There are three levels of entry, from coarsest to finest:

```python
# 1. Whole pipeline. Same work as the CLI; writes the same files into results_dir.
from pocketmapper.pocketmapper import PocketMapper
PocketMapper().search(query="4Q5J:B_F", target="4Q5J:A_E", results_dir="./out")

# 2. Settings-driven components. resolve_paths() is REQUIRED before use (see gotcha below).
from pocketmapper.pocketmapper import Settings
from pocketmapper.qt_processor import QTProcessor
settings = Settings(query="4Q5J:B_F", target="4Q5J:A_E").resolve_paths()
query_df, target_df = QTProcessor(settings=settings).process_qt_cmdline_input()

# 3. Standalone components. These take no Settings and are the easiest pieces to reuse.
from pocketmapper.structure_fetcher import StructureFetcher
from pocketmapper.lib_struct import parse_pocket_from_struct
fetcher = StructureFetcher()
fetcher.set_output_directory("./structs"); fetcher.update_cache()
fetcher.fetch_structures([{"struct_type": "pdb", "struct_info": "4Q5J"}])
pocket = parse_pocket_from_struct("./structs/4Q5J.cif.gz", "B", [100, 101, 102])
```

Level 3 covers `StructureFetcher`, `StructurePreprocessor`, `PisaDownloader`, `PisaParser`,
`SequenceAligner`, `StructureAligner`, `PocketCalculator`, and the `lib`/`lib_struct` functions — all
constructed with no arguments and driven by explicit paths.

Things to know when consuming it this way:

- **`pocketmapper/__init__.py` only exports `main` and `__version__`.** Submodules are reachable as
  `pocketmapper.lib` etc. only as a side effect of `pocketmapper.pocketmapper` importing them — always use
  explicit `from pocketmapper.<module> import <name>` rather than relying on that.
- **Always call `Settings(...).resolve_paths()`** before handing a `Settings` to any component. Unresolved
  derived paths are `None`, and the failure is an opaque `TypeError: expected str, bytes or os.PathLike
  object, not NoneType` from deep inside `os.path.join`, not a helpful error. `search()` does this for you;
  direct component use does not.
- **`search()` has global side effects**: it calls `logging.config.dictConfig`, which reconfigures the *root*
  logger and will stomp on a host application's logging setup. It also creates directories and, at the end,
  `shutil.rmtree`s `query_dir`/`target_dir`/`foldseek_tmp_dir` (`_delete_tmp`, marked `# TODO this is unsafe`)
  — so do not point those at a directory holding anything you want to keep.
- **Failures raise `PocketMapperError` rather than exiting.** This is what makes embedding viable; catch it
  instead of the `sys.exit(1)` that `main()` turns it into. Preserve this when adding error paths.
- Results come back through files, not return values — `search()` returns `None`, so read
  `pocket_comparison.tsv` / `alignment.tsv` from `results_dir` (paths available on the `Settings` object).

## Conventions

- **Logging.** The root formatter is `"%(levelname)s: %(stage)s - %(msg)s"`, so *every* log call must pass
  `extra={"stage": "..."}` (usually a local `stage` dict or `self._log_extra`) or the record fails to format.
  Verbosity is numeric: 4=DEBUG, 3=INFO (default), 2=WARNING, else ERROR.
- **Errors.** Log with `logging.critical(...)` then `raise PocketMapperError(...)` (from `pocketmapper.exceptions`);
  `main()` catches it and exits 1. Do not call `exit()`/`sys.exit()` inside modules — that was deliberately removed.
- **Settings.** `Settings` is a frozen-by-convention dataclass updated via `dataclasses.replace`. Resolution order is
  defaults → `--settings` JSON file → explicit CLI args → `resolve_paths()` (which only fills derived paths still `None`,
  so a settings file can pin any individual path). New options go on the dataclass *and* in the `cli_overrides` dict
  in `_configure_workflow` *and* in the `search()` signature.
- **Fetcher/preprocessor API.** `StructureFetcher` and `StructurePreprocessor` follow
  `set_output_directory()` → `update_cache()` → `fetch_*`/`preprocess_records()`. Caching is a plain `os.listdir`
  snapshot, so `update_cache()` must be called after the output dir is set and before work begins.
- **Structure parsing** is gemmi throughout (`.cif.gz` on disk); Biopython is used only for pairwise alignment
  and SVD superposition.

## Repo layout notes

- `lib.py` is the legacy grab-bag; `compare_pockets`, the BLOSUM62 matrix reader, and the similarity scorers live
  there. Several functions in it (`pdb_preprocessing_gemmi`, `calculate_pockets`, `pocket_overlap`,
  `download_pisa_info`) are superseded by the class-based modules and are no longer called by `search()`.
- `blosum62.bla` and `human_domains/` (a bundled Foldseek DB) ship as package data — see
  `[tool.setuptools.package-data]` in `pyproject.toml`.
- `tests/e2e/fixtures/` holds the batch input files and a local `4Q5J.cif.gz`. Cases run with that directory
  as their working directory, because `testfile.txt` refers to `4Q5J.cif.gz` by relative path — which is also
  what makes it a local-file-input test. Keep that relative reference if you edit the fixtures.
- `build/` and `dist/` are stale checked-in artifacts containing an older version of the package. Ignore them;
  never edit `build/lib/pocketmapper/`.
