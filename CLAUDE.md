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
tests/e2e/run_e2e.sh --list          # the 16 cases and their tags
tests/e2e/run_e2e.sh -t core         # ~35s, no human_domains searches
tests/e2e/run_e2e.sh -o /tmp/pm_e2e  # everything, results under a chosen directory
tests/e2e/run_e2e.sh test_7          # one case by name
```

Cases are grouped by what they exercise and numbered in that order: 1–7 structure-vs-structure pairs
(all `core`), 8–12 human_domains DB targets, 13–14 the larger Foldseek DB targets, 15–16 the local aligner.
Blank lines separate the groups in the `CASES` heredoc and are skipped by the runner (a `#` comment there
would *not* be). Adding a case means inserting it in its group and renumbering what follows, so prefer
describing a case by what it does rather than pinning to its number.

Each case shells out to the real CLI and hits live wwPDB / AlphaFold / PDBe PISA — there are no mocks, so the
suite needs network access. It asserts exit status plus the presence (and, where a pair is known to produce
hits, the non-emptiness) of `pocket_comparison.tsv`.

**Always run against the existing cache at `tests/e2e/e2e_results/pocketmapper_cache`** — pass
`POCKETMAPPER_E2E_CACHE="$PWD/tests/e2e/e2e_results/pocketmapper_cache"`, or `--cache-dir` for a direct
`pocketmapper search`. It is several GB of already-downloaded structures, PISA responses and Foldseek DBs, and
a cold cache is dramatically slower: PISA is fetched per entry behind a rate-limiting sleep, so a full-PDB run
means thousands of calls at ~3/s. Note the runner's `-o` default is `$PWD/e2e_results`, *not* relative to the
script — running it from the repo root silently creates a second, empty cache at the top level and re-downloads
everything. To deliberately test cold-cache behaviour, point at a throwaway directory *inside*
`tests/e2e/e2e_results/` rather than at /tmp, so the downloads are still reusable.

Cases needing resources that aren't present are reported as SKIP rather than failing: the `needs-pdb-fsdb`
case wants `POCKETMAPPER_PDB_FSDB` pointed at a prebuilt Foldseek PDB database, and `needs-pdb-download`
downloads the full PDB Foldseek DB (2GB download, 7GB unzipped) so it only runs when named explicitly.

Foldseek is an optional external binary (`conda install -c conda-forge -c bioconda foldseek`); it is
invoked via `subprocess.run` and is required for `--foldseek True` and for any `foldseek_db` target.

Whether a case uses Foldseek is decided by `--foldseek` in its own `args` field, not by the runner — cases
tagged `local` omit it to exercise the BLOSUM62 aligner, and only the Foldseek ones are skipped when the
binary is missing. Keep it that way: when every case forced `--foldseek`, the suite could not see the local
branch at all, which is how it shipped broken.

## Pipeline

`main()` → `fire.Fire(PocketMapper())`, so **every public method on `PocketMapper` is a CLI subcommand**
(hence the leading underscores on all internals — that is deliberate, to keep them out of fire's help).
`search()` in `pocketmapper/pocketmapper.py` is the whole workflow and reads top to bottom:

1. `_configure_workflow` → `Settings` (see below), creates directories, dumps `job_settings.json`, configures logging.
2. `_configure_query_target` → `QTProcessor` parses `--query`/`--target` into two DataFrames of `QTRecord`s
   (one `process_qt_cmdline_input` call per side).
3. `_fetch_missing_structures` (or `_fetch_missing_fsdb`) → downloads mmCIF from wwPDB / AlphaFold into `structure_dir`.
4. `_alignment` → either `_foldseek_preprocessing` + `_foldseek_alignment`, or `_local_alignment`. Both write `alignment.tsv`.
5. `_get_pockets` → union of `_retrieve_pisa_pockets` | `_retrieve_passthrough_pockets` | `_retrieve_vdw_pockets`.
6. `_compare_pockets_based_on_alignment` → `pocket_comparison.compare_pockets` → `pocket_comparison.tsv`.
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

A local-file entry like `4Q5J.cif.gz:B_F` resolves to `vdw`, not `pisa` — `B_F` matches the vdw regex and PISA
is PDB-only. That is how the mixed-input fixtures reach the vdw code.

### Pocket dict shape

Every pocket method returns the same nested dict, produced or extended by `lib_struct.parse_pocket_from_struct`:
top-level `res_auth_ids` (list of author seqids as strings), `ca_sequence`, `pocket_exists`, `has_coords`, plus
one entry per residue keyed by the **string** author seqid holding `res_code`, `res_code_single`, `seq_pos`
(0-based index among CA-bearing residues — this is what maps into the alignment) and `ca_coords`.
Residues without a CA atom get `seq_pos = -1` and are excluded, because Foldseek only sees CA-bearing residues.

### Foldseek-DB targets

When the target is a bundled Foldseek DB, `self.fsdb_target` is set and several branches change: no target
structures are fetched or preprocessed, and `_align_structs` reconstructs target PDBs from the DB via
`foldseek createsubdb` + `convert2pdb`. A `foldseek_db` target without `--foldseek True` is a hard error.

What the target "pocket" is then depends on which DB it is, and the two cases are genuinely different:

- **A PDB DB** (`pdb`, or a local prebuilt one) has hits that are real PDB chains, so they get real PISA
  pockets. `_expand_fsdb_pdb_targets` runs first in `_get_pockets`: it reads the hit names out of
  `alignment.tsv`, resolves each through `lib.parse_foldseek_pdb_entry_name`, asks PISA which chains each hit
  chain touches, and appends one ordinary `pisa` record per interface to `_target_df` — so
  `_retrieve_pisa_pockets` then handles them like any other pisa entry and no separate pocket code exists.
  `self._fsdb_pdb_target` is set and `compare_pockets` is called with `alphafold=False`, so every `pocket_2_*`
  column is populated. **Hits with no usable PISA data are dropped**, not compared against a stand-in.
- **Any other DB** (`human_domains`, built from AlphaFold models) keeps the old behaviour:
  `compare_pockets(alphafold=True)` synthesises a whole-chain "pocket" per hit and the `pocket_2_*` columns
  stay empty.

Three things about the PDB path are load-bearing:

**The generated records carry the Foldseek entry name as their `preprocess_name`**, not the one `QTProcessor`
derives. That field is the alignment join key — it is what links these pockets back to their alignment rows
(via `preproc_to_ids`) and to their Foldseek transforms (`foldseek_transform` looks up `u`/`t` by it). The
records are otherwise built by `QTProcessor.parse_individual_qt` on a synthesised `"<PDB>:<chain>_<partner>"`
string, so that method is now part of the pipeline and not only of CLI parsing — its output shape is depended
on in two places.

**The assembly id is deliberately discarded.** `4q5j-assembly1_B` and `4q5j-assembly2_B` both resolve to
`4Q5J:B_F`, so one `pocket_id` can sit behind two `preprocess_name`s. The pocket is computed once, and
`compare_pockets`'s `existing_calcs` set means only the first assembly's alignment row is scored — so the
transform used is whichever assembly Foldseek reported first. `_align_structs` de-duplicates on `pocket_id`
when mapping back to entry names for exactly this reason; without that its `.loc` lookup returns extra rows
and superposes the same structure twice.

**Pockets come from the wwPDB asymmetric unit while Foldseek's `tseq` comes from the assembly.** These agree
in the ordinary case (verified: a 4Q5J self-comparison through a PDB-named DB gives `overlap_count ==
pocket_len`, identity 1.0, RMSD ~1e-14), and `compare_pockets`'s 0.8 sequence-identity guard catches them when
they don't — a populated `incorrect_mapping.json` is the signal that an entry's assembly and AU numbering have
diverged.

**There is no cap on how many hits get enriched**, by choice. A full-PDB search is genuinely large — `4Q5J:B_F`
against the bundled `pdb` DB returns ~4,970 hits across ~3,620 entries, and PISA is fetched per entry with a
rate-limiting sleep, so the first run takes hours. The interface cache makes reruns cheap, and
`_expand_fsdb_pdb_targets` logs both counts before starting so the wait is legible. Add a cap here if that
ever becomes untenable.

## Invariants

The contracts that hold the pipeline together. Breaking one of these generally produces silently wrong
output rather than an error, which is what makes them worth stating. The Foldseek-DB path carries three
more of its own — see that section.

### Identifiers and join keys

**`seq_pos` is the value everything hinges on.** It is the residue's index among the CA-bearing residues *of
its own chain*, and it is what maps a pocket residue into the alignment. Any pocket method that computes it
differently from `lib_struct.parse_pocket_from_struct` will silently produce zero overlap rather than an
error. When adding a pocket method, check a pocket against itself: self-comparison must yield
`overlap_count == pocket_len`.

**`preprocess_name` is the join key.** It is `<basename>_<chain><md5-of-that>` (e.g. `4Q5J_B_<hash>`), computed
once in `QTProcessor.parse_individual_qt`. Alignments are keyed by it; pockets are keyed by `pocket_id`;
`_compare_pockets_based_on_alignment` builds `preproc_to_ids` to bridge the two. One `preprocess_name` can map
to several `pocket_id`s (same chain, different pockets).

**Aligned structures are named by `lib.safe_filename(query_id)`, not by the `pocket_id` itself.** A
`pocket_id` is raw user input, so it may be a path or carry a long residue list; `safe_filename` keeps only
the basename, sanitises it and appends an md5 of the full original. So `aligned_structures/*.pdb` filenames
are not directly greppable for an input string — match on the `MOLECULE` records inside instead.

### Table schemas

**The alignment table's column order is a positional contract.** `SequenceAligner.align_records` builds the
exact same 18 columns that the Foldseek `--format-output` flag requests
(`query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,lddt,qaln,taln,u,t,qseq,tseq`),
and `pocket_comparison.compare_pockets` reads them **by index** (`row[0]`…`row[17]`). Changing or reordering
columns in one producer without the other, or without updating the indices in `compare_pockets`, breaks
silently.

**The comparison table has a fixed schema.** `pocket_comparison.POCKET_COMPARISON_COLUMNS` declares every
column `compare_pockets` can produce, and the result is reindexed onto it before returning, so all 33 exist
on every run — rows that stop early (no overlap, no coordinates, an alphafold/foldseek-db target with no
pocket 2) leave the later fields empty rather than dropping them. Add a new output field to that list as well
as to the row dict; a column produced but not declared is kept and logged as a warning rather than silently
dropped, so the list cannot quietly drift.

### Structural alignment (step 7)

**`_align_structs` only superposes targets that actually overlap the query.** It filters on
`overlap_count > 0` before ranking, because a target sharing no pocket residues has no common residue set to
superpose on and empty overlap metrics that would sort arbitrarily. Queries with no overlapping target are
skipped with a log line, and a run where nothing overlaps returns early instead of proceeding.

**It does not work at all on the local-aligner path.** `SequenceAligner` writes `"-"` for the `u`/`t`
transforms, which `foldseek_transform` catches and logs per record, so a non-Foldseek run still completes but
its aligned PDB contains only the query.

## Using it as a library

Nothing in the package depends on being launched from a terminal — fire only wraps `PocketMapper` at the
`main()` boundary. There are two levels of entry, coarse and fine:

```python
# 1. Whole pipeline. Same work as the CLI; writes the same files into results_dir.
from pocketmapper.pocketmapper import PocketMapper
PocketMapper().search(query="4Q5J:B_F", target="4Q5J:A_E", results_dir="./out")

# 2. Individual components. None of them take a Settings -- they are driven by explicit values.
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
functions. **`Settings` stays inside `pocketmapper.py`** — it is resolved there and unpacked at each call
site into the individual values a component needs. Keep it that way when adding components: a component
that reaches into a `Settings` cannot be used without building one, and hides which fields it depends on.

Things to know when consuming it this way:

- **`pocketmapper/__init__.py` only exports `main` and `__version__`.** Submodules are reachable as
  `pocketmapper.lib` etc. only as a side effect of `pocketmapper.pocketmapper` importing them — always use
  explicit `from pocketmapper.<module> import <name>` rather than relying on that.
- **Always call `Settings(...).resolve_paths()`** if you build a `Settings` yourself to source paths from.
  Unresolved derived paths are `None`, and passing one on gives an opaque `TypeError: expected str, bytes
  or os.PathLike object, not NoneType` from deep inside `os.path.join`, not a helpful error. `search()`
  does this for you.
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
  The one remaining `exit()` is `_check_help_search`, which prints `HELP_MESSAGE` and quits when `--help` is
  passed; it is a CLI affordance, so `PocketMapper().search(help=True)` will kill a host process — a library
  caller should never pass it.
- **Settings.** `Settings` is a frozen-by-convention dataclass updated via `dataclasses.replace`. Resolution order is
  defaults → `--settings` JSON file → explicit CLI args → `resolve_paths()` (which only fills derived paths still `None`,
  so a settings file can pin any individual path). New options go on the dataclass *and* in the `cli_overrides` dict
  in `_configure_workflow` *and* in the `search()` signature — and, if they are user-facing, in `HELP_MESSAGE`
  (`constants.py`) and the README's Options list, neither of which is generated from the dataclass.
- **Fetcher/preprocessor API.** `StructureFetcher` and `StructurePreprocessor` follow
  `set_output_directory()` → `update_cache()` → `fetch_*`/`preprocess_records()`. Caching is a plain `os.listdir`
  snapshot, so `update_cache()` must be called after the output dir is set and before work begins.
- **Structure parsing** is gemmi throughout (`.cif.gz` on disk); Biopython is used only for pairwise alignment
  and SVD superposition.

## Repo layout notes

- `lib.py` holds only generic, stateless helpers — `jsonify_dict`, `safe_filename`, the BLOSUM62 matrix
  reader and the similarity scorers. Nothing in it knows about `Settings`, the pipeline or the pocket dict
  shape; keep it that way, and put workflow logic in a component module instead. It used to be a grab-bag:
  the superseded copies of the preprocessing, pocket-calculation and PISA-download logic were deleted in
  favour of the class-based modules, and `compare_pockets` moved out to `pocket_comparison.py`.
- `PocketCalculator.atp_pocket_overlap` is uncalled but deliberately retained for planned ATP-pocket work —
  leave it in place rather than pruning it as dead code.
- `pocket_comparison.py` owns pipeline step 6 — `POCKET_COMPARISON_COLUMNS`, `compare_pockets`, and the
  Foldseek column-order contract they depend on.
- `constants.py` holds `SINGLE_AA_CODE` (the one three-to-one letter table — duplicates elsewhere were
  deleted; it maps the modified residues `SEP`/`TPO`/`PTR`/`MSE` and every caller defaults unknowns to `"X"`)
  and `HELP_MESSAGE`, the `--help` text, which also documents the settings-file-only "Advanced Options".
- `blosum62.bla` and `human_domains/` (a bundled Foldseek DB) ship as package data — see
  `[tool.setuptools.package-data]` in `pyproject.toml`.
- `tests/e2e/fixtures/` holds the batch input files and a local `4Q5J.cif.gz`. Cases run with that directory
  as their working directory, because `testfile.txt` refers to `4Q5J.cif.gz` by relative path — which is also
  what makes it a local-file-input test. Keep that relative reference if you edit the fixtures.
- `build/` and `dist/` are stale checked-in artifacts containing an older version of the package. Ignore them;
  never edit `build/lib/pocketmapper/`.
