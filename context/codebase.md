# Codebase review

A review of structure, layout and consistency as of commit `b13dfd1` (2026-10-09). It covers every module
in `pocketmapper/` (~8.9k lines), the config and CI files, and `context/`. Line numbers are against
that commit.

Each finding carries a source tag:
- **[probed]**: reproduced read-only in `pocketmapper_dev`.
- **[grep]**: confirmed by search.
- **[read]**: found by reading the code.

## Summary

The pipeline is mature and heavily documented. The steps follow a clear two-layer design, and
`search` equals the chained steps by construction. The main debts are:
- per-step plumbing copied six times
- the CA-walk (`seq_pos`) invariant implemented in four places
- stateless classes beside plain functions, with no rule for which to use
- settings validated in several places
- three crashes on plausible input
- some drift between the docs and the code

## Layout

```
pocketmapper/
  cli.py                  argparse front end; the only argv/exit-code module
  pocketmapper.py         PocketMapper.search: builds Settings, calls every step
  settings.py             Settings, job-file layering, path resolution, resolve_* validators
  constants.py            schemas, defaults, value sets, CLI epilogs
  lib.py                  "generic" helpers (logging, temp dirs, scopes) + domain helpers
  entries.py              JSON I/O, failed_entries.json, QTRecord helpers
  exceptions.py           PocketMapperError
  qt_processor.py         input grammar -> QTRecord
  foldseek.py             binary wrapper, bundled DB paths
  sequence_aligner.py     --aligner seq
  structure_preprocessor.py  single-chain cut for Foldseek
  structure_aligner.py    superposition + PDB writing
  pocket_comparison.py    compare_pockets, POCKET_COMPARISON_COLUMNS
  steps/                  parse, fetch_structures, align, pockets, compare, superpose
  pockets/                Pocket shape, fetcher/dispatch, structure/pisa/vdw builders
  downloads/              lib_download, structure_downloader, pisa_downloader
  human_domains/          bundled Foldseek DB + offset table
  blosum62.bla
```

- **Uneven grouping.** `pockets/` and `downloads/` became subpackages. The alignment and
  superposition components (`sequence_aligner`, `structure_preprocessor`, `structure_aligner`,
  `foldseek`, `pocket_comparison`) are still flat. This does no harm, but there is no stated rule.
- **Steps import each other.**
  - `align`, `pockets` and `superpose` import `fetch_structures.fetch_missing_*`.
  - Every step imports `parse.parse_job_entries` or one of its siblings.
  - So the shared helpers live inside step modules rather than in a shared module.
- **`entries.py` holds three concerns:**
  - generic JSON I/O (`require_file`, `read_json`, `write_json`)
  - the failed-entries log
  - QTRecord helpers

  The name covers only the last two.
- **`lib.py` is not generic.** Its docstring says "Nothing here knows about the pipeline"
  (`lib.py:4`). Yet `lib.py:355-570` holds:
  - BLOSUM scoring
  - Foldseek PDB entry-name parsing
  - `fsdb_pocket_mode`
  - `split_chain_info`
  - Foldseek one-letter codes
  - UniProt offset maps

## Strengths

- **Two-layer steps.** The entry function does settings, logging and temp. The core function takes
  explicit values. `search` calls the same entry functions.
- **Errors.** There is one error type, and `cli()` is the only exit point.
- **Logging discipline.** Every logger sits under the package logger, which has a `NullHandler`. The
  handlers are scoped and removed in `finally`.
- **Declared contracts.** Table schemas are declared, the pockets file is versioned, and an
  undeclared column is flagged, not dropped.
- **Safety guards.** `temp_dir` emptying and deletion are guarded by `is_within`. Scopes nest and are
  released in `finally`.
- **Downloads.** There are two entry points by workload, `.part` writes, and a whitelist of retryable
  statuses.
- **CI.** The `compat` matrix guards the 3.10 floor separately from lint.
- **Documentation.** Docstrings are thorough and mostly accurate. Non-obvious decisions are
  explained at the site.

## Structural duplication

- **Step entry boilerplate, ×6** (`steps/*.py`). Every entry function does the same sequence:
  1. take a signature of ~20 parameters
  2. build an `arguments` dict
  3. loop over `require_setting`
  4. call `resolve_paths`
  5. open `run_scope` + `log_to_file`
  6. call `dump_settings` if it is the outermost scope

  The Args entries for the seven cache and path parameters are copied into six steps and `search`.
  - Suggestion: one context manager, e.g. `settings.step_scope(command, job_file, arguments,
    required)`, yielding `values`.
  - Keep the signatures explicit, because the CLI parity check relies on them.
- **CA-residue walk, ×4.** This is the `seq_pos` invariant, the costliest one to get wrong. Each
  site walks the chain differently:

  | Site | Walk | CA test |
  |---|---|---|
  | `pockets/structure.py:72` | `chain.first_conformer()` | `get_ca()`, which also requires element C |
  | `pockets/vdw.py:60` | `get_polymer().first_conformer()` | `"CA" in res`, by name |
  | `sequence_aligner.py:91` | `get_polymer().first_conformer()` | `"CA" in res`, by name |
  | `pockets/vdw.py:136` (`atp_pocket_overlap`, kept deliberately) | `get_polymer()`, no `first_conformer` | `"CA" in res` |

  The last site breaks the documented "every chain walk reads `first_conformer()`" invariant.
  Suggestion: one `ca_residues(chain)` helper used by all four.
- **Near-identical parse helpers.** `parse_side` and `parse_structure_side` (`steps/parse.py:357-426`)
  differ only in the `QTProcessor` method they call.
- **Three TSV writers of the same shape:** parse's `entries_path`, `write_structures_table` and
  `write_pockets_table`. Each does `make_dir` → DataFrame → `to_csv` → log.
- **Re-implemented dedup.** The `PocketFetcher` dedup (`pockets/pocket_fetcher.py:191-197`) is
  `entries.unique_by` over again.
- **Pairwise-identical validators** in `settings.py`:
  - `resolve_aligner` and `resolve_pisa_source`: a choice validator.
  - `resolve_delete_tmp` and `resolve_fetch_missing`: a 0/1 flag validator.
- **Log-and-raise, ×45.** `logger.critical(msg)` followed by `raise PocketMapperError(msg)`. Some
  pairs already disagree (`qt_processor.py:160-161`). Suggestion: a `fail(msg, log_extra)` helper
  that logs and returns the exception.
- **Duplicated "foldseek required" message**, verbatim in `pocketmapper.py:282-285` and
  `steps/align.py:141-144`.

## Consistency

- **Classes vs functions.** These classes hold nothing but `log_extra`:
  - `SequenceAligner`
  - `StructurePreprocessor`
  - `StructureAligner`
  - `PocketFetcher`
  - `PisaParser`
  - `VdWCalculator`

  `PocketMapper` holds nothing at all. Steps and pocket builders are plain functions. Pick one rule.
  Converting the classes would change the library API.
- **Validation is scattered.**
  - **Pocket methods:** the accepted values live on a `QTProcessor` instance
    (`accepted_pocket_methods`). Every other value set is in `constants`. The CLI help hardcodes the
    same list.
  - **Grammar regexes:** these are instance attributes. `steps/pockets.py:359` reaches into
    `qtprocessor.pocket_methods["pisa"][0]`.
  - **`pisa_source`:** validated in `settings` and again in `PisaDownloader.__init__`.
  - **`align_struct_method`:** validated in `settings` and again in `StructureAligner.align_structs`.
- **Log stages.** There are 27 distinct labels, with mixed tense and case:
  - "Determine Query/Target Types"
  - "Processing Inputs"
  - "Retrieving passthrough Pockets" next to "Retrieving PISA Pockets"
  - "End"

  Few of them are step names. **[grep]**
- **Progress bars.** Twelve `tqdm` bars write to stderr whatever the verbosity. That breaks the
  logging contract, under which a library caller sees nothing unless it attaches a handler. **[grep]**
- **JSON writing has four styles:**
  - `entries.write_json`: indent 4, makes the dir.
  - `settings.dump_settings`: indent 4, uses `os.makedirs`.
  - compare's reports (`steps/compare.py:206,216`): raw `json.dump`, no indent.
  - `pocket_fetcher`: compact, and the caller makes the dir.
- **Resource lookup has two methods.**
  - `BLOSUM_PATH` is built from `__file__` (`steps/compare.py:32`).
  - `human_domains` uses `importlib.resources` (`foldseek.py:145`).
- **Step numbering has three schemes.**
  - `pocketmapper.py:8-17` numbers parse as 2.
  - The `steps/*` docstrings call parse "Step 1".
  - `resolve_settings` carries `# 4b`–`# 5` from an older scheme, with no 1–4a
    (`pocketmapper.py:277-304`).
- **Naming.**
  - `QT`/`qt` throughout. `process_qt_cmdline_input` also serves job files and library calls.
  - `p1`/`p2` in `pocket_comparison` against query/target everywhere else.
  - `pocket_comparison.superpose` (a function) vs `steps.superpose` (a step).
  - `pocketmapper.pocketmapper`.
  - camelCase `replaceNonCommonResidues` (`sequence_aligner.py:26`).
- **Deviations from `coding-standards.md`.**
  - `typing.NamedTuple` annotations on `MappedPocket` and `Context` (`pocket_comparison.py:80,492`).
    The same module also uses `collections.namedtuple` (`AlignmentRow`).
  - `open(...).read()` without `with` and an unused `max_score` (`lib.py:411-443`).
  - `isinstance(qt_input, type(None))` (`qt_processor.py:159`).
  - A stale "Exiting." in a message that does not exit (`qt_processor.py:160`).
  - `raise` without `from error` (`pocket_comparison.py:241`).
- **Docstring wrap glitches.**
  - An overlong Raises line in `steps/fetch_structures.py:111`.
  - A broken Raises paragraph in `steps/pockets.py:115-118`.
- **CLI-only text in `constants.py`.** It holds the CLI epilogs, but its docstring scopes it to
  tables that several modules must agree on.

## Defects

1. **A blank line in an entries file aborts the run.** **[probed]**
   - `read_entries` keeps `""` as an entry (`qt_processor.py:240`).
   - `determine_struct_type("")` then resolves the path to `work_dir`, finds a directory, and raises
     `Directory input is not currently supported: ` (`qt_processor.py:442-444`).
   - The run fails with a misleading message instead of skipping the line.
2. **A vdw pair with no contacts crashes compare.** **[probed]**
   - `pocket_overlap` sets `pocket_exists = True` even with zero residues (`pockets/vdw.py:99`).
     4Q5J `A_F`, `B_E` and `E_F` each give 0.
   - `compare_pocket_pair` then divides by a `query_len` of 0 (`pocket_comparison.py:457`) and raises
     `ZeroDivisionError`. That is not a `PocketMapperError`, so the CLI prints a traceback.
   - pisa and `parse_pocket_from_struct` set `pocket_exists` only when residues are found.
3. **A missing chain raises a bare `ValueError`.** **[probed for gemmi; path found by reading]**
   - gemmi raises `ValueError: Chain Z not found` from `model[name]`.
   - Two sites index this way and don't catch it: `pockets/vdw.py:60` and `sequence_aligner.py:91`.
   - Parse never checks that a chain exists, so `--aligner seq` or a vdw entry naming an absent chain
     ends in a traceback.
   - The Foldseek preprocessing path and `parse_pocket_from_struct` handle the same case and skip the
     entry.
4. **The `mismatch` column of `alignment.tsv` means different things by aligner.** **[read]**
   The seq aligner writes a fraction (`sequence_aligner.py:131`); Foldseek writes a count. Nothing
   downstream reads it, but the column is user-facing.
5. **The PISA interface cache is written without `.part`.** **[read]**
   - `parse_assemblies` writes `<pdb_code>.json` directly (`downloads/pisa_downloader.py:425-428`).
   - `download_missing_interfaces` trusts any `*.json` it finds (line 132).
   - An interrupted write leaves a truncated file that later runs trust, against the `.part` guarantee
     stated in the module and in `downloads/CLAUDE.md`.
6. **The per-method pocket "cache" is write-only.** **[grep]**
   - `fetch_pockets` writes `<pocket_dir>/<method>_pockets.json` (`pocket_fetcher.py:208`), but
     `load_pockets` has no caller.
   - Pockets are rebuilt every run; only the PISA downloads under `pocket_dir/pisa/` are cached.

## Doc drift

- **`project-overview.md` "Declared schemas"** places `POCKET_COMPARISON_COLUMNS` in `constants`.
  It is defined in `pocket_comparison.py:46`.
- **`project-overview.md` "Logging and errors"** says "Declared stages are Title Case step names".
  That is not true (see Consistency).
- **`--pocket_dir` help and `Settings` docs** call it a "cache of parsed pockets". See defect 6.
- **The `lib.py` module docstring** claims to be pipeline-agnostic (see Layout).

## Testing

- There are no unit tests. The e2e suite (71 cases) hits live wwPDB and PDBe.
- The parity checks in `project-overview.md` "Settings" are run by hand:
  - each step's signature == its subparser's dests
  - `Settings` fields == `search`'s signature == its `arguments` keys
  - `SETTING_DEFAULTS` ⊆ `JOB_KEYS`
- These are cheap, offline pytest targets:
  - the parity checks above
  - the `QTProcessor` grammar (inference vs forced methods, rejections, blank lines)
  - `settings.layer_settings` / `resolve_paths`
  - `entries.unlisted` (dedup)
  - `compare_pockets` on synthetic rows, which is how defect 2 was reproduced

## Hand-maintained sync points

- A new option touches 4–5 places: signature, `arguments`, `OPTIONS`/`COMMANDS`, README and
  `SETTING_DEFAULTS`.
- The supported Python range is stated in five places.
- The README Options section duplicates the CLI help strings.

All three would be covered by the parity test above.

## Recommendations (in priority order)

1. **Fix defects 1–3:**
   - skip blank entry lines
   - set vdw `pocket_exists` from the residues it found
   - check that the chain exists at parse time, or catch it in the builders and the seq aligner and
     report the entry as failed
2. **Unify the CA walk** into one helper.
3. **Fix defect 5**: write the interface file through `.part`. **Fix defect 6**: either read the
   per-method files or stop writing them, and correct the docs.
4. **Add `step_scope`** to remove the step-entry boilerplate.
5. **Make the grammar and validation single-source:**
   - move the pocket-method values and grammar regexes to module level or `constants`
   - validate each setting in one place
6. **Add a small offline pytest** for the parity checks and the grammar.
7. **Settle classes vs functions** for stateless components.
8. **Housekeeping:**
   - normalise the stage names
   - gate `tqdm` on verbosity (`disable=`)
   - unify JSON writing and resource lookup
   - fix the naming and docstring nits above
9. **Correct the doc drift** listed above.
