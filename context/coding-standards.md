# Coding standards

A rule earns a line here only when breaking it fails silently rather than loudly. What a docstring
or comment at the code site already says appears as a pointer to that site, never as a second copy.
Generic conventions first, then what is specific to this pipeline.

## Generic

- **Comments explain why, and what breaks otherwise.** Density is deliberately high — see
  `constants.py`, where nearly every constant carries a rationale above it. A comment that restates
  the code is noise; one recording a measurement, a deliberate omission or a load-bearing ordering
  earns its place.
- **Docstrings are Google-style** — `Args:` / `Returns:` / `Raises:`. Match that rather than the
  reST outlier in `structure_preprocessor.py` or the section-less prose in `lib.py`, and keep
  `Returns:` honest — several already claim `None` and return a value.
- **No type annotations on functions.** No function in the package carries them. Annotate
  dataclass and `NamedTuple` fields only, in `str | None` form (`Settings`, `QTRecord`).
- **`os.path`, never `pathlib`; f-strings, never `%` or `.format()`** — inside log calls too.
  Currently no exceptions to either.
- **Reuse before adding.** `lib.py` and `lib_struct.py` already hold the generic helpers.
- **Formatting is not a standard here.** black and flake8 run on pre-commit; their settings live in
  `pyproject.toml` and `.flake8`.

## Specific

### Logging

Every log call must pass `extra={"stage": "..."}` — usually a local `stage` dict or
`self._log_extra` — or the record fails to format against the root formatter
`"%(levelname)s: %(stage)s - %(msg)s"` (set in `PocketMapper.__init__`). Every call currently does;
nothing enforces it.

### Errors

`logging.critical(...)` then `raise PocketMapperError(...)` (from `pocketmapper.exceptions`);
`main()` catches it and exits 1. Raising rather than exiting is what makes the package embeddable —
preserve it when adding error paths.

**No `exit()`/`sys.exit()` inside modules** — deliberately removed, which no code comment can show.
The one survivor is the CLI affordance `_check_help_search`, so `PocketMapper().search(help=True)`
will kill a host process — library callers must not pass it.

### Settings

A new option goes in **five hand-maintained places**: the `Settings` dataclass, the `cli_overrides`
dict in `_configure_workflow`, the `search()` signature, `HELP_MESSAGE` (`constants.py`) and the
README's Options list. Neither of the last two is generated from the dataclass. Miss one and the
option is silently ignored.

Resolution order, and the tri-state `foldseek` and `align_struct_method` settings, are documented
where they are resolved — the `Settings` docstring and the `# 4b.` / `# 4c.` comments in
`_configure_workflow`, which also give the reasons those call sites are load-bearing. Keep new
resolution logic there, not at the call sites consuming the resolved value.

### Components

- `StructureFetcher` and `StructurePreprocessor` follow `set_output_directory()` → `update_cache()`
  → `fetch_*`/`preprocess_records()`. The ordering is required and stated nowhere in the code.
- **That cache never hits.** `update_cache()` stores `os.listdir` basenames, but `fetch_alphafold`,
  `fetch_mmcif` and `preprocess_records` test a full output path for membership in it, so every
  structure is refetched and re-split on every run. Fix the comparison rather than documenting
  around it if this becomes the bottleneck.
- **Structure parsing is gemmi throughout** (`.cif.gz` on disk). Biopython is used only for
  pairwise alignment (`sequence_aligner.py`) and SVD superposition (`pocket_comparison.py`).
- Components never take a `Settings` — see "As a library".

### Extending the pipeline

Each of these fails silently rather than raising. Most are documented at the code site; the list is
here so you know to go and look.

- **A new pocket method** must compute `seq_pos` exactly as `lib_struct.parse_pocket_from_struct`
  does — the index among the CA-bearing residues of its own chain. A different convention yields
  zero overlap with no error. Check it by comparing a pocket against itself: `overlap_count ==
  pocket_len`.
- **A new alignment column** goes into `constants.ALIGNMENT_COLUMNS`, never one producer alone — see
  the note above the constant.
- **A new comparison output field** goes into `pocket_comparison.POCKET_COMPARISON_COLUMNS` as well
  as the row dict; an undeclared column is kept and warned about, not dropped.
- **`compare_pockets` must not write to a pocket dict** — return a `_MappedPocket` instead
  (`pocket_comparison.py` module docstring).
- **Never hand a raw `p2_to_p1_u`/`p2_to_p1_t` cell to gemmi** — go through
  `pocket_comparison.parse_pocket_transform`, whose docstring carries the measured evidence.
- **Changing step 6 without changing behaviour**: there are no unit tests, so capture
  `compare_pockets`' arguments from a real run and diff old against new output.

### As a library

fire wraps `PocketMapper` only at the `main()` boundary, so nothing here needs a terminal. Two
levels of entry:

```python
# 1. Whole pipeline. Same work as the CLI; writes the same files into results_dir.
from pocketmapper.pocketmapper import PocketMapper
PocketMapper().search(query="4Q5J:B_F", target="4Q5J:A_E", results_dir="./out")

# 2. Any component -- driven by explicit values, never a Settings. QTProcessor parses one side per
#    call; `name` labels that side in logs and errors.
from pocketmapper.qt_processor import QTProcessor
qtp = QTProcessor(
    structure_dir="./structs", foldseek_preprocessed_structure_dir="./preproc", fsdb_dir="./fsdb"
)
query_df = qtp.process_qt_cmdline_input("4Q5J:B_F", name="query")
```

- **`Settings` stays inside `pocketmapper.py`**, unpacked at each call site into the values a
  component needs. A component reaching into a `Settings` can't be used without building one, and
  hides which fields it depends on.
- **`pocketmapper/__init__.py` only exports `main` and `__version__`.** Submodules are reachable as
  `pocketmapper.lib` etc. only as a side effect of `pocketmapper.pocketmapper` importing them —
  always use explicit `from pocketmapper.<module> import <name>`.
- **Always call `Settings(...).resolve_paths()`** if you build one yourself, or an unset derived
  path gives an opaque `TypeError: expected str, bytes or os.PathLike object, not NoneType` from
  deep inside `os.path.join`. `search()` does this for you.
- **`search()` has global side effects**: `logging.config.dictConfig` reconfigures the *root* logger
  and stomps on a host app's logging setup, and `_delete_tmp` `shutil.rmtree`s
  `query_dir`/`target_dir`/`foldseek_tmp_dir` at the end (marked `# TODO this is unsafe`).
- Results come back through files — `search()` returns `None`, so read `pocket_comparison.tsv` /
  `alignment.tsv` from `results_dir` (paths available on `Settings`).
