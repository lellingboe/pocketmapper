# Coding standards

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
- **The `align_struct_method` setting is tri-value**, and `_resolve_align_struct_method` collapses `"auto"` to
  `"pocket"` or `"foldseek"` so nothing downstream sees it. Same load-bearing call site as `_resolve_foldseek`, and
  immediately after it: it reads the resolved `foldseek` bool, needs logging configured, and must precede the
  settings dump so `job_settings.json` records the resolved value. Both of its errors are pure config, so an
  impossible combination fails before anything is downloaded.
- **Fetcher/preprocessor API.** `StructureFetcher` and `StructurePreprocessor` follow `set_output_directory()`
  → `update_cache()` → `fetch_*`/`preprocess_records()`. Caching is a plain `os.listdir` snapshot, so
  `update_cache()` must be called after the output dir is set and before work begins.
- **Structure parsing** is gemmi throughout (`.cif.gz` on disk); Biopython only for pairwise alignment and SVD
  superposition.

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
