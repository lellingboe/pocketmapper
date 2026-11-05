# Using PocketMapper as a library

`PocketMapper().search(...)` runs the same workflow as the CLI and writes the same files; the result
tables are written under `results_dir`, and `search()` returns the resolved settings. Each step command
is a function, `pocketmapper.<step>.<step>.<command>` (such as `pocketmapper.align.align.align`;
`fetch_structures` is `pocketmapper.fetch.fetch.fetch_structures`), taking its CLI options as keyword
arguments; its `job_file` also takes a dict, such as the settings `search()` returns. The individual
components (`parse.qt_processor`, `lib.structure_manager`, `pockets.pisa_downloader`,
`align.sequence_aligner`, `superpose.structure_aligner`, `pockets.pocket_manager`, ...) are each usable
on their own. Note that `search()` deletes its temporary directories on the way out.

The options are those of the CLI; see the [command-line reference](options.md).

## Logging

Logging goes through the `pocketmapper` logger and never touches the root logger, so your own logging
setup decides what is shown. `search()` still writes `info.log` and sets the `pocketmapper` logger's level
from `verbosity` for the length of the call. For the CLI's console format, use
`pocketmapper.lib.logs.format_handler(logging.StreamHandler())` and add the handler to
`logging.getLogger("pocketmapper")`.

## Deferring superposition

The last step is the one you can defer: run with `align_count=0` and no aligned structures are written,
then run `pocketmapper.superpose.superpose.superpose(job_file=settings)`. For finer control, use
a `StructureAligner` built on a `StructureManager` (which finds each structure from the run's cache dirs
and `local_root`) and hand its `align_structs` the run's records, read from its `entries.tsv` with
`pocketmapper.lib.entries.entry_sides(pocketmapper.lib.entries.read_entries(path))`, and its result
files: it takes `query_ids` and `target_ids` to superpose only part of the run, and `overwrite=False` to
leave alone what it has already written. (Against a Foldseek database target it also needs `fsdb_path`,
`manager.record_path(pocketmapper.lib.entries.fsdb_record(sides["target"]))`; `superpose.superpose`
handles that.)

```python
from pocketmapper import PocketMapper
from pocketmapper.lib.entries import entry_sides, read_entries
from pocketmapper.lib.structure_manager import StructureManager
from pocketmapper.superpose.structure_aligner import StructureAligner

settings = PocketMapper().search("4Q5J:A_E", "4Q5J:B_F", align_count=0, results_dir="results")
sides = entry_sides(read_entries(settings["entries_path"]))
manager = StructureManager(
    pdb_dir=settings["pdb_dir"],
    alphafold_dir=settings["alphafold_dir"],
    fsdb_dir=settings["fsdb_dir"],
    local_root=settings["local_root"],
)

StructureAligner(manager).align_structs(
    query_records=sides["query"],
    target_records=sides["target"],
    pocket_comparison=settings["pocket_comparison_path"],
    alignment=settings["alignment_path"],
    out_dir=settings["aligned_structure_dir"],
    method=settings["align_struct_method"],
    align_count=10,
    query_ids=["4Q5J:A_E"],
)
```
