# parse

Cross-module facts about `parse/`.

## Where it is stated

- The grammar: the `qt_processor` module docstring and the tables in `QTProcessor.__init__`; `docs/input.md`
  "Input format" is the user-facing copy.
- The step and its component: the `parse.parse` and `qt_processor` module docstrings.
- Entries files vs single entries: `QTProcessor.read_entries`.
- The residue-list pattern: its comment in `QTProcessor.__init__` and `test_residue_pattern_linear`; keep new
  patterns free of nested quantifiers.
- Settings (`query_pocket_method`, `target_pocket_method`): `context/settings.md`.

## What other packages rely on

- **Parse is the only step that parses**, and the only one that opens local inputs other than structures: it
  reads entries files and hashes local structure files. Every later step reads its records from
  `entries.tsv`, so a change here reaches an existing run only when parse and the later steps are rerun.
- **`QTRecord`'s fields are `lib.entries.RECORD_FIELDS`, in the same order.** Neither is derived from the
  other (`lib` cannot import `parse`). A new field goes in both places; `test_qt_processor` checks the order.
- **`PocketManager.expand_pdb_hits` also calls `QTProcessor`.** It builds one without `local_root`, reads
  `pocket_methods["pisa"][0]` to skip chain pairs the grammar cannot spell, calls
  `parse_individual_qt(..., pocket_method="pisa")`, and replaces `preprocess_name` with the hit name. A
  change to that table's layout or to that method's signature breaks the PDB-database path
  (`test_pocket_manager.py` offline; `pdb_fsdb` and `huge` end to end).
- **`preprocess_name` names the domain chain, not the pocket** (invariant in `project-overview.md`). So
  `4Q5J:B`, `4Q5J:B_F` and `4Q5J:B:10,11` share one name and are aligned once. It also names the files in
  `foldseek_preprocessed_structure_dir`, so a new naming scheme leaves the old cached files unused. Its hash
  covers the local file parse resolved under `local_root`; a later step run with another `local_root` can
  load a different file of the same name, undetected.

## A new pocket method needs

- a row in `pocket_methods`;
- a place in each `struct_type_pocket_methods` tuple that should offer it (tuple order is precedence, and
  `whole_chain` stays first);
- a builder in `pockets.pocket_manager.POCKET_BUILDERS` (and an entry in `PocketManager`'s `builder_options`
  if it takes options);
- the hand-written method lists: the `Option.help` of `query_pocket_method` in `lib/options.py`
  (`test_pocket_method_help_lists_accepted_methods` checks it), the `PocketMapper.search` docstring,
  `docs/options.md` (Options table) and `docs/input.md` (both "Input format" tables).

`accepted_pocket_methods` reads the table, so validation needs no other change.

Longer chain ids: see `pockets/CLAUDE.md` "Single-character chains".
