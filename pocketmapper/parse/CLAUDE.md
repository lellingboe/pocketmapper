# parse

Facts about `parse/` that no single docstring states, so they load only when working here. The input
grammar is the two tables in `QTProcessor.__init__`; README "Input format" is the user-facing copy of it.

## Where it sits

- **Step 1, and the only code that turns entry strings into records.** It uses no network and reads no
  structure; it opens a local file only to hash it. It writes `entries.tsv` afresh. Every later step reads
  its records from that file (`lib.entries`) and never parses again, so a change here reaches an existing
  run only when parse and the later steps are rerun.
- **Two layers.** `parse.py` is the step: settings, log, the entries table and the checks that need both
  sides (a Foldseek-DB query is dropped, a DB target must be the only target, each side needs a record).
  `qt_processor.py` is the component. It treats both sides the same way and knows nothing of settings
  or the table.

## What other packages rely on

- **`QTRecord`'s fields are `lib.entries.RECORD_FIELDS`, in the same order.** Neither is derived from the
  other. A new field goes in both places; `test_qt_processor` checks the order.
- **`pockets.pocket_manager.PocketManager.expand_pdb_hits` also calls `QTProcessor`.** It builds one
  without `local_root`, reads `pocket_methods["pisa"][0]` to skip chain pairs the grammar cannot spell, calls
  `parse_individual_qt(..., pocket_method="pisa")`, and replaces `preprocess_name` with the Foldseek hit
  name. A change to that table's layout or to that method's signature breaks the PDB-database path. Only
  the `pdb_fsdb` and `huge` tests (`tests/integration/test_databases.py`) cover that path.
- **`preprocess_name` names the domain chain, not the pocket** (invariant in `project-overview.md`). So
  `4Q5J:B`, `4Q5J:B_F` and `4Q5J:B:10,11` share one name and are aligned once. It also names the files in
  `foldseek_preprocessed_structure_dir`, so a new naming scheme leaves the old cached files unused and
  every chain is preprocessed again.

## Changing the grammar

- **A new pocket method needs:** a row in `pocket_methods`; a place in each `struct_type_pocket_methods`
  tuple that should offer it (tuple order is precedence, and `whole_chain` stays first); a builder in
  `pockets.pocket_manager.POCKET_BUILDERS`; and the hand-written method lists in `cli.py`
  (`--query_pocket_method` help) and the README (Options table, both "Input format" tables).
  `accepted_pocket_methods` reads the table, so validation needs no other change.
- **A residue list is spelled `\d+(,\d+)*,?`.** The older `(\d+,?)+` matches the same strings but
  backtracks exponentially when a long list has a bad last character; keep new patterns free of nested
  quantifiers (`test_residue_pattern_linear`).
- **A top-level input that is an existing file is a file of entries**, unless `is_structure_file` says
  its name is a structure's (`STRUCTURE_EXTENSIONS`). Inside an entries file, any existing file is a
  `local_file` entry.
- **Every pattern allows only single-character chains.** To allow longer chain ids, also check
  `lib.residues.split_chain_info` and the PDB-database skip described in `pockets/CLAUDE.md`.
- `context/settings.md` covers `query_pocket_method` / `target_pocket_method` as settings; read it before
  you change their defaults or how they are layered.

## Tests

`tests/unit/parse/` runs offline and without foldseek: `test_qt_processor.py` covers the grammar tables,
path resolution and record shape; `test_parse.py` covers the entries table and the checks that need both
sides.
