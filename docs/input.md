# Entries and targets

How to write a query or target entry, how PocketMapper decides where each pocket comes from, and what
searching a Foldseek database gives you.

## Input format

Query and target entries are colon-separated:

```text
STRUCTURE[:CHAIN[:RESIDUES]]
```

- **STRUCTURE** — a 4-character PDB ID (`4Q5J`), a UniProt accession (`P04637`, fetched from
  AlphaFold), or a path to a local structure file. IDs and accessions are case-insensitive.
- **CHAIN** — the chain the pocket belongs to. For the interface-based methods this is two chains joined
  by an underscore (`B_F`), where the first is the chain carrying the pocket and the second is its
  binding partner. Optional; omitted, it defaults to chain `A`.
- **RESIDUES** — optional comma-separated author residue ids (`10,11,12`). An id may carry a one-letter
  insertion code (`100,100A,100B`), which is case-sensitive; `100`, `100A` and `100B` are three
  residues. Omitting the whole pocket part makes the entry an open search — see below.

How the pocket residues are derived is inferred from the shape of the entry:

| Example | Method | Meaning |
| --- | --- | --- |
| `4Q5J:B_F` | pisa | Interface residues of chain B against chain F, taken from PDBe PISA, or, with `--pisa_bsa_threshold`, residues buried past it |
| `4Q5J:A:10,11,12` | passthrough | Exactly the residues listed, on chain A (reported in ascending order) |
| `4Q5J:A_B:10,11,12` | vdw | Van der Waals contact residues on chain A against chain B. The residue list is ignored: on a PDB entry the trailing `:` part only selects vdw over pisa (or force it with `--query_pocket_method vdw`) |
| `P04637:A:10,11,12` | passthrough | AlphaFold model for the accession, residues listed |
| `4Q5J:B` | whole_chain | Open search — the whole of chain B is treated as the pocket |
| `4Q5J` | whole_chain | Open search over the default chain, `A` |
| `human_domains` | — | Bundled Foldseek database (valid as a target only) |

Which methods an entry can reach depends on its structure type:

| Structure type | Methods available |
| --- | --- |
| PDB entry | `pisa`, `passthrough`, `vdw`, `whole_chain` |
| Local file | `passthrough`, `vdw`, `whole_chain` |
| AlphaFold accession | `passthrough`, `whole_chain` |

PISA is only available for PDB entries, so a local file with an interface-style chain spec
(`4Q5J.cif.gz:B_F`) resolves to `vdw`, not `pisa`. An AlphaFold model is a single chain, so the
methods that need a binding partner are not offered for one at all. A passthrough entry needs an
explicit residue list — without one the entry is an open search instead.

The inferred method (`auto`, the default) can be overridden with `--query_pocket_method` /
`--target_pocket_method`, which does not relax any of the above: the entry must still spell what the
method reads, and the method must be available for that structure type. An entry that cannot satisfy
both is skipped with a warning naming what was missing, before anything is downloaded, and the run
fails only if it leaves a side with no usable entries — so one bad line in a batch file does not stop
the rest. An unrecognised method name is rejected outright.

A pocket's residue ids are checked against the structure, whatever the pocket method: a residue the chain
does not have, or has without a CA atom, is dropped from the pocket, with one warning per pocket listing
them. An entry left with none of its residues is skipped (`pocket_not_built`) rather than compared, and the
rest of the run continues. A residue id listed more than once is collapsed to one, also with a warning.

### Batch files

For batch runs, pass a path to a file containing one such entry per line instead of a single entry;
blank lines are ignored. A path ending in `.cif`, `.mmcif`, `.pdb` or `.ent` (optionally `.gz`) is read
as a structure entry, not as a file of entries. Query and target files are read independently, and every
query is compared against every target.

## Closed vs open search

The two shapes of entry ask different questions.

A **closed search** names a pocket on both sides (`4Q5J:B_F` against `4Q5J:A_E`) and asks *how do these two
pockets compare?* Every column of `pocket_comparison.tsv` is populated: both pockets are described, the
overlap is scored, and `jaccard_index` sizes the shared residues against their union.

An **open search** names a structure but no pocket (`4Q5J:B`, or just `4Q5J`) and asks *does the query
pocket resemble anything on this chain at all?* The whole chain becomes the pocket, so the answer is
carried by `overlap_count` and `query_overlap_ids` — how many of the query pocket's residues the target
chain covers, and which ones.

Because there is no pocket on the target to describe, the descriptive and length-normalised columns
(`target_res_ids`, `target_len`, `target_seq`, `target_pct_aln`, `jaccard_index`) are left empty —
the same shape a Foldseek-database row has. `jaccard_index` in particular needs a second pocket to size
against, and a whole chain's length would swamp the union it normalises by. The overlap itself is still
fully reported: `target_overlap_ids` gives the author residue ids the query pocket maps onto, and
because a named structure has real coordinates the superposition columns (`rmsd`, `ca_dists`, the
transforms) are populated too, which a Foldseek-database row cannot offer.

Open and closed targets can be mixed freely in one batch file.

## Databases

A target can name a Foldseek database instead of a structure, which searches the query pocket against
every entry in it. Both bundled databases require the foldseek binary, and
`--align_struct_method pocket` is rejected against either (there is no second pocket to fit to).

**`human_domains`** ships inside the package — no download. Its entries are structural domains carved out
of human UniProt sequences. There is no interface to compute on a domain, so the target side behaves like
an open search: no `target_*` descriptors, no `jaccard_index`, no superposition columns. The database
ships an offset table, so `target_overlap_ids` is reported in **UniProt residue numbers** rather than
positions within the domain. The query side (`query_overlap_ids`) is always author residue ids.

**`pdb`** is downloaded on first use (`foldseek databases PDB`) into `<cache_dir>/fsdb/pdb`. Hits are
real PDB chains, so PocketMapper fetches PISA data for each one and compares against a real interface
pocket; hits with no usable PISA data are dropped rather than compared against a stand-in. **This is
slow the first time**: `4Q5J:B_F` returns roughly 4,970 hits across 3,620 entries. The default
`--pisa_source ftp` downloads their PISA files concurrently; with `--pisa_source api` each assembly is a
separate rate-limited request, so the first run takes hours. Both counts are logged before the fetching
starts, and reruns are fast from the interface cache.

A Foldseek database you build yourself also works as a target (pass its path, with
`--target_pocket_method foldseek_db`), but only the bundled `human_domains` has an offset table, so its
target residue ids are 0-indexed positions within the entry. That is logged at INFO, because the same
column then means different things on different runs.
