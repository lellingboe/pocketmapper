# PocketMapper

PocketMapper is a command-line tool to compare the binding surfaces of protein-protein interactions.
Give it two pockets, or a pocket and a set of targets, and it reports which residues the surfaces share,
how well those residues are conserved, and how closely the pockets superpose.

> **Status:** early development. Output names and locations may change between releases.

## Installation

PocketMapper supports **Python 3.10 to 3.14** and is published on
[PyPI](https://pypi.org/project/pocketmapper/).

```bash
pip install pocketmapper
```

Installation of [Foldseek](https://github.com/steineggerlab/foldseek) is recommended, an external binary that pip does
not install:

```bash
# With conda
conda install -c conda-forge -c bioconda foldseek

# Or with Homebrew
brew tap brewsci/bio
brew trust brewsci/bio
brew install foldseek
```

Precompiled Foldseek binaries are also available at <https://dev.mmseqs.com/foldseek/>.

Foldseek is the default aligner. To run without it, pass `--aligner seq`, which uses a built-in
BLOSUM62 sequence aligner instead. Searching a Foldseek database target (`human_domains`, `pdb`) always
needs the binary.

Check the installed version with `pocketmapper --version`.

## Quick start

```bash
pocketmapper search 4Q5J:B_F 4Q5J:A_E --results_dir ./out
```

This compares the PISA interface of chain B against chain F of PDB entry 4Q5J with the interface of
chain A against chain E. The main results in `./out`:

| File | Contents |
| --- | --- |
| `pocket_comparison.tsv` | One row per pocket pair: shared residues, overlap scores, RMSD. |
| `alignment.tsv` | The chain alignments the comparison is built on. |
| `aligned_structures/` | The top targets superposed onto each query, as PDB files. |

Downloads are cached in `pocketmapper_cache/`, so a second run over the same structures is fast. See
[Outputs](https://github.com/lellingboe/pocketmapper/blob/main/docs/outputs.md) for every file and column.

## How it works

A `search` runs six steps:

1. **parse**: reads the query and target entries and decides how each pocket is derived.
2. **fetch_structures**: downloads the mmCIF files (or AlphaFold models) it needs.
3. **align**: aligns the query and target chains with Foldseek or the sequence aligner.
4. **pockets**: gets the pocket residues from PDBe PISA, van der Waals contacts, or your list.
5. **compare**: maps both pockets onto the alignment and scores their overlap.
6. **superpose**: writes the top targets superposed onto each query.

Each step is also a command of its own, so you can run, inspect and rerun steps one at a time. See
[Running one step at a time](https://github.com/lellingboe/pocketmapper/blob/main/docs/options.md#running-one-step-at-a-time).

## Specifying entries

Each query and target is written as:

```text
STRUCTURE[:CHAIN[:RESIDUES]]
```

| Example | Pocket |
| --- | --- |
| `4Q5J:B_F` | Interface of chain B against chain F, from PDBe PISA |
| `4Q5J:A_B:10,11,12` | Van der Waals contacts of chain A against chain B |
| `4Q5J:A:10,11,12` | Exactly the residues listed, on chain A |
| `P04637:A:10,11,12` | The same, on the AlphaFold model for a UniProt accession |
| `4Q5J:B` | Open search: the whole of chain B |
| `path/to/model.cif` | A local structure file (chain A, open search) |
| `human_domains`, `pdb` | A Foldseek database (target only) |

An entry that names a pocket on both sides is a *closed search*: how do these two pockets compare? An
entry with no pocket is an *open search*: does the query pocket match anything on this chain?

To compare many entries, pass a file with one entry per line instead. Every query is compared against
every target.

The full rules (which methods each structure type allows, residue checks, and what the Foldseek
databases return) are in
[Entries and targets](https://github.com/lellingboe/pocketmapper/blob/main/docs/input.md).

## Examples

```bash
# One pair, aligned with Foldseek
pocketmapper search 4Q5J:B_F 4Q5J:A_E --results_dir ./out_fs

# The same pair with the built-in sequence aligner (no Foldseek needed)
pocketmapper search 4Q5J:B_F 4Q5J:A_E --aligner seq --results_dir ./out_local

# Open search: is this pocket anywhere on chain A of 4Q5J?
pocketmapper search 4Q5J:B_F 4Q5J:A --results_dir ./out_open

# Search the bundled Foldseek database of human domains
pocketmapper search 4Q5J:B_F human_domains --results_dir ./out_hd

# Batch mode: one entry per line, settings from a job file
pocketmapper search queries.txt targets.txt --job_file job.json
```

Run `pocketmapper search --help` for the full option list.

## Documentation

- [Command-line reference](https://github.com/lellingboe/pocketmapper/blob/main/docs/options.md): every
  option, the step commands, job files and the choice of aligner.
- [Entries and targets](https://github.com/lellingboe/pocketmapper/blob/main/docs/input.md): the entry
  format, pocket methods, open vs closed search, and Foldseek databases.
- [Outputs](https://github.com/lellingboe/pocketmapper/blob/main/docs/outputs.md): every result file,
  the cache, and every `pocket_comparison.tsv` column.
- [Python library](https://github.com/lellingboe/pocketmapper/blob/main/docs/library.md): calling
  PocketMapper from Python.

## Contact / Authors

PocketMapper is developed by Lachlan Ellingboe (<Lachlan.Ellingboe@icr.ac.uk>).

Source, issues and feature requests: <https://github.com/lellingboe/pocketmapper>
