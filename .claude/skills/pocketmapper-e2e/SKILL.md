---
name: pocketmapper-e2e
description: Run PocketMapper's end-to-end test suite (tests/e2e/run_e2e.sh) against the existing warm cache, add new cases to it, and keep CI's pinned cases intact. Use this skill whenever the user wants to run tests, run the suite, check that a change still works, verify or smoke-test the pipeline, add or edit a test case, mentions run_e2e.sh or e2e_results, or asks "did I break anything" after touching pocketmapper code — even if they never say the words "e2e" or "end-to-end". There are no unit tests in this repo, so any request to test anything means this suite.
---

# PocketMapper e2e suite

`tests/e2e/run_e2e.sh` is the only test suite in this repo. Each of its 23 cases shells out to
the real `pocketmapper` CLI against live wwPDB, AlphaFold and PDBe PISA — no mocks, network
required — and asserts exit status plus the presence (and, where hits are expected,
non-emptiness) of `pocket_comparison.tsv`.

## Always run against the warm cache

`tests/e2e/e2e_results/pocketmapper_cache` holds ~5 GB of fetched structures, PISA responses
and Foldseek databases. A cold cache is not merely slower: PISA is fetched one entry at a time
behind a rate-limiting sleep, so a cold run re-downloads thousands of responses at roughly
3/s. The runner's default out-dir is `$PWD/e2e_results` — resolved against the *current
directory*, not the script — so invoking it from the repo root without `-o` silently creates a
second, empty cache and starts over.

From the repo root, this is the invocation:

```bash
POCKETMAPPER_E2E_CACHE="$PWD/tests/e2e/e2e_results/pocketmapper_cache" \
  tests/e2e/run_e2e.sh -o "$PWD/tests/e2e/e2e_results" -t core
```

`-c` is the flag equivalent of `POCKETMAPPER_E2E_CACHE`; either is fine, but pass one of them
every time. For a direct `pocketmapper search` outside the runner, pass the same path as
`--cache_dir`.

To deliberately test cold-cache behaviour, point `-o` at a throwaway directory *inside*
`tests/e2e/e2e_results/` so the downloads it makes stay reusable next time.

Everything under `tests/e2e/e2e_results/` is gitignored — never commit results, logs or cache.

## Running

The script must be executed, never sourced: it sets `set -u` and calls `exit`, which under a
login shell kills the terminal.

```bash
tests/e2e/run_e2e.sh --list          # the 23 cases with their tags and descriptions
tests/e2e/run_e2e.sh -n -t core      # dry run: print the commands, touch nothing
tests/e2e/run_e2e.sh -t core         # ~35s warm; every group except human_domains searches
tests/e2e/run_e2e.sh test_core_7     # one or more cases by name
tests/e2e/run_e2e.sh -k              # keep existing per-case results instead of clearing
tests/e2e/run_e2e.sh -v 2            # PocketMapper verbosity (4=DEBUG default .. 1=ERROR)
```

Tags: `core` (quick), `local` (the BLOSUM62 aligner, runs without foldseek), `human_domains`,
`slow`, and `huge` (excluded unless the case is named explicitly). The group prefix in a case
name tracks what it exercises, not its tag column — `test_open_*` and five of the six
`test_local_*` cases carry `core` too.

Missing resources SKIP rather than fail: `needs-pdb-fsdb` wants `POCKETMAPPER_PDB_FSDB` pointed
at a prebuilt Foldseek PDB database; `needs-pdb-download` downloads the full PDB Foldseek DB
(2 GB download, 7 GB unzipped), so it is also tagged `huge` and runs only when named.

Case names run in groups: `test_core_*` structure-vs-structure pairs, `test_open_*` open
whole-chain targets, `test_domains_*` `human_domains` DB targets, `test_fsdb_*` larger Foldseek
DB targets, `test_local_*` the local aligner.

### What to run for a given change

Running everything means the `human_domains` searches, which are slow. Pick by what was
touched:

| Changed | Run |
| --- | --- |
| Alignment, pockets, comparison, superposition | `-t core` |
| Local BLOSUM62 aligner (`SequenceAligner`) | `-t local` |
| A pocket method, or `pocket_comparison.py` | `-t core`, then `test_domains_1` |
| The Foldseek-DB path | `test_fsdb_1` (needs `POCKETMAPPER_PDB_FSDB` set) |
| `Settings`, CLI plumbing, option validation | `test_local_5 test_local_6` — the two rejection cases |
| Anything you're unsure of | `-t core` first; it's the cheap signal |

## Reading the result

Each case prints `RUN`/`PASS`/`FAIL`/`SKIP` and the script ends with a `passed: N failed: N
skipped: N` line, exiting 1 only if something failed. Per-case output lands in
`<out-dir>/<name>/` with the full CLI log at `<out-dir>/<name>.log` — read the log, not just
the summary, when a case fails.

**A skip is not a pass.** Without the `foldseek` binary on PATH, 17 of the 23 cases skip and
the run still exits 0. Report the skip count alongside the pass count, and say what was
skipped and why. Foldseek is the CLI's default aligner; install it with
`conda install -c conda-forge -c bioconda foldseek` if the user wants full coverage.

The `expect` field in each case decides what is asserted:

- `rows` — must exit 0 **and** write at least one data row to `pocket_comparison.tsv`.
- `ok` — must exit 0 and write the file; zero hits is a legitimate outcome for that pair.
- `fail` — must exit non-zero (a rejected option combination). Nothing is asserted about output.

## Adding a case

Cases live in the `CASES` heredoc at the top of `run_e2e.sh`, five pipe-separated fields:

```
name | tags | expect | description | args
```

`args` is passed to `pocketmapper search` verbatim. `--verbosity`, `--cache_dir` and
`--results_dir` are appended by the runner — don't put them in a case. Cases run with
`tests/e2e/fixtures/` as their working directory, because `testfile.txt` refers to
`4Q5J.cif.gz` by a relative path; keep that relative reference if you edit the fixtures.

Things that fail silently rather than loudly if you get them wrong:

- **Append to a group's tail; never insert mid-group.** Inserting renumbers everything below
  it in that group, and CI pins cases *by name* (see below).
- **Blank lines separate groups and are skipped by the runner; a `#` line inside the heredoc
  is not** — it would be parsed as a case with a garbage name. Put annotations in the header
  comment above the heredoc, where the rest of this convention is already documented.
- **Foldseek is assumed.** A case is skipped when the binary is missing unless its `args`
  contain the literal string `--foldseek False` — that exact substring is what the gate
  matches, checked before the catch-all. Tag such a case `local`. Keep at least one case on
  the local branch: when every case ran Foldseek, the suite couldn't see that branch, which is
  how it once shipped broken.
- **`args` is word-split on spaces.** No quoted arguments, and no spaces inside a residue list.
- **`@PDB_FSDB@` expands to `$POCKETMAPPER_PDB_FSDB`**; pair it with the `needs-pdb-fsdb` tag
  so the case skips cleanly when that isn't set.
- **Choose `expect` honestly.** Use `rows` only when hits are genuinely guaranteed for that
  pair — an over-optimistic `rows` turns into a flaky failure that hides real regressions.

Then verify: `run_e2e.sh -n <name>` to check the command it builds, then run it by name
against the warm cache and confirm it passes for the reason you expect (check the row count
and the log, not just the `PASS`).

## CI

`.github/workflows/test_and_deploy.yml` gates on lint, then runs exactly two cases —
`test_core_1` and `test_local_1`, one per aligner — on Linux and macOS. Windows is only
covered by a dependency-install job, because foldseek has no Windows build.

- Those two names are **pinned in the workflow**. They sit at the head of their groups so
  appending doesn't move them, but any mid-group insert or rename silently repoints CI at a
  different test. Check the workflow after touching case names.
- A new case does **not** run in CI unless you edit the workflow. Adding one there spends
  runner minutes on a live-network job with a 20-minute timeout and one retry — worth it for a
  genuinely new code path, not for a variation on an existing pair.
- Lint runs first and blocks e2e, so run `black ./ --check -l 120 && flake8` locally before
  pushing.

## Don't

- Don't run the suite from the repo root without `-o` and a cache path — that's the expensive
  mistake this skill exists to prevent.
- Don't add mocks or fixtures that stub the network. These are deliberately live smoke tests
  of the whole pipeline; a mocked case asserts nothing about the thing that breaks.
- Don't edit `build/lib/pocketmapper/` — it's a stale checked-in artifact of an older version.
