"""
Integration tests for the step commands through the CLI, chained and alone, offline against the seeded cache
unless marked `network`.

`parse` takes `--results_dir` and `--cache_dir`; every later step takes `--job_file <results>/parse_settings.json`.
"""

import json
import os
import shutil

import pytest

from pocketmapper.lib.entries import read_entries


def results_file(results_dir, name):
    """The path of `name` under `results_dir`."""
    return os.path.join(results_dir, name)


def pocket_residues(pockets_path):
    """Pocket id -> its residue auth ids, from a pockets.json."""
    with open(pockets_path) as f:
        return {pocket_id: pocket["res_auth_ids"] for pocket_id, pocket in json.load(f)["pockets"].items()}


def parse_args(*entries, results_dir, cache_dir):
    """The argv of a `parse` over `entries` into `results_dir`, with `cache_dir`."""
    return ("parse", *entries, "--results_dir", results_dir, "--cache_dir", cache_dir)


# test_steps_6 (with its prefixes, test_steps_2-5), test_steps_19, and test_steps_8's cd
def test_steps_6_19_8cd(tmp_path, seed_cache, run_cli, comparison_rows):
    """Each step, chained off parse's settings and run from another cwd, writes its outputs."""
    results_dir = os.path.join(str(tmp_path), "results")
    elsewhere = os.path.join(str(tmp_path), "elsewhere")
    os.makedirs(elsewhere)
    job_file = results_file(results_dir, "parse_settings.json")

    assert run_cli(*parse_args("4Q5J:A_E", "4Q5J:B_F", results_dir=results_dir, cache_dir=seed_cache)) == 0
    assert os.path.isfile(job_file)
    assert os.path.isfile(results_file(results_dir, "entries.tsv"))

    assert run_cli("fetch_structures", "--job_file", job_file, cwd=elsewhere) == 0

    assert run_cli("align", "--aligner", "seq", "--job_file", job_file, cwd=elsewhere) == 0
    assert os.path.isfile(results_file(results_dir, "alignment.tsv"))

    pockets_tsv = results_file(results_dir, "pockets.tsv")
    assert run_cli("pockets", "--job_file", job_file, "--pockets_tsv_path", pockets_tsv, cwd=elsewhere) == 0
    assert os.path.isfile(results_file(results_dir, "pockets.json"))
    assert os.path.isfile(pockets_tsv)

    assert run_cli("compare", "--job_file", job_file, cwd=elsewhere) == 0
    assert len(comparison_rows(results_dir)) > 0

    assert run_cli("superpose", "--job_file", job_file, cwd=elsewhere) == 0
    aligned_dir = results_file(results_dir, "aligned_structures")
    assert os.path.isdir(aligned_dir) and os.listdir(aligned_dir)

    # Nothing was written relative to the cwd the steps ran from
    assert os.listdir(elsewhere) == []


@pytest.mark.parametrize(
    "aligner_args",
    [
        pytest.param(["--aligner", "seq"], id="seq"),
        pytest.param([], id="foldseek", marks=pytest.mark.foldseek),
    ],
)
def test_steps_7(aligner_args, tmp_path, seed_cache, run_cli, comparison_rows, same_outputs):
    """Chained steps give the comparison and alignment tables search gives (PISA batch vs itself)."""
    search_dir = os.path.join(str(tmp_path), "search")
    steps_dir = os.path.join(str(tmp_path), "steps")
    job_file = results_file(steps_dir, "parse_settings.json")

    search_argv = ("search", "pdb_pisa_in.txt", "pdb_pisa_in.txt", *aligner_args)
    assert run_cli(*search_argv, "--results_dir", search_dir, "--cache_dir", seed_cache) == 0
    assert run_cli(*parse_args("pdb_pisa_in.txt", "pdb_pisa_in.txt", results_dir=steps_dir, cache_dir=seed_cache)) == 0
    assert run_cli("fetch_structures", "--job_file", job_file) == 0
    assert run_cli("align", *aligner_args, "--job_file", job_file) == 0
    for step in ("pockets", "compare", "superpose"):
        assert run_cli(step, "--job_file", job_file) == 0, step

    assert len(comparison_rows(steps_dir)) > 0
    assert same_outputs(search_dir, steps_dir)


def test_steps_15(tmp_path, seed_cache, run_cli, comparison_rows):
    """compare reads a search's job_settings.json, and --pocket_comparison_path moves only its output."""
    search_dir = os.path.join(str(tmp_path), "search")
    other_path = os.path.join(str(tmp_path), "other", "pocket_comparison.tsv")
    search_argv = ("search", "4Q5J:A_E", "4Q5J:B_F", "--aligner", "seq", "--align_count", "0")
    assert run_cli(*search_argv, "--results_dir", search_dir, "--cache_dir", seed_cache) == 0

    # search wrote its own comparison; compare must leave it as it was
    default_path = results_file(search_dir, "pocket_comparison.tsv")
    before = os.stat(default_path).st_mtime_ns
    with open(default_path) as f:
        content = f.read()

    job_file = results_file(search_dir, "job_settings.json")
    assert run_cli("compare", "--job_file", job_file, "--pocket_comparison_path", other_path) == 0

    assert len(comparison_rows(os.path.dirname(other_path))) > 0
    assert os.stat(default_path).st_mtime_ns == before
    with open(default_path) as f:
        assert f.read() == content


def test_steps_21(tmp_path, seed_cache, run_cli):
    """fetch_structures follows the entries table parse moved with --entries_path."""
    results_dir = os.path.join(str(tmp_path), "results")
    moved_path = os.path.join(str(tmp_path), "moved", "entries.tsv")
    job_file = results_file(results_dir, "parse_settings.json")

    argv = parse_args("4Q5J:A_E", "4Q5J:B_F", results_dir=results_dir, cache_dir=seed_cache)
    assert run_cli(*argv, "--entries_path", moved_path) == 0
    with open(job_file) as f:
        assert json.load(f)["entries_path"] == moved_path

    # fetch_structures fails on a missing entries table, so exit 0 means it read the moved one
    assert run_cli("fetch_structures", "--job_file", job_file) == 0
    assert os.path.isfile(moved_path)
    assert not os.path.exists(results_file(results_dir, "entries.tsv"))


def test_steps_23(tmp_path, seed_cache, run_cli, comparison_rows, failed):
    """Rerunning pockets clears its own marks and marks the unbuildable pocket again, once."""
    results_dir = os.path.join(str(tmp_path), "results")
    job_file = results_file(results_dir, "parse_settings.json")
    entries_path = results_file(results_dir, "entries.tsv")
    pocket_id = "4Q5J:A:9999"

    # invalid_residues.txt line 2 is deliberately absent from the chain
    assert run_cli(*parse_args("invalid_residues.txt", "4Q5J:B_F", results_dir=results_dir, cache_dir=seed_cache)) == 0
    assert run_cli("fetch_structures", "--job_file", job_file) == 0
    assert run_cli("align", "--aligner", "seq", "--job_file", job_file) == 0
    assert run_cli("pockets", "--job_file", job_file) == 0
    assert failed(entries_path, "pocket_not_built", pocket_id)
    assert run_cli("pockets", "--job_file", job_file) == 0
    assert run_cli("compare", "--job_file", job_file) == 0

    assert failed(entries_path, "pocket_not_built", pocket_id)
    marked = [row for row in read_entries(entries_path) if row["failed_reason"]]
    assert [(row["pocket_id"], row["failed_step"]) for row in marked] == [(pocket_id, "pockets")]
    assert len(comparison_rows(results_dir)) > 0


def test_steps_24(tmp_path, seed_cache, run_cli):
    """pockets takes --pisa_bsa_threshold directly, and it changes the pisa pockets it builds."""
    results_dir = os.path.join(str(tmp_path), "results")
    job_file = results_file(results_dir, "parse_settings.json")
    bonds_path = results_file(results_dir, "pockets.json")
    bsa_path = results_file(results_dir, "pockets_bsa.json")

    assert run_cli(*parse_args("4Q5J:B_F", "4Q5J:A_E", results_dir=results_dir, cache_dir=seed_cache)) == 0
    assert run_cli("fetch_structures", "--job_file", job_file) == 0
    assert run_cli("pockets", "--job_file", job_file) == 0
    assert run_cli("pockets", "--job_file", job_file, "--pisa_bsa_threshold", "25", "--pockets_path", bsa_path) == 0

    # Without a threshold pisa pockets come from interface bonds, so the buried-surface ones differ
    bonds = pocket_residues(bonds_path)
    bsa = pocket_residues(bsa_path)
    assert set(bsa) == set(bonds) == {"4Q5J:B_F", "4Q5J:A_E"}
    assert bsa != bonds


# Step under test, the steps run before it (online), and its own arguments besides --offline 1
OFFLINE_STEPS = [
    pytest.param("fetch_structures", [], [], id="fetch_structures"),
    pytest.param("align", [], ["--aligner", "seq", "--fetch_missing", "1"], id="align"),
    pytest.param("pockets", [], ["--fetch_missing", "1"], id="pockets"),
    pytest.param(
        "superpose",
        [["align", "--aligner", "seq"], ["pockets"], ["compare"]],
        ["--fetch_missing", "1"],
        id="superpose",
    ),
]


@pytest.mark.parametrize("step, before, step_args", OFFLINE_STEPS)
def test_offline_step(step, before, step_args, tmp_path, seed_cache, run_cli, failed):
    """With --offline 1 a step downloads nothing, even with --fetch_missing 1, and marks the uncached entry."""
    results_dir = os.path.join(str(tmp_path), "results")
    job_file = results_file(results_dir, "parse_settings.json")
    entries_path = results_file(results_dir, "entries.tsv")
    uncached = "P24941:A:160"
    query_path = os.path.join(str(tmp_path), "query.txt")
    with open(query_path, "w") as f:
        f.write(f"4Q5J:B_F\n{uncached}\n")

    assert run_cli(*parse_args(query_path, "4Q5J:A_E", results_dir=results_dir, cache_dir=seed_cache)) == 0
    for argv in before:
        assert run_cli(*argv, "--job_file", job_file) == 0, argv[0]

    # Uncached from here on; 4Q5J stays. The no_network guard fails the test on any download attempt.
    os.remove(os.path.join(seed_cache, "alphafold_structures", "P24941.cif.gz"))
    assert run_cli(step, *step_args, "--offline", "1", "--job_file", job_file) == 0

    marked = [(row["pocket_id"], row["failed_step"]) for row in read_entries(entries_path) if row["failed_reason"]]
    assert marked == [(uncached, step)]
    assert failed(entries_path, "structure_not_found", uncached)


@pytest.mark.parametrize("command", ["search", "fetch_structures"])
def test_offline_refresh(command, tmp_path, seed_cache, run_cli):
    """--offline 1 with --refresh_structures 1 exits 1, before search dumps its settings."""
    results_dir = os.path.join(str(tmp_path), "results")
    flags = ("--offline", "1", "--refresh_structures", "1")
    if command == "search":
        argv = ("search", "4Q5J:B_F", "4Q5J:A_E", "--aligner", "seq", "--results_dir", results_dir)
        assert run_cli(*argv, "--cache_dir", seed_cache, *flags) == 1
        assert not os.path.exists(results_file(results_dir, "job_settings.json"))
    else:
        assert run_cli(*parse_args("4Q5J:B_F", "4Q5J:A_E", results_dir=results_dir, cache_dir=seed_cache)) == 0
        job_file = results_file(results_dir, "parse_settings.json")
        assert run_cli("fetch_structures", "--job_file", job_file, *flags) == 1


@pytest.mark.network
def test_steps_18(tmp_path, seed_cache, run_cli):
    """fetch_structures downloads a PDB and an AlphaFold structure into the cache paths parse moved."""
    results_dir = os.path.join(str(tmp_path), "results")
    pdb_dir = os.path.join(results_dir, "pdb")
    alphafold_dir = os.path.join(results_dir, "af")
    job_file = results_file(results_dir, "parse_settings.json")

    argv = parse_args("4Q5J:A_E", "P24941:A:160", results_dir=results_dir, cache_dir=seed_cache)
    assert run_cli(*argv, "--pdb_dir", pdb_dir, "--alphafold_dir", alphafold_dir) == 0
    assert run_cli("fetch_structures", "--job_file", job_file) == 0

    assert os.path.isfile(os.path.join(pdb_dir, "4Q5J.cif.gz"))
    assert os.path.isfile(os.path.join(alphafold_dir, "P24941.cif.gz"))


@pytest.mark.network
def test_steps_20(tmp_path, seed_cache, run_cli):
    """pockets --fetch_missing 1 downloads a structure missing from the moved pdb_dir, then builds the pockets."""
    results_dir = os.path.join(str(tmp_path), "results")
    pdb_dir = os.path.join(str(tmp_path), "pdb")
    job_file = results_file(results_dir, "parse_settings.json")

    # The seed supplies the PISA interfaces, so only the structure is downloaded
    argv = parse_args("4Q5J:B_F", "4Q5J:A_E", results_dir=results_dir, cache_dir=seed_cache)
    assert run_cli(*argv, "--pdb_dir", pdb_dir) == 0
    assert run_cli("pockets", "--job_file", job_file, "--fetch_missing", "1") == 0

    assert os.path.isfile(os.path.join(pdb_dir, "4Q5J.cif.gz"))
    assert os.path.isfile(results_file(results_dir, "pockets.json"))


LOCAL_POCKET = "sub/local.cif.gz:A:1101,1104,1105,1108,1125"
PDB_POCKET = "4Q5J:A:1101,1104,1105,1108,1125"


def unfailed_ids(entries_path):
    """The pocket ids of the entries table's rows that are not marked failed."""
    return {row["pocket_id"] for row in read_entries(entries_path) if not row["failed_reason"]}


@pytest.mark.parametrize("use_dump", [True, False], ids=["dump", "explicit_paths"])
def test_local_root_later_step(use_dump, tmp_path, seed_cache, run_cli):
    """A relative local path parsed from `local_root` resolves in a later step run from another cwd."""
    local_root = str(tmp_path / "work")
    elsewhere = str(tmp_path / "elsewhere")
    os.makedirs(os.path.join(local_root, "sub"))
    os.makedirs(elsewhere)
    shutil.copyfile(
        os.path.join(seed_cache, "pdb_structures", "4Q5J.cif.gz"), os.path.join(local_root, "sub", "local.cif.gz")
    )
    results_dir = os.path.join(str(tmp_path), "results")
    job_file = results_file(results_dir, "parse_settings.json")
    entries_path = results_file(results_dir, "entries.tsv")
    pockets_path = results_file(results_dir, "pockets.json")

    assert (
        run_cli(*parse_args(LOCAL_POCKET, LOCAL_POCKET, results_dir=results_dir, cache_dir=seed_cache), cwd=local_root)
        == 0
    )
    if use_dump:
        # The dump holds local_root, so the step finds the file from another cwd
        argv = ("pockets", "--job_file", job_file, "--offline", "1")
    else:
        argv = (
            "pockets",
            "--entries", entries_path,
            "--alignment", results_file(results_dir, "alignment.tsv"),
            "--pockets_path", pockets_path,
            "--pdb_dir", os.path.join(seed_cache, "pdb_structures"),
            "--alphafold_dir", os.path.join(seed_cache, "alphafold_structures"),
            "--pocket_dir", os.path.join(seed_cache, "pockets"),
            "--fsdb_dir", os.path.join(seed_cache, "foldseek_dbs"),
            "--local_root", local_root,
            "--offline", "1",
        )  # fmt: skip
    assert run_cli(*argv, cwd=elsewhere) == 0

    assert unfailed_ids(entries_path) == {LOCAL_POCKET}
    assert set(pocket_residues(pockets_path)) == {LOCAL_POCKET}
    assert os.listdir(elsewhere) == []


def test_pdb_dir_changed_after_parse(tmp_path, seed_cache, run_cli):
    """A later step locates a PDB structure in its own `--pdb_dir`, not the one parse was given."""
    results_dir = os.path.join(str(tmp_path), "results")
    job_file = results_file(results_dir, "parse_settings.json")
    entries_path = results_file(results_dir, "entries.tsv")
    new_pdb_dir = str(tmp_path / "new_pdb")
    os.makedirs(new_pdb_dir)
    shutil.move(os.path.join(seed_cache, "pdb_structures", "4Q5J.cif.gz"), os.path.join(new_pdb_dir, "4Q5J.cif.gz"))

    assert run_cli(*parse_args(PDB_POCKET, PDB_POCKET, results_dir=results_dir, cache_dir=seed_cache)) == 0
    # The old dir no longer holds the structure, so only the new dir can satisfy the step
    argv = ("pockets", "--job_file", job_file, "--offline", "1")
    assert run_cli(*argv) == 0
    assert unfailed_ids(entries_path) == set()
    assert run_cli(*argv, "--pdb_dir", new_pdb_dir) == 0
    assert unfailed_ids(entries_path) == {PDB_POCKET}
