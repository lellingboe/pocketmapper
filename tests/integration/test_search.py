"""
Integration tests for single `search` runs through the CLI.

Every test passes --results_dir and --cache_dir. Unless marked `network`, a test runs offline against its own copy
of the seed cache. Each test id is the e2e case it replaces, except the self-comparison, which is new.
"""

import json
import os
import shutil

import pytest

from pocketmapper.lib.entries import read_entries

# The PISA interface pair most cases compare
QUERY_ID = "4Q5J:A_E"
TARGET_ID = "4Q5J:B_F"

# Each row: extra search arguments, then the aligner and align_struct_method job_settings.json must record.
# The id is the old e2e case name.
SEARCH_CASES = [
    pytest.param([], "foldseek", "foldseek", id="test_core_1", marks=pytest.mark.foldseek),
    pytest.param(
        ["--align_struct_method", "pocket"], "foldseek", "pocket", id="test_core_8", marks=pytest.mark.foldseek
    ),
    pytest.param(["--aligner", "seq"], "seq", "pocket", id="test_local_1"),
    pytest.param(["--aligner", "seq", "--align_struct_method", "pocket"], "seq", "pocket", id="test_local_4"),
]


def run_search(run_cli, results_dir, cache_dir, *args, cwd=None):
    """Run `search` with `args` into `results_dir` against `cache_dir`, assert it exits 0 and return `results_dir`."""
    code = run_cli("search", *args, "--results_dir", results_dir, "--cache_dir", cache_dir, cwd=cwd)
    assert code == 0, f"search {' '.join(args)} exited {code}"
    return results_dir


def read_json(path):
    """The JSON object in the file at `path`."""
    with open(path) as f:
        return json.load(f)


def find_row(rows, query, target):
    """The single comparison row for `query` against `target`."""
    found = [row for row in rows if row["query"] == query and row["target"] == target]
    assert len(found) == 1, f"expected one {query} vs {target} row, got {len(found)}"
    return found[0]


def unfailed(entries_path):
    """Assert that no entry of `entries_path` is marked failed, and return its rows."""
    entries = read_entries(entries_path)
    marked = [(row["pocket_id"], row["failed_reason"]) for row in entries if row["failed_reason"] is not None]
    assert not marked, f"entries marked failed: {marked}"
    return entries


def assert_self_rows_complete(rows):
    """Assert every row comparing a pocket with itself overlaps on every residue, the `seq_pos` check."""
    self_rows = [row for row in rows if row["query"] == row["target"]]
    assert self_rows, "no self-comparison rows"
    for row in self_rows:
        assert int(row["overlap_count"]) == int(row["query_len"]) == int(row["target_len"]), row["query"]
        assert row["query_overlap_ids"] == row["query_res_ids"], row["query"]


@pytest.mark.parametrize("args, aligner, method", SEARCH_CASES)
def test_search(args, aligner, method, tmp_path, seed_cache, run_cli, comparison_rows):
    """The 4Q5J interface pair overlaps, and the top target is superposed with the recorded method."""
    results_dir = run_search(run_cli, str(tmp_path / "results"), seed_cache, QUERY_ID, TARGET_ID, *args)
    row = find_row(comparison_rows(results_dir), QUERY_ID, TARGET_ID)
    assert int(row["overlap_count"]) > 0
    # Foldseek scores its alignment; the sequence aligner writes "-" for both
    scored = row["evalue"] != "-" and row["lddt"] != "-"
    assert scored == (aligner == "foldseek")

    settings = read_json(os.path.join(results_dir, "job_settings.json"))
    assert (settings["aligner"], settings["align_struct_method"]) == (aligner, method)
    aligned = os.listdir(os.path.join(results_dir, "aligned_structures"))
    assert len(aligned) == 1 and aligned[0].endswith(".pdb")
    with open(os.path.join(results_dir, "info.log")) as f:
        assert f"structural alignment of target structures on the {method}" in f.read()


def test_self_comparison(tmp_path, seed_cache, run_cli, comparison_rows):
    """A pocket compared with itself overlaps on every residue: `seq_pos` agrees with the alignment."""
    results_dir = run_search(run_cli, str(tmp_path / "results"), seed_cache, QUERY_ID, QUERY_ID, "--aligner", "seq")
    rows = comparison_rows(results_dir)
    row = find_row(rows, QUERY_ID, QUERY_ID)
    assert int(row["query_len"]) > 1
    assert float(row["jaccard_index"]) == 1.0
    assert_self_rows_complete(rows)


def test_local_2(tmp_path, seed_cache, run_cli, comparison_rows, distinct_queries):
    """Mixed inputs (PDB, local mmCIF, AlphaFold) all compare, and each matches itself on every residue."""
    results_dir = run_search(
        run_cli, str(tmp_path / "results"), seed_cache, "testfile.txt", "testfile.txt", "--aligner", "seq"
    )
    entries = unfailed(os.path.join(results_dir, "entries.tsv"))
    assert {row["struct_type"] for row in entries} == {"pdb", "local_file", "alphafold"}
    assert len(distinct_queries(results_dir)) == 3
    assert_self_rows_complete(comparison_rows(results_dir))


def test_local_3(tmp_path, seed_cache, run_cli, comparison_rows):
    """An open whole-chain target has no pocket of its own and overlaps every aligned query residue."""
    results_dir = run_search(run_cli, str(tmp_path / "results"), seed_cache, QUERY_ID, "4Q5J:B", "--aligner", "seq")
    row = find_row(comparison_rows(results_dir), QUERY_ID, "4Q5J:B")
    assert row["target_len"] == "" and row["jaccard_index"] == ""
    aligned_count = round(float(row["query_pct_aln"]) * int(row["query_len"]))
    assert int(row["overlap_count"]) == aligned_count > 0


def test_local_7(tmp_path, seed_cache, run_cli, distinct_queries):
    """Same-named local files in different directories get distinct names, so neither query is lost."""
    work_dir = str(tmp_path / "work")
    # Different contents under one basename, as the old same_name/ fixture had: a 4Q5J copy and a P24941 copy
    sources = {
        "a": os.path.join(seed_cache, "pdb_structures", "4Q5J.cif.gz"),
        "b": os.path.join(seed_cache, "alphafold_structures", "P24941.cif.gz"),
    }
    for name, source in sources.items():
        os.makedirs(os.path.join(work_dir, name))
        shutil.copyfile(source, os.path.join(work_dir, name, "structure.cif.gz"))
    with open(os.path.join(work_dir, "same_name.txt"), "w") as f:
        f.write("a/structure.cif.gz:A:1101,1104,1105,1108,1125\nb/structure.cif.gz:A:160,161,162,163,164,165\n")

    results_dir = run_search(
        run_cli,
        str(tmp_path / "results"),
        seed_cache,
        "same_name.txt",
        "same_name.txt",
        "--aligner",
        "seq",
        cwd=work_dir,
    )
    entries = unfailed(os.path.join(results_dir, "entries.tsv"))
    assert len({row["preprocess_name"] for row in entries}) == 2
    assert len(distinct_queries(results_dir)) == 2


def test_core_4(tmp_path, seed_cache, run_cli, comparison_rows):
    """A PISA interface against a single-residue AlphaFold pocket runs and drops neither entry; no overlap is fine."""
    results_dir = run_search(
        run_cli, str(tmp_path / "results"), seed_cache, TARGET_ID, "P24941:A:160", "--aligner", "seq"
    )
    unfailed(os.path.join(results_dir, "entries.tsv"))
    # Zero rows is a legitimate outcome, but the table must be written
    comparison_rows(results_dir)


def test_core_6(tmp_path, seed_cache, run_cli, comparison_rows):
    """Two AlphaFold passthrough pockets overlap."""
    query = "P06493:A:160,161,162,163,164,165"
    target = "P24941:A:160,161,162,163,164,165"
    results_dir = run_search(run_cli, str(tmp_path / "results"), seed_cache, query, target, "--aligner", "seq")
    entries = unfailed(os.path.join(results_dir, "entries.tsv"))
    assert {row["struct_type"] for row in entries} == {"alphafold"}
    row = find_row(comparison_rows(results_dir), query, target)
    assert int(row["overlap_count"]) > 0


def test_core_7(tmp_path, seed_cache, run_cli, comparison_rows, distinct_queries):
    """Two pockets on one query chain (pisa and passthrough) are each compared, with their own residues."""
    pockets = {"4Q5J:A_E", "4Q5J:A:1101,1104,1105,1108,1125"}
    results_dir = run_search(
        run_cli, str(tmp_path / "results"), seed_cache, "multi_pocket_chain.txt", TARGET_ID, "--aligner", "seq"
    )
    entries = unfailed(os.path.join(results_dir, "entries.tsv"))
    queries = [row for row in entries if row["side"] == "query"]
    assert {row["pocket_id"] for row in queries} == pockets
    assert len({row["preprocess_name"] for row in queries}) == 1, "the two pockets should share one chain"
    assert distinct_queries(results_dir) == pockets
    rows = comparison_rows(results_dir)
    assert len({find_row(rows, pocket, TARGET_ID)["query_res_ids"] for pocket in pockets}) == 2


def test_settings_1(tmp_path, seed_cache, run_cli):
    """Every path option set on the command line is used, and the default location for it is not."""
    root = str(tmp_path)
    results_dir = os.path.join(root, "results")
    # The cache dirs move out of cache_dir, so a run that ignored the options would find nothing there
    cache_paths = {"pdb_dir": "pdb_structures", "alphafold_dir": "alphafold_structures", "pocket_dir": "pockets"}
    custom = {}
    for key, name in cache_paths.items():
        custom[key] = os.path.join(root, f"custom_{name}")
        shutil.move(os.path.join(seed_cache, name), custom[key])
    custom |= {
        "alignment_path": os.path.join(results_dir, "custom_alignment.tsv"),
        "aligned_structure_dir": os.path.join(results_dir, "custom_aligned"),
        "job_settings_path": os.path.join(results_dir, "custom_settings.json"),
        "log_path": os.path.join(results_dir, "custom.log"),
        "temp_dir": os.path.join(results_dir, "custom_temp"),
    }
    options = [word for key, path in custom.items() for word in (f"--{key}", path)]
    run_search(run_cli, results_dir, seed_cache, QUERY_ID, TARGET_ID, "--aligner", "seq", *options)

    settings = read_json(custom["job_settings_path"])
    for key, path in custom.items():
        if key != "job_settings_path":
            assert settings[key] == path, key
    for key in ("alignment_path", "log_path"):
        assert os.path.getsize(custom[key]) > 0, key
    assert os.listdir(custom["aligned_structure_dir"])
    # The PISA data is only under the custom pocket_dir and the network is blocked, so a built pocket shows it was used
    assert read_json(settings["pockets_path"])["pockets"][QUERY_ID] is not None
    entries = read_entries(settings["entries_path"])
    assert all(os.path.dirname(row["struct_path"]) == custom["pdb_dir"] for row in entries)

    # Nothing at the defaults. temp_dir is deleted on success, so only its default's absence is observable.
    for name in cache_paths.values():
        assert not os.path.exists(os.path.join(seed_cache, name)), name
    for name in ("alignment.tsv", "aligned_structures", "job_settings.json", "info.log", "tmp"):
        assert not os.path.exists(os.path.join(results_dir, name)), name


@pytest.mark.network
@pytest.mark.foldseek
def test_core_3(tmp_path, network_cache, run_cli, comparison_rows, distinct_queries):
    """Chain IDs differing only in case (4DX9 a_b and A_B) stay two pockets with their own residues."""
    pockets = {"4DX9:a_b", "4DX9:A_B"}
    results_dir = run_search(run_cli, str(tmp_path / "results"), network_cache, "4DX9.txt", "4DX9.txt")
    entries = unfailed(os.path.join(results_dir, "entries.tsv"))
    assert len({row["preprocess_name"] for row in entries}) == 2
    assert distinct_queries(results_dir) == pockets
    rows = comparison_rows(results_dir)
    assert (
        find_row(rows, "4DX9:a_b", "4DX9:a_b")["query_res_ids"]
        != find_row(rows, "4DX9:A_B", "4DX9:A_B")["query_res_ids"]
    )
    assert_self_rows_complete(rows)


@pytest.mark.network
@pytest.mark.foldseek
def test_core_10(tmp_path, network_cache, run_cli, comparison_rows):
    """A microheterogeneous residue (4Z0Y:A 252, HS8 and HIS) leaves the chain walk aligned with Foldseek's."""
    results_dir = run_search(run_cli, str(tmp_path / "results"), network_cache, "4Z0Y:A_E", "4Z0Y:C_G")
    unfailed(os.path.join(results_dir, "entries.tsv"))
    row = find_row(comparison_rows(results_dir), "4Z0Y:A_E", "4Z0Y:C_G")
    # Chains A and C are the same protein, so every overlapping residue has the same letter on both
    assert int(row["overlap_count"]) > 0
    assert float(row["overlap_identity"]) == 1.0
