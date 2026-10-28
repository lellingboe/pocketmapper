"""
Integration tests for Foldseek-database targets through the CLI: the bundled `human_domains` and a PDB database.

`human_domains` ships in the package and its target structures are extracted from it, so those cases run offline
from the seed cache. A PDB database's hits are expanded into PISA records, which downloads PISA data and a
structure for every hit entry, so those cases are also marked `network` and share the network cache.
"""

import os

import pytest

# The CDK2 pocket of test_domains_4, on 1B38 chain A
CDK2_RESIDUES = (
    "8,9,10,11,12,13,14,15,16,17,18,19,20,30,31,32,33,34,35,47,48,49,50,51,52,53,54,55,56,57,58,59,61,62,63,64,65,"
    "66,67,68,69,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,117,118,119,120,121,122,123,124,125,126,127,"
    "128,129,130,131,132,133,134,135,143,144,145,146,147,148,149"
)

# Each row: the search arguments, before --results_dir and --cache_dir, and the cache fixture. Run from fixtures/,
# which testfile.txt needs for its relative 4Q5J.cif.gz. The id is the old e2e case name.
DOMAINS_CASES = [
    pytest.param(["4Q5J:B_F", "human_domains"], "seed_cache", id="test_domains_1"),
    pytest.param(["testfile.txt", "human_domains"], "seed_cache", id="test_domains_2"),
    pytest.param(
        [f"1B38:A:{CDK2_RESIDUES}", "human_domains"],
        "network_cache",
        id="test_domains_4",
        marks=pytest.mark.network,
    ),
]

# The steps after parse, in order; each reads parse's settings
LATER_STEPS = ["fetch_structures", "align", "pockets", "compare", "superpose"]


def pdb_fsdb_path():
    """The PDB Foldseek database named by `POCKETMAPPER_PDB_FSDB`, which the `pdb_fsdb` marker checks exists."""
    return os.path.abspath(os.environ["POCKETMAPPER_PDB_FSDB"])


def run_chain(run_cli, parse_args, cache_dir, results_dir):
    """Run parse from fixtures/, then every later step from `results_dir` on parse's settings; assert each exits 0."""
    code = run_cli("parse", *parse_args, "--results_dir", results_dir, "--cache_dir", cache_dir)
    assert code == 0, "parse failed"
    settings_path = os.path.join(results_dir, "parse_settings.json")
    for step in LATER_STEPS:
        code = run_cli(step, "--job_file", settings_path, cwd=results_dir)
        assert code == 0, f"{step} failed"


@pytest.mark.foldseek
@pytest.mark.parametrize("args, cache", DOMAINS_CASES)
def test_domains(args, cache, request, tmp_path, run_cli, comparison_rows):
    """A search against human_domains exits 0 and writes at least one comparison row."""
    cache_dir = request.getfixturevalue(cache)
    results_dir = os.path.join(str(tmp_path), "results")
    code = run_cli("search", *args, "--results_dir", results_dir, "--cache_dir", cache_dir)
    assert code == 0
    assert len(comparison_rows(results_dir)) > 0


@pytest.mark.foldseek
def test_steps_8(tmp_path, seed_cache, run_cli, same_outputs):
    """Chained steps against human_domains, run from the results dir after parse, give what search gives."""
    search_dir = os.path.join(str(tmp_path), "results", "search")
    results_dir = os.path.join(str(tmp_path), "results")
    args = ["testfile.txt", "human_domains"]
    assert run_cli("search", *args, "--results_dir", search_dir, "--cache_dir", seed_cache) == 0
    run_chain(run_cli, args, seed_cache, results_dir)
    assert same_outputs(results_dir, search_dir)


@pytest.mark.foldseek
@pytest.mark.slow
@pytest.mark.pdb_fsdb
@pytest.mark.network
def test_fsdb_1(tmp_path, network_cache, run_cli, comparison_rows):
    """A PISA interface against a local PDB Foldseek database exits 0 and writes at least one comparison row."""
    results_dir = os.path.join(str(tmp_path), "results")
    args = ["4Q5J:B_F", pdb_fsdb_path(), "--target_pocket_method", "foldseek_db"]
    assert run_cli("search", *args, "--results_dir", results_dir, "--cache_dir", network_cache) == 0
    assert len(comparison_rows(results_dir)) > 0


@pytest.mark.foldseek
@pytest.mark.slow
@pytest.mark.pdb_fsdb
@pytest.mark.network
def test_steps_17(tmp_path, network_cache, run_cli, same_outputs):
    """Chained steps against a PDB Foldseek database give what search gives, and write the aligned structures."""
    search_dir = os.path.join(str(tmp_path), "results", "search")
    results_dir = os.path.join(str(tmp_path), "results")
    args = ["4Q5J:B_F", pdb_fsdb_path(), "--target_pocket_method", "foldseek_db"]
    assert run_cli("search", *args, "--results_dir", search_dir, "--cache_dir", network_cache) == 0
    run_chain(run_cli, args, network_cache, results_dir)
    assert same_outputs(results_dir, search_dir)
    aligned_dir = os.path.join(results_dir, "aligned_structures")
    assert os.path.isdir(aligned_dir) and os.listdir(aligned_dir), f"no aligned structures in {aligned_dir}"


@pytest.mark.foldseek
@pytest.mark.huge
@pytest.mark.network
def test_fsdb_2(tmp_path, network_cache, run_cli, comparison_rows):
    """A PISA interface against the downloaded full-PDB Foldseek database writes at least one comparison row."""
    results_dir = os.path.join(str(tmp_path), "results")
    code = run_cli("search", "4Q5J:A_E", "pdb", "--results_dir", results_dir, "--cache_dir", network_cache)
    assert code == 0
    assert len(comparison_rows(results_dir)) > 0
