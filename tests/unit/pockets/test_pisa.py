"""
Unit tests for `pocketmapper.pockets.pisa`: PISA pockets by bond and by buried surface, on the seeded 4Q5J
interfaces, reading the interface cache, and the PISA downloader: with downloads off, and online against a fake
server, what it caches and when it rebuilds an interface file.
"""

import json
import os

import pytest

from pocketmapper.lib.structure_manager import StructureManager
from pocketmapper.pockets import pisa_downloader
from pocketmapper.pockets.pisa import load_interfaces, pisa_pockets, pisa_pockets_from_records
from pocketmapper.pockets.pisa_downloader import PisaDownloader, write_json

# The start of the one warning `download_missing_interfaces` logs with downloads off
SKIPPED_WARNING = "Downloads are off"

# (pdb_list, cached files removed first, count of entries skipped). 4Q5J is the seeded entry; 1ABC and
# 2XYZ are not cached. With only its interface file removed, 4Q5J still has its summary and assemblies.
OFFLINE_CASES = [
    pytest.param(["4Q5J"], [], 0, id="all_cached"),
    pytest.param(["4Q5J", "1abc"], [], 1, id="mixed"),
    pytest.param(["1ABC", "2XYZ", "1abc"], [], 2, id="none_cached"),
    pytest.param(["4Q5J"], [os.path.join("interface_pairs", "4q5j.json")], 1, id="interfaces_missing"),
]

# A threshold that removes residues from both seeded pockets
BSA_THRESHOLD = 25

# (pocket_id, residues by bond, buried at 0, buried past BSA_THRESHOLD), counted by hand from the seed
POCKET_SIZES = [
    pytest.param("4Q5J:A_E", 11, 19, 14, id="4Q5J:A_E"),
    pytest.param("4Q5J:B_F", 7, 18, 14, id="4Q5J:B_F"),
]


@pytest.fixture
def interface_dir(pocket_dir):
    """The seeded PISA interface files, in the test's own pocket directory."""
    return os.path.join(pocket_dir, "pisa", "interface_pairs")


@pytest.fixture
def records(processor):
    """The 4Q5J query records, built by the parser."""
    records = []
    for entry in ("4Q5J:A_E", "4Q5J:B_F"):
        ((_, record, reason),) = processor.process_qt_cmdline_input(entry, "query")
        assert reason is None
        records.append(record)
    return records


def residues(records, interface_dir, bsa_threshold=None):
    """pocket_id -> the set of residue ids `pisa_pockets_from_records` gives under `bsa_threshold`."""
    pockets = pisa_pockets_from_records(records, interface_dir, bsa_threshold=bsa_threshold)
    return {pocket_id: set(pocket.res_auth_ids) for pocket_id, pocket in pockets.items()}


@pytest.mark.parametrize("pocket_id, bonded, buried, past_threshold", POCKET_SIZES)
def test_bsa_threshold_narrows_pocket(
    network_attempts, records, interface_dir, pocket_id, bonded, buried, past_threshold
):
    """A threshold keeps a strict subset of the buried residues, and with 0 every bonded residue is buried."""
    by_bond = residues(records, interface_dir)[pocket_id]
    at_zero = residues(records, interface_dir, 0)[pocket_id]
    narrowed = residues(records, interface_dir, BSA_THRESHOLD)[pocket_id]
    assert (len(by_bond), len(at_zero), len(narrowed)) == (bonded, buried, past_threshold)
    assert narrowed < at_zero
    assert by_bond <= at_zero


def tree(root):
    """Every file under `root`, relative to it."""
    return {os.path.relpath(os.path.join(d, f), root) for d, _, files in os.walk(root) for f in files}


@pytest.mark.parametrize("pdb_list, removed, skipped", OFFLINE_CASES)
def test_downloader_offline(no_downloads, warnings, pocket_dir, pdb_list, removed, skipped):
    """With downloads off, entries not in `interface_pairs/` are skipped with one warning, and nothing is written."""
    pisa_dir = os.path.join(pocket_dir, "pisa")
    for path in removed:
        os.remove(os.path.join(pisa_dir, path))
    before = tree(pisa_dir)

    failures = PisaDownloader(pisa_dir, download=False).download_missing_interfaces(pdb_list)

    assert failures == {}
    assert tree(pisa_dir) == before
    assert not os.path.exists(os.path.join(pisa_dir, "errors.json"))
    skipped_warnings = [message for message in warnings() if message.startswith(SKIPPED_WARNING)]
    if skipped:
        assert len(skipped_warnings) == 1
        assert f" {skipped} " in skipped_warnings[0]
    else:
        assert skipped_warnings == []


def test_pisa_pockets_offline(no_downloads, processor, pdb_dir, pocket_dir):
    """With downloads off, a cached entry's pocket is built and an uncached entry's is left out."""
    records = []
    for entry in ("4Q5J:A_E", "1ABC:A_B"):
        ((_, record, reason),) = processor.process_qt_cmdline_input(entry, "query")
        assert reason is None
        records.append(record)

    structures = StructureManager(pdb_dir, download=False)
    pockets = pisa_pockets(records, structures, PisaDownloader(os.path.join(pocket_dir, "pisa"), download=False))

    assert list(pockets) == ["4Q5J:A_E"]
    assert pockets["4Q5J:A_E"].pocket_exists
    assert not os.path.exists(os.path.join(pocket_dir, "pisa", "errors.json"))


# Reading the interface cache and rebuilding it online

# The seeded entry's assembly ids, as its summary lists them
SEED_ASSEMBLIES = ["1", "2"]

# A synthetic entry, served only by the `fake_pisa` fixture
CODE = "9zzz"

# (file content, what `load_interfaces` returns), "seed" standing for the seeded entry's interfaces
LOAD_CASES = [
    pytest.param({"assemblies": SEED_ASSEMBLIES, "interfaces": "seed"}, "seed", id="new_shape"),
    pytest.param({"assemblies": ["1"], "interfaces": {}}, {}, id="no_interfaces"),
    pytest.param("seed", None, id="old_shape"),
    pytest.param(None, None, id="missing"),
]

# (assembly ids the summary lists, chain ids of each assembly's interfaces): entries with no usable interface
NEGATIVE_CASES = [
    pytest.param(["1", "2"], [[["A", "B", "C"]], [["A", "B", "C"]]], id="all_skipped"),
    pytest.param(["1"], [[]], id="no_interfaces"),
    pytest.param([], [], id="no_assemblies"),
]

# How the seeded 4Q5J interface file is left before an online run, which must rebuild it
STALE_CASES = [
    pytest.param("old_shape", id="old_shape"),
    pytest.param("incomplete", id="incomplete"),
]


def summary(code, asm_ids):
    """A PDBe entry summary naming `asm_ids` as the entry's assemblies."""
    return {code: [{"assemblies": [{"assembly_id": asm_id} for asm_id in asm_ids]}]}


def assembly(code, interfaces):
    """A PISA assembly file with one interface per list of chain ids in `interfaces`."""
    entries = [
        {"interface_id": i, "molecules": [{"chain_id": chain_id} for chain_id in chain_ids]}
        for i, chain_ids in enumerate(interfaces, start=1)
    ]
    return {code: {"assembly": {"interfaces": entries}}}


class FakePisa:
    """What the fake server holds, and every request it was sent."""

    def __init__(self):
        """An empty server that has been sent no request."""
        self.summaries = {}  # pdb_code -> summary file content
        self.assemblies = {}  # "<pdb_code>_<assembly_id>" -> assembly file content
        self.requests = []

    def serve_summaries(self, url, out_fpath, data=None, **kwargs):
        """Write the summaries held for the requested codes, as one batch response; False if it holds none."""
        codes = data.decode().split(",")
        self.requests.append(f"summaries:{','.join(codes)}")
        response = {}
        for code in codes:
            response.update(self.summaries.get(code, {}))
        if not response:
            return False
        write_json(response, out_fpath)
        return True

    def serve_assembly(self, url, out_fpath, **kwargs):
        """Write the assembly `url` names if it is held, and report whether it was."""
        self.requests.append(url)
        if url not in self.assemblies:
            return False
        write_json(self.assemblies[url], out_fpath)
        return True


@pytest.fixture
def fake_pisa(monkeypatch, network_attempts):
    """A fake PISA server that summaries and assemblies are requested from, with the network blocked."""
    server = FakePisa()
    monkeypatch.setattr(pisa_downloader, "download_api", server.serve_summaries)
    # An assembly's url is its "<pdb_code>_<assembly_id>" key in the server
    sources = {
        source: (lambda code, asm: f"{code}_{asm}", server.serve_assembly) for source in pisa_downloader.PISA_SOURCES
    }
    monkeypatch.setattr(pisa_downloader, "ASSEMBLY_SOURCES", sources)
    return server


def read_json(path):
    """The parsed content of the JSON file at `path`."""
    with open(path) as f:
        return json.load(f)


@pytest.fixture
def seed_entry(pocket_dir):
    """The seeded 4Q5J interface file's path in the test's pocket directory, and its parsed content."""
    path = os.path.join(pocket_dir, "pisa", "interface_pairs", "4q5j.json")
    return path, read_json(path)


@pytest.mark.parametrize("content, expected", LOAD_CASES)
def test_load_interfaces(network_attempts, interface_dir, seed_entry, content, expected):
    """A new-shape file gives its `"interfaces"`; an old-shape or missing one counts as not cached."""
    path, seed = seed_entry

    def resolve(value):
        """`value` with "seed" standing for the seeded entry's interfaces."""
        if value == "seed":
            return seed["interfaces"]
        if isinstance(value, dict):
            return {key: resolve(item) for key, item in value.items()}
        return value

    if content is None:
        os.remove(path)
    else:
        write_json(resolve(content), path)
    assert load_interfaces("4Q5J", interface_dir) == resolve(expected)


@pytest.mark.parametrize("asm_ids, interfaces", NEGATIVE_CASES)
def test_no_usable_interface_cached(fake_pisa, tmp_path, asm_ids, interfaces):
    """An entry with no usable interface gets a file listing every assembly, and a rerun sends no request."""
    pisa_dir = os.path.join(str(tmp_path), "pisa")
    fake_pisa.summaries[CODE] = summary(CODE, asm_ids)
    for asm_id, asm_interfaces in zip(asm_ids, interfaces, strict=True):
        fake_pisa.assemblies[f"{CODE}_{asm_id}"] = assembly(CODE, asm_interfaces)

    PisaDownloader(pisa_dir).download_missing_interfaces([CODE.upper()])
    assert read_json(os.path.join(pisa_dir, "interface_pairs", f"{CODE}.json")) == {
        "assemblies": asm_ids,
        "interfaces": {},
    }
    assert len(fake_pisa.requests) == 1 + len(asm_ids)

    fake_pisa.requests.clear()
    assert PisaDownloader(pisa_dir).download_missing_interfaces([CODE]) == {}
    assert fake_pisa.requests == []


def test_partial_download_completed(fake_pisa, tmp_path):
    """An assembly that failed to download is left out of the file, and the next run fetches only it."""
    pisa_dir = os.path.join(str(tmp_path), "pisa")
    interface_path = os.path.join(pisa_dir, "interface_pairs", f"{CODE}.json")
    fake_pisa.summaries[CODE] = summary(CODE, ["1", "2", "3"])
    fake_pisa.assemblies[f"{CODE}_1"] = assembly(CODE, [["A", "B"]])
    fake_pisa.assemblies[f"{CODE}_3"] = assembly(CODE, [["C", "D"]])

    failures = PisaDownloader(pisa_dir).download_missing_interfaces([CODE])
    assert failures == {"assembly_downloading": [f"{CODE}_2"]}
    assert read_json(interface_path)["assemblies"] == ["1", "3"]

    fake_pisa.requests.clear()
    fake_pisa.assemblies[f"{CODE}_2"] = assembly(CODE, [["A", "C"]])
    assert PisaDownloader(pisa_dir).download_missing_interfaces([CODE]) == {}
    assert fake_pisa.requests == [f"{CODE}_2"]
    entry = read_json(interface_path)
    assert entry["assemblies"] == ["1", "2", "3"]
    assert list(entry["interfaces"]) == ["A,B", "A,C", "C,D"]


def make_stale(path, seed, how):
    """Rewrite the interface file at `path` in the older flat shape, or listing only the first seeded assembly."""
    if how == "old_shape":
        write_json(seed["interfaces"], path)
    else:
        write_json({"assemblies": SEED_ASSEMBLIES[:1], "interfaces": seed["interfaces"]}, path)


@pytest.mark.parametrize("how", STALE_CASES)
def test_stale_file_rebuilt_online(fake_pisa, pocket_dir, seed_entry, how):
    """Online, an old-shape or incomplete interface file is rebuilt from the cached summary and assemblies."""
    path, seed = seed_entry
    make_stale(path, seed, how)

    assert PisaDownloader(os.path.join(pocket_dir, "pisa")).download_missing_interfaces(["4Q5J"]) == {}
    assert fake_pisa.requests == []
    assert read_json(path) == seed


@pytest.mark.parametrize("how", STALE_CASES)
def test_stale_file_kept_offline(no_downloads, warnings, pocket_dir, seed_entry, how):
    """Offline, an existing interface file counts as cached whatever its shape, and nothing is written."""
    path, seed = seed_entry
    make_stale(path, seed, how)
    with open(path) as f:
        before = f.read()
    pisa_dir = os.path.join(pocket_dir, "pisa")
    files = tree(pisa_dir)

    assert PisaDownloader(pisa_dir, download=False).download_missing_interfaces(["4Q5J"]) == {}
    with open(path) as f:
        assert f.read() == before
    assert tree(pisa_dir) == files
    assert not any(message.startswith(SKIPPED_WARNING) for message in warnings())
