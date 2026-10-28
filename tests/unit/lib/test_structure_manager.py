"""
Unit tests for `lib.structure_manager`: cache paths, fetching from a `file://` source, and loading with validation.
"""

import logging
import os
from urllib.parse import urljoin
from urllib.request import pathname2url

import gemmi
import pytest

from pocketmapper.lib import structure_manager
from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.residues import ca_residues
from pocketmapper.lib.structure_manager import StructureManager, load_structure

FIXTURE = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "..", "..", "integration", "fixtures", "4Q5J.cif.gz")
)
SENTINEL = b"sentinel"


def read_bytes(path):
    """The contents of `path`, or None if it does not exist."""
    if not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        return f.read()


def has_ca_residues(st):
    """Whether any chain of the first model of `st` has a CA residue."""
    return any(ca_residues(chain) for chain in st[0])


@pytest.fixture
def file_source(monkeypatch):
    """Point the `pdb` URL builder at a `file://` URL of the 4Q5J fixture, whatever the code."""
    url = urljoin("file:", pathname2url(FIXTURE))
    monkeypatch.setitem(structure_manager.STRUCTURE_SOURCES, "pdb", (lambda pdb_code: url, None))


@pytest.fixture
def pdb_record(work_dir):
    """A 4Q5J record whose `struct_path` is under `work_dir`, not yet written."""
    return {"struct_type": "pdb", "struct_info": "4Q5J", "struct_path": os.path.join(work_dir, "pdb", "4Q5J.cif.gz")}


@pytest.mark.parametrize(
    "kwargs, struct_type, struct_info, expected",
    [
        pytest.param({"pdb_dir": "pdb"}, "pdb", "4q5j", os.path.join("pdb", "4Q5J.cif.gz"), id="pdb_lower"),
        pytest.param({"pdb_dir": "pdb"}, "pdb", "4Q5J", os.path.join("pdb", "4Q5J.cif.gz"), id="pdb_upper"),
        pytest.param(
            {"alphafold_dir": "af"}, "alphafold", "p06493", os.path.join("af", "P06493.cif.gz"), id="alphafold"
        ),
        pytest.param({}, "local_file", os.path.join("in", "my.pdb"), os.path.join("in", "my.pdb"), id="local_file"),
    ],
)
def test_path(kwargs, struct_type, struct_info, expected):
    """`path` builds `<dir>/<ID upper-cased>.cif.gz`, and leaves a local file's path unchanged."""
    assert StructureManager(**kwargs).path(struct_type, struct_info) == expected


@pytest.mark.parametrize(
    "kwargs, struct_type, match",
    [
        pytest.param({"pdb_dir": "pdb"}, "foldseek_db", "Unknown structure type", id="foldseek_db"),
        pytest.param({"pdb_dir": "pdb"}, "nonsense", "Unknown structure type", id="unknown"),
        pytest.param({"alphafold_dir": "af"}, "pdb", "needs pdb_dir", id="no_pdb_dir"),
        pytest.param({"pdb_dir": "pdb"}, "alphafold", "needs alphafold_dir", id="no_alphafold_dir"),
    ],
)
def test_path_rejects(kwargs, struct_type, match):
    """`path` raises for an unknown type, or when the directory the type needs was not given."""
    with pytest.raises(PocketMapperError, match=match):
        StructureManager(**kwargs).path(struct_type, "4Q5J")


@pytest.mark.parametrize(
    "kwargs, existing, available, contents",
    [
        pytest.param({}, None, True, "fixture", id="missing_downloaded"),
        pytest.param({}, SENTINEL, True, SENTINEL, id="cached_kept"),
        pytest.param({"refresh": True}, SENTINEL, True, "fixture", id="refresh_overwrites"),
        pytest.param({"download": False}, None, False, None, id="offline_missing"),
        pytest.param({"download": False}, SENTINEL, True, SENTINEL, id="offline_cached"),
    ],
)
def test_fetch(file_source, pdb_record, kwargs, existing, available, contents):
    """`fetch` downloads a missing file, keeps a cached one unless refreshing, and sends nothing when offline."""
    if existing is not None:
        os.makedirs(os.path.dirname(pdb_record["struct_path"]))
        with open(pdb_record["struct_path"], "wb") as f:
            f.write(existing)
    found = StructureManager(**kwargs).fetch([pdb_record])
    assert found == {"4Q5J": available}
    expected = read_bytes(FIXTURE) if contents == "fixture" else contents
    assert read_bytes(pdb_record["struct_path"]) == expected


@pytest.mark.parametrize("exists", [True, False])
def test_fetch_local_file(work_dir, exists):
    """A local file is reported available only if it exists."""
    struct_path = FIXTURE if exists else os.path.join(work_dir, "missing.cif")
    record = {"struct_type": "local_file", "struct_info": "my_structure", "struct_path": struct_path}
    assert StructureManager().fetch([record]) == {"my_structure": exists}


def test_fetch_unknown_type(work_dir, caplog):
    """An unknown type is reported unavailable, with a warning."""
    record = {"struct_type": "nonsense", "struct_info": "X", "struct_path": os.path.join(work_dir, "X")}
    with caplog.at_level(logging.WARNING):
        assert StructureManager().fetch([record]) == {"X": False}
    assert "Unknown structure type nonsense" in caplog.text


def test_fetch_dedupes(file_source, pdb_record, work_dir, monkeypatch):
    """Records sharing a `struct_info` give one download job, to the first record's `struct_path`."""
    calls = []
    real_download_many = structure_manager.download_many

    def recording_download_many(jobs, **kwargs):
        """Record the jobs, then run the real `download_many`."""
        calls.append(dict(jobs))
        return real_download_many(jobs, **kwargs)

    monkeypatch.setattr(structure_manager, "download_many", recording_download_many)
    duplicate = dict(pdb_record, struct_path=os.path.join(work_dir, "other", "4Q5J.cif.gz"))
    found = StructureManager().fetch([pdb_record, duplicate])
    assert found == {"4Q5J": True}
    assert [[job.out_fpath for job in jobs.values()] for jobs in calls] == [[pdb_record["struct_path"]]]


def test_refresh_needs_download():
    """`refresh=True` with `download=False` is rejected."""
    with pytest.raises(PocketMapperError, match="cannot be refreshed"):
        StructureManager(download=False, refresh=True)


@pytest.fixture
def pdb_copy(work_dir):
    """A factory that writes the fixture to `name` under `work_dir` in PDB format and returns its path."""

    def write(name, options=None):
        """Write the fixture as PDB to `name`, with `options` if given, and return the path."""
        path = os.path.join(work_dir, name)
        st = gemmi.read_structure(FIXTURE)
        if options is None:
            st.write_pdb(path)
        else:
            st.write_pdb(path, options)
        return path

    return write


def test_load_structure_mmcif():
    """The mmCIF fixture loads with CA residues."""
    st = load_structure(FIXTURE)
    assert st is not None and has_ca_residues(st)


def test_load_structure_pdb_format(pdb_copy):
    """A PDB-format copy of the fixture loads with CA residues."""
    st = load_structure(pdb_copy("4Q5J.pdb"))
    assert st is not None and has_ca_residues(st)


@pytest.mark.parametrize("entities, has_entities", [(True, True), (False, False)])
def test_load_structure_entities(pdb_copy, entities, has_entities):
    """`entities=False` skips `setup_entities()`, on a PDB copy without SEQRES, which has no entities as read."""
    path = pdb_copy("4Q5J_noseqres.pdb", gemmi.PdbWriteOptions(seqres_records=False))
    st = load_structure(path, entities=entities)
    assert (len(st.entities) > 0) == has_entities


def truncated_fixture():
    """The first half of the gzipped fixture."""
    data = read_bytes(FIXTURE)
    return data[: len(data) // 2]


@pytest.mark.parametrize(
    "name, content, message",
    [
        pytest.param("missing.cif", None, "does not exist", id="missing"),
        pytest.param("truncated.cif.gz", truncated_fixture(), "cannot be read", id="truncated"),
        pytest.param(
            "malformed.cif", "data_x\nloop_\n_atom_site.id\n'unterminated\n", "cannot be read", id="malformed"
        ),
        pytest.param("empty.cif", "", "cannot be read", id="empty"),
        pytest.param("unknown.xyz", "hello\n", "cannot be read", id="unknown_extension"),
        pytest.param("no_model.cif", "data_x\n_cell.length_a 10\n", "has no model", id="no_model"),
        pytest.param("no_chain.pdb", "not a structure\n", "has no chain", id="no_chain"),
    ],
)
def test_load_structure_rejects(write_file, work_dir, caplog, name, content, message):
    """A missing, unreadable, model-less or chain-less file gives None and a warning."""
    path = os.path.join(work_dir, name) if content is None else write_file(name, content)
    with caplog.at_level(logging.WARNING):
        assert load_structure(path) is None
    assert message in caplog.text


def test_get(file_source, work_dir):
    """`get` goes from no file to a loaded structure in one call, caching the file."""
    pdb_dir = os.path.join(work_dir, "pdb")
    st = StructureManager(pdb_dir=pdb_dir).get("pdb", "4q5j")
    assert st is not None and has_ca_residues(st)
    assert os.path.exists(os.path.join(pdb_dir, "4Q5J.cif.gz"))
