"""
Unit tests for `pocketmapper.parse.qt_processor`: the input grammar, path resolution and `QTRecord` building.
"""

import builtins
import hashlib
import logging
import os
import re

import pytest

from pocketmapper.lib.entries import RECORD_FIELDS
from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.foldseek import BUNDLED_HUMAN_DOMAINS_DB
from pocketmapper.parse.qt_processor import QTProcessor

STRUCT_IDS = {"pdb": "4Q5J", "alphafold": "P24941", "local_file": "struct.cif"}

# (pocket info, {struct_type: inferred method}). The entry is "<id><pocket info>" for each struct_type's id.
INFERENCE_TABLE = [
    ("", {"pdb": "whole_chain", "alphafold": "whole_chain", "local_file": "whole_chain"}),
    (":", {"pdb": "whole_chain", "alphafold": "whole_chain", "local_file": "whole_chain"}),
    (":A", {"pdb": "whole_chain", "alphafold": "whole_chain", "local_file": "whole_chain"}),
    (":A:", {"pdb": "whole_chain", "alphafold": "whole_chain", "local_file": "whole_chain"}),
    (":A_B", {"pdb": "pisa", "alphafold": None, "local_file": "vdw"}),
    (":A:10,11", {"pdb": "passthrough", "alphafold": "passthrough", "local_file": "passthrough"}),
    (":A_B:10,11", {"pdb": "vdw", "alphafold": None, "local_file": "vdw"}),
    (":A_B:", {"pdb": "vdw", "alphafold": None, "local_file": "vdw"}),
    (":AB", {"pdb": None, "alphafold": None, "local_file": None}),
    (":A:x", {"pdb": None, "alphafold": None, "local_file": None}),
]

INFERENCE_CASES = [
    pytest.param(f"{STRUCT_IDS[struct_type]}{info}", struct_type, expected, id=f"{struct_type}-{info or 'bare'}")
    for info, by_type in INFERENCE_TABLE
    for struct_type, expected in by_type.items()
]


def md5(data):
    """The md5 hex digest of `data` (str or bytes)."""
    return hashlib.md5(data if isinstance(data, bytes) else data.encode()).hexdigest()


# Construction and paths


def test_work_dir_defaults_to_cwd(tmp_path, monkeypatch):
    """With no `work_dir`, the working directory is used and the cache directories resolve against it."""
    monkeypatch.chdir(tmp_path)
    processor = QTProcessor(pdb_dir="pdb", alphafold_dir="af", fsdb_dir="fsdb")
    assert processor.work_dir == os.getcwd()
    assert os.path.samefile(processor.work_dir, str(tmp_path))
    assert processor.pdb_dir == os.path.join(os.getcwd(), "pdb")
    assert processor.alphafold_dir == os.path.join(os.getcwd(), "af")


def test_relative_cache_dirs_resolve_against_work_dir(processor, work_dir):
    """Relative `pdb_dir`/`alphafold_dir` are made absolute against `work_dir`."""
    assert processor.pdb_dir == os.path.join(work_dir, "pdb")
    assert processor.alphafold_dir == os.path.join(work_dir, "af")


def test_absolute_cache_dirs_are_kept(work_dir):
    """Absolute `pdb_dir`/`alphafold_dir` are kept as given."""
    pdb_dir = os.path.join(work_dir, "cache", "pdb")
    alphafold_dir = os.path.join(work_dir, "cache", "af")
    processor = QTProcessor(pdb_dir, alphafold_dir, "fsdb", work_dir=os.path.join(work_dir, "work"))
    assert processor.pdb_dir == pdb_dir
    assert processor.alphafold_dir == alphafold_dir


def test_work_path_normalises(processor, work_dir):
    """`work_path` resolves against `work_dir` and normalises `..`."""
    assert processor.work_path(os.path.join("a", "..", "b")) == os.path.join(work_dir, "b")


def test_work_path_keeps_absolute(processor, work_dir):
    """`work_path` keeps an absolute path."""
    path = os.path.join(work_dir, "elsewhere", "file.cif")
    assert processor.work_path(path) == path


def test_accepted_pocket_methods(processor):
    """The accepted methods, in the order the error message lists them."""
    assert processor.accepted_pocket_methods() == ("auto", "whole_chain", "pisa", "passthrough", "vdw", "foldseek_db")


# determine_struct_type


@pytest.mark.parametrize(
    "struct_str, expected",
    [
        ("4Q5J", "pdb"),
        ("4q5j", "pdb"),
        ("1abc", "pdb"),
        ("P24941", "alphafold"),
        ("A0A023GPI8", "alphafold"),
        ("struct.cif", "local_file"),
        ("XYZ", None),
        ("4Q5JX", None),
        ("missing.cif", None),
    ],
)
def test_determine_struct_type(processor, write_file, struct_str, expected):
    """Each identifier is classified by pattern, then by an existing file; anything else gives a reason."""
    write_file("struct.cif", b"any bytes")
    struct_type, reason = processor.determine_struct_type(struct_str)
    assert struct_type == expected
    if expected is None:
        assert struct_str in reason
    else:
        assert reason is None


def test_accession_beats_file(processor, write_file):
    """A file named like a PDB id is still read as a PDB id."""
    write_file("4Q5J", b"any bytes")
    assert processor.determine_struct_type("4Q5J") == ("pdb", None)


def test_directory_struct_raises(processor, work_dir):
    """A directory is not a supported structure."""
    os.makedirs(os.path.join(work_dir, "somedir"))
    with pytest.raises(PocketMapperError, match="Directory input"):
        processor.determine_struct_type("somedir")


# determine_ref_struct_path


@pytest.mark.parametrize(
    "struct_info, struct_type, expected_parts",
    [
        ("4Q5J", "pdb", ("pdb", "4Q5J.cif.gz")),
        ("P24941", "alphafold", ("af", "P24941.cif.gz")),
        ("struct.cif", "local_file", ("struct.cif",)),
        (os.path.join("sub", "struct.cif"), "local_file", ("sub", "struct.cif")),
    ],
)
def test_determine_ref_struct_path(processor, work_dir, struct_info, struct_type, expected_parts):
    """Fetched structures live in their cache directory; a local file's path is made absolute."""
    assert processor.determine_ref_struct_path(struct_info, struct_type) == os.path.join(work_dir, *expected_parts)


def test_determine_ref_struct_path_unknown_type(processor):
    """An unknown struct_type raises."""
    with pytest.raises(PocketMapperError, match="Unknown structure type"):
        processor.determine_ref_struct_path("4Q5J", "nonsense")


# pocket_info


@pytest.mark.parametrize("entry, expected", [("4Q5J", ""), ("4Q5J:A", "A"), ("4Q5J:A:1,2", "A:1,2")])
def test_pocket_info(processor, entry, expected):
    """The pocket info is everything after the first ':'."""
    assert processor.pocket_info(entry) == expected


# determine_pocket_method


@pytest.mark.parametrize("entry, struct_type, expected", INFERENCE_CASES)
def test_determine_pocket_method(processor, entry, struct_type, expected):
    """The first supported method whose pattern matches the pocket info is inferred."""
    assert processor.determine_pocket_method(entry, struct_type) == expected


def test_determine_pocket_method_unknown_struct_type(processor):
    """An unknown struct_type supports no method."""
    assert processor.determine_pocket_method("4Q5J:A", "nonsense") is None


# validate_pocket_method


@pytest.mark.parametrize(
    "entry, method, struct_type, supported",
    [
        ("P24941:A_B", "pisa", "alphafold", "whole_chain, passthrough"),
        ("struct.cif:A_B", "pisa", "local_file", "whole_chain, passthrough, vdw"),
        ("4Q5J:A", "bogus", "pdb", "whole_chain, pisa, passthrough, vdw"),
    ],
)
def test_validate_unsupported_method(processor, entry, method, struct_type, supported):
    """A method the struct_type does not support gives a reason listing the supported ones, in order."""
    reason = processor.validate_pocket_method(entry, method, struct_type)
    assert reason.endswith(f"{struct_type} entries support: {supported}")


@pytest.mark.parametrize(
    "entry, method, struct_type",
    [
        ("4Q5J:A", "passthrough", "pdb"),
        ("4Q5J:A", "pisa", "pdb"),
        ("4Q5J:A", "vdw", "pdb"),
        ("struct.cif:A", "vdw", "local_file"),
    ],
)
def test_validate_pattern_mismatch(processor, entry, method, struct_type):
    """A supported method whose pattern does not match gives a reason with what the method needs."""
    reason = processor.validate_pocket_method(entry, method, struct_type)
    assert reason is not None
    assert processor.pocket_methods[method][1] in reason


@pytest.mark.parametrize(
    "entry, method, struct_type",
    [
        ("4Q5J:A", "whole_chain", "pdb"),
        ("4Q5J:A_B", "vdw", "pdb"),
        ("P24941:A:1,2", "passthrough", "alphafold"),
        ("struct.cif:A_B", "vdw", "local_file"),
    ],
)
def test_validate_match(processor, entry, method, struct_type):
    """A supported method whose pattern matches is valid."""
    assert processor.validate_pocket_method(entry, method, struct_type) is None


@pytest.mark.parametrize("entry, struct_type, expected", [c for c in INFERENCE_CASES if c.values[2] is not None])
def test_inferred_method_validates(processor, entry, struct_type, expected):
    """Every method inference gives passes validation, so the two tables agree."""
    assert processor.validate_pocket_method(entry, expected, struct_type) is None


# parse_residue_info


@pytest.mark.parametrize("residue_info", [None, ""])
def test_parse_residue_info_missing(processor, residue_info):
    """An absent residue list gives a reason."""
    result, reason = processor.parse_residue_info("4Q5J:A", residue_info)
    assert result is None
    assert "No residue ids" in reason


@pytest.mark.parametrize(
    "residue_info, expected, repeated",
    [
        ("10,11,12", "10,11,12", None),
        ("12,10,11", "12,10,11", None),
        ("07,7", "7", "7"),
        ("5,5,5,6,6", "5,6", "5,6"),
    ],
)
def test_parse_residue_info_canonical(processor, caplog, residue_info, expected, repeated):
    """Ids are made canonical and deduplicated in typed order, with one warning naming each repeat once."""
    caplog.set_level(logging.WARNING, logger="pocketmapper")
    result, reason = processor.parse_residue_info(f"4Q5J:A:{residue_info}", residue_info)
    assert (result, reason) == (expected, None)
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    if repeated is None:
        assert warnings == []
    else:
        assert len(warnings) == 1
        assert f"Residue id(s) {repeated} listed more than once" in warnings[0]


@pytest.mark.parametrize(
    "residue_info, bad_id",
    [("0", "0"), ("-1", "-1"), ("1.5", "1.5"), ("a", "a"), ("1,,2", ""), ("²", "²")],
)
def test_parse_residue_info_rejected(processor, residue_info, bad_id):
    """Anything but a positive decimal integer is rejected, naming the bad id."""
    result, reason = processor.parse_residue_info(f"4Q5J:A:{residue_info}", residue_info)
    assert result is None
    assert f"Residue id '{bad_id}'" in reason


# build_record / parse_individual_qt


def test_build_record_defaults(processor, work_dir):
    """A bare PDB id is an open search over the default chain."""
    record, reason = processor.parse_individual_qt("4Q5J", "auto")
    assert reason is None
    assert record.pocket_id == "4Q5J"
    assert record.struct_info == "4Q5J"
    assert record.chain_info == "A"
    assert record.residue_info is None
    assert record.struct_type == "pdb"
    assert record.struct_path == os.path.join(work_dir, "pdb", "4Q5J.cif.gz")
    assert record.pocket_method == "whole_chain"


def test_pocket_id_verbatim_residues_canonical(processor):
    """`pocket_id` keeps the entry as typed while `residue_info` is canonical."""
    record, reason = processor.parse_individual_qt("4Q5J:A:07,7", "auto")
    assert reason is None
    assert record.pocket_id == "4Q5J:A:07,7"
    assert record.residue_info == "7"
    assert record.pocket_method == "passthrough"


def test_preprocess_name_format(processor):
    """The name is `<stem>_<domain_chain>_<md5 of stem_chain>` for a fetched structure."""
    record, _ = processor.parse_individual_qt("4Q5J:B", "auto")
    assert re.fullmatch(r"4Q5J_B_[0-9a-f]{32}", record.preprocess_name)
    assert record.preprocess_name == f"4Q5J_B_{md5('4Q5J_B')}"


def test_preprocess_name_format_local(processor, write_file):
    """A local file's name hashes its contents too."""
    write_file("struct.cif", b"structure bytes")
    record, _ = processor.parse_individual_qt("struct.cif:A", "auto")
    assert record.preprocess_name == f"struct_A_{md5('struct_A_' + md5(b'structure bytes'))}"


def test_preprocess_name_uses_domain_chain(processor):
    """Two pockets on one chain share a preprocess name."""
    vdw, _ = processor.parse_individual_qt("4Q5J:A_E", "auto")
    passthrough, _ = processor.parse_individual_qt("4Q5J:A:1,2", "auto")
    assert vdw.preprocess_name == passthrough.preprocess_name


def test_preprocess_name_stem_stops_at_first_dot(processor, write_file):
    """The stem is the basename up to its first '.'."""
    write_file("4Q5J.cif.gz", b"structure bytes")
    record, _ = processor.parse_individual_qt("4Q5J.cif.gz:B_F", "auto")
    assert record.struct_type == "local_file"
    assert record.preprocess_name.startswith("4Q5J_B_")


def test_preprocess_name_local_contents(processor, write_file):
    """Same bytes share a name across paths; other bytes under the same basename, or an edit, change it."""

    def name(entry):
        """The preprocess name built for `entry`."""
        record, reason = processor.parse_individual_qt(entry, "auto")
        assert reason is None
        return record.preprocess_name

    write_file(os.path.join("dir1", "struct.cif"), b"first")
    write_file(os.path.join("dir2", "struct.cif"), b"first")
    write_file(os.path.join("dir3", "struct.cif"), b"second")
    original = name("dir1/struct.cif:A")
    assert name("dir2/struct.cif:A") == original
    assert name("dir3/struct.cif:A") != original

    write_file(os.path.join("dir1", "struct.cif"), b"edited")
    assert name("dir1/struct.cif:A") != original


def test_forced_method_mismatch(processor):
    """A forced method the entry cannot supply gives a reason."""
    record, reason = processor.parse_individual_qt("4Q5J:A", "passthrough")
    assert record is None
    assert reason is not None


def test_forced_method_match(processor):
    """A forced method the entry supplies gives a record with that method, not the one inference picks (pisa)."""
    record, reason = processor.parse_individual_qt("4Q5J:A_B", "vdw")
    assert reason is None
    assert record.pocket_method == "vdw"


@pytest.mark.parametrize(
    "entry, expected_path",
    [("human_domains", BUNDLED_HUMAN_DOMAINS_DB), ("pdb", None)],
)
def test_bundled_foldseek_db(processor, work_dir, entry, expected_path):
    """A bundled DB name is a DB record even under "auto", with no chain or preprocess name."""
    record, reason = processor.parse_individual_qt(entry, "auto")
    assert reason is None
    assert record.struct_type == "foldseek_db"
    assert record.struct_path == (expected_path or os.path.join(work_dir, "fsdb", "pdb"))
    assert record.chain_info is None
    assert record.preprocess_name is None


def test_forced_foldseek_db_relative_path(processor, work_dir):
    """A forced foldseek_db entry's path resolves against `work_dir`."""
    record, reason = processor.parse_individual_qt(os.path.join("dbs", "mydb"), "foldseek_db")
    assert reason is None
    assert record.struct_type == "foldseek_db"
    assert record.struct_path == os.path.join(work_dir, "dbs", "mydb")


def test_unparseable_entry(processor):
    """An entry with no recognisable structure gives a reason."""
    record, reason = processor.parse_individual_qt("XYZ:A", "auto")
    assert record is None
    assert "XYZ" in reason


def test_directory_entry_raises(processor, work_dir):
    """An entry naming a directory raises."""
    os.makedirs(os.path.join(work_dir, "somedir"))
    with pytest.raises(PocketMapperError, match="Directory input"):
        processor.parse_individual_qt("somedir:A", "auto")


# process_qt_cmdline_input


def test_missing_input_raises(processor):
    """No input raises."""
    with pytest.raises(PocketMapperError, match="query input is required"):
        processor.process_qt_cmdline_input(None, "query")


def test_unknown_method_raises(processor):
    """An unknown pocket method raises, listing every accepted one."""
    with pytest.raises(PocketMapperError) as excinfo:
        processor.process_qt_cmdline_input("4Q5J:A", "query", pocket_method="psia")
    assert "'psia'" in str(excinfo.value)
    assert "Choose one of: auto, whole_chain, pisa, passthrough, vdw, foldseek_db." in str(excinfo.value)


def test_single_entry_record_fields(processor):
    """A single entry gives one record dict holding exactly RECORD_FIELDS, in order."""
    records, rejected = processor.process_qt_cmdline_input("4Q5J:A", "query")
    assert rejected == []
    assert len(records) == 1
    assert tuple(records[0]) == RECORD_FIELDS


def test_entries_file(processor, write_file):
    """An entries file gives its records in file order, stripped, with a bad line rejected."""
    write_file(os.path.join("inputs", "entries.txt"), "4Q5J:A\n  P24941:A:1,2  \nXYZ:A\n4Q5J:B_F\n")
    records, rejected = processor.process_qt_cmdline_input(os.path.join("inputs", "entries.txt"), "target")
    assert [r["pocket_id"] for r in records] == ["4Q5J:A", "P24941:A:1,2", "4Q5J:B_F"]
    assert len(rejected) == 1
    entry, reason = rejected[0]
    assert entry == "XYZ:A"
    assert "XYZ" in reason


def test_entries_file_forced_method(processor, write_file):
    """A forced method rejects only the lines that cannot use it, and keeps the rest in file order."""
    path = write_file("entries.txt", "4Q5J:B_F\n4Q5J:A\n4Q5J:A_E\n")
    records, rejected = processor.process_qt_cmdline_input(path, "query", pocket_method="pisa")
    assert [(r["pocket_id"], r["pocket_method"]) for r in records] == [("4Q5J:B_F", "pisa"), ("4Q5J:A_E", "pisa")]
    assert len(rejected) == 1
    entry, reason = rejected[0]
    assert entry == "4Q5J:A"
    assert processor.pocket_methods["pisa"][1] in reason


def test_foldseek_db_file_not_read(processor, write_file, work_dir):
    """Under foldseek_db, an existing file is a DB, not an entries file."""
    write_file("mydb", "4Q5J:A\n4Q5J:B\n")
    records, rejected = processor.process_qt_cmdline_input("mydb", "target", pocket_method="foldseek_db")
    assert rejected == []
    assert len(records) == 1
    assert records[0]["pocket_id"] == "mydb"
    assert records[0]["struct_type"] == "foldseek_db"
    assert records[0]["struct_path"] == os.path.join(work_dir, "mydb")


def test_log_stage_names_side(processor):
    """The log stage names the side being processed."""
    processor.process_qt_cmdline_input("4Q5J:A", "query")
    assert processor.log_extra["stage"] == "Processing query"


def test_unreadable_entries_file_raises(processor, write_file, monkeypatch):
    """An entries file that cannot be read raises."""
    path = write_file("entries.txt", "4Q5J:A\n")
    real_open = builtins.open

    def failing_open(file, *args, **kwargs):
        """Raise OSError for the entries file, else open normally."""
        if file == path:
            raise OSError("permission denied")
        return real_open(file, *args, **kwargs)

    monkeypatch.setattr(builtins, "open", failing_open)
    with pytest.raises(PocketMapperError, match="Problem reading the file entries.txt"):
        processor.process_qt_cmdline_input("entries.txt", "query")
