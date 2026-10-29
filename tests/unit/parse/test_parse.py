"""
Unit tests for the parse step: `parse_side`, `parse_entries` and the `parse` entry function.
"""

import json
import logging
import os

import pytest

from pocketmapper.lib.constants import PACKAGE_LOGGER
from pocketmapper.lib.entries import ENTRY_COLUMNS, RECORD_FIELDS, entry_sides, read_entries
from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.scopes import OUTER_COMMAND, run_scope
from pocketmapper.parse.parse import parse, parse_entries, parse_side


@pytest.fixture(autouse=True)
def no_leaked_state():
    """Check that a test leaves the package logger's handlers and level, and the run scope, as it found them."""
    package_logger = logging.getLogger(PACKAGE_LOGGER)
    handlers = list(package_logger.handlers)
    level = package_logger.level
    yield
    assert package_logger.handlers == handlers
    assert package_logger.level == level
    assert OUTER_COMMAND == []


@pytest.fixture
def warnings(caplog):
    """A function returning the messages of the package's warnings captured so far."""
    caplog.set_level(logging.WARNING, logger=PACKAGE_LOGGER)

    def messages():
        """The messages of the captured WARNING records logged under the package logger."""
        return [
            record.getMessage()
            for record in caplog.records
            if record.levelno == logging.WARNING and record.name.startswith(PACKAGE_LOGGER)
        ]

    return messages


@pytest.fixture
def run_parse(cache_dirs, local_root):
    """A function calling `parse` with the cache directories and working directory under `tmp_path`."""

    def run(**kwargs):
        """Call `parse` with `cache_dir`, the cache paths and `local_root` under `tmp_path`, plus `kwargs`."""
        return parse(cache_dir=os.path.join(local_root, "cache"), local_root=local_root, **cache_dirs, **kwargs)

    return run


# parse_side


@pytest.mark.parametrize(
    "entries, expected",
    [
        ("4Q5J:A", [("4Q5J:A", "4Q5J:A")]),
        (["4Q5J:A", "P24941"], [("4Q5J:A", "4Q5J:A"), ("P24941", "P24941")]),
    ],
)
def test_parse_side_sources(local_root, entries, expected):
    """A string and a list of sources each give a row per entry, with `source` set to the input it came from."""
    rows = parse_side(entries, "query", "auto", local_root)
    assert [(row["pocket_id"], row["source"]) for row in rows] == expected
    assert all(row["side"] == "query" for row in rows)


def test_parse_side_order(local_root, write_file, monkeypatch):
    """Rows follow source order, and within a source the order of its entries, parsed or rejected."""
    monkeypatch.chdir(local_root)
    write_file("first.txt", "XYZ\n4Q5J:A\n1abc:B\n")
    write_file("second.txt", "4Q5J:AB\nP24941\n")
    rows = parse_side(["first.txt", "second.txt", "2ER6"], "target", "auto", local_root)
    assert [(row["pocket_id"], row["source"], row["failed_reason"]) for row in rows] == [
        ("XYZ", "first.txt", "invalid_entry"),
        ("4Q5J:A", "first.txt", None),
        ("1abc:B", "first.txt", None),
        ("4Q5J:AB", "second.txt", "invalid_entry"),
        ("P24941", "second.txt", None),
        ("2ER6", "2ER6", None),
    ]


def test_parse_side_rejected_row(local_root):
    """A rejected row is marked dropped by parse as `invalid_entry`, holding only `pocket_id` of the record fields."""
    (row,) = parse_side("XYZ", "query", "auto", local_root)
    assert row["failed_step"] == "parse"
    assert row["failed_reason"] == "invalid_entry"
    assert row["failed_detail"] == "Could not determine structure type for XYZ"
    assert row["pocket_id"] == "XYZ"
    assert all(row[key] is None for key in RECORD_FIELDS if key != "pocket_id")


def test_parse_side_columns(local_root, write_file):
    """Every row, parsed or rejected, has exactly the entries table's columns, in order."""
    write_file("struct.cif", b"any bytes")
    rows = parse_side(["4Q5J:A_E", "struct.cif:A:10,11", "human_domains", "XYZ"], "target", "auto", local_root)
    assert len(rows) == 4
    assert all(tuple(row) == ENTRY_COLUMNS for row in rows)


# parse_entries


def test_parse_entries_normal(local_root):
    """A valid search returns the rows and the sides `entry_sides` builds from them."""
    sides, rows = parse_entries(["4Q5J:A_E", "XYZ"], "4Q5J:B_F", "auto", "auto", local_root)
    assert sides == entry_sides(rows)
    assert [record["pocket_id"] for record in sides["query"]] == ["4Q5J:A_E"]
    assert [record["pocket_id"] for record in sides["target"]] == ["4Q5J:B_F"]
    assert [row["pocket_id"] for row in rows] == ["4Q5J:A_E", "XYZ", "4Q5J:B_F"]


def test_parse_entries_fsdb_query_dropped(local_root):
    """A Foldseek-database query entry is marked `invalid_entry`, and the other query entries are kept."""
    sides, rows = parse_entries(["human_domains", "4Q5J:A_E"], "4Q5J:B_F", "auto", "auto", local_root)
    (row,) = [row for row in rows if row["pocket_id"] == "human_domains"]
    assert row["failed_step"] == "parse"
    assert row["failed_reason"] == "invalid_entry"
    assert row["failed_detail"] == "A Foldseek database cannot be a query entry: human_domains"
    assert [record["pocket_id"] for record in sides["query"]] == ["4Q5J:A_E"]


@pytest.mark.parametrize(
    "query, target, match",
    [
        ("human_domains", "4Q5J:B_F", "^no valid query entries$"),
        ("XYZ", "4Q5J:B_F", "^no valid query entries$"),
        ("4Q5J:A_E", "XYZ", "^no valid target entries$"),
        ("XYZ", "ABC", "^no valid query entries; no valid target entries$"),
    ],
)
def test_parse_entries_no_valid_side(local_root, query, target, match):
    """A side left with no valid entry raises, one error naming every such side."""
    with pytest.raises(PocketMapperError, match=match):
        parse_entries(query, target, "auto", "auto", local_root)


@pytest.mark.parametrize(
    "target, target_pocket_method",
    [
        (["4Q5J:B_F", "human_domains"], "auto"),
        (["human_domains", "4Q5J:B_F"], "auto"),
        (["user_db", "other_db"], "foldseek_db"),
    ],
)
def test_parse_entries_fsdb_target_not_alone(local_root, target, target_pocket_method):
    """A Foldseek-database target beside any other kept target entry raises."""
    with pytest.raises(PocketMapperError, match="must be the only target entry, but 2 target entries"):
        parse_entries("4Q5J:A_E", target, "auto", target_pocket_method, local_root)


@pytest.mark.parametrize(
    "target, target_pocket_method",
    [("human_domains", "auto"), ("user_db", "foldseek_db"), (["human_domains", "XYZ"], "auto")],
)
def test_parse_entries_fsdb_target_alone(local_root, target, target_pocket_method):
    """
    A Foldseek-database target on its own passes.

    A dropped entry beside it does not count.
    """
    sides, _ = parse_entries("4Q5J:A_E", target, "auto", target_pocket_method, local_root)
    (record,) = sides["target"]
    assert record["struct_type"] == "foldseek_db"


@pytest.mark.parametrize(
    "query, target, raises",
    [
        (["4Q5J:A_E", "XYZ", "human_domains"], "4Q5J:B_F", False),
        ("XYZ", "4Q5J:B_F", True),
        ("4Q5J:A_E", ["human_domains", "4Q5J:B_F"], True),
    ],
)
def test_parse_entries_writes_table(local_root, query, target, raises):
    """The entries table is written with every row, before the sides are checked, so even when that check raises."""
    entries_path = os.path.join(local_root, "out", "entries.tsv")
    if raises:
        with pytest.raises(PocketMapperError):
            parse_entries(query, target, "auto", "auto", local_root, entries_path)
        rows = parse_side(query, "query", "auto", local_root)
        rows += parse_side(target, "target", "auto", local_root)
    else:
        _, rows = parse_entries(query, target, "auto", "auto", local_root, entries_path)
    assert read_entries(entries_path) == rows


def test_parse_entries_no_table(local_root):
    """Without an `entries_path`, nothing is written."""
    parse_entries("4Q5J:A_E", "4Q5J:B_F", "auto", "auto", local_root)
    assert os.listdir(local_root) == []


def test_parse_entries_warns_per_dropped_row(local_root, warnings):
    """Each dropped row logs one warning, its detail ending `; skipping this entry`."""
    _, rows = parse_entries(["XYZ", "human_domains", "4Q5J:A_E"], ["4Q5J:AB", "4Q5J:B_F"], "auto", "auto", local_root)
    dropped = [row for row in rows if row["failed_step"] is not None]
    assert [row["pocket_id"] for row in dropped] == ["XYZ", "human_domains", "4Q5J:AB"]
    assert warnings() == [f"{row['failed_detail']}; skipping this entry" for row in dropped]


# parse


def test_parse_results_dir(run_parse, local_root):
    """With a `results_dir`, parse writes the entries table, its settings and its log, and returns the table."""
    results_dir = os.path.join(local_root, "results")
    table = run_parse(query="4Q5J:A_E", target=["4Q5J:B_F", "XYZ"], results_dir=results_dir)
    assert sorted(os.listdir(results_dir)) == ["entries.tsv", "info.log", "parse_settings.json"]
    assert list(table.columns) == list(ENTRY_COLUMNS)
    assert list(table["pocket_id"]) == ["4Q5J:A_E", "4Q5J:B_F", "XYZ"]
    assert len(read_entries(os.path.join(results_dir, "entries.tsv"))) == 3
    with open(os.path.join(results_dir, "parse_settings.json")) as f:
        settings = json.load(f)
    assert settings["query"] == "4Q5J:A_E"
    assert settings["entries_path"] == os.path.join(results_dir, "entries.tsv")
    assert "job_settings_path" not in settings


def test_parse_no_results_dir(run_parse, local_root, write_file):
    """Without a `results_dir`, parse writes no file."""
    write_file("struct.cif", b"any bytes")
    table = run_parse(query="struct.cif:A_B", target="4Q5J:B_F")
    assert list(table.columns) == list(ENTRY_COLUMNS)
    assert os.listdir(local_root) == ["struct.cif"]


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"target": "4Q5J:B_F"}, "No query given"),
        ({"query": "4Q5J:A_E"}, "No target given"),
        ({"query": "4Q5J:A_E", "target": "4Q5J:B_F", "job_file": {"query": "1abc"}}, "query is set both"),
        ({"query": "4Q5J:A_E", "target": "4Q5J:B_F", "job_file": {"target": "1abc"}}, "target is set both"),
    ],
)
def test_parse_query_target_errors(run_parse, local_root, kwargs, match):
    """A query or target given neither way, or both in the job file and as an argument, raises before writing."""
    with pytest.raises(PocketMapperError, match=match):
        run_parse(results_dir=os.path.join(local_root, "results"), **kwargs)
    assert os.listdir(local_root) == []


def test_parse_job_file_dict(run_parse):
    """A job file given as a dict supplies the query and target."""
    table = run_parse(job_file={"query": "4Q5J:A_E", "target": "4Q5J:B_F"})
    assert list(zip(table["side"], table["pocket_id"])) == [("query", "4Q5J:A_E"), ("target", "4Q5J:B_F")]


@pytest.mark.parametrize(
    "kwargs, expected",
    [
        ({"query_pocket_method": "vdw"}, [("query", "vdw"), ("target", "pisa")]),
        ({"target_pocket_method": "vdw"}, [("query", "pisa"), ("target", "vdw")]),
    ],
)
def test_parse_pocket_methods(run_parse, kwargs, expected):
    """Each side's forced pocket method reaches only that side's entries."""
    table = run_parse(query="4Q5J:A_B", target="4Q5J:A_B", **kwargs)
    assert list(zip(table["side"], table["pocket_method"])) == expected


def test_parse_nested_in_search(run_parse, local_root):
    """Inside an open search scope, parse writes its entries table and log but no settings dump."""
    results_dir = os.path.join(local_root, "results")
    with run_scope("search"):
        run_parse(query="4Q5J:A_E", target="4Q5J:B_F", results_dir=results_dir)
    assert sorted(os.listdir(results_dir)) == ["entries.tsv", "info.log"]


def test_read_entries_ignores_struct_path_column(local_root):
    """An entries table written with the former `struct_path` column still reads, the column left out."""
    _, rows = parse_entries("4Q5J:A_E", "4Q5J:B_F", "auto", "auto", local_root)
    path = os.path.join(local_root, "old_entries.tsv")
    with open(path, "w") as f:
        f.write("\t".join(list(ENTRY_COLUMNS) + ["struct_path"]) + "\n")
        for row in rows:
            f.write("\t".join(row[column] or "" for column in ENTRY_COLUMNS) + "\t/old/path.cif.gz\n")
    assert read_entries(path) == rows
