"""
Unit tests for `cli`: each command's options against the function it dispatches to, the job keys against
the registry, and an invalid option failing before any network access.
"""

import argparse
import inspect
import logging
import os
import socket

import pytest

from pocketmapper import cli
from pocketmapper.align.align import align
from pocketmapper.compare.compare import compare
from pocketmapper.fetch.fetch import fetch_structures
from pocketmapper.lib.constants import PACKAGE_LOGGER
from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.options import (
    COMMANDS,
    JOB_KEYS,
    OPTIONS,
    SETTING_DEFAULTS,
    STEP_ONLY_KEYS,
    command_keys,
    command_options,
    job_key,
)
from pocketmapper.lib.scopes import HELD_TEMP_DIRS, OUTER_COMMAND
from pocketmapper.parse.parse import parse
from pocketmapper.pockets.pockets import pockets
from pocketmapper.superpose.superpose import superpose


class NetworkAttempt(RuntimeError):
    """Raised in place of a connection or name lookup the test forbids."""


def subparsers():
    """Command -> its subparser, from `cli.build_parser()`."""
    parser = cli.build_parser()
    actions = [action for action in parser._actions if isinstance(action, argparse._SubParsersAction)]
    assert len(actions) == 1
    return actions[0].choices


def option_dests(subparser):
    """The dests of a subparser's options and positionals, without `help`."""
    return [action.dest for action in subparser._actions if not isinstance(action, argparse._HelpAction)]


def parameters(function):
    """The parameter names of `function`."""
    return list(inspect.signature(function).parameters)


def test_every_command_dispatched():
    """The parser, the registry's commands, the descriptions and the dispatch table name the same commands."""
    assert set(subparsers()) == set(COMMANDS) == set(cli.COMMAND_DESCRIPTIONS) == set(cli.DISPATCH)


@pytest.mark.parametrize("command", list(COMMANDS))
def test_options_match_signature(command):
    """A command's options are exactly the parameters of the function it dispatches to, each once."""
    dests = option_dests(subparsers()[command])
    assert len(dests) == len(set(dests)), f"duplicate options in {command}"
    assert set(dests) == set(parameters(cli.DISPATCH[command]))


@pytest.mark.parametrize("command", list(COMMANDS))
def test_entry_signature_matches_registry(command):
    """An entry function takes exactly its command's registry options, `job_file` and `self` aside."""
    names = [name for name in parameters(cli.DISPATCH[command]) if name not in ("job_file", "self")]
    assert sorted(names) == sorted(option.name for option in command_options(command) if option.name != "job_file")


def test_job_keys():
    """`JOB_KEYS` is search's keys plus `STEP_ONLY_KEYS`, each key once, and covers every default."""
    assert sorted(JOB_KEYS) == sorted(command_keys("search") + STEP_ONLY_KEYS)
    assert len(JOB_KEYS) == len(set(JOB_KEYS))
    assert set(SETTING_DEFAULTS) <= set(JOB_KEYS)
    # Every option's job key is a job key, except `job_file`, which has none
    assert {job_key(option) for option in OPTIONS.values()} - {None} == set(JOB_KEYS)


@pytest.mark.parametrize(
    "command, expected",
    [
        ("fetch_structures", {"pdb_dir", "alphafold_dir", "fsdb_dir", "local_root"}),
        ("align", {"pdb_dir", "alphafold_dir", "fsdb_dir", "local_root"}),
        ("superpose", {"pdb_dir", "alphafold_dir", "fsdb_dir", "local_root"}),
        ("compare", {"fsdb_dir", "local_root"}),
        ("pockets", {"pdb_dir", "alphafold_dir", "fsdb_dir", "local_root"}),
    ],
)
def test_step_path_options(command, expected):
    """A step takes the cache directories and `local_root` it needs to locate a record's structure."""
    assert expected <= set(option_dests(subparsers()[command]))


@pytest.mark.parametrize("command", list(COMMANDS))
def test_options_default_none(command):
    """Every option defaults to None, so the job file and SETTING_DEFAULTS apply under it."""
    subparser = subparsers()[command]
    for action in subparser._actions:
        if not isinstance(action, argparse._HelpAction):
            assert action.default is None, f"{command} --{action.dest} has a default"


def test_bogus_aligner_fails_before_network(local_root, monkeypatch, capsys):
    """`search --aligner bogus` exits 1 with the reason, having tried no connection or name lookup."""
    attempts = []

    def blocked_connect(self, address):
        """Record a connection attempt and refuse it."""
        attempts.append(("connect", address))
        raise NetworkAttempt(f"connect to {address}")

    def blocked_getaddrinfo(host, *args, **kwargs):
        """Record a name lookup and refuse it."""
        attempts.append(("getaddrinfo", host))
        raise NetworkAttempt(f"getaddrinfo for {host}")

    monkeypatch.setattr(socket.socket, "connect", blocked_connect)
    monkeypatch.setattr(socket, "getaddrinfo", blocked_getaddrinfo)
    package_logger = logging.getLogger(PACKAGE_LOGGER)
    handlers = list(package_logger.handlers)
    level = package_logger.level
    results_dir = os.path.join(local_root, "results")
    cache_dir = os.path.join(local_root, "cache")
    argv = ["search", "4Q5J:A_E", "4Q5J:B_F", "--aligner", "bogus", "--results_dir", results_dir]

    with pytest.raises(SystemExit) as exit_info:
        cli.cli(argv + ["--cache_dir", cache_dir])

    assert exit_info.value.code == 1
    assert attempts == []
    assert "Unknown aligner 'bogus'" in capsys.readouterr().out
    # Nothing was fetched or parsed: the cache was never made
    assert not os.path.exists(cache_dir)
    assert package_logger.handlers == handlers
    assert package_logger.level == level
    assert OUTER_COMMAND == []
    assert not HELD_TEMP_DIRS


@pytest.mark.parametrize(
    "step, setting",
    [
        (parse, {"verbosity": "loud"}),
        (fetch_structures, {"threads": 0}),
        (fetch_structures, {"delete_tmp": 2}),
        (align, {"aligner": "bogus"}),
        (align, {"threads": 0}),
        (pockets, {"pisa_source": "bogus"}),
        (compare, {"verbosity": "loud"}),
        (superpose, {"threads": 0}),
        (superpose, {"align_count": -1}),
        (superpose, {"align_struct_method": "bogus"}),
    ],
)
def test_step_validates_before_reading_entries(tmp_path, step, setting):
    """A step rejects a bad setting before it reads the entries table, even for an empty alignment."""
    alignment = tmp_path / "alignment.tsv"
    alignment.write_text("u\n")
    comparison = tmp_path / "comparison.tsv"
    comparison.write_text("")
    job = {
        "query": "4Q5J:A",
        "target": "4Q5J:B",
        "cache_dir": str(tmp_path / "cache"),
        "entries_path": str(tmp_path / "missing_entries.tsv"),
        "alignment_path": str(alignment),
        "pocket_comparison_path": str(comparison),
        "pockets_path": str(tmp_path / "pockets.json"),
        "aligned_structure_dir": str(tmp_path / "aligned"),
        "temp_dir": str(tmp_path / "tmp"),
        "pdb_dir": str(tmp_path / "pdb"),
        "alphafold_dir": str(tmp_path / "af"),
        "fsdb_dir": str(tmp_path / "fsdb"),
        "pocket_dir": str(tmp_path / "pockets"),
        "foldseek_preprocessed_structure_dir": str(tmp_path / "pre"),
    } | setting
    key = next(iter(setting))
    with pytest.raises(PocketMapperError, match=key):
        step(job_file=job)
