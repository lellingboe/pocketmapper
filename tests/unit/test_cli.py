"""
Unit tests for `cli`: each command's options against the function it dispatches to, the job keys against
`Settings`, and an invalid option failing before any network access.
"""

import argparse
import inspect
import logging
import os
import socket
from dataclasses import fields

import pytest

from pocketmapper import cli
from pocketmapper import search as search_module
from pocketmapper.lib.constants import PACKAGE_LOGGER
from pocketmapper.lib.scopes import HELD_TEMP_DIRS, OUTER_COMMAND
from pocketmapper.lib.settings import JOB_KEYS, SETTING_DEFAULTS, STEP_ONLY_KEYS, Settings

SETTINGS_FIELDS = [field.name for field in fields(Settings)]


class NetworkAttempt(RuntimeError):
    """Raised in place of a connection or name lookup the test forbids."""


class ArgumentsCaptured(Exception):
    """Raised to stop `search` once its `arguments` dict has been recorded."""


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
    """The parser, `COMMANDS` and the dispatch table name the same commands."""
    assert set(subparsers()) == set(cli.COMMANDS) == set(cli.DISPATCH)


@pytest.mark.parametrize("command", list(cli.COMMANDS))
def test_options_match_signature(command):
    """A command's options are exactly the parameters of the function it dispatches to, each once."""
    dests = option_dests(subparsers()[command])
    assert len(dests) == len(set(dests)), f"duplicate options in {command}"
    assert set(dests) == set(parameters(cli.DISPATCH[command]))


def test_settings_match_search_signature():
    """`Settings` has a field for every parameter of `search` but `job_file`, and no other."""
    assert sorted(SETTINGS_FIELDS) == sorted(name for name in parameters(cli.DISPATCH["search"]) if name != "job_file")


def test_settings_match_search_arguments(monkeypatch):
    """The `arguments` dict `search` layers holds exactly the `Settings` fields."""
    captured = {}

    def capture(job_file, arguments):
        """Record `arguments` and stop the run."""
        captured.update(arguments)
        raise ArgumentsCaptured

    monkeypatch.setattr(search_module, "layer_settings", capture)
    with pytest.raises(ArgumentsCaptured):
        search_module.PocketMapper().search()
    assert sorted(captured) == sorted(SETTINGS_FIELDS)


def test_job_keys():
    """`JOB_KEYS` is the `Settings` fields plus `STEP_ONLY_KEYS`, and every `SETTING_DEFAULTS` key is a job key."""
    assert list(JOB_KEYS) == SETTINGS_FIELDS + list(STEP_ONLY_KEYS)
    assert set(SETTING_DEFAULTS) <= set(JOB_KEYS)


@pytest.mark.parametrize("command", list(cli.COMMANDS))
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
