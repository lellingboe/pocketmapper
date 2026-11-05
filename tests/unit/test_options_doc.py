"""
Checks that docs/options.md documents exactly the options each command takes.
"""

import os
import re

import pytest

from pocketmapper.lib.options import COMMANDS, command_options

OPTIONS_DOC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "docs", "options.md")

# Documented in prose rather than in a step's row: every step takes these, and every step after parse --entries
SHARED_STEP_FLAGS = ("--job_file", "--verbosity", "--log_path")
# Flags the registry does not hold
IGNORED_FLAGS = ("--help",)


def read_options_doc():
    """
    Read docs/options.md.

    Returns:
        str: Its text.
    """
    with open(OPTIONS_DOC, encoding="utf-8") as handle:
        return handle.read()


def section(text, heading):
    """
    The text under a heading, up to the next heading of the same or a higher level.

    Args:
        text (str): The options doc.
        heading (str): The heading line, e.g. "#### Cache options".

    Returns:
        str: The section body.
    """
    level = len(heading.split()[0])
    match = re.search(rf"^{re.escape(heading)}\n(.*?)(?=^#{{1,{level}}} |\Z)", text, re.S | re.M)
    assert match, f"options doc has no heading {heading!r}"
    return match.group(1)


def row_flags(cells):
    """
    The documented long flags and positionals in some text.

    Args:
        cells (str): A table row or cell.

    Returns:
        set of str: Each `--flag`, and `query` or `target` for `QUERY` or `TARGET`.
    """
    found = set(re.findall(r"`(--\w+)`", cells))
    found |= {word.lower() for word in re.findall(r"`(QUERY|TARGET)`", cells)}
    return found - set(IGNORED_FLAGS)


def table_flags(body):
    """
    The flags in the first column of every table row in some text.

    Args:
        body (str): A section of the options doc.

    Returns:
        set of str: The documented flags.
    """
    found = set()
    for line in body.splitlines():
        if line.startswith("|"):
            found |= row_flags(line.split("|")[1])
    return found


def documented_flags(command, text):
    """
    The flags the options doc documents for a command, with group links and shared prose expanded.

    Args:
        command (str): One of COMMANDS.
        text (str): The options doc.

    Returns:
        set of str: The documented flags.
    """
    if command == "search":
        return (
            table_flags(section(text, "### Arguments")) | table_flags(section(text, "### Options")) | group_flags(text)
        )
    step_table = section(text, "### Running one step at a time")
    # The step table has two: what each step does, then the options it takes
    options_table = step_table[step_table.index("| Command | Options |") :]
    row = next(line for line in options_table.splitlines() if line.startswith(f"| `{command}` |"))
    found = row_flags(row)
    # A link such as [cache options](#cache-options) stands for that group's table
    for anchor in re.findall(r"\(#([\w-]+-options)\)", row):
        found |= table_flags(section(text, f"#### {anchor.replace('-', ' ').capitalize()}"))
    found |= set(SHARED_STEP_FLAGS)
    if command != "parse":
        found.add("--entries")
    return found


def group_flags(text):
    """
    The flags in every `#### ... options` group.

    Args:
        text (str): The options doc.

    Returns:
        set of str: The documented flags.
    """
    found = set()
    for heading in re.findall(r"^(#### \w[\w ]* options)$", text, re.M):
        found |= table_flags(section(text, heading))
    return found


def registry_flags(command):
    """
    The flags a command takes, by the registry.

    Args:
        command (str): One of COMMANDS.

    Returns:
        set of str: Each long flag, and the positionals by name.
    """
    found = set()
    for option in command_options(command):
        flags = option.flags or (f"--{option.name}",)
        found |= {flag for flag in flags if flag.startswith("--") or not flag.startswith("-")}
    return found


@pytest.mark.parametrize("command", COMMANDS)
def test_options_doc_matches_registry(command):
    """docs/options.md documents exactly the flags `command_options(command)` has."""
    documented = documented_flags(command, read_options_doc())
    expected = registry_flags(command)
    assert documented == expected, (
        f"options doc vs registry for {command}: "
        f"missing from options doc {sorted(expected - documented)}, "
        f"not in registry {sorted(documented - expected)}"
    )
