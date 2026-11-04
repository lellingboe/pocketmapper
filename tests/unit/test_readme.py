"""
Checks that the README documents exactly the options each command takes.
"""

import os
import re

import pytest

from pocketmapper.lib.options import COMMANDS, command_options

README = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "README.md")

# Documented in prose rather than in a step's row: every step takes these, and every step after parse --entries
SHARED_STEP_FLAGS = ("--job_file", "--verbosity", "--log_path")
# Flags the registry does not hold
IGNORED_FLAGS = ("--help",)


def read_readme():
    """
    Read the README.

    Returns:
        str: Its text.
    """
    with open(README, encoding="utf-8") as handle:
        return handle.read()


def section(text, heading):
    """
    The text under a heading, up to the next heading of the same or a higher level.

    Args:
        text (str): The README.
        heading (str): The heading line, e.g. "#### Cache options".

    Returns:
        str: The section body.
    """
    level = len(heading.split()[0])
    match = re.search(rf"^{re.escape(heading)}\n(.*?)(?=^#{{1,{level}}} |\Z)", text, re.S | re.M)
    assert match, f"README has no heading {heading!r}"
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
        body (str): A README section.

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
    The flags the README documents for a command, with group links and shared prose expanded.

    Args:
        command (str): One of COMMANDS.
        text (str): The README.

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
        text (str): The README.

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
def test_readme_options_match_registry(command):
    """The README documents exactly the flags `command_options(command)` has."""
    documented = documented_flags(command, read_readme())
    expected = registry_flags(command)
    assert documented == expected, (
        f"README vs registry for {command}: "
        f"missing from README {sorted(expected - documented)}, "
        f"not in registry {sorted(documented - expected)}"
    )
