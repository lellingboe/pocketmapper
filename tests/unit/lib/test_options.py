"""
Unit tests for `lib.options`: the registry's consistency, the tables derived from it, `option_help`, and
each validator factory.
"""

import re

import pytest

from pocketmapper.lib import options
from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.options import (
    CACHE_PATH_DEFAULTS,
    COMMANDS,
    GROUPS,
    JOB_KEYS,
    OPTIONS,
    RESULTS_PATH_DEFAULTS,
    SETTING_DEFAULTS,
    STEP_ONLY_KEYS,
    STEP_ONLY_PATH_KEYS,
    command_keys,
    command_options,
    job_key,
    option_help,
)

VALIDATED = [name for name, option in OPTIONS.items() if option.validate is not None]


def raises(message):
    """A `pytest.raises` context for a PocketMapperError whose message starts with `message`, matched literally."""
    return pytest.raises(PocketMapperError, match=f"^{re.escape(message)}")


def validate(name, value):
    """Run option `name`'s validator on `value`."""
    return OPTIONS[name].validate(value)


# The registry


def test_options_well_formed():
    """Every option names known commands and a known group, and is keyed by its name."""
    for name, option in OPTIONS.items():
        assert option.name == name
        assert option.commands and set(option.commands) <= set(COMMANDS), name
        assert option.group in GROUPS, name
        assert set(option.command_help or {}) <= set(option.commands), name


def test_under_unique():
    """No two options derive the same path, and each derives it from cache_dir or results_dir."""
    unders = [option.under for option in OPTIONS.values() if option.under]
    assert len(unders) == len(set(unders))
    assert {under[0] for under in unders} == {"cache_dir", "results_dir"}


def test_input_keys_are_job_keys():
    """An input option's key is the name of an option that sets it, so it is a job key."""
    for option in OPTIONS.values():
        if option.key:
            assert option.key in JOB_KEYS, option.name
            assert OPTIONS[option.key].key == "", option.name


@pytest.mark.parametrize("command", COMMANDS)
def test_command_options_nonempty(command):
    """Every command takes `job_file` and at least one setting."""
    names = [option.name for option in command_options(command)]
    assert "job_file" in names
    assert command_keys(command)


@pytest.mark.parametrize("command", COMMANDS)
def test_command_keys(command):
    """A command's job keys are its options' keys, each once, without `job_file`."""
    keys = command_keys(command)
    assert len(keys) == len(set(keys))
    assert set(keys) == {job_key(option) for option in command_options(command)} - {None}
    assert set(keys) <= set(JOB_KEYS)


# The tables derived from it, as they were hand-maintained before the registry (74ec5b9)


def test_derived_tables():
    """The derived defaults, path defaults and job keys hold what the hand-written ones did."""
    assert SETTING_DEFAULTS == {
        "cache_dir": "pocketmapper_cache",
        "verbosity": 3,
        "aligner": "foldseek",
        "align_count": 10,
        "align_struct_method": "auto",
        "query_pocket_method": "auto",
        "target_pocket_method": "auto",
        "delete_tmp": 1,
        "pisa_source": "ftp",
        "fetch_missing": 0,
        "offline": 0,
        "refresh_structures": 0,
    }
    assert CACHE_PATH_DEFAULTS == {
        "pdb_dir": "pdb_structures",
        "alphafold_dir": "alphafold_structures",
        "pocket_dir": "pockets",
        "foldseek_preprocessed_structure_dir": "foldseek_preprocessed_structures",
        "fsdb_dir": "fsdb",
    }
    assert RESULTS_PATH_DEFAULTS == {
        "temp_dir": "tmp",
        "aligned_structure_dir": "aligned_structures",
        "alignment_path": "alignment.tsv",
        "pocket_comparison_path": "pocket_comparison.tsv",
        "pockets_path": "pockets.json",
        "entries_path": "entries.tsv",
        "log_path": "info.log",
    }
    assert set(STEP_ONLY_KEYS) == {"pockets_tsv_path", "fetch_missing"}
    assert STEP_ONLY_PATH_KEYS == ("pockets_tsv_path",)
    assert set(JOB_KEYS) == set(command_keys("search")) | set(STEP_ONLY_KEYS)
    assert len(JOB_KEYS) == 32


# option_help


@pytest.mark.parametrize(
    "name, command, default",
    [
        ("verbosity", "align", "3"),
        ("cache_dir", "search", "pocketmapper_cache"),
        ("pdb_dir", "search", "<cache_dir>/pdb_structures"),
        ("log_path", "parse", "<results_dir>/info.log"),
        ("alignment_path", "search", "<results_dir>/alignment.tsv"),
        ("alignment_path", "align", "the job file's alignment_path"),
        ("temp_dir", "fetch_structures", "the job file's temp_dir"),
        ("log_path", "compare", "the job file's log_path"),
        ("pdb_dir", "superpose", "the job file's pdb_dir"),
        ("entries", "pockets", "the job file's entries_path"),
        ("threads", "align", "one per available core"),
        ("local_root", "superpose", "the job file's, else the current directory"),
        ("results_dir", "search", "pocketmapper_results_<YYMMDD_HHMMSS>"),
        ("results_dir", "parse", "none"),
        ("job_settings_path", "search", "<results_dir>/job_settings.json"),
        ("job_settings_path", "parse", "<results_dir>/parse_settings.json"),
    ],
)
def test_option_help_default(name, command, default):
    """The default text follows the static default, then the derived path, then the job file, then `default_help`."""
    assert option_help(OPTIONS[name], command).endswith(f" (default: {default})")


def test_option_help_command_help():
    """A command in `command_help` gets its own wording; the others get `help`."""
    option = OPTIONS["results_dir"]
    assert option_help(option, "search") == "Where results are written. (default: pocketmapper_results_<YYMMDD_HHMMSS>)"
    assert option_help(option, "parse").startswith("Where the log, the settings and the entries table are written.")


@pytest.mark.parametrize("name", ["query", "target"])
def test_option_help_no_default(name):
    """An option with no default to state gets its help alone."""
    assert option_help(OPTIONS[name], "search") == OPTIONS[name].help


# Validators: every one


@pytest.mark.parametrize("name", VALIDATED)
@pytest.mark.parametrize("value", [True, False])
def test_validate_rejects_bool(name, value):
    """No validator takes a bool, though True == 1 and False == 0."""
    with pytest.raises(PocketMapperError, match=f"^{re.escape(name)}|^Unknown {re.escape(name)}"):
        validate(name, value)


# choice: aligner, pisa_source, align_struct_method


@pytest.mark.parametrize(
    "name, value, expected",
    [
        ("aligner", "foldseek", "foldseek"),
        ("aligner", "seq", "seq"),
        ("aligner", "SEQ", "seq"),
        ("aligner", "Foldseek", "foldseek"),
        ("pisa_source", "ftp", "ftp"),
        ("pisa_source", "api", "api"),
        ("pisa_source", "FTP", "ftp"),
        ("pisa_source", "Api", "api"),
        ("align_struct_method", "AUTO", "auto"),
        ("align_struct_method", "Pocket", "pocket"),
    ],
)
def test_choice_accepts(name, value, expected):
    """A known value, in any case, is returned lowercased."""
    assert validate(name, value) == expected


@pytest.mark.parametrize(
    "name, value, choices",
    [
        ("aligner", "bogus", "foldseek, seq"),
        ("aligner", "", "foldseek, seq"),
        ("aligner", None, "foldseek, seq"),
        ("aligner", 1, "foldseek, seq"),
        ("pisa_source", "http", "ftp, api"),
        ("pisa_source", "", "ftp, api"),
        ("pisa_source", None, "ftp, api"),
        ("pisa_source", 1, "ftp, api"),
        ("align_struct_method", "bogus", "auto, pocket, foldseek"),
    ],
)
def test_choice_rejects(name, value, choices):
    """An unknown value, or one that is not a string, is rejected by name."""
    with raises(f"Unknown {name} {value!r}. Choose one of: {choices}."):
        validate(name, value)


# positive_int_or_cores: threads


@pytest.mark.parametrize("threads", [1, 2, 64])
def test_threads_accepts(threads):
    """A positive integer is returned unchanged."""
    assert validate("threads", threads) == threads


@pytest.mark.parametrize("cpu_count, expected", [(6, 6), (None, 1)])
def test_threads_unset(monkeypatch, cpu_count, expected):
    """None gives one thread per available core, or 1 when the core count is unknown."""
    monkeypatch.setattr(options.os, "cpu_count", lambda: cpu_count)
    assert validate("threads", None) == expected


@pytest.mark.parametrize("threads", [0, -1, True, False, 2.0, "2"])
def test_threads_rejects(threads):
    """Zero, a negative, a bool, a float or a string is rejected."""
    with raises(f"threads must be a positive integer, got {threads!r}."):
        validate("threads", threads)


# flag: delete_tmp, fetch_missing, offline, refresh_structures

FLAGS = ["delete_tmp", "fetch_missing", "offline", "refresh_structures"]


@pytest.mark.parametrize("name", FLAGS)
@pytest.mark.parametrize("value", [0, 1])
def test_flag_accepts(name, value):
    """0 and 1 are returned unchanged."""
    assert validate(name, value) == value


@pytest.mark.parametrize("name", FLAGS)
@pytest.mark.parametrize("value", [2, -1, True, False, "1", None])
def test_flag_rejects(name, value):
    """Anything but the int 0 or 1, a bool included, is rejected by setting name."""
    with raises(f"{name} must be 1 or 0, got {value!r}."):
        validate(name, value)


# percent_below_100: pisa_bsa_threshold


@pytest.mark.parametrize("threshold, expected", [(None, None), (0, 0.0), (50, 50.0), (25.5, 25.5), (99.9, 99.9)])
def test_pisa_bsa_threshold_accepts(threshold, expected):
    """None stays None; a number in [0, 100) is returned as a float."""
    result = validate("pisa_bsa_threshold", threshold)
    assert result == expected
    assert result is None or isinstance(result, float)


@pytest.mark.parametrize("threshold", [100, 100.0, 150, -1, -0.1, float("nan"), float("inf"), True, "50"])
def test_pisa_bsa_threshold_rejects(threshold):
    """100 and above, a negative, NaN, a bool or a string is rejected."""
    with raises(f"pisa_bsa_threshold must be a number at least 0 and below 100, got {threshold!r}."):
        validate("pisa_bsa_threshold", threshold)


# non_negative_int: align_count; integer: verbosity


@pytest.mark.parametrize("name, value", [("align_count", 0), ("align_count", 10), ("verbosity", 0), ("verbosity", -1)])
def test_int_accepts(name, value):
    """An integer in range is returned unchanged."""
    assert validate(name, value) == value


@pytest.mark.parametrize(
    "name, value, message",
    [
        ("align_count", -1, "align_count must be a non-negative integer, got -1."),
        ("align_count", 1.0, "align_count must be a non-negative integer, got 1.0."),
        ("align_count", None, "align_count must be a non-negative integer, got None."),
        ("align_count", "3", "align_count must be a non-negative integer, got '3'."),
        ("verbosity", 3.0, "verbosity must be an integer, got 3.0."),
        ("verbosity", None, "verbosity must be an integer, got None."),
        ("verbosity", "3", "verbosity must be an integer, got '3'."),
    ],
)
def test_int_rejects(name, value, message):
    """A float, None or a string is rejected, and so is a negative count."""
    with raises(message):
        validate(name, value)
