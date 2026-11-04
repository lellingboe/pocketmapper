"""
Unit tests for `lib.settings`: `resolve_align_struct_method`, the offline/refresh conflict, `layer_command`,
`validate_options`, `layer_settings`, `require_setting`, `read_job_file` and the path resolution. The
single-setting validators are tested in `test_options.py`.
"""

import json
import os
import re

import pytest

from pocketmapper.lib.constants import DEFAULT_ALIGN_COUNT, DEFAULT_ALIGNER
from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.options import command_options
from pocketmapper.lib.settings import (
    JOB_KEYS,
    SETTING_DEFAULTS,
    check_offline_refresh,
    layer_command,
    layer_settings,
    read_job_file,
    require_setting,
    resolve_align_struct_method,
    resolve_auto_align_struct_method,
    resolve_paths,
    resolve_step_paths,
    validate_options,
)
from pocketmapper.search import PocketMapper


def raises(message):
    """A `pytest.raises` context for a PocketMapperError whose message starts with `message`, matched literally."""
    return pytest.raises(PocketMapperError, match=f"^{re.escape(message)}")


@pytest.fixture
def job_source(request, write_file):
    """A function turning a job dict into the job file form `request.param` names: "dict" or "path"."""

    def make(job):
        """`job` itself, or the path of a JSON file holding it."""
        if request.param == "dict":
            return job
        return write_file("job.json", json.dumps(job))

    return make


# resolve_align_struct_method


@pytest.mark.parametrize(
    "method, aligner, expected",
    [
        ("auto", "foldseek", "foldseek"),
        ("auto", "seq", "pocket"),
        ("AUTO", "seq", "pocket"),
        ("pocket", "foldseek", "pocket"),
        ("pocket", "seq", "pocket"),
        ("foldseek", "foldseek", "foldseek"),
        ("Foldseek", "foldseek", "foldseek"),
    ],
)
def test_resolve_align_struct_method_accepts(method, aligner, expected):
    """A known method resolves to "pocket" or "foldseek", "auto" following the aligner."""
    assert resolve_align_struct_method(method, aligner) == expected


@pytest.mark.parametrize(
    "method, aligner, message",
    [
        ("foldseek", "seq", "align_struct_method 'foldseek' needs Foldseek's whole-chain transform"),
        ("FOLDSEEK", "seq", "align_struct_method 'foldseek' needs Foldseek's whole-chain transform"),
        ("bogus", "seq", "Unknown align_struct_method 'bogus'. Choose one of: auto, pocket, foldseek."),
        ("bogus", "foldseek", "Unknown align_struct_method 'bogus'."),
        (None, "foldseek", "Unknown align_struct_method None."),
        (1, "seq", "Unknown align_struct_method 1."),
    ],
)
def test_resolve_align_struct_method_rejects(method, aligner, message):
    """An unknown method is rejected, and so is "foldseek" with the seq aligner, which has no chain transform."""
    with raises(message):
        resolve_align_struct_method(method, aligner)


def test_resolve_auto_align_struct_method_takes_validated_value():
    """The resolution half leaves a concrete method as it is and follows the aligner for "auto"."""
    assert resolve_auto_align_struct_method("pocket", "foldseek") == "pocket"
    assert resolve_auto_align_struct_method("auto", "foldseek") == "foldseek"


# check_offline_refresh, alone and in `PocketMapper.resolve_settings`

OFFLINE_REFRESH_CONFLICT = "offline and refresh_structures cannot both be 1"


@pytest.mark.parametrize("offline, refresh_structures", [(0, 0), (1, 0), (0, 1)])
def test_check_offline_refresh_accepts(offline, refresh_structures):
    """Any pair but both 1 passes."""
    check_offline_refresh(offline, refresh_structures)


def test_check_offline_refresh_rejects():
    """Refreshing structures with downloads turned off is rejected."""
    with raises(OFFLINE_REFRESH_CONFLICT):
        check_offline_refresh(1, 1)


def test_resolve_settings_offline_refresh_conflict():
    """`resolve_settings` rejects offline and refresh_structures both 1 before the settings are dumped."""
    values = {
        "aligner": "seq",
        "align_struct_method": "auto",
        "threads": 1,
        "delete_tmp": 1,
        "pisa_source": "ftp",
        "pisa_bsa_threshold": None,
        "offline": 1,
        "refresh_structures": 1,
    }
    with raises(OFFLINE_REFRESH_CONFLICT):
        PocketMapper().resolve_settings(values)


# layer_command and validate_options


def parameter_names(command):
    """The parameter names `command`'s entry function must have, `job_file` aside."""
    return [option.name for option in command_options(command) if option.name != "job_file"]


def test_layer_command_maps_names_to_keys():
    """`job_file` and `self` are dropped, and an input option sets the job key it reads."""
    arguments = dict.fromkeys(parameter_names("compare")) | {"self": object(), "job_file": None}
    arguments["entries"] = "entries.tsv"
    values = layer_command("compare", None, arguments)
    assert values["entries_path"] == "entries.tsv"
    assert "entries" not in values
    assert "self" not in values


@pytest.mark.parametrize(
    "drop, add, message",
    [
        ("align_count", None, "missing align_count; not registered none."),
        (None, "bogus", "missing none; not registered bogus."),
    ],
)
def test_layer_command_rejects_drifted_signature(drop, add, message):
    """A parameter missing from the registry's options, or one it lacks, is rejected on the first call."""
    arguments = {name: None for name in parameter_names("search") if name != drop}
    if add is not None:
        arguments[add] = None
    with raises(f"The search parameters do not match its options in lib.options.OPTIONS: {message}"):
        layer_command("search", None, arguments | {"job_file": None})


@pytest.mark.parametrize(
    "key, value, message",
    [
        ("align_count", -1, "align_count must be a non-negative integer, got -1."),
        ("align_count", "10", "align_count must be a non-negative integer, got '10'."),
        ("verbosity", "3", "verbosity must be an integer, got '3'."),
        ("verbosity", True, "verbosity must be an integer, got True."),
    ],
)
def test_validate_options_rejects(key, value, message):
    """`validate_options` applies the validators of `align_count` and `verbosity`, which nothing checked before."""
    values = layer_settings(None, {key: value})
    with raises(message):
        validate_options(values, "search")


def test_validate_options_normalises_present_keys():
    """A present setting is normalised, an absent one is not added, and the input is left as it was."""
    values = {"aligner": "SEQ", "threads": 2}
    assert validate_options(values, "search") == {"aligner": "seq", "threads": 2}
    assert values["aligner"] == "SEQ"


# layer_settings and require_setting


def test_layer_settings_defaults():
    """With no job file and no arguments, every job key is set: its SETTING_DEFAULTS entry, else None."""
    values = layer_settings(None, {})
    assert set(values) == set(JOB_KEYS)
    for key in JOB_KEYS:
        assert values[key] == SETTING_DEFAULTS.get(key)


@pytest.mark.parametrize("job_source", ["dict", "path"], indirect=True)
def test_layer_settings_priority(job_source):
    """An argument beats the job file, which beats SETTING_DEFAULTS; an argument of None is unset."""
    job = job_source({"align_count": 3, "verbosity": 2, "threads": 4})
    values = layer_settings(job, {"align_count": 5, "verbosity": None, "aligner": None})
    assert values["align_count"] == 5
    assert values["verbosity"] == 2
    assert values["threads"] == 4
    assert values["aligner"] == DEFAULT_ALIGNER


@pytest.mark.parametrize("job_source", ["dict", "path"], indirect=True)
def test_layer_settings_job_file_over_default(job_source):
    """A job file value replaces a SETTING_DEFAULTS one when no argument is given."""
    values = layer_settings(job_source({"align_count": DEFAULT_ALIGN_COUNT + 1}), {"align_count": None})
    assert values["align_count"] == DEFAULT_ALIGN_COUNT + 1


@pytest.mark.parametrize("job_source", ["dict", "path"], indirect=True)
def test_layer_settings_drops_job_settings_path(job_source):
    """`job_settings_path` is never read from a job file, only from the arguments."""
    job = job_source({"job_settings_path": "from_job.json"})
    assert layer_settings(job, {})["job_settings_path"] is None
    assert layer_settings(job, {"job_settings_path": "from_arg.json"})["job_settings_path"] == "from_arg.json"


@pytest.mark.parametrize("job_source", ["dict", "path"], indirect=True)
def test_layer_settings_query_target_from_job_file(job_source):
    """A job file alone can supply query and target, and `require_setting` then accepts both."""
    job = job_source({"query": "4Q5J:A_E", "target": "4Q5J:B_F"})
    values = layer_settings(job, {"query": None, "target": None})
    assert (values["query"], values["target"]) == ("4Q5J:A_E", "4Q5J:B_F")
    for key in ("query", "target"):
        require_setting(values, key)


@pytest.mark.parametrize("job_source", ["dict", "path"], indirect=True)
@pytest.mark.parametrize("key", ["query", "target"])
def test_layer_settings_query_target_given_twice(job_source, key):
    """Query or target set both in the job file and as an argument is rejected, naming it."""
    job = job_source({"query": "4Q5J:A_E", "target": "4Q5J:B_F"})
    with raises(f"{key} is set both in the job file and as an argument; give it only once."):
        layer_settings(job, {key: "4Q5J:A_E"})


@pytest.mark.parametrize("key", ["query", "target"])
def test_require_setting_missing(key):
    """A setting with no default that neither the arguments nor the job file gave is rejected, naming it."""
    values = layer_settings({"query": None, "target": None}, {})
    with raises(f"No {key} given; pass it as an argument or set it in the job file."):
        require_setting(values, key)


# read_job_file


def test_read_job_file_none_and_dict():
    """None reads as an empty job; a dict is copied, not returned itself."""
    job = {"align_count": 3}
    assert read_job_file(None) == {}
    copy = read_job_file(job)
    assert copy == job
    assert copy is not job


def test_read_job_file_missing(local_root):
    """A path that does not exist is rejected, naming it."""
    path = os.path.join(local_root, "missing.json")
    with raises(f"Job file not found: {path}"):
        read_job_file(path)


@pytest.mark.parametrize(
    "content, message",
    [
        ("not json", "Error reading job file: {path}. Is it in JSON format?"),
        ("", "Error reading job file: {path}. Is it in JSON format?"),
        ('["align_count", 3]', 'Job file {path} must hold a JSON object of {{"option": value}}.'),
        ('{"bogus": 1}', "Unknown setting(s) in {path}: bogus"),
        ('{"pocket_method": "pisa"}', "Unknown setting(s) in {path}: pocket_method. No longer used: pocket_method"),
        (
            '{"work_dir": "."}',
            "Unknown setting(s) in {path}: work_dir. No longer used: work_dir "
            "(renamed to local_root; set local_root in the job file instead)",
        ),
    ],
)
def test_read_job_file_rejects(write_file, content, message):
    """A file that is not JSON, not a JSON object, or names an unknown or removed key is rejected."""
    path = write_file("job.json", content)
    with raises(message.format(path=path)):
        read_job_file(path)


def test_read_job_file_unknown_key_in_dict():
    """An unknown key in a job dict is rejected the same way as in a file."""
    with raises("Unknown setting(s) in the job settings: bogus"):
        read_job_file({"bogus": 1})


# resolve_paths, resolve_step_paths


@pytest.mark.parametrize("local_root_given", [False, True])
def test_resolve_paths_against_cwd(local_root, tmp_path, monkeypatch, local_root_given):
    """Every path resolves against the working directory whatever `local_root` is; `local_root` is absolute."""
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    monkeypatch.chdir(cwd)
    root = os.path.join(local_root, "root")
    values = layer_settings(None, {"local_root": root if local_root_given else None, "results_dir": "out"})
    values["pockets_tsv_path"] = "pockets.tsv"
    values["fsdb_dir"] = os.path.join("rel", "fsdb")
    resolved = resolve_paths(values, "parse")
    assert resolved["local_root"] == (root if local_root_given else str(cwd.resolve()))
    for key in (
        "results_dir",
        "cache_dir",
        "pdb_dir",
        "fsdb_dir",
        "entries_path",
        "job_settings_path",
        "pockets_tsv_path",
    ):
        assert resolved[key].startswith(str(cwd.resolve())), key
    assert resolved["fsdb_dir"] == os.path.join(str(cwd.resolve()), "rel", "fsdb")


def test_resolve_paths_relative_local_root(tmp_path, monkeypatch):
    """A relative `local_root` is made absolute against the working directory."""
    monkeypatch.chdir(tmp_path)
    resolved = resolve_paths(layer_settings(None, {"local_root": "data"}), "parse")
    assert resolved["local_root"] == os.path.join(str(tmp_path.resolve()), "data")


@pytest.mark.parametrize("local_root_given", [None, "data"])
def test_resolve_step_paths_local_root(tmp_path, monkeypatch, local_root_given):
    """An optional `local_root` is made absolute against the working directory, and is the cwd when unset."""
    monkeypatch.chdir(tmp_path)
    cwd = os.getcwd()
    values = layer_settings(None, {"local_root": local_root_given, "fsdb_dir": "fsdb"})
    resolved = resolve_step_paths(values, ("fsdb_dir",), ("local_root",))
    assert resolved["fsdb_dir"] == os.path.join(cwd, "fsdb")
    assert resolved["local_root"] == (os.path.join(cwd, "data") if local_root_given else cwd)
