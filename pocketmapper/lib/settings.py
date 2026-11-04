"""
Run configuration: the job file layered under the arguments, the validation of a command's options, the
paths, the settings dump, and the checks that span several settings.

Every command layers its settings the same way, in `layer_command` and `layer_settings`: an argument given
beats the job file, which beats `SETTING_DEFAULTS`. Search and parse then fill in the paths still unset under
`cache_dir` and `results_dir` (`resolve_paths`); every other step takes its paths resolved, as its job
file or arguments name them (`resolve_step_paths`). A job file may hold any of `JOB_KEYS`. Those, the
defaults and the single-setting validators that `validate_options` applies all come from the registry in
`lib.options`.

Every validator and `resolve_*` takes the value as given -- from the command line, a job file or a library
call, so any JSON value -- and returns it checked and normalised, or logs a critical and raises.
"""

import json
import logging
import os
from datetime import datetime

from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.options import (
    CACHE_PATH_DEFAULTS,
    JOB_KEYS,
    OPTIONS,
    RESULTS_PATH_DEFAULTS,
    SETTING_DEFAULTS,
    STEP_ONLY_PATH_KEYS,
    command_options,
    job_key,
)

logger = logging.getLogger(__name__)


# Job keys no longer accepted -> what replaced them, for the error a stale job file gets
REMOVED_JOB_KEYS = {
    key: "steps read their records from parse's entries.tsv; rerun parse for a new job file"
    for key in ("query_records_path", "target_records_path", "entries", "struct_type", "out_dir", "structures_tsv_path")
} | {
    "failed_entries_path": "entries.tsv records the entries each step dropped; rerun parse for a new job file",
    "pocket_method": "pockets builds the pocket method parse chose for each entry; rerun parse for a new job file",
    "work_dir": "renamed to local_root; set local_root in the job file instead",
}

# The registry's single-setting validators, by the names the commands call
resolve_aligner = OPTIONS["aligner"].validate
resolve_threads = OPTIONS["threads"].validate
resolve_delete_tmp = OPTIONS["delete_tmp"].validate
resolve_fetch_missing = OPTIONS["fetch_missing"].validate
resolve_offline = OPTIONS["offline"].validate
resolve_refresh_structures = OPTIONS["refresh_structures"].validate
resolve_pisa_source = OPTIONS["pisa_source"].validate
resolve_pisa_bsa_threshold = OPTIONS["pisa_bsa_threshold"].validate

# The settings dump search writes; parse writes parse_settings.json
SEARCH_SETTINGS_NAME = "job_settings.json"


def default_results_dir():
    """
    A fresh results directory name for a run that was given none.

    Returns:
        str: pocketmapper_results_<YYMMDD_HHMMSS>, relative to the working directory.
    """
    return f"pocketmapper_results_{datetime.now().strftime('%y%m%d_%H%M%S')}"


def read_job_file(job_file):
    """
    Read a job file: job key -> value.

    Args:
        job_file (str, dict or None): Path to a JSON job file, the same already loaded, or None for none.

    Returns:
        dict: A copy of the job's values; empty for None.

    Raises:
        PocketMapperError: If the file is missing, unreadable or not a JSON object, or the job names a
            key not in JOB_KEYS.
    """
    log_extra = {"stage": "Configuring Settings"}

    if job_file is None:
        return {}
    if isinstance(job_file, dict):
        job = dict(job_file)
        source = "the job settings"
    else:
        if not os.path.isfile(job_file):
            logger.critical(f"Job file not found: {job_file}", extra=log_extra)
            raise PocketMapperError(f"Job file not found: {job_file}")
        try:
            with open(job_file) as f:
                job = json.load(f)
        except Exception as e:
            logger.critical(f"Error reading job file: {job_file}. Is it in JSON format?", extra=log_extra)
            raise PocketMapperError(f"Error reading job file: {job_file}. Is it in JSON format?") from e
        if not isinstance(job, dict):
            msg = f'Job file {job_file} must hold a JSON object of {{"option": value}}.'
            logger.critical(msg, extra=log_extra)
            raise PocketMapperError(msg)
        source = job_file

    unknown = [key for key in job if key not in JOB_KEYS]
    if unknown:
        removed = [key for key in unknown if key in REMOVED_JOB_KEYS]
        msg = f"Unknown setting(s) in {source}: {', '.join(unknown)}"
        if removed:
            msg += ". No longer used: " + "; ".join(f"{key} ({REMOVED_JOB_KEYS[key]})" for key in removed)
        logger.critical(msg, extra=log_extra)
        raise PocketMapperError(msg)
    return job


def layer_settings(job_file, arguments):
    """
    Layer the arguments over the job file over SETTING_DEFAULTS.

    Args:
        job_file (str, dict or None): As `read_job_file`.
        arguments (dict): Job key -> the value passed, None for unset. Holds only the settings the
            caller takes as arguments.

    Returns:
        dict: Every job key -> its value: the argument if not None, else the job file's, else its
            SETTING_DEFAULTS entry, else None. `job_settings_path` comes only from the arguments.
            Nothing is checked or derived.

    Raises:
        PocketMapperError: If the job file cannot be read, or query or target is both an argument and
            in the job file.
    """
    log_extra = {"stage": "Configuring Settings"}

    job = read_job_file(job_file)
    # Never read from a job file: a step run from another command's dump would overwrite it
    job.pop("job_settings_path", None)
    # query and target must come from exactly one of the two: the one silently overridden would name
    # a search other than the one run
    for key in ("query", "target"):
        if arguments.get(key) is not None and job.get(key) is not None:
            msg = f"{key} is set both in the job file and as an argument; give it only once."
            logger.critical(msg, extra=log_extra)
            raise PocketMapperError(msg)

    values = dict.fromkeys(JOB_KEYS)
    values.update(SETTING_DEFAULTS)
    values.update(job)
    values.update({key: value for key, value in arguments.items() if value is not None})
    return values


def layer_command(command, job_file, arguments):
    """
    Layer an entry function's arguments over its job file, checking them against the command's options.

    Must be called as the entry function's first statement, with `locals()`, so `arguments` holds its
    parameters and nothing else.

    Args:
        command (str): One of `lib.options.COMMANDS`.
        job_file (str, dict or None): As `read_job_file`.
        arguments (dict): Parameter name -> the value passed, None for unset: the entry function's
            `locals()`. `job_file` and `self` are dropped.

    Returns:
        dict: As `layer_settings`, each parameter set under its job key (`entries` -> `entries_path`).

    Raises:
        PocketMapperError: If the parameters are not exactly the command's options but `job_file`, or as
            `layer_settings`.
    """
    log_extra = {"stage": "Configuring Settings"}

    names = {name: value for name, value in arguments.items() if name not in ("job_file", "self")}
    expected = [option.name for option in command_options(command) if option.name != "job_file"]
    missing = [name for name in expected if name not in names]
    extra = [name for name in names if name not in expected]
    if missing or extra:
        msg = (
            f"The {command} parameters do not match its options in lib.options.OPTIONS: "
            f"missing {', '.join(missing) or 'none'}; not registered {', '.join(extra) or 'none'}."
        )
        logger.critical(msg, extra=log_extra)
        raise PocketMapperError(msg)
    # None is left out, so an unset input option cannot hide another value of the key it shares
    keyed = {job_key(OPTIONS[name]): value for name, value in names.items() if value is not None}
    return layer_settings(job_file, keyed)


def validate_options(values, command):
    """
    Check and normalise each of a command's settings that has a validator, in registry order.

    Logs, as a validator may (the core count `threads` picks for None), so run it with the command's log
    open.

    Args:
        values (dict): Job key -> value, from `layer_command`. A key it lacks is not checked.
        command (str): One of `lib.options.COMMANDS`.

    Returns:
        dict: A copy of `values`, each checked setting replaced by its validated value.

    Raises:
        PocketMapperError: On the first setting its validator rejects.
    """
    values = dict(values)
    for option in command_options(command):
        key = job_key(option)
        if option.validate is not None and key in values:
            values[key] = option.validate(values[key])
    return values


def require_setting(values, key):
    """
    Check that a setting with no default was given.

    Args:
        values (dict): Job key -> value, from `layer_settings`.
        key (str): The setting.

    Returns:
        None

    Raises:
        PocketMapperError: If `values[key]` is None.
    """
    if values[key] is None:
        msg = f"No {key} given; pass it as an argument or set it in the job file."
        logger.critical(msg, extra={"stage": "Configuring Settings"})
        raise PocketMapperError(msg)


def resolve_paths(values, command):
    """
    Fill in `local_root`, `results_dir` and any derived path left unset, and make every path absolute.

    For search and parse, the commands that take `cache_dir` and `results_dir`.

    Args:
        values (dict): Job key -> value, with `cache_dir` set.
        command (str): The command these settings are for, which names its settings dump:
            SEARCH_SETTINGS_NAME for "search", else <command>_settings.json.

    Returns:
        dict: A copy of `values` with every path absolute. Every path resolves against the working
            directory, `local_root` included, which defaults to it. Paths already set are kept. For
            search, `results_dir` defaults to a timestamped name; for any other command it stays
            unset, and so does every results path not given, but `temp_dir`, which then defaults to
            <cache_dir>/tmp. The rest default to locations under `cache_dir` or `results_dir`.
    """
    values = dict(values)
    values["local_root"] = absolute_path(values["local_root"] if values["local_root"] is not None else os.getcwd())
    if values["results_dir"] is None and command == "search":
        values["results_dir"] = default_results_dir()
    results_dir = values["results_dir"] = absolute_path(values["results_dir"])
    values["cache_dir"] = absolute_path(values["cache_dir"])
    for key in CACHE_PATH_DEFAULTS:
        values[key] = absolute_path(cache_path(values["cache_dir"], key, values[key]))
    for key in RESULTS_PATH_DEFAULTS:
        if values[key] is None and results_dir is not None:
            values[key] = results_path(results_dir, key)
        values[key] = absolute_path(values[key])
    if values["job_settings_path"] is None and results_dir is not None:
        name = SEARCH_SETTINGS_NAME if command == "search" else f"{command}_settings.json"
        values["job_settings_path"] = os.path.join(results_dir, name)
    values["job_settings_path"] = absolute_path(values["job_settings_path"])
    if values["temp_dir"] is None:
        values["temp_dir"] = os.path.join(values["cache_dir"], "tmp")
    for key in STEP_ONLY_PATH_KEYS:
        values[key] = absolute_path(values[key])
    return values


def dump_settings(values):
    """
    Write a command's settings to its `job_settings_path`, as a job file any later command can take.

    Args:
        values (dict): Job key -> value. Written without `job_settings_path`, which a job file never
            sets. Nothing is written when `job_settings_path` is None.

    Returns:
        None

    Raises:
        PocketMapperError: If the file cannot be written.
    """
    log_extra = {"stage": "Configuring Settings"}

    path = values["job_settings_path"]
    if path is None:
        return
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            json.dump({key: value for key, value in values.items() if key != "job_settings_path"}, f, indent=4)
    except OSError as e:
        msg = f"Could not write the settings to {path}: {e}"
        logger.critical(msg, extra=log_extra)
        raise PocketMapperError(msg) from e
    logger.info(f"Settings written to {path}", extra=log_extra)


def absolute_path(path):
    """
    Resolve a path against the working directory.

    Args:
        path (str or None): The path; an absolute one is kept.

    Returns:
        str: The absolute, normalised path, or None for none.
    """
    return os.path.abspath(path) if path is not None else None


def resolve_step_paths(values, required, optional=()):
    """
    Make a step's paths absolute, checking that each one it needs is set.

    A step derives no path from `cache_dir` or `results_dir`: each comes resolved from its job file, such
    as parse's settings, or from its arguments. A relative one resolves against the working directory.
    An optional `local_root` left unset becomes the working directory.

    Args:
        values (dict): Job key -> value, from `layer_settings`.
        required (tuple): The path settings the step needs.
        optional (tuple, optional): The path settings it can do without. Defaults to none.

    Returns:
        dict: A copy of `values` with those paths absolute.

    Raises:
        PocketMapperError: If a required path is not set.
    """
    values = dict(values)
    for key in required:
        require_setting(values, key)
    for key in required + optional:
        if values[key] is not None:
            values[key] = os.path.abspath(values[key])
    if "local_root" in optional and values["local_root"] is None:
        values["local_root"] = os.getcwd()
    return values


def cache_path(cache_dir, key, path=None):
    """
    A cache path setting, or its default under `cache_dir` when unset.

    Args:
        cache_dir (str): The cache root.
        key (str): A key of CACHE_PATH_DEFAULTS, e.g. "pdb_dir".
        path (str, optional): The value given. Defaults to None, for unset.

    Returns:
        str: `path`, or the default.
    """
    return path if path is not None else os.path.join(cache_dir, CACHE_PATH_DEFAULTS[key])


def results_path(results_dir, key, path=None):
    """
    A results path setting, or its default under `results_dir` when unset.

    Args:
        results_dir (str): The results root.
        key (str): A key of RESULTS_PATH_DEFAULTS, e.g. "alignment_path".
        path (str, optional): The value given. Defaults to None, for unset.

    Returns:
        str: `path`, or the default.
    """
    return path if path is not None else os.path.join(results_dir, RESULTS_PATH_DEFAULTS[key])


def resolve_auto_align_struct_method(align_struct_method, aligner):
    """
    Turn a validated `align_struct_method` into "pocket" or "foldseek" for the aligner in use.

    Args:
        align_struct_method (str): "pocket", "foldseek" or "auto", already validated.
        aligner (str): The aligner that produced the alignment, already validated.

    Returns:
        str: "pocket" or "foldseek". "auto" gives "foldseek" with the foldseek aligner and
            "pocket" with the seq aligner, which produces no chain transform.

    Raises:
        PocketMapperError: If "foldseek" was asked for with the "seq" aligner.
    """
    log_extra = {"stage": "Configuring Settings"}

    method = align_struct_method
    if method == "auto":
        method = "foldseek" if aligner == "foldseek" else "pocket"
        logger.info(
            f"align_struct_method 'auto' resolved to '{method}' (aligner '{aligner}' is in use)",
            extra=log_extra,
        )
    elif method == "foldseek" and aligner != "foldseek":
        msg = (
            "align_struct_method 'foldseek' needs Foldseek's whole-chain transform, but this run "
            "uses the local BLOSUM62 aligner, which does not produce one. Use "
            "--align_struct_method pocket, or --aligner foldseek."
        )
        logger.critical(msg, extra=log_extra)
        raise PocketMapperError(msg)

    return method


def resolve_align_struct_method(align_struct_method, aligner):
    """
    Validate the tri-value `align_struct_method` setting and turn it into "pocket" or "foldseek".

    Args:
        align_struct_method (str): The requested method: "pocket", "foldseek" or "auto".
        aligner (str): The aligner that produced the alignment, already validated.

    Returns:
        str: As `resolve_auto_align_struct_method`.

    Raises:
        PocketMapperError: If the value is not one of ALIGN_STRUCT_METHODS, or as
            `resolve_auto_align_struct_method`.
    """
    return resolve_auto_align_struct_method(OPTIONS["align_struct_method"].validate(align_struct_method), aligner)


def check_offline_refresh(offline, refresh_structures):
    """
    Check that structures are not to be refreshed with downloads turned off.

    Args:
        offline (int): The `offline` setting, already validated.
        refresh_structures (int): The `refresh_structures` setting, already validated.

    Returns:
        None

    Raises:
        PocketMapperError: If both are 1.
    """
    log_extra = {"stage": "Configuring Settings"}

    if offline and refresh_structures:
        msg = "offline and refresh_structures cannot both be 1: refreshing structures needs downloads."
        logger.critical(msg, extra=log_extra)
        raise PocketMapperError(msg)


def check_fsdb_aligner(aligner):
    """
    Check that a Foldseek-database target is searched with the foldseek aligner.

    Args:
        aligner (str): The aligner, already validated.

    Returns:
        None

    Raises:
        PocketMapperError: If `aligner` is not "foldseek".
    """
    log_extra = {"stage": "Configuring Settings"}

    if aligner != "foldseek":
        msg = "A Foldseek database target requires --aligner foldseek."
        logger.critical(msg, extra=log_extra)
        raise PocketMapperError(msg)


def check_fsdb_align_struct_method(align_struct_method):
    """
    Check that a Foldseek-database target is not superposed on the pocket.

    Args:
        align_struct_method (str): The method as given, or resolved.

    Returns:
        None

    Raises:
        PocketMapperError: If `align_struct_method` is "pocket".
    """
    log_extra = {"stage": "Configuring Settings"}

    # Which kind of database it is is not known until its hits are read, so both kinds are rejected
    if align_struct_method == "pocket":
        msg = (
            "align_struct_method 'pocket' is not available against a Foldseek database "
            "target. A human_domains-style hit has no coordinates to superpose at all, and "
            "a PDB database's structures are assemblies while its pockets come from the "
            "wwPDB asymmetric unit, so a pocket fit would be applied in the wrong frame. "
            "Use --align_struct_method foldseek."
        )
        logger.critical(msg, extra=log_extra)
        raise PocketMapperError(msg)
