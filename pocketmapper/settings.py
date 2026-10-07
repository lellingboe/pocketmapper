"""
Run configuration: the `Settings` record, the job file layered under the arguments, and the
validators that resolve each setting.

Every command layers its settings the same way, in `layer_settings`: an argument given beats the
job file, which beats `SETTING_DEFAULTS`; `resolve_paths` then fills in the paths still unset. A job
file may hold any of `JOB_KEYS`: the `Settings` fields, which are search's, plus the options only
some steps take.

Every `resolve_*` takes the value as given -- from the command line, a job file or a library call, so
any JSON value -- and returns it checked and normalised, or logs a critical and raises.
`require_foldseek` is separate from `resolve_aligner`: some steps need the binary whatever aligner was
chosen.
"""

import json
import logging
import os
from dataclasses import dataclass
from dataclasses import fields
from datetime import datetime

from pocketmapper.constants import ALIGN_STRUCT_METHODS
from pocketmapper.constants import ALIGNERS
from pocketmapper.constants import DEFAULT_ALIGN_COUNT
from pocketmapper.constants import DEFAULT_ALIGN_STRUCT_METHOD
from pocketmapper.constants import DEFAULT_ALIGNER
from pocketmapper.constants import DEFAULT_CACHE_DIR
from pocketmapper.constants import DEFAULT_DELETE_TMP
from pocketmapper.constants import DEFAULT_FETCH_MISSING
from pocketmapper.constants import DEFAULT_PISA_SOURCE
from pocketmapper.constants import DEFAULT_POCKET_METHOD
from pocketmapper.constants import DEFAULT_STRUCT_TYPE
from pocketmapper.constants import DEFAULT_VERBOSITY
from pocketmapper.constants import DELETE_TMP_VALUES
from pocketmapper.constants import FETCH_MISSING_VALUES
from pocketmapper.constants import FOLDSEEK_INSTALL_HINT
from pocketmapper.constants import PISA_SOURCES
from pocketmapper.exceptions import PocketMapperError
from pocketmapper.foldseek import check_foldseek

logger = logging.getLogger(__name__)


@dataclass
class Settings:
    """
    Fully resolved PocketMapper run configuration.

    A record of what a run actually used, built once by `resolve_settings` and dumped to
    job_settings.json. It has no defaults: each field arrives from the arguments to `search()`, the
    job file, SETTING_DEFAULTS or run-time resolution, in that priority order. No field is optional.
    """

    query: str
    target: str
    cache_dir: str
    results_dir: str
    # Absolute. Entries, a user Foldseek database and every relative path setting resolve against it.
    work_dir: str
    # A forced pocket method, or "auto" to infer one per entry. Unlike align_struct_method, "auto"
    # is kept here: it is resolved for each entry, not once for the run.
    query_pocket_method: str
    target_pocket_method: str
    # "foldseek" or "seq" (the local BLOSUM62 sequence aligner).
    aligner: str
    align_count: int
    # "pocket" or "foldseek"; "auto" is resolved to one of those before a Settings is built.
    align_struct_method: str
    verbosity: int
    threads: int
    # 1 deletes temp_dir at the end of the run; 0 keeps it for inspection.
    delete_tmp: int
    # "ftp" or "api": where PISA assembly interfaces are fetched from.
    pisa_source: str
    pdb_dir: str
    alphafold_dir: str
    pocket_dir: str
    foldseek_preprocessed_structure_dir: str
    temp_dir: str
    aligned_structure_dir: str
    alignment_path: str
    pocket_comparison_path: str
    pockets_path: str
    failed_entries_path: str
    job_settings_path: str
    log_path: str
    fsdb_dir: str


# What a job file may hold, in Settings order: every Settings field, plus the options only some steps
# take. `layer_settings` seeds every key.
STEP_ONLY_KEYS = (
    "entries",
    "struct_type",
    "out_dir",
    "structures_tsv_path",
    "pocket_method",
    "pockets_tsv_path",
    "fetch_missing",
)

# Step-only paths, resolved against work_dir like every path setting
STEP_ONLY_PATH_KEYS = ("out_dir", "structures_tsv_path", "pockets_tsv_path")
JOB_KEYS = tuple(field.name for field in fields(Settings)) + STEP_ONLY_KEYS

# Job keys no longer accepted -> what replaced them, for the error a stale job file gets
REMOVED_JOB_KEYS = {
    "query_records_path": "steps re-derive records from query and target; rerun parse for a new job file",
    "target_records_path": "steps re-derive records from query and target; rerun parse for a new job file",
}

# Setting -> its static default. The others default to None: resolved at run time, or required.
SETTING_DEFAULTS = {
    "cache_dir": DEFAULT_CACHE_DIR,
    "verbosity": DEFAULT_VERBOSITY,
    "aligner": DEFAULT_ALIGNER,
    "align_count": DEFAULT_ALIGN_COUNT,
    "align_struct_method": DEFAULT_ALIGN_STRUCT_METHOD,
    "query_pocket_method": DEFAULT_POCKET_METHOD,
    "target_pocket_method": DEFAULT_POCKET_METHOD,
    "pocket_method": DEFAULT_POCKET_METHOD,
    "delete_tmp": DEFAULT_DELETE_TMP,
    "pisa_source": DEFAULT_PISA_SOURCE,
    "fetch_missing": DEFAULT_FETCH_MISSING,
    "struct_type": DEFAULT_STRUCT_TYPE,
}

# Path setting -> its default, relative to cache_dir
CACHE_PATH_DEFAULTS = {
    "pdb_dir": "pdb_structures",
    "alphafold_dir": "alphafold_structures",
    "pocket_dir": "pockets",
    "foldseek_preprocessed_structure_dir": "foldseek_preprocessed_structures",
    "fsdb_dir": "fsdb",
}

# Path setting -> its default, relative to results_dir
RESULTS_PATH_DEFAULTS = {
    "temp_dir": "tmp",
    "aligned_structure_dir": "aligned_structures",
    "alignment_path": "alignment.tsv",
    "pocket_comparison_path": "pocket_comparison.tsv",
    "pockets_path": "pockets.json",
    "failed_entries_path": "failed_entries.json",
    "log_path": "info.log",
}

# The settings dump search writes; each step writes <command>_settings.json
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
    Fill in `work_dir`, `results_dir` and any derived path left unset, and make every path absolute.

    Args:
        values (dict): Job key -> value, with `cache_dir` set.
        command (str): The command these settings are for, which names its settings dump:
            SEARCH_SETTINGS_NAME for "search", else <command>_settings.json.

    Returns:
        dict: A copy of `values` with every path absolute. `work_dir` defaults to the working
            directory, and relative paths resolve against it. Paths already set are kept. For
            search, `results_dir` defaults to a timestamped name; for any other command it stays
            unset, and so does every results path not given, but `temp_dir`, which then defaults to
            <cache_dir>/tmp. The rest default to locations under `cache_dir` or `results_dir`.
    """
    values = dict(values)
    work_dir = values["work_dir"] = os.path.abspath(
        values["work_dir"] if values["work_dir"] is not None else os.getcwd()
    )
    if values["results_dir"] is None and command == "search":
        values["results_dir"] = default_results_dir()
    results_dir = values["results_dir"] = work_path(work_dir, values["results_dir"])
    values["cache_dir"] = work_path(work_dir, values["cache_dir"])
    for key in CACHE_PATH_DEFAULTS:
        values[key] = work_path(work_dir, cache_path(values["cache_dir"], key, values[key]))
    for key in RESULTS_PATH_DEFAULTS:
        if values[key] is None and results_dir is not None:
            values[key] = results_path(results_dir, key)
        values[key] = work_path(work_dir, values[key])
    if values["job_settings_path"] is None and results_dir is not None:
        name = SEARCH_SETTINGS_NAME if command == "search" else f"{command}_settings.json"
        values["job_settings_path"] = os.path.join(results_dir, name)
    values["job_settings_path"] = work_path(work_dir, values["job_settings_path"])
    if values["temp_dir"] is None:
        values["temp_dir"] = os.path.join(values["cache_dir"], "tmp")
    for key in STEP_ONLY_PATH_KEYS:
        values[key] = work_path(work_dir, values[key])
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


def work_path(work_dir, path):
    """
    Resolve a path against the working directory a run's settings name.

    Args:
        work_dir (str): The run's absolute `work_dir`.
        path (str or None): The path; an absolute one is kept.

    Returns:
        str: The absolute, normalised path, or None for none.
    """
    return os.path.normpath(os.path.join(work_dir, path)) if path is not None else None


def input_path(values, given, key):
    """
    The file a step reads: the one given, else the one a path setting names.

    Args:
        values (dict): Job key -> value, from `resolve_paths`.
        given (str or None): The input as passed to the step, None for unset. Resolved against
            `work_dir`.
        key (str): The path setting it defaults to, e.g. "alignment_path".

    Returns:
        str: The absolute path.
    """
    return work_path(values["work_dir"], given) if given is not None else values[key]


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


def resolve_aligner(aligner):
    """
    Validate the `aligner` setting.

    Args:
        aligner (str): The requested aligner, "foldseek" or "seq", in any case.

    Returns:
        str: The aligner, lowercased.

    Raises:
        PocketMapperError: If the value is not one of ALIGNERS.
    """
    log_extra = {"stage": "Configuring Settings"}

    # A job file can hold any JSON value, not only a str
    normalised = aligner.lower() if isinstance(aligner, str) else aligner
    if normalised not in ALIGNERS:
        msg = f"Unknown aligner {aligner!r}. Choose one of: {', '.join(ALIGNERS)}."
        logger.critical(msg, extra=log_extra)
        raise PocketMapperError(msg)
    return normalised


def require_foldseek(needed_by, alternative=None):
    """
    Check that the foldseek binary can run, by running it.

    Args:
        needed_by (str): What needs the binary, opening the error message, e.g. "The foldseek aligner
            was selected".
        alternative (str, optional): A way to avoid needing it, appended to the error message.
            Defaults to None, which suggests none.

    Returns:
        None

    Raises:
        PocketMapperError: If foldseek cannot be run.
    """
    log_extra = {"stage": "Configuring Settings"}

    if check_foldseek():
        return
    msg = (
        f"{needed_by} but 'foldseek' could not be run; it is either not on PATH or not executable "
        f"(run with --verbosity 4 for the reason). {FOLDSEEK_INSTALL_HINT}"
    )
    if alternative:
        msg += f" {alternative}"
    logger.critical(msg, extra=log_extra)
    raise PocketMapperError(msg)


def resolve_align_struct_method(align_struct_method, aligner):
    """
    Turn the tri-value `align_struct_method` setting into "pocket" or "foldseek".

    Args:
        align_struct_method (str): The requested method: "pocket", "foldseek" or "auto".
        aligner (str): The aligner that produced the alignment, already validated.

    Returns:
        str: "pocket" or "foldseek". "auto" gives "foldseek" with the foldseek aligner and
            "pocket" with the seq aligner, which produces no chain transform.

    Raises:
        PocketMapperError: If the value is not one of ALIGN_STRUCT_METHODS, or "foldseek" was
            asked for with the "seq" aligner.
    """
    log_extra = {"stage": "Configuring Settings"}

    # A job file can hold any JSON value, not only a str
    method = align_struct_method.lower() if isinstance(align_struct_method, str) else align_struct_method
    if method not in ALIGN_STRUCT_METHODS:
        msg = (
            f"Unknown align_struct_method {align_struct_method!r}. "
            f"Choose one of: {', '.join(ALIGN_STRUCT_METHODS)}."
        )
        logger.critical(msg, extra=log_extra)
        raise PocketMapperError(msg)

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


def resolve_threads(threads):
    """
    Turn an unset `threads` setting into a concrete core count.

    Args:
        threads (int or None): The requested core count.

    Returns:
        int: `threads`, or one per available core when it is None.

    Raises:
        PocketMapperError: If `threads` is not a positive integer.
    """
    log_extra = {"stage": "Configuring Settings"}

    if threads is None:
        # Not os.process_cpu_count() (3.13+) or os.sched_getaffinity (Linux only); None on an
        # exotic platform
        threads = os.cpu_count() or 1
        logger.info(f"threads unset, using one per available core ({threads})", extra=log_extra)
        return threads

    # A job file can hold any JSON value, and a bool passes isinstance(threads, int)
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        msg = f"threads must be a positive integer, got {threads!r}."
        logger.critical(msg, extra=log_extra)
        raise PocketMapperError(msg)

    return threads


def resolve_delete_tmp(delete_tmp):
    """
    Validate the `delete_tmp` setting.

    Args:
        delete_tmp (int): 1 to delete temp_dir at the end of the run, 0 to keep it.

    Returns:
        int: The value, unchanged.

    Raises:
        PocketMapperError: If the value is not one of DELETE_TMP_VALUES, including a bool.
    """
    log_extra = {"stage": "Configuring Settings"}

    # A job file can hold any JSON value, and True == 1
    if isinstance(delete_tmp, bool) or delete_tmp not in DELETE_TMP_VALUES:
        msg = f"delete_tmp must be 1 or 0, got {delete_tmp!r}."
        logger.critical(msg, extra=log_extra)
        raise PocketMapperError(msg)

    return delete_tmp


def resolve_fetch_missing(fetch_missing):
    """
    Validate the `fetch_missing` setting.

    Args:
        fetch_missing (int): 1 to download an entry structure missing from the cache, 0 to skip the
            entry.

    Returns:
        int: The value, unchanged.

    Raises:
        PocketMapperError: If the value is not one of FETCH_MISSING_VALUES, including a bool.
    """
    log_extra = {"stage": "Configuring Settings"}

    # A job file can hold any JSON value, and True == 1
    if isinstance(fetch_missing, bool) or fetch_missing not in FETCH_MISSING_VALUES:
        msg = f"fetch_missing must be 1 or 0, got {fetch_missing!r}."
        logger.critical(msg, extra=log_extra)
        raise PocketMapperError(msg)

    return fetch_missing


def resolve_pisa_source(pisa_source):
    """
    Validate the `pisa_source` setting.

    Args:
        pisa_source (str): "ftp" or "api", in any case.

    Returns:
        str: The source, lowercased.

    Raises:
        PocketMapperError: If the value is not one of PISA_SOURCES.
    """
    log_extra = {"stage": "Configuring Settings"}

    # A job file can hold any JSON value, not only a str
    normalised = pisa_source.lower() if isinstance(pisa_source, str) else pisa_source
    if normalised not in PISA_SOURCES:
        msg = f"Unknown pisa_source {pisa_source!r}. Choose one of: {', '.join(PISA_SOURCES)}."
        logger.critical(msg, extra=log_extra)
        raise PocketMapperError(msg)

    return normalised


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
