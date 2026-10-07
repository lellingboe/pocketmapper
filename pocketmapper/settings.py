"""
Run configuration: the `Settings` record and the validators that resolve each setting.

Every `resolve_*` takes the value as given -- from the command line, a job file or a library call, so
any JSON value -- and returns it checked and normalised, or logs a critical and raises.
`require_foldseek` is separate from `resolve_aligner`: some steps need the binary whatever aligner was
chosen.
"""

import logging
import os
from dataclasses import dataclass
from datetime import datetime

from pocketmapper.constants import ALIGN_STRUCT_METHODS
from pocketmapper.constants import ALIGNERS
from pocketmapper.constants import DELETE_TMP_VALUES
from pocketmapper.constants import FOLDSEEK_INSTALL_HINT
from pocketmapper.constants import PISA_SOURCES
from pocketmapper.exceptions import PocketMapperError
from pocketmapper.foldseek import check_foldseek

logger = logging.getLogger(__name__)


@dataclass
class Settings:
    """
    Fully resolved PocketMapper run configuration.

    A record of what a run actually used, built once at the end of `configure_workflow` and dumped
    to job_settings.json. It has no defaults: each field arrives from the job file, the arguments to
    `search()` or run-time resolution, in that priority order. No field is optional.
    """

    query: str
    target: str
    cache_dir: str
    results_dir: str
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
    query_records_path: str
    target_records_path: str
    pockets_path: str
    failed_entries_path: str
    job_settings_path: str
    log_path: str
    fsdb_dir: str


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
    "query_records_path": "query_records.json",
    "target_records_path": "target_records.json",
    "pockets_path": "pockets.json",
    "failed_entries_path": "failed_entries.json",
    "job_settings_path": "job_settings.json",
    "log_path": "info.log",
}


def default_results_dir():
    """
    A fresh results directory name for a run that was given none.

    Returns:
        str: pocketmapper_results_<YYMMDD_HHMMSS>, relative to the working directory.
    """
    return f"pocketmapper_results_{datetime.now().strftime('%y%m%d_%H%M%S')}"


def resolve_paths(values):
    """
    Fill in `results_dir` and any derived path left unset.

    Args:
        values (dict): Settings field name -> value, with `cache_dir` set.

    Returns:
        dict: A copy of `values` with every path set. Paths already set are kept; `results_dir`
            defaults to a timestamped name, the rest to locations under `cache_dir` or `results_dir`.
    """
    values = dict(values)
    if values["results_dir"] is None:
        values["results_dir"] = default_results_dir()
    for key in CACHE_PATH_DEFAULTS:
        values[key] = cache_path(values["cache_dir"], key, values[key])
    for key in RESULTS_PATH_DEFAULTS:
        values[key] = results_path(values["results_dir"], key, values[key])
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
