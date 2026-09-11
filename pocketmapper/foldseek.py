"""
Everything PocketMapper does with the Foldseek binary and its bundled databases.

Foldseek is an optional external dependency. `run_foldseek` is the single point at which the package
runs a subcommand, so the binary name, the debug log and the failure convention are written once
(`check_foldseek`'s probe is the one other invocation, and reports rather than raises). Callers build
their own argument lists: the subcommands share no shape worth abstracting over.
"""

import logging
import os
import subprocess
from importlib.resources import files

from pocketmapper.exceptions import PocketMapperError

# The external binary, resolved off PATH. Optional, and never bundled.
FOLDSEEK_BINARY = "foldseek"


def check_foldseek():
    """
    Report whether the Foldseek binary is installed and runnable.

    Probes by actually running `foldseek -h`, which exits 0 without touching any input, rather than
    only resolving the name off PATH: a binary that is present but not executable, built for another
    architecture, or on a noexec mount resolves fine and then fails at the first real subcommand.
    Foldseek's output is discarded and nothing is raised; the reason for a False is logged at debug.

    Returns:
        bool: True if `foldseek -h` ran and exited 0.
    """
    try:
        subprocess.run(
            [FOLDSEEK_BINARY, "-h"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except (OSError, subprocess.CalledProcessError) as e:
        logging.debug(f"'{FOLDSEEK_BINARY}' is not runnable: {e}")
        return False
    return True


def run_foldseek(args, log_extra=None):
    """
    Run one Foldseek subcommand, logging the command first and raising on failure.

    Output is not captured, so Foldseek writes its own progress to stdout/stderr. Verbosity is
    Foldseek's to control, through whatever flags `args` carries.

    Args:
        args (list): The subcommand and its arguments, without the binary name; it is prepended here.
        log_extra (dict, optional): Logging `extra`, e.g. `{"stage": "Foldseek Alignment"}`. Defaults
            to None, which logs these records under this function's own name.

    Returns:
        None

    Raises:
        PocketMapperError: If Foldseek is not callable, or exits non-zero.
    """
    cmd = [FOLDSEEK_BINARY] + list(args)
    cmd_str = " ".join([str(x) for x in cmd])
    logging.debug(f"Running Foldseek with command: {cmd_str}", extra=log_extra)
    try:
        subprocess.run(cmd, check=True)
    except OSError as e:
        # Missing, not executable, or otherwise unable to exec -- all the same failure to the caller.
        msg = f"Foldseek is not callable, so '{cmd_str}' could not be run: {e}"
        logging.critical(msg, extra=log_extra)
        raise PocketMapperError(msg) from e
    except subprocess.CalledProcessError as e:
        msg = f"Foldseek exited with code {e.returncode} running '{cmd_str}'"
        logging.critical(msg, extra=log_extra)
        raise PocketMapperError(msg) from e


def bundled_human_domains_path(filename):
    """
    Resolve a file shipped in the package's `human_domains` directory.

    Anchored on the `pocketmapper` package rather than on `human_domains` itself: that directory holds
    only Foldseek DB files and has no __init__.py, so passing it to `files()` resolves a namespace
    package, which importlib.resources only learned to handle in 3.12. Going through the parent regular
    package works on every supported version.

    Args:
        filename (str): Name of the file within the `human_domains` directory.

    Returns:
        str: Absolute path to that file.
    """
    return str(files("pocketmapper").joinpath("human_domains", filename))


# The bundled human-domains Foldseek DB. Versioned here and nowhere else -- bump it on a DB refresh.
BUNDLED_HUMAN_DOMAINS_DB = bundled_human_domains_path("human_v3_20260901")

# The UniProt coordinates of every entry in the DB above, shipped beside it and keyed by the same
# entry names. Refresh the two together: an entry the table has lost aborts the run.
BUNDLED_HUMAN_DOMAINS_OFFSET_TABLE = bundled_human_domains_path("offset_table.tsv")


def bundled_foldseek_dbs(fsdb_dir):
    """
    The Foldseek database names accepted in place of a structure, mapped to what is known about each.

    Takes `fsdb_dir` rather than being a constant because `pdb` is downloaded into it on first use, so
    its path is not known until the cache directory is settled. `human_domains` ships inside the
    package and is fixed at import.

    Args:
        fsdb_dir (str): Cache directory for downloaded Foldseek databases.

    Returns:
        dict: DB name -> {"db_path": path to the database, "offset_path": path to the UniProt offset
            table shipped beside it, or None for a database that ships none}.
    """
    return {
        "human_domains": {
            "db_path": BUNDLED_HUMAN_DOMAINS_DB,
            "offset_path": BUNDLED_HUMAN_DOMAINS_OFFSET_TABLE,
        },
        "pdb": {
            "db_path": os.path.join(fsdb_dir, "pdb"),
            "offset_path": None,  # PDB hits carry real author ids; nothing to renumber
        },
    }
