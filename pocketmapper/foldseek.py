"""
Everything PocketMapper does with the Foldseek binary and its bundled databases.

Foldseek is an optional external dependency. `run_foldseek` is the single point at which the package
shells out to it, so the binary name, the debug log and the failure convention are written once.
Callers build their own argument lists: the subcommands share no shape worth abstracting over.
"""

import logging
import os
import shutil
import subprocess
from importlib.resources import files

from pocketmapper.exceptions import PocketMapperError

# The external binary, resolved off PATH. Optional, and never bundled.
FOLDSEEK_BINARY = "foldseek"


def check_foldseek():
    """
    Report whether the Foldseek binary is on PATH.

    Returns:
        bool: True if `foldseek` is callable.
    """
    return shutil.which(FOLDSEEK_BINARY) is not None


def run_foldseek(args, stage):
    """
    Run one Foldseek subcommand, logging the command first and raising on failure.

    Output is not captured, so Foldseek writes its own progress to stdout/stderr. Verbosity is
    Foldseek's to control, through whatever flags `args` carries.

    Args:
        args (list): The subcommand and its arguments, without the binary name; it is prepended here.
        stage (dict): Logging `extra`, e.g. `{"stage": "Foldseek Alignment"}`. Required, because the
            root log format interpolates a `stage` key and a record without one fails to format.

    Returns:
        None

    Raises:
        PocketMapperError: If Foldseek is not callable, or exits non-zero.
    """
    cmd = [FOLDSEEK_BINARY] + list(args)
    cmd_str = " ".join([str(x) for x in cmd])
    logging.debug(f"Running Foldseek with command: {cmd_str}", extra=stage)
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError as e:
        msg = f"Foldseek is not callable, so '{cmd_str}' could not be run: {e}"
        logging.critical(msg, extra=stage)
        raise PocketMapperError(msg) from e
    except subprocess.CalledProcessError as e:
        msg = f"Foldseek exited with code {e.returncode} running '{cmd_str}'"
        logging.critical(msg, extra=stage)
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


def bundled_human_domains_offset_table(struct_path):
    """
    The offset table shipped beside the bundled human-domains DB, or None for any other DB.

    Only the bundled DB ships a table, so any other database gets None rather than a table that does
    not describe it. Compared on the resolved path, so naming the bundled DB by its path and naming it
    "human_domains" give the same answer.

    Args:
        struct_path (str): The Foldseek DB path in use.

    Returns:
        str | None: Path to the table, or None when `struct_path` is not the bundled DB.
    """
    if struct_path != BUNDLED_HUMAN_DOMAINS_DB:
        return None
    return BUNDLED_HUMAN_DOMAINS_OFFSET_TABLE


def bundled_foldseek_dbs(fsdb_dir):
    """
    The Foldseek database names accepted in place of a structure, mapped to their paths.

    Takes `fsdb_dir` rather than being a constant because `pdb` is downloaded into it on first use, so
    its path is not known until the cache directory is settled. `human_domains` ships inside the
    package and is fixed at import.

    Args:
        fsdb_dir (str): Cache directory for downloaded Foldseek databases.

    Returns:
        dict: DB name -> path.
    """
    return {
        "human_domains": BUNDLED_HUMAN_DOMAINS_DB,
        "pdb": os.path.join(fsdb_dir, "pdb"),
    }
