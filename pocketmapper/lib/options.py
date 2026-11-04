"""
The option registry: every option once, with its commands, help, default and validator.

`OPTIONS` is the single source the per-command lists, the static defaults, the job keys and the path
defaults are derived from. Its order is the `--help` order within each group of `GROUPS`, and the order
of `JOB_KEYS`. A job key is the option's `key`, or its name when that is empty; `job_file` has none.

Each `validate` takes the value as given -- from the command line, a job file or a library call, so any
JSON value -- and returns it checked and normalised, or logs a critical and raises.

This module imports only `lib.constants` and `lib.exceptions` from the package, so every other module
can import it.
"""

import logging
import os
from typing import NamedTuple

from pocketmapper.lib.constants import (
    ALIGN_STRUCT_METHODS,
    ALIGNERS,
    DEFAULT_ALIGN_COUNT,
    DEFAULT_ALIGN_STRUCT_METHOD,
    DEFAULT_ALIGNER,
    DEFAULT_CACHE_DIR,
    DEFAULT_DELETE_TMP,
    DEFAULT_FETCH_MISSING,
    DEFAULT_OFFLINE,
    DEFAULT_PISA_SOURCE,
    DEFAULT_POCKET_METHOD,
    DEFAULT_REFRESH_STRUCTURES,
    DEFAULT_VERBOSITY,
    FLAG_VALUES,
    PISA_SOURCES,
)
from pocketmapper.lib.exceptions import PocketMapperError

logger = logging.getLogger(__name__)

# The commands, in `--help` order; every one but search is a step
STEPS = ("parse", "fetch_structures", "align", "pockets", "compare", "superpose")
COMMANDS = STEPS + ("search",)

# The `--help` groups, in the order they are printed. "main" options are the parser's own; argparse prints
# them first whatever the order.
GROUPS = ("main", "in", "aligned structure", "cache", "out", "temp", "advanced")

# The metavars of a path option
PATH_METAVARS = ("PATH", "DIR")


class Option(NamedTuple):
    """
    One option: a parameter of each command in `commands`, its `--help` entry and its job key.
    """

    # Parameter name and argparse dest
    name: str
    # The commands that take it, a subset of COMMANDS
    commands: tuple
    # Help without the "(default: ...)" part, which `option_help` adds
    help: str
    # Static default -> SETTING_DEFAULTS; None for one resolved at run time, or none
    default: object = None
    # value -> value, logs a critical and raises PocketMapperError; None for no check here
    validate: object = None
    # Job key when it differs from `name`: the name of the output option an input option reads ("entries" ->
    # "entries_path"); "" for `name`; None for no job key
    key: str = ""
    # One of GROUPS
    group: str = "main"
    # Every flag, e.g. ("-a", "--aligner"); ("query",) is a positional; () for ("--<name>",)
    flags: tuple = ()
    # argparse type, e.g. int or float; None for str
    type: object = None
    metavar: str = "STR"
    # (root setting, path under it): the default of a path derived from cache_dir or results_dir
    under: tuple = ()
    # Default text when there is no static or derived default, e.g. "one per available core"
    default_help: str = ""
    # Command -> (help, default text), where both differ from `help` and the default rules
    command_help: dict = None


def reject(message):
    """
    Log a bad setting as critical and raise.

    Args:
        message (str): What is wrong.

    Raises:
        PocketMapperError: Always, with `message`.
    """
    logger.critical(message, extra={"stage": "Configuring Settings"})
    raise PocketMapperError(message)


def flag(name):
    """
    A validator for a setting that is 1 or 0.

    Args:
        name (str): The setting, for the message.

    Returns:
        function: value -> the value, unchanged.
    """

    def validate(value):
        """
        Check that `value` is 1 or 0.

        Args:
            value (object): The value given.

        Returns:
            int: The value, unchanged.

        Raises:
            PocketMapperError: If the value is not one of FLAG_VALUES, including a bool.
        """
        # A job file can hold any JSON value, and True == 1
        if isinstance(value, bool) or value not in FLAG_VALUES:
            reject(f"{name} must be 1 or 0, got {value!r}.")
        return value

    return validate


def choice(name, values):
    """
    A validator for a setting that is one of a few strings, in any case.

    Args:
        name (str): The setting, for the message.
        values (tuple): The accepted values, lowercase.

    Returns:
        function: value -> the value, lowercased.
    """

    def validate(value):
        """
        Check that `value` is one of the accepted values, in any case.

        Args:
            value (object): The value given.

        Returns:
            str: The value, lowercased.

        Raises:
            PocketMapperError: If the value is not one of them, or not a string.
        """
        # A job file can hold any JSON value, not only a str
        normalised = value.lower() if isinstance(value, str) else value
        if normalised not in values:
            reject(f"Unknown {name} {value!r}. Choose one of: {', '.join(values)}.")
        return normalised

    return validate


def positive_int_or_cores(name):
    """
    A validator for a core count, which None sets to one per available core.

    Args:
        name (str): The setting, for the messages.

    Returns:
        function: value -> a positive int.
    """

    def validate(value):
        """
        Check that `value` is a positive integer, or turn None into the available core count.

        Logs the core count it picks for None.

        Args:
            value (object): The value given.

        Returns:
            int: The value, or one per available core when it is None.

        Raises:
            PocketMapperError: If the value is not a positive integer, a bool included.
        """
        if value is None:
            # Not os.process_cpu_count() (3.13+) or os.sched_getaffinity (Linux only); None on an
            # exotic platform
            cores = os.cpu_count() or 1
            logger.info(
                f"{name} unset, using one per available core ({cores})", extra={"stage": "Configuring Settings"}
            )
            return cores
        # A job file can hold any JSON value, and a bool passes isinstance(value, int)
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            reject(f"{name} must be a positive integer, got {value!r}.")
        return value

    return validate


def percent_below_100(name):
    """
    A validator for an optional percentage, at least 0 and below 100.

    Args:
        name (str): The setting, for the message.

    Returns:
        function: value -> a float, or None.
    """

    def validate(value):
        """
        Check that `value` is None or a number in [0, 100).

        Args:
            value (object): The value given.

        Returns:
            float or None: The value as a float, or None.

        Raises:
            PocketMapperError: If the value is not a number, is a bool, or is outside [0, 100), NaN included.
        """
        if value is None:
            return None
        # A job file can hold any JSON value, and a bool passes isinstance(value, int). The range check is
        # negated so NaN, which argparse's float and json both accept, fails it too.
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not (0 <= value < 100):
            reject(f"{name} must be a number at least 0 and below 100, got {value!r}.")
        return float(value)

    return validate


def non_negative_int(name):
    """
    A validator for a count that may be 0.

    Args:
        name (str): The setting, for the message.

    Returns:
        function: value -> the value, unchanged.
    """

    def validate(value):
        """
        Check that `value` is an integer of at least 0.

        Args:
            value (object): The value given.

        Returns:
            int: The value, unchanged.

        Raises:
            PocketMapperError: If the value is not an integer, is a bool, or is negative.
        """
        # A job file can hold any JSON value, and a bool passes isinstance(value, int)
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            reject(f"{name} must be a non-negative integer, got {value!r}.")
        return value

    return validate


def integer(name):
    """
    A validator for any integer.

    Args:
        name (str): The setting, for the message.

    Returns:
        function: value -> the value, unchanged.
    """

    def validate(value):
        """
        Check that `value` is an integer.

        Args:
            value (object): The value given.

        Returns:
            int: The value, unchanged.

        Raises:
            PocketMapperError: If the value is not an integer, or is a bool.
        """
        # A job file can hold any JSON value, and a bool passes isinstance(value, int)
        if isinstance(value, bool) or not isinstance(value, int):
            reject(f"{name} must be an integer, got {value!r}.")
        return value

    return validate


# The steps that read parse's entries table
LATER_STEPS = STEPS[1:]

OPTIONS = {
    option.name: option
    for option in (
        # main
        Option(
            "query",
            ("parse", "search"),
            "Query entry, or a file with one entry per line. STRUCT[:CHAIN[:RESIDUES]], e.g. 4Q5J:B_F. "
            "Required unless the job file sets query.",
            flags=("query",),
            metavar="QUERY",
        ),
        Option(
            "target",
            ("parse", "search"),
            "Target entry, a file with one entry per line, or a Foldseek DB name: human_domains, pdb. "
            "Required unless the job file sets target.",
            flags=("target",),
            metavar="TARGET",
        ),
        Option(
            "job_file",
            COMMANDS,
            "JSON file of {\"option\": value}, e.g. a run's job_settings.json or parse's parse_settings.json. "
            "Arguments override it, but it sets every path, so pass a specific path option to move one. A step "
            "other than parse takes every path it uses from here unless given.",
            key=None,
            flags=("-j", "--job_file"),
            metavar="PATH",
            default_help="none",
        ),
        Option(
            "verbosity",
            COMMANDS,
            "Log level: 4=DEBUG, 3=INFO, 2=WARNING, else ERROR.",
            default=DEFAULT_VERBOSITY,
            validate=integer("verbosity"),
            flags=("-v", "--verbosity"),
            type=int,
            metavar="INT",
        ),
        Option(
            "aligner",
            ("align", "search"),
            "Chain aligner: foldseek (needs the binary) or seq (built-in BLOSUM62 sequence aligner).",
            default=DEFAULT_ALIGNER,
            validate=choice("aligner", ALIGNERS),
            flags=("-a", "--aligner"),
        ),
        # Pocket methods are checked by parse.qt_processor.QTProcessor, which owns their values. "auto" is
        # kept as the value: it is resolved for each entry, not once for the run.
        Option(
            "query_pocket_method",
            ("parse", "search"),
            "Query pocket method: auto (infer it from each entry), pisa, passthrough, vdw or whole_chain.",
            default=DEFAULT_POCKET_METHOD,
            flags=("-q", "--query_pocket_method"),
        ),
        Option(
            "target_pocket_method",
            ("parse", "search"),
            "As --query_pocket_method, for targets; also accepts foldseek_db.",
            default=DEFAULT_POCKET_METHOD,
            flags=("-t", "--target_pocket_method"),
        ),
        Option(
            "threads",
            ("fetch_structures", "align", "superpose", "search"),
            "Cap on the cores Foldseek uses.",
            validate=positive_int_or_cores("threads"),
            flags=("-T", "--threads"),
            type=int,
            metavar="INT",
            default_help="one per available core",
        ),
        Option(
            "fetch_missing",
            ("align", "pockets", "superpose"),
            "1 downloads an entry structure missing from the cache; 0 skips the entry, as fetch_structures has "
            "not fetched it. --offline 1 overrides it.",
            default=DEFAULT_FETCH_MISSING,
            validate=flag("fetch_missing"),
            type=int,
            metavar="INT",
        ),
        # in: the files a step reads, each layered over the job key of the file a step writes
        Option(
            "entries",
            LATER_STEPS,
            "The entries table to read.",
            key="entries_path",
            group="in",
            metavar="PATH",
        ),
        Option(
            "pocket_comparison",
            ("superpose",),
            "The pocket comparison table to read.",
            key="pocket_comparison_path",
            group="in",
            metavar="PATH",
        ),
        Option(
            "alignment",
            ("pockets", "compare", "superpose"),
            "The alignment table to read.",
            key="alignment_path",
            group="in",
            metavar="PATH",
        ),
        Option(
            "pockets",
            ("compare", "superpose"),
            "The pockets file to read.",
            key="pockets_path",
            group="in",
            metavar="PATH",
        ),
        # Absolute once resolved. Only the relative paths of local structure files and user Foldseek
        # databases resolve against it; every path setting resolves against the working directory.
        Option(
            "local_root",
            COMMANDS,
            "Directory that relative local structure files and Foldseek databases resolve against. Other paths "
            "resolve against the current directory. A run's settings dump records it, so a later step run from "
            "elsewhere resolves them the same way.",
            group="in",
            metavar="DIR",
            default_help="the job file's, else the current directory",
        ),
        # aligned structure
        Option(
            "align_count",
            ("superpose", "search"),
            "How many top-scoring targets to superpose onto each query; 0 disables.",
            default=DEFAULT_ALIGN_COUNT,
            validate=non_negative_int("align_count"),
            group="aligned structure",
            type=int,
            metavar="INT",
        ),
        # "auto" is resolved to "pocket" or "foldseek" once for the run, from the aligner
        Option(
            "align_struct_method",
            ("superpose", "search"),
            "Which transform superposes a target onto its query: auto, pocket or foldseek.",
            default=DEFAULT_ALIGN_STRUCT_METHOD,
            validate=choice("align_struct_method", ALIGN_STRUCT_METHODS),
            group="aligned structure",
        ),
        # cache
        Option(
            "cache_dir",
            ("parse", "search"),
            "Where structures and PISA responses are cached.",
            default=DEFAULT_CACHE_DIR,
            group="cache",
            metavar="DIR",
        ),
        Option(
            "pdb_dir",
            ("parse", "fetch_structures", "align", "pockets", "superpose", "search"),
            "Cache of fetched PDB structures.",
            group="cache",
            metavar="DIR",
            under=("cache_dir", "pdb_structures"),
        ),
        Option(
            "alphafold_dir",
            ("parse", "fetch_structures", "align", "pockets", "superpose", "search"),
            "Cache of fetched AlphaFold structures.",
            group="cache",
            metavar="DIR",
            under=("cache_dir", "alphafold_structures"),
        ),
        Option(
            "pocket_dir",
            ("parse", "pockets", "search"),
            "Cache of PISA data.",
            group="cache",
            metavar="DIR",
            under=("cache_dir", "pockets"),
        ),
        Option(
            "foldseek_preprocessed_structure_dir",
            ("parse", "align", "search"),
            "Cache of the single-chain structures Foldseek is given.",
            group="cache",
            metavar="DIR",
            under=("cache_dir", "foldseek_preprocessed_structures"),
        ),
        Option(
            "fsdb_dir",
            COMMANDS,
            "Cache of bundled Foldseek databases.",
            group="cache",
            metavar="DIR",
            under=("cache_dir", "fsdb"),
        ),
        # out
        Option(
            "results_dir",
            ("parse", "search"),
            "Where results are written.",
            group="out",
            metavar="DIR",
            default_help="pocketmapper_results_<YYMMDD_HHMMSS>",
            command_help={
                "parse": (
                    "Where the log, the settings and the entries table are written. Without it, none is "
                    "written unless its own path is given.",
                    "none",
                ),
            },
        ),
        Option(
            "aligned_structure_dir",
            ("superpose", "search"),
            "Where superposed structures for the top hits are written.",
            group="out",
            metavar="DIR",
            under=("results_dir", "aligned_structures"),
        ),
        Option(
            "alignment_path",
            ("align", "search"),
            "Where to write the alignment table.",
            group="out",
            metavar="PATH",
            under=("results_dir", "alignment.tsv"),
        ),
        Option(
            "pocket_comparison_path",
            ("compare", "search"),
            "Where to write the pocket comparison table.",
            group="out",
            metavar="PATH",
            under=("results_dir", "pocket_comparison.tsv"),
        ),
        Option(
            "pockets_path",
            ("pockets", "search"),
            "Where to write the pockets.",
            group="out",
            metavar="PATH",
            under=("results_dir", "pockets.json"),
        ),
        Option(
            "pockets_tsv_path",
            ("pockets",),
            "Where to write the table of every pocket: pocket_id, chain, res_auth_ids, method.",
            group="out",
            metavar="PATH",
            default_help="none",
        ),
        Option(
            "entries_path",
            ("parse", "search"),
            "Where to write the entries table: how each entry resolved, and which step dropped it and why.",
            group="out",
            metavar="PATH",
            under=("results_dir", "entries.tsv"),
        ),
        # Never read from a job file. Its default is named for the command by lib.settings.resolve_paths.
        Option(
            "job_settings_path",
            ("parse", "search"),
            "Where to write the run's resolved settings.",
            group="out",
            metavar="PATH",
            default_help="<results_dir>/job_settings.json",
            command_help={
                "parse": (
                    "Where to write the settings, as a job file for the later steps; never read from a job file.",
                    "<results_dir>/parse_settings.json",
                ),
            },
        ),
        Option(
            "log_path",
            COMMANDS,
            "The run log, appended to.",
            group="out",
            metavar="PATH",
            under=("results_dir", "info.log"),
        ),
        # temp
        Option(
            "temp_dir",
            ("fetch_structures", "align", "search"),
            "Scratch space, emptied before use and deleted at the end of a successful run.",
            group="temp",
            metavar="DIR",
            under=("results_dir", "tmp"),
        ),
        Option(
            "delete_tmp",
            ("fetch_structures", "align", "search"),
            "1 deletes --temp_dir at the end of a successful run; 0 keeps it.",
            default=DEFAULT_DELETE_TMP,
            validate=flag("delete_tmp"),
            group="temp",
            type=int,
            metavar="INT",
        ),
        # advanced
        Option(
            "pisa_source",
            ("pockets", "search"),
            "Where PISA interfaces are fetched from: ftp (EBI FTP server) or api (paced PDBe API).",
            default=DEFAULT_PISA_SOURCE,
            validate=choice("pisa_source", PISA_SOURCES),
            group="advanced",
        ),
        Option(
            "pisa_bsa_threshold",
            ("pockets", "search"),
            "Build PISA pockets from buried surface: a residue joins when more than this percent of its "
            "accessible surface is buried by the partner chain.",
            validate=percent_below_100("pisa_bsa_threshold"),
            group="advanced",
            type=float,
            metavar="FLOAT",
            default_help="unset, interface bonds",
        ),
        Option(
            "offline",
            ("fetch_structures", "align", "pockets", "superpose", "search"),
            "1 downloads nothing (no structures, PISA interfaces or Foldseek databases) and skips what is not "
            "cached; 0 downloads what the cache lacks.",
            default=DEFAULT_OFFLINE,
            validate=flag("offline"),
            group="advanced",
            type=int,
            metavar="INT",
        ),
        Option(
            "refresh_structures",
            ("fetch_structures", "search"),
            "1 re-downloads PDB and AlphaFold structures in fetch_structures even when cached; 0 keeps them. "
            "Not with --offline 1.",
            default=DEFAULT_REFRESH_STRUCTURES,
            validate=flag("refresh_structures"),
            group="advanced",
            type=int,
            metavar="INT",
        ),
    )
}


def job_key(option):
    """
    The job key an option sets.

    Args:
        option (Option): The option.

    Returns:
        str or None: Its `key`, else its name; None for an option with no job key (`job_file`).
    """
    return option.name if option.key == "" else option.key


def command_options(command):
    """
    The options a command takes, in registry order.

    Args:
        command (str): One of COMMANDS.

    Returns:
        list of Option: Its options, `job_file` included.
    """
    return [option for option in OPTIONS.values() if command in option.commands]


def command_keys(command):
    """
    The job keys a command's options set, in registry order.

    Args:
        command (str): One of COMMANDS.

    Returns:
        tuple: Each job key once; an input option and the output option it reads share one.
    """
    return tuple(dict.fromkeys(job_key(option) for option in command_options(command) if job_key(option)))


def option_help(option, command):
    """
    An option's `--help` text for one command, its default appended.

    The default text is, in order: the static default; for a command that takes the root a path is
    derived from, the path under it; for a derived path or an input file in any other command, the job
    file's value; else `default_help`. A command in `command_help` has its own help and default text.

    Args:
        option (Option): The option.
        command (str): One of `option.commands`.

    Returns:
        str: The help, ending "(default: ...)" unless the option has no default to state.
    """
    text, default_help = (option.command_help or {}).get(command, (option.help, option.default_help))
    if option.default is not None:
        default = str(option.default)
    elif option.under and command in OPTIONS[option.under[0]].commands:
        default = f"<{option.under[0]}>/{option.under[1]}"
    elif option.under or option.key:
        default = f"the job file's {job_key(option)}"
    else:
        default = default_help
    return f"{text} (default: {default})" if default else text


# What a job file may hold and `lib.settings.layer_settings` seeds, in registry order: the options that set
# their own name. An input option sets the key of the output option it reads.
JOB_KEYS = tuple(name for name, option in OPTIONS.items() if option.key == "")

# Setting -> its static default. The others default to None: resolved at run time, or required.
SETTING_DEFAULTS = {key: OPTIONS[key].default for key in JOB_KEYS if OPTIONS[key].default is not None}

# Path setting -> its default, relative to cache_dir or results_dir
CACHE_PATH_DEFAULTS = {key: OPTIONS[key].under[1] for key in JOB_KEYS if OPTIONS[key].under[:1] == ("cache_dir",)}
RESULTS_PATH_DEFAULTS = {key: OPTIONS[key].under[1] for key in JOB_KEYS if OPTIONS[key].under[:1] == ("results_dir",)}

# The job keys search does not take, and the paths among them, which resolve against the working directory
STEP_ONLY_KEYS = tuple(key for key in JOB_KEYS if key not in command_keys("search"))
STEP_ONLY_PATH_KEYS = tuple(key for key in STEP_ONLY_KEYS if OPTIONS[key].metavar in PATH_METAVARS)
