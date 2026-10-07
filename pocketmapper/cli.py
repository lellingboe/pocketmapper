"""
Command-line front end: argparse over `PocketMapper.search` and the step entry functions in
`pocketmapper.steps`.

This is the only module that knows about `sys.argv`, terminals or exit codes. Everything here exists
to turn a command line into one function's keyword arguments and nothing more: each subcommand's
options are exactly its function's parameters, and `cli()` passes the parsed arguments straight
through. The pipeline itself stays importable without a terminal.

`OPTIONS` holds every option once; `COMMANDS` lists which options each subcommand takes, grouped
for `--help`. Three parsing details are load-bearing, each for a reason the code alone would not show:

- `search`'s and `parse`'s query and target are optional positionals; there are no
  `--query`/`--target` options. A job file may supply them instead, and `settings.layer_settings`
  and `settings.require_setting` require each from exactly one of the two.
- Every option defaults to None, for unset: an argument given beats the job file, and the defaults
  (`settings.SETTING_DEFAULTS`, or resolved at run time) apply only after it, so a real default here
  would hide every job-file value. The help still states each default.
- No `choices=` anywhere. The same values arrive from the job file and from library calls, which
  never pass through this parser, so validation lives downstream where every path reaches it.

Author: Lachlan Ellingboe
"""

import argparse
import logging
import sys

from pocketmapper.constants import CLI_COMMAND_EPILOGS
from pocketmapper.constants import DEFAULT_ALIGN_COUNT
from pocketmapper.constants import DEFAULT_ALIGN_STRUCT_METHOD
from pocketmapper.constants import DEFAULT_ALIGNER
from pocketmapper.constants import DEFAULT_CACHE_DIR
from pocketmapper.constants import DEFAULT_DELETE_TMP
from pocketmapper.constants import DEFAULT_FETCH_MISSING
from pocketmapper.constants import DEFAULT_PISA_SOURCE
from pocketmapper.constants import DEFAULT_POCKET_METHOD
from pocketmapper.constants import DEFAULT_VERBOSITY
from pocketmapper.constants import PACKAGE_LOGGER
from pocketmapper.exceptions import PocketMapperError
from pocketmapper.lib import format_handler
from pocketmapper.pocketmapper import PocketMapper
from pocketmapper.steps.align import align
from pocketmapper.steps.compare import compare
from pocketmapper.steps.fetch_structures import fetch_structures
from pocketmapper.steps.parse import parse
from pocketmapper.steps.pockets import pockets
from pocketmapper.steps.superpose import superpose


def input_file(flag, name, default):
    """
    The spec of an option naming a file a command reads.

    Args:
        flag (str): The option, e.g. "--alignment".
        name (str): What the file holds, for the help.
        default (str): The default, for the help.

    Returns:
        tuple: (flags, add_argument keyword arguments).
    """
    return [flag], dict(default=None, metavar="PATH", help=f"The {name} to read. (default: {default})")


def output_file(flag, name, default):
    """
    The spec of an option naming a file a command writes.

    Args:
        flag (str): The option, e.g. "--alignment_path".
        name (str): What the file holds, for the help.
        default (str): The default, for the help.

    Returns:
        tuple: (flags, add_argument keyword arguments).
    """
    return [flag], dict(default=None, metavar="PATH", help=f"Where to write the {name}. (default: {default})")


# Option name -> (flags, add_argument keyword arguments). The name is the dest, except for a variant
# of an option whose wording differs between commands, named after it with a suffix.
OPTIONS = {
    "query_optional": (
        ["query"],
        dict(
            nargs="?",
            default=None,
            metavar="QUERY",
            help="Query entry, or a file with one entry per line. STRUCT[:CHAIN[:RESIDUES]], e.g. 4Q5J:B_F. "
            "Required unless the job file sets query.",
        ),
    ),
    "target_optional": (
        ["target"],
        dict(
            nargs="?",
            default=None,
            metavar="TARGET",
            help="Target entry, a file with one entry per line, or a Foldseek DB name: human_domains, pdb. "
            "Required unless the job file sets target.",
        ),
    ),
    "job_file": (
        ["-j", "--job_file"],
        dict(
            default=None,
            metavar="PATH",
            help="JSON file of {\"option\": value}, e.g. a run's job_settings.json or a step's "
            "<command>_settings.json. Arguments override it, but it sets every path, so pass a specific path "
            "option to move one. (default: none)",
        ),
    ),
    "work_dir": (
        ["--work_dir"],
        dict(
            default=None,
            metavar="DIR",
            help="Directory that entries, entries files and relative paths resolve against. A run's "
            "settings dump records it, so a later step run from elsewhere resolves them the same way. "
            "(default: the current directory)",
        ),
    ),
    "verbosity": (
        ["-v", "--verbosity"],
        dict(
            type=int,
            default=None,
            metavar="INT",
            help=f"Log level: 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. (default: {DEFAULT_VERBOSITY})",
        ),
    ),
    "aligner": (
        ["-a", "--aligner"],
        dict(
            default=None,
            metavar="STR",
            help="Chain aligner: foldseek (needs the binary) or seq (built-in BLOSUM62 sequence aligner). "
            f"(default: {DEFAULT_ALIGNER})",
        ),
    ),
    "query_pocket_method": (
        ["-q", "--query_pocket_method"],
        dict(
            default=None,
            metavar="STR",
            help="Query pocket method: auto (infer it from each entry), pisa, passthrough, vdw or whole_chain. "
            f"(default: {DEFAULT_POCKET_METHOD})",
        ),
    ),
    "target_pocket_method": (
        ["-t", "--target_pocket_method"],
        dict(
            default=None,
            metavar="STR",
            help="As --query_pocket_method, for targets; also accepts foldseek_db. "
            f"(default: {DEFAULT_POCKET_METHOD})",
        ),
    ),
    "threads": (
        ["-T", "--threads"],
        dict(
            type=int,
            default=None,
            metavar="INT",
            help="Cap on the cores Foldseek uses. (default: one per available core)",
        ),
    ),
    "align_count": (
        ["--align_count"],
        dict(
            type=int,
            default=None,
            metavar="INT",
            help="How many top-scoring targets to superpose onto each query; 0 disables. "
            f"(default: {DEFAULT_ALIGN_COUNT})",
        ),
    ),
    "align_struct_method": (
        ["--align_struct_method"],
        dict(
            default=None,
            metavar="STR",
            help="Which transform superposes a target onto its query: auto, pocket or foldseek. "
            f"(default: {DEFAULT_ALIGN_STRUCT_METHOD})",
        ),
    ),
    "cache_dir": (
        ["--cache_dir"],
        dict(
            default=None,
            metavar="DIR",
            help=f"Where structures, pockets and PISA responses are cached. (default: {DEFAULT_CACHE_DIR})",
        ),
    ),
    "pdb_dir": (
        ["--pdb_dir"],
        dict(
            default=None, metavar="DIR", help="Cache of fetched PDB structures. (default: <cache_dir>/pdb_structures)"
        ),
    ),
    "alphafold_dir": (
        ["--alphafold_dir"],
        dict(
            default=None,
            metavar="DIR",
            help="Cache of fetched AlphaFold structures. (default: <cache_dir>/alphafold_structures)",
        ),
    ),
    "pocket_dir": (
        ["--pocket_dir"],
        dict(default=None, metavar="DIR", help="Cache of parsed pockets. (default: <cache_dir>/pockets)"),
    ),
    "foldseek_preprocessed_structure_dir": (
        ["--foldseek_preprocessed_structure_dir"],
        dict(
            default=None,
            metavar="DIR",
            help="Cache of the single-chain structures Foldseek is given. "
            "(default: <cache_dir>/foldseek_preprocessed_structures)",
        ),
    ),
    "fsdb_dir": (
        ["--fsdb_dir"],
        dict(default=None, metavar="DIR", help="Cache of bundled Foldseek databases. (default: <cache_dir>/fsdb)"),
    ),
    "results_dir": (
        ["--results_dir"],
        dict(
            default=None,
            metavar="DIR",
            help="Where results are written. (default: pocketmapper_results_<YYMMDD_HHMMSS>)",
        ),
    ),
    "results_dir_optional": (
        ["--results_dir"],
        dict(
            default=None,
            metavar="DIR",
            help="Where the log, the settings and the failed entries are written. Without it, none is "
            "written unless its own path is given. (default: none)",
        ),
    ),
    "results_dir_required": (
        ["--results_dir"],
        dict(
            default=None,
            metavar="DIR",
            help="The results directory; inputs and outputs default to files in it. Required here or in the "
            "job file.",
        ),
    ),
    "alignment": input_file("--alignment", "alignment table", "<results_dir>/alignment.tsv"),
    "pockets": input_file("--pockets", "pockets file", "<results_dir>/pockets.json"),
    "pocket_comparison": input_file(
        "--pocket_comparison", "pocket comparison table", "<results_dir>/pocket_comparison.tsv"
    ),
    "aligned_structure_dir": (
        ["--aligned_structure_dir"],
        dict(
            default=None,
            metavar="DIR",
            help="Where superposed structures for the top hits are written. (default: <results_dir>/aligned_structures)",
        ),
    ),
    "alignment_path": output_file("--alignment_path", "alignment table", "<results_dir>/alignment.tsv"),
    "pockets_path": output_file("--pockets_path", "pockets", "<results_dir>/pockets.json"),
    "pocket_comparison_path": output_file(
        "--pocket_comparison_path", "pocket comparison table", "<results_dir>/pocket_comparison.tsv"
    ),
    "failed_entries_path": (
        ["--failed_entries_path"],
        dict(
            default=None,
            metavar="PATH",
            help="Where the entries dropped along the way are listed, with the reason. "
            "(default: <results_dir>/failed_entries.json)",
        ),
    ),
    "job_settings_path": output_file(
        "--job_settings_path", "run's resolved settings", "<results_dir>/job_settings.json"
    ),
    "job_settings_path_step": output_file(
        "--job_settings_path",
        "settings, as a job file for later steps; never read from a job file",
        "<results_dir>/<command>_settings.json",
    ),
    "log_path": (
        ["--log_path"],
        dict(default=None, metavar="PATH", help="The run log, appended to. (default: <results_dir>/info.log)"),
    ),
    "temp_dir": (
        ["--temp_dir"],
        dict(
            default=None,
            metavar="DIR",
            help="Scratch space, emptied before use and deleted at the end. (default: <results_dir>/tmp)",
        ),
    ),
    "delete_tmp": (
        ["--delete_tmp"],
        dict(
            type=int,
            default=None,
            metavar="INT",
            help=f"1 deletes --temp_dir at the end; 0 keeps it. (default: {DEFAULT_DELETE_TMP})",
        ),
    ),
    "fetch_missing": (
        ["--fetch_missing"],
        dict(
            type=int,
            default=None,
            metavar="INT",
            help="1 downloads an entry structure missing from the cache; 0 skips the entry, as fetch_structures "
            f"has not fetched it. (default: {DEFAULT_FETCH_MISSING})",
        ),
    ),
    "pisa_source": (
        ["--pisa_source"],
        dict(
            default=None,
            metavar="STR",
            help="Where PISA interfaces are fetched from: ftp (EBI FTP server) or api (paced PDBe API). "
            f"(default: {DEFAULT_PISA_SOURCE})",
        ),
    ),
}

# Command -> (description, [(--help group title, or None for the main options, [option names])]).
# Groups follow the options' lifetime rather than their kind: in, aligned structure, cache, out,
# temp, advanced. argparse prints groups after the main options, in the order declared.
COMMANDS = {
    "parse": (
        "Check how each query and target entry parses, before anything is fetched. No network.",
        [
            (
                None,
                [
                    "query_optional",
                    "target_optional",
                    "job_file",
                    "verbosity",
                    "query_pocket_method",
                    "target_pocket_method",
                ],
            ),
            ("in options", ["work_dir"]),
            (
                "cache options",
                [
                    "cache_dir",
                    "pdb_dir",
                    "alphafold_dir",
                    "pocket_dir",
                    "foldseek_preprocessed_structure_dir",
                    "fsdb_dir",
                ],
            ),
            ("out options", ["results_dir_optional", "failed_entries_path", "job_settings_path_step", "log_path"]),
        ],
    ),
    "fetch_structures": (
        "Download the structures and Foldseek database the query and target entries need.",
        [
            (None, ["job_file", "verbosity", "threads"]),
            ("in options", ["work_dir"]),
            (
                "cache options",
                [
                    "cache_dir",
                    "pdb_dir",
                    "alphafold_dir",
                    "pocket_dir",
                    "foldseek_preprocessed_structure_dir",
                    "fsdb_dir",
                ],
            ),
            (
                "out options",
                ["results_dir_optional", "failed_entries_path", "job_settings_path_step", "log_path"],
            ),
            ("temp options", ["temp_dir", "delete_tmp"]),
        ],
    ),
    "align": (
        "Align the query chains against the target chains.",
        [
            (None, ["job_file", "verbosity", "aligner", "threads", "fetch_missing"]),
            ("in options", ["work_dir"]),
            (
                "cache options",
                [
                    "cache_dir",
                    "pdb_dir",
                    "alphafold_dir",
                    "pocket_dir",
                    "foldseek_preprocessed_structure_dir",
                    "fsdb_dir",
                ],
            ),
            (
                "out options",
                [
                    "results_dir_required",
                    "alignment_path",
                    "failed_entries_path",
                    "job_settings_path_step",
                    "log_path",
                ],
            ),
            ("temp options", ["temp_dir", "delete_tmp"]),
        ],
    ),
    "pockets": (
        "Build the pocket of every query and target entry.",
        [
            (None, ["job_file", "verbosity", "fetch_missing"]),
            ("in options", ["work_dir", "alignment"]),
            (
                "cache options",
                [
                    "cache_dir",
                    "pdb_dir",
                    "alphafold_dir",
                    "pocket_dir",
                    "foldseek_preprocessed_structure_dir",
                    "fsdb_dir",
                ],
            ),
            (
                "out options",
                [
                    "results_dir_required",
                    "pockets_path",
                    "failed_entries_path",
                    "job_settings_path_step",
                    "log_path",
                ],
            ),
            ("advanced options", ["pisa_source"]),
        ],
    ),
    "compare": (
        "Compare the pockets of every aligned query/target pair.",
        [
            (None, ["job_file", "verbosity"]),
            ("in options", ["work_dir", "alignment", "pockets"]),
            (
                "cache options",
                [
                    "cache_dir",
                    "pdb_dir",
                    "alphafold_dir",
                    "pocket_dir",
                    "foldseek_preprocessed_structure_dir",
                    "fsdb_dir",
                ],
            ),
            (
                "out options",
                [
                    "results_dir_required",
                    "pocket_comparison_path",
                    "failed_entries_path",
                    "job_settings_path_step",
                    "log_path",
                ],
            ),
        ],
    ),
    "superpose": (
        "Superpose the top targets of each query onto it.",
        [
            (None, ["job_file", "verbosity", "threads", "fetch_missing"]),
            ("in options", ["work_dir", "pocket_comparison", "alignment", "pockets"]),
            ("aligned structure options", ["align_count", "align_struct_method"]),
            (
                "cache options",
                [
                    "cache_dir",
                    "pdb_dir",
                    "alphafold_dir",
                    "pocket_dir",
                    "foldseek_preprocessed_structure_dir",
                    "fsdb_dir",
                ],
            ),
            (
                "out options",
                [
                    "results_dir_required",
                    "aligned_structure_dir",
                    "failed_entries_path",
                    "job_settings_path_step",
                    "log_path",
                ],
            ),
        ],
    ),
    "search": (
        "Run the full search workflow: every step above, in order.",
        [
            (
                None,
                [
                    "query_optional",
                    "target_optional",
                    "job_file",
                    "verbosity",
                    "aligner",
                    "query_pocket_method",
                    "target_pocket_method",
                    "threads",
                ],
            ),
            ("aligned structure options", ["align_count", "align_struct_method"]),
            ("in options", ["work_dir"]),
            (
                "cache options",
                [
                    "cache_dir",
                    "pdb_dir",
                    "alphafold_dir",
                    "pocket_dir",
                    "foldseek_preprocessed_structure_dir",
                    "fsdb_dir",
                ],
            ),
            (
                "out options",
                [
                    "results_dir",
                    "aligned_structure_dir",
                    "alignment_path",
                    "pocket_comparison_path",
                    "pockets_path",
                    "failed_entries_path",
                    "job_settings_path",
                    "log_path",
                ],
            ),
            ("temp options", ["temp_dir", "delete_tmp"]),
            ("advanced options", ["pisa_source"]),
        ],
    ),
}


def build_parser():
    """
    Build the top-level parser and one subparser per command in `COMMANDS`.

    Returns:
        argparse.ArgumentParser: The configured parser.
    """
    parser = argparse.ArgumentParser(
        prog="pocketmapper",
        description="PocketMapper - compare the binding surfaces of protein-protein interactions.",
    )
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")

    for command, (description, groups) in COMMANDS.items():
        subparser = subparsers.add_parser(
            command,
            help=description,
            description=description,
            epilog=CLI_COMMAND_EPILOGS[command],
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )
        for title, names in groups:
            group = subparser if title is None else subparser.add_argument_group(title)
            for name in names:
                flags, kwargs = OPTIONS[name]
                group.add_argument(*flags, **kwargs)

    return parser


def cli(argv=None):
    """
    Console-script entry point.

    Adds a stdout handler to the `pocketmapper` logger for the length of the call. Calls `sys.exit(1)`
    on a `PocketMapperError`.

    Args:
        argv (list, optional): Argument list to parse. Defaults to sys.argv[1:].

    Returns:
        None: Exits 1 on a pipeline error; results are written to the run's results_dir.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    # No subcommand is not an error: print the help and exit 0
    if args.command is None:
        parser.print_help()
        return

    # Command -> the function its options are the parameters of
    dispatch = {
        "parse": parse,
        "fetch_structures": fetch_structures,
        "align": align,
        "pockets": pockets,
        "compare": compare,
        "superpose": superpose,
        "search": PocketMapper().search,
    }
    kwargs = {dest: value for dest, value in vars(args).items() if dest != "command"}

    package_logger = logging.getLogger(PACKAGE_LOGGER)
    handler = format_handler(logging.StreamHandler(sys.stdout))
    package_logger.addHandler(handler)
    try:
        dispatch[args.command](**kwargs)
    except PocketMapperError:
        # Already logged at the raise site
        sys.exit(1)
    finally:
        package_logger.removeHandler(handler)
