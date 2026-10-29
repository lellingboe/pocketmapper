"""
Command-line front end: argparse over `PocketMapper.search` and the step entry functions in
`pocketmapper.<step>.<step>`.

This is the only module that knows about `sys.argv`, terminals or exit codes. Everything here exists
to turn a command line into one function's keyword arguments and nothing more: each subcommand's
options are exactly its function's parameters, and `cli()` passes the parsed arguments straight
through. The pipeline itself stays importable without a terminal.

`OPTIONS` holds every option once; `COMMANDS` lists which options each subcommand takes, grouped
for `--help`; `DISPATCH` maps each subcommand to its function. Three parsing details are load-bearing,
each for a reason the code alone would not show:

- `search`'s and `parse`'s query and target are optional positionals; there are no
  `--query`/`--target` options. A job file may supply them instead, and `lib.settings.layer_settings`
  and `lib.settings.require_setting` require each from exactly one of the two.
- Every option defaults to None, for unset: an argument given beats the job file, and the defaults
  (`lib.settings.SETTING_DEFAULTS`, or resolved at run time) apply only after it, so a real default here
  would hide every job-file value. The help still states each default.
- No `choices=` anywhere. The same values arrive from the job file and from library calls, which
  never pass through this parser, so validation lives downstream where every path reaches it.

Author: Lachlan Ellingboe
"""

import argparse
import logging
import sys

from pocketmapper import __version__
from pocketmapper.align.align import align
from pocketmapper.compare.compare import compare
from pocketmapper.fetch.fetch import fetch_structures
from pocketmapper.lib.constants import (
    CLI_COMMAND_EPILOGS,
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
    PACKAGE_LOGGER,
)
from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.logs import format_handler
from pocketmapper.parse.parse import parse
from pocketmapper.pockets.pockets import pockets
from pocketmapper.search import PocketMapper
from pocketmapper.superpose.superpose import superpose


def input_file(flag, name, key):
    """
    The spec of an option naming a file a step reads.

    Args:
        flag (str): The option, e.g. "--alignment".
        name (str): What the file holds, for the help.
        key (str): The job key it defaults to, e.g. "alignment_path", for the help.

    Returns:
        tuple: (flags, add_argument keyword arguments).
    """
    return [flag], dict(default=None, metavar="PATH", help=f"The {name} to read. (default: the job file's {key})")


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
            help="JSON file of {\"option\": value}, e.g. a run's job_settings.json or parse's "
            "parse_settings.json. Arguments override it, but it sets every path, so pass a specific path "
            "option to move one. A step other than parse takes every path it uses from here unless given. "
            "(default: none)",
        ),
    ),
    "local_root": (
        ["--local_root"],
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
            help=f"As --query_pocket_method, for targets; also accepts foldseek_db. (default: {DEFAULT_POCKET_METHOD})",
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
            help=f"Where structures and PISA responses are cached. (default: {DEFAULT_CACHE_DIR})",
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
        dict(default=None, metavar="DIR", help="Cache of PISA data. (default: <cache_dir>/pockets)"),
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
            help="Where the log, the settings and the entries table are written. Without it, none is "
            "written unless its own path is given. (default: none)",
        ),
    ),
    "entries": input_file("--entries", "entries table", "entries_path"),
    "alignment": input_file("--alignment", "alignment table", "alignment_path"),
    "pockets": input_file("--pockets", "pockets file", "pockets_path"),
    "pocket_comparison": input_file("--pocket_comparison", "pocket comparison table", "pocket_comparison_path"),
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
    "job_settings_path": output_file(
        "--job_settings_path", "run's resolved settings", "<results_dir>/job_settings.json"
    ),
    "job_settings_path_parse": output_file(
        "--job_settings_path",
        "settings, as a job file for the later steps; never read from a job file",
        "<results_dir>/parse_settings.json",
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
            help="Scratch space, emptied before use and deleted at the end of a successful run. "
            "(default: <results_dir>/tmp)",
        ),
    ),
    "delete_tmp": (
        ["--delete_tmp"],
        dict(
            type=int,
            default=None,
            metavar="INT",
            help=f"1 deletes --temp_dir at the end of a successful run; 0 keeps it. (default: {DEFAULT_DELETE_TMP})",
        ),
    ),
    "entries_path": output_file(
        "--entries_path",
        "entries table: how each entry resolved, and which step dropped it and why",
        "<results_dir>/entries.tsv",
    ),
    "pockets_tsv_path": output_file(
        "--pockets_tsv_path", "table of every pocket: pocket_id, chain, res_auth_ids, method", "none"
    ),
    "fetch_missing": (
        ["--fetch_missing"],
        dict(
            type=int,
            default=None,
            metavar="INT",
            help="1 downloads an entry structure missing from the cache; 0 skips the entry, as fetch_structures "
            f"has not fetched it. --offline 1 overrides it. (default: {DEFAULT_FETCH_MISSING})",
        ),
    ),
    "offline": (
        ["--offline"],
        dict(
            type=int,
            default=None,
            metavar="INT",
            help="1 downloads nothing (no structures, PISA interfaces or Foldseek databases) and skips "
            f"what is not cached; 0 downloads what the cache lacks. (default: {DEFAULT_OFFLINE})",
        ),
    ),
    "refresh_structures": (
        ["--refresh_structures"],
        dict(
            type=int,
            default=None,
            metavar="INT",
            help="1 re-downloads PDB and AlphaFold structures in fetch_structures even when cached; 0 keeps "
            f"them. Not with --offline 1. (default: {DEFAULT_REFRESH_STRUCTURES})",
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
    "pisa_bsa_threshold": (
        ["--pisa_bsa_threshold"],
        dict(
            type=float,
            default=None,
            metavar="FLOAT",
            help="Build PISA pockets from buried surface: a residue joins when more than this percent of its "
            "accessible surface is buried by the partner chain. (default: unset, interface bonds)",
        ),
    ),
}

# Command -> (description, [(--help group title, or None for the main options, [option names])]).
# Groups follow the options' lifetime rather than their kind: in, aligned structure, cache, out,
# temp, advanced. argparse prints groups after the main options, in the order declared.
COMMANDS = {
    "parse": (
        "Check how each query and target entry parses, before anything is fetched, and print a table of "
        "them. No network.",
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
            ("in options", ["local_root"]),
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
                ["results_dir_optional", "entries_path", "job_settings_path_parse", "log_path"],
            ),
        ],
    ),
    "fetch_structures": (
        "Download the structures and Foldseek databases the entries table names.",
        [
            (None, ["job_file", "verbosity", "threads"]),
            ("in options", ["entries"]),
            ("out options", ["log_path"]),
            ("temp options", ["temp_dir", "delete_tmp"]),
            ("advanced options", ["offline", "refresh_structures"]),
        ],
    ),
    "align": (
        "Align the query chains against the target chains.",
        [
            (None, ["job_file", "verbosity", "aligner", "threads", "fetch_missing"]),
            ("in options", ["entries"]),
            ("cache options", ["foldseek_preprocessed_structure_dir"]),
            ("out options", ["alignment_path", "log_path"]),
            ("temp options", ["temp_dir", "delete_tmp"]),
            ("advanced options", ["offline"]),
        ],
    ),
    "pockets": (
        "Build the pocket of every query and target entry.",
        [
            (None, ["job_file", "verbosity", "fetch_missing"]),
            ("in options", ["entries", "alignment"]),
            ("cache options", ["pdb_dir", "alphafold_dir", "pocket_dir", "fsdb_dir"]),
            ("out options", ["pockets_path", "pockets_tsv_path", "log_path"]),
            ("advanced options", ["pisa_source", "pisa_bsa_threshold", "offline"]),
        ],
    ),
    "compare": (
        "Compare the pockets of every aligned query/target pair.",
        [
            (None, ["job_file", "verbosity"]),
            ("in options", ["entries", "alignment", "pockets"]),
            ("out options", ["pocket_comparison_path", "log_path"]),
        ],
    ),
    "superpose": (
        "Superpose the top targets of each query onto it.",
        [
            (None, ["job_file", "verbosity", "threads", "fetch_missing"]),
            ("in options", ["entries", "pocket_comparison", "alignment", "pockets"]),
            ("aligned structure options", ["align_count", "align_struct_method"]),
            ("out options", ["aligned_structure_dir", "log_path"]),
            ("advanced options", ["offline"]),
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
            ("in options", ["local_root"]),
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
                    "results_dir",
                    "aligned_structure_dir",
                    "alignment_path",
                    "pocket_comparison_path",
                    "pockets_path",
                    "entries_path",
                    "job_settings_path",
                    "log_path",
                ],
            ),
            ("temp options", ["temp_dir", "delete_tmp"]),
            ("advanced options", ["pisa_source", "pisa_bsa_threshold", "offline", "refresh_structures"]),
        ],
    ),
}

# Command -> the function its options are the parameters of
DISPATCH = {
    "parse": parse,
    "fetch_structures": fetch_structures,
    "align": align,
    "pockets": pockets,
    "compare": compare,
    "superpose": superpose,
    "search": PocketMapper().search,
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
    parser.add_argument("-V", "--version", action="version", version=f"pocketmapper {__version__}")
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

    Adds a stdout handler to the `pocketmapper` logger for the length of the call. Prints the table
    `parse` returns once it has run. Calls `sys.exit(1)` on a `PocketMapperError`.

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

    kwargs = {dest: value for dest, value in vars(args).items() if dest != "command"}

    # Before dispatch: layer_settings can log a critical on a bad job file before any log is open, and
    # the inherited root WARNING level lets it through
    package_logger = logging.getLogger(PACKAGE_LOGGER)
    handler = format_handler(logging.StreamHandler(sys.stdout))
    package_logger.addHandler(handler)
    try:
        result = DISPATCH[args.command](**kwargs)
        # After the run, so it does not interleave with the log, which goes to stdout too
        if args.command == "parse":
            result.to_csv(sys.stdout, sep="\t", index=False)
    except PocketMapperError:
        # Already logged at the raise site
        sys.exit(1)
    finally:
        package_logger.removeHandler(handler)
