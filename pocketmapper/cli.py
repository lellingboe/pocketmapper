"""
Command-line front end: argparse over `PocketMapper.search`.

This is the only module that knows about `sys.argv`, terminals or exit codes. `search` is the whole
CLI -- one subcommand, because it is the one public method on `PocketMapper` -- and everything here
exists to turn a command line into its keyword arguments and nothing more. The pipeline itself stays
importable without a terminal; see `pocketmapper.pocketmapper`.

Three parsing details are load-bearing, each for a reason the code alone would not show:

- Query and target are optional positionals; there are no `--query`/`--target` options. A job file
  may supply them instead, and `configure_workflow` requires each from exactly one of the two.
- Defaults are real values, shared with `PocketMapper.search` through `constants`. The job file is
  layered on top of the parsed arguments, so a default here never hides a job-file value. Options
  whose default depends on the run default to None and are resolved downstream.
- No `choices=` anywhere. The same values arrive from the job file, which never passes through
  this parser, so validation lives downstream where both paths reach it.

Author: Lachlan Ellingboe
"""

import argparse
import sys

from pocketmapper.constants import CLI_SEARCH_EPILOG
from pocketmapper.constants import DEFAULT_ALIGN_COUNT
from pocketmapper.constants import DEFAULT_ALIGN_STRUCT_METHOD
from pocketmapper.constants import DEFAULT_ALIGNER
from pocketmapper.constants import DEFAULT_CACHE_DIR
from pocketmapper.constants import DEFAULT_DELETE_TMP
from pocketmapper.constants import DEFAULT_VERBOSITY
from pocketmapper.exceptions import PocketMapperError
from pocketmapper.pocketmapper import PocketMapper

# Accepted spellings for a boolean option value. Fire used to `ast.literal_eval` the token, so
# the True/False spellings must keep working exactly as written.
TRUE_VALUES = ("true", "1", "yes")
FALSE_VALUES = ("false", "0", "no")


def bool_arg(value):
    """
    Parse a boolean option value from the command line.

    Args:
        value (str): The token following the option.

    Returns:
        bool: The parsed value.

    Raises:
        argparse.ArgumentTypeError: If the token is not a recognised boolean spelling.
    """
    lowered = value.strip().lower()
    if lowered in TRUE_VALUES:
        return True
    if lowered in FALSE_VALUES:
        return False
    raise argparse.ArgumentTypeError(f"expected True or False, got {value!r}")


def build_parser():
    """
    Build the top-level parser and its one `search` subcommand.

    Option help here is the whole per-option reference -- it is what `search --help` prints, so it
    must stay in step with the `Settings` dataclass and the README's Options tables. Every `Settings`
    field is reachable from here; the job file is a convenience, never the only route to one.
    Only what argparse cannot generate (the examples) lives in `CLI_SEARCH_EPILOG`, which hangs off
    the `search` subparser alone: the bare `pocketmapper --help` is the subcommand list and nothing
    more.

    Returns:
        argparse.ArgumentParser: The configured parser.
    """
    parser = argparse.ArgumentParser(
        prog="pocketmapper",
        description="PocketMapper - compare the binding surfaces of protein-protein interactions.",
    )
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")

    search = subparsers.add_parser(
        "search",
        help="Run the full search workflow.",
        description="Run the full search workflow.",
        epilog=CLI_SEARCH_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    search.add_argument(
        "query",
        nargs="?",
        default=None,
        metavar="QUERY",
        help="Query entry, or a file with one entry per line. STRUCT[:CHAIN[:RESIDUES]], e.g. 4Q5J:B_F. "
        "Required unless the job file sets query.",
    )
    search.add_argument(
        "target",
        nargs="?",
        default=None,
        metavar="TARGET",
        help="Target entry, a file with one entry per line, or a Foldseek DB name: human_domains, pdb. "
        "Required unless the job file sets target.",
    )
    search.add_argument(
        "--job_file",
        default=None,
        metavar="PATH",
        help='JSON file of {"option": value}, query and target included; it overrides CLI args. ' "(default: none)",
    )
    search.add_argument(
        "--verbosity",
        type=int,
        default=DEFAULT_VERBOSITY,
        metavar="INT",
        help=f"Log level: 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. (default: {DEFAULT_VERBOSITY})",
    )
    search.add_argument(
        "--aligner",
        default=DEFAULT_ALIGNER,
        metavar="STR",
        help="Chain aligner: foldseek (needs the binary) or seq (built-in BLOSUM62 sequence aligner). "
        f"(default: {DEFAULT_ALIGNER})",
    )
    search.add_argument(
        "--query_pocket_method",
        default=None,
        metavar="STR",
        help="Force the query pocket method rather than inferring it from the entry: "
        "pisa, passthrough, vdw or whole_chain. (default: unset)",
    )
    search.add_argument(
        "--target_pocket_method",
        default=None,
        metavar="STR",
        help="As --query_pocket_method, for targets; also accepts foldseek_db. (default: unset)",
    )
    search.add_argument(
        "--threads",
        type=int,
        default=None,
        metavar="INT",
        help="Cap on the cores Foldseek uses. (default: one per available core)",
    )

    # Grouped by lifetime rather than by kind: the twelve path options roughly double the option
    # count, and leaving them in one list would bury --aligner and --query_pocket_method among
    # them. argparse prints groups after the main options, in the order declared.
    aligned_structure_options = search.add_argument_group(
        "aligned structure options",
    )
    aligned_structure_options.add_argument(
        "--align_count",
        type=int,
        default=DEFAULT_ALIGN_COUNT,
        metavar="INT",
        help=f"How many top-scoring targets to superpose onto each query; 0 disables. "
        f"(default: {DEFAULT_ALIGN_COUNT})",
    )
    aligned_structure_options.add_argument(
        "--align_struct_method",
        default=DEFAULT_ALIGN_STRUCT_METHOD,
        metavar="STR",
        help=f"Which transform superposes a target onto its query: auto, pocket or foldseek. "
        f"(default: {DEFAULT_ALIGN_STRUCT_METHOD})",
    )

    cache_paths = search.add_argument_group(
        "cache options",
    )
    cache_paths.add_argument(
        "--cache_dir",
        default=DEFAULT_CACHE_DIR,
        metavar="DIR",
        help=f"Where structures, pockets and PISA responses are cached. (default: {DEFAULT_CACHE_DIR})",
    )
    cache_paths.add_argument(
        "--structure_dir",
        default=None,
        metavar="DIR",
        help="Cache of fetched reference structures. (default: <cache_dir>/ref_structures)",
    )
    cache_paths.add_argument(
        "--pocket_dir",
        default=None,
        metavar="DIR",
        help="Cache of parsed pockets. (default: <cache_dir>/pockets)",
    )
    cache_paths.add_argument(
        "--foldseek_preprocessed_structure_dir",
        default=None,
        metavar="DIR",
        help="Cache of the single-chain structures Foldseek is given. "
        "(default: <cache_dir>/foldseek_preprocessed_structures)",
    )
    cache_paths.add_argument(
        "--fsdb_dir",
        default=None,
        metavar="DIR",
        help="Cache of bundled Foldseek databases. (default: <cache_dir>/fsdb)",
    )

    out_paths = search.add_argument_group(
        "out options",
    )
    out_paths.add_argument(
        "--results_dir",
        default=None,
        metavar="DIR",
        help="Where results are written. (default: pocketmapper_results_<YYMMDD_HHMMSS>)",
    )
    out_paths.add_argument(
        "--aligned_structure_dir",
        default=None,
        metavar="DIR",
        help="Where superposed structures for the top hits are written. (default: <results_dir>/aligned_structures)",
    )
    out_paths.add_argument(
        "--alignment_path",
        default=None,
        metavar="PATH",
        help="Where the alignment table is written. (default: <results_dir>/alignment.tsv)",
    )
    out_paths.add_argument(
        "--pocket_comparison_path",
        default=None,
        metavar="PATH",
        help="Where the pocket comparison table is written. (default: <results_dir>/pocket_comparison.tsv)",
    )
    out_paths.add_argument(
        "--job_settings_path",
        default=None,
        metavar="PATH",
        help="Where this run's resolved settings are dumped. (default: <results_dir>/job_settings.json)",
    )
    out_paths.add_argument(
        "--log_path",
        default=None,
        metavar="PATH",
        help="Where the run log is written. (default: <results_dir>/info.log)",
    )

    temp_paths = search.add_argument_group(
        "temp options",
    )
    temp_paths.add_argument(
        "--temp_dir",
        default=None,
        metavar="DIR",
        help="Per-run scratch, emptied before use and deleted at the end. (default: <results_dir>/tmp)",
    )
    temp_paths.add_argument(
        "--delete_tmp",
        nargs="?",
        const=True,
        type=bool_arg,
        default=DEFAULT_DELETE_TMP,
        metavar="BOOL",
        help=f"Delete --temp_dir at the end of the run; False keeps it. (default: {DEFAULT_DELETE_TMP})",
    )

    return parser


def cli(argv=None):
    """
    Console-script entry point.

    Exits 1 on `PocketMapperError`, which has already been logged with full stage context at the raise
    site. Modules raise rather than calling `exit()` precisely so the package stays embeddable -- keep
    it that way when adding error paths, and keep every exit in this module.

    Args:
        argv (list, optional): Argument list to parse. Defaults to sys.argv[1:].

    Returns:
        None: Exits 1 on a pipeline error; results are written to the run's results_dir.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    # No subcommand is not an error: fire printed its group help and exited 0 here, and CI's compat
    # job runs the bare entry point to prove the install works.
    if args.command is None:
        parser.print_help()
        return

    try:
        PocketMapper().search(
            query=args.query,
            target=args.target,
            job_file=args.job_file,
            cache_dir=args.cache_dir,
            results_dir=args.results_dir,
            verbosity=args.verbosity,
            threads=args.threads,
            aligner=args.aligner,
            align_count=args.align_count,
            align_struct_method=args.align_struct_method,
            query_pocket_method=args.query_pocket_method,
            target_pocket_method=args.target_pocket_method,
            delete_tmp=args.delete_tmp,
            structure_dir=args.structure_dir,
            pocket_dir=args.pocket_dir,
            foldseek_preprocessed_structure_dir=args.foldseek_preprocessed_structure_dir,
            temp_dir=args.temp_dir,
            aligned_structure_dir=args.aligned_structure_dir,
            alignment_path=args.alignment_path,
            pocket_comparison_path=args.pocket_comparison_path,
            job_settings_path=args.job_settings_path,
            log_path=args.log_path,
            fsdb_dir=args.fsdb_dir,
        )
    except PocketMapperError:
        # Already logged with full stage context at the raise site.
        sys.exit(1)
