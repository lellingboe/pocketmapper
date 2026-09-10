"""
Command-line front end: argparse over `PocketMapper.search`.

This is the only module that knows about `sys.argv`, terminals or exit codes. `search` is the whole
CLI -- one subcommand, because it is the one public method on `PocketMapper` -- and everything here
exists to turn a command line into its keyword arguments and nothing more. The pipeline itself stays
importable without a terminal; see `pocketmapper.pocketmapper`.

Three parsing details are load-bearing, each for a reason the code alone would not show:

- Query and target are positional and required; there are no `--query`/`--target` options. The
  `Settings` fields of those names stay for library callers, but argparse always supplies both here,
  so a settings file's `query`/`target` can never win on the CLI path.
- Every option defaults to None, never to a `Settings` default. `_configure_workflow` layers the JSON
  settings file under the CLI arguments by testing `is not None`, so a non-None default here would
  make the settings file unoverridable.
- No `choices=` anywhere. The same values arrive from the settings file, which never passes through
  this parser, so validation lives downstream where both paths reach it.

Author: Lachlan Ellingboe
"""

import argparse
import sys

from pocketmapper.constants import CLI_SEARCH_EPILOG
from pocketmapper.exceptions import PocketMapperError
from pocketmapper.pocketmapper import PocketMapper

# Accepted spellings for a boolean option value. Fire used to `ast.literal_eval` the token, so
# `--foldseek False` is what every local-aligner e2e case passes and what the runner's skip gate
# greps for; the True/False spellings must keep working exactly as written.
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
    field is reachable from here; the settings JSON is a convenience, never the only route to one.
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
        metavar="QUERY",
        help="Query entry, or a file with one entry per line. STRUCT[:CHAIN[:RESIDUES]], e.g. 4Q5J:B_F.",
    )
    search.add_argument(
        "target",
        metavar="TARGET",
        help="Target entry, a file with one entry per line, or a Foldseek DB name: human_domains, pdb.",
    )
    search.add_argument(
        "--settings",
        default=None,
        metavar="PATH",
        help='JSON file of {"option": value}; CLI args override it. (default: none)',
    )
    search.add_argument(
        "--verbosity",
        type=int,
        default=None,
        metavar="INT",
        help="Log level: 4=DEBUG, 3=INFO, 2=WARNING, else ERROR. (default: 3)",
    )
    # nargs="?" with const=True is what makes the bare `--foldseek` mean True while `--foldseek False`
    # still parses, matching what fire did and what the e2e cases pass.
    search.add_argument(
        "--foldseek",
        nargs="?",
        const=True,
        type=bool_arg,
        default=None,
        metavar="BOOL",
        help="Require the Foldseek aligner (True) or forbid it (False); unset auto-detects the binary. "
        "(default: unset)",
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

    # Grouped by lifetime rather than by kind: the twelve path options roughly double the option
    # count, and leaving them in one list would bury --foldseek and --query_pocket_method among
    # them. argparse prints groups after the main options, in the order declared.
    aligned_structure_options = search.add_argument_group(
        "aligned structure options",
    )
    aligned_structure_options.add_argument(
        "--align_count",
        type=int,
        default=None,
        metavar="INT",
        help="How many top-scoring targets to superpose onto each query; 0 disables. (default: 10)",
    )
    aligned_structure_options.add_argument(
        "--align_struct_method",
        default=None,
        metavar="STR",
        help="Which transform superposes a target onto its query: auto, pocket or foldseek. (default: auto)",
    )

    cache_paths = search.add_argument_group(
        "cache options",
    )
    cache_paths.add_argument(
        "--cache_dir",
        default=None,
        metavar="DIR",
        help="Where structures, pockets and PISA responses are cached. (default: pocketmapper_cache)",
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
        "--delete_tmp",
        nargs="?",
        const=True,
        type=bool_arg,
        default=None,
        metavar="BOOL",
        help="Delete the directories below at the end of the run; False keeps them. (default: True)",
    )
    temp_paths.add_argument(
        "--query_dir",
        default=None,
        metavar="DIR",
        help="Per-run query structures, deleted at the end of the run. (default: <results_dir>/query_structures)",
    )
    temp_paths.add_argument(
        "--target_dir",
        default=None,
        metavar="DIR",
        help="Per-run target structures, deleted at the end of the run. (default: <results_dir>/target_structures)",
    )
    # Grouped by lifetime, not by where its default comes from: it is scratch that a Foldseek run
    # deletes, so it belongs here rather than under the caches that survive the run -- and it is the
    # one option in this group whose default hangs off --cache_dir.
    temp_paths.add_argument(
        "--foldseek_tmp_dir",
        default=None,
        metavar="DIR",
        help="Foldseek's scratch directory, deleted after a Foldseek run. (default: <cache_dir>/foldseek_tmp)",
    )

    return parser


def main(argv=None):
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
            settings=args.settings,
            cache_dir=args.cache_dir,
            results_dir=args.results_dir,
            verbosity=args.verbosity,
            foldseek=args.foldseek,
            align_count=args.align_count,
            align_struct_method=args.align_struct_method,
            query_pocket_method=args.query_pocket_method,
            target_pocket_method=args.target_pocket_method,
            delete_tmp=args.delete_tmp,
            structure_dir=args.structure_dir,
            pocket_dir=args.pocket_dir,
            foldseek_tmp_dir=args.foldseek_tmp_dir,
            foldseek_preprocessed_structure_dir=args.foldseek_preprocessed_structure_dir,
            query_dir=args.query_dir,
            target_dir=args.target_dir,
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


if __name__ == "__main__":
    main()
