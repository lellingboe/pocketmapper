"""
Command-line front end: argparse over `PocketMapper.search`.

This is the only module that knows about `sys.argv`, terminals or exit codes. `search` is the whole
CLI -- one subcommand, because it is the one public method on `PocketMapper` -- and everything here
exists to turn a command line into its keyword arguments and nothing more. The pipeline itself stays
importable without a terminal; see `pocketmapper.pocketmapper`.

Three parsing details are load-bearing, each for a reason the code alone would not show:

- Query and target are accepted both positionally and as `--query`/`--target`. Nearly every e2e case
  uses the positional form, which predates this module.
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
_TRUE_VALUES = ("true", "1", "yes")
_FALSE_VALUES = ("false", "0", "no")


def _bool_arg(value):
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
    if lowered in _TRUE_VALUES:
        return True
    if lowered in _FALSE_VALUES:
        return False
    raise argparse.ArgumentTypeError(f"expected True or False, got {value!r}")


def _merge_positional(parser, named, positional, name):
    """
    Collapse the positional and `--name` spellings of one argument into a single value.

    Both spellings are supported, so they need separate argparse dests -- a `nargs="?"` positional
    sharing a dest with an option overwrites that option with its own default whenever it is absent,
    silently discarding `--query`.

    Args:
        parser (argparse.ArgumentParser): Parser to raise the usage error through.
        named (str or None): Value from `--name`.
        positional (str or None): Value from the positional slot.
        name (str): Option name, for the error message.

    Returns:
        str or None: Whichever was supplied, or None if neither was.
    """
    if named is not None and positional is not None:
        parser.error(f"{name} given both positionally ({positional!r}) and as --{name} ({named!r})")
    return named if named is not None else positional


def _build_parser():
    """
    Build the top-level parser and its one `search` subcommand.

    Option help here is the whole per-option reference -- it is what `search --help` prints, so it
    must stay in step with the `Settings` dataclass and the README's Options table. Only what argparse
    cannot generate (the settings-file-only paths, and the examples) lives in `CLI_SEARCH_EPILOG`,
    which hangs off the `search` subparser alone: the bare `pocketmapper --help` is the subcommand
    list and nothing more.

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
        description="PocketMapper - compare the binding surfaces of protein-protein interactions.",
        epilog=CLI_SEARCH_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    search.add_argument(
        "query_pos",
        nargs="?",
        default=None,
        metavar="QUERY",
        help="Positional spelling of --query.",
    )
    search.add_argument(
        "target_pos",
        nargs="?",
        default=None,
        metavar="TARGET",
        help="Positional spelling of --target.",
    )
    search.add_argument(
        "--query",
        default=None,
        metavar="STR",
        help="Query entry, or a file with one entry per line. STRUCT[:CHAIN[:RESIDUES]], e.g. 4Q5J:B_F. (required)",
    )
    search.add_argument(
        "--target",
        default=None,
        metavar="STR",
        help="Target entry, a file with one entry per line, or a Foldseek DB name: human_domains, pdb. (required)",
    )
    search.add_argument(
        "--settings",
        default=None,
        metavar="PATH",
        help='JSON file of {"option": value}; CLI args override it. (default: none)',
    )
    search.add_argument(
        "--cache_dir",
        default=None,
        metavar="DIR",
        help="Where structures, pockets and PISA responses are cached. (default: pocketmapper_cache)",
    )
    search.add_argument(
        "--results_dir",
        default=None,
        metavar="DIR",
        help="Where results are written. (default: pocketmapper_results_<YYMMDD_HHMMSS>)",
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
        type=_bool_arg,
        default=None,
        metavar="BOOL",
        help="Require the Foldseek aligner (True) or forbid it (False); unset auto-detects the binary. "
        "(default: unset)",
    )
    search.add_argument(
        "--align_count",
        type=int,
        default=None,
        metavar="INT",
        help="How many top-scoring targets to superpose onto each query; 0 disables. (default: 10)",
    )
    search.add_argument(
        "--align_struct_method",
        default=None,
        metavar="STR",
        help="Which transform superposes a target onto its query: auto, pocket or foldseek. (default: auto)",
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
    # So main() can report a usage error against the subparser the argument actually belongs to;
    # the top-level parser's usage line names only COMMAND, which is no help to someone who
    # mistyped a search option.
    search.set_defaults(subparser=search)
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
    parser = _build_parser()
    args = parser.parse_args(argv)

    # No subcommand is not an error: fire printed its group help and exited 0 here, and CI's compat
    # job runs the bare entry point to prove the install works.
    if args.command is None:
        parser.print_help()
        return

    query = _merge_positional(args.subparser, args.query, args.query_pos, "query")
    target = _merge_positional(args.subparser, args.target, args.target_pos, "target")

    try:
        PocketMapper().search(
            query=query,
            target=target,
            settings=args.settings,
            cache_dir=args.cache_dir,
            results_dir=args.results_dir,
            verbosity=args.verbosity,
            foldseek=args.foldseek,
            align_count=args.align_count,
            align_struct_method=args.align_struct_method,
            query_pocket_method=args.query_pocket_method,
            target_pocket_method=args.target_pocket_method,
        )
    except PocketMapperError:
        # Already logged with full stage context at the raise site.
        sys.exit(1)


if __name__ == "__main__":
    main()
