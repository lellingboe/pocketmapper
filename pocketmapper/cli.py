"""
Command-line front end: argparse over `PocketMapper.search` and the step entry functions in
`pocketmapper.<step>.<step>`.

This is the only module that knows about `sys.argv`, terminals or exit codes. Everything here exists
to turn a command line into one function's keyword arguments and nothing more: each subcommand's
options are exactly its function's parameters, and `cli()` passes the parsed arguments straight
through. The pipeline itself stays importable without a terminal.

`lib.options` holds every option once, with its commands and `--help` group; `build_parser` reads it, and
`DISPATCH` maps each subcommand to its function. Three parsing details are load-bearing,
each for a reason the code alone would not show:

- `search`'s and `parse`'s query and target are optional positionals; there are no
  `--query`/`--target` options. A job file may supply them instead, and `lib.settings.layer_command`
  and `lib.settings.require_setting` require each from exactly one of the two.
- Every option defaults to None, for unset: an argument given beats the job file, and the defaults
  (`lib.options.SETTING_DEFAULTS`, or resolved at run time) apply only after it, so a real default here
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
    PACKAGE_LOGGER,
)
from pocketmapper.lib.exceptions import PocketMapperError
from pocketmapper.lib.logs import format_handler
from pocketmapper.lib.options import COMMANDS, GROUPS, command_options, option_help
from pocketmapper.parse.parse import parse
from pocketmapper.pockets.pockets import pockets
from pocketmapper.search import PocketMapper
from pocketmapper.superpose.superpose import superpose

# Command -> its description, for `--help`
COMMAND_DESCRIPTIONS = {
    "parse": "Check how each query and target entry parses, before anything is fetched, and print a table of "
    "them. No network.",
    "fetch_structures": "Download the structures and Foldseek databases the entries table names.",
    "align": "Align the query chains against the target chains.",
    "pockets": "Build the pocket of every query and target entry.",
    "compare": "Compare the pockets of every aligned query/target pair.",
    "superpose": "Superpose the top targets of each query onto it.",
    "search": "Run the full search workflow: every step above, in order.",
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


def add_option(group, option, command):
    """
    Add one registry option to an argparse group, defaulting to None.

    Args:
        group (argparse._ActionsContainer): The parser or argument group to add it to.
        option (lib.options.Option): The option.
        command (str): The command whose help text it gets.
    """
    kwargs = dict(default=None, metavar=option.metavar, help=option_help(option, command))
    if option.type is not None:
        kwargs["type"] = option.type
    flags = option.flags or (f"--{option.name}",)
    # A positional is the only option whose flag has no dash; it is optional, since a job file may supply it
    if not flags[0].startswith("-"):
        kwargs["nargs"] = "?"
    group.add_argument(*flags, **kwargs)


def build_parser():
    """
    Build the top-level parser and one subparser per command in `lib.options.COMMANDS`.

    Returns:
        argparse.ArgumentParser: The configured parser.
    """
    parser = argparse.ArgumentParser(
        prog="pocketmapper",
        description="PocketMapper - compare the binding surfaces of protein-protein interactions.",
    )
    parser.add_argument("-V", "--version", action="version", version=f"pocketmapper {__version__}")
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")

    for command in COMMANDS:
        description = COMMAND_DESCRIPTIONS[command]
        subparser = subparsers.add_parser(
            command,
            help=description,
            description=description,
            epilog=CLI_COMMAND_EPILOGS[command],
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )
        options = command_options(command)
        for title in GROUPS:
            members = [option for option in options if option.group == title]
            if not members:
                continue
            group = subparser if title == "main" else subparser.add_argument_group(f"{title} options")
            for option in members:
                add_option(group, option, command)

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

    # Before dispatch: layering can log a critical on a bad job file before any log is open, and
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
