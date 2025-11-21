"""Entry point for the `mahler` command-line interface."""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence

from .commands import Command, available_commands


def build_parser(commands: Sequence[Command]) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mahler",
        description=(
            "MAHLER: Metadynamics-Anchored Hybrid Learning for Engineering off-Rates"
        ),
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable verbose logging for troubleshooting.",
    )
    subparsers = parser.add_subparsers(dest="command", metavar="<command>")

    for command in commands:
        subparser = subparsers.add_parser(command.name, help=command.help)
        command.configure_parser(subparser)
        subparser.set_defaults(handler=command)

    return parser


def dispatch(args: argparse.Namespace) -> int:
    handler: Command | None = getattr(args, "handler", None)
    if handler is None:
        # argparse will print usage when return code non-zero
        return 1
    return handler.run(args)


def main(argv: Sequence[str] | None = None) -> int:
    commands = available_commands()
    parser = build_parser(commands)
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 0
    # Keep third-party loggers quiet and only elevate MAHLER logs when debugging.
    logging.basicConfig(level=logging.WARNING)
    mahler_logger = logging.getLogger("mahler")
    mahler_logger.setLevel(logging.DEBUG if args.debug else logging.INFO)
    return dispatch(args)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
