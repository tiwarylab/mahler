"""`mahler ifmetad` command stub."""

from __future__ import annotations

import argparse


class IfMetadCommand:
    """Stub for running infrequent metadynamics simulations."""

    name = "ifmetad"
    help = "Run infrequent metadynamics."

    @staticmethod
    def configure_parser(parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--config",
            help="Path to infrequent metadynamics configuration file.",
        )

    @staticmethod
    def run(args: argparse.Namespace) -> int:
        raise NotImplementedError("`mahler ifmetad` is not implemented yet.")
