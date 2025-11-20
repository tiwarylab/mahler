"""`tempo rave` command stub."""

from __future__ import annotations

import argparse


class RaveCommand:
    """Stub for running the RAVE protocol to obtain a latent space."""

    name = "rave"
    help = "Run the RAVE protocol and acquire a latent space."

    @staticmethod
    def configure_parser(parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--config",
            help="Path to RAVE configuration file.",
        )

    @staticmethod
    def run(args: argparse.Namespace) -> int:
        raise NotImplementedError("`tempo rave` is not implemented yet.")
