"""`mahler reweight` command stub."""

from __future__ import annotations

import argparse


class ReweightCommand:
    """Stub for reweighting trajectories."""

    name = "reweight"
    help = "Reweight trajectories."

    @staticmethod
    def configure_parser(parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--config",
            help="Path to reweighting configuration file.",
        )

    @staticmethod
    def run(args: argparse.Namespace) -> int:
        raise NotImplementedError("`mahler reweight` is not implemented yet.")
