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
            "--traj",
            "-t",
            required=True,
            help="Directory containing trajectory files.",
        )
        parser.add_argument(
            "--colvar",
            "-c",
            required=True,
            help="Directory containing corresponding colvar files.",
        )
        parser.add_argument(
            "--fasta",
            "-f",
            required=True,
            help="FASTA file with sequences to score.",
        )
        parser.add_argument(
            "--output",
            "-o",
            required=True,
            help="Directory where analysis outputs will be written.",
        )
        parser.add_argument(
            "--bootstrap",
            "-n",
            type=int,
            default=0,
            help="Number of bootstrap samples to draw (set to 0 disable).",
        )
        parser.add_argument(
            "--cache",
            "-k",
            required=False,
            help="Directory where score cache files will be read/written. If none is provided, defaults to the output directory.",
        )

    @staticmethod
    def run(args: argparse.Namespace) -> int:
        from .reweight import execute
        if args.cache is not None:
            cache_directory = args.cache
        else:
            cache_directory = args.output
        return execute(
            traj_directory=args.traj,
            colvar_directory=args.colvar,
            sequence_fasta=args.fasta,
            score_directory=cache_directory,
            output_directory=args.output,
            n_bootstrap=args.bootstrap
        )
