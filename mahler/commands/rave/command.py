"""`mahler rave` command stub."""

from __future__ import annotations

import argparse


class RaveCommand:
    """Stub for running the RAVE protocol to obtain a latent space."""

    name = "rave"
    help = "Run the RAVE protocol and acquire a latent space."

    @staticmethod
    def configure_parser(parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--colvar",
            "-i",
            required=True,
            metavar="COLVAR_DIR",
            help="Directory containing COLVAR files to select from.",
        )
        parser.add_argument(
            "--output",
            "-o",
            required=True,
            metavar="OUTPUT_DIR",
            help="Directory where filtered COLVAR files and SPIB model will be written.",
        )
        parser.add_argument(
            "--max-colvar",
            type=int,
            default=30,
            metavar="MAX",
            help="Maximum number of collective variables to keep after AMINO selection.",
        )
        parser.add_argument(
            "--topology",
            default=None,
            metavar="TOPOLOGY",
            help="Optional topology file used to generate Chimera selection commands.",
        )
        parser.add_argument(
            "--suffix",
            default="dat",
            metavar="SUFFIX",
            help="File suffix for both input COLVAR files and filtered outputs (default: dat).",
        )

    @staticmethod
    def run(args: argparse.Namespace) -> int:
        
        from .rave import execute
        return execute(
            file_path=args.colvar,
            output_path=args.output,
            max_colvar=args.max_colvar,
            topology=args.topology,
            suffix=args.suffix,
        )

