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
            metavar="<Colvar_dir>",
            help="Directory containing COLVAR files to select from.",
        )
        parser.add_argument(
            "--output",
            "-o",
            required=True,
            metavar="<output_dir>",
            help="Directory where filtered COLVAR files and SPIB model will be written.",
        )
        parser.add_argument(
            "--antigen_chains",
            "-ag",
            default="A",
            metavar="CHAINS",
            help="Comma-separated chain IDs treated as antigen (letters).",
        )
        parser.add_argument(
            "--antibody_chains",
            "-ab",
            default="B,C",
            metavar="CHAINS",
            help="Comma-separated chain IDs treated as antibody (letters).",
        )
        parser.add_argument(
            "--structure",
            "-s",
            required=True,
            metavar="<PDB_file>",
            help="Starting structure PDB file.",
        )
        parser.add_argument(
            "--max-colvar",
            "-n",
            type=int,
            default=30,
            metavar="<max_n>",
            help="Maximum number of collective variables to keep after AMINO selection.",
        )
        parser.add_argument(
            "--suffix",
            default="dat",
            metavar="<suffix>",
            help="File suffix for both input COLVAR files and filtered outputs (default: dat).",
        )

    @staticmethod
    def run(args: argparse.Namespace) -> int:

        def parse_chain_ids(value: str) -> tuple[str, ...]:
            return tuple(part.strip() for part in value.split(",") if part.strip())
        
        from .rave import execute
        return execute(
            file_path=args.colvar,
            output_path=args.output,
            ag_chains=parse_chain_ids(args.antigen_chains),
            ab_chains=parse_chain_ids(args.antibody_chains),
            max_colvar=args.max_colvar,
            topology=args.structure,
            suffix=args.suffix,
        )

