"""`tempo fold` command stub."""

from __future__ import annotations

import argparse


class FoldCommand:
    """Stub for running AlphaFold to acquire representative structures."""

    name = "fold"
    help = "Run AlphaFold to acquire representative structures."

    @staticmethod
    def configure_parser(parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--sequence", 
            required=True,
            help="Path to the input sequence file."
        )
        parser.add_argument(
            "--template", 
            default=None,
            help="Path to the template structure directory."
        )
        parser.add_argument(
            "--output", 
            default="structures/",
            help="Path for storing generated structures.",
        )

    @staticmethod
    def run(args: argparse.Namespace) -> int:
        from .alphafold import predict, post_process
        predict(args.sequence, args.template)
        post_process("alphafold", args.output)
