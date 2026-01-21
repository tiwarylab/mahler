"""`mahler ifmetad` command stub."""

from __future__ import annotations

import argparse
from pathlib import Path

class IfMetadCommand:
    """Stub for running infrequent metadynamics simulations."""

    name = "ifmetad"
    help = "Run infrequent metadynamics."

    @staticmethod
    def configure_parser(parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--pdb",
            "-i",
            required=True,
            metavar="PDB",
            help="One or more PDB files describing prepared structures.",
        )
        parser.add_argument(
            "--out",
            "-o",
            metavar="DIR",
            help="Output directory to save traj/colvar produced by the MD run.",
        )
        parser.add_argument(
            "--plumed",
            "-p",
            metavar="PLUMED_SCRIPT",
            help="Path to the PLUMED script for infrequent metadynamics.",
        )
        parser.add_argument(
            "--time", 
            metavar="TIME",
            type=int,
            default=50,
            help="Total simulation time in nanoseconds (default: 50 ns).",
        )
        parser.add_argument(
            "--xtc-freq",
            metavar="XTC_FREQ",
            type=float,
            default=10.0,
            help="Frequency (in ps) to write frames to the XTC trajectory (default: 10.0 ps).",
        )
        parser.add_argument(
            "--final-pdb",
            metavar="FINAL_PDB",
            default=None,
            help="Path to save the final PDB structure after simulation.",
        )
        parser.add_argument(
            "--temperature", "-t",
            metavar="TEMP",
            type=float,
            default=298.15,
            help="Simulation temperature in Kelvin (default: 298.15 K).",
        )

    @staticmethod
    def run(args: argparse.Namespace) -> int:
        from .ifmetad import execute
        return execute(
            pdb_file=Path(args.pdb),
            xtc_file=Path(args.out),
            colvar_file=Path(args.out),
            plumed=Path(args.plumed),
            time_ns=args.time,
            xtc_freq_ps=args.xtc_freq,
            final_pdb=Path(args.final_pdb) if args.final_pdb is not None else None,
            temperature=args.temperature,
        )