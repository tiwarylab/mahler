"""`mahler mdrun` command stub."""

from __future__ import annotations

import argparse
import numpy as np

class MDRunCommand:
    """Stub for running MD production trajectories."""

    name = "mdrun"
    help = "Execute MD production simulations after preparation."

    @staticmethod
    def configure_parser(parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--pdb",
            required=True,
            metavar="PDB",
            help="One or more PDB files describing prepared structures.",
        )
        parser.add_argument(
            "--index",
            metavar="INDEX",
            default=None,
            help="Path to the NumPy index (.npy) file used for selections.",
        )
        parser.add_argument(
            "--colvar",
            metavar="COLVAR",
            default=None,
            help="Output file path for the generated COLVAR data. Leave empty to use <pdb_name>.dat.",
        )
        parser.add_argument(
            "--traj",
            default=None,
            metavar="TRAJ",
            help="Output trajectory file produced by the MD run. No trajectory will be saved if not provided.",
        )
        parser.add_argument(
            "--time", "-t",
            metavar="TIME",
            type=int,
            default=50,
            help="Total simulation time in nanoseconds (default: 50 ns).",
        )
        parser.add_argument(
            "--xtc-freq",
            metavar="XTC_FREQ",
            type=float,
            default=1.0,
            help="Frequency (in ps) to write frames to the XTC trajectory (default: 1.0 ps).",
        )
        parser.add_argument(
            "--checkpnt",
            metavar="CHECKPNT",
            default=None,
            help="Path to save the simulation checkpoint file.",
        )
        parser.add_argument(
            "--final-pdb",
            metavar="FINAL_PDB",
            default=None,
            help="Path to save the final PDB structure after simulation.",
        )

    @staticmethod
    def run(args: argparse.Namespace) -> int:
        from .ubs import execute
        if args.index is None:
            args.index = None
        else:
            args.index = np.load(args.index)
            if args.index.shape[1] != 2:
                raise ValueError("Index file must have shape (N, 2) for pairs of atoms.")
        return execute(
            pdb_file=args.pdb,
            index=args.index,
            xtc_file=args.traj,
            colvar_file=args.colvar,
            time_ns=args.time,
            xtc_freq_ps=args.xtc_freq,
            checkpnt_file=args.checkpnt,
            final_pdb=args.final_pdb,
        )