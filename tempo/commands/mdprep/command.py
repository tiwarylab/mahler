"""`tempo mdprep` command stub."""

from __future__ import annotations

import argparse

from ..base import CommandResult
from .feature import execute


class MDPrepCommand:
    """Stub for preparing MD clusters before downstream workflows."""

    name = "mdprep"
    help = "Prepare MD structures by clustering trajectories."

    @staticmethod
    def configure_parser(parser: argparse.ArgumentParser) -> None:
        parser.description = (
            "Cluster MD-ready structures and create feature selections from a reference trajectory."
        )
        parser.add_argument("--prefix", required=True, help="Path to structures or trajectory.")
        parser.add_argument(
            "--output_dir",
            required=True,
            help="Directory to write feature selection artifacts and clustering output.",
        )
        parser.add_argument(
            "--ref",
            default=None,
            help="Reference PDB used for contact detection (defaults to first frame).",
        )
        parser.add_argument(
            "--steric_clash_cutoff",
            type=float,
            default=1.0,
            help="Minimum non-bonded distance for steric clash filtering (Angstrom).",
        )
        parser.add_argument(
            "--antigen_chains",
            default="A",
            metavar="CHAINS",
            help="Comma-separated chain IDs treated as antigen (letters).",
        )
        parser.add_argument(
            "--antibody_chains",
            default="B,C",
            metavar="CHAINS",
            help="Comma-separated chain IDs treated as antibody (letters).",
        )
        parser.add_argument(
            "--exclude_cb",
            action="store_true",
            help="Exclude CB atoms from feature selections (default: include).",
        )
        parser.add_argument(
            "--n_features",
            type=int,
            default=200,
            help="Number of highest-variance features to retain.",
        )
        parser.add_argument(
            "--n_clusters",
            type=int,
            default=15,
            help="Number of cluster centroids to compute from PCA space.",
        )
        parser.add_argument(
            "--random_seed",
            type=int,
            default=42,
            help="Random seed for clustering reproducibility.",
        )
        parser.add_argument(
            "--qval_cutoff",
            type=float,
            default=0.85,
            help="Minimum Q-value fraction to keep structures before clustering.",
        )
        parser.add_argument(
            "--json",
            default=None,
            help="Path to JSON file containing arguments (overrides CLI flags).",
        )

    @staticmethod
    def run(args: argparse.Namespace) -> int:
        if args.json:
            import json
            from pathlib import Path

            json_config = json.loads(Path(args.json).read_text())
            execute(**json_config)
            return 0

        def parse_chain_ids(value: str) -> tuple[str, ...]:
            return tuple(part.strip() for part in value.split(",") if part.strip())

        execute(
            prefix=args.prefix,
            output_dir=args.output_dir,
            ref=args.ref,
            steric_clash_cutoff=args.steric_clash_cutoff,
            ag_chains=parse_chain_ids(args.antigen_chains),
            ab_chains=parse_chain_ids(args.antibody_chains),
            include_cb=not args.exclude_cb,
            n_features=args.n_features,
            n_clusters=args.n_clusters,
            random_seed=args.random_seed,
            qval_cutoff=args.qval_cutoff,
        )
        return 0
