from __future__ import annotations

from collections.abc import Mapping, Sequence
from os import PathLike
from pathlib import Path

import json
import logging

import mdtraj as md
import numpy as np
import tqdm.rich as tqdm
from cvtoolkit import Colvar

from mahler.proteinmpnn import ScoreMPNN
from mahler.utils.fasta import parse_fasta
from mahler.utils.poisson import bootstrap, fit_poisson
from mahler.utils.topology import check_chains
from mahler.utils.os import list_files

__all__ = ["execute", "score_sequence", "reweighted_time"]

LOGGER = logging.getLogger("mahler.reweight")


Pathish = PathLike[str] | str

def execute(
    traj_directory: Pathish,
    colvar_directory: Pathish,
    sequence_fasta: Pathish,
    output_directory: Pathish,
    score_directory: Pathish | None = None,
    n_bootstrap: int = 1000,
    topology_file: Pathish | None = None,
) -> int:
    """Run the reweighting pipeline for every trajectory/colvar pair."""

    files_mapping = _map_files(traj_directory, colvar_directory)
    if not files_mapping:
        msg = "No trajectory files were found in the provided directory."
        raise FileNotFoundError(msg)

    output_path = Path(output_directory)
    output_path.mkdir(parents=True, exist_ok=True)
    if score_directory is not None:
        score_path = Path(score_directory)
        score_path.mkdir(parents=True, exist_ok=True)

    sequence_dict = parse_fasta(sequence_fasta)
    sequences = list(sequence_dict.values())
    seq_keys = list(sequence_dict.keys())

    all_times: list[float] = []
    all_rw_times: list[np.ndarray] = []

    for traj_file, colvar_file in files_mapping.items():
        LOGGER.info(
            "Processing trajectory file %s with colvar file %s.",
            traj_file,
            colvar_file,
        )
        cache = score_path / f"{traj_file.stem}.npy" if score_directory else None
        time, rw_time = reweighted_time(
            sequences=sequences,
            traj_file=traj_file,
            colvar_file=colvar_file,
            cache=cache,
            stride=500,
            frames_per_batch=10,
            topology_file=topology_file,
        )
        all_times.append(time)
        all_rw_times.append(rw_time)

    rw_array = np.stack(all_rw_times, axis=0)
    with (output_path / "reweighted_times.json").open("w", encoding="utf-8") as fh:
        json.dump({"wt": all_times, **_serialize_sequences(rw_array, seq_keys)}, fh, indent=2)

    wt_array = np.asarray(all_times)
    if n_bootstrap > 0:
        output = {
            "wt": bootstrap(wt_array, n_samples=n_bootstrap).tolist(),
            **{
                key: bootstrap(rw_array[:, i], n_samples=n_bootstrap).tolist()
                for i, key in enumerate(seq_keys)
            },
        }
    else:
        output = {
            "wt": float(fit_poisson(wt_array)[0]),
            **{
                key: float(fit_poisson(rw_array[:, i])[0]) for i, key in enumerate(seq_keys)
            },
        }

    with (output_path / "poisson_fit.json").open("w", encoding="utf-8") as fh:
        json.dump(output, fh, indent=2)

    return 0

    
def _map_files(
        traj_directory: Pathish, 
        colvar_directory: Pathish
) -> dict[Path, Path]:
    """
    Check completeness of directories and map trajectory files to colvar files.

    :returns: dict mapping trajectory file paths to colvar file paths
    :raises: FileNotFoundError. If either directory does not exist or is not a directory.
    :raises: ValueError. If there is not exactly one colvar file for each trajectory file.
    """

    # check completeness of directories
    traj_directory = Path(traj_directory)
    colvar_directory = Path(colvar_directory)
    if not traj_directory.is_dir():
        raise FileNotFoundError(f"Trajectory directory {traj_directory} does not exist or is not a directory.")
    if not colvar_directory.is_dir():
        raise FileNotFoundError(f"Colvar directory {colvar_directory} does not exist or is not a directory.")

    traj_files = list_files(traj_directory)
    print(traj_files)
    files_mapping: dict[Path, Path] = {}
    for tf in traj_files:
        prefix = tf.stem
        cf = list(colvar_directory.glob(f"{prefix}.*"))
        if len(cf) != 1:
            raise ValueError(f"Cannot find exactly one colvar file for trajectory file {tf}.")
        else:
            files_mapping[tf] = cf[0]
    return files_mapping


def _serialize_sequences(array: np.ndarray, keys: Sequence[str]) -> Mapping[str, list[float]]:
    """Return a JSON-friendly mapping of sequence names to per-trajectory values."""

    if array.shape[1] != len(keys):
        msg = "Number of sequence keys does not match array columns."
        raise ValueError(msg)
    return {key: array[:, i].tolist() for i, key in enumerate(keys)}

def score_sequence(
    sequences: Sequence[str],
    traj_file: Pathish,
    cache: Pathish | None = None,
    decoding_order: Sequence[str] | None = None,
    frames_per_batch: int = 1,
    topology_file: Pathish | None = None,
) -> np.ndarray:
    """Score every trajectory frame against the provided sequences."""

    if frames_per_batch <= 0:
        msg = f"frames_per_batch must be positive, received {frames_per_batch}."
        raise ValueError(msg)

    try:
        traj: md.Trajectory = md.load(str(traj_file), top=topology_file)
    except Exception as e:
        msg = f"Failed to load trajectory file {traj_file} with topology {topology_file}: {e}"
        raise RuntimeError(msg) from e

    if cache is not None:
        # check if cache exists and is non-empty
        if Path(cache).is_file() and Path(cache).stat().st_size > 0:
            result = np.load(str(cache), allow_pickle=False)
            LOGGER.info(f"Loaded scores from cache file {cache}.")
            result -= result[:, [0]]
            return result
        else:
            LOGGER.warning(f"Cache file {cache} does not exist or is empty. Proceeding to score trajectory.")

    scorer = ScoreMPNN()
    scores: list[np.ndarray] = []
    chains = check_chains(traj)
    with tqdm.tqdm(total=len(traj), desc="Scoring frames") as progress:
        for start in range(0, len(traj), frames_per_batch):
            end = min(start + frames_per_batch, len(traj))
            chunk = traj.slice(slice(start, end), copy=False)
            chunk_scores = scorer.score(
                sequences,
                chunk,
                chains,
                decoding_order=decoding_order,
            ).T
            scores.append(chunk_scores)
            progress.update(end - start)
    result = np.concatenate(scores, axis=0)
    if cache is not None:
        np.save(str(cache), result, allow_pickle=False)
        LOGGER.info(f"Saved scores to cache file {cache}.")

    # Normalize with respect to the first sequence for each frame.
    result -= result[:, [0]]
    return result


def reweighted_time(
    sequences: Sequence[str],
    traj_file: Pathish,
    colvar_file: Pathish,
    cache: Pathish | None = None,
    stride: int = 500,
    decoding_order: Sequence[str] | None = None,
    frames_per_batch: int = 1,
    topology_file: Pathish | None = None,
) -> tuple[float, np.ndarray]:
    """Return WT and per-sequence reweighted times for a trajectory."""
    if stride <= 0:
        msg = f"Stride must be positive, received {stride}."
        raise ValueError(msg)

    colvar = Colvar.from_file(str(colvar_file))
    dt = float(colvar.time[1] - colvar.time[0])

    beta = 1.0 / (1.987204259e-3 * 310)  # kcal/mol/K
    acc = np.exp(beta * colvar["metad.bias"])

    result = score_sequence(
        sequences,
        traj_file,
        cache=cache,
        decoding_order=decoding_order,
        frames_per_batch=frames_per_batch,
        topology_file=topology_file,
    )

    buffed_result = np.repeat(result, stride, axis=0)

    size = min(buffed_result.shape[0], len(acc))
    if size != len(acc):
        raise RuntimeWarning(
            "The trajectory is too short for the colvar data. Is the stride correct?"
            f" Trajectory length: {result.shape[0]}, stride: {stride}, colvar length: {len(acc)}."
        )

    time = float(np.sum(acc[:size]) * dt)
    reweighted_time = np.sum(
        acc[:size, None] * np.exp(-buffed_result[:size]),
        axis=0,
    ) * dt

    return time, reweighted_time
