from af2rave.amino import AMINO
from af2rave.spib import SPIBProcess
from mahler.utils.os import file_exists

from os import PathLike
from pathlib import Path
import numpy as np
import json

import logging
LOGGER = logging.getLogger("mahler.rave")


def execute(
        file_path: PathLike[str], 
        output_path: PathLike[str],
        ag_chains: list[str],
        ab_chains: list[str],
        max_colvar: int = 30,
        topology: PathLike[str] | None = None,
        suffix: str = "dat"
) -> int:
    
    colvar_files = Path(file_path)
    if colvar_files.is_dir():
        colvar_files = [str(p) for p in colvar_files.glob(f"*.{suffix}") if p.is_file()]
        LOGGER.info(f"Found {len(colvar_files)} colvar files.")
    else:
        raise ValueError(f"{file_path} is not a directory.")
    if len(colvar_files) == 0:
        print(f"No colvar files found matching: {file_path}")
        return -1

    # now check if results are already present
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    all_outputs_exist = all(
        file_exists(_get_output_filename(c, output_path, suffix)) 
        for c in colvar_files
    )
    if all_outputs_exist and file_exists(output_path / "index.npy"):
        LOGGER.info("Looks like AMINO has been finished. Skipping AMINO selection.")

        with open(output_path / "amino_result.json", "r") as f:
            result = json.load(f)        
            print(result)
        with open(output_path / "amino_result.txt", "w") as f:
            f.write("\n".join(_parse_result(topology, result)))
        LOGGER.info(f"AMINO selected {len(result)} CVs.")
    else:
        if not (output_path / "amino_result.json").exists():
            LOGGER.info("Running AMINO to select important CVs.")
            a = AMINO.from_file(colvar_files, n=max_colvar)
            result = a.result
            with open(output_path / "amino_result.json", "w") as f:
                json.dump(result, f, indent=4)
        else:
            LOGGER.info("AMINO result file already exists. Loading existing results.")
            with open(output_path / "amino_result.json", "r") as f:
                result = json.load(f)        
        for c in colvar_files:
            _select(c, output_path, result)
        LOGGER.info(f"AMINO selected {len(result)} CVs.")
        with open(output_path / "amino_result.txt", "w") as f:
            f.write("\n".join(_parse_result(topology, result)))
        index = _parse_index_from_result(result)
        np.save(output_path / "index.npy", index)

    # SPIB
    if (output_path / "spib_model.pkl").exists():
        LOGGER.info("SPIB model already exists. Skipping SPIB training.")
    else:
        spib = SPIBProcess(
            traj = list(output_path.glob(f"*.{suffix}")),
            init="tica:50"
        )
        result = spib.run(time_lag=100, lr_scheduler_gamma=0.9)
        result.to_file(str(output_path / "spib_model.pkl"))

    from .plumed import print_plumed
    print_plumed(
        input_file=topology,
        ag_chains=ag_chains,
        ab_chains=ab_chains,
        spib_model=output_path / "spib_model.pkl",
        index_file=output_path / "index.npy",
        out=output_path / "plumed.dat"
    )

    return 0


def _parse_index_from_result(result) -> np.ndarray:
    """Parse AMINO result to obtain selected residue indices."""
    indices = list()
    for r in result:
        _, i, j = r.split('_')
        indices += [(int(i), int(j))]
    indices = np.array(indices, dtype=int)
    return indices

def _do_amino(
        colvar_files: list[PathLike[str]],
        max_colvar: int,
) -> list[str]:
    """Run AMINO on the provided colvar files to select important CVs."""
    a = AMINO.from_file(colvar_files, n=max_colvar)
    return a.result


def _parse_result(topology: PathLike[str], result) -> list[int]:
    """Parse AMINO result to obtain selected residue indices."""

    import mdtraj as md
    from af2rave.feature import chimera_representation
    
    top = md.load(topology).topology
    commands = []
    for r in result:
        _, i, j = r.split('_')
        i, j = int(i), int(j)
        rep_i = chimera_representation(top, i)
        rep_j = chimera_representation(top, j)
        commands += [f"distance {rep_i} {rep_j}"]

        if not rep_i.endswith("CA"):
            commands += [f"show {rep_i} a"]
        if not rep_j.endswith("CA"):
            commands += [f"show {rep_j} a"]
    commands = sorted(set(commands))
    return commands


def _select(
        input_cv: PathLike[str],
        output_path: PathLike[str],
        result: list[str]
) -> None:
    """Select columns from colvar file based on AMINO result."""

    from cvtoolkit import Colvar

    cv = Colvar.from_file(input_cv)
    cv_chosen = cv.choose(result)
    cv_chosen.write(str(_get_output_filename(input_cv, output_path)))

def _get_output_filename(
        input_cv: PathLike[str],
        output_path: PathLike[str],
        suffix: str = "dat"
) -> PathLike[str]:
    """Generate output filename for selected colvar file."""
    filename = Path(input_cv).name
    output_filename = Path(output_path) / filename
    return output_filename.with_suffix(f".{suffix}")