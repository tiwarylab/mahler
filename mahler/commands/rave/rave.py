from af2rave.amino import AMINO
from af2rave.spib import SPIBProcess

from os import PathLike
from pathlib import Path
import json

import logging
LOGGER = logging.getLogger("mahler.rave")


def execute(
        file_path: PathLike[str], 
        output_path: PathLike[str],
        max_colvar: int = 30,
        topology: PathLike[str] | None = None,
        suffix: str = "dat"
) -> int:
    
    colvar_files = Path(file_path)
    if colvar_files.is_dir():
        colvar_files = [str(p) for p in colvar_files.glob(f"*.{suffix}") if p.is_file()]
    else:
        raise ValueError(f"{file_path} is not a directory.")
    if len(colvar_files) == 0:
        print(f"No colvar files found matching: {file_path}")
        return -1

    # now check if results are already present
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    if all((_get_output_filename(c, output_path, suffix)).exists() for c in colvar_files):
        LOGGER.info("All selected colvar files already exist. Skipping AMINO selection.")
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

    if topology is not None:
        _parse_result(topology, result)

    # SPIB
    if (output_path / "spib_model.pkl").exists():
        LOGGER.info("SPIB model already exists. Skipping SPIB training.")
    else:
        spib = SPIBProcess(
            output_path.glob(f"*.{suffix}"),
            init="tica:50"
        )
        result = spib.run(time_lag=100, lr_scheduler_gamma=0.9)
        result.to_file(str(output_path / "spib_model.pkl"))

    return 0


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
    from af2rave.feature import chimera_representation, atom_name
    
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