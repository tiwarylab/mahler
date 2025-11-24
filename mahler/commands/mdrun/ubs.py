from os import PathLike
from collections.abc import Sequence
from pathlib import Path

import mdtraj as md
import numpy as np
import openmm.app as app
from af2rave.simulation import UnbiasedSimulation

import logging
LOGGER = logging.getLogger("mahler.mdrun")


def execute(
    pdb_file: PathLike[str],
    index: Sequence[int],
    xtc_file: PathLike[str] = None,
    colvar_file: PathLike[str] = None,
    time_ns: int = 50,
    xtc_freq_ps: float = 1.0,
    checkpnt_file: PathLike[str] | None = None,
    final_pdb: PathLike[str] | None = None,
) -> int:
    """Run an unbiased molecular dynamics simulation and persist checkpoints."""

    step_ps = 0.002
    steps = int(time_ns * 1000 / step_ps)

    # setting up the trajectory reporter
    if xtc_file is not None:
        report_interval = int(xtc_freq_ps / step_ps)
        if report_interval <= 500:
            LOGGER.warning(f"XTC report interval is {report_interval} steps ({xtc_freq_ps} ps), which may be too frequent.")
        xtc_rep = app.xtcreporter.XTCReporter(
            xtc_file,
            reportInterval=report_interval,
            atomSubset=_find_protein_subset(pdb_file),
        )
    else:
        LOGGER.warning("No trajectory file provided; no trajectory will be saved.")
        xtc_rep = None

    if index is None:
        raise NotImplementedError("mahler.mdrun requires an index of atom pairs.")
    
    # default values for colvar reporting
    if colvar_file is None:
        colvar_file = Path(pdb_file).with_suffix(".dat")

    ubs = UnbiasedSimulation(
        pdb_file,
        list_of_index=np.asarray(index, dtype=int),
        xtc_reporter=xtc_rep,
        cv_file=colvar_file,
        cv_freq=500,   # I would want to fix this number to avoid further problems.
    )

    ubs.run(steps)
    if final_pdb:
        ubs.save_pdb(final_pdb)
    if checkpnt_file:
        ubs.save_checkpoint(checkpnt_file)
    
    return 0


def _find_protein_subset(pdb_file: PathLike[str]) -> np.ndarray:
    """Return atom indices corresponding to protein residues."""
    traj = md.load_pdb(pdb_file)
    protein_atoms = traj.topology.select("protein")
    return protein_atoms
