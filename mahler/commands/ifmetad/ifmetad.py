try:
    from openmmplumed import PlumedForce
except ImportError:
    raise ImportError(
        "The 'ifmetad' command requires the OpenMM-Plumed plugin. "
        "Please install it following the instructions at "
        "https://github.com/openmm/openmm-plumed"
    )


import openmm.app as app
from pathlib import Path
from af2rave.simulation import UnbiasedSimulation
from mahler.utils.topology import find_protein_subset
from mahler.utils.os import file_exists
import logging
LOGGER = logging.getLogger("mahler.ifmetad")

from os import PathLike


def execute(
        pdb_file: Path,
        xtc_file: Path,
        colvar_file: Path,
        plumed: PathLike[str],
        time_ns: int = 50,
        xtc_freq_ps: float = 10.0,
        final_pdb: PathLike[str] = None,
        temperature: float = 298.15,
) -> int:

    step_ps = 0.002
    steps = int(time_ns * 1000 / step_ps)

    # setting up the trajectory reporter
    report_interval = int(xtc_freq_ps / step_ps)
    if xtc_file.is_dir():
        i, xtc_file = _find_next_name(xtc_file, pdb_file, "xtc")
    xtc_rep = app.xtcreporter.XTCReporter(
        str(xtc_file),
        reportInterval=report_interval,
        atomSubset=find_protein_subset(pdb_file),
    )
    LOGGER.info(f"XTC trajectory will be saved to {xtc_file}.")

    # setting up PLUMED colvar redirection
    if colvar_file.is_dir():
        colvar_file = colvar_file / (f"{pdb_file.stem}_{i:02d}.dat")

    if not file_exists(plumed, non_empty=True):
        LOGGER.error(f"PLUMED script {plumed} does not exist or is empty.")
        return 1
    with open(plumed, "r") as f:
        plumed_script = f.read()
        plumed_script = plumed_script.replace("COLVAR", str(colvar_file))
        plumed_script = plumed_script.replace("HILLS", str(colvar_file.with_suffix(".hills")))
    plumed_force = PlumedForce(plumed_script)
    LOGGER.info(f"Using PLUMED script from {plumed}.")

    forcefield = app.ForceField("amber14/protein.ff14SB.xml", "amber14/tip3p.xml")

    ubs = UnbiasedSimulation(
        str(pdb_file),
        xtc_reporter=xtc_rep,
        forcefield=forcefield,
        custom_forces=[plumed_force],
        temp=temperature
    )

    ubs.run(steps)
    if final_pdb:
        ubs.save_pdb(final_pdb)
    return 0

def _find_next_name(
        base_dir: Path,
        pdb_file: Path,
        postfix: str = "xtc",
):
    for i in range(1, 100):
        candidate = base_dir / f"{pdb_file.stem}_{i:02d}.{postfix}"
        if not file_exists(candidate, non_empty=False):
            # calling dibs on this name to avoid race conditions over other replicas
            candidate.touch()
            return i, candidate
    raise FileExistsError(
        f"Could not find a free name for {pdb_file.stem} in {base_dir} after 99 attempts."
    )