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
        out_path: Path,
        replica_id: int,
        plumed: PathLike[str],
        time_ns: int,
        xtc_freq_ps: float,
        temperature: float,
) -> int:

    # setting up the trajectory reporter
    step_ps = 0.002
    steps = int(time_ns * 1000 / step_ps)
    report_interval = int(xtc_freq_ps / step_ps)

    # setting up output paths
    if not out_path.is_dir() and out_path.exists():
        LOGGER.error(f"Output path {out_path} exists and is not a directory.")
        return 1
    out_path.mkdir(parents=True, exist_ok=True)
    
    prefix = f"{pdb_file.stem}_r.{replica_id}"
    xtc_file = out_path / f"{prefix}.xtc"
    colvar_file = out_path / f"{prefix}.dat"
    hills_file = out_path / f"{prefix}.hills"
    final_pdb = out_path / f"{prefix}.pdb"
    plumed_dat = out_path / f"{prefix}.plumed"

    xtc_rep = app.xtcreporter.XTCReporter(
        str(xtc_file),
        reportInterval=report_interval,
        atomSubset=find_protein_subset(pdb_file),
    )
    LOGGER.info(f"XTC trajectory will be saved to {xtc_file}.")
    LOGGER.info(f"PLUMED COLVAR will be saved to {colvar_file}.")

    if not file_exists(plumed, non_empty=True):
        LOGGER.error(f"PLUMED script {plumed} does not exist or is empty.")
        return 1
    with open(plumed, "r") as f:
        plumed_script = f.read()
        plumed_script = plumed_script.replace("COLVAR", str(colvar_file))
        plumed_script = plumed_script.replace("HILLS", str(hills_file))
    plumed_force = PlumedForce(plumed_script)
    with open(plumed_dat, "w") as f:
        f.write(plumed_script)
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
    ubs.save_pdb(str(final_pdb))
    return 0
