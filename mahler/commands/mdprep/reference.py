# This matches the reference pdb to out AF2 PDB indexing

import mdtraj as md
import numpy as np
from af2rave.simulation.utils import TopologyMap
from mahler.utils.os import file_exists 

from os import PathLike
import logging
logger = logging.getLogger("mahler.mdprep.index") 

def generate_reference(
        old_pdb: PathLike[str],
        af2_pdb: PathLike[str],
        output_pdb: PathLike[str],
) -> np.ndarray[int]:
    """
    Generate a reference PDB file that matches the indexing of an AF2 PDB file.
    Returns the sulfide bond pairs in AF2 topology.
    
    :param old_pdb: Path to the reference PDB file.
    :param af2_pdb: Path to the AF2 PDB file.
    :param output_pdb: Path to the output PDB file.
    :return: np.ndarray of sulfide bond pairs in AF2 topology.
    :rtype: np.ndarray[int]
    :raises FileNotFoundError: If input files do not exist or are empty.
    :raises ValueError: If atom index mapping fails.
    """
    
    if not file_exists(old_pdb):
        raise FileNotFoundError(f"The reference PDB file {old_pdb} does not exist or is empty.")
    if not file_exists(af2_pdb):
        raise FileNotFoundError(f"The AF2 PDB file {af2_pdb} does not exist or is empty.")
    if file_exists(output_pdb):
        logger.warning(f"The output PDB file {output_pdb} already exists and will be overwritten.")

    ref_traj = md.load(old_pdb)
    af2_traj = md.load(af2_pdb)
    n_atoms_af2 = af2_traj.n_atoms

    topmap = TopologyMap(ref_traj.topology, af2_traj.topology)
    try:
        idx = topmap.map_atom_index(np.arange(n_atoms_af2))
    except ValueError as e:
        raise ValueError(f"Error mapping atom indices between reference and AF2 PDBs: {e}")

    # for every atom in the af2 topology, find the corresponding atom in the ref topology
    xyz = np.zeros_like(af2_traj.xyz)
    for i_af2, i_ref in enumerate(idx):
        xyz[:, i_af2, :] = ref_traj.xyz[0, i_ref, :]

    new_traj = md.Trajectory(xyz, af2_traj.topology)
    new_traj.save_pdb(output_pdb)

    # now get the sulfide bonds
    sulfide_pairs = []
    for b in ref_traj.topology.bonds:
        if b[0].name == "SG" and b[1].name == "SG":
            sulfide_pairs.append((b[0].residue.index, b[1].residue.index))

    return np.array(sulfide_pairs, dtype=int)