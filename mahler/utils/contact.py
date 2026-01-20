from typing import TextIO
import os

import mdtraj as md
import numpy as np

import logging
LOGGER = logging.getLogger(__name__)

class NativeContact:

    '''
    Native contact class to compute the native contacts of a protein.
    The native contacts are defined as the pairs of atoms that are within a
    certain distance from each other in the reference structure.
    
    The native contacts are computed using the following formula:
    frac_of_contacts = 1 / (1 + exp(beta * (dist - lamda * dist_ref)))
    where dist is the distance between the atoms in the trajectory,
    dist_ref is the distance between the atoms in the reference structure,
    lamda is the distance cutoff, and beta is the steepness of the function.

    Arguments
    ---------
    reference : md.Trajectory
        The reference trajectory to compute the native contacts from.
    selection : list[str]
        The selection of atoms to compute the native contacts from.
        The selection should be in the format of mdtraj.
    cutoff : float
        The distance cutoff to compute the native contacts.
        The default value is 0.45 nm.
    '''

    def __init__(self, reference: md.Trajectory, selection: list[str], cutoff: float = 0.45):
        
        self._native = reference
        self._top: md.Topology = reference.topology
        self._selection = selection

        LOGGER.debug(f"Computing native contacts for {selection[0]} and {selection[1]}")
        pairs = self._top.select_pairs(selection[0], selection[1])
        dist = md.compute_distances(self._native, pairs, periodic=False)[0]

        self._native_contacts = pairs[dist < cutoff]
        self._dist = dist[dist < cutoff]
        self._inv_S = 1.0 / self._native_contacts.shape[0]

    @property
    def contacts(self) -> np.ndarray:
        """Return the native contacts as an (n_contacts, 2) array."""
        return self._native_contacts
    
    @property
    def contacts_atoms(self) -> np.ndarray:
        """Return the sorted, unique atom indices (0-based) involved in native contacts."""
        # Flatten once and use np.unique for deterministic ordering
        return np.unique(self.contacts)
    
    @property
    def contacts_residues(self) -> dict[str, list[int]]:
        """Return residue indices (0-based) per chain that participate in native contacts."""
        residues: dict[str, list[int]] = {}
        for atom_index in self.contacts_atoms:
            residue = self._top.atom(int(atom_index)).residue
            chain = residue.chain.chain_id
            if chain is None:
                raise RuntimeError(f"Residue {residue} has no chain designation.")
            residues.setdefault(chain, []).append(residue.index)

        # Remove duplicates per chain and keep indices sorted for determinism
        for chain, residue_list in residues.items():
            residues[chain] = sorted(set(residue_list))
        return residues
    
    @property
    def contacts_ca(self) -> dict[str, list[int]]:
        """Return CA atom indices (0-based) per chain represented in native contacts."""

        # Use sets for O(1) membership when checking which residues contribute to contacts
        residues_by_chain = {
            chain: set(residue_indices) for chain, residue_indices in self.contacts_residues.items()
        }
        ca_per_chain: dict[str, list[int]] = {}

        for atom_index in self._top.select("name CA"):
            atom = self._top.atom(int(atom_index))
            residue = atom.residue
            chain = residue.chain.chain_id

            if chain in residues_by_chain and residue.index in residues_by_chain[chain]:
                ca_per_chain.setdefault(chain, []).append(int(atom_index))

        for chain, indices in ca_per_chain.items():
            ca_per_chain[chain] = sorted(indices)

        return ca_per_chain
    

    @property
    def ref_ifd(self) -> np.ndarray:
        """Return the inter-group distance (IFD) in the reference structure."""
        return self.compute_ifd(self._native)

    def compute_ifd(self, traj: md.Trajectory) -> np.ndarray:

        # all the CA atoms involving in contacts
        ca_idx = []
        for chain in self.contacts_ca.values():
            ca_idx += chain
        
        # intersect with selection
        ca_idx_sel0 = sorted(list(set(traj.top.select(self._selection[0])).intersection(ca_idx)))
        ca_idx_sel1 = sorted(list(set(traj.top.select(self._selection[1])).intersection(ca_idx)))

        # compute distances
        ag_com = np.array(traj.xyz[:, ca_idx_sel0, :].mean(axis=1))
        ab_com = np.array(traj.xyz[:, ca_idx_sel1, :].mean(axis=1))
        new_dist = np.linalg.norm(ag_com - ab_com, axis=1) * 10  # Convert to Angstroms
        return new_dist

    def compute(self, traj: md.Trajectory, lamda: float = 1.8, beta: float = 50) -> np.ndarray:
        '''
        Compute the fraction of native contacts of the trajectory.
    
        Arguments
        ---------
        traj : md.Trajectory
            The trajectory to compute the native contacts from.
        lamda : float
            The distance cutoff to compute the native contacts.
            The default value is 1.8.
        beta : float
            The steepness of the function to compute the native contacts.
            The default value is 50.0/ns. When using Angstrom divide by 10.
        '''

        dist = md.compute_distances(traj, self._native_contacts, periodic=False)
        
        with np.errstate(over='ignore'):
            # Compute the native contacts
            kernel = np.exp(beta * (dist - lamda * self._dist))
        frac_of_contacts = np.sum(self._inv_S / (1.0 + kernel), axis=1)
        
        return frac_of_contacts
    
    def print_plumed(self,
            header: bool = True,
            filename: str | os.PathLike[str] | TextIO | None = None) -> None:
        '''
        Print the plumed input file for the native contacts.
        Provide a filename to write to, or None to print to stdout.
        '''

        f: TextIO | None
        should_close = False
        if filename is None:
            f = None
        elif hasattr(filename, "write"):
            f = filename  # already a TextIO-like object
        else:
            f = open(filename, "w")
            should_close = True

        if header:
            print("UNITS LENGTH=A", file=f)  
        print("cmap: CONTACTMAP ...", file=f)
        for i, (a, b) in enumerate(self._native_contacts):
            print((f"    ATOMS{i+1}={a+1},{b+1} "
                   f"SWITCH{i+1}={{Q R_0=0.01 BETA=5.0 LAMBDA=1.8 REF={self._dist[i]*10:.4f}}} " 
                   f"WEIGHT{i+1}={self._inv_S:.7f}"
            ), file=f)
        print("    SUM\n...", file=f)
        if header:
            print("PRINT ARG=cmap FILE=COLVAR.dat", file=f)

        if should_close and f is not None:
            f.close()
