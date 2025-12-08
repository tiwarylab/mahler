import mdtraj as md
import logging

from collections.abc import Sequence

LOGGER = logging.getLogger("mahler.topology")

def check_chains(traj: md.Trajectory) -> str:
    """Check the chains existing in the trajectory and fix problems."""
    top: md.Topology = traj.topology
    idx_protein = top.select("protein")

    top_pro = top.subset(idx_protein)
    pro_chain_ids = [chain.chain_id.strip() for chain in top_pro.chains]
    if "" not in pro_chain_ids:             # we are good to go
        chains = " ".join(pro_chain_ids)
        LOGGER.info(f"Protein chains found: {chains}")
        return chains

    pro_chain_ids = []
    LOGGER.warning("Some protein chains have no chain IDs. Renaming them...")
    for i, c in enumerate(top.chains):
        # check if this chain is a protein chain
        atom_indices = [atom.index for atom in c.atoms]
        if all(idx in idx_protein for idx in atom_indices):
            c.chain_id = chr(ord('A') + i)
            pro_chain_ids.append(c.chain_id)
            LOGGER.info(f"Renaming protein chain {c.index} to '{c.chain_id}'.")
        else:
            LOGGER.info(f"Skipping non-protein chain {c.index}.")
    return " ".join(pro_chain_ids)

def chain_idx_from_chain_id(
        top: md.Topology,
        protein_only: bool = True,
) -> list[int]:
    """Get the atom indices for a given chain ID."""

    if protein_only:
        idx_protein = top.select("protein")
        top_pro = top.subset(idx_protein)
    else:
        top_pro = top

    return {
        c.chain_id: c.index for c in top_pro.chains
    }
