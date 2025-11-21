from __future__ import annotations

from collections.abc import Sequence
import logging
from pathlib import Path

from af2rave.feature import FeatureSelection
from af2rave.simulation import SimulationBox
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import pickle

from tempo.utils.contact import NativeContact

LOGGER = logging.getLogger(__name__)


def execute(
        prefix: str,
        output_dir: str,
        ref: str | None = None,
        steric_clash_cutoff: float = 1.0,
        ag_chains: Sequence[str] = ("A",),
        ab_chains: Sequence[str] = ("B", "C"),
        include_cb: bool = True,
        n_features: int = 200,
        n_clusters: int = 15,
        random_seed: int = 42,
        qval_cutoff: float = 0.85
):
    LOGGER.info("Starting mdprep workflow for prefix '%s'", prefix)

    if prefix == output_dir:
        raise ValueError("Input and output directories cannot be the same.")
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    fs = _load(prefix, ref, steric_clash_cutoff)
    LOGGER.debug("Loaded feature selection with %d structures", len(fs.pdb_name))

    # Find chain_id (0-based) from one-letter chain names in PDB
    antigen_chains = _find_chains(fs, ag_chains)
    antibody_chains = _find_chains(fs, ab_chains)
    LOGGER.debug("Antigen chains resolved to indices: %s", antigen_chains)
    LOGGER.debug("Antibody chains resolved to indices: %s", antibody_chains)

    # Build native contacts and feature selections
    nc = _build_native_contacts(fs, antigen_chains, antibody_chains)
    selection_string = _get_feature_selection_string(
        fs, nc.contacts_residues, include_cb, antigen_chains, antibody_chains
    )
    LOGGER.debug("Constructed %d native contacts", len(nc.contacts))
    feature_names = _select_features(fs, selection_string, n_features, output_dir)
    LOGGER.info("Selected top %d features", len(feature_names))

    # remove structures with q-val too low
    qval = _filter_by_qval(fs, nc, qval_cutoff)
    LOGGER.info("Filtered structures; %d remain above qval %.2f", len(qval), qval_cutoff)

    # cluster into centers
    center_id = _cluster(fs, feature_names, n_clusters, 
                         output_dir=output_dir, random_seed=random_seed, qval=qval)
    LOGGER.info("Identified %d cluster centers", len(center_id))

    # Save PDBs
    _box_structures(fs, feature_names, center_id, output_dir)
    LOGGER.info("Boxed structures written to %s", output_dir)

    
def _load(prefix: str, ref: str | None = None, steric_clash_cutoff: float = 1.0) -> FeatureSelection:
    """Return a FeatureSelection, loading/saving a cache to avoid reprocessing."""
    cache_path = Path(prefix, "feature_selection.pkl")
    if cache_path.exists():
        LOGGER.debug("Loading FeatureSelection cache from %s", cache_path)
        with cache_path.open("rb") as handle:
            cached = pickle.load(handle)
        if isinstance(cached, FeatureSelection) and len(cached.pdb_name) > 0:
            return cached
        LOGGER.warning("Cached FeatureSelection invalid; rebuilding from prefix %s", prefix)

    LOGGER.info("Building FeatureSelection from prefix %s", prefix)
    fs = FeatureSelection(prefix, ref_pdb=ref)
    fs.apply_filter(fs.steric_clash_filter(steric_clash_cutoff))
    try:
        with cache_path.open("wb") as handle:
            pickle.dump(fs, handle)
        LOGGER.debug("FeatureSelection cached to %s", cache_path)
    except OSError as exc:
        LOGGER.warning("Unable to write FeatureSelection cache %s: %s", cache_path, exc)
    return fs


def _find_chains(fs: FeatureSelection, letters: Sequence[str]):

    top = fs.top
    result = []
    for l in letters:
        chain_id = None
        for c in top.chains:
            if c.chain_id == l:
                chain_id = c.index
                break
        if chain_id is None:
            raise ValueError(f"Chain ID {l} not found in topology.")
        else:
            result.append(chain_id)
    return result


def _build_native_contacts(
        fs: FeatureSelection,
        antigen_chains: Sequence[int] = (0,),
        antibody_chains: Sequence[int] = (1, 2)
) -> NativeContact:
    """Construct a NativeContact object limited to the heavy atoms of the provided chains."""
    antigen_sel = " ".join(str(chain) for chain in antigen_chains)
    antibody_sel = " ".join(str(chain) for chain in antibody_chains)

    nc = NativeContact(
        fs._ref,
        [
            f"chainid {antigen_sel} and mass > 1.1",
            f"chainid {antibody_sel} and mass > 1.1",
        ],
    )

    return nc

def _get_feature_selection_string(
    fs: FeatureSelection,
    native_contact_residues: dict[str, Sequence[int]],
    include_cb: bool = True,
    antigen_chains: Sequence[int] = (0,),
    antibody_chains: Sequence[int] = (1, 2),
) -> tuple[str, str]:
    """Build MDTraj selection strings for antigen/antibody contacts derived from a reference."""
    chain_index_to_id = {chain.index: chain.chain_id for chain in fs._ref.topology.chains}

    def format_chain_clause(chain_index: int) -> str | None:
        chain_id = chain_index_to_id.get(chain_index)
        if chain_id is None:
            return None
        residues = native_contact_residues.get(chain_id)
        if not residues:
            return None
        residue_list = " ".join(str(residue) for residue in residues)
        return f"(chainid {chain_index} and resid {residue_list})"

    def combine_chain_clauses(chain_indices: Sequence[int]) -> str:
        clauses = [clause for idx in chain_indices if (clause := format_chain_clause(idx))]
        if not clauses:
            raise ValueError(
                "Native-contact selection produced no residues for chains "
                f"{', '.join(map(str, chain_indices))}."
            )
        return " or ".join(clauses)

    prefix = "name CA CB" if include_cb else "name CA"
    selection_ag = f"{prefix} and ({combine_chain_clauses(antigen_chains)})"
    selection_ab = f"{prefix} and ({combine_chain_clauses(antibody_chains)})"
    return selection_ag, selection_ab


def _select_features(
        fs: FeatureSelection,
        selection: tuple[str, str],
        n_features: int = 200,
        output_dir: str | None = None
    ) -> list[str]:

    LOGGER.debug(selection)
    names, cv = fs.rank_feature(selection)

    if output_dir:

        fig, ax = plt.subplots(figsize=(3, 2), tight_layout=True, dpi=300)
        counts, _, _ = ax.hist(cv, bins=100, label=f"{len(names)}-dim", density=True)
        ax.vlines(cv[n_features], 0, np.max(counts) * 0.6, 
                  color="red", linestyle="--", linewidth=1, 
                  label=f"Antigen={cv[n_features]:.3f}"
        )

        ax.set_ylim(0, np.max(counts) * 1.2)
        ax.set_xlabel("Coefficient of Variation")
        ax.set_ylabel("Density")
        ax.set_xlim(0, 0.5)
        ax.legend()
        fig.savefig(output_dir + "/feature_selection.png")

    return names[:n_features]


def _filter_by_qval(
        fs: FeatureSelection,
        nc: NativeContact,
        qval_cutoff: float = 0.85,
) -> list[str]:

    pdb_names = fs.pdb_name
    qval = nc.compute(fs.traj)
    valid_id = [i for i, q in enumerate(qval) if q > qval_cutoff]
    filter = [pdb_names[i] for i in valid_id]
    fs.apply_filter(filter)
    return qval[valid_id]

def _cluster(
        fs: FeatureSelection,
        names: Sequence[str],
        n_clusters: int = 15,
        dimensions: int = 8,
        random_seed: int = 42,
        output_dir: str | None = None,
        **plot_kwargs
):
    '''
    Perform PCA to reduce the feature space to a lower-dimensional subspace
    and do K-centroid clustering in it to get representative structures.
    '''

    _, pca_result = _run_pca_in_subspace(fs, names, n_components=dimensions)
    center_id = _get_pca_centroids(pca_result, n_clusters, random_seed)

    if output_dir:
        fig, ax = plt.subplots(figsize=(4,3), dpi=300)
        ax.set_box_aspect(1)

        scatter_kwargs = {"s": 3, 
                         "c": plot_kwargs.get("qval", "k"),
                         "cmap": "gnuplot", 
                         "norm": "linear", 
                         "edgecolors": "none",
                         "vmin": 0
                         }

        cbar = ax.scatter(pca_result[:,0], pca_result[:,1], **scatter_kwargs)
        cb_ax = fig.add_axes([.85,.124,.02,.754])
        fig.colorbar(cbar, orientation='vertical', cax=cb_ax)
        cb_ax.set_ylabel("Q-val from reference structure", loc="center")

        # cluster center in PC space
        pca_cc = pca_result[center_id]

        cc_kwargs = {"s": 10, 
                    "edgecolor": "black",
                    "marker": "H",
                    "facecolor": "None", 
                    "linewidth": 0.5}
        ax.scatter(pca_cc[:,0], pca_cc[:,1], **cc_kwargs)
        for c in center_id:
            ax.text(pca_result[c,0], pca_result[c,1], fs.pdb_name[c].split("/")[-1], fontsize=3)
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        fig.savefig(output_dir + "/pca_cluster.png")

    return center_id


def _run_pca_in_subspace(
        fs: FeatureSelection, 
        feature_names: Sequence[str], 
        n_components: int = 2, 
        **kwargs
) -> tuple[PCA, np.ndarray]:

    # Extract time series data from features
    z = np.array([fs.features[fn] for fn in feature_names], dtype=np.float64).T
    LOGGER.debug("Feature data shape before PCA: %s", z.shape)

    # Ensure there are enough features to compute the requested components
    if z.shape[1] < n_components:
        raise ValueError(f"Number of components ({n_components}) cannot exceed available features ({z.shape[1]}).")

    # Perform PCA
    pca = PCA(n_components=n_components, **kwargs)
    transformed_data = pca.fit_transform(z)

    return pca, transformed_data


def _get_pca_centroids(coord, n_clusters=15, random_seed=42) -> np.ndarray:

    np.random.seed(random_seed)
    random_indices = np.random.choice(coord.shape[0], n_clusters, replace=False)
    initial_centroids = coord[random_indices]

    kmeans = KMeans(n_clusters=n_clusters, init=initial_centroids, n_init=1)
    kmeans.fit(coord)

    cid = np.asarray([
        np.linalg.norm(coord - cc, axis=1).argmin() 
        for cc in kmeans.cluster_centers_
    ])

    return cid


def _box_structures(
        fs: FeatureSelection,
        feature_names: Sequence[str],
        center_id: Sequence[int],
        output_dir: str
):

    pdb_names = [fs.pdb_name[i] for i in center_id]
    atom_index = [fs.atom_pairs[n] for n in feature_names]

    for p in pdb_names:

        prefix = p.split("/")[-1]
        
        box = SimulationBox(p)
        box.create_box()
        new_index = box.map_atom_index(atom_index)
        np.save(f"{output_dir}/index.npy", new_index)

        box.save_pdb(f"{output_dir}/{prefix}")
        LOGGER.info(f"Wrote boxed structure {prefix}")
