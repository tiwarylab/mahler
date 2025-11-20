from __future__ import annotations

from pathlib import Path

import mdtraj as md
import pytest

pytest.importorskip("af2rave.feature")
from af2rave.feature import FeatureSelection

from tempo.commands.mdprep import feature
from tempo.utils.contact import NativeContact


@pytest.fixture(scope="module")
def example_pdb_path() -> Path:
    return Path(__file__).with_name("example.pdb")


@pytest.fixture(scope="module")
def feature_selection(example_pdb_path: Path) -> FeatureSelection:
    return feature._load(str(example_pdb_path), ref=str(example_pdb_path), steric_clash_cutoff=0.5)


def test_load_returns_real_feature_selection(example_pdb_path: Path, feature_selection: FeatureSelection) -> None:
    assert isinstance(feature_selection, FeatureSelection)
    assert feature_selection.ref_pdb == str(example_pdb_path)
    assert feature_selection.pdb_name == [str(example_pdb_path)]
    assert len(feature_selection.traj) == 1


def test_get_feature_selection_string_matches_native_contacts(feature_selection: FeatureSelection) -> None:
    antigen = (0,)
    antibody = (1, 2)

    selection_ag, selection_ab = feature._get_feature_selection_string(
        feature_selection, include_cb=False, antigen_chains=antigen, antibody_chains=antibody
    )

    nc = NativeContact(
        feature_selection._ref,
        [
            "chainid 0 and mass > 1.1",
            "chainid 1 2 and mass > 1.1",
        ],
    )
    residues = nc.contacts_residues
    chain_index_to_id = {chain.index: chain.chain_id for chain in feature_selection._ref.topology.chains}

    def build_expected(indices: tuple[int, ...]) -> str:
        clauses = []
        for idx in indices:
            chain_id = chain_index_to_id[idx]
            resid_list = " ".join(str(r) for r in residues[chain_id])
            clauses.append(f"(chainid {idx} and resid {resid_list})")
        return f"name CA and ({' or '.join(clauses)})"

    assert selection_ag == build_expected(antigen)
    assert selection_ab == build_expected(antibody)
