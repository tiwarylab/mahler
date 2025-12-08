from __future__ import annotations

from pathlib import Path

import mdtraj as md
import numpy as np
import pytest

from mahler.utils.contact import NativeContact


@pytest.fixture(scope="module")
def example_traj() -> md.Trajectory:
    pdb_path = Path(__file__).with_name("example.pdb")
    return md.load(pdb_path)


@pytest.fixture(scope="module")
def native_contact(example_traj: md.Trajectory) -> NativeContact:
    selection = ["chainid 0 and mass > 1.1", "chainid 1 2 and mass > 1.1"]
    return NativeContact(example_traj, selection)


def test_contacts_property_matches_reference(native_contact: NativeContact, example_traj: md.Trajectory) -> None:
    assert native_contact.contacts.shape == (301, 2)
    np.testing.assert_array_equal(
        native_contact.contacts[:5],
        np.array(
            [
                [224, 3019],
                [226, 3019],
                [228, 3019],
                [229, 3019],
                [229, 3021],
            ],
            dtype=int,
        ),
    )

    frac = native_contact.compute(example_traj)
    assert frac.shape == (1,)
    assert np.allclose(frac[0], 1.0, rtol=0, atol=1e-4)


def test_contacts_atoms_property_returns_sorted_unique_atoms(native_contact: NativeContact) -> None:
    atoms = native_contact.contacts_atoms
    assert atoms.shape == (183,)
    assert np.all(atoms[:-1] <= atoms[1:])

    expected = np.unique(native_contact.contacts.flatten())
    np.testing.assert_array_equal(atoms, expected)


def test_contacts_residues_property_groups_atoms_by_chain(native_contact: NativeContact, example_traj: md.Trajectory) -> None:
    residues = native_contact.contacts_residues
    assert set(residues) == {"A", "B", "C"}

    # Ensure each chain entry is sorted and deduplicated
    for residue_list in residues.values():
        assert residue_list == sorted(residue_list)
        assert len(residue_list) == len(set(residue_list))

    expected: dict[str, set[int]] = {}
    for atom_idx in native_contact.contacts_atoms:
        atom = example_traj.topology.atom(int(atom_idx))
        chain_id = atom.residue.chain.chain_id
        expected.setdefault(chain_id, set()).add(atom.residue.index)

    expected_sorted = {chain: sorted(indices) for chain, indices in expected.items()}
    assert residues == expected_sorted


def test_contacts_ca_property_returns_sorted_numpy_arrays(native_contact: NativeContact, example_traj: md.Trajectory) -> None:
    ca_map = native_contact.contacts_ca
    assert set(ca_map) == {"A", "B", "C"}

    for chain, indices in ca_map.items():
        assert isinstance(indices, list)
        assert all(isinstance(idx, int) for idx in indices)
        assert all(indices[i] <= indices[i + 1] for i in range(len(indices) - 1))

    # Validate against residues participating in native contacts
    expected: dict[str, list[int]] = {chain: [] for chain in native_contact.contacts_residues}
    for chain, residue_indices in native_contact.contacts_residues.items():
        for residue_idx in residue_indices:
            residue = example_traj.topology.residue(residue_idx)
            ca_atom = next((atom.index for atom in residue.atoms if atom.name == "CA"), None)
            if ca_atom is not None:
                expected[chain].append(ca_atom)

    expected_sorted = {chain: sorted(indices) for chain, indices in expected.items()}
    assert ca_map == expected_sorted
