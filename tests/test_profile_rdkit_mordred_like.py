from __future__ import annotations

from rdkit import Chem

from descriptors.mordred_rdkit_registry import SUPPORTED_MORDRED_2D_DESCRIPTORS
from scripts import profile_rdkit_mordred_like


def test_descriptor_groups_cover_supported_descriptors_once():
    grouped = [
        descriptor
        for descriptors in profile_rdkit_mordred_like.RDKIT_DESCRIPTOR_GROUPS.values()
        for descriptor in descriptors
    ]

    assert sorted(grouped) == sorted(SUPPORTED_MORDRED_2D_DESCRIPTORS)
    assert len(grouped) == len(set(grouped))


def test_grouped_rdkit_timing_calculates_descriptor_values():
    mol = Chem.MolFromSmiles("CCO")

    results = profile_rdkit_mordred_like._time_rdkit_groups([("ethanol", mol)], 1)

    assert sum(result.descriptor_count for result in results) == len(
        SUPPORTED_MORDRED_2D_DESCRIPTORS
    )
    assert all(result.elapsed >= 0 for result in results)
