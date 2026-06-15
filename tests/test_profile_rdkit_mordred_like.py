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


def test_markdown_report_includes_configuration_and_timing_table():
    config = profile_rdkit_mordred_like.ProfileConfig(
        molecule_count=2,
        repeat=3,
        descriptors_per_molecule=5,
    )
    results = (
        profile_rdkit_mordred_like.TimingResult("rdkit", 30, 0.25),
        profile_rdkit_mordred_like.TimingResult("rdkit/chi", 6, 0.125),
    )

    report = profile_rdkit_mordred_like._format_markdown_report(config, results)

    assert report.startswith("# RDKit Mordred-like descriptor profile")
    assert "| setting | value |" in report
    assert "| molecules | 2 |" in report
    assert "| repeats | 3 |" in report
    assert "| descriptors per molecule | 5 |" in report
    assert (
        "| calculator | descriptor values | elapsed seconds | "
        "descriptor values/second |"
    ) in report
    assert "| rdkit | 30 | 0.250000 | 120.00 |" in report
    assert "| rdkit/chi | 6 | 0.125000 | 48.00 |" in report


def test_markdown_report_handles_zero_elapsed_time():
    config = profile_rdkit_mordred_like.ProfileConfig(
        molecule_count=1,
        repeat=1,
        descriptors_per_molecule=1,
    )
    results = (profile_rdkit_mordred_like.TimingResult("rdkit", 1, 0.0),)

    report = profile_rdkit_mordred_like._format_markdown_report(config, results)

    assert "| rdkit | 1 | 0.000000 | inf |" in report
