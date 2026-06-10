from __future__ import annotations

import ast
import json
import math
from pathlib import Path

import pytest
from rdkit import Chem

from descriptors import calc_rdkit_mordred_like_2d
from descriptors import rdkit_mordred_like
from descriptors.mordred_rdkit_registry import (
    EXACT_NAME_RDKIT_DESCRIPTORS,
    MORDRED_RDKIT_ALIASES,
    SUPPORTED_MORDRED_2D_DESCRIPTORS,
)
from tests.mordred_reference import calc_mordred_2d, mordred_2d_descriptor_names

ABS_TOL = 1e-8
REL_TOL = 1e-6
EXPECTED_MORDRED_MISSING_VALUES = {
    ("ammonium", "PetitjeanIndex"): "nan",
    ("ammonium", "RotRatio"): "nan",
    ("ammonium", "TopoShapeIndex"): "nan",
    ("ammonium", "Xp-0d"): "numeric",
    ("ammonium", "mZagreb1"): "nan",
    ("methane", "PetitjeanIndex"): "nan",
    ("methane", "RotRatio"): "nan",
    ("methane", "TopoShapeIndex"): "nan",
    ("methane", "Xp-0d"): "numeric",
    ("methane", "mZagreb1"): "nan",
    ("tetrahalo_methane", "Vabc"): "nan",
}


def _load_validation_molecules():
    panel_path = Path(__file__).with_name("smiles_panel.smi")
    molecules = []
    for line in panel_path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        smiles, name = line.split(maxsplit=1)
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"invalid validation SMILES for {name}: {smiles}"
        molecules.append((name, mol))
    return molecules


def test_production_modules_do_not_import_mordred():
    package_root = Path(__file__).resolve().parents[1] / "descriptors"
    for path in package_root.rglob("*.py"):
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported = [node.module]
            else:
                continue
            assert not any(
                name == "mordred"
                or name.startswith("mordred.")
                or name == "mordredcommunity"
                or name.startswith("mordredcommunity.")
                for name in imported
            ), f"production module {path} imports a Mordred package"


def test_supported_descriptor_list_is_locked():
    expected_path = Path(__file__).with_name("expected_supported.json")
    expected = tuple(json.loads(expected_path.read_text()))
    assert SUPPORTED_MORDRED_2D_DESCRIPTORS == expected


def test_alias_registry_covers_renamed_descriptors():
    renamed = set(SUPPORTED_MORDRED_2D_DESCRIPTORS) - set(EXACT_NAME_RDKIT_DESCRIPTORS)
    assert set(MORDRED_RDKIT_ALIASES) == renamed


def test_supported_descriptors_exist_in_mordred():
    mordred_names = mordred_2d_descriptor_names()
    missing = set(SUPPORTED_MORDRED_2D_DESCRIPTORS) - mordred_names
    assert not missing


def test_calculator_rejects_none():
    with pytest.raises(ValueError, match="RDKit Mol"):
        calc_rdkit_mordred_like_2d(None)


def test_calculator_reuses_expensive_context_values(monkeypatch):
    counts = {
        "distance": 0,
        "adjacency": 0,
        "rings": 0,
        "ring_filters": 0,
        "fused": 0,
        "hydrogen": 0,
        "kekulized": 0,
    }
    context_cls = rdkit_mordred_like._DescriptorContext

    def count_calls(method_name, count_name):
        original = getattr(context_cls, method_name)

        def counted(self):
            counts[count_name] += 1
            return original(self)

        monkeypatch.setattr(context_cls, method_name, counted)

    count_calls("_compute_distance_matrix", "distance")
    count_calls("_compute_adjacency_matrix", "adjacency")
    count_calls("_compute_ring_atom_sets", "rings")
    count_calls("_compute_fused_ring_systems", "fused")
    count_calls("_compute_implicit_hydrogen_count", "hydrogen")
    count_calls("_compute_kekulized_mol", "kekulized")
    original_ring_filter = rdkit_mordred_like._ring_matches_filters

    def counted_ring_filter(*args, **kwargs):
        counts["ring_filters"] += 1
        return original_ring_filter(*args, **kwargs)

    monkeypatch.setattr(
        rdkit_mordred_like,
        "_ring_matches_filters",
        counted_ring_filter,
    )

    mol = Chem.MolFromSmiles("c1ccc2ccccc2c1")
    values = calc_rdkit_mordred_like_2d(mol)

    assert values["WPath"] > 0
    assert values["nFRing"] == 1
    assert counts == {
        "distance": 1,
        "adjacency": 1,
        "rings": 1,
        "ring_filters": 0,
        "fused": 1,
        "hydrogen": 1,
        "kekulized": 1,
    }


@pytest.mark.parametrize(("name", "mol"), _load_validation_molecules())
def test_rdkit_values_match_mordred_reference(name, mol):
    rdkit_values = calc_rdkit_mordred_like_2d(mol)
    reference_values = calc_mordred_2d(mol, SUPPORTED_MORDRED_2D_DESCRIPTORS)

    assert set(rdkit_values) == set(SUPPORTED_MORDRED_2D_DESCRIPTORS)
    for descriptor_name in SUPPORTED_MORDRED_2D_DESCRIPTORS:
        rdkit_value = rdkit_values[descriptor_name]
        reference_value = reference_values[descriptor_name]
        assert isinstance(rdkit_value, int | float)
        if not isinstance(reference_value, int | float):
            expected = EXPECTED_MORDRED_MISSING_VALUES.get((name, descriptor_name))
            assert expected is not None, (
                f"{descriptor_name} returned undocumented non-numeric Mordred "
                f"value for {name}: {reference_value!r}"
            )
            if expected == "numeric":
                assert not math.isnan(float(rdkit_value)), (
                    f"{descriptor_name} should provide an RDKit numeric value "
                    f"for {name} when Mordred is missing"
                )
            elif expected == "nan":
                assert math.isnan(float(rdkit_value)), (
                    f"{descriptor_name} should return NaN for {name} "
                    f"when both implementations are undefined"
                )
            else:
                raise AssertionError(f"unknown missing-value expectation: {expected}")
            continue
        assert math.isclose(
            rdkit_value,
            reference_value,
            abs_tol=ABS_TOL,
            rel_tol=REL_TOL,
        ), f"{descriptor_name} mismatch for {name}: RDKit={rdkit_value!r}, Mordred={reference_value!r}"
