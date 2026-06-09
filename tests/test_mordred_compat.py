from __future__ import annotations

import ast
import json
import math
from pathlib import Path

import pytest
from rdkit import Chem

from descriptors import calc_rdkit_mordred_like_2d
from descriptors.mordred_rdkit_registry import (
    MORDRED_RDKIT_ALIASES,
    SUPPORTED_MORDRED_2D_DESCRIPTORS,
)
from tests.mordred_reference import calc_mordred_2d, mordred_2d_descriptor_names

ABS_TOL = 1e-8
REL_TOL = 1e-6


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
    renamed = set(SUPPORTED_MORDRED_2D_DESCRIPTORS) - {"BertzCT"}
    assert set(MORDRED_RDKIT_ALIASES) == renamed


def test_supported_descriptors_exist_in_mordred():
    mordred_names = mordred_2d_descriptor_names()
    missing = set(SUPPORTED_MORDRED_2D_DESCRIPTORS) - mordred_names
    assert not missing


def test_calculator_rejects_none():
    with pytest.raises(ValueError, match="RDKit Mol"):
        calc_rdkit_mordred_like_2d(None)


@pytest.mark.parametrize(("name", "mol"), _load_validation_molecules())
def test_rdkit_values_match_mordred_reference(name, mol):
    rdkit_values = calc_rdkit_mordred_like_2d(mol)
    reference_values = calc_mordred_2d(mol, SUPPORTED_MORDRED_2D_DESCRIPTORS)

    assert set(rdkit_values) == set(SUPPORTED_MORDRED_2D_DESCRIPTORS)
    for descriptor_name in SUPPORTED_MORDRED_2D_DESCRIPTORS:
        rdkit_value = rdkit_values[descriptor_name]
        reference_value = reference_values[descriptor_name]
        assert isinstance(rdkit_value, int | float)
        assert isinstance(reference_value, int | float), (
            f"{descriptor_name} returned non-numeric Mordred value "
            f"for {name}: {reference_value!r}"
        )
        assert math.isclose(
            rdkit_value,
            reference_value,
            abs_tol=ABS_TOL,
            rel_tol=REL_TOL,
        ), f"{descriptor_name} mismatch for {name}: RDKit={rdkit_value!r}, Mordred={reference_value!r}"
