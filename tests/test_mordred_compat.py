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
# Validation-panel cases where Mordred returns a missing value but the RDKit-only
# calculator deliberately provides a meaningful numeric value. These are the only
# Mordred-missing cases that must be documented: when Mordred is missing and RDKit
# is also NaN, both implementations agree the descriptor is undefined and there is
# no oracle value to check, so those cases are accepted without enumeration.
EXPECTED_RDKIT_IMPROVEMENTS: frozenset[tuple[str, str]] = frozenset()


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


def test_unsupported_descriptor_list_is_locked():
    unsupported_path = Path(__file__).with_name("unsupported_mordred.json")
    expected = tuple(json.loads(unsupported_path.read_text()))
    all_mordred = mordred_2d_descriptor_names()
    unsupported = tuple(sorted(n for n in all_mordred if n not in set(SUPPORTED_MORDRED_2D_DESCRIPTORS)))
    assert unsupported == expected


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


def test_calculator_names_subset_matches_full_run():
    mol = Chem.MolFromSmiles("c1ccc2ccccc2c1")
    full = calc_rdkit_mordred_like_2d(mol)
    requested = ["Xp-2d", "ATS0Z", "MW", "nFRing"]
    subset = calc_rdkit_mordred_like_2d(mol, names=requested)
    assert list(subset) == requested
    for name in requested:
        assert subset[name] == full[name]


def test_calculator_rejects_unknown_names():
    mol = Chem.MolFromSmiles("CCO")
    with pytest.raises(KeyError, match="unsupported descriptor name"):
        calc_rdkit_mordred_like_2d(mol, names=["MW", "not_a_descriptor"])


def test_calculator_names_subset_skips_unrelated_families(monkeypatch):
    # Requesting only chi descriptors must not trigger the autocorrelation
    # family computation, proving subset evaluation is demand-driven.
    context_cls = rdkit_mordred_like._DescriptorContext
    original = context_cls.autocorrelation_values.func
    calls = {"n": 0}

    def spy(self):
        calls["n"] += 1
        return original(self)

    # A plain property is enough to observe access without the cached_property
    # __set_name__ binding; caching is irrelevant when the expected count is 0.
    monkeypatch.setattr(context_cls, "autocorrelation_values", property(spy))

    mol = Chem.MolFromSmiles("CCO")
    calc_rdkit_mordred_like_2d(mol, names=["Xp-2d", "Xp-3d"])
    assert calls["n"] == 0


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
            # Mordred is missing for this molecule/descriptor.
            if (name, descriptor_name) in EXPECTED_RDKIT_IMPROVEMENTS:
                assert not math.isnan(float(rdkit_value)), (
                    f"{descriptor_name} is documented as an RDKit improvement for "
                    f"{name} but returned NaN while Mordred is missing"
                )
            else:
                assert math.isnan(float(rdkit_value)), (
                    f"{descriptor_name} is numeric for {name} while Mordred is "
                    f"missing; if this is a deliberate RDKit improvement, add "
                    f"({name!r}, {descriptor_name!r}) to EXPECTED_RDKIT_IMPROVEMENTS"
                )
            continue
        assert math.isclose(
            rdkit_value,
            reference_value,
            abs_tol=ABS_TOL,
            rel_tol=REL_TOL,
        ), f"{descriptor_name} mismatch for {name}: RDKit={rdkit_value!r}, Mordred={reference_value!r}"
