"""Test-only Mordred oracle for descriptor compatibility checks."""

from __future__ import annotations

from functools import lru_cache
from typing import Iterable


def _load_mordred():
    try:
        from mordredcommunity import Calculator, descriptors
    except ModuleNotFoundError:
        from mordred import Calculator, descriptors

    return Calculator, descriptors


@lru_cache(maxsize=1)
def _all_2d_descriptors_by_name():
    Calculator, descriptors = _load_mordred()
    calc = Calculator(descriptors, ignore_3D=True)
    return {str(descriptor): descriptor for descriptor in calc.descriptors}


def mordred_2d_descriptor_names() -> set[str]:
    return set(_all_2d_descriptors_by_name())


def calc_mordred_2d(mol, names: Iterable[str]) -> dict[str, float | int]:
    Calculator, _ = _load_mordred()
    descriptors_by_name = _all_2d_descriptors_by_name()
    selected = [descriptors_by_name[name] for name in names]
    result = Calculator(selected, ignore_3D=True)(mol).asdict()
    return {str(name): value for name, value in result.items()}

