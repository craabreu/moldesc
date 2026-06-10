"""Profile RDKit-only and Mordred descriptor calculators."""

from __future__ import annotations

import argparse
from collections.abc import Iterable
from typing import NamedTuple
from pathlib import Path
import sys
from time import perf_counter

from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from descriptors import calc_rdkit_mordred_like_2d, rdkit_mordred_like
from descriptors.mordred_rdkit_registry import (
    AUTOCORRELATION_DESCRIPTORS,
    BCUT_DESCRIPTORS,
    CHI_DESCRIPTORS,
    CONSTITUTIONAL_DESCRIPTORS,
    ESTATE_ATOM_TYPE_DESCRIPTORS,
    TOPOLOGICAL_CHARGE_DESCRIPTORS,
    EXACT_NAME_RDKIT_DESCRIPTORS,
    GRAPH_TOPOLOGY_DESCRIPTORS,
    PATH_COUNT_DESCRIPTORS,
    PHYSICAL_PROPERTY_DESCRIPTORS,
    RING_COUNT_DESCRIPTORS,
    SMALL_GRAPH_FORMULA_DESCRIPTORS,
    SPECTRAL_DESCRIPTORS,
    SUPPORTED_MORDRED_2D_DESCRIPTORS,
    WALK_COUNT_DESCRIPTORS,
)


class TimingResult(NamedTuple):
    label: str
    descriptor_count: int
    elapsed: float


_NAMED_DESCRIPTOR_GROUPS: dict[str, tuple[str, ...]] = {
    "exact_rdkit": EXACT_NAME_RDKIT_DESCRIPTORS,
    "estate_atom_types": ESTATE_ATOM_TYPE_DESCRIPTORS,
    "ring_counts": RING_COUNT_DESCRIPTORS,
    "graph_topology": GRAPH_TOPOLOGY_DESCRIPTORS,
    "path_counts": PATH_COUNT_DESCRIPTORS,
    "walk_counts": WALK_COUNT_DESCRIPTORS,
    "autocorrelation": AUTOCORRELATION_DESCRIPTORS,
    "bcut": BCUT_DESCRIPTORS,
    "constitutional": CONSTITUTIONAL_DESCRIPTORS,
    "topological_charge": TOPOLOGICAL_CHARGE_DESCRIPTORS,
    "chi": CHI_DESCRIPTORS,
    "spectral": SPECTRAL_DESCRIPTORS,
    "small_graph_formula": SMALL_GRAPH_FORMULA_DESCRIPTORS,
    "physical_properties": PHYSICAL_PROPERTY_DESCRIPTORS,
}

_GROUPED_DESCRIPTOR_NAMES = {
    descriptor
    for descriptors in _NAMED_DESCRIPTOR_GROUPS.values()
    for descriptor in descriptors
}

RDKIT_DESCRIPTOR_GROUPS: dict[str, tuple[str, ...]] = {
    **_NAMED_DESCRIPTOR_GROUPS,
    "scalar_aliases": tuple(
        name
        for name in SUPPORTED_MORDRED_2D_DESCRIPTORS
        if name not in _GROUPED_DESCRIPTOR_NAMES
    ),
}


def _load_molecules(panel_path: Path):
    molecules = []
    for line in panel_path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        smiles, name = line.split(maxsplit=1)
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            msg = f"invalid SMILES for {name}: {smiles}"
            raise ValueError(msg)
        molecules.append((name, mol))
    return molecules


def _load_mordred_calculator():
    try:
        from mordredcommunity import Calculator, descriptors
    except ModuleNotFoundError:
        from mordred import Calculator, descriptors

    all_descriptors = {
        str(descriptor): descriptor
        for descriptor in Calculator(descriptors, ignore_3D=True).descriptors
    }
    selected = [
        all_descriptors[name]
        for name in SUPPORTED_MORDRED_2D_DESCRIPTORS
    ]
    return Calculator(selected, ignore_3D=True)


def _time_rdkit(molecules, repeat: int) -> tuple[int, float]:
    descriptor_count = 0
    start = perf_counter()
    for _ in range(repeat):
        for _, mol in molecules:
            descriptor_count += len(calc_rdkit_mordred_like_2d(mol))
    return descriptor_count, perf_counter() - start


def _time_rdkit_descriptor_names(
    molecules,
    repeat: int,
    descriptor_names: Iterable[str],
) -> tuple[int, float]:
    names = tuple(descriptor_names)
    functions = tuple(rdkit_mordred_like._DESCRIPTOR_FUNCTIONS[name] for name in names)
    descriptor_count = 0
    start = perf_counter()
    for _ in range(repeat):
        for _, mol in molecules:
            context = rdkit_mordred_like._DescriptorContext(mol)
            for function in functions:
                function(context)
            descriptor_count += len(functions)
    return descriptor_count, perf_counter() - start


def _time_rdkit_groups(molecules, repeat: int) -> tuple[TimingResult, ...]:
    return tuple(
        TimingResult(
            group_name,
            *_time_rdkit_descriptor_names(molecules, repeat, descriptor_names),
        )
        for group_name, descriptor_names in RDKIT_DESCRIPTOR_GROUPS.items()
    )


def _time_mordred(molecules, repeat: int) -> tuple[int, float]:
    calculator = _load_mordred_calculator()
    descriptor_count = 0
    start = perf_counter()
    for _ in range(repeat):
        for _, mol in molecules:
            descriptor_count += len(calculator(mol).asdict())
    return descriptor_count, perf_counter() - start


def _print_result(label: str, descriptor_count: int, elapsed: float) -> None:
    print(f"{label} descriptor values: {descriptor_count}")
    print(f"{label} elapsed seconds: {elapsed:.6f}")
    print(f"{label} descriptor values/second: {descriptor_count / elapsed:.2f}")


def _print_group_results(results: tuple[TimingResult, ...]) -> None:
    print("rdkit grouped timings:")
    for result in results:
        _print_result(f"rdkit {result.label}", result.descriptor_count, result.elapsed)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeat", type=int, default=100)
    parser.add_argument(
        "--backend",
        choices=("rdkit", "mordred", "both"),
        default="rdkit",
    )
    parser.add_argument(
        "--panel",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "tests" / "smiles_panel.smi",
    )
    parser.add_argument(
        "--grouped",
        action="store_true",
        help="also time RDKit descriptors by implementation family",
    )
    args = parser.parse_args()

    molecules = _load_molecules(args.panel)
    print(f"molecules: {len(molecules)}")
    print(f"repeats: {args.repeat}")
    print(f"descriptors per molecule: {len(SUPPORTED_MORDRED_2D_DESCRIPTORS)}")

    if args.backend in {"rdkit", "both"}:
        descriptor_count, elapsed = _time_rdkit(molecules, args.repeat)
        _print_result("rdkit", descriptor_count, elapsed)
        if args.grouped:
            _print_group_results(_time_rdkit_groups(molecules, args.repeat))

    if args.backend in {"mordred", "both"}:
        descriptor_count, elapsed = _time_mordred(molecules, args.repeat)
        _print_result("mordred", descriptor_count, elapsed)


if __name__ == "__main__":
    main()
