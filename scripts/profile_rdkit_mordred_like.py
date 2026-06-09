"""Profile RDKit-only and Mordred descriptor calculators."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from time import perf_counter

from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from descriptors import calc_rdkit_mordred_like_2d
from descriptors.mordred_rdkit_registry import SUPPORTED_MORDRED_2D_DESCRIPTORS


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
    args = parser.parse_args()

    molecules = _load_molecules(args.panel)
    print(f"molecules: {len(molecules)}")
    print(f"repeats: {args.repeat}")
    print(f"descriptors per molecule: {len(SUPPORTED_MORDRED_2D_DESCRIPTORS)}")

    if args.backend in {"rdkit", "both"}:
        descriptor_count, elapsed = _time_rdkit(molecules, args.repeat)
        _print_result("rdkit", descriptor_count, elapsed)

    if args.backend in {"mordred", "both"}:
        descriptor_count, elapsed = _time_mordred(molecules, args.repeat)
        _print_result("mordred", descriptor_count, elapsed)


if __name__ == "__main__":
    main()
