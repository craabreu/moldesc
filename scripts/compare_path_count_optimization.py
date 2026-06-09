"""Compare the current path-count implementation with the legacy algorithm."""

from __future__ import annotations

import argparse
from collections.abc import Callable
from pathlib import Path
import statistics
import sys
from time import perf_counter

from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from descriptors import calc_rdkit_mordred_like_2d
from descriptors import rdkit_mordred_like


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


def _legacy_bond_ids_to_atom_ids(ctx, path) -> tuple[int, ...]:
    path_iter = iter(path)

    try:
        atom0_from, atom0_to = ctx.bond_atom_pairs[next(path_iter)]
    except StopIteration:
        return ()

    try:
        atom1_from, atom1_to = ctx.bond_atom_pairs[next(path_iter)]
    except StopIteration:
        return atom0_from, atom0_to

    if atom0_from in [atom1_from, atom1_to]:
        atoms = [atom0_to, atom0_from]
        current = atom1_from if atom0_from == atom1_to else atom1_to
    else:
        atoms = [atom0_from, atom0_to]
        current = atom1_from if atom0_to == atom1_to else atom1_to

    for bond_index in path_iter:
        atom_from, atom_to = ctx.bond_atom_pairs[bond_index]
        atoms.append(current)

        if atom_from == current:
            current = atom_to
        else:
            current = atom_from

    atoms.append(current)
    return tuple(atoms)


def _legacy_compute_path_count(ctx, order: int) -> tuple[int, float]:
    path_count = 0
    pi_path_count = 0.0

    for path in Chem.FindAllPathsOfLengthN(ctx.mol, order):
        atom_ids = set()
        previous = None
        pi_weight = 1.0

        for atom_index in _legacy_bond_ids_to_atom_ids(ctx, path):
            if atom_index in atom_ids:
                break

            atom_ids.add(atom_index)

            if previous is not None:
                bond = ctx.mol.GetBondBetweenAtoms(previous, atom_index)
                pi_weight *= bond.GetBondTypeAsDouble()

            previous = atom_index
        else:
            path_count += 1
            pi_path_count += pi_weight

    return path_count, pi_path_count


def _time_calculator(molecules, repeat: int) -> tuple[int, float]:
    descriptor_count = 0
    start = perf_counter()
    for _ in range(repeat):
        for _, mol in molecules:
            descriptor_count += len(calc_rdkit_mordred_like_2d(mol))
    return descriptor_count, perf_counter() - start


def _run_mode(
    label: str,
    molecules,
    repeat: int,
    path_counter: Callable | None,
) -> tuple[int, float]:
    context_cls = rdkit_mordred_like._DescriptorContext
    original = context_cls._compute_path_count
    if path_counter is not None:
        context_cls._compute_path_count = path_counter
    try:
        descriptor_count, elapsed = _time_calculator(molecules, repeat)
    finally:
        context_cls._compute_path_count = original

    print(f"{label} descriptor values: {descriptor_count}")
    print(f"{label} elapsed seconds: {elapsed:.6f}")
    print(f"{label} descriptor values/second: {descriptor_count / elapsed:.2f}")
    return descriptor_count, elapsed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeat", type=int, default=100)
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument(
        "--panel",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "tests" / "smiles_panel.smi",
    )
    args = parser.parse_args()

    molecules = _load_molecules(args.panel)
    print(f"molecules: {len(molecules)}")
    print(f"repeats: {args.repeat}")
    print(f"rounds: {args.rounds}")

    timings: dict[str, list[float]] = {"current": [], "legacy": []}
    for round_index in range(args.rounds):
        print(f"round: {round_index + 1}")
        _, current_elapsed = _run_mode("current", molecules, args.repeat, None)
        _, legacy_elapsed = _run_mode(
            "legacy",
            molecules,
            args.repeat,
            _legacy_compute_path_count,
        )
        timings["current"].append(current_elapsed)
        timings["legacy"].append(legacy_elapsed)

    current_median = statistics.median(timings["current"])
    legacy_median = statistics.median(timings["legacy"])
    improvement = (legacy_median - current_median) / legacy_median * 100
    print(f"current median seconds: {current_median:.6f}")
    print(f"legacy median seconds: {legacy_median:.6f}")
    print(f"median improvement percent: {improvement:.2f}")


if __name__ == "__main__":
    main()
