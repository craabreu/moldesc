"""Mordred-compatible 2D descriptors computed with RDKit only."""

from __future__ import annotations

from collections.abc import Callable
import re

from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

from .mordred_rdkit_registry import (
    RING_COUNT_DESCRIPTORS,
    SUPPORTED_MORDRED_2D_DESCRIPTORS,
)

DescriptorValue = float | int
DescriptorFunction = Callable[[Chem.Mol], DescriptorValue]

_RING_COUNT_PATTERN = re.compile(r"^n(?:(G12|\d+))?(F)?([aA])?(H)?Ring$")
_HALOGEN_ATOMIC_NUMBERS = {9, 17, 35, 53}
_ATOM_SYMBOLS_BY_DESCRIPTOR = {
    "nB": "B",
    "nC": "C",
    "nN": "N",
    "nO": "O",
    "nS": "S",
    "nP": "P",
    "nF": "F",
    "nCl": "Cl",
    "nBr": "Br",
    "nI": "I",
}


def _total_atom_count_including_hydrogen(mol: Chem.Mol) -> int:
    return mol.GetNumAtoms() + sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())


def _average_molecular_weight(mol: Chem.Mol) -> float:
    atom_count = _total_atom_count_including_hydrogen(mol)
    if atom_count == 0:
        return float("nan")
    return Descriptors.ExactMolWt(mol) / atom_count


def _distance_matrix(mol: Chem.Mol):
    return Chem.GetDistanceMatrix(mol, force=True)


def _adjacency_valences(mol: Chem.Mol) -> list[float]:
    matrix = Chem.GetAdjacencyMatrix(mol, useBO=False, force=True)
    return [float(value) for value in matrix.sum(axis=0)]


def _diameter(mol: Chem.Mol) -> int:
    matrix = _distance_matrix(mol)
    if matrix.size == 0:
        return 0
    return int(matrix.max())


def _radius(mol: Chem.Mol) -> int:
    matrix = _distance_matrix(mol)
    if matrix.size == 0:
        return 0
    return int(matrix.max(axis=0).min())


def _topological_shape_index(mol: Chem.Mol) -> float:
    radius = _radius(mol)
    return (_diameter(mol) - radius) / radius


def _petitjean_index(mol: Chem.Mol) -> float:
    diameter = _diameter(mol)
    return (_diameter(mol) - _radius(mol)) / diameter


def _wiener_path_index(mol: Chem.Mol) -> int:
    return int(0.5 * _distance_matrix(mol).sum())


def _wiener_polarity_index(mol: Chem.Mol) -> int:
    return int(0.5 * (_distance_matrix(mol) == 3).sum())


def _zagreb_index_1(mol: Chem.Mol) -> float:
    return sum(valence**2 for valence in _adjacency_valences(mol))


def _zagreb_index_2(mol: Chem.Mol) -> float:
    valences = _adjacency_valences(mol)
    return float(
        sum(
            valences[bond.GetBeginAtomIdx()] * valences[bond.GetEndAtomIdx()]
            for bond in mol.GetBonds()
        )
    )


def _modified_zagreb_index_1(mol: Chem.Mol) -> float:
    return sum(valence**-2 for valence in _adjacency_valences(mol))


def _modified_zagreb_index_2(mol: Chem.Mol) -> float:
    valences = _adjacency_valences(mol)
    return float(
        sum(
            (valences[bond.GetBeginAtomIdx()] * valences[bond.GetEndAtomIdx()]) ** -1
            for bond in mol.GetBonds()
        )
    )


def _hydrogen_atom_count(mol: Chem.Mol) -> int:
    return sum(
        1 if atom.GetAtomicNum() == 1 else atom.GetTotalNumHs()
        for atom in mol.GetAtoms()
    )


def _atom_count_by_symbol(symbol: str) -> DescriptorFunction:
    return lambda mol: sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == symbol)


def _halogen_atom_count(mol: Chem.Mol) -> int:
    return sum(
        1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in _HALOGEN_ATOMIC_NUMBERS
    )


def _aromatic_atom_count(mol: Chem.Mol) -> int:
    return sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())


def _is_single_bond(bond: rdchem.Bond) -> bool:
    return bond.GetBondType() == rdchem.BondType.SINGLE


def _is_double_bond(bond: rdchem.Bond) -> bool:
    return bond.GetBondType() == rdchem.BondType.DOUBLE


def _is_triple_bond(bond: rdchem.Bond) -> bool:
    return bond.GetBondType() == rdchem.BondType.TRIPLE


def _is_aromatic_bond(bond: rdchem.Bond) -> bool:
    return bond.GetIsAromatic()


def _implicit_hydrogen_bond_count(mol: Chem.Mol) -> int:
    return sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())


def _bond_count(mol: Chem.Mol) -> int:
    return mol.GetNumBonds() + _implicit_hydrogen_bond_count(mol)


def _bond_count_by_predicate(
    mol: Chem.Mol,
    predicate: Callable[[rdchem.Bond], bool],
    *,
    include_implicit_hydrogen_bonds: bool = False,
) -> int:
    count = sum(1 for bond in mol.GetBonds() if predicate(bond))
    if include_implicit_hydrogen_bonds:
        count += _implicit_hydrogen_bond_count(mol)
    return count


def _multiple_bond_count(mol: Chem.Mol) -> int:
    return _bond_count_by_predicate(
        mol,
        lambda bond: _is_double_bond(bond)
        or _is_triple_bond(bond)
        or _is_aromatic_bond(bond),
    )


def _kekulized_bond_count(
    mol: Chem.Mol,
    bond_type: rdchem.BondType,
    *,
    include_implicit_hydrogen_bonds: bool = False,
) -> int:
    kekule_mol = Chem.Mol(mol)
    Chem.Kekulize(kekule_mol, clearAromaticFlags=True)
    count = sum(1 for bond in kekule_mol.GetBonds() if bond.GetBondType() == bond_type)
    if include_implicit_hydrogen_bonds:
        count += _implicit_hydrogen_bond_count(mol)
    return count


def _ring_matches_filters(
    mol: Chem.Mol,
    atoms: set[int],
    size: str | None,
    aromaticity: str | None,
    hetero: str | None,
) -> bool:
    if size == "G12":
        if len(atoms) < 12:
            return False
    elif size is not None and len(atoms) != int(size):
        return False

    if aromaticity is not None:
        is_aromatic = all(mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in atoms)
        if aromaticity == "a" and not is_aromatic:
            return False
        if aromaticity == "A" and is_aromatic:
            return False

    if hetero is not None and not any(
        mol.GetAtomWithIdx(atom).GetAtomicNum() != 6 for atom in atoms
    ):
        return False

    return True


def _fused_ring_systems(mol: Chem.Mol) -> list[set[int]]:
    rings = [set(ring) for ring in mol.GetRingInfo().AtomRings()]
    parent = list(range(len(rings)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    for left, left_ring in enumerate(rings):
        for right in range(left + 1, len(rings)):
            if len(left_ring & rings[right]) >= 2:
                union(left, right)

    component_atoms: dict[int, set[int]] = {}
    component_ring_counts: dict[int, int] = {}
    for index, ring in enumerate(rings):
        root = find(index)
        component_atoms.setdefault(root, set()).update(ring)
        component_ring_counts[root] = component_ring_counts.get(root, 0) + 1

    return [
        atoms
        for root, atoms in component_atoms.items()
        if component_ring_counts[root] >= 2
    ]


def _ring_count_descriptor(name: str) -> DescriptorFunction:
    match = _RING_COUNT_PATTERN.match(name)
    if match is None:
        msg = f"unsupported ring count descriptor: {name}"
        raise ValueError(msg)

    size, fused, aromaticity, hetero = match.groups()

    def calc(mol: Chem.Mol) -> int:
        if fused is not None:
            atom_sets = _fused_ring_systems(mol)
        else:
            atom_sets = [set(ring) for ring in mol.GetRingInfo().AtomRings()]

        return sum(
            1
            for atoms in atom_sets
            if _ring_matches_filters(mol, atoms, size, aromaticity, hetero)
        )

    return calc


_DESCRIPTOR_FUNCTIONS: dict[str, DescriptorFunction] = {
    "AMW": _average_molecular_weight,
    "BertzCT": Descriptors.BertzCT,
    "Diameter": _diameter,
    "FCSP3": rdMolDescriptors.CalcFractionCSP3,
    "MW": Descriptors.ExactMolWt,
    "PetitjeanIndex": _petitjean_index,
    "Radius": _radius,
    "SMR": Crippen.MolMR,
    "SLogP": Crippen.MolLogP,
    "TopoPSA": lambda mol: rdMolDescriptors.CalcTPSA(mol, includeSandP=True),
    "TopoPSA(NO)": rdMolDescriptors.CalcTPSA,
    "TopoShapeIndex": _topological_shape_index,
    "WPath": _wiener_path_index,
    "WPol": _wiener_polarity_index,
    "Zagreb1": _zagreb_index_1,
    "Zagreb2": _zagreb_index_2,
    "mZagreb1": _modified_zagreb_index_1,
    "mZagreb2": _modified_zagreb_index_2,
    "nAromAtom": _aromatic_atom_count,
    "nAromBond": lambda mol: _bond_count_by_predicate(mol, _is_aromatic_bond),
    "nAtom": _total_atom_count_including_hydrogen,
    "nBridgehead": rdMolDescriptors.CalcNumBridgeheadAtoms,
    "nBonds": _bond_count,
    "nBondsA": lambda mol: _bond_count_by_predicate(mol, _is_aromatic_bond),
    "nBondsD": lambda mol: _bond_count_by_predicate(mol, _is_double_bond),
    "nBondsKD": lambda mol: _kekulized_bond_count(mol, rdchem.BondType.DOUBLE),
    "nBondsKS": lambda mol: _kekulized_bond_count(
        mol,
        rdchem.BondType.SINGLE,
        include_implicit_hydrogen_bonds=True,
    ),
    "nBondsM": _multiple_bond_count,
    "nBondsO": Chem.Mol.GetNumBonds,
    "nBondsS": lambda mol: _bond_count_by_predicate(
        mol,
        _is_single_bond,
        include_implicit_hydrogen_bonds=True,
    ),
    "nBondsT": lambda mol: _bond_count_by_predicate(mol, _is_triple_bond),
    "nHBAcc": rdMolDescriptors.CalcNumHBA,
    "nHBDon": rdMolDescriptors.CalcNumHBD,
    "nH": _hydrogen_atom_count,
    "nHeavyAtom": Descriptors.HeavyAtomCount,
    "nHetero": rdMolDescriptors.CalcNumHeteroatoms,
    "nRot": rdMolDescriptors.CalcNumRotatableBonds,
    "nSpiro": rdMolDescriptors.CalcNumSpiroAtoms,
    "nX": _halogen_atom_count,
}

for _name, _symbol in _ATOM_SYMBOLS_BY_DESCRIPTOR.items():
    _DESCRIPTOR_FUNCTIONS[_name] = _atom_count_by_symbol(_symbol)

for _name in RING_COUNT_DESCRIPTORS:
    _DESCRIPTOR_FUNCTIONS[_name] = _ring_count_descriptor(_name)

for _name in SUPPORTED_MORDRED_2D_DESCRIPTORS:
    if _name not in _DESCRIPTOR_FUNCTIONS and hasattr(Descriptors, _name):
        _DESCRIPTOR_FUNCTIONS[_name] = getattr(Descriptors, _name)


def calc_rdkit_mordred_like_2d(mol: Chem.Mol) -> dict[str, DescriptorValue]:
    """Return Mordred-compatible 2D descriptors computed using RDKit only."""

    if mol is None:
        msg = "mol must be an RDKit Mol, not None"
        raise ValueError(msg)

    return {name: _DESCRIPTOR_FUNCTIONS[name](mol) for name in SUPPORTED_MORDRED_2D_DESCRIPTORS}
