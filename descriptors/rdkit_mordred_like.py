"""Mordred-compatible 2D descriptors computed with RDKit only."""

from __future__ import annotations

from collections.abc import Callable
from functools import cached_property
import math
import re

from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

from .mordred_rdkit_registry import (
    PATH_COUNT_DESCRIPTORS,
    RING_COUNT_DESCRIPTORS,
    SUPPORTED_MORDRED_2D_DESCRIPTORS,
    WALK_COUNT_DESCRIPTORS,
)

DescriptorValue = float | int
DescriptorFunction = Callable[["_DescriptorContext"], DescriptorValue]

_RING_COUNT_PATTERN = re.compile(r"^n(?:(G12|\d+))?(F)?([aA])?(H)?Ring$")
_PATH_COUNT_PATTERN = re.compile(r"^(T)?(?:(pi)PC|MPC)(\d+)$")
_WALK_COUNT_PATTERN = re.compile(r"^(T)?(?:(M)WC|(SR)W)(\d+)$")
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


class _DescriptorContext:
    """Per-molecule cache for descriptor calculations."""

    def __init__(self, mol: Chem.Mol) -> None:
        self.mol = mol
        self._adjacency_power_cache = {}
        self._path_count_cache = {}

    @cached_property
    def atoms(self) -> tuple[rdchem.Atom, ...]:
        return tuple(self.mol.GetAtoms())

    @cached_property
    def bonds(self) -> tuple[rdchem.Bond, ...]:
        return tuple(self.mol.GetBonds())

    @cached_property
    def bond_atom_pairs(self) -> tuple[tuple[int, int], ...]:
        return tuple(
            (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in self.bonds
        )

    @cached_property
    def bond_orders(self) -> tuple[float, ...]:
        return tuple(bond.GetBondTypeAsDouble() for bond in self.bonds)

    @cached_property
    def distance_matrix(self):
        return self._compute_distance_matrix()

    def _compute_distance_matrix(self):
        return Chem.GetDistanceMatrix(self.mol, force=True)

    @cached_property
    def adjacency_matrix(self):
        return self._compute_adjacency_matrix()

    def _compute_adjacency_matrix(self):
        return Chem.GetAdjacencyMatrix(self.mol, useBO=False, force=True)

    @cached_property
    def adjacency_valences(self) -> tuple[float, ...]:
        return tuple(float(value) for value in self.adjacency_matrix.sum(axis=0))

    def adjacency_power(self, order: int):
        if order not in self._adjacency_power_cache:
            if order == 1:
                self._adjacency_power_cache[order] = self.adjacency_matrix
            else:
                self._adjacency_power_cache[order] = self.adjacency_power(
                    order - 1
                ).dot(self.adjacency_matrix)
        return self._adjacency_power_cache[order]

    def path_count(self, order: int) -> tuple[int, float]:
        if order not in self._path_count_cache:
            self._path_count_cache[order] = self._compute_path_count(order)
        return self._path_count_cache[order]

    def _path_bond_ids_to_atom_ids_and_pi_weight(
        self,
        path,
    ) -> tuple[tuple[int, ...], float]:
        path_iter = iter(path)
        pi_weight = 1.0

        try:
            bond_index = next(path_iter)
        except StopIteration:
            return (), pi_weight

        pi_weight *= self.bond_orders[bond_index]
        atom0_from, atom0_to = self.bond_atom_pairs[bond_index]

        try:
            bond_index = next(path_iter)
        except StopIteration:
            return (atom0_from, atom0_to), pi_weight

        pi_weight *= self.bond_orders[bond_index]
        atom1_from, atom1_to = self.bond_atom_pairs[bond_index]

        if atom0_from in [atom1_from, atom1_to]:
            atoms = [atom0_to, atom0_from]
            current = atom1_from if atom0_from == atom1_to else atom1_to
        else:
            atoms = [atom0_from, atom0_to]
            current = atom1_from if atom0_to == atom1_to else atom1_to

        for bond_index in path_iter:
            atom_from, atom_to = self.bond_atom_pairs[bond_index]
            pi_weight *= self.bond_orders[bond_index]
            atoms.append(current)

            if atom_from == current:
                current = atom_to
            else:
                current = atom_from

        atoms.append(current)
        return tuple(atoms), pi_weight

    def _compute_path_count(self, order: int) -> tuple[int, float]:
        path_count = 0
        pi_path_count = 0.0

        for path in Chem.FindAllPathsOfLengthN(self.mol, order):
            atom_ids = set()

            atom_path, pi_weight = self._path_bond_ids_to_atom_ids_and_pi_weight(path)
            for atom_index in atom_path:
                if atom_index in atom_ids:
                    break

                atom_ids.add(atom_index)
            else:
                path_count += 1
                pi_path_count += pi_weight

        return path_count, pi_path_count

    @cached_property
    def ring_atom_sets(self) -> tuple[set[int], ...]:
        return self._compute_ring_atom_sets()

    def _compute_ring_atom_sets(self) -> tuple[set[int], ...]:
        return tuple(set(ring) for ring in self.mol.GetRingInfo().AtomRings())

    @cached_property
    def fused_ring_systems(self) -> tuple[set[int], ...]:
        return self._compute_fused_ring_systems()

    def _compute_fused_ring_systems(self) -> tuple[set[int], ...]:
        rings = self.ring_atom_sets
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

        return tuple(
            atoms
            for root, atoms in component_atoms.items()
            if component_ring_counts[root] >= 2
        )

    @cached_property
    def implicit_hydrogen_count(self) -> int:
        return self._compute_implicit_hydrogen_count()

    def _compute_implicit_hydrogen_count(self) -> int:
        return sum(atom.GetTotalNumHs() for atom in self.atoms)

    @cached_property
    def kekulized_mol(self) -> Chem.Mol:
        return self._compute_kekulized_mol()

    def _compute_kekulized_mol(self) -> Chem.Mol:
        kekule_mol = Chem.Mol(self.mol)
        Chem.Kekulize(kekule_mol, clearAromaticFlags=True)
        return kekule_mol

    @cached_property
    def exact_molecular_weight(self) -> float:
        return Descriptors.ExactMolWt(self.mol)

    @cached_property
    def total_atom_count_including_hydrogen(self) -> int:
        return self.mol.GetNumAtoms() + self.implicit_hydrogen_count

    @cached_property
    def atom_symbol_counts(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for atom in self.atoms:
            symbol = atom.GetSymbol()
            counts[symbol] = counts.get(symbol, 0) + 1
        return counts

    @cached_property
    def aromatic_atom_count(self) -> int:
        return sum(1 for atom in self.atoms if atom.GetIsAromatic())

    @cached_property
    def halogen_count(self) -> int:
        return sum(
            1
            for atom in self.atoms
            if atom.GetAtomicNum() in _HALOGEN_ATOMIC_NUMBERS
        )


def _average_molecular_weight(ctx: _DescriptorContext) -> float:
    atom_count = ctx.total_atom_count_including_hydrogen
    if atom_count == 0:
        return float("nan")
    return ctx.exact_molecular_weight / atom_count


def _diameter(ctx: _DescriptorContext) -> int:
    if ctx.distance_matrix.size == 0:
        return 0
    return int(ctx.distance_matrix.max())


def _radius(ctx: _DescriptorContext) -> int:
    if ctx.distance_matrix.size == 0:
        return 0
    return int(ctx.distance_matrix.max(axis=0).min())


def _topological_shape_index(ctx: _DescriptorContext) -> float:
    radius = _radius(ctx)
    return (_diameter(ctx) - radius) / radius


def _petitjean_index(ctx: _DescriptorContext) -> float:
    diameter = _diameter(ctx)
    return (_diameter(ctx) - _radius(ctx)) / diameter


def _wiener_path_index(ctx: _DescriptorContext) -> int:
    return int(0.5 * ctx.distance_matrix.sum())


def _wiener_polarity_index(ctx: _DescriptorContext) -> int:
    return int(0.5 * (ctx.distance_matrix == 3).sum())


def _zagreb_index_1(ctx: _DescriptorContext) -> float:
    return sum(valence**2 for valence in ctx.adjacency_valences)


def _zagreb_index_2(ctx: _DescriptorContext) -> float:
    valences = ctx.adjacency_valences
    return float(
        sum(
            valences[bond.GetBeginAtomIdx()] * valences[bond.GetEndAtomIdx()]
            for bond in ctx.bonds
        )
    )


def _modified_zagreb_index_1(ctx: _DescriptorContext) -> float:
    return sum(valence**-2 for valence in ctx.adjacency_valences)


def _modified_zagreb_index_2(ctx: _DescriptorContext) -> float:
    valences = ctx.adjacency_valences
    return float(
        sum(
            (valences[bond.GetBeginAtomIdx()] * valences[bond.GetEndAtomIdx()]) ** -1
            for bond in ctx.bonds
        )
    )


def _hydrogen_atom_count(ctx: _DescriptorContext) -> int:
    return sum(
        1 if atom.GetAtomicNum() == 1 else atom.GetTotalNumHs()
        for atom in ctx.atoms
    )


def _atom_count_by_symbol(symbol: str) -> DescriptorFunction:
    return lambda ctx: ctx.atom_symbol_counts.get(symbol, 0)


def _halogen_atom_count(ctx: _DescriptorContext) -> int:
    return ctx.halogen_count


def _aromatic_atom_count(ctx: _DescriptorContext) -> int:
    return ctx.aromatic_atom_count


def _is_single_bond(bond: rdchem.Bond) -> bool:
    return bond.GetBondType() == rdchem.BondType.SINGLE


def _is_double_bond(bond: rdchem.Bond) -> bool:
    return bond.GetBondType() == rdchem.BondType.DOUBLE


def _is_triple_bond(bond: rdchem.Bond) -> bool:
    return bond.GetBondType() == rdchem.BondType.TRIPLE


def _is_aromatic_bond(bond: rdchem.Bond) -> bool:
    return bond.GetIsAromatic()


def _implicit_hydrogen_bond_count(ctx: _DescriptorContext) -> int:
    return ctx.implicit_hydrogen_count


def _bond_count(ctx: _DescriptorContext) -> int:
    return len(ctx.bonds) + ctx.implicit_hydrogen_count


def _bond_count_by_predicate(
    ctx: _DescriptorContext,
    predicate: Callable[[rdchem.Bond], bool],
    *,
    include_implicit_hydrogen_bonds: bool = False,
) -> int:
    count = sum(1 for bond in ctx.bonds if predicate(bond))
    if include_implicit_hydrogen_bonds:
        count += ctx.implicit_hydrogen_count
    return count


def _multiple_bond_count(ctx: _DescriptorContext) -> int:
    return _bond_count_by_predicate(
        ctx,
        lambda bond: _is_double_bond(bond)
        or _is_triple_bond(bond)
        or _is_aromatic_bond(bond),
    )


def _kekulized_bond_count(
    ctx: _DescriptorContext,
    bond_type: rdchem.BondType,
    *,
    include_implicit_hydrogen_bonds: bool = False,
) -> int:
    count = sum(
        1 for bond in ctx.kekulized_mol.GetBonds() if bond.GetBondType() == bond_type
    )
    if include_implicit_hydrogen_bonds:
        count += ctx.implicit_hydrogen_count
    return count


def _ring_matches_filters(
    ctx: _DescriptorContext,
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
        is_aromatic = all(ctx.mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in atoms)
        if aromaticity == "a" and not is_aromatic:
            return False
        if aromaticity == "A" and is_aromatic:
            return False

    if hetero is not None and not any(
        ctx.mol.GetAtomWithIdx(atom).GetAtomicNum() != 6 for atom in atoms
    ):
        return False

    return True


def _ring_count_descriptor(name: str) -> DescriptorFunction:
    match = _RING_COUNT_PATTERN.match(name)
    if match is None:
        msg = f"unsupported ring count descriptor: {name}"
        raise ValueError(msg)

    size, fused, aromaticity, hetero = match.groups()

    def calc(ctx: _DescriptorContext) -> int:
        if fused is not None:
            atom_sets = ctx.fused_ring_systems
        else:
            atom_sets = ctx.ring_atom_sets

        return sum(
            1
            for atoms in atom_sets
            if _ring_matches_filters(ctx, atoms, size, aromaticity, hetero)
        )

    return calc


def _path_count_descriptor(name: str) -> DescriptorFunction:
    match = _PATH_COUNT_PATTERN.match(name)
    if match is None:
        msg = f"unsupported path count descriptor: {name}"
        raise ValueError(msg)

    total_prefix, pi_prefix, order_text = match.groups()
    order = int(order_text)
    use_pi = pi_prefix is not None
    total = total_prefix is not None

    def raw_value_for_order(ctx: _DescriptorContext, path_order: int) -> int | float:
        if path_order == 0:
            return ctx.mol.GetNumAtoms()
        path_count, pi_path_count = ctx.path_count(path_order)
        return pi_path_count if use_pi else path_count

    def value_for_order(ctx: _DescriptorContext, path_order: int) -> int | float:
        value = raw_value_for_order(ctx, path_order)
        if use_pi:
            return math.log(value + 1)
        return value

    def calc(ctx: _DescriptorContext) -> int | float:
        if total:
            value = sum(
                raw_value_for_order(ctx, path_order)
                for path_order in range(0, order + 1)
            )
            if use_pi:
                return math.log(value + 1)
            return value
        return value_for_order(ctx, order)

    return calc


def _walk_count_descriptor(name: str) -> DescriptorFunction:
    match = _WALK_COUNT_PATTERN.match(name)
    if match is None:
        msg = f"unsupported walk count descriptor: {name}"
        raise ValueError(msg)

    total_prefix, molecular_walk, self_returning_walk, order_text = match.groups()
    order = int(order_text)
    total = total_prefix is not None
    self_returning = self_returning_walk is not None

    def value_for_order(ctx: _DescriptorContext, walk_order: int) -> float:
        matrix_power = ctx.adjacency_power(walk_order)
        if self_returning:
            return math.log(matrix_power.trace() + 1)
        if walk_order == 1:
            return 0.5 * matrix_power.sum()
        return math.log(matrix_power.sum() + 1)

    def calc(ctx: _DescriptorContext) -> float:
        if total:
            values = [float(ctx.mol.GetNumAtoms())]
            start = 2 if self_returning else 1
            values.extend(
                value_for_order(ctx, walk_order)
                for walk_order in range(start, order + 1)
            )
            return sum(values)
        return value_for_order(ctx, order)

    return calc


def _rdkit_descriptor(
    function: Callable[[Chem.Mol], DescriptorValue],
) -> DescriptorFunction:
    return lambda ctx: function(ctx.mol)


_DESCRIPTOR_FUNCTIONS: dict[str, DescriptorFunction] = {
    "AMW": _average_molecular_weight,
    "BertzCT": _rdkit_descriptor(Descriptors.BertzCT),
    "Diameter": _diameter,
    "FCSP3": _rdkit_descriptor(rdMolDescriptors.CalcFractionCSP3),
    "MW": lambda ctx: ctx.exact_molecular_weight,
    "NaaNH": _rdkit_descriptor(Descriptors.fr_Ar_NH),
    "NaaO": _rdkit_descriptor(Descriptors.fr_furan),
    "NddsN": _rdkit_descriptor(Descriptors.fr_nitro),
    "NsNH2": _rdkit_descriptor(Descriptors.fr_NH2),
    "NsSH": _rdkit_descriptor(Descriptors.fr_SH),
    "NssO": _rdkit_descriptor(Descriptors.fr_ether),
    "NssS": _rdkit_descriptor(Descriptors.fr_sulfide),
    "NtN": _rdkit_descriptor(Descriptors.fr_nitrile),
    "NtsC": _rdkit_descriptor(Descriptors.fr_nitrile),
    "PetitjeanIndex": _petitjean_index,
    "Radius": _radius,
    "SMR": _rdkit_descriptor(Crippen.MolMR),
    "SLogP": _rdkit_descriptor(Crippen.MolLogP),
    "TopoPSA": lambda ctx: rdMolDescriptors.CalcTPSA(
        ctx.mol,
        includeSandP=True,
    ),
    "TopoPSA(NO)": _rdkit_descriptor(rdMolDescriptors.CalcTPSA),
    "TopoShapeIndex": _topological_shape_index,
    "WPath": _wiener_path_index,
    "WPol": _wiener_polarity_index,
    "Zagreb1": _zagreb_index_1,
    "Zagreb2": _zagreb_index_2,
    "mZagreb1": _modified_zagreb_index_1,
    "mZagreb2": _modified_zagreb_index_2,
    "nAromAtom": _aromatic_atom_count,
    "nAromBond": lambda ctx: _bond_count_by_predicate(ctx, _is_aromatic_bond),
    "nAtom": lambda ctx: ctx.total_atom_count_including_hydrogen,
    "nBridgehead": _rdkit_descriptor(rdMolDescriptors.CalcNumBridgeheadAtoms),
    "nBonds": _bond_count,
    "nBondsA": lambda ctx: _bond_count_by_predicate(ctx, _is_aromatic_bond),
    "nBondsD": lambda ctx: _bond_count_by_predicate(ctx, _is_double_bond),
    "nBondsKD": lambda ctx: _kekulized_bond_count(ctx, rdchem.BondType.DOUBLE),
    "nBondsKS": lambda ctx: _kekulized_bond_count(
        ctx,
        rdchem.BondType.SINGLE,
        include_implicit_hydrogen_bonds=True,
    ),
    "nBondsM": _multiple_bond_count,
    "nBondsO": lambda ctx: len(ctx.bonds),
    "nBondsS": lambda ctx: _bond_count_by_predicate(
        ctx,
        _is_single_bond,
        include_implicit_hydrogen_bonds=True,
    ),
    "nBondsT": lambda ctx: _bond_count_by_predicate(ctx, _is_triple_bond),
    "nHBAcc": _rdkit_descriptor(rdMolDescriptors.CalcNumHBA),
    "nHBDon": _rdkit_descriptor(rdMolDescriptors.CalcNumHBD),
    "nH": _hydrogen_atom_count,
    "nHeavyAtom": _rdkit_descriptor(Descriptors.HeavyAtomCount),
    "nHetero": _rdkit_descriptor(rdMolDescriptors.CalcNumHeteroatoms),
    "nRot": _rdkit_descriptor(rdMolDescriptors.CalcNumRotatableBonds),
    "nSpiro": _rdkit_descriptor(rdMolDescriptors.CalcNumSpiroAtoms),
    "nX": _halogen_atom_count,
}

for _name, _symbol in _ATOM_SYMBOLS_BY_DESCRIPTOR.items():
    _DESCRIPTOR_FUNCTIONS[_name] = _atom_count_by_symbol(_symbol)

for _name in RING_COUNT_DESCRIPTORS:
    _DESCRIPTOR_FUNCTIONS[_name] = _ring_count_descriptor(_name)

for _name in PATH_COUNT_DESCRIPTORS:
    _DESCRIPTOR_FUNCTIONS[_name] = _path_count_descriptor(_name)

for _name in WALK_COUNT_DESCRIPTORS:
    _DESCRIPTOR_FUNCTIONS[_name] = _walk_count_descriptor(_name)

for _name in SUPPORTED_MORDRED_2D_DESCRIPTORS:
    if _name not in _DESCRIPTOR_FUNCTIONS and hasattr(Descriptors, _name):
        _DESCRIPTOR_FUNCTIONS[_name] = _rdkit_descriptor(getattr(Descriptors, _name))


def calc_rdkit_mordred_like_2d(mol: Chem.Mol) -> dict[str, DescriptorValue]:
    """Return Mordred-compatible 2D descriptors computed using RDKit only."""

    if mol is None:
        msg = "mol must be an RDKit Mol, not None"
        raise ValueError(msg)

    context = _DescriptorContext(mol)
    return {
        name: _DESCRIPTOR_FUNCTIONS[name](context)
        for name in SUPPORTED_MORDRED_2D_DESCRIPTORS
    }
