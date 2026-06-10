"""Mordred-name 2D descriptors computed with RDKit only."""

from __future__ import annotations

from collections.abc import Callable
from functools import cached_property
import math
import re
from typing import NamedTuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from rdkit.Chem.EState import AtomTypes

from .mordred_rdkit_registry import (
    AUTOCORRELATION_M_DESCRIPTORS,
    AUTOCORRELATION_Z_DESCRIPTORS,
    BCUT_Z_DESCRIPTORS,
    ESTATE_ATOM_TYPE_DESCRIPTORS,
    PATH_COUNT_DESCRIPTORS,
    PHYSICAL_PROPERTY_DESCRIPTORS,
    RING_COUNT_DESCRIPTORS,
    SMALL_GRAPH_FORMULA_DESCRIPTORS,
    SUPPORTED_MORDRED_2D_DESCRIPTORS,
    WALK_COUNT_DESCRIPTORS,
)

DescriptorValue = float | int
DescriptorFunction = Callable[["_DescriptorContext"], DescriptorValue]


class _RingMetadata(NamedTuple):
    atoms: frozenset[int]
    size: int
    is_aromatic: bool
    has_hetero: bool

_RING_COUNT_PATTERN = re.compile(r"^n(?:(G12|\d+))?(F)?([aA])?(H)?Ring$")
_PATH_COUNT_PATTERN = re.compile(r"^(T)?(?:(pi)PC|MPC)(\d+)$")
_WALK_COUNT_PATTERN = re.compile(r"^(T)?(?:(M)WC|(SR)W)(\d+)$")
_AUTOCORRELATION_PATTERN = re.compile(
    r"^(AATSC|AATS|ATSC|ATS|MATS|GATS)(\d+)(Z|m)$"
)
_HALOGEN_ATOMIC_NUMBERS = {9, 17, 35, 53}
_ACID_GROUP_SMARTS = (
    "[O;H1]-[C,S,P]=O",
    "[*;-;!$(*~[*;+])]",
    "[NH](S(=O)=O)C(F)(F)F",
    "n1nnnc1",
)
_BASE_GROUP_SMARTS = (
    "[NH2]-[CX4]",
    "[NH](-[CX4])-[CX4]",
    "N(-[CX4])(-[CX4])-[CX4]",
    "[*;+;!$(*~[*;-])]",
    "N=C-N",
    "N-C=N",
)
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
_POLARIZABILITY_94_BY_ATOMIC_NUM = {
    1: 0.666793,
    2: 0.2050522,
    3: 24.33,
    4: 5.6,
    5: 3.03,
    6: 1.67,
    7: 1.1,
    8: 0.802,
    9: 0.557,
    10: 0.39432,
    11: 24.11,
    12: 10.6,
    13: 6.8,
    14: 5.53,
    15: 3.63,
    16: 2.9,
    17: 2.18,
    18: 1.6411,
    19: 43.06,
    20: 22.8,
    21: 17.8,
    22: 14.6,
    23: 12.4,
    24: 11.6,
    25: 9.4,
    26: 8.4,
    27: 7.5,
    28: 6.8,
    29: 6.2,
    30: 5.75,
    31: 8.12,
    32: 5.84,
    33: 4.31,
    34: 3.77,
    35: 3.05,
    36: 2.4844,
    37: 47.24,
    38: 23.5,
    39: 22.7,
    40: 17.9,
    41: 15.7,
    42: 12.8,
    43: 11.4,
    44: 9.6,
    45: 8.6,
    46: 4.8,
    47: 6.78,
    48: 7.36,
    49: 10.2,
    50: 7.84,
    51: 6.6,
    52: 5.5,
    53: 5.35,
    54: 4.044,
    55: 59.42,
    56: 39.7,
    57: 31.1,
    58: 29.6,
    59: 28.2,
    60: 31.4,
    61: 30.1,
    62: 28.8,
    63: 27.7,
    64: 23.5,
    65: 25.5,
    66: 24.5,
    67: 23.6,
    68: 22.7,
    69: 21.8,
    70: 20.9,
    71: 21.9,
    72: 16.2,
    73: 13.1,
    74: 11.1,
    75: 9.7,
    76: 8.5,
    77: 7.6,
    78: 6.5,
    79: 5.8,
    80: 5.02,
    81: 7.6,
    82: 7.01,
    83: 7.4,
    84: 6.8,
    85: 6.0,
    86: 5.3,
    87: 48.6,
    88: 38.3,
    89: 32.1,
    90: 32.1,
    91: 25.4,
    92: 24.9,
    93: 24.8,
    94: 24.5,
    95: 23.3,
    96: 23.0,
    97: 22.7,
    98: 20.5,
    99: 19.7,
    100: 23.8,
    101: 18.2,
    102: 16.4,
    112: 4.06,
    114: 4.59,
}
_MCGOWAN_VOLUME_BY_ATOMIC_NUM = {
    1: 8.71,
    2: 6.75,
    3: 22.23,
    4: 20.27,
    5: 18.31,
    6: 16.35,
    7: 14.39,
    8: 12.43,
    9: 10.47,
    10: 8.51,
    11: 32.71,
    12: 30.75,
    13: 28.79,
    14: 26.83,
    15: 24.87,
    16: 22.91,
    17: 20.95,
    18: 18.99,
    19: 51.89,
    20: 50.28,
    21: 48.68,
    22: 47.07,
    23: 45.47,
    24: 43.86,
    25: 42.26,
    26: 40.65,
    27: 39.05,
    28: 37.44,
    29: 35.84,
    30: 34.23,
    31: 32.63,
    32: 31.02,
    33: 29.42,
    34: 27.81,
    35: 26.21,
    36: 24.6,
    37: 60.22,
    38: 58.61,
    39: 57.01,
    40: 55.4,
    41: 53.8,
    42: 52.19,
    43: 50.59,
    44: 48.98,
    45: 47.38,
    46: 45.77,
    47: 44.17,
    48: 42.56,
    49: 40.96,
    50: 39.35,
    51: 37.75,
    52: 36.14,
    53: 34.54,
    54: 32.93,
    55: 77.25,
    56: 76.0,
    57: 74.75,
    58: 73.49,
    59: 72.24,
    60: 70.99,
    61: 69.74,
    62: 68.49,
    63: 67.23,
    64: 65.98,
    65: 64.73,
    66: 63.48,
    67: 62.23,
    68: 60.97,
    69: 59.72,
    70: 58.47,
    71: 57.22,
    72: 55.97,
    73: 54.71,
    74: 53.46,
    75: 52.21,
    76: 50.96,
    77: 49.71,
    78: 48.45,
    79: 47.2,
    80: 45.95,
    81: 44.7,
    82: 43.45,
    83: 42.19,
    84: 40.94,
    85: 39.69,
    86: 38.44,
    87: 75.59,
    88: 74.34,
    89: 73.09,
    90: 71.83,
    91: 70.58,
    92: 69.33,
    93: 68.08,
    94: 66.83,
    95: 65.57,
    96: 64.32,
    97: 63.07,
    98: 61.82,
    99: 60.57,
    100: 59.31,
    101: 58.06,
    102: 56.81,
    103: 55.56,
}
# Standard atomic weights as used by Mordred's mass property. These differ from
# RDKit's GetAtomicWeight beyond tolerance for several elements (e.g. S, Cl, B),
# so the values are replicated here rather than read from RDKit.
_MASS_BY_ATOMIC_NUM = {
    1: 1.008,
    2: 4.002602,
    3: 6.94,
    4: 9.012182,
    5: 10.81,
    6: 12.011,
    7: 14.007,
    8: 15.999,
    9: 18.9984032,
    10: 20.1797,
    11: 22.98976928,
    12: 24.305,
    13: 26.9815386,
    14: 28.085,
    15: 30.973762,
    16: 32.06,
    17: 35.45,
    18: 39.948,
    19: 39.0983,
    20: 40.078,
    21: 44.955912,
    22: 47.867,
    23: 50.9415,
    24: 51.9961,
    25: 54.938045,
    26: 55.845,
    27: 58.933195,
    28: 58.6934,
    29: 63.546,
    30: 65.38,
    31: 69.723,
    32: 72.63,
    33: 74.9216,
    34: 78.96,
    35: 79.904,
    36: 83.798,
    37: 85.4678,
    38: 87.62,
    39: 88.90585,
    40: 91.224,
    41: 92.90638,
    42: 95.96,
    43: 98.0,
    44: 101.07,
    45: 102.9055,
    46: 106.42,
    47: 107.8682,
    48: 112.411,
    49: 114.818,
    50: 118.71,
    51: 121.76,
    52: 127.6,
    53: 126.90447,
    54: 131.293,
    55: 132.9054519,
    56: 137.327,
    57: 138.90547,
    58: 140.116,
    59: 140.90765,
    60: 144.242,
    61: 145.0,
    62: 150.36,
    63: 151.964,
    64: 157.25,
    65: 158.92535,
    66: 162.5,
    67: 164.93032,
    68: 167.259,
    69: 168.93421,
    70: 173.054,
    71: 174.9668,
    72: 178.49,
    73: 180.94788,
    74: 183.84,
    75: 186.207,
    76: 190.23,
    77: 192.217,
    78: 195.084,
    79: 196.966569,
    80: 200.59,
    81: 204.38,
    82: 207.2,
    83: 208.9804,
    84: 210.0,
    85: 210.0,
    86: 222.0,
    87: 223.0,
    88: 226.0,
    89: 227.0,
    90: 232.03806,
    91: 231.03588,
    92: 238.02891,
    93: 237.0,
    94: 244.0,
    95: 243.0,
    96: 247.0,
    97: 247.0,
    98: 251.0,
    99: 252.0,
    100: 257.0,
    101: 258.0,
    102: 259.0,
    103: 262.0,
    104: 261.0,
    105: 262.0,
    106: 266.0,
    107: 264.0,
    108: 269.0,
    109: 268.0,
    110: 271.0,
}
_BONDI_RADII_BY_ATOMIC_NUM = {
    1: 1.20,
    5: 2.13,
    6: 1.70,
    7: 1.55,
    8: 1.52,
    9: 1.47,
    14: 2.10,
    15: 1.80,
    16: 1.80,
    17: 1.75,
    33: 1.85,
    34: 1.90,
    35: 1.85,
}
_VABC_ATOM_CONTRIBUTION_BY_ATOMIC_NUM = {
    atomic_num: 4.0 / 3.0 * math.pi * radius**3
    for atomic_num, radius in _BONDI_RADII_BY_ATOMIC_NUM.items()
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
    def ring_metadata(self) -> tuple[_RingMetadata, ...]:
        return self._compute_ring_metadata(self.ring_atom_sets)

    @cached_property
    def fused_ring_metadata(self) -> tuple[_RingMetadata, ...]:
        return self._compute_ring_metadata(self.fused_ring_systems)

    def _compute_ring_metadata(
        self,
        ring_sets: tuple[set[int], ...],
    ) -> tuple[_RingMetadata, ...]:
        metadata = []
        for atoms in ring_sets:
            atom_tuple = tuple(self.mol.GetAtomWithIdx(atom) for atom in atoms)
            metadata.append(
                _RingMetadata(
                    atoms=frozenset(atoms),
                    size=len(atoms),
                    is_aromatic=all(atom.GetIsAromatic() for atom in atom_tuple),
                    has_hetero=any(atom.GetAtomicNum() != 6 for atom in atom_tuple),
                )
            )
        return tuple(metadata)

    @cached_property
    def framework_linker_atoms(self) -> set[int]:
        return self._compute_framework_linker_atoms()

    def _compute_framework_linker_atoms(self) -> set[int]:
        rings = self.ring_atom_sets
        if len(rings) < 2:
            return set()

        atom_to_node: dict[int, tuple[str, int]] = {}
        ring_nodes: set[tuple[str, int]] = set()
        for ring_index, ring in enumerate(rings):
            ring_node = ("R", ring_index)
            ring_nodes.add(ring_node)
            for atom_index in ring:
                atom_to_node[atom_index] = ring_node

        graph: dict[tuple[str, int], set[tuple[str, int]]] = {}
        for bond in self.bonds:
            begin_index = bond.GetBeginAtomIdx()
            end_index = bond.GetEndAtomIdx()
            begin = atom_to_node.get(begin_index, ("A", begin_index))
            end = atom_to_node.get(end_index, ("A", end_index))
            if begin == end:
                continue
            graph.setdefault(begin, set()).add(end)
            graph.setdefault(end, set()).add(begin)

        linkers: set[int] = set()
        ring_node_list = list(ring_nodes)
        for left_index, left in enumerate(ring_node_list):
            for right in ring_node_list[left_index + 1:]:
                path = _shortest_path(graph, left, right)
                if path is not None:
                    linkers.update(
                        index for node_type, index in path if node_type == "A"
                    )
        return linkers

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
    def explicit_hydrogen_mol(self) -> Chem.Mol:
        return Chem.AddHs(self.mol)

    @cached_property
    def autocorrelation_distance_matrix(self):
        return Chem.GetDistanceMatrix(self.explicit_hydrogen_mol, force=True)

    @cached_property
    def autocorrelation_atomic_numbers(self) -> tuple[float, ...]:
        return tuple(
            float(atom.GetAtomicNum()) for atom in self.explicit_hydrogen_mol.GetAtoms()
        )

    @cached_property
    def autocorrelation_atomic_masses(self) -> tuple[float, ...]:
        return tuple(
            _atomic_property_value(_MASS_BY_ATOMIC_NUM, atom.GetAtomicNum())
            for atom in self.explicit_hydrogen_mol.GetAtoms()
        )

    @cached_property
    def autocorrelation_property_vectors(self) -> dict[str, tuple[float, ...]]:
        """Per-atom property vectors keyed by Mordred property suffix."""

        return {
            "Z": self.autocorrelation_atomic_numbers,
            "m": self.autocorrelation_atomic_masses,
        }

    @cached_property
    def autocorrelation_values(self) -> dict[str, float]:
        """Moreau-Broto/Moran/Geary autocorrelations for every property and lag.

        All property vectors are stacked and processed together so the shared
        per-lag graph-distance work (the adjacency-at-distance matrix and its
        quadratic forms) is computed once per lag instead of once per property.
        """

        vectors = self.autocorrelation_property_vectors
        suffixes = tuple(vectors)
        property_matrix = np.array([vectors[s] for s in suffixes], dtype=float)
        atom_count = property_matrix.shape[1]

        centered_matrix = property_matrix - property_matrix.mean(axis=1, keepdims=True)
        centered_square_sums = (centered_matrix**2).sum(axis=1)
        distance_matrix = np.asarray(self.autocorrelation_distance_matrix)

        results: dict[str, float] = {}
        for order in range(0, 9):
            if order == 0:
                pair_count = atom_count
                ats = (property_matrix**2).sum(axis=1)
                atsc = centered_square_sums
                degrees = None
            else:
                adjacency = (distance_matrix == order).astype(float)
                pair_count = int(adjacency.sum() // 2)
                ats = 0.5 * ((property_matrix @ adjacency) * property_matrix).sum(axis=1)
                atsc = 0.5 * ((centered_matrix @ adjacency) * centered_matrix).sum(axis=1)
                degrees = adjacency.sum(axis=1)

            for index, suffix in enumerate(suffixes):
                results[f"ATS{order}{suffix}"] = float(ats[index])
                results[f"ATSC{order}{suffix}"] = float(atsc[index])
                results[f"AATS{order}{suffix}"] = (
                    float(ats[index] / pair_count) if pair_count else float("nan")
                )
                results[f"AATSC{order}{suffix}"] = (
                    float(atsc[index] / pair_count) if pair_count else float("nan")
                )

            if order == 0:
                continue

            # Geary numerator, vectorized across properties. With B the
            # adjacency-at-distance matrix and w a property vector,
            #   sum_ij B_ij (w_i - w_j)^2 = 2 (w^2 . deg) - 2 (w^T B w)
            # and w^T B w = 2 * ATS, so the doubly-counted sum is
            #   2 (w^2 . deg) - 4 * ATS, which Mordred divides by 4 * pair_count.
            weighted_degree = (property_matrix**2) @ degrees
            geary_numerators = (
                (2.0 * weighted_degree - 4.0 * ats) / (4.0 * pair_count)
                if pair_count
                else None
            )

            for index, suffix in enumerate(suffixes):
                centered_square_sum = centered_square_sums[index]
                aatsc = atsc[index] / pair_count if pair_count else float("nan")
                results[f"MATS{order}{suffix}"] = (
                    float(atom_count * aatsc / centered_square_sum)
                    if centered_square_sum
                    else float("nan")
                )

                geary_denominator = (
                    centered_square_sum / (atom_count - 1)
                    if atom_count > 1
                    else float("nan")
                )
                if pair_count and geary_denominator and not math.isnan(geary_denominator):
                    results[f"GATS{order}{suffix}"] = float(
                        geary_numerators[index] / geary_denominator
                    )
                else:
                    results[f"GATS{order}{suffix}"] = float("nan")

        return results

    @cached_property
    def bcut_z_values(self) -> dict[str, float]:
        burden_matrix = 0.001 * np.ones((len(self.atoms), len(self.atoms)))
        for bond in self.bonds:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            begin_index = begin_atom.GetIdx()
            end_index = end_atom.GetIdx()
            weight = bond.GetBondTypeAsDouble() / 10.0
            if begin_atom.GetDegree() == 1 or end_atom.GetDegree() == 1:
                weight += 0.01

            burden_matrix[begin_index, end_index] = weight
            burden_matrix[end_index, begin_index] = weight

        for atom in self.atoms:
            burden_matrix[atom.GetIdx(), atom.GetIdx()] = atom.GetAtomicNum()

        eigenvalues = np.linalg.eig(burden_matrix)[0]
        if np.iscomplexobj(eigenvalues):
            eigenvalues = eigenvalues.real

        sorted_eigenvalues = np.sort(eigenvalues)[-1::-1]
        return {
            "BCUTZ-1h": float(sorted_eigenvalues[0]),
            "BCUTZ-1l": float(sorted_eigenvalues[-1]),
        }

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

    @cached_property
    def estate_atom_type_counts(self) -> dict[str, int]:
        counts = {name: 0 for name in ESTATE_ATOM_TYPE_DESCRIPTORS}
        for atom_types in AtomTypes.TypeAtoms(self.mol):
            for atom_type in atom_types:
                name = f"N{atom_type}"
                if name in counts:
                    counts[name] += 1
        return counts


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
    if radius == 0:
        return float("nan")
    return (_diameter(ctx) - radius) / radius


def _petitjean_index(ctx: _DescriptorContext) -> float:
    diameter = _diameter(ctx)
    if diameter == 0:
        return float("nan")
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
    if any(valence == 0 for valence in ctx.adjacency_valences):
        return float("nan")
    return sum(valence**-2 for valence in ctx.adjacency_valences)


def _modified_zagreb_index_2(ctx: _DescriptorContext) -> float:
    valences = ctx.adjacency_valences
    return float(
        sum(
            (valences[bond.GetBeginAtomIdx()] * valences[bond.GetEndAtomIdx()]) ** -1
            for bond in ctx.bonds
        )
    )


def _shortest_path(
    graph: dict[tuple[str, int], set[tuple[str, int]]],
    start: tuple[str, int],
    end: tuple[str, int],
) -> tuple[tuple[str, int], ...] | None:
    queue = [(start, (start,))]
    seen = {start}
    for node, path in queue:
        if node == end:
            return path
        for neighbor in graph.get(node, set()):
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append((neighbor, (*path, neighbor)))
    return None


def _atom_bond_connectivity_index(ctx: _DescriptorContext) -> float:
    valences = ctx.adjacency_valences
    value = 0.0
    for bond in ctx.bonds:
        begin_valence = valences[bond.GetBeginAtomIdx()]
        end_valence = valences[bond.GetEndAtomIdx()]
        denominator = begin_valence * end_valence
        if denominator:
            value += math.sqrt((begin_valence + end_valence - 2) / denominator)
    return value


def _graovac_ghorbani_index(ctx: _DescriptorContext) -> float:
    distance_matrix = ctx.distance_matrix
    value = 0.0
    for bond in ctx.bonds:
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        begin_count = int(np.sum(distance_matrix[begin, :] < distance_matrix[end, :]))
        end_count = int(np.sum(distance_matrix[end, :] < distance_matrix[begin, :]))
        denominator = begin_count * end_count
        if denominator:
            value += math.sqrt((begin_count + end_count - 2) / denominator)
    return value


def _eccentric_connectivity_index(ctx: _DescriptorContext) -> int:
    if ctx.distance_matrix.size == 0:
        return 0
    return int(
        sum(
            valence * eccentricity
            for valence, eccentricity in zip(
                ctx.adjacency_valences,
                ctx.distance_matrix.max(axis=0),
                strict=True,
            )
        )
    )


def _fragment_complexity(ctx: _DescriptorContext) -> float:
    atom_count = ctx.mol.GetNumAtoms()
    bond_count = ctx.mol.GetNumBonds()
    hetero_atom_count = sum(1 for atom in ctx.atoms if atom.GetAtomicNum() != 6)
    return abs(bond_count**2 - atom_count**2 + atom_count) + hetero_atom_count / 100


def _framework_molecular_fraction(ctx: _DescriptorContext) -> float:
    atom_count = ctx.total_atom_count_including_hydrogen
    if atom_count == 0:
        return float("nan")
    framework_atoms = set().union(*ctx.ring_atom_sets) if ctx.ring_atom_sets else set()
    framework_atoms.update(ctx.framework_linker_atoms)
    return len(framework_atoms) / atom_count


def _atomic_property_value(table: dict[int, float], atomic_num: int) -> float:
    return table.get(atomic_num, float("nan"))


def _atomic_polarizability(ctx: _DescriptorContext) -> float:
    return sum(
        _atomic_property_value(_POLARIZABILITY_94_BY_ATOMIC_NUM, atom.GetAtomicNum())
        for atom in ctx.explicit_hydrogen_mol.GetAtoms()
    )


def _bond_polarizability(ctx: _DescriptorContext) -> float:
    value = 0.0
    for bond in ctx.explicit_hydrogen_mol.GetBonds():
        begin = bond.GetBeginAtom().GetAtomicNum()
        end = bond.GetEndAtom().GetAtomicNum()
        begin_pol = _atomic_property_value(_POLARIZABILITY_94_BY_ATOMIC_NUM, begin)
        end_pol = _atomic_property_value(_POLARIZABILITY_94_BY_ATOMIC_NUM, end)
        value += abs(begin_pol - end_pol)
    return value


def _mcgowan_volume(ctx: _DescriptorContext) -> float:
    mol = ctx.explicit_hydrogen_mol
    atom_sum = sum(
        _atomic_property_value(_MCGOWAN_VOLUME_BY_ATOMIC_NUM, atom.GetAtomicNum())
        for atom in mol.GetAtoms()
    )
    return atom_sum - mol.GetNumBonds() * 6.56


def _vabc_volume(ctx: _DescriptorContext) -> float:
    mol = ctx.explicit_hydrogen_mol
    atom_contribution = sum(
        _atomic_property_value(
            _VABC_ATOM_CONTRIBUTION_BY_ATOMIC_NUM,
            atom.GetAtomicNum(),
        )
        for atom in mol.GetAtoms()
    )
    aromatic_ring_count = _ring_metadata_count(
        ctx,
        size=None,
        fused=None,
        aromaticity="a",
        hetero=None,
    )
    aliphatic_ring_count = _ring_metadata_count(
        ctx,
        size=None,
        fused=None,
        aromaticity="A",
        hetero=None,
    )
    return (
        atom_contribution
        - 5.92 * mol.GetNumBonds()
        - 14.7 * aromatic_ring_count
        - 3.8 * aliphatic_ring_count
    )


def _hydrogen_atom_count(ctx: _DescriptorContext) -> int:
    return sum(
        1 if atom.GetAtomicNum() == 1 else atom.GetTotalNumHs()
        for atom in ctx.atoms
    )


def _rotatable_bond_ratio(ctx: _DescriptorContext) -> float:
    bond_count = len(ctx.bonds)
    if bond_count == 0:
        return float("nan")
    return rdMolDescriptors.CalcNumRotatableBonds(ctx.mol) / bond_count


def _autocorrelation_descriptor(name: str) -> DescriptorFunction:
    if _AUTOCORRELATION_PATTERN.match(name) is None:
        msg = f"unsupported autocorrelation descriptor: {name}"
        raise ValueError(msg)

    return lambda ctx: ctx.autocorrelation_values[name]


def _bcut_z_descriptor(name: str) -> DescriptorFunction:
    return lambda ctx: ctx.bcut_z_values[name]


def _atom_count_by_symbol(symbol: str) -> DescriptorFunction:
    return lambda ctx: ctx.atom_symbol_counts.get(symbol, 0)


def _halogen_atom_count(ctx: _DescriptorContext) -> int:
    return ctx.halogen_count


def _aromatic_atom_count(ctx: _DescriptorContext) -> int:
    return ctx.aromatic_atom_count


def _estate_atom_type_count(name: str) -> DescriptorFunction:
    return lambda ctx: ctx.estate_atom_type_counts[name]


def _smarts_count_descriptor(smarts: tuple[str, ...]) -> DescriptorFunction:
    pattern = Chem.MolFromSmarts("[" + ",".join(f"$({value})" for value in smarts) + "]")
    if pattern is None:
        msg = f"invalid SMARTS pattern set: {smarts!r}"
        raise ValueError(msg)

    return lambda ctx: len(ctx.mol.GetSubstructMatches(pattern))


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


def _ring_metadata_matches_filters(
    metadata: _RingMetadata,
    size: str | None,
    aromaticity: str | None,
    hetero: str | None,
) -> bool:
    if size == "G12":
        if metadata.size < 12:
            return False
    elif size is not None and metadata.size != int(size):
        return False

    if aromaticity == "a" and not metadata.is_aromatic:
        return False
    if aromaticity == "A" and metadata.is_aromatic:
        return False

    if hetero is not None and not metadata.has_hetero:
        return False

    return True


def _ring_metadata_count(
    ctx: _DescriptorContext,
    size: str | None,
    fused: str | None,
    aromaticity: str | None,
    hetero: str | None,
) -> int:
    metadata_values = ctx.fused_ring_metadata if fused is not None else ctx.ring_metadata
    return sum(
        1
        for metadata in metadata_values
        if _ring_metadata_matches_filters(metadata, size, aromaticity, hetero)
    )


def _ring_count_descriptor(name: str) -> DescriptorFunction:
    match = _RING_COUNT_PATTERN.match(name)
    if match is None:
        msg = f"unsupported ring count descriptor: {name}"
        raise ValueError(msg)

    size, fused, aromaticity, hetero = match.groups()

    def calc(ctx: _DescriptorContext) -> int:
        return _ring_metadata_count(ctx, size, fused, aromaticity, hetero)

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
    "ABC": _atom_bond_connectivity_index,
    "ABCGG": _graovac_ghorbani_index,
    "AMW": _average_molecular_weight,
    "BertzCT": _rdkit_descriptor(Descriptors.BertzCT),
    "Diameter": _diameter,
    "ECIndex": _eccentric_connectivity_index,
    "FCSP3": _rdkit_descriptor(rdMolDescriptors.CalcFractionCSP3),
    "MW": lambda ctx: ctx.exact_molecular_weight,
    "PetitjeanIndex": _petitjean_index,
    "Radius": _radius,
    "RotRatio": _rotatable_bond_ratio,
    "SMR": _rdkit_descriptor(Crippen.MolMR),
    "SLogP": _rdkit_descriptor(Crippen.MolLogP),
    "TopoPSA": lambda ctx: rdMolDescriptors.CalcTPSA(
        ctx.mol,
        includeSandP=True,
    ),
    "TopoPSA(NO)": _rdkit_descriptor(rdMolDescriptors.CalcTPSA),
    "TopoShapeIndex": _topological_shape_index,
    "VMcGowan": _mcgowan_volume,
    "Vabc": _vabc_volume,
    "WPath": _wiener_path_index,
    "WPol": _wiener_polarity_index,
    "Xp-0d": _rdkit_descriptor(Descriptors.Chi0),
    "Xp-1d": _rdkit_descriptor(Descriptors.Chi1),
    "Zagreb1": _zagreb_index_1,
    "Zagreb2": _zagreb_index_2,
    "apol": _atomic_polarizability,
    "bpol": _bond_polarizability,
    "fMF": _framework_molecular_fraction,
    "fragCpx": _fragment_complexity,
    "mZagreb1": _modified_zagreb_index_1,
    "mZagreb2": _modified_zagreb_index_2,
    "nAcid": _smarts_count_descriptor(_ACID_GROUP_SMARTS),
    "nAromAtom": _aromatic_atom_count,
    "nAromBond": lambda ctx: _bond_count_by_predicate(ctx, _is_aromatic_bond),
    "nAtom": lambda ctx: ctx.total_atom_count_including_hydrogen,
    "nBase": _smarts_count_descriptor(_BASE_GROUP_SMARTS),
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

for _name in ESTATE_ATOM_TYPE_DESCRIPTORS:
    _DESCRIPTOR_FUNCTIONS[_name] = _estate_atom_type_count(_name)

for _name in RING_COUNT_DESCRIPTORS:
    _DESCRIPTOR_FUNCTIONS[_name] = _ring_count_descriptor(_name)

for _name in PATH_COUNT_DESCRIPTORS:
    _DESCRIPTOR_FUNCTIONS[_name] = _path_count_descriptor(_name)

for _name in WALK_COUNT_DESCRIPTORS:
    _DESCRIPTOR_FUNCTIONS[_name] = _walk_count_descriptor(_name)

for _name in (*AUTOCORRELATION_Z_DESCRIPTORS, *AUTOCORRELATION_M_DESCRIPTORS):
    _DESCRIPTOR_FUNCTIONS[_name] = _autocorrelation_descriptor(_name)

for _name in BCUT_Z_DESCRIPTORS:
    _DESCRIPTOR_FUNCTIONS[_name] = _bcut_z_descriptor(_name)

for _name in SMALL_GRAPH_FORMULA_DESCRIPTORS:
    if _name not in _DESCRIPTOR_FUNCTIONS:
        msg = f"missing small graph/formula implementation: {_name}"
        raise RuntimeError(msg)

for _name in PHYSICAL_PROPERTY_DESCRIPTORS:
    if _name not in _DESCRIPTOR_FUNCTIONS:
        msg = f"missing physical-property implementation: {_name}"
        raise RuntimeError(msg)

for _name in SUPPORTED_MORDRED_2D_DESCRIPTORS:
    if _name not in _DESCRIPTOR_FUNCTIONS and hasattr(Descriptors, _name):
        _DESCRIPTOR_FUNCTIONS[_name] = _rdkit_descriptor(getattr(Descriptors, _name))


def calc_rdkit_mordred_like_2d(mol: Chem.Mol) -> dict[str, DescriptorValue]:
    """Return Mordred-name 2D descriptors computed using RDKit only."""

    if mol is None:
        msg = "mol must be an RDKit Mol, not None"
        raise ValueError(msg)

    context = _DescriptorContext(mol)
    return {
        name: _DESCRIPTOR_FUNCTIONS[name](context)
        for name in SUPPORTED_MORDRED_2D_DESCRIPTORS
    }
