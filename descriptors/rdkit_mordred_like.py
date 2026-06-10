"""Mordred-name 2D descriptors computed with RDKit only."""

from __future__ import annotations

from collections.abc import Callable
from functools import cached_property
import csv
import math
from pathlib import Path
import re
from typing import NamedTuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors, rdPartialCharges
from rdkit.Chem.EState import AtomTypes

from .mordred_rdkit_registry import (
    AUTOCORRELATION_DESCRIPTORS,
    BCUT_DESCRIPTORS,
    ESTATE_ATOM_TYPE_DESCRIPTORS,
    PATH_COUNT_DESCRIPTORS,
    CONSTITUTIONAL_DESCRIPTORS,
    PHYSICAL_PROPERTY_DESCRIPTORS,
    TOPOLOGICAL_CHARGE_DESCRIPTORS,
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
    r"^(AATSC|AATS|ATSC|ATS|MATS|GATS)(\d+)(are|se|pe|dv|Z|m|v|p|i|d|s|c)$"
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


def _load_atomic_property_tables() -> dict[str, dict[int, float]]:
    """Load per-element property tables from the bundled CSV data file.

    Each column becomes a ``{atomic_number: value}`` dict; blank cells are
    treated as missing so lookups fall back to NaN via
    :func:`_atomic_property_value`.
    """

    path = Path(__file__).with_name("atomic_properties.csv")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        columns = [name for name in reader.fieldnames if name != "atomic_number"]
        tables: dict[str, dict[int, float]] = {name: {} for name in columns}
        for row in reader:
            atomic_num = int(row["atomic_number"])
            for name in columns:
                cell = row[name]
                if cell != "":
                    tables[name][atomic_num] = float(cell)
    return tables


_ATOMIC_PROPERTY_TABLES = _load_atomic_property_tables()
_MASS_BY_ATOMIC_NUM = _ATOMIC_PROPERTY_TABLES["mass"]
_VDW_VOLUME_BY_ATOMIC_NUM = _ATOMIC_PROPERTY_TABLES["vdw_volume"]
_SANDERSON_EN_BY_ATOMIC_NUM = _ATOMIC_PROPERTY_TABLES["sanderson_en"]
_PAULING_EN_BY_ATOMIC_NUM = _ATOMIC_PROPERTY_TABLES["pauling_en"]
_ALLRED_ROCOW_EN_BY_ATOMIC_NUM = _ATOMIC_PROPERTY_TABLES["allred_rocow_en"]
_POLARIZABILITY_94_BY_ATOMIC_NUM = _ATOMIC_PROPERTY_TABLES["polarizability"]
_IONIZATION_POTENTIAL_BY_ATOMIC_NUM = _ATOMIC_PROPERTY_TABLES["ionization_potential"]
_MCGOWAN_VOLUME_BY_ATOMIC_NUM = _ATOMIC_PROPERTY_TABLES["mcgowan_volume"]
_PERIOD_BY_ATOMIC_NUM = _ATOMIC_PROPERTY_TABLES["period"]
_BONDI_RADII_BY_ATOMIC_NUM = _ATOMIC_PROPERTY_TABLES["bondi_radius"]

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

    def _element_property_vector(
        self, table: dict[int, float]
    ) -> tuple[float, ...]:
        return tuple(
            _atomic_property_value(table, atom.GetAtomicNum())
            for atom in self.explicit_hydrogen_mol.GetAtoms()
        )

    @cached_property
    def gasteiger_charges(self) -> tuple[float, ...]:
        # ComputeGasteigerCharges annotates atoms in place, so work on a copy to
        # avoid mutating the cached explicit-hydrogen molecule used elsewhere.
        mol = Chem.Mol(self.explicit_hydrogen_mol)
        rdPartialCharges.ComputeGasteigerCharges(mol)
        return tuple(
            float(atom.GetDoubleProp("_GasteigerCharge")) for atom in mol.GetAtoms()
        )

    @cached_property
    def autocorrelation_property_vectors(self) -> dict[str, tuple[float, ...]]:
        """Per-atom property vectors keyed by Mordred property suffix."""

        atoms = list(self.explicit_hydrogen_mol.GetAtoms())
        return {
            "Z": tuple(float(atom.GetAtomicNum()) for atom in atoms),
            "m": self._element_property_vector(_MASS_BY_ATOMIC_NUM),
            "v": self._element_property_vector(_VDW_VOLUME_BY_ATOMIC_NUM),
            "se": self._element_property_vector(_SANDERSON_EN_BY_ATOMIC_NUM),
            "pe": self._element_property_vector(_PAULING_EN_BY_ATOMIC_NUM),
            "are": self._element_property_vector(_ALLRED_ROCOW_EN_BY_ATOMIC_NUM),
            "p": self._element_property_vector(_POLARIZABILITY_94_BY_ATOMIC_NUM),
            "i": self._element_property_vector(_IONIZATION_POTENTIAL_BY_ATOMIC_NUM),
            "d": tuple(float(_sigma_electron_count(atom)) for atom in atoms),
            "dv": tuple(_valence_electron_count(atom) for atom in atoms),
            "s": tuple(_intrinsic_state(atom) for atom in atoms),
            "c": self.gasteiger_charges,
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
    def _burden_matrix_base(self) -> np.ndarray:
        n = len(self.atoms)
        B = 0.001 * np.ones((n, n))
        for bond in self.bonds:
            i = bond.GetBeginAtom().GetIdx()
            j = bond.GetEndAtom().GetIdx()
            w = bond.GetBondTypeAsDouble() / 10.0
            if bond.GetBeginAtom().GetDegree() == 1 or bond.GetEndAtom().GetDegree() == 1:
                w += 0.01
            B[i, j] = w
            B[j, i] = w
        return B

    @cached_property
    def _bcut_gasteiger_diagonal(self) -> list[float]:
        # BCUT uses the heavy-atom mol; _GasteigerHCharge holds the implicit-H contribution.
        mol = Chem.Mol(self.mol)
        rdPartialCharges.ComputeGasteigerCharges(mol)
        result = []
        for atom in mol.GetAtoms():
            q = atom.GetDoubleProp("_GasteigerCharge")
            if atom.HasProp("_GasteigerHCharge"):
                q += atom.GetDoubleProp("_GasteigerHCharge")
            result.append(q)
        return result

    @cached_property
    def bcut_values(self) -> dict[str, float]:
        atoms = self.atoms
        n = len(atoms)
        props = ("Z", "m", "v", "se", "pe", "are", "p", "i", "d", "dv", "s", "c")
        diag_matrix = np.array(
            [
                [float(a.GetAtomicNum()) for a in atoms],
                [_MASS_BY_ATOMIC_NUM.get(a.GetAtomicNum(), float("nan")) for a in atoms],
                [_VDW_VOLUME_BY_ATOMIC_NUM.get(a.GetAtomicNum(), float("nan")) for a in atoms],
                [_SANDERSON_EN_BY_ATOMIC_NUM.get(a.GetAtomicNum(), float("nan")) for a in atoms],
                [_PAULING_EN_BY_ATOMIC_NUM.get(a.GetAtomicNum(), float("nan")) for a in atoms],
                [_ALLRED_ROCOW_EN_BY_ATOMIC_NUM.get(a.GetAtomicNum(), float("nan")) for a in atoms],
                [_POLARIZABILITY_94_BY_ATOMIC_NUM.get(a.GetAtomicNum(), float("nan")) for a in atoms],
                [_IONIZATION_POTENTIAL_BY_ATOMIC_NUM.get(a.GetAtomicNum(), float("nan")) for a in atoms],
                [float(_sigma_electron_count(a)) for a in atoms],
                [_valence_electron_count(a) for a in atoms],
                [_intrinsic_state(a) for a in atoms],
                self._bcut_gasteiger_diagonal,
            ],
            dtype=float,
        )  # shape (12, n)

        # Build a batch of 12 Burden matrices differing only in their diagonals.
        nan_rows = np.isnan(diag_matrix).any(axis=1)
        safe_diag = diag_matrix.copy()
        safe_diag[nan_rows] = 0.0  # sentinel to avoid LinAlgError; results overwritten below

        batch = np.broadcast_to(self._burden_matrix_base, (len(props), n, n)).copy()
        idx = np.arange(n)
        batch[:, idx, idx] = safe_diag

        ev = np.linalg.eig(batch)[0]  # (12, n), one call for all properties
        if np.iscomplexobj(ev):
            ev = ev.real
        sorted_ev = np.sort(ev, axis=1)[:, ::-1]  # descending per row
        sorted_ev[nan_rows] = float("nan")

        results: dict[str, float] = {}
        for k, prop in enumerate(props):
            results[f"BCUT{prop}-1h"] = float(sorted_ev[k, 0])
            results[f"BCUT{prop}-1l"] = float(sorted_ev[k, -1])
        return results

    @cached_property
    def constitutional_values(self) -> dict[str, float]:
        """Constitutional sums S_p = Σ(p_i/p_C) and means M_p = S_p/A over explicit-H atoms."""
        atoms = list(self.explicit_hydrogen_mol.GetAtoms())
        n = len(atoms)
        carbon = 6
        prop_funs: list[tuple[str, dict[int, float]]] = [
            ("Z",   {z: float(z) for z in range(1, 119)}),
            ("m",   _MASS_BY_ATOMIC_NUM),
            ("v",   _VDW_VOLUME_BY_ATOMIC_NUM),
            ("se",  _SANDERSON_EN_BY_ATOMIC_NUM),
            ("pe",  _PAULING_EN_BY_ATOMIC_NUM),
            ("are", _ALLRED_ROCOW_EN_BY_ATOMIC_NUM),
            ("p",   _POLARIZABILITY_94_BY_ATOMIC_NUM),
            ("i",   _IONIZATION_POTENTIAL_BY_ATOMIC_NUM),
        ]
        results: dict[str, float] = {}
        for suffix, table in prop_funs:
            carbon_val = table[carbon]
            vals = [table.get(a.GetAtomicNum(), float("nan")) / carbon_val for a in atoms]
            s = sum(vals) if not any(math.isnan(v) for v in vals) else float("nan")
            results[f"S{suffix}"] = s
            results[f"M{suffix}"] = s / n if not math.isnan(s) else float("nan")
        return results

    @cached_property
    def topological_charge_values(self) -> dict[str, float]:
        """GGI/JGI/JGT topological charge descriptors (heavy-atom mol)."""
        D = self.distance_matrix.astype(float)
        A = self.adjacency_matrix.astype(float)
        n = len(self.atoms)

        D2 = D.copy()
        D2[D2 != 0] **= -2
        np.fill_diagonal(D2, 0)
        M = A @ D2
        CT = M - M.T  # antisymmetric charge-term matrix

        D_lower = D * np.tri(n)
        D_lower[D_lower == 0] = np.inf

        results: dict[str, float] = {}
        jgt = 0.0
        for k in range(1, 11):
            f = D_lower == k
            ct_k = CT[f]
            count_k = len(ct_k)
            results[f"GGI{k}"] = float(np.abs(ct_k).sum())
            jgi_k = float((np.abs(ct_k) / count_k).sum()) if count_k > 0 else 0.0
            results[f"JGI{k}"] = jgi_k
            jgt += jgi_k
        results["JGT10"] = jgt
        return results

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


_PERIODIC_TABLE = Chem.GetPeriodicTable()


def _sigma_electron_count(atom: Chem.Atom) -> int:
    """Number of non-hydrogen neighbors (Mordred ``d`` property)."""

    return sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() != 1)


def _valence_electron_count(atom: Chem.Atom) -> float:
    """Kier-Hall valence-electron count (Mordred ``dv`` property)."""

    atomic_num = atom.GetAtomicNum()
    if atomic_num == 1:
        return 0.0
    formal_charge = atom.GetFormalCharge()
    outer = _PERIODIC_TABLE.GetNOuterElecs(atomic_num) - formal_charge
    core = atomic_num - formal_charge
    hydrogens = atom.GetTotalNumHs() + sum(
        1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 1
    )
    return (outer - hydrogens) / (core - outer - 1)


def _intrinsic_state(atom: Chem.Atom) -> float:
    """Electrotopological intrinsic state (Mordred ``s`` property)."""

    sigma = _sigma_electron_count(atom)
    if sigma == 0:
        return float("nan")
    period = _PERIOD_BY_ATOMIC_NUM.get(atom.GetAtomicNum(), float("nan"))
    valence = _valence_electron_count(atom)
    return ((2.0 / period) ** 2 * valence + 1) / sigma


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


def _bcut_descriptor(name: str) -> DescriptorFunction:
    return lambda ctx: ctx.bcut_values[name]


def _constitutional_descriptor(name: str) -> DescriptorFunction:
    return lambda ctx: ctx.constitutional_values[name]


def _topological_charge_descriptor(name: str) -> DescriptorFunction:
    return lambda ctx: ctx.topological_charge_values[name]


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

for _name in AUTOCORRELATION_DESCRIPTORS:
    _DESCRIPTOR_FUNCTIONS[_name] = _autocorrelation_descriptor(_name)

for _name in BCUT_DESCRIPTORS:
    _DESCRIPTOR_FUNCTIONS[_name] = _bcut_descriptor(_name)

for _name in CONSTITUTIONAL_DESCRIPTORS:
    _DESCRIPTOR_FUNCTIONS[_name] = _constitutional_descriptor(_name)

for _name in TOPOLOGICAL_CHARGE_DESCRIPTORS:
    _DESCRIPTOR_FUNCTIONS[_name] = _topological_charge_descriptor(_name)

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
