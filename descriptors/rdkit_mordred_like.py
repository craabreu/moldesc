"""Mordred-compatible 2D descriptors computed with RDKit only."""

from __future__ import annotations

from collections.abc import Callable

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

from .mordred_rdkit_registry import SUPPORTED_MORDRED_2D_DESCRIPTORS

DescriptorValue = float | int
DescriptorFunction = Callable[[Chem.Mol], DescriptorValue]


def _total_atom_count_including_hydrogen(mol: Chem.Mol) -> int:
    return mol.GetNumAtoms() + sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())


def _average_molecular_weight(mol: Chem.Mol) -> float:
    atom_count = _total_atom_count_including_hydrogen(mol)
    if atom_count == 0:
        return float("nan")
    return Descriptors.ExactMolWt(mol) / atom_count


_DESCRIPTOR_FUNCTIONS: dict[str, DescriptorFunction] = {
    "AMW": _average_molecular_weight,
    "BertzCT": Descriptors.BertzCT,
    "FCSP3": rdMolDescriptors.CalcFractionCSP3,
    "MW": Descriptors.ExactMolWt,
    "SMR": Crippen.MolMR,
    "SLogP": Crippen.MolLogP,
    "TopoPSA": lambda mol: rdMolDescriptors.CalcTPSA(mol, includeSandP=True),
    "TopoPSA(NO)": rdMolDescriptors.CalcTPSA,
    "nBridgehead": rdMolDescriptors.CalcNumBridgeheadAtoms,
    "nHBAcc": rdMolDescriptors.CalcNumHBA,
    "nHBDon": rdMolDescriptors.CalcNumHBD,
    "nHeavyAtom": Descriptors.HeavyAtomCount,
    "nHetero": rdMolDescriptors.CalcNumHeteroatoms,
    "nRot": rdMolDescriptors.CalcNumRotatableBonds,
    "nSpiro": rdMolDescriptors.CalcNumSpiroAtoms,
}

for _name in SUPPORTED_MORDRED_2D_DESCRIPTORS:
    if _name not in _DESCRIPTOR_FUNCTIONS and hasattr(Descriptors, _name):
        _DESCRIPTOR_FUNCTIONS[_name] = getattr(Descriptors, _name)


def calc_rdkit_mordred_like_2d(mol: Chem.Mol) -> dict[str, DescriptorValue]:
    """Return Mordred-compatible 2D descriptors computed using RDKit only."""

    if mol is None:
        msg = "mol must be an RDKit Mol, not None"
        raise ValueError(msg)

    return {name: _DESCRIPTOR_FUNCTIONS[name](mol) for name in SUPPORTED_MORDRED_2D_DESCRIPTORS}
