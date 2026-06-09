"""Registry for Mordred-compatible descriptors implemented with RDKit."""

from __future__ import annotations

MORDRED_RDKIT_ALIASES: dict[str, str] = {
    "MW": "ExactMolWt",
    "AMW": "ExactMolWt / total atom count including hydrogens",
    "nHeavyAtom": "HeavyAtomCount",
    "nHetero": "CalcNumHeteroatoms",
    "SLogP": "MolLogP",
    "SMR": "MolMR",
    "TopoPSA(NO)": "CalcTPSA",
    "TopoPSA": "CalcTPSA(includeSandP=True)",
    "nHBAcc": "CalcNumHBA",
    "nHBDon": "CalcNumHBD",
    "nRot": "CalcNumRotatableBonds",
    "nSpiro": "CalcNumSpiroAtoms",
    "nBridgehead": "CalcNumBridgeheadAtoms",
    "FCSP3": "CalcFractionCSP3",
}

EXACT_NAME_RDKIT_DESCRIPTORS: tuple[str, ...] = (
    "BertzCT",
    "LabuteASA",
    *(f"PEOE_VSA{i}" for i in range(1, 14)),
    *(f"SMR_VSA{i}" for i in range(1, 10)),
    *(f"SlogP_VSA{i}" for i in range(1, 12)),
    *(f"EState_VSA{i}" for i in range(1, 11)),
    *(f"VSA_EState{i}" for i in range(1, 10)),
)

SUPPORTED_MORDRED_2D_DESCRIPTORS: tuple[str, ...] = tuple(
    sorted((*MORDRED_RDKIT_ALIASES, *EXACT_NAME_RDKIT_DESCRIPTORS))
)
