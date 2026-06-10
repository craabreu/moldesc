"""Registry for Mordred-name descriptors implemented with RDKit."""

from __future__ import annotations

MORDRED_RDKIT_ALIASES: dict[str, str] = {
    "MW": "ExactMolWt",
    "AMW": "ExactMolWt / total atom count including hydrogens",
    "nAtom": "atom count including implicit hydrogens",
    "nAromAtom": "aromatic atom count",
    "nH": "hydrogen atom count including implicit hydrogens",
    "nB": "boron atom count",
    "nC": "carbon atom count",
    "nN": "nitrogen atom count",
    "nO": "oxygen atom count",
    "nS": "sulfur atom count",
    "nP": "phosphorus atom count",
    "nF": "fluorine atom count",
    "nCl": "chlorine atom count",
    "nBr": "bromine atom count",
    "nI": "iodine atom count",
    "nX": "halogen atom count",
    "nHeavyAtom": "HeavyAtomCount",
    "nHetero": "CalcNumHeteroatoms",
    "nBonds": "bond count including implicit hydrogen bonds",
    "nBondsO": "ordinary bond count excluding implicit hydrogen bonds",
    "nBondsS": "single bond count including implicit hydrogen bonds",
    "nBondsD": "double bond count",
    "nBondsT": "triple bond count",
    "nBondsA": "aromatic bond count",
    "nBondsM": "multiple bond count",
    "nBondsKS": "kekulized single bond count including implicit hydrogen bonds",
    "nBondsKD": "kekulized double bond count",
    "nAromBond": "aromatic bond count",
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
    "nAcid": "Mordred AcidicGroupCount SMARTS",
    "nBase": "Mordred BasicGroupCount SMARTS",
    "RotRatio": "CalcNumRotatableBonds / heavy bond count",
    "Xp-0d": "Chi0",
    "Xp-1d": "Chi1",
}

ESTATE_ATOM_TYPE_DESCRIPTORS: tuple[str, ...] = (
    "NsLi",
    "NssBe",
    "NssssBe",
    "NssBH",
    "NsssB",
    "NssssB",
    "NsCH3",
    "NdCH2",
    "NssCH2",
    "NtCH",
    "NdsCH",
    "NaaCH",
    "NsssCH",
    "NddC",
    "NtsC",
    "NdssC",
    "NaasC",
    "NaaaC",
    "NssssC",
    "NsNH3",
    "NsNH2",
    "NssNH2",
    "NdNH",
    "NssNH",
    "NaaNH",
    "NtN",
    "NsssNH",
    "NdsN",
    "NaaN",
    "NsssN",
    "NddsN",
    "NaasN",
    "NssssN",
    "NsOH",
    "NdO",
    "NssO",
    "NaaO",
    "NsF",
    "NsSiH3",
    "NssSiH2",
    "NsssSiH",
    "NssssSi",
    "NsPH2",
    "NssPH",
    "NsssP",
    "NdsssP",
    "NsssssP",
    "NsSH",
    "NdS",
    "NssS",
    "NaaS",
    "NdssS",
    "NddssS",
    "NsCl",
    "NsGeH3",
    "NssGeH2",
    "NsssGeH",
    "NssssGe",
    "NsAsH2",
    "NssAsH",
    "NsssAs",
    "NsssdAs",
    "NsssssAs",
    "NsSeH",
    "NdSe",
    "NssSe",
    "NaaSe",
    "NdssSe",
    "NddssSe",
    "NsBr",
    "NsSnH3",
    "NssSnH2",
    "NsssSnH",
    "NssssSn",
    "NsI",
    "NsPbH3",
    "NssPbH2",
    "NsssPbH",
    "NssssPb",
)

MORDRED_RDKIT_ALIASES.update(
    {name: "RDKit EState atom type count" for name in ESTATE_ATOM_TYPE_DESCRIPTORS}
)

_RING_SIZE_PREFIXES = ("", *(str(i) for i in range(3, 13)), "G12")
_FUSED_RING_SIZE_PREFIXES = ("", *(str(i) for i in range(4, 13)), "G12")
_RING_CLASS_PREFIXES = ("", "a", "A")

RING_COUNT_DESCRIPTORS: tuple[str, ...] = tuple(
    f"n{size}{ring_class}{hetero}Ring"
    for size in _RING_SIZE_PREFIXES
    for ring_class in _RING_CLASS_PREFIXES
    for hetero in ("", "H")
) + tuple(
    f"n{size}F{ring_class}{hetero}Ring"
    for size in _FUSED_RING_SIZE_PREFIXES
    for ring_class in _RING_CLASS_PREFIXES
    for hetero in ("", "H")
)

MORDRED_RDKIT_ALIASES.update(
    {name: "RDKit ring info count" for name in RING_COUNT_DESCRIPTORS}
)

GRAPH_TOPOLOGY_DESCRIPTORS: tuple[str, ...] = (
    "Diameter",
    "Radius",
    "TopoShapeIndex",
    "PetitjeanIndex",
    "WPath",
    "WPol",
    "Zagreb1",
    "Zagreb2",
    "mZagreb1",
    "mZagreb2",
)

MORDRED_RDKIT_ALIASES.update(
    {
        name: "RDKit graph distance/adjacency calculation"
        for name in GRAPH_TOPOLOGY_DESCRIPTORS
    }
)

PATH_COUNT_DESCRIPTORS: tuple[str, ...] = (
    *(f"MPC{i}" for i in range(2, 11)),
    "TMPC10",
    *(f"piPC{i}" for i in range(1, 11)),
    "TpiPC10",
)

WALK_COUNT_DESCRIPTORS: tuple[str, ...] = (
    *(f"MWC{i:02d}" for i in range(1, 11)),
    "TMWC10",
    *(f"SRW{i:02d}" for i in range(2, 11)),
    "TSRW10",
)

# Moreau-Broto (ATS/AATS), centered (ATSC/AATSC), Moran (MATS) and Geary (GATS)
# autocorrelation families. Uncentered families are undefined for the Gasteiger
# charge property, so it only appears in the centered families.
_AUTOCORRELATION_FULL_FAMILIES: tuple[str, ...] = (
    *(f"ATS{i}" for i in range(0, 9)),
    *(f"ATSC{i}" for i in range(0, 9)),
    *(f"AATS{i}" for i in range(0, 9)),
    *(f"AATSC{i}" for i in range(0, 9)),
    *(f"MATS{i}" for i in range(1, 9)),
    *(f"GATS{i}" for i in range(1, 9)),
)
_AUTOCORRELATION_CENTERED_FAMILIES: tuple[str, ...] = (
    *(f"ATSC{i}" for i in range(0, 9)),
    *(f"AATSC{i}" for i in range(0, 9)),
    *(f"MATS{i}" for i in range(1, 9)),
    *(f"GATS{i}" for i in range(1, 9)),
)

# Property suffix -> the families it appears in, plus a human-readable label.
AUTOCORRELATION_PROPERTIES: dict[str, tuple[tuple[str, ...], str]] = {
    "Z": (_AUTOCORRELATION_FULL_FAMILIES, "atomic-number"),
    "m": (_AUTOCORRELATION_FULL_FAMILIES, "atomic-mass"),
    "v": (_AUTOCORRELATION_FULL_FAMILIES, "van der Waals volume"),
    "se": (_AUTOCORRELATION_FULL_FAMILIES, "Sanderson electronegativity"),
    "pe": (_AUTOCORRELATION_FULL_FAMILIES, "Pauling electronegativity"),
    "are": (_AUTOCORRELATION_FULL_FAMILIES, "Allred-Rocow electronegativity"),
    "p": (_AUTOCORRELATION_FULL_FAMILIES, "polarizability"),
    "i": (_AUTOCORRELATION_FULL_FAMILIES, "ionization potential"),
    "d": (_AUTOCORRELATION_FULL_FAMILIES, "sigma-electron count"),
    "dv": (_AUTOCORRELATION_FULL_FAMILIES, "valence-electron count"),
    "s": (_AUTOCORRELATION_FULL_FAMILIES, "intrinsic state"),
    "c": (_AUTOCORRELATION_CENTERED_FAMILIES, "Gasteiger charge"),
}

AUTOCORRELATION_DESCRIPTORS: tuple[str, ...] = tuple(
    f"{family}{suffix}"
    for suffix, (families, _label) in AUTOCORRELATION_PROPERTIES.items()
    for family in families
)

_BCUT_PROPERTIES: tuple[str, ...] = (
    "Z", "m", "v", "se", "pe", "are", "p", "i", "d", "dv", "s", "c",
)

BCUT_DESCRIPTORS: tuple[str, ...] = tuple(
    f"BCUT{prop}{suffix}"
    for prop in _BCUT_PROPERTIES
    for suffix in ("-1h", "-1l")
)

SMALL_GRAPH_FORMULA_DESCRIPTORS: tuple[str, ...] = (
    "ABC",
    "ABCGG",
    "ECIndex",
    "fragCpx",
    "fMF",
)

PHYSICAL_PROPERTY_DESCRIPTORS: tuple[str, ...] = (
    "VMcGowan",
    "Vabc",
    "apol",
    "bpol",
)

# Constitutional property sums (S*) and means (M*), normalized to carbon.
# Formula: S_p = Σ(p_i / p_C) over all atoms including explicit H;
# M_p = S_p / total atom count (including H).
_CONSTITUTIONAL_PROPERTIES: tuple[str, ...] = (
    "Z", "m", "v", "se", "pe", "are", "p", "i",
)

CONSTITUTIONAL_DESCRIPTORS: tuple[str, ...] = tuple(
    f"{prefix}{prop}"
    for prefix in ("S", "M")
    for prop in _CONSTITUTIONAL_PROPERTIES
)

TOPOLOGICAL_CHARGE_DESCRIPTORS: tuple[str, ...] = (
    *(f"GGI{k}" for k in range(1, 11)),
    *(f"JGI{k}" for k in range(1, 11)),
    "JGT10",
)

MORDRED_RDKIT_ALIASES.update(
    {name: "RDKit path count calculation" for name in PATH_COUNT_DESCRIPTORS}
)
MORDRED_RDKIT_ALIASES.update(
    {name: "RDKit adjacency walk count calculation" for name in WALK_COUNT_DESCRIPTORS}
)
for _suffix, (_families, _label) in AUTOCORRELATION_PROPERTIES.items():
    MORDRED_RDKIT_ALIASES.update(
        {
            f"{family}{_suffix}": f"RDKit {_label} autocorrelation calculation"
            for family in _families
        }
    )
MORDRED_RDKIT_ALIASES.update(
    {
        name: "RDKit Burden matrix eigenvalue calculation"
        for name in BCUT_DESCRIPTORS
    }
)
MORDRED_RDKIT_ALIASES.update(
    {
        name: "RDKit graph formula calculation"
        for name in SMALL_GRAPH_FORMULA_DESCRIPTORS
    }
)
MORDRED_RDKIT_ALIASES.update(
    {
        name: "RDKit explicit-hydrogen physical property calculation"
        for name in PHYSICAL_PROPERTY_DESCRIPTORS
    }
)
MORDRED_RDKIT_ALIASES.update(
    {
        name: "RDKit carbon-normalized constitutional sum/mean calculation"
        for name in CONSTITUTIONAL_DESCRIPTORS
    }
)
MORDRED_RDKIT_ALIASES.update(
    {
        name: "RDKit charge-term matrix topological charge calculation"
        for name in TOPOLOGICAL_CHARGE_DESCRIPTORS
    }
)

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
