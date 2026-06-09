# descriptors

RDKit-only Mordred-compatible 2D descriptor subset.

Production code in `descriptors/` does not import `mordred` or
`mordredcommunity`. Mordred is used only by tests as a numerical reference.

## Supported Mordred-Compatible Descriptors

The locked supported list is in `tests/expected_supported.json` and is exposed
in code as `SUPPORTED_MORDRED_2D_DESCRIPTORS`.

Exact Mordred/RDKit name match:

- `BertzCT`
- `LabuteASA`
- `PEOE_VSA1` through `PEOE_VSA13`
- `SMR_VSA1` through `SMR_VSA9`
- `SlogP_VSA1` through `SlogP_VSA11`
- `EState_VSA1` through `EState_VSA10`
- `VSA_EState1` through `VSA_EState9`

Mordred names implemented via RDKit aliases:

- `AMW` -> `ExactMolWt / total atom count including hydrogens`
- `FCSP3` -> `CalcFractionCSP3`
- `MW` -> `ExactMolWt`
- `SMR` -> `MolMR`
- `SLogP` -> `MolLogP`
- `TopoPSA` -> `CalcTPSA(includeSandP=True)`
- `TopoPSA(NO)` -> `CalcTPSA`
- `nBridgehead` -> `CalcNumBridgeheadAtoms`
- `nHBAcc` -> `CalcNumHBA`
- `nHBDon` -> `CalcNumHBD`
- `nHeavyAtom` -> `HeavyAtomCount`
- `nHetero` -> `CalcNumHeteroatoms`
- `nRot` -> `CalcNumRotatableBonds`
- `nSpiro` -> `CalcNumSpiroAtoms`

Mordred names implemented with explicit RDKit-only helpers:

- Atom counts: `nAtom`, `nAromAtom`, `nH`, `nB`, `nC`, `nN`, `nO`,
  `nS`, `nP`, `nF`, `nCl`, `nBr`, `nI`, `nX`
- Bond counts: `nBonds`, `nBondsO`, `nBondsS`, `nBondsD`, `nBondsT`,
  `nBondsA`, `nBondsM`, `nBondsKS`, `nBondsKD`, `nAromBond`
- Ring counts: the full Mordred `n*Ring` family, including size-specific,
  aromatic/aliphatic, hetero, and fused-ring-system variants
- Graph topology: `Diameter`, `Radius`, `TopoShapeIndex`, `PetitjeanIndex`,
  `WPath`, `WPol`, `Zagreb1`, `Zagreb2`, `mZagreb1`, `mZagreb2`

`BalabanJ` is intentionally not included. RDKit and Mordred values did not
match on the validation panel. RDKit `TPSA` is also not included under that
name because Mordred exposes the compatible descriptors as `TopoPSA` and
`TopoPSA(NO)`. `RotRatio` is not included because Mordred can return missing
values for zero-heavy-edge molecules and this package does not define a
production missing-value policy yet. RDKit-native `Chi*`, `Kappa*`, `fr_*`,
`AUTOCORR2D`, and `BCUT2D_*` descriptors are intentionally excluded from the
Mordred-compatible output unless separately validated.

## Usage

```python
from rdkit import Chem
from descriptors import calc_rdkit_mordred_like_2d

mol = Chem.MolFromSmiles("CCO")
values = calc_rdkit_mordred_like_2d(mol)
```

## Tests

Run:

```bash
python -m pytest -q
```
