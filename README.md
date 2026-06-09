# descriptors

RDKit-only Mordred-name 2D descriptor subset.

Production code in `descriptors/` does not import `mordred` or
`mordredcommunity`. Mordred is used only by tests as a numerical reference.

## Compatibility Contract

The locked supported list is in `tests/expected_supported.json` and is exposed
in code as `SUPPORTED_MORDRED_2D_DESCRIPTORS`.

Supported descriptors use Mordred names and are validated against Mordred on the
test panel. When Mordred returns a numeric value, the RDKit-only value must
match within tolerance. When Mordred returns a missing value, this package may
return either `NaN` for an undefined descriptor or a documented RDKit numeric
value when the RDKit implementation is meaningful and deterministic. Returning
missing where Mordred returns numeric is a compatibility failure.

## Supported Descriptor Groups

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
- Rotatable bond ratio: `RotRatio`
- Path counts: `MPC2` through `MPC10`, `TMPC10`, `piPC1` through
  `piPC10`, `TpiPC10`
- Walk counts: `MWC01` through `MWC10`, `TMWC10`, `SRW02` through
  `SRW10`, `TSRW10`
- EState atom-type counts: the full Mordred `N...` family, including common
  organic types such as `NsCH3`, `NaaCH`, `NssssC`, `NsOH`, `NdO`, `NssO`,
  `NaaN`, `NddsN`, `NsSH`, `NssS`, `NsF`, `NsCl`, `NsBr`, and `NsI`
- Acid/base group counts: `nAcid` and `nBase`, implemented from the Mordred
  SMARTS definitions with RDKit substructure matching
- Chi topological indices: `Xp-0d` and `Xp-1d`, implemented as validated RDKit
  `Chi0` and `Chi1` aliases

`BalabanJ` is intentionally not included. RDKit and Mordred values did not
match on the validation panel. RDKit `TPSA` is also not included under that
name because Mordred exposes the compatible descriptors as `TopoPSA` and
`TopoPSA(NO)`. Other RDKit-native `Chi*`, `Kappa*`, `fr_*`, autocorrelation,
and `BCUT2D_*` descriptors are intentionally excluded from the Mordred-name
output unless separately validated. Future functional-group expansion should
start from Mordred descriptor names and targeted counterexamples, not from bulk
RDKit `fr_*` helpers.

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
