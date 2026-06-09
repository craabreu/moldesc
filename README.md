# descriptors

RDKit-only Mordred-compatible 2D descriptor subset.

Production code in `descriptors/` does not import `mordred` or
`mordredcommunity`. Mordred is used only by tests as a numerical reference.

## Supported Mordred-Compatible Descriptors

The locked supported list is in `tests/expected_supported.json` and is exposed
in code as `SUPPORTED_MORDRED_2D_DESCRIPTORS`.

Exact Mordred/RDKit name match:

- `BertzCT`

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

`BalabanJ` is intentionally not included. RDKit and Mordred values did not
match on the validation panel.

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

