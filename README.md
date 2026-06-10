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
- EState atom-type descriptors: the full Mordred `N...` (count), `MAX...`
  (max EState index), and `MIN...` (min EState index) families for all 79
  atom types, including common organic types such as `sCH3`, `aaCH`, `ssssC`,
  `sOH`, `dO`, `ssO`, `aaN`, `ddsN`, `sSH`, `ssS`, `sF`, `sCl`, `sBr`, and
  `sI`; `MAX*`/`MIN*` return NaN when no atom of that type is present
- Acid/base group counts: `nAcid` and `nBase`, implemented from the Mordred
  SMARTS definitions with RDKit substructure matching
- Chi connectivity indices: the full Kier-Hall family ported from Mordred's `Chi.py`
  — path (`Xp-0d` through `Xp-7dv`), averaged path (`AXp-0d` through `AXp-7dv`),
  cluster (`Xc-3d` through `Xc-6dv`), chain (`Xch-3d` through `Xch-7dv`), and
  path-cluster (`Xpc-4d` through `Xpc-6dv`) sub-families, each weighted by sigma
  electrons (`d`) or valence electrons (`dv`); all values computed via
  `FindAllSubgraphsOfLengthN` + a DFS subgraph classifier; returns NaN when any
  participating atom has property ≤ 0 (matching Mordred behavior), and NaN for
  averaged variants when no subgraphs exist (matching Mordred's zero-division
  policy)
- Autocorrelation descriptors: the Moreau-Broto (`ATS`/`AATS`, lags 0-8),
  centered (`ATSC`/`AATSC`, lags 0-8), Moran (`MATS`, lags 1-8) and Geary
  (`GATS`, lags 1-8) families, each weighted by an atomic property and suffixed
  with its Mordred property code:
  - `Z` atomic number, `m` mass, `v` van der Waals volume, `se`/`pe`/`are`
    Sanderson/Pauling/Allred-Rocow electronegativity, `p` polarizability,
    `i` ionization potential, `d` sigma-electron count, `dv` valence-electron
    count, `s` intrinsic state, and `c` Gasteiger charge
  - the charge property `c` only has centered families (`ATSC`/`AATSC`/`MATS`/
    `GATS`), matching Mordred
  - the averaged/normalized variants return `NaN` for molecules with no atom
    pairs at that graph distance (undefined in both Mordred and RDKit)
- BCUT descriptors: `BCUT*-1h` and `BCUT*-1l` for all 12 Mordred property
  suffixes (`Z`, `m`, `v`, `se`, `pe`, `are`, `p`, `i`, `d`, `dv`, `s`, `c`),
  implemented via RDKit Burden matrix eigenvalues; Gasteiger charge (`c`) uses
  the heavy-atom mol with implicit-H contribution (`_GasteigerHCharge`),
  matching Mordred's behavior
- Small graph/formula descriptors: `ABC`, `ABCGG`, `ECIndex`, `fragCpx`, and
  `fMF`
- Physical-property table descriptors: `apol`, `bpol`, `VMcGowan`, and `Vabc`
- Constitutional property sums and means: `SZ`/`MZ`, `Sm`/`Mm`, `Sv`/`Mv`,
  `Sse`/`Mse`, `Spe`/`Mpe`, `Sare`/`Mare`, `Sp`/`Mp`, `Si`/`Mi` — each
  computed as `Σ(p_i/p_carbon)` over all atoms including explicit hydrogens
- Topological charge descriptors: `GGI1`–`GGI10` (raw), `JGI1`–`JGI10`
  (mean), and `JGT10` (global), computed from the antisymmetric charge-term
  matrix `CT = A·D⁻² − (A·D⁻²)ᵀ`
- Matrix-spectral descriptors: `SpAbs_*`, `SpMax_*`, `SpDiam_*`, `SpAD_*`,
  `SpMAD_*`, `LogEE_*`, `VE1_*`/`VE2_*`/`VE3_*`, `VR1_*`/`VR2_*`/`VR3_*`,
  and `SM1_*` (where applicable), each computed over four matrix types:
  - `_A` — unweighted heavy-atom adjacency matrix
  - `_D` — heavy-atom topological distance matrix (already available from RDKit)
  - `_Dt` — detour matrix (longest simple paths), computed by DFS from every
    source atom; undefined for disconnected molecules and single-atom molecules
    (`SM1_Dt` only, matching Mordred's behavior)
  - `_DzZ`, `_Dzm`, `_Dzv`, `_Dzse`, `_Dzpe`, `_Dzare`, `_Dzp`, `_Dzi` —
    Barysz matrices weighted by the eight table-based atomic properties
    (atomic number, mass, vdW volume, Sanderson/Pauling/Allred-Rocow
    electronegativities, polarizability, ionization potential); off-diagonal
    entries are Floyd-Warshall shortest paths with bond weights
    `C²/(P[i]·P[j]·π_ij)`, diagonal `1 − C/P[i]` (C = carbon reference)
  - `SM1` is defined only for `_Dt` and `_Dz*` (not `_A` or `_D`, where the
    diagonal is always zero so the trace is always zero)

`BalabanJ` is intentionally not included. RDKit and Mordred values did not
match on the validation panel. RDKit `TPSA` is also not included under that
name because Mordred exposes the compatible descriptors as `TopoPSA` and
`TopoPSA(NO)`. `Vabc` returns documented `NaN` for atoms outside its Bondi
radius table, such as iodine in the validation panel. RDKit-native `Kappa*`,
`fr_*`, and `BCUT2D_*` descriptors are intentionally excluded from the
Mordred-name output unless separately validated. Future functional-group
expansion should start from Mordred descriptor names and targeted
counterexamples, not from bulk RDKit `fr_*` helpers.

## Usage

```python
from rdkit import Chem
from descriptors import calc_rdkit_mordred_like_2d

mol = Chem.MolFromSmiles("CCO")
values = calc_rdkit_mordred_like_2d(mol)
```

Pass `names` to compute only a subset. Because each descriptor family is a
separately cached calculation on the shared per-molecule context, requesting a
subset skips the work of any family no requested descriptor touches:

```python
subset = calc_rdkit_mordred_like_2d(mol, names=["MW", "Xp-2d", "ATS0Z"])
```

## Tests

Run:

```bash
python -m pytest -q
```
