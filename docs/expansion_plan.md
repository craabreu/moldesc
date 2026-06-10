# Descriptor Expansion Plan

## Summary

Expand only the validated Mordred-name RDKit-only output. Keep the
public API unchanged:

```python
calc_rdkit_mordred_like_2d(mol) -> dict[str, float | int]
```

Every supported descriptor must be RDKit-only in production code, present in
Mordred, included in `tests/expected_supported.json`, and validated against the
Mordred test oracle on the validation panel. When Mordred returns numeric
values, RDKit values must match within tolerance. When Mordred returns missing
values, RDKit may return either a documented numeric improvement or `NaN` for a
descriptor that is undefined in both implementations.

## Current State

- Supported Mordred-name descriptors: 394.
- Validation molecules: 65.
- Compatibility-checked panel cases: 25,610.
- Exact-name RDKit/Mordred overlap is exhausted except `BalabanJ`, which fails
  compatibility and must remain unsupported.
- `[NH4+]` and methane are included in the validation panel to lock documented
  behavior when Mordred returns missing values for zero-heavy-edge molecules.

## Expansion Order

1. Completed: add graph-topology descriptors that can be computed directly
   from RDKit graph distance and adjacency data:
   - `Diameter`
   - `Radius`
   - `TopoShapeIndex`
   - `PetitjeanIndex`
   - `WPath`
   - `WPol`
   - `Zagreb1`
   - `Zagreb2`
   - `mZagreb1`
   - `mZagreb2`

2. Completed: add path and walk count families:
   - `MPC*`
   - `MWC*`
   - `piPC*`, `TMPC10`, `TpiPC10`, `SRW*`, `TMWC10`, and `TSRW10`

3. Completed: investigate small functional-group aliases only with targeted counterexample
   molecules:
   - nitrile, nitro, ether, thiol, sulfide, and aromatic heteroatom aliases.
   - Added only aliases that passed targeted counterexamples: `NaaNH`, `NaaO`,
     `NddsN`, `NsNH2`, `NsSH`, `NssO`, `NssS`, `NtN`, and `NtsC`.

4. Completed: expand Mordred EState atom-type count descriptors:
   - Added the full `N...` EState atom-type family using RDKit
     `AtomTypes.TypeAtoms`.
   - This replaces the earlier small hand-picked functional aliases with a
     general atom-type implementation.
   - All 79 `N...` descriptors match Mordred on the validation panel.

5. Completed: investigate `fr_*`-style functional group counts as narrow,
   counterexample-driven batches:
   - Prefer exact Mordred atom-type descriptors when a `fr_*` helper only
     happens to overlap on a simple molecule.
   - Added `nAcid` and `nBase` using the Mordred AcidBase SMARTS definitions
     with RDKit substructure matching.
   - Added targeted panel molecules for carboxylates, tetrazole, aliphatic
     amines, guanidine, sulfonamide, and aromatic/amide counterexamples.
   - Decision for future expansion: do not bulk-add RDKit `fr_*` helpers.
     Only add functional-group descriptors when a Mordred descriptor name has
     a proven RDKit-only implementation and targeted counterexamples.

6. Completed: investigate topological index families separately:
   - Added `Xp-1d` as a validated alias for RDKit `Chi1`.
   - Added acyclic alkane panel molecules to cover simple Chi path behavior.
   - `Kier1`, `Kier2`, and `Kier3` remain excluded because RDKit `Kappa*`
     values do not match Mordred `Kier*` values when Mordred returns numeric
     values.
   - `BalabanJ` remains excluded because RDKit and Mordred values differ on
     many validation molecules.

7. Completed: define missing-value policy:
   - Mordred numeric values must match RDKit numeric values within tolerance.
   - Mordred missing values are allowed only when explicitly documented in the
     compatibility test expectations.
   - RDKit numeric values are allowed for documented Mordred-missing cases when
     the RDKit implementation is meaningful and deterministic.
   - RDKit `NaN` is allowed for documented cases that are undefined in both
     implementations.
   - Added `RotRatio` and `Xp-0d`; reintroduced methane and `[NH4+]` to the
     validation panel.

8. Completed: investigate autocorrelation families separately:
   - Added a shared RDKit-only implementation for atomic-number (`Z`)
     autocorrelation descriptors using explicit hydrogens and graph distances.
   - Added `ATS0Z` through `ATS8Z` and `ATSC0Z` through `ATSC8Z`.
   - Added low-lag averaged and normalized descriptors with no missing cases on
     the validation panel: `AATS0Z` through `AATS2Z`, `AATSC0Z` through
     `AATSC2Z`, `MATS1Z`, `MATS2Z`, `GATS1Z`, and `GATS2Z`.
   - Deferred higher-lag averaged `Z` descriptors and non-`Z` properties because
     they require broader documented missing-value expectations and property
     vector validation.

9. Next: investigate BCUT descriptors separately:
   - Mordred `BCUT*` names do not directly match RDKit `BCUT2D_*` names.
   - Treat them as unsupported until a descriptor-by-descriptor numerical
     match is demonstrated.

## Implementation Rules

- Production code must not import `mordred` or `mordredcommunity`.
- Add explicit descriptor groups to `mordred_rdkit_registry.py`.
- Add small RDKit-only helper functions to `rdkit_mordred_like.py`.
- Update `tests/expected_supported.json` only after compatibility tests pass.
- Update the README with newly supported groups and exclusions.

## Test Plan

- Use Mordred through `tests/mordred_reference.py` only.
- Verify every supported descriptor:
  - exists in Mordred,
  - is produced by the RDKit-only calculator,
  - matches Mordred on the validation molecule panel,
  - remains in the locked supported list.
- Run:

```bash
python -m pytest -q
```
