# Descriptor Expansion Plan

## Summary

Expand only the validated Mordred-compatible RDKit-only output. Keep the
public API unchanged:

```python
calc_rdkit_mordred_like_2d(mol) -> dict[str, float | int]
```

Every supported descriptor must be RDKit-only in production code, present in
Mordred, included in `tests/expected_supported.json`, and numerically validated
against the Mordred test oracle on the validation panel.

## Current State

- Supported Mordred-compatible descriptors: 361.
- Validation molecules: 52.
- Compatibility-checked numerical values: 18,772.
- Exact-name RDKit/Mordred overlap is exhausted except `BalabanJ`, which fails
  compatibility and must remain unsupported.
- `[NH4+]` is intentionally excluded from the validation panel for now because
  Mordred returns missing values for some graph-ratio descriptors on
  zero-heavy-edge molecules, and production missing-value behavior is not
  defined yet.

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

5. Next: investigate `fr_*`-style functional group counts as narrow,
   counterexample-driven batches:
   - Prefer exact Mordred atom-type descriptors when a `fr_*` helper only
     happens to overlap on a simple molecule.
   - Add panel molecules for each accepted functional group.
   - Keep failures and near misses documented outside the supported list.

6. Next: investigate topological index families separately:
   - `Chi*`
   - `Kappa*`
   - `BalabanJ` remains excluded unless a Mordred-equivalent RDKit-only
     implementation is proven.

7. Next: investigate autocorrelation families separately:
   - `ATS*`
   - `AATS*`
   - `ATSC*`
   - `AATSC*`
   - `MATS*`
   - `GATS*`
   - These are large families, so add them only after a shared implementation
     and panel coverage are in place.

8. Next: investigate BCUT descriptors separately:
   - Mordred `BCUT*` names do not directly match RDKit `BCUT2D_*` names.
   - Treat them as unsupported until a descriptor-by-descriptor numerical
     match is demonstrated.

9. Next: define missing-value policy before adding descriptors such as
   `RotRatio` or before reintroducing `[NH4+]` to the validation panel.

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
