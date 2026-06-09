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

- Supported Mordred-compatible descriptors: 291.
- Validation molecules: 52.
- Compatibility-checked numerical values: 15,132.
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

4. Next: treat `Chi*`, `Kappa*`, `fr_*`, `AUTOCORR2D`, and `BCUT2D_*` as separate
   investigations. Do not add them in bulk to the Mordred-compatible output.

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
