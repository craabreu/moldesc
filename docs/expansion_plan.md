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

- Supported Mordred-name descriptors: 1580.
- Validation molecules: 65.
- Exact-name RDKit/Mordred overlap is exhausted except `BalabanJ`, which fails
  compatibility and must remain unsupported.
- `[NH4+]` and methane are included in the validation panel to lock documented
  behavior when Mordred returns missing values for zero-heavy-edge molecules.
- The public entry point accepts an optional `names=` subset; because each
  descriptor family is a separately cached calculation on the shared per-molecule
  context, requesting a subset skips any family no requested descriptor touches.

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
     values. `descriptastorus` wraps the same RDKit `Kappa*` functions and
     does not help. No external library provides a compatible implementation.
   - `BalabanJ` remains excluded because RDKit and Mordred values differ on
     many validation molecules.
   - RDKit `fr_*` functional-group counters (85 descriptors) are not Mordred
     descriptors at all; Mordred has no `fr_*` family. They cannot be added
     as Mordred-name descriptors and are permanently excluded from this scope.

7. Completed: define missing-value policy:
   - Mordred numeric values must match RDKit numeric values within tolerance.
   - When Mordred is missing and RDKit deliberately provides a meaningful,
     deterministic numeric value, the case must be documented in the
     `EXPECTED_RDKIT_IMPROVEMENTS` set in the compatibility tests.
   - When Mordred is missing and RDKit also returns `NaN`, both implementations
     agree the descriptor is undefined; there is no oracle value to check, so the
     case is accepted without per-case documentation.
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

9. Completed: investigate BCUT descriptors separately:
   - RDKit `BCUT2D_*` descriptors do not numerically match Mordred `BCUT*`
     names on the validation panel and are not treated as aliases.
   - Added `BCUTZ-1h` and `BCUTZ-1l` with a RDKit-only Burden matrix
     implementation using atomic numbers on the diagonal.
   - Deferred other `BCUT*` properties until Mordred atomic-property tables or
     equivalent RDKit-only property vectors are implemented and validated.

10. Completed: expand higher-lag atomic-number autocorrelation descriptors:

    - Exposed the full lag 0-8 averaged/normalized family: `AATS0Z`..`AATS8Z`,
      `AATSC0Z`..`AATSC8Z`, `MATS1Z`..`MATS8Z`, and `GATS1Z`..`GATS8Z`. The
      underlying order-sums were already computed for all lags, so this only
      removed lag gates.
    - The averaged/normalized autocorrelation family (`AATS*`, `AATSC*`,
      `MATS*`, `GATS*`) divides by the atom-pair count at a graph distance, so it
      is undefined (Mordred Missing, RDKit `NaN`) when a molecule has no atom
      pairs at that distance. Under the missing-value policy (step 7) these
      both-`NaN` cases need no enumeration; the raw `ATS*`/`ATSC*` sums stay
      numeric and are checked against Mordred as usual.

11. Completed: implement non-`Z` autocorrelation property vectors:

    - Generalized the autocorrelation machinery into a single vectorized
      `autocorrelation_values` that stacks every property vector and computes all
      properties together per lag. The shared per-lag graph-distance work and the
      Geary numerator run once per lag in numpy, so adding properties scales well.
    - Added every Mordred autocorrelation property, all validated against the
      oracle on the panel with zero mismatches:
      - per-element table properties: `m` (mass), `v` (van der Waals volume),
        `se`/`pe`/`are` (Sanderson/Pauling/Allred-Rocow electronegativity),
        `p` (polarizability), `i` (ionization potential).
      - environment-dependent properties computed per atom with RDKit: `d`
        (sigma electrons), `dv` (valence electrons), `s` (intrinsic state), and
        `c` (Gasteiger charge via `rdPartialCharges.ComputeGasteigerCharges`).
    - `c` only appears in the centered families (`ATSC`/`AATSC`/`MATS`/`GATS`),
      matching Mordred.
    - Per-element tables now live in the bundled `atomic_properties.csv` data
      file rather than hard-coded dicts.

12. Completed: expand remaining Mordred `BCUT*` properties:

    Added `BCUT*-1h` and `BCUT*-1l` for all 12 Mordred property suffixes. The
    shared Burden matrix off-diagonal (bond weights) is built once; the diagonal
    is swapped per property before solving the eigenvalue problem.

    - Table properties (`m`, `v`, `se`, `pe`, `are`, `p`, `i`) reuse the same
      per-element tables already loaded for autocorrelation.
    - Environment-dependent properties (`d`, `dv`, `s`) reuse the same per-atom
      helper functions (`_sigma_electron_count`, `_valence_electron_count`,
      `_intrinsic_state`).
    - Gasteiger charge (`c`) is computed on the heavy-atom mol (matching
      Mordred's `explicit_hydrogens = False`); the implicit-hydrogen contribution
      is captured via `_GasteigerHCharge`, which Mordred adds to `_GasteigerCharge`
      to get the effective per-atom charge.
    - Zero mismatches on the full validation panel; 4 both-NaN cases (molecules
      with atoms outside the property table, e.g. iodine-containing compounds for
      some properties) auto-accepted under the missing-value policy.

13. Completed: add small standalone graph/formula descriptors:

    Added descriptors that need only RDKit atom, bond, ring, and distance data:
    `ABC`, `ABCGG`, `ECIndex`, `fragCpx`, and `fMF`. These are useful
    non-alias additions that avoid new property-table infrastructure.

14. Completed: add simple physical-property table descriptors:

    Added `apol`, `bpol`, `VMcGowan`, and `Vabc` using explicit RDKit-only
    atomic constants and explicit-hydrogen molecules to match Mordred behavior.
    `Vabc` returns documented `NaN` for atoms outside its Bondi radius table,
    such as iodine in the validation panel.

15. Completed: add constitutional property sums and means:

    Added `SZ`, `MZ`, `Sm`, `Mm`, `Sv`, `Mv`, `Sse`, `Mse`, `Spe`, `Mpe`,
    `Sare`, `Mare`, `Sp`, `Mp`, `Si`, `Mi` — 16 descriptors covering all 8
    table-property variants in both sum and mean form.

    Formula: `S_p = Σ(p_i / p_C)` summed over all atoms including explicit
    hydrogens, normalized to the carbon reference value; `M_p = S_p / A` where
    A is total atom count including H. Reuses the same per-element property
    tables already loaded for autocorrelation and BCUT. Zero mismatches on the
    full validation panel.

16. Completed: add topological charge descriptors:

    Added `GGI1`–`GGI10` (raw sum), `JGI1`–`JGI10` (mean per distance), and
    `JGT10` (global sum of JGI1–JGI10) — 21 descriptors.

    Algorithm: build the antisymmetric charge-term matrix `CT = (A @ D⁻²) − (A @ D⁻²)ᵀ`
    from the heavy-atom adjacency and distance matrices; for each lag k extract
    lower-triangle pairs at distance k and sum `|CT_ij|` (GGI) or `|CT_ij| /
    count_k` (JGI). When no atom pairs exist at distance k the result is 0 (empty
    sum), not NaN — matching Mordred's `np.abs([]).sum() == 0` behavior. Zero
    mismatches on the full validation panel.

17. Completed: implement the full Mordred connectivity-index (Xp-*/Xc-*/Xch-*/Xpc-*/AXp-*) family:

    Added all 56 Kier-Hall chi descriptors via a clean-room port of Mordred's
    `Chi.py`. Zero mismatches on the full validation panel.

    Algorithm:
    - `FindAllSubgraphsOfLengthN(mol, order)` enumerates connected subgraphs.
    - Each subgraph is classified by DFS: a back-edge detected during traversal
      marks it as `chain`; otherwise degrees determine `path` (all ≤ 2),
      `path_cluster` (has 2 and ≥ 3), or `cluster` (only ≤ 1 and ≥ 3).
    - For each subgraph of the relevant type, compute `Σ prod(P[node])^{-0.5}`.
    - Property `d`: sigma electrons = heavy-atom neighbor count (`_sigma_electron_count`).
    - Property `dv`: Kier-Hall valence electrons = `(Zv − h) / (Z − Zv − 1)`,
      formal-charge-adjusted, already implemented as `_valence_electron_count`.
    - Returns NaN when any node's property ≤ 0 (matches Mordred's fail-on-zero-product rule).
    - Averaged (`AXp-*`) returns NaN when subgraph count = 0 (matches Mordred's
      `ZeroDivisionError → Missing` policy).

    RDKit's built-in `Chi0v` / `Chi1v` were NOT used as aliases because they
    apply a slightly different valence formula that diverges from Mordred on
    charged atoms (e.g. nitro groups, quaternary ammonium, carboxylates).
    All `Xp-*` / `AXp-*` descriptors are computed uniformly via `chi_values`.

    Prior `Xp-0d` and `Xp-1d` aliases to RDKit `Chi0` / `Chi1` were replaced
    by the same `chi_values` path; the two `EXPECTED_RDKIT_IMPROVEMENTS` entries
    for ammonium/methane were removed (now both return NaN, matching Mordred).

    Descriptors added (+54 vs prior state):
    - `Xp-0dv`..`Xp-7dv`, `Xp-2d`..`Xp-7d` (14 new path variants)
    - `AXp-0d`..`AXp-7d`, `AXp-0dv`..`AXp-7dv` (16 averaged path)
    - `Xc-3d`..`Xc-6d`, `Xc-3dv`..`Xc-6dv` (8 cluster)
    - `Xch-3d`..`Xch-7d`, `Xch-3dv`..`Xch-7dv` (10 chain)
    - `Xpc-4d`..`Xpc-6d`, `Xpc-4dv`..`Xpc-6dv` (6 path-cluster)

18. Completed: implement the full matrix-spectral family (141 descriptors):

    Added `SpAbs_*`, `SpMax_*`, `SpDiam_*`, `SpAD_*`, `SpMAD_*`, `LogEE_*`,
    `VE1_*`/`VE2_*`/`VE3_*`, `VR1_*`/`VR2_*`/`VR3_*`, and `SM1_*` (where
    applicable) across four matrix types. Zero mismatches on the full panel.

    Matrices and their RDKit source:
    - `_A` — `Chem.GetAdjacencyMatrix` (already cached for walk counts etc.)
    - `_D` — `Chem.GetDistanceMatrix` (already cached for topological indices)
    - `_Dt` — longest-simple-path matrix via DFS from each source atom.
      The DFS marks visited nodes and backtracks, accumulating the maximum
      distance reached at each target node across all simple paths; the result
      matrix is symmetrized with `np.maximum(D, D.T)`.  Returns NaN for all
      descriptors when the molecule is disconnected (Mordred `require_connected`),
      and returns NaN for `SM1_Dt` when n==1 (Mordred returns `np.int64(0)` for
      that edge case, which Python's `isinstance(v, int|float)` treats as missing,
      so the test expects NaN from our side).
    - `_DzX` (8 variants) — Barysz matrix: direct bond weights
      `C²/(P[i]·P[j]·π_ij)` (C = carbon reference; π_ij = `GetBondTypeAsDouble`)
      then Floyd-Warshall shortest paths; diagonal filled with `1 − C/P[i]`
      afterward.  Implemented with pure-numpy O(n³) Floyd-Warshall.

    Spectral aggregates (`_fill_spectral` helper, shared across all matrices):
    - SpAbs = Σ|λ|, SpMax = max λ, SpDiam = max λ − min λ
    - SpAD = Σ|λ − mean(λ)|, SpMAD = SpAD/N
    - LogEE = log(Σexp(λ) + 1) via log-sum-exp (numerically stable)
    - SM1 = trace(matrix) (no eigendecomposition needed)
    - VE1 = Σ|v_i|, VE2 = VE1/N, VE3 = log(0.1·N·VE1)
      (v = leading eigenvector; VE3 → NaN when VE1 = 0)
    - VR1 = Σ (v_i·v_j)^{−½} over bonds (Randić-like), VR2 = VR1/N,
      VR3 = log(0.1·N·VR1); VR1 → NaN if any bond product ≤ 0;
      VR3 → NaN if VR1 = 0 (matches Mordred's ZeroDivisionError→Missing)

    SM1 is defined only for `_Dt` and `_Dz*`; for `_A` and `_D` the diagonal
    is always zero so the trace is always zero and Mordred removed SM1_A/SM1_D
    in v1.1.0.

    All 65 panel molecules × 141 descriptors checked against Mordred;
    12 both-NaN cases (VR3_* for ammonium and methane, SM1_Dt for ammonium and
    methane), auto-accepted under the missing-value policy.

19. Completed: implement EState atom-type max/min/sum descriptors (237 descriptors):

    Added `MAX*`, `MIN*` (79 each), and `S*` (79) — the `max`, `min`, and `sum`
    aggregation variants of the same `AtomTypeEState` family already implemented
    for `N*` counts.

    Algorithm: one `estate_atom_type_agg` cached property pairs
    `AtomTypes.TypeAtoms(mol)` with `EStateIndices(mol)` (both on the heavy-atom
    mol) to group EState index values by atom type. From the grouped dict:
    - `MAX*` / `MIN*`: `max`/`min(values)`, NaN when no atom of that type is
      present (matching Mordred's `rethrow_na(ValueError)` policy).
    - `S*`: `sum(values)`, 0 when absent (matching Mordred's `sum([]) == 0`
      behavior — the `sum` aggregation never returns NaN).

    The 79 base types are derived from the existing `ESTATE_ATOM_TYPE_DESCRIPTORS`
    tuple (strip the `N` prefix), so the type list stays in one canonical place.
    All three aggregations share a single TypeAtoms + EStateIndices pass.

    Numeric oracle comparisons: 334 (MAX/MIN, most cells both-NaN since rare
    atom types absent from the panel) + 5,135 (S*, every cell numeric since
    sum never NaN). Zero mismatches; 9,936 both-NaN cases auto-accepted.

20. Completed: implement information content descriptors (42 descriptors):

    Added `IC0`–`IC5`, `TIC0`–`TIC5`, `SIC0`–`SIC5`, `BIC0`–`BIC5`,
    `CIC0`–`CIC5`, `MIC0`–`MIC5`, and `ZMIC0`–`ZMIC5` — Shannon entropy
    over atom-neighborhood equivalence classes at orders 0–5.

    Algorithm (matching Mordred exactly):
    - Uses explicit-hydrogen kekulized mol — Mordred's `explicit_hydrogens=True`
      and `kekulize=True` flags; atom count A includes H.
    - For order 0: atom code = atomic number.
    - For order m > 0: BFS tree rooted at each atom expanded m levels; each
      atom's code is the sorted tuple of all root-to-leaf trails, where each
      trail element is `(bond_type, (atomic_num, degree))`. Two atoms have the
      same code iff their neighborhoods are structurally identical up to depth m.
    - `IC_m = -Σ(g/A · log₂(g/A))` over equivalence-class sizes g.
    - `TIC_m = A · IC_m`
    - `SIC_m = IC_m / log₂(A)` — NaN when A=1 (single-atom mol, log₂(1)=0)
    - `BIC_m = IC_m / log₂(Σπ_b)` — NaN when no bonds (Σπ_b=0)
    - `CIC_m = log₂(A) − IC_m`
    - `MIC_m` — IC weighted by representative-atom atomic mass per group
    - `ZMIC_m` — IC weighted by group_size × representative-atom atomic_num

    2,730 numeric oracle comparisons (all cells numeric; no NaN from IC family).
    Zero mismatches on the 65-molecule panel.

21. Completed: implement extended topochemical atom (ETA/AETA, 45 descriptors)
    and molecular distance edge (MDE, 19 descriptors):

    **ETA/AETA (45):** `ETA_alpha`, `AETA_alpha`, `ETA_beta_s/ns_d/ns/beta`,
    their `AETA_*` averages, `ETA_dBeta`, `ETA_eta/eta_L/eta_R/eta_RL`,
    their `AETA_*` averages, `ETA_eta_F/FL/B/BR` and averages, `ETA_epsilon_1`–`5`,
    `ETA_dEpsilon_A–D`, `ETA_dAlpha_A/B`, `ETA_psi_1`, `ETA_dPsi_A/B`,
    `ETA_shape_p/y/x`.

    Algorithm:
    - All 45 require a connected molecule (`EtaBase.require_connected = True`);
      disconnected mols return NaN for all 45.
    - Per-atom properties on the kekulized heavy-atom mol (aromatic flags
      preserved via `Chem.Kekulize(mol)`, default `clearAromaticFlags=False`):
      - `alpha_i = (Z − Zv) / (Zv · (period − 1))` for period > 1, else 0
      - `eps_i = 0.3 · Zv − alpha_i`
      - `beta_sigma_i = Σ_{j∈heavy_nbrs} (0.5 if |eps_i−eps_j|≤0.3 else 0.75)`
      - `beta_ns_i = Σ_{bonds} f·y` where f=2 for triple bonds else 1;
        y=2 for aromatic, 1.5 if |Δeps|>0.3, else 1.0; f·y=0 for SINGLE bonds
      - `beta_d_i = 0.5` if non-ring, non-aromatic, with lone pairs, adjacent
        to aromatic atom; else 0
      - `gamma_i = alpha_i / (beta_sigma_i + beta_ns_i + beta_d_i)`, NaN if
        denominator = 0
    - `ETA_beta_s = Σ beta_sigma_i / 2` (halved to avoid double-counting)
    - `ETA_beta_ns_d = Σ beta_d_i` (per-atom, no double-counting)
    - `ETA_beta_ns = Σ (beta_ns_i / 2 + beta_d_i)`
    - `ETA_beta = ETA_beta_s + ETA_beta_ns`
    - `ETA_eta = Σ_{i<j, r>0} sqrt(gamma_i · gamma_j) / r` (heavy-atom D)
    - Reference mol (all-C, all-SINGLE): same topology → same distance matrix;
      only gamma values change. `ETA_eta_R` uses original D with reference-mol
      gammas (matching Mordred's `DistanceMatrix(explicit_hydrogens=False)`
      dependency shared between original and reference `EtaCompositeIndex`).
    - `ETA_dAlpha_A = max((alpha − alpha_R) / n, 0)`
    - `ETA_dAlpha_B = max((alpha_R − alpha) / n, 0)`
    - `ETA_eta_B = eta_NL − eta_RL` where `eta_NL = 1.0` (n=2), or
      `sqrt(2) + 0.5·(n−3)` (n≥3); NaN for n≤1
    - `ETA_eta_BR = eta_NL − eta_RL + 0.086 · ring_count`
    - Epsilon variants:
      - type 1: mean eps over all atoms with explicit H
      - type 2: mean eps over heavy atoms only
      - type 3: mean eps over all atoms in reference mol with H (all-C, all-SINGLE + AddHs)
      - type 4: mean eps over all atoms in saturated mol (C-C bonds → SINGLE; other bonds kept; + AddHs)
      - type 5: mean eps over heavy atoms + H bonded to heteroatoms
    - `ETA_dEpsilon_A = eps1 − eps3`, B: eps1−eps4, C: eps3−eps4, D: eps2−eps5
    - `ETA_psi_1 = alpha / (n · eps2)`;
      `ETA_dPsi_A = max(0.714 − psi1, 0)`, `ETA_dPsi_B = max(psi1 − 0.714, 0)`
    - `ETA_shape_p/y/x = Σ(alpha_i for degree=1/3/4) / alpha_total`

    4,160 oracle cells (65 × 64), zero mismatches.

    **MDE (19):** `MDEC-{v1}{v2}` (10), `MDEN-{v1}{v2}` (6), `MDEO-{v1}{v2}` (3).
    Algorithm: collect atom-index pairs (i,j) with i<j where both atoms share
    the same element and matching heavy-atom degrees v1/v2; `n / exp(Σlog(D[i,j]) / n)`,
    NaN when n=0.

## Remaining Families (future expansion)

33 Mordred 2D descriptors remain unsupported (`tests/unsupported_mordred.json`).
All are 2D-computable in principle. Continue the validated-increment pattern.
Current unsupported list is in `tests/unsupported_mordred.json`.

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
