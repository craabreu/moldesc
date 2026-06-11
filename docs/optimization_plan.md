# Optimization Plan

Performance plan for the RDKit-only descriptor calculator. This file records both
the **completed** optimizations (with the reasoning and the rejected experiments,
so they are not re-attempted) and a **backlog** of future candidates ranked by the
current profiling baseline. The compatibility contract is fixed (every supported
descriptor must still match the Mordred oracle within tolerance), so every
optimization here is **behavior-preserving**: it changes how a value is computed,
never what the value is.

## Strategy

1. **Profile before optimizing.** Use `scripts/profile_rdkit_mordred_like.py
   --grouped` to rank descriptor families, but read it by *total elapsed per
   group*, not throughput-per-descriptor — a small family can still dominate
   wall-clock. Then drill in with `cProfile` (sort by `tottime`) to separate
   Python-level cost from the C++/RDKit calls underneath.
2. **Prefer provably-equivalent algebraic shortcuts over micro-tuning.** The
   biggest wins come from recognizing that an expensive intermediate result is
   never actually needed, or that a quantity has a cheaper closed form. Verify
   the equivalence numerically against the current implementation across the
   full panel and every parameter (order/lag/property) *before* replacing code,
   then rely on the oracle test suite as the safety net.
3. **Batch into numpy when the per-item work is uniform.** Stacking many small
   matrices or vectors and issuing one vectorized call amortizes Python/dispatch
   overhead.
4. **Resist speculative optimization.** If a fast-path is not demonstrably
   helping on the validation panel, the added branch/complexity is not worth it.
   Keep the simplest code that hits the measured target.

## Architecture

The calculator builds one `_DescriptorContext` per molecule and computes each
descriptor through a flat `name → function` map. Expensive intermediates
(distance/adjacency matrices, the explicit-hydrogen molecule, per-property atom
vectors) are `cached_property` values on that shared context, so work is done
once per molecule and reused across every family that needs it. Profiling
confirmed there is no meaningful cross-family redundancy to harvest and no global
restructuring worth doing: the grouped benchmark *appears* to show redundant
setup only because it rebuilds a fresh context per group, which the real entry
point never does.

One consequence of the cached-family design is that **subset evaluation is
naturally cheap** — touching only one family's descriptors never triggers any
other family's computation. `calc_rdkit_mordred_like_2d(mol, names=...)` exposes
this: requesting only the 56 chi descriptors costs ~23% of a full run, because
autocorrelation, BCUT, walk counts, etc. are simply never accessed. This is the
right lever for callers that need a handful of descriptors, and it is a thin
wrapper over the existing dispatch rather than a new code path.

## Completed optimizations

### Information content: build the BFS tree once, collect all orders in one walk

The IC family (`IC`/`TIC`/`SIC`/`BIC`/`CIC`/`MIC`/`ZMIC`, orders 0-5) groups atoms
by a canonical code of their depth-*m* BFS subtree. Mordred's `BFSTree` builds that
code per atom *per order* — for order *m* it resets and expands the tree *m* times,
then walks it to regenerate every root-to-leaf trail. The original port mirrored
this exactly, so each atom's tree was rebuilt and re-walked six times.

The key invariant: `_ic_expand_tree` only ever *adds* a deeper level and never
mutates a shallower one, so the depth-5 tree truncated to depth *m* equals the
independently built order-*m* tree. The tree is therefore built **once per root**
(to depth 5) and every order's trails are gathered in a **single DFS**
(`_ic_root_codes` / `_ic_walk` in `rdkit_mordred_like.py`): each tree position is a
frontier leaf of the order-*d* tree where *d* is its depth, so it contributes its
trail to `order_trails[d]`; a *natural* leaf (no children in the full tree) stays a
leaf at every deeper order, so it also contributes to orders *d+1..5*. Per-atom
masses and atomic numbers are precomputed into lists instead of re-fetched from the
molecule inside each group's Shannon-entropy sum.

A Morgan/Weisfeiler-Lehman integer-refinement rewrite was **prototyped and
rejected**: Mordred's BFS tree duplicates an atom across cycles (two sibling
frontier nodes can both claim the same unvisited neighbor) and updates its
`visited` set asymmetrically *during* sibling expansion, so the codes are
genuinely order-dependent rather than a clean neighborhood hash. The WL prototype
matched on acyclic molecules but diverged on every ring system in the panel
(furan, pyrrole, spiro, norbornane). Recorded so it is not re-attempted — the
speedup has to come from removing redundant work in the *same* traversal, not from
swapping the traversal.

Verified bit-exact against the current implementation (0 / 4285 per-(root,order)
codes) and against the Mordred oracle (0 / 2730 cells) before landing. IC isolated
time dropped ~47.7 → ~24.5 ms/panel-pass.

### chi/Kier: one shared subgraph enumeration

`kier_values` (`Kier1`/`Kier2`/`Kier3`) needs the order-2 and order-3 *path*
subgraph counts, which it obtained by re-running `FindAllSubgraphsOfLengthN` plus
the chi subgraph classifier for orders 2-3 — duplicating the orders 2-7 sweep
`chi_values` already performs. The enumeration is now a shared
`chi_subgraph_accumulators` cached property returning `(sums, nan_flag, counts)`;
`chi_values` consumes the sums/counts for its connectivity indices and
`kier_values` reads `counts[("path", 2)]` / `counts[("path", 3)]`, so the subgraph
sweep runs once per molecule.

Both changes together: full panel ~255 → ~228 ms/pass (commit `cae5a5f`).

### Autocorrelation: fused per-lag numpy pass

All autocorrelation properties share the same per-lag graph-distance work (the
adjacency-at-distance matrix and its quadratic forms). Rather than re-running
that O(n²) work once per property, `autocorrelation_values` stacks every
property vector into a `P×n` matrix and computes all properties together per
lag. The Geary numerator uses the identity
`Σᵢⱼ Bᵢⱼ(wᵢ−wⱼ)² = 2(w²·deg) − 4·ATS` so it never materializes the pairwise
difference matrix. Adding properties now scales with a matrix width, not with
repeated graph traversals — autocorrelation is the cheapest compute-intensive
group per descriptor as a result.

Two follow-up passes cut this family a further ~40% (1.51s → 0.91s), once
profiling showed the molecules are small enough that *Python-level* overhead, not
the numpy work, dominated:

- **Build the property vectors in one atom walk.** The element-table properties
  (mass, vdW volume, electronegativities, polarizability, ionization potential)
  previously re-iterated the molecule once per table. `autocorrelation_property_vectors`
  now reads atomic numbers once and indexes each table from that list, and it
  computes the sigma/valence counts a single time and reuses them to derive the
  intrinsic state (`_intrinsic_state_from`) instead of recomputing both inside it.
  Vector building dropped ~56% (0.68s → 0.30s).
- **Precompute descriptor keys; vectorize MATS/GATS.** The per-lag loop formatted
  ~600 f-string keys and called `float()` per element every molecule, even though
  the keys depend only on `(family, order, suffix)`. They are now built once into
  `_AC_KEYS`; each family's row is produced as a single numpy array (MATS/GATS via
  masked `np.where`, preserving the earlier scalar arithmetic order so values stay
  bit-for-bit identical) and assigned with one `zip(keys, row.tolist())`. Verified
  bit-identical to the prior output across the full panel (65 molecules × 624
  descriptors) before landing.

### BCUT: one batched eigensolve instead of twelve

The 12 Burden matrices for the 12 BCUT properties differ only in their diagonal;
the off-diagonal bond-weight structure is identical. `bcut_values` builds the
shared off-diagonal once, broadcasts it into a `(12, n, n)` batch, swaps each
property's diagonal in, and calls `np.linalg.eig` **once** on the batch instead
of 12 times. NaN-diagonal rows (atoms outside a property table) are replaced
with a `0.0` sentinel before the solve — LAPACK raises `LinAlgError` on NaN
input — and overwritten with NaN after sorting. This lifted the BCUT group from
~85k to ~138k descriptors/sec.

### Matrix-spectral: batch every matrix through one eigensolve

The spectral family computes eigenvalue/eigenvector aggregates over 11 matrices
per molecule — adjacency, distance, detour, and eight Barysz variants. The
first implementation built and diagonalized each matrix separately; profiling
ranked it the single most expensive group (2.64s on the panel), with the cost
split across three Python-level bottlenecks: eight separate Floyd-Warshall
loops, eleven separate `eigh` dispatches, and per-matrix aggregate extraction.

All three share the BCUT lesson — every matrix for one molecule is the same
`n×n` shape, so they batch:

- **One batched Floyd-Warshall for the eight Barysz matrices.** They differ only
  in their edge weights (`C²/(P[i]·P[j]·π_ij)`), not their structure, so
  `_compute_barysz_matrices` stacks them into a `(k, n, n)` tensor and relaxes
  all of them against each pivot at once: `w[:, :, p:p+1] + w[:, p:p+1, :]`.
- **One batched `eigh`.** Adjacency, distance, detour, and the Barysz batch are
  stacked into a single `(k, n, n)` array and diagonalized in one call;
  `eigh` returns `(k, n)` ascending eigenvalues and `(k, n, n)` eigenvectors,
  so the leading pair is always at index `n−1`. (NaN Barysz matrices — atoms
  outside a property table — are dropped from the batch *before* the solve and
  filled with NaN directly, since LAPACK rejects NaN input.)
- **Vectorized aggregates.** SpAbs/SpMax/SpDiam/SpAD/SpMAD/LogEE and VE/VR are
  computed across the whole batch with array ops; even VR1's per-bond product is
  a single fancy-indexed `lead[:, ai] * lead[:, aj]` with a masked `np.where`
  for the "≤ 0 ⇒ NaN" rule. The per-suffix result dict is then filled from
  `.tolist()` conversions (one C-level conversion per aggregate) rather than
  extracting ~12·k numpy scalars individually — the same per-element-`float()`
  cost the autocorrelation pass eliminated.

Verified bit-for-bit identical (max abs diff 0.0) to the per-matrix code across
the full panel before landing. The group dropped from 2.64s to ~1.43s (−46%).
The remaining cost is the Barysz bond-weight construction (a scatter assignment
that does not vectorize cleanly) and the detour-matrix DFS (inherently
recursive); both were left alone as the simplest code meeting the target.

### Path counts: "simple path ⟺ N+1 distinct atoms"

`MPC*`/`piPC*` count *self-avoiding* paths and sum their bond-order products.
The original code, for each path returned by `FindAllPathsOfLengthN`, did two
passes: first reconstruct the ordered atom sequence (branchy bond-to-atom
connectivity logic), then re-iterate building a `set` to test that no atom
repeats. cProfile showed the cost was almost entirely this Python per-path work,
not the C++ enumeration — and the ordered sequence was only ever used for the
distinctness check.

Key invariant: RDKit paths never repeat a *bond*, so a path of `N` bonds is an
atom-simple path **iff it touches exactly `N+1` distinct atoms** (revisiting an
atom on a bond-distinct trail would close a cycle, raising the edge/vertex
count). That collapses both passes into one: union the bond endpoints into a
set, multiply the bond orders, and accept the path when `len(atoms) == order+1`.
Verified to reproduce identical counts and π-weights across the full panel × all
orders before replacing the code. Result: ~11% faster on the group and ~40 lines
of fragile connectivity reconstruction deleted.

A ring-free fast-path (every bond-trail in an acyclic graph is automatically
simple, so the set check can be skipped) was prototyped and **rejected**: the
validation panel is ring-heavy, so cyclic molecules dominate the path
enumeration and the branch showed no measurable benefit while adding complexity.
Recorded here so it is not re-attempted without a workload that actually
warrants it.

A later pass confirmed this group is now **at its floor**, so no further change
was made. Measured split of the ~0.88s group: ~0.46s is the RDKit C++
enumeration itself (orders 1–10 are ten separate `FindAllPathsOfLengthN` calls
per molecule — there is no range-returning API to fold them into one), leaving
only ~0.42s of Python post-processing. Two attempts to cut that lost:

- **Vectorize per order.** Within one order every path has exactly `order`
  bonds, so the paths form a rectangular `(P, order)` matrix and
  `np.array(paths)` builds it in C. The whole order then reduces to fancy-indexed
  endpoint lookup + per-row sort for the distinct-atom count and a `prod` for the
  π-weight. It is **slower** (1.08s vs 0.88s): path sets per (molecule, order)
  are small, so numpy's fixed per-call overhead dominates. A clear case of the
  "batch only when per-item work is uniform *and large*" caveat.
- **Tighter Python loop.** Replacing the two `set.add` calls with one
  `set.update(pair)` per bond, and the inline π-multiply with `math.prod`, both
  regressed slightly. The existing two-`add`-plus-inline-multiply loop is already
  the fastest pure-Python form.

Recorded so the vectorization is not re-attempted on this (small-molecule)
workload.

### Chi connectivity: classify subgraphs by degree, not by DFS

The Kier-Hall chi family (`Xp-*`/`Xc-*`/`Xch-*`/`Xpc-*`) classifies every
connected bond-subgraph from `FindAllSubgraphsOfLengthN` as path, cluster,
path-cluster, or chain. The original port mirrored Mordred's `Chi.py`: build a
neighbor adjacency dict per subgraph and run a recursive DFS to detect a
back-edge (chain) and collect the degree set. cProfile showed the DFS plus
neighbor-dict construction were the dominant cost (~1.5s of the group, with
~1M recursive `_dfs` calls).

The same distinct-node-count invariant from path counts applies: a subgraph has
exactly `order` bonds (`E == order`) and is connected, so it contains a cycle
**iff `E ≥ V`**, i.e. iff its distinct-node count is `≤ order`. That is the chain
test — no traversal. For the remaining (tree) case, `V == order + 1`, the type
is fixed by the in-subgraph degree multiset alone: all degrees ≤ 2 → path; a
degree-2 node present alongside a branch point → path-cluster; only terminal and
branch nodes (no degree 2) → cluster. So one pass counting bond-endpoint degrees
replaces the whole DFS. Verified to reproduce Mordred's classifier on all 3601
panel subgraphs (zero mismatches) before replacing the code. Result: the chi
group dropped ~33% (1.82s → 1.22s) and the full panel ~10%.

## Future optimization plan

### Baseline to work from

Current grouped profile (`scripts/profile_rdkit_mordred_like.py --grouped
--repeat 60`, full panel, ~13.5s total), ranked by elapsed:

| group | elapsed (s) | state |
| --- | --- | --- |
| scalar_aliases | 2.89 | opaque grab-bag — split it (A) |
| chi | 1.63 | optimized once; likely near floor (D) |
| information_content | 1.50 | just halved; residual (E) |
| eta | 1.46 | **never optimized** (C) |
| autocorrelation | 1.27 | heavily optimized; near floor |
| path_counts | 1.06 | optimized; confirmed at floor |
| spectral | 0.84 | optimized (batched eigensolve) |

Read the grouped numbers with the Architecture caveat in mind: each group is timed
with a *fresh* `_DescriptorContext`, so `scalar_aliases` is inflated by setup
(Gasteiger charges, distance/detour matrices) that the single-context entry point
amortizes. Confirm a candidate against the real `calc_rdkit_mordred_like_2d` run
before trusting a group delta.

Every candidate below must follow the same discipline as the completed work:
**numerically diff the new output against the current implementation across the
full panel and every parameter (order/lag/property) before replacing code**, then
run `conda run -n open3d python -m pytest tests/ -q` as the final oracle gate.
Record any rejected experiment with its workload caveat.

### A. Enabler — split the `scalar_aliases` profiling bucket ✓ done

Added `carbon_types`, `drug_likeness`, `graph_scalars` named groups to
`scripts/profile_rdkit_mordred_like.py`. What remains in `scalar_aliases` is only
the genuinely trivial counts (atom/bond/ring counts, `MW`, `TopoPSA`, etc.).

### B. Cache Crippen logP/MR and HBA/HBD; reuse in drug-likeness filters ✓ done

Added `mol_log_p`, `mol_mr`, `num_hba`, `num_hbd` cached properties on
`_DescriptorContext` and rewired `SMR`, `SLogP`, `nHBAcc`, `nHBDon`, `_lipinski`,
`_ghose_filter` to read them. Eliminated up to 3× redundant `MolLogP` and 2×
redundant `MolMR`/`CalcNumHBA`/`CalcNumHBD` calls per molecule.

### C. ETA family — cache reference/saturated mols ✓ done

`cProfile` showed `_eta_build_reference_mol` and `_eta_atom_properties` were each
called **2× per molecule** (once directly, once via `_eta_reference_mol_with_h`
which unconditionally re-called `_eta_build_reference_mol`). Additionally, the
kekulized heavy-atom mol was rebuilt inside `eta_values` on every call even though
`eta_values` is already a `cached_property`. Fixed by adding four `cached_property`
entries on `_DescriptorContext`: `_eta_kekulized_mol`, `_eta_reference_mol`,
`_eta_reference_mol_with_h` (reuses `_eta_reference_mol` via `Chem.AddHs`),
`_eta_saturated_mol`. ETA isolated: ~29 → ~22 ms/pass (~25%); full panel ~228 →
~225 ms. All 77 tests pass.

### D. chi accumulation (1.63s) — measure first, likely near floor

Past the DFS→degree-count rewrite. The remaining cost is the C++
`FindAllSubgraphsOfLengthN` sweep (one call per order, orders 2-7 — there is no
range-returning API to fold them) plus the per-subgraph Python product loop in
`chi_subgraph_accumulators`. A measurement pass could test vectorizing the
endpoint-value products per order (subgraphs of one order form a rectangular
node-index matrix), but heed the path-count lesson already recorded above: numpy
*lost* there because per-(molecule, order) sets are small and fixed call overhead
dominated. Treat as "measure, expect to leave alone."

### E. IC residual (1.50s, post-halving) — lower priority

After the single-walk rewrite the residual is dominated by building and sorting the
trail tuples. A possible further step: intern each trail to an integer id per order
so the grouping step compares ints instead of nested tuples. Confirm with
`cProfile` that the sort/compare actually dominates before attempting — it may
already be at a reasonable floor, in which case record that and stop.
