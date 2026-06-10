# Optimization Notes

Performance notes for the RDKit-only descriptor calculator. The compatibility
contract is fixed (every supported descriptor must still match the Mordred
oracle within tolerance), so every optimization here is **behavior-preserving**:
it changes how a value is computed, never what the value is.

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

## Worked examples

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
