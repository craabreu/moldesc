# AGENTS.md

## Project goal

This repository implements a **RDKit-only subset of Mordred 2D descriptors**.

The objective is to compute as many Mordred-name 2D descriptors as possible using RDKit in production code, without importing `mordred` or `mordredcommunity` outside of tests, benchmark scripts, or development utilities.

`mordredcommunity` may be used as a reference/oracle in tests to validate numerical equivalence when Mordred returns numeric values and to document allowed RDKit behavior when Mordred returns missing values.

## Non-negotiable rules

- Production/library code must not import `mordred` or `mordredcommunity`.
- RDKit is the only chemistry backend allowed in production descriptor calculations.
- Do not claim a descriptor is supported unless it is covered by tests against `mordredcommunity`.
- Preserve Mordred descriptor names for compatible descriptors whenever possible.
- For renamed RDKit equivalents, document the mapping explicitly.
- Do not silently add approximate descriptor families to the Mordred-name output.
- Approximate or RDKit-native-only descriptors must go into a separate namespace, for example `RDKit_*`.
- Prefer small, testable changes over broad refactors.

## Definition of a supported descriptor

A descriptor is considered supported only if all conditions are true:

1. It is a 2D descriptor in Mordred/Mordred-community.
2. It can be computed from RDKit only.
3. It has either:
   - the same name as Mordred, or
   - a documented alias mapping from Mordred name to RDKit implementation.
4. It numerically matches Mordred-community within accepted tolerances whenever Mordred returns numeric values.
5. Any case where Mordred returns a missing value is explicitly documented in compatibility test expectations as either:
   - RDKit numeric value accepted because the RDKit implementation is meaningful and deterministic, or
   - RDKit `NaN` accepted because the descriptor is undefined in both implementations.
6. It is included in the locked supported-descriptor list.
7. It has unit-test coverage.

## Expected package structure

Use this structure unless the repository already has a different established layout:

```text
your_package/
  rdkit_mordred_like.py          # production RDKit-only descriptor code
  mordred_rdkit_registry.py      # supported names and alias mappings

tests/
  mordred_reference.py           # test-only Mordred-community oracle
  test_mordred_compat.py         # numerical compatibility tests
  smiles_panel.smi               # chemically diverse validation panel
  expected_supported.json        # locked supported descriptor names

scripts/
  discover_rdkit_mordred_overlap.py
  write_supported_list.py
```

If the actual package has a different name, adapt imports accordingly without changing the separation between production code and test-only Mordred reference code.

## Production implementation guidelines

Production descriptor code should expose a function similar to:

```python
def calc_rdkit_mordred_like_2d(mol):
    """Return Mordred-name 2D descriptors computed using RDKit only."""
```

Expected behavior:

- Input is an RDKit `Mol`.
- Output is a plain `dict[str, float | int]`.
- Keys are Mordred descriptor names.
- Values should be numeric whenever possible.
- Raise a clear `ValueError` for invalid inputs such as `mol is None`.
- Do not mutate the input molecule unless explicitly documented.

Recommended conservative descriptor groups:

- Molecular weight and atom-count descriptors:
  - `MW`
  - `AMW`
  - `nHeavyAtom`
  - `nHetero`
- Crippen descriptors:
  - `SLogP`
  - `SMR`
- H-bonding, polarity, and rotors:
  - `TopoPSA(NO)`
  - `TPSA`
  - `nHBAcc`
  - `nHBDon`
  - `nRot`
- Ring/topology descriptors:
  - `nSpiro`
  - `nBridgehead`
  - `FCSP3`
  - `BalabanJ`
  - `BertzCT`
- Surface-area families:
  - `LabuteASA`
  - `PEOE_VSA1` through `PEOE_VSA13`
  - `SMR_VSA1` through `SMR_VSA9`
  - `SlogP_VSA1` through `SlogP_VSA11`
  - `EState_VSA1` through `EState_VSA10`
  - `VSA_EState1` through `VSA_EState9`

These are starting points, not automatic guarantees. Keep only descriptors that pass compatibility tests.

## Alias registry

Maintain explicit mappings for renamed equivalents. Example:

```python
MORDRED_RDKIT_ALIASES = {
    "MW": "ExactMolWt",
    "SLogP": "MolLogP",
    "SMR": "MolMR",
    "nHBAcc": "NumHAcceptors",
    "nHBDon": "NumHDonors",
    "nRot": "NumRotatableBonds",
    "nHeavyAtom": "HeavyAtomCount",
    "FCSP3": "FractionCSP3",
}
```

When adding a new alias:

1. Add the mapping to the registry.
2. Add or update the RDKit-only implementation.
3. Add the descriptor to the supported list only after tests pass.
4. Document any tolerance or edge-case behavior.

## Testing requirements

Run tests after changes:

```bash
pytest -q
```

For compatibility debugging, use:

```bash
pytest -q -s tests/test_mordred_compat.py
```

Tests should verify:

- Production modules do not import `mordred` or `mordredcommunity`.
- Every supported descriptor is present in Mordred-community.
- Every supported descriptor is produced by the RDKit-only calculator.
- RDKit-only values match Mordred-community numeric values on the validation panel.
- Mordred missing values are handled only through explicit test expectations.
- The supported descriptor list remains stable unless intentionally updated.

Use `math.isclose` or `numpy.isclose` for numerical comparisons. Suggested starting tolerances:

```python
ABS_TOL = 1e-8
REL_TOL = 1e-6
```

Loosen tolerances only with a comment explaining why.

## Validation molecule panel

The validation panel should include chemically diverse molecules, not only trivial examples.

Include examples covering:

- aliphatic molecules
- aromatic molecules
- heteroaromatics
- acids
- amides
- amines
- charged species
- nitro groups
- heterocycles
- fused rings
- bridged or spiro systems when relevant

A minimal starter panel may include:

```text
CCO ethanol
CC(=O)O acetic_acid
CCN(CC)CC triethylamine
c1ccccc1 benzene
c1ccncc1 pyridine
O=C(O)c1ccccc1 benzoic_acid
CCOC(=O)c1ccccc1 ethyl_benzoate
C1CCCCC1 cyclohexane
CC(C)(C)c1ccc(O)cc1 tert_butyl_phenol
O=C(N)c1ccccc1 benzamide
CCS(=O)(=O)C sulfone
C[N+](C)(C)C tetramethylammonium
O=[N+]([O-])c1ccccc1 nitrobenzene
C1CCC2CCCCC2C1 decalin
C1OC1 epoxide
C1COCCO1 dioxane
```

Expand this panel when adding descriptors that are sensitive to rings, charges, hydrogen handling, aromaticity, or specific functional groups.

## Discovery workflow

When trying to add more descriptors:

1. Use a discovery script to compare Mordred-community and RDKit outputs across the validation panel.
2. Identify exact-name overlaps first.
3. Add renamed aliases only when there is a clear RDKit equivalent.
4. Add the RDKit implementation.
5. Run the compatibility tests.
6. Add passing descriptors to `expected_supported.json`.
7. Keep failing or approximate descriptors out of the supported Mordred-name output.

Do not add large descriptor families wholesale unless each descriptor has been validated.

## Handling approximate or partial descriptor families

Some RDKit descriptor families overlap conceptually with Mordred but may not be Mordred-equivalent. Examples include:

- `BCUT2D_*`
- `Chi*`
- `Kappa*`
- `fr_*` functional group counts
- `AUTOCORR2D`

These may be useful, but keep them separate unless they pass Mordred-name compatibility tests.

Preferred naming for unsupported RDKit-native descriptors:

```text
RDKit_BCUT2D_MWHI
RDKit_Chi0v
RDKit_fr_Ar_N
RDKit_AUTOCORR2D_001
```

Do not place these in the supported Mordred-name descriptor dictionary unless validated.

## Code style

- Keep descriptor code explicit and readable.
- Prefer registries and small helper functions over clever introspection.
- Avoid broad exception swallowing.
- If RDKit raises an exception for a descriptor, return a controlled missing value only if the project already has a missing-value policy.
- Keep descriptor names stable.
- Avoid changing molecule sanitization or hydrogen handling globally.
- Add comments for any descriptor where RDKit/Mordred equivalence is non-obvious.

## Missing values and errors

Use a consistent policy for descriptors that cannot be computed.

Preferred behavior:

- Invalid input molecule: raise `ValueError`.
- Descriptor calculation failure: return a numeric value when RDKit provides a meaningful deterministic value, return `float("nan")` for documented undefined cases, or raise a descriptor-specific error for unsupported/unknown failures.
- Do not return Mordred error objects from production code.
- Do not silently coerce non-numeric outputs to zero.

Document every Mordred-missing validation-panel case in the compatibility test expectations.

## Documentation requirements

When adding or removing supported descriptors, update the relevant documentation or generated report.

The report should distinguish:

1. Exact Mordred/RDKit name matches.
2. Mordred names implemented via RDKit aliases.
3. RDKit-native descriptors that are useful but not in the supported Mordred-name set.
4. Unsupported Mordred descriptors.

Avoid statements like "RDKit supports Mordred descriptor X" unless the test suite proves numeric compatibility where Mordred is numeric and documents behavior where Mordred is missing.

## Dependency policy

Production dependencies:

- `rdkit`
- standard library
- existing package dependencies already accepted by the project

Development/test dependencies may include:

- `mordredcommunity`
- `pytest`
- `pandas` or `numpy` for discovery/reporting scripts

Do not introduce new runtime dependencies unless necessary.

## Commands

Use the commands appropriate to the repository. If no project-specific commands exist, default to:

```bash
pytest -q
```

For formatting or linting, inspect the repository before adding assumptions. Do not introduce a new formatter or linter configuration unless explicitly requested.

## Pull request / patch expectations

Each patch should include:

- RDKit-only implementation change.
- Test or test update.
- Supported descriptor list update, if applicable.
- Short explanation of which descriptors were added and how they were validated.

Avoid unrelated cleanup in the same patch.
