# General Cell Canonicalization

This note describes the requirements for canonicalizing general input cell vectors in CP2K. It is a
design checklist for future implementation work, not a user-facing promise that all cases are
handled already.

CP2K accepts cell vectors through `CELL%A`, `CELL%B`, `CELL%C`, `CELL%ABC`, `CELL%ALPHA_BETA_GAMMA`,
and several external cell or coordinate file formats. Internally, many algorithms assume the
conventional orientation where vector `A` lies on the Cartesian x axis and vector `B` lies in the xy
plane. If a user provides a rotated but otherwise equivalent cell, merely rewriting `CELL` into that
canonical orientation is not enough. Every input quantity whose meaning depends on the Cartesian
frame must be transformed consistently.

## Invariants

A complete implementation must preserve the following invariants.

- Scalar observables such as total energy, pressure, and electronic charge are unchanged by
  canonicalization.
- Cartesian vectors such as positions, velocities, forces, electric fields, and constraint
  directions are rotated into the canonical frame.
- Fractional lattice coordinates and reciprocal `B_VECTOR` coordinates keep their numeric values.
- Cartesian reciprocal-space vectors, including explicit k-points and band-path points, are
  transformed as reciprocal vectors.
- Rank-2 Cartesian tensors such as stress are transformed consistently with the chosen tensor
  convention.
- Restart files written after canonicalization restart into the same physical state, and existing
  restart files are either interpreted unambiguously or rejected with a clear message.

## Transform

Let `H_in` be the cell matrix that maps fractional coordinates to Cartesian coordinates in the input
frame, and let `H_can` be the canonical cell matrix for the same lattice:

```text
r_in  = H_in  s
r_can = H_can s
```

The Cartesian transform from the input frame to the canonical frame is therefore

```text
T = H_can inv(H_in)
r_can = T r_in
```

Cartesian real-space vectors use `T`. Cartesian reciprocal-space vectors must preserve phases,
`q_in . r_in = q_can . r_can`, and therefore use

```text
q_can = transpose(inv(T)) q_in
```

For a pure canonicalization that preserves the metric, `T` is an orthogonal rotation or
rotoinversion, but the implementation should still use the general formulas above.

## Policy

Canonicalization must be explicit and auditable. A safe staged policy is:

- Add a single transform state after the final input cell has been read, after `MULTIPLE_UNIT_CELL`
  and requested cell symmetry have been applied.
- Keep both `H_in` and `H_can` available to all input readers that may parse cell-dependent data.
- Do not silently canonicalize when an enabled input section contains a Cartesian quantity that has
  no implemented transform.
- Initially expose the behavior behind an opt-in keyword or an internal guarded mode. Automatic
  canonicalization should only become the default after the test matrix below is covered.
- Print one concise summary whenever canonicalization is applied, including the original cell,
  canonical cell, determinant sign, and maximum metric deviation.

## Input Coverage

The following inputs must be covered before canonicalization can be considered complete.

| Area                      | Required handling                                                                                                        | Minimal regression coverage                                                                |
| ------------------------- | ------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------ |
| `CELL`                    | Preserve lengths, angles, volume, periodicity, symmetry, and `CELL_REF` semantics.                                       | Canonical and rotated cells give identical scalar observables.                             |
| `COORD`                   | Keep scaled coordinates unchanged; transform Cartesian coordinates with `T`.                                             | Paired scaled and Cartesian coordinates produce the same energy and rotated forces.        |
| `SHELL_COORD`             | Apply the same scaled or Cartesian policy as atom coordinates.                                                           | A core-shell force-field input is invariant under a rotated cell.                          |
| `CORE_COORD`              | Apply the same scaled or Cartesian policy as atom coordinates.                                                           | A core-shell restart or static calculation is invariant.                                   |
| External coordinate files | Treat CP2K, XYZ, extended XYZ, PDB, CIF, and related file formats according to their native coordinate convention.       | Each supported format has at least one rotated-cell pair, or is rejected when unsupported. |
| `KPOINTS GENERAL`         | Keep `B_VECTOR` values unchanged; transform `CART_BOHR` and `CART_ANGSTROM` values as reciprocal vectors.                | Explicit general k-points match the equivalent canonical input.                            |
| Band paths                | Apply the same reciprocal-space policy as `KPOINTS GENERAL`.                                                             | Band endpoints and eigenvalues are invariant for `B_VECTOR` and Cartesian paths.           |
| Electric and other fields | Transform Cartesian field vectors with `T`; define the behavior of Berry-phase and periodic-field directions explicitly. | Forces and polarization are invariant for fields along nontrivial directions.              |
| Constraints               | Transform Cartesian constraint vectors, normals, fixed directions, and collective-variable directions.                   | Constrained optimizations or single-step MD runs are invariant.                            |
| Velocities                | Transform atomic, shell, and core velocities with `T`.                                                                   | Kinetic energy and first MD step are invariant.                                            |
| Stress                    | Transform stress tensors with the documented Cartesian tensor convention.                                                | Stress components rotate correctly; scalar pressure is invariant.                          |
| Restarts                  | Write canonicalized restarts with enough metadata to avoid frame ambiguity.                                              | Restarted runs match uninterrupted runs and old incompatible restarts fail clearly.        |

## Implementation Checkpoints

1. Introduce a small canonicalization state that stores `H_in`, `H_can`, `T`, the reciprocal
   transform, and whether the transformation changes handedness.
1. Thread that state through topology, coordinate, k-point, band-structure, field, constraint,
   velocity, stress, and restart setup code.
1. Apply transforms at input-boundary routines, so physics kernels continue to see one consistent
   internal representation.
1. Add guard rails that detect unhandled Cartesian input sections and either disable
   canonicalization before data are consumed or abort before the calculation starts.
1. Add paired regression tests for every row in the coverage table.
1. Only after the full matrix is covered, consider enabling canonicalization automatically for
   general cells.

## Non-goals

This work does not replace local bug fixes that avoid unsafe symmetry paths for noncanonical cells.
Those fixes remain valuable while the broader input canonicalization infrastructure is designed and
tested.
