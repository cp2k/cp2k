# K-Points

Periodic Quickstep calculations can sample the Brillouin zone through the `&DFT%KPOINTS` section.
The most common setup is a Monkhorst-Pack mesh:

```text
&DFT
  &KPOINTS
    SCHEME MONKHORST-PACK 6 6 6
    WAVEFUNCTIONS COMPLEX
  &END KPOINTS
&END DFT
```

`SCHEME GAMMA` uses only the Gamma point. `SCHEME MONKHORST-PACK` and `SCHEME MACDONALD` generate
regular meshes. `SCHEME GENERAL` uses explicitly listed k-points:

```text
&KPOINTS
  SCHEME GENERAL
  KPOINT 0.0 0.0 0.0 1.0
  KPOINT 0.5 0.0 0.0 1.0
&END KPOINTS
```

K-points are given in reciprocal lattice-vector coordinates by default (`UNITS B_VECTOR`). Cartesian
coordinates can be selected with `UNITS CART_BOHR` or `UNITS CART_ANGSTROM`.

## Symmetry Reduction

Atomic symmetry can reduce the number of k-points that have to be solved explicitly:

```text
&KPOINTS
  SCHEME MONKHORST-PACK 8 8 8
  SYMMETRY ON
  SYMMETRY_BACKEND SPGLIB
  WAVEFUNCTIONS COMPLEX
&END KPOINTS
```

`SYMMETRY_BACKEND` selects the backend that provides and applies atom and k-point symmetry
operations:

- `K290`: the established CP2K default.
- `SPGLIB`: use the space-group operations returned by SPGLIB, including fractional translations.

`SYMMETRY_REDUCTION_METHOD` selects the method used to build the irreducible k-point set. If
`SYMMETRY_BACKEND` is set explicitly and `SYMMETRY_REDUCTION_METHOD` is omitted, the reduction
method follows the backend. `SYMMETRY_REDUCTION_METHOD SPGLIB` with `SYMMETRY_BACKEND K290` can be
used as a comparison mode: SPGLIB proposes the k-point orbits, while only mappings represented by
the K290 backend are used for the actual transformations.

`INVERSION_SYMMETRY_ONLY ON` restricts the reduction to time-reversal/inversion symmetry.
`FULL_GRID ON` disables symmetry reduction while still generating the regular mesh.

For general atomic k-point symmetry, use complex wavefunctions. Real wavefunctions are only valid
for Gamma and special k-points with real Bloch phases.

## Explicit K-Point Sets

`SCHEME GENERAL` can be used with symmetry reduction when the explicit k-point set is equally
weighted and closed under the requested symmetry operations:

```text
&KPOINTS
  SCHEME GENERAL
  SYMMETRY ON
  SYMMETRY_BACKEND SPGLIB
  KPOINT -0.25 -0.25 -0.25 1.0
  KPOINT -0.25 -0.25  0.25 1.0
  KPOINT -0.25  0.25 -0.25 1.0
  KPOINT -0.25  0.25  0.25 1.0
  KPOINT  0.25 -0.25 -0.25 1.0
  KPOINT  0.25 -0.25  0.25 1.0
  KPOINT  0.25  0.25 -0.25 1.0
  KPOINT  0.25  0.25  0.25 1.0
&END KPOINTS
```

The same closure requirement applies to K290, SPGLIB, and mixed SPGLIB-reduction/K290-backend
setups. For band paths or other intentionally nonuniform explicit lists, keep `SYMMETRY OFF`,
because the order and weights of the points are part of the requested property.

## Cell Convention

CP2K's k-point machinery assumes the standard CP2K cell convention: vector `A` lies along the
Cartesian X axis and vector `B` lies in the XY plane. Using `ABC` together with `ALPHA_BETA_GAMMA`,
or reading a CIF file, lets CP2K construct the cell in that convention from orientation-independent
lattice parameters.

## Moving Geometries

For moving geometries (`GEO_OPT`, `CELL_OPT`, `MD`, and related run types), atomic k-point symmetry
is rebuilt from the current cell and coordinates instead of reusing operations from the initial
geometry. This applies to regular Monkhorst-Pack/MacDonald meshes and to closed `GENERAL` k-point
sets. If the current geometry no longer supports a symmetry operation, the reduced set changes
accordingly or the explicit `GENERAL` set is rejected if it is no longer symmetry-closed.

`KEEP_SPACE_GROUP T` is the safest way to combine geometry optimization with full atomic k-point
symmetry. `KEEP_SYMMETRY T` constrains the cell metric, but does not by itself keep atoms on
space-group-related positions.

## Wannier90

The Wannier90 interface can either keep its historical k-point path controlled by
`DFT%PRINT%WANNIER90%MP_GRID`, or use the full SCF k-point mesh via
`DFT%PRINT%WANNIER90%KPOINTS_SOURCE SCF`. The SCF source supports Gamma, Monkhorst-Pack, MacDonald,
and explicit `SCHEME GENERAL` k-point meshes. If the SCF calculation used K290 or SPGLIB symmetry
reduction, the corresponding unreduced mesh is regenerated or recovered for the Wannier90 export,
because Wannier90 expects the full mesh and its nearest-neighbour connectivity. The Wannier90
interface still diagonalizes the full export mesh instead of reconstructing full-mesh MO
coefficients from the irreducible SCF orbitals.

## Related Pages

- <https://www.cp2k.org/faq:kpoints>
- [](hartree-fock/ri_kpoints)
