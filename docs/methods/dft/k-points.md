# K-Points

Periodic Quickstep calculations approximate Brillouin-zone integrals with a finite, weighted set of
k-points. The sampling set used for the self-consistent-field (SCF) calculation is defined in
[&DFT%KPOINTS](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS).

This page introduces k-point sampling in CP2K, explains the available sampling schemes, and
summarizes the experimental atomic-symmetry reduction path. A band-structure path is a different
object from the SCF integration mesh; it is defined through
[&DFT%PRINT%BAND_STRUCTURE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.BAND_STRUCTURE), not through an
arbitrary `SCHEME GENERAL` list.

## Why k-points are needed

### Brillouin-zone integration

For a periodic system, Bloch's theorem labels the one-electron states by a crystal momentum
$\mathbf{k}$ in the Brillouin zone. Quantities such as the electronic density, total energy, and
occupations contain integrals over this zone. In a numerical calculation, CP2K replaces such an
integral by a weighted sum over a finite set of k-points,

$$
  \frac{1}{\Omega_\mathrm{BZ}}\int_\mathrm{BZ} f(\mathbf{k})\,d\mathbf{k}
  \approx \sum_{\mathbf{k}} w_{\mathbf{k}} f(\mathbf{k}),
$$

where $w_{\mathbf{k}}$ are normalized k-point weights.

A Gamma-only calculation samples only $\mathbf{k}=0$. It is often appropriate for isolated systems,
large supercells, or other cases where the Brillouin zone is sufficiently small. Smaller primitive
cells, metals, and systems with strongly dispersive bands usually need a converged k-point mesh.

### K-points in CP2K

For every sampled k-point, CP2K solves a k-dependent Kohn--Sham problem and combines the resulting
quantities with the k-point weights. The `&KPOINTS` section therefore describes the sampling used
during the SCF calculation, rather than a post-processing path through selected high-symmetry
points.

Omitting `&KPOINTS` gives the usual Gamma-only calculation (`SCHEME NONE`). `SCHEME GAMMA` instead
creates an explicit one-point k-point set at Gamma. The physical sampling is the same in most cases,
but the two inputs use different implementation paths.

Complex wavefunctions are the default for k-point calculations; see
[WAVEFUNCTIONS](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS.WAVEFUNCTIONS). Real wavefunctions are only valid
for Gamma and special k-points whose Bloch phases can be represented as real. Use complex
wavefunctions for a general mesh or for atomic k-point symmetry reduction.

## Choosing and converging a mesh

Although rules of thumb can provide a useful starting point, reliable results require converging the
k-point mesh for the specific system and property of interest. Total energies, forces, stresses,
metallic occupations, density of states, and band edges can converge at different rates. Increase
the mesh density until the relevant quantity no longer changes at the accuracy required for the
calculation.

For slabs, wires, and other low-dimensional systems, sample the periodic directions and normally use
one k-point in a non-periodic or vacuum direction. Enlarging a real-space supercell reduces the
Brillouin zone and can reduce the required k-point density, but does not by itself remove
finite-size effects.

```{important}
Electronic smearing and DOS broadening do not replace k-point convergence. In particular, a smooth
DOS obtained from a sparse mesh may still be physically unconverged. See
[](../electronic_structure/dos.md#broadening-k-points-and-gaps).
```

For a conventional band structure, first converge the SCF calculation on an appropriate integration
mesh. Then use [&BAND_STRUCTURE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.BAND_STRUCTURE) and
[&KPOINT_SET](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.BAND_STRUCTURE.KPOINT_SET) to define the path.

## Sampling schemes

[SCHEME](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS.SCHEME) selects one of the schemes described below.
Regular meshes are evaluated as full meshes by default: atomic symmetry reduction is only requested
when [SYMMETRY](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS.SYMMETRY) is explicitly enabled.

### Gamma-only sampling

For the conventional Gamma-only calculation, omit `&KPOINTS` entirely:

```text
&DFT
  ...
&END DFT
```

An explicit Gamma-point set can be requested when a k-point calculation path is required by the
workflow:

```text
&DFT
  &KPOINTS
    SCHEME GAMMA
  &END KPOINTS
&END DFT
```

### Monkhorst--Pack meshes

A Monkhorst--Pack mesh is the usual regular sampling scheme for periodic calculations:

```text
&DFT
  &KPOINTS
    SCHEME MONKHORST-PACK 6 6 6
  &END KPOINTS
&END DFT
```

The three integers specify the mesh dimensions along the reciprocal lattice vectors. Use
[GAMMA_CENTERED](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS.GAMMA_CENTERED) to generate a Gamma-centered
variant:

```text
&KPOINTS
  SCHEME MONKHORST-PACK 6 6 6
  GAMMA_CENTERED T
&END KPOINTS
```

Gamma centering is supported for Monkhorst--Pack meshes. It is most useful when an even number of
subdivisions is used and the mesh is required to include Gamma point.

### MacDonald meshes

A MacDonald mesh specifies both the mesh dimensions and an explicit shift:

```text
&KPOINTS
  SCHEME MACDONALD 4 4 4 0.25 0.25 0.25
&END KPOINTS
```

The first three values define the mesh dimensions; the final three define the shift.

### Explicit k-point sets

`SCHEME GENERAL` accepts an explicitly supplied weighted set of k-points:

```text
&KPOINTS
  SCHEME GENERAL
  KPOINT 0.0 0.0 0.0 1.0
  KPOINT 0.5 0.0 0.0 1.0
&END KPOINTS
```

Each [KPOINT](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS.KPOINT) line contains three coordinates and one
weight. CP2K normalizes the supplied weights internally.

By default, [UNITS](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS.UNITS) is `B_VECTOR`, so the coordinates are
expressed in reciprocal-lattice-vector coordinates. Cartesian coordinates can instead be selected
with `CART_BOHR` or `CART_ANGSTROM`; their units are $2\pi/\mathrm{Bohr}$ and $2\pi/\mathrm{\AA}$,
respectively.

```{note}
`SCHEME GENERAL` defines an integration set for the SCF calculation. It is not the usual interface
for a high-symmetry band path. Use `&DFT%PRINT%BAND_STRUCTURE` for that purpose.
```

### Parallelization over k-points

[PARALLEL_GROUP_SIZE](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS.PARALLEL_GROUP_SIZE) controls how MPI
processes are grouped for a k-point calculation. Its value is the number of MPI processes assigned
to one k-point group. The group size must divide the total number of MPI processes, and the
resulting number of groups must divide the number of k-points.

The default `-1` selects the smallest valid number of processes per group. `0` uses all processes
for each k-point, while a positive value requests that exact group size. This setting is a
parallelization choice and does not change the physical k-point mesh.

## Related workflows and output

The k-point mesh used for SCF can also be used by several electronic-structure workflows. Their
additional requirements are documented separately:

- [](../electronic_structure/dos) explains DOS and PDOS output. DOS/PDOS commonly need a denser mesh
  than a geometry optimization.
- [](../electronic_structure/molecular_orbitals.md#k-point-mo-output-mokp) documents the `.mokp`
  k-point MO output. Molden output is not available for k-point calculations.
- [](hartree-fock/ri_kpoints) documents RI-HFX with k-point sampling and includes a band-structure
  example.
- [&WANNIER90](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90) is an experimental interface. With
  [KPOINTS_SOURCE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.KPOINTS_SOURCE) set to `SCF`, it can
  use the full SCF mesh from `&DFT%KPOINTS`; Wannier90 export requires a complete mesh.
- [](../optimization/geometry_and_cell_opt) describes geometry and cell optimization. Their
  interaction with experimental atomic k-point symmetry is discussed below.

## K-point symmetry reduction

For regular Monkhorst--Pack and MacDonald meshes, CP2K distinguishes two levels of k-point
reduction:

1. **k-space inversion (time-reversal) reduction**, which pairs $\mathbf{k}$ and $-\mathbf{k}$ and
   is used by default for regular meshes; and
1. **atomic (space-group) symmetry reduction**, which uses additional operations that map the
   current periodic structure onto itself.

The [SYMMETRY](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS.SYMMETRY) keyword controls the second level. It is
off by default; this does **not** disable the default time-reversal reduction for regular meshes.

### Time-reversal reduction

For a regular Monkhorst--Pack or MacDonald mesh, CP2K normally combines inversion-related
$\mathbf{k}$ and $-\mathbf{k}$ points. This is the standard reduction path when both `SYMMETRY F`
(the default) and `FULL_GRID F` (the default) are used:

```text
&KPOINTS
  SCHEME MONKHORST-PACK 8 8 8
&END KPOINTS
```

The reduction can also be requested explicitly with
[INVERSION_SYMMETRY_ONLY](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS.INVERSION_SYMMETRY_ONLY). This is
useful when `SYMMETRY T` is present in a shared input template, but only time-reversal reduction is
desired:

```text
&KPOINTS
  SCHEME MONKHORST-PACK 8 8 8
  SYMMETRY T
  INVERSION_SYMMETRY_ONLY T
&END KPOINTS
```

`SCHEME GENERAL` preserves the supplied list by default. When inversion-only reduction is requested
for an explicit list, each $\mathbf{k}$/$-\mathbf{k}$ pair must be present with equal weights.

To calculate every point of a regular mesh explicitly, disable atomic symmetry and request the full
mesh:

```text
&KPOINTS
  SCHEME MONKHORST-PACK 8 8 8
  SYMMETRY F
  FULL_GRID T
&END KPOINTS
```

```{note}
For Monkhorst--Pack and MacDonald meshes, `FULL_GRID T` together with `SYMMETRY T` disables atomic
symmetry reduction but retains k-space inversion (time-reversal) reduction. Use both `SYMMETRY F`
and `FULL_GRID T` to obtain a strict full-mesh reference calculation.
```

### Atomic (space-group) symmetry reduction

```{warning}
Atomic k-point symmetry reduction is experimental. Validate the energy, forces, stress, and any
other target property against an equivalent full-mesh calculation before using it for production.
```

Atomic symmetry reduction further groups regular-grid k-points that are related by operations of the
current atomic structure. Enable it with `SYMMETRY T`:

```text
&KPOINTS
  SCHEME MONKHORST-PACK 8 8 8
  SYMMETRY T
  WAVEFUNCTIONS COMPLEX
&END KPOINTS
```

This combines the default time-reversal reduction with the additional atomic symmetry operations.
Complex wavefunctions are required for general atomic symmetry operations with nontrivial Bloch
phases.

#### Compatible sampling sets

Atomic symmetry reduction applies to regular Monkhorst--Pack and MacDonald meshes. It can also be
used with `SCHEME GENERAL`, provided that all explicit weights are equal and that the complete set
is closed under every requested symmetry operation. A nonuniform `GENERAL` list, including a band
path, should keep `SYMMETRY F`.

#### Cell requirements

For regular Monkhorst--Pack and MacDonald meshes, full atomic reduction currently requires a cell
matrix in the standard CP2K lower-triangular convention. If a full atomic reduction is requested for
a non-orthogonal cell or for a cell matrix outside that convention, CP2K warns and falls back to
[INVERSION_SYMMETRY_ONLY](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS.INVERSION_SYMMETRY_ONLY).

Defining the cell through [ABC](#CP2K_INPUT.FORCE_EVAL.SUBSYS.CELL.ABC) and
[ALPHA_BETA_GAMMA](#CP2K_INPUT.FORCE_EVAL.SUBSYS.CELL.ALPHA_BETA_GAMMA), or reading a suitable CIF
structure, lets CP2K construct the standard cell orientation from orientation-independent lattice
parameters.

#### Symmetry backends

[K290](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS.SYMMETRY_BACKEND) is the established default backend. The
optional `SPGLIB` backend uses symmetry operations returned by spglib, including fractional
translations:

```text
&KPOINTS
  SCHEME MONKHORST-PACK 8 8 8
  SYMMETRY T
  SYMMETRY_BACKEND SPGLIB
  WAVEFUNCTIONS COMPLEX
&END KPOINTS
```

This option requires CP2K to be built with
[spglib](../../technologies/libraries.md#spglib-crystal-symmetries-tools). If `SYMMETRY_BACKEND` is
specified and
[SYMMETRY_REDUCTION_METHOD](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS.SYMMETRY_REDUCTION_METHOD) is
omitted, the reduction method follows the selected backend.

`SYMMETRY_REDUCTION_METHOD SPGLIB` together with `SYMMETRY_BACKEND K290` is a comparison mode:
SPGLIB proposes the k-point orbits, while K290 operations are used for the actual transformations.
It is primarily useful for validation and development rather than as a default production setup.

#### Moving geometries

For `GEO_OPT`, `CELL_OPT`, molecular dynamics, and related calculations, CP2K determines atomic
k-point symmetry from the current cell and coordinates rather than assuming that the initial
operations remain valid. The irreducible k-point set can consequently change as the geometry
evolves. A `SCHEME GENERAL` list must remain symmetry-closed at every step or it is rejected.

When a geometry or cell optimization is intended to preserve the full space group, use the relevant
[KEEP_SPACE_GROUP](#CP2K_INPUT.MOTION.GEO_OPT.KEEP_SPACE_GROUP) setting. For `CELL_OPT`, see also
the discussion of [KEEP_SYMMETRY](#CP2K_INPUT.MOTION.CELL_OPT.KEEP_SYMMETRY) in
[](../optimization/geometry_and_cell_opt.md#constraints-cell-degrees-of-freedom-and-symmetry).

#### Validation and troubleshooting

For every new system or workflow, compare an atomic-symmetry-reduced calculation with the strict
full-mesh reference:

```text
&KPOINTS
  SCHEME MONKHORST-PACK 8 8 8
  SYMMETRY F
  FULL_GRID T
  WAVEFUNCTIONS COMPLEX
&END KPOINTS
```

[&DFT%PRINT%KPOINTS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.KPOINTS) prints k-point information and is
useful for checking the generated set. For further diagnostics, use
[VERBOSE](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS.VERBOSE) and, where necessary, adjust
[EPS_SYMMETRY](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS.EPS_SYMMETRY). The
[DEBUG_FULL_KPOINT_SYMMETRY](#CP2K_INPUT.FORCE_EVAL.DFT.KPOINTS.DEBUG_FULL_KPOINT_SYMMETRY) option
is intended for expert finite-difference debugging.
