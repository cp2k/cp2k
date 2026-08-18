# Spatial distribution functions

`spatial_distribution.x` computes orientationally resolved, three-dimensional spatial distribution
functions (SDFs) from CP2K XYZ trajectories. It is a generalized Fortran implementation of the
water-specific `rot_anal_H2O.c` workflow originally written by Thomas D. Kuehne: neither the cell
geometry nor an O-H-H atom order is hard coded.

For every selected reference triplet, the program applies periodic boundary conditions, constructs a
local orthonormal frame, rotates target-atom displacements into that frame, and accumulates them on
a Cartesian grid. It writes one Gaussian cube per target selection and an XYZ file containing the
average reference geometry in the local frame.

## Output and scope

For every `TARGET`, the program writes `<prefix>-<name>.cube`. It also writes
`<prefix>-reference.xyz` with the averaged reference triplet. These files contain all accumulated
three-dimensional information.

The separate radial check table produced by the legacy water program is intentionally not repeated.
It was a validation projection of the same histogram rather than an independent result, and a
spherical average or coordination integral can be derived from the generated cube without coupling
this tool to a particular radial binning convention.

## Build and run

The tool is standalone and requires a Fortran 2008 compiler. It does not have to be linked to CP2K.

```console
make
./spatial_distribution.x water.in
```

See [`water.in.example`](water.in.example) for a complete bulk-water setup. Coordinates and cells
are read in angstrom; Gaussian cube coordinates are written in bohr, as required by the cube format.

## Input

The input is line based. Empty lines and text following `#` or `!` are ignored. Paths are resolved
relative to the input file. A path containing spaces can be enclosed in single or double quotes.

| Keyword                                      | Meaning                                                                                                                                           |
| -------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| `TRAJECTORY file`                            | CP2K multi-frame XYZ trajectory.                                                                                                                  |
| `CELL_FILE file`                             | CP2K `&MOTION / &PRINT / &CELL` output. Records are matched to the `i =` step in each XYZ comment.                                                |
| `CELL ax ay az bx by bz cx cy cz`            | Static cell instead of `CELL_FILE`; A, B, and C are cell vectors in angstrom.                                                                     |
| `REFERENCE origin arm1 arm2`                 | One-based atom indices defining one local frame. Repeat as needed.                                                                                |
| `REFERENCE_PATTERN o a1 a2 os a1s a2s count` | Generate `count` references with independent strides for the origin and both arms. This covers both interleaved and species-grouped trajectories. |
| `FRAME_MODE BISECTOR`                        | Put the normalized arm bisector on local z and arm2-minus-arm1 on local x. This reproduces the water convention.                                  |
| `FRAME_MODE TRIAD`                           | Put origin-to-arm1 on local z and the perpendicular part of origin-to-arm2 on local x. This works for a general non-collinear triplet.            |
| `TARGET name labels`                         | Create an SDF for exact XYZ labels. `labels` may be comma-separated; `*` selects all atoms.                                                       |
| `BOUNDS xmin xmax ymin ymax zmin zmax`       | Histogram edges in angstrom. Default: -5 to 5 on every axis.                                                                                      |
| `GRID nx ny nz`                              | Number of histogram bins. Default: 81 on every axis.                                                                                              |
| `CUTOFF radius`                              | Optional spherical cutoff in angstrom. Zero disables it.                                                                                          |
| `PERIODIC axes`                              | Periodic lattice directions: `XYZ`, any subset, or `NONE`. Default: `XYZ`.                                                                        |
| `EXCLUDE_REFERENCE logical`                  | Exclude the current origin and its two arms from target histograms. Default: true.                                                                |
| `FIRST_FRAME n`, `LAST_FRAME n`, `STRIDE n`  | Select trajectory frames. `LAST_FRAME 0` means all remaining frames.                                                                              |
| `OUTPUT_PREFIX path`                         | Prefix for cube and reference-geometry files. Default: `sdf`.                                                                                     |

Use either `CELL_FILE`, `CELL`, or neither. With neither, each frame must be in CP2K's extended-XYZ
format and contain a `Lattice="..."` field.

The trajectory must keep the same atom ordering in every frame. References are listed explicitly
because an XYZ trajectory contains no bond topology. This makes the tool applicable to water,
solvated ions, molecular liquids, surfaces, and other environments without guessing molecules from
distances.

CP2K can write atomic-kind names instead of element symbols with
`MOTION / PRINT / TRAJECTORY / PRINT_ATOM_KIND T`. This is useful when a target must distinguish
chemically different atoms of the same element.

For a bent molecule such as water, swapping `arm1` and `arm2` reverses the local x and y directions.
For linear reference geometries, use `TRIAD` with a third atom that is not collinear with the first
arm.

## Normalization

For target selection `s`, the output in voxel `v` is

```text
g_s(v) = N_s(v) / {Delta_V * sum_f [N_ref(f) * N_s(f) / V(f)]}.
```

Here, `N_s(v)` is the accumulated number of target observations, `N_ref(f)` is the number of
references, `N_s(f)` the number of selected target atoms, and `V(f)` the instantaneous cell volume.
The per-frame density weighting is important for NPT trajectories. The SDF is dimensionless and
approaches one in a homogeneous bulk region. Excluding atoms belonging to the current reference does
not change the bulk-density denominator, matching the usual RDF/SDF self-exclusion convention.

The minimum-image displacement is found by a bounded lattice search, so skewed triclinic cells and
partial periodicity do not rely on orthorhombic assumptions. Histogram values are written in the
standard cube ordering (z is the fastest index), with voxel centers rather than duplicated end
points on the bounds.

## Tests

```console
make test
```

The tests cover the water bisector convention with CP2K cell-history matching, rotation invariance,
and a non-water `TRIAD` reference in a triclinic cell.
