# Wannier90 interface

CP2K can generate the input and matrix files required by [Wannier90](https://wannier.org/) through
[&WANNIER90](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90). The interface is experimental. It prepares
the k-point mesh, eigenvalues, and overlap matrices from a periodic Quickstep calculation; the
subsequent construction and use of Wannier functions are performed by Wannier90.

Wannier-function construction requires a complete, uniformly weighted k-point mesh with its
nearest-neighbour connectivity. A high-symmetry band path is not a suitable input mesh. See
[](../dft/k-points) for k-point sampling and convergence.

## Basic workflow

1. Run a periodic, diagonalization-based SCF calculation with a converged k-point mesh and enough
   bands for the intended Wannierization.
1. Enable `&DFT%PRINT%WANNIER90`. CP2K writes the Wannier90 input and data files.
1. Add the Wannier90 settings specific to the calculation, such as projections, disentanglement
   windows, or post-processing options, to the generated `.win` file, then run Wannier90 with the
   same seed name.

For example, the following exports a four-function Wannierization using the k-point mesh already
used by the SCF calculation:

```text
&FORCE_EVAL
  &DFT
    &KPOINTS
      SCHEME MONKHORST-PACK 6 6 6
    &END KPOINTS
    &PRINT
      &WANNIER90
        SEED_NAME silicon
        KPOINTS_SOURCE SCF
        WANNIER_FUNCTIONS 4
        ADDED_MOS 4
      &END WANNIER90
    &END PRINT
  &END DFT
&END FORCE_EVAL
```

[WANNIER_FUNCTIONS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.WANNIER_FUNCTIONS) sets the number of
Wannier functions. [ADDED_MOS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.ADDED_MOS) provides
additional bands for the export, and
[EXCLUDE_BANDS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.EXCLUDE_BANDS) can remove selected bands
from it. Choose the exported band window and the subsequent Wannier90 settings for the particular
material and target property.

## In-process localization

When CP2K is compiled with `CP2K_USE_WANNIER90=ON`, the optional
[LIBRARY](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.LIBRARY) subsection runs Wannier90 v4 directly,
without starting a separate executable. See [](../../technologies/libraries) for build instructions.
The ordinary file-export interface remains the default and does not require this dependency.

```text
&WANNIER90 ON
  KPOINTS_SOURCE SCF
  WANNIER_FUNCTIONS 1
  &LIBRARY ON
    NUM_ITER 1000
    CONV_WINDOW 5
    CONV_TOL 1.e-10
    WRITE_INPUTS T
  &END LIBRARY
&END WANNIER90
```

The library obtains the complete mesh through the existing export path, including reconstruction of
symmetry-reduced SCF orbitals. Neighbour connectivity is supplied by Wannier90. Overlap blocks are
distributed over the CP2K communicator by k-point. The reported centres and quadratic spreads are
converted to bohr and bohr squared; `CONV_TOL` uses Wannier90's angstrom-squared convention. The
optimizer log is written to `SEED_NAME.library.wout`. A calculation that reaches `NUM_ITER` without
satisfying the convergence criterion is rejected, even if the library reports no runtime error. In
Wannier90 4.0.2 this requires checking its explicit convergence report, since exhaustion of the
iteration budget is not reflected in the API return code.

`WRITE_INPUTS T` also writes matching `.win`, `.amn`, `.mmn`, and `.eig` files for comparison with
an external Wannier90 run. By default these additional files are not written. The library uses
`INITIAL_PROJECTIONS AO_SCDM` by default: CP2K projects the exported MOs onto its AO basis with the
complex overlap metric and selects one fixed set of AO trials over the complete mesh using pivoted
QR. This is an AO-based, SCDM-inspired construction, not sampling of the density matrix on a
real-space grid. The projections use the same Bloch gauge as the overlap matrices. Per-k-point rank
checks reject a trial set that does not span the target space.

`INITIAL_PROJECTIONS IDENTITY` is available for diagnostics but retains the arbitrary MO gauge. It
can converge to different local minima for symmetry-reconstructed and newly diagonalized MOs. The
initial library path requires a spin-unpolarized calculation with as many Wannier functions as
exported bands, without `EXCLUDE_BANDS`. Even with AO trials, convergence to a stationary
localization result does not guarantee the global minimum of the spread functional.

## Generated files

With `SEED_NAME silicon`, CP2K writes the following Wannier90 files:

- `silicon.win`: a starting Wannier90 input file containing the cell, atomic positions, exported
  band count, and k-point mesh;
- `silicon.mmn`: overlap matrices between neighbouring k-points;
- `silicon.eig`: eigenvalues for the exported bands; and
- `silicon.amn`: an identity projection matrix, only when `USE_BLOCH_PHASES T` is used.

CP2K regenerates these files when the calculation is run. Preserve a separate copy of a completed
Wannier90 input file, or add project-specific settings after the CP2K export has finished.

## Selecting the k-point source

[KPOINTS_SOURCE](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.KPOINTS_SOURCE) selects how CP2K builds
the export mesh.

### Use the SCF mesh

`KPOINTS_SOURCE SCF` uses the k-point mesh from `&DFT%KPOINTS`. It supports explicit Gamma,
Monkhorst--Pack, MacDonald, and equally weighted `GENERAL` meshes, and is the preferred choice when
the Wannier90 export should match the SCF calculation:

```text
&DFT
  &KPOINTS
    SCHEME MONKHORST-PACK 6 6 6
  &END KPOINTS
  &PRINT
    &WANNIER90
      KPOINTS_SOURCE SCF
      ...
    &END WANNIER90
  &END PRINT
&END DFT
```

`KPOINTS_SOURCE SCF` requires an active `&DFT%KPOINTS` section. For an explicit `GENERAL` mesh, the
points must have equal weights and must form a complete mesh from which Wannier90 connectivity can
be constructed. If CP2K cannot infer the mesh dimensions, set
[MP_GRID](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.MP_GRID) explicitly.

When the SCF calculation used k-point symmetry reduction, CP2K regenerates the corresponding full
mesh for the export. Wannier90 needs that full mesh even though the SCF calculation solved only its
irreducible subset.

### Use a separate Monkhorst--Pack mesh

`KPOINTS_SOURCE MP_GRID` is the historical default. It builds a full Gamma-centred uniform mesh from
[MP_GRID](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.MP_GRID), independently of the SCF k-point
setup. For a like-for-like comparison with `KPOINTS_SOURCE SCF`, use the same dimensions and
`GAMMA_CENTERED T` in `&DFT%KPOINTS`:

```text
&PRINT
  &WANNIER90
    KPOINTS_SOURCE MP_GRID
    MP_GRID 6 6 6
    ...
  &END WANNIER90
&END PRINT
```

This path performs the necessary full-mesh diagonalizations for the Wannier90 files. It can be used
when the SCF calculation is Gamma-only or when the export mesh intentionally differs from the SCF
mesh, but using a separately chosen mesh requires its own convergence assessment.

### Explicit overlap loops and external topology analysis

`KPOINTS_SOURCE NNKP` reads arbitrary fractional k-points and directed connections from `NNKP_FILE`.
The Wannier90-format file must contain `real_lattice`, `recip_lattice`, `kpoints` and `nnkpts`
blocks. The reciprocal-vector shift in each connection specifies the periodic closure; it must not
be discarded when a loop crosses a Brillouin-zone boundary.

```text
&WANNIER90
  KPOINTS_SOURCE NNKP
  NNKP_FILE loop.nnkp
  SEED_NAME loop
  WILSON_LOOP T
&END WANNIER90
```

CP2K keeps the converged SCF potential fixed while diagonalizing the requested post-SCF k-points.
The SCF mesh is independent of these overlap loops and must be converged separately. The exported
`.mmn` matrices contain the Gaussian-basis metric and periodic phase/image information:

`M(k,b) = C(k)^dagger O(k,b) C(k+b)`.

Plain Euclidean overlaps of AO coefficient vectors are not suitable. Directed cross-k overlaps use
ordered, nonsymmetric AO pair matrices. `SPIN_CHANNEL` selects a single collinear channel; UKS
channels are never concatenated in one `.mmn` file. In SOC mode, `EXCLUDE_BANDS` selects spinor
bands instead.

The NNKP/MMN file interface can be used by an external
[Z2Pack overlap-system adapter](https://z2pack.greschd.ch/en/latest/reference/other_systems.html).
Such an adapter supplies each requested closed loop in an NNKP file, runs CP2K with an unchanged SCF
setup, and returns the corresponding ordered MMN matrices to Z2Pack. The Python adapter is
maintained separately from CP2K; neither the native calculation nor CP2K's tests require Z2Pack. No
Wannier90 library or Wannier fit is needed for explicit overlap loops.

### Native Wilson loops and Z2 analysis

`WILSON_LOOP T` calculates Wilson eigenphases from the SVD-unitarized links. With
`KPOINTS_SOURCE WILSON`, CP2K generates a surface internally and doubles its longitudinal and
transverse resolution until the convergence checks pass or `WILSON_MAX_REFINEMENT` is reached.
`Z2 T` additionally calculates largest-gap crossing parity on a time-reversal half-plane.

The following illustrates an eight-electron SOC system with five scalar bands (ten spinors), of
which eight are retained. Adjust the band count and exclusions for the actual system:

```text
&WANNIER90
  KPOINTS_SOURCE WILSON
  SOC T
  Z2 T
  TIME_REVERSAL T
  EXCLUDE_BANDS 9 10
  WILSON_ORIGIN 0 0 0
  WILSON_DIRECTION 1 0 0
  WILSON_TRANSVERSE 0 0.5 0
  WILSON_MESH 4 3
  WILSON_MAX_REFINEMENT 2
&END WANNIER90
```

SOC uses restricted SCF followed by second-variational pseudopotential SOC, not self-consistent
noncollinear DFT. SOC-capable pseudopotentials are required. Converge the scalar unoccupied space
(`SCF/ADDED_MOS` and `WANNIER90/ADDED_MOS`), basis, cutoffs and SCF mesh independently.

`TIME_REVERSAL T` is a user assertion about the Hamiltonian, not an automatic symmetry proof. Native
Z2 requires an even, lowest-energy occupied spinor subspace with one state per electron and an
excluded conduction band. The supported native Z2 surfaces use two distinct reciprocal coordinate
axes, with winding one along the loop and one half transversely. A single plane provides one 2D
invariant, not all four strong/weak 3D indices. Arbitrary loops remain available for Wilson phases
and external analysis.

Checks cover singular links, sampled band gaps, boundary Kramers pairs, WCC changes under joint
refinement, adjacent-line movement and gap separation, and parity stability. A coarse mesh can still
miss a gap closing or rapid evolution between samples. Repeat with finer starting meshes and tighter
tolerances. Gapless graphene does not have a well-defined insulating Z2 invariant without specifying
and resolving a gap-opening Hamiltonian.

The `.wilson` file contains the loop index, minimum link singular value, and sorted hybrid Wannier
centres (WCC) in `[0,1)`. Printed Berry phases are in radians. Refinements overwrite the seed output
with the final mesh; the log retains diagnostics from all levels.

The regular CP2K test runner includes `topology_wilson_unittest` (known trivial/nontrivial BHZ
models, Chern models, gauge/reversal invariance and failure checks) and the four short helium/neon
inputs in `tests/QS/regtest-topology`. These are mathematical and smoke tests, not
material-convergence benchmarks. Larger DFT/SOC, adaptive-surface and MPI-scaling validations are
separate manual work.

### Native first Chern number and spectral gaps

`CHERN T` evaluates the determinant Wilson-phase winding on a full closed surface. For example:

```text
&WANNIER90
  KPOINTS_SOURCE WILSON
  CHERN T
  WILSON_DIRECTION 1 0 0
  WILSON_TRANSVERSE 0 1 0
  WILSON_MESH 16 17
  WILSON_MAX_REFINEMENT 3
  REQUIRE_GLOBAL_GAP T
&END WANNIER90
```

Set `EXCLUDE_BANDS` for the intended isolated subspace, with at least one computed band above it
when testing a gap. `CHERN` requires integer transverse winding, not a Z2 half-surface, and cannot
be combined with `Z2 T`. Scalar, single-collinear-channel and second-variational SOC states are
supported without imposing time reversal. The sign follows Z2Pack's increasing-transverse-coordinate
Wilson winding. No magnetic field, SOC term or other change to the Hamiltonian is implied.

Checks require endpoint WCC closure, resolved determinant-phase steps, stable integer winding, and
WCC convergence under joint refinement. The native unit test includes known C1=0,+1,-1 two-band
models, orientation reversal, direct-sum additivity and rejection of invalid surfaces. A helium
regression exercises `CHERN`, gap checking and state export through the regular CP2K test runner.

For a selected lowest-band prefix, output distinguishes the minimum sampled direct separation from
the sampled indirect gap `min(E[N+1]) - max(E[N])`. An isolated band bundle can have positive direct
separation but a negative indirect gap. Optional `REQUIRE_GLOBAL_GAP T` rejects a missing common
spectral interval above the prefix; it requires Wilson analysis and an excluded band above the
selected states. No point-dependent Fermi shifts are applied. Finite sampling cannot rule out a
missed bulk gap closing. Across separately self-consistent geometries, a common energy reference
must be justified before interpreting an indirect gap.

### Gaussian state snapshots for phason analysis

The native [`topology_phasons` postprocessor](phason-topology.md) evaluates physical cross-geometry
links, mixed-parameter C1 and four-parameter C2 from these snapshots, without Z2Pack or Python.
External consumers may use the same file interface.

For explicit `NNKP` or `WILSON` points, `STATE_EXPORT T` writes `SEED_NAME.topology`. The versioned
text output can be large and is disabled by default. It provides physical AO states and basis
metadata for cross-geometry analysis, not a Wannier fit. The export itself does not evaluate an
invariant; the native postprocessor or an external consumer performs that analysis. Existing `.mmn`
files describe cross-k links at a fixed geometry only.

Version 1 uses atomic units and contains, in order:

1. The header `CP2K_TOPOLOGY_STATE 1` and dimensions: atom count, AO count, selected state count,
   k-point count, spinor component count, total computed band count and collinear channel.
1. One-based selected band indices, the three direct lattice vectors, and three periodicity flags.
1. For each atom: index, kind index, number of Gaussian sets, AO count and canonical periodic
   centre. Each set stores its first atom-local AO index, spherical AO count, primitive count,
   Cartesian count, minimum angular momentum and screening radius. Each primitive stores its
   exponent and radius, followed by Cartesian powers and full spherical contraction coefficients.
1. For each point: index, fractional reciprocal coordinates, all computed eigenvalues, and selected
   complex AO coefficients in column-major order. Each complex value is a real/imaginary pair. SOC
   spinor components are stacked by AO; a scalar export contains only its selected channel.

The required moving-basis link is `C_a^dagger O_ab C_b`, with the cross-geometry Gaussian operator,
not a Euclidean coefficient overlap or an independently orthogonalized basis identification. A
physical phason family needs consistent cell, orbital/atom count, selected rank, spin convention and
self-consistent branch, together with explicit endpoint geometry and subspace sewing. Dropping and
adding atoms in a finite patch is not a fixed-rank cycle. External software must check metric
normalization, sampled isolation and link singular values, then demonstrate mesh convergence. The
snapshot format alone does not establish a topological invariant for a material.

## Reusing SCF orbitals

With `KPOINTS_SOURCE SCF`, [REUSE_SCF_MOS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.REUSE_SCF_MOS)
is enabled by default. CP2K reuses the SCF orbital coefficients directly when the SCF mesh is
already complete. It can also reconstruct some time-reversal and atomic-symmetry-related points from
a symmetry-reduced SCF mesh.

```{note}
SCF orbital reuse requires all relevant symmetry k-point data to be available within one k-point
parallel group. If the SCF calculation distributes k-points over multiple groups, CP2K cannot
currently collect the orbitals across those groups and instead falls back to a full-mesh
diagonalization for the Wannier90 files.

Set `PARALLEL_GROUP_SIZE 0` to keep all MPI processes in one k-point group and enable the reuse
path. This disables parallelization over k-points, however, so it is a compatibility setting for
SCF MO reuse rather than a general performance recommendation.
```

When the SCF k-point data are available, CP2K can reuse orbitals directly from a complete SCF mesh
or reconstruct missing points from a symmetry-reduced mesh. A symmetry-reconstructed export window
must contain complete degenerate subspaces. If the window cuts through a degenerate subspace, CP2K
falls back to a full-mesh diagonalization before writing the Wannier90 files. This preserves a
well-defined exported subspace.

[VALIDATE_REUSE_SCF_MOS](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.VALIDATE_REUSE_SCF_MOS) builds a
full-mesh reference and compares it with the reconstructed orbitals. It is expensive and intended
for development and diagnostic use, not routine production calculations.

## Bloch phases and projections

[USE_BLOCH_PHASES](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.USE_BLOCH_PHASES) applies the CP2K
Bloch-phase gauge and writes an identity `.amn` projection file. It is valid only when
`WANNIER_FUNCTIONS` equals the number of exported bands. Disentanglement calculations, or any case
with fewer Wannier functions than exported bands, still require explicit Wannier90 projections.

## Limitations

The CP2K Wannier90 interface is experimental. In particular:

- for Wannier-function construction, use a complete k-point mesh rather than a band path;
- verify the convergence of the SCF and export meshes for the target quantity;
- inspect CP2K output when exporting from a symmetry-reduced SCF mesh, since CP2K may reconstruct
  the missing orbitals or fall back to full-mesh diagonalization; and
- consult the Wannier90 documentation for localization, projection, disentanglement, interpolation,
  and post-processing settings that are not controlled by CP2K.
