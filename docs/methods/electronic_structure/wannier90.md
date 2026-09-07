# Wannier90 interface

CP2K can generate the input and matrix files required by [Wannier90](https://wannier.org/) through
[&WANNIER90](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90). The interface is experimental. It prepares
the k-point mesh, eigenvalues, and overlap matrices from a periodic Quickstep calculation; the
subsequent construction and use of Wannier functions are performed by Wannier90.

Wannier90 requires a complete, uniformly weighted k-point mesh with its nearest-neighbour
connectivity. A high-symmetry band path is not a suitable input mesh. See [](../dft/k-points) for
k-point sampling and convergence.

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

`EXCLUDE_BANDS` refers to the original, one-based MO indices, including the additional bands. CP2K
removes these states consistently from the eigenvalues, overlap matrices, and library projections
before passing data to Wannier90. The retained bands are renumbered from one in ascending original
order. Repeated indices are ignored, and out-of-range indices are rejected. A companion
`SEED_NAME_band_indices.dat` file lists the exported band index and its original MO index. The
generated files are already filtered: do not add the same exclusion list to the `.win` file. Library
localization rejects an exclusion that cuts a degenerate band group inside the outer energy window
(adjacent eigenvalues within `1.e-8` hartree), since selecting individual states there would depend
on the arbitrary MO gauge. Retain or remove the entire group instead.

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
checks detect trial sets that lose rank at individual k-points. CP2K then tries other QR anchors and
column exchanges that reduce the total rank defect, without changing the band space. If no full-rank
set is found, the calculation stops. Numerically tied QR columns are selected in their original
candidate order, so roundoff in equivalent MO gauges does not arbitrarily select a different
physical trial.

`INITIAL_PROJECTIONS AO_HYBRID` forms tetrahedral s/p combinations from every radial s-shell and
p-shell pair on each atom before selecting the projections. AOs not used in a hybrid group remain
available, including d and higher angular momenta. These candidates are constructed from the full AO
projection bank, not by rotating an already selected `.amn` matrix. They can provide different,
better localized starting functions for degenerate valence and semicore spaces. This is an
alternative start, not an automatic guarantee of the lowest minimum: compare converged spreads and
centres with the AO-based start and converge the localization settings for the desired observable.
Equal spreads alone do not establish equality of individual centres in a nearly flat localization
minimum.

`INITIAL_PROJECTIONS AUTO` compares the AO start, the all-atom hybrid start and mixed starts with
hybrids on one eligible atomic kind at a time. Each family is minimized independently with the
requested `NUM_CG_STEPS` and with pure steepest descent (`NUM_CG_STEPS 0`), avoiding duplicate runs
when zero was already requested. This checks sensitivity to the optimization path as well as to the
projections. The default `NUM_CG_STEPS 5` matches Wannier90's conjugate-gradient reset interval.

For more than one Wannier function, AUTO also tests two fixed unitary mixtures of the best original
trial set. These apply successive real Givens rotations of +30 and -30 degrees to adjacent trial
columns, identically on the entire mesh. Both start from the same original projections, not from the
localized orbitals. The band subspace, projection singular values, raw overlaps and eigenvalues are
unchanged. Each mixed start is minimized with the requested reset interval and with
`MAX(20, NUM_CG_STEPS)`, skipping duplicate settings. This adds at most four minimizations and
reduces the observed sensitivity to the minimization path in degenerate band spaces. It does not use
random noise or alter the convergence criteria. A one-function calculation skips these mixtures,
which would have no effect.

`NUM_PRINT_CYCLES 10` limits the iteration log volume without changing the optimization or its
convergence test. Set it to one to retain every iteration; initial and final states and the explicit
convergence report are always written. All trials use the same `NUM_ITER`, `CONV_WINDOW`, and
`CONV_TOL`; an unconverged or rank-deficient trial cannot win. CP2K stops if no admissible trial
converges.

Trials run sequentially in fresh Wannier90 instances sharing copies of the original input matrices,
not the overlaps modified by an earlier minimization. CP2K retains the optimized state with the
smallest converged total spread. This is a finite candidate comparison, not a guarantee of a
globally optimal or unique set of Wannier functions.

With `WRITE_INPUTS T`, each executed trial writes `SEED_NAME.trial-N.win` and
`SEED_NAME.trial-N.amn`; its log is `SEED_NAME.trial-N.library.wout`. The raw `.mmn` and `.eig`
matrices are shared by all trials. For an external comparison of a particular trial, use its `.win`
and `.amn` together with those raw matrices under a common seed name in a separate directory. The
canonical `.win`, `.amn`, and `.library.wout` are replaced by the winning trial's settings, initial
projections, and log. Writing these files does not reset the optimized library matrices.

`INITIAL_PROJECTIONS IDENTITY` is available for diagnostics but retains the arbitrary MO gauge. It
can converge to different local minima for symmetry-reconstructed and newly diagonalized MOs. The
library handles collinear spin channels independently. Even with AO trials, convergence to a
stationary localization result does not guarantee the global minimum of the spread functional.

CP2K reports both the total spread and its gauge-invariant contribution. The latter depends on the
selected band subspace, but not on the unitary rotations used to localize it. It therefore helps
distinguish a change of subspace from different local minima within the same subspace.

### Spin-polarized calculations

[SPIN_CHANNEL](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.SPIN_CHANNEL) selects `BOTH` (the
default), `ALPHA`, or `BETA`. For two-spin calculations, CP2K writes independent files with the seed
suffixes `_up` and `_down`. Each channel uses its own eigenvalues, MO coefficients, projections,
Hamiltonian and library instance. The ordinary file exporter uses the same separation; it does not
concatenate the two spin channels into a single `.eig` or `.mmn` file. With one spin channel, the
original seed name is unchanged, and selecting `BETA` is an input error.

Specify `WANNIER_FUNCTIONS` once to use the same count in both channels, or twice to select the
alpha and beta counts in that order. The order is independent of `SPIN_CHANNEL`, so the following
exports only the beta channel with three Wannier functions:

```text
&WANNIER90 ON
  SEED_NAME magnetic_crystal
  KPOINTS_SOURCE SCF
  SPIN_CHANNEL BETA
  WANNIER_FUNCTIONS 5
  WANNIER_FUNCTIONS 3
  &LIBRARY ON
    INITIAL_PROJECTIONS AUTO
  &END LIBRARY
&END WANNIER90
```

Band exclusions and library settings, including energy windows, apply to every selected channel.
They must be valid for each channel individually. One spin can require disentanglement while the
other has equal band and Wannier counts. The printed total and invariant spreads belong to the
preceding spin channel, not to a sum over spins. This is collinear spin support, not a spinor or
spin-orbit-coupled Wannierization.

### Disentanglement

When the number of exported bands is larger than `WANNIER_FUNCTIONS`, the library first minimizes
the gauge-invariant spread to select a connected subspace, then localizes that subspace. Use
`ADDED_MOS` to include extra bands in the export. The AO projection matrix may be rectangular; the
same initial candidates and raw band-space overlaps are used for the external reference.

`DIS_NUM_ITER`, `DIS_CONV_WINDOW`, `DIS_CONV_TOL`, and `DIS_MIX_RATIO` control this first stage.
Both disentanglement and localization must converge for an AUTO candidate to be admissible.
`DIS_WIN_MIN` and `DIS_WIN_MAX` optionally restrict the outer energy window. `DIS_FROZ_MAX` enables
a frozen inner window; `DIS_FROZ_MIN` is optional. Energies are in eV unless an explicit CP2K unit
is provided, and refer to the exported eigenvalues, not to a shifted Fermi-energy zero. These
windows must leave at least `WANNIER_FUNCTIONS` bands at every k-point and cannot freeze more than
that many states. With equal band and Wannier counts, energy windows are rejected. AO candidates are
projected into the allowed outer window before their selection. If the frozen states already fill
the entire target subspace, only those states enter the projection selection; irrelevant unoccupied
bands cannot steer the choice of initial orbitals in this case.

### Optimized outputs

`WRITE_HR T` writes the converged real-space Hamiltonian (`SEED_NAME_hr.dat`, in eV) and the native
Wigner-Seitz displacement file (`SEED_NAME_wsvec.dat`) through Wannier90's public postprocessing
API. `WRITE_U_MATRICES T` writes `SEED_NAME_u.mat`; with disentanglement it also writes
`SEED_NAME_u_dis.mat` and their product `SEED_NAME_v.mat`. These options are independent of
`WRITE_INPUTS` and default to false. In AUTO mode, only the selected converged state produces these
canonical outputs; later losing trials cannot overwrite them. When `WRITE_U_MATRICES T` is used for
a complete-space rerun with the same seed, CP2K removes obsolete `_u_dis.mat` and `_v.mat` files
from a previous disentangled calculation. Outputs whose write options are disabled are not
refreshed.

The matrix formats and conventions match external Wannier90. In particular, Wannier90 v4 stores the
rows of `u_dis` and `v` in **packed outer-window order**, with zero padding beyond the window. Do
not multiply these rows directly by the original SCF coefficients when the outer window excludes
lower bands. CP2K additionally writes `SEED_NAME_band_map.dat`: after a comment and the number of
k-points and exported bands, each row lists the one-based k-point index, exported band index, matrix
row (zero for an excluded band), and eigenvalue in eV. With no disentanglement, this is the identity
band mapping for `u`. No `.chk` restart file is produced by these options. When `EXCLUDE_BANDS` is
used, compose this packed-window mapping with `_band_indices.dat` to obtain the original MO indices.

Real-space interpolation requires a suitable, converged k-point mesh. Wannier90 warns when the mesh
lacks Gamma; generating an `_hr.dat` file alone does not validate off-mesh interpolation.

## Generated files

In the ordinary file-export mode, `SEED_NAME silicon` produces the following Wannier90 files:

- `silicon.win`: a starting Wannier90 input file containing the cell, atomic positions, exported
  band count, and k-point mesh;
- `silicon.mmn`: overlap matrices between neighbouring k-points;
- `silicon.eig`: eigenvalues for the exported bands;
- `silicon.amn`: an identity projection matrix, only when `USE_BLOCH_PHASES T` is used; and
- `silicon_band_indices.dat`: the mapping from exported bands to original MO indices, including the
  identity mapping when no bands are excluded.

With `LIBRARY ON`, the additional input files require `WRITE_INPUTS T`, and `.amn` contains the
selected `INITIAL_PROJECTIONS` rather than necessarily an identity matrix. Optimized outputs are
controlled separately by `WRITE_HR` and `WRITE_U_MATRICES`.

CP2K regenerates the enabled files when the calculation is run. Preserve a separate copy of a
completed Wannier90 input file, or add project-specific settings after the CP2K export has finished.

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

The atom/AO transformation uses the same lattice-periodic Bloch gauge as the SCF matrices. Atom cell
shifts contribute phases; folding a k-point by a reciprocal lattice vector does not add an extra
atom-position phase to the SCF coefficients. This is distinct from the optional export gauge
described below.

A finite integration grid can weakly break crystal symmetry even when the AO overlap metric is
preserved. The reference validation can then reject the reconstructed eigenvalues and select
full-mesh diagonalization. Converge both `CUTOFF` and `REL_CUTOFF` as well as the SCF threshold when
testing quantitative reconstruction; SCF convergence alone does not remove grid errors.

## Bloch phases and projections

[USE_BLOCH_PHASES](#CP2K_INPUT.FORCE_EVAL.DFT.PRINT.WANNIER90.USE_BLOCH_PHASES) applies the CP2K
Bloch-phase gauge. In the ordinary file-export mode, it also writes identity `.amn` projections and
requires `WANNIER_FUNCTIONS` to equal the number of exported bands. Disentanglement in that mode
still requires explicit Wannier90 projections.

With `LIBRARY ON`, `INITIAL_PROJECTIONS` determines the trial matrix instead. The physical AO and
hybrid projections use the same gauge as the overlaps, including when there are more bands than
Wannier functions. The equal-count restriction therefore does not apply to the library path.

## Limitations

The CP2K Wannier90 interface is experimental. In particular:

- use a complete k-point mesh rather than a band path;
- verify the convergence of the SCF and export meshes for the target quantity;
- inspect CP2K output when exporting from a symmetry-reduced SCF mesh, since CP2K may reconstruct
  the missing orbitals or fall back to full-mesh diagonalization; and
- consult the Wannier90 documentation for localization, projection, disentanglement, interpolation,
  and post-processing settings that are not controlled by CP2K.
