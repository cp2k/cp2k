# Modular trajectory analysis tools

This directory contains five standalone Fortran programs for post-processing CP2K trajectories:

- `trajana.x` performs geometry, hysteretic-event, RDF, VACF, hydrogen-bond, scattering, and
  collective-current analyses.
- `trajconvert.x` wraps, unwraps, converts, and reduces trajectories to group centers.
- `twopt.x` evaluates conventional two-phase thermodynamics directly from CP2K trajectories.
- `twopt3d.x` maps local solvent density and intermolecular 2PT entropy in three dimensions.
- `fresean.x` extracts frequency-selective anharmonic vibrational modes from molecular dynamics.

The programs share the same XYZ reader, triclinic-cell geometry, frame selection, atom selection,
and group definitions. They do not use CP2K source modules and can therefore be built independently.
The legacy `tools/trajana.f90` remains available while workflows are migrated.

## Build and test

```console
make
make test
```

`make check` performs a bounds-checked build with compiler warnings treated as errors and then runs
the reference tests. A different compiler or flags can be provided on the command line, for example:

```console
make FC=ifx FFLAGS="-O2"
```

## Trajectory and cell input

The trajectory must use XYZ format. Coordinates are interpreted as angstrom; velocity units are
preserved. The number and ordering of atoms must remain constant. Input and output default to
standard input and standard output and can be changed with `--input` and `--output`.

Cell information can be supplied in one of three ways:

```console
--cell "ax ay az bx by bz cx cy cz"
--cell-file PROJECT-1.cell
```

or through `Lattice="..."` in an extended-XYZ comment. A CP2K `.cell` file is read one data record
per XYZ frame and may contain a changing triclinic cell. The default periodicity is `XYZ`;
`--periodic XY`, `XZ`, `YZ`, `X`, `Y`, `Z`, or `NONE` is also accepted. RDF normalization currently
requires `XYZ`.

All modes accept `--first N`, `--last N`, and `--stride N`. Atom selections use one of:

```text
all
name:O,H
index:1,4-12,20
```

## Geometry analysis

The `geometry` mode is the migration path for the original `trajana` distance, angle, and torsion
calculations:

```console
./trajana.x geometry \
  --input PROJECT-pos-1.xyz \
  --actions examples/geometry.in \
  --cell-file PROJECT-1.cell \
  --output geometry.dat
```

The action-file syntax remains `DISTANCE i j`, `ANGLE i j k`, and `TORSION i j k l`; only the first
letter is significant. Geometry uses the minimum image when a cell is available. Distances are
written in angstrom and angles in degrees.

## Hysteretic events in geometric observables

The `events` mode detects threshold-crossing episodes in the same distance, angle, or torsion
actions used by `geometry`. Separate entry and exit values suppress rapid recrossings near a single
threshold. For example, C--O elongation episodes that start at 1.70 angstrom and end below 1.65
angstrom are obtained with:

```console
./trajana.x events \
  --input PROJECT-pos-1.xyz \
  --actions co-bonds.in \
  --cell-file PROJECT-1.cell \
  --dt-fs 0.5 \
  --direction above \
  --entry 1.70 \
  --exit 1.65 \
  --output co-events.dat \
  --summary co-event-summary.dat
```

For `--direction above`, the exit value must be smaller than the entry value. With
`--direction below`, the inequalities and required threshold order are reversed. All actions in one
run use the same direction and thresholds; run observables with different units or definitions
separately. Distances use angstrom, while angles and torsions use degrees. Periodic distances and
angles use the same triclinic minimum-image implementation as `geometry`.

The event output reports the first and last active sampled frames, sampled residence time, extremum,
and whether the observable returned past the exit threshold. An event that remains active at the end
of the selected trajectory is written with `returned_past_exit=0` and counted as open in the
summary. Each active sample contributes the selected-frame time step to the residence time,
including a one-sample event. Times start at zero for the first selected frame, and `--stride`
multiplies the input `--dt-fs` automatically. The summary reports `first_entry_ps=-1` when no event
occurred. Event counts and residence times are descriptive trajectory statistics; probabilities and
rates require an independently sampled trajectory ensemble.

## Radial distribution function

```console
./trajana.x rdf \
  --input PROJECT-pos-1.xyz \
  --cell-file PROJECT-1.cell \
  --select-a name:O \
  --select-b name:H \
  --bins 400 \
  --rmax 7.5 \
  --output rdf-oh.dat
```

The histogram uses all distinct ordered pairs from selections A and B. For each frame it is
normalized by the exact shell volume and instantaneous cell volume. This also gives the correct
normalization when the selections overlap. The third output column is the running mean number of B
atoms around an A atom. If `--rmax` is omitted, half the shortest periodic cell height of the first
selected frame is used. Larger cutoffs are rejected.

## Velocity autocorrelation and spectrum

The `vacf` mode expects an XYZ-formatted CP2K velocity trajectory:

```console
./trajana.x vacf \
  --input PROJECT-vel-1.xyz \
  --select name:O \
  --dt-fs 0.5 \
  --max-lag 10000 \
  --window hann \
  --output vacf-o.dat \
  --spectrum spectrum-o.dat \
  --peaks peaks-o.dat
```

The VACF is

```text
C(t) = average_i,t0 [v_i(t0) dot v_i(t0+t)]
```

and uses every available time origin with the unbiased `N-t` denominator. Autocorrelations and
spectra are evaluated with a dependency-free radix-2 FFT, not the proprietary ESSL routine used by
the old programs. `--remove-mean` removes the time-average of every Cartesian velocity component
before the calculation. `--window` affects only the optional spectrum and can be `none` or `hann`.
The spectrum is not mass weighted and does not include a quantum correction. Its frequency axis is
written simultaneously as wavenumber, THz, and meV. `--peaks` reports local maxima above
`--peak-threshold` (default: 0.05 of the largest nonzero-frequency peak), together with an
interpolated FWHM and the corresponding exponential-decay lifetime. Run the mode separately with
`--select name:H` and `--select name:O` for species-resolved spectra.

`--dt-fs` is the time between input frames. When `--stride` is used, the output time step is
multiplied by the stride automatically.

## Frequency-selective anharmonic modes

`fresean.x` implements the vibrational-analysis part of FREquency-SElective ANharmonic (FRESEAN)
mode analysis directly for CP2K trajectories. It forms the mass-weighted velocity cross-correlation
matrix at every sampled frequency, convolves it with a Gaussian spectral window, and diagonalizes
the resulting real symmetric matrix. Its eigenvectors identify collective motions contributing at a
chosen frequency without a harmonic or quasi-harmonic approximation.

```console
./fresean.x \
  --velocity PROJECT-vel-1.xyz \
  --position PROJECT-pos-1.xyz \
  --reference reference.xyz \
  --mass-file masses.dat \
  --dt-fs 0.5 \
  --correlation-frames 500 \
  --sigma-cm 10 \
  --frequency-cm 50 \
  --mode-count 10 \
  --output fresean-eigenvalues.dat \
  --mode-file fresean-modes.xyz \
  --mode-spectrum fresean-mode-spectra.dat \
  --mode-timeseries fresean-mode-timeseries.dat
```

The position and velocity files must contain synchronized frames with identical atom ordering.
Supplying `--position` enables the mass-weighted rotational fit used by FRESEAN: each selected
structure is aligned to `--reference`, and the same rotation is applied to its velocities. The first
selected position frame is the default reference. Molecules must be whole before fitting; use
`trajconvert.x unwrap` first when periodic crossings occur. Without positions, velocities are
analyzed in their input orientation, which is appropriate only when overall rotation has already
been removed or is intentionally part of the analysis.

The mass file uses the same `LABEL MASS[g/mol]` format as `twopt.x`. CP2K's default `bohr/au_time`
velocity unit is assumed; `--velocity-unit angstrom/ps` and `angstrom/fs` are also available. The
reported kinetic temperature is a sampling diagnostic inferred from the selected mass-weighted
velocities. Use `--constraints N` when constraints reduce the selected degrees of freedom; a
fractional value is accepted for a subsystem whose center-of-mass constraint is shared with the
surrounding system. `--remove-mean` is optional because solute translation can be a physically
meaningful zero-frequency contribution.

The main output contains the total VDoS, the largest eigenvalues, and their retained fraction at
every frequency. It is normalized so that summing the total positive-frequency VDoS gives the
effective number of selected degrees of freedom. `--mode-file` writes the leading eigenvectors at
the closest sampled frequency to `--frequency-cm`; components are divided by the square root of the
atomic mass so that they are suitable as displacement vectors. `--mode-spectrum` evaluates
`q^T C(omega) q` for those fixed modes over the full spectrum.

With synchronized positions, `--mode-timeseries` projects every aligned frame onto the fixed
mass-weighted eigenvectors selected at `--frequency-cm`. For mode vector `e_k`, it writes

```text
q_k(t)    = e_k^T M^(1/2) [R_aligned(t) - R_reference]
qdot_k(t) = e_k^T M^(1/2) V_aligned(t)
E_k(t)    = 1/2 [qdot_k(t)^2 + omega^2 q_k(t)^2].
```

The time mean of each projected coordinate is removed before evaluating the energy. Coordinates,
velocities, harmonic energy proxies in kJ/mol, energy-equivalent quanta, and instantaneous fractions
normalized over the extracted modes are written for every selected frame. This output requires a
positive `--frequency-cm` and `--position`; the latter ensures that coordinates are projected
directly instead of being reconstructed by integrating velocities. The energies use the selected
frequency for every extracted vector and are diagnostics rather than exact quantum populations.
Individual vectors inside a degenerate subspace are not unique, although sums over a complete
degenerate subspace are invariant.

The output header also reports the first `1/e` time of the unbiased autocorrelation of the
normalized complex phase vector `omega*q_k + i*qdot_k`. A value of `-1` means that no crossing
occurred within half of the selected trajectory. Under continuous external or thermostatted driving
this is an effective driven phase-correlation time, not a field-free
intramolecular-vibrational-redistribution lifetime.

The exact discrete procedure follows the equations and normalization demonstrated in the
[FRESEAN tutorial](https://github.com/HeydenLabASU-collab/FRESEAN-tutorial), including the symmetric
correlation interval of length `2*correlation-frames-1` and Gaussian convolution. The arbitrary-size
FFT and symmetric eigensolver are dependency-free Fortran implementations, so FFTW, GSL, GROMACS,
and Python are not required. At least twice as many selected trajectory frames as correlation points
are required. Matrix storage scales as `correlation-frames*(3*N_selected)^2`; use `--select` and
verify convergence with the trajectory length, correlation time, Gaussian width, and selection.

The method was introduced by Sauer and Heyden,
[*J. Chem. Theory Comput.* **19**, 5481-5490 (2023)](https://doi.org/10.1021/acs.jctc.2c01309). The
independent implementation was cross-checked against the vibrational matrix/eigenvector path in the
[reference repository](https://github.com/HeydenLabASU-collab/FRESEAN-metadynamics). Its
protein-specific coarse graining, mode-to-COLVAR conversion, PLUMED input, DCCM, resampling,
reweighting, and metadynamics workflows are intentionally outside this trajectory-analysis tool. The
optional generalized-normal-mode joint diagonalization is a separate approximation and is not used.
The tutorial's clustering example is deliberately not automated because its authors describe the
cutoff-based procedure as an exploratory example rather than a general analysis recipe.

## Two-phase thermodynamics

`twopt.x` is an independent Fortran implementation of the conventional Lin-Blanco-Goddard 2PT
method. It reads the XYZ velocity trajectory written by CP2K, constructs the mass-weighted velocity
density of states (DoS), separates its diffusive hard-sphere and harmonic solid parts, and reports
entropy, Helmholtz free energy, internal energy, heat capacity, zero-point energy, diffusivity,
fluidicity, and packing fraction. It does not invoke or copy the legacy C++ program.

For an atomic liquid:

```console
./twopt.x \
  --velocity PROJECT-vel-1.xyz \
  --mass-file masses.dat \
  --dt-fs 2.0 \
  --temperature 300 \
  --cell-file PROJECT-1.cell \
  --output 2pt-thermo.dat \
  --spectrum 2pt-spectrum.dat \
  --vacf 2pt-vacf.dat
```

The mass file maps trajectory labels to masses in g/mol:

```text
Ar 39.948
```

CP2K prints velocities in `bohr/au_time` by default, and this is therefore the default assumed by
`twopt.x`. If `&MOTION / &PRINT / &VELOCITIES` contains a different `UNIT`, pass
`--velocity-unit angstrom/ps` or `--velocity-unit angstrom/fs`. `--velocity-scale` is an additional
multiplier for uncommon input conventions. The tool converts velocities internally to angstrom/ps
and reports the kinetic temperatures of all resolved channels as a consistency diagnostic.

Molecular 2PT additionally requires synchronized positions and one molecule per line in a group
file. It is not restricted to water:

```text
M1 1 2 3
M2 4 5 6
```

```console
./twopt.x \
  --velocity PROJECT-vel-1.xyz \
  --position PROJECT-pos-1.xyz \
  --groups molecules.groups \
  --mass-file masses.dat \
  --cell-file PROJECT-1.cell \
  --dt-fs 2.0 \
  --temperature 300 \
  --rotational-symmetry 2 \
  --constraints 2 \
  --output 2pt-water.dat \
  --spectrum 2pt-water-spectrum.dat
```

For every frame, the molecular velocity is decomposed into center-of-mass translation, rigid-body
rotation, and internal vibration. Molecules split across periodic boundaries are reconstructed with
the triclinic minimum image. The inertia tensor is diagonalized without an external linear-algebra
dependency; linear molecules are detected from its null principal moment. Translation and rotation
receive independent self-consistent fluidicities, while internal motion is harmonic. Every atom must
occur in exactly one group. `--rotational-symmetry` currently applies to all molecular groups, so a
mixed system with different symmetry numbers should be analyzed in separate homogeneous runs or with
a more specialized mixture implementation.

The DoS output obeys the sum rules for the translational, rotational, and unconstrained internal
degrees of freedom. The default `--window none` matches the standard 2PT definition and preserves
the zero-frequency DoS used for diffusion; `--window hann` is available for sensitivity checks but
changes that quantity. System center-of-mass drift is removed frame by frame unless
`--keep-system-drift` is specified. The `lin2003` hard-sphere entropy convention is the default for
compatibility with established 2PT results; `--entropy-convention rigorous` omits the additional
compressibility term used by that convention.

The reported energy and free energy are spectral contributions unless the mean classical MD energy
of the simulation box is supplied with `--energy-kj-mol`. This adds the usual reference
`E_MD-E_classical,DoS`. Likewise, `--classical-cv-j-mol-k` accepts a heat capacity obtained from
energy fluctuations and applies only the quantum correction from the DoS. Without these optional
values, the entropy, DoS partition, ZPE, and diffusion remain usable, but the total absolute energy,
free energy, and liquid heat capacity must not be over-interpreted.

The zero-frequency DoS is particularly sensitive to sampling. Production use should normally start
from at least several thousand equilibrated frames and demonstrate convergence with trajectory
length, print interval, window choice, and system size. The implementation follows the method and
output definitions described by Lin, Blanco, and Goddard (*J. Chem. Phys.* **119**, 11792, 2003),
Lin, Maiti, and Goddard (*J. Phys. Chem. B* **114**, 8191, 2010), and Pascal, Lin, and Goddard
(*PCCP* **13**, 169, 2011). The legacy and current `py-xPT` implementations in
[`atlas-nano/codes`](https://github.com/atlas-nano/codes) were used as cross-checking references.
The 3PT memory-cage refinement, finite-size diffusivity corrections, and LAMMPS-specific DMA
workflow in that repository are intentionally outside the scope of this compact CP2K trajectory
tool.

## Spatially resolved two-phase thermodynamics

`twopt3d.x` implements the trajectory-based part of the 3D-2PT method for CP2K. It assigns each
solvent molecule to the voxel containing its center of mass at the correlation time origin,
accumulates the forward translational and rotational velocity autocorrelation functions there, and
converts their symmetrized spectra into local 2PT entropies. The resulting Gaussian CUBE maps
contain the molecular number density, translational entropy, rotational entropy, total
intermolecular entropy, and translational and rotational fluidicities.

The solvent is defined with the same group and mass files as molecular `twopt.x`. All groups must
describe the same molecular species in the same atom order, but the species is not restricted to
water:

```console
./twopt3d.x \
  --velocity PROJECT-vel-1.xyz \
  --position PROJECT-pos-1.xyz \
  --groups solvent.groups \
  --mass-file masses.dat \
  --cell-file PROJECT-1.cell \
  --dt-fs 0.5 \
  --temperature 300 \
  --grid "80 80 80" \
  --spacing 0.5 \
  --origin "-20 -20 -20" \
  --correlation-frames 500 \
  --minimum-origins 50 \
  --rotational-symmetry 2 \
  --output-prefix solvation
```

For a mobile solute, `--align-select index:1-42` translates and rotates every solvent origin into a
solute-fixed frame. `--reference reference.xyz` is optional; the first selected position frame is
otherwise used. Alignment atoms and solvent molecules must be whole. A cell is used to reconstruct
split solvent molecules and to take the minimum-image displacement from the aligned solute. At least
three non-collinear alignment atoms are required so that all three rotational axes are defined.
Without alignment, molecular centers are wrapped into the supplied periodic cell and the grid
remains in the laboratory frame, which is useful for interfaces and crystals.

The local translational signal is `sqrt(M) V_COM`. The rotational signal follows the current
[`trvdos`](https://github.com/HeydenLabASU-collab/trvdos) implementation: angular momentum is
projected onto the instantaneous principal axes and divided by the square root of the corresponding
principal moment. Axis signs are followed continuously. This generalizes the fixed water axes and
hard-coded water moments of the original [`3D-2PT`](https://github.com/HeydenLabASU/3D-2PT) program.
Linear molecules receive two rotational degrees of freedom and use the laboratory-frame
angular-momentum vector to avoid an arbitrary basis inside their degenerate perpendicular
eigenspace; nonlinear molecules receive three. The current thermodynamic map is intended for rigid
solvent models. Internal vibrations and torsions of flexible complex solvents are not included in
the reported entropy and should be analyzed separately until a topology-aware internal-coordinate
extension is available.

Each voxel spectrum is normalized to its molecular degrees of freedom before applying the same 2PT
partition and quantum harmonic weights as `twopt.x`. Finite local sampling can make a symmetrized
VACF spectrum slightly negative. Such bins are projected to zero before normalization; the removed
spectral fraction is written to `*-negative-vdos-translation.cube` and
`*-negative-vdos-rotation.cube`. `*-origins.cube` contains the independent time-origin count used
for the density and the `--minimum-origins` mask. Optional `--vacf` and `--spectrum` files expose
the voxel-resolved intermediate data for convergence checks. Production maps require convergence
with trajectory length, correlation time, voxel size, origin count, and the optional lag `--window`.

With `--bulk-entropy VALUE`, where `VALUE` is the translational plus rotational bulk-solvent entropy
in J/(mol K), the program additionally writes the local excess entropy, `-T Delta S` per solvent
molecule, and its number-density-weighted contribution. The summary contains the spatial integrals
of excess entropy and `-T Delta S`. The reference can be obtained from a separate bulk calculation
with `twopt.x` (sum of its translation and rotation channels) or from a converged bulk-like region.
Voxels below `--minimum-origins` are excluded from these integrals.

The maps use the standard Gaussian CUBE format and can be combined with existing tools such as
CP2K's `tools/cubecruncher` or [`cubeWorks`](https://github.com/HeydenLabASU/cubeWorks). Generic
CUBE arithmetic, filtering, and peak finding are not duplicated in `twopt3d.x`; thermodynamic
density weighting is performed directly so that its units and sampling mask remain explicit.

The original package also computes solute-solvent and solvent-solvent Lennard-Jones/Coulomb energy
maps from a classical pairwise force field and combines them with entropy into enthalpy and free
energy. `twopt3d.x` deliberately does not label such maps as generally available from CP2K: an
arbitrary DFT or many-body CP2K potential has no unique per-solvent pair-energy decomposition. The
entropy and density maps are therefore methodologically portable, while enthalpy and solvation free
energy require a separate, explicitly justified energy decomposition for the chosen Hamiltonian.

The spatial conditioning and thermodynamic definitions follow Persson et al.,
[*J. Chem. Theory Comput.* **13**, 4467-4481 (2017)](https://doi.org/10.1021/acs.jctc.7b00184).

## Hydrogen-bond network

The `hbond` mode is a tested reimplementation of the scientifically useful part of
`science_count_nn_zone.f`; it does not reuse that fixed-size F77 code. Molecules or donor groups are
defined by a group file. The first index is the donor/acceptor heavy atom and all following indices
are its donor hydrogens:

```text
W1 1 2 3
W2 4 5 6
```

Run the analysis with:

```console
./trajana.x hbond \
  --input PROJECT-pos-1.xyz \
  --cell-file PROJECT-1.cell \
  --groups examples/water.groups \
  --output hbond-frames.dat \
  --summary hbond-populations.dat
```

For every donor O-H and acceptor O pair the default Wernet-type criterion is

```text
r_OO < 3.3 angstrom - 0.00044 angstrom/degree^2 * theta^2
theta <= 45 degrees
r_OH <= 1.25 angstrom
```

where `theta` is the O-H...O deviation from linear donation as measured at the donor. The constants
can be changed with `--roo-zero`, `--angle-coefficient`, `--angle-max`, and `--oh-max`. The frame
output reports the number of bonds and the free donor-H fraction. The optional summary contains the
population of each `(accepted, donated)` group state. Bifurcated bonds are counted individually.

## Scattering functions and collective currents

`dynStruct.f90` was not ported because its complex exponential was reduced to `exp(real(i*x)) = 1`
and its atom-pair loop did not use the second atom. The `dsf` mode implements the coherent density
correlation from its definition:

```text
rho(q,t) = sum_j exp(i q dot r_j(t))
F(q,tau) = average_q,t [rho(q,t+tau) rho(q,t)*] / N
```

Cartesian q vectors in inverse angstrom are provided one per line. All vectors in a file must have
the same magnitude; list symmetry-equivalent vectors to average one reciprocal-space shell. Use one
run per magnitude when constructing a dispersion curve.

The basic coherent calculation is:

```console
./trajana.x dsf \
  --input PROJECT-pos-1.xyz \
  --select name:O \
  --q-file examples/q-vectors.dat \
  --dt-fs 0.5 \
  --output intermediate-scattering.dat \
  --spectrum dynamic-structure-factor.dat
```

The first output contains the real and imaginary intermediate scattering function and its normalized
real part. The optional spectrum is the averaged periodogram of the density modes and also writes
`omega^2*S(q,omega)/q^2`. The latter can be compared directly with a velocity spectrum in the
small-q limit. `--remove-mean` removes elastic mean-mode components, and `--window none|hann`
controls spectral leakage. Frequency axes are written as wavenumber, THz, and meV.

The self (incoherent) intermediate scattering function

```text
F_self(q,tau) = average_q,j,t [exp(i q dot (r_j(t+tau)-r_j(t)))]
```

is requested from the same position trajectory with:

```console
./trajana.x dsf \
  --input PROJECT-pos-1.xyz \
  --select name:H \
  --q-file examples/q-vectors.dat \
  --dt-fs 0.5 \
  --output coherent-h.dat \
  --self-output incoherent-h.dat \
  --self-spectrum incoherent-spectrum-h.dat
```

For the longitudinal and two-polarization-averaged transverse molecular-current correlations, add a
synchronized XYZ velocity trajectory:

```console
./trajana.x dsf \
  --input water-com-pos.xyz \
  --velocity water-com-vel.xyz \
  --select all \
  --q-file q-shell-1.dat \
  --dt-fs 0.5 \
  --window hann \
  --output coherent.dat \
  --spectrum coherent-spectrum.dat \
  --current currents.dat \
  --current-spectrum current-spectra.dat \
  --summary moments.dat \
  --peaks collective-peaks.dat
```

This evaluates

```text
J(q,t) = sum_j v_j(t) exp(i q dot r_j(t))
C_L(q,tau) = average [J_L(q,t+tau) J_L(q,t)*] / N
C_T(q,tau) = average over the two transverse polarizations / N
```

The current-spectrum output also contains `omega^2*S(q,omega)/q^2`, which should agree with the
longitudinal spectrum when positions and velocities use consistent units. If velocities are not in
angstrom/fs, convert them with `--velocity-scale FACTOR`. The summary contains `S(q)`, both current
values at zero lag, the density frequency `omega0^2=q^2*C_L(0)/S(q)`, the second current-spectrum
moments, and their high-frequency velocities. With `--mass-density RHO` in g/cm3 it additionally
writes `rho*c^2` in GPa. These moment-derived quantities should only be interpreted after
convergence with respect to trajectory length, time step, window, and q.

`--peaks` can be used by VACF and `dsf`. For collective spectra it reports q, peak frequency,
energy, phase velocity, FWHM, and the lifetime convention `tau=1/(pi*FWHM_frequency)`. Multiple
longitudinal peaks are retained rather than collapsed into a single sound branch.

For molecular-center currents, first generate matched center-of-mass position and velocity files.
The same group and mass files must be used for both:

```console
./trajconvert.x center --input PROJECT-pos-1.xyz --groups water.groups \
  --mode mass --mass-file masses.dat --cell-file PROJECT-1.cell --output water-com-pos.xyz
./trajconvert.x center --input PROJECT-vel-1.xyz --groups water.groups \
  --mode mass --mass-file masses.dat --vectors --output water-com-vel.xyz
```

### Coverage of the water analyses

The implemented trajectory observables cover the reusable analyses in the two classic water studies:

- O-O, O-H, and H-H pair distributions and running coordination numbers: `rdf` with the
  corresponding selections.
- H and O velocity autocorrelations and vibrational spectra: `vacf`, including meV output and peak
  positions.
- Proton incoherent scattering: `dsf --select name:H --self-output ... --self-spectrum ...`.
- Center-of-mass static/dynamic structure factors and longitudinal/transverse current correlations:
  `dsf` on center trajectories, including spectral peaks, widths, lifetimes, dispersion velocities,
  and second moments.

The force-field parameter tables, free-molecule normal-mode calculations, and the historical
memory-function fit are not trajectory observables and are intentionally not reproduced. Modern
trajectories can be transformed directly; extrapolating a short correlation with a chosen memory
kernel would add a model assumption rather than another measured observable.

## Trajectory transformations

Wrap atoms or convert them to fractional coordinates with:

```console
./trajconvert.x wrap --input in.xyz --cell-file PROJECT-1.cell --output wrapped.xyz
./trajconvert.x fractional --input in.xyz --cell-file PROJECT-1.cell --output fractional.xyz
```

Atomwise unwrapping follows crossings between consecutive frames:

```console
./trajconvert.x unwrap --input wrapped.xyz --cell-file PROJECT-1.cell --output continuous.xyz
```

As with any frame-to-frame unwrap, each fractional displacement must be less than half a periodic
cell length between stored frames.

Group centers use the same group format as the hydrogen-bond analysis:

```console
./trajconvert.x center \
  --input PROJECT-pos-1.xyz \
  --cell-file PROJECT-1.cell \
  --groups examples/water.groups \
  --mode mass \
  --mass-file examples/masses.dat \
  --wrap-output \
  --output water-com.xyz
```

`--mode geometry` assigns every atom unit weight. `--mode mass` requires an explicit mass for every
XYZ label, avoiding guesses for custom CP2K kind names. Groups crossing a periodic boundary are
reconstructed around their first atom before their center is calculated. For velocity or other
vector trajectories, use `--vectors` to disable this positional reconstruction.

## License

These sources are part of CP2K and are distributed under the same `GPL-2.0-or-later` license. The
external programs and publications cited above were used as methodological references and
cross-checks; no third-party source code is included in this tool.
