# Modular trajectory analysis tools

This directory contains two standalone Fortran programs for post-processing CP2K trajectories:

- `trajana.x` performs geometry, RDF, VACF, hydrogen-bond, scattering, and collective-current
  analyses.
- `trajconvert.x` wraps, unwraps, converts, and reduces trajectories to group centers.

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
