# Modular trajectory analysis tools

This directory contains two standalone Fortran programs for post-processing CP2K trajectories:

- `trajana.x` performs geometry, RDF, VACF, hydrogen-bond, and dynamic structure-factor analyses.
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
  --spectrum spectrum-o.dat
```

The VACF is

```text
C(t) = average_i,t0 [v_i(t0) dot v_i(t0+t)]
```

and uses every available time origin with the unbiased `N-t` denominator. Autocorrelations and
spectra are evaluated with a dependency-free radix-2 FFT, not the proprietary ESSL routine used by
the old programs. `--remove-mean` removes the time-average of every Cartesian velocity component
before the calculation. `--window` affects only the optional spectrum and can be `none` or `hann`.
The spectrum is not mass weighted and does not include a quantum correction.

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

## Dynamic structure factor

`dynStruct.f90` was not ported because its complex exponential was reduced to `exp(real(i*x)) = 1`
and its atom-pair loop did not use the second atom. The `dsf` mode implements the coherent density
correlation from its definition:

```text
rho(q,t) = sum_j exp(i q dot r_j(t))
F(q,tau) = average_q,t [rho(q,t+tau) rho(q,t)*] / N
```

Cartesian q vectors in inverse angstrom are provided one per line:

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
real part. The optional spectrum is the averaged periodogram of the density modes. Include
symmetry-equivalent `q` vectors in the file when a shell average is desired. `--remove-mean` removes
the elastic mean-density component, and `--window none|hann` controls spectral leakage.

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
