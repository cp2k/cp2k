# Gaussian fits of OpenMX PAO basis sets

`data/BASIS_OMX` contains approximate Gaussian representations of numerical pseudo-atomic orbitals
(PAOs) from the [OpenMX 2019 database](https://www.openmx-square.org/vps_pao2019/). The basis names
are `OMX-FIT-<OpenMX selection>`, for example `OMX-FIT-H7.0-s2p1`. No changes to CP2K's basis reader
or electronic solvers are required.

## Scope and limitations

The shipped library contains **279 selections covering all 81 elements** in the pinned catalog. The
catalog records 268 PAO files and 1,415 shell selections; the other 1,136 selections fail one or
more numerical fit thresholds and are not shipped. This is not every possible combination of radial
functions. Published shell-count patterns are expanded over the linked cutoffs within each
element/valence family; this expansion is a construction policy, not an OpenMX recommendation for
every resulting combination.

This is **not a lossless format conversion**. Compactly supported numerical PAOs are approximated by
Gaussian contractions with nonzero tails. Angular momentum and radial multiplicity, hence the AO
dimension of each selection, are preserved; radial shapes, overlaps and matrix elements are
approximate. The converter uses default CP2K primitive and contraction normalization. Do not change
normalization options when using these fits.

Passing the fit thresholds does **not** establish production accuracy for energies, forces, stress,
transferability, or multi-center overlap conditioning. The fits use up to 64 Gaussian primitives per
angular channel and are not optimized for CP2K performance. Benchmark against established CP2K bases
before using them in an application.

No pseudopotentials, pseudopotential fits, density-matrix converters or restart writers are
included. Select the pseudopotential and its valence electron count explicitly. Using a GTH
potential with these bases does not reproduce the OpenMX Hamiltonian. In particular, source
semicore/open-core families must not be assumed compatible with a GTH potential just because the
element matches. A raw OpenMX density matrix cannot be used unchanged: AO ordering, normalization,
changed radial functions and the target overlap metric matter. The permutation recorded in the JSON
is only angular-order metadata, not a complete density-matrix transformation.

## Using the library

For the H basis, for example:

```text
&DFT
  BASIS_SET_FILE_NAME BASIS_OMX
  POTENTIAL_FILE_NAME GTH_POTENTIALS
  ...
&END DFT
&SUBSYS
  ...
  &KIND H
    BASIS_SET OMX-FIT-H7.0-s2p1
    POTENTIAL GTH-PBE-q1
  &END KIND
&END SUBSYS
```

`examples/H2O_basis_smoke.inp` is a runnable PBE/GTH test using only shipped entries: O7.0-s1p1 and
H7.0-s2p1. It has 14 AOs and eight valence electrons. Run it from `tools/openmx/examples` so the
relative data paths resolve. It demonstrates successful basis use, not chemical accuracy or
agreement with an OpenMX calculation.

## Numerical method and acceptance

The fitter reads the radial PAO tables and uses cubic Hermite interpolation with three-point
interior derivatives, a regular `r**l` inner extrapolation, and zero beyond the outer PAO grid. The
inner extrapolation is an approximation to OpenMX's boundary handling.

Each angular channel is fitted with shared even-tempered Gaussian exponents and independent radial
contractions. An SVD least-squares solve determines the coefficients; an outer optimization varies
the exponent endpoints. The objective combines radial L2 and gradient errors with derivative weight
0.2 bohr. Gaussian tails are included. A higher-order quadrature independently checks the normalized
exported functions.

The database builder tries 24, 40 and 64 primitives as needed. Every exported selection must satisfy
all of these bounds, channel by channel:

- Normalized radial L2 error at most 0.01.
- Absolute error in each radial orbital's kinetic expectation at most 0.001 hartree.
- Maximum same-center overlap error at most 0.001.
- Analytic versus numerical Gaussian normalization error at most 1e-6.

These are orbital-fit criteria, not total-energy error bounds. Diagnostics also record tail leakage,
source radial norms, SVD rank and exponent-search convergence. The basis file is accompanied by
`BASIS_OMX.json`, including source hashes/URLs, fit settings, full diagnostics for exported
selections, and rejection reasons for omitted selections. Cached fits retain detailed diagnostics
for omitted selections as well.

## Rebuilding and fitting additional selections

Requires Python 3.10 or newer, NumPy and SciPy:

```sh
python3 -m pip install -r tools/openmx/requirements.txt
python3 tools/openmx/collect_database.py \
  --cache /path/to/openmx-data \
  --manifest tools/openmx/openmx2019_sources.json --fetch-only
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python3 tools/openmx/build_database.py \
  --manifest tools/openmx/openmx2019_sources.json \
  --pao-dir /path/to/openmx-data/pao \
  --fit-cache /path/to/openmx-data/fits \
  --output /path/to/new/BASIS_OMX --workers 4
```

Use a new output path. Existing outputs require explicit `--overwrite`. The downloader verifies PAO
SHA-256 hashes against the pinned manifest. Generation is offline once PAOs are present. Fits are
cached with source/code hashes and numerical settings; the build is deterministic within a numerical
environment, but NumPy/SciPy/BLAS differences can change final digits or threshold decisions.
Recheck diagnostics and validation after regeneration.

To discover a new catalog, omit `--fetch-only` and choose a new manifest path; review it rather than
silently replacing the pinned catalog. `build_database.py --include-inaccurate` explicitly exports
numerically normalizable fits that fail accuracy limits, with per-entry warnings. It was **not**
used for the shipped library.

The lower-level converter can fit an explicit selection, including one not in the catalog:

```sh
python3 tools/openmx/convert_basis.py H7.0-s2p1 \
  --pao-dir /path/to/openmx-data/pao --output /path/to/new/BASIS_H
```

Unlike the database builder, this single-fit interface enforces only its radial L2 limit and the
normalization check by default; inspect the other diagnostics before accepting its output.
`--allow-inaccurate` relaxes only the radial L2 rejection and does not certify accuracy.

## Tests

Offline tests cover the parser, synthetic radial fits, real-harmonic ordering, source pinning, cache
reuse, quality gates, and every shipped basis entry:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python3 -m unittest discover -s tools/openmx -p 'test_*.py' -v
```

Set `CP2K_EXE` to the absolute executable path to additionally run the H2O SCF test. It checks
convergence, 14 AOs, eight electrons, `Tr(P S)`, positive overlap, unit AO norms, and kinetic
diagonal elements against fit diagnostics. No OpenMX installation or downloaded PAOs are needed for
these tests.

To load **every** shipped entry through the actual CP2K basis reader:

```sh
python3 tools/openmx/validate_database.py data/BASIS_OMX \
  --cp2k /absolute/path/to/cp2k.psmp \
  --output /path/to/new/basis-validation.json
```

This uses ghost kinds and `RUN_TYPE NONE`; it tests parsing and initialization without choosing
physical pseudopotentials. Without `--cp2k`, only independent Gaussian normalization and
basis/report consistency are checked.

## Provenance and attribution

The OpenMX 2019 database credits T. Ozaki and H. Kawai; its authors retain copyright and distribute
the data under the GNU General Public Licence without warranty. These derived Gaussian fits retain
the source attribution and SHA-256 provenance. The original PAO files are not vendored. The
conversion tools are GPL-2.0-or-later, like CP2K.

The source database cites T. Ozaki, Phys. Rev. B **67**, 155108 (2003), and T. Ozaki and H. Kino,
Phys. Rev. B **69**, 195113 (2004).
