# Scalar NC-UPF validation

These checks cover the native scalar norm-conserving UPF 2.0.1 operators. They do not establish
USPP, PAW, spin-orbit UPF, or complete SSSP transferability support. SSSP release v2.0 is distinct
from the UPF format version.

## Reproducible regression and integral checks

Build `cp2k-bin`, `atom_upf_unittest`, `upf_projector_integrals_unittest`,
`upf_local_integrals_unittest`, `gapw_background_unittest`, and `gapw_1c_basis_unittest` with the
usual CP2K CMake configuration. Run each unit executable on one and two MPI ranks. The UPF integral
executables optionally accept UPF filenames for additional finite-integral checks. The SSSP
download, hash verification, and import-audit commands are given in [README.md](README.md).

Run the registered energy/force tests with the default timeout:

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_STACKSIZE=128M
ulimit -s unlimited
python3 tests/do_regtest.py --restrictdir QS/regtest-kind \
  --restrictdir QS/regtest-upf --mpiranks 2 --ompthreads 1 \
  --maxtasks 4 /path/to/build/bin psmp
```

Validation used GNU Fortran 15.2, MPI, OpenBLAS, FFTW, DBCSR 2.9.1, and a CPU build with bounds
checking. LIBXC and LIBINT2 were disabled; PBE and PADE used CP2K's internal functionals.
Calculations ran on Terok. The new SBr2 regression cases took 23, 57, and 38 seconds for GPW,
GAPW_XC, and GAPW, respectively, and passed all nine matchers. The 12 existing KIND cases also
passed.

Independent analytic Gaussian references check both local and nonlocal operators, their center
derivatives, and required moments. The local matrix, contracted derivative, and RI maximum scaled
discrepancies were 1.65e-13, 2.28e-14, and 2.40e-15. The coarse-EPS_PPL negative control rejects
truncation of the radial integration mesh. The one-center basis negative control rejects discarding
admissible diffuse exponents after an earlier candidate was rejected.

All 82 NC entries in the four pinned SSSP collections (58 distinct files) passed the import and
finite-integral audit. Among them, 40 contain NLCC data; their maximum relative radial L2
core-density fit error was 4.54e-7. These checks do not replace self-consistent benchmarks across
all elements and materials.

## Finite-difference controls

For the larger SBr2 controls, start from the registered inputs, restore a 10-Angstrom cubic cell,
translate all atoms by (1, 1, 1) Angstrom, use 600/80-Ry cutoffs, and tighten EPS_SCF to 1e-10. Keep
the unmodified UPFs, UZH-MOLOPT DZVP bases, and the GAPW hard radius and atomic quadratures. Compare
the sulfur z force with total energies at z displacements of plus/minus 0.0005 Angstrom. For yy
stress, set `STRESS_TENSOR ANALYTICAL`, scale both the y cell vector and all y coordinates by 1
plus/minus 1e-4, and use minus the energy strain derivative divided by the original volume. Compare
final FORCE_EVAL energies, not the intermediate SCF energy field.

| Mode    |  Center energy (Ha) | Force discrepancy (Ha/bohr) | yy stress discrepancy (bar) |
| ------- | ------------------: | --------------------------: | --------------------------: |
| GPW     | -42.022299953158296 |                     7.77e-9 |                    0.001118 |
| GAPW_XC |  -42.04622862017529 |                     7.54e-9 |                    0.001444 |
| GAPW    |  -42.04311512017074 |                     6.61e-9 |                    0.001134 |

These 15 converged calculations check derivatives of each discrete energy. Their differences across
modes do not establish a converged density splitting. Separate charged, spin-polarized controls and
a reciprocal-space Gaussian reference verify the fixed-charge periodic GAPW background repair. The
charged molecular controls predate the native local-potential replacement and used the former fitted
local potential. Broader response, CNEO, and planar-countercharge validation remains outside these
checks.

## Independent local-operator and basis controls

The existing `tests/QS/regtest-kind/H.pbe-hgh.UPF` table has an analytic HGH reference with Z=1,
local radius 0.2 bohr, and local coefficients (-4.178900438, 0.724463313). The two radial functions
agree to 1.5e-14 Ha pointwise. H2 with positions (0,0,0) and (0.72,0,0) Angstrom in an 8-Angstrom
cell, aug-TZV2P-GTH, PADE, and 1000/100-Ry grids gives a GPW energy difference of 7.62e-8 Ha between
analytic and tabulated potentials.

For the GAPW compensation control use FORCE_PAW, GAPW_1C_BASIS ORB, EPS_DEFAULT=1e-12,
EPSRHO0=1e-12, EPSSVD=1e-12, LMAXN0=LMAXN1=6, 400 radial and 590 angular points, and EPS_SCF=1e-10
with the default diagonalization solver. Set ALPHA0_HARD to the respective Gaussian core exponent:
12.5 for analytic HGH and 9.876543209876543 for UPF. With EPSFIT=1e-10 the all-soft limit recovers
each GPW energy within 2.1e-13 Ha. With EPSFIT=1e-4 and HARD_EXP_RADIUS=1.2, a nonzero one-center XC
correction of about -0.016295 Ha is retained. Analytic HGH and UPF energies differ by 2.93e-9 Ha,
while both remain about 9.1 microhartree above their GPW limits. This distinguishes operator
agreement from decomposition convergence. Matching the compensation exponent is a diagnostic
cancellation, not a general parameter recommendation. Extended one-center bases can be
ill-conditioned; unsuccessful OT controls are not counted as accuracy results.

An enlarged diagnostic SBr2 basis with 14 geometric exponents from 0.04 to 24 bohr^-2 in each s
through g channel, canonical overlap truncation at 1e-6, 1000/80-Ry grids, EPS_DEFAULT=1e-14,
CHOLESKY OFF, and EPS_SCF=1e-10 gives -42.255480473001 Ha. Original-UPF Quantum ESPRESSO at 120/960
Ry gives -42.255645780 Ha. The difference is 0.165307 mHa; the sulfur z force differs by 1.0874e-5
Ha/bohr. These are finite-basis comparisons, not universal error bounds or optimized new bases.
Existing UZH-MOLOPT bases remain usable starting points with matching valence counts and explicit
convergence checks.
