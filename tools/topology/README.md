# CP2K overlaps, Wilson loops, and Z2Pack

Development implementation; not yet a production-validated topology facility. The changes are in
CP2K's `DFT/PRINT/WANNIER90` path. No Wannier90 library or Wannier fit is needed for explicit
overlap loops. Z2Pack remains an independent Python dependency, not a CP2K build dependency.

## Two calculation routes

1. `KPOINTS_SOURCE NNKP` reads arbitrary fractional k-points and directed connections (including the
   reciprocal-vector closure) from a Wannier90 `.nnkp` file. `CP2KSystem` implements Z2Pack's
   `OverlapSystem` and supplies those overlaps to Z2Pack's line/surface calculations.
1. `WILSON_LOOP T` evaluates those loops inside CP2K. `KPOINTS_SOURCE WILSON` generates a surface
   internally and doubles both resolutions until the WCC converge. `Z2 T` additionally evaluates
   largest-gap crossing parity.

Both routes use the actual Gaussian-basis Berry operator:

`M(k,b) = C(k)^dagger O(k,b) C(k+b)`.

The AO matrix contains the nonorthogonal-basis metric and phase information. The explicit-loop path
uses ordered, nonsymmetric atom-pair matrices and
`O(k,b) = sum_R exp(i*(k+b).R) <mu,0|exp(-i*b.r)|nu,R>`, with the same periodic atom images as
CP2K's k-point Hamiltonian. Plain Euclidean overlaps of AO coefficient vectors are not correct. The
legacy regular-grid exporter is unchanged; its Hermitian-operator shortcuts are not a reference for
cross-k links. The SCF potential is kept fixed while the requested post-SCF k-points are
diagonalized. The adapter runs the same SCF problem for every requested line; its SCF mesh must
never be replaced by the current Wilson loop.

## Z2Pack interface

Install `numpy`, `scipy`, `z2pack`, and optionally `pytest`. Put this directory on `PYTHONPATH`.
Start from `examples/helium.inp` and use:

```python
import numpy as np
import z2pack
from cp2k_z2pack import CP2KSystem

system = CP2KSystem(
    input_file="examples/helium.inp",
    lattice=4 * np.eye(3),  # direct lattice vectors as ROWS, in Angstrom
    command=["/path/to/patched/cp2k.psmp"],
    workdir="helium-loops",
    num_bands=1,
    env={"OMP_NUM_THREADS": "2"},
)
result = z2pack.line.run(system=system, line=lambda t: [t, 0, 0])
print(result.wcc)
```

The CP2K input must contain `KPOINTS_SOURCE NNKP`, `NNKP_FILE loop.nnkp`, and `SEED_NAME loop`. Set
`WILSON_LOOP T` to also calculate native phases and check the sampled band separation. Choose
`EXCLUDE_BANDS` explicitly to select the intended subspace. `SPIN_CHANNEL` selects one collinear
channel; UKS channels are never concatenated in a single `.mmn` file. In SOC mode, exclusions select
spinor bands instead.

`input_files` copies named dependencies, including restart files, to every run. `command` is an
argument list (for example with `mpiexec`), never shell text. Each line calculation has its own
directory with the input, requested points, CP2K log, `.eig`, `.mmn`, and singular-value
diagnostics. Failed runs are retained. The final point must equal the first modulo an integer
reciprocal vector.

By default the adapter polar-unitarizes each overlap using SVD, matching the native kernel.
`polar=False` exposes raw CP2K overlaps to Z2Pack. Comparing both at increasingly fine sampling is
an additional useful numerical check. Agreement of two analyses of the same overlaps does **not**
independently validate the underlying AO operator or the DFT Hamiltonian.

For a spinful time-reversal insulator, use `z2pack.surface.run` on a half surface such as
`surface=lambda s,t: [t, s/2, 0]`, and check every line and surface convergence criterion before
calling `z2pack.invariant.z2(result)`.

## Native calculation

`examples/neon-soc.inp` is a validated smoke test for a trivial closed-shell SOC insulator.
`examples/stanene-soc.inp` tests a bonded two-atom system with full-basis second-variational SOC.
The syntax is:

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
&END
```

The example has 8 electrons and 5 scalar bands: 10 spinors, of which 8 are retained. The band count
must be adapted for every material. Second-variational SOC uses CP2K's pseudopotential SOC integrals
and a restricted SCF. Converge the number of scalar unoccupied orbitals (`SCF/ADDED_MOS` plus the
additional `WANNIER90/ADDED_MOS`) and all ordinary DFT cutoffs. SOC-capable pseudopotentials are
required. This is not self-consistent noncollinear DFT.

`TIME_REVERSAL T` is a user assertion about the Hamiltonian, not an automatic symmetry proof. Native
Z2 requires an even, lowest-energy occupied spinor subspace with one state per electron and an
excluded conduction band. Presently native Z2 surfaces use two distinct reciprocal coordinate axes,
with winding one along the loop and one half transversely. Arbitrary loops remain available for
Wilson phases and for the external Z2Pack route.

The native checks cover singular links, sampled band gaps, boundary Kramers pairs, WCC changes under
joint refinement, adjacent-line movement and gap separation, and parity stability. A coarse mesh can
still miss a gap closing or rapid evolution between samples. Repeat with finer starting meshes and
tighter tolerances. A single time-reversal plane gives one 2D invariant, not all four 3D strong/weak
indices. Gapless graphene does not have a well-defined insulating Z2 invariant without specifying
and resolving a gap-opening Hamiltonian.

The `.wilson` file contains loop index, minimum link singular value, and sorted WCC in `[0,1)`.
Printed Berry phases are in radians. Refinements overwrite the seed output with the final mesh; the
log retains diagnostics from all levels.

## Validation

The mandatory Misc CI check runs only the short model and adapter tests, including known Z2=0 and
Z2=1 models. The regular CP2K regression suite contains only the small helium and neon smoke tests.
The more expensive DFT/SOC references, stanene surfaces, adaptive Z2Pack surfaces, and MPI sweeps
below are manual validation tools; they are deliberately not added to automatic CI.

For the model tests, `FC` selects a GNU Fortran compiler and `WILSON_LAPACK_FLAGS` overrides the
default `-llapack -lblas` link flags. The tests fail, rather than silently skip, if the compiler is
missing. Manual checks are available as follows; the stanene calculations can take substantially
longer than the smoke tests:

```sh
python -m pytest -q tools/topology/test_topology.py
OMP_NUM_THREADS=1 python tools/topology/validate_cp2k.py /path/cp2k.psmp /path/results --soc
python tools/topology/validate_overlap_geometry.py /path/cp2k.psmp /path/geometry
python tools/topology/validate_multicentre_overlap.py /path/cp2k.psmp /path/multicentre
python tools/topology/validate_soc_reference.py /path/cp2k.psmp /path/neon-reference
python tools/topology/validate_soc_reference.py /path/cp2k.psmp /path/stanene-reference --stanene
python tools/topology/validate_stanene.py /path/cp2k.psmp /path/stanene
python tools/topology/validate_z2pack_stanene.py /path/cp2k.psmp /path/stanene/stanene-RESTART.kp /path/stanene-z2pack
python tools/topology/validate_parallel.py /path/cp2k.psmp /path/parallel --stanene
```

The unit tests compile the same Fortran kernel as CP2K, compare known trivial and nontrivial BHZ
models with Z2Pack, and check non-Abelian gauge changes, reversed loops, singular links, boundary
degeneracy, unresolved surfaces, nonorthogonal-cell input, and `.mmn` ordering/closure validation.

The end-to-end script compares a helium DFT loop against Z2Pack and tests native surface refinement.
With `--soc`, it compares native/Z2Pack SOC surfaces for neon. The old/new regular-connection
comparison is diagnostic only: it exposes a difference in the legacy AO shortcut, not a criterion
for correctness. Separate tests enforce cross-k adjoint symmetry, zero-step normalization, loop
reversal, shifted atomic centres in a skew cell, and full-basis SOC spectra against CP2K's existing
band-structure implementation. The Stanen scripts test native refinement, independent Z2Pack
analysis of the same matrices (raw and polar), and a separate adaptive Z2Pack-driven DFT surface.
The parallel script compares neon across 1x1, 2x1, 2x2, and 4x1 MPI/OpenMP layouts with a full
scalar SOC basis, runs the adaptive Z2Pack interface on two MPI ranks, and optionally validates the
full nontrivial stanene surface on two ranks. It requires `mpiexec` and sufficient execution slots
for four ranks. All inputs and logs are retained; use fresh result directories for new runs.

The fast helium/neon tests are registered in `tests/QS/regtest-topology`. The SOC implementation
also fixes an uninitialized screened projector buffer in `core_ppnl.F`; the Hermiticity tolerance
has not been relaxed or hidden by post-hoc symmetrization.

Still required before production material predictions: basis/virtual-space, cutoff and SCF-mesh
convergence, broader CP2K regression coverage, and scaling tests. This build is validated with the
MPI/OpenMP layouts listed above; 18/18 grid and topology regressions pass with two ranks and two
threads. The previous GCC 16 bounds-check crash in replicated grid transfers is fixed using local
pointer aliases, with bounds checking still enabled; see STATUS.md. Multiple k-point groups are not
used for the explicit export. Successful BHZ tests alone would validate only the mathematics, not
the DFT SOC implementation.

## References

- [Z2Pack overlap interface](https://z2pack.greschd.ch/en/latest/reference/other_systems.html)
- [Gresch et al., Phys. Rev. B 95, 075146 (2017)](https://doi.org/10.1103/PhysRevB.95.075146)
- [Soluyanov and Vanderbilt, Phys. Rev. B 83, 235401 (2011)](https://doi.org/10.1103/PhysRevB.83.235401)
