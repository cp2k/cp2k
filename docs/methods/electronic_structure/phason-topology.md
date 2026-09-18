# Native topology across sliding geometries

The `topology_phasons` executable evaluates physical cross-geometry overlaps and first/second Chern
pairings directly from CP2K `STATE_EXPORT` snapshots. It uses CP2K's Gaussian `cossin` integrals and
native LAPACK kernels. Neither Python, Z2Pack nor Wannier90 is a runtime dependency. Build with the
CMake target `topology_phasons`; the executable suffix follows the CP2K build, for example
`topology_phasons.psmp`.

The SCF calculations and the topology analysis are separate. Generate converged electronic states at
the prescribed geometries, export the requested post-SCF k-points, then analyze these files.
Changing the analysis mesh must not change the underlying SCF potential at a fixed geometry. The
selected bands, spin convention, cell and atom count must remain consistent along the family.

## Workflow

1. Choose a physically closed parameter family, a selected isolated subspace and explicit atom maps
   at its periodic boundaries. Two parameters give C1. Four parameters, for example
   `(kx, ky, phi_x, phi_y)`, give the second Chern character pairing C2. A varying-atom-count twist
   scan is not such a family.

1. Enable `STATE_EXPORT T` in `DFT/PRINT/WANNIER90` with `NNKP` or `WILSON` points. Keep the
   exported version-1 files, including all computed eigenvalues.

1. Write a mesh manifest as specified below and run, using one process:

   ```sh
   OMP_NUM_THREADS=4 topology_phasons.psmp mesh coarse.mesh > coarse.log
   ```

1. Repeat with every axis refined by an integer factor of at least two. To check a pair of meshes
   against an absolute invariant tolerance:

   ```sh
   OMP_NUM_THREADS=4 topology_phasons.psmp converge coarse.mesh fine.mesh 0.01
   ```

   This requires compatible dimensions, ranks, reciprocal periods and initial physical subspaces,
   and checks both the change and the finer estimate's distance to an integer. The estimates
   themselves are never rounded. Supply the same physical family on both meshes; matching metadata
   alone cannot prove that all intermediate SCF states represent the intended branch. Additional
   refinements, shifted grids and ordinary electronic-structure convergence remain important.
   Numerical convergence does not certify a bulk spectral gap.

## Physical mesh manifest, version 1

The following records appear in order. There are no inline comments:

1. `CP2K_PHASON_MESH 1`
1. Parameter dimension (2 or 4), selected rank.
1. Number of intervals along each axis, each at least 3.
1. Positive metric tolerance, direct-gap tolerance (hartree), minimum link singular value, maximum
   principal plaquette phase (radians, less than pi). Recommended starting values:
   `1e-7 1e-7 1e-8 1.5707963267948966`.
1. Three fractional reciprocal coordinates for the origin.
1. One three-component reciprocal-period vector per parameter axis. Components must be integers; a
   pure displacement axis has vector `0 0 0`. The geometry files, not these vectors, define the
   corresponding displacement.
1. Atom count and logical common-energy-reference assertion, e.g. `1 F`.
1. One atom permutation per axis, using one-based indices. End atom i corresponds to start atom
   permutation(i), modulo periodic lattice translations. Each map must be bijective and the declared
   maps must commute.
1. One frame record for every index `0..n_1, ..., 0..n_d`, including endpoints, with **the first
   index varying fastest**. A record contains a quoted snapshot filename, its one-based k-point
   index and an additive energy offset in hartree. Use absolute paths for portable invocation from
   another working directory.

For a mixed `(kx, phi_x)` mesh, the reciprocal-period rows are `1 0 0` and `0 0 0`. At an endpoint
kx=1, the manifest may reference the stored kx=0 state: the native evaluator retains the reciprocal
seam phase. At a displacement endpoint, supply the corresponding endpoint geometry/state, not an
arbitrary copy chosen only to force closure.

At every used frame the code checks AO-metric normalization and sampled selected/excluded-band
separation. At a seam it verifies the atom mapping, contracted basis, geometry modulo allowed
lattice translations and closure of the physical selected subspace. The endpoint polar sewing matrix
is included in the link product. Independent eigenvector gauges are allowed.

`common_energy_reference=T` asserts that the supplied offsets put spectra on a justified common
scale. For a lowest-band prefix the code then reports `min(E[N+1])-max(E[N])`. Otherwise it reports
the direct isolation and explicitly omits the indirect-gap interpretation. It never independently
zeros each Fermi energy. The global-gap assertion is not required to define an isolated band bundle.

## Output and conventions

- `CHERN_RAW`: unrounded C1 or C2, selected by the parameter dimension.
- `MINIMUM_LINK_SINGULAR_VALUE`: conditioning before polar unitarization.
- `MAXIMUM_PLAQUETTE_PHASE`: phase-resolution diagnostic.
- Metric, seam, direct-gap and, when justified, indirect-gap diagnostics.
- `REFINEMENT_CHANGE` and `INTEGER_RESIDUAL` in convergence mode. An unsuccessful check terminates
  with nonzero status rather than printing a converged invariant.

The physical link is `C_a^dagger O_ab C_b`. `O_ab` includes displaced Gaussian centres, periodic
images, the nonorthogonal metric and Bloch/reciprocal-seam phases. Both SOC spinor components are
summed. No Hamiltonian, SOC term or magnetic field is added by the postprocessor.

For links directed from x to x+mu the plaquette is
`U_mu(x) U_nu(x+mu) U_mu(x+nu)^dagger U_nu(x)^dagger`. C1 is minus the summed principal determinant
phases divided by `2*pi`, matching CP2K's native Wilson winding convention. In four dimensions the
full matrix logarithms F are retained:

`C2 = -sum Tr(F12 F34 - F13 F24 + F14 F23)/(4*pi^2)`.

F is anti-Hermitian and includes the plaquette area. This is the physics second Chern **character
pairing**, not Z2 and not generally the second Chern class. The estimator is not exactly quantized
on a finite mesh. Unresolved phases, singular links and nonunitary plaquettes are rejected.

## Reference links and individual overlaps

For independent model checks, `topology_phasons.psmp links model.links` reads:
`CP2K_TOPOLOGY_LINKS 1`, dimension/rank, mesh intervals, the same four tolerances, then complex raw
matrices (real/imaginary pairs) in vertex, direction, column-major-matrix order. Vertices have no
repeated endpoints and the first mesh coordinate varies fastest. Physical sewing must already be
included in these supplied links. `converge-links coarse.links fine.links tolerance` compares two
such model meshes. This mode cannot validate a model's basis or geometry.

`topology_phasons.psmp overlap left.topology 1 right.topology 1 0 0 0` prints a single physical
overlap. The final three integers shift the right reciprocal coordinate. Frame normalization and
sampled isolation are checked first.

## Memory, parallelism and validation

Screened atom-row contractions use OpenMP and never construct the full dense AO operator. At most a
few selected-state frames and rank-sized working matrices are retained. Link matrices are streamed
through an anonymous scratch file and plaquette curvatures are processed one vertex at a time, not
stored as a dense four-dimensional tensor. Scratch disk use scales as
`16 * rank^2 * dimension * product(mesh)` bytes. Thread-local rank-sized matrices and text snapshot
parsing can still be expensive for large occupied spaces. Use one process, not `mpiexec -n >1`; this
auxiliary program does not distribute its file traversal over MPI ranks. No persistent link cache is
trusted or reused.

The registered native unit tests cover moving normalized Gaussians, adjoint and seam-gauge checks,
singular/underresolved links, C1 determinant phases, C2 gauge covariance and refinement of the
rank-two 4D Dirac model with exact C2=-1. The existing Wilson/C1/Z2 tests remain unchanged. Material
predictions require system-specific basis, SCF, mesh and size convergence in addition to these
tests.

Method references: [Fukui, Hatsugai and Suzuki (2005)](https://doi.org/10.1143/JPSJ.74.1674),
[Mochol-Grzelak et al. (2018)](https://arxiv.org/abs/1803.07003), and
[Rosa, Ruzzene and Prodan (2021)](https://doi.org/10.1038/s42005-021-00630-3).
