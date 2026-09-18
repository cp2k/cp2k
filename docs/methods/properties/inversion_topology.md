# Native inversion representations and symmetry indicators

CP2K can evaluate inversion representations directly from Gaussian-basis Bloch
states, including second-variational spin-orbit coupling (SOC). The analysis does
not require Z2Pack, Wannier90 execution, or maximally localized Wannier functions.
The input resides in `DFT/PRINT/WANNIER90`, alongside the native Wilson-loop
analysis and overlap export.

This implementation covers the **inversion subgroup**, with spinful
time-reversal symmetry for the topological classification. It is not an
automatic irreducible-representation or elementary-band-representation (EBR)
catalogue for all space groups. Additional crystal symmetries are not used to
classify the bands.

## Input

After a converged complex-k-point GPW calculation, request:

```text
&WANNIER90
  KPOINTS_SOURCE TRIM
  SOC T
  TIME_REVERSAL T
  INVERSION_TQC T
  TQC_DIMENSION 3
  PARITY_ORIGIN 0.0 0.0 0.0
  ! EXCLUDE_BANDS must exclude all unoccupied spinor bands.
&END WANNIER90
```

`KPOINTS_SOURCE TRIM` generates the eight time-reversal-invariant momenta (TRIM)
in three dimensions. `TQC_DIMENSION 2` instead generates the four TRIM of the
fractional `kz=0` plane. The SCF potential is held fixed during the additional
diagonalizations. A slab still uses a three-dimensionally periodic cell.

`PARITY_ORIGIN` is the inversion center in fractional direct-lattice coordinates
in the input cell setting. It is not inferred from the structure. For example,
inversion about `(1/2,1/2,1/2)` has translation `(1,1,1)` relative to inversion
about the origin. These operations can have different individual parities at
zone-boundary TRIM, so the setting must accompany reported symmetry data.

`INVERSION_TQC` requires an even number of occupied spinor states, selected as
the lowest consecutive bands. The number must equal the electron count. Compute
excluded bands as well so that a selected/excluded spectral separation can be
checked. Converge the scalar-state basis used for second-variational SOC.
`TIME_REVERSAL T` asserts physical time-reversal symmetry of the Hamiltonian.
The code additionally checks its action on the selected subspace at each TRIM.

`PARITY T` without `INVERSION_TQC` is also available for isolated scalar or
spinor subspaces. It reports even and odd inversion-irrep multiplicities,
counting individual states rather than Kramers pairs. The `TRIM` point source
implicitly enables this diagnostic. Alternatively, `NNKP` or `WILSON` point
sources can be used, provided they include the required TRIM. The `TRIM` source
itself does not define Wilson loops and cannot be combined with `WILSON_LOOP T`.

`PARITY_TOLERANCE` defaults to `1e-6`. It controls dimensionless representation
residuals and the atom-mapping distance in Bohr. `PARITY_ENERGY_TOL` defaults to
`1e-6` Hartree and controls projected symmetry/eigenvalue commutators.
`WILSON_GAP_TOL` sets the minimum selected/excluded separation at sampled points.
Fractional occupations or a metallic occupied manifold do not define the
insulating classification described here.

## Metric-aware inversion representations

For a basis orbital of angular momentum $l$ on atom $a$, let

$$2\mathbf c-\mathbf r_a=\mathbf r_b+\mathbf L$$

with fractional coordinates and integer lattice translation $\mathbf L$.
Inversion maps that orbital to the same Gaussian function on atom $b$ with
multiplier $(-1)^l\exp(2\pi i\mathbf k\cdot\mathbf L)$. Atom mapping requires
identical CP2K atomic kinds and a unique inversion partner. Spatial inversion
does not rotate spin. The resulting sparse coefficient-space operator is $P$.

For selected Bloch coefficients $C$ with $C^\dagger S C=I$, the representation is

$$D(P)=C^\dagger S P C.$$

For spinors, $S$ is repeated on the two spin blocks. CP2K checks metric
covariance, subspace normalization, Hermiticity and unitarity of $D(P)$, and
commutation with the selected eigenvalues. This treats arbitrary mixtures
within degenerate eigenspaces without assigning a parity to a gauge-dependent
individual vector. A native character-decomposition kernel yields the even and
odd multiplicities from $\chi(E)$ and $\chi(P)$. That kernel also accepts a
complete unitary character table supplied per group element, but automatic
general little-group tables and AO rotation representations are not provided.

For spinful time reversal, the sewing matrix of $\Theta=i\sigma_y K$ is checked
for unitarity, skew symmetry and consistency with the eigenvalues. Each parity
multiplicity must be even before converting from states to Kramers pairs.

## Indicators and atomic signatures

Let $n^-_\kappa$ be the number of odd occupied Kramers pairs at TRIM
$\kappa\in\{0,1\}^d$, where $\mathbf k=\kappa/2$ and the first reciprocal
coordinate is the fastest-changing bit. The Fu--Kane parity index is

$$\nu=\sum_\kappa n^-_\kappa\pmod 2.$$

In three dimensions this is the strong index. The weak index $\nu_i$ is the
same sum restricted to $\kappa_i=1$. CP2K defines its inversion indicator as

$$z_4=\sum_\kappa n^-_\kappa\pmod 4.$$

The sign convention is stated explicitly because other conventions exchange
$z_4=1$ and $z_4=3$. In two dimensions only the plane $\mathbb Z_2$ index is
reported. A bulk insulating gap and physical time-reversal symmetry are
prerequisites for interpreting these as insulating topological indices.

The inversion-subgroup EBR catalogue is generated analytically. One local
Kramers pair of parity $s=\pm1$ at an inversion center $\mathbf a/2$ has
character $s(-1)^{\mathbf a\cdot\kappa}$. For $N$ occupied Kramers pairs,
define $t_\kappa=N-2n^-_\kappa$. Its discrete Fourier coefficients are

$$d_\mathbf a=2^{-d}\sum_\kappa
  (-1)^{\mathbf a\cdot\kappa}t_\kappa.$$

The signature belongs to the **integer lattice of atomic signatures** if all
$d_\mathbf a$ are integers and $N-\sum_\mathbf a d_\mathbf a$ is even. It admits
a **nonnegative atomic decomposition** if additionally
$N\geq\sum_\mathbf a|d_\mathbf a|$. These tests use exact integer arithmetic.
A nonnegative decomposition is printed as even/odd local Kramers-pair counts
at each center relative to `PARITY_ORIGIN`. A general inversion-related pair
of sites has a signature already represented by the sum of even and odd
centered EBRs, so it adds no new signature generator.

An atomic-compatible signature is not a proof of trivial topology. An integer
but not nonnegative decomposition is an EBR-signature obstruction; calling it
fragile requires excluding stable topology invisible to these symmetry data.
A nonzero stable indicator is likewise conditional on band isolation and the
assumed symmetries. TRIM checks alone do not demonstrate a gap throughout the
Brillouin zone. Wilson loops provide complementary information.

## Tests and implementation boundaries

`topology_symmetry_unittest` tests character decomposition including a
two-dimensional irrep, metric-aware scalar/spinor inversion, degenerate-state
gauge changes, rejection of broken inversion and time reversal, all individual
inversion-center EBRs, strong and inversion-$z_4$ examples, and a signed-only
signature. It independently enumerates all rank-two atomic sums and compares
the solver for all 81 two-dimensional and 6,561 three-dimensional candidate
parity signatures, including indicator/atomic-lattice consistency.

`QS/regtest-topology/neon-tqc.inp` and `neon-tqc-plane.inp` exercise the complete
SCF-to-parity/EBR path with SOC in three and two dimensions and with different
inversion origins. They are atomic-limit regression tests, not tests of a
topologically nontrivial material.

`stanene-tqc.inp` uses the fixed buckled-stanene benchmark geometry and SOC
basis to obtain the nontrivial plane parity index, complementing Wilson-loop
validation. These end-to-end tests have been run with two MPI ranks and two
OpenMP threads per rank. Stanene was also checked with one rank and one thread.
The existing single-rank, multiple-thread local k-point diagonalization path
can produce non-normalized states with the tested GNU 16/OpenBLAS 0.3.33/DBCSR
2.10.0 build. The new metric checks reject those states. Use one OpenMP thread
per rank for single-rank runs until that independent threading issue is
resolved; this extension does not claim to fix it.

Automatic space-group identification, nonsymmorphic little-group sewing
matrices, general double-group irreps, compatibility graphs, magnetic groups
and a complete space-group EBR catalogue are outside this implementation.

## References

- L. Fu and C. L. Kane, *Topological insulators with inversion symmetry*,
  [Phys. Rev. B **76**, 045302 (2007)](https://doi.org/10.1103/PhysRevB.76.045302).
- B. Bradlyn et al., *Topological quantum chemistry*,
  [Nature **547**, 298--305 (2017)](https://doi.org/10.1038/nature23268).
- H. C. Po, A. Vishwanath and H. Watanabe, *Symmetry-based indicators of band
  topology in the 230 space groups*,
  [Nature Communications **8**, 50 (2017)](https://doi.org/10.1038/s41467-017-00133-2).
- J. Cano and B. Bradlyn, *Band representations and topological quantum chemistry*,
  [Annual Review of Condensed Matter Physics **12**, 225--246 (2021)](https://arxiv.org/abs/2006.04890).
