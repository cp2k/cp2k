# Spectral localizers

`FORCE_EVAL / PROPERTIES / SPECTRAL_LOCALIZER` evaluates a local topological index and protection
gap from the converged GPW/GAPW Hamiltonian and AO overlap. It is independent of Wannier90 and Kubo
transport. No additional library is required.

## Finite systems

For `CELL PERIODIC NONE` and `POISSON PERIODIC NONE`, the operators are the analytic Cartesian AO
moments. For the ordered plane $(X,Y)$, the covariant localizer is

$$
L=\begin{pmatrix}
H-ES & \kappa[(X-xS)-i(Y-yS)]\\
\kappa[(X-xS)+i(Y-yS)] & -(H-ES)
\end{pmatrix},\qquad B=\operatorname{diag}(S,S).
$$

Here $H$, $S$, $X$ and $Y$ are AO matrices, $(x,y)$ and $E$ are the query position and energy, and
$\kappa$ has units of energy per length. For class A, the index is half the signature of $L$. The
gap is the smallest absolute generalized eigenvalue of $(L,B)$, not a pivot magnitude. LAPACK
pivoted LDL factorization is checked against the eigenvalue result.

```text
&PROPERTIES
  &SPECTRAL_LOCALIZER
    FORMULATION FINITE
    INVARIANT CHERN
    PLANE XY
    POSITION [bohr] 0 0 0
    ENERGY [hartree] 0.0
    KAPPA [hartree*bohr^-1] 0.01
  &END
&END
```

Choose the energy and scale for the system. Repeated positions and lists of energies and scales
produce scans. A real scalar Hamiltonian is not expected to have a nonzero class-A Chern index.

## Spin-orbit coupling and periodic systems

`SOC T` adds GTH spin-orbit coupling after a restricted SCF calculation in the complete AO spinor
space. It needs SOC pseudopotential parameters and is not self-consistent noncollinear DFT.
`INVARIANT Z2` additionally requires `TIME_REVERSAL T`. The code checks this symmetry and uses an
orientation-sensitive Pfaffian sign in two dimensions. `DIMENSION 3` uses the class-AII chiral-block
determinant construction and all three position directions.

`FORMULATION PERIODIC` uses analytic sine/cosine AO integrals on a complete Born-von-Karman torus.
Set `MP_GRID` to its replication factors and use `ETA` (energy units) instead of `KAPPA`.
Nonperiodic directions must have size one. For two dimensions, `PLANE` must agree with the two
periodic directions in both `CELL` and `POISSON`. Three dimensions require `PERIODIC XYZ`. The torus
is built at the frozen SCF potential, independently of the SCF mesh and its symmetry reduction.
Finite Cartesian coordinates are not substituted for periodic position operators.

## Numerical controls

`SOLVER DENSE` is the reference implementation. `MAX_AO` bounds the complete scalar AO space
including torus cells. Spinors double that space, and dense localizers require quadratic storage.
`EPS_METRIC`, `MATRIX_TOLERANCE` and `EPS_GAP` control metric conditioning, operator checks and gap
resolution. Unresolved gaps do not yield an index. `SPECTRAL_FLATTENING` optionally replaces the
shifted Hamiltonian by its metric-covariant matrix sign after checking an electronic gap. Its scale
and gap are distinct from the final localizer gap.

`DFT / PRINT / AO_MATRICES / POSITION` exports the finite AO moments in the same ordering as overlap
and Hamiltonian matrices. `SOC` exports the three real antisymmetric GTH components, whose physical
orbital operators include a factor of $i$. Both exports are restricted to isolated systems.

These are finite-basis local diagnostics. Check stability against basis, system size, query position
and scale before interpreting a material index. Broader examples are available in
[TopologicalCP2K](https://github.com/DCM-Uni-Paderborn/TopologicalCP2K).

## Related AO transport operators

`KUBO_TRANSPORT / CURRENT_OPERATOR PROJECTED_AO` uses the finite projected Hamiltonian commutator
with analytic AO positions. `BLOCH` uses analytic Hamiltonian/overlap derivatives and the AO
connection on an independent full `MP_GRID`. `SOC` adds post-SCF GTH spin-orbit coupling after a
restricted SCF. `SYMMETRY` reconstructs property eigenframes and checks current covariance.

`HALL_RESPONSE` additionally evaluates the antisymmetric DC charge response with the same positive
scalar dissipation. It introduces neither magnetic order nor microscopic scattering, and a
time-reversal-symmetric Hamiltonian has zero net charge Hall response. `MAX_AO` and `MAX_MEMORY_MB`
bound dense property storage. The atom-embedding current remains the default.
