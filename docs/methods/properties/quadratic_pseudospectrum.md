# Quadratic pseudospectrum

`FORCE_EVAL / PROPERTIES / QUADRATIC_PSEUDOSPECTRUM` finds states approximately localized in both
position and energy from a converged GPW/GAPW calculation. It reports no topological index and
requires no additional library.

For Hermitian covariant AO observables $A_j$, queries $\lambda_j$, positive weights $w_j$ and
overlap $S$, the generalized eigenproblem is

$$
Q c = q S c,\qquad
Q=\sum_j w_j^2(A_j-\lambda_j S)S^{-1}(A_j-\lambda_j S).
$$

The lowest eigenpairs describe joint approximate states. The output includes $\sqrt{q}$, separate
energy/position residuals and expectation values. States satisfy $c^\dagger S c=1$. Residuals
concern projected finite-basis operators, not exact continuum variances.

```text
&PROPERTIES
  &QUADRATIC_PSEUDOSPECTRUM
    FORMULATION FINITE
    POSITION [bohr] 0 0 0
    ENERGY [hartree] 0.0
    KAPPA [hartree*bohr^-1] 0.01
    NSTATES 2
    &STATES
    &END
  &END
&END
```

## Observables

- `FINITE` uses analytic Cartesian moments. It requires matching `CELL/POISSON PERIODIC NONE`
  without explicit k-points.
- `PERIODIC` uses sine/cosine coordinates on the full `MP_GRID` torus at the frozen SCF potential.
  Only periodic directions are localized. Nonperiodic mesh dimensions must be one.
- `TRANSLATION` uses analytic translated-Gaussian overlaps for $U(a)\psi(r)=\psi(r+a)$. Specify
  repeated `TRANSLATION_VECTOR` and `WAVE_VECTOR` queries and positive `TRANSLATION_SCALE` instead
  of `KAPPA`. The target phase is $\exp(i k\cdot a)$ without an additional $2\pi$. This path
  includes the full continuum translation norm and reports leakage outside the AO span, rather than
  treating compressed AO translation matrices as unitary.

The periodic paths construct a complete analysis torus, independently of SCF symmetry reduction.
`SOC T` uses post-SCF GTH spin-orbit coupling after a restricted calculation with suitable
pseudopotentials.

## Output and limits

`SOLVER DENSE` uses a metric-aware dense reference calculation. `MAX_AO` bounds the complete scalar
AO space, including periodic cells. `EPS_METRIC` rejects ill-conditioned overlap matrices, and
`EPS_EIGEN` bounds eigenpair residuals. `NSTATES` selects the number of states.

The `STATES` print section writes complex AO coefficients. Spin-up precedes spin-down for SOC, with
cell-major ordering inside each spin block on a torus. Individual vectors in a degenerate subspace
are gauge dependent. Compare subspaces, residuals and observables instead of individual
coefficients.

Choose energy, position and scales for the physical question and check basis and size dependence.
Broader examples are available in
[TopologicalCP2K](https://github.com/DCM-Uni-Paderborn/TopologicalCP2K).
