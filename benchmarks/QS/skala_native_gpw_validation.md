# Native GPW/SKALA validation notes

These checks document the current validation scope of the experimental CP2K-native GPW/SKALA grid
path. They are not part of the standard regression suite because several finite-difference force
checks are too expensive for routine CI, but they are useful for reproducing the branch-level
validation.

Scope:

- Gamma-only `METHOD GPW` with GTH pseudopotentials.
- `GAUXC%FUNCTIONAL PBE`, `GAUXC%MODEL SKALA`, and `GAUXC%NATIVE_GRID T`.
- Energy, VXC matrix construction, and nuclear gradients.
- No k-points, no stress tensor, no ADMM, no ROKS, no NLCC pseudopotentials, and no native-grid
  GAPW.

Representative force finite-difference checks from the Spark validation build:

| System | Spin | MPI x OpenMP | Cutoff/rel_cutoff [Ry] | Displacement [bohr] | Force difference [Ha/bohr] |
| ------ | ---- | ------------ | ---------------------- | ------------------- | -------------------------- |
| H2     | RKS  | 2 x 2        | 150/30                 | 1.0e-4              | 4.23e-6                    |
| OH     | UKS  | 2 x 2        | 250/40                 | 2.0e-4              | 7.6e-7                     |
| NH3    | RKS  | 2 x 2        | 250/40                 | 2.0e-4              | 5.77e-6                    |
| H2O    | RKS  | 2 x 2        | 250/40                 | 2.0e-4              | 3.00e-5                    |

The standard `QS/regtest-gauxc` block covers the compact H2 force check and the energy/feature
diagnostics. The longer OH, NH3, and H2O checks should be kept as external validation unless their
runtime can be reduced substantially.

GAPW status:

- The all-electron GAPW/SKALA path remains the molecular GauXC AO-matrix path.
- The native GPW grid path intentionally aborts for `METHOD GAPW` because a valid extension needs a
  consistent reconstructed all-electron density, gradient, and kinetic-energy density on the SKALA
  feature grid.
- `METHOD GAPW_XC` remains outside the current scope.
