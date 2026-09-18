# Fractional occupation density

Fractional occupation number weighted density (FOD), introduced by Grimme and Hansen
{ref}`Grimme2015`, provides a qualitative indication of static electron correlation and its spatial
distribution. It is not a definitive test for multireference character and does not replace a
multireference calculation.

For each spin channel CP2K uses the converged Fermi-Dirac occupations and chemical potential:

$$
\rho_{\mathrm{FOD}}(\mathbf r)=\sum_{i\sigma}w_{i\sigma}|\psi_{i\sigma}(\mathbf r)|^2,
\qquad
w_{i\sigma}=\begin{cases}
f_{\max}-f_{i\sigma}, & \epsilon_{i\sigma}\leq\mu_\sigma,\\
f_{i\sigma}, & \epsilon_{i\sigma}>\mu_\sigma.
\end{cases}
$$

Here `f_max` is two for a restricted calculation and one for each unrestricted spin channel. At the
chemical potential both expressions give half the maximum occupation. The reported
`N_FOD (orbital sum)` is the sum of these weights, with no real-space quadrature error. The SCF
wavefunction and density are not changed by the analysis.

Enable `FORCE_EVAL/DFT/PRINT/FOD` explicitly. This initial implementation supports GPW at the Gamma
point, with positive-temperature `SCF/SMEAR/METHOD FERMI_DIRAC`. Other smearing schemes, GAPW and
explicit k-point calculations are rejected rather than producing a different quantity.

```text
&SCF
  ADDED_MOS 20
  &DIAGONALIZATION
  &END
  &SMEAR
    METHOD FERMI_DIRAC
    ELECTRONIC_TEMPERATURE [K] 5000
  &END
&END
&PRINT
  &FOD
    CUBE T
    STRIDE 1
  &END
&END
```

This is an input fragment, not a complete calculation. Choose `ADDED_MOS` for the system and check
convergence of `N_FOD` with the number of available orbitals. A still-occupied last orbital produces
a warning. The reference molecular protocol uses TPSS/def2-TZVP and 5000 K. Changing the functional,
temperature, basis or pseudopotential changes the diagnostic; a PBE/MOLOPT calculation is not a
numerical reproduction of that protocol. Metallic thermal occupations must not be interpreted
directly as molecular multireference character.

`CUBE` defaults to false: scalar-only analysis avoids copying orbitals or allocating density grids.
With `CUBE T`, CP2K additionally prints the full-grid integral and writes the nonnegative FOD
density in electrons/bohr cubed. Converge the density grid when comparing that integral to the
orbital sum. `STRIDE` subsamples only the output cube; use `STRIDE 1` to integrate the cube without
additional sampling error. The cube writer is collective and works in serial and MPI.
