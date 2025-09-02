# Constrained Nuclear-Electronic Orbital DFT

## Introduction

Constrained nuclear-electronic orbital density functional theory (CNEO-DFT)[^1] treats some or all
nuclei quantum mechanically while maintaining well-defined molecular geometries/crystal structures.
The method applies position constraints to quantum nuclei, creating effective potential energy
surfaces that incorporate nuclear quantum effects, particularly zero-point energies due to quantum
delocalization. Standard geometry optimization, vibrational analysis and classical molecular
dynamics can then be performed on these effective potential energy surfaces in place of conventional
Born-Oppenheimer surfaces.

## Theory

CNEO-DFT applies position constraints to quantum nuclear densities within the multicomponent DFT
framework. The energy functional depends on both the electronic density
$\rho^{\text{e}}(\mathbf{r})$ and nuclear densities $\rho^{\text{n}}_I(\mathbf{r})$:

$$
\begin{multline}
E[ \rho^{\text{e}}, \{ \rho^{\text{n}}_I \} ] = T_{\text{s}}^{\text{e}} [ \rho^{\text{e}} ] + \sum_I T_{\text{s}}^{\text{n}, I} [ \rho^{\text{n}}_I ] + \int \text{d} \mathbf{r} V_{\text{ext}} (\mathbf{r}) \left[ \rho^{\text{e}} (\mathbf{r}) - \sum_I Z_I \rho^{\text{n}}_I (\mathbf{r}) \right] \\ + E_{\text{H}} [ \rho^{\text{e}}, \{ \rho^{\text{n}}_I \} ] + E_{\text{xc}}^{\text{e}} [ \rho^{\text{e}} ] + E_{\text{c}} [ \rho^{\text{e}}, \{ \rho^{\text{n}}_I \} ]
\end{multline}
$$

The position constraint for each quantum nucleus $I$ is:

$$
\langle \mathbf{r} \rangle_I = \int \mathbf{r} \rho^{\text{n}}_I (\mathbf{r}) \text{d} \mathbf{r} = \mathbf{R}_I
$$

where $\mathbf{R}_I$ is the position expectation value corresponding to classical molecular/crystal
geometries. This constraint is enforced through Lagrange multipliers $\mathbf{f}_I$ in the CNEO
nuclear Kohn-Sham equation:

$$
\left[ -\frac{1}{2M_I}\nabla^2 + v_{\text{eff}}^{\text{n}, I} + \mathbf{f}_I \cdot (\mathbf{r} - \mathbf{R}_I) \right] \phi_i^{\text{n}, I} = \varepsilon_i^{\text{n}, I} \phi_i^{\text{n}, I}
$$

Self-consistent solution at each nuclear geometry generates a point on the CNEO effective energy
surface, which can be interpreted as a constrained minimized energy surface (CMES). The CMES
theoretical framework[^2] provides the theoretical justification for performing dynamics
propagations on these effective energy surfaces.

Analytic gradients of the CNEO energy with respect to quantum nuclear position expectation values as
well as with respect to classical nuclear positions provide the forces needed for geometry
optimization, vibrational analysis, and molecular dynamics simulations.

## How to Use Periodic CNEO-DFT

### Prerequisites

[Periodic CNEO-DFT](#Chen2025) requires GAPW for atoms containing quantum nuclei. Users should be
familiar with GAPW setup in CP2K.

### Basic Input Structure

In addition to electronic basis set files, nuclear basis set files also need to be specified. For
example, using `NUCLEAR_BASIS_SETS` for PB series protonic basis sets.[^3]

```
&DFT
  BASIS_SET_FILE_NAME NUCLEAR_BASIS_SETS  ! Nuclear basis
  &QS
    METHOD GAPW               ! GAPW is required for CNEO
  &END QS
&END DFT
```

Users can also provide other basis files with custom nuclear basis functions.

### Defining Quantum Nuclei

Atoms with quantum nuclei are specified by setting `POTENTIAL CNEO`, which activates the CNEO-DFT
treatment for those atoms. These atoms require both electronic and nuclear basis sets: the
electronic basis must be suitable for all-electron calculations to properly describe
electron-nuclear interactions in the core region, while nuclear basis functions describe the quantum
nuclear orbitals:

```
&KIND H
  BASIS_SET Ahlrichs-def2-TZVP   ! All-electron electronic basis
  BASIS_SET NUC PB4-D            ! Nuclear basis
  POTENTIAL CNEO                 ! This enables CNEO-DFT calculation
&END KIND
```

### Classical Nuclei

Classical nuclei use standard setup with either pseudopotentials or all-electron treatments. Users
may choose to skip GAPW treatment for atoms with classical nuclei by setting `GPW_TYPE .TRUE.` as
only atoms with quantum nuclei strictly require GAPW.

### Nuclear Basis Sets

PB series[^3] (PB4-D to PB6-H) protonic basis sets are provided in `NUCLEAR_BASIS_SETS`. It also
contains two versions with scaled exponents as examples for quantum nuclei other than protons:

- **PB4-D** to **PB6-H**: Optimized basis functions for protons (with spherical harmonics)
- **PB4-D_D**: Scaled exponents of PB4-D for deuterium
- **PB4-D_Mu**: Scaled exponents of PB4-D for muonium (for testing purposes)

Although scaled versions for basis sets other than `PB4-D` are not provided, users can manually
scale the exponents by $\sqrt{m_{\text{deuteron}}/m_{\text{proton}}}$ if they wish to use those
basis sets to study deuterium. However, this simple modification may not be optimal. For quantum
nuclei with general masses and charges, basis development remains very limited, and one may need to
construct basis functions on their own (e.g., even-tempered basis).

### Mass Specification

For hydrogen atoms, the program automatically assumes $^1\text{H}$ isotope and uses the
corresponding nuclear mass (proton mass = 1.007825 u minus electron mass). For isotopes or other
quantum nuclei, specify the neutral atom mass explicitly:

```
&KIND H
  BASIS_SET Ahlrichs-def2-TZVP
  BASIS_SET NUC PB4-D_D          ! Nuclear basis suitable for the mass
  POTENTIAL CNEO
  MASS 2.01410177811             ! Neutral deuterium atom mass
&END KIND
```

The program automatically subtracts the electron mass to obtain the nuclear mass. For elements
beyond hydrogen, use pure isotope masses rather than average atomic masses. Note that the mean-field
treatment of electron-nuclear interactions has been benchmarked primarily for hydrogen systems,[^4]
though preliminary full-quantum studies without correlation exist.[^5] The lack of optimized nuclear
basis sets is another possible issue.

### K-Point Sampling

Electronic degrees of freedom support standard k-point sampling for extended systems. K-point
sampling affects only the electronic subsystem, while quantum nuclei, due to their highly-localized
picture, remain described by localized basis functions and are independent of the k-point grid.

## When to Use CNEO-DFT

**Hydrogen-containing systems** where proton quantum effects significantly influence molecular
properties. Key applications include:

- **Vibrational spectroscopy** for improved frequency predictions that capture anharmonic
  contributions from quantum nuclear motion
- **Isotope effect studies** to understand how nuclear mass changes (H/D) affect molecular
  properties through altered zero-point energies and quantum delocalization
- **Surface chemistry** involving adsorption, diffusion, or transfer processes where nuclear quantum
  effects alter thermodynamic and kinetic behavior
- **Reaction dynamics** where zero-point energy and nuclear shallow tunneling influence barrier
  heights and reaction rates

[^1]: https://doi.org/10.1063/1.5143371

[^2]: https://doi.org/10.1021/acs.jpclett.2c02905

[^3]: https://doi.org/10.1063/5.0009233

[^4]: https://doi.org/10.1063/5.0243086

[^5]: https://doi.org/10.1063/5.0014001
