# Electrostatics and Poisson Solvers

The electrostatic boundary conditions are part of the physical model, not just a numerical setting.
They determine which periodic images interact, whether a charged calculation needs a compensating
charge, and which total-energy differences are meaningful. Set them explicitly before converging the
cell size.

## Periodicity settings

Two settings have to be distinguished:

- [CELL/PERIODIC](#CP2K_INPUT.FORCE_EVAL.SUBSYS.CELL.PERIODIC) controls the periodicity used for the
  atomic system, including geometry and pair lists.
- [POISSON/PERIODIC](#CP2K_INPUT.FORCE_EVAL.DFT.POISSON.PERIODIC) selects the periodic directions of
  the electrostatic Green's function.

They should normally describe the same physical system. CP2K permits different values because some
specialized workflows need them, but a mismatch can make different parts of a calculation use
different images. Do not use a mismatch merely to suppress interactions in a vacuum direction.

The generalized [IMPLICIT](#CP2K_INPUT.FORCE_EVAL.DFT.POISSON.IMPLICIT) solver has an additional,
independent [BOUNDARY_CONDITIONS](#CP2K_INPUT.FORCE_EVAL.DFT.POISSON.IMPLICIT.BOUNDARY_CONDITIONS)
keyword. Its default is `PERIODIC`, so `POISSON_SOLVER IMPLICIT` does not become a 1D or 2D solver
solely because `POISSON/PERIODIC` was set to `X`, `XY`, or another reduced periodicity. Configure
the `IMPLICIT` subsection explicitly when nonperiodic, mixed, or dielectric boundary conditions are
required.

For example, a slab that is periodic in the $x$ and $y$ directions can use the analytic 2D Green's
function as follows:

```
&SUBSYS
  &CELL
    ABC 10.0 10.0 20.0
    PERIODIC XY
  &END CELL
  # ...
&END SUBSYS

&DFT
  &POISSON
    PERIODIC XY
    POISSON_SOLVER ANALYTIC
  &END POISSON
  # ...
&END DFT
```

## Solver overview

Choose [POISSON_SOLVER](#CP2K_INPUT.FORCE_EVAL.DFT.POISSON.POISSON_SOLVER) from the intended
periodicity and boundary conditions. The input reference is the authoritative list of currently
supported combinations.

| Solver      | Supported periodicity         | Main considerations                                                                                                                                                                                                                                                                                 |
| ----------- | ----------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `PERIODIC`  | 3D                            | Standard fully periodic reciprocal-space solution.                                                                                                                                                                                                                                                  |
| `ANALYTIC`  | 0D, 1D, 2D                    | Analytic reciprocal-space Green's functions; convergence can be slow.                                                                                                                                                                                                                               |
| `MT`        | 0D, 2D                        | Martyna--Tuckerman image decoupling. The cell must be at least twice the extent of the complete electronic density. [](#Martyna1999)                                                                                                                                                                |
| `MULTIPOLE` | 0D                            | Fits the total charge with atom-centred Gaussians to decouple periodic images. [](#Bl%C3%B6chl1995)                                                                                                                                                                                                 |
| `WAVELET`   | 0D, 2D, 3D                    | Requires FFTSG-compatible grid dimensions. The default [EXTENDED_FFT_LENGTHS](#CP2K_INPUT.GLOBAL.EXTENDED_FFT_LENGTHS) setting, `.FALSE.`, guarantees this independently of the selected FFT library. The density must decay to zero at nonperiodic cell faces. [](#Genovese2006) [](#Genovese2007) |
| `IMPLICIT`  | Generalized 0D--3D conditions | Supports periodic, Neumann, and mixed Dirichlet/periodic conditions and spatially varying dielectrics. Configure its boundary conditions explicitly. [](#BaniHashemian2016)                                                                                                                         |

```{warning}
`MT` can return completely wrong results when the electronic density reaches beyond half the cell.
The relevant extent includes diffuse density tails, not only the nuclear coordinates.
```

For an isolated, localized molecule or cluster, `MT`, `WAVELET`, or `MULTIPOLE` are common 0D
choices; `ANALYTIC` and `IMPLICIT` provide further options. A slab normally needs a genuine 2D
solver or explicitly constructed mixed boundary conditions. For a wire, use a solver that supports
1D periodicity. Fully periodic bulk calculations normally use `PERIODIC`.

## Systems with a net charge

The electrostatic energy of a non-neutral periodic cell is not defined without an additional
convention. With the 3D `PERIODIC` solver, CP2K uses a homogeneous compensating background and omits
the singular zero reciprocal-space component. The electronic density, ionic charge, and background
therefore form one electrostatic problem; separate electron--background or ion--background energies
are not independently meaningful.

Increasing the cell volume reduces the background density but also changes the interaction with
periodic images. Consequently, a charged 3D calculation has a cell-size- and cell-shape-dependent
finite-size error that usually decays only algebraically. CP2K does not infer a universal physical
correction merely from [CHARGE](#CP2K_INPUT.FORCE_EVAL.DFT.CHARGE). Compare energies only between
calculations that use the same boundary-condition and potential-zero convention, or apply a
correction appropriate to the physical model.

An isolated 0D solver can give a finite energy for a localized charged cluster without a periodic
background, provided that the density is contained in the numerical cell. The situation differs for
partially periodic systems: a net charge per unit length or area produces a long-ranged field.
Without an explicit countercharge, electrode, electrolyte, or other compensating model, the energy
of a charged 1D wire or 2D slab generally has no vacuum-independent isolated limit. Increasing the
vacuum can then change the reported energy rather than converge it.

```{important}
Total energies from different Poisson solvers are not automatically on the same reference for a
charged system. A constant shift of the electrostatic potential is irrelevant for a neutral cell
but changes the energy assigned to a cell with nonzero net charge. Solver-to-solver agreement is
therefore not, by itself, a valid correctness test.
```

Implicit solvent, planar countercharges, and electrode models can make a charged interface a
well-defined physical problem, but they describe different systems. Select such a model from the
experimental or thermodynamic boundary conditions rather than as a numerical correction applied
afterwards.

## Converging vacuum and cell size

Use a convergence series that preserves the physical model:

1. Choose the periodic directions, Poisson solver, charge-compensation model, and potential
   reference first.
1. Converge [CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) and
   [REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF) separately.
1. For a 0D system, enlarge every nonperiodic direction. For a slab or wire, enlarge only the
   nonperiodic direction or directions. Keep the structure centred and keep the basis, functional,
   SCF thresholds, and periodic cell vectors fixed.
1. Record energy differences relative to the largest cell together with forces and, where relevant,
   the planar-averaged potential and charge density. Check that the density is negligible at every
   nonperiodic boundary required by the solver.
1. Repeat at a higher grid cutoff if changes in the FFT grid could be confused with a vacuum-size
   trend.

For a neutral localized system, the total energy should reach a plateau once image interactions,
density truncation, and grid changes are below the target accuracy. For a charged 3D supercell,
analyse the finite-size trend under the homogeneous-background convention. For a charged 1D or 2D
system, do not interpret the absence of a plateau as a numerical failure until the electrostatic
model contains the physically required compensation.

Diffuse basis functions can create states concentrated in the vacuum, especially for charged
surfaces. Inspect orbitals and density when adding vacuum or diffuse functions: changing the cell
may otherwise change the electronic state as well as its electrostatic finite-size error.
