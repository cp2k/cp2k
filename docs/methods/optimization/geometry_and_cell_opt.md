# Geometry and cell optimization

Geometry optimization relaxes atomic positions by repeatedly evaluating the energy and forces and
moving the atoms toward a stationary point on the potential-energy surface. It is commonly used to
remove artificial forces from a starting structure, obtain a local minimum, prepare structures for
molecular dynamics or property calculations, and, with specialised settings, search for transition
states.

Cell optimization is the related variable-cell problem: the atomic positions and the simulation cell
are relaxed together so that the structure is compatible with a target external pressure or stress
condition.

In CP2K, fixed-cell geometry optimization is selected with
[`RUN_TYPE GEO_OPT`](#CP2K_INPUT.GLOBAL.RUN_TYPE) and controlled by
[`MOTION / GEO_OPT`](#CP2K_INPUT.MOTION.GEO_OPT). If the cell volume or shape should also be
relaxed, use [`RUN_TYPE CELL_OPT`](#CP2K_INPUT.GLOBAL.RUN_TYPE) and
[`MOTION / CELL_OPT`](#CP2K_INPUT.MOTION.CELL_OPT).

```{note}
Some keywords, such as `OPTIMIZER`, `MAX_ITER`, `MAX_DR`, `RMS_DR`, `MAX_FORCE`, and
`RMS_FORCE`, are used by both [`MOTION / GEO_OPT`](#CP2K_INPUT.MOTION.GEO_OPT) and
[`MOTION / CELL_OPT`](#CP2K_INPUT.MOTION.CELL_OPT) with the same meaning. They are linked only once
below unless the variable-cell case needs an additional setting.
```

## GEO_OPT or CELL_OPT?

Use fixed-cell `GEO_OPT` when the lattice vectors are known, intentionally fixed, or irrelevant to
the problem. This is the usual choice for molecules and clusters in a large box, slabs with a fixed
substrate or vacuum region, calculations on a previously optimized crystal cell, and workflows where
only internal atomic coordinates should relax.

Use `CELL_OPT` when the equilibrium cell volume, shape, pressure, or residual stress is part of the
question. It is most common for bulk crystals, strained solids, two-dimensional materials with
relaxed in-plane lattice constants, and structures prepared for fixed-cell molecular dynamics. A
cell optimization is not a replacement for finite-temperature pressure sampling; if the desired
quantity is a thermal average at finite temperature, use an appropriate NPT molecular-dynamics
workflow instead.

`CELL_OPT` has different operating modes selected by [`TYPE`](#CP2K_INPUT.MOTION.CELL_OPT.TYPE):

- `DIRECT_CELL_OPT` optimizes atomic positions and cell vectors in the same optimization loop. The
  stress tensor is computed at every step and this is the usual direct variable-cell optimization
  mode.
- `GEO_OPT` performs an inner geometry optimization between cell-optimization steps. In this mode,
  the [`MOTION / GEO_OPT`](#CP2K_INPUT.MOTION.GEO_OPT) section must also be defined. It can be
  useful when the atoms should be well relaxed before each cell update, but it is more expensive.
- `MD` uses an MD run to compute the stress tensor used for cell optimization. This is a specialised
  workflow and is not the normal choice for zero-temperature structural relaxation.

## Basic setup

Every optimization needs a force-evaluation method in [`FORCE_EVAL`](#CP2K_INPUT.FORCE_EVAL). For a
fixed-cell optimization, set `RUN_TYPE GEO_OPT` and provide a
[`GEO_OPT`](#CP2K_INPUT.MOTION.GEO_OPT) block:

```none
&GLOBAL
  PROJECT my_geo_opt
  RUN_TYPE GEO_OPT
&END GLOBAL

&FORCE_EVAL
  ! Define the method that provides energies and forces.
  ! This can be Quickstep DFT, a force field, a machine-learning potential, QM/MM, ...
&END FORCE_EVAL

&MOTION
  &GEO_OPT
    TYPE MINIMIZATION
    OPTIMIZER BFGS
    MAX_ITER 200
    MAX_DR    3.0E-03
    RMS_DR    1.5E-03
    MAX_FORCE 4.5E-04
    RMS_FORCE 3.0E-04
  &END GEO_OPT
&END MOTION
```

For a variable-cell optimization, use `RUN_TYPE CELL_OPT` and provide a
[`CELL_OPT`](#CP2K_INPUT.MOTION.CELL_OPT) block. Set the intended
[`EXTERNAL_PRESSURE`](#CP2K_INPUT.MOTION.CELL_OPT.EXTERNAL_PRESSURE) explicitly rather than relying
on defaults:

```none
&GLOBAL
  PROJECT my_cell_opt
  RUN_TYPE CELL_OPT
&END GLOBAL

&FORCE_EVAL
  ! Define the method that provides energies, forces, and stress.
&END FORCE_EVAL

&MOTION
  &CELL_OPT
    TYPE DIRECT_CELL_OPT
    OPTIMIZER BFGS
    MAX_ITER 200
    EXTERNAL_PRESSURE 1.01325E+00
    PRESSURE_TOLERANCE 100.0
    MAX_DR    3.0E-03
    RMS_DR    1.5E-03
    MAX_FORCE 4.5E-04
    RMS_FORCE 3.0E-04
  &END CELL_OPT
&END MOTION
```

For `GEO_OPT`, [`TYPE`](#CP2K_INPUT.MOTION.GEO_OPT.TYPE) selects the optimization target.
`MINIMIZATION` searches for a local minimum and is the normal choice. `TRANSITION_STATE` activates
transition-state search machinery; it requires additional method-specific settings in
[`MOTION / GEO_OPT / TRANSITION_STATE`](#CP2K_INPUT.MOTION.GEO_OPT.TRANSITION_STATE) and should not
be treated as an ordinary minimization.

[`MAX_ITER`](#CP2K_INPUT.MOTION.GEO_OPT.MAX_ITER) limits the number of optimization iterations. One
optimization iteration can require more than one force evaluation, depending on the optimizer, line
search, and the selected `CELL_OPT` mode.

## Starting structure and cell

Start from a chemically and physically sensible structure. Severe close contacts, unrealistic
coordination environments, or a badly prepared cell can lead to unstable forces, SCF failures,
unphysical atom displacements, or optimizer steps that are difficult to recover from.

For `CELL_OPT`, the initial cell matters as much as the initial coordinates. The starting volume and
shape should be close enough to the expected structure that the pressure and stress are not
dominated by preparation artefacts. For slabs, two-dimensional materials, or systems with vacuum, do
not relax the vacuum direction unless this is physically intended; constrain the appropriate cell
components or use a fixed-cell optimization after choosing the desired cell.

For variable-cell optimizations, a conventional cell representation is recommended for readability
and easier troubleshooting, while non-standard cell orientations can usually still be optimized.

## Choosing an optimizer

The [`OPTIMIZER`](#CP2K_INPUT.MOTION.GEO_OPT.OPTIMIZER) keyword selects the algorithm used to update
the geometry and, for `CELL_OPT`, the cell. The most useful optimizer-specific controls are the
trust radius for quasi-Newton methods and the line-search settings for conjugate gradients.

- `BFGS` is the default and is usually efficient for small and medium-sized systems. It builds an
  approximate Hessian and can converge rapidly when the starting structure is already reasonable. If
  the optimization takes overly large steps, oscillates, or becomes unstable, the
  [`TRUST_RADIUS`](#CP2K_INPUT.MOTION.GEO_OPT.BFGS.TRUST_RADIUS) in the `BFGS` subsection is often
  the first parameter to reduce. A smaller trust radius makes the optimization more conservative; an
  excessively small value can make convergence slow.
- `LBFGS` is the limited-memory variant and is usually more suitable for large systems where storing
  a full approximate Hessian would be expensive. It also has a
  [`TRUST_RADIUS`](#CP2K_INPUT.MOTION.GEO_OPT.LBFGS.TRUST_RADIUS) control. This can be useful for
  difficult relaxations, although the default behaviour is often adequate for well-prepared large
  systems.
- `CG` is a robust conjugate-gradient minimizer. It is often a good fallback for poor starting
  structures, noisy forces or stresses, or systems where quasi-Newton steps are unstable. Each CG
  optimization step involves a line search along the search direction, which can require several
  force evaluations and make `CG` more expensive than `BFGS`.

For `CG`, the one-dimensional minimization along each search direction is controlled by
[`CG / LINE_SEARCH`](#CP2K_INPUT.MOTION.GEO_OPT.CG.LINE_SEARCH). The
[`TYPE`](#CP2K_INPUT.MOTION.GEO_OPT.CG.LINE_SEARCH.TYPE) keyword selects the line-search algorithm:

- `2PNT` is the cheapest option and extrapolates from two points. It can be efficient when the force
  evaluations are smooth and the optimization path is well behaved. The maximum line-search step can
  be limited with
  [`MAX_ALLOWED_STEP`](#CP2K_INPUT.MOTION.GEO_OPT.CG.LINE_SEARCH.2PNT.MAX_ALLOWED_STEP).
- `GOLD` performs a golden-section/Brent-style line search. It is more expensive but more robust and
  is a safer default for difficult minimizations. The initial bracketing step is controlled by
  [`INITIAL_STEP`](#CP2K_INPUT.MOTION.GEO_OPT.CG.LINE_SEARCH.GOLD.INITIAL_STEP), which may need to
  be reduced for structures with close contacts or unstable early steps.
- `FIT` uses a fitted one-dimensional model. It is the most robust against numerical noise, but also
  the most expensive, and is mainly useful for difficult cases where cheaper line searches behave
  poorly.

If an optimization takes unphysically large steps or becomes unstable, first check the starting
structure or cell and the force/stress quality. Then consider reducing the `BFGS`/`LBFGS` trust
radius, constraining problematic cell degrees of freedom, switching to `CG`, or using a short
preliminary optimization with looser settings before tightening the convergence criteria.

## Convergence criteria

The most important atomic convergence criteria are
[`MAX_FORCE`](#CP2K_INPUT.MOTION.GEO_OPT.MAX_FORCE),
[`RMS_FORCE`](#CP2K_INPUT.MOTION.GEO_OPT.RMS_FORCE), [`MAX_DR`](#CP2K_INPUT.MOTION.GEO_OPT.MAX_DR),
and [`RMS_DR`](#CP2K_INPUT.MOTION.GEO_OPT.RMS_DR). Force criteria are given in Hartree/Bohr, while
displacement criteria are given in Bohr. For `CELL_OPT`, the pressure target and pressure
convergence are controlled by [`EXTERNAL_PRESSURE`](#CP2K_INPUT.MOTION.CELL_OPT.EXTERNAL_PRESSURE)
and [`PRESSURE_TOLERANCE`](#CP2K_INPUT.MOTION.CELL_OPT.PRESSURE_TOLERANCE). The optimization is
considered converged only when all active criteria are satisfied.

When `CELL_OPT` uses `TYPE GEO_OPT`, both the inner `GEO_OPT` criteria and the outer `CELL_OPT`
criteria should be chosen consistently.

A typical convergence report for a variable-cell optimization looks like:

```none
 OPT| **************************************************************************
 OPT| Step number                                                              1
 OPT| Optimization method                                                   BFGS
 OPT| Total energy [hartree]                                      -45.5348295392
 OPT| Internal pressure [bar]                                   47536.1640430316
 OPT| Effective energy change [hartree]                            -0.0002840779
 OPT| Predicted energy change [hartree]                            -0.0001679733
 OPT| Step size                                                     0.0105821341
 OPT| Trust radius                                                  0.3779452266
 OPT| Decrease in energy                                                     YES
 OPT|
 OPT| Maximum step size                                             0.0105821341
 OPT| Convergence limit for maximum step size                       0.0030000000
 OPT| Maximum step size is converged                                          NO
 OPT|
 OPT| RMS step size                                                 0.0033463738
 OPT| Convergence limit for RMS step size                           0.0015000000
 OPT| RMS step size is converged                                              NO
 OPT|
 OPT| Maximum gradient                                              0.0073208554
 OPT| Convergence limit for maximum gradient                        0.0004500000
 OPT| Maximum gradient is converged                                           NO
 OPT|
 OPT| RMS gradient                                                  0.0023150598
 OPT| Convergence limit for RMS gradient                            0.0003000000
 OPT| RMS gradient is converged                                               NO
 OPT|
 OPT| Pressure deviation [bar]                                  47535.1507930316
 OPT| Pressure tolerance [bar]                                    100.0000000000
 OPT| Pressure is converged                                                   NO
 OPT| **************************************************************************
```

In the output, forces are reported as gradients.

Geometry and cell optimization normally try to lower the energy, but convergence is determined by
the active force, displacement, and pressure criteria rather than by the total energy alone.
Intermediate steps may sometimes increase the energy. A very small energy increase at the final step
is usually acceptable if all active convergence criteria are satisfied. In some cases, however, it
may indicate that the structure is not a real minimum, which can be checked by vibrational analysis
where imaginary frequencies reveal unstable modes.

## Constraints, cell degrees of freedom, and symmetry

Atomic constraints are defined in [`MOTION / CONSTRAINT`](#CP2K_INPUT.MOTION.CONSTRAINT). A common
example is [`FIXED_ATOMS`](#CP2K_INPUT.MOTION.CONSTRAINT.FIXED_ATOMS), which keeps selected atoms or
Cartesian components fixed:

```none
&MOTION
  &CONSTRAINT
    &FIXED_ATOMS
      COMPONENTS_TO_FIX XYZ
      LIST 1 2 3
    &END FIXED_ATOMS
  &END CONSTRAINT
&END MOTION
```

[`COMPONENTS_TO_FIX`](#CP2K_INPUT.MOTION.CONSTRAINT.FIXED_ATOMS.COMPONENTS_TO_FIX) selects the
Cartesian components to constrain, and [`LIST`](#CP2K_INPUT.MOTION.CONSTRAINT.FIXED_ATOMS.LIST)
contains atom indices in the order used in [`COORD`](#CP2K_INPUT.FORCE_EVAL.SUBSYS.COORD).
Constraints are useful, for example, for fixing the bottom layers of a slab, holding part of a
support structure, or imposing a chemically intended restriction. Check them carefully before
production calculations: accidental constraints are a common reason for apparently converged but
physically wrong structures.

For `CELL_OPT`, also decide which cell degrees of freedom should be optimized. Useful keywords
include [`CONSTRAINT`](#CP2K_INPUT.MOTION.CELL_OPT.CONSTRAINT) for fixing selected cell components,
[`KEEP_ANGLES`](#CP2K_INPUT.MOTION.CELL_OPT.KEEP_ANGLES) for keeping cell angles fixed,
[`KEEP_VOLUME`](#CP2K_INPUT.MOTION.CELL_OPT.KEEP_VOLUME) for fixed-volume cell-shape relaxation, and
[`KEEP_SYMMETRY`](#CP2K_INPUT.MOTION.CELL_OPT.KEEP_SYMMETRY) or
[`KEEP_SPACE_GROUP`](#CP2K_INPUT.MOTION.CELL_OPT.KEEP_SPACE_GROUP) for symmetry-constrained cell
optimization. For surfaces or layered systems, it is often better to relax only the physically
meaningful cell directions rather than allowing the full cell to change.

Symmetry constraints should be used only when the intended minimum is known to preserve that
symmetry; otherwise they can prevent the structure from relaxing to a lower-symmetry minimum. For
periodic fixed-cell `GEO_OPT`, see [`KEEP_SPACE_GROUP`](#CP2K_INPUT.MOTION.GEO_OPT.KEEP_SPACE_GROUP)
and related symmetry keywords.

## Electronic-structure settings

### SCF quality

For ab-initio geometry or cell optimization, every optimization step requires an
electronic-structure calculation. The force quality depends on the SCF convergence, grid settings,
basis set, pseudopotentials, and method-specific thresholds. For `CELL_OPT`, the stress tensor and
pressure must also be accurate enough for reliable cell updates. Important Quickstep settings
include [`EPS_DEFAULT`](#CP2K_INPUT.FORCE_EVAL.DFT.QS.EPS_DEFAULT),
[`CUTOFF`](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF),
[`REL_CUTOFF`](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF), and
[`EPS_SCF`](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.EPS_SCF).

Optimization is driven by forces, and `CELL_OPT` is additionally driven by stress, rather than by
total energies alone. Therefore, do not choose `EPS_SCF`, `CUTOFF`, `REL_CUTOFF`, or `EPS_DEFAULT`
only from a final static-energy test. These settings usually need to be at least as strict as, and
often stricter than, those used for routine single-point energy calculations, because noisy or
poorly converged forces and stresses can slow down the optimization, cause oscillations, produce
line-search failures, or lead to an unreliable structure or cell.

### Extrapolation

Consecutive optimization steps are often close to each other, so extrapolating the electronic
initial guess can greatly reduce the number of SCF iterations. The relevant Quickstep keywords are
[`EXTRAPOLATION`](#CP2K_INPUT.FORCE_EVAL.DFT.QS.EXTRAPOLATION) and
[`EXTRAPOLATION_ORDER`](#CP2K_INPUT.FORCE_EVAL.DFT.QS.EXTRAPOLATION_ORDER).

`GEXT_PROJ` is usually a strong first choice when applicable, while `ASPC` is also highly
recommended as a robust general option. Other high-order schemes such as `PS` and `GEXT_PROJ_QTR`
can also be worth considering. Simpler previous-step guesses, especially `USE_PREV_WF`, are reliable
fallbacks when a more conservative initial guess is preferred, and `USE_GUESS` can be used when a
fresh initial guess is needed.

Special care is needed with density-matrix-based guesses such as `USE_PREV_P` and `LINEAR_P`,
especially when the cell or neighbour lists change. This is particularly relevant for `CELL_OPT`,
but it can also appear in fixed-cell workflows after restarts or large geometry changes. CP2K may
print a warning such as

```none
*** WARNING in qs_wf_history_methods.F:849 :: Change in cell ***
*** neighborlist: might affect quality of initial guess      ***
```

This indicates that the previous density matrix may not provide a high-quality initial guess. If
this warning appears together with abnormal changes in the geometry, cell, total energy, forces,
stress, or SCF behaviour, stop using that density-matrix-based method and switch to a safer
alternative.

## Output files and restarts

The main output file contains the optimization progress and convergence checks. CP2K also writes
geometry and restart files whose names are based on
[`PROJECT_NAME`](#CP2K_INPUT.GLOBAL.PROJECT_NAME) (or `PROJECT`):

- `project-pos-1.xyz` contains the sequence of geometries visited during the optimization. The last
  frame is the latest optimized geometry.
- Cell information and stress output, when requested, are controlled by print keys such as
  [`MOTION / PRINT / CELL`](#CP2K_INPUT.MOTION.PRINT.CELL) and
  [`MOTION / PRINT / STRESS`](#CP2K_INPUT.MOTION.PRINT.STRESS).
- `project-1.restart` is a CP2K restart input containing the latest geometry, cell, and relevant
  settings.
- `project-1.restart.bak-*` are backup restart files from previous optimization steps.

If a job stops before convergence, continue from the latest restart file, for example:

```none
cp2k -i project-1.restart -o project-restart.out
```

For expensive electronic-structure calculations, wavefunction restart files can reduce the cost of
continuation or follow-up calculations; see
[`FORCE_EVAL / DFT / SCF / PRINT / RESTART`](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.PRINT.RESTART). For BFGS
optimizations, Hessian-restart options are available in the `BFGS` subsection of the active
optimization driver.

## Practical checklist

Before trusting an optimized structure or cell, it's suggested to check that:

- All active force and displacement convergence criteria are satisfied;
- For `CELL_OPT`, the pressure criterion is satisfied and the final cell is physically sensible;
- The optimized cell corresponds to the intended problem: fixed cell for `GEO_OPT`, relaxed cell for
  `CELL_OPT`, and constrained directions where needed;
- No unintended constraints, fixed atoms, cell constraints, or symmetry restrictions were active;
- The SCF procedure converged reliably in the last few optimization steps;
- The final structure and cell were re-used consistently for subsequent calculations.
