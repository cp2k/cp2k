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
[MOTION/GEO_OPT](#CP2K_INPUT.MOTION.GEO_OPT). If the cell volume or shape should also be relaxed,
use [`RUN_TYPE CELL_OPT`](#CP2K_INPUT.GLOBAL.RUN_TYPE) and
[MOTION/CELL_OPT](#CP2K_INPUT.MOTION.CELL_OPT).

```{note}
Some keywords, such as `OPTIMIZER`, `MAX_ITER`, `MAX_DR`, `RMS_DR`, `MAX_FORCE`, and
`RMS_FORCE`, are used by both [MOTION/GEO_OPT](#CP2K_INPUT.MOTION.GEO_OPT) and
[MOTION/CELL_OPT](#CP2K_INPUT.MOTION.CELL_OPT) with the same meaning. They are linked only once
below unless the variable-cell case needs an additional setting.
```

## GEO_OPT or CELL_OPT?

Use fixed-cell `GEO_OPT` when the lattice vectors are known, intentionally fixed, or irrelevant to
the problem. This is the usual choice for molecules and clusters in a large box, slabs with a fixed
substrate or vacuum region, calculations on a previously optimized crystal cell, and workflows where
only internal atomic coordinates should relax.

Use `CELL_OPT` when the equilibrium cell volume, shape, pressure, or residual stress is part of the
question. It is most common for bulk crystals, strained solids, two-dimensional materials with
relaxed in-plane lattice constants, and structures prepared for fixed-cell molecular dynamics.

```{note}
An optimization is not a replacement for finite-temperature pressure sampling; if the desired
quantity is a thermal average at finite temperature, use a molecular-dynamics workflow with
appropriate ensemble instead.

For the usual electronic-structure methods based on the Born-Oppenheimer approximation
(which, by neglecting nuclear motion, provides the concept of "potential-energy surface"
in the first place), there is no involvement of the temperature (kinetic energy of the
nuclei) in an optimization, and the energy that is being optimized does not include
components like zero-point energy and thermal corrections (harmonic or anharmonic).
The notion "optimization at 0 K" is really a misnomer in this regard.
```

## Basic setup

Every optimization needs a force-evaluation method in [FORCE_EVAL](#CP2K_INPUT.FORCE_EVAL). For a
fixed-cell optimization, set `RUN_TYPE GEO_OPT` and provide a [GEO_OPT](#CP2K_INPUT.MOTION.GEO_OPT)
block:

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
    MAX_ITER  200
    MAX_DR    3.0E-03
    RMS_DR    1.5E-03
    MAX_FORCE 4.5E-04
    RMS_FORCE 3.0E-04
  &END GEO_OPT
&END MOTION
```

For a variable-cell optimization, use `RUN_TYPE CELL_OPT` and provide a
[CELL_OPT](#CP2K_INPUT.MOTION.CELL_OPT) block. Set the intended
[EXTERNAL_PRESSURE](#CP2K_INPUT.MOTION.CELL_OPT.EXTERNAL_PRESSURE) explicitly rather than relying on
defaults, and an additional keyword [STRESS_TENSOR](#CP2K_INPUT.FORCE_EVAL.STRESS_TENSOR) makes
sense as well:

```none
&GLOBAL
  PROJECT my_cell_opt
  RUN_TYPE CELL_OPT
&END GLOBAL

&FORCE_EVAL
  ! Define the method that provides energies, forces, and stress.
  STRESS_TENSOR ANALYTICAL ! Compute full stress tensor analytically
&END FORCE_EVAL

&MOTION
  &CELL_OPT
    OPTIMIZER BFGS
    MAX_ITER  200
    EXTERNAL_PRESSURE 1.01325E+00
    PRESSURE_TOLERANCE 100.0
    MAX_DR    3.0E-03
    RMS_DR    1.5E-03
    MAX_FORCE 4.5E-04
    RMS_FORCE 3.0E-04
  &END CELL_OPT
&END MOTION
```

```{note}
The availability of analytical stress tensor depends on the method for evaluating energy and force.
Some electronic-structure methods may only support numerical stress tensor from finite differences,
which would take quite more time at each iteration of optimization. This can be confirmed by setting
[FORCE_EVAL/PRINT/STRESS_TENSOR](#CP2K_INPUT.FORCE_EVAL.PRINT.STRESS_TENSOR) and seeing if numerical
stress is mentioned in the output.
```

For `GEO_OPT`, [TYPE](#CP2K_INPUT.MOTION.GEO_OPT.TYPE) selects the optimization target.
`MINIMIZATION` searches for a local minimum and is the normal choice. `TRANSITION_STATE` activates
transition-state search machinery; it requires additional method-specific settings in
[MOTION/GEO_OPT/TRANSITION_STATE](#CP2K_INPUT.MOTION.GEO_OPT.TRANSITION_STATE) and should not be
treated as an ordinary minimization.

[MAX_ITER](#CP2K_INPUT.MOTION.GEO_OPT.MAX_ITER) limits the number of optimization iterations. One
optimization iteration can require more than one force evaluation, depending on the optimizer, line
search, and the selected `CELL_OPT` mode.

An optimization may terminate if one of the following events takes place:

- The convergence criteria (see below) are satisfied before reaching `MAX_ITER`;
- The `MAX_ITER` is reached, regardless of whether the convergence criteria are satisfied or not;
- Some external control for the program is activated, such as hitting the global walltime as
  specified by [WALLTIME](#CP2K_INPUT.GLOBAL.WALLTIME), or discovering a file in the working
  directory with filename `EXIT` or `EXIT_GEO` that is created by user in case of emergency.

## Starting structure and cell

Start from a chemically and physically sensible structure. Severe close contacts, unrealistic
coordination environments, or a badly prepared cell can lead to unstable forces, SCF failures,
unphysical atom displacements, or optimizer steps that are difficult to recover from. For instance,
"amorphous cell" structures generated by randomly packing molecules and/or atoms in a box can be
extremely challenging for optimization, especially if the chemical bonds between atoms are
significantly scrambled. Same for snapshots from a high-temperature molecular-dynamics simulation in
vaporized or molten conditions with lots of broken chemical bonds.

For `CELL_OPT`, the initial cell matters as much as the initial coordinates. The starting volume and
shape should be close enough to the expected structure and density, such that the pressure and
stress are not dominated by preparation artefacts.

A visualization of the structure and cell in modelling programs, with the box and the neighboring
periodic images displayed, will be very helpful; neither large vacuous gaps nor crowded cluster of
atoms should occur near the boundary of the box on the directions consistent with the periodicity,
and conversely, sufficient vacuum space on the non-periodic directions is crucial for eliminating
unwanted interactions across the boundary.

For surface slabs, two- or one-dimensional materials, or systems with vacuum, do not relax the
vacuum direction unless physically intended; it is better to constrain the appropriate cell
components or use a fixed-cell optimization after choosing the desired cell.

```{warning}
**Do not use experimental structure blindly.**

If available, **computational** materials databases are the most recommended avenue
to obtain structures that are "computation-ready", or even better, already optimized
with some electronic-structure methods. On the other hand, structures that are from
**experimental** characterization are frequently not "computation-ready", and thus
should not be subject to optimization without careful validation in pre-processing.
This can be prominent for `cif` and `pdb` structures determined by powder or single-
crystal XRD which can be affected by sample quality and thermal motion.

- Watch out for crystallographic disorder and atoms with low resolution or fractional
  occupation: using the superposition of all atoms as if every occupancy is 1.00 is
  highly likely to introduce contacting or even overlapping atoms.
- Beware of composition: the atomic structure may not match the intended macroscopic,
  charge-neutral chemical formula, owing to missing or duplicated hydrogen atoms,
  small counter ions, solvent or ligand molecules.

Possible resolutions vary from simple manual editing in the modelling stage, to
utilization of supercells and enumeration of special quasirandom structures (common
for materials with dopants), and to more rigorous XRD refinement and application of
quantum crystallography methods.
```

```{note}
For variable-cell optimizations, the [CELL_OPT](#CP2K_INPUT.MOTION.CELL_OPT) section
clarifies that a cell representation with the first vector A along the X-axis and
the second vector B on the XY plane is expected. This is to ensure that the cell
vectors are updated correctly in response to the external pressure in the current
implementation; even when an unconventional vector definition could sometimes feel
easy (for instance the primitive rhombohedral cell of the face-centered cubic lattice),
an appropriate linear transformation to the cell and the atomic coordinates should
still be applied to create the model for input. For more details, see online
[github discussion](https://github.com/cp2k/cp2k/issues/3384).
```

## Choosing an optimizer

The [OPTIMIZER](#CP2K_INPUT.MOTION.GEO_OPT.OPTIMIZER) keyword selects the algorithm used to update
the geometry and, for `CELL_OPT`, the cell. The most useful optimizer-specific controls are the
trust radius for quasi-Newton methods and the line-search settings for conjugate gradients.

- `BFGS` is the default and is usually efficient for small and medium-sized systems. It builds an
  approximate Hessian and can converge rapidly when the starting structure is already reasonable. If
  the optimization takes overly large steps, oscillates, or becomes unstable, the
  [TRUST_RADIUS](#CP2K_INPUT.MOTION.GEO_OPT.BFGS.TRUST_RADIUS) in the `BFGS` subsection is often the
  first parameter to reduce. A smaller trust radius makes the optimization more conservative; an
  excessively small value can make convergence slow. The number of required iterations can often be
  reduced by enabling the [USE_MODEL_HESSIAN](#CP2K_INPUT.MOTION.GEO_OPT.BFGS.USE_MODEL_HESSIAN)
  option.
- `LBFGS` is the limited-memory variant and is usually more suitable for large systems where storing
  a full approximate Hessian would be expensive. It also has a
  [TRUST_RADIUS](#CP2K_INPUT.MOTION.GEO_OPT.LBFGS.TRUST_RADIUS) control. This can be useful for
  difficult relaxations, although the default behaviour is often adequate for well-prepared large
  systems.
- `CG` is a robust conjugate-gradient minimizer. It is often a good fallback for poor starting
  structures, noisy forces or stresses, or systems where quasi-Newton steps are unstable. Each CG
  optimization step involves a line search along the search direction, which can require several
  force evaluations and make `CG` more expensive than `BFGS`.

For `CG`, the one-dimensional minimization along each search direction is controlled by
[CG/LINE_SEARCH](#CP2K_INPUT.MOTION.GEO_OPT.CG.LINE_SEARCH). The
[TYPE](#CP2K_INPUT.MOTION.GEO_OPT.CG.LINE_SEARCH.TYPE) keyword selects the line-search algorithm:

- `2PNT` is the cheapest option and extrapolates from two points. It can be efficient when the force
  evaluations are smooth and the optimization path is well behaved. The maximum line-search step can
  be limited with
  [MAX_ALLOWED_STEP](#CP2K_INPUT.MOTION.GEO_OPT.CG.LINE_SEARCH.2PNT.MAX_ALLOWED_STEP).
- `GOLD` performs a golden-section/Brent-style line search. It is more expensive but more robust and
  is a safer default for difficult minimizations. The initial bracketing step is controlled by
  [INITIAL_STEP](#CP2K_INPUT.MOTION.GEO_OPT.CG.LINE_SEARCH.GOLD.INITIAL_STEP), which may need to be
  reduced for structures with close contacts or unstable early steps.
- `FIT` uses a fitted one-dimensional model. It is the most robust against numerical noise, but also
  the most expensive, and is mainly useful for difficult cases where cheaper line searches behave
  poorly.

If an optimization takes unphysically large steps or becomes unstable, first check the starting
structure or cell and the force/stress quality. Then consider reducing the `BFGS`/`LBFGS` trust
radius, constraining problematic cell degrees of freedom, switching to `CG`, or using a short
preliminary optimization with looser settings before tightening the convergence criteria.

## Convergence criteria

The most important atomic convergence criteria are
[MAX_FORCE](#CP2K_INPUT.MOTION.GEO_OPT.MAX_FORCE),
[RMS_FORCE](#CP2K_INPUT.MOTION.GEO_OPT.RMS_FORCE), [MAX_DR](#CP2K_INPUT.MOTION.GEO_OPT.MAX_DR), and
[RMS_DR](#CP2K_INPUT.MOTION.GEO_OPT.RMS_DR). Force criteria are given in unit of Hartree/Bohr, while
displacement criteria are given in Bohr. In addition, for `CELL_OPT`, the pressure target and
pressure convergence in bar are specified by
[EXTERNAL_PRESSURE](#CP2K_INPUT.MOTION.CELL_OPT.EXTERNAL_PRESSURE) and
[PRESSURE_TOLERANCE](#CP2K_INPUT.MOTION.CELL_OPT.PRESSURE_TOLERANCE) respectively.

A typical convergence report for a variable-cell optimization, as controlled by
[PROGRAM_RUN_INFO](#CP2K_INPUT.MOTION.CELL_OPT.PRINT.PROGRAM_RUN_INFO), looks like the following
where forces are reported as gradients. For a fixed-cell optimization, the
[PROGRAM_RUN_INFO](#CP2K_INPUT.MOTION.GEO_OPT.PRINT.PROGRAM_RUN_INFO) does not mention the pressure
but is otherwise the same.

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

Satisfying all of these criteria within `MAX_ITER` steps is a sufficient condition for the
convergence of optimization. With the `BFGS` and `CG` optimizers, this is also the necessary
condition. With the `LBFGS` optimizer, however, this is not a necessary condition because `LBFGS`
has its own extra pair of criteria named
[WANTED_PROJ_GRADIENT](#CP2K_INPUT.MOTION.GEO_OPT.LBFGS.WANTED_PROJ_GRADIENT) and
[WANTED_REL_F_ERROR](#CP2K_INPUT.MOTION.GEO_OPT.LBFGS.WANTED_REL_F_ERROR) which may override the
aforementioned general ones. It is possible to see an optimization task with the `LBFGS` optimizer
finish with the message below, even when one or more of the general criteria have not been met yet;
try restarting the optimization task with tightened criteria or alternative optimizers until
satisfactory convergence.

```none
 ***********************************************
 * Specific L-BFGS convergence criteria         
 * WANTED_PROJ_GRADIENT and WANTED_REL_F_ERROR  
 * satisfied .... run CONVERGED!                
 ***********************************************
```

The `LBFGS` optimizer has its own [PRINT_LEVEL](#CP2K_INPUT.MOTION.GEO_OPT.LBFGS.PRINT_LEVEL) which
controls the verbosity of details printed to the main output. In addition, a separate `iterate.dat`
file is also created by `LBFGS` optimizer where the columns `projg` and `f` can be checked against
the `WANTED_PROJ_GRADIENT` and `WANTED_REL_F_ERROR` criteria. Omitting some headers, the main body
of an `iterate.dat` excerpt is reproduced below.

```none
   it   nf  nseg  nact  sub  itls  stepl    tstep     projg        f
    0    1     -     -   -     -     -        -     4.632D-02 -1.082D+00
    1    2     1     0  ---    0  1.0D+00  6.6D-02  5.016D-02 -1.087D+00
    2    4     1     0  ---    1  1.9D+00  1.3D-01  5.707D-02 -1.097D+00
    3    6     1     0  ---    1  1.7D+00  1.3D-01  6.212D-02 -1.108D+00
    4    8     1     0  ---    1  1.5D+00  1.3D-01  6.413D-02 -1.120D+00
    5   10     1     0  ---    1  1.5D+00  1.3D-01  6.112D-02 -1.132D+00
    6   11     3     2  con    0  1.0D+00  1.3D-01  4.979D-02 -1.143D+00
    7   12     3     2  con    0  1.0D+00  1.3D-01  2.440D-02 -1.150D+00
    8   14     1     0  con    1  5.1D-01  6.6D-02  3.685D-03 -1.152D+00
    9   15     1     0  con    0  1.0D+00  1.2D-02  7.569D-04 -1.152D+00
```

Geometry and cell optimization normally try to lower the energy, but convergence is determined by
the active force, displacement, and pressure criteria rather than by the total energy alone.
Intermediate steps may sometimes increase the energy. A very small energy increase at the final step
is usually acceptable if all active convergence criteria are satisfied. In some cases, however, it
may indicate that the structure is not a real minimum, which can be checked by vibrational analysis
where imaginary frequencies reveal unstable modes.

Generally speaking, the usual convergence criteria may be satisfied within dozens or a few hundreds
of iterations. The default `MAX_ITER` value of 200 is mostly sufficient and may be raised to about
300 ~ 400 to deal with difficult cases like massive flexible structures. Running a single-point
calculation with `RUN_TYPE ENERGY_FORCE` beforehand, and then multiplying the elapsed time by
expected iterations, provides an estimate for the total time cost of an optimization task. To
monitor a long optimization, do not focus on these numeric criteria alone; it is recommended to
download the geometry (trajectory) file and observe the evolution of the structure and cell in a
visualizer program occasionally, and if anything goes wrong, halt the program in time.

## Constraints, cell degrees of freedom, and symmetry

Atomic constraints are defined in [MOTION/CONSTRAINT](#CP2K_INPUT.MOTION.CONSTRAINT). A common
example is [FIXED_ATOMS](#CP2K_INPUT.MOTION.CONSTRAINT.FIXED_ATOMS), which keeps selected atoms or
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

[COMPONENTS_TO_FIX](#CP2K_INPUT.MOTION.CONSTRAINT.FIXED_ATOMS.COMPONENTS_TO_FIX) selects the
Cartesian components to constrain, and [LIST](#CP2K_INPUT.MOTION.CONSTRAINT.FIXED_ATOMS.LIST)
contains atom indices in the order used in [COORD](#CP2K_INPUT.FORCE_EVAL.SUBSYS.COORD). Constraints
are useful, for example, for fixing the bottom layers of a slab, holding part of a support
structure, or imposing a chemically intended restriction. Check them carefully before production
calculations: accidental constraints are a common reason for apparently converged but physically
wrong structures. For example, in a model of surface or micropore adsorption (physisorption or
chemisorption), the substrate should be allowed to relax locally rather than fully fixed to account
for the realistic interaction between the substrate (adsorbent) and the adsorbate.

For `CELL_OPT`, also decide which cell degrees of freedom should be optimized. Useful keywords
include [CONSTRAINT](#CP2K_INPUT.MOTION.CELL_OPT.CONSTRAINT) for fixing selected cell components,
[KEEP_ANGLES](#CP2K_INPUT.MOTION.CELL_OPT.KEEP_ANGLES) for keeping cell angles fixed,
[KEEP_VOLUME](#CP2K_INPUT.MOTION.CELL_OPT.KEEP_VOLUME) for fixed-volume cell-shape relaxation, and
[KEEP_SYMMETRY](#CP2K_INPUT.MOTION.CELL_OPT.KEEP_SYMMETRY) or
[KEEP_SPACE_GROUP](#CP2K_INPUT.MOTION.CELL_OPT.KEEP_SPACE_GROUP) for symmetry-constrained cell
optimization. For surfaces or layered systems, it is often better to relax only the physically
meaningful cell directions rather than allowing the full cell to change.

Symmetry constraints should be used only when the intended minimum is known to preserve that
symmetry; otherwise they can prevent the structure from relaxing to a lower-symmetry minimum. For
periodic fixed-cell `GEO_OPT`, see [KEEP_SPACE_GROUP](#CP2K_INPUT.MOTION.GEO_OPT.KEEP_SPACE_GROUP)
and related symmetry keywords.

```{note}
The `KEEP_SYMMETRY` keyword should always be used in conjunction with the cell symmetry
specified by [FORCE_EVAL/SUBSYS/CELL/SYMMETRY](#CP2K_INPUT.FORCE_EVAL.SUBSYS.CELL.SYMMETRY).
These keywords act only on the cell vectors and their gradients in a variable-cell
optimization, and they both assume a conventional cell where A is along the X-axis and
B is on the XY plane as noted in "Starting structure and cell" above.

To allow for the use of `KEEP_SPACE_GROUP`, the `spglib` library should be installed
and linked to the CP2K build in order to detect and preserve the space group. Usually
`KEEP_SPACE_GROUP` is used together with `KEEP_SYMMETRY`.
```

## Electronic-structure settings

Density functional theory (DFT) is an electronic-structure (wavefunction) method routinely used for
optimization. This section elaborate on the relevant aspects, assuming basic knowledge about the
method which can be found at [](../dft/index.md).

### SCF quality

Every optimization step requires an electronic-structure calculation, and the force quality depends
on the SCF convergence, grid settings, basis set, pseudopotentials, and method-specific thresholds.
For `CELL_OPT`, the stress tensor and pressure must also be accurate enough for reliable cell
updates. (The SCF convergence for each iteration is not to be confused with the convergence of
geometry in terms of stepsize, gradient and pressure across iterations in the optimization process.)

Optimization is driven by forces, and `CELL_OPT` is additionally driven by stress, rather than by
total energies alone. Therefore, simply adjusting settings and tuning parameters from a final
static-energy test may not be sufficient. They usually need to be at least as strict as, and often
stricter than, those used for routine single-point energy calculations, because noisy or poorly
converged forces and stresses can slow down the optimization, cause oscillations, produce
line-search failures, or lead to an unreliable structure or cell.

Important Quickstep settings include [CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF),
[REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF),
[EPS_DEFAULT](#CP2K_INPUT.FORCE_EVAL.DFT.QS.EPS_DEFAULT), and
[EPS_SCF](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.EPS_SCF).

```{note}
The mesh spacing of the plane-wave and real-space grid for electronic-structure
calculations is determined both by the `CUTOFF` setting and by the dimensions of
the cell. A variable-cell optimization may witness significant variations of the
cell, which will also affect the number of grid points and introduce artificial
discontinuities to the energy, force, stress and other properties even if the
change in structure is small.

It is recommended for better stability to set up a reference cell in the section
[FORCE_EVAL/SUBSYS/CELL/CELL_REF](#CP2K_INPUT.FORCE_EVAL.SUBSYS.CELL.CELL_REF),
with the dimension of the reference cell larger than the expected upper bound within
reach during the whole variable-cell optimization process, which will be kept constant
so that the number of grid points is held fixed and the change in energy more smooth.
```

### Extrapolation

Consecutive optimization steps are often close to each other, so extrapolating the electronic
initial guess can greatly reduce the number of SCF iterations. The relevant Quickstep keywords are
[EXTRAPOLATION](#CP2K_INPUT.FORCE_EVAL.DFT.QS.EXTRAPOLATION) and
[EXTRAPOLATION_ORDER](#CP2K_INPUT.FORCE_EVAL.DFT.QS.EXTRAPOLATION_ORDER).

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

```{note}
The applicability of extrapolation schemes depends on versions. To exemplify, extrapolation
was formerly available for gamma-only calculations, but the latest versions of CP2K has
started supporting extrapolation with full k-point sampling. However, when k-points are
requested to be reduced by symmetry, the dynamic refreshing of symmetry of k-points alongside
geometry updates could make extrapolation unavailable and automatically fall back to `USE_GUESS`.
```

```{danger}
**Be responsible, and do not ignore SCF convergence failure blindly.**

By default, a failure in SCF convergence aborts the program with a message like:
`SCF run NOT converged. To continue the calculation regardless, please set the keyword
IGNORE_CONVERGENCE_FAILURE.` Setting
[IGNORE_CONVERGENCE_FAILURE](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.IGNORE_CONVERGENCE_FAILURE)
to `.TRUE.` turns it into a warning `SCF run NOT converged` that allows for optimization
to continue. However, bad SCF convergence leads to unreliable energy, force, and stress
on the current step, introducing error to the updated structure and the extrapolated
wavefunction on the next step. An optimization process with most or all steps failing
to converge SCF cycles does not yield trustworthy results in the end.

Therefore, `IGNORE_CONVERGENCE_FAILURE` should only be considered as a **last resort**
out of desperation, rather than a universal remedy used on a regular basis or even as
the default. There are dozens of measures available towards SCF convergence on top of
a realistic, chemically sensible model structure and appropriate net charge, spin
multiplicity, atomic magnetization, and k-point sampling; please try achieving SCF
convergence on the starting structure in a single-point calculation in the first place.
```

## Output files and restarts

The main output file contains the optimization progress and convergence checks. Geometry and restart
files are also written with names based on the global
[PROJECT_NAME](#CP2K_INPUT.GLOBAL.PROJECT_NAME) (or `PROJECT`) as well as relevant printkeys under
the [MOTION/PRINT](#CP2K_INPUT.MOTION.PRINT) section.

- `project-pos-1.xyz` contains the sequence of geometries visited during the optimization. This is
  controlled by the [TRAJECTORY](#CP2K_INPUT.MOTION.PRINT.TRAJECTORY) section with the default
  format of XYZ (or XMOL). Note that the original XYZ specification does not convey periodicity, and
  not all of the visualizers can parse the additional cell information in the comment line of XYZ
  (e.g. in extended XYZ specification); switching to other formats using the `FORMAT` keyword may be
  needed for visualization.
- `project-FINAL-1_{iter}.cif` and `project-FINAL-1_{iter}.xyz` is a CIF and an extended XYZ file
  for the final structure, where `{iter}` is the number of iterations (steps) during optimization
  until reaching convergence or hitting `MAX_ITER`. These are controlled by
  [FINAL_STRUCTURE](#CP2K_INPUT.MOTION.PRINT.FINAL_STRUCTURE).
- `project-1.cell` contains the cell vectors during optimization, as controlled by
  [CELL](#CP2K_INPUT.MOTION.PRINT.CELL).
- `project-1.stress` contains the stress tensors during optimization, as controlled by
  [STRESS](#CP2K_INPUT.MOTION.PRINT.STRESS).
- `project-1.restart` is a CP2K restart input containing the latest geometry, cell, and relevant
  settings, as controlled by [RESTART](#CP2K_INPUT.MOTION.PRINT.RESTART). When `BACKUP_COPIES` under
  this section is non-zero, backup restart files from previous optimization steps will also be
  preserved as `project-1.restart.bak-*`.

The restart file can be used directly as a new input file like

```none
cp2k -i project-1.restart -o project-restart.out
```

Alternatively, the [EXT_RESTART](#CP2K_INPUT.EXT_RESTART) section can be used for finer restart
controls. These will be convenient in case an optimization does not converge and a new optimization
starting from the latest structure is needed.

For expensive electronic-structure calculations, wavefunction restart files can reduce the cost of
continuation or follow-up calculations; see
[FORCE_EVAL/DFT/SCF/PRINT/RESTART](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.PRINT.RESTART). For BFGS
optimizations, Hessian-restart options are available in the `BFGS` subsection of the active
optimization driver.

## Practical checklist

Before trusting an optimized structure or cell, it is suggested to check that:

- All active force and displacement convergence criteria are satisfied;
- For `CELL_OPT`, the pressure criterion is satisfied and the final cell is physically sensible;
- The optimized cell corresponds to the intended problem: fixed cell for `GEO_OPT`, relaxed cell for
  `CELL_OPT`, and constrained directions where needed;
- No unintended constraints, fixed atoms, cell constraints, or symmetry restrictions were active;
- The SCF procedure converged reliably in the last few optimization steps;
- The final structure and cell were re-used consistently for subsequent calculations.
