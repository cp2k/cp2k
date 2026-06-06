# Molecular Dynamics

Molecular dynamics (MD) simulates how atoms and molecules move in time. It is widely used to sample
finite-temperature behaviour, compute thermodynamic averages, and analyse dynamical properties such
as diffusion, vibrations, structural fluctuations, and reaction events. MD trajectories are also a
common starting point for more advanced simulation and analysis techniques.

In CP2K, the forces used for MD are provided by a selected `FORCE_EVAL` method. The same MD driver
can therefore be used with classical force fields, machine-learning potentials, QM/MM,
semi-empirical methods, and Born-Oppenheimer ab-initio molecular dynamics (AIMD).

```{youtube} cGDiVSWpXLc
---
url_parameters: ?start=2
align: center
privacy_mode:
---
```

## Basic setup

Select an MD calculation with [`RUN_TYPE MD`](#CP2K_INPUT.GLOBAL.RUN_TYPE). The main MD controls are
in [MOTION/MD](#CP2K_INPUT.MOTION.MD): [ENSEMBLE](#CP2K_INPUT.MOTION.MD.ENSEMBLE) selects the
propagator, [STEPS](#CP2K_INPUT.MOTION.MD.STEPS) or [MAX_STEPS](#CP2K_INPUT.MOTION.MD.MAX_STEPS)
sets the length of the run, [TIMESTEP](#CP2K_INPUT.MOTION.MD.TIMESTEP) sets the integration step,
and [TEMPERATURE](#CP2K_INPUT.MOTION.MD.TEMPERATURE) is used for velocity initialisation and as the
constant-temperature target.

A minimal fixed-energy MD block is:

```none
&GLOBAL
  RUN_TYPE MD
  PROJECT_NAME my_md
&END GLOBAL

&MOTION
  &MD
    ENSEMBLE NVE
    STEPS 10000
    TIMESTEP 0.5
    TEMPERATURE 300
  &END MD
&END MOTION
```

The force method is defined independently in [FORCE_EVAL](#CP2K_INPUT.FORCE_EVAL).

## Preliminary considerations

Universal to every atomistic simulation and every program, a number of crucial points must be
carefully considered before setting up a molecular dynamics task and investing copious amounts of
computational resources. As a starter, this section discusses preliminary factors in the choice of
the spatiotemporal scale of the trajectory, periodic boundary condition. force method, and the
initial structure.

### Space and time scale

Even when taking the ergodic hypothesis as granted, statistical analysis on target properties at
thermodynamic equilibrium has to be performed on top of sufficient sampling. The (auto)correlation
length and (auto)correlation time provide estimates for the space and time scale where measurements
of some property retain a memory of past states and are not fully statistically independent. It is
important to learn about their typical values in the same or similar systems, and take adequate
multiples as the size of the system and the length of the trajectory.

On the other hand, non-equilibrium irreversible processes such as phase transitions and chemical
reactions are also observable at some time and space resolution. Experimentally it is determined by
the sensitivity of instrumental characterization and analysis methods, with certain precise
definitions of the states before and after the process. Unfortunately, many experiments take place
in real life in a macroscopic time interval that exceeds the capability of simulation; for instance,
a chemical reaction taking days, hours, or even just seconds to happen may not be replicatable in
AIMD simulations with trajectory on the picosecond or nanosecond scale starting with only a bunch of
reactants randomly packed together, however thermodynamically favorable it is. Again, it is better
to survey literatures with relevant kinetic data and have a clear understanding of the typical free
energy barriers and reaction rates before proceeding to simulations. Even static thermochemistry
calculations for the profile of energy landscape via optimization of minima and transition states,
vibrational analysis and evaluation of single-point energy would be helpful. To simulate difficult
reactions by molecular dynamics in a short trajectory, advanced _enhanced sampling_ methods such as
[](./metadynamics) may be necessary.

### Periodic boundary condition

**Periodic boundary condition** ({term}`PBC`) may or may not be used in the molecular dynamics
simulation, and there are caveats for both ends.

On one hand, very small periodic cells for production MD should be avoided unless the finite-size
effects are known to be acceptable. Small cells produce artificial correlations through periodic
images, exaggerate statistical temperature and pressure fluctuations, and can make NPT cell dynamics
unstable or unphysical. Sometimes using a supercell of a crystalline structure is preferred if its
primitive cell is too small.

On the other hand, without PBC in MD, a common risk is that strong repulsion or very weak attraction
causes the system to expand and effectively diminishes desired interactions between molecules. A
possible remedy is to define an external potential or a constraint that act as a soft wall to
reflect or push stray molecules back to the center, but this is artificial and arbitrary.

### Force method

To propagate the trajectory in a molecular dynamics simulation, the energy and force is evaluated on
each timestep to provide acceleration via Newton's second law. The choice of the force method
determines the description of chemical bonding and/or weak interactions, the configuration space
that can be sampled, and the eventual computational cost. Careful selection of function forms and
parameters as well as thorough validation against the target property is vital to any application;
in particular, the part of "cross" interactions between two or more components described with more
than one unified force method requires very well-founded justification.

The strength and limitation of each force method is decided by design and parametrization. For
instance, classical force fields in molecular mechanics usually model chemical bonds by harmonic
potentials on the stretching, bending and rigid dihedral terms (and periodic torsion potentials on
the flexible dihedral terms). The simplicity of function forms allows for speedy evaluation and the
conformational fluctuation or conversion can be captured well, but the bonding pattern is fixed to
near equilibrium and bond breaking and forming in chemical reactions cannot be captured. To describe
the bond rearrangement and the underlying electronic structure, AIMD based on quantum chemical
methods becomes requisite even though it is slower and more costly. The middle grounds have recently
becoming filled with emerging methods like the machine-learning potentials.

### Initial structure

Because most of the points apply here as well, please first refer to the documentation about the
[starting structure and cell](../optimization/geometry_and_cell_opt.md#starting-structure-and-cell)
of a geometry optimization.

In general, start from a physically reasonable structure with sensible bond lengths, intermolecular
distances, density, and cell parameters. A majority of the force methods demonstrate extremely
strong repulsion at very small nuclei distances, and thus molecular dynamics simulation employing
these methods is not a reliable way to fix severe close contacts, unrealistic densities, or badly
prepared cells, because such problems can cause unstable forces, SCF failures, or immediate heating
and atom ejection.

In fact, when starting from scratch, a geometry optimization preceding molecular dynamics simulation
is highly recommended. During geometry optimization, the structure is relaxed and the maximum force
on atoms is significantly lowered, which puts an upper limit to the initial acceleration in the
subsequent molecular dynamics simulation and prevents rapid atomic motion and extreme temperature
rise. With electronic-structure methods, the geometry optimization also offers a chance to confirm
SCF convergence and estimate the time consumption on each step. The only "downside" is that the
temperature may start from or drop to a relatively low value at early simulation due to the small
forces and slow motion, and may take a little longer to heat up to the desired temperature; see the
section on equilibration below for how to address this.

## Initial velocities

If velocities are not provided by a restart file or an external trajectory, CP2K can initialise
atomic velocities from [TEMPERATURE](#CP2K_INPUT.MOTION.MD.TEMPERATURE). The
[INITIALIZATION_METHOD](#CP2K_INPUT.MOTION.MD.INITIALIZATION_METHOD) keyword controls how the
initial positions and velocities are generated.

For isolated molecules, clusters, slabs, or other systems where global translation or rotation is
not part of the intended dynamics, it can be useful to remove centre-of-mass translation and, when
meaningful, overall angular motion; see [COMVEL_TOL](#CP2K_INPUT.MOTION.MD.COMVEL_TOL),
[ANGVEL_TOL](#CP2K_INPUT.MOTION.MD.ANGVEL_TOL), [ANGVEL_ZERO](#CP2K_INPUT.MOTION.MD.ANGVEL_ZERO),
and related keywords. For periodic bulk systems, small total-momentum drift is usually less
important.

## Choosing an ensemble

The most common choices are `NVE`, `NVT`, `NPT_I`, `NPT_F`, and `LANGEVIN`.

- `NVE` is microcanonical dynamics. Use it to test the time step, force convergence, grid settings,
  and extrapolation through energy conservation. It is also suitable for production when isolated
  Hamiltonian dynamics is desired.
- `NVT` is canonical sampling at fixed cell. It is usually the most convenient ensemble for
  equilibration and for production when the volume is known or fixed.
- `NPT_I` uses isotropic cell fluctuations. It is appropriate when only the cell volume should
  fluctuate, for example for liquids, gases, or approximately isotropic solids.
- `NPT_F` uses a flexible cell. Use it when the cell shape and anisotropic stress relaxation are
  part of the physical question, for example for solids or interfaces.
- `LANGEVIN` is stochastic constant-temperature dynamics. It is useful for robust equilibration,
  dissipative dynamics, or models where a stochastic heat bath is intentional.

Other supported values of [ENSEMBLE](#CP2K_INPUT.MOTION.MD.ENSEMBLE) are more specialised: `NPE_I`,
`NPE_F`, `ISOKIN`, `REFTRAJ`, `MSST`, `MSST_DAMPED`, `HYDROSTATICSHOCK`, `NVT_ADIABATIC`, and
`NPT_IA`. See the input reference for details.

Note that NPT simulations are sensitive to system size and barostat settings; small cells can show
large pressure fluctuations and artificial cell oscillations. For small cells, consider fixed-cell
`NVT` after a separate cell optimisation, or use a larger supercell before production `NPT`.

## Thermostats

Constant-temperature ensembles use [MOTION/MD/THERMOSTAT](#CP2K_INPUT.MOTION.MD.THERMOSTAT). The
input keyword default for [TYPE](#CP2K_INPUT.MOTION.MD.THERMOSTAT.TYPE) is currently `NOSE`, but for
routine `NVT` equilibration and production `CSVR` is often the safer practical default.

- `CSVR`: recommended general-purpose choice. It gives canonical sampling through velocity rescaling
  and is robust for both equilibration and production.
- `NOSE`: Nose-Hoover chain thermostat. It can be appropriate for well-equilibrated systems, but it
  is less robust when the system is far from equilibrium or very small, and should not be the
  default recommendation for new users.
- `GLE`: generalised Langevin thermostat. Use it for specialised coloured-noise applications or
  workflows that need a fitted GLE kernel.
- `AD_LANGEVIN`: adaptive Langevin thermostat. Use it when adaptive stochastic thermostatting is
  specifically desired.

A practical `NVT` setup is:

```none
&MOTION
  &MD
    ENSEMBLE NVT
    TEMPERATURE 300
    TIMESTEP 0.5
    &THERMOSTAT
      TYPE CSVR
      &CSVR
        TIMECON 100.0
      &END CSVR
    &END THERMOSTAT
  &END MD
&END MOTION
```

The [CSVR/TIMECON](#CP2K_INPUT.MOTION.MD.THERMOSTAT.CSVR.TIMECON) (unit: fs) value controls the
coupling strength. Useful starting values are `50` for rapid heating or cooling, `100` as a simple
robust default, and `200` or larger for gentler production sampling after equilibration.

Do not use crude velocity rescaling as a production thermostat.
[TEMP_TOL](#CP2K_INPUT.MOTION.MD.TEMP_TOL) can trigger velocity rescaling when the instantaneous
temperature deviates from the target, but this keyword is obsolescent; a CSVR thermostat with a
short time constant is usually the better early- equilibration replacement.

## Barostats and NPT simulations

NPT simulations additionally use [MOTION/MD/BAROSTAT](#CP2K_INPUT.MOTION.MD.BAROSTAT). The pressure
target is set by [PRESSURE](#CP2K_INPUT.MOTION.MD.BAROSTAT.PRESSURE), and the barostat time scale by
[TIMECON](#CP2K_INPUT.MOTION.MD.BAROSTAT.TIMECON). For `NPT_F`,
[VIRIAL](#CP2K_INPUT.MOTION.MD.BAROSTAT.VIRIAL) can restrict which Cartesian components of the
virial are used for cell relaxation.

Use `NPT_I` when isotropic volume fluctuations are sufficient. Use `NPT_F` only when the cell shape
should respond to anisotropic stress. Full-cell NPT is powerful but can amplify problems from small
cells, poor initial structures, noisy stresses, or too aggressive barostat coupling. A robust
workflow is to equilibrate first in `NVT`, then run a short NPT test, inspect the cell trajectory
and pressure oscillations, and only then start production sampling. If the cell volume or shape
oscillates strongly, use a larger barostat time constant, a larger supercell, or fixed-cell `NVT` if
pressure sampling is not essential.

## Time step and validation

The time step must resolve the fastest relevant nuclear motion. For AIMD with light atoms, 0.5 ~ 1.0
fs is a common starting range. Use a smaller value for high temperature, weak bonds, reactive
events, poor starting structures, or hydrogen-rich systems. Constraints, deuteration, classical
force fields, or multiple-time-step schemes may allow larger time steps, but always validate the
choice.

A sensitive validation test is a short `NVE` trajectory. The total-energy oscillation and long-time
drift should be small relative to the kinetic energy or target temperature. Energy conservation
alone is not always sufficient for heterogeneous systems; also inspect the region or property of
interest when only part of the system matters.

## Electronic-structure settings for AIMD

Each AIMD step requires an electronic-structure calculation. Important Quickstep settings include
[EPS_SCF](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.EPS_SCF),
[EPS_DEFAULT](#CP2K_INPUT.FORCE_EVAL.DFT.QS.EPS_DEFAULT), and the plane-wave grid parameters in
[MGRID](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID), especially
[CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.CUTOFF) and
[REL_CUTOFF](#CP2K_INPUT.FORCE_EVAL.DFT.MGRID.REL_CUTOFF).

AIMD settings do not always need to match high-precision static energies or final property
calculations. For long trajectories, stable forces and negligible drift are often more important
than over-converging every SCF cycle. For example, if one element in your system use `CUTOFF 400` in
optimization and/or energy tasks, then in AIMD you can set to `CUTOFF 300`. However, too loose
settings can create noisy forces, poor energy conservation, artificial heating, or wrong dynamics.

## Extrapolation of the electronic initial guess

In a normal MD trajectory, consecutive geometries are close to each other. CP2K can use previous
electronic solutions to build a better initial guess for the next step through
[EXTRAPOLATION](#CP2K_INPUT.FORCE_EVAL.DFT.QS.EXTRAPOLATION) and
[EXTRAPOLATION_ORDER](#CP2K_INPUT.FORCE_EVAL.DFT.QS.EXTRAPOLATION_ORDER). A good extrapolation keeps
the electronic state continuous and usually reduces the number of SCF iterations at each MD step,
which is often crucial for efficiency.

For standard Born-Oppenheimer MD, `ASPC` is the default and generally the best first choice. It is
similar in spirit to `PS`, but is designed for MD stability rather than only initial-guess accuracy.
`EXTRAPOLATION_ORDER 3` is a good starting point for `ASPC` and `PS`.

`GEXT_PROJ` and `GEXT_PROJ_QTR` are also worth considering. They can use higher extrapolation
orders, typically 4 ~ 10. `GEXT_PROJ_QTR` is better suited for long-term dynamical stability in
principle.

## Equilibration

Equilibration should remove artefacts of the initial structure and velocity distribution before any
production statistics are collected. It does not have to be performed in `NVT`. A good equilibration
setup is usually the same as, or close to, the intended production ensemble: use `NVT` for
fixed-cell canonical sampling, `NPT_I` or `NPT_F` when pressure and cell relaxation are part of the
target workflow, and `LANGEVIN` when stochastic damping is desired. A typical protocol is:

- Prepare a reasonable starting structure;
- Set the target temperature directly for ordinary heating, or use staged heating/cooling when an
  abrupt temperature change is likely to destabilise the trajectory or when the thermal path is
  important;
- Equilibrate in the chosen target-like ensemble with suitable thermostat and, when needed, barostat
  time constants;
- Continue with adjusted thermostat or barostat parameters for gentler production sampling;
- Discard the equilibration part before computing averages.

For heating with `CSVR`, a practical two-stage strategy is to use a relatively short thermostat time
constant first, for example `TIMECON 50`, to bring the system close to the target temperature
quickly, and then switch to a gentler production setting such as `TIMECON 200`. This is useful when
the initial and target temperatures are not too far apart and the structure is already reasonable,
but a short equilibration stage is still desired before collecting statistics.

If a separate equilibration stage is not necessary, or if a simple robust setup is preferred, one
can also run the whole trajectory with an intermediate coupling strength, for example `TIMECON 100`.
This is often sufficient for ordinary room-temperature or moderately high-temperature simulations,
and the equilibration and production parts may then differ only in which initial frames are
discarded from the analysis.

If the target temperature is very high, directly heating to that temperature can make the simulation
likely to crash. In such cases, staged heating is recommended: increase
[TEMPERATURE](#CP2K_INPUT.MOTION.MD.TEMPERATURE) through several short equilibration stages until
the final target temperature is reached. This is especially useful for AIMD at several thousand
kelvin, reactive systems, poor initial structures, or systems with close contacts. The thermostat
`TIMECON` still controls the coupling strength within each stage, but it does not replace the need
for a gradual temperature schedule when an abrupt jump would destabilise the trajectory.

For cooling or annealing, it is usually best to lower the target temperature in stages, especially
when the final structure depends on the cooling path. One-step cooling can also be used, but the
thermostat coupling should be weak enough, with a sufficiently large `TIMECON` so that the
temperature decreases slowly rather than being quenched abruptly. The appropriate schedule depends
on the system and on whether the goal is sampling, relaxation, or deliberate annealing.

For stochastic equilibration, use the `LANGEVIN` ensemble and adjust
[GAMMA](#CP2K_INPUT.MOTION.MD.LANGEVIN.GAMMA): smaller values are closer to NVE-like dynamics and
equilibrate more slowly; larger values thermostat more aggressively.

Temperature is an instantaneous kinetic estimator and is expected to fluctuate. A single snapshot at
the target temperature is not evidence of equilibration. Check structural, energetic, and
property-specific observables. For heterogeneous systems, do not rely only on total energy or global
temperature; also inspect the region of interest.

## Trajectories, energies, and restarts

The trajectory is controlled by [MOTION/PRINT/TRAJECTORY](#CP2K_INPUT.MOTION.PRINT.TRAJECTORY). By
default, positions are written to files such as `project-pos-1.xyz`, where `project` is the value of
[PROJECT_NAME](#CP2K_INPUT.GLOBAL.PROJECT_NAME). Velocities, forces, cell vectors, stresses, and
restart files are controlled by other print keys in [MOTION/PRINT](#CP2K_INPUT.MOTION.PRINT).

The file `project-1.ener` is often the first file to inspect. It contains columns such as

```none
md_step, time[fs], e_kin[hartree], temp[K], e_pot[hartree], e_tot[hartree], elapsed_time_per_step[s]
```

For thermostatted or barostatted simulations, also inspect the corresponding conserved quantity or
cell output when available. For production runs, write restart files frequently enough that the
trajectory can be continued without losing significant sampling time.

## Validation checklist

Before running a long production trajectory, it is suggested to check that:

- The planned length of trajectory and reserved computational time is long enough to overcome
  (auto)correlations and/or observe phenomena;
- The force method and key parameters is validated on the particular system and is suitable for
  describing the target property;
- The structure is physically reasonable and the shortest interatomic distances are sensible;
- The periodic cell is large enough for the target property;
- The SCF procedure converges reliably on representative snapshots;
- The extrapolation method reduces SCF iterations without causing instability;
- The time step gives acceptable energy conservation in a short test;
- Thermostat and barostat time constants do not dominate the physical dynamics of interest;
- NPT cell fluctuations are reasonable and not dominated by small-cell oscillations;
- Output frequencies are sufficient for the target analysis but not so high that I/O dominates.
