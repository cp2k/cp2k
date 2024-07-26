# Molecular Dynamics

Molecular dynamics is a good method to perform thermodynamical averages, and to look at dynamical
properties. It is also a good starting point for many other more advanced techniques.

During an MD one expects the density to change more or less continuously, thus it is possible to use
the solutions of the previous steps to create the new initial guess for the density (and
wavefunctions for OT). Indeed cp2k can use different
[EXTRAPOLATION](#CP2K_INPUT.FORCE_EVAL.DFT.QS.EXTRAPOLATION) techniques.

For MD a good extrapolation is `PS` with an
[EXTRAPOLATION_ORDER](#CP2K_INPUT.FORCE_EVAL.DFT.QS.EXTRAPOLATION_ORDER) of 3. If you are doing
something where from one step to the other there is little continuity, (geometry optimization,
montecarlo), then something with little continuity like `LINEAR_PS` or just the previous density
(`USE_PREV_P`) might be better. If you change the particle number then pick `USE_GUESS` that
restarts from scratch with the restart method choosen in
[SCF_GUESS](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.SCF_GUESS). The `ASPC` extrapolation is used in
connection with langevin and a reduced SCF convergence (just on SCF step) and history restart to
simulate big systems.

```{youtube} cGDiVSWpXLc
---
url_parameters: ?start=2
align: center
privacy_mode:
---
```

## Trajectory

With the printkey [TRAJECTORY](#CP2K_INPUT.MOTION.PRINT.TRAJECTORY) you can control the output of
the trajectory. The name of the file is by default `projectname-pos-1.xyz` and projectname is
whatever you have written in [PROJECT_NAME](#CP2K_INPUT.GLOBAL.PROJECT_NAME). In the prinkey you can
also change the format from xyz to something else.

An interesting file to check is the `projectname-1.ener` file, in it you can find the following
columns:

```none
md_step, time[fs], e_kin [hartree], temp[K], e_pot [hartree], e_tot [hartree], cpu_time_per_step [s]
```

You should always check it and look at how the system equilibrates.

The `.ener` file and other interesting trajectory files are all controlled with from the
[PRINT](#CP2K_INPUT.MOTION.PRINT) section.

## MD Convergence

If the MD has to be trusted then one should be sure that the trajectory can be trusted. Actually a
simple, and very sensitive test that there are no big technical errors is to perform an NVE
trajectory and look at the energy conservation. The energy conservation has normally two features,
short time oscillations (that are larger when the system is not yet equilibrated) and a long range
drift. If you express these in *K*, then you can compare them with the temperature that you are
simulating at. Another (equivalent) possibility is to express them as fraction of the kinetic
energy. The oscillation and drift (normally per *ps*, but it also depends on how many picoseconds
you want to simulate, and if you want an NVE trajectory) should be small with respect to the kinetic
energy (1% or less is a good value, but obviously it depends on the accuracy that you want to
achieve, more might be acceptable, or less needed).

To improve the energy conservation one can either improve the forces with
[EPS_SCF](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.EPS_SCF) and
[EPS_DEFAULT](#CP2K_INPUT.FORCE_EVAL.DFT.QS.EPS_DEFAULT), improve $\tilde n$ increasing the cutoff,
and reduce the timestep. For the timestep a good value is normally around 1 fs, but if you are
simulating high temperature or special system you might need to reduce it. Normally having good
forces and larger timestep is advantageous. Obviously the timestep cannot be larger than the period
of the highest oscillation of the system (typically H atoms, that is the reason why sometimes D
atoms are used).

Conserving the total energy well is not enough if the system is very heterogeneous and the
interesting part is small. In that case there is the risk that even a large error on it might pass
unnoticed. If this is to be excepted (for example the atom with the sharpest gaussian basis set is
in the interesting part) then checking that part of the system (oscillations, kinetic energy,...) is
advisable.

## Equilibration

For the result to be meaningful the system should be equilibrated, and not in a very unlikely state.
If one is doing shooting or transition path sampling then this is less true, but still the system
should not be in a very high-energy state.

So at the beginning you should build your system in some meaningful way, or from classical
simulations. To have a better starting point you can optimize the energy (if your system is small),
you anneal it, but it is not always needed. Then you can equilibrate it at a given temperature.

To equilibrate a system one can use a crude velocity rescale when it is too far away from the goal
temperature (as established by [TEMP_TOL](#CP2K_INPUT.MOTION.MD.TEMP_TOL)). Doing an NVE simulation
with [TEMP_TOL](#CP2K_INPUT.MOTION.MD.TEMP_TOL) is better way to equilibrate than using the NVT
ensamble that uses the Nose-Hoover chain thermostats and might give spurios effects if you are far
from equilibrium, as it tryes to conserve some extended quantity.the thermostat can store energy
that he will give back at some later point. It should be taken into account that the temperature is
not constant, but does oscillate in the NVE ensamble, these oscillations are connected with the
specific heat and are inversely proportional with sqrt(N) where N is the system size.

A faster way to equilibrate the system is to use `LANGEVIN` as
[ENSEMBLE](#CP2K_INPUT.MOTION.MD.ENSEMBLE), and use a [GAMMA](#CP2K_INPUT.MOTION.MD.LANGEVIN.GAMMA)
of 0.001 \[1/fs\] (or larger if you want to equilibrate faster). Langevin introduces a viscous
dissipative force and a random forces that are in equilibrioum at the given temperature. For
equilibration purposes (and not to simulate removed degrees of freedom, like a solvent), the smaller
gamma the longer it takes to equilibrate, and the closer the trajectory is to an NVE trajectory, the
larger gamma the more the "environment thermostat" influences the system. Also With langevin if your
system is in a bad initial configuration it is a good idea to set a
[TEMP_TOL](#CP2K_INPUT.MOTION.MD.TEMP_TOL), or first anneal a little the system.
