# Langevin Dynamics

In this tutorial, we are going to show the reader how to perform Langevin molecular dynamics for a
sub set of atoms in the simulation cell, with the rest of the atoms undergoing Born-Oppenheimer
molecular dynamics. We assume the reader has already got the basic knowhow of performing molecular
dynamics using CP2K.

To be able to perform this calculation, you must have CP2K version 2.5 or above.

We will use a simple 64 atoms face centred cubic bulk Si as an example. The system will start from a
relaxed ground state structure (i.e. a geometry optimisation calculation has first been performed
already), with the first atom will be displaced slightly from its optimal position, which then
kick-starts the molecular dynamics. The initial atomic velocities are set to be zero. The first 24
atoms in the system will be performing NVT Langevin dynamics, while the rest will be performing NVE
Born-Oppenheimer dynamics.

The example files can be downloaded from
[here](https://github.com/cp2k/cp2k-examples/tree/master/langevin_regions). The calculation were
carried out using CP2K version 2.5.

## Input Flags

All we need to do to perform this calculation is to add/modify a few flags in the
[MD](#CP2K_INPUT.MOTION.MD) subsection in the main CP2K input file.

The relevant sections for the example are listed below:

In the [MD](#CP2K_INPUT.MOTION.MD) subsection, we need to set the
[ENSEMBLE](#CP2K_INPUT.MOTION.MD.ENSEMBLE) keyword to `LANGEVIN`

```none
ENSEMBLE LANGEVIN
```

This ensures we will be doing Langevin molecular dynamics. The method of applying mixed NVE and NVT
dynamics **only works** if Langevin MD is switched on.

The important thing to do next is to define the thermal regions, which controls whether each region
will be performing NVT Langevin MD or NVE Born-Oppenheimer MD. In our example, we have inside the
[MD](#CP2K_INPUT.MOTION.MD) subsection the following:

```none
&THERMAL_REGION
  DO_LANGEVIN_DEFAULT F
  &DEFINE_REGION
    TEMPERATURE $temp
    DO_LANGEVIN T
    LIST  1..24
  &END DEFINE_REGION
  &DEFINE_REGION
    # TEMPERATURE $temp
    DO_LANGEVIN F
    LIST  25..64
  &END DEFINE_REGION
  &PRINT
    &LANGEVIN_REGIONS
    &END LANGEVIN_REGIONS
  &END PRINT
&END THERMAL_REGION
```

The [DO_LANGEVIN_DEFAULT](#CP2K_INPUT.MOTION.MD.THERMAL_REGION.DO_LANGEVIN_DEFAULT) keyword defines
if Langevin MD is to be performed for all atoms that are **outside** the regions defined in
[THEMAL_REGION](#CP2K_INPUT.MOTION.MD.THERMAL_REGION); the default value is `F`, which means any
atoms not included in the thermal regions will undergo NVE MD by default. If the value is set to
`T`, then any atoms not included in the thermal regions will undergo Langevin MD by default.

Each of the subsections

```none
&DEFINE_REGION
  TEMPERATURE $temp
  DO_LANGEVIN T
  LIST  1..24
&END DEFINE_REGION
```

defines a thermal region. The subsections may be repeated an arbitrary $N$ number of times for $N$
thermal regions. Inside, the
[DO_LANGEVIN](#CP2K_INPUT.MOTION.MD.THERMAL_REGION.DEFINE_REGION.DO_LANGEVIN) keyword defines if the
atoms defined in the region is to undergo NVT Langevin MD (`F`), or NVE MD (`F`); the
[LIST](#CP2K_INPUT.MOTION.MD.THERMAL_REGION.DEFINE_REGION.LIST) keyword defines the list of atoms in
the particular thermal region;
[TEMPERATURE](#CP2K_INPUT.MOTION.MD.THERMAL_REGION.DEFINE_REGION.TEMPERATURE) defines the target
temperature for the region, which is only taken into account if
[DO_LANGEVIN](#CP2K_INPUT.MOTION.MD.THERMAL_REGION.DEFINE_REGION.DO_LANGEVIN) is set to `T`. Note
that in this example, temperature is set by referring to an input preprocessor variable `$temp`,
whose value (500 K) is defined at the top of the main input file.

By default, [DO_LANGEVIN](#CP2K_INPUT.MOTION.MD.THERMAL_REGION.DEFINE_REGION.DO_LANGEVIN) is set to
`T`, however, this will only be taken into account if the [ENSEMBLE](#CP2K_INPUT.MOTION.MD.ENSEMBLE)
keyword in the `MD` subsection is set to `LANGEVIN`, therefore it will not effect the definition of
the thermal region for molecular dynamics using ensembles other than `LANGEVIN`.

In our example, we have defined two regions. The first region contains atoms 1 to 24, undergoing NVT
Langevin MD with target temperature of 500 K, and the second region contains atoms 25 to 64,
undergoing NVE Born-Oppenheimer MD. Note that since
[DO_LANGEVIN_DEFAULT](#CP2K_INPUT.MOTION.MD.THERMAL_REGION.DO_LANGEVIN_DEFAULT) is set to `F` (by
default), in principle, we do not have to define the second region. It is defined here in this
example to show how these regions can be defined.

The [TEMPERATURE](#CP2K_INPUT.MOTION.MD.TEMPERATURE) keyword in [MD](#CP2K_INPUT.MOTION.MD)
subsection defines the target temperature of the molecular dynamics for all atoms **left out of**
the defined thermal regions. If the [THERMAL_REGION](#CP2K_INPUT.MOTION.MD.THERMAL_REGION)
subsection is not explicitly present in the input file, then CP2K assumes all atoms in the
simulation cell undergoes Langevin MD. In this case this
[TEMPERATURE](#CP2K_INPUT.MOTION.MD.TEMPERATURE) keyword defines the target temperature for the
entire system. If the [THERMAL_REGION](#CP2K_INPUT.MOTION.MD.THERMAL_REGION) subsection is
explicitly present in the input file, then this [TEMPERATURE](#CP2K_INPUT.MOTION.MD.TEMPERATURE)
keyword sets the **default** target temperature for the atoms if they undergo Langevin MD. This
value is overridden by the
[TEMPERATURE](#CP2K_INPUT.MOTION.MD.THERMAL_REGION.DEFINE_REGION.TEMPERATURE) keywords in each
[DEFINE_REGION](#CP2K_INPUT.MOTION.MD.THERMAL_REGION.DEFINE_REGION) subsections.

Information on the NVT and NVE regions may be printed out by using the
[PRINT](#CP2K_INPUT.MOTION.MD.THERMAL_REGION.PRINT) subsection in the
[THERMAL_REGION](#CP2K_INPUT.MOTION.MD.THERMAL_REGION) subsection:

```none
&PRINT
  &LANGEVIN_REGIONS
  &END LANGEVIN_REGIONS
&END PRINT
```

Simply add the subsection
[LANGEVIN_REGIONS](#CP2K_INPUT.MOTION.MD.THERMAL_REGION.PRINT.LANGEVIN_REGIONS) to
[PRINT](#CP2K_INPUT.MOTION.MD.THERMAL_REGION.PRINT). The region information will be written in an
output file with suffix: `lgv_regions`.

## Importance of Initial Velocity To Consistency Of Calculations

If the initial velocities of the atoms are not explicitly defined in the input, CP2K will randomise
the atomic velocities to give an initial temperature corresponding to the target temperature defined
by [TEMPERATURE](#CP2K_INPUT.MOTION.MD.TEMPERATURE) keyword in [MD](#CP2K_INPUT.MOTION.MD)
subsection.

Due to the stochastic nature of Langevin MD, the trajectories of the atoms generated as a result of
the pseudo-random number generators will be dependent on the initial velocities of the atoms.
Therefore, if you are to perform to two calculations with the same physical thermal regions setup,
but do not specify the initial velocities, there is a chance that the velocities and energies at
each MD step can be different for the two calculations. This can arise, in the cases where the
setups are physically the same, but computationally different in terms of input setup: for example:
setting [DO_LANGEVIN_DEFAULT](#CP2K_INPUT.MOTION.MD.THERMAL_REGION.DO_LANGEVIN_DEFAULT) to `T` and
define a NVE region with atoms 25 to 64
([DO_LANGEVIN](#CP2K_INPUT.MOTION.MD.THERMAL_REGION.DEFINE_REGION.DO_LANGEVIN) set to `F`), is
physically the same as setting
[DO_LANGEVIN_DEFAULT](#CP2K_INPUT.MOTION.MD.THERMAL_REGION.DO_LANGEVIN_DEFAULT) to `F`, and define a
NVT region with atoms 1 to 24
([DO_LANGEVIN](#CP2K_INPUT.MOTION.MD.THERMAL_REGION.DEFINE_REGION.DO_LANGEVIN) set to `T`). However,
this difference in the input may cause a different randomised set of initial velocities at the start
of the calculations, and make the two calculations not matching exactly step-by-step. Of course, the
overall physical results will still be the same.
