# Metadynamics

## Introduction

**Metadynamics** is a method that allows the acceleration of rare events and estimation of the free
energy of a system undergoing conformational transitions.

In the following exercise, we will explore the dynamic equilibrium between formic acid and water
molecules at the rutile (110) $TiO_2$ surface. For more information about the keywords used in the
input files, refer to [METADYN](#CP2K_INPUT.MOTION.FREE_ENERGY.METADYN).

The tasks to be performed are:

- Set up and run preliminary simulations to understand the dynamics of formic acid and water on
  $TiO_2$ obtained from DFT-based Born-Oppenheimer MD simulations;
- Metadynamics simulation to trigger the replacement of an adsorbed formate/formic acid by a water
  molecule by following changes of two different collective variables.

All relevant files can be downloaded
[here](https://github.com/cp2k/cp2k-examples/tree/master/metadynamics).

```{youtube} t8jxq4e0Qtw
---
url_parameters: ?start=1839
align: center
privacy_mode:
---
```

## First task: dynamics of formic acid and water molecules on rutile (110)

In the proposed example, formic acid molecules are initially adsorbed onto the rutile surface as
formate and $H^+$, with some water molecules co-adsorbed and some others in the vicinity of the
system. $HCOO^-$ is bound either in a bridged, bidentate way to surface titanium atoms or in a
monodentate manner, with only one of its oxygen atoms bonded to a lattice Ti. The corresponding
proton is attached to a lattice oxygen atom forming a hydroxyl.

The model we are using is fully periodic, and enough space must be left above the free water
molecules to avoid interactions with the images along the z-direction. Here, we added a ~20 Å
vacuum.

```
    &CELL
      ABC 19.659 17.806 33.110 
      ALPHA_BETA_GAMMA 90 90 90
      PERIODIC XYZ
    &END CELL
```

The first step is to perform simple MD simulations for a few picoseconds at a constant temperature
to monitor possible rearrangements of the adsorbates. In this case, the initial equilibration of the
whole system is obtained by a short MD run of 10 ps, at 300 K. The coordinates of the surface after
MD are given in
[snapshot_MD-300K.xyz](https://github.com/cp2k/cp2k-examples/blob/master/metadynamics/snapshot_MD-300K.xyz).
By visualizing the trajectory produced by the MD run, we notice there are no considerable
rearrangements of the adsorbates apart from simple fluctuations of the H-bond environment.

The next step is to set up a few collective variables (CV) that can be later used for the
metadynamics (MTD) simulations. These CVs have to be carefully chosen and must describe the relevant
configuration changes one aims to describe. It might also be useful to study the typical behavior of
the selected CVs along an unbiased MD run by performing preliminary runs to monitor the dynamics of
these variables. To do so, an MTD simulation needs to be set up without the addition of any type of
bias.

```
  &FREE_ENERGY
    &METADYN
      DO_HILLS .FALSE.
…
    &END METADYN
&END FREE_ENERGY
```

The evolution of the chosen variables is then monitored while the system explores the configurations
around the initial structure, i.e., belonging to the initial minimum on the free energy surface
(FES). Understanding the typical fluctuation amplitudes of the CVs is important for two reasons: i)
one learns which variations of the CV can occur spontaneously and would not need to be biased, and
in which cases the CV cannot move without activation, and ii) one would be able to appropriately set
up the width of the Gaussian hills that will build up the potential along the biased MTD simulation.

To save computational resources, the previously obtained trajectory can be used by simply setting up
a [REFTRAJ](#CP2K_INPUT.MOTION.MD.REFTRAJ) simulation together with the definition of the relevant
CVs. In this way, the trajectory will be read, and the instantaneous value of the selected CVs will
be computed at every timestep.

```
  &MD
    ENSEMBLE REFTRAJ
    STEPS 20596
    &REFTRAJ
      TRAJ_FILE_NAME trajectory_MD-300K.xyz
      STRIDE 1
      EVAL_ENERGY_FORCES .FALSE.
    &END REFTRAJ
  &END MD
```

The input
[TiO2_reftraj.inp](https://github.com/cp2k/cp2k-examples/blob/master/metadynamics/TiO2_reftraj.inp)
was prepared with the definition of two CVs. The first one is the coordination number (CN) between
one of the oxygen atoms of one specific monodentate formate and the lattice titanium atom to which
it is bound and describes the interaction between these two atoms. This CV should be approximately 0
if this formate leaves and moves far from the surface, and 1 when the molecule is close to the
surface.

```
    &COLVAR
      &COORDINATION
          ATOMS_FROM 453 ! oxygen from formate 
          ATOMS_TO 286 ! index of Ti kind
          R0 [angstrom] 2.9
          ND 12 
          NN 8
      &END COORDINATION
    &END COLVAR
```

NN and ND determine the curvature of the function used to compute the CN, and $R_0$ is the reference
O-Ti distance,
$CN_{O-Ti} = \frac{1}{N_O} \sum_{i O} \sum_{j Ti} \frac{1-(\frac{r_{ij}}{R_0})^{NN}}{1-(\frac{r_{ij}}{R_0})^{ND}}$
.

The second CV describes the interaction between the oxygen atom of one specific water molecule and
the same titanium atom used in the description of the other CV. We also use the CN between these
species to define this collective variable, and its value should also lie between zero and one.

```
    &COLVAR
      &COORDINATION
        ATOMS_FROM 501 ! oxygen from water 
        ATOMS_TO 286 ! index of Ti kind 
        R0 [angstrom] 2.9
        ND 12 
        NN 8
      &END COORDINATION
    &END COLVAR
```

Although no bias is added at this time, for each defined COLVAR an MTD variable is initialized. The
[PRINT/COLVAR](#CP2K_INPUT.MOTION.FREE_ENERGY.METADYN.PRINT.COLVAR) section controls the printing of
the COLVAR output file, which contains, among other information, the instantaneous values of the CVs
(second and third columns) at the indicated time, in fs (first column).

By plotting the instantaneous values of the CVs along the 10 ps MD run, the amplitude of the
equilibrium fluctuations can be evaluated and then used to set up the size of the Gaussian hills
that build up the biasing potential during the MTD simulation. With the NN and ND values of 8 and
12, respectively, the first CV fluctuates close to 1, with an amplitude smaller than 0.2, whereas
the second one fluctuates close to 0, with an amplitude of approximately 0.05. This indicates that
this specific formate has some freedom to move away from the surface, but without fully desorbing.
As a consequence, the free water molecule is not able to get close to the titanium site, which
explains the low fluctuation amplitude in the second CV.

## Second task: metadynamics of the dynamic equilibrium between formic acid and water on rutile (110)

The MTD simulation employs the above described CVs, and the input file can be found in
[TiO2_metadyn.inp](https://github.com/cp2k/cp2k-examples/blob/master/metadynamics/TiO2_metadyn.inp).
The [MOTION/FREE_ENERGY/METADYN](#CP2K_INPUT.MOTION.FREE_ENERGY.METADYN) input section was modified
to activate the MTD algorithm.

```
  &FREE_ENERGY
    &METADYN
      DO_HILLS .TRUE. 
      NT_HILLS 60
      WW 0.5E-03 

      &METAVAR 
        COLVAR 1
        SCALE 0.05 
      &END METAVAR

      &METAVAR
        COLVAR 2
        SCALE 0.05
      &END METAVAR

      &PRINT
        &COLVAR
           COMMON_ITERATION_LEVELS 3
           &EACH
             MD 1
           &END EACH
        &END COLVAR

        &HILLS
           COMMON_ITERATION_LEVELS 3
           &EACH
             MD 1
           &END EACH
        &END HILLS
      &END PRINT
    &END METADYN
  &END FREE_ENERGY
```

One Gaussian hill is deposited every NT_HILLS timesteps, with the height of the hill given by WW, in
Hartree. Together with the width of the Gaussian hills, given by SCALE, these parameters determine
the accuracy of the description of the FES through the MTD biasing potential. Since each variable
has, in principle, different dimensions and dynamics, the shape of the hills filling up the $N_{CV}$
-dimensional configurations space, as defined by the chosen CVs, is not isotropic. The parameter
SCALE, associated with the i-th MTD variable, determines the amplitude of the Gaussian in the i-th
space-direction of the $N_{CV}$ -dimensional configuration space. All three parameters (WW,
NT_HILLS, and SCALE) can be modified along the same MTD run by simply restarting the simulations
employing different values in the input file. This is a very useful feature in case the dynamics of
one or more variables changes after some event has occurred.

For an efficient and accurate exploration of the configurations space, it is important that the
added hills are not too large, otherwise, important features of the topography of the FES might not
be sufficiently well-resolved, or even the MTD-trajectory could follow the wrong, not minimum energy
path. On the other hand, using too small hills might require excessively long simulation times.

The printing of the HILLS file is controlled by
[PRINT/HILLS](#CP2K_INPUT.MOTION.FREE_ENERGY.METADYN.PRINT.HILLS), with the information given in the
following order: timestep, coordinates of the center in the CV space for CV1, coordinates of the
center in the CV space for CV2, width of CV1, width of CV2, and height of the hill.

The resulting trajectory is about 34 ps long and shows the desorption of one formate as formic acid,
with subsequent adsorption of a water molecule to the freed titanium site. Additionally, in this
case, we added a [reflective wall](#CP2K_INPUT.MOTION.FREE_ENERGY.METADYN.METAVAR.WALL) after around
17 ps of simulation, which prevents the desorbed formic acid from moving too far from the surface.
Other quantities can be monitored from the
[PRINT/COLVAR](#CP2K_INPUT.MOTION.FREE_ENERGY.METADYN.PRINT.COLVAR) output file, in the following
order: timestep, the instantaneous value of CV1, the instantaneous value of CV2, the instantaneous
gradient of the bias potential computed wrt CV1, the instantaneous gradient of the bias potential
computed wrt CV2, instantaneous gradient wrt CV1 of wall potentials (if present), instantaneous
gradient wrt CV2 of wall potentials (if present), the instantaneous value of the bias potential,
instantaneous values of the wall potentials (if present).

Finally, by using this information, the FES can be reconstructed with the help of the graph script
available with cp2k, in the *tools* directory. The command line, in this case, would look as
follows:

```
graph.psmp -i HILLS -stride 10 -ndim 2 -ndw 1 2 -cp2k -integrated_fes 
```

Where HILLS is the last restart file printed during the metadynamics simulation.
