# Geometry Optimisation

## Introduction

This tutorial is designed to illustrate how to relax the structure of a system (without changing the
cell dimensions) using CP2K. We use the relaxation of a water ($H_2O$) molecule as an example.

The example files can be downloaded from
[here](https://github.com/cp2k/cp2k-examples/tree/master/geometry_optimisation). The calculation was
carried out with CP2K version 2.4.

It should be noted that before running the geometry optimisation, the reader should have already
know how to perform a simple Kohn-Sham Density Functional Theory energy and force calculation (this
is covered in tutorial
[Calculating Energy and Forces using QUICKSTEP](https://www.cp2k.org/howto:static_calculation)), and
they should also know how how to find a sufficient grid cutoff for the static energy calculations
(this is covered in tutorial
[Converging the CUTOFF and REL_CUTOFF](https://www.cp2k.org/howto:converging_cutoff)).

DIIS (direct inversion in the iterative subspace or direct inversion of the iterative subspace),
also known as Pulay mixing, is an extrapolation technique. DIIS was developed by Peter Pulay in the
field of computational quantum chemistry with the intent to accelerate and stabilize the convergence
of the Hartree–Fock self-consistent field method.

```{note}
The [ALPHA](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.MIXING.ALPHA) and
[NBUFFER](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.MIXING.NBUFFER) parameters might have to be reduced to
avoid a [bad condition number](https://github.com/cp2k/cp2k/issues/2360).
```

## Input Files

The input file for a geometry calculation is shown below:

```none
&GLOBAL
  PROJECT H2O
  RUN_TYPE GEO_OPT
  PRINT_LEVEL LOW
&END GLOBAL
&FORCE_EVAL
  METHOD QS
  &SUBSYS
    &CELL
      ABC 12.4138 12.4138 12.4138
    &END CELL
    &COORD
      O      12.235322       1.376642      10.869880
      H      12.415139       2.233125      11.257611
      H      11.922476       1.573799       9.986994
    &END COORD
    &KIND H
      BASIS_SET DZVP-GTH-PADE
      POTENTIAL GTH-PADE-q1
    &END KIND
    &KIND O
      BASIS_SET DZVP-GTH-PADE
      POTENTIAL GTH-PADE-q6
    &END KIND
  &END SUBSYS
  &DFT
    BASIS_SET_FILE_NAME ./BASIS_SET
    POTENTIAL_FILE_NAME ./POTENTIAL
    &QS
      EPS_DEFAULT 1.0E-7
    &END QS
    &MGRID
      CUTOFF 200
      NGRIDS 4
      REL_CUTOFF 30
    &END MGRID
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-05
      MAX_SCF 200
      &DIAGONALIZATION T
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING T
        ALPHA 0.5
        METHOD PULAY_MIXING
        NPULAY 5
      &END MIXING
      &PRINT
        &RESTART OFF
        &END RESTART
      &END PRINT
    &END SCF
    &XC
      &XC_FUNCTIONAL PADE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
&END FORCE_EVAL
&MOTION
  &GEO_OPT
    TYPE MINIMIZATION
    MAX_DR    1.0E-03
    MAX_FORCE 1.0E-03
    RMS_DR    1.0E-03
    RMS_FORCE 1.0E-03
    MAX_ITER 200
    OPTIMIZER CG
    &CG
      MAX_STEEP_STEPS  0
      RESTART_LIMIT 9.0E-01
    &END CG
  &END GEO_OPT
  &CONSTRAINT
    &FIXED_ATOMS
      COMPONENTS_TO_FIX XYZ
      LIST 1
    &END FIXED_ATOMS
  &END CONSTRAINT
&END MOTION
```

The reader should already be familiar with the [GLOBAL](#CP2K_INPUT.GLOBAL) and
[FORCE_EVAL](#CP2K_INPUT.FORCE_EVAL) sections. For geometry optimisation calculations, we must set
[RUN_TYPE](#CP2K_INPUT.GLOBAL.RUN_TYPE) in [GLOBAL](#CP2K_INPUT.GLOBAL) section to `GEO_OPT`:

```none
RUN_TYPE GEO_OPT
```

In this example, we note that we have chosen diagonalisation of the Kohn-Sham Hamiltonian for the
evaluation of wavefunctions, and used Pulay mixing for the self-consistency loops. 5 histories are
used for Pulay mixing.

The important section for geometry optimisation settings are contained in subsection
[GEO_OPT](#CP2K_INPUT.MOTION.GEO_OPT) of [MOTION](#CP2K_INPUT.MOTION) section. Note that
[GEO_OPT](#CP2K_INPUT.MOTION.GEO_OPT) subsection only applies to the calculation where the cell
dimensions do not change. Calculations which allows the relaxation of the cell are covered in a
separate tutorial.

```none
&GEO_OPT
  TYPE MINIMIZATION
  MAX_DR    1.0E-03
  MAX_FORCE 1.0E-03
  RMS_DR    1.0E-03
  RMS_FORCE 1.0E-03
  MAX_ITER 200
  OPTIMIZER CG
  &CG
    MAX_STEEP_STEPS  0
    RESTART_LIMIT 9.0E-01
  &END CG
&END GEO_OPT
```

The [TYPE](#CP2K_INPUT.MOTION.GEO_OPT.TYPE) keyword sets whether the geometry optimisation is for
finding the local minima (`MINIMIZATION`) or for finding the saddle point transition state
(`TRANSITION_STATE`). The keywords [MAX_DR](#CP2K_INPUT.MOTION.GEO_OPT.MAX_DR),
[MAX_FORCE](#CP2K_INPUT.MOTION.GEO_OPT.MAX_FORCE), [RMS_DR](#CP2K_INPUT.MOTION.GEO_OPT.RMS_DR) and
[RMS_FORCE](#CP2K_INPUT.MOTION.GEO_OPT.RMS_FORCE) set the criteria of whether an optimised geometry
is reached. [MAX_DR](#CP2K_INPUT.MOTION.GEO_OPT.MAX_DR) and
[RMS_DR](#CP2K_INPUT.MOTION.GEO_OPT.RMS_DR) (in Bohr) are the tolerance on the maximum and
root-mean-square of atomic displacements from the previous geometry optimisation iteration;
[MAX_FORCE](#CP2K_INPUT.MOTION.GEO_OPT.MAX_FORCE) and
[RMS_FORCE](#CP2K_INPUT.MOTION.GEO_OPT.RMS_FORCE) (in Bohr/Hartree) are the tolerance on the maximum
and root-mean-square of atomic forces. The geometry is considered to be optimised **only when all
four criteria are satisfied**. The keyword [MAX_ITER](#CP2K_INPUT.MOTION.GEO_OPT.MAX_ITER) sets the
maximum number of geometry optimisation iterations.
[OPTIMIZER](#CP2K_INPUT.MOTION.GEO_OPT.OPTIMIZER) sets the algorithm for finding the stationary
points; in this example we have chosen the conjugate gradients (`CG`) method.

The [CG](#CP2K_INPUT.MOTION.GEO_OPT.CG) subsection sets options for the conjugate gradients
algorithm. In this case, we have configured it so that no steepest descent steps are to be performed
before the start of the conjugate gradients algorithm; and the CG algorithm should be reset (and one
steepest descent step is performed) if the cosine of the angles between two consecutive searching
directions is less than 0.9.

```none
&CONSTRAINT
  &FIXED_ATOMS
    COMPONENTS_TO_FIX XYZ
    LIST 1
  &END FIXED_ATOMS
&END CONSTRAINT
```

We can add constraints to atomic movements by using the [CONSTRAINT](#CP2K_INPUT.MOTION.CONSTRAINT)
subsection in [MOTION](#CP2K_INPUT.MOTION) section. In this example, we choose to fix particular
atoms using the [FIXED_ATOMS](#CP2K_INPUT.MOTION.CONSTRAINT.FIXED_ATOMS) subsection. The keyword
[COMPONENTS_TO_FIX](#CP2K_INPUT.MOTION.CONSTRAINT.FIXED_ATOMS.COMPONENTS_TO_FIX) sets which of the
`X` `Y` `Z` directions are to be fixed, and in this case, the atoms will be completely pinned in all
directions (`XYZ`). The list of atoms to be constrained are given by the
[LIST](#CP2K_INPUT.MOTION.CONSTRAINT.FIXED_ATOMS.LIST) keyword:

```none
LIST 1 2 3 ... N
```

The numbers to the right of [LIST](#CP2K_INPUT.MOTION.CONSTRAINT.FIXED_ATOMS.LIST) are the list of
atomic indices, and correspond to the order (from top to bottom) of the atoms given in the
[COORD](#CP2K_INPUT.FORCE_EVAL.SUBSYS.COORD) subsection of [SUBSYS](#CP2K_INPUT.FORCE_EVAL.SUBSYS)
(of [FORCE_EVAL](#CP2K_INPUT.FORCE_EVAL)). In our example, we have fixed the oxygen atom during
geometry optimisation, so that the water molecule will not move around while its structure is being
relaxed.

## Results

The example is run using the serial version of the CP2K binaries:

```none
cp2k.sopt -o H2O.out H2O.inp &
```

After the job has finished, you should obtain the following files:

- `H2O.out`
- `H2O-pos-1.xyz`
- `H2O-1.restart`
- `H2O-1.restart.bak-1`
- `H2O-1.restart.bak-2`
- `H2O-1.restart.bak-3`

Again, the file `H2O.out` contains the main output of the job. `H2O-pos-1.xyz` contains the trace of
atomic coordinates at each geometry optimisation step in the `xyz` file format. The last set of
atomic coordinates corresponds to the relaxed structure. `H2O-1.restart` is a CP2K input file,
similar to `H2O.inp`, which contains the latest atomic coordinates of the water molecule. Should the
job die for some reason, you can continue the job using the latest atomic coordinates by using
command:

```none
cp2k.sopt -o H2O.out H2O-1.restart &
```

You can of course also use `H2O-1.restart` as a template for writing an input for further
calculations using the relaxed atomic structures.

The files `H2O-1.restart.bak-*` are backup restart files with atomic coordinates obtained from the
previous 1, 2 and 3 geometry optimisation iterations. `H2O-1.restart.bak-1` should be the same as
`H2O-1.restart`.

In the main output file `H2O.out`, at the end of each geometry optimisation step, we will have the
following information:

```none
--------  Informations at step =     1 ------------
 Optimization Method        =                   SD
 Total Energy               =       -17.1643447508
 Real energy change         =        -0.0006776683
 Decrease in energy         =                  YES
 Used time                  =               90.837

 Convergence check :
 Max. step size             =         0.0336570168
 Conv. limit for step size  =         0.0010000000
 Convergence in step size   =                   NO
 RMS step size              =         0.0168136889
 Conv. limit for RMS step   =         0.0010000000
 Convergence in RMS step    =                   NO
 Max. gradient              =         0.0182785685
 Conv. limit for gradients  =         0.0010000000
 Conv. for gradients        =                   NO
 RMS gradient               =         0.0091312361
 Conv. limit for RMS grad.  =         0.0010000000
 Conv. for gradients        =                   NO
---------------------------------------------------
```

The above output segment states that at the end of geometry optimisation step 1, the total energy of
the system is -17.1643447508 (Ha) and none of the criteria for optimised geometry has been reached.
The iteration therefore will carry on, until all criteria becomes `YES`.

At the end of geometry optimisation, one should obtain something like:

```none
--------  Informations at step =    11 ------------
 Optimization Method        =                   SD
 Total Energy               =       -17.1646204766
 Real energy change         =        -0.0000000529
 Decrease in energy         =                  YES
 Used time                  =               49.893

 Convergence check :
 Max. step size             =         0.0003393150
 Conv. limit for step size  =         0.0010000000
 Convergence in step size   =                  YES
 RMS step size              =         0.0001493298
 Conv. limit for RMS step   =         0.0010000000
 Convergence in RMS step    =                  YES
 Max. gradient              =         0.0001787448
 Conv. limit for gradients  =         0.0010000000
 Conv. in gradients         =                  YES
 RMS gradient               =         0.0000786642
 Conv. limit for RMS grad.  =         0.0010000000
 Conv. in RMS gradients     =                  YES
---------------------------------------------------
```

which clearly shows all criteria have been satisfied.

The final Kohn-Sham energies can be obtained at the end of the output:

```none
*******************************************************************************
***                    GEOMETRY OPTIMIZATION COMPLETED                      ***
*******************************************************************************

                   Reevaluating energy at the minimum

Number of electrons:                                                          8
Number of occupied orbitals:                                                  4
Number of molecular orbitals:                                                 4

Number of orbital functions:                                                 23
Number of independent orbital functions:                                     23

 Parameters for the always stable predictor-corrector (ASPC) method:

 ASPC order: 3

 B(1) =   3.000000
 B(2) =  -3.428571
 B(3) =   1.928571
 B(4) =  -0.571429
 B(5) =   0.071429

Extrapolation method: ASPC


SCF WAVEFUNCTION OPTIMIZATION

 Step     Update method      Time    Convergence         Total energy    Change
 ------------------------------------------------------------------------------
    1 Pulay/Diag. 0.50E+00    0.5     0.00005615       -17.1646204762 -1.72E+01
    2 Pulay/Diag. 0.50E+00    1.0     0.00000563       -17.1646347711 -1.43E-05

 *** SCF run converged in     2 steps ***


 Electronic density on regular grids:         -8.0000016293       -0.0000016293
 Core density on regular grids:                7.9999992554       -0.0000007446
 Total charge density on r-space grids:       -0.0000023739
 Total charge density g-space grids:          -0.0000023739

 Overlap energy of the core charge distribution:               0.00000004555422
 Self energy of the core charge distribution:                -43.83289054591484
 Core Hamiltonian energy:                                     12.82175605770555
 Hartree energy:                                              17.97395116120845
 Exchange-correlation energy:                                 -4.12745148966141

 Total energy:                                               -17.16463477110803

ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):              -17.164634771108034
```
