# Monte Carlo

In the most common set up of Gibbs Ensemble Monte Carlo (GEMC) simulation, 2 boxes are utilized to
represent vapor and liquid phases. In order to equilibrate the system, different types of moves are
used:

1. Translations, rotations, and conformational changes
1. Volume exchanges
1. Particle swaps.

The particles are swapped between boxes to equilibrate the chemical potential, volume moves
equilibrate pressure, and the rest of the moves within a box are performed to maintain thermal
equilibrium. The main advantage of the GEMC simulation is that coexisting phases can be simulated
without a physical interface using a unified partition function. Thus, we used GEMC simulations to
determined vapor-liquid coexistence curves for a system.

## Files required to run GEMC

In order to run GEMC certain input files are needed; the two main input files for each box
<path:GEMC_NVT_box1.inp>, <path:GEMC_NVT_box2.inp>, a topology file <path:topology_atoms_WAT.psf>
for the particular component, and a bias file that contain information of the approximate potential
for that component <path:bias_template.inp>. For example, we have provided here the sample files for
64 water molecules.

## Sample of input files

Starting with the input `GEMC_NVT_box1.inp` file first, we note that the location of the basis set
file, as well as the potential file, is declared (in this case these two files are present in the
current working directory). [FORCE_EVAL](#CP2K_INPUT.FORCE_EVAL) initializes the parameters needed
to calculate the energy and forces to describe your system. The `Quickstep` module is used in order
to use of electronic structure methods.

```none
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME GTH_BASIS_SETS
    POTENTIAL_FILE_NAME POTENTIAL
    ...
```

The SCF section in the code below generates atomic density.

```none
&SCF
  SCF_GUESS ATOMIC
&END SCF
```

One can increase the SCF iterations by including the
[MAX_SCF](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.MAX_SCF) keyword. Also,
[EPS_SCF](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.EPS_SCF) declares the expected SCF convergence. Continuing
down the code of the input file, the section shown below details which functional we intend on
using, in our case BLYP. We also use the DFTD2 dispersion correction. Additionally, the
[XC_DERIV](#CP2K_INPUT.ATOM.METHOD.XC.XC_GRID.XC_DERIV) keyword specifies about method used to
compute derivatives.

```none
&XC
  &XC_FUNCTIONAL BLYP
  &END XC_FUNCTIONAL
  &VDW_POTENTIAL
    POTENTIAL_TYPE PAIR_POTENTIAL
    &PAIR_POTENTIAL
      R_CUTOFF 40.0
      TYPE DFTD2
      REFERENCE_FUNCTIONAL BLYP
    &END PAIR_POTENTIAL
  &END VDW_POTENTIAL
  &XC_GRID
    XC_DERIV SPLINE2
    XC_SMOOTH_RHO NONE
  &END XC_GRID
&END XC
```

The cell and cell ref of the box 1 in angstorms.

```none
&CELL
  ABC 13.7151207699 13.7151207699 13.7151207699
  &CELL_REF
    ABC 13.7151207699 13.7151207699 13.7151207699
  &END CELL_REF
&END CELL
```

After the above section of code, there is a listing of coordinates for each atom. After this, the
section

```none
&KIND H
  BASIS_SET TZV2P-GTH
  POTENTIAL GTH-BLYP-q1
&END KIND
```

declares the basis set(TZV2P) intended to be used for the simulation. The following code

```none
&GLOBAL
  PROJECT H2O_MC
  RUN_TYPE MC
  PRINT_LEVEL LOW
&END GLOBAL
```

is intended to described the type of run. Consequently, [RUN_TYPE](#CP2K_INPUT.GLOBAL.RUN_TYPE)
should be set to `MC`, as this is what we are running.

```none
&MOTION
  &MC
    ENSEMBLE GEMC_NVT
    TEMPERATURE 398.0
    IPRINT 1
    LBIAS yes
    LSTOP yes
    NMOVES 8
    NSWAPMOVES 640
    NSTEP 5
    PRESSURE 1.013
    RESTART no
    BOX2_FILE_NAME GEMC_NVT_box1.inp
    RESTART_FILE_NAME mc_restart_2
    ...
  &END MC
&END MOTION
```

The [ENSEMBLE](#CP2K_INPUT.MOTION.MC.ENSEMBLE) `GEMC_NVT` illustrates the particular type of
simulation. It should be noted that the condition [LBIAS](#CP2K_INPUT.MOTION.MC.LBIAS) `YES` must be
true, as we pre sample moves with a classical potential. The [LSTOP](#CP2K_INPUT.MOTION.MC.LSTOP)
keyword determines whether the simulation increment will be in cycles (no), or in steps (yes). The
[NSTEP](#CP2K_INPUT.MOTION.MC.NSTEP) keyword gives the number of MC cycles in a particular
simulation run, and should be adjusted according to the length of the simulation. The line
[RESTART](#CP2K_INPUT.MOTION.MC.RESTART) `NO` should only be set to 'no' for the initial run, then
switched to 'yes' after the first simulation run is complete. The keyword
[BOX2_FILE_NAME](#CP2K_INPUT.MOTION.MC.BOX2_FILE_NAME) gives the file name of the input for Box2 and
uses it as a reference such that the two input files are read together. `GEMC_NVT_box2.inp` has a
similar line that references Box 1, for example: `BOX2_FILE_NAME GEMC_NVT_box1.inp`.

## Sample of output files

Sample output files for 64 $H_2O$ molecules are provided below. The input file for Box 1 and Box2
has already been explained above. We additionally include a sample output file,
<path:GEMC_NVT_box1.out>. Information used to calculate the density at the end of the run is towards
the end of the file and looks like the following:

```none
  |                   BOX 1                      |
  ------------------------------------------------

  ********************************************************************************
  Average Energy [Hartrees]:                                        -1085.04223205
  Average number of molecules:                                         63.00000000
  Average Volume [angstroms**3]:                                       2579.876452
  --------------------------------------------------------------------------------
  Quickstep Moves                           Attempted       Accepted       Percent
                                                    4              1        25.000
  --------------------------------------------------------------------------------
  --------------------------------------------------------------------------------
  Move Data for Molecule Type     1
  --------------------------------------------------------------------------------
  Conformational Moves                      Attempted       Accepted       Percent
                                                    9              2        22.222
  Bond Changes                              Attempted       Accepted       Percent
                                                    5              2        40.000
  --------------------------------------------------------------------------------
  Angle Changes                             Attempted       Accepted       Percent
                                                    4              0         0.000
  --------------------------------------------------------------------------------
  Conformational Moves Rejected BecauseBox Was Empty:     0
  -------------------------------------------------------------------------------
  Translation Moves                         Attempted       Accepted       Percent
                                                   12              1         8.333
  --------------------------------------------------------------------------------
  Rotation Moves                            Attempted       Accepted       Percent
                                                   11              1         9.091
  --------------------------------------------------------------------------------
  Biased Move Data
  --------------------------------------------------------------------------------
  Bond Changes                              Attempted       Accepted       Percent
                                                    5              5       100.000
  --------------------------------------------------------------------------------
  Angle Changes                             Attempted       Accepted       Percent
                                                    4              4       100.000
  --------------------------------------------------------------------------------
  Translation Moves                         Attempted       Accepted       Percent
                                                   12              8        66.667
  --------------------------------------------------------------------------------
  Rotation Moves                            Attempted       Accepted       Percent
                                                   11              5        45.455
  --------------------------------------------------------------------------------
  ********************************************************************************
  ------------------------------------------------
  |                   BOX 2                      |
  ------------------------------------------------
  ********************************************************************************
  Average Energy [Hartrees]:                                          -17.20378753
  Average number of molecules:                                          1.00000000
  Average Volume [angstroms**3]:                                       2579.876452

  --------------------------------------------------------------------------------
  --------------------------------------------------------------------------------
  Move Data for Molecule Type     1
  --------------------------------------------------------------------------------
  Swap Moves into this box                  Attempted       Empty          Percent
                                                    1              0         0.000
                    Growths                 Attempted       Sucessful      Percent
                                                    1              1       100.000
                      Total                 Attempted       Accepted       Percent
                                                    1              0         0.000
  -------------------------------------------------------------------------------
  Biased Move Data
  --------------------------------------------------------------------------------
  ********************************************************************************
```

Additionally, other relevant information is included in this file, such as the percentage of
accepted moves.
