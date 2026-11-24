! -*- cp2k -*-
! common subsystem for ch4-ch4 vdw cp2k test
&FORCE_EVAL
  METHOD ${METHOD}
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME POTENTIAL
    &MGRID
      CUTOFF 50
      NGRIDS 5
    &END MGRID
    &POISSON
      PERIODIC NONE
      POISSON_SOLVER WAVELET
    &END POISSON
    &QS
      METHOD GPW
    &END QS
    &SCF
      EPS_SCF ${EPS_SCF}
      MAX_SCF 20
      SCF_GUESS ${SCF_GUESS}
      &OT ON
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF
        EPS_SCF ${EPS_SCF}
        MAX_SCF 5
      &END OUTER_SCF
    &END SCF
    &XC
      &XC_FUNCTIONAL ${DFT}
      &END XC_FUNCTIONAL
      &XC_GRID
        XC_DERIV NN50_SMOOTH
        XC_SMOOTH_RHO NN50
      &END XC_GRID
    &END XC
  &END DFT
  &MIXED
    MIXING_TYPE GENMIX
    NGROUPS 1
    &GENERIC
      MIXING_FUNCTION E1+E2
      VARIABLES E1 E2
    &END GENERIC
  &END MIXED
  &MM
    &FORCEFIELD
      PARMTYPE OFF
      &CHARGE
        ATOM H
        CHARGE 0.0
      &END CHARGE
      &CHARGE
        ATOM C
        CHARGE 0.0
      &END CHARGE
      &NONBONDED
        &GENPOT
          ATOMS C C
          FUNCTION (CSIX/(R**6))*((1.0 - EXP(-1.0*D*((R/RZERO)**7)))**4)
          PARAMETERS CSIX RZERO D
          RCUT [angstrom] 14.9
          RMAX [angstrom] 15.0
          RMIN [angstrom] 0.2
          VALUES -216.9850 6.57936 2.9688
          VARIABLES R
        &END GENPOT
        &GENPOT
          ATOMS C H
          FUNCTION 0.0
          RCUT 14.9
          VARIABLES R
        &END GENPOT
        &GENPOT
          ATOMS H H
          FUNCTION 0.0
          RCUT 14.9
          VARIABLES R
        &END GENPOT
      &END NONBONDED
      &SPLINE
        EMAX_SPLINE 10.0
      &END SPLINE
    &END FORCEFIELD
    &POISSON
      &EWALD
        EWALD_TYPE NONE
      &END EWALD
    &END POISSON
  &END MM
  &PRINT
    &FORCES ON
    &END FORCES
  &END PRINT
  &SUBSYS
    &CELL
      ABC [angstrom] 10.0 10.0 10.0
      PERIODIC NONE
    &END CELL
    &KIND H
      BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-${DFT}-q1
    &END KIND
    &KIND C
      BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-${DFT}-q4
    &END KIND
    &PRINT
      &ATOMIC_COORDINATES
      &END ATOMIC_COORDINATES
    &END PRINT
    &TOPOLOGY
      COORD_FILE_FORMAT XYZ
      COORD_FILE_NAME ../sample_xyz/ch4-ch4-in.xyz
      &CENTER_COORDINATES
      &END CENTER_COORDINATES
    &END TOPOLOGY
  &END SUBSYS
&END FORCE_EVAL
