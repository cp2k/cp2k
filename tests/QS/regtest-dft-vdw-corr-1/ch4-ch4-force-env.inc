! -*- cp2k -*-
! common subsystem for ch4-ch4 vdw cp2k test
&FORCE_EVAL
  METHOD ${METHOD}
  &PRINT
    &FORCES ON
    &END
  &END
  &MIXED
    MIXING_TYPE GENMIX
    NGROUPS  1
    &GENERIC
      MIXING_FUNCTION E1+E2
      VARIABLES E1 E2
    &END
  &END
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME POTENTIAL
    &MGRID
      NGRIDS 5
      CUTOFF 50
    &END MGRID
    &QS
      METHOD GPW
    &END QS
    &POISSON
      PERIODIC NONE
      POISSON_SOLVER WAVELET
    &END POISSON
    &SCF
      SCF_GUESS ${SCF_GUESS}
      MAX_SCF 20
      EPS_SCF ${EPS_SCF}
      &OT ON
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF
        MAX_SCF 5
        EPS_SCF ${EPS_SCF}
      &END OUTER_SCF
    &END SCF
    &XC
      &XC_FUNCTIONAL ${DFT}
      &END XC_FUNCTIONAL
      &XC_GRID
        XC_SMOOTH_RHO NN50
        XC_DERIV NN50_SMOOTH
      &END
    &END XC
  &END DFT
  &MM
    &POISSON
      &EWALD
        EWALD_TYPE NONE
      &END
    &END
    &FORCEFIELD
      PARMTYPE OFF
      &SPLINE
        EMAX_SPLINE 10.0
      &END
      &CHARGE
        ATOM H
        CHARGE 0.0
      &END
      &CHARGE
        ATOM C
        CHARGE 0.0
      &END
      &NONBONDED
        &GENPOT
          ATOMS        C  C
          RCUT [angstrom]        14.9
          RMIN [angstrom]        0.2
          RMAX [angstrom]        15.0
          VARIABLES    R
          PARAMETERS   CSIX      RZERO    D    
          VALUES      -216.9850 6.57936 2.9688  
          FUNCTION    (CSIX/(R**6))*((1.0 - EXP(-1.0*D*((R/RZERO)**7)))**4)
        &END
        &GENPOT
          ATOMS        C  H
          RCUT         14.9
          VARIABLES    R
          FUNCTION     0.0
        &END
        &GENPOT
          ATOMS        H  H
          RCUT         14.9
          VARIABLES    R
          FUNCTION     0.0
        &END
      &END
    &END
  &END MM
  &SUBSYS
    &CELL
      ABC [angstrom] 10.0 10.0 10.0
      PERIODIC NONE
    &END
    &PRINT
      &ATOMIC_COORDINATES
      &END
    &END
    &TOPOLOGY
      COORD_FILE_FORMAT XYZ
      COORD_FILE_NAME ../sample_xyz/ch4-ch4-in.xyz
      &CENTER_COORDINATES
      &END
    &END TOPOLOGY
    &KIND H
      BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-${DFT}-q1
    &END
    &KIND C
      BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-${DFT}-q4
    &END
  &END SUBSYS
&END
