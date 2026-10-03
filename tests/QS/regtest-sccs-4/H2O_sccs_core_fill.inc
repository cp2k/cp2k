# SCCS with CORE_FILL for H2O with the O nucleus exactly on a grid point (origin).
# Without CORE_FILL the nearly vanishing GTH pseudo-valence density at the O nucleus
# turns that grid point into solvent ("solvent bubble"), which gives a large spurious
# positive polarisation energy (E_total ~ -15.70 instead of ~ -17.23 Ha) and forces
# that are inconsistent with the energy.
&GLOBAL
  PRINT_LEVEL low
  PROJECT ${PROJECT}
  RUN_TYPE ${RUN}
&END GLOBAL

&DEBUG
  CHECK_ATOM_FORCE 1 Y
  DEBUG_FORCES yes
  DEBUG_STRESS_TENSOR no
  DX 0.002
  EPS_NO_ERROR_CHECK 1.0E-5
  STOP_ON_MISMATCH off
&END DEBUG

&FORCE_EVAL
  METHOD Quickstep
  # STRESS_TENSOR analytical
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &MGRID
      CUTOFF 240
      NGRIDS 4
      &RS_GRID
        DISTRIBUTION_TYPE distributed
      &END RS_GRID
    &END MGRID
    &PRINT
      &SCCS ${SCCS}
        &DENSITY_GRADIENT off
          &EACH
            QS_SCF 0
          &END EACH
        &END DENSITY_GRADIENT
        &DIELECTRIC_FUNCTION off
          &EACH
            QS_SCF 0
          &END EACH
        &END DIELECTRIC_FUNCTION
        &EACH
          QS_SCF 1
        &END EACH
        &POLARISATION_POTENTIAL off
          &EACH
            QS_SCF 0
          &END EACH
        &END POLARISATION_POTENTIAL
      &END SCCS
    &END PRINT
    &QS
      EPS_DEFAULT 1.0E-8
    &END QS
    &SCCS ${SCCS}
      ALPHA [N*m^-1] 0.0
      BETA [kbar] 0.0
      DELTA_RHO 2.0E-5
      DERIVATIVE_METHOD ${DM}
      DIELECTRIC_CONSTANT 78.36
      EPS_SCCS ${EPS_SCCS}
      GAMMA [mN/m] 0.0
      @IF ${OT}
      EPS_SCF 0.03
      @ENDIF
      @IF ${TD}
      EPS_SCF 0.3
      @ENDIF
      MAX_ITER 100
      METHOD ${METHOD}
      CORE_FILL on
      CORE_FILL_SPREAD [bohr] 0.5
      @IF ${KIND_SPREAD}
      &CORE_FILL_KIND
        H 0.5
        O 0.4
      &END CORE_FILL_KIND
      @ENDIF
      MIXING 0.6
      &ANDREUSSI
        RHO_MAX 0.001
        RHO_MIN 0.0001
      &END ANDREUSSI
      &FATTEBERT-GYGI
        BETA 1.3
        RHO_ZERO 0.0004
      &END FATTEBERT-GYGI
    &END SCCS
    &SCF
      @IF ${TD}
      EPS_SCF ${EPS_SCF}
      MAX_SCF 300
      SCF_GUESS atomic
      &MIXING on
        ALPHA 0.2
      &END MIXING
      @ENDIF
      @IF ${OT}
      EPS_SCF 1.0E-4
      MAX_SCF 31
      SCF_GUESS atomic
      &OT on
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF on
        EPS_SCF 1.0E-4
        MAX_SCF 50
      &END OUTER_SCF
      @ENDIF
      &PRINT
        &RESTART
          BACKUP_COPIES 0
          &EACH
            QS_SCF 10
          &END EACH
        &END RESTART
      &END PRINT
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 8.0 8.0 8.0
    &END CELL
    &COORD
      O         0.0000000000        0.0000000000        0.0000000000
      H         0.7600000000        0.5860000000        0.0000000000
      H        -0.7600000000        0.5860000000        0.0000000000
    &END COORD
    &KIND H
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE-q1
    &END KIND
    &KIND O
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE-q6
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
