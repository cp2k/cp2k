&QS
  METHOD GPW
  EPS_DEFAULT 1.0E-12
  MAP_CONSISTENT
  EXTRAPOLATION ASPC
  EXTRAPOLATION_ORDER 3
  &CDFT
    TYPE_OF_CONSTRAINT BECKE
    &OUTER_SCF ON
      EPS_SCF 1.0e-0
      TYPE BECKE_CONSTRAINT
      OPTIMIZER BISECT
      BISECT_TRUST_COUNT 8
      EXTRAPOLATION_ORDER 2
      MAX_SCF 0
      STEP_SIZE -0.1
    &END
  &END CDFT
  &BECKE_RESTRAINT
    @IF ( ${BECKE_ADJUST_SIZE} == TRUE )
      ! Defaults to false
      ADJUST_SIZE     TRUE
      ATOMIC_RADII 0.630 0.320  
    @ENDIF
    @IF ( ${BECKE_ATOMIC_CHARGES} == TRUE )
      ! Defaults to false
      ATOMIC_CHARGES  TRUE
    @ENDIF
    STRENGTH        ${BECKE_STR}
    TARGET          ${BECKE_TARGET}
    @IF ( ${BECKE_CUTOFF_ELEMENT} == TRUE )
      CUTOFF_TYPE     ELEMENT
      ! Note these values are in angstrom
      ELEMENT_CUTOFF  2.0 2.0
    @ENDIF
    @IF ( ${BECKE_GLOBAL_CUTOFF} == TRUE )
      CUTOFF_TYPE     GLOBAL
      GLOBAL_CUTOFF   2.0
    @ENDIF
    @IF ( ${BECKE_IN_MEMORY} == TRUE )
      ! Defaults to false
      IN_MEMORY       TRUE
    @ENDIF
    @IF ( ${BECKE_CAVITY_CONFINE} == TRUE )
      ! Defaults to NONE
      CONFINE         CAVITY
      EPS_CAVITY      1.0E-6
      CAVITY_SHAPE    ${BECKE_CAVITY_SHAPE}
      ! For shape single
      CAVITY_RADIUS   1.3
      CAVITY_USE_BOHR FALSE    
    @ENDIF
    @IF ( ${BECKE_SHOULD_SKIP} == TRUE )
      ! Defaults to false
      SHOULD_SKIP     TRUE
    @ENDIF
    @IF ( ${BECKE_CAVITY_PRINT} == TRUE )
      ! Defaults to false      
      CAVITY_PRINT    TRUE
    @ENDIF
    ATOMS           1..3
    COEFF           1 1 1
    CONSTRAINT_TYPE TOTAL
    &PROGRAM_RUN_INFO ON
      &EACH
        QS_SCF 1
      &END EACH
      COMMON_ITERATION_LEVELS 2
      ADD_LAST NUMERIC
      FILENAME ./${PROJECT_NAME}
    &END PROGRAM_RUN_INFO
  &END BECKE_RESTRAINT
&END QS
