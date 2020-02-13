&QS
  METHOD GPW
  EPS_DEFAULT 1.0E-12
  EXTRAPOLATION ASPC
  EXTRAPOLATION_ORDER 3
  &CDFT
    TYPE_OF_CONSTRAINT BECKE
    @IF ( ${BECKE_ATOMIC_CHARGES} == TRUE )
      ! Defaults to false
      ATOMIC_CHARGES  TRUE
    @ENDIF
    STRENGTH        ${BECKE_STR}
    ! The constraint target: sum_i coeff_i * N_i
    ! where N_i is the number of VALENCE electrons on i
    TARGET          ${BECKE_TARGET}
    @IF ( ${BECKE_FRAGMENT} == TRUE )
      &ATOM_GROUP
        ATOMS           1
        COEFF           1
        FRAGMENT_CONSTRAINT
      &END ATOM_GROUP
      FRAGMENT_A_FILE_NAME He+-noconstraint-ELECTRON_DENSITY-1_0.cube
      FRAGMENT_B_FILE_NAME H-noconstraint-ELECTRON_DENSITY-1_0.cube
    @ENDIF
    @IF ( ${BECKE_FRAGMENT} == FALSE )
      &ATOM_GROUP
        ATOMS           1..2
        COEFF           1 -1
      &END ATOM_GROUP
    @ENDIF
    &OUTER_SCF ON
      EPS_SCF 1.0E-2
      TYPE CDFT_CONSTRAINT
      OPTIMIZER ${OPTIMIZER}
      BISECT_TRUST_COUNT 8
      EXTRAPOLATION_ORDER 2
      MAX_SCF  ${MAX_SCF}
      STEP_SIZE -1.0
      &CDFT_OPT
        JACOBIAN_TYPE ${JACOBIAN_TYPE}
        BROYDEN_TYPE ${BROYDEN_TYPE}
        JACOBIAN_STEP 1.0E-4
        JACOBIAN_FREQ ${JACOBIAN_FREQ_1} ${JACOBIAN_FREQ_2}
        @IF ( ${RESTART_JACOBIAN} == TRUE )
          JACOBIAN_RESTART TRUE
          JACOBIAN_VECTOR -0.91773106138385885
        @ENDIF
        @IF ( ${SCF_CONTINUE_LS} == TRUE )
          ! Continue line search until norm no longer decreases or
          ! max steps are reached
          CONTINUE_LS TRUE
          ! We use more LS steps and a smaller update factor to
          ! demonstrate the quick exit conditions of CONTINUE_LS
          MAX_LS 8
          FACTOR_LS 0.8
        @ENDIF
      &END CDFT_OPT
    &END
    &BECKE_CONSTRAINT
      @IF ( ${BECKE_ADJUST_SIZE} == TRUE )
        ! Defaults to false
        ADJUST_SIZE     TRUE
        ATOMIC_RADII 0.460 0.320
      @ENDIF
      @IF ( ${BECKE_CUTOFF_ELEMENT} == TRUE )
        CUTOFF_TYPE     ELEMENT
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
        ! Defaults to FALSE
        CAVITY_CONFINE  TRUE
        EPS_CAVITY      1.0E-6
        CAVITY_SHAPE    ${BECKE_CAVITY_SHAPE}
        ! For shape single
        CAVITY_RADIUS   1.3
        CAVITY_USE_BOHR        FALSE
      @ENDIF
      @IF ( ${BECKE_SHOULD_SKIP} == TRUE )
        ! Defaults to false
        SHOULD_SKIP     TRUE
      @ENDIF
      @IF ( ${BECKE_CAVITY_PRINT} == TRUE )
        ! Defaults to false
        CAVITY_PRINT    TRUE
      @ENDIF
    &END BECKE_CONSTRAINT
    &PROGRAM_RUN_INFO ON
      &EACH
        QS_SCF 1
      &END EACH
      COMMON_ITERATION_LEVELS 2
      ADD_LAST NUMERIC
      FILENAME ./${PROJECT_NAME}
    &END PROGRAM_RUN_INFO
  &END CDFT
&END QS
