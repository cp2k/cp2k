    # CP2K_SECTION_PATH: FORCE_EVAL/DFT
    &QS
      EPS_DEFAULT 1.0E-12
      EXTRAPOLATION ASPC
      EXTRAPOLATION_ORDER 3
      METHOD GPW
      &CDFT
        STRENGTH ${HIRSHFELD_STR}
        TARGET ${HIRSHFELD_TARGET}
        TYPE_OF_CONSTRAINT HIRSHFELD
        &ATOM_GROUP
          ATOMS 1 2
          COEFF 1 -1
        &END ATOM_GROUP
        &HIRSHFELD_CONSTRAINT
          SHAPE_FUNCTION DENSITY
        &END HIRSHFELD_CONSTRAINT
        &OUTER_SCF ON
          EPS_SCF 1.0E-3
          MAX_SCF 50
          OPTIMIZER NEWTON_LS
          TYPE CDFT_CONSTRAINT
          &CDFT_OPT ON
            CONTINUE_LS TRUE
          &END CDFT_OPT
        &END OUTER_SCF
        &PROGRAM_RUN_INFO ON
          ! Print weight function to cube file
          @IF ( ${PRINT_WEIGHT_FUNCTION} == TRUE )
          &WEIGHT_FUNCTION
          &END WEIGHT_FUNCTION
          @ENDIF
          &EACH
            QS_SCF 1
          &END EACH
          COMMON_ITERATION_LEVELS 2
          ADD_LAST NUMERIC
          FILENAME ./${PROJECT_NAME}
        &END PROGRAM_RUN_INFO
      &END CDFT
    &END QS
