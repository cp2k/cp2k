&GLOBAL
  PROJECT ${seed}
  PRINT_LEVEL LOW
  RUN_TYPE ENERGY
&END GLOBAL
&FORCE_EVAL
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    &KPOINTS
      SCHEME MONKHORST-PACK 2 2 2
      GAMMA_CENTERED T
      SYMMETRY T
      SYMMETRY_BACKEND SPGLIB
      EPS_SYMMETRY 1.e-8
      PARALLEL_GROUP_SIZE 0
      WAVEFUNCTIONS COMPLEX
    &END KPOINTS
    &QS
      METHOD GPW
      EPS_DEFAULT 1.e-12
    &END QS
    &MGRID
      CUTOFF 100
      REL_CUTOFF 30
    &END MGRID
    &SCF
      @IF ( ${disentangle} == T )
        ADDED_MOS 4
      @ENDIF
      EPS_SCF 1.e-9
      EPS_EIGVAL 1.e-10
      MAX_SCF 150
      SCF_GUESS ATOMIC
      CHOLESKY OFF
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.3
      &END MIXING
      &PRINT
        &RESTART OFF
        &END RESTART
      &END PRINT
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
    &PRINT
      &WANNIER90 ON
        KPOINTS_SOURCE SCF
        SEED_NAME ${seed}
        WANNIER_FUNCTIONS ${num_wann}
        @IF ( ${exclude} /= 0 )
          EXCLUDE_BANDS ${exclude}
        @ENDIF
        REUSE_SCF_MOS ${reuse}
        VALIDATE_REUSE_SCF_MOS ${reuse}
        USE_BLOCH_PHASES T
        &LIBRARY ON
          @IF ( ${disentangle} == T )
            DIS_FROZ_MAX -11.0
            DIS_NUM_ITER 2000
            DIS_CONV_TOL 1.e-12
          @ENDIF
          INITIAL_PROJECTIONS ${projections}
          NUM_ITER ${num_iter}
          CONV_TOL ${conv_tol}
          NUM_PRINT_CYCLES 100
        &END LIBRARY
      &END WANNIER90
    &END PRINT
  &END DFT
  &SUBSYS
    &CELL
      A 0 ${half_a} ${half_a}
      B ${half_a} 0 ${half_a}
      C ${half_a} ${half_a} 0
    &END CELL
    &COORD
      SCALED
      ${cation} 0 0 0
      ${anion} 0.5 0.5 0.5
    &END COORD
    &KIND ${cation}
      BASIS_SET TZV2P-MOLOPT-GGA-GTH-q${q_cation}
      POTENTIAL GTH-GGA-q${q_cation}
    &END KIND
    &KIND ${anion}
      BASIS_SET TZV2P-MOLOPT-GGA-GTH-q7
      POTENTIAL GTH-GGA-q7
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
