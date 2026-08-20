&GLOBAL
  PRINT_LEVEL LOW
  PROJECT ${PROJ-H2O-256_C-NNP_MD-bitexact}
  RUN_TYPE MD
&END GLOBAL

&MOTION
  &MD
    ENSEMBLE NVT
    STEPS 5
    TEMPERATURE 300
    TIMESTEP 0.5
    &PRINT
      FORCE_LAST
      &ENERGY
        &EACH
          MD 50
        &END EACH
      &END ENERGY
    &END PRINT
    &THERMOSTAT
      REGION GLOBAL
      TYPE CSVR
      &CSVR
        TIMECON 30.00
      &END CSVR
    &END THERMOSTAT
  &END MD
  &PRINT
    &TRAJECTORY
      &EACH
        MD 500
      &END EACH
    &END TRAJECTORY
  &END PRINT
&END MOTION

&FORCE_EVAL
  METHOD NNP
  &NNP
    GPU_PRECISION ${PREC-DOUBLE}
    NNP_INPUT_FILE_NAME NNP/bulkH2O-jcp2020-cnnp/nnp-1/input.nn
    SCALE_FILE_NAME NNP/bulkH2O-jcp2020-cnnp/nnp-1/scaling.data
    USE_GPU ${GPU-OFF}
    &MODEL
      WEIGHTS NNP/bulkH2O-jcp2020-cnnp/nnp-1/weights
    &END MODEL
    &MODEL
      WEIGHTS NNP/bulkH2O-jcp2020-cnnp/nnp-2/weights
    &END MODEL
    &MODEL
      WEIGHTS NNP/bulkH2O-jcp2020-cnnp/nnp-3/weights
    &END MODEL
    &MODEL
      WEIGHTS NNP/bulkH2O-jcp2020-cnnp/nnp-4/weights
    &END MODEL
    &MODEL
      WEIGHTS NNP/bulkH2O-jcp2020-cnnp/nnp-5/weights
    &END MODEL
    &MODEL
      WEIGHTS NNP/bulkH2O-jcp2020-cnnp/nnp-6/weights
    &END MODEL
    &MODEL
      WEIGHTS NNP/bulkH2O-jcp2020-cnnp/nnp-7/weights
    &END MODEL
    &MODEL
      WEIGHTS NNP/bulkH2O-jcp2020-cnnp/nnp-8/weights
    &END MODEL
    &PRINT
      &ENERGIES
        &EACH
          MD 100
        &END EACH
      &END ENERGIES
    &END PRINT
  &END NNP
  &SUBSYS
    &CELL
      ABC [angstrom] 24.84000000 24.84000000 12.42000000
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_FORMAT XYZ
      COORD_FILE_NAME ../coords/H2O-256.xyz
    &END TOPOLOGY
  &END SUBSYS
&END FORCE_EVAL
