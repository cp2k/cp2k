&GLOBAL
  PRINT_LEVEL LOW
  PROJECT H2O-64_NNP_MD
  RUN_TYPE MD
&END GLOBAL

&MOTION
  &MD
    ENSEMBLE NVT
    STEPS 10
    TEMPERATURE 300
    TIMESTEP 0.5
    &PRINT
      FORCE_LAST
      &ENERGY
        &EACH
          MD 1
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
        MD 1
      &END EACH
    &END TRAJECTORY
  &END PRINT
&END MOTION

&FORCE_EVAL
  METHOD NNP
  &NNP
    NNP_INPUT_FILE_NAME NNP/bulkH2O-jcp2020-cnnp/nnp-1/input.nn
    SCALE_FILE_NAME NNP/bulkH2O-jcp2020-cnnp/nnp-1/scaling.data
    USE_GPU ${GPU-OFF}
    &MODEL
      WEIGHTS NNP/bulkH2O-jcp2020-cnnp/nnp-1/weights
    &END MODEL
    &PRINT
      &ENERGIES
        &EACH
          MD 1
        &END EACH
      &END ENERGIES
      &FORCES
        &EACH
          MD 1
        &END EACH
      &END FORCES
    &END PRINT
  &END NNP
  &SUBSYS
    &CELL
      ABC [angstrom] 12.42 12.42 12.42
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_FORMAT XYZ
      COORD_FILE_NAME ../coords/H2O-64.xyz
    &END TOPOLOGY
  &END SUBSYS
&END FORCE_EVAL
