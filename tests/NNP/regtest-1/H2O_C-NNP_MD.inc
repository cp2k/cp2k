&GLOBAL
  PRINT_LEVEL LOW
  PROJECT H2O_C-NNP_MD
  RUN_TYPE MD
&END GLOBAL

&MOTION
  &MD
    ENSEMBLE NVE
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
          MD 1
        &END EACH
      &END ENERGIES
      &EXTRAPOLATION
        &EACH
          MD 1
        &END EACH
      &END EXTRAPOLATION
    &END PRINT
  &END NNP
  &SUBSYS
    &CELL
      ABC [angstrom] 12.42 12.42 12.42
      PERIODIC NONE
    &END CELL
    &TOPOLOGY
      COORD_FILE_FORMAT XYZ
      COORD_FILE_NAME ../coords/H2O-1.xyz
    &END TOPOLOGY
  &END SUBSYS
&END FORCE_EVAL
