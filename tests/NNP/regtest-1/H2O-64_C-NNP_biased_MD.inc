&GLOBAL
  PRINT_LEVEL LOW
  PROJECT H2O-64_C-NNP_biased_MD
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
    &BIAS
      SIGMA_0 [hartree] 0.0001
      K_B 5000
      ALIGN_NNP_ENERGIES 3.7512433673534811E+03 3.7512444099257964E+03 3.7512439185529488E+03 3.7512431538151554E+03 \
      3.7512424733561511E+03    3.7512436390725293E+03    3.7512420236900825E+03    3.7512430756277122E+03
      &PRINT
        &BIAS_ENERGY
          &EACH
            MD 1
          &END EACH
        &END BIAS_ENERGY
      &END PRINT
    &END BIAS
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
      &FORCES
        &EACH
          MD 1
        &END EACH
      &END FORCES
      &FORCES_SIGMA
        &EACH
          MD 1
        &END EACH
      &END FORCES_SIGMA
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
