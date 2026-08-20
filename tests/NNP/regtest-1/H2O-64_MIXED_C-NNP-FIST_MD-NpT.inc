&GLOBAL
  PRINT_LEVEL LOW
  PROJECT ${PROJ-H2O-64_MIXED_C-NNP-FIST_MD-NpT-cpu}
  RUN_TYPE MD
&END GLOBAL

&MOTION
  &MD
    ENSEMBLE NPT_F
    STEPS 2
    TEMPERATURE 300
    TIMESTEP 0.5
    &BAROSTAT
      PRESSURE 1
    &END BAROSTAT
    &PRINT
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
&END MOTION

&MULTIPLE_FORCE_EVALS
  FORCE_EVAL_ORDER 2 3
  MULTIPLE_SUBSYS F
&END MULTIPLE_FORCE_EVALS

! Committee NNP summed with a FIST point-charge baseline. Both halves supply an
! analytic virial, so this pins the mixed force environment's stress
! aggregation on whichever code path USE_GPU selects (regtest-1 runs the
! CPU path, regtest-gpu the device).
!
! Every force environment carries its own STRESS_TENSOR on purpose:
! pv_availability is derived per force environment and is NOT inherited from the
! MIXED level, so dropping it from the NNP section below would leave the NNP
! virial at zero while the MIXED level still printed a complete-looking stress
! tensor holding only the FIST term. Constant-pressure runs abort on the
! omission; NVT and NVE do not.
&FORCE_EVAL
  METHOD MIXED
  STRESS_TENSOR ANALYTICAL
  &MIXED
    MIXING_TYPE GENMIX
    &GENERIC
      MIXING_FUNCTION a+b
      VARIABLES a b
    &END GENERIC
  &END MIXED
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

&FORCE_EVAL
  METHOD NNP
  STRESS_TENSOR ANALYTICAL
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
  &END NNP
&END FORCE_EVAL

&FORCE_EVAL
  METHOD FIST
  STRESS_TENSOR ANALYTICAL
  &MM
    &FORCEFIELD
      IGNORE_MISSING_CRITICAL_PARAMS
      &CHARGE
        ATOM O
        CHARGE -0.834
      &END CHARGE
      &CHARGE
        ATOM H
        CHARGE 0.417
      &END CHARGE
    &END FORCEFIELD
    &POISSON
      &EWALD
        ALPHA 0.35
        EWALD_TYPE SPME
        GMAX 25
        O_SPLINE 6
      &END EWALD
    &END POISSON
  &END MM
&END FORCE_EVAL
