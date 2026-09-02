# Dilute electrolyte: 62 waters + 1 Na + 1 Cl in the 12.42 A cube. With a
# single ion pair 10.47 A apart (minimum image, cutoff 12 bohr = 6.35 A) the
# heterogeneous angular symmetry-function groups (Na; O,Cl), (Cl; O,Na),
# (Na; Cl,Cl), ... have an empty first and/or second neighbour list on every
# step. This pins the empty-hetero-angular-group handling on the GPU path: a
# group skipped at symmetry-function time is also skipped at force-scatter
# time, since scattering a slot the descriptor stage never wrote would read
# stale device buffers and produce wrong forces. The regtest-1 run pins the
# reference trajectory on the CPU path; the regtest-gpu twin runs this same
# body on the device.
&GLOBAL
  PRINT_LEVEL LOW
  PROJECT NaCl-dilute-hetero
  RUN_TYPE MD
&END GLOBAL

&MOTION
  &MD
    ENSEMBLE NVT
    STEPS 3
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
&END MOTION

&FORCE_EVAL
  METHOD NNP
  &NNP
    NNP_INPUT_FILE_NAME NNP/NaClaq-cnnp/train_001/input.nn
    SCALE_FILE_NAME NNP/NaClaq-cnnp/train_001/scaling.data
    USE_GPU ${GPU-OFF}
    &MODEL
      WEIGHTS NNP/NaClaq-cnnp/train_001/weights
    &END MODEL
  &END NNP
  &SUBSYS
    &CELL
      ABC [angstrom] 12.42 12.42 12.42
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_FORMAT XYZ
      COORD_FILE_NAME ../coords/NaCl-dilute.xyz
    &END TOPOLOGY
  &END SUBSYS
&END FORCE_EVAL
