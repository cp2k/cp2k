&MOTION
  &MD
    ENSEMBLE ${ENSEMBLE}
    STEPS ${STEPS}
    TEMPERATURE 150.0
    TIMESTEP ${TIMESTEP}
    &PRINT ON
      &ENERGY
        FILENAME =${BASE_NAME}-${ID}.ener
        &EACH
          MD 1
        &END EACH
      &END ENERGY
    &END PRINT
    &THERMOSTAT
      REGION MASSIVE
      &NOSE
        LENGTH 6
        MTS 2
        TIMECON [wavenumber_t] 2300
        YOSHIDA 3
      &END NOSE
    &END THERMOSTAT
  &END MD
  &PRINT
    &FORCES OFF
      FILENAME =${BASE_NAME}-${ID}_force.xyz
      FORMAT XYZ
      &EACH
        MD 2
      &END EACH
    &END FORCES
    &RESTART ON
      FILENAME =${BASE_NAME}-${ID}.restart
      &EACH
        MD 1000
      &END EACH
    &END RESTART
    &TRAJECTORY ON
      FILENAME =${BASE_NAME}-${ID}.xyz
      FORMAT XYZ
      &EACH
        MD 2
      &END EACH
    &END TRAJECTORY
    &VELOCITIES OFF
      FILENAME =${BASE_NAME}-${ID}_vel.xyz
      FORMAT XYZ
      &EACH
        MD 2
      &END EACH
    &END VELOCITIES
  &END PRINT
&END MOTION
