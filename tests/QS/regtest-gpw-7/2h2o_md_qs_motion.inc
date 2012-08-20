&MOTION
  &MD
    ENSEMBLE ${ENSEMBLE}
    STEPS ${STEPS}
    TIMESTEP ${TIMESTEP}
    TEMPERATURE 150.0
    &THERMOSTAT
      REGION MASSIVE 
      &NOSE
        LENGTH 6
        YOSHIDA 3
        TIMECON  [wavenumber_t] 2300
        MTS 2
      &END NOSE
    &END
    &PRINT ON
      &ENERGY
        &EACH
          MD 1
        &END EACH
        FILENAME =${BASE_NAME}-${ID}.ener
      &END ENERGY
    &END PRINT
  &END MD
  &PRINT
    &TRAJECTORY ON
      &EACH
        MD 2
      &END EACH
      FILENAME =${BASE_NAME}-${ID}.xyz 
      FORMAT XYZ
    &END TRAJECTORY
    &VELOCITIES OFF
      &EACH
        MD 2
      &END EACH
      FILENAME =${BASE_NAME}-${ID}_vel.xyz 
      FORMAT XYZ
    &END VELOCITIES
    &FORCES OFF
      &EACH
        MD 2
      &END EACH
      FILENAME =${BASE_NAME}-${ID}_force.xyz 
      FORMAT XYZ
    &END FORCES
    &RESTART ON
      &EACH
        MD 1000
      &END EACH
      FILENAME =${BASE_NAME}-${ID}.restart 
    &END RESTART
  &END PRINT
&END MOTION
