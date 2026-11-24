  &MM
    &FORCEFIELD
      PARMTYPE CHM
      PARM_FILE_NAME ../sample_pot/deca_ala.pot
    &END FORCEFIELD
    &POISSON
      &EWALD
        ALPHA .36
        EWALD_TYPE ewald
        GMAX 29
      &END EWALD
    &END POISSON
    &PRINT
      &FF_INFO off
        SPLINE_DATA
        SPLINE_INFO
      &END FF_INFO
    &END PRINT
  &END MM
