&FORCE_EVAL
  METHOD FIST
  STRESS_TENSOR ANALYTICAL
  &MM
    &FORCEFIELD
      &CHARGE
        ATOM Ar
        CHARGE 0.0
      &END CHARGE
      &NONBONDED
        &LENNARD-JONES
          ATOMS Ar Ar
          EPSILON [hartree] 0.01
          RCUT [angstrom] 7.0
          SIGMA [angstrom] 3.4
        &END LENNARD-JONES
      &END NONBONDED
    &END FORCEFIELD
    &POISSON
      PERIODIC ${PERIODIC}
      &EWALD
        EWALD_TYPE NONE
      &END EWALD
    &END POISSON
  &END MM
  &SUBSYS
    &CELL
      ABC ${CELL_A} 18.0 18.0
      PERIODIC ${PERIODIC}
    &END CELL
    &COORD
      Ar 0.0 0.0 0.0
      Ar ${R} 0.0 0.0
    &END COORD
    &TOPOLOGY
      CONNECTIVITY OFF
    &END TOPOLOGY
  &END SUBSYS
&END FORCE_EVAL
