  METHOD FIST
  &MM
    &FORCEFIELD
      &CHARGE
        ATOM Ar
        CHARGE 0.0
      &END CHARGE
      &NONBONDED
        &LENNARD-JONES
          ATOMS Ar Ar
          EPSILON [hartree] 0.001
          RCUT [angstrom] 7.0
          SIGMA [angstrom] 2.4
        &END LENNARD-JONES
      &END NONBONDED
    &END FORCEFIELD
    &POISSON
      PERIODIC X
      &EWALD
        EWALD_TYPE NONE
      &END EWALD
    &END POISSON
  &END MM
  &SUBSYS
    &CELL
      ABC 18.0 18.0 18.0
      PERIODIC X
    &END CELL
    &COORD
      Ar -1.5 0.0 0.0
      Ar 19.5 0.0 0.0
    &END COORD
    &TOPOLOGY
      CONNECTIVITY OFF
    &END TOPOLOGY
  &END SUBSYS
