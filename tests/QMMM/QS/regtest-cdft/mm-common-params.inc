&MM
  &FORCEFIELD
     IGNORE_MISSING_CRITICAL_PARAMS
     &SPLINE
       RCUT_NB         1.2500000000000000E+01
       EMAX_SPLINE     1.0000000000000000E+02
     &END SPLINE
     &NONBONDED
       &LENNARD-JONES
         ATOMS OW OW
         EPSILON [kjmol]  0.63639
         SIGMA [nm]     0.315057422683
       &END LENNARD-JONES
       &LENNARD-JONES
         ATOMS OW HW
         EPSILON [kjmol]    0.3499745204439889
         SIGMA [nm]       0.17752938756375 
       &END LENNARD-JONES
       &LENNARD-JONES
         ATOMS HW HW
         EPSILON [kjmol]   0.192464
         SIGMA [nm]       0.0400013524445
       &END LENNARD-JONES
     &END NONBONDED
     &CHARGE
       ATOM OW
       CHARGE -0.834
     &END CHARGE
     &CHARGE
       ATOM HW
       CHARGE 0.417
     &END CHARGE
     &BOND
       ATOMS OW HW
       KIND  HARMONIC
       K  [kjmol*nm^-2]    376560.0
       R0 [nm]       0.09572 
     &END BOND
     &BOND
       ! Dummy
       ATOMS HW HW
       KIND  HARMONIC
       K     0.0000000000000000E+00
       R0    0.09572 
     &END BOND
     &BEND
       ATOMS HW OW HW
       KIND  HARMONIC
       K [kjmol*rad^-2] 460.24
       THETA0 [deg]     104.52 
     &END BEND
     &BEND
      ! Dummy
       ATOMS OW HW HW
       KIND  HARMONIC
       K     0.0000000000000000E+00
       THETA0     1.9111355309337910E+00
     &END BEND
  &END FORCEFIELD
  &POISSON
    &EWALD
      EWALD_TYPE spme
      ALPHA 0.44
      GMAX 30
    &END EWALD
  &END POISSON
  &PRINT 
    &FF_INFO ON
    &END FF_INFO
  &END PRINT
&END MM