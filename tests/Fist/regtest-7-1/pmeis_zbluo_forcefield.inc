    &FORCEFIELD
      MULTIPLE_POTENTIAL
      &CHARGE
        ATOM U
        CHARGE 4.0
      &END
      &CHARGE
        ATOM O
        CHARGE -2.0
      &END
      &SHELL U
        CORE_CHARGE   -2.840000
        SHELL_CHARGE 6.840000
        MASS_FRACTION  0.01
        SPRING  [eV*angstrom^-2] 171.556
        MAX_DISTANCE [angstrom] 0.8
      &END SHELL
      &SHELL O
        CORE_CHARGE 1.186267
        SHELL_CHARGE -3.186267
        MASS_FRACTION 0.1
        SPRING [eV*angstrom^-2] 70.824
        MAX_DISTANCE [angstrom] 0.8
      &END SHELL
      &NONBONDED
        &WILLIAMS
          atoms U  O
          A [eV] 844.41
          B [angstrom^-1] 2.35159661652
          C [eV*angstrom^6] 0.0
          RCUT 10.4
    RMIN [angstrom] 0.8
        &END WILLIAMS
  &GENPOT
          atoms U O
    FUNCTION P/X*(a1*EXP(-3.2*X/AU)+a2*EXP(-0.9423*X/AU)+a3*EXP(-0.4029*X/AU)+a4*EXP(-0.2016*X/AU)) 
    PARAMETERS P AU a1 a2 a3 a4
    VALUES  51.11239406426095344911  .37659841887685972746 0.1818 0.5099 0.2802  0.02817
    VARIABLES X
    RMAX [angstrom]  0.4
          RCUT 10.4
  &END
  &GENPOT
          atoms U O
    FUNCTION a0 + a1*X + a2*X**2 + a3 *X**3 + a4* X**4 + a5 * X**5
    PARAMETERS a0 a1 a2 a3 a4  a5
    VALUES 103.31962895    -182.45077437       32.96257294      140.94822916     -111.18252029       24.77706269 
    VARIABLES X
    RMIN [angstrom]  0.4
    RMAX [angstrom]  0.8
          RCUT 10.4
  &END
        &BUCK4RANGES
          atoms O O
          A [eV] 20908.03
          B [angstrom^-1] 7.71604938272
          C [eV*angstrom^6] 229.04
          R1 1.17
          R2 1.62
          R3 2.84
          POLY1 [eV] 4204.472427467161  -7959.112702279119   6021.074988330987  -2270.766871550818 426.319753380798 -31.850801875744
          POLY2 [eV] 6.734634669569  -6.692698089020 1.656245757939  -0.122636741736
          RCUT 10.4
        &END BUCK4RANGES
        &WILLIAMS
          atoms U U
          A [eV] 0.0
          B [angstrom^-1] 0.0
          C [eV*angstrom^6] 0.0
          RCUT 10.4
        &END WILLIAMS
      &END NONBONDED
      &SPLINE
        R0_NB   0.05
        EMAX_SPLINE 100
      &END
    &END FORCEFIELD
