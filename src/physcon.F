!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/physcon [1.0] *
!!
!!   NAME
!!     physcon
!!
!!   FUNCTION
!!      Definition of physical constants.
!!
!!      List of variables:
!!
!!      a_bohr    : Bohr radius [m].
!!      a_fine    : Fine-structure constant.
!!      angstrom  : Conversion factor [Bohr] -> [Angstrom]
!!      bohr      : Conversion factor [Angstrom] -> [Bohr].
!!      c_light   : Speed of light in vacuum [m/s].
!!      e_charge  : Elementary charge [C].
!!      e_mass    : Electron mass [kg].
!!      a_mass    : Atomic mass unit[kg].
!!      esu       : Conversion factors [Hartree] -> [esu].
!!      evolt     : Conversion factor [Hartree] -> [eV].
!!      h_planck  : Planck constant [J*s].
!!      hz        : Conversion factor [Hartree] -> [Hz].
!!      kcalmol   : Conversion factor [Hartree] -> [kcal/mol].
!!      kjmol     : Conversion factor [Hartree] -> [kJ/mol].
!!      joule     : Conversion factor [Hartree] -> [J].
!!      mhz       : Conversion factor [Hartree] -> [MHz].
!!      mu_perm   : Permeability of vacuum [N/A**2].
!!      boltzmann : Boltzman constant [J/Kelvin].
!!      n_avogadro: Avogrado constant [1/mol].
!!      rydberg   : Rydberg constant [1/m].
!!      vibfac    : Conversion factor [Hartree/Bohr**2] -> [1/cm].
!!      wavenum   : Conversion factor [Hartree] -> [1/cm] (wave numbers).
!!
!!   AUTHOR
!!     Matthias Krack
!!
!!   MODIFICATION HISTORY
!!     Adapted for use in CP2K (JGH)
!!
!!   SOURCE
!******************************************************************************
!

MODULE physcon

  USE kinds, ONLY : dbl

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: a_bohr, a_fine, angstrom, bohr, c_light, e_charge, e_mass, &
       evolt, h_planck, hz, joule, kcalmol, kjmol, mhz, mu_perm, n_avogadro, &
       rydberg, vibfac, wavenum, massunit, init_physcon, a_mass, boltzmann, &
       write_physcon

  INTEGER, PARAMETER :: maxmom = 3

  REAL ( dbl ) :: a_bohr, a_fine, angstrom, bohr, c_light, e_charge, e_mass, &
       evolt, h_planck, hz, joule, kcalmol, kjmol, mhz, mu_perm, n_avogadro, &
       rydberg, vibfac, wavenum, massunit, a_mass, boltzmann

  REAL ( dbl ), DIMENSION (maxmom) :: esu

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** physcon/init_physcon [1.0] *
!!
!!   NAME
!!     init_physcon
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   CREATION DATE
!!     14-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!     Fundamental physical constants (SI units)
!!     Literature: The 1986 CODATA Recommended Values of the
!!                 Fundamental Physical Constants,
!!                 E. R. Cohen, B. N. Taylor,
!!                 J. Res. Natl. Bur. Stand. 92(2), 1 (1987)
!!   SOURCE
!******************************************************************************

SUBROUTINE init_physcon()

  USE mathconstants, ONLY : pi
  USE periodic_table, ONLY : init_periodic_table

  IMPLICIT NONE

! Locals
  INTEGER :: i

!-----------------------------------------------------------------------------!

!   *** Speed of light in vacuum [m/s] ***
  c_light = 2.99792458E+8_dbl

!   *** Permeability of vacuum [N/A**2] ***
  mu_perm = 4.0_dbl * pi * 1.0E-7_dbl

!   *** Planck constant [J*s] ***
  h_planck = 6.6260755E-34_dbl

!   *** Elementary charge [C] ***
  e_charge = 1.60217733E-19_dbl

!   *** Electron mass [kg] ***
  e_mass = 9.1093897E-31_dbl

!   *** Atomic mass unit [kg] ***
  a_mass = 1.6605402E-27_dbl

!   *** Avogradro constant [1/mol] ***
  n_avogadro = 6.0221367E23_dbl

!   *** Fine-structure constant ***
  a_fine = 0.5_dbl * mu_perm * c_light * e_charge ** 2 / h_planck

!   *** Rydberg constant [1/m] ***
  rydberg = 0.5_dbl * e_mass * c_light * a_fine ** 2 / h_planck

!   *** Bohr radius [m] ***
  a_bohr = a_fine / ( 4.0_dbl * pi * rydberg )

!   *** Conversion factors ***

!   *** [Bohr] -> [Angstrom] ***
  angstrom = 1.0E+10_dbl * a_bohr

!   *** [Angstrom] -> [Bohr] ***
  bohr = 1.0_dbl / angstrom

!   *** [Hartree] -> [J] ***
  joule = 2.0_dbl * rydberg * h_planck * c_light

!   *** [Hartree] -> [kJ/mol] ***
  kjmol = 0.001_dbl * joule * n_avogadro

!   *** [Hartree] -> [kcal/mol] ***
  kcalmol = kjmol / 4.184_dbl

!   *** [Hartree] -> [eV] ***
  evolt = joule / e_charge

!   *** [Hartree] -> [Hz] ***
  hz = joule / h_planck
  mhz = 1.0E-6_dbl * hz

!   *** [Hartree] -> [1/cm] (wave numbers) ***
  wavenum = 0.02_dbl * rydberg

!   *** [Hartree/Bohr**2] -> [1/cm] (wave numbers) ***
  vibfac = 5.0_dbl * SQRT ( kjmol ) / ( pi * a_bohr * c_light )

!   *** [Hartree] -> [esu] (electrostatic units) ***
  esu(1) = 1.0E+21_dbl * a_bohr * c_light * e_charge
  DO i = 2, maxmom
     esu ( i ) = esu ( i - 1 )/ bohr
  END DO

!   *** Electron mass [amu] ***
  massunit = a_mass / e_mass

!   *** Boltzmann constant [Joule/Kelvin] ***
  boltzmann = 1.380658E-23_dbl

!   *** Periodic Table ***
  CALL init_periodic_table()

END SUBROUTINE init_physcon

!******************************************************************************

SUBROUTINE write_physcon ( iunit )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: iunit

!------------------------------------------------------------------------------

  WRITE (iunit,'(9x,A )' ) &
       ' *** Fundamental physical constants (SI units)             *** '
  WRITE (iunit,'(9x,A )' ) &
       ' *** Literature: The 1986 CODATA Recommended Values of the *** '
  WRITE (iunit,'(9x,A )' ) &
       ' ***             Fundamental Physical Constants,           *** '
  WRITE (iunit,'(9x,A )' ) &
       ' ***             E. R. Cohen, B. N. Taylor,                *** '
  WRITE (iunit,'(9x,A )' ) &
       ' ***             J. Res. Natl. Bur. Stand. 92(2), 1 (1987) *** '
  WRITE (iunit,'(A )' )
  WRITE (iunit,'(A,T61,G20.14 )' ) ' Bohr radius [m].', a_bohr
  WRITE (iunit,'(A,T61,G20.14 )' ) ' Fine-structure constant.', a_fine
  WRITE (iunit,'(A,T61,F20.14 )' ) &
       ' Conversion factor [Bohr] -> [Angstrom]', angstrom
  WRITE (iunit,'(A,T61,F20.14 )' ) &
       ' Conversion factor [Angstrom] -> [Bohr].', bohr
  WRITE (iunit,'(A,T61,F20.0 )' ) ' Speed of light in vacuum [m/s].', &
       c_light
  WRITE (iunit,'(A,T61,G20.9 )' ) ' Elementary charge [C].', e_charge
  WRITE (iunit,'(A,T61,G20.9 )' ) ' Atomic mass unit[kg].', a_mass
  WRITE (iunit,'(A,T61,G20.9 )' ) ' Electron mass [kg].', e_mass
  WRITE (iunit,'(A,T61,F20.5 )' ) ' Electron mass [amu].', massunit
  WRITE (iunit,'(A,T61,F20.14 )' ) &
       ' Conversion factors [Hartree] -> [esu].', esu(1)
  WRITE (iunit,'(A,T61,F20.14 )' ) ' Conversion factor [Hartree] -> [eV].' &
       , evolt
  WRITE (iunit,'(A,T61,G20.8 )' ) ' Planck constant [J*s].', h_planck
  WRITE (iunit,'(A,T61,G20.8 )' ) ' Boltzmann constant [J/Kelvin].', &
       boltzmann
  WRITE (iunit,'(A,T61,G20.14 )' ) ' Conversion factor [Hartree] -> [Hz].' &
       , hz
  WRITE (iunit,'(A,T61,F20.14 )' ) &
       ' Conversion factor [Hartree] -> [kcal/mol].', kcalmol
  WRITE (iunit,'(A,T61,F20.14 )' ) &
       ' Conversion factor [Hartree] -> [kJ/mol].', kjmol
  WRITE (iunit,'(A,T61,G20.14 )' ) ' Conversion factor [Hartree] -> [J].', &
       joule
  WRITE (iunit,'(A,T61,F20.4 )' ) ' Conversion factor [Hartree] -> [MHz].' &
       , mhz
  WRITE (iunit,'(A,T61,G20.14 )' ) ' Permeability of vacuum [N/A**2].', &
       mu_perm
  WRITE (iunit,'(A,T61,E20.9 )' ) ' Avogrado constant [1/mol].', &
       n_avogadro
  WRITE (iunit,'(A,T61,F20.8 )' ) ' Rydberg constant [1/m].', rydberg
  WRITE (iunit,'(A,T61,F20.10 )' ) &
       ' Conversion factor [Hartree/Bohr**2] -> [1/cm].', vibfac
  WRITE (iunit,'(A,T61,F20.10 )' ) &
       ' Conversion factor [Hartree] -> [1/cm] (wave numbers).', wavenum
  WRITE ( iunit, '()' )

END SUBROUTINE write_physcon

!!*****
!******************************************************************************

END MODULE physcon

!******************************************************************************
