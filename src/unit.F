!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/unit [1.0] *
!!
!!   NAME
!!     unit
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE unit

  USE kinds, ONLY : dbl

  USE mathconstants, ONLY : pi
  USE physcon, ONLY : e_charge, angstrom, a_bohr, e_mass, joule, a_mass, &
       boltzmann, massunit, bohr
  USE termination, ONLY : stop_program

  PRIVATE
  PUBLIC :: set_units, unit_convert_type

  TYPE unit_convert_type
     REAL ( dbl ) :: econv, mconv, lconv, tconv, pconv, eps0
     CHARACTER ( LEN = 10 ) :: e_label, l_label, vol_label,  &
          temp_label, pres_label,angl_label
  END TYPE unit_convert_type

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** unit/set_units [1.0] *
!!
!!   NAME
!!     unit
!!
!!   FUNCTION
!!     routine which defines internal units and conversion factors
!!     Internal units for the code are as follows:
!!     mass is in amu; energy is in kelvin; length is in angstroms
!!     IT IS ASSUMED THAT
!!     input time is in femtoseconds; output pressure in bar
!!     mass is in amu ALWAYS
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE set_units ( unit_char, units )

  IMPLICIT NONE

! Arguments
  CHARACTER ( LEN = * ), INTENT ( IN ) :: unit_char
  TYPE ( unit_convert_type ), INTENT ( OUT ) :: units

!------------------------------------------------------------------------------

  units % temp_label = '[K]'
  units % pres_label = '[bar]'
  units % angl_label = '[deg]'

  IF ( unit_char == 'KELVIN' ) THEN

! mass is amu; energy is kelvin; length is angstroms
     units % mconv = 1.0_dbl
     units % econv = 1.0_dbl
     units % lconv = 1.0_dbl
     units % tconv = 1.0E-15_dbl * angstrom / a_bohr &
          / ( SQRT ( a_mass / boltzmann ) )
     units % pconv = boltzmann * 1.0E+25_dbl
     units % eps0 = boltzmann * a_bohr / ( angstrom * e_charge ** 2 )
     units % l_label = '[A]'
     units % vol_label = '[A^3]'
     units % e_label = '[K]'

! gcunit is group code unit

  ELSE IF ( unit_char == 'GCUNIT' .OR. unit_char == 'ATOMIC' ) THEN

! mass is amu; energy is hartree; length is bohr
     units % mconv = massunit
     units % econv = boltzmann / joule
     units % lconv = bohr
     units % tconv = 1.0E-15_dbl / ( a_bohr * SQRT ( e_mass / joule ) )
     units % pconv = joule / ( a_bohr ** 3 * 1.0E+05_dbl )
     units % eps0 = joule * a_bohr / e_charge ** 2
     units % l_label = '[bohr]'
     units % vol_label = '[bohr^3]'
     units % e_label = '[hartree]'
  ELSE

     CALL stop_program ( 'unit', 'unknown unit type' )

  END IF

END SUBROUTINE set_units

!!*****
!******************************************************************************

END MODULE unit

!******************************************************************************
