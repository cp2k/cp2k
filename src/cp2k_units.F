!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/cp2k_units [1.0] *
!!
!!   NAME
!!     cp2k_units
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

MODULE cp2k_units
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : pi
  USE physcon, ONLY : e_charge, angstrom, a_bohr, e_mass, joule, a_mass, &
       c_light, boltzmann, kcalmol, massunit, bohr, evolt, h_bar, permittivity
  USE termination, ONLY : stop_memory

  PRIVATE

  TYPE unit_convert_type
 !    PRIVATE  later
     REAL ( dbl ) :: mass, fs, kcal, kelvin, bar, wavenum, angstrom, eps0
     CHARACTER ( LEN = 10 ) :: e_label, l_label, vol_label,  &
                               temp_label, pres_label,       &
                               angl_label
  END TYPE unit_convert_type
  
  TYPE ( unit_convert_type ), POINTER :: units
  PUBLIC :: convert_to_cp2k_units, unit_convert_type, get_cp2k_units, time_convert, init_cp2k_units

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** cp2k_units/init_cp2k_units [1.0] *
!!
!!   NAME
!!     init_cp2k_units
!!
!!   FUNCTION
!!     routine which defines internal units and conversion factors
!!     internal units of CP2K are in atomic_units
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE init_cp2k_units ()

  IMPLICIT NONE
  INTEGER :: isos

!------------------------------------------------------------------------------

  ALLOCATE ( units, stat = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'cp2k_units', &
       'units', 1 )

  units % temp_label = '[K]'
  units % pres_label = '[bar]'
  units % angl_label = '[deg]'
  units % l_label = '[bohr]'
  units % vol_label = '[bohr^3]'
  units % e_label = '[hartree]'

! amu to mass AU
  units % mass = massunit
! fs to time AU
  units % fs = joule / 1.e15_dbl /  h_bar
! 1/cm to time AU
  units % wavenum = 1._dbl / 100._dbl / c_light * h_bar / joule
! kcal/mole to energy AU
  units % kcal = 1._dbl / kcalmol
! kelvin to energy AU
  units % kelvin = boltzmann / joule
! angstrom to length AU
  units % angstrom = bohr
! bar to pressure AU
  units % bar = joule / ( a_bohr ** 3 * 1.0E+05_dbl )
! permittivity to AU
  units % eps0 = joule * a_bohr / e_charge ** 2

END SUBROUTINE init_cp2k_units

!******************************************************************************
!!****** cp2k_units/get_cp2k_units [1.0] *
!!
!!   NAME
!!     get_cp2k_units
!!
!!   FUNCTION
!!     routine which defines internal units and conversion factors
!!     internal units of CP2K are in atomic_units
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE get_cp2k_units ( unit_struc )
  IMPLICIT NONE

  TYPE ( unit_convert_type ), POINTER :: unit_struc
  unit_struc => units

END SUBROUTINE get_cp2k_units

!******************************************************************************
!!****** units/time_convert [1.0] *
!!
!!   NAME
!!     time_convert
!!
!!   FUNCTION
!!     Converts time in AU to the internal units of CP2K
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE time_convert ( time )
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( INOUT ) :: time
  time = time / units % fs * 1.e-15_dbl / ( a_bohr * SQRT ( e_mass / joule ) )
END SUBROUTINE time_convert 


END MODULE cp2k_units

!******************************************************************************
