!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000,2004  CP2K developers group                            !
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
  USE kinds,                           ONLY: dbl
  USE physcon,                         ONLY: a_bohr,&
                                             bohr,&
                                             boltzmann,&
                                             e_charge,&
                                             femtoseconds,&
                                             h_bar,&
                                             hz,&
                                             joule,&
                                             kcalmol,&
                                             massunit,&
                                             permittivity,&
                                             wavenum
  PRIVATE

  TYPE unit_convert_type
 !    PRIVATE  later
     REAL ( dbl ) :: mass, fs, kcal, kelvin, bar, wavenum, angstrom, eps0
     CHARACTER ( LEN = 10 ) :: e_label, l_label, vol_label,  &
                               temp_label, pres_label,       &
                               angl_label
  END TYPE unit_convert_type
  
  TYPE ( unit_convert_type ), TARGET, SAVE :: units
  PUBLIC :: unit_convert_type, get_cp2k_units, init_cp2k_units

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

  units % temp_label = '[K]'
  units % pres_label = '[bar]'
  units % angl_label = '[deg]'
  units % l_label = '[bohr]'
  units % vol_label = '[bohr^3]'
  units % e_label = '[hartree]'

! amu to mass AU
  units % mass = massunit
! fs to time AU
  units % fs = 1._dbl/femtoseconds
! 1/cm to time AU
  units % wavenum =  joule * wavenum / hz / h_bar
! kcal/mole to energy AU
  units % kcal = 1._dbl / kcalmol
! kelvin to energy AU
  units % kelvin = boltzmann / joule
! angstrom to length AU
  units % angstrom = bohr
! bar to pressure AU
  units % bar = joule / ( a_bohr ** 3 * 1.0E+05_dbl )
! permittivity to AU
  units % eps0 = permittivity * joule * a_bohr / e_charge ** 2

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

    TYPE(unit_convert_type), POINTER         :: unit_struc

  unit_struc => units

END SUBROUTINE get_cp2k_units

END MODULE cp2k_units

!******************************************************************************
