!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/pair_potential_types [1.0] *
!!
!!   NAME
!!     pair_potential_types
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

MODULE pair_potential_types

  USE kinds,                           ONLY: dbl

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: pair_potential_type
  PUBLIC :: lj_pot_type, williams_pot_type, goodwin_pot_type

  TYPE lj_pot_type
    REAL ( dbl ) :: epsilon
    REAL ( dbl ) :: sigma6
    REAL ( dbl ) :: sigma12
  END TYPE lj_pot_type

  TYPE williams_pot_type
    REAL ( dbl ) :: a
    REAL ( dbl ) :: b
    REAL ( dbl ) :: c
  END TYPE williams_pot_type

  TYPE goodwin_pot_type
    REAL ( dbl ) :: vr0
    REAL ( dbl ) :: m, mc
    REAL ( dbl ) :: d, dc
  END TYPE goodwin_pot_type

  TYPE pair_potential_type
    CHARACTER ( LEN = 40 ) :: type
    TYPE (lj_pot_type ) :: lj
    TYPE (williams_pot_type ) :: willis
    TYPE (goodwin_pot_type ) :: goodwin
    REAL ( dbl ) :: rcutsq
    REAL ( dbl ) :: energy_cutoff
    REAL ( dbl ) :: e_cutoff_coul
    REAL ( dbl ) :: lrc
  END TYPE pair_potential_type

!!*****
!******************************************************************************

END MODULE pair_potential_types
