!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!***** cp2k/extended_system_types [1.0] *
!!
!!   NAME
!!     extended_system_types
!!
!!   FUNCTION
!!   Lumps all possible extended system variables into one
!!   type for easy access and passing
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE extended_system_types

  USE kinds, ONLY : dbl
  USE nose, ONLY : lnhc_parameters_type, npt_info_type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: extended_system_type

!  PUBLIC :: npt_info, nhc_part, nhc_baro, nhc_coef

  TYPE extended_system_type
     TYPE ( lnhc_parameters_type ), POINTER :: nhc_part  ( : )
     TYPE ( lnhc_parameters_type ), POINTER :: nhc_baro ( : )
     TYPE ( lnhc_parameters_type ), POINTER :: nhc_coef 
     TYPE ( npt_info_type ), POINTER :: npt_info ( :, : )
  END TYPE extended_system_type

END MODULE extended_system_types

!!*****
!******************************************************************************
