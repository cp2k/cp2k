!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/structure_factor_types [1.0] *
!!
!!   NAME
!!     structure_factor_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!!
!******************************************************************************

MODULE structure_factor_types

  USE kinds,                           ONLY: dp

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: structure_factor_type

  TYPE structure_factor_type
     COMPLEX (KIND=dp), DIMENSION ( :, : ), POINTER :: ex, ey, ez
     INTEGER, DIMENSION ( :, : ), POINTER :: centre
     INTEGER :: lb ( 3 )
  END TYPE structure_factor_type

END MODULE structure_factor_types

!!*****
!******************************************************************************
