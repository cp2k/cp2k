!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/ewald_parameters_types [1.0] *
!!
!!   NAME
!!     ewald_parameters_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE ewald_parameters_types

  USE kinds, ONLY : dbl

  IMPLICIT NONE

  TYPE ewald_parameters_type
    CHARACTER( LEN = 12 ) :: ewald_type  ! type of ewald
    CHARACTER( LEN = 12 ) :: unit_type   !  unit_type
    REAL ( dbl ) :: alpha               ! ewald alpha
    REAL ( dbl ) :: eps0                ! permittivity constant
    REAL ( dbl ) :: epsilon             ! tolerance for small grid (PME)
    INTEGER   :: gmax                ! max Miller index
    INTEGER   :: ns_max              ! # grid points for small grid (PME)
    INTEGER   :: o_spline            ! order of spline (SPME)
    INTEGER   :: gtot                ! total number of g-vectors
    INTEGER   :: group               ! parallel group
  END TYPE ewald_parameters_type

  PRIVATE
  PUBLIC :: ewald_parameters_type

END MODULE ewald_parameters_types

!******************************************************************************
!!*****
