!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/dg_types [1.0] *
!!
!!   NAME
!!     dg_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!******************************************************************************
MODULE dg_types

  USE coefficient_lists,               ONLY: coef_list_type
  USE dg_rho0s,                        ONLY: dg_rho0_type
  USE particle_lists,                  ONLY: particle_list_type

  IMPLICIT NONE

! Information on the assignment function for Ewald
  TYPE dg_type
    INTEGER :: grid_index
    INTEGER :: nparts
    INTEGER :: ncoefs
    TYPE ( dg_rho0_type ) :: dg_rho0
    TYPE ( particle_list_type ), POINTER :: plist
    TYPE ( coef_list_type ), POINTER :: clist
  END TYPE dg_type

  PRIVATE
  PUBLIC :: dg_type

!-----------------------------------------------------------------------------!

END MODULE dg_types

!-----------------------------------------------------------------------------!
!!*****
