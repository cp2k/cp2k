!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/spme [1.0] *
!!
!!   NAME
!!     spme
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

MODULE spme

  USE kinds, ONLY: dbl
  USE dg_types, ONLY : dg_type
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE pw_grid_types, ONLY : pw_grid_type
  USE simulation_cell, ONLY : cell_type
  USE termination, ONLY : stop_program, stop_memory
  USE timings, ONLY: timeset, timestop

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: spme_evaluate

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE spme_evaluate ( dg, fg_coulomb, vg_coulomb, pv_g, box, &
     grid_spme, ewald_param )

  IMPLICIT NONE

! Arguments
  TYPE ( dg_type ), INTENT ( IN ), TARGET :: dg
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( pw_grid_type ), INTENT ( IN ), TARGET :: grid_spme
  TYPE ( ewald_parameters_type ), INTENT ( IN ) :: ewald_param
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( OUT ) :: fg_coulomb
  REAL ( dbl ), INTENT ( OUT ) :: vg_coulomb
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( OUT ) :: pv_g

! Locals
  INTEGER :: handle
  
!------------------------------------------------------------------------------
  
  CALL timeset ( 'SPME', 'I', ' ', handle )
  
  CALL timestop ( 0.0_dbl, handle )

END SUBROUTINE spme_evaluate

!******************************************************************************

END MODULE spme

!******************************************************************************
