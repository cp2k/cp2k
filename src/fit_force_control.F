!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/force_control [1.0] *
!!
!!   NAME
!!     fit_force_control
!!
!!   FUNCTION
!!     Interface routine for the fit calculations
!!
!!   SOURCE
!******************************************************************************

MODULE fit_force_control

  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE fist_force, ONLY : fist_force_control => force_control
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE md, ONLY : simulation_parameters_type, thermodynamic_type
  USE structure_types, ONLY : structure_type
  USE pair_potential, ONLY : potentialparm_type 
  USE fit_parameters_types, ONLY : fit_parameter_type
  USE termination, ONLY : stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: fit_force

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** force_control/force [1.0] *
!!
!!   NAME
!!     force
!!
!******************************************************************************

SUBROUTINE fit_force ( struc, potparm, thermo, ewald_param, &
                     fit_param, globenv)
  
  IMPLICIT NONE

! Arguments
  TYPE ( structure_type ), INTENT ( INOUT ) :: struc
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE (potentialparm_type), DIMENSION (:,:), INTENT (INOUT) :: potparm
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE (fit_parameter_type), INTENT(INOUT) :: fit_param
  TYPE ( global_environment_type), INTENT ( IN ) :: globenv
!------------------------------------------------------------------------------
!locals
  LOGICAL, PARAMETER :: box_change = .FALSE.

  thermo % ptens_aval = .FALSE.

  SELECT CASE ( fit_param % fit_type )
  CASE DEFAULT
     CALL stop_program ( 'fit_force', 'not implemented' )

  CASE ( 'AMOEBA_RI' )
     CALL fist_force_control ( struc % molecule, struc % pnode, struc % part, &
          struc % box, struc % box_ref, thermo, potparm, ewald_param,  &
          box_change, globenv, struc%ll_data (1) )

  CASE ( 'AMOEBA_POL' )

     CALL stop_program ( 'fit_force', 'not implemented' )
     
  END SELECT

END SUBROUTINE fit_force

!!*****
!******************************************************************************

END MODULE fit_force_control

!******************************************************************************
