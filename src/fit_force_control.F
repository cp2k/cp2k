!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/force_env_types [1.0] *
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

  USE coefficient_types, ONLY : coeff_type, PLANEWAVES, ATOMICORBITALS
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE fist_force, ONLY : fist_force_control => force_env_types
  USE pol_force_ao, ONLY : pol_force_control_ao => force_env_types
  USE pol_force_pw, ONLY : pol_force_control_pw => force_env_types
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE md, ONLY : simulation_parameters_type, thermodynamic_type
  USE structure_types, ONLY : structure_type, interaction_type
  USE fit_parameters_types, ONLY : fit_parameter_type
  USE termination, ONLY : stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: fit_force

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** force_env_types/force [1.0] *
!!
!!   NAME
!!     force
!!
!******************************************************************************

SUBROUTINE fit_force ( struc, inter, thermo, simpar, ewald_param, &
                     box_change, fit_param, globenv)
  
  IMPLICIT NONE

! Arguments
  TYPE ( structure_type ), INTENT ( INOUT ) :: struc
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE ( interaction_type ), INTENT ( INOUT ) :: inter
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE (fit_parameter_type), INTENT(INOUT) :: fit_param
  TYPE ( global_environment_type), INTENT ( IN ) :: globenv
  LOGICAL, INTENT (in) :: box_change

!------------------------------------------------------------------------------
!locals
  thermo % ptens_aval = .FALSE.

  SELECT CASE ( fit_param % fit_type )
  CASE DEFAULT
     CALL stop_program ( 'fit_force', 'not implemented' )

  CASE ( 'AMOEBA_RI' )
     CALL fist_force_control ( struc % molecule, struc % pnode, struc % part, &
          struc % box, struc % box_ref, thermo, inter %potparm, ewald_param,  &
          box_change, globenv, struc%ll_data (1) )

  CASE ( 'FORCES_POL' )
     SELECT CASE ( struc % coef_pos ( 1 ) % in_use )
     CASE ( ATOMICORBITALS )
       CALL pol_force_control_ao ( struc % molecule, struc % pnode, struc % part, &
          struc % box, struc % box_ref, struc % drho_basis_info, &
          struc % rho0_basis_info, struc % coef_pos (1), struc % coef_vel (1), &
          struc % coef_force (1), thermo, inter%potparm, inter%empparm,  &
          ewald_param, box_change, globenv, struc%ll_data )
     CASE ( PLANEWAVES )
       CALL pol_force_control_pw ( struc % molecule, struc % pnode, struc % part, &
          struc % box, struc % box_ref, struc % drho_basis_info, &
          struc % rho0_basis_info, struc % coef_pos (1), struc % coef_vel (1), &
          struc % coef_force (1), struc % projector ( 1 ), struc % pw_grid ( 1 ), &
          thermo, inter%potparm, inter%empparm, ewald_param, box_change, globenv, &
          struc%ll_data, simpar )
     END SELECT

  CASE ( 'AMOEBA_POL' )

     CALL stop_program ( 'fit_force', 'not implemented' )

  END SELECT

END SUBROUTINE fit_force

!!*****
!******************************************************************************

END MODULE fit_force_control

!******************************************************************************
