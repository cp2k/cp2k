!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/force_control [1.0] *
!!
!!   NAME
!!     force_control
!!
!!   FUNCTION
!!     Interface routine for the force calculations
!!
!!   AUTHOR
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     cjm, FEB-20-2001: pass variable box_ref
!!
!!   SOURCE
!******************************************************************************

MODULE force_control

  USE brillouin, ONLY : kpoint_type
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE dft_structure_control, ONLY : dft_force_control
  USE dft_types, ONLY : dft_control_type
  USE fist_force, ONLY : fist_force_control => force_control
  USE global_types, ONLY : global_environment_type
  USE tbmd_force, ONLY : tbmd_force_control => force_control
  USE kinds, ONLY : dbl
  USE md, ONLY : simulation_parameters_type, thermodynamic_type
  USE structure_types, ONLY : structure_type, interaction_type
  USE termination, ONLY : stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: force

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** force_control/force [1.0] *
!!
!!   NAME
!!     force
!!
!!   FUNCTION
!!     Interface routine for force calculations
!!
!!   AUTHOR
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     Harald Forbert (Dec-2000): Changes for multiple linked lists
!!     JGH (21-12-2000): Flag to test if pressure is available
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE force ( struc, inter, thermo, simpar, ewald_param, globenv, &
     dft_control, kp )
  
  IMPLICIT NONE

! Arguments
  TYPE ( structure_type ), INTENT ( INOUT ) :: struc
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE ( interaction_type ), INTENT ( IN ) :: inter
  TYPE ( global_environment_type), INTENT ( IN ) :: globenv
  TYPE ( dft_control_type ), INTENT ( IN ) :: dft_control
  TYPE ( kpoint_type ), INTENT ( IN ) :: kp
  
!------------------------------------------------------------------------------

  thermo % ptens_aval = .FALSE.

  SELECT CASE ( simpar % program )
  CASE DEFAULT
     CALL stop_program ( 'force', 'not implemented' )

  CASE ( 'FIST' )
     CALL fist_force_control ( struc % molecule, struc % pnode, struc % part, &
          struc % box, struc % box_ref, thermo, inter%potparm, ewald_param,  &
          simpar % ensemble, globenv, struc%ll_data )

  CASE ( 'TBMD' )
     CALL tbmd_force_control ( struc % molecule, struc % pnode, struc % part, &
          struc % box, struc % box_ref, thermo, inter%potparm, ewald_param, &
          simpar % ensemble, globenv , struc%ll_data )
     
  CASE ( 'WAVE' )
     CALL dft_force_control ( struc % molecule, struc % pnode, struc % part, &
          struc % box, thermo, inter % potparm, ewald_param, &
          simpar % ensemble, struc % ll_data, dft_control, kp, &
          dftpar = globenv, wavepar = globenv )
     
  END SELECT

END SUBROUTINE force

!!*****
!******************************************************************************

END MODULE force_control

!******************************************************************************
