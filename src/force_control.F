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
  USE pol_force, ONLY : pol_force_control => force_control
  USE global_types, ONLY : global_environment_type
  USE tbmd_force, ONLY : tbmd_force_control => force_control
  USE kinds, ONLY : dbl
  USE md, ONLY : simulation_parameters_type, thermodynamic_type
  USE structure_types, ONLY : structure_type, interaction_type
  USE termination, ONLY : stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: force

  INTERFACE force
     MODULE PROCEDURE force_1, force_n
  END INTERFACE
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

SUBROUTINE force_1 ( struc, inter, thermo, simpar, ewald_param, box_change, &
     globenv, dft_control, kp )
  
  IMPLICIT NONE

! Arguments
  TYPE ( structure_type ), INTENT ( INOUT ) :: struc
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  LOGICAL, INTENT ( IN ) :: box_change
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
          box_change, globenv, struc%ll_data (1) )

  CASE ( 'POL' )
     CALL pol_force_control ( struc % molecule, struc % pnode, struc % part, &
          struc % box, struc % box_ref, struc % drho_basis_info, &
	  struc % rho0_basis_info, struc % coef_pos (1), struc % coef_vel (1), &
	  struc % coef_force (1), thermo, inter%potparm, ewald_param,  &
          simpar % ensemble, globenv, struc%ll_data )

  CASE ( 'TBMD' )
     CALL tbmd_force_control ( struc % molecule, struc % pnode, struc % part, &
          struc % box, struc % box_ref, thermo, inter%potparm, ewald_param, &
          box_change, globenv , struc%ll_data (1) )
     
  CASE ( 'WAVE' )
     CALL dft_force_control ( struc % molecule, struc % pnode, struc % part, &
          struc % box, thermo, inter % potparm, ewald_param, &
          box_change, struc % ll_data (1), dft_control, kp, &
          dftpar = globenv, wavepar = globenv )
     
  END SELECT

END SUBROUTINE force_1

!!*****
!******************************************************************************

SUBROUTINE force_n ( struc, inter, thermo, simpar, ewald_param, box_change, &
     globenv, dft_control, kp )
  
  IMPLICIT NONE

! Arguments
  TYPE ( structure_type ), DIMENSION ( : ), INTENT ( INOUT ) :: struc
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  LOGICAL, INTENT ( IN ) :: box_change
  TYPE ( interaction_type ), INTENT ( IN ) :: inter
  TYPE ( global_environment_type), INTENT ( IN ) :: globenv
  TYPE ( dft_control_type ), INTENT ( IN ) :: dft_control
  TYPE ( kpoint_type ), INTENT ( IN ) :: kp
  
!------------------------------------------------------------------------------
! Local
  INTEGER :: i, n
  REAL (dbl) :: e_pot
  
!------------------------------------------------------------------------------

  thermo % ptens_aval = .FALSE.

  n = SIZE ( struc )

  e_pot = 0.0_dbl
  
  DO i = 1, n
    SELECT CASE ( simpar % program )
    CASE DEFAULT
       CALL stop_program ( 'force', 'not implemented' )

    CASE ( 'FIST' )
       CALL fist_force_control ( struc ( i ) % molecule, struc ( i ) % pnode,  &
         struc ( i ) % part, struc ( i ) % box, struc ( i ) % box_ref,  &
         thermo, inter%potparm, ewald_param,  box_change, globenv,  &
         struc ( i ) % ll_data ( 1 ) )

    CASE ( 'POL' )
       CALL stop_program ( 'force', 'POL PIMD not implemented' )

    CASE ( 'TBMD' )
       CALL stop_program ( 'force', 'TBMD PIMD not implemented' )

    CASE ( 'WAVE' )
       CALL stop_program ( 'force', 'WAVE PIMD not implemented' )
     
    END SELECT

    e_pot = e_pot + thermo % pot

  END DO

  thermo % pot = e_pot

END SUBROUTINE force_n
!******************************************************************************
END MODULE force_control
!******************************************************************************
