!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE force_control

  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE fist_force, ONLY : fist_force_control => force_control
  USE fist_global, ONLY : fistpar
  USE tbmd_force, ONLY : tbmd_force_control => force_control
  USE tbmd_global, ONLY : tbmdpar
  USE kinds, ONLY : dbl
  USE md, ONLY : simulation_parameters_type, thermodynamic_type
  USE pair_potential, ONLY : potentialparm_type
  USE stop_program, ONLY : stop_prg
  USE structure_types, ONLY : structure_type
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: force
  
CONTAINS

!******************************************************************************

SUBROUTINE force ( struc, potparm, thermo, simpar, ewald_param )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( structure_type ), INTENT ( INOUT ) :: struc
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE ( potentialparm_type ), DIMENSION ( :, : ), INTENT ( IN ) :: potparm
  
!------------------------------------------------------------------------------
  
  SELECT CASE ( simpar % program )
  CASE DEFAULT
     CALL stop_prg ( 'force', 'not implemented' )
     
  CASE ( 'FIST' )
     CALL fist_force_control ( struc % molecule, struc % pnode, struc % part, &
          struc % box, thermo, potparm, ewald_param, simpar % ensemble, &
          fistpar )
     
  CASE ( 'TBMD' )
     CALL tbmd_force_control ( struc % molecule, struc % pnode, struc % part, &
          struc % box, thermo, potparm, ewald_param, simpar % ensemble, &
          tbmdpar )
     
  END SELECT
  
END SUBROUTINE force

!******************************************************************************

END MODULE force_control
