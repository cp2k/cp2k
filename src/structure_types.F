!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/structure_types [1.0] *
!!
!!   NAME
!!     structure_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     Harald Forbert (Dec-2000): Add ll-data to structure_type
!!     CJM 14-April-2001 : added coeff_type to structure_type
!!
!!   SOURCE
!******************************************************************************

MODULE structure_types

  USE kinds, ONLY : dbl
  USE atomic_kinds, ONLY :  kind_info_type
  USE coefficient_types, ONLY : coeff_type
  USE empirical_parameters, ONLY : empirical_parameter_type
  USE linklist_types, ONLY : linklist_internal_data_type
  USE molecule_types, ONLY : molecule_structure_type, particle_node_type
  USE pair_potential, ONLY : potentialparm_type
  USE particle_types, ONLY : particle_type
  USE pw_grid_types, ONLY : pw_grid_type
  USE simulation_cell, ONLY : cell_type
  USE tbmd_types, ONLY : tbatom_param_type, tb_hopping_type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: structure_type, interaction_type

  PUBLIC :: cell_type, particle_type, particle_node_type
  PUBLIC :: molecule_structure_type
  PUBLIC :: pimd_data_type
  PUBLIC :: potentialparm_type, tbatom_param_type, tb_hopping_type

  TYPE pimd_data_type
     REAL ( dbl ), DIMENSION ( : , : ), POINTER :: u
     REAL ( dbl ) :: lbd
     REAL ( dbl ) :: slbd
     REAL ( dbl ) :: meff
     REAL ( dbl ) :: ismeff
     REAL ( dbl ), DIMENSION ( 2 , 2 ) :: harmonic_prop
  END TYPE pimd_data_type

  TYPE structure_type
     CHARACTER ( LEN = 30 ) :: name
     TYPE ( cell_type ) :: box
     TYPE ( cell_type ) :: box_ref
     TYPE ( particle_type ), POINTER :: part ( : )
     TYPE ( particle_node_type ), POINTER :: pnode ( : )
     TYPE ( molecule_structure_type ), POINTER :: molecule ( : )
     TYPE ( linklist_internal_data_type), POINTER :: ll_data ( : )
     TYPE ( pimd_data_type ) :: pimd_data
     TYPE ( kind_info_type ), POINTER :: rho0_basis_info ( : )
     TYPE ( kind_info_type ), POINTER :: drho_basis_info ( : )
     TYPE ( coeff_type ), POINTER, DIMENSION ( : )  :: coef_pos         
     TYPE ( coeff_type ), POINTER, DIMENSION ( : )  :: coef_vel         
     TYPE ( coeff_type ), POINTER, DIMENSION ( : )  :: coef_force         
     TYPE ( coeff_type ), POINTER, DIMENSION ( : )  :: projector
     TYPE ( pw_grid_type ), POINTER, DIMENSION ( : )  :: pw_grid
  END TYPE structure_type

  TYPE interaction_type
     TYPE ( potentialparm_type ), DIMENSION ( :, : ), POINTER :: potparm
     TYPE ( tbatom_param_type ), DIMENSION ( : ), POINTER :: tbatom
     TYPE ( tb_hopping_type ), DIMENSION ( :, : ), POINTER :: tbhop
     TYPE ( empirical_parameter_type ), DIMENSION (:), POINTER :: empparm
  END TYPE interaction_type

END MODULE structure_types

!!*****
!******************************************************************************
