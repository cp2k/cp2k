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
!!
!!   SOURCE
!******************************************************************************

MODULE structure_types
  
  USE kinds, ONLY : dbl
  USE molecule_types, ONLY : molecule_structure_type, particle_node_type
  USE particle_types, ONLY : particle_type
  USE simulation_cell, ONLY : cell_type
  USE pair_potential, ONLY : potentialparm_type
  USE tbmd_types, ONLY : tbatom_param_type, tb_hopping_type
  USE linklist_types, ONLY : linklist_internal_data_type
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: structure_type, interaction_type

  PUBLIC :: cell_type, particle_type, particle_node_type
  PUBLIC :: molecule_structure_type
  PUBLIC :: potentialparm_type, tbatom_param_type, tb_hopping_type
  
  TYPE structure_type
     CHARACTER ( LEN = 30 ) :: name
     TYPE ( cell_type ) :: box
     TYPE ( particle_type ), POINTER :: part ( : )
     TYPE ( particle_node_type ), POINTER :: pnode ( : )
     TYPE ( molecule_structure_type ), POINTER :: molecule ( : )
     TYPE ( linklist_internal_data_type) :: ll_data
  END TYPE structure_type
  
  TYPE interaction_type
     TYPE ( potentialparm_type ), DIMENSION ( :, : ), POINTER :: potparm
     TYPE ( tbatom_param_type ), DIMENSION ( : ), POINTER :: tbatom
     TYPE ( tb_hopping_type ), DIMENSION ( :, : ), POINTER :: tbhop
  END TYPE interaction_type

END MODULE structure_types

!!*****
!******************************************************************************
