!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE structure_types
  
  USE kinds, ONLY : dbl
  USE molecule_types, ONLY : molecule_structure_type, particle_node_type
  USE particle_types, ONLY : particle_type
  USE simulation_cell, ONLY : cell_type
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: structure_type
  PUBLIC :: cell_type, particle_type, particle_node_type
  PUBLIC :: molecule_structure_type
  
  TYPE structure_type
     CHARACTER ( LEN = 30 ) :: name
     TYPE (cell_type ) :: box
     TYPE (particle_type ), POINTER :: part ( : )
     TYPE (particle_node_type ), POINTER :: pnode ( : )
     TYPE (molecule_structure_type ), POINTER :: molecule ( : )
  END TYPE structure_type
  
END MODULE structure_types
