!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/linklist_cell_types [1.0] *
!!
!!   NAME
!!     linklist_cell_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE linklist_cell_types

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: cell_neighborlist_type, neighbor_cells_type
  PUBLIC :: cell_atoms_type, cell_atoms_linklist_type

  TYPE cell_atoms_linklist_type
     INTEGER :: atom
     LOGICAL :: used
     TYPE ( cell_atoms_linklist_type ), POINTER :: next
  END TYPE cell_atoms_linklist_type

  TYPE cell_atoms_type
     TYPE ( cell_atoms_linklist_type ), POINTER :: ll
  END TYPE cell_atoms_type

  TYPE neighbor_cells_type
     INTEGER, DIMENSION ( :, : ), POINTER :: index
  END TYPE neighbor_cells_type

  TYPE cell_neighborlist_type
     INTEGER, DIMENSION ( : ), POINTER :: natoms
     TYPE ( neighbor_cells_type ), DIMENSION ( :, : ), POINTER :: neighbor
     TYPE ( cell_atoms_type ), DIMENSION ( : ), POINTER :: ptype
  END TYPE cell_neighborlist_type

END MODULE linklist_cell_types

!******************************************************************************
!!*****
