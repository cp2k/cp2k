!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/particle_types [1.0] *
!!
!!   NAME
!!     particle_types
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

MODULE particle_types

  USE kinds, ONLY : dbl
  USE linklist_cell_types, ONLY : cell_neighborlist_type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: particle_type, particle_prop_type

  TYPE particle_prop_type
     REAL ( dbl ) :: charge
     REAL ( dbl ) :: mass
     INTEGER    :: ptype
  END TYPE particle_prop_type

  TYPE particle_type
     INTEGER :: iatom
     INTEGER ::  grid_id 
     INTEGER ::  node 
     REAL ( dbl ), DIMENSION (3) :: r
     REAL ( dbl ), DIMENSION (3) :: f
     REAL ( dbl ), DIMENSION (3) :: v
     TYPE (particle_prop_type ), POINTER :: prop
     TYPE ( cell_neighborlist_type ), POINTER :: cell
     TYPE ( cell_neighborlist_type ), POINTER :: cell_ol
     INTEGER, DIMENSION(:), POINTER :: coef_list 
  END TYPE particle_type

END MODULE particle_types

!******************************************************************************
!!*****
