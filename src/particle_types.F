!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

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
     REAL ( dbl ), DIMENSION (3) :: r
     REAL ( dbl ), DIMENSION (3) :: f
     REAL ( dbl ), DIMENSION (3) :: v
     TYPE (particle_prop_type ), POINTER :: prop
     TYPE ( cell_neighborlist_type ), POINTER :: cell
  END TYPE particle_type
  
END MODULE particle_types
