!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/linklist_types [1.0] *
!!
!!   NAME
!!     linklist_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     HAF
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE linklist_types

  USE linklist_cell_types, ONLY : cell_neighborlist_type, &
       cell_atoms_linklist_type, &
       neighbor_cells_type
  USE molecule_types, ONLY : linklist_neighbor
  USE linklist_utilities, ONLY : startlist_images_type
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: linklist_internal_data_type

  TYPE linklist_internal_data_type
     TYPE ( cell_neighborlist_type ), DIMENSION (:,:,:), POINTER :: pp_cell_ll
     TYPE ( cell_atoms_linklist_type ), DIMENSION (:,:,:,:), POINTER  :: &
         pp_startlist_cell
     TYPE ( neighbor_cells_type ), DIMENSION ( :, : ), POINTER :: pp_neighbor
     TYPE ( linklist_neighbor ), POINTER :: pp_startlist_nl
     TYPE ( startlist_images_type ), POINTER :: pp_startlist_im ( : )
     INTEGER, DIMENSION ( 3 ) :: pp_cell
     INTEGER, POINTER :: pp_images(:,:,:)
     INTEGER, POINTER :: pp_ncell ( : )
     TYPE ( global_environment_type ) :: globenv_ll
     INTEGER :: natom_types
     INTEGER :: list_type
     INTEGER :: counter, last_update, num_update
     INTEGER :: print_level
     REAL ( dbl ), DIMENSION ( :, : ), POINTER :: rlist_cut
     REAL ( dbl ), DIMENSION ( :, : ), POINTER :: rlist_cutsq
     REAL ( dbl ) :: verlet_skin
     REAL ( dbl ) :: lup, aup
     REAL ( dbl ), DIMENSION ( :, : ), POINTER :: r_last_update
  END TYPE linklist_internal_data_type

END MODULE linklist_types

!******************************************************************************
!!*****
