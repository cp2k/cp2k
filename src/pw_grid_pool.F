!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pw_grid_pool [1.0] *
!!
!!   NAME
!!     pw_grid_pool
!!
!!   FUNCTION
!!     Manages a pool of grids (to be used for example as tmp objects),
!!     but can also be used to intantiate grids that are never given back.
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE pw_grid_pool
  USE cp_log_handling, ONLY: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log
  USE cp_error_handling, ONLY: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger
  USE kinds, ONLY: wp=>dp
  USE cp_para_types, ONLY: cp_para_env_type
  USE timings, ONLY: timeset, timestop
  USE cp_linked_list_pw, ONLY: cp_sll_pw_type
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='pw_grid_pool'

!***
!****************************************************************************

  !!****s* pw_grid_pool/pw_grid_pool_type [1.0] *
  !!
  !!   NAME
  !!    pw_grid_pool_type
  !!
  !!   FUNCTION
  !!     Manages a pool of grids (to be used for example as tmp objects),
  !!     but can also be used to intantiate grids that are never given back.
  !!     
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     -
  !!
  !!   AUTHOR
  !!     Fawzi Mohamed
  !!
  !!   MODIFICATION HISTORY
  !!     08.2002 created [fawzi]
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE pw_grid_pool_type
     INTEGER :: ref_count, id_nr
     TYPE(cp_sll_pw_type), POINTER :: real1d_grids, real3d_grids,&
          complex1d_grids, complex3d_grids
  END TYPE pw_grid_pool_type
  !!***
  !****************************************************************************

CONTAINS


END MODULE pw_grid_pool
