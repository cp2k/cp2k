!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_rho_types [1.0] *
!!
!!   NAME
!!     qs_rho_types
!!
!!   FUNCTION
!!     superstucture that hold various representations of the density and 
!!     keeps track of which ones are valid
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
module qs_rho_types
  use cp_log_handling, only: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log
  use cp_error_handling, only: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger
  use kinds, only: wp=>dp
  use cp_para_types, only: cp_para_env_type
  use timings, only: timeset, timestop
  use coefficient_types,        only: coeff_type, coeff_deallocate
  use cube_utils,               only: cube_info_type
  use global_types, only: global_environment_type,DEBUG
  use l_utils,                  only: l_info_type
  use gaussian_gridlevels,      only: gridlevel_info_type
  use sparse_matrix_types,      only: real_matrix_p_type
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='qs_rho_types'
  integer, private, save :: last_rho_id_nr=0

  public :: qs_rho_type

!***
!****************************************************************************

  !!****s* qs_rho_types/qs_rho_type [1.0] *
  !!
  !!   NAME
  !!     qs_rho_type
  !!
  !!   FUNCTION
  !!     keeps the density in various representations, keeping track of
  !!     which ones are valid.
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     most attributes are array with either lda or lsd_alpha,lsd_beta.
  !!     - rho_ao: the filterd rho in the localized atom basis (to have rho(r)
  !!       the filtred matrix is enough, but rho(r,r') is lost).
  !!     - rho_r: grids with rho in the real space
  !!     - rho_g: grids with rho in the g space
  !!     - rho_ao_valid, rho_g_valid, rho_r_valid: if the corresponding
  !!       component is valid
  !!     - ref_count: the reference count, when it becomes 0 the type
  !!       is deallocated.
  !!
  !!   AUTHOR
  !!     Fawzi Mohamed
  !!
  !!   MODIFICATION HISTORY
  !!     08.2002 created [fawzi]
  !!
  !!   SOURCE
  !***************************************************************************
  type qs_rho_type
     type(real_matrix_p_type), dimension(:), pointer :: rho_ao
     type(coeff_type), dimension(:), pointer :: rho_g, rho_r
     logical :: rho_ao_valid, &
                rho_g_valid, &
                rho_r_valid
     integer :: ref_count, &
                id_nr, &
                updates_sice_rebuild, &
                rebuild_each
  end type qs_rho_type
  !!***
  !****************************************************************************

end module qs_rho_types

