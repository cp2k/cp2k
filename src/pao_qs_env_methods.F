!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_qs_env_methods [1.0] *
!!
!!   NAME
!!     pao_qs_env_methods
!!
!!   FUNCTION
!!     pao related methods that work on the qs_env methods
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
module pao_qs_env_methods
  use cp_log_handling, only: cp_to_string, cp_failure_level, cp_warning_level,&
       cp_note_level
  use cp_error_handling, only: cp_debug, cp_error_type, cp_assert, &
       cp_error_message
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use pao_types, only: pao_env_type, pao_env_get, pao_glob_angles_type
  use pao_bas_projection_types, only: pao_bas_projection_type
  use qs_environment_types, only: qs_environment_type, get_qs_env
  use atomic_kind_types, only: atomic_kind_type, get_atomic_kind,&
       get_atomic_kind_set
  use pao_glob_angles_methods, only: pao_g_ang_did_change
  use cp_sparse_matrix, only: cp_sparse_matrix_p_type, cp_sparse_matrix_type,&
       cp_sp_dealloc, cp_sp_create, cp_sp_release
  use cp_b_matrix_structure, only: cp_b_matrix_struct_type, cp_bmstruct_get,&
       cp_bmstruct_create, cp_bmstruct_release, cp_b_dims_type
  use cp_matrix_utils, only: cp_sm_set_dist_sparsity_from_sm
  use sparse_matrix_types, only: real_matrix_p_type, real_matrix_type
  use particle_types, only: particle_type
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='pao_qs_env_methods'

  public :: cp_get_b_projection, cp_did_change
  ! underlying routines
  public :: pao_qs_env_get_b_projection, pao_qs_env_get, pao_qs_env_did_change
  public :: pao_bas_projection_type ! for the intel compiler
  interface cp_get_b_projection
     module procedure pao_qs_env_get_b_projection
  end interface
  interface cp_did_change
     module procedure pao_qs_env_did_change
  end interface
!***
!****************************************************************************
contains

!!****f* pao_qs_env_methods/pao_qs_env_get_b_projection [1.0] *
!!
!!   NAME
!!     pao_qs_env_get_b_projection
!!
!!   SYNOPSIS
!!     Function pao_qs_env_get_b_projection(qs_env, atom_nr, error)&
!!         Result(res)
!!       Type(pao_bas_projection_type), Pointer:: res
!!       Type(qs_environment_type), Intent (IN):: qs_env
!!       Integer, Intent (IN):: atom_nr
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_qs_env_get_b_projection
!!
!!   FUNCTION
!!     returns the base projection projection on atom_nr
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     qs_env: the qs environement you get the data from
!!     atom_nr: the number of the atom
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     - 1.2002 created [fawzi]
!!     - 19.04.02 used new projection
!!     - 23.04.02 moved to qs
!!
!!*** **********************************************************************
  function pao_qs_env_get_b_projection(qs_env, atom_nr, error) result(res)
    type(pao_bas_projection_type), pointer ::res
    type(qs_environment_type), intent(in) :: qs_env
    integer, intent(in) :: atom_nr
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='pao_qs_env_get_b_projection',&
        routineP=moduleN//':'//routineN
    type(particle_type), dimension(:), pointer :: particle_set
    failure=.false.

    call get_qs_env(qs_env, particle_set=particle_set)
    call get_atomic_kind(particle_set(atom_nr)%atomic_kind,pao_bas_proj=res)
  end function pao_qs_env_get_b_projection
!***************************************************************************

!!****f* pao_qs_env_methods/pao_qs_env_get [1.0] *
!!
!!   NAME
!!     pao_qs_env_get
!!
!!   SYNOPSIS
!!     Subroutine pao_qs_env_get(qs_env, global_env, pao_env,&
!!         pao_max_min_bas, pao_max_full_bas, pao_max_excl_bas,&
!!         pao_max_l_angles, pao_tot_min_bas, pao_tot_full_bas,&
!!         pao_tot_excl_bas, pao_min_density_m, pao_min_bas_m_struct,&
!!         pao_full_bas_m_struct, pao_angles_att,&
!!         pao_w_subspace_pos_precision, pao_w_in_subspace_precision,&
!!         error)
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(pao_env_type), Pointer, Optional:: pao_env
!!       Integer, Intent (OUT), Optional:: pao_max_min_bas,&
!!         pao_max_full_bas, pao_max_excl_bas, pao_max_l_angles,&
!!         pao_tot_min_bas, pao_tot_full_bas, pao_tot_excl_bas
!!       Real(Kind=wp), Optional, Intent (OUT)::&
!!         pao_w_subspace_pos_precision
!!       Real(Kind=wp), Optional, Intent (OUT)::&
!!         pao_w_in_subspace_precision
!!       Type(pao_glob_angles_type), Optional, Pointer:: pao_angles_att
!!       Type(cp_sparse_matrix_p_type), Dimension(:), Optional, Pointer::&
!!         pao_min_density_m
!!       Type(cp_b_matrix_struct_type), Pointer, Optional::&
!!         pao_min_bas_m_struct, pao_full_bas_m_struct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_qs_env_get
!!
!!   FUNCTION
!!     returns pao related info from the qs_environment_type
!!
!!   NOTES
!!     use this function rather than get_atomic_kind_set: some values
!!     might be cached in pao_env in the future
!!
!!   INPUTS
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
subroutine pao_qs_env_get(qs_env, global_env, pao_env, pao_max_min_bas,&
     pao_max_full_bas, pao_max_excl_bas, pao_max_l_angles,&
     pao_tot_min_bas, pao_tot_full_bas, pao_tot_excl_bas,&
     pao_min_density_m, pao_min_bas_m_struct, pao_full_bas_m_struct,&
     pao_angles_att, pao_w_subspace_pos_precision,&
     pao_w_in_subspace_precision, error)
  type(qs_environment_type), intent(in), target :: qs_env
  type(global_environment_type), intent(in), target :: global_env
  type(pao_env_type), pointer, optional :: pao_env
  integer, intent(out), optional :: pao_max_min_bas,&
     pao_max_full_bas, pao_max_excl_bas, pao_max_l_angles,&
     pao_tot_min_bas, pao_tot_full_bas, pao_tot_excl_bas
    real(kind=wp), optional, intent(out) :: pao_w_subspace_pos_precision
    real(kind=wp), optional, intent(out) :: pao_w_in_subspace_precision
    type(pao_glob_angles_type), optional, pointer :: pao_angles_att
    type(cp_sparse_matrix_p_type), dimension(:), optional, pointer ::&
         pao_min_density_m
    type(cp_b_matrix_struct_type), pointer, optional :: pao_min_bas_m_struct, &
         pao_full_bas_m_struct
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='pao_qs_env_get',&
        routineP=moduleN//':'//routineN
  type(atomic_kind_type), dimension(:), pointer :: atomic_kind_set
  type(pao_env_type), pointer :: p_env
  failure=.false.
  
  call get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, pao_env=p_env)
  CPPrecondition(associated(atomic_kind_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(associated(p_env),cp_failure_level,routineP,error,failure)
  call get_atomic_kind_set(atomic_kind_set,pao_max_min_bas=pao_max_min_bas,&
       pao_max_full_bas=pao_max_full_bas,&
       pao_max_excl_bas=pao_max_excl_bas, pao_tot_min_bas=pao_tot_min_bas,&
       pao_tot_full_bas=pao_tot_full_bas, pao_tot_excl_bas=pao_tot_excl_bas,&
       pao_max_l_angles=pao_max_l_angles)
  call pao_env_get(p_env,&
       wanted_subspace_pos_precision=pao_w_subspace_pos_precision,&
       wanted_in_subspace_precision=pao_w_in_subspace_precision,&
       angles_att=pao_angles_att,min_density_m=pao_min_density_m,&
       min_bas_m_struct=pao_min_bas_m_struct,&
       full_bas_m_struct=pao_full_bas_m_struct, error=error)
  if (present(pao_env)) pao_env => p_env
end subroutine pao_qs_env_get
!***************************************************************************

!!****f* pao_qs_env_methods/pao_qs_env_did_change [1.0] *
!!
!!   NAME
!!     pao_qs_env_did_change
!!
!!   SYNOPSIS
!!     Subroutine pao_qs_env_did_change(qs_env, global_env, h_changed,&
!!         s_changed, angles_changed, full_reset, error)
!!       Type(qs_environment_type), Intent (INOUT):: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Logical, Intent (IN), Optional:: h_changed, s_changed,&
!!         angles_changed, full_reset
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_qs_env_did_change
!!
!!   FUNCTION
!!     invalidates the cache when some data changes
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env: the environement that has changed
!!     - s_changed: if the overlap has changed (implies h_changed
!!     - angles_changed: if the angles that describe tha actual subspace 
!!       have changed
!!     - full_reset: if the structures should be regenerated from scratch
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
subroutine pao_qs_env_did_change(qs_env,global_env,&
     h_changed, s_changed,&
     angles_changed, full_reset, error)
  type (qs_environment_type), intent(inout) :: qs_env
  type(global_environment_type), intent(in), target :: global_env
  logical, intent(in), optional :: h_changed, s_changed, &
       angles_changed, full_reset
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='pao_qs_env_did_change',&
       routineP=moduleN//':'//routineN
  logical :: s_c, full_r
  integer :: i
  type(pao_env_type), pointer :: pao_env
  type(real_matrix_p_type), pointer :: overlap_rmp
  type(cp_b_dims_type), pointer :: min_dims, full_dims
  type(cp_b_matrix_struct_type), pointer :: old_min_struct, old_full_struct
  failure=.false.
  s_c=.false.; full_r=.false.
  
  call get_qs_env(qs_env,pao_env=pao_env)
  CPPrecondition(associated(pao_env),cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     CPPrecondition(pao_env%initialized,cp_failure_level,routineP,error,failure)
  end if
  if (.not. failure) then
     if (present(s_changed)) s_c=s_changed
     if (present(full_reset)) full_r=full_reset
     if (s_c.or.full_r) then
!MK        call get_qs_env(qs_env, s=overlap_rmp)
        ! invalidate min_density_m
        pao_env%min_d_m_valid=.false.
        do i=1,size(pao_env%min_density_m)
           call cp_sp_release(pao_env%min_density_m(i)%matrix,error=error)
        end do
        ! invalidate matrix structure
        call cp_bmstruct_get(pao_env%min_bas_m_struct,&
             global_env=global_env, b_dims=min_dims, error=error)
        call cp_bmstruct_get(pao_env%full_bas_m_struct,&
             global_env=global_env, b_dims=full_dims, error=error)
        old_min_struct => pao_env%min_bas_m_struct
        old_full_struct => pao_env%full_bas_m_struct
        call cp_bmstruct_create(pao_env%min_bas_m_struct,&
               symmetric=.true.,&
               global_env=global_env,&
               b_dims=min_dims,&
               error=error)
        call cp_sm_set_dist_sparsity_from_sm(pao_env%min_bas_m_struct,&
             global_env=global_env,sample_matrix=overlap_rmp%matrix,&
             error=error)
        call cp_bmstruct_release(old_min_struct,error=error)
        call cp_bmstruct_create(pao_env%full_bas_m_struct,&
             symmetric=.true.,&
             global_env=global_env,&
             b_dims=full_dims,&
             d_struct=pao_env%min_bas_m_struct%distribution,&
             sparsity=pao_env%min_bas_m_struct%sparsity,&
             error=error)
        call cp_bmstruct_release(old_full_struct,error=error)
        
        do i=1,size(pao_env%min_density_m)
           call cp_sp_create(pao_env%min_density_m(i)%matrix,&
                global_env=global_env,&
                matrix_structure=pao_env%min_bas_m_struct,error=error)
        end do
     end if
     call pao_g_ang_did_change(pao_env%angles_att,qs_env=qs_env,&
          global_env=global_env, h_changed=h_changed,&
          s_changed=s_changed, angles_changed=angles_changed,&
          full_reset=full_reset, error=error)
  end if
end subroutine pao_qs_env_did_change
!***************************************************************************

end module pao_qs_env_methods
