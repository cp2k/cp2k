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
  use coefficient_types,        only: coeff_array_type,coeff_type,&
       coeff_deallocate
  use cube_utils,               only: cube_info_type
  use global_types, only: global_environment_type,DEBUG
  use qs_collocate_density,     only: calculate_rho_elec
  use sparse_matrix_types,      only: real_matrix_type,&
       real_matrix_p_type
  use cp_block_matrix,   only: cp_block_matrix_p_type, &
       cp_block_matrix_type,cp_sp_array_dealloc
  use l_utils,                  only: l_info_type
  use gaussian_gridlevels,      only: gridlevel_info_type
  use qs_environment_types, only: qs_environment_type, get_qs_env
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='qs_rho_types'
  integer, private, save :: last_rho_id_nr=0

  public :: qs_rho_type
  public :: qs_rho_retain, qs_rho_release, qs_rho_create, qs_rho_did_change,&
       qs_rho_get, qs_rho_update_rho_r_g_new_ao

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
     type(cp_block_matrix_p_type), dimension(:), pointer :: rho_ao
     type(coeff_array_type), dimension(:), pointer :: rho_g, rho_r
     logical :: rho_ao_valid, rho_g_valid, rho_r_valid
     integer :: ref_count, id_nr, updates_sice_rebuild, rebuild_each
  end type qs_rho_type
  !!***
  !****************************************************************************

contains

!!****f* qs_rho_types/qs_rho_create [1.0] *
!!
!!   NAME
!!     qs_rho_create
!!
!!   FUNCTION
!!     allocates and initializes a rho type
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - rho_struct: the type to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!     - *_ptr: will be deallocated by this type
!!     - rho_ao_valid: defaults to true, the others to false
!!     for the explanation see the attributes of qs_rho_type
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
subroutine qs_rho_create(rho_struct,rho_ao_ptr,rho_r_ptr, rho_g_ptr,&
     rho_ao_valid,rho_r_valid,rho_g_valid, rebuild_each,error)
  type(qs_rho_type), pointer :: rho_struct
  type(cp_block_matrix_p_type), dimension(:), pointer :: rho_ao_ptr
  type(coeff_array_type), dimension(:), pointer :: rho_g_ptr,&
       rho_r_ptr
  logical, intent(in),optional :: rho_ao_valid, rho_g_valid,&
       rho_r_valid  
  type(cp_error_type), optional, intent(inout) :: error
  integer, intent(in), optional :: rebuild_each
  
  logical :: failure
  integer :: stat
  character(len=*), parameter :: routineN='qs_rho_create',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(associated(rho_ao_ptr),cp_failure_level,routineP,error,failure)
  CPPrecondition(associated(rho_r_ptr),cp_failure_level,routineP,error,failure)
  CPPrecondition(associated(rho_g_ptr),cp_failure_level,routineP,error,failure)
  if (.not.failure) then
     allocate(rho_struct,stat=stat)
     CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
  end if
  if (.not. failure) then
     last_rho_id_nr=last_rho_id_nr+1
     rho_struct%id_nr=last_rho_id_nr
     rho_struct%ref_count=1
     rho_struct%rebuild_each=5
     rho_struct%rho_ao => rho_ao_ptr
     rho_struct%rho_r => rho_r_ptr
     rho_struct%rho_g => rho_g_ptr
     rho_struct%rho_ao_valid = .true.
     rho_struct%rho_r_valid = .false.
     rho_struct%rho_g_valid = .false.
     if (present(rebuild_each)) rho_struct%rebuild_each=rebuild_each
     if (present(rho_ao_valid)) rho_struct%rho_ao_valid=rho_ao_valid
     if (present(rho_r_valid)) rho_struct%rho_r_valid=rho_r_valid
     if (present(rho_g_valid)) rho_struct%rho_g_valid=rho_g_valid
     call cp_assert(rho_ao_valid.or.rho_r_valid.or.rho_g_valid,&
          cp_failure_level,cp_assertion_failed,fromWhere=routineP,&
          message="at least one representation must be valid "//&
          CPSourceFileRef,&
          error=error,failure=failure)
  end if
end subroutine qs_rho_create
!***************************************************************************

!!****f* qs_rho_types/qs_rho_retain [1.0] *
!!
!!   NAME
!!     qs_rho_retain
!!
!!   FUNCTION
!!     retains a rho_struct by increasing the reference count by one
!!     (to be called when you what to keep a shared copy)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - rho_struct: the structure to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
subroutine qs_rho_retain(rho_struct,error)
  type(qs_rho_type), pointer :: rho_struct
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='qs_rho_retain',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(associated(rho_struct),cp_failure_level,routineP,error,failure)
  if (.not.failure) then
     CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
     rho_struct%ref_count=rho_struct%ref_count+1
  end if
end subroutine qs_rho_retain
!***************************************************************************

!!****f* qs_rho_types/qs_rho_release [1.0] *
!!
!!   NAME
!!     qs_rho_release
!!
!!   FUNCTION
!!     releases a rho_struct by decreasing the reference count by one
!!     and deallocating if it reaces 0 (to be called when you don't whant
!!     anymore a shared copy)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - rho_struct: the structure to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
subroutine qs_rho_release(rho_struct,error)
  type(qs_rho_type), pointer :: rho_struct
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='qs_rho_release',&
        routineP=moduleN//':'//routineN
  integer :: stat, i, j
  failure=.false.
  
  CPPrecondition(associated(rho_struct),cp_failure_level,routineP,error,failure)
  if (.not.failure) then
     CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
     rho_struct%ref_count=rho_struct%ref_count-1
     if (rho_struct%ref_count<1) then
        call cp_sp_array_dealloc(rho_struct%rho_ao,error=error)
        if (associated(rho_struct%rho_r)) then
           do i=1,size(rho_struct%rho_r)
              do j=1,size(rho_struct%rho_r(i)%coeff_array)
                 call coeff_deallocate(rho_struct%rho_r(i)%coeff_array(j))
              end do
              deallocate(rho_struct%rho_r(i)%coeff_array,stat=stat)
              CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
           end do
           deallocate(rho_struct%rho_r,stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        end if
        if (associated(rho_struct%rho_g)) then
           do i=1,size(rho_struct%rho_g)
              do j=1,size(rho_struct%rho_g(i)%coeff_array)
                 call coeff_deallocate(rho_struct%rho_g(i)%coeff_array(j))
              end do
              deallocate(rho_struct%rho_g(i)%coeff_array,stat=stat)
              CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
           end do
           deallocate(rho_struct%rho_g,stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        end if
        deallocate(rho_struct,stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
     end if
  end if
  nullify(rho_struct)
end subroutine qs_rho_release
!***************************************************************************

!!****f* qs_rho_types/qs_rho_did_change [1.0] *
!!
!!   NAME
!!     qs_rho_did_change
!!
!!   FUNCTION
!!     this function should be called wehn something changes to update the
!!     internal status 
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - rho_struct: the density structure that changed
!!     - rho_ao_valid,rho_r_valid, rho_g_valid: actual valid representations
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
subroutine qs_rho_did_change(rho_struct,rho_ao_valid,rho_r_valid, &
     rho_g_valid, error)
  type(qs_rho_type), pointer :: rho_struct
  logical, intent(in), optional :: rho_ao_valid, rho_r_valid, rho_g_valid
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='qs_rho_did_change',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(associated(rho_struct),cp_failure_level,routineP,error,failure)
  CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     if (present(rho_ao_valid)) rho_struct%rho_ao_valid=rho_ao_valid
     if (present(rho_r_valid)) rho_struct%rho_r_valid=rho_r_valid
     if (present(rho_g_valid)) rho_struct%rho_g_valid=rho_g_valid
     call cp_assert(rho_ao_valid.or.rho_r_valid.or.rho_g_valid,&
          cp_failure_level,cp_assertion_failed,fromWhere=routineP,&
          message="at least one representation must be valid "//&
          CPSourceFileRef,&
          error=error,failure=failure)
  end if
end subroutine qs_rho_did_change
!***************************************************************************

!!****f* qs_rho_types/qs_rho_get [1.0] *
!!
!!   NAME
!!     qs_rho_get
!!
!!   FUNCTION
!!     returns info about the density described by this object.
!!     If some representation is not available an error is issued
!!
!!   NOTES
!!     as the construction is costly a get method that gives an error
!!     if not valid (like this at the moment) can be useful
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
subroutine qs_rho_get(rho_struct, rho_ao, rho_r, rho_g, rho_ao_valid,&
     rho_r_valid, rho_g_valid, rebuild_each, error)
  type(qs_rho_type), pointer :: rho_struct
  type(cp_block_matrix_p_type), dimension(:), pointer,optional :: rho_ao
  type(coeff_array_type), dimension(:), pointer,optional :: rho_g, rho_r
  logical, intent(out), optional :: rho_ao_valid, rho_g_valid, rho_r_valid
  integer, intent(out), optional :: rebuild_each
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='qs_rho_get',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(associated(rho_struct),cp_failure_level,routineP,error,failure)
  CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     if (present(rho_ao)) then
        if (.not.rho_struct%rho_ao_valid) then
           call cp_unimplemented_error(routineP,"to do",error=error)
           nullify(rho_ao)
        else
           rho_ao => rho_struct%rho_ao
        end if
     end if
     if (present(rho_r)) then
        if (.not.rho_struct%rho_r_valid) then
           call cp_unimplemented_error(routineP,"to do",error=error)
           nullify(rho_r)
        else
           rho_r => rho_struct%rho_r
        end if
     end if
     if (present(rho_g)) then
        if (.not.rho_struct%rho_g_valid) then
           call cp_unimplemented_error(routineP,"to do",error=error)
           nullify(rho_g)
        else
           rho_g => rho_struct%rho_g
        end if
     end if
     if (present(rho_ao_valid)) rho_ao_valid=rho_struct%rho_ao_valid
     if (present(rho_r_valid)) rho_r_valid=rho_struct%rho_r_valid
     if (present(rho_g_valid)) rho_g_valid=rho_struct%rho_g_valid
     if (present(rebuild_each)) rebuild_each=rho_struct%rebuild_each
  end if
end subroutine qs_rho_get
!***************************************************************************

!!****f* qs_rho_types/qs_rho_update_rho_r_g_new_ao [1.0] *
!!
!!   NAME
!!     qs_rho_update_rho_r_g_new_ao
!!
!!   FUNCTION
!!     updates rho_r and rho_g to the new (given) ao. Does NOT modify
!!     rho_ao
!!
!!   NOTES
!!     ugly interface, too many arguments
!!
!!   INPUTS
!!     - rho_struct: the rho structure that should be changed
!!     - rho_ao_new: the new value for rho_ao
!!     - total_rho_rspace: returns the integrated charge in r space
!!     - rebuild: if false builds rho_r by just updating the parts
!!       of rho that changed with respect to the old density
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!     - tmp_rspace: is a temporary grid like rho_struct%rho_r(1)coeff_array
!!       at the moment on output it will contain exactly the same as
!!       it, but try not rely on it
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
subroutine qs_rho_update_rho_r_g_new_ao(rho_struct, rho_ao_new,gridlevel_info,&
     cube_info,l_info,tmp_rspace,qs_env,global_env,total_rho_rspace,&
     rebuild,error)
  type(qs_rho_type), pointer :: rho_struct
  type(cube_info_type), dimension(:), pointer :: cube_info
  type(global_environment_type), intent(in) :: global_env
  type(gridlevel_info_type), intent(in) :: gridlevel_info
  type(l_info_type), intent(in) :: l_info
  type(qs_environment_type), intent(in) :: qs_env
  type(cp_block_matrix_p_type),dimension(:), intent(in) :: rho_ao_new
  real(wp), intent(out),optional :: total_rho_rspace
  logical, intent(in),optional :: rebuild
  type(coeff_type), dimension(:), pointer :: tmp_rspace
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure,r_build
  character(len=*), parameter :: routineN='qs_rho_update_rho_r_g_new_ao',&
        routineP=moduleN//':'//routineN
  real(kind=wp) :: t_rho_r
  failure=.false.
  
  CPPrecondition(associated(rho_struct),cp_failure_level,routineP,error,failure)
  CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     r_build=rho_struct%updates_sice_rebuild > rho_struct%rebuild_each
     if (present(rebuild)) r_build=rebuild
     if (.not.(rho_struct%rho_ao_valid.and.rho_struct%rho_r_valid)) then
        r_build=.true.
     end if
     if (r_build) then
        rho_struct%updates_sice_rebuild=0
     else
        rho_struct%updates_sice_rebuild=rho_struct%updates_sice_rebuild+1
     end if
     call calculate_rho_elec(rho_ao_new(1)%matrix%sm,&
          rho_struct%rho_ao(1)%matrix%sm,r_build,tmp_rspace,&
          rho_struct%rho_r(1)%coeff_array,rho_struct%rho_g(1)%coeff_array,&
          gridlevel_info,cube_info,l_info,t_rho_r,&
          qs_env,global_env)
     rho_struct%rho_ao_valid=.false.
     rho_struct%rho_r_valid=.true.
     rho_struct%rho_g_valid=.true.
     if (present(total_rho_rspace)) total_rho_rspace=t_rho_r
  end if
end subroutine qs_rho_update_rho_r_g_new_ao
!***************************************************************************

end module qs_rho_types

