!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_l_angles_methods [1.0] *
!!
!!   NAME
!!     pao_l_angles_methods
!!
!!   FUNCTION
!!     methods to access and modify the handle to the local angles
!!
!!   NOTES
!!     get methods with inout could be changed to in by using a ptr var 
!!     to modify the subcomponents
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
module pao_l_angles_methods
  use cp_log_handling, only: cp_to_string, cp_failure_level, cp_warning_level
  use cp_error_handling, only: cp_error_type, cp_assert, cp_error_message,&
       cp_assertion_failed, cp_debug, cp_precondition_failed, cp_error_init,&
       cp_dealloc_ref, cp_internal_error, cp_error_get_logger
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use pao_types, only: pao_env_type, pao_glob_angles_type, pao_env_get
  use pao_glob_angles_methods, only: pao_g_ang_get, cp_validate
  use sparse_matrix_types, only: real_matrix_type, replicate_matrix, &
       copy_matrix, real_matrix_p_type,&
       allocate_matrix, deallocate_matrix, get_matrix_info, get_block_node,&
       add_block_node, put_block_node
  use cp_sparse_matrix, only: cp_sparse_matrix_type, cp_sparse_matrix_p_type,&
       cp_get_matrix
  use pao_bas_projection_types, only: pao_bas_projection_type, cp_proj_get
  use pao_plain_rot, only: cp_left_compose_p_rot
  use qs_environment_types, only: qs_environment_type, get_qs_env
  use atomic_kind_types, only:atomic_kind_type, get_atomic_kind
  use particle_types, only: particle_type
  use pao_qs_env_methods, only: pao_qs_env_get
  use mathlib, only: symmetrize_matrix, power_matrix
  use cp_array_utils, only: cp_2d_r_output
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='pao_env_methods'
! types
  public :: pao_local_angles_type

! core procedures
  public :: cp_init, cp_dealloc_ref, cp_set, cp_get, &
       cp_did_change

! core functions
  public :: cp_validate, cp_valid

! special procedures
  public :: cp_set_angle

! special functions
  public :: cp_get_angle, &
       cp_get_size_min_bas, cp_get_size_full_bas, cp_get_size_excl_bas, &
       cp_get_p_indexes, cp_get_p_index, cp_get_non_p_index,&
       cp_get_non_p_indexes, cp_get_unitary_t,&
       cp_get_g_index

  public :: pao_l_compose_u ! to do remove from here and plain rot tests

! underlying functions
  public :: pao_l_ang_init, pao_l_ang_dealloc_ref,&
       pao_l_ang_validate, pao_l_ang_set,&
       pao_l_ang_get, pao_l_ang_did_change,&
       pao_l_ang_get_angle,&
       pao_l_ang_set_angle, pao_l_ang_get_size_min_bas,&
       pao_l_ang_get_size_full_bas, pao_l_ang_get_size_excl_bas,&
       pao_l_ang_get_p_indexes, pao_l_ang_get_p_index, &
       pao_l_ang_get_non_p_indexes, pao_l_ang_get_non_p_index,&
       pao_l_ang_get_unitary_t, pao_l_ang_set_unitary_t, &
       pao_l_ang_get_g_index

! ============== interfaces ===========

! +++ begin of core methods +++

!! initializes the first, argument.
!! Support optional initial values that depend on the 
!! first argument 
  interface cp_init
     module procedure pao_l_ang_init
  end interface

!! Deallocates the memory allocated by the first argument.
!! Supports error=error
  interface cp_dealloc_ref
     module procedure pao_l_ang_dealloc_ref
  end interface

!! returns true if the content of the first argumenst is valid.
!! writes the errors ar warnings, full validation (long)
  interface cp_validate
     module procedure pao_l_ang_validate
  end interface

!! returns true if the first argument is valid. Only a fast minimal 
!! validation is performed (no dangling pointers)
  interface cp_valid
     module procedure pao_l_angles_valid
  end interface

!! sets the value of various optonal attributes of the firs argument.
!! the optional attributes depend on the type of the first
  interface cp_set
     module procedure pao_l_ang_set
  end interface

!! gets values from the first argument via optional keyword arguments
!! that depend on the type of the first argument
  interface cp_get
     module procedure pao_l_ang_get
  end interface

!! tells the object that its internal values have changes and its cached
!! values may be invalid
  interface cp_did_change
     module procedure pao_l_ang_did_change
  end interface

! +++ end of the core methods +++

!! returns the angle in the first arg identified by the second (min_bas)
!! and the third arg (full_bas) 
  interface cp_get_angle
     module procedure pao_l_ang_get_angle
  end interface

!! sets the angle in the first arg identified by the second (min_bas)
!! and the third arg (full_bas) to the value in the fourth (value)
  interface cp_set_angle
     module procedure pao_l_ang_set_angle
  end interface

!! returns the size of the minimal basis
  interface cp_get_size_min_bas
     module procedure pao_l_ang_get_size_min_bas
  end interface

!! returns the size of the full basis
  interface cp_get_size_full_bas
     module procedure pao_l_ang_get_size_full_bas
  end interface

!! returns the size of the excluded basis
  interface cp_get_size_excl_bas
     module procedure pao_l_ang_get_size_excl_bas
  end interface

!! returns the projection indexes
  interface cp_get_p_indexes
     module procedure pao_l_ang_get_p_indexes
  end interface

!! returns the indexes of the excluded basis
  interface cp_get_non_p_indexes
     module procedure pao_l_ang_get_non_p_indexes
  end interface

!! returns an index of the projection
  interface cp_get_p_index
     module procedure pao_l_ang_get_p_index
  end interface

!! returns an index of the excluded basis
  interface cp_get_non_p_index
     module procedure pao_l_ang_get_non_p_index
  end interface

!! gets the unitary transformation
  interface cp_get_unitary_t
     module procedure pao_l_ang_get_unitary_t
  end interface

!! sets the unitary transformation
  interface cp_set_unitary_t
     module procedure pao_l_ang_set_unitary_t
  end interface

!! returns the global indexes
interface cp_get_g_index
   module procedure pao_l_ang_get_g_index
end interface

!!***
!****************************************************************************

!!****s* pao_types/pao_local_angles_type [1.0] *
!!
!!   NAME
!!     pao_local_angles_type
!!
!!   FUNCTION
!!     This structure hold the angles that regard one atom.
!!     it should be used to access the angles when you want an "atom view"
!!     (for example to build the unitary transformation).
!!
!!   NOTES
!!      Do not use the global structure, the idea is that parallelization
!!      (distributing various atoms to various processors) should be easier
!!      using this structure.
!!
!!   ATTRIBUTES
!!     initialized: internal flag to (maybe) find errors that come
!!          from not having initialized the structure. is true if 
!!          this type has been initialized
!!     atom_nr : the atom this object focalizes on
!!     glob_arg: a pointer to the global angles
!!     tmp_m: tmeporary array, used for internal calculations
!!      
!!     readonly (virtual) attributes:
!!
!!     n_min_bas: the size of the minimal basis on this atom.
!!     n_full_bas: the size of the full basis on this atom.
!!     pao_env: the pao environment this local angles lives in.
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  type pao_local_angles_type
     private
     logical :: initialized
     integer :: atom_nr, begins_at, size_min_bas
     real(kind=wp), dimension(:,:), pointer :: tmp_m
     type(qs_environment_type), pointer :: qs_env
     type(global_environment_type), pointer :: global_env
     type(pao_bas_projection_type), pointer :: projection
     type(pao_glob_angles_type), pointer :: glob_angles
  end type pao_local_angles_type
!!***
!****************************************************************************

contains

! ====== core methods ======

!!****f* pao_l_angles_methods/pao_l_ang_init [1.0] *
!!
!!   NAME
!!     pao_l_ang_init
!!
!!   SYNOPSIS
!!     Subroutine pao_l_ang_init(local_angles, qs_env, global_env, atom_nr,&
!!         glob_angles, error)
!!       Type(pao_local_angles_type), Intent (OUT):: local_angles
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Integer, Intent (IN):: atom_nr
!!       Type(pao_glob_angles_type), Target:: glob_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_l_ang_init
!!
!!   FUNCTION
!!     Initializes the pao angles structure
!!
!!   NOTES
!!     add source arg (i.e. copy)?
!!
!!   INPUTS
!!     local_angles: the pao_local angles to be initialized
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling
!!
!!     For the other arguments see pao_l_ang_set
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine pao_l_ang_init(local_angles, qs_env, global_env, atom_nr,&
       glob_angles, error)
    type(pao_local_angles_type), intent(out):: local_angles
    type(qs_environment_type), intent(in), target :: qs_env
    type(global_environment_type), intent(in), target :: global_env
    integer, intent(in) :: atom_nr
    type(pao_glob_angles_type), target :: glob_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_l_angles_methods:pao_l_ang_init'
    failure=.false.

    local_angles%initialized=.true.
    local_angles%atom_nr=-1
    local_angles%global_env => global_env
    local_angles%qs_env => qs_env
    local_angles%glob_angles => glob_angles
    nullify(local_angles%projection, local_angles%tmp_m)
    call pao_l_ang_set(local_angles, atom_nr, error=error)
  end subroutine pao_l_ang_init
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_dealloc_ref [1.0] *
!!
!!   NAME
!!     pao_l_ang_dealloc_ref
!!
!!   SYNOPSIS
!!     Subroutine pao_l_ang_dealloc_ref(local_angles, error)
!!       Type(pao_local_angles_type), Intent (INOUT):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_l_ang_dealloc_ref
!!
!!   FUNCTION
!!     deallocates the memory that the given local angles structure
!!     has allocated. Does not deallocate the local angles structure itself
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles structure to be deallocated
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine pao_l_ang_dealloc_ref(local_angles, error)
    real(kind=wp) ::res
    type(pao_local_angles_type), intent(inout) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: stat
    character(len=*), parameter :: routineP='pao_l_angles_methods:pao_l_ang_dealloc_ref'
    failure=.false.

    CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       if (associated(local_angles%tmp_m)) then
          deallocate(local_angles%tmp_m,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       end if
       nullify(local_angles%glob_angles, local_angles%projection, &
            local_angles%qs_env, local_angles%global_env, local_angles%tmp_m)
       local_angles%initialized=.false.
    end if failureIf
  end subroutine pao_l_ang_dealloc_ref
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_validate [1.0] *
!!
!!   NAME
!!     pao_l_ang_validate
!!
!!   SYNOPSIS
!!     Function pao_l_ang_validate(local_angles, error) Result(res)
!!       Logical:: res
!!       Type(pao_local_angles_type), Intent (IN):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_l_ang_validate
!!
!!   FUNCTION
!!     checks if the local angles object is valid (returns .true. if it is)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles object to check
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function pao_l_ang_validate(local_angles, error) result(res)
    logical ::res
    type(pao_local_angles_type), intent(in) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_l_angles_methods:pao_l_ang_validate'
    type(pao_env_type), pointer :: pao_env
    integer :: i, beginAt
    failure=.false.
    nullify(pao_env)

    CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
    CPAssert(associated(local_angles%glob_angles),cp_warning_level,routineP,error,failure)
    CPAssert(associated(local_angles%qs_env),cp_warning_level,routineP,error,failure)
    CPAssert(associated(local_angles%global_env),cp_warning_level,routineP,error,failure)
    if (cp_debug .and. .not. failure) then
       call cp_assert(cp_validate(local_angles%glob_angles,&
            qs_env=local_angles%qs_env,global_env=local_angles%global_env,&
            error=error),cp_warning_level,cp_assertion_failed,routineP,&
            "referenced pao_glob_angles is invalid, in file "//&
            CPSourceFileRef,&
            error,failure)
    end if
    if (.not.failure) then
       CPAssert(local_angles%atom_nr>0,cp_warning_level,routineP,error,failure)
!       CPAssert(local_angles%atom_nr<=cp_get_n_atom(pao_env),cp_warning_level,routineP,error,failure)
       if (.not.failure) then
          CPAssert(local_angles%begins_at>=0,cp_warning_level,routineP,error,failure)
          call cp_assert(local_angles%begins_at < &
               size(local_angles%glob_angles%angles),&
               cp_warning_level,cp_assertion_failed,routineP,&
               "begins_at out of range in "//&
               CPSourceFileRef,&
               error,failure)
          call get_qs_env(local_angles%qs_env, pao_env=pao_env)
          beginAt=pao_env%angles_begin_at(local_angles%atom_nr)
          CPAssert(local_angles%begins_at==beginAt,cp_warning_level,routineP,error,failure)
       end if
    end if
    res=.not.failure
  end function pao_l_ang_validate
!***************************************************************************


!!****f* pao_l_angles_methods/pao_l_ang_set [1.0] *
!!
!!   NAME
!!     pao_l_ang_set
!!
!!   SYNOPSIS
!!     Subroutine pao_l_ang_set(local_angles, atom_nr, error)
!!       Type(pao_local_angles_type), Intent (INOUT):: local_angles
!!       Integer, Intent (IN):: atom_nr
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_l_ang_set
!!
!!   FUNCTION
!!     sets the values of the various attributes of the local angles
!!     structure.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the object to change
!!     atom_nr: the numbre of the atom this object focalizes on
!!     glob_angles: the global angles object that this object references
!!          SHARED.
!!     source: a local angles object that is copied into this(other arguments
!!             have the precedence for setting the values, the values
!!             of source are taken as default, if the value is not directly
!!             set)
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine pao_l_ang_set(local_angles, atom_nr,error)
    type(pao_local_angles_type), intent(inout):: local_angles
    integer, intent(in) :: atom_nr
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_l_angles_methods:pao_l_ang_set'
    type(pao_env_type), pointer :: pao_env
    type(particle_type), dimension(:), pointer :: particle_set
    type(cp_error_type) :: iError
    failure=.false.
    nullify(pao_env)
    call cp_error_init(iError, template_error=error)

    CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       local_angles%atom_nr=atom_nr
       CPAssert(atom_nr>0,cp_warning_level,routineP,error,failure)
       if (local_angles%atom_nr>0) then
          call cp_assert(associated(local_angles%glob_angles),&
               cp_failure_level,cp_precondition_failed,routineP,&
               "global_angles not associated in "//&
               CPSourceFileRef,&
               error,failure)
          call get_qs_env(local_angles%qs_env,pao_env=pao_env,&
               particle_set=particle_set)
          CPPrecondition(associated(pao_env),cp_failure_level,routineP,error,failure)
          CPPrecondition(associated(particle_set),cp_failure_level,routineP,error,failure)
          if (.not.failure) then
             local_angles%begins_at= pao_env%angles_begin_at(atom_nr)
             call get_atomic_kind(particle_set(atom_nr)%atomic_kind,&
                  pao_bas_proj=local_angles%projection)
             CPPostcondition(associated(local_angles%projection),cp_failure_level,routineP,error,failure)
          end if
          if (.not.failure) then
             call cp_proj_get(local_angles%projection,&
                  size_min_bas=local_angles%size_min_bas, error=error)
          end if
       else
          nullify(local_angles%projection)
       end if
    end if failureIf
    if (failure) then
       nullify(local_angles%projection)
    end if
    call cp_dealloc_ref(iError,error=error)
  end subroutine pao_l_ang_set
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get [1.0] *
!!
!!   NAME
!!     pao_l_ang_get
!!
!!   SYNOPSIS
!!     Subroutine pao_l_ang_get(local_angles, atom_nr, pao_env, unitary_t,&
!!         cache_valid, size_min_bas, size_full_bas, size_excl_bas,&
!!         proj_indexes, non_proj_indexes, projection, angles,&
!!         nui_injection, atomic_ortho, error)
!!       Type(pao_local_angles_type), Intent (INOUT):: local_angles
!!       Integer, Optional, Intent (OUT):: atom_nr
!!       Integer, Optional, Intent (OUT):: size_excl_bas, size_full_bas,&
!!         size_min_bas
!!       Integer, Optional, Dimension(:), Pointer:: proj_indexes,&
!!         non_proj_indexes
!!       Type(pao_env_type), Optional, Pointer:: pao_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Real(Kind=wp), Dimension(:,:), Pointer, Optional:: unitary_t,&
!!         atomic_ortho, nui_injection
!!       Logical, Intent (OUT), Optional:: cache_valid
!!       Type(pao_bas_projection_type), Pointer, Optional:: projection
!!       Real(Kind=wp), Dimension(:), Pointer, Optional:: angles
!!     End Subroutine pao_l_ang_get
!!
!!   FUNCTION
!!     gets various attributes of the pao local angles structure
!!
!!   NOTES
!!     I try to keep the glob_angles private.
!!     I have not added al local angles matrix
!!
!!   INPUTS
!!     local_angles: the local angles to get the info from.
!!     atom_nr: the atom this local angles is focalized on.
!!     n_min_bas: the size of the minimal basis on this atom.
!!     n_full_bas: the size of the full basis on this atom.
!!     n_excl_bas: the size of the excluded basis on this atom
!!            (n_full_bas-n_min_bas).
!!     pao_env: the pao environment this local angles lives in.
!!     angles: the local angles in a min_bas first order
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine pao_l_ang_get(local_angles, atom_nr, &
       pao_env, unitary_t,cache_valid,&
       size_min_bas, size_full_bas, &
       size_excl_bas, proj_indexes, non_proj_indexes,&
       projection, angles, NUi_injection, atomic_ortho, error)
    type(pao_local_angles_type), intent(inout) :: local_angles
    integer, optional, intent(out) :: atom_nr
    integer, optional, intent(out) :: size_excl_bas,size_full_bas,&
         size_min_bas
    integer, optional, dimension(:), pointer :: proj_indexes,&
         non_proj_indexes
    type(pao_env_type), optional, pointer :: pao_env
    type(cp_error_type), optional, intent(inout) :: error
    real(kind=wp), dimension(:,:), pointer, optional :: unitary_t,&
         atomic_ortho, NUi_injection
    logical, intent(out), optional :: cache_valid
    type(pao_bas_projection_type), pointer, optional :: projection
    real(kind=wp), dimension(:), pointer, optional :: angles

    logical :: failure
    character(len=*), parameter :: routineP='pao_l_angles_methods:pao_l_ang_get'
    type(pao_env_type), pointer :: my_pao_env
    failure=.false.
    nullify(my_pao_env)

    CPAssert(local_angles%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       if (present(atom_nr)) then
          atom_nr=local_angles%atom_nr
       end if
       if (present(pao_env)) then
          call get_qs_env(local_angles%qs_env, pao_env=pao_env)
       end if
       if (present(projection)) projection => local_angles%projection
       CPPrecondition(associated(local_angles%projection),cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
       call cp_proj_get(local_angles%projection, size_min_bas=size_min_bas,&
            size_full_bas=size_full_bas, size_excl_bas=size_excl_bas,&
            proj_indexes=proj_indexes, non_proj_indexes=non_proj_indexes,&
            error=error)
       if (present(unitary_t)) then
          unitary_t => pao_l_ang_get_unitary_t(local_angles,error=error)
       end if
       if (present(cache_valid)) then
          if (associated(local_angles%glob_angles%ut_calculated).and.&
               local_angles%atom_nr>0) then
             cache_valid=local_angles%glob_angles%ut_calculated(&
                  local_angles%atom_nr)
          else
             cache_valid=.false.
          end if
       end if
       if (present(angles)) then
          nullify(angles)
          ! seems to give problems on some platforms
          CPAssert(.false.,cp_warning_level,routineP,error,failure)
          call get_qs_env(local_angles%qs_env, pao_env=my_pao_env)
          angles => local_angles%glob_angles%angles(local_angles%begins_at:&
               (my_pao_env%angles_begin_at(local_angles%atom_nr+1)-1))
       end if
       if (present(NUi_injection)) then
          NUi_injection => pao_l_ang_get_NUi_injection(local_angles,error=error)
       end if
       if (present(atomic_ortho)) then
          atomic_ortho => pao_l_ang_get_atomic_ortho(local_angles,error=error)
       end if
    end if
  end subroutine pao_l_ang_get
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_angles_valid [1.0] *
!!
!!   NAME
!!     pao_l_angles_valid
!!
!!   FUNCTION
!!     performs minimal validation on the object
!!
!!   NOTES
!!     do not check pao_env?
!!
!!   INPUTS
!!     local_angles: the object to validate
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function pao_l_angles_valid(local_angles, error) result(res)
    logical ::res
    type(pao_local_angles_type), intent(in) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_l_angles_methods:pao_l_angles_valid'
    failure=.false.

    CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
    CPAssert(associated(local_angles%glob_angles),cp_warning_level,routineP,error,failure)
    res=.not.failure
  end function pao_l_angles_valid
!***************************************************************************

! ====== special methods =====

!!****f* pao_l_angles_methods/pao_l_ang_get_p_indexes [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_p_indexes
!!
!!   SYNOPSIS
!!     Function pao_l_ang_get_p_indexes(local_angles, error) Result(res)
!!       Integer, Dimension(:), Pointer:: res
!!       Type(pao_local_angles_type), Intent (IN):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_l_ang_get_p_indexes
!!
!!   FUNCTION
!!     returns the indexes of the projection
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles you get values from
!!     i: number of the index you want to get
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     1.2002 created [fawzi]
!!     4.2002 use bas_projection [fawzi]
!!
!!*** **********************************************************************
  function pao_l_ang_get_p_indexes(local_angles,error) result(res)
    integer, dimension(:), pointer ::res
    type(pao_local_angles_type), intent(in) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='pao_l_ang_get_p_indexes',&
         routineP=moduleN//':'//routineN
    failure=.false.
    if (cp_debug .or. debug_this_module) then
       CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(local_angles%projection),cp_failure_level,routineP,error,failure)
    end if

    if (.not. failure) then
       call cp_proj_get(local_angles%projection, proj_indexes=res,error=error)
    else
       nullify(res)
    end if
  end function pao_l_ang_get_p_indexes
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get_non_p_indexes [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_non_p_indexes
!!
!!   SYNOPSIS
!!     Function pao_l_ang_get_non_p_indexes(local_angles, error)&
!!         Result(res)
!!       Integer, Dimension(:), Pointer:: res
!!       Type(pao_local_angles_type), Intent (IN):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_l_ang_get_non_p_indexes
!!
!!   FUNCTION
!!     returns the indexes of the excluded basis
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles you get values from
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function pao_l_ang_get_non_p_indexes(local_angles,error) result(res)
    integer, dimension(:), pointer ::res
    type(pao_local_angles_type), intent(in) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='pao_l_ang_get_p_indexes',&
         routineP=moduleN//':'//routineN
    failure=.false.
    if (cp_debug .or. debug_this_module) then
       CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(local_angles%projection),cp_failure_level,routineP,error,failure)
    end if

    if (.not. failure) then
       call cp_proj_get(local_angles%projection, non_proj_indexes=res,error=error)
    else
       nullify(res)
    end if
  end function pao_l_ang_get_non_p_indexes
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get_p_index [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_p_index
!!
!!   SYNOPSIS
!!     Function pao_l_ang_get_p_index(local_angles, min_bas, error)&
!!         Result(res)
!!       Integer:: res
!!       Type(pao_local_angles_type), Intent (IN):: local_angles
!!       Integer, Intent (IN):: min_bas
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_l_ang_get_p_index
!!
!!   FUNCTION
!!     returns the index min_bas of the projection
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles you get values from
!!     min_bas: the index of the minimal basis (1:size_min_bas)
!!     i: number of the index you want to get
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     1.2002 created [fawzi]
!!     4.2002 use bas_projection [fawzi]
!!
!!*** **********************************************************************
  function pao_l_ang_get_p_index(local_angles, min_bas, error) result(res)
    integer ::res
    type(pao_local_angles_type), intent(in) :: local_angles
    integer, intent(in) :: min_bas
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='pao_l_ang_get_p_index',&
         routineP=moduleN//':'//routineN
    integer, dimension(:), pointer :: p
    failure=.false.
    if (cp_debug .or. debug_this_module) then
       CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(local_angles%projection),cp_failure_level,routineP,error,failure)
    end if

    if (.not. failure) then
       call cp_proj_get(local_angles%projection, proj_indexes=p,error=error)
       res=p(min_bas)
    else
       res=-1
    end if
  end function pao_l_ang_get_p_index
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get_non_p_index [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_non_p_index
!!
!!   SYNOPSIS
!!     Function pao_l_ang_get_non_p_index(local_angles, excl_bas, error)&
!!         Result(res)
!!       Integer:: res
!!       Integer, Intent (IN):: excl_bas
!!       Type(pao_local_angles_type), Intent (IN):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_l_ang_get_non_p_index
!!
!!   FUNCTION
!!     returns the index excl_bas of the excluded basis
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - local_angles: the local angles you get values from
!!     - excl_bas: the local index on the excluded basis (1:size_excl_bas)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     1.2002 created [fawzi]
!!     4.2002 use bas projection [fawzi]
!!
!!*** **********************************************************************
  function pao_l_ang_get_non_p_index(local_angles,excl_bas,error) result(res)
    integer ::res
    integer, intent(in) :: excl_bas
    type(pao_local_angles_type), intent(in) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='pao_l_ang_get_p_index',&
         routineP=moduleN//':'//routineN
    integer, dimension(:), pointer :: np
    failure=.false.
    if (cp_debug .or. debug_this_module) then
       CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(local_angles%projection),cp_failure_level,routineP,error,failure)
    end if

    if (.not. failure) then
       call cp_proj_get(local_angles%projection, non_proj_indexes=np,error=error)
       res=np(excl_bas)
    else
       res=-1
    end if
  end function pao_l_ang_get_non_p_index
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get_angle [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_angle
!!
!!   SYNOPSIS
!!     Function pao_l_ang_get_angle(local_angles, min_bas, excl_bas, error)&
!!         Result(res)
!!       Real(Kind=wp):: res
!!       Type(pao_local_angles_type), Intent (IN):: local_angles
!!       Integer, Intent (IN):: min_bas, excl_bas
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_l_ang_get_angle
!!
!!   FUNCTION
!!     returns the requested angle
!!
!!   NOTES
!!     should be fast, and hopefully inlined
!!
!!   INPUTS
!!     local_angles: the local angles object you want to get the angles from
!!     min_bas: the local index of the minimal basis of which you want 
!!            the angle
!!     excl_bas: the local index of the excluded basis of which you want the
!!            angle
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function pao_l_ang_get_angle(local_angles, min_bas, excl_bas, &
       error) result(res)
    real(kind=wp) ::res
    type(pao_local_angles_type), intent(in) :: local_angles
    integer, intent(in) :: min_bas, excl_bas
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_l_angles_methods:pao_l_ang_get_angle'
    integer :: size_excl_bas
    type(pao_glob_angles_type), pointer :: glob_angles
    failure=.false.
    glob_angles => local_angles%glob_angles

    debugIf: if (debug_this_module.and.cp_debug) then
       CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(glob_angles),cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          CPPrecondition(associated(local_angles%projection),cp_failure_level,routineP,error,failure)
       end if
       if (.not.failure) then
          call cp_proj_get(local_angles%projection,size_excl_bas=size_excl_bas,&
               error=error)
          call cp_assert(min_bas>0.and.min_bas <= &
               local_angles%size_min_bas,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "min_bas out of range in "//&
               CPSourceFileRef,&
               error,failure)
          call cp_assert(excl_bas>0.and.excl_bas <= &
               size_excl_bas,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "excl_bas out of range in "//&
               CPSourceFileRef,&
               error,failure)
       end if
    end if debugIf
    failureIf: if (.not. failure) then
       res=glob_angles%angles(local_angles%begins_at + &
            local_angles%size_min_bas*(excl_bas-1) + &
            min_bas)
    else failureIf
       res=-1
    end if failureIf
  end function pao_l_ang_get_angle
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_set_angle [1.0] *
!!
!!   NAME
!!     pao_l_ang_set_angle
!!
!!   SYNOPSIS
!!     Subroutine pao_l_ang_set_angle(local_angles, min_bas, excl_bas,&
!!         value, error)
!!       Type(pao_local_angles_type), Intent (IN):: local_angles
!!       Integer, Intent (IN):: min_bas, excl_bas
!!       Real(Kind=wp), Intent (IN):: value
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_l_ang_set_angle
!!
!!   FUNCTION
!!     returns the requested angle
!!
!!   NOTES
!!     should be fast, and hopefully inlined
!!
!!   INPUTS
!!     local_angles: the local angles object you want to set the angles from
!!     min_bas: the local index of the minimal basis of which you want 
!!            the angle
!!     full_bas: the local index of the full basis of which you want the
!!            angle
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine pao_l_ang_set_angle(local_angles, min_bas, excl_bas, &
       value,error)
    type(pao_local_angles_type), intent(in) :: local_angles !actually the substructure globalangles is modifed
    integer, intent(in) :: min_bas, excl_bas
    real(kind=wp), intent(in) :: value
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='pao_l_ang_set_angle',&
         routineP=moduleN//':'//routineN
    type(pao_glob_angles_type), pointer :: glob_angles
    integer :: size_excl_bas
    failure=.false.
    glob_angles => local_angles%glob_angles

    debugIf: if (debug_this_module.and.cp_debug) then
       CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(glob_angles),cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          call cp_assert(min_bas>0.and.min_bas <= &
               local_angles%size_min_bas,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "min_bas out of range in "//&
               CPSourceFileRef,&
               error,failure)
          call cp_proj_get(local_angles%projection,size_excl_bas=size_excl_bas,&
               error=error)
          call cp_assert(excl_bas>0.and.excl_bas <= &
               size_excl_bas,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "excl_bas out of range in "//&
               CPSourceFileRef,&
               error,failure)
       end if
    end if debugIf
    failureIf: if (.not. failure) then
       glob_angles%angles(local_angles%begins_at + &
            local_angles%size_min_bas*(excl_bas-1) + &
            min_bas) = value
    end if failureIf
    nullify(glob_angles)
  end subroutine pao_l_ang_set_angle
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get_size_excl_bas [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_size_excl_bas
!!
!!   SYNOPSIS
!!     Function pao_l_ang_get_size_excl_bas(local_angles, error)&
!!         Result(res)
!!       Integer:: res
!!       Type(pao_local_angles_type), Intent (IN):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_l_ang_get_size_excl_bas
!!
!!   FUNCTION
!!     returns the size of the excluded basis
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles object you get the information from
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function pao_l_ang_get_size_excl_bas(local_angles, error) result(res)
    integer ::res
    type(pao_local_angles_type), intent(in) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_l_angles_methods:pao_l_ang_get_size_excl_bas'
    failure=.false.

    CPPrecondition(local_angles%initialized, cp_failure_level,routineP,error,failure)
    CPPrecondition(associated(local_angles%projection),cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       call cp_proj_get(local_angles%projection,size_excl_bas=res,error=error)
    else
       res=-1
    end if
  end function pao_l_ang_get_size_excl_bas
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get_size_min_bas [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_size_min_bas
!!
!!   SYNOPSIS
!!     Function pao_l_ang_get_size_min_bas(local_angles, error) Result(res)
!!       Integer:: res
!!       Type(pao_local_angles_type), Intent (IN):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_l_ang_get_size_min_bas
!!
!!   FUNCTION
!!     returns the size of the minimal basis
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles object you get the information from
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function pao_l_ang_get_size_min_bas(local_angles, error) result(res)
    integer ::res
    type(pao_local_angles_type), intent(in) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_l_angles_methods:pao_l_ang_get_size_min_bas'
    failure=.false.

    CPPrecondition(local_angles%initialized, cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       res=local_angles%size_min_bas
    else
       res=-1
    end if
  end function pao_l_ang_get_size_min_bas
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get_size_full_bas [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_size_full_bas
!!
!!   SYNOPSIS
!!     Function pao_l_ang_get_size_full_bas(local_angles, error)&
!!         Result(res)
!!       Integer:: res
!!       Type(pao_local_angles_type), Intent (IN):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_l_ang_get_size_full_bas
!!
!!   FUNCTION
!!     returns the size of the full basis
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles object you get the information from
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function pao_l_ang_get_size_full_bas(local_angles, error) result(res)
    integer ::res
    type(pao_local_angles_type), intent(in) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_l_angles_methods:pao_l_ang_get_size_full_bas'
    failure=.false.

    CPPrecondition(local_angles%initialized, cp_failure_level,routineP,error,failure)
    CPPrecondition(associated(local_angles%projection),cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       call cp_proj_get(local_angles%projection,size_full_bas=res,error=error)
    else
       res=-1
    end if
  end function pao_l_ang_get_size_full_bas
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_did_change [1.0] *
!!
!!   NAME
!!     pao_l_ang_did_change
!!
!!   SYNOPSIS
!!     Subroutine pao_l_ang_did_change(local_angles, error)
!!       Type(pao_local_angles_type), Intent (INOUT):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_l_ang_did_change
!!
!!   FUNCTION
!!     invalidate the cached unitary transformation of this atom.
!!     to be called when the angles change
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles object you get the info from
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
  subroutine pao_l_ang_did_change(local_angles,error)
! changes in local_angles%glob_angles
    type(pao_local_angles_type), intent(inout) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_l_angles_methods:pao_l_ang_did_change'
    logical, dimension(:), pointer :: ut_calculated
    failure=.false.

    CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(local_angles%atom_nr>0,cp_warning_level,routineP,error,failure)
    if (.not. failure) then
       if (associated(local_angles%glob_angles%ut_calculated)) then
          local_angles%glob_angles%ut_calculated(local_angles%atom_nr)=.false.
       end if
       if (associated(local_angles%glob_angles%NUi_injection_calculated)) then
          local_angles%glob_angles%NUi_injection_calculated&
               (local_angles%atom_nr)=.false.
       end if
    end if
  end subroutine pao_l_ang_did_change
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_set_unitary_t [1.0] *
!!
!!   NAME
!!     pao_l_ang_set_unitary_t
!!
!!   SYNOPSIS
!!     Subroutine pao_l_ang_set_unitary_t(local_angles, unitary_t, error)
!!       Type(pao_local_angles_type), Intent (INOUT):: local_angles
!!       Real(Kind=wp), Dimension(:,:), Intent (IN), Target:: unitary_t
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_l_ang_set_unitary_t
!!
!!   FUNCTION
!!     changes the cached unitary transformation
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles you want to change
!!     unitary_t: the new unitary trasformation
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
  subroutine pao_l_ang_set_unitary_t(local_angles,unitary_t,error)
! changes in glob_angles
    type(pao_local_angles_type), intent(inout) :: local_angles
    real(kind=wp), dimension(:,:), intent(in), target :: unitary_t
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_l_angles_methods:pao_l_ang_set_unitary_t'
    type(cp_sparse_matrix_type), pointer :: full_u_t
    real(kind=wp), dimension(:,:), pointer :: unitary_t_ptr
    failure=.false.

    CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(local_angles%atom_nr>0,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       call pao_g_ang_get(local_angles%glob_angles,qs_env=local_angles%qs_env,&
            global_env=local_angles%global_env, unitary_t=full_u_t,&
            error=error)
       CPPrecondition(associated(full_u_t),cp_warning_level,routineP,error,failure)
       call cp_assert(associated(local_angles%glob_angles%ut_calculated),&
            cp_failure_level,cp_assertion_failed,routineP,&
            "ut_calculated not allocated in "//&
            CPSourceFileRef,&
            error,failure)
       if (.not.failure) then
          local_angles%glob_angles%ut_calculated(local_angles%atom_nr)=.true.
          unitary_t_ptr => unitary_t
          call put_block_node(cp_get_matrix(full_u_t),&
               block_row=local_angles%atom_nr,&
               block_col=local_angles%atom_nr, block=unitary_t_ptr)
       end if
    end if
  end subroutine pao_l_ang_set_unitary_t
!***************************************************************************

!!****f* pao_plain_rot/pao_l_compose_u [1.0] *
!!
!!   NAME
!!     pao_l_compose_u
!!
!!   SYNOPSIS
!!     Subroutine pao_l_compose_u(matrix, local_angles, error)
!!       Real(Kind=wp), Dimension(:,:), Intent (INOUT):: matrix
!!       Type(pao_local_angles_type), Intent (INOUT):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_l_compose_u
!!
!!   FUNCTION
!!     calculates U*matrix where U is the unitary transformation defined by
!!     the given angles. the sequence is: smallbasis then full basis in 
!!     (local) growing order: u(1,1)*u(2,1)*...*u(size_small_bas,1)*...*
!!      u(size_small_bas,size_full_bas)*matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix: the matrix that will be modified
!!     local_angles: the angles that define U
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine pao_l_compose_u(matrix,local_angles,error)
    real(kind=wp), dimension(:,:), intent(inout) :: matrix
    type(pao_local_angles_type), intent(inout) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_plain_rot:pao_l_compose_u'
    integer :: n_full_bas, i_min_bas, i_excl_bas
    integer, dimension(:), pointer :: proj_indexes, non_proj_indexes
    failure=.false.

    call pao_l_ang_get(local_angles,&
         size_full_bas=n_full_bas, proj_indexes=proj_indexes,&
         non_proj_indexes=non_proj_indexes)
    CPPrecondition(size(proj_indexes)>0,cp_warning_level,routineP,error,failure)
    CPPrecondition(size(non_proj_indexes)>=0,cp_warning_level,routineP,error,failure)
    CPPrecondition(size(matrix,1)==n_full_bas,cp_warning_level,routineP,error,failure)
    CPPrecondition(size(matrix,2)==n_full_bas,cp_warning_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       do i_excl_bas=size(non_proj_indexes),1,-1
          do i_min_bas=size(proj_indexes),1,-1
             call cp_left_compose_p_rot(matrix=matrix,&
                  angle=cp_get_angle(local_angles,min_bas=i_min_bas,&
                  excl_bas=i_excl_bas,error=error),&
                  i=proj_indexes(i_min_bas),&
                  j=non_proj_indexes(i_excl_bas),&
                  error=error)
          end do
       end do
    end if failureIf
  end subroutine pao_l_compose_u
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get_unitary_t [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_unitary_t
!!
!!   SYNOPSIS
!!     Function pao_l_ang_get_unitary_t(local_angles, error) Result(res)
!!       Type(pao_local_angles_type), Intent (INOUT):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Real(Kind=wp), Dimension(:,:), Pointer:: res
!!     End Function pao_l_ang_get_unitary_t
!!
!!   FUNCTION
!!     returns the unitary trasformation defined by these angles
!!
!!   NOTES
!!     the trasformation is cached
!!
!!   INPUTS
!!     local_angles: the local angles you get the information from
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
  function pao_l_ang_get_unitary_t(local_angles,error) result(res)
    type(pao_local_angles_type), intent(inout) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error
    real(kind=wp), dimension(:,:), pointer :: res

    logical :: failure
    character(len=*), parameter :: routineP='pao_l_angles_methods:pao_l_ang_get_unitary_t'
    integer :: i
    type(cp_sparse_matrix_type), pointer :: unitary_t
    failure=.false.
    nullify(res)

    if (cp_debug .and. debug_this_module) then
       CPPrecondition(cp_valid(local_angles),cp_failure_level,routineP,error,failure)
    end if
    if (.not. failure) then
       call pao_g_ang_get(local_angles%glob_angles,qs_env=local_angles%qs_env,&
            global_env=local_angles%global_env, unitary_t=unitary_t,&
            error=error)
       CPPrecondition(associated(unitary_t),cp_failure_level,routineP,error,failure)
       call cp_assert(associated(local_angles%glob_angles%ut_calculated),&
            cp_failure_level,cp_assertion_failed,routineP,&
            "ut_calculated not associated in "//&
            CPSourceFileRef,&
            error,failure)
       if (.not.failure) then
          if (local_angles%glob_angles%ut_calculated(local_angles%atom_nr)) then
             call get_block_node(cp_get_matrix(unitary_t),&
                  block=res,block_row=local_angles%atom_nr,&
                  block_col=local_angles%atom_nr)
             CPPostcondition(associated(res),cp_failure_level,routineP,error,failure)
          else
             call get_block_node(cp_get_matrix(unitary_t),&
                  block=res,block_row=local_angles%atom_nr,&
                  block_col=local_angles%atom_nr)
             if (.not.associated(res)) then
                call add_block_node(cp_get_matrix(unitary_t),&
                     local_angles%atom_nr,local_angles%atom_nr,block=res)
             end if
             CPPostcondition(associated(res),cp_failure_level,routineP,error,failure)
             if (.not.failure) then
                CPPrecondition(size(res,1)==size(res,2),cp_failure_level,routineP,error,failure)
             end if
             if (.not.failure) then
                res=0.0_wp
                forall (i=1:size(res,1)) res(i,i)=1.0_wp
                call pao_l_compose_u(res,local_angles,error=error)
                local_angles%glob_angles%ut_calculated(local_angles%atom_nr)=.true.
             end if
          end if
       end if
    end if
    if (failure) nullify(res)
  end function pao_l_ang_get_unitary_t
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get_g_index [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_g_index
!!
!!   SYNOPSIS
!!     Function pao_l_ang_get_g_index(local_angles, min_bas, excl_bas,&
!!         error) Result(res)
!!       Integer:: res
!!       Type(pao_local_angles_type), Intent (INOUT):: local_angles
!!       Integer, Intent (IN):: min_bas, excl_bas
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_l_ang_get_g_index
!!
!!   FUNCTION
!!     returns the global angles index for the angle described with
!!     the given local indexes. 
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the actual local angles object
!!     min_bas: the local minimal basis index
!!     excl_bas: the local excluded basis index
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
function pao_l_ang_get_g_index(local_angles, min_bas, excl_bas, error)&
     result(res)
  integer :: res
  type(pao_local_angles_type), intent(inout) :: local_angles
  integer, intent(in) :: min_bas, excl_bas
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='pao_l_ang_get_g_index',&
       routineP=moduleN//':'//routineN
  failure=.false.
  
  if (cp_debug .and. debug_this_module) then
     CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
     CPPrecondition(local_angles%atom_nr>0,cp_failure_level,routineP,error,failure)
  end if
  if (.not. failure) then
     res = local_angles%begins_at+local_angles%size_min_bas*(excl_bas-1)&
          + min_bas
  else
     res=-1
  end if
end function pao_l_ang_get_g_index
!***************************************************************************

!!****f* pao_l_angles_methods [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_atomic_ortho
!!
!!   FUNCTION
!!     returns the atomic orthogonalization matrix, calculating it
!!     if necessary.
!!
!!   NOTES
!!     cache tmp_m (but submatrix indexing on most compilers is not optimized)
!!
!!   INPUTS
!!     - local_angles: the local_angles object that defines the actual atom
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     6.2002 created [fawzi]
!!
!!*** **********************************************************************
function pao_l_ang_get_atomic_ortho(local_angles,error) result(res)
  type(pao_local_angles_type), intent(inout) :: local_angles
  type(cp_error_type), optional, intent(inout) :: error
  real(kind=wp), dimension(:,:), pointer :: res

  logical :: failure
  character(len=*), parameter :: routineN='pao_l_ang_get_atomic_ortho',&
       routineP=moduleN//':'//routineN
  integer :: stat
  integer, save :: iter=0
  type(cp_sparse_matrix_type), pointer :: atomic_ortho
  type(pao_env_type), pointer :: pao_env
  type(real_matrix_p_type), dimension(:), pointer :: s
  real(kind=wp), dimension(:,:), pointer :: block_val, tmp_m
  real(kind=wp) :: eps_eigval
  failure=.false.
  nullify(res, atomic_ortho, block_val, tmp_m, pao_env, s)

  CPPrecondition(cp_valid(local_angles),cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     call pao_qs_env_get(local_angles%qs_env,&
          global_env=local_angles%global_env, pao_env=pao_env,&
          pao_atomic_ortho=atomic_ortho,error=error)
     CPPrecondition(associated(atomic_ortho),cp_failure_level,routineP,error,failure)
     CPPrecondition(associated(pao_env),cp_failure_level,routineP,error,failure)
     if (.not.failure) then
        call cp_assert(associated(pao_env%atomic_ortho_valid),&
             cp_failure_level,cp_assertion_failed,routineP,&
             "atomic_ortho_valid not associated in "//&
             CPSourceFileRef,&
             error,failure)
     end if
  end if
  if (.not.failure) then
     if (pao_env%atomic_ortho_valid(local_angles%atom_nr)) then
        call get_block_node(cp_get_matrix(atomic_ortho),&
             block=res,block_row=local_angles%atom_nr,&
             block_col=local_angles%atom_nr)
        CPPostcondition(associated(res),cp_failure_level,routineP,error,failure)
     else
        call get_qs_env(local_angles%qs_env,s=s)
        call get_block_node(cp_get_matrix(atomic_ortho),&
             block=res,block_row=local_angles%atom_nr,&
             block_col=local_angles%atom_nr)
        if (.not.associated(res)) then
           call add_block_node(cp_get_matrix(atomic_ortho),&
                local_angles%atom_nr,local_angles%atom_nr,block=res)
        end if
        CPPostcondition(associated(res),cp_failure_level,routineP,error,failure)
        if (.not.failure) then
           CPPrecondition(size(res,1)==size(res,2),cp_failure_level,routineP,error,failure)
        end if
        if (.not.failure) then
           allocate(tmp_m(size(res,1),size(res,2)),stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        end if
        if (.not.failure) then
           call get_block_node(s(1)%matrix,block_row=local_angles%atom_nr,&
                block_col=local_angles%atom_nr,block=block_val)
           CPPrecondition(size(block_val,1)==size(res,1),cp_failure_level,routineP,error,failure)
           CPPrecondition(size(block_val,2)==size(res,2),cp_failure_level,routineP,error,failure)
           if (.not.failure) then
              tmp_m=block_val
              eps_eigval= 10 * epsilon(0.0_wp) ! to do: use the scf value
              call power_matrix(tmp_m,res,-0.5_wp,threshold=eps_eigval)
              call symmetrize_matrix(res,"upper_to_lower") ! remove?
              if (cp_debug.and.debug_this_module) then
                 iter=iter+1
                 call cp_2d_r_output(cp_error_get_logger(error),&
                      outputName="atomic_ortho",&
                      fromWhere=routineP, iter=iter,array=res,&
                      comment="atom_nr="//cp_to_string(local_angles%atom_nr),&
                      error=error)
              end if
           end if
           deallocate(tmp_m,stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        end if
     end if
  end if
  if (failure) nullify(res)
end function pao_l_ang_get_atomic_ortho
!***************************************************************************

!!****f* pao_l_angles_methods/pao_l_ang_get_NUi_injection [1.0] *
!!
!!   NAME
!!     pao_l_ang_get_NUi_injection
!!
!!   FUNCTION
!!     returns the injection from the small basis to the full basis,
!!     calculating it if necessary (N ortho, U unitary, i basis injection)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - local_angles: the actual local angles object that defines
!!       the actual atom 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
function pao_l_ang_get_NUi_injection(local_angles, error) result(res)
  real(kind=wp), dimension(:,:), pointer ::res
  type(pao_local_angles_type), intent(inout) :: local_angles
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='pao_l_ang_get_NUi_injection',&
        routineP=moduleN//':'//routineN
  integer :: stat, max_full_bas, i
  integer, save :: iter=0
  integer, dimension(:), pointer :: proj_indexes
  type(cp_sparse_matrix_type), pointer :: NUi_injection
  type(pao_env_type), pointer :: pao_env
  real(kind=wp), dimension(:,:), pointer :: n_m,u_m
  failure=.false.
  nullify(res, NUi_injection, n_m, u_m, pao_env)

  CPPrecondition(cp_valid(local_angles),cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     call pao_qs_env_get(local_angles%qs_env,&
          global_env=local_angles%global_env, pao_env=pao_env,&
          pao_max_full_bas=max_full_bas,&
          error=error)
     CPPrecondition(associated(NUi_injection),cp_failure_level,routineP,error,failure)
     CPPrecondition(associated(pao_env),cp_failure_level,routineP,error,failure)
     if (.not.failure) then
        call cp_assert(associated(local_angles%glob_angles%NUi_injection_calculated),&
             cp_failure_level,cp_assertion_failed,routineP,&
             "NUi_injection_calculated not associated in "//&
             CPSourceFileRef,&
             error,failure)
     end if
  end if
  if (.not.failure) then
     if (local_angles%glob_angles%NUi_injection_calculated&
          (local_angles%atom_nr)) then
        call get_block_node(cp_get_matrix(NUi_injection),&
             block=res,block_row=local_angles%atom_nr,&
             block_col=local_angles%atom_nr)
        CPPostcondition(associated(res),cp_failure_level,routineP,error,failure)
     else
        call get_block_node(cp_get_matrix(NUi_injection),&
             block=res,block_row=local_angles%atom_nr,&
             block_col=local_angles%atom_nr)
        if (.not.associated(res)) then
           call add_block_node(cp_get_matrix(NUi_injection),&
                local_angles%atom_nr,local_angles%atom_nr,block=res)
        end if
        CPPostcondition(associated(res),cp_failure_level,routineP,error,failure)
        if (.not.failure) then
           CPPrecondition(size(res,1)==size(res,2),cp_failure_level,routineP,error,failure)
        end if
        if (.not.failure) then
           if (associated(local_angles%tmp_m)) then
              call cp_assert(size(local_angles%tmp_m,1) >= size(res,1) .and.&
                   size(local_angles%tmp_m,2) >= size(res,2),cp_failure_level,&
                   cp_internal_error, routineP, &
                   "local_angles%tmp_m not big enough, in"//&
                   CPSourceFileRef,&
                   error=error,failure=failure)
           else
              allocate(local_angles%tmp_m(max_full_bas,max_full_bas),stat=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           end if
        end if
        u_m => pao_l_ang_get_unitary_t(local_angles, error=error)
        CPPrecondition(associated(u_m),cp_failure_level,routineP,error,failure)
        n_m => pao_l_ang_get_atomic_ortho(local_angles, error=error)
        CPPrecondition(associated(n_m),cp_failure_level,routineP,error,failure)
        if (.not.failure) then
           CPPrecondition(size(u_m,1)==size(n_m,2),cp_failure_level,routineP,error,failure)
        end if
        if (.not.failure) then
           call dgemm('N','N',size(n_m,1),size(u_m,2),size(u_m,2),1.0_wp,&
                n_m,size(n_m,1),u_m,size(u_m,1),0.0_wp,&
                local_angles%tmp_m,size(local_angles%tmp_m,1))
           call cp_proj_get(local_angles%projection,proj_indexes=proj_indexes,&
                error=error)
           CPPrecondition(associated(proj_indexes),cp_failure_level,routineP,error,failure)
           if (.not.failure) then
              do i=1,size(proj_indexes)
                 res(:,i)=local_angles%tmp_m(:,proj_indexes(i))
              end do
              if (cp_debug.and.debug_this_module) then
                 iter=iter+1
                 call cp_2d_r_output(cp_error_get_logger(error),&
                      outputName="atomic_NUi",&
                      fromWhere=routineP, iter=iter,array=res,&
                      comment="atom_nr="//cp_to_string(local_angles%atom_nr),&
                      error=error)
              end if
           end if
        end if
     end if
  end if
  if (failure) nullify(res)
end function pao_l_ang_get_NUi_injection
!***************************************************************************


end module pao_l_angles_methods
