!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_pao_l_angles_methods [1.0] *
!!
!!   NAME
!!     cp_pao_l_angles_methods
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
module cp_pao_l_angles_methods
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use cp_pao_types
  use cp_pao_glob_angles_methods
  use cp_pao_env_methods
  use matrix_types, only: real_matrix_type, replicate_matrix, copy_matrix,&
       allocate_matrix, deallocate_matrix, get_matrix_info, get_block_node,&
       add_block_node, put_block_node
  use cp_pao_matrix_methods
  use cp_pao_plain_rot
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_pao_env_methods'

! core procedures
  public :: cp_init, cp_dealloc_ref, cp_set, cp_get, &
       cp_did_change

! core functions
  public :: cp_validate, cp_valid

! special procedures
  public :: cp_set_angle

! special functions
  public :: cp_get_pao_env, cp_get_angle, &
       cp_get_min_bas_size, cp_get_full_bas_size, cp_get_excl_bas_size, &
       cp_get_p_index, cp_set_p_index, &
       cp_get_non_p_index, cp_get_unitary_t,cp_get_glob_angles,&
       cp_get_g_index

  public :: cp_pao_l_compose_u ! to do remove from here and plain rot tests

! underlying functions
  public :: cp_pao_l_ang_init, cp_pao_l_ang_dealloc_ref,&
       cp_pao_l_ang_validate, cp_pao_l_ang_set,&
       cp_pao_l_ang_get, cp_pao_l_ang_did_change,&
       cp_pao_l_ang_get_pao_env, cp_pao_l_ang_get_angle,&
       cp_pao_l_ang_set_angle, cp_pao_l_ang_get_min_bas_size,&
       cp_pao_l_ang_get_full_bas_size, cp_pao_l_ang_get_excl_bas_size,&
       cp_pao_l_ang_get_p_index, cp_pao_l_ang_set_p_index,&
       cp_pao_l_ang_get_non_p_index,cp_pao_l_ang_get_unitary_t,&
       cp_pao_l_ang_set_unitary_t, cp_pao_l_ang_get_glob_angles,&
       cp_pao_l_ang_get_g_index

! ============== interfaces ===========

! +++ begin of core methods +++

!! initializes the first, argument.
!! Support optional initial values that depend on the 
!! first argument 
  interface cp_init
     module procedure cp_pao_l_ang_init
  end interface

!! Deallocates the memory allocated by the first argument.
!! Supports error=error
  interface cp_dealloc_ref
     module procedure cp_pao_l_ang_dealloc_ref
  end interface

!! returns true if the content of the first argumenst is valid.
!! writes the errors ar warnings, full validation (long)
  interface cp_validate
     module procedure cp_pao_l_ang_validate
  end interface

!! returns true if the first argument is valid. Only a fast minimal 
!! validation is performed (no dangling pointers)
  interface cp_valid
     module procedure cp_pao_l_angles_valid
  end interface

!! sets the value of various optonal attributes of the firs argument.
!! the optional attributes depend on the type of the first
  interface cp_set
     module procedure cp_pao_l_ang_set
  end interface

!! gets values from the first argument via optional keyword arguments
!! that depend on the type of the first argument
  interface cp_get
     module procedure cp_pao_l_ang_get
  end interface

!! tells the object that its internal values have changes and its cached
!! values may be invalid
  interface cp_did_change
     module procedure cp_pao_l_ang_did_change
  end interface

! +++ end of the core methods +++

!! returns the pao environement the first object lives in
  interface cp_get_pao_env
     module procedure cp_pao_l_ang_get_pao_env
  end interface

!! returns the angle in the first arg identified by the second (min_bas)
!! and the third arg (full_bas) 
  interface cp_get_angle
     module procedure cp_pao_l_ang_get_angle
  end interface

!! sets the angle in the first arg identified by the second (min_bas)
!! and the third arg (full_bas) to the value in the fourth (value)
  interface cp_set_angle
     module procedure cp_pao_l_ang_set_angle
  end interface

!! returns the size of the minimal basis
  interface cp_get_min_bas_size
     module procedure cp_pao_l_ang_get_min_bas_size
  end interface

!! returns the size of the full basis
  interface cp_get_full_bas_size
     module procedure cp_pao_l_ang_get_full_bas_size
  end interface

!! returns the size of the excluded basis
  interface cp_get_excl_bas_size
     module procedure cp_pao_l_ang_get_excl_bas_size
  end interface

!! returns the projection index
  interface cp_get_p_index
     module procedure cp_pao_l_ang_get_p_index
  end interface

!! sets the projection index
  interface cp_set_p_index
     module procedure cp_pao_l_ang_set_p_index
  end interface

!! returns the index of the excluded basis
  interface cp_get_non_p_index
     module procedure cp_pao_l_ang_get_non_p_index
  end interface

!! gets the unitary transformation
  interface cp_get_unitary_t
     module procedure cp_pao_l_ang_get_unitary_t
  end interface

!! sets the unitary transformation
  interface cp_set_unitary_t
     module procedure cp_pao_l_ang_set_unitary_t
  end interface

!! gets the undelying glob_angles object
  interface cp_get_glob_angles
     module procedure cp_pao_l_ang_get_glob_angles
  end interface

!! returns the global index
interface cp_get_g_index
   module procedure cp_pao_l_ang_get_g_index
end interface

!!***
!****************************************************************************

contains

! ====== core methods ======

!!****f* cp_pao_l_angles_methods/cp_pao_l_ang_init [1.0] *
!!
!!   NAME
!!     cp_pao_l_ang_init
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_l_ang_init(local_angles, atom_nr, glob_angles,&
!!         source, error)
!!       Type(cp_pao_local_angles_type), Intent (OUT):: local_angles
!!       Integer, Optional, Intent (IN):: atom_nr
!!       Type(cp_pao_glob_angles_type), Optional, Target:: glob_angles
!!       Type(cp_pao_local_angles_type), Intent (IN), Optional:: source
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_l_ang_init
!!
!!   FUNCTION
!!     Initializes the pao angles structure
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the pao_local angles to be initialized
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling
!!
!!     For the other arguments see cp_pao_l_ang_set
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_pao_l_ang_init(local_angles, atom_nr, glob_angles, &
       source, error)
    type(cp_pao_local_angles_type), intent(out):: local_angles
    integer, optional, intent(in) :: atom_nr
    type(cp_pao_glob_angles_type), optional, target :: glob_angles
    type(cp_pao_local_angles_type), intent(in), optional :: source
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_l_angles_methods:cp_pao_l_ang_init'
    failure=.false.

    local_angles%initialized=.true.
    local_angles%atom_nr=-1
    nullify(local_angles%glob_angles)
    call cp_assert((present(glob_angles).and.present(atom_nr)).or.&
         present(source), cp_warning_level, cp_precondition_failed,&
         routineP,"You did not fully initialize the local_angles",&
         error,failure)
    call cp_pao_l_ang_set(local_angles, atom_nr, glob_angles, &
         source, error=error)
  end subroutine cp_pao_l_ang_init
!***************************************************************************

!!****f* cp_pao_l_angles_methods/cp_pao_l_ang_dealloc_ref [1.0] *
!!
!!   NAME
!!     cp_pao_l_ang_dealloc_ref
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_l_ang_dealloc_ref(local_angles, error)
!!       Type(cp_pao_local_angles_type), Intent (INOUT):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_l_ang_dealloc_ref
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
  subroutine cp_pao_l_ang_dealloc_ref(local_angles, error)
    real(kind=wp) ::res
    type(cp_pao_local_angles_type), intent(inout) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_l_angles_methods:cp_pao_l_ang_dealloc_ref'
    failure=.false.

    CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       nullify(local_angles%glob_angles)
    end if failureIf
  end subroutine cp_pao_l_ang_dealloc_ref
!***************************************************************************

!!****f* cp_pao_l_angles_methods/cp_pao_l_ang_validate [1.0] *
!!
!!   NAME
!!     cp_pao_l_ang_validate
!!
!!   SYNOPSIS
!!     Function cp_pao_l_ang_validate(local_angles, error) Result(res)
!!       Logical:: res
!!       Type(cp_pao_local_angles_type), Intent (IN):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_l_ang_validate
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
  function cp_pao_l_ang_validate(local_angles, error) result(res)
    logical ::res
    type(cp_pao_local_angles_type), intent(in) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_l_angles_methods:cp_pao_l_ang_validate'
    type(cp_pao_env_type), pointer :: pao_env
    integer :: i, beginAt
    failure=.false.
    nullify(pao_env)

    CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
    CPAssert(associated(local_angles%glob_angles),cp_warning_level,routineP,error,failure)
    if (.not. failure) then
       call cp_assert(cp_validate(local_angles%glob_angles,&
            error=error),cp_warning_level,cp_assertion_failed,routineP,&
            "referenced pao_glob_angles is invalid, in file "//&
            CPSourceFileRef(),error,failure)
       pao_env => local_angles%glob_angles%pao_env
! pao_env has been checked by glob angles validation
    end if
    if (.not.failure) then
       CPAssert(local_angles%atom_nr>0,cp_warning_level,routineP,error,failure)
       CPAssert(local_angles%atom_nr<=cp_get_n_atom(pao_env),cp_warning_level,routineP,error,failure)
       if (.not.failure) then
          CPAssert(local_angles%begins_at>=0,cp_warning_level,routineP,error,failure)
          call cp_assert(local_angles%begins_at < &
               size(local_angles%glob_angles%angles),&
               cp_warning_level,cp_assertion_failed,routineP,&
               "begins_at out of range in "//CPSourceFileRef(),&
               error,failure)
          beginAt=sum((pao_env%nr_full_bas(1:(local_angles%atom_nr-1))-&
               pao_env%nr_min_bas(1:(local_angles%atom_nr-1)))*&
               pao_env%nr_min_bas(1:(local_angles%atom_nr-1)))
          CPAssert(local_angles%begins_at==beginAt,cp_warning_level,routineP,error,failure)
       end if
    end if
    res=.not.failure
  end function cp_pao_l_ang_validate
!***************************************************************************


!!****f* cp_pao_l_angles_methods/cp_pao_l_ang_set [1.0] *
!!
!!   NAME
!!     cp_pao_l_ang_set
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_l_ang_set(local_angles, atom_nr, glob_angles,&
!!         source, error)
!!       Type(cp_pao_local_angles_type), Intent (INOUT):: local_angles
!!       Integer, Optional, Intent (IN):: atom_nr
!!       Type(cp_pao_glob_angles_type), Optional, Target:: glob_angles
!!       Type(cp_pao_local_angles_type), Intent (IN), Optional:: source
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_l_ang_set
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
  subroutine cp_pao_l_ang_set(local_angles, atom_nr, glob_angles,&
       source,error)
    type(cp_pao_local_angles_type), intent(inout):: local_angles
    integer, optional, intent(in) :: atom_nr
    type(cp_pao_glob_angles_type), optional, target :: glob_angles 
    type(cp_pao_local_angles_type), intent(in), optional :: source
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_l_angles_methods:cp_pao_l_ang_set'
    type(cp_pao_env_type), pointer :: pao_env
    type(cp_error_type) :: iError
    failure=.false.
    nullify(pao_env)
    call cp_init_error(iError, template_error=error)

    CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       if (present(source)) then
          CPPrecondition(source%initialized,cp_failure_level,routineP,error,failure)
          if (.not.failure) then
             local_angles%atom_nr=source%atom_nr
             local_angles%begins_at=source%begins_at
             local_angles%glob_angles => source%glob_angles
          end if
       end if
       if (present(glob_angles)) then
          local_angles%glob_angles => glob_angles
       end if
       if (present(atom_nr)) then
          local_angles%atom_nr=atom_nr
          debugIf: if (cp_debug) then
             CPAssert(atom_nr>0,cp_warning_level,routineP,error,failure)
             if (associated(local_angles%glob_angles)) then
                if (associated(local_angles%glob_angles%pao_env)) then
                   call cp_assert(atom_nr <= &
                        cp_get_n_atom(local_angles%glob_angles%pao_env),&
                        cp_warning_level,cp_assertion_failed,routineP,&
                        "atom_nr too big in "//CPSourceFileRef(),&
                        error,failure)
                end if
             end if
          end if debugIf
       end if
! updates begins_at
       if (associated(local_angles%glob_angles) .and. &
            local_angles%atom_nr>0) then
          pao_env => local_angles%glob_angles%pao_env
          if (associated(pao_env)) then
             if (cp_debug) then
                call cp_assert(cp_validate(pao_env, &
                     test_pos_att=.false., error=iError),&
                     cp_warning_level,cp_assertion_failed,routineP,&
                     "the environement should be valid when set in "//&
                     CPSourceFileRef(),error,failure)
                call cp_assert(local_angles%atom_nr <= cp_get_n_atom(pao_env),&
                     cp_warning_level,cp_assertion_failed,routineP,&
                     "the atom_nr should be either negative or valid in "//&
                     CPSourceFileRef(),error,failure)
             end if
             if (.not.failure) then
                local_angles%begins_at= &
                     sum(pao_env%nr_min_bas(1:(local_angles%atom_nr-1))*&
                     (pao_env%nr_full_bas(1:(local_angles%atom_nr-1))-&
                     pao_env%nr_min_bas(1:(local_angles%atom_nr-1))))
             end if
          end if
       end if
    end if failureIf
    call cp_dealloc_ref(iError,error=error)
  end subroutine cp_pao_l_ang_set
!***************************************************************************

!!****f* cp_pao_l_angles_methods/cp_pao_l_ang_get [1.0] *
!!
!!   NAME
!!     cp_pao_l_ang_get
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_l_ang_get(local_angles, atom_nr, n_min_bas,&
!!         n_full_bas, pao_env, n_excl_bas, unitary_t, cache_valid, error)
!!       Type(cp_pao_local_angles_type), Intent (INOUT):: local_angles
!!       Integer, Optional, Intent (OUT):: atom_nr, n_min_bas, n_full_bas,&
!!         n_excl_bas
!!       Type(cp_pao_env_type), Optional, Pointer:: pao_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Real(Kind=wp), Dimension(:,:), Pointer, Optional:: unitary_t
!!       Logical, Intent (OUT), Optional:: cache_valid
!!     End Subroutine cp_pao_l_ang_get
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
!!     pao_env: the pao environement this local angles lives in.
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
  subroutine cp_pao_l_ang_get(local_angles, atom_nr, n_min_bas, &
       n_full_bas, pao_env, n_excl_bas,unitary_t,cache_valid, error)
    type(cp_pao_local_angles_type), intent(inout) :: local_angles
    integer, optional, intent(out) :: atom_nr, n_min_bas, n_full_bas,&
         n_excl_bas
    type(cp_pao_env_type), optional, pointer :: pao_env
    type(cp_error_type), optional, intent(inout) :: error
    real(kind=wp), dimension(:,:), pointer, optional :: unitary_t
    logical, intent(out), optional :: cache_valid

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_l_angles_methods:cp_pao_l_ang_get'
    type(cp_pao_env_type), pointer :: my_pao_env
    failure=.false.
    nullify(my_pao_env)

    CPAssert(local_angles%initialized,cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       if (present(atom_nr)) then
          atom_nr=local_angles%atom_nr
       end if
       if (associated(local_angles%glob_angles)) then
          my_pao_env => local_angles%glob_angles%pao_env
       end if
       if (present(pao_env)) then
          pao_env => my_pao_env
       end if
       if (present(n_min_bas)) then
          CPAssert(local_angles%atom_nr>0,cp_failure_level,routineP,error,failure)
! also check upper bound?
          n_min_bas=my_pao_env%nr_min_bas(local_angles%atom_nr)
       end if
       if (present(n_full_bas)) then
          CPAssert(local_angles%atom_nr>0,cp_failure_level,routineP,error,failure)
          n_full_bas=my_pao_env%nr_full_bas(local_angles%atom_nr)
       end if
       if (present(n_excl_bas)) then
          CPAssert(local_angles%atom_nr>0,cp_failure_level,routineP,error,failure)
          n_excl_bas=my_pao_env%nr_full_bas(local_angles%atom_nr)-&
               my_pao_env%nr_min_bas(local_angles%atom_nr)
       end if
       if (present(unitary_t)) then
          unitary_t => cp_pao_l_ang_get_unitary_t(local_angles,error=error)
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
       nullify(my_pao_env)
    end if failureIf
  end subroutine cp_pao_l_ang_get
!***************************************************************************

!!****f* cp_pao_l_angles_methods/cp_pao_l_angles_valid [1.0] *
!!
!!   NAME
!!     cp_pao_l_angles_valid
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
  function cp_pao_l_angles_valid(local_angles, error) result(res)
    logical ::res
    type(cp_pao_local_angles_type), intent(in) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_l_angles_methods:cp_pao_l_angles_valid'
    failure=.false.

    CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
    CPAssert(associated(local_angles%glob_angles),cp_warning_level,routineP,error,failure)
    if (cp_debug.and.debug_this_module) then
       if (.not. failure) then
          if (.not.associated(local_angles%glob_angles%pao_env)) then
             CPAssert(.false.,cp_warning_level,routineP,error,failure)
          else if (.not.(associated(local_angles%glob_angles%pao_env%nr_full_bas)&
               .and.associated(local_angles%glob_angles%pao_env%nr_min_bas))) then
             failure=cp_valid(local_angles%glob_angles%pao_env, error=error)
          end if
       end if
    end if
    res=.not.failure
  end function cp_pao_l_angles_valid
!***************************************************************************

! ====== special methods =====

!!****f* cp_pao_l_angles_methods/cp_pao_l_ang_get_pao_env [1.0] *
!!
!!   NAME
!!     cp_pao_l_ang_get_pao_env
!!
!!   SYNOPSIS
!!     Function cp_pao_l_ang_get_pao_env(local_angles, error) Result(res)
!!       Type(cp_pao_env_type), Pointer:: res
!!       Type(cp_pao_local_angles_type), Intent (IN):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_l_ang_get_pao_env
!!
!!   FUNCTION
!!     returns the global environement the local angles object lives in
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles object you are interested in
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
  function cp_pao_l_ang_get_pao_env(local_angles, error) result(res)
    type(cp_pao_env_type), pointer ::res
    type(cp_pao_local_angles_type), intent(in) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_l_angles_methods:cp_pao_l_ang_get_pao_env'
    failure=.false.
    CPPrecondition(local_angles%initialized, cp_failure_level,routineP,error,failure)
    CPPrecondition(associated(local_angles%glob_angles),cp_warning_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       res => local_angles%glob_angles%pao_env
    else failureIf
       nullify(res)
    end if failureIf
  end function cp_pao_l_ang_get_pao_env
!***************************************************************************

!!****f* cp_pao_l_angles_methods/cp_pao_l_ang_get_p_index [1.0] *
!!
!!   NAME
!!     cp_pao_l_ang_get_p_index
!!
!!   SYNOPSIS
!!     Function cp_pao_l_ang_get_p_index(local_angles, i, error)&
!!         Result(res)
!!       Integer:: res
!!       Type(cp_pao_local_angles_type), Intent (IN):: local_angles
!!       Integer, Intent (IN):: i
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_l_ang_get_p_index
!!
!!   FUNCTION
!!     returns the i-te index of the projection
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
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function cp_pao_l_ang_get_p_index(local_angles, i,error) result(res)
    integer ::res
    type(cp_pao_local_angles_type), intent(in) :: local_angles
    integer, intent(in) :: i
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_l_angles_methods:cp_pao_l_ang_get_p_index'
    type(cp_pao_env_type), pointer :: pao_env
    failure=.false.
    pao_env => cp_get_pao_env(local_angles)
    if (cp_debug .or. debug_this_module) then
       CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(pao_env),cp_failure_level,routineP,error,failure)
    end if

    failureIf: if (.not. failure) then
       res=cp_get_p_index(pao_env,atom_nr=local_angles%atom_nr,&
            i=i,error=error)
    else failureIf
       res=-1
    end if failureIf
  end function cp_pao_l_ang_get_p_index
!***************************************************************************

!!****f* cp_pao_l_angles_methods/cp_pao_l_ang_set_p_index [1.0] *
!!
!!   NAME
!!     cp_pao_l_ang_set_p_index
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_l_ang_set_p_index(local_angles, i, value, error)
!!       Type(cp_pao_local_angles_type), Intent (INOUT):: local_angles
!!       Integer, Intent (IN):: i, value
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_l_ang_set_p_index
!!
!!   FUNCTION
!!     set the value of the i-te index of the projection
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles you set values in
!!     i: number of the index you want to set
!!     value: the new nalue of the i-te index
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
  subroutine cp_pao_l_ang_set_p_index(local_angles,i,value,error)
    type(cp_pao_local_angles_type), intent(inout) :: local_angles ! actually only pao_env is changed
    integer, intent(in) :: i,value
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_l_angles_methods:cp_pao_l_ang_set_p_index'
    type(cp_pao_env_type), pointer :: pao_env
    failure=.false.
    pao_env => cp_get_pao_env(local_angles)
    if (cp_debug .or. debug_this_module) then
       CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(pao_env),cp_failure_level,routineP,error,failure)
    end if

    failureIf: if (.not. failure) then
       call cp_set_p_index(pao_env,atom_nr=local_angles%atom_nr,&
            i=i,value=value,error=error)
    end if failureIf
  end subroutine cp_pao_l_ang_set_p_index
!***************************************************************************

!!****f* cp_pao_l_angles_methods/cp_pao_l_ang_get_non_p_index [1.0] *
!!
!!   NAME
!!     cp_pao_l_ang_get_non_p_index
!!
!!   SYNOPSIS
!!     Function cp_pao_l_ang_get_non_p_index(local_angles, i, error)&
!!         Result(res)
!!       Integer:: res
!!       Type(cp_pao_local_angles_type), Intent (IN):: local_angles
!!       Integer, Intent (IN):: i
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_l_ang_get_non_p_index
!!
!!   FUNCTION
!!     returns the i-te index of the excluded basis
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
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function cp_pao_l_ang_get_non_p_index(local_angles, i,error) result(res)
    integer ::res
    type(cp_pao_local_angles_type), intent(in) :: local_angles
    integer, intent(in) :: i
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_l_angles_methods:cp_pao_l_ang_get_non_p_index'
    type(cp_pao_env_type), pointer :: pao_env
    failure=.false.
    pao_env => cp_get_pao_env(local_angles)
    if (cp_debug .or. debug_this_module) then
       CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(pao_env),cp_failure_level,routineP,error,failure)
    end if

    failureIf: if (.not. failure) then
       res=cp_get_non_p_index(pao_env,atom_nr=local_angles%atom_nr,&
            i=i,error=error)
    else failureIf
       res=-1
    end if failureIf
  end function cp_pao_l_ang_get_non_p_index
!***************************************************************************

!!****f* cp_pao_l_angles_methods/cp_pao_l_ang_get_angle [1.0] *
!!
!!   NAME
!!     cp_pao_l_ang_get_angle
!!
!!   SYNOPSIS
!!     Function cp_pao_l_ang_get_angle(local_angles, min_bas, excl_bas,&
!!         error) Result(res)
!!       Real(Kind=wp):: res
!!       Type(cp_pao_local_angles_type), Intent (IN):: local_angles
!!       Integer, Intent (IN):: min_bas, excl_bas
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_l_ang_get_angle
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
  function cp_pao_l_ang_get_angle(local_angles, min_bas, excl_bas, &
       error) result(res)
    real(kind=wp) ::res
    type(cp_pao_local_angles_type), intent(in) :: local_angles
    integer, intent(in) :: min_bas, excl_bas
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_l_angles_methods:cp_pao_l_ang_get_angle'
    type(cp_pao_glob_angles_type), pointer :: glob_angles
    failure=.false.
    glob_angles => local_angles%glob_angles

    debugIf: if (debug_this_module.and.cp_debug) then
       CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(glob_angles),cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          CPPrecondition(associated(glob_angles%pao_env),cp_failure_level,routineP,error,failure)
       end if
       if (.not.failure) then
          call cp_assert(associated(glob_angles%pao_env%nr_min_bas),&
               cp_failure_level,cp_assertion_failed,routineP,&
               "min bas not initialized in "//CPSourceFileRef(),&
               error,failure)
          call cp_assert(associated(glob_angles%pao_env%nr_full_bas),&
               cp_failure_level,cp_assertion_failed,routineP,&
               "full bas not initialized in "//CPSourceFileRef(),&
               error,failure)
       end if
       if (.not.failure) then
          call cp_assert(min_bas>0.and.min_bas <= &
               glob_angles%pao_env%nr_min_bas(local_angles%atom_nr),&
               cp_failure_level,cp_assertion_failed,routineP,&
               "min_bas out of range in "//CPSourceFileRef(),&
               error,failure)
          call cp_assert(excl_bas>0.and.excl_bas <= &
               glob_angles%pao_env%nr_full_bas(local_angles%atom_nr)-&
               glob_angles%pao_env%nr_min_bas(local_angles%atom_nr),&
               cp_failure_level,cp_assertion_failed,routineP,&
               "excl_bas out of range in "//CPSourceFileRef(),&
               error,failure)
       end if
! also check begins_at???
    end if debugIf
    failureIf: if (.not. failure) then
       res=glob_angles%angles(local_angles%begins_at + &
            glob_angles%pao_env%nr_min_bas(local_angles%atom_nr)*(min_bas-1) + &
            excl_bas)
    else failureIf
       res=-1
    end if failureIf
  end function cp_pao_l_ang_get_angle
!***************************************************************************

!!****f* cp_pao_l_angles_methods/cp_pao_l_ang_set_angle [1.0] *
!!
!!   NAME
!!     cp_pao_l_ang_set_angle
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_l_ang_set_angle(local_angles, min_bas, excl_bas,&
!!         value, error)
!!       Type(cp_pao_local_angles_type), Intent (IN):: local_angles
!!       Integer, Intent (IN):: min_bas, excl_bas
!!       Real(Kind=wp), Intent (IN):: value
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_l_ang_set_angle
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
  subroutine cp_pao_l_ang_set_angle(local_angles, min_bas, excl_bas, &
       value,error)
    type(cp_pao_local_angles_type), intent(in) :: local_angles !actually the substructure globalangles is modifed
    integer, intent(in) :: min_bas, excl_bas
    real(kind=wp), intent(in) :: value
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_l_angles_methods:cp_pao_l_ang_set_angle'
    type(cp_pao_glob_angles_type), pointer :: glob_angles
    failure=.false.
    glob_angles => local_angles%glob_angles

    debugIf: if (debug_this_module.and.cp_debug) then
       CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(glob_angles),cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          CPPrecondition(associated(glob_angles%pao_env),cp_failure_level,routineP,error,failure)
       end if
       if (.not.failure) then
          call cp_assert(associated(local_angles%glob_angles%pao_env%nr_min_bas),&
               cp_failure_level,cp_assertion_failed,routineP,&
               "env min bas not initialized in "//CPSourceFileRef(),&
               error,failure)
          call cp_assert(associated(local_angles%glob_angles%pao_env%nr_full_bas),&
               cp_failure_level,cp_assertion_failed,routineP,&
               "env full bas not initialized in "//CPSourceFileRef(),&
               error,failure)
       end if
       if (.not.failure) then
          call cp_assert(min_bas>0.and.min_bas <= &
               glob_angles%pao_env%nr_min_bas(local_angles%atom_nr),&
               cp_failure_level,cp_assertion_failed,routineP,&
               "min_bas out of range in "//CPSourceFileRef(),&
               error,failure)
          call cp_assert(excl_bas>0.and.excl_bas <= &
               glob_angles%pao_env%nr_full_bas(local_angles%atom_nr)-&
               glob_angles%pao_env%nr_min_bas(local_angles%atom_nr),&
               cp_failure_level,cp_assertion_failed,routineP,&
               "excl_bas out of range in "//CPSourceFileRef(),&
               error,failure)
       end if
! also check begins_at???
    end if debugIf
    failureIf: if (.not. failure) then
       glob_angles%angles(local_angles%begins_at + &
            glob_angles%pao_env%nr_min_bas(local_angles%atom_nr)*(min_bas-1) + &
            excl_bas) = value
    end if failureIf
    nullify(glob_angles)
  end subroutine cp_pao_l_ang_set_angle
!***************************************************************************

!!****f* cp_pao_l_angles_methods/cp_pao_l_ang_get_excl_bas_size [1.0] *
!!
!!   NAME
!!     cp_pao_l_ang_get_excl_bas_size
!!
!!   SYNOPSIS
!!     Function cp_pao_l_ang_get_excl_bas_size(local_angles, error)&
!!         Result(res)
!!       Integer:: res
!!       Type(cp_pao_local_angles_type), Intent (IN):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_l_ang_get_excl_bas_size
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
  function cp_pao_l_ang_get_excl_bas_size(local_angles, error) result(res)
    integer ::res
    type(cp_pao_local_angles_type), intent(in) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_l_angles_methods:cp_pao_l_ang_get_excl_bas_size'
    type(cp_pao_env_type), pointer :: pao_env
    failure=.false.

    CPPrecondition(local_angles%initialized, cp_failure_level,routineP,error,failure)
    pao_env => cp_get_pao_env(local_angles)
    CPPrecondition(associated(pao_env),cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       res=cp_pao_env_get_excl_bas_size(pao_env,atom_nr=local_angles%atom_nr,&
            error=error)
    else
       res=-1
    end if
  end function cp_pao_l_ang_get_excl_bas_size
!***************************************************************************

!!****f* cp_pao_l_angles_methods/cp_pao_l_ang_get_min_bas_size [1.0] *
!!
!!   NAME
!!     cp_pao_l_ang_get_min_bas_size
!!
!!   SYNOPSIS
!!     Function cp_pao_l_ang_get_min_bas_size(local_angles, error)&
!!         Result(res)
!!       Integer:: res
!!       Type(cp_pao_local_angles_type), Intent (IN):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_l_ang_get_min_bas_size
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
  function cp_pao_l_ang_get_min_bas_size(local_angles, error) result(res)
    integer ::res
    type(cp_pao_local_angles_type), intent(in) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_l_angles_methods:cp_pao_l_ang_get_min_bas_size'
    type(cp_pao_env_type), pointer :: pao_env
    failure=.false.

    CPPrecondition(local_angles%initialized, cp_failure_level,routineP,error,failure)
    pao_env => cp_get_pao_env(local_angles)
    CPPrecondition(associated(pao_env),cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       res=cp_get_min_bas_size(pao_env,atom_nr=local_angles%atom_nr,&
            error=error)
    else
       res=-1
    end if
  end function cp_pao_l_ang_get_min_bas_size
!***************************************************************************

!!****f* cp_pao_l_angles_methods/cp_pao_l_ang_get_full_bas_size [1.0] *
!!
!!   NAME
!!     cp_pao_l_ang_get_full_bas_size
!!
!!   SYNOPSIS
!!     Function cp_pao_l_ang_get_full_bas_size(local_angles, error)&
!!         Result(res)
!!       Integer:: res
!!       Type(cp_pao_local_angles_type), Intent (IN):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_l_ang_get_full_bas_size
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
  function cp_pao_l_ang_get_full_bas_size(local_angles, error) result(res)
    integer ::res
    type(cp_pao_local_angles_type), intent(in) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_l_angles_methods:cp_pao_l_ang_get_full_bas_size'
    type(cp_pao_env_type), pointer :: pao_env
    failure=.false.

    CPPrecondition(local_angles%initialized, cp_failure_level,routineP,error,failure)
    pao_env => cp_get_pao_env(local_angles)
    CPPrecondition(associated(pao_env),cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       res=cp_get_full_bas_size(pao_env,atom_nr=local_angles%atom_nr,&
            error=error)
    else
       res=-1
    end if
  end function cp_pao_l_ang_get_full_bas_size
!***************************************************************************

!!****f* cp_pao_l_angles_methods/cp_pao_l_ang_get_glob_angles [1.0] *
!!
!!   NAME
!!     cp_pao_l_ang_get_glob_angles
!!
!!   SYNOPSIS
!!     Function cp_pao_l_ang_get_glob_angles(local_angles, error)&
!!         Result(res)
!!       Type(cp_pao_local_angles_type), Intent (IN):: local_angles
!!       Type(cp_pao_glob_angles_type), Pointer:: res
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_l_ang_get_glob_angles
!!
!!   FUNCTION
!!     returns the global angles object that this object accesses
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     local_angles: the local angles you get the info from
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
  function cp_pao_l_ang_get_glob_angles(local_angles,error) result(res)
    type(cp_pao_local_angles_type), intent(in) :: local_angles
    type(cp_pao_glob_angles_type), pointer ::res
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_l_angles_methods:cp_pao_l_ang_get_glob_angles'
    failure=.false.

    CPPrecondition(cp_valid(local_angles),cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       res => local_angles%glob_angles
    else
       nullify(res)
    end if
  end function cp_pao_l_ang_get_glob_angles
!***************************************************************************


!!****f* cp_pao_l_angles_methods/cp_pao_l_ang_did_change [1.0] *
!!
!!   NAME
!!     cp_pao_l_ang_did_change
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_l_ang_did_change(local_angles, error)
!!       Type(cp_pao_local_angles_type), Intent (INOUT):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_l_ang_did_change
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
  subroutine cp_pao_l_ang_did_change(local_angles,error)
! changes in local_angles%glob_angles
    type(cp_pao_local_angles_type), intent(inout) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_l_angles_methods:cp_pao_l_ang_did_change'
    logical, dimension(:), pointer :: ut_calculated
    failure=.false.

    CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(local_angles%atom_nr>0,cp_warning_level,routineP,error,failure)
    if (.not. failure) then
       if (associated(local_angles%glob_angles%ut_calculated)) then
          local_angles%glob_angles%ut_calculated(local_angles%atom_nr)=.false.
       end if
    end if
  end subroutine cp_pao_l_ang_did_change
!***************************************************************************

!!****f* cp_pao_l_angles_methods/cp_pao_l_ang_set_unitary_t [1.0] *
!!
!!   NAME
!!     cp_pao_l_ang_set_unitary_t
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_l_ang_set_unitary_t(local_angles, unitary_t,&
!!         error)
!!       Type(cp_pao_local_angles_type), Intent (INOUT):: local_angles
!!       Real(Kind=wp), Dimension(:,:), Intent (IN), Target:: unitary_t
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_l_ang_set_unitary_t
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
  subroutine cp_pao_l_ang_set_unitary_t(local_angles,unitary_t,error)
! changes in glob_angles
    type(cp_pao_local_angles_type), intent(inout) :: local_angles
    real(kind=wp), dimension(:,:), intent(in), target :: unitary_t
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_l_angles_methods:cp_pao_l_ang_set_unitary_t'
    type(cp_pao_matrix_type), pointer :: full_u_t
    real(kind=wp), dimension(:,:), pointer :: unitary_t_ptr
    failure=.false.

    CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(local_angles%atom_nr>0,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       call cp_get(local_angles%glob_angles,unitary_t=full_u_t,&
            error=error)
       CPPrecondition(associated(full_u_t),cp_warning_level,routineP,error,failure)
       call cp_assert(associated(local_angles%glob_angles%ut_calculated),&
            cp_failure_level,cp_assertion_failed,routineP,&
            "ut_calculated not allocated in "//&
            CPSourceFileRef(),error,failure)
       if (.not.failure) then
          local_angles%glob_angles%ut_calculated(local_angles%atom_nr)=.true.
          unitary_t_ptr => unitary_t
          call put_block_node(cp_get_matrix(full_u_t),&
               block_row=local_angles%atom_nr,&
               block_col=local_angles%atom_nr, block=unitary_t_ptr)
       end if
    end if
  end subroutine cp_pao_l_ang_set_unitary_t
!***************************************************************************

!!****f* cp_pao_plain_rot/cp_pao_l_compose_u [1.0] *
!!
!!   NAME
!!     cp_pao_l_compose_u
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_l_compose_u(matrix, local_angles, error)
!!       Real(Kind=wp), Dimension(:,:), Intent (INOUT):: matrix
!!       Type(cp_pao_local_angles_type), Intent (INOUT):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_l_compose_u
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
  subroutine cp_pao_l_compose_u(matrix,local_angles,error)
    real(kind=wp), dimension(:,:), intent(inout) :: matrix
    type(cp_pao_local_angles_type), intent(inout) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_plain_rot:cp_pao_l_compose_u'
    integer :: n_min_bas, n_excl_bas, i_min_bas,i_excl_bas,n_full_bas
    failure=.false.

    call cp_get(local_angles,n_excl_bas=n_excl_bas,n_min_bas=n_min_bas,&
         n_full_bas=n_full_bas)
    CPPrecondition(n_min_bas>0,cp_warning_level,routineP,error,failure)
    CPPrecondition(n_excl_bas>=0,cp_warning_level,routineP,error,failure)
    CPPrecondition(size(matrix,1)==n_full_bas,cp_warning_level,routineP,error,failure)
    CPPrecondition(size(matrix,2)==n_full_bas,cp_warning_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       do i_excl_bas=n_excl_bas,1,-1
          do i_min_bas=n_min_bas,1,-1
             call cp_left_compose_p_rot(matrix=matrix,&
                  angle=cp_get_angle(local_angles,min_bas=i_min_bas,&
                  excl_bas=i_excl_bas,error=error),&
                  i=cp_get_p_index(local_angles,i=i_min_bas,error=error),&
                  j=cp_get_non_p_index(local_angles,i=i_excl_bas,error=error),&
                  error=error)
          end do
       end do
    end if failureIf
  end subroutine cp_pao_l_compose_u
!***************************************************************************

!!****f* cp_pao_l_angles_methods/cp_pao_l_ang_get_unitary_t [1.0] *
!!
!!   NAME
!!     cp_pao_l_ang_get_unitary_t
!!
!!   SYNOPSIS
!!     Function cp_pao_l_ang_get_unitary_t(local_angles, error) Result(res)
!!       Type(cp_pao_local_angles_type), Intent (INOUT):: local_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Real(Kind=wp), Dimension(:,:), Pointer:: res
!!     End Function cp_pao_l_ang_get_unitary_t
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
  function cp_pao_l_ang_get_unitary_t(local_angles,error) result(res)
    type(cp_pao_local_angles_type), intent(inout) :: local_angles
    type(cp_error_type), optional, intent(inout) :: error
    real(kind=wp), dimension(:,:), pointer :: res

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_l_angles_methods:cp_pao_l_ang_get_unitary_t'
    integer :: i
    type(cp_pao_matrix_type), pointer :: unitary_t
    failure=.false.
    nullify(res)

    CPPrecondition(cp_valid(local_angles),cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       call cp_get(local_angles%glob_angles,unitary_t=unitary_t,&
            error=error)
       CPPrecondition(associated(unitary_t),cp_failure_level,routineP,error,failure)
       call cp_assert(associated(local_angles%glob_angles%ut_calculated),&
            cp_failure_level,cp_assertion_failed,routineP,&
            "ut_calculated not associated in "//&
            CPSourceFileRef(),error,failure)
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
                call cp_pao_l_compose_u(res,local_angles,error=error)
                local_angles%glob_angles%ut_calculated(local_angles%atom_nr)=.true.
             end if
          end if
       end if
    end if
    if (failure) nullify(res)
  end function cp_pao_l_ang_get_unitary_t
!***************************************************************************

!!****f* cp_pao_l_angles_methods/cp_pao_l_ang_get_g_index [1.0] *
!!
!!   NAME
!!     cp_pao_l_ang_get_g_index
!!
!!   SYNOPSIS
!!     Function cp_pao_l_ang_get_g_index(local_angles, min_bas, excl_bas,&
!!         error) Result(res)
!!       Integer:: res
!!       Type(cp_pao_local_angles_type), Intent (INOUT):: local_angles
!!       Integer, Intent (IN):: min_bas, excl_bas
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_l_ang_get_g_index
!!
!!   FUNCTION
!!     returns the global angles index for the angle described with
!!     the given local indices. 
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
function cp_pao_l_ang_get_g_index(local_angles, min_bas, excl_bas, error)&
     result(res)
  integer :: res
  type(cp_pao_local_angles_type), intent(inout) :: local_angles
  integer, intent(in) :: min_bas, excl_bas
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  integer :: handle
  character(len=*), parameter :: routineN='cp_pao_l_ang_get_g_index',&
	routineP=moduleN//':'//routineN
  failure=.false.
  
  !call timeset(routineN//','//moduleN,'I',"",handle)
  if (cp_debug .and. debug_this_module) then
     CPPrecondition(local_angles%initialized,cp_failure_level,routineP,error,failure)
     CPPrecondition(local_angles%atom_nr>0,cp_failure_level,routineP,error,failure)
  end if
  if (.not. failure) then
     res = local_angles%glob_angles% pao_env% angles_begin_at&
          (local_angles%atom_nr)+&
          local_angles%glob_angles%pao_env%nr_min_bas(local_angles%atom_nr)*&
          (min_bas-1) + excl_bas
  else
     res=-1
  end if
  !call timestop(0.0_wp,handle)
end function cp_pao_l_ang_get_g_index
!***************************************************************************

end module cp_pao_l_angles_methods
