!----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations        !
!   Copyright (C) 2001  CP2K developers group                                !
!----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_glob_angles_methods [1.0] *
!!
!!   NAME
!!     pao_glob_angles_methods
!!
!!   FUNCTION
!!     metodhs to access and modify the angles that define the unitary
!!     transformation
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!****************************************************************************
module pao_glob_angles_methods
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use pao_types
  use cp_sparse_matrix
  use cp_b_matrix_structure
  use pao_utils
  use pao_qs_env_methods
  use qs_environment_types
  use cp_array_utils!, only: cp_guarantee_size
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='pao_glob_angles_methods'

! core procedures
  public :: cp_init, cp_dealloc_ref, cp_set, cp_get, cp_create, cp_dealloc,&
       cp_did_change

! core functions
  public :: cp_validate, cp_valid

! the underlying functions
  public :: pao_g_ang_init, pao_g_ang_dealloc_ref, &
       pao_g_ang_create, pao_g_ang_dealloc,&
       pao_g_ang_validate, pao_g_ang_valid,&
       pao_g_ang_set, pao_g_ang_get, pao_g_ang_did_change

! ============== interfaces ===========

! +++ begin of core methods +++

!! initializes the first, argument.
!! Support optional initial values that depend on the 
!! first argument 
  interface cp_init
     module procedure pao_g_ang_init
  end interface

!! Deallocates the memory allocated by the first argument.
!! Supports error=error
  interface cp_dealloc_ref
     module procedure pao_g_ang_dealloc_ref
  end interface

!! Allocates and initialize a new object returning a pointer to it in 
!! the first argument
  interface cp_create
     module procedure pao_g_ang_create
  end interface

!! Deallocates the memory used by the first argument, and the object
!! itself
  interface cp_dealloc
     module procedure pao_g_ang_dealloc
  end interface

!! returns true if the content of the first argumenst is valid.
!! writes the errors ar warnings, full validation (long)
  interface cp_validate
     module procedure pao_g_ang_validate
  end interface

!! returns true if the first argument is valid. Only a fast minimal 
!! validation is performed (no dangling pointers)
  interface cp_valid
     module procedure pao_g_ang_valid
  end interface

!! sets the value of various optonal attributes of the firs argument.
!! the optional attributes depend on the type of the first
  interface cp_set
     module procedure pao_g_ang_set
  end interface

!! gets values from the first argument via optional keyword arguments
!! that depend on the type of the first argument
  interface cp_get
     module procedure pao_g_ang_get
  end interface

!! tells the object that its internal values have changes and its cached
!! values may be invalid
  interface cp_did_change
     module procedure pao_g_ang_did_change
  end interface

! +++ end of the core methods +++

!!***
!****************************************************************************

contains

! ======= core methods ======

!!****f* pao_glob_angles_methods/pao_g_ang_init [1.0] *
!!
!!   NAME
!!     pao_g_ang_init
!!
!!   SYNOPSIS
!!     Subroutine pao_g_ang_init(glob_angles, pao_env, angles, source,&
!!         ut_calculated, unitary_t, min_overlap_m, min_hamiltonian_m,&
!!         error)
!!       Type(pao_glob_angles_type), Intent (OUT):: glob_angles
!!       Type(pao_env_type), Intent (IN), Target:: pao_env
!!       Real(Kind=wp), Optional, Dimension(:), Intent (IN):: angles
!!       Type(pao_glob_angles_type), Intent (IN), Optional:: source
!!       Logical, Dimension(:), Intent (IN), Optional:: ut_calculated
!!       Type(pao_matrix_type), Intent (IN), Optional:: unitary_t,&
!!         min_hamiltonian_m, min_overlap_m
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_g_ang_init
!!
!!   FUNCTION
!!     initializes the global angles
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: the global angles to be initialized
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
  subroutine pao_g_ang_init(glob_angles, qs_env, global_env, angles, source,&
       ut_calculated, unitary_t, min_overlap_m, min_hamiltonian_m, error)
    type(pao_glob_angles_type), intent(out) :: glob_angles
    type(qs_environment_type), intent(in), target :: qs_env
    type(global_environment_type), intent(in), target :: global_env
    real(kind=wp), optional, dimension(:), intent(in) :: angles
    type(pao_glob_angles_type), intent(in), optional :: source
    logical, dimension(:), intent(in), optional :: ut_calculated
    type(cp_sparse_matrix_type), intent(in), optional, target :: unitary_t,&
         min_overlap_m
    type(cp_sparse_matrix_p_type), intent(in), dimension(:), optional ::&
         min_hamiltonian_m
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_glob_angles_methods:pao_g_ang_init'
    failure=.false.

    nullify(glob_angles%angles, glob_angles%unitary_t,&
         glob_angles%ut_calculated, glob_angles%min_overlap_m, &
         glob_angles%min_hamiltonian_m)
    glob_angles%m_s_m_valid=.false.
    glob_angles%m_h_m_valid=.false.
    glob_angles%initialized = .true.
    if (.not.failure) then
       call pao_g_ang_set(glob_angles,qs_env=qs_env, global_env=global_env,&
            angles=angles,&
            source=source, ut_calculated=ut_calculated,&
            unitary_t=unitary_t,min_hamiltonian_m=min_hamiltonian_m,&
            min_overlap_m=min_overlap_m, error=error)
    end if
  end subroutine pao_g_ang_init
!***************************************************************************

!!****f* pao_glob_angles_methods/pao_g_ang_dealloc_ref [1.0] *
!!
!!   NAME
!!     pao_g_ang_dealloc_ref
!!
!!   SYNOPSIS
!!     Subroutine pao_g_ang_dealloc_ref(glob_angles, error)
!!       Type(pao_glob_angles_type), Intent (INOUT):: glob_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_g_ang_dealloc_ref
!!
!!   FUNCTION
!!     deallocates the memory allocated by a glob_angles object.
!!     does NOT deallocate the object itself
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_glob_angles: the object to release
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
  subroutine pao_g_ang_dealloc_ref(glob_angles,global_env, error)
    type(pao_glob_angles_type), intent(inout) :: glob_angles
    type(global_environment_type), intent(in), target :: global_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_glob_angles_methods:pao_g_ang_dealloc_ref'
    integer ::stat
    failure=.false.

    CPPrecondition(glob_angles%initialized, cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       if (associated(glob_angles%angles)) then
          deallocate(glob_angles%angles, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       end if
       if (associated(glob_angles%ut_calculated)) then
          deallocate(glob_angles%ut_calculated,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       end if
       call cp_dealloc(glob_angles%unitary_t,error=error)
       call cp_dealloc(glob_angles%min_hamiltonian_m,error=error)
       call cp_dealloc(glob_angles%min_overlap_m,error=error)
    end if failureIf
  end subroutine pao_g_ang_dealloc_ref
!***************************************************************************

!!****f* pao_glob_angles_methods/pao_g_ang_dealloc [1.0] *
!!
!!   NAME
!!     pao_g_ang_dealloc
!!
!!   SYNOPSIS
!!     Subroutine pao_g_ang_dealloc(glob_angles, error)
!!       Type(pao_glob_angles_type), Pointer:: glob_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_g_ang_dealloc
!!
!!   FUNCTION
!!     deallocated the memory allocated by this object and the oject itself
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: the object to be deallocated
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
  subroutine pao_g_ang_dealloc(glob_angles,global_env,error)
    type(pao_glob_angles_type), pointer :: glob_angles
    type(global_environment_type), intent(in), target :: global_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_glob_angles_methods:pao_g_ang_dealloc'
    integer :: stat
    failure=.false.

    if (associated(glob_angles)) then
       CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
       if (.not. failure) then
          call pao_g_ang_dealloc_ref(glob_angles, global_env=global_env,&
               error=error)
       end if
       deallocate(glob_angles,stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    end if
  end subroutine pao_g_ang_dealloc
!***************************************************************************

!!****f* pao_glob_angles_methods/pao_g_ang_create [1.0] *
!!
!!   NAME
!!     pao_g_ang_create
!!
!!   SYNOPSIS
!!     Subroutine pao_g_ang_create(glob_angles, pao_env, angles, source,&
!!         ut_calculated, unitary_t, min_overlap_m, min_hamiltonian_m,&
!!         error)
!!       Type(pao_glob_angles_type), Pointer:: glob_angles
!!       Type(pao_env_type), Intent (IN), Target:: pao_env
!!       Real(Kind=wp), Optional, Dimension(:), Intent (IN):: angles
!!       Type(pao_glob_angles_type), Intent (IN), Optional:: source
!!       Logical, Dimension(:), Intent (IN):: ut_calculated
!!       Type(pao_matrix_type), Intent (IN):: unitary_t,&
!!         min_hamiltonian_m, min_overlap_m
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_g_ang_create
!!
!!   FUNCTION
!!     creates a new glob_angles object
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: will be a pointer to the newly created object
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!     for the others inputs see pao_g_ang_init
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine pao_g_ang_create(glob_angles, qs_env, global_env, angles, source,&
       ut_calculated, unitary_t, min_overlap_m, min_hamiltonian_m, error)
    type(pao_glob_angles_type), pointer :: glob_angles
    type(qs_environment_type), intent(in), target :: qs_env
    type(global_environment_type), intent(in), target :: global_env
    real(kind=wp), optional, dimension(:), intent(in) :: angles
    type(pao_glob_angles_type), intent(in), optional :: source
    logical, dimension(:), intent(in) :: ut_calculated
    type(cp_sparse_matrix_type), intent(in), target :: unitary_t,&
         min_overlap_m
    type(cp_sparse_matrix_p_type), intent(in), dimension(:), optional ::&
         min_hamiltonian_m
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_glob_angles_methods:pao_g_ang_create'
    integer :: stat
    failure=.false.

    allocate(glob_angles, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       call pao_g_ang_init(glob_angles, qs_env, global_env, angles, source,&
            ut_calculated, unitary_t, min_overlap_m, min_hamiltonian_m, error)
    end if
  end subroutine pao_g_ang_create
!***************************************************************************

!!****f* pao_glob_angles_methods/pao_g_ang_validate [1.0] *
!!
!!   NAME
!!     pao_g_ang_validate
!!
!!   SYNOPSIS
!!     Function pao_g_ang_validate(glob_angles, error) Result(res)
!!       Logical:: res
!!       Type(pao_glob_angles_type), Intent (IN):: glob_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_g_ang_validate
!!
!!   FUNCTION
!!     returns true if the angles are valid
!!
!!   NOTES
!!     Writes out all the failures as warning, use error to change
!!     this behaviour
!!
!!   INPUTS
!!     glob_angles: the gangles to be tested
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
  function pao_g_ang_validate(glob_angles,qs_env,global_env, error) result(res)
    logical ::res
    type(pao_glob_angles_type), intent(in) :: glob_angles
    type(qs_environment_type), intent(in),target :: qs_env
    type(global_environment_type), intent(in), target :: global_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_glob_angles_methods:pao_g_ang_validate'
    type(pao_env_type), pointer :: pao_env
    failure=.false.
    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       CPAssert(associated(glob_angles%angles),cp_warning_level,routineP,error,failure)
       call get_qs_env(qs_env,pao_env=pao_env)
       CPAssert(associated(pao_env),cp_warning_level,routineP,error,failure)
       if (.not.failure) then
          call cp_assert(pao_env%angles_begin_at &
               (size(pao_env%angles_begin_at))-1== &
               size(glob_angles%angles),cp_warning_level,&
               cp_assertion_failed, routineP,&
               "number of angles different from last begin_at in "//&
               CPSourceFileRef,&
               error,failure)
       end if
!!FM       if (associated(glob_angles%ut_calculated)) then
!!FM          call cp_assert(size(glob_angles%ut_calculated)==&
!!FM               glob_angles%pao_env%n_atoms, cp_warning_level,&
!!FM               cp_assertion_failed, routineP,&
!!FM               "cached invalid ut_calculated in "//&
!!FM               CPSourceFileRef,&
!!FM               error,failure)
!!FM       end if
!!FM       if (associated(glob_angles%unitary_t)) then
!!FM          call cp_assert(cp_validate(glob_angles%unitary_t,&
!!FM               error=error), cp_warning_level,&
!!FM               cp_assertion_failed, routineP,&
!!FM               "cached invalid unitary_t in "//&
!!FM               CPSourceFileRef,&
!!FM               error,failure)
!!FM          CPAssert(associated(glob_angles%ut_calculated),cp_warning_level,routineP,error,failure)
!!FM       end if
       if (associated(glob_angles%min_hamiltonian_m)) then
!!FM          call cp_assert(cp_validate(glob_angles%min_hamiltonian_m,&
!!FM               error=error), cp_warning_level,&
!!FM               cp_assertion_failed, routineP,&
!!FM               "cached invalid min_hamiltonian_m in "//&
!!FM               CPSourceFileRef,&
!!FM               error,failure)
       else
          CPAssert(.not.glob_angles%m_h_m_valid,cp_warning_level,routineP,error,failure)
       end if
       if (associated(glob_angles%min_overlap_m)) then
!!FM          call cp_assert(cp_validate(glob_angles%min_overlap_m,&
!!FM               error=error), cp_warning_level,&
!!FM               cp_assertion_failed, routineP,&
!!FM               "cached invalid min_overlap_m in "//&
!!FM               CPSourceFileRef,&
!!FM               error,failure)
       else
          CPAssert(.not.glob_angles%m_s_m_valid,cp_warning_level,routineP,error,failure)
       end if
    end if failureIf
    res=.not.failure
  end function pao_g_ang_validate
!***************************************************************************

!!****f* pao_glob_angles_methods/pao_g_ang_set [1.0] *
!!
!!   NAME
!!     pao_g_ang_set
!!
!!   SYNOPSIS
!!     Subroutine pao_g_ang_set(glob_angles, pao_env, angles, source,&
!!         ut_calculated, unitary_t, min_overlap_m, min_hamiltonian_m,&
!!         m_s_m_valid, m_h_m_valid, error)
!!       Type(pao_glob_angles_type), Intent (INOUT):: glob_angles
!!       Type(pao_env_type), Optional, Intent (IN), Target:: pao_env
!!       Real(Kind=wp), Optional, Dimension(:), Intent (IN):: angles
!!       Type(pao_glob_angles_type), Optional, Intent (IN):: source
!!       Type(pao_matrix_type), Optional, Intent (IN)::&
!!         min_hamiltonian_m, min_overlap_m, unitary_t
!!       Logical, Dimension(:), Optional, Intent (IN):: ut_calculated
!!       Logical, Intent (IN), Optional:: m_h_m_valid, m_s_m_valid
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_g_ang_set
!!
!!   FUNCTION
!!     can change the value of the various components of pao_glob_angles
!!
!!   NOTES
!!     you should almost never need this function
!!
!!   INPUTS
!!     glob_angles: the structure wose values are set
!!     pao_env: the pao env that this glob_angles refers to (SHARED)
!!     angles: the angles themselves as a long 1d array (copied)
!!     source: a glob angles that is copied into this (other arguments
!!             have the precedence for setting the values, the values
!!             of source are taken as default, if the value is not directly
!!             set)
!!     ut_calculated: flags to see which blocka have been calculated
!!     unitary_t: the cached unitary matrix. Be careful: only the 
!!             blocks that have been calculated are valid.
!!             Every processor calculates oly the blocks he needs.
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
  subroutine pao_g_ang_set(glob_angles, qs_env, global_env, angles,source,&
       ut_calculated, unitary_t,min_overlap_m, min_hamiltonian_m,&
       m_s_m_valid, m_h_m_valid, error)
    type(pao_glob_angles_type), intent(inout) :: glob_angles
    type(qs_environment_type), intent(in), target :: qs_env
    type(global_environment_type), intent(in), target :: global_env
    real(kind=wp), optional, dimension(:), intent(in) :: angles
    type(pao_glob_angles_type), optional, intent(in) :: source
    type(cp_sparse_matrix_type), optional, intent(in), target :: &
         min_overlap_m, unitary_t
    type(cp_sparse_matrix_p_type), dimension(:), intent(in),&
         optional :: min_hamiltonian_m
    logical, dimension(:), optional, intent(in) :: ut_calculated
    logical, intent(in), optional :: m_h_m_valid, m_s_m_valid
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_glob_angles_methods:pao_g_ang_set'
    type(pao_env_type), pointer :: pao_env
    integer :: stat
    type(cp_error_type) :: iError
    call cp_error_init(iError,template_error=error)
    failure=.false.

    call get_qs_env(qs_env, pao_env=pao_env)
    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(associated(pao_env),cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       if (present(source)) then
          CPPrecondition(source%initialized,cp_failure_level,routineP,error,failure)
          if (.not.failure) then
             call cp_set(iError, print_level=max(cp_failure_level,&
                  cp_error_get_print_level(iError)),&
                  stop_level=max(cp_failure_level,&
                  cp_error_get_stop_level(iError)))
             if (associated(source%angles)) then
                call cp_1d_r_guarantee_size(glob_angles%angles,size(source%angles),&
                     error=iError)
                call cp_assert(cp_error_get_level(iError)<cp_failure_level,&
                     cp_error_get_level(iError),cp_assertion_failed,routineP,&
                     "sub routine cp_guarantee_size error in"//&
                     CPSourceFileRef,&
                     error,failure)
                glob_angles%angles=source%angles
             else if (associated(glob_angles%angles)) then
                deallocate(glob_angles%angles, stat=stat)
                CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
             end if
             call pao_g_angles_set_unitary_t(glob_angles,&
                  ut_calculated=source%ut_calculated,&
                  new_value=source%unitary_t,&
                  error=error)
             call pao_g_angles_set_min_s_m(glob_angles,&
                  new_value=source%min_overlap_m,&
                  error=error)
             call pao_g_angles_set_min_h_m(glob_angles,&
                  new_value=source%min_hamiltonian_m,&
                  error=error)
             glob_angles%m_h_m_valid=source%m_h_m_valid
             glob_angles%m_s_m_valid=source%m_s_m_valid
             call cp_dealloc_ref(iError,error)
             call cp_init(iError, template_error=error)
          end if
       end if
       if (present(angles)) then
          call cp_guarantee_size(glob_angles%angles,size(angles),&
               error=iError)
          call cp_assert(cp_error_get_level(iError)<cp_failure_level,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "error in cp_guarantee_size called in "//&
               CPSourceFileRef,&
               error,failure)
          if (.not. failure) then
             glob_angles%angles = angles
          end if
       end if
       call pao_g_angles_set_unitary_t(glob_angles,&
            ut_calculated=ut_calculated,&
            new_value=unitary_t,&
            error=error)
       if (present(min_overlap_m)) then
          call pao_g_angles_set_min_s_m(glob_angles,&
               new_value=min_overlap_m,&
               error=error)
       end if
       if (present(min_hamiltonian_m)) then
          call pao_g_angles_set_min_h_m(glob_angles,&
               new_value=min_hamiltonian_m,&
               error=error)
       end if
       if (.not.failure) then
          call cp_guarantee_size(glob_angles%angles,&
               pao_env%angles_begin_at(size(pao_env%angles_begin_at))-1,&
               error=error)
       end if
    end if failureIf
    call cp_dealloc_ref(iError,error=error)
  end subroutine pao_g_ang_set
!***************************************************************************

!!****f* pao_glob_angles_methods/pao_g_ang_get [1.0] *
!!
!!   NAME
!!     pao_g_ang_get
!!
!!   SYNOPSIS
!!     Subroutine pao_g_ang_get(glob_angles, pao_env, angles,&
!!         m_s_m_valid, m_h_m_valid, ut_calculated, unitary_t, error)
!!       Type(pao_glob_angles_type), Intent (INOUT):: glob_angles
!!       Type(pao_env_type), Optional, Pointer:: pao_env
!!       Real(Kind=wp), Optional, Dimension(:), Pointer:: angles
!!       Logical, Intent (OUT), Optional:: m_h_m_valid, m_s_m_valid
!!       Logical, Dimension(:), Pointer, Optional:: ut_calculated
!!       Type(pao_matrix_type), Pointer, Optional:: unitary_t
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_g_ang_get
!!
!!   FUNCTION
!!     retrive the various attributes of the glob angles object
!!
!!   NOTES
!!     try not to use the angles attribute
!!
!!   INPUTS
!!     glob_angles: the object you get the info from
!!     pao_env: the pao environement this object lives in SHARED
!!     angles: a big 1D array with all the angles that define the unitary
!!            transformation (try not to use this)
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling
!!     ut_calculated: which block diagonal blocks of the unitary trasformation
!!            are valid (allocated if necessary)
!!     unitary_t: unitary trasformation. Be careful, only the blocks flagged
!!            in ut_calculated are valid (allocated if necessary)
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine pao_g_ang_get(glob_angles, qs_env, global_env, angles,&
       m_s_m_valid, m_h_m_valid, ut_calculated, unitary_t, &
       error)
    type(pao_glob_angles_type), intent(inout) :: glob_angles
    type(qs_environment_type), intent(in), target :: qs_env
    type(global_environment_type), intent(in), target :: global_env
    real(kind=wp), optional, dimension(:), pointer :: angles
    logical, intent(out), optional :: m_h_m_valid, m_s_m_valid
    logical, dimension(:), pointer, optional :: ut_calculated
    type(cp_sparse_matrix_type), pointer, optional :: unitary_t
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_glob_angles_methods:pao_g_ang_get'
    type(cp_error_type) :: iError
    type(pao_env_type), pointer :: pao_env
    type(cp_b_matrix_struct_type), pointer :: full_bas_m_struct
    integer :: stat, n_atoms
    call cp_init(iError,template_error=error)
    failure=.false.

    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    call get_qs_env(qs_env,natoms=n_atoms)
    failureIf: if (.not. failure) then
       if (present(angles)) angles => glob_angles%angles
       if (present(m_s_m_valid)) m_s_m_valid = glob_angles%m_s_m_valid
       if (present(m_h_m_valid)) m_h_m_valid=glob_angles%m_h_m_valid
       if (present(unitary_t)) then
          call get_qs_env(qs_env, pao_env=pao_env)
          call cp_get(pao_env, full_bas_m_struct=full_bas_m_struct)
          if (.not.associated(glob_angles%unitary_t)) then
             call cp_create(glob_angles%unitary_t,global_env=global_env,&
                  matrix_structure=full_bas_m_struct,&
                  error=error)
          end if
          if (.not.associated(glob_angles%ut_calculated)) then
             allocate(glob_angles%ut_calculated(n_atoms),&
                  stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             if (stat==0) glob_angles%ut_calculated=.false.
          end if
          unitary_t => glob_angles%unitary_t
       end if
       if (present(ut_calculated).and.&
            .not.associated(glob_angles%ut_calculated)) then
          allocate(glob_angles%ut_calculated(n_atoms),&
               stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          if (stat==0) glob_angles%ut_calculated=.false.
          ut_calculated => glob_angles%ut_calculated
       end if
    end if failureIf
  end subroutine pao_g_ang_get
!***************************************************************************

!!****f* pao_glob_angles_methods/pao_g_ang_valid [1.0] *
!!
!!   NAME
!!     pao_g_ang_valid
!!
!!   SYNOPSIS
!!     Function pao_g_ang_valid(glob_angles, error) Result(res)
!!       Logical:: res
!!       Type(pao_glob_angles_type), Intent (IN):: glob_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function pao_g_ang_valid
!!
!!   FUNCTION
!!     minimal validation on the glob angles object (no dangling pointers)
!!
!!   NOTES
!!     do not check pao_env?
!!
!!   INPUTS
!!     - glob_angles: the globa angles object to test
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function pao_g_ang_valid(glob_angles, qs_env, global_env, error) result(res)
    logical ::res
    type(pao_glob_angles_type), intent(in) :: glob_angles
    type(qs_environment_type), intent(in), target :: qs_env
    type(global_environment_type), intent(in), target :: global_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_glob_angles_methods:pao_g_ang_valid'
    failure=.false.

    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       if (.not.(associated(glob_angles%angles).and.&
            (.not.associated(glob_angles%unitary_t).or.&
            associated(glob_angles%ut_calculated)) .and.&
            (.not.glob_angles%m_s_m_valid .or.&
            associated(glob_angles%min_overlap_m)).and.&
            (.not.glob_angles%m_h_m_valid .or.&
            associated(glob_angles%min_hamiltonian_m))&
            )) then
          CPAssert(associated(glob_angles%angles),cp_warning_level,routineP,error,failure)
          if (associated(glob_angles%unitary_t)) then
             CPAssert(associated(glob_angles%ut_calculated),cp_warning_level,routineP,error,failure)
          end if
          if (glob_angles%m_s_m_valid) then
             CPAssert(associated(glob_angles%min_overlap_m),cp_warning_level,routineP,error,failure)
          end if
          if (glob_angles%m_h_m_valid) then
             CPAssert(associated(glob_angles%min_hamiltonian_m),cp_warning_level,routineP,error,failure)
          end if
       end if
    end if
    res=.not.failure
  end function pao_g_ang_valid
!***************************************************************************

! ===== special methods ====

!!****f* pao_glob_angles_methods/pao_g_angles_set_unitary_t [1.0] *
!!
!!   NAME
!!     pao_g_angles_set_unitary_t
!!
!!   FUNCTION
!!     sets the unitary transformation
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: the structure wose values are set
!!     new_value: the new cached unitary matrix.
!!     ut_calculated: which blocks are valid
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
  subroutine pao_g_angles_set_unitary_t(glob_angles,&
       new_value, ut_calculated,error)
    type(pao_glob_angles_type), intent(inout) :: glob_angles
    type(cp_sparse_matrix_type), target, optional, intent(in) :: new_value
    logical, dimension(:), intent(in), optional, target :: ut_calculated
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_glob_angles_methods:pao_g_angles_set_unitary_t'
    logical, dimension(:),pointer :: ut_c
    type(cp_sparse_matrix_type), pointer :: n_value
    type (cp_error_type) :: iError
    call cp_init(iError, template_error=error, &
         print_level=max(cp_failure_level, cp_error_get_print_level(error)))
    failure=.false.

    if (present(ut_calculated)) then
       ut_c => ut_calculated
       if (associated(ut_c)) then
          call cp_guarantee_size(glob_angles%ut_calculated,&
               size(ut_calculated),error=iError)
          call cp_error_propagate_error(iError,routineP,&
               CPSourceFileRef,&
               error,failure,cp_failure_level)
          glob_angles%ut_calculated=ut_calculated
       else if (associated(glob_angles%ut_calculated)) then
          deallocate(glob_angles%ut_calculated)
       end if
    end if
    if (present(new_value)) then
       n_value => new_value
       if (associated(n_value)) then
          call cp_retain(n_value,error=error)
          call cp_release(glob_angles%unitary_t, error=error)
          glob_angles%unitary_t => new_value
       else if (associated(glob_angles%unitary_t)) then
! do not dealloc ??
          call cp_dealloc(glob_angles%unitary_t,error=error)
       end if
    end if
  end subroutine pao_g_angles_set_unitary_t
!***************************************************************************

!!****f* pao_glob_angles_methods/pao_g_angles_set_min_s_m [1.0] *
!!
!!   NAME
!!     pao_g_angles_set_min_s_m
!!
!!   FUNCTION
!!     sets the unitary transformation
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: the structure wose values are set
!!     min_overlap_m: the overlap matrix in the minimal basis
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
  subroutine pao_g_angles_set_min_s_m(glob_angles, new_value,error)
    type(pao_glob_angles_type), intent(inout) :: glob_angles
    type(cp_sparse_matrix_type), target :: new_value
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='pao_glob_angles_methods:pao_g_angles_set_min_s_m'
    type(cp_sparse_matrix_type), pointer :: n_value
    failure=.false.
    n_value => new_value

    if (associated(n_value)) then
       glob_angles%m_s_m_valid=.true.
       call cp_retain(new_value, error=error)
       call cp_release(glob_angles%min_overlap_m,error=error)
       glob_angles%min_overlap_m => new_value
    else if (associated(glob_angles%min_overlap_m)) then
! do not dealloc ??
       call cp_dealloc(glob_angles%min_overlap_m,error=error)
    end if
  end subroutine pao_g_angles_set_min_s_m
!***************************************************************************

!!****f* pao_glob_angles_methods/pao_g_angles_set_min_h_m [1.0] *
!!
!!   NAME
!!     pao_g_angles_set_min_h_m
!!
!!   FUNCTION
!!     sets the unitary transformation
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: the structure wose values are set
!!     min_hamiltonian_m: the hamiltonian matrix in the minimal basis
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
  subroutine pao_g_angles_set_min_h_m(glob_angles, new_value,error)
    type(pao_glob_angles_type), intent(inout) :: glob_angles
    type(cp_sparse_matrix_p_type), dimension(:), target :: new_value
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='pao_g_angles_set_min_h_m',&
         routineP=moduleN//':'//routineN
    type(cp_sparse_matrix_p_type), dimension(:), pointer :: n_value
    integer :: i, stat
    type(cp_sparse_matrix_p_type), dimension(:), pointer :: tmp_array
    failure=.false.
    n_value => new_value

    if (associated(n_value)) then
       glob_angles%m_h_m_valid=.true.
       allocate(tmp_array(size(new_value)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          do i=1,size(new_value)
             call cp_retain(new_value(i)%matrix, error=error)
             tmp_array(i)%matrix => new_value(i)%matrix
          end do
          call cp_dealloc(glob_angles%min_hamiltonian_m, error=error)
          glob_angles%min_hamiltonian_m => tmp_array
       end if
    else if (associated(glob_angles%min_hamiltonian_m)) then
       ! do not dealloc ??
       call cp_dealloc(glob_angles%min_hamiltonian_m,error=error)
    end if
  end subroutine pao_g_angles_set_min_h_m
!***************************************************************************

!!****f* pao_glob_angles_methods/pao_g_ang_did_change [1.0] *
!!
!!   NAME
!!     pao_g_ang_did_change
!!
!!   SYNOPSIS
!!     Subroutine pao_g_ang_did_change(glob_angles, h_changed, s_changed,&
!!         angles_changed, full_reset, error)
!!       Type(pao_glob_angles_type), Intent (INOUT):: glob_angles
!!       Logical, Intent (IN), Optional:: h_changed, s_changed,&
!!         angles_changed, full_reset
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_g_ang_did_change
!!
!!   FUNCTION
!!     Invalidate all the cached values.
!!     Should be called when the value of the angles change.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: the globa angles object that should reset his cache
!!     full_reset: if true deallocates all the memory of the cache
!!            nedeed when the matrix structure changes (i.e. when the 
!!            atoms move). defaults to false
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
  subroutine pao_g_ang_did_change(glob_angles,h_changed, s_changed,&
       angles_changed, full_reset, error)
    type(pao_glob_angles_type), intent(inout) :: glob_angles
    logical, intent(in), optional :: h_changed, s_changed, &
         angles_changed, full_reset
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: handle, stat
    character(len=*), parameter :: routineN='pao_g_ang_did_change',&
         routineP=moduleN//':'//routineN
    logical :: h_c, s_c, a_c, full_r
    !call timeset(routineN,'I',moduleN,handle)
    failure=.false.
    h_c=.true. ; s_c=.true. ; a_c=.true. ; full_r=.false.

    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       if (present(h_changed)) h_c=h_changed
       if (present(s_changed)) s_c=s_changed
       if (present(angles_changed)) a_c=angles_changed
       if (present(full_reset)) full_r=full_reset
       glob_angles%m_h_m_valid=glob_angles%m_h_m_valid.and.&
            .not.(h_c.or.s_c.or.a_c.or.full_r)
       glob_angles%m_s_m_valid=glob_angles%m_s_m_valid.and.&
            .not.(s_c.or.a_c.or.full_r)
       if (associated(glob_angles%ut_calculated).and.&
            (a_c.or.full_r)) then
          glob_angles%ut_calculated=.false.
       end if
       if (full_r.or.s_c) then
          call cp_dealloc(glob_angles%unitary_t,error=error)
          call cp_dealloc(glob_angles%min_hamiltonian_m,error=error)
          call cp_dealloc(glob_angles%min_overlap_m)
          if (associated(glob_angles%ut_calculated)) then
             ! do not dealloc??
             deallocate(glob_angles%ut_calculated, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          end if
       end if
    end if
    !call timestop(0.0_wp,handle)
  end subroutine pao_g_ang_did_change
!***************************************************************************

end module pao_glob_angles_methods
