!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_pao_glob_angles_methods [1.0] *
!!
!!   NAME
!!     cp_pao_glob_angles_methods
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
module cp_pao_glob_angles_methods
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use cp_pao_types
  use cp_pao_matrix_methods
  use cp_pao_utils
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_pao_glob_angles_methods'

! core procedures
  public :: cp_init, cp_dealloc_ref, cp_set, cp_get, cp_create, cp_dealloc,&
       cp_did_change

! core functions
  public :: cp_validate, cp_valid

! special functions
  public :: cp_get_pao_env

! the underlying functions
  public :: cp_pao_g_ang_init, cp_pao_g_ang_dealloc_ref, &
       cp_pao_g_ang_create, cp_pao_g_ang_dealloc,&
       cp_pao_g_ang_validate, cp_pao_g_ang_valid,&
       cp_pao_g_ang_set, cp_pao_g_ang_get, cp_pao_g_ang_did_change,&
       cp_pao_g_ang_get_pao_env

! ============== interfaces ===========

! +++ begin of core methods +++

!! initializes the first, argument.
!! Support optional initial values that depend on the 
!! first argument 
  interface cp_init
     module procedure cp_pao_g_ang_init
  end interface

!! Deallocates the memory allocated by the first argument.
!! Supports error=error
  interface cp_dealloc_ref
     module procedure cp_pao_g_ang_dealloc_ref
  end interface

!! Allocates and initialize a new object returning a pointer to it in 
!! the first argument
  interface cp_create
     module procedure cp_pao_g_ang_create
  end interface

!! Deallocates the memory used by the first argument, and the object
!! itself
  interface cp_dealloc
     module procedure cp_pao_g_ang_dealloc
  end interface

!! returns true if the content of the first argumenst is valid.
!! writes the errors ar warnings, full validation (long)
  interface cp_validate
     module procedure cp_pao_g_ang_validate
  end interface

!! returns true if the first argument is valid. Only a fast minimal 
!! validation is performed (no dangling pointers)
  interface cp_valid
     module procedure cp_pao_g_ang_valid
  end interface

!! sets the value of various optonal attributes of the firs argument.
!! the optional attributes depend on the type of the first
  interface cp_set
     module procedure cp_pao_g_ang_set
  end interface

!! gets values from the first argument via optional keyword arguments
!! that depend on the type of the first argument
  interface cp_get
     module procedure cp_pao_g_ang_get
  end interface

!! tells the object that its internal values have changes and its cached
!! values may be invalid
  interface cp_did_change
     module procedure cp_pao_g_ang_did_change
  end interface

! +++ end of the core methods +++

!! returns the pao environement the first object lives in
  interface cp_get_pao_env
     module procedure cp_pao_g_ang_get_pao_env
  end interface

!!***
!****************************************************************************

contains

! ======= core methods ======

!!****f* cp_pao_glob_angles_methods/cp_pao_g_ang_init [1.0] *
!!
!!   NAME
!!     cp_pao_g_ang_init
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_g_ang_init(glob_angles, pao_env, angles, source,&
!!         ut_calculated, unitary_t, min_overlap_m, min_hamiltonian_m,&
!!         error)
!!       Type(cp_pao_glob_angles_type), Intent (OUT):: glob_angles
!!       Type(cp_pao_env_type), Intent (IN), Target:: pao_env
!!       Real(Kind=wp), Optional, Dimension(:), Intent (IN):: angles
!!       Type(cp_pao_glob_angles_type), Intent (IN), Optional:: source
!!       Logical, Dimension(:), Intent (IN), Optional:: ut_calculated
!!       Type(cp_pao_matrix_type), Intent (IN), Optional:: unitary_t,&
!!         min_hamiltonian_m, min_overlap_m
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_g_ang_init
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
  subroutine cp_pao_g_ang_init(glob_angles, pao_env, angles, source,&
       ut_calculated, unitary_t, min_overlap_m, min_hamiltonian_m, error)
    type(cp_pao_glob_angles_type), intent(out) :: glob_angles
    type(cp_pao_env_type), intent(in), target :: pao_env
    real(kind=wp), optional, dimension(:), intent(in) :: angles
    type(cp_pao_glob_angles_type), intent(in), optional :: source
    logical, dimension(:), intent(in), optional :: ut_calculated
    type(cp_pao_matrix_type), intent(in), optional :: unitary_t,&
         min_hamiltonian_m, min_overlap_m
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_glob_angles_methods:cp_pao_g_ang_init'
    failure=.false.

    nullify(glob_angles%pao_env, glob_angles%angles, glob_angles%unitary_t,&
         glob_angles%ut_calculated, glob_angles%min_overlap_m, &
         glob_angles%min_hamiltonian_m)
    glob_angles%m_s_m_valid=.false.
    glob_angles%m_h_m_valid=.false.
    glob_angles%pao_env => pao_env
    glob_angles%initialized = .true.
    CPPrecondition(associated(glob_angles%pao_env),cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       call cp_pao_g_ang_set(glob_angles, angles=angles,&
            source=source, ut_calculated=ut_calculated,&
            unitary_t=unitary_t,min_hamiltonian_m=min_hamiltonian_m,&
            min_overlap_m=min_overlap_m, error=error)
    end if
  end subroutine cp_pao_g_ang_init
!***************************************************************************

!!****f* cp_pao_glob_angles_methods/cp_pao_g_ang_dealloc_ref [1.0] *
!!
!!   NAME
!!     cp_pao_g_ang_dealloc_ref
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_g_ang_dealloc_ref(glob_angles, error)
!!       Type(cp_pao_glob_angles_type), Intent (INOUT):: glob_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_g_ang_dealloc_ref
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
  subroutine cp_pao_g_ang_dealloc_ref(glob_angles, error)
    type(cp_pao_glob_angles_type), intent(inout) :: glob_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_glob_angles_methods:cp_pao_g_ang_dealloc_ref'
    integer ::stat
    failure=.false.

    CPPrecondition(glob_angles%initialized, cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       nullify(glob_angles%pao_env)
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
  end subroutine cp_pao_g_ang_dealloc_ref
!***************************************************************************

!!****f* cp_pao_glob_angles_methods/cp_pao_g_ang_dealloc [1.0] *
!!
!!   NAME
!!     cp_pao_g_ang_dealloc
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_g_ang_dealloc(glob_angles, error)
!!       Type(cp_pao_glob_angles_type), Pointer:: glob_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_g_ang_dealloc
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
  subroutine cp_pao_g_ang_dealloc(glob_angles,error)
    type(cp_pao_glob_angles_type), pointer :: glob_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_glob_angles_methods:cp_pao_g_ang_dealloc'
    integer :: stat
    failure=.false.

    if (associated(glob_angles)) then
       CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
       if (.not. failure) then
          call cp_pao_g_ang_dealloc_ref(glob_angles,error=error)
       end if
       deallocate(glob_angles,stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    end if
  end subroutine cp_pao_g_ang_dealloc
!***************************************************************************

!!****f* cp_pao_glob_angles_methods/cp_pao_g_ang_create [1.0] *
!!
!!   NAME
!!     cp_pao_g_ang_create
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_g_ang_create(glob_angles, pao_env, angles, source,&
!!         ut_calculated, unitary_t, min_overlap_m, min_hamiltonian_m,&
!!         error)
!!       Type(cp_pao_glob_angles_type), Pointer:: glob_angles
!!       Type(cp_pao_env_type), Intent (IN), Target:: pao_env
!!       Real(Kind=wp), Optional, Dimension(:), Intent (IN):: angles
!!       Type(cp_pao_glob_angles_type), Intent (IN), Optional:: source
!!       Logical, Dimension(:), Intent (IN):: ut_calculated
!!       Type(cp_pao_matrix_type), Intent (IN):: unitary_t,&
!!         min_hamiltonian_m, min_overlap_m
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_g_ang_create
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
!!     for the others inputs see cp_pao_g_ang_init
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_pao_g_ang_create(glob_angles, pao_env, angles, source,&
       ut_calculated, unitary_t, min_overlap_m, min_hamiltonian_m, error)
    type(cp_pao_glob_angles_type), pointer :: glob_angles
    type(cp_pao_env_type), intent(in), target :: pao_env
    real(kind=wp), optional, dimension(:), intent(in) :: angles
    type(cp_pao_glob_angles_type), intent(in), optional :: source
    logical, dimension(:), intent(in) :: ut_calculated
    type(cp_pao_matrix_type), intent(in) :: unitary_t, min_hamiltonian_m,&
         min_overlap_m
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_glob_angles_methods:cp_pao_g_ang_create'
    integer :: stat
    failure=.false.

    allocate(glob_angles, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       call cp_pao_g_ang_init(glob_angles, pao_env, angles, source,&
            ut_calculated, unitary_t, min_overlap_m, min_hamiltonian_m, error)
    end if
  end subroutine cp_pao_g_ang_create
!***************************************************************************

!!****f* cp_pao_glob_angles_methods/cp_pao_g_ang_validate [1.0] *
!!
!!   NAME
!!     cp_pao_g_ang_validate
!!
!!   SYNOPSIS
!!     Function cp_pao_g_ang_validate(glob_angles, error) Result(res)
!!       Logical:: res
!!       Type(cp_pao_glob_angles_type), Intent (IN):: glob_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_g_ang_validate
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
  function cp_pao_g_ang_validate(glob_angles,error) result(res)
    logical ::res
    type(cp_pao_glob_angles_type), intent(in) :: glob_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_glob_angles_methods:cp_pao_g_ang_validate'
    failure=.false.
    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       CPAssert(associated(glob_angles%angles),cp_warning_level,routineP,error,failure)
       CPAssert(associated(glob_angles%pao_env),cp_warning_level,routineP,error,failure)
       if (.not.failure) then
          call cp_assert(cp_valid(glob_angles%pao_env,&
               error=error),cp_warning_level,&
               cp_assertion_failed,routineP,&
               "validation of referred globAngles failed in "//&
               CPSourceFileRef,&
               error,failure)
       end if
       if (.not.failure) then
          call cp_assert(glob_angles%pao_env%angles_begin_at &
               (glob_angles%pao_env%n_atoms+1)== &
               size(glob_angles%angles),cp_warning_level,&
               cp_assertion_failed, routineP,&
               "number of angles different from last begin_at in "//&
               CPSourceFileRef,&
               error,failure)
          if (dot_product(glob_angles%pao_env%nr_min_bas,&
               glob_angles%pao_env%nr_full_bas-glob_angles%pao_env%nr_min_bas)&
               /= size(glob_angles%angles)) then
             call cp_assert(.false., cp_warning_level, cp_assertion_failed,&
                  routineP,'wrong number of angles in cp_pao_glob_angles_type',&
                  error,failure)
          end if
       end if
       if (associated(glob_angles%ut_calculated)) then
          call cp_assert(size(glob_angles%ut_calculated)==&
               glob_angles%pao_env%n_atoms, cp_warning_level,&
               cp_assertion_failed, routineP,&
               "cached invalid ut_calculated in "//&
               CPSourceFileRef,&
               error,failure)
       end if
       if (associated(glob_angles%unitary_t)) then
          call cp_assert(cp_validate(glob_angles%unitary_t,&
               error=error), cp_warning_level,&
               cp_assertion_failed, routineP,&
               "cached invalid unitary_t in "//&
               CPSourceFileRef,&
               error,failure)
          CPAssert(associated(glob_angles%ut_calculated),cp_warning_level,routineP,error,failure)
       end if
       if (associated(glob_angles%min_hamiltonian_m)) then
          call cp_assert(cp_validate(glob_angles%min_hamiltonian_m,&
               error=error), cp_warning_level,&
               cp_assertion_failed, routineP,&
               "cached invalid min_hamiltonian_m in "//&
               CPSourceFileRef,&
               error,failure)
       else
          CPAssert(.not.glob_angles%m_h_m_valid,cp_warning_level,routineP,error,failure)
       end if
       if (associated(glob_angles%min_overlap_m)) then
          call cp_assert(cp_validate(glob_angles%min_overlap_m,&
               error=error), cp_warning_level,&
               cp_assertion_failed, routineP,&
               "cached invalid min_overlap_m in "//&
               CPSourceFileRef,&
               error,failure)
       else
          CPAssert(.not.glob_angles%m_s_m_valid,cp_warning_level,routineP,error,failure)
       end if
    end if failureIf
    res=.not.failure
  end function cp_pao_g_ang_validate
!***************************************************************************

!!****f* cp_pao_glob_angles_methods/cp_pao_g_ang_set [1.0] *
!!
!!   NAME
!!     cp_pao_g_ang_set
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_g_ang_set(glob_angles, pao_env, angles, source,&
!!         ut_calculated, unitary_t, min_overlap_m, min_hamiltonian_m,&
!!         m_s_m_valid, m_h_m_valid, error)
!!       Type(cp_pao_glob_angles_type), Intent (INOUT):: glob_angles
!!       Type(cp_pao_env_type), Optional, Intent (IN), Target:: pao_env
!!       Real(Kind=wp), Optional, Dimension(:), Intent (IN):: angles
!!       Type(cp_pao_glob_angles_type), Optional, Intent (IN):: source
!!       Type(cp_pao_matrix_type), Optional, Intent (IN)::&
!!         min_hamiltonian_m, min_overlap_m, unitary_t
!!       Logical, Dimension(:), Optional, Intent (IN):: ut_calculated
!!       Logical, Intent (IN), Optional:: m_h_m_valid, m_s_m_valid
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_g_ang_set
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
  subroutine cp_pao_g_ang_set(glob_angles, pao_env, angles,source,&
       ut_calculated, unitary_t,min_overlap_m, min_hamiltonian_m,&
       m_s_m_valid, m_h_m_valid, error)
    type(cp_pao_glob_angles_type), intent(inout) :: glob_angles
    type(cp_pao_env_type), optional, intent(in), target :: pao_env
    real(kind=wp), optional, dimension(:), intent(in) :: angles
    type(cp_pao_glob_angles_type), optional, intent(in) :: source
    type(cp_pao_matrix_type), optional, intent(in) :: min_hamiltonian_m,&
         min_overlap_m, unitary_t
    logical, dimension(:), optional, intent(in) :: ut_calculated
    logical, intent(in), optional :: m_h_m_valid, m_s_m_valid
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_glob_angles_methods:cp_pao_g_ang_set'
    integer :: stat
    type(cp_error_type) :: iError
    call cp_init_error(iError,template_error=error)
    failure=.false.

    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       if (present(source)) then
          CPPrecondition(source%initialized,cp_failure_level,routineP,error,failure)
          if (.not.failure) then
             glob_angles%pao_env => source%pao_env
             call cp_set(iError, print_level=max(cp_failure_level,&
                  cp_get_error_print_level(iError)),&
                  stop_level=max(cp_failure_level,&
                  cp_get_error_stop_level(iError)))
             if (associated(source%angles)) then
                call cp_guarantee_size(glob_angles%angles,size(source%angles),&
                     error=iError)
                call cp_assert(cp_get_error_level(iError)<cp_failure_level,&
                     cp_get_error_level(iError),cp_assertion_failed,routineP,&
                     "sub routine cp_guarantee_size error in"//&
                     CPSourceFileRef,&
                     error,failure)
                glob_angles%angles=source%angles
             else if (associated(glob_angles%angles)) then
                deallocate(glob_angles%angles, stat=stat)
                CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
             end if
             call cp_pao_g_angles_set_unitary_t(glob_angles,&
                  ut_calculated=source%ut_calculated,&
                  new_value=source%unitary_t,&
                  error=error)
             call cp_pao_g_angles_set_min_s_m(glob_angles,&
                  new_value=source%min_overlap_m,&
                  error=error)
             call cp_pao_g_angles_set_min_h_m(glob_angles,&
                  new_value=source%min_hamiltonian_m,&
                  error=error)
             glob_angles%m_h_m_valid=source%m_h_m_valid
             glob_angles%m_s_m_valid=source%m_s_m_valid
             call cp_dealloc_ref(iError,error)
             call cp_init(iError, template_error=error)
          end if
       end if
       if (present(pao_env)) &
            glob_angles%pao_env => pao_env
       if (present(angles)) then
          call cp_guarantee_size(glob_angles%angles,size(angles),&
               error=iError)
          call cp_assert(cp_get_error_level(iError)<cp_failure_level,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "error in cp_guarantee_size called in "//&
               CPSourceFileRef,&
               error,failure)
          if (.not. failure) then
             glob_angles%angles = angles
          end if
       end if
       call cp_pao_g_angles_set_unitary_t(glob_angles,&
            ut_calculated=ut_calculated,&
            new_value=unitary_t,&
            error=error)
       if (present(min_overlap_m)) then
          call cp_pao_g_angles_set_min_s_m(glob_angles,&
               new_value=min_overlap_m,&
               error=error)
       end if
       if (present(min_hamiltonian_m)) then
          call cp_pao_g_angles_set_min_h_m(glob_angles,&
               new_value=min_hamiltonian_m,&
               error=error)
       end if
       if (cp_debug.and.debug_this_module) then
          if (associated(glob_angles%pao_env)) then
             call cp_assert(cp_valid(glob_angles%pao_env, &
                  error=error), cp_failure_level, &
                  cp_precondition_failed, routineP, &
                  "if set the pao_env must be valid in"//&
                  CPSourceFileRef,&
                  error,failure)
          end if
       end if
       if (.not.failure) then
          call cp_guarantee_size(glob_angles%angles,&
               dot_product(glob_angles%pao_env%nr_min_bas, &
               glob_angles%pao_env%nr_full_bas-&
               glob_angles%pao_env%nr_min_bas), &
               error=error)
       end if
    end if failureIf
    call cp_dealloc_ref(iError,error=error)
  end subroutine cp_pao_g_ang_set
!***************************************************************************

!!****f* cp_pao_glob_angles_methods/cp_pao_g_ang_get [1.0] *
!!
!!   NAME
!!     cp_pao_g_ang_get
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_g_ang_get(glob_angles, pao_env, angles,&
!!         m_s_m_valid, m_h_m_valid, ut_calculated, unitary_t, error)
!!       Type(cp_pao_glob_angles_type), Intent (INOUT):: glob_angles
!!       Type(cp_pao_env_type), Optional, Pointer:: pao_env
!!       Real(Kind=wp), Optional, Dimension(:), Pointer:: angles
!!       Logical, Intent (OUT), Optional:: m_h_m_valid, m_s_m_valid
!!       Logical, Dimension(:), Pointer, Optional:: ut_calculated
!!       Type(cp_pao_matrix_type), Pointer, Optional:: unitary_t
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_g_ang_get
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
  subroutine cp_pao_g_ang_get(glob_angles, pao_env, angles,&
       m_s_m_valid, m_h_m_valid, ut_calculated, unitary_t, &
       error)
    type(cp_pao_glob_angles_type), intent(inout) :: glob_angles
    type(cp_pao_env_type), optional, pointer :: pao_env
    real(kind=wp), optional, dimension(:), pointer :: angles
    logical, intent(out), optional :: m_h_m_valid, m_s_m_valid
    logical, dimension(:), pointer, optional :: ut_calculated
    type(cp_pao_matrix_type), pointer, optional :: unitary_t
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_glob_angles_methods:cp_pao_g_ang_get'
    type(cp_error_type) :: iError
    integer :: stat
    call cp_init(iError,template_error=error)
    failure=.false.

    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       if (present(pao_env)) pao_env => glob_angles%pao_env
       if (present(angles)) angles => glob_angles%angles
       if (present(m_s_m_valid)) m_s_m_valid = glob_angles%m_s_m_valid
       if (present(m_h_m_valid)) m_h_m_valid=glob_angles%m_h_m_valid
       if (present(unitary_t)) then
          if (.not.associated(glob_angles%unitary_t)) then
             call cp_create(glob_angles%unitary_t, min_bas_rows=.false.,&
                  min_bas_cols=.true., pao_env=glob_angles%pao_env,&
                  error=error)
          end if
          if (.not.associated(glob_angles%ut_calculated)) then
             allocate(glob_angles%ut_calculated(glob_angles%pao_env%n_atoms),&
                  stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             if (stat==0) glob_angles%ut_calculated=.false.
          end if
       end if
       if (present(ut_calculated).and.&
            .not.associated(glob_angles%ut_calculated)) then
          allocate(glob_angles%ut_calculated(glob_angles%pao_env%n_atoms),&
               stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          if (stat==0) glob_angles%ut_calculated=.false.
       end if
    end if failureIf
  end subroutine cp_pao_g_ang_get
!***************************************************************************

!!****f* cp_pao_glob_angles_methods/cp_pao_g_ang_valid [1.0] *
!!
!!   NAME
!!     cp_pao_g_ang_valid
!!
!!   SYNOPSIS
!!     Function cp_pao_g_ang_valid(glob_angles, error) Result(res)
!!       Logical:: res
!!       Type(cp_pao_glob_angles_type), Intent (IN):: glob_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_g_ang_valid
!!
!!   FUNCTION
!!     minimal validation on the glob angles object (no dangling pointers)
!!
!!   NOTES
!!     do not check pao_env?
!!
!!   INPUTS
!!     glob_angles: the globa angles object to test
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
  function cp_pao_g_ang_valid(glob_angles, error) result(res)
    logical ::res
    type(cp_pao_glob_angles_type), intent(in) :: glob_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_glob_angles_methods:cp_pao_g_ang_valid'
    failure=.false.

    CPPrecondition(glob_angles%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       if (.not.(associated(glob_angles%pao_env).and.&
            associated(glob_angles%angles).and.&
            (.not.associated(glob_angles%unitary_t).or.&
            associated(glob_angles%ut_calculated)) .and.&
            (.not.glob_angles%m_s_m_valid .or.&
            associated(glob_angles%min_overlap_m)).and.&
            (.not.glob_angles%m_h_m_valid .or.&
            associated(glob_angles%min_hamiltonian_m))&
            )) then
          CPAssert(associated(glob_angles%pao_env),cp_warning_level,routineP,error,failure)
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
       else if(cp_debug.and.debug_this_module) then
          if (.not.(associated(glob_angles%pao_env%nr_full_bas).and.&
               associated(glob_angles%pao_env%nr_min_bas))) then
             failure=cp_valid(glob_angles%pao_env, error=error)
          end if
       end if
    end if
    res=.not.failure
  end function cp_pao_g_ang_valid
!***************************************************************************

! ===== special methods ====

!!****f* cp_pao_glob_angles_methods/cp_pao_g_angles_set_unitary_t [1.0] *
!!
!!   NAME
!!     cp_pao_g_angles_set_unitary_t
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
  subroutine cp_pao_g_angles_set_unitary_t(glob_angles, new_value, &
       ut_calculated,error)
    type(cp_pao_glob_angles_type), intent(inout) :: glob_angles
    type(cp_pao_matrix_type), target, optional, intent(in) :: new_value
    logical, dimension(:), intent(in), optional, target :: ut_calculated
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_glob_angles_methods:cp_pao_g_angles_set_unitary_t'
    logical, dimension(:),pointer :: ut_c
    type(cp_pao_matrix_type), pointer :: n_value
    type (cp_error_type) :: iError
    call cp_init(iError, template_error=error, &
         print_level=max(cp_failure_level, cp_get_error_print_level(error)))
    failure=.false.

    if (present(ut_calculated)) then
       ut_c => ut_calculated
       if (associated(ut_c)) then
          call cp_guarantee_size(glob_angles%ut_calculated,&
               size(ut_calculated),error=iError)
          call cp_propagate_error(iError,routineP,&
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
          if (.not.associated(glob_angles%unitary_t)) then
             call cp_create(glob_angles%unitary_t, &
                  pao_env=glob_angles%pao_env,&
                  source=n_value,&
                  error=error)
          else
             call cp_set(glob_angles%unitary_t,&
                  source=n_value,&
                  error=error)
          end if
       else if (associated(glob_angles%unitary_t)) then
! do not dealloc ??
          call cp_dealloc(glob_angles%unitary_t,error=error)
       end if
    end if
  end subroutine cp_pao_g_angles_set_unitary_t
!***************************************************************************

!!****f* cp_pao_glob_angles_methods/cp_pao_g_angles_set_min_s_m [1.0] *
!!
!!   NAME
!!     cp_pao_g_angles_set_min_s_m
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
  subroutine cp_pao_g_angles_set_min_s_m(glob_angles, new_value,error)
    type(cp_pao_glob_angles_type), intent(inout) :: glob_angles
    type(cp_pao_matrix_type), target :: new_value
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_glob_angles_methods:cp_pao_g_angles_set_min_s_m'
    type(cp_pao_matrix_type), pointer :: n_value
    failure=.false.
    n_value => new_value

    if (associated(n_value)) then
       glob_angles%m_s_m_valid=.true.
       if (.not.associated(glob_angles%min_overlap_m)) then
          call cp_create(glob_angles%min_overlap_m, &
               pao_env=glob_angles%pao_env,&
               source=n_value,&
               error=error)
       else
          call cp_set(glob_angles%min_overlap_m,&
               source=n_value,&
               error=error)
       end if
    else if (associated(glob_angles%min_overlap_m)) then
! do not dealloc ??
       call cp_dealloc(glob_angles%min_overlap_m,error=error)
    end if
  end subroutine cp_pao_g_angles_set_min_s_m
!***************************************************************************

!!****f* cp_pao_glob_angles_methods/cp_pao_g_angles_set_min_h_m [1.0] *
!!
!!   NAME
!!     cp_pao_g_angles_set_min_h_m
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
  subroutine cp_pao_g_angles_set_min_h_m(glob_angles, new_value,error)
    type(cp_pao_glob_angles_type), intent(inout) :: glob_angles
    type(cp_pao_matrix_type), target :: new_value
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_glob_angles_methods:cp_pao_g_angles_set_min_h_m'
    type(cp_pao_matrix_type), pointer :: n_value
    failure=.false.
    n_value => new_value

    if (associated(n_value)) then
       glob_angles%m_h_m_valid=.true.
       if (.not.associated(glob_angles%min_hamiltonian_m)) then
          call cp_create(glob_angles%min_hamiltonian_m, &
               pao_env=glob_angles%pao_env,&
               source=n_value,&
               error=error)
       else
          call cp_set(glob_angles%min_hamiltonian_m,&
               source=n_value,&
               error=error)
       end if
    else if (associated(glob_angles%min_hamiltonian_m)) then
! do not dealloc ??
       call cp_dealloc(glob_angles%min_hamiltonian_m,error=error)
    end if
  end subroutine cp_pao_g_angles_set_min_h_m
!***************************************************************************

!!****f* cp_pao_glob_angles_methods/cp_pao_g_ang_get_pao_env [1.0] *
!!
!!   NAME
!!     cp_pao_g_ang_get_pao_env
!!
!!   SYNOPSIS
!!     Function cp_pao_g_ang_get_pao_env(glob_angles, error) Result(res)
!!       Type(cp_pao_env_type), Pointer:: res
!!       Type(cp_pao_glob_angles_type), Intent (IN):: glob_angles
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_pao_g_ang_get_pao_env
!!
!!   FUNCTION
!!     returns the global environement the global angles object lives in
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_angles: the global angles object you are interested in
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
  function cp_pao_g_ang_get_pao_env(glob_angles, error) result(res)
    type(cp_pao_env_type), pointer ::res
    type(cp_pao_glob_angles_type), intent(in) :: glob_angles
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_pao_glob_angles_methods:cp_pao_g_ang_get_pao_env'
    failure=.false.
    CPPrecondition(glob_angles%initialized, cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       res => glob_angles%pao_env
    else failureIf
       nullify(res)
    end if failureIf
  end function cp_pao_g_ang_get_pao_env
!***************************************************************************

!!****f* cp_pao_glob_angles_methods/cp_pao_g_ang_did_change [1.0] *
!!
!!   NAME
!!     cp_pao_g_ang_did_change
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_g_ang_did_change(glob_angles, h_changed, s_changed,&
!!         angles_changed, full_reset, error)
!!       Type(cp_pao_glob_angles_type), Intent (INOUT):: glob_angles
!!       Logical, Intent (IN), Optional:: h_changed, s_changed,&
!!         angles_changed, full_reset
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_g_ang_did_change
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
  subroutine cp_pao_g_ang_did_change(glob_angles,h_changed, s_changed,&
       angles_changed, full_reset, error)
    type(cp_pao_glob_angles_type), intent(inout) :: glob_angles
    logical, intent(in), optional :: h_changed, s_changed, &
         angles_changed, full_reset
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: handle, stat
    character(len=*), parameter :: routineN='cp_pao_g_ang_did_change',&
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
  end subroutine cp_pao_g_ang_did_change
!***************************************************************************

end module cp_pao_glob_angles_methods
