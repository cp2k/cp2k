!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_pao_utils [1.0] *
!!
!!   NAME
!!     cp_pao_utils
!!
!!   FUNCTION
!!     all kinds of utilities
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
module cp_pao_utils
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp => dp
  use global_types, only: global_environment_type
  use matrix_types
  use message_passing, only: mp_sum
  implicit none
  private
  logical, private, parameter :: debug_this_module=.true.
  character(len=*), parameter :: moduleN='cp_pao_utils'

! =========== public objects ==========
! types
  public :: cp_1d_i_array

! procedures
  public :: cp_guarantee_size, cp_get_random_integer


! ========== util types =========
  type cp_1d_i_array
     integer, dimension(:), pointer :: array
  end type cp_1d_i_array

! ============== interfaces ===========

!! if the pointer does not point to an array of the requested size 
!! allocates it anew
  interface cp_guarantee_size
     module procedure cp_guarantee_size_1d_int, cp_guarantee_size_1d_real, &
          cp_guarantee_size_1d_logical
  end interface


!!***
!****************************************************************************


contains

!!****f* cp_pao_utils/cp_get_random_integer [1.0] *
!!
!!   NAME
!!     cp_get_random_integer
!!
!!   SYNOPSIS
!!     Function cp_get_random_integer(from, to, error) Result(res)
!!       Integer:: res
!!       Integer, Optional, Intent (IN):: from, to
!!       Type(cp_error_type), Optional, Intent (IN):: error
!!     End Function cp_get_random_integer
!!
!!   FUNCTION
!!     This function returns a random integer in the requested range
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     from: minumum random number defaults to 1
!!     to: maximum random number defaults to 100
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling, but does NOT report
!!            the error in this variable (to make it pure)
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!****************************************************************************
  !pure
  function cp_get_random_integer(from,to,error) result(res)
    integer :: res
    integer, optional, intent(in) :: from, to
    type(cp_error_type), optional, intent(in) :: error

    logical :: failure
    character(len=*), parameter :: routineP='util:cp_get_random_integer'
    real :: rand
    integer :: fromi, toi
    type(cp_error_type) :: iError
    failure=.FALSE.
!call cp_init(iError,template_error=error)
    if (.not. present(from)) then
       fromi=1
    else
       fromi=from
    end if
    if (.not. present(to)) then
       toi=100
    else
       toi=to
    end if
!CPPrecondition(fromi<=toi,cp_failure_level,routineP,iError,failure)
    if (.not.fromi<=toi) failure=.true.
    if (.NOT. failure) then
       call random_number(rand)
       res=fromi+floor(rand*real(toi-fromi+1))
!CPPostcondition(res<=toi.and.res>=fromi,cp_warning_level,routineP,iError,failure)
       if (res>toi) res=toi
       if (res<fromi) res=fromi
    else
       res=-1
    end if
!call cp_dealloc_ref(iError,error=error)
  end function cp_get_random_integer
!!***
!****************************************************************************


! ========== memory utilities =========

!!****f* cp_pao_utils/cp_guarantee_size_1d_int [1.0] *
!!
!!   NAME
!!     cp_guarantee_size_1d_int
!!
!!   FUNCTION
!!     If the size of the array is changes reallocate it.
!!     Issues a warning when the size changes (but not on allocation
!!     and deallocation).
!!     
!!     The data is NOT preserved (if you want to preserve the data see
!!     the realloc in the module memory_utilities)
!!
!!   NOTES
!!     this is a different behaviour than the realloc in the module
!!     memory_utilities. It is quite low level
!!
!!   INPUTS
!!     array: the array to reallocate if necessary
!!     n: the wanted size
!!     init: the initial value of the elements when resized
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
!!   SOURCE
!***************************************************************************
  subroutine cp_guarantee_size_1d_int(array, n, init, error)
    integer, dimension(:), pointer :: array
    integer, intent(in) :: n
    integer, intent(in), optional :: init
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: stat
    character(len=*), parameter :: routineP='cp_pao_utils:cp_guarantee_size_1d_int'
    failure=.false.

    CPPrecondition(n>=0,cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       if (associated(array)) then
          if (size(array) /= n) then
             CPErrorMessage(cp_warning_level,routineP,'size has changed',error)
             deallocate(array, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          end if
       end if
       if (.not.associated(array)) then
          allocate(array(n), stat=stat)
          CPPostcondition(stat==0, cp_failure_level,routineP,error,failure)
          if (present(init).and..not.failure) array=init
       end if
    end if failureIf
  end subroutine cp_guarantee_size_1d_int
!!***
!***************************************************************************

!!****f* cp_pao_utils/cp_guarantee_size_1d_logical [1.0] *
!!
!!   NAME
!!     cp_guarantee_size_1d_logical
!!
!!   FUNCTION
!!     If the size of the array is changes reallocate it.
!!     Issues a warning when the size changes (but not on allocation
!!     and deallocation).
!!     
!!     The data is NOT preserved (if you want to preserve the data see
!!     the realloc in the module memory_utilities)
!!
!!   NOTES
!!     this is a different behaviour than the realloc in the module
!!     memory_utilities. It is quite low level
!!
!!   INPUTS
!!     array: the array to reallocate if necessary
!!     n: the wanted size
!!     init: the initial value of the elements when resized
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
!!   SOURCE
!***************************************************************************
  subroutine cp_guarantee_size_1d_logical(array, n, init, error)
    logical, dimension(:), pointer :: array
    integer, intent(in) :: n
    logical, intent(in), optional :: init
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: stat
    character(len=*), parameter :: routineP='cp_pao_utils:cp_guarantee_size_1d_logical'
    failure=.false.

    CPPrecondition(n>=0,cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       if (associated(array)) then
          if (size(array) /= n) then
             CPErrorMessage(cp_warning_level,routineP,'size has changed',error)
             deallocate(array, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          end if
       end if
       if (.not.associated(array)) then
          allocate(array(n), stat=stat)
          CPPostcondition(stat==0, cp_failure_level,routineP,error,failure)
          if (present(init).and..not.failure) array=init
       end if
    end if failureIf
  end subroutine cp_guarantee_size_1d_logical
!!***
!***************************************************************************

!!****f* cp_pao_utils/cp_guarantee_size_1d_real [1.0] *
!!
!!   NAME
!!     cp_guarantee_size_1d_real
!!
!!   FUNCTION
!!     If the size of the array is changes reallocate it.
!!     If warning_on_change is present and true issues a warning 
!!     whe the size changes (but not on allocation and deallocation).
!!     
!!     The data is NOT preserved (if you want to preserve the data see
!!     the realloc in the module memory_utilities)
!!
!!   NOTES
!!     this is a different behaviour than the realloc in the module
!!     memory_utilities. It is quite low level
!!
!!   INPUTS
!!     array: the array to reallocate if necessary
!!     n: the wanted size
!!     init: the initial value of the elements when resized
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
!!   SOURCE
!***************************************************************************
  subroutine cp_guarantee_size_1d_real(array, n, init, error)
    real(kind=wp), dimension(:), pointer :: array
    integer, intent(in) :: n
    real(kind=wp), optional, intent(in) :: init
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: stat
    character(len=*), parameter :: routineP='cp_pao_utils:cp_guarantee_size_1d_real'
    failure=.false.

    CPPrecondition(n>=0,cp_failure_level,routineP,error,failure)
    failureIf: if (.not. failure) then
       if (associated(array)) then
          if (size(array) /= n) then
             CPErrorMessage(cp_warning_level,routineP,'size has changed',error)
             deallocate(array, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          end if
       end if
       if (.not.associated(array)) then
          allocate(array(n), stat=stat)
          CPPostcondition(stat==0, cp_failure_level,routineP,error,failure)
          if (present(init).and..not.failure) array=init
       end if
    end if failureIf
  end subroutine cp_guarantee_size_1d_real
!!***
!***************************************************************************

end module cp_pao_utils

