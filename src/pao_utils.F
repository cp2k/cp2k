!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_utils [1.0] *
!!
!!   NAME
!!     pao_utils
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
module pao_utils
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp => dp
  use global_types, only: global_environment_type
  use sparse_matrix_types
  use message_passing, only: mp_sum
  implicit none
  private
  logical, private, parameter :: debug_this_module=.true.
  character(len=*), parameter :: moduleN='pao_utils'

! =========== public objects ==========
! types

! procedures
  public :: cp_get_random_integer


! ============== interfaces ===========

!!FM!! if the pointer does not point to an array of the requested size 
!!FM!! allocates it anew
!!FM  interface cp_guarantee_size
!!FM     module procedure cp_guarantee_size_1d_int, cp_guarantee_size_1d_real, &
!!FM          cp_guarantee_size_1d_logical
!!FM  end interface


!!***
!****************************************************************************


contains

!!****f* pao_utils/cp_get_random_integer [1.0] *
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

!!FM! ========== memory utilities =========
!!FM
!!FM!!****f* pao_utils/cp_guarantee_size_1d_int [1.0] *
!!FM!!
!!FM!!   NAME
!!FM!!     cp_guarantee_size_1d_int
!!FM!!
!!FM!!   FUNCTION
!!FM!!     If the size of the array is changes reallocate it.
!!FM!!     Issues a warning when the size changes (but not on allocation
!!FM!!     and deallocation).
!!FM!!     
!!FM!!     The data is NOT preserved (if you want to preserve the data see
!!FM!!     the realloc in the module memory_utilities)
!!FM!!
!!FM!!   NOTES
!!FM!!     this is a different behaviour than the realloc in the module
!!FM!!     memory_utilities. It is quite low level
!!FM!!
!!FM!!   INPUTS
!!FM!!     array: the array to reallocate if necessary
!!FM!!     n: the wanted size
!!FM!!     init: the initial value of the elements when resized
!!FM!!     error: variable to control error logging, stopping,... 
!!FM!!            see module cp_error_handling 
!!FM!!
!!FM!!   AUTHOR
!!FM!!     @author Fawzi Mohamed
!!FM!!     @version 1.2002
!!FM!!
!!FM!!   MODIFICATION HISTORY
!!FM!!     none
!!FM!!
!!FM!!   SOURCE
!!FM!***************************************************************************
!!FM  subroutine cp_guarantee_size_1d_int(array, n, init, error)
!!FM    integer, dimension(:), pointer :: array
!!FM    integer, intent(in) :: n
!!FM    integer, intent(in), optional :: init
!!FM    type(cp_error_type), optional, intent(inout) :: error
!!FM
!!FM    logical :: failure
!!FM    integer :: stat
!!FM    character(len=*), parameter :: routineP='pao_utils:cp_guarantee_size_1d_int'
!!FM    failure=.false.
!!FM
!!FM    CPPrecondition(n>=0,cp_failure_level,routineP,error,failure)
!!FM    failureIf: if (.not. failure) then
!!FM       if (associated(array)) then
!!FM          if (size(array) /= n) then
!!FM             CPErrorMessage(cp_warning_level,routineP,'size has changed',error)
!!FM             deallocate(array, stat=stat)
!!FM             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
!!FM          end if
!!FM       end if
!!FM       if (.not.associated(array)) then
!!FM          allocate(array(n), stat=stat)
!!FM          CPPostcondition(stat==0, cp_failure_level,routineP,error,failure)
!!FM          if (present(init).and..not.failure) array=init
!!FM       end if
!!FM    end if failureIf
!!FM  end subroutine cp_guarantee_size_1d_int
!!FM!!***
!!FM!***************************************************************************
!!FM
!!FM!!****f* pao_utils/cp_guarantee_size_1d_logical [1.0] *
!!FM!!
!!FM!!   NAME
!!FM!!     cp_guarantee_size_1d_logical
!!FM!!
!!FM!!   FUNCTION
!!FM!!     If the size of the array is changes reallocate it.
!!FM!!     Issues a warning when the size changes (but not on allocation
!!FM!!     and deallocation).
!!FM!!     
!!FM!!     The data is NOT preserved (if you want to preserve the data see
!!FM!!     the realloc in the module memory_utilities)
!!FM!!
!!FM!!   NOTES
!!FM!!     this is a different behaviour than the realloc in the module
!!FM!!     memory_utilities. It is quite low level
!!FM!!
!!FM!!   INPUTS
!!FM!!     array: the array to reallocate if necessary
!!FM!!     n: the wanted size
!!FM!!     init: the initial value of the elements when resized
!!FM!!     error: variable to control error logging, stopping,... 
!!FM!!            see module cp_error_handling 
!!FM!!
!!FM!!   AUTHOR
!!FM!!     @author Fawzi Mohamed
!!FM!!     @version 1.2002
!!FM!!
!!FM!!   MODIFICATION HISTORY
!!FM!!     none
!!FM!!
!!FM!!   SOURCE
!!FM!***************************************************************************
!!FM  subroutine cp_guarantee_size_1d_logical(array, n, init, error)
!!FM    logical, dimension(:), pointer :: array
!!FM    integer, intent(in) :: n
!!FM    logical, intent(in), optional :: init
!!FM    type(cp_error_type), optional, intent(inout) :: error
!!FM
!!FM    logical :: failure
!!FM    integer :: stat
!!FM    character(len=*), parameter :: routineP='pao_utils:cp_guarantee_size_1d_logical'
!!FM    failure=.false.
!!FM
!!FM    CPPrecondition(n>=0,cp_failure_level,routineP,error,failure)
!!FM    failureIf: if (.not. failure) then
!!FM       if (associated(array)) then
!!FM          if (size(array) /= n) then
!!FM             CPErrorMessage(cp_warning_level,routineP,'size has changed',error)
!!FM             deallocate(array, stat=stat)
!!FM             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
!!FM          end if
!!FM       end if
!!FM       if (.not.associated(array)) then
!!FM          allocate(array(n), stat=stat)
!!FM          CPPostcondition(stat==0, cp_failure_level,routineP,error,failure)
!!FM          if (present(init).and..not.failure) array=init
!!FM       end if
!!FM    end if failureIf
!!FM  end subroutine cp_guarantee_size_1d_logical
!!FM!!***
!!FM!***************************************************************************
!!FM
!!FM!!****f* pao_utils/cp_guarantee_size_1d_real [1.0] *
!!FM!!
!!FM!!   NAME
!!FM!!     cp_guarantee_size_1d_real
!!FM!!
!!FM!!   FUNCTION
!!FM!!     If the size of the array is changes reallocate it.
!!FM!!     If warning_on_change is present and true issues a warning 
!!FM!!     whe the size changes (but not on allocation and deallocation).
!!FM!!     
!!FM!!     The data is NOT preserved (if you want to preserve the data see
!!FM!!     the realloc in the module memory_utilities)
!!FM!!
!!FM!!   NOTES
!!FM!!     this is a different behaviour than the realloc in the module
!!FM!!     memory_utilities. It is quite low level
!!FM!!
!!FM!!   INPUTS
!!FM!!     array: the array to reallocate if necessary
!!FM!!     n: the wanted size
!!FM!!     init: the initial value of the elements when resized
!!FM!!     error: variable to control error logging, stopping,... 
!!FM!!            see module cp_error_handling 
!!FM!!
!!FM!!   AUTHOR
!!FM!!     @author Fawzi Mohamed
!!FM!!     @version 1.2002
!!FM!!
!!FM!!   MODIFICATION HISTORY
!!FM!!     none
!!FM!!
!!FM!!   SOURCE
!!FM!***************************************************************************
!!FM  subroutine cp_guarantee_size_1d_real(array, n, init, error)
!!FM    real(kind=wp), dimension(:), pointer :: array
!!FM    integer, intent(in) :: n
!!FM    real(kind=wp), optional, intent(in) :: init
!!FM    type(cp_error_type), optional, intent(inout) :: error
!!FM
!!FM    logical :: failure
!!FM    integer :: stat
!!FM    character(len=*), parameter :: routineP='pao_utils:cp_guarantee_size_1d_real'
!!FM    failure=.false.
!!FM
!!FM    CPPrecondition(n>=0,cp_failure_level,routineP,error,failure)
!!FM    failureIf: if (.not. failure) then
!!FM       if (associated(array)) then
!!FM          if (size(array) /= n) then
!!FM             CPErrorMessage(cp_warning_level,routineP,'size has changed',error)
!!FM             deallocate(array, stat=stat)
!!FM             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
!!FM          end if
!!FM       end if
!!FM       if (.not.associated(array)) then
!!FM          allocate(array(n), stat=stat)
!!FM          CPPostcondition(stat==0, cp_failure_level,routineP,error,failure)
!!FM          if (present(init).and..not.failure) array=init
!!FM       end if
!!FM    end if failureIf
!!FM  end subroutine cp_guarantee_size_1d_real
!!FM!!***
!!FM!***************************************************************************

end module pao_utils

