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
  use cp_log_handling, only: cp_to_string, cp_failure_level, cp_warning_level
  use cp_error_handling, only: cp_error_type, cp_assert, cp_error_message,&
       cp_assertion_failed
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



end module pao_utils

