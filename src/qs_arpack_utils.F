!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****m* cp2k/qs_arpack_utils *
!!
!!   NAME
!!     qs_arpack_utils
!!
!!   FUNCTION
!!     Utility functions for the arpack modules
!!
!!   SYNOPSIS
!!     -
!!
!!   MODIFICATION HISTORY
!!     10-2002, TCH, created
!!
!!****
module qs_arpack_utils

  use global_types, only : global_environment_type
  use termination,  only : stop_program

  implicit none

  private
  public :: arpack_confirm

contains

! *****************************************************************************

!!****f* cp2k/m/s *
!!
!!   NAME
!!     arpack_confirm
!!
!!   FUNCTION
!!     Checks if arpack is present (defined). If optional
!!     is_present is given, stores the information there. Otherwise
!!     stops the program.
!!
!!   ARGUMENTS
!!     - is_present : logical, optional. stores the answer
!!
!!   NOTES
!!     -
!!
!!   MODIFICATION HISTORY
!!     10-2002, TCH, extracted from qs_arpack_diag_sy
!!
!!
!!****
  subroutine arpack_confirm(globenv,is_present)
    
    implicit none
    
    ! arguments
    type(global_environment_type), intent(IN) :: globenv
    logical, optional, intent(INOUT)          :: is_present

#if defined(__ARPACK)
    if (present(is_present)) then
       is_present=.true.
    endif
#else
    if (present(is_present)) then
       is_present=.false.
    else
       call stop_program("arpack_confirm",& 
            "(p)arpack library not present but needed",globenv)
    endif
#endif
    
  end subroutine arpack_confirm
          
! *****************************************************************************

end module qs_arpack_utils
