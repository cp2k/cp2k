!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/xc_derivative_desc [1.0] *
!!
!!   NAME
!!     xc_derivative_desc
!!
!!   FUNCTION
!!     Module with functions to handle derivative descriptors.
!!
!!   NOTES
!!     -
!!
!!****
module xc_derivative_desc

  use cp_log_handling,   only: cp_fatal_level, cp_failure_level, cp_warning_level, &
                               cp_note_level, cp_to_string, cp_logger_type
  use cp_error_handling, only: cp_debug, cp_error_type, cp_error_init, cp_error_dealloc_ref, &
                               cp_error_message, cp_assert, cp_assertion_failed, &
                               cp_internal_error, cp_a_l, cp_error_check
  use kinds,             only: wp=>dp
  use global_types,      only: global_environment_type
  use timings,           only: timeset, timestop

  implicit none

  private

  integer, parameter :: &
       MAX_LABEL_LENGTH           = 10, &
       MAX_DERIVATIVE             = 4,  &
       MAX_DERIVATIVE_DESC_LENGTH = (MAX_LABEL_LENGTH+2)*MAX_DERIVATIVE

  logical, parameter :: debug_this_module=.false.

  character(len=*), parameter :: moduleN='xc_derivative_desc'

  public :: MAX_DERIVATIVE_DESC_LENGTH, &
            standardize_dervitative_desc, &
            create_split_derivative_desc

!****************************************************************************

contains

!!****f* xc_derivative_types/standardize_dervitative_desc *
!!
!!   NAME
!!     standardize_dervitative_desc
!!
!!   FUNCTION
!!     returns a unique string for a given derivative string
!!
!!   ARGUMENTS
!!     - deriv_desc: input string which describes the derivative and
!!                   should be standardized
!!     - error     : variable to control error logging, stopping,... 
!!                   see module cp_error_handling 
!!
!!****
 function standardize_dervitative_desc(deriv_desc, error) result(res)
    character(len=*), intent(in) :: deriv_desc
    type(cp_error_type), intent(inout), optional :: error
    character(len=MAX_DERIVATIVE_DESC_LENGTH) :: res

    logical :: failure
    character(len=*), parameter :: routineN='standardize_dervitative_desc',&
         routineP=moduleN//':'//routineN
    character(len=MAX_LABEL_LENGTH), dimension(:), pointer :: deriv_array

    integer :: i, pos, l_label
    logical :: ordered
    character(len=MAX_LABEL_LENGTH) :: tmp

    deriv_array => create_split_derivative_desc(deriv_desc, error=error)

    ordered=.false.
    do while(.not.ordered)
       ordered=.true.
       do i=1,size(deriv_array)-1
          if (deriv_array(i)>deriv_array(i+1)) then
             tmp=deriv_array(i+1)
             deriv_array(i+1)=deriv_array(i)
             deriv_array(i)=tmp
             ordered=.false.
          end if
       end do
    end do
    
    res=""
    pos=1
    do i=1,size(deriv_array)
       l_label=len_trim(deriv_array(i))
       res(pos:pos+l_label+1)='('//deriv_array(i)(1:l_label)//')'
       pos=pos+l_label+2
    end do

  end function standardize_dervitative_desc

  !!****f* xc_derivative_types/create_split_derivative_desc [1.0] *
  !!
  !!   NAME
  !!     create_split_derivative_desc
  !!
  !!   FUNCTION
  !!     returns an array of 1st order derivative descriptions
  !!
  !!   ARGUMENTS
  !!     - deriv_desc: input string which describes the derivative
  !!     - error     : variable to control error logging, stopping,... 
  !!                   see module cp_error_handling 
  !!
  !!*** *********************************************************************
  function create_split_derivative_desc(deriv_desc, error) result(res)
    implicit none
    character(len=*), intent(in) :: deriv_desc
    type(cp_error_type), intent(inout), optional :: error
    character(len=MAX_LABEL_LENGTH), dimension(:), pointer :: res

    logical :: failure
    character(len=*), parameter :: routineN='create_split_derivative_desc',&
         routineP=moduleN//':'//routineN
    integer :: nderiv, i, j, stat

    nderiv=0
    do i=1,len(deriv_desc)
       if (deriv_desc(i:i)=='(') then
          nderiv=nderiv+1
       end if
    end do

    allocate(res(nderiv),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    nderiv=0
    do i=1,len(deriv_desc)
       if (deriv_desc(i:i)=='(') then
          nderiv=nderiv+1
          do j=i+1,len(deriv_desc)
             if (deriv_desc(j:j)==')') exit
          end do
!tc: should we do a check on the derivative descriptions?
          res(nderiv)=deriv_desc(i+1:j-1)
       end if
    end do
    
  end function create_split_derivative_desc

! *****************************************************************************

end module xc_derivative_desc
