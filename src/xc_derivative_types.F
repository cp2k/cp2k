!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****m* cp2k/xc_derivative_types *
!!
!!   NAME
!!     xc_derivative_types
!!
!!   FUNCTION
!!     Provides types for the management of the xc-functionals and
!!     their derivatives.
!!
!!   NOTES
!!     
!!     
!!****
module xc_derivative_types

  use kinds,              only: wp => dp
  use cp_error_handling,  only: cp_error_type, cp_a_l, cp_assert, cp_error_get_logger
  use cp_log_handling,    only: cp_logger_type, cp_failure_level, cp_warning_level, &
                                cp_to_string, cp_logger_get_default_unit_nr
  use xc_derivative_desc, only: MAX_DERIVATIVE_DESC_LENGTH, &
                                standardize_dervitative_desc

  implicit none

  private

  character(len=*), parameter, private :: moduleN='xc_derivative_types'

!!****t* cp2k/xc_derivative_type  *
!!
!!   NAME
!!     xc_derivative_type
!!
!!   FUNCTION
!!     A node of a linked list contains the actual derivative together
!!     with a string describing it.
!!
!!****
type xc_derivative_type
   integer                                   :: ref_count, id_nr
   character(len=MAX_DERIVATIVE_DESC_LENGTH) :: derivative_desc
   real(kind=wp), dimension(:), pointer      :: derivative
   type(xc_derivative_type), pointer         :: next
end type xc_derivative_type

!!****t* cp2k/xc_derivative_set_type  *
!!
!!   NAME
!!     xc_derivative_set_type
!!
!!   FUNCTION
!!     A derivative set contains the different derivatives of a xc-functional
!!     in form of a linked list
!!
!!****
type xc_derivative_set_type
   integer :: ref_count, id_nr, deriv_size
   type(xc_derivative_type), pointer :: first
end type xc_derivative_set_type


contains

! *****************************************************************************

!!****f* xc_derivative_types/xc_dset_get_derivative [1.0] *
!!
!!   NAME
!!     xc_dset_get_derivative
!!
!!   FUNCTION
!!     returns the requested xc_derivative
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - deriv_set: the set where to search for the derivative
!!     - deriv_desc: the description of the derivative you want to have
!!     - allocate_deriv: if the derivative should be allocated when not present
!!       Defaults to false.
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi & Thomas
!!
!!*** **********************************************************************
function xc_dset_get_derivative(deriv_set, deriv_desc, allocate_deriv, error)&
     result(res)
  type(xc_derivative_set_type), pointer :: deriv_set
  character(len=*), intent(in) :: deriv_desc
  logical, intent(in), optional :: allocate_deriv
  type(cp_error_type), optional, intent(inout) :: error
  type(xc_derivative_type), pointer :: res
  
  logical :: failure, my_allocate_deriv
  character(len=*), parameter :: routineN='xc_dset_get_derivative',&
        routineP=moduleN//':'//routineN
  character(len=MAX_DERIVATIVE_DESC_LENGTH) :: std_deriv_desc
  type(xc_derivative_type), pointer :: pos

  failure=.false.
  CPPrecondition(associated(deriv_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
  my_allocate_deriv=.false.
  if (present(allocate_deriv)) my_allocate_deriv=allocate_deriv
  nullify(res)
  std_deriv_desc=standardize_dervitative_desc(deriv_desc)
  if (.not. failure) then
     pos => deriv_set%first
     do while(associated(pos))
        if (pos%derivative_desc==std_deriv_desc) then
           res => pos
           exit
        end if
        pos => pos%next
     end do
     if (.not.associated(res).and.my_allocate_deriv) then
        call xc_derivative_create(res,std_deriv_desc,deriv_set%deriv_size, error=error)
        res%next => deriv_set%first
        deriv_set%first => res
     end if
  end if
end function xc_dset_get_derivative
!***************************************************************************

!!****f* xc_derivative_types/xc_derivative_create [1.0] *
!!
!!   NAME
!!     xc_derivative_create
!!
!!   FUNCTION
!!     allocates and initializes a derivative type
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - xc_deriv: the object to create
!!     - deriv_desc: the derivative description
!!     - deriv_size: the size of the (local) point in the derivative array
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2003 created [fawzi]
!!
!!*** **********************************************************************
subroutine xc_derivative_create(xc_deriv, deriv_desc, deriv_size, error)
  type(xc_derivative_type), pointer :: xc_deriv
  character(len=*), intent(in) :: deriv_desc
  integer, intent(in) :: deriv_size
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='f',&
        routineP=moduleN//':'//routineN

  failure=.false.
  
  if (.not. failure) then
     
  end if
end subroutine xc_derivative_create
!***************************************************************************


end module xc_derivative_types

