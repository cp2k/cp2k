!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_sll_int_tests [1.0] *
!!
!!   NAME
!!     cp_sll_int_tests
!!
!!   FUNCTION
!!     do some tests on a linked list of a scalar numeric type (int)
!!
!!   NOTES
!!     this is a template to test both the pointer and non pointer version
!!     of the list. DO NOT modify the .F files, modify the .template file
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
module cp_sll_int_tests
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use cp_linked_list_int_ptr
  use cp_array_utils, only: cp_output
  use cp_linked_list_int
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character, private, parameter :: moduleN='cp_sll_int_tests'

  public :: cp_test_int_sll

!!***
!****************************************************************************
contains

!!****f* cp_sll_int_tests/cp_test_int_sll [1.0] *
!!
!!   NAME
!!     cp_test_int_sll
!!
!!   SYNOPSIS
!!     Function cp_test_int_sll(glob_env, error) Result(res)
!!       Logical:: res
!!       Type(global_environment_type), Target:: glob_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_test_int_sll
!!
!!   FUNCTION
!!     tests the integer single linked list
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     glob_env: global environement
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
function cp_test_int_sll(glob_env, error) result(res)
  logical :: res
  type(global_environment_type), target :: glob_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_test_int_sll',&
       routineP=moduleN//routineN
  integer, parameter :: n=6
  integer, dimension(n) :: array
  integer :: i
  type(cp_sll_int_type), pointer :: reverseList, iterator
  integer :: el
  failure=.false.
  nullify(reverseList,iterator)
  
  forall (i=1:size(array)) array(i)=i*i
  call cp_create(reverseList,error=error)
  do i=1,size(array)
     call cp_insert(reverseList,array(i),error=error)
  end do
  CPAssert(cp_get_length(reverseList)==n,cp_warning_level,routineP,error,failure)
  if (.not. failure) then
     iterator => reverseList
     do i=1,size(array)
        call cp_assert(cp_get_element_at(reverseList,i)==array(n-i+1),&
             cp_warning_level,cp_assertion_failed,routineP,&
             "reverse list ("//cp_to_string(i)//")="//cp_to_string(&
             cp_get_element_at(reverseList,i))//" /= "//cp_to_string(array(n-i+1))&
             //"array("//cp_to_string(n-i+1)//") in "//&
             CPSourceFileRef,&
             error,failure)
        if (.not.cp_next(iterator,el_att=el,error=error)) then
           CPPrecondition(.false.,cp_warning_level,routineP,error,failure)
        end if
        if (.not.failure) then
           call cp_assert(el==array(n-i+1),&
                cp_warning_level,cp_assertion_failed,routineP,&
                "iterator /= array(n-i+1) with i="//cp_to_string(i)//&
                " in "//CPSourceFileRef,&
                error,failure)
        end if
     end do
     CPAssert(.not.cp_next(iterator),cp_warning_level,routineP,error,failure)
  end if
  call cp_dealloc(reverseList,error=error)
  if (debug_this_module.and.cp_debug) then
     if (failure) then
        call cp_error_message(cp_warning_level, routineP,&
             "int single linked list tests: tests FAILED",error)
        iterator=>reverseList
        print *,"list="
        do
           if (.not.cp_next(iterator,el_att=el,error=error)) exit
           print *, "  ",el,","
        end do
        print *,"array=",array
     else
        call cp_error_message(cp_warning_level, routineP,&
             "int single linked list tests: SUCESS",error)
     end if
  end if
  res=.not.failure
end function cp_test_int_sll
!!***
!***************************************************************************

! template def put here so that line numbers in template and derived 
! files are almost the same
! [template(el_type,plain_el_type,name_el_type)]
! ARGS:
!  name_el_type = "int"
!  plain_el_type = "integer"
!  el_type = "integer"


end module cp_sll_int_tests
