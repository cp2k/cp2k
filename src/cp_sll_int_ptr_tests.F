!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_sll_int_ptr_tests [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_tests
!!
!!   FUNCTION
!!     do some tests on a linked list of a scalar numeric type (int_ptr)
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
MODULE cp_sll_int_ptr_tests
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_linked_list_int_ptr
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_sll_int_ptr_tests'

  PUBLIC :: cp_test_int_ptr_sll

!!***
!****************************************************************************
CONTAINS

!!****f* cp_sll_int_ptr_tests/cp_test_int_ptr_sll [1.0] *
!!
!!   NAME
!!     cp_test_int_ptr_sll
!!
!!   SYNOPSIS
!!     Function cp_test_int_ptr_sll(glob_env, error) Result(res)
!!       Logical:: res
!!       Type(global_environment_type), Target:: glob_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_test_int_ptr_sll
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
FUNCTION cp_test_int_ptr_sll(glob_env, error) RESULT(res)
  LOGICAL :: res
  TYPE(global_environment_type), TARGET :: glob_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_test_int_ptr_sll',&
       routineP=moduleN//routineN
  INTEGER, PARAMETER :: n=6
  INTEGER, DIMENSION(n) :: array
  INTEGER :: i
  TYPE(cp_sll_int_ptr_type), POINTER :: reverseList, iterator
  INTEGER, POINTER :: el
  failure=.FALSE.
  NULLIFY(reverseList,iterator)
  
  FORALL (i=1:SIZE(array)) array(i)=i*i
  CALL cp_create(reverseList,error=error)
  DO i=1,SIZE(array)
     CALL cp_insert(reverseList,array(i),error=error)
  END DO
  CPAssert(cp_get_length(reverseList)==n,cp_warning_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     iterator => reverseList
     DO i=1,SIZE(array)
        CALL cp_assert(cp_get_element_at(reverseList,i)==array(n-i+1),&
             cp_warning_level,cp_assertion_failed,routineP,&
             "reverse list ("//cp_to_string(i)//")="//cp_to_string(&
             cp_get_element_at(reverseList,i))//" /= "//cp_to_string(array(n-i+1))&
             //"array("//cp_to_string(n-i+1)//") in "//&
             CPSourceFileRef,&
             error,failure)
        IF (.not.cp_next(iterator,el_att=el,error=error)) THEN
           CPPrecondition(.FALSE.,cp_warning_level,routineP,error,failure)
        END IF
        IF (.not.failure) THEN
           CALL cp_assert(el==array(n-i+1),&
                cp_warning_level,cp_assertion_failed,routineP,&
                "iterator /= array(n-i+1) with i="//cp_to_string(i)//&
                " in "//&
                CPSourceFileRef,&
                error,failure)
        END IF
     END DO
     CPAssert(.not.cp_next(iterator),cp_warning_level,routineP,error,failure)
  END IF
  CALL cp_dealloc(reverseList,error=error)
  IF (debug_this_module.and.cp_debug) THEN
     IF (failure) THEN
        CALL cp_error_message(cp_warning_level, routineP,&
             "int_ptr single linked list tests: tests FAILED",error)
        iterator=>reverseList
        PRINT *,"list="
        DO
           IF (.not.cp_next(iterator,el_att=el,error=error)) EXIT
           PRINT *, "  ",el,","
        END DO
        PRINT *,"array=",array
     ELSE
        CALL cp_error_message(cp_warning_level, routineP,&
             "int_ptr single linked list tests: SUCESS",error)
     END IF
  END IF
  res=.not.failure
END FUNCTION cp_test_int_ptr_sll
!!***
!***************************************************************************

! template def put here so that line numbers in template and derived 
! files are almost the same
! [template(el_type,plain_el_type,name_el_type)]
! ARGS:
!  el_type = "integer, pointer"
!  name_el_type = "int_ptr"
!  plain_el_type = "integer"


END MODULE cp_sll_int_ptr_tests
