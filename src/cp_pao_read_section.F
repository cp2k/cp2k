!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_pao_read_section [1.0] *
!!
!!   NAME
!!     cp_pao_read_section
!!
!!   FUNCTION
!!     reads the pao section and initializes the environement with it
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     3.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
module cp_pao_read_section
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use qs_parser, only: start_parser, stop_parser, test_object, read_object
  use string_utilities, only: uppercase
  use cp_linked_list_int
  use cp_pao_types, only: cp_pao_env_type
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_pao_read_section'

  public :: cp_pao_env_init_from_file, cp_test_pao

!***
!****************************************************************************
contains

!!****f* cp_pao_read_section/cp_pao_env_init_from_file [1.0] *
!!
!!   NAME
!!     cp_pao_env_init_from_file
!!
!!   SYNOPSIS
!!     Subroutine cp_pao_env_init_from_file(pao_env, global_env,&
!!         section_label, error)
!!       Type(cp_pao_env_type), Intent (OUT):: pao_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Character(Len=*), Intent (IN):: section_label
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_pao_env_init_from_file
!!
!!   FUNCTION
!!     initialize the pao environement from a section file
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     pao_env: the environement to be initialized
!!     global_env: the global environement
!!     start_section: the string that starts the section
!!     stop_section: the string that stops the section
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     3.2002 created [fawzi]
!!
!!*** **********************************************************************
  subroutine cp_pao_env_init_from_file(pao_env, global_env, section_label,&
       error)
    type(cp_pao_env_type), intent(out) :: pao_env
    type(global_environment_type), intent(in), target :: global_env
    character(len=*), intent(in) :: section_label
    type(cp_error_type), optional, intent(inout) :: error

!!$    logical :: failure, found
!!$    real(kind=wp) :: wanted_subspace_pos_precision,&
!!$         wanted_in_subspace_precision
!!$    character(len=*), parameter :: routineN='cp_pao_env_init_from_file',&
!!$         routineP=moduleN//':'//routineN
!!$    type(cp_sll_int_type), pointer :: min_bas_sll
!!$    character(len=40) :: keyword,label
!!$    character         :: section_character, comment_character
!!$    integer :: min_b_index,i,j
!!$    integer, parameter :: max_kinds=200, max_bas_size=1000
!!$    type(cp_error_type) :: newErr
!!$    type(cp_logger_type) :: fileLogger
!!$    failure=.false.
!!$    wanted_subspace_pos_precision=1.0e-4
!!$    wanted_in_subspace_precision=1.0e-4
!!$
!!$!    call get_special_input_character(comment_character,section_character)
!!$    call start_input_session(global_env%input_file_name,global_env)
!!$    call search(start_section,.true.,found)
!!$    if (found) then
!!$       do
!!$          call read_object(keyword,skip_lines=1)
!!$          call uppercase(keyword)
!!$
!!$          select case (trim(keyword))
!!$          case ("SUBSPACE_POS_PRECISION")
!!$             call read_object(wanted_subspace_pos_precision)
!!$          case ("IN_SUBSPACE_PRECISION")
!!$             call read_object(wanted_in_subspace_precision)
!!$          case ("DO_TESTS")
!!$             call cp_init(fileLogger,filename="pao_tests.log")
!!$             call cp_init(newErr,logger=fileLogger)
!!$             call cp_log(cp_error_get_logger(error),cp_note_level,&
!!$                  routineP," performing pao tests")
!!$             if (.not.cp_test_pao(global_env,newErr)) then
!!$                call cp_error_message(cp_warning_level,routineP,&
!!$                     '+++++ some of the pao tests FAILED ++++',error)
!!$             else
!!$                call cp_error_message(cp_note_level,routineP,&
!!$                     ' +++++ all the pao tests were SUCCESSFULL ++++',error)
!!$             end if
!!$             call cp_dealloc_ref(newErr)
!!$             call cp_dealloc_ref(fileLogger)
!!$          case ("&MIN_BAS")
!!$             do i=1,max_kinds
!!$                call read_object(label,skip_lines=1)
!!$                call uppercase(label)
!!$                if (trim(label)==stop_section) exit
!!$                CPInvariant(label(1:1) /= section_character,cp_failure_level,routineP,error,failure)
!!$                if (failure) exit
!!$                do j=1,max_bas_size
!!$                   if (.not.cp_int_try_read(min_b_index)) exit
!!$                   call cp_insert(min_bas_sll,min_b_index,error=error)
!!$                end do
!!$                call cp_sll_int_output(cp_error_get_logger(error),&
!!$                     'minbas_index',&
!!$                     routineP,i,min_bas_sll,&
!!$                     comment="label="//label,error=error)
!!$                call cp_dealloc(min_bas_sll,error=error)
!!$             end do
!!$          case DEFAULT
!!$             CPPostcondition(keyword == stop_section,cp_failure_level,routineP,error,failure)
!!$             exit
!!$          end select
!!$          
!!$       end do
!!$    else
!!$       call cp_log(cp_error_get_logger(error),cp_warning_level,routineP,&
!!$            "no pao section found")
!!$    end if
!!$    call finish_input_session
  end subroutine cp_pao_env_init_from_file

!!****f* cp_pao_read_section/cp_test_pao [1.0] *
!!
!!   NAME
!!     cp_test_pao
!!
!!   SYNOPSIS
!!     Function cp_test_pao(global_env, error) Result(res)
!!       Use cp_pao_plain_rot_tests, Only: cp_pao_test_plain_rot
!!       Use cp_pao_types_tests, Only: cp_test_pao_types
!!       Use cp_sll_int_tests, Only: cp_test_int_sll
!!       Use cp_sll_int_ptr_tests, Only: cp_test_int_ptr_sll
!!       Use cp_lbfgs_optimizer_tests, Only: cp_lbfgs_e_rosen_test
!!       Type(global_environment_type):: global_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Logical:: res
!!     End Function cp_test_pao
!!
!!   FUNCTION
!!     does all the automatic tests of pao (and more or less related)
!!     routines.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     3.2002 noved from fd_tests to here [fawzi]
!!
!!*** **********************************************************************
function cp_test_pao(global_env,error) result(res)  
  use cp_pao_plain_rot_tests, only: cp_pao_test_plain_rot
  use cp_pao_types_tests, only: cp_test_pao_types
  use cp_sll_int_tests, only: cp_test_int_sll
  use cp_sll_int_ptr_tests, only: cp_test_int_ptr_sll
  use cp_lbfgs_optimizer_tests, only: cp_lbfgs_e_rosen_test

  type(global_environment_type) :: global_env
  type(cp_error_type), optional, intent(inout) :: error
  logical :: res
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_test_pao',&
        routineP=moduleN//':'//routineN
  failure=.false.

  call cp_error_message(cp_note_level,routineP,&
       ' +++++ began pao tests ++++',error)
  call cp_pao_test_plain_rot(global_env,failure,error=error)
  call cp_test_pao_types(global_env,failure,error=error)
  if (.not.cp_test_int_sll(global_env,error)) failure=.true.
  if (.not.cp_test_int_ptr_sll(global_env,error)) failure=.true.
  if (.not.cp_lbfgs_e_rosen_test(global_env,error)) failure=.true.
  if (failure) then
     call cp_error_message(cp_warning_level,routineP,&
          '+++++ some of the pao tests FAILED ++++',error)
  else
     call cp_error_message(cp_note_level,routineP,&
          ' +++++ all the pao tests were SUCCESSFULL ++++',error)
  end if
  res=.not.failure
end function cp_test_pao
!***************************************************************************

end module cp_pao_read_section
