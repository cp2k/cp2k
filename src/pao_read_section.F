!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_read_section [1.0] *
!!
!!   NAME
!!     pao_read_section
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
module pao_read_section
  use cp_log_handling, only: cp_to_string, cp_failure_level, cp_warning_level,&
       cp_log, cp_logger_type, cp_note_level, cp_logger_create, &
       cp_logger_release
  use cp_error_handling, only: cp_error_type, cp_assert, cp_error_message,&
       cp_assertion_failed, cp_debug, cp_precondition_failed, cp_error_init,&
       cp_error_dealloc_ref, cp_error_get_logger, cp_error_get_print_level,&
       cp_internal_error
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use qs_parser, only: start_parser, finish_parser, test_object, read_object
  use qs_environment_types, only: qs_environment_type, get_qs_env
  use string_utilities, only: uppercase
  use cp_linked_list_int, only: cp_sll_int_type, cp_create, cp_dealloc,&
       cp_next, cp_get
  use pao_types, only: pao_env_type
  use pao_bas_projection_types, only: pao_bas_projection_type, cp_get
  use atomic_kind_types, only: atomic_kind_type, get_atomic_kind
  use pao_env_methods, only: pao_env_init
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='pao_read_section'

  public :: pao_env_init_from_file, cp_test_pao

!***
!****************************************************************************
contains

!!****f* pao_read_section/pao_env_init_from_file [1.0] *
!!
!!   NAME
!!     pao_env_init_from_file
!!
!!   SYNOPSIS
!!     Subroutine pao_env_init_from_file(qs_env, global_env, section_label,&
!!         error)
!!       Type(qs_environment_type), Intent (INOUT):: qs_env
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Character(Len=*), Intent (IN), Optional:: section_label
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pao_env_init_from_file
!!
!!   FUNCTION
!!     initialize the pao environement from a section file
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env: the qs_environement that contains the pao environement
!!       to be initialized
!!     - global_env: the global environement
!!     - section_label: the name of the section (defaults to "PAO")
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     3.2002 created [fawzi]
!!
!!*** **********************************************************************
  subroutine pao_env_init_from_file(qs_env, global_env, section_label,&
       error)
    type(qs_environment_type), intent(inout) :: qs_env
    type(global_environment_type), intent(in), target :: global_env
    character(len=*), intent(in), optional :: section_label
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure, found
    real(kind=wp) :: wanted_subspace_pos_precision,&
         wanted_in_subspace_precision
    character(len=*), parameter :: routineN='pao_env_init_from_file',&
         routineP=moduleN//':'//routineN, default_section_label="PAO"
    character(len=40) :: keyword,label
    character(len=3) :: obj_kind
    integer :: min_b_index,i,j
    integer, parameter :: max_kinds=200, max_bas_size=1000
    type(cp_error_type) :: newErr
    type(cp_logger_type) :: fileLogger
    type(pao_env_type), pointer :: pao_env
    type(atomic_kind_type), dimension(:), pointer :: atomic_kind_set
    type(pao_bas_projection_type), pointer :: pao_bas_proj
    character(len=60) :: name
    failure=.false.
    wanted_subspace_pos_precision=1.0e-4
    wanted_in_subspace_precision=1.0e-4

    call get_qs_env(qs_env, pao_env=pao_env)
    if (present(section_label)) then
       call start_parser(file_name=global_env%input_file_name,&
            globenv=global_env,&
            start_section_label=section_label,&
            section_found=found)
    else
       call start_parser(file_name=global_env%input_file_name,&
            globenv=global_env,&
            start_section_label=default_section_label,&
            section_found=found)
    end if
    
    if (found) then
       do
          obj_kind=test_object(newline=.TRUE.)
          select case (obj_kind)
          case ("EOS")
             exit
          case ("EOF")
             call cp_assert(.false.,cp_warning_level,cp_assertion_failed,&
                  routineP,&
                  "file ended before the pao section was ended"//&
                  CPSourceFileRef,&
                  error,failure)
             exit
          end select
          
          call read_object(keyword,skip_lines=1)
          call uppercase(keyword)

          select case (trim(keyword))
          case ("SUBSPACE_POS_PRECISION")
             call read_object(wanted_subspace_pos_precision)
          case ("IN_SUBSPACE_PRECISION")
             call read_object(wanted_in_subspace_precision)
          case ("DO_TESTS") ! to do: new logger
             !call cp_logger_create(fileLogger,filename="pao_tests.log")
             call cp_error_init(newErr)!,logger=fileLogger)
             call cp_log(cp_error_get_logger(error),cp_note_level,&
                  routineP," performing pao tests")
             ! if it crashes while testing use error instread of newErr
             ! to log to sdtout (often output to a file is cached and lost)
             if (.not.cp_test_pao(global_env,newErr)) then
                call cp_error_message(cp_warning_level,routineP,&
                     '+++++ some of the pao tests FAILED ++++',error)
             else
                call cp_error_message(cp_note_level,routineP,&
                     ' +++++ all the pao tests were SUCCESSFULL ++++',error)
             end if
             call cp_error_dealloc_ref(newErr)
             !call cp_logger_dealloc_ref(fileLogger)
          case DEFAULT
             call cp_error_message(cp_warning_level,routineP,&
                  "unknown keyword '"//keyword//"' in pao section",&
                  error)
             exit
          end select
       end do

       ! initializes pao_env
       call pao_env_init(pao_env, qs_env=qs_env, global_env=global_env,&
            wanted_subspace_pos_precision=wanted_subspace_pos_precision,&
            wanted_in_subspace_precision=wanted_in_subspace_precision,&
            error=error)
    else
       call cp_log(cp_error_get_logger(error),cp_warning_level,routineP,&
            "no pao section found")
    end if
    call finish_parser()
  end subroutine pao_env_init_from_file

!!****f* pao_read_section/cp_test_pao [1.0] *
!!
!!   NAME
!!     cp_test_pao
!!
!!   SYNOPSIS
!!     Function cp_test_pao(global_env, error) Result(res)
!!       Use pao_plain_rot_tests, Only: pao_test_plain_rot
!!       Use pao_types_tests, Only: cp_test_pao_types
!!       Use cp_sll_int_tests, Only: cp_test_int_sll
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
!!     - global_env: the global environement
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     3.2002 moved from fd_tests to here [fawzi]
!!
!!*** **********************************************************************
function cp_test_pao(global_env,error) result(res)  
  use pao_plain_rot_tests, only: pao_test_plain_rot
  use pao_types_tests, only: cp_test_pao_types
  use cp_sll_int_tests, only: cp_test_int_sll
  !use cp_sll_int_ptr_tests, only: cp_test_int_ptr_sll
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
  call pao_test_plain_rot(global_env,failure,error=error)
  call cp_test_pao_types(global_env,failure,error=error)
  if (.not.cp_test_int_sll(global_env,error)) failure=.true.
  !if (.not.cp_test_int_ptr_sll(global_env,error)) failure=.true.
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

end module pao_read_section
