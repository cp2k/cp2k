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
MODULE pao_read_section
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_lbfgs_optimizer_tests,        ONLY: cp_lbfgs_e_rosen_test
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_log,&
                                             cp_logger_create,&
                                             cp_logger_release,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_sll_int_tests,                ONLY: cp_test_int_sll
!use cp_sll_int_ptr_tests, only: cp_test_int_ptr_sll

  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE pao_env_methods,                 ONLY: pao_env_init
  USE pao_plain_rot_tests,             ONLY: pao_test_plain_rot
  USE pao_types,                       ONLY: pao_env_type
  USE pao_types_tests,                 ONLY: cp_test_pao_types
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_parser,                       ONLY: finish_parser,&
                                             read_object,&
                                             start_parser,&
                                             test_object
  USE string_utilities,                ONLY: uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='pao_read_section'

  PUBLIC :: pao_env_init_from_file, cp_test_pao

!***
!****************************************************************************
CONTAINS

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
  SUBROUTINE pao_env_init_from_file(qs_env, global_env, section_label,&
       error)
    TYPE(qs_environment_type), INTENT(inout) :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: section_label
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: default_section_label = "PAO", &
      routineN = 'pao_env_init_from_file', routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: max_bas_size = 1000, &
                                                max_kinds = 200

    CHARACTER(len=3)                         :: obj_kind
    CHARACTER(len=40)                        :: keyword
    LOGICAL                                  :: failure, found
    REAL(KIND=dp)                            :: wanted_in_subspace_precision, &
                                                wanted_subspace_pos_precision
    TYPE(cp_error_type)                      :: newErr
    TYPE(cp_logger_type), POINTER            :: errorLogger, fileLogger
    TYPE(pao_env_type), POINTER              :: pao_env

    failure=.FALSE.
    errorLogger => cp_error_get_logger(error)
    wanted_subspace_pos_precision=1.0d-4
    wanted_in_subspace_precision=1.0d-4

    CALL get_qs_env(qs_env, pao_env=pao_env)
    IF (PRESENT(section_label)) THEN
       CALL start_parser(file_name=global_env%input_file_name,&
            globenv=global_env,&
            start_section_label=section_label,&
            section_found=found)
    ELSE
       CALL start_parser(file_name=global_env%input_file_name,&
            globenv=global_env,&
            start_section_label=default_section_label,&
            section_found=found)
    END IF
    
    IF (found) THEN
       DO
          obj_kind=test_object(newline=.TRUE.)
          SELECT CASE (obj_kind)
          CASE ("EOS")
             EXIT
          CASE ("EOF")
             CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,&
                  routineP,&
                  "file ended before the pao section was ended"//&
                  CPSourceFileRef,&
                  error,failure)
             EXIT
          END SELECT
          
          CALL read_object(keyword,skip_lines=1)
          CALL uppercase(keyword)

          SELECT CASE (TRIM(keyword))
          CASE ("SUBSPACE_POS_PRECISION")
             CALL read_object(wanted_subspace_pos_precision)
          CASE ("IN_SUBSPACE_PRECISION")
             CALL read_object(wanted_in_subspace_precision)
          CASE ("DO_TESTS") ! to do: new logger
             CALL cp_logger_create(fileLogger,global_filename="pao_tests.log",&
                  local_filename="pao_tests.log",&
                  template_logger=errorLogger)
             CALL cp_error_init(newErr,logger=fileLogger)
             CALL cp_log(errorLogger,cp_note_level,&
                  routineP," performing pao tests")
             ! if it crashes while testing use error instread of newErr
             ! to log to sdtout (often output to a file is cached and lost)
             IF (.NOT.cp_test_pao(global_env,newErr)) THEN
                CALL cp_error_message(cp_warning_level,routineP,&
                     '+++++ some of the pao tests FAILED ++++',error)
             ELSE
                CALL cp_error_message(cp_note_level,routineP,&
                     ' +++++ all the pao tests were SUCCESSFULL ++++',error)
             END IF
             CALL cp_error_dealloc_ref(newErr)
             CALL cp_logger_release(fileLogger)
          CASE DEFAULT
             CALL cp_error_message(cp_warning_level,routineP,&
                  "unknown keyword '"//keyword//"' in pao section",&
                  error)
             EXIT
          END SELECT
       END DO

       ! initializes pao_env
       CALL pao_env_init(pao_env, qs_env=qs_env, global_env=global_env,&
            wanted_subspace_pos_precision=wanted_subspace_pos_precision,&
            wanted_in_subspace_precision=wanted_in_subspace_precision,&
            error=error)
    ELSE
       CALL cp_log(errorLogger,cp_warning_level,routineP,&
            "no pao section found")
    END IF
    CALL finish_parser()
  END SUBROUTINE pao_env_init_from_file

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
FUNCTION cp_test_pao(global_env,error) RESULT(res)  
    TYPE(global_environment_type)            :: global_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_test_pao', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CALL cp_error_message(cp_note_level,routineP,&
       ' +++++ began pao tests ++++',error)
  CALL pao_test_plain_rot(global_env,failure,error=error)
  CALL cp_test_pao_types(global_env,failure,error=error)
  IF (.NOT.cp_test_int_sll(global_env,error)) failure=.TRUE.
  !if (.not.cp_test_int_ptr_sll(global_env,error)) failure=.true.
  IF (.NOT.cp_lbfgs_e_rosen_test(global_env,error)) failure=.TRUE.
  IF (failure) THEN
     CALL cp_error_message(cp_warning_level,routineP,&
          '+++++ some of the pao tests FAILED ++++',error)
  ELSE
     CALL cp_error_message(cp_note_level,routineP,&
          ' +++++ all the pao tests were SUCCESSFULL ++++',error)
  END IF
  res=.NOT.failure
END FUNCTION cp_test_pao
!***************************************************************************

END MODULE pao_read_section
