!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2004  CP2K developers group
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"

MODULE farming_methods
  USE farming_types,                   ONLY: farming_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length
  USE cp_log_handling,                 ONLY: cp_logger_type, &
                                             cp_failure_level, &
                                             cp_to_string
  USE cp_error_handling,               ONLY: cp_error_type, &
                                             cp_a_l, &
                                             cp_assert, &
                                             cp_error_init, &
                                             cp_error_dealloc_ref, &
                                             cp_assertion_failed 
  USE machine,                         ONLY: m_chdir,m_getcwd
  USE qs_parser,                       ONLY: finish_parser,&
                                             read_object,&
                                             start_parser,&
                                             stop_parser,&
                                             test_object,&
                                             open_file,&
                                             close_file
  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: farming_parse_input

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "farming_runs"

CONTAINS

  SUBROUTINE farming_parse_input(farming_env,globenv,error)
    TYPE(farming_env_type), POINTER          :: farming_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'farming_parse_input', &
      routineP = moduleN//':'//routineN
    CHARACTER(LEN=3)                         :: test_result
    INTEGER :: stat,i 
    LOGICAL :: found,failure
    CHARACTER(LEN=default_string_length) :: dum
    failure=.FALSE.

    ! count the number of jobs first
    farming_env % njobs = 0
    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="JOBS",section_found=found)

    IF (found) THEN
       DO WHILE (test_object(newline=.TRUE.) /= "EOS")
          farming_env % njobs = farming_env % njobs + 1
          CALL read_object(dum,newline=.TRUE.)
       ENDDO
    ELSE
       failure=.true.
       CALL stop_parser(routine=routineP,message="section &JOBS not found in "//globenv%input_file_name)      
    ENDIF
    CALL finish_parser()

    ! get the job info
    ALLOCATE(farming_env%Job_cwd(farming_env % njobs),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(farming_env%Job_input(farming_env % njobs),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
        CALL start_parser(file_name=globenv%input_file_name,&
                          globenv=globenv,&
                          start_section_label="JOBS",section_found=found)
        i=0
        DO WHILE (test_object(newline=.TRUE.) /= "EOS")
           i=i+1
           IF (i<=farming_env % njobs) THEN
              CALL read_object(farming_env%Job_cwd(i),newline=.TRUE.)           
              CALL read_object(farming_env%Job_input(i))           
           ELSE
              failure=.true.
              CALL stop_parser(routine=routineP,message="inconsistent reading in the JOBS section")
           ENDIF
        ENDDO
        CALL finish_parser()
    ENDIF

  END SUBROUTINE

END MODULE farming_methods
!!*****
