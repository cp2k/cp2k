!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief - writing and printing the files, trajectory (pos, cell, dipoles) as 
!>        well as restart files
!>        - usually just the Markov Chain elements are regarded, the elements 
!>        beside this trajectory are neglected
!>        - futrthermore (by option) just the accepted configurations 
!>          are print out to reduce the file sizes
!> \par History
!>      12.2012 created [Mandes Schönherr]
!> \author Mandes
! *****************************************************************************

MODULE tmc_file_io
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE physcon,                         ONLY: au2a => angstrom
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE tmc_analysis_types,              ONLY: tmc_analysis_env
  USE tmc_calculations,                ONLY: get_cell_scaling,&
                                             get_scaled_cell
  USE tmc_move_types,                  ONLY: nr_mv_types
  USE tmc_stati,                       ONLY: TMC_STATUS_FAILED,&
                                             TMC_STATUS_OK,&
                                             TMC_STATUS_WAIT_FOR_NEW_TASK,&
                                             tmc_default_restart_in_file_name,&
                                             tmc_default_restart_out_file_name,&
                                             tmc_default_trajectory_file_name
  USE tmc_tree_types,                  ONLY: elem_array_type,&
                                             tree_type
  USE tmc_types,                       ONLY: tmc_env_type,&
                                             tmc_param_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tmc_file_io'

  ! filename manipulation
  PUBLIC :: expand_file_name_char, expand_file_name_temp, expand_file_name_int
  ! read/write restart file
  PUBLIC :: print_restart_file, read_restart_file
  ! write the configuration
  PUBLIC :: write_result_list_element
  PUBLIC :: write_element_in_file
  PUBLIC :: write_dipoles_in_file
  ! analysis read
  PUBLIC :: analyse_files_open, read_element_from_file, analyse_files_close

CONTAINS

!------------------------------------------------------------------------------
! routines for manipulating the file name
!------------------------------------------------------------------------------
! *****************************************************************************
!> \brief placing a character string at the end of a file name 
!>        (instead of the ending) 
!> \param file_name original file name 
!> \param extra string to be added before the file extention
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param result_file_name the new filename
!> \author Mandes 11.2012
! *****************************************************************************
  FUNCTION expand_file_name_ending(file_name, extra, error) RESULT(result_file_name)
    CHARACTER(LEN=*)                         :: file_name, extra
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(LEN=default_path_length)       :: result_file_name

    CHARACTER(LEN=*), PARAMETER :: routineN = 'expand_file_name_ending', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ind
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(file_name.NE."",cp_failure_level,routineP,error,failure)

    ind = INDEX(file_name,".",BACK = .TRUE.)
    IF(.NOT.ind.EQ.0) THEN
       WRITE(result_file_name,*)file_name(1:ind-1),".",&
             TRIM(ADJUSTL(extra))
    ELSE
       WRITE(result_file_name,*)TRIM(file_name),".",extra
    END IF
    result_file_name=TRIM(ADJUSTL(result_file_name))
    CPPostcondition(result_file_name.NE."",cp_failure_level,routineP,error,failure)
  END FUNCTION expand_file_name_ending

! *****************************************************************************
!> \brief placing a character string at the end of a file name 
!>        (before the file extention) 
!> \param file_name original file name 
!> \param extra string to be added before the file extention
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param result_file_name the new filename
!> \author Mandes 11.2012
! *****************************************************************************
  FUNCTION expand_file_name_char(file_name, extra, error) RESULT(result_file_name)
    CHARACTER(LEN=*)                         :: file_name, extra
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(LEN=default_path_length)       :: result_file_name

    CHARACTER(LEN=*), PARAMETER :: routineN = 'expand_file_name_char', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ind
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(file_name.NE."",cp_failure_level,routineP,error,failure)

    ind = INDEX(file_name,".",BACK = .TRUE.)
    IF(.NOT.ind.EQ.0) THEN
       WRITE(result_file_name,*)file_name(1:ind-1),"_",&
             TRIM(ADJUSTL(extra)),file_name(ind:LEN_TRIM(file_name))
    ELSE
       WRITE(result_file_name,*)TRIM(file_name),"_",extra
    END IF
    result_file_name=TRIM(ADJUSTL(result_file_name))
    CPPostcondition(result_file_name.NE."",cp_failure_level,routineP,error,failure)
  END FUNCTION expand_file_name_char

! *****************************************************************************
!> \brief placing the temperature at the end of a file name 
!>        (before the file extention) 
!> \param file_name original file name 
!> \param rvalue temperature to be added
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param result_file_name the new filename
!> \author Mandes 11.2012
! *****************************************************************************
  FUNCTION expand_file_name_temp(file_name,rvalue, error) RESULT(result_file_name)
    CHARACTER(LEN=*)                         :: file_name
    REAL(KIND=dp)                            :: rvalue
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(LEN=default_path_length)       :: result_file_name

    CHARACTER(LEN=*), PARAMETER :: routineN = 'expand_file_name_temp', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=18)                        :: rval_to_string
    INTEGER                                  :: ind
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(file_name.NE."",cp_failure_level,routineP,error,failure)

    rval_to_string = ""

    WRITE(rval_to_string,"(F16.2)")rvalue
    ind = INDEX(file_name,".",BACK = .TRUE.)
    IF(.NOT.ind.EQ.0) THEN
       WRITE(result_file_name,*)file_name(1:ind-1),"_T", &
              TRIM(ADJUSTL(rval_to_string)),file_name(ind:LEN_TRIM(file_name))
    ELSE
       IF(LEN(file_name).EQ.0)THEN
          WRITE(result_file_name,*)TRIM(file_name),"T", TRIM(ADJUSTL(rval_to_string)),&
            file_name(ind:LEN_TRIM(file_name))
       ELSE
          WRITE(result_file_name,*)TRIM(file_name),"_T", TRIM(ADJUSTL(rval_to_string))
       END IF
    END IF
    result_file_name=TRIM(ADJUSTL(result_file_name))
    CPPostcondition(result_file_name.NE."",cp_failure_level,routineP,error,failure)
  END FUNCTION expand_file_name_temp

! *****************************************************************************
!> \brief placing an integer at the end of a file name 
!>        (before the file extention) 
!> \param file_name original file name 
!> \param ivalue number to be added
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param result_file_name the new filename
!> \author Mandes 11.2012
! *****************************************************************************
  FUNCTION expand_file_name_int(file_name,ivalue, error) RESULT(result_file_name)
    CHARACTER(LEN=*)                         :: file_name
    INTEGER                                  :: ivalue
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(LEN=default_path_length)       :: result_file_name

    CHARACTER(LEN=*), PARAMETER :: routineN = 'expand_file_name_int', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=18)                        :: rval_to_string
    INTEGER                                  :: ind
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(file_name.NE."",cp_failure_level,routineP,error,failure)

    rval_to_string = ""

    WRITE(rval_to_string,*)ivalue
    ind = INDEX(file_name,".",BACK = .TRUE.)
    IF(.NOT.ind.EQ.0) THEN
       WRITE(result_file_name,*)file_name(1:ind-1),"_", &
              TRIM(ADJUSTL(rval_to_string)),file_name(ind:LEN_TRIM(file_name))
    ELSE
       IF(LEN(file_name).EQ.0)THEN
          WRITE(result_file_name,*)TRIM(file_name),"", TRIM(ADJUSTL(rval_to_string)),&
            file_name(ind:LEN_TRIM(file_name))
       ELSE
          WRITE(result_file_name,*)TRIM(file_name),"_", TRIM(ADJUSTL(rval_to_string)),&
            file_name(ind:LEN_TRIM(file_name))
       END IF
    END IF
    result_file_name=TRIM(ADJUSTL(result_file_name))
    CPPostcondition(result_file_name.NE."",cp_failure_level,routineP,error,failure)
  END FUNCTION expand_file_name_int

!------------------------------------------------------------------------------
! routines for reading and writing RESTART file
!------------------------------------------------------------------------------
! *****************************************************************************
!> \brief prints out the TMC restart files with all last configurations and 
!>        counters etc.
!> \param tmc_env ...
!> \param job_counts the counters for counting the submitted different job types
!> \param timings ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param tmv_env the tmc environment, storing result lists and counters an in 
!>        temperatures UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE print_restart_file(tmc_env, job_counts, timings, error)
    TYPE(tmc_env_type), POINTER              :: tmc_env
    INTEGER, DIMENSION(:)                    :: job_counts
    REAL(KIND=dp), DIMENSION(4)              :: timings
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_restart_file', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: c_tmp, file_name
    INTEGER                                  :: f_unit, i
    LOGICAL                                  :: failure

    failure = .FALSE.
    c_tmp=""
    CPPrecondition(ASSOCIATED(tmc_env),cp_failure_level,routineP,error, failure)
    CPPrecondition(ASSOCIATED(tmc_env%m_env),cp_failure_level,routineP,error, failure)
    CPPrecondition(ASSOCIATED(tmc_env%params),cp_failure_level,routineP,error, failure)
    CPPrecondition(ASSOCIATED(tmc_env%m_env%gt_act),cp_failure_level,routineP,error, failure)

    IF(failure) RETURN

    WRITE(c_tmp, FMT='(I9.9)') tmc_env%m_env%result_count(0)
    file_name=TRIM(expand_file_name_char(&
                        file_name=tmc_default_restart_out_file_name, &
                        extra=c_tmp, error=error)) 
    CALL open_file(file_name=file_name, file_status="REPLACE", &
                   file_action="WRITE", file_form="UNFORMATTED", &
                   unit_number=f_unit)
    WRITE(f_unit) SIZE(tmc_env%params%Temp)
    WRITE(f_unit) tmc_env%params%Temp(:), &
              tmc_env%m_env%gt_act%nr, &
              tmc_env%m_env%gt_act%rng_seed, &
              tmc_env%m_env%gt_act%rnd_nr, &
              tmc_env%m_env%gt_act%prob_acc, &
              tmc_env%m_env%gt_act%mv_conf, &
              tmc_env%m_env%gt_act%mv_next_conf, &
              tmc_env%m_env%result_count(0:), &
              tmc_env%params%move_types%mv_weight, &
              tmc_env%params%move_types%acc_count, &
              tmc_env%params%move_types%mv_count, &
              tmc_env%params%move_types%subbox_acc_count, &
              tmc_env%params%move_types%subbox_count, &
              tmc_env%params%cell%hmat, &
              job_counts, &
              timings
    DO i=1, SIZE(tmc_env%params%Temp)
       WRITE(f_unit) tmc_env%m_env%result_list(i)%elem%nr, &
                 tmc_env%m_env%result_list(i)%elem%rng_seed, &
                 tmc_env%m_env%result_list(i)%elem%pos, &
                 tmc_env%m_env%result_list(i)%elem%vel, &
                 tmc_env%m_env%result_list(i)%elem%box_scale, &
                 tmc_env%m_env%result_list(i)%elem%potential, &
                 tmc_env%m_env%result_list(i)%elem%e_pot_approx, &
                 tmc_env%m_env%result_list(i)%elem%ekin, &
                 tmc_env%m_env%result_list(i)%elem%ekin_before_md, &
                 tmc_env%m_env%result_list(i)%elem%temp_created
    END DO
    CALL close_file(unit_number=f_unit)
    ! write the file, where the restart file name is written in
    CALL open_file(file_name=tmc_default_restart_in_file_name, &
                   file_action="WRITE", file_status="REPLACE", &
                   unit_number=f_unit)
    WRITE(f_unit,*)TRIM(file_name)
    CALL close_file(unit_number=f_unit)
  END SUBROUTINE print_restart_file

! *****************************************************************************
!> \brief reads the TMC restart file with all last configurations and 
!>        counters etc.
!> \param tmc_env ...
!> \param job_counts the counters for counting the submitted different job types
!> \param timings ...
!> \param file_name the restart file name
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param tmv_env the tmc environment, storing result lists and counters an in 
!>        temperatures UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE read_restart_file(tmc_env, job_counts, timings, file_name, error)
    TYPE(tmc_env_type), POINTER              :: tmc_env
    INTEGER, DIMENSION(:)                    :: job_counts
    REAL(KIND=dp), DIMENSION(4)              :: timings
    CHARACTER(LEN=*)                         :: file_name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_restart_file', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: file_ptr, i, temp_size
    LOGICAL                                  :: failure, flag
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: tmp_temp
    REAL(KIND=dp), DIMENSION(nr_mv_types)    :: mv_weight_tmp

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(tmc_env),cp_failure_level,routineP,error, failure)
    CPPrecondition(ASSOCIATED(tmc_env%m_env),cp_failure_level,routineP,error, failure)
    CPPrecondition(ASSOCIATED(tmc_env%params),cp_failure_level,routineP,error, failure)
    CPPrecondition(ASSOCIATED(tmc_env%m_env%gt_act),cp_failure_level,routineP,error, failure)
    IF(failure) RETURN

    IF(file_name.EQ.tmc_default_restart_in_file_name) THEN
      INQUIRE(FILE=tmc_default_restart_in_file_name, EXIST=flag)
      CPPrecondition(flag,cp_failure_level,routineP,error, failure)
      CALL open_file(file_name=tmc_default_restart_in_file_name, file_status="OLD", &
                   file_action="READ", unit_number=file_ptr)
      READ(file_ptr,*) file_name
      CALL close_file(unit_number=file_ptr)
    END IF

    CALL open_file(file_name=file_name, file_status="OLD", file_form="UNFORMATTED", &
                   file_action="READ", unit_number=file_ptr)
    READ(file_ptr) temp_size
    CALL cp_assert(temp_size.EQ. SIZE(tmc_env%params%Temp),&
                   cp_failure_level,cp_assertion_failed,routineP,&
                   "the actual specified temperatures does not "//&
                   "fit in amount with the one from restart file ",&
                   error, failure=flag)
    ALLOCATE(tmp_temp(temp_size))
    READ(file_ptr) tmp_temp(:), &
                   tmc_env%m_env%gt_act%nr, &
                   tmc_env%m_env%gt_act%rng_seed, &
                   tmc_env%m_env%gt_act%rnd_nr, &
                   tmc_env%m_env%gt_act%prob_acc, &
                   tmc_env%m_env%gt_act%mv_conf, &         !
                   tmc_env%m_env%gt_act%mv_next_conf, &    !
                   tmc_env%m_env%result_count(0:), &
                   mv_weight_tmp, &                    !
                   tmc_env%params%move_types%acc_count, &
                   tmc_env%params%move_types%mv_count, &
                   tmc_env%params%move_types%subbox_acc_count, &
                   tmc_env%params%move_types%subbox_count, &          !
                   tmc_env%params%cell%hmat, &
                   job_counts, &
                   timings

    CALL cp_assert(ALL(ABS(tmc_env%params%Temp(:)-tmp_temp(:)) .LT. 0.005),&
                   cp_failure_level,cp_assertion_failed,routineP,&
                   "the temperatures differ from the previous calculation. There were the following temperatures used:",&
                   error, failure=flag)
    CALL cp_assert(ALL(mv_weight_tmp(:).EQ.tmc_env%params%move_types%mv_weight(:)),&
                   cp_warning_level,cp_assertion_failed,routineP,&
                   "The amount of mv types differs between the original and the restart run. ",&
                   error, failure=flag)

    DO i=1, SIZE(tmc_env%params%Temp)
       tmc_env%m_env%gt_act%conf(i)%elem => tmc_env%m_env%result_list(i)%elem
       READ(file_ptr) tmc_env%m_env%result_list(i)%elem%nr, &
                      tmc_env%m_env%result_list(i)%elem%rng_seed, &
                      tmc_env%m_env%result_list(i)%elem%pos, &
                      tmc_env%m_env%result_list(i)%elem%vel, &
                      tmc_env%m_env%result_list(i)%elem%box_scale, &
                      tmc_env%m_env%result_list(i)%elem%potential, &
                      tmc_env%m_env%result_list(i)%elem%e_pot_approx, &
                      tmc_env%m_env%result_list(i)%elem%ekin, &
                      tmc_env%m_env%result_list(i)%elem%ekin_before_md, &
                      tmc_env%m_env%result_list(i)%elem%temp_created
    END DO
    CALL close_file(unit_number=file_ptr)
  END SUBROUTINE read_restart_file

  !----------------------------------------------------------------------------
  ! printing configuration in file
  !----------------------------------------------------------------------------

! *****************************************************************************
!> \brief select the correct configuration to print out the 
!>        (coordinates, forces, cell ...)
!> \param result_list list of configurations for each temperature
!> \param result_count list with number of Markov Chain number 
!>          for each teperature (index 0 for global tree)
!> \param conf_updated index of the updated (modified element)
!> \param accepted acceptance flag
!> \param tmc_params TMC environment parameters
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE write_result_list_element(result_list, result_count,conf_updated,&
                                       accepted, tmc_params, error)
    TYPE(elem_array_type), DIMENSION(:), &
      POINTER                                :: result_list
    INTEGER, DIMENSION(:), POINTER           :: result_count
    INTEGER                                  :: conf_updated
    LOGICAL, INTENT(IN)                      :: accepted
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_result_list_element', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: file_name
    INTEGER                                  :: handle, i
    LOGICAL                                  :: failure

    failure   = .FALSE.
    file_name = ""

    CPPrecondition(ASSOCIATED(result_list),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(result_count),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params%Temp),cp_failure_level,routineP,error,failure)
    CPPrecondition(conf_updated.GE.0,cp_failure_level,routineP,error,failure)
    CPPrecondition(conf_updated.LE.SIZE(tmc_params%Temp),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    IF(conf_updated.EQ.0)THEN
      ! for debugging print every configuration of every temperature
      DO i=1, SIZE(tmc_params%Temp)
        WRITE(file_name,*)"every_step_",TRIM(tmc_default_trajectory_file_name)
        CALL write_element_in_file(elem=result_list(i)%elem, &
               tmc_params=tmc_params, conf_nr=result_count(0),&
               file_name=expand_file_name_temp(file_name=file_name,rvalue=tmc_params%Temp(i), error=error), &
               error=error)
      END DO
    ELSE
      IF((.NOT.tmc_params%print_only_diff_conf).OR.&
         (tmc_params%print_only_diff_conf.AND.accepted))THEN
        CALL write_element_in_file(elem=result_list(conf_updated)%elem, &
               tmc_params=tmc_params, conf_nr=result_count(conf_updated),&
               file_name=expand_file_name_temp(file_name=TRIM(tmc_default_trajectory_file_name),&
                           rvalue=tmc_params%Temp(conf_updated), error=error),&
               error=error)
      END IF
    END IF
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE write_result_list_element

! *****************************************************************************
!> \brief writes the trajectory element in a file from sub tree element
!> \param elem actual tree element to be printed out
!> \param tmc_params TMC environment parameters
!> \param temp_index ...
!> \param accepted ...
!> \param file_name ...
!> \param conf_nr Markov chain element number
!> \param conf_info whole header line
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param filename file name will be extended by type of file (pos, cell,...) UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE write_element_in_file(elem, tmc_params, temp_index, accepted, file_name, conf_nr, &
                                   conf_info, error)
    TYPE(tree_type), POINTER                 :: elem
    TYPE(tmc_param_type), POINTER            :: tmc_params
    INTEGER, OPTIONAL                        :: temp_index
    LOGICAL, OPTIONAL                        :: accepted
    CHARACTER(LEN=*), OPTIONAL               :: file_name
    INTEGER, OPTIONAL                        :: conf_nr
    CHARACTER(LEN=*), OPTIONAL               :: conf_info
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_element_in_file', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: file_name_act, tmp_name
    CHARACTER(LEN=default_string_length)     :: header
    INTEGER                                  :: file_ptr, handle, i, nr_atoms
    LOGICAL                                  :: failure, file_exists, print_it
    REAL(KIND=dp)                            :: vol
    REAL(KIND=dp), DIMENSION(3, 3)           :: hmat_scaled

    failure   = .FALSE.
    file_name_act = ""
    tmp_name  = ""
    header    = ""
    print_it  = .TRUE.

    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params%atoms),cp_failure_level,routineP,error,failure)
    CPPrecondition(PRESENT(conf_nr).OR.PRESENT(conf_info),cp_failure_level,routineP,error,failure)

    IF(print_it) THEN
      ! start the timing
      CALL timeset(routineN,handle)
  
      ! set default file name
      IF(PRESENT(file_name)) THEN
        CPPrecondition(file_name.NE."",cp_failure_level,routineP,error,failure)
        file_name_act = file_name
      ELSE
        CPPrecondition(ASSOCIATED(tmc_params%Temp),cp_failure_level,routineP,error,failure)
        CPPrecondition(PRESENT(temp_index),cp_failure_level,routineP,error,failure)
        file_name_act = expand_file_name_temp(file_name=tmc_default_trajectory_file_name, &
                                          rvalue=tmc_params%Temp(temp_index), error=error)
      END IF
  
      nr_atoms=SIZE(elem%pos)/tmc_params%dim_per_elem
 
      ! set header (for coordinate or force file)
      IF(tmc_params%print_trajectory.OR.tmc_params%print_forces) THEN
        IF(PRESENT(conf_info)) THEN
          WRITE(header,*) TRIM(ADJUSTL(conf_info))
        ELSE
          !WRITE(header,FMT="(A,I8,A,F20.10)") " i = ", conf_nr,", E = ", elem%potential
          WRITE(header,FMT="(A,I8,A,F20.10,F20.10,A,I8,I8)")"i =",conf_nr," ,E =",&
            elem%potential, elem%ekin," st elem",elem%sub_tree_nr,elem%nr
        END IF
      END IF
  
      ! write the coordinates
      IF(tmc_params%print_trajectory) THEN
        tmp_name = expand_file_name_ending(file_name_act,"xyz",error)
        CALL open_file(file_name=tmp_name, file_status="UNKNOWN", &
                       file_action="WRITE", file_position="APPEND", &
                       unit_number=file_ptr)
        WRITE(file_ptr,FMT="(I8)") nr_atoms
        WRITE(file_ptr,*) TRIM(header)
        DO i=1, SIZE(elem%pos), tmc_params%dim_per_elem
           WRITE(file_ptr,FMT="(A4,1X,1000F20.10)") &
             TRIM(tmc_params%atoms((i-1)/tmc_params%dim_per_elem+1)%name), &
             elem%pos(i:i+tmc_params%dim_per_elem-1)*au2a
        END DO
        CALL close_file(unit_number=file_ptr)
      END IF
  
      ! write the forces
      IF(tmc_params%print_forces) THEN
         tmp_name = expand_file_name_ending(file_name_act,"frc",error)
         CALL open_file(file_name=tmp_name, file_status="UNKNOWN", &
                       file_action="WRITE", file_position="APPEND", &
                       unit_number=file_ptr)
         WRITE(file_ptr,FMT="(I8)") nr_atoms
         WRITE(file_ptr,*) TRIM(header)
         DO i=1, SIZE(elem%pos), tmc_params%dim_per_elem
            WRITE(file_ptr,FMT="(A4,1X,1000F20.10)") &
              TRIM(tmc_params%atoms((i-1)/tmc_params%dim_per_elem+1)%name), &
              elem%frc(i:i+tmc_params%dim_per_elem-1)
         END DO
         CALL close_file(unit_number=file_ptr)
      END IF

      ! write the cell dipoles
      IF(tmc_params%print_dipole) THEN
         CALL write_dipoles_in_file(file_name=file_name_act, &
                conf_nr=conf_nr, dip=elem%dipole, error=error)
      END IF

      ! write the cell file
      IF(tmc_params%print_cell) THEN
         tmp_name = expand_file_name_ending(file_name_act,"cell",error)
         ! header
         INQUIRE(FILE=tmp_name, EXIST=file_exists)   ! file_exists will be TRUE if the file exist
         IF(.NOT.file_exists) THEN
            CALL open_file(file_name=tmp_name, file_status="NEW", &
                           file_action="WRITE", unit_number=file_ptr)
            WRITE (file_ptr,FMT='(A,9(7X,A2," [Angstrom]"),6X,A)')&
               "# MC step ","Ax","Ay","Az","Bx","By","Bz","Cx","Cy","Cz",&
               "Volume [Angstrom^3]"
         ELSE
            CALL open_file(file_name=tmp_name, file_status="OLD", &
                           file_action="WRITE", file_position="APPEND", &
                           unit_number=file_ptr)
         END IF
         CALL get_scaled_cell(cell=tmc_params%cell, &
                box_scale=elem%box_scale, scaled_hmat=hmat_scaled, &
                vol=vol, error=error)
         WRITE (file_ptr,FMT="(I8,9(1X,F19.10),1X,F24.10)") conf_nr, &
           hmat_scaled(:,:)*au2a, vol* au2a**3 
           !TODO better cell output e.g. using cell_types routine
         CALL close_file(unit_number=file_ptr)
      END IF

      ! write the different energies
      IF(tmc_params%print_energies) THEN
         tmp_name = expand_file_name_ending(file_name_act,"ener",error)
         ! header
         INQUIRE(FILE=tmp_name, EXIST=file_exists)   ! file_exists will be TRUE if the file exist
         IF(.NOT.file_exists) THEN
            CALL open_file(file_name=tmp_name, file_status="NEW", &
                           file_action="WRITE", unit_number=file_ptr)
            WRITE (file_ptr,FMT='(A,4A20)')&
               "# MC step "," exact "," approx "," last SCF "," kinetic "
         ELSE
           CALL open_file(file_name=tmp_name, file_status="OLD", &
                          file_action="WRITE", file_position="APPEND", &
                          unit_number=file_ptr)
         END IF
         WRITE(file_ptr,FMT="(I8,14F20.10)")conf_nr, elem%potential, elem%e_pot_approx, &
           elem%scf_energies(MOD(elem%scf_energies_count,4)+1), elem%ekin
         CALL close_file(unit_number=file_ptr)
      END IF

      ! end the timing
      CALL timestop(handle)
    END IF
  END SUBROUTINE write_element_in_file

! *****************************************************************************
!> \brief writes the cell dipoles in dipole trajectory file
!> \param file_name ...
!> \param conf_nr ...
!> \param dip ...
!> \param file_ext ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE write_dipoles_in_file(file_name, conf_nr, dip, file_ext, error)
    CHARACTER(LEN=default_path_length)       :: file_name
    INTEGER                                  :: conf_nr
    REAL(KIND=dp), DIMENSION(:), POINTER     :: dip
    CHARACTER(LEN=*), INTENT(in), OPTIONAL   :: file_ext
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_dipoles_in_file', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: file_name_tmp
    INTEGER                                  :: file_ptr
    LOGICAL                                  :: failure, file_exists

    CPPrecondition(ASSOCIATED(dip),cp_failure_level,routineP,error,failure)

    IF(PRESENT(file_ext)) THEN
      CPPrecondition(file_ext.NE."",cp_failure_level,routineP,error,failure)
      file_name_tmp = expand_file_name_ending(file_name,TRIM(file_ext),error)
    ELSE
      file_name_tmp = expand_file_name_ending(file_name,"dip",error)
    END IF
    INQUIRE(FILE=file_name_tmp, EXIST=file_exists)
    IF(.NOT.file_exists) THEN
      CALL open_file(file_name=file_name_tmp, file_status="NEW", &
                     file_action="WRITE", unit_number=file_ptr)
      WRITE(file_ptr, FMT='(A8,10A20)')"# conf_nr", "dip_x [C Angstrom]", &
                                 "dip_y [C Angstrom]",  "dip_z [C Angstrom]"
    ELSE
      CALL open_file(file_name=file_name_tmp, file_status="OLD", &
                     file_action="WRITE", file_position="APPEND", &
                     unit_number=file_ptr)
    END IF
    WRITE(file_ptr,FMT="(I8,10F20.10)")conf_nr, dip(:)
    CALL close_file(unit_number=file_ptr)
  END SUBROUTINE write_dipoles_in_file

  !----------------------------------------------------------------------------
  ! read configuration from file
  !----------------------------------------------------------------------------

! *****************************************************************************
!> \brief read the trajectory element from a file from sub tree element
!> \param elem actual tree element to be printed out
!> \param tmc_ana TMC analysis environment parameters
!> \param conf_nr Markov chain element number 
!>        (input the old number and read only if conf nr from file is greater
!> \param stat ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 03.2013
! *****************************************************************************
  SUBROUTINE read_element_from_file(elem, tmc_ana, conf_nr, stat, error)
    TYPE(tree_type), POINTER                 :: elem
    TYPE(tmc_analysis_env), POINTER          :: tmc_ana
    INTEGER                                  :: conf_nr, stat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_element_from_file', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: conf_nr_old, handle, i_tmp
    LOGICAL                                  :: failure, &
                                                files_conf_missmatch, flag

    failure   = .FALSE.
    stat = TMC_STATUS_OK
    conf_nr_old = conf_nr
    files_conf_missmatch = .FALSE.

    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_ana),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_ana%atoms),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    ! read the coordinates
    IF(tmc_ana%id_traj.GT.0) THEN
      i_tmp = conf_nr_old
      CALL read_pos_from_file(elem=elem, tmc_ana=tmc_ana, stat=stat, &
                              conf_nr=i_tmp, error=error)
      IF(stat.EQ.TMC_STATUS_WAIT_FOR_NEW_TASK) THEN
        CALL cp_assert(.FALSE.,&
             cp_warning_level,cp_assertion_failed,routineP,&
             'end of position file reached at line '//&
             cp_to_string(REAL(tmc_ana%lc_traj,KIND=dp))//", last element "//&
             cp_to_string(tmc_ana%last_elem%nr),&
             error, failure=flag, only_ionode=.TRUE.)
      ELSE
        CPPostcondition(i_tmp.GT.conf_nr_old,cp_failure_level,routineP,error,failure)
        conf_nr = i_tmp
        elem%nr = i_tmp
      END IF
    END IF

    ! read the forces
    ! TODO if necessary

    ! read the dipoles file
    IF(tmc_ana%id_dip.GT.0 .AND. stat.EQ.TMC_STATUS_OK) THEN
      i_tmp = conf_nr_old
      search_conf_dip: DO
        CALL read_dipole_from_file(elem=elem, tmc_ana=tmc_ana, stat=stat, &
                                   conf_nr=i_tmp, error=error)
        IF(stat.EQ.TMC_STATUS_WAIT_FOR_NEW_TASK) THEN 
          CALL cp_assert(.FALSE.,&
                 cp_warning_level,cp_assertion_failed,routineP,&
                 'end of dipole file reached at line'//&
                 cp_to_string(REAL(tmc_ana%lc_dip,KIND=dp)),&
                 error, failure=flag, only_ionode=.TRUE.)
          EXIT search_conf_dip
        END IF
        ! check consitence with pos file
        IF(tmc_ana%id_traj.GT.0) THEN
          IF(i_tmp .EQ. conf_nr) THEN
            files_conf_missmatch = .FALSE.
            EXIT search_conf_dip
          ELSE
            ! the configuration numbering differ from the position file,
            !  but we keep on searching for the correct configuration
            files_conf_missmatch = .TRUE.
          END IF
        ! if no pos file, just take the next conf
        ELSE IF(i_tmp.GT.conf_nr_old) THEN
          conf_nr = i_tmp
          elem%nr = i_tmp
          EXIT search_conf_dip
        END IF
      END DO search_conf_dip
    END IF

    ! read the cell file
    IF(tmc_ana%id_cell.GT.0 .AND. stat.EQ.TMC_STATUS_OK) THEN
      search_conf_cell: DO
        CALL read_cell_from_file(elem=elem, tmc_ana=tmc_ana, stat=stat, &
                                 conf_nr=i_tmp, error=error)
        IF(stat.EQ.TMC_STATUS_WAIT_FOR_NEW_TASK) THEN
          CALL cp_assert(.FALSE.,&
                 cp_warning_level,cp_assertion_failed,routineP,&
                 'end of cell file reached at line at line'//&
                 cp_to_string(REAL(tmc_ana%lc_cell, KIND=dp)),&
                 error, failure=flag, only_ionode=.TRUE.)
          EXIT search_conf_cell
        END IF
        ! check consitence with pos file
        IF(tmc_ana%id_traj.GT.0) THEN
          IF(i_tmp .EQ. conf_nr) THEN
            files_conf_missmatch = .FALSE.
            EXIT search_conf_cell
          ELSE 
            ! the configuration numbering differ from the position file,
            !  but we keep on searching for the correct configuration
            files_conf_missmatch = .TRUE.
          END IF
        ! if no pos file, just take the next conf
        ELSE IF(i_tmp.GT.conf_nr_old) THEN
          conf_nr = i_tmp
          elem%nr = i_tmp
          EXIT search_conf_cell
        END IF
      END DO search_conf_cell
      
    END IF

    ! write the different energies
    ! TODO if necessary

    CALL cp_assert(.NOT.files_conf_missmatch,&
           cp_warning_level,cp_assertion_failed,routineP,&
           'there is a missmatch in the configuration numbering. '//&
           "Read number of lines (pos|cell|dip)"//&
           cp_to_string(tmc_ana%lc_traj)//"|"//&
           cp_to_string(tmc_ana%lc_cell)//"|"//&
           cp_to_string(tmc_ana%lc_dip), &
           error, failure=flag, only_ionode=.TRUE.)

    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE read_element_from_file

! *****************************************************************************
!> \brief search for the next configurational position in file
!> \param elem actual tree element to be read
!> \param tmc_ana ...
!> \param stat ...
!> \param conf_nr Markov chain element number 
!>        (input the old number and read only if conf nr from file is greater
!> \param header_info ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 03.2013
! *****************************************************************************
  SUBROUTINE read_pos_from_file(elem, tmc_ana, stat, conf_nr, header_info, error)
    TYPE(tree_type), POINTER                 :: elem
    TYPE(tmc_analysis_env), POINTER          :: tmc_ana
    INTEGER                                  :: stat, conf_nr
    CHARACTER(LEN=*), OPTIONAL               :: header_info
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_pos_from_file', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: c_tmp
    INTEGER                                  :: handle, i, i_tmp, status
    LOGICAL                                  :: failure

    failure = .FALSE.
    stat = TMC_STATUS_FAILED

    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(elem%pos),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_ana),cp_failure_level,routineP,error,failure)
    CPPrecondition(tmc_ana%id_traj.GT.0,cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    search_next_conf:DO
      c_tmp(:) = " "
      READ(tmc_ana%id_traj,'(A)',IOSTAT=status) c_tmp(:)
      tmc_ana%lc_traj = tmc_ana%lc_traj +1
      CALL cp_assert(status.LE.0,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "configuration header read error at line: "//&
               cp_to_string(tmc_ana%lc_traj)//": "//c_tmp, error)
      IF(status .LT. 0) THEN ! end of file reached
        stat = TMC_STATUS_WAIT_FOR_NEW_TASK
        EXIT search_next_conf
      END IF
      IF(INDEX(c_tmp,"=").GT.0)THEN
        READ(c_tmp(INDEX(c_tmp,"=")+1:),*,IOSTAT=status) i_tmp ! read the configuration number
        CALL cp_assert(status.EQ.0,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "configuration header read error (for conf nr) at line: "//&
               cp_to_string(tmc_ana%lc_traj), error)
        IF(i_tmp .GT. conf_nr) THEN
          ! TODO we could also read the energy ...
          conf_nr = i_tmp
          IF(PRESENT(header_info)) header_info = c_tmp
          stat = TMC_STATUS_OK
          EXIT search_next_conf
        END IF
      END IF
    END DO search_next_conf

    IF(stat.EQ. TMC_STATUS_OK) THEN
      pos_loop:DO i=1, SIZE(elem%pos), tmc_ana%dim_per_elem
        READ(tmc_ana%id_traj,FMT="(A4,1X,1000F20.10)", IOSTAT=status) &
             c_tmp, elem%pos(i:i+tmc_ana%dim_per_elem-1)
        tmc_ana%lc_traj = tmc_ana%lc_traj +1
        IF(status.NE.0) THEN
          CALL cp_assert(.FALSE.,&
                 cp_failure_level,cp_assertion_failed,routineP,&
                 "configuration pos read error at line: "//&
                 cp_to_string(tmc_ana%lc_traj), error)
        END IF
      END DO pos_loop
      elem%pos(:) = elem%pos(:)/au2a
    END IF

    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE read_pos_from_file

! *****************************************************************************
!> \brief search for the dipole entry
!> \param elem actual tree element to be read
!> \param tmc_ana ...
!> \param stat ...
!> \param conf_nr Markov chain element number 
!>        (input the old number and read only if conf nr from file is greater
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 03.2013
! *****************************************************************************
  SUBROUTINE read_dipole_from_file(elem, tmc_ana, stat, conf_nr, error)
    TYPE(tree_type), POINTER                 :: elem
    TYPE(tmc_analysis_env), POINTER          :: tmc_ana
    INTEGER                                  :: stat, conf_nr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_dipole_from_file', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=250)                       :: c_tmp
    INTEGER                                  :: handle, status
    LOGICAL                                  :: failure

    failure = .FALSE.
    stat = TMC_STATUS_FAILED

    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(elem%dipole),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_ana),cp_failure_level,routineP,error,failure)
    CPPrecondition(tmc_ana%id_dip.GT.0,cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)
    READ(tmc_ana%id_dip,FMT="(A)", IOSTAT=status) c_tmp
    tmc_ana%lc_dip = tmc_ana%lc_dip +1
    IF(status.EQ.0) THEN
      ! skip the initial line (header)
      IF(INDEX(c_tmp,"#").GT.0) THEN
        READ(tmc_ana%id_dip,FMT="(A)", IOSTAT=status) c_tmp
        tmc_ana%lc_dip = tmc_ana%lc_dip +1
      END IF
    END IF
    IF(status.EQ.0) THEN
      READ(c_tmp,FMT="(I8,10F20.10)", IOSTAT=status) &
        conf_nr, elem%dipole(:)
    END IF
    IF(status.EQ.0) THEN ! success
      stat = TMC_STATUS_OK
    ELSE IF(status.LT.0) THEN ! end of file reached
      stat = TMC_STATUS_WAIT_FOR_NEW_TASK
    ELSE
      CALL cp_assert(status.EQ.0,&
             cp_warning_level,cp_assertion_failed,routineP,&
             "configuration dipole read error at line: "//cp_to_string(tmc_ana%lc_dip), error)
      stat = TMC_STATUS_FAILED
    END IF

    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE read_dipole_from_file

! *****************************************************************************
!> \brief search for the cell entry
!> \param elem actual tree element to be read
!> \param tmc_ana ...
!> \param stat ...
!> \param conf_nr Markov chain element number 
!>        (input the old number and read only if conf nr from file is greater
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 03.2013
! *****************************************************************************
  SUBROUTINE read_cell_from_file(elem, tmc_ana, stat, conf_nr, error)
    TYPE(tree_type), POINTER                 :: elem
    TYPE(tmc_analysis_env), POINTER          :: tmc_ana
    INTEGER                                  :: stat, conf_nr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_cell_from_file', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=250)                       :: c_tmp
    INTEGER                                  :: handle, status
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: r_tmp
    REAL(KIND=dp), DIMENSION(3, 3)           :: hmat

    failure = .FALSE.
    stat = TMC_STATUS_FAILED

    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_ana),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_ana%cell),cp_failure_level,routineP,error,failure)
    CPPrecondition(tmc_ana%id_cell.GT.0,cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    READ(tmc_ana%id_cell,FMT="(A)", IOSTAT=status) c_tmp
    tmc_ana%lc_cell = tmc_ana%lc_cell +1
    
    IF(status.EQ.0) THEN
      ! skip the initial line (header)
      IF(INDEX(c_tmp,"#").GT.0) THEN
        READ(tmc_ana%id_cell,FMT="(A)", IOSTAT=status) c_tmp
        tmc_ana%lc_cell = tmc_ana%lc_cell +1
      END IF
    END IF
    IF(status.EQ.0) THEN
        READ(c_tmp,FMT="(I8,9(1X,F19.10),1X,F24.10)", IOSTAT=status) conf_nr, &
             hmat(:,:), r_tmp
    END IF
    IF(status.LT.0) THEN ! end of file reached
      stat = TMC_STATUS_WAIT_FOR_NEW_TASK
    ELSE IF(status.GT.0) THEN
      CALL cp_assert(status.EQ.0,&
           cp_failure_level,cp_assertion_failed,routineP,&
           "configuration cell read error at line: "//cp_to_string(tmc_ana%lc_cell), error)
      stat = TMC_STATUS_FAILED
    ELSE
      IF(elem%nr.LT.0) elem%nr=conf_nr
      hmat(:,:) = hmat(:,:)/au2a
      ! get the box scaling
      CALL get_cell_scaling(cell=tmc_ana%cell, scaled_hmat=hmat, &
                            box_scale=elem%box_scale, error=error)
      stat = TMC_STATUS_OK
    END IF
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE read_cell_from_file

  !----------------------------------------------------------------------------
  ! get the configurations from file and calc 
  !----------------------------------------------------------------------------

! *****************************************************************************
!> \brief opens the files for reading configurations data to analyze
!> \param tmc_ana ...
!> \param stat ...
!> \param dir_ind ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param 
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE analyse_files_open(tmc_ana, stat, dir_ind, error)
    TYPE(tmc_analysis_env), POINTER          :: tmc_ana
    INTEGER                                  :: stat
    INTEGER, OPTIONAL                        :: dir_ind
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'analyse_files_open', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: dir_name, file_name_act, &
                                                file_name_temp
    INTEGER                                  :: handle
    LOGICAL                                  :: failure, file_exists

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(tmc_ana),cp_failure_level,routineP,error,failure)

    stat = TMC_STATUS_WAIT_FOR_NEW_TASK

    ! start the timing
    CALL timeset(routineN,handle)

    IF(PRESENT(dir_ind)) THEN
      CPPrecondition(ASSOCIATED(tmc_ana%dirs),cp_failure_level,routineP,error,failure)
      CPPrecondition(dir_ind.GT.0,cp_failure_level,routineP,error,failure)
      CPPrecondition(dir_ind.LE.SIZE(tmc_ana%dirs),cp_failure_level,routineP,error,failure)

      IF(INDEX(tmc_ana%dirs(dir_ind),"/", BACK=.TRUE.).EQ.&
         LEN_TRIM(tmc_ana%dirs(dir_ind))) THEN
        dir_name = TRIM(tmc_ana%dirs(dir_ind))
      ELSE
        dir_name = TRIM(tmc_ana%dirs(dir_ind))//"/"
      END IF
    ELSE
      dir_name = "./"
    END IF

    ! open the files 
    file_name_temp = expand_file_name_temp(&
                       file_name=tmc_default_trajectory_file_name, &
                       rvalue=tmc_ana%temperature, error=error)
    ! position file
    IF(tmc_ana%costum_pos_file_name.NE."") THEN
      file_name_act = TRIM(dir_name)//tmc_ana%costum_pos_file_name
    ELSE
      file_name_act = TRIM(dir_name)//&
                      expand_file_name_ending(file_name_temp,"xyz",error)
    END IF
    INQUIRE(FILE=file_name_act, EXIST=file_exists)
    IF(file_exists) THEN
      CALL open_file(file_name=file_name_act, file_status="OLD", &
                   file_action="READ", unit_number=tmc_ana%id_traj)
      WRITE(tmc_ana%io_unit,FMT='(T2,A,"| ",A,T41,A40)')"TMC_ANA",&
            "read xyz file",TRIM(file_name_act)
    END IF

    ! cell file
    IF(tmc_ana%costum_cell_file_name.NE."") THEN
      file_name_act = TRIM(dir_name)//tmc_ana%costum_cell_file_name
    ELSE 
      file_name_act = TRIM(dir_name)//&
                      expand_file_name_ending(file_name_temp,"cell",error)
    END IF
    INQUIRE(FILE=file_name_act, EXIST=file_exists)
    IF(file_exists) THEN
      CALL open_file(file_name=file_name_act, file_status="OLD", &
                     file_action="READ", unit_number=tmc_ana%id_cell)
      WRITE(tmc_ana%io_unit,FMT='(T2,A,"| ",A,T41,A40)')"TMC_ANA",&
            "read cell file",TRIM(file_name_act)
    END IF

    ! dipole file
    IF(tmc_ana%costum_dip_file_name.NE."") THEN
      file_name_act = TRIM(dir_name)//tmc_ana%costum_dip_file_name
    ELSE 
      file_name_act = TRIM(dir_name)//&
                      expand_file_name_ending(file_name_temp,"dip",error)
    END IF
    INQUIRE(FILE=file_name_act, EXIST=file_exists)
    IF(file_exists) THEN
      CALL open_file(file_name=file_name_act, file_status="OLD", &
                   file_action="READ", unit_number=tmc_ana%id_dip)
      WRITE(tmc_ana%io_unit,FMT='(T2,A,"| ",A,T41,A40)')"TMC_ANA",&
            "read dip file",TRIM(file_name_act)
    END IF

    IF(tmc_ana%id_traj.GT.0 .OR. tmc_ana%id_cell.GT.0 .OR. &
                   tmc_ana%id_dip.GT.0) THEN
      stat = TMC_STATUS_OK
    ELSE
      CALL cp_assert(.FALSE.,&
           cp_warning_level,cp_assertion_failed,routineP,&
           "There is no file to open for temperature "//cp_to_string(tmc_ana%temperature)//&
           "K in directory "//TRIM(dir_name),&
           error)
    END IF
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE analyse_files_open

! *****************************************************************************
!> \brief close the files for reading configurations data to analyze
!> \param tmc_ana ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param 
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE analyse_files_close(tmc_ana, error)
    TYPE(tmc_analysis_env), POINTER          :: tmc_ana
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'analyse_files_close', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(tmc_ana),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    ! position file
    IF(tmc_ana%id_traj.GT.0)  CALL close_file(unit_number=tmc_ana%id_traj) 

    ! cell file
    IF(tmc_ana%id_cell.GT.0) CALL close_file(unit_number=tmc_ana%id_cell)

    ! dipole file
    IF(tmc_ana%id_dip.GT.0) CALL close_file(unit_number=tmc_ana%id_dip)

    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE analyse_files_close

END MODULE tmc_file_io
