!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!


! *****************************************************************************
!> \brief Master's routines for the swarm-framework
!> \author Ole Schuett
! *****************************************************************************
MODULE swarm_master
  USE cp_external_control,             ONLY: external_control
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE glbopt_master,                   ONLY: glbopt_master_finalize,&
                                             glbopt_master_init,&
                                             glbopt_master_steer,&
                                             glbopt_master_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: swarm_do_glbopt
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length
  USE swarm_message,                   ONLY: swarm_message_add,&
                                             swarm_message_equal,&
                                             swarm_message_file_read,&
                                             swarm_message_file_write,&
                                             swarm_message_free,&
                                             swarm_message_get,&
                                             swarm_message_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'swarm_master'


 PUBLIC :: swarm_master_type
 PUBLIC :: swarm_master_init, swarm_master_finalize
 PUBLIC :: swarm_master_steer

 TYPE swarm_message_p_type
   TYPE(swarm_message_type), POINTER                   :: p => Null()
 END TYPE swarm_message_p_type


 TYPE swarm_master_type
   PRIVATE
   INTEGER                                             :: behavior = -1
   TYPE(glbopt_master_type), POINTER                   :: glbopt => Null()
   !possibly more behaviors ...
   INTEGER                                             :: iw = 0
   INTEGER                                             :: i_iteration = 0
   INTEGER                                             :: max_iter = 0
   LOGICAL                                             :: should_stop = .FALSE.
   INTEGER                                             :: n_workers = -1
   INTEGER                                             :: comlog_unit
   TYPE(cp_error_type)                                 :: error
   TYPE(section_vals_type), POINTER                    :: swarm_section => Null()
   TYPE(cp_para_env_type), POINTER                     :: para_env => Null()
   TYPE(swarm_message_p_type), DIMENSION(:), POINTER   :: queued_commands => Null()
   TYPE(global_environment_type), POINTER              :: globenv => Null()
   LOGICAL                                             :: ignore_last_iteration = .FALSE.
   INTEGER                                             :: n_waiting = 0
 END TYPE swarm_master_type


 CONTAINS


! *****************************************************************************
!> \brief Initializes the swarm master
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_master_init(master, para_env, globenv, root_section, input_path, n_workers, error)
    TYPE(swarm_master_type)                  :: master
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: root_section
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    INTEGER, INTENT(IN)                      :: n_workers
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(cp_logger_type), POINTER            :: logger

    master%swarm_section => section_vals_get_subs_vals(root_section, "SWARM", error=error)

    logger => cp_error_get_logger(error)
    master%n_workers = n_workers
    master%error = error
    master%para_env => para_env
    master%globenv => globenv
    ALLOCATE(master%queued_commands(master%n_workers))
    master%iw = cp_print_key_unit_nr(logger, master%swarm_section,&
          "PRINT%MASTER_RUN_INFO",extension=".masterLog",error=error)

    CALL section_vals_val_get(master%swarm_section,"BEHAVIOR", i_val=master%behavior, error=error)

    ! uses logger%iter_info%project_name to construct filename
    master%comlog_unit = cp_print_key_unit_nr(logger,master%swarm_section,"PRINT%COMMUNICATION_LOG",&
            !middle_name="comlog", extension=".xyz", &
            extension=".comlog", &
            file_action="WRITE", file_position="REWIND", error=error)

    CALL section_vals_val_get(master%swarm_section,"MAX_ITER", i_val=master%max_iter, error=error)

    SELECT CASE(master%behavior)
      CASE(swarm_do_glbopt)
         ALLOCATE(master%glbopt)
         CALL glbopt_master_init(master%glbopt, para_env, root_section, input_path, n_workers, master%iw, error)
      CASE DEFAULT
         STOP "swarm_master_init: got unkown behavior"
    END SELECT

    CALL replay_comlog(master)
 END SUBROUTINE swarm_master_init


! *****************************************************************************
!> \brief Helper routine for swarm_master_init, restarts a calculation
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE replay_comlog(master)
    TYPE(swarm_master_type)                  :: master

    CHARACTER(LEN=default_path_length)       :: filename
    CHARACTER(LEN=default_string_length)     :: command_log
    INTEGER                                  :: handle, i, worker_id
    LOGICAL                                  :: at_end, explicit
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(swarm_message_type)                 :: cmd_log, report_log
    TYPE(swarm_message_type), &
      DIMENSION(master%n_workers)            :: last_commands
    TYPE(swarm_message_type), POINTER        :: cmd_now

    NULLIFY(parser)

    ! Initialize parser for trajectory
    CALL section_vals_val_get(master%swarm_section,"REPLAY_COMMUNICATION_LOG",&
       c_val=filename, explicit=explicit, error=master%error)

    IF(.NOT. explicit) RETURN
    IF(master%iw>0) WRITE(master%iw,'(A,A)') &
       " SWARM| Starting replay of communication-log: ", TRIM(filename)

    CALL timeset("swarm_master_replay_comlog", handle)
    CALL parser_create(parser, filename, para_env=master%para_env, error=master%error)

    at_end = .FALSE.
    DO
      CALL swarm_message_file_read(report_log, parser, at_end, master%error)
      IF(at_end) EXIT

      CALL swarm_message_file_read(cmd_log, parser, at_end, master%error)
      IF(at_end) EXIT

      ALLOCATE(cmd_now)
      CALL swarm_master_steer(master, report_log, cmd_now)

      !TODO: maybe we should just exit the loop instead of stopping?
      CALL swarm_message_get(cmd_log, "command", command_log)
      IF(TRIM(command_log) /= "shutdown") THEN
        IF(.NOT. commands_equal(cmd_now, cmd_log, master%iw)) STOP "wrong behaviour"
      END IF

      CALL swarm_message_free(cmd_log)
      CALL swarm_message_free(report_log)
      CALL swarm_message_get(cmd_now, "worker_id", worker_id)
      CALL swarm_message_free(last_commands(worker_id))
      last_commands(worker_id) = cmd_now
      DEALLOCATE(cmd_now)
    END DO

    CALL swarm_message_free(report_log) !don't worry about double-frees
    CALL swarm_message_free(cmd_log)

    IF(master%iw>0) WRITE(master%iw,'(A,A)') &
       " SWARM| Reached end of communication log. Queueing last commands."

    DO i=1, master%n_workers
       ALLOCATE(master%queued_commands(i)%p)
       master%queued_commands(i)%p = last_commands(i)
    END DO

    CALL parser_release(parser, error=master%error)
    CALL timestop(handle)
 END SUBROUTINE replay_comlog


! *****************************************************************************
!> \brief Helper routine for replay_comlog, compares two commands
!> \author Ole Schuett
! *****************************************************************************
 FUNCTION commands_equal(cmd1, cmd2, iw) RESULT(res)
    TYPE(swarm_message_type)                 :: cmd1, cmd2
    INTEGER                                  :: iw
    LOGICAL                                  :: res

    res = swarm_message_equal(cmd1, cmd2)
    IF(.NOT. res .AND. iw>0) THEN
       WRITE(iw, *) "Command 1:"
       CALL swarm_message_file_write(cmd1, iw)
       WRITE(iw, *) "Command 2:"
       CALL swarm_message_file_write(cmd2, iw)
    END IF
 END FUNCTION commands_equal


! *****************************************************************************
!> \brief Central steering routine of the swarm master
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_master_steer(master, report, cmd)
    TYPE(swarm_master_type), INTENT(INOUT)   :: master
    TYPE(swarm_message_type), INTENT(IN)     :: report
    TYPE(swarm_message_type), INTENT(OUT)    :: cmd

    CHARACTER(len=default_string_length)     :: command, status
    INTEGER                                  :: handle, worker_id
    LOGICAL                                  :: should_stop

    should_stop = .FALSE.

    CALL timeset("swarm_master_steer", handle)

    ! First check if there are queued commands for this worker
    CALL swarm_message_get(report, "worker_id", worker_id)

    IF(ASSOCIATED(master%queued_commands(worker_id)%p)) THEN
       cmd = master%queued_commands(worker_id)%p
       DEALLOCATE(master%queued_commands(worker_id)%p)
       IF(master%iw>0) WRITE (master%iw,'(A,A,A,I9,1X,A)') ' SWARM| ',&
          REPEAT("*", 9), " Sending out queued command to worker: ",&
          worker_id,REPEAT("*", 9)
       CALL timestop(handle)
       RETURN
    END IF


    IF(.NOT. master%ignore_last_iteration) THEN
       ! There are no queued commands. Do the normal processing.
       master%i_iteration = master%i_iteration + 1

       IF(master%iw>0) WRITE (master%iw,'(A,A,1X,I8,A,A)') ' SWARM| ',REPEAT("*", 15),&
           master%i_iteration,' Master / Worker Communication  ',REPEAT("*", 15)
    ENDIF

    IF(master%i_iteration>=master%max_iter .AND. .NOT. master%should_stop) THEN
       IF(master%iw>0)  WRITE (master%iw, '(A)')  " SWARM| Reached MAX_ITER. Quitting."
       master%should_stop = .TRUE.
    ENDIF

    IF(.NOT. master%should_stop) THEN
      CALL external_control(master%should_stop,"SWARM",master%globenv,error=master%error)
      IF(master%should_stop .AND. master%iw>0) &
         WRITE (master%iw, *)  " SWARM| Received stop from external_control. Quitting."
    END IF

    !IF(unit > 0) &

    IF(master%should_stop) THEN
       CALL swarm_message_add(cmd, "command", "shutdown")
       IF(master%iw>0) WRITE(master%iw,'(1X,A,T71,I10)') &
         "SWARM| Sending shutdown command to worker",worker_id
    ELSE
       SELECT CASE(master%behavior)
         CASE(swarm_do_glbopt)
            CALL glbopt_master_steer(master%glbopt, report, cmd, should_stop)
         CASE DEFAULT
            STOP "swarm_master_steer: got unkown behavior"
       END SELECT

       IF(should_stop) THEN
          CALL swarm_message_free(cmd)
          CALL swarm_message_add(cmd, "command", "shutdown") !overwrite command
          IF(master%iw>0) WRITE(master%iw,'(1X,A,T71,I10)') &
              "SWARM| Sending shutdown command to worker",worker_id
           master%should_stop = .TRUE.
       ENDIF
    ENDIF

    CALL swarm_message_add(cmd, "worker_id", worker_id)

    ! Don't pollute comlog with "continue waiting"-commands.
    CALL swarm_message_get(report, "status", status)
    CALL swarm_message_get(cmd, "command", command)
    IF(TRIM(status)=="wait_done") master%n_waiting = master%n_waiting -1
    IF(TRIM(command)=="wait")     master%n_waiting = master%n_waiting +1
    IF(master%n_waiting < 0) STOP "master%n_waiting < 0"
    IF(TRIM(status)/="wait_done" .OR. TRIM(command)/="wait") THEN
       CALL swarm_message_file_write(report, master%comlog_unit)
       CALL swarm_message_file_write(cmd, master%comlog_unit)
       IF(master%n_waiting>0 .AND. master%iw>0) WRITE(master%iw,'(1X,A,T71,I10)') &
         "SWARM| Number of waiting workers:", master%n_waiting
       master%ignore_last_iteration = .FALSE.
    ELSE
       master%ignore_last_iteration = .TRUE.
    END IF
    CALL timestop(handle)
 END SUBROUTINE swarm_master_steer


! *****************************************************************************
!> \brief Finalizes the swarm master
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_master_finalize(master)
    TYPE(swarm_master_type)                  :: master

    TYPE(cp_logger_type), POINTER            :: logger

    IF(master%iw > 0) THEN
       WRITE(master%iw,"(1X,A,T71,I10)") "SWARM| Total number of iterations ", master%i_iteration
       WRITE(master%iw,"(A)") " SWARM| Shutting down the master."
    ENDIF

    SELECT CASE(master%behavior)
      CASE(swarm_do_glbopt)
         CALL glbopt_master_finalize(master%glbopt)
         DEALLOCATE(master%glbopt)
      CASE DEFAULT
         STOP "swarm_master_finalize: got unkown behavior"
    END SELECT

    DEALLOCATE(master%queued_commands)

    logger => cp_error_get_logger(master%error)
    CALL cp_print_key_finished_output(master%iw, logger,&
            master%swarm_section, "PRINT%MASTER_RUN_INFO", error=master%error)
    CALL cp_print_key_finished_output(master%comlog_unit, logger,&
            master%swarm_section, "PRINT%COMMUNICATION_LOG", error=master%error)

    !CALL rm_timer_env() !pops the top-most timer
 END SUBROUTINE swarm_master_finalize


END MODULE swarm_master

