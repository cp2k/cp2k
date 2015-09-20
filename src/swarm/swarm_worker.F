!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Workers's routines for the swarm-framework
!> \author Ole Schuett
! *****************************************************************************
MODULE swarm_worker
  USE cp_output_handling,              ONLY: cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE glbopt_worker,                   ONLY: glbopt_worker_execute,&
                                             glbopt_worker_finalize,&
                                             glbopt_worker_init,&
                                             glbopt_worker_type
  USE input_constants,                 ONLY: swarm_do_glbopt
  USE input_section_types,             ONLY: section_type,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length
  USE swarm_message,                   ONLY: swarm_message_add,&
                                             swarm_message_get,&
                                             swarm_message_haskey,&
                                             swarm_message_type
#include "../common/cp_common_uses.f90"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'swarm_worker'

 PUBLIC :: swarm_worker_init, swarm_worker_finalize
 PUBLIC :: swarm_worker_execute
 PUBLIC :: swarm_worker_type

  TYPE swarm_worker_type
   PRIVATE
   INTEGER                                  :: id = -1
   INTEGER                                  :: iw = -1
   INTEGER                                  :: behavior = -1
   TYPE(glbopt_worker_type), POINTER        :: glbopt => Null()
   !possbily more behaviors...
 END TYPE swarm_worker_type


 CONTAINS


! *****************************************************************************
!> \brief Initializes a swarm worker
!> \param worker ...
!> \param para_env ...
!> \param input_declaration ...
!> \param root_section ...
!> \param input_path ...
!> \param worker_id ...
!> \author Ole Schuett
! *****************************************************************************
   SUBROUTINE swarm_worker_init(worker, para_env, input_declaration, root_section,&
                 input_path, worker_id)
    TYPE(swarm_worker_type), INTENT(INOUT)   :: worker
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_type), POINTER              :: input_declaration
    TYPE(section_vals_type), POINTER         :: root_section
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    INTEGER, INTENT(in)                      :: worker_id

    TYPE(cp_logger_type), POINTER            :: logger

    worker%id = worker_id

    ! getting an output unit for logging
    logger => cp_get_default_logger()
    worker%iw = cp_print_key_unit_nr(logger,root_section,&
          "SWARM%PRINT%WORKER_RUN_INFO",extension=".workerLog")

    CALL section_vals_val_get(root_section,"SWARM%BEHAVIOR", i_val=worker%behavior)

    SELECT CASE(worker%behavior)
      CASE(swarm_do_glbopt)
         ALLOCATE(worker%glbopt)
         CALL glbopt_worker_init(worker%glbopt, input_declaration, para_env,&
                      root_section, input_path, worker_id, worker%iw)
      CASE DEFAULT
         STOP "swarm_worker_init: got unknown behavior"
    END SELECT

  END SUBROUTINE swarm_worker_init


! *****************************************************************************
!> \brief Central execute routine of the swarm worker
!> \param worker ...
!> \param cmd ...
!> \param report ...
!> \param should_stop ...
!> \author Ole Schuett
! *****************************************************************************
   SUBROUTINE swarm_worker_execute(worker, cmd, report, should_stop)
    TYPE(swarm_worker_type), INTENT(INOUT)   :: worker
    TYPE(swarm_message_type), INTENT(IN)     :: cmd
    TYPE(swarm_message_type), INTENT(OUT)    :: report
    LOGICAL, INTENT(INOUT)                   :: should_stop

    CHARACTER(LEN=default_string_length)     :: command

     CALL swarm_message_get(cmd, "command", command)
     CALL swarm_message_add(report, "worker_id", worker%id)

     IF(TRIM(command) == "shutdown") THEN
        IF(worker%iw>0) WRITE(worker%iw,*) "SWARM| Received shutdown command, quitting."
        should_stop = .TRUE.
     ELSE IF(TRIM(command) == "wait") THEN  !only needed for serial driver
        CALL swarm_message_add(report, "status", "wait_done")
     ELSE
        SELECT CASE(worker%behavior)
           CASE(swarm_do_glbopt)
             CALL glbopt_worker_execute(worker%glbopt, cmd, report)
           CASE DEFAULT
              STOP "swarm_worker_execute: got unknown behavior"
        END SELECT
     ENDIF

     IF(.NOT. swarm_message_haskey(report, "status")) &
        CALL swarm_message_add(report, "status", "ok")

   END SUBROUTINE swarm_worker_execute


! *****************************************************************************
!> \brief Finalizes a swarm worker
!> \param worker ...
!> \author Ole Schuett
! *****************************************************************************
   SUBROUTINE swarm_worker_finalize(worker)
    TYPE(swarm_worker_type), INTENT(INOUT)   :: worker

     SELECT CASE(worker%behavior)
      CASE(swarm_do_glbopt)
         CALL glbopt_worker_finalize(worker%glbopt)
         DEALLOCATE(worker%glbopt)
      CASE DEFAULT
         STOP "swarm_worker_finalize: got unknown behavior"
    END SELECT

   END SUBROUTINE swarm_worker_finalize


END MODULE swarm_worker

