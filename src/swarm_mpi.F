!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Handles the MPI communication of the swarm framework.
!> \author Ole Schuett
! *****************************************************************************
MODULE swarm_mpi
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_iter_types,                   ONLY: cp_iteration_info_create,&
                                             cp_iteration_info_release,&
                                             cp_iteration_info_type
  USE cp_para_env,                     ONLY: cp_para_env_create,&
                                             cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length
  USE machine,                         ONLY: default_output_unit
  USE message_passing,                 ONLY: mp_any_source,&
                                             mp_bcast,&
                                             mp_comm_free,&
                                             mp_comm_split,&
                                             mp_comm_split_direct,&
                                             mp_environ,&
                                             mp_sum,&
                                             mp_sync
  USE swarm_message,                   ONLY: swarm_message_get,&
                                             swarm_message_mpi_bcast,&
                                             swarm_message_mpi_recv,&
                                             swarm_message_mpi_send,&
                                             swarm_message_type
#include "./common/cp_common_uses.f90"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'swarm_mpi'

 PUBLIC :: swarm_mpi_type, swarm_mpi_init, swarm_mpi_finalize
 PUBLIC :: swarm_mpi_send_report, swarm_mpi_recv_report
 PUBLIC :: swarm_mpi_send_command, swarm_mpi_recv_command


 TYPE swarm_mpi_type
    TYPE(cp_para_env_type), POINTER          :: world => Null()
    TYPE(cp_para_env_type), POINTER          :: worker  => Null()
    TYPE(cp_para_env_type), POINTER          :: master  => Null()
    INTEGER, DIMENSION(:), ALLOCATABLE       :: wid2group
    CHARACTER(LEN=default_path_length)       :: master_output_path = ""
 END TYPE swarm_mpi_type

 CONTAINS


! *****************************************************************************
!> \brief Initialize MPI communicators for a swarm run.
!> \param swarm_mpi ...
!> \param world_para_env ...
!> \param root_section ...
!> \param n_workers ...
!> \param worker_id ...
!> \param iw ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
   SUBROUTINE swarm_mpi_init(swarm_mpi, world_para_env, root_section, n_workers, worker_id, iw, error)
    TYPE(swarm_mpi_type)                     :: swarm_mpi
    TYPE(cp_para_env_type), POINTER          :: world_para_env
    TYPE(section_vals_type), POINTER         :: root_section
    INTEGER, INTENT(IN)                      :: n_workers
    INTEGER, INTENT(OUT)                     :: worker_id
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: n_groups_created, &
                                                pe_per_worker, subgroup, &
                                                subgroup_rank, subgroup_size, &
                                                worker_group
    LOGICAL                                  :: im_the_master
    INTEGER, DIMENSION(:), POINTER           :: group_distribution_p
    INTEGER, &
      DIMENSION(0:world_para_env%num_pe-2), &
      TARGET                                 :: group_distribution

! ====== Setup of MPI-Groups ======

    worker_id = -1
    swarm_mpi%world => world_para_env

    IF (MOD(swarm_mpi%world%num_pe-1, n_workers) /= 0) &
       STOP "number of processors-1 is not divisible by n_workers."
    IF (swarm_mpi%world%num_pe < n_workers + 1) &
       STOP "There are not enough processes for n_workers + 1. Aborting."

    pe_per_worker = (swarm_mpi%world%num_pe-1)/n_workers

    IF(iw>0) THEN
       WRITE(iw,'(A,45X,I8)') " SWARM| Number of mpi ranks", swarm_mpi%world%num_pe
       WRITE(iw,'(A,47X,I8)') " SWARM| Number of workers", n_workers
    ENDIF

    ! the last task becomes the master. Preseves node-alignment of other tasks.
    im_the_master = (swarm_mpi%world%mepos == swarm_mpi%world%num_pe-1)

    ! First split split para_env%group into a master- and a workers-groups...
    IF (im_the_master) THEN
       CALL mp_comm_split_direct(swarm_mpi%world%group, subgroup, 1)
       CALL mp_environ(subgroup_size, subgroup_rank, subgroup)
       IF(subgroup_size/=1) STOP "swarm: mp_comm_split_direct failed (master)"
       CALL cp_para_env_create(swarm_mpi%master, group=subgroup, error=error)
       !WRITE (*,*) "this is a master ", swarm_mpi%master%mepos, swarm_mpi%master%num_pe
    ELSE
       CALL mp_comm_split_direct(swarm_mpi%world%group, subgroup, 2)
       CALL mp_environ(subgroup_size, subgroup_rank, subgroup)
       !WRITE (*,*) "Hello, this is a Worker - there are ",subgroup_size, " of us."
       IF(subgroup_size/=swarm_mpi%world%num_pe-1) STOP "swarm: mp_comm_split_direct failed (worker)"
    ENDIF

    ALLOCATE(swarm_mpi%wid2group(n_workers))
    swarm_mpi%wid2group = 0

    IF(.NOT. im_the_master) THEN
       ! ...then split workers-group into n_workers groups - one for each worker.
       group_distribution_p => group_distribution
       CALL mp_comm_split(subgroup, worker_group, n_groups_created, group_distribution_p, n_subgroups=n_workers)
       worker_id = group_distribution(subgroup_rank) + 1 ! shall start by 1
       IF(n_groups_created/=n_workers) STOP "swarm: mp_comm_split failed."
       CALL cp_para_env_create(swarm_mpi%worker, group=worker_group, error=error)
       CALL mp_comm_free(subgroup)

       !WRITE (*,*) "this is worker ", worker_id, swarm_mpi%worker%mepos, swarm_mpi%worker%num_pe

       ! collect world-ranks of each worker groups rank-0 node
       IF(swarm_mpi%worker%mepos == 0) &
          swarm_mpi%wid2group(worker_id) = swarm_mpi%world%mepos

    ENDIF

    CALL mp_sum(swarm_mpi%wid2group, swarm_mpi%world%group)
    !WRITE (*,*), "wid2group table: ",swarm_mpi%wid2group


    CALL logger_init_master(swarm_mpi, error)
    CALL logger_init_worker(swarm_mpi, root_section, worker_id, error)
  END SUBROUTINE swarm_mpi_init


! *****************************************************************************
!> \brief Helper routine for swarm_mpi_init, configures the master's logger.
!> \param swarm_mpi ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE logger_init_master(swarm_mpi, error)
    TYPE(swarm_mpi_type)                     :: swarm_mpi
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: output_unit
    TYPE(cp_logger_type), POINTER            :: logger

! broadcast master_output_path to all ranks

    IF(swarm_mpi%world%mepos == swarm_mpi%world%source) THEN
       logger => cp_error_get_logger(error)
       output_unit = logger%default_local_unit_nr
       swarm_mpi%master_output_path = output_unit2path(output_unit)
       IF(output_unit/=default_output_unit) &
         CLOSE(output_unit)
    END IF

    CALL mp_bcast(swarm_mpi%master_output_path, swarm_mpi%world%source, swarm_mpi%world%group)

    IF(ASSOCIATED(swarm_mpi%master)) &
       CALL error_add_new_logger(swarm_mpi%master, swarm_mpi%master_output_path, error=error)
  END SUBROUTINE logger_init_master


! *****************************************************************************
!> \brief Helper routine for logger_init_master, inquires filename for given unit.
!> \param output_unit ...
!> \retval output_path ...
!> \author Ole Schuett
! *****************************************************************************
  FUNCTION output_unit2path(output_unit) RESULT(output_path)
    INTEGER, INTENT(IN)                      :: output_unit
    CHARACTER(LEN=default_path_length)       :: output_path

     output_path = "__STD_OUT__"
     IF(output_unit /= default_output_unit)&
         INQUIRE(unit=output_unit, name=output_path)
  END FUNCTION output_unit2path


! *****************************************************************************
!> \brief Helper routine for swarm_mpi_init, configures the workers's logger.
!> \param swarm_mpi ...
!> \param root_section ...
!> \param worker_id ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE logger_init_worker(swarm_mpi, root_section, worker_id, error)
    TYPE(swarm_mpi_type)                     :: swarm_mpi
    TYPE(section_vals_type), POINTER         :: root_section
    INTEGER                                  :: worker_id
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=default_path_length)       :: output_path
    CHARACTER(len=default_string_length)     :: new_project_name, &
                                                project_name, worker_name
    TYPE(cp_iteration_info_type), POINTER    :: new_iter_info
    TYPE(cp_logger_type), POINTER            :: old_logger

    NULLIFY(old_logger, new_iter_info)
    IF(ASSOCIATED(swarm_mpi%worker)) THEN
       old_logger => cp_error_get_logger(error)
       project_name = old_logger%iter_info%project_name
       IF(worker_id > 99999) &
          STOP "logger_init_worker: Did not expect so many workers."
       WRITE(worker_name,"(A,I5.5)") 'WORKER', worker_id
       IF(LEN_TRIM(project_name)+1+LEN_TRIM(worker_name)>default_string_length) &
          STOP "logger_init_worker: project name too long"
       output_path = TRIM(project_name) // "-"// TRIM(worker_name) // ".out"
       new_project_name = TRIM(project_name) // "-"// TRIM(worker_name)
       CALL section_vals_val_set(root_section, "GLOBAL%PROJECT_NAME", c_val=new_project_name,error=error)
       CALL cp_iteration_info_create(new_iter_info, new_project_name)
       CALL error_add_new_logger(swarm_mpi%worker, output_path, new_iter_info, error)
       CALL cp_iteration_info_release(new_iter_info)
    ENDIF
  END SUBROUTINE logger_init_worker


! *****************************************************************************
!> \brief Helper routine for logger_init_master and logger_init_worker
!> \param para_env ...
!> \param output_path ...
!> \param iter_info ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE error_add_new_logger(para_env, output_path, iter_info, error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    CHARACTER(LEN=default_path_length)       :: output_path
    TYPE(cp_iteration_info_type), OPTIONAL, &
      POINTER                                :: iter_info
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: output_unit
    TYPE(cp_logger_type), POINTER            :: new_logger, old_logger

    NULLIFY(new_logger, old_logger)
    output_unit = -1
    IF(para_env%source == para_env%mepos) THEN
      ! open output_unit according to output_path
      output_unit = default_output_unit
      IF (output_path /= "__STD_OUT__") &
        CALL open_file(file_name=output_path,file_status="UNKNOWN",&
           file_action="WRITE", file_position="APPEND", unit_number=output_unit)
    ENDIF

    old_logger => cp_error_get_logger(error)
    CALL cp_logger_create(new_logger, para_env=para_env,&
       default_global_unit_nr=output_unit, close_global_unit_on_dealloc=.FALSE.,&
       template_logger=old_logger, iter_info=iter_info)

    CALL cp_add_default_logger(new_logger)
    CALL cp_error_set(error, logger=new_logger)
    CALL cp_logger_release(new_logger)
  END SUBROUTINE error_add_new_logger


! *****************************************************************************
!> \brief Finalizes the MPI communicators of a swarm run.
!> \param swarm_mpi ...
!> \param root_section ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_mpi_finalize(swarm_mpi, root_section, error)
    TYPE(swarm_mpi_type)                     :: swarm_mpi
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CALL mp_sync(swarm_mpi%world%group)
    CALL logger_finalize(swarm_mpi, root_section, error)

    IF(ASSOCIATED(swarm_mpi%worker)) CALL cp_para_env_release(swarm_mpi%worker, error)
    IF(ASSOCIATED(swarm_mpi%master)) CALL cp_para_env_release(swarm_mpi%master, error)
    NULLIFY(swarm_mpi%worker, swarm_mpi%master)
    DEALLOCATE(swarm_mpi%wid2group)
  END SUBROUTINE swarm_mpi_finalize


! *****************************************************************************
!> \brief Helper routine for swarm_mpi_finalize, restores the original loggers
!> \param swarm_mpi ...
!> \param root_section ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE logger_finalize(swarm_mpi, root_section, error)
    TYPE(swarm_mpi_type)                     :: swarm_mpi
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: output_unit
    LOGICAL                                  :: I_am_rank0
    TYPE(cp_logger_type), POINTER            :: logger, old_logger

    NULLIFY(logger, old_logger)
    logger => cp_error_get_logger(error)
    output_unit = logger%default_local_unit_nr
    IF(output_unit>0 .AND. output_unit/=default_output_unit)&
        CALL close_file(output_unit)

    CALL cp_rm_default_logger() !pops the top-most logger
    old_logger => cp_error_get_logger(error)
    CALL cp_error_set(error, logger=old_logger)

    ! restore GLOBAL%PROJECT_NAME
    CALL section_vals_val_set(root_section, "GLOBAL%PROJECT_NAME", &
     c_val=old_logger%iter_info%project_name, error=error)

    CALL mp_sync(swarm_mpi%world%group)

    ! do this only on master's rank 0
    I_am_rank0 = swarm_mpi%world%mepos==swarm_mpi%world%source
    IF(I_am_rank0 .AND. output_unit/=default_output_unit) THEN
       output_unit = old_logger%default_local_unit_nr
       OPEN(unit=output_unit, file=swarm_mpi%master_output_path, &
            status="UNKNOWN", action="WRITE", position="APPEND")
    ENDIF
  END SUBROUTINE logger_finalize



! *****************************************************************************
!> \brief Sends a report via MPI
!> \param swarm_mpi ...
!> \param report ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_mpi_send_report(swarm_mpi, report)
    TYPE(swarm_mpi_type)                     :: swarm_mpi
    TYPE(swarm_message_type)                 :: report

    INTEGER                                  :: dest, tag

! Only rank-0 of worker group sends it's report

     IF(swarm_mpi%worker%source /= swarm_mpi%worker%mepos) RETURN

     dest = swarm_mpi%world%num_pe-1
     tag = 42
     CALL swarm_message_mpi_send(report, group=swarm_mpi%world%group, dest=dest, tag=tag)

  END SUBROUTINE swarm_mpi_send_report


! *****************************************************************************
!> \brief Receives a report via MPI
!> \param swarm_mpi ...
!> \param report ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_mpi_recv_report(swarm_mpi, report)
    TYPE(swarm_mpi_type)                     :: swarm_mpi
    TYPE(swarm_message_type), INTENT(OUT)    :: report

    INTEGER                                  :: src, tag

     tag=42
     src=mp_any_source

     CALL swarm_message_mpi_recv(report, group=swarm_mpi%world%group, src=src, tag=tag)

  END SUBROUTINE swarm_mpi_recv_report


! *****************************************************************************
!> \brief Sends a command via MPI
!> \param swarm_mpi ...
!> \param cmd ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_mpi_send_command(swarm_mpi, cmd)
    TYPE(swarm_mpi_type)                     :: swarm_mpi
    TYPE(swarm_message_type)                 :: cmd

    INTEGER                                  :: dest, tag, worker_id

     CALL swarm_message_get(cmd, "worker_id", worker_id)
     tag = 42
     dest = swarm_mpi%wid2group(worker_id)

     CALL swarm_message_mpi_send(cmd, group=swarm_mpi%world%group, dest=dest, tag=tag)

  END SUBROUTINE swarm_mpi_send_command


! *****************************************************************************
!> \brief Receives a command via MPI and broadcasts it within a worker.
!> \param swarm_mpi ...
!> \param cmd ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_mpi_recv_command(swarm_mpi, cmd)
    TYPE(swarm_mpi_type)                     :: swarm_mpi
    TYPE(swarm_message_type), INTENT(OUT)    :: cmd

    INTEGER                                  :: src, tag

! This is a two step communication schema.
! First: The rank-0 of the worker groups receives the command from the master.

     IF(swarm_mpi%worker%source == swarm_mpi%worker%mepos) THEN
        src = swarm_mpi%world%num_pe-1!
        tag = 42
        CALL swarm_message_mpi_recv(cmd, group=swarm_mpi%world%group, src=src, tag=tag)

     ENDIF

!     ! Second: The command is broadcasted within the worker group.
      CALL swarm_message_mpi_bcast(cmd, src=swarm_mpi%worker%source, group=swarm_mpi%worker%group)

  END SUBROUTINE swarm_mpi_recv_command


END MODULE swarm_mpi

