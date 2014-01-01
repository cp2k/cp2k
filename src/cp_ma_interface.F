!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   MACHINE ARCH interface in CP2K
!> \author  Christiane Pousa Ribeiro
!> \date    2011
!> \version 0.1
!>
!> <b>Modification history:</b>
!> - Created 2011
! *****************************************************************************
MODULE cp_ma_interface
  USE cp2k_info,                       ONLY: r_pid
  USE cp_error_handling
  USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit,&
                                             cp_logger_type,&
                                             cp_logger_would_log,&
                                             cp_note_level,&
                                             cp_warning_level
  USE cp_output_handling
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dbcsr_config,                    ONLY: dbcsr_get_conf_use_comm_thread
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE ma_affinity
  USE ma_config
  USE ma_process_mapping,              ONLY: ma_mpi_reordering
  USE ma_topology,                     ONLY: ma_destroy_topology
  USE machine_architecture,            ONLY: &
       ma_finalize_machine, ma_get_nnodes, ma_hwloc, ma_init_machine, &
       ma_int_hwloc, ma_int_libnuma, ma_int_none, ma_interface, ma_libnuma, &
       ma_show_machine_branch, ma_show_machine_full, ma_show_topology
  USE machine_architecture_types
  USE message_passing,                 ONLY: mp_environ
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

!$ USE OMP_LIB

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_ma_interface'

  PRIVATE

  ! Interface to libma
  PUBLIC :: cp_ma_config, cp_ma_init_lib, cp_ma_finalize_lib
  PUBLIC :: cp_ma_run_on, cp_ma_thread_run_on
  PUBLIC :: cp_ma_current_thread_run
  PUBLIC :: cp_ma_set_neighbors, cp_ma_get_neighbors
  PUBLIC :: cp_ma_thread_sched, cp_ma_verify_place
  PUBLIC :: cp_ma_mpi_sched
  PUBLIC :: cp_ma_mempol
  PUBLIC :: cp_ma_print_machine
  PUBLIC :: cp_ma_get_cuDev
  PUBLIC :: cp_ma_net_affinity
  PUBLIC :: cp_ma_mpi_reorder_strategy
  PUBLIC :: cp_ma_mpi_reordering
  PUBLIC :: cp_ma_print_strategy
  PUBLIC :: cp_ma_set_mpi_reordering
  PUBLIC :: cp_ma_default_affinity
  PUBLIC :: has_ma, has_ma_topology

CONTAINS

! *****************************************************************************
!> \brief Configures options for Machine Architecture
! *****************************************************************************
  SUBROUTINE cp_ma_config(root_section, error)
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_config', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=1)                         :: mpi_reorder, mpi_sched, &
                                                use_sched
    INTEGER, DIMENSION(:), POINTER           :: mem, proc, thr
    LOGICAL :: comm_thread, print_branch, print_full, print_proc, &
      print_resume, print_thread, print_thread_cur
    TYPE(section_vals_type), POINTER         :: ma_section

    NULLIFY(proc,mem,thr)

    ma_section => section_vals_get_subs_vals(root_section,&
      "GLOBAL%MACHINE_ARCH",error=error)
    CALL section_vals_val_get(ma_section,&
         "PRINT_BRANCH", l_val=print_branch, error=error)
    CALL section_vals_val_get(ma_section,&
         "PRINT_RESUME", l_val=print_resume, error=error)
    CALL section_vals_val_get(ma_section,&
         "PRINT_FULL", l_val=print_full, error=error)
    CALL section_vals_val_get(ma_section,&
         "PRINT_PROC", l_val=print_proc, error=error)
    CALL section_vals_val_get(ma_section,&
         "PRINT_THREAD", l_val=print_thread, error=error)
    CALL section_vals_val_get(ma_section,&
         "PRINT_THREAD_CUR", l_val=print_thread_cur, error=error)
    CALL section_vals_val_get(ma_section,&
         "SCHED_MPI", c_val=mpi_sched, error=error)
    CALL section_vals_val_get(ma_section,&
         "SCHED_THREAD", c_val=use_sched, error=error)
    CALL section_vals_val_get(ma_section,&
         "MEM_POL", c_val=use_mempol, error=error)
    CALL section_vals_val_get(ma_section,&
         "MPI_REORDERING", c_val=mpi_reorder, error=error)  
    CALL section_vals_val_get(ma_section,"PROCESS",&
            i_vals=proc, error=error)
    CALL section_vals_val_get(ma_section,"MEMORY",&
                    i_vals=mem, error=error)
    CALL section_vals_val_get(ma_section,"THREAD",&
                    i_vals=thr, error=error)

    CALL ma_set_conf_print_branch (print_branch)
    CALL ma_set_conf_print_resume (print_resume)
    CALL ma_set_conf_print_full (print_full)
    CALL ma_set_conf_print_proc (print_proc)
    CALL ma_set_conf_print_thread (print_thread)
    CALL ma_set_conf_print_thread_cur (print_thread_cur)
    CALL ma_set_conf_mempol (use_mempol)
    CALL ma_set_conf_sched (use_sched)
    CALL ma_set_conf_mpiSched (mpi_sched)
    CALL ma_set_conf_mpi_reordering(mpi_reorder)
    comm_thread = dbcsr_get_conf_use_comm_thread()
    CALL ma_set_conf_comm_thread(comm_thread)
    CALL ma_set_all_affinty(proc, thr, mem)

  END SUBROUTINE cp_ma_config

! *****************************************************************************
!> \brief Set the network card affinity for a MPI
!> \var proc is mpi for MPI process and threads for OpenMP threads
!> \var id is the MPI rank or thread openmp id
! NOTE: ONLY work on scenarios of one MPI per compute node
! *****************************************************************************
SUBROUTINE cp_ma_net_affinity(proc,id)
    INTEGER                                  :: proc, id

   CALL ma_set_net_affinity (proc,id)

END SUBROUTINE cp_ma_net_affinity

! *****************************************************************************
!> \brief Get the cuda devices for the MPIs
! *****************************************************************************
  SUBROUTINE cp_ma_get_cuDev(cuDev, myRank)
    INTEGER                                  :: cuDev, myRank

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_get_cuDev', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_num
    TYPE(cp_error_type)                      :: error
    TYPE(cp_logger_type), POINTER            :: logger

    logger => cp_error_get_logger(error)

    IF ( has_ma ) THEN
            unit_num = cp_logger_get_default_io_unit(logger)
            cuDev = ma_set_gpu_affinity(myRank)

            IF (unit_num > 0 .AND. cp_logger_would_log(logger, cp_note_level)) THEN
              IF (ma_MPI_nGpus()) THEN
                WRITE(unit_num,"(T2,A)")"WARNING: Number of MPI greater than GPUs"
              ENDIF
            ENDIF
    ELSE
        CALL cp_assert (.FALSE.,&
        cp_warning_level, cp_unimplemented_error_nr, routineP,&
        "Need hwloc/libnuma library to query machine architecture.")
    ENDIF

  END SUBROUTINE  cp_ma_get_cuDev

! *****************************************************************************
!> \brief Print the configuration options for Machine Architecture
! *****************************************************************************
  SUBROUTINE cp_ma_print_config(unit_nr, error)
    INTEGER, INTENT(IN), OPTIONAL            :: unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_print_config', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=1)                         :: sched_thread
    INTEGER                                  :: pol_sched, unit_num
    TYPE(cp_logger_type), POINTER            :: logger

    sched_thread = ma_get_conf_sched()

    logger => cp_error_get_logger(error)

    IF (has_ma) THEN
      IF (sched_thread .EQ. 'L' .OR. sched_thread .EQ. 'l') THEN
           pol_sched = linear
      ELSE IF (sched_thread .EQ. 'S' .OR. sched_thread .EQ. 's') THEN
           pol_sched = scatter
      ELSE IF (sched_thread .EQ. 'G' .OR. sched_thread .EQ. 'g') THEN
           pol_sched = group
      ELSE IF (sched_thread .EQ. 'M' .OR. sched_thread .EQ. 'm') THEN
           pol_sched = manual
      ELSE
           pol_sched = nosched
      ENDIF

    IF(PRESENT (unit_nr)) THEN
      unit_num = unit_nr
    ELSE
      unit_num = cp_logger_get_default_io_unit(logger)
    ENDIF
     IF (unit_num > 0 .AND. cp_logger_would_log(logger, cp_note_level)) THEN
      SELECT CASE (pol_sched)
      CASE (nosched)
          WRITE(unit_num,'()')
          WRITE(unit_num,'(T2,A)') "SCHED | No scheduling"
      CASE (linear)
          WRITE(unit_num,'()')
          WRITE(unit_num,'(T2,A)') "SCHED | Applying linear scheduling"
      CASE (scatter)
          WRITE(unit_num,'()')
          WRITE(unit_num,'(T2,A)') "SCHED | Applying scatter scheduling"
      CASE (group)
          WRITE(unit_num,'()')
          WRITE(unit_num,'(T2,A)') "SCHED | Applying group scheduling"
      CASE DEFAULT
          WRITE(unit_num,'()')
          WRITE(unit_num,'(T2,A)') "SCHED | No valid scheduling"
      END SELECT
     ENDIF
   ELSE
        CALL cp_assert (.FALSE.,&
        cp_warning_level, cp_unimplemented_error_nr, routineP,&
        "Need hwloc/libnuma library to query machine architecture.")
   ENDIF
   END SUBROUTINE cp_ma_print_config

! *****************************************************************************
!> \brief Prints the Machine Architecture
! *****************************************************************************
  SUBROUTINE cp_ma_print_machine(unit_nr, error)
    INTEGER, INTENT(IN), OPTIONAL            :: unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_print_machine', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_num
    LOGICAL                                  :: print_branch, print_full, &
                                                print_resume
    TYPE(cp_logger_type), POINTER            :: logger

    print_resume = ma_get_conf_print_resume()
    print_full = ma_get_conf_print_full()
    print_branch = ma_get_conf_print_branch()

   logger => cp_error_get_logger(error)

   IF(PRESENT (unit_nr)) THEN
     unit_num = unit_nr
   ELSE
     unit_num = cp_logger_get_default_io_unit(logger)
   ENDIF
   IF (unit_num > 0 .AND. cp_logger_would_log(logger, cp_note_level)) THEN

        SELECT CASE (ma_interface)
        CASE (ma_int_hwloc)
          IF (print_resume) CALL ma_show_topology(unit_num)
          IF (print_branch) CALL ma_show_machine_branch(unit_num)
          IF (print_full)   CALL ma_show_machine_full(unit_num)
        CASE (ma_int_libnuma)
          IF (print_resume) CALL ma_show_topology(unit_num)
        CASE (ma_int_none)
          IF (print_resume .OR. print_branch .OR. print_full) THEN
           CALL cp_assert (.FALSE.,cp_warning_level, cp_unimplemented_error_nr,&
                routineP, "Need hwloc/libNuma library to query machine architecture.")
          ENDIF
        END SELECT
   ENDIF
  END SUBROUTINE cp_ma_print_machine


! *****************************************************************************
!> \brief Prints where the process is running - Also print the memory policy
!> \brief associated with the process
! *****************************************************************************
  SUBROUTINE cp_ma_run_on(para_env, error, unit_nr)
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER, OPTIONAL                        :: unit_nr

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_run_on', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_num
    LOGICAL                                  :: print_mem, print_proc
    TYPE(cp_logger_type), POINTER            :: logger

    print_proc = ma_get_conf_print_proc()
    print_mem = (ma_get_nnodes() .GT. 0 .AND. ma_get_conf_mempol() .NE. 'N')

    logger => cp_error_get_logger(error)

    IF(PRESENT (unit_nr)) THEN
     unit_num = unit_nr
    ELSE
     unit_num = cp_logger_get_default_io_unit(logger)
    ENDIF

    IF ( has_ma ) THEN
      CALL ma_print_proc_affinity(print_proc, print_mem, unit_num)
    ELSE
     CALL cp_assert (.NOT. print_proc,&
        cp_warning_level, cp_unimplemented_error_nr, routineP,&
        "Need hwloc/libnuma library to query machine architecture.")
    ENDIF
  END SUBROUTINE cp_ma_run_on

! *****************************************************************************
!> \brief Get the number of neighbors for a process. Neighbors is the process
!> \brief running on the same machine
! *****************************************************************************
  FUNCTION cp_ma_get_neighbors(error) RESULT (nr_neighbors)
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER                                  :: nr_neighbors

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_get_neighbors', &
      routineP = moduleN//':'//routineN

   nr_neighbors = ma_get_neighbors()

  END FUNCTION cp_ma_get_neighbors

! *****************************************************************************
!> \brief Set the number of neighbors for a process. Neighbors is the process
!> \brief running on the same machine
! *****************************************************************************
  SUBROUTINE cp_ma_set_neighbors(para_env, error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_set_neighbors', &
      routineP = moduleN//':'//routineN

    TYPE(ma_mp_type)                         :: ma_env

    ma_env%numproc = para_env%num_pe
    ma_env%mp_group = para_env%group
    ma_env%myid = r_pid
    ma_env%myproc = para_env%mepos

    IF ( has_ma ) THEN
        CALL ma_set_neighbors(ma_env)
    ELSE
       CALL cp_assert (.FALSE.,&
        cp_warning_level, cp_unimplemented_error_nr, routineP,&
        "Need hwloc/libnuma library to query machine architecture.")
     ENDIF
  END SUBROUTINE cp_ma_set_neighbors

! *****************************************************************************
!> \brief Prints where the threads are running
! Note: Should be used with the cp_ma_sched and outside of a parallel OMP
! *****************************************************************************
  SUBROUTINE cp_ma_thread_run_on(unit_nr, error)
    INTEGER, INTENT(IN), OPTIONAL            :: unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_thread_run_on', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=1)                         :: sched_thread
    INTEGER                                  :: unit_num
    LOGICAL                                  :: print_thread
    TYPE(cp_logger_type), POINTER            :: logger

    print_thread = ma_get_conf_print_thread()
    sched_thread = ma_get_conf_sched()

    logger => cp_error_get_logger(error)

   IF(PRESENT (unit_nr)) THEN
     unit_num = unit_nr
   ELSE
     unit_num = cp_logger_get_default_io_unit(logger)
   ENDIF

  IF (ma_interface .NE. ma_int_none) THEN
     IF (print_thread) THEN
      CALL ma_thread_running_on(unit_num)
     ENDIF
  ELSE
     CALL cp_assert (.NOT. print_thread,&
       cp_warning_level, cp_unimplemented_error_nr, routineP,&
       "Need hwloc library to query machine architecture.")
  ENDIF
  END SUBROUTINE cp_ma_thread_run_on

! *****************************************************************************
!> \brief Prints where the threads are running
! Note: must be used within a parallel OMP region
! *****************************************************************************
  SUBROUTINE cp_ma_current_thread_run(id, unit_nr, error)
    INTEGER                                  :: id
    INTEGER, INTENT(IN), OPTIONAL            :: unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_current_thread_run', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_num
    LOGICAL                                  :: print_thread_cur
    TYPE(cp_logger_type), POINTER            :: logger

    print_thread_cur = ma_get_conf_print_thread_cur()

    logger => cp_error_get_logger(error)

   IF(PRESENT (unit_nr)) THEN
     unit_num = unit_nr
   ELSE
     unit_num = cp_logger_get_default_io_unit(logger)
   ENDIF

   IF (unit_num > 0 .AND. cp_logger_would_log(logger, cp_note_level)) THEN
    IF (ma_interface .NE. ma_int_none) THEN
     IF (print_thread_cur) THEN
        CALL ma_current_thread_run(id, unit_num)
     ENDIF
    ELSE
       CALL cp_assert (.NOT. print_thread,&
       cp_warning_level, cp_unimplemented_error_nr, routineP,&
       "Need hwloc library to query machine architecture.")
    ENDIF
   ENDIF
  END SUBROUTINE cp_ma_current_thread_run

! *****************************************************************************
!> \brief Apply a memory binding strategy
! Note: set the configuration keywords in the input file
! *****************************************************************************
  SUBROUTINE cp_ma_mempol(error)
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_mempol', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=1)                         :: mempol
    INTEGER                                  :: mem_pol

    mempol = ma_get_conf_mempol()

    IF (ma_interface .NE. ma_int_none) THEN
      IF (mempol .EQ. 'D' .OR. mempol .EQ. 'd') THEN
           mem_pol = os
      ELSE IF (mempol .EQ. 'L' .OR. mempol .EQ. 'l') THEN
           mem_pol = local
      ELSE IF (mempol .EQ. 'I' .OR. mempol .EQ. 'i') THEN
           mem_pol = interleave
      ELSE IF (mempol .EQ. 'M' .OR. mempol .EQ. 'm') THEN
           mem_pol = manual
      ELSE
           mem_pol = none_pol
      ENDIF

    IF (mem_pol .NE. none_pol) THEN
       isdefault = .FALSE.

       CALL ma_set_mempol(mem_pol)
    ENDIF

    ELSE
      CALL cp_assert (.FALSE.,&
       cp_warning_level, cp_unimplemented_error_nr, routineP,&
       "Need hwloc or libnuma library to query machine architecture.")
    ENDIF
  END SUBROUTINE cp_ma_mempol

! *****************************************************************************
!> \brief Get or apply a process scheduling strategy
! Note: set the configuration keywords in the input file
! *****************************************************************************
  SUBROUTINE cp_ma_mpi_sched(error)
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_mpi_sched', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=1)                         :: sched_mpi
    INTEGER                                  :: pol_sched
    LOGICAL                                  :: comm_thread

    comm_thread = ma_get_conf_comm_thread()
    sched_mpi = ma_get_conf_mpiSched()

    IF (ma_interface .NE. ma_int_none) THEN
      IF (sched_mpi .EQ. 'L' .OR. sched_mpi .EQ. 'l') THEN
           pol_sched = linear
      ELSE IF (sched_mpi .EQ. 'S' .OR. sched_mpi .EQ. 's') THEN
           pol_sched = scatter
      ELSE IF (sched_mpi .EQ. 'G' .OR. sched_mpi .EQ. 'g') THEN
           pol_sched = group
      ELSE IF (sched_mpi .EQ. 'M' .OR. sched_mpi .EQ. 'm') THEN
           pol_sched = manual
      ELSE IF (sched_mpi .EQ. 'D' ) THEN
           pol_sched = def
      ELSE
           pol_sched = nosched
      ENDIF
      ! confirm the inputs for manual mapping
      IF (pol_sched .EQ. manual .AND. &
         .NOT. ma_valid_inputs(ma_get_neighbors())) THEN
         CALL stop_program(routineN,moduleN,__LINE__,&
              "Need to specify Processors and memory for all MPI")
      ELSE
        IF (pol_sched .NE. def )  CALL ma_sched_mpi(pol_sched)
      END IF
    ELSE
      CALL cp_assert (.FALSE.,&
       cp_warning_level, cp_unimplemented_error_nr, routineP,&
       "Need hwloc/libnuma library to query machine architecture.")
    ENDIF
  END SUBROUTINE cp_ma_mpi_sched



! *****************************************************************************
!> \brief Get or apply a thread scheduling strategy
! Note: set the configuration keywords in the input file
! *****************************************************************************
  SUBROUTINE cp_ma_thread_sched(id, error)
    INTEGER, INTENT(IN)                      :: id
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_thread_sched', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=1)                         :: sched_thread
    INTEGER                                  :: pol_sched

    sched_thread = ma_get_conf_sched()

    IF (ma_interface .NE. ma_int_none) THEN
      IF (sched_thread .EQ. 'L' .OR. sched_thread .EQ. 'l') THEN
           pol_sched = linear
      ELSE IF (sched_thread .EQ. 'S' .OR. sched_thread .EQ. 's') THEN
           pol_sched = scatter
      ELSE IF (sched_thread .EQ. 'G' .OR. sched_thread .EQ. 'g') THEN
           pol_sched = group
      ELSE IF (sched_thread .EQ. 'M' .OR. sched_thread .EQ. 'm') THEN
           pol_sched = manual
      ELSE IF (sched_thread .EQ. 'D' ) THEN
           pol_sched = def
      ELSE
           pol_sched = nosched
      ENDIF
      IF (pol_sched .NE. def) CALL ma_sched_threads(id, pol_sched)
    ELSE
      CALL cp_assert (.FALSE.,&
       cp_warning_level, cp_unimplemented_error_nr, routineP,&
       "Need hwloc/libnuma library to query machine architecture.")
   ENDIF

  END SUBROUTINE cp_ma_thread_sched

! *****************************************************************************
!> \brief Verify the thread scheduling
! *****************************************************************************
  SUBROUTINE cp_ma_verify_place(id, unit_nr, error)
    INTEGER, INTENT(IN)                      :: id
    INTEGER, INTENT(IN), OPTIONAL            :: unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_verify_place', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=1)                         :: sched_thread
    INTEGER                                  :: unit_num
    LOGICAL                                  :: print_thread
    TYPE(cp_logger_type), POINTER            :: logger

    print_thread = ma_get_conf_print_thread()
    sched_thread = ma_get_conf_sched()

    logger => cp_error_get_logger(error)

    IF(PRESENT (unit_nr)) THEN
      unit_num = unit_nr
    ELSE
      unit_num = cp_logger_get_default_io_unit(logger)
    ENDIF
    IF (unit_num > 0 .AND. cp_logger_would_log(logger, cp_note_level)) THEN
    IF (ma_interface .NE. ma_int_none) THEN
       IF (sched_thread .EQ. 'L' .OR. sched_thread .EQ. 'l' .OR. &
           sched_thread .EQ. 'S' .OR. sched_thread .EQ. 's' .OR. &
           sched_thread .EQ. 'G' .OR. sched_thread .EQ. 'g' .OR. &
           sched_thread .EQ. 'N' .OR. sched_thread .EQ. 'n') THEN
         IF (print_thread) THEN
           CALL ma_verify_place(id, unit_num)
         ENDIF
       ENDIF
     ELSE
      CALL cp_assert (.NOT. print_thread,&
       cp_warning_level, cp_unimplemented_error_nr, routineP,&
       "Need hwloc library to query machine architecture.")
     ENDIF
    ENDIF
  END SUBROUTINE cp_ma_verify_place

  SUBROUTINE cp_ma_set_mpi_reordering(mpi_reordering)
    CHARACTER(len=1), INTENT(IN)             :: mpi_reordering

    CALL ma_set_conf_mpi_reordering(mpi_reordering)

  END SUBROUTINE cp_ma_set_mpi_reordering

  SUBROUTINE cp_ma_print_strategy(unit_nr, error)
    INTEGER, INTENT(IN), OPTIONAL            :: unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_print_strategy', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: mpi_reorder
    INTEGER                                  :: unit_num
    TYPE(cp_logger_type), POINTER            :: logger

    mpi_reorder = ma_get_conf_mpi_reordering()

    logger => cp_error_get_logger(error)

    IF(PRESENT (unit_nr)) THEN
      unit_num = unit_nr
    ELSE
      unit_num = cp_logger_get_default_io_unit(logger)
    ENDIF

    IF (unit_num > 0 .AND. cp_logger_would_log(logger, cp_note_level)) THEN
      IF (mpi_reorder .EQ. 'H' .OR. mpi_reorder .EQ. 'h') THEN
        WRITE(unit_num,'(A)') " MPI REORDERING| Hilbert curve"
      ELSE IF (mpi_reorder .EQ. 'P' .OR. mpi_reorder .EQ. 'p') THEN
        WRITE(unit_num,'(A)') " MPI REORDERING| Peano curve"
      ELSE IF (mpi_reorder .EQ. 'D' .OR. mpi_reorder .EQ. 'd') THEN
        WRITE(unit_num,'(A)') " MPI REORDERING| Default"
      ELSE IF (mpi_reorder .EQ. 'R' .OR. mpi_reorder .EQ. 'r') THEN
        WRITE(unit_num,'(A)') " MPI REORDERING| Round-robin"
      ELSE IF (mpi_reorder .EQ. 'F' .OR. mpi_reorder .EQ. 'f') THEN
        WRITE(unit_num,'(A)') " MPI REORDERING| Hilbert-Peano"
      ELSE IF (mpi_reorder .EQ. 'C' .OR. mpi_reorder .EQ. 'c') THEN
        WRITE(unit_num,'(A)') " MPI REORDERING| Cannon Heuristic"
      ELSE IF (mpi_reorder .EQ. 'O' .OR. mpi_reorder .EQ. 'o') THEN
        WRITE(unit_num,'(A)') " MPI REORDERING| Own"
      ELSE
        WRITE(unit_num,'(A)') " MPI REORDERING| No strategy selected"
      ENDIF
        WRITE(unit_num,'(A)')
    ENDIF

  END SUBROUTINE cp_ma_print_strategy

  SUBROUTINE cp_ma_mpi_reorder_strategy(mapping_method,strategy)
    INTEGER, INTENT(IN)                      :: mapping_method
    CHARACTER, INTENT(OUT)                   :: strategy

   SELECT CASE(mapping_method)
     CASE(1)
       strategy = 'D'
     CASE(2)
       strategy = 'O'
     CASE(3)
      strategy = 'H'
     CASE(4)
      strategy = 'P'
     CASE(5)
      strategy = 'R'
     CASE(6)
      strategy = 'F'
     CASE(7)
      strategy = 'C'
     CASE default
      strategy = 'X'
   END SELECT
  END SUBROUTINE cp_ma_mpi_reorder_strategy

  SUBROUTINE cp_ma_mpi_reordering (mp_comm,mp_reorder,error)
    INTEGER, INTENT(INOUT)                   :: mp_comm
    CHARACTER(len=1), INTENT(IN), OPTIONAL   :: mp_reorder
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_mpi_reordering', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=1)                         :: mpi_reorder
    INTEGER                                  :: handle, reorder

    CALL timeset(routineN,handle)

    IF(.NOT. PRESENT(mp_reorder)) THEN
       mpi_reorder = ma_get_conf_mpi_reordering()
    ELSE
       mpi_reorder = mp_reorder
    ENDIF

    IF (has_ma_topology) THEN
      IF (mpi_reorder .EQ. 'H' .OR. mpi_reorder .EQ. 'h') THEN
        reorder = hilbert 
      ELSE IF (mpi_reorder .EQ. 'P' .OR. mpi_reorder .EQ. 'p') THEN
        reorder = peano
      ELSE IF (mpi_reorder .EQ. 'N' .OR. mpi_reorder .EQ. 'n') THEN
        reorder = snake      
      ELSE IF (mpi_reorder .EQ. 'D' .OR. mpi_reorder .EQ. 'd') THEN
        reorder = packed
      ELSE IF (mpi_reorder .EQ. 'R' .OR. mpi_reorder .EQ. 'r') THEN
        reorder = round_robin
      ELSE IF (mpi_reorder .EQ. 'F' .OR. mpi_reorder .EQ. 'f') THEN
        reorder = hilbert_peano
      ELSE IF (mpi_reorder .EQ. 'S' .OR. mpi_reorder .EQ. 's') THEN
        reorder = switch
      ELSE IF (mpi_reorder .EQ. 'C' .OR. mpi_reorder .EQ. 'c') THEN
        reorder = cannon
      ELSE IF (mpi_reorder .EQ. 'O' .OR. mpi_reorder .EQ. 'o') THEN
        reorder = own
      ELSE
        reorder = none_order   
      ENDIF
     ELSE
         WRITE(*,*) "No support for process mapping"  
         reorder = packed
     ENDIF
     IF (reorder .NE. none_order) THEN
       CALL ma_mpi_reordering(reorder,mp_comm)
     ENDIF

     CALL timestop(handle) 
  END SUBROUTINE cp_ma_mpi_reordering


! *****************************************************************************
!> \brief Initializes Machine Architecture
! *****************************************************************************
  SUBROUTINE cp_ma_init_lib (init_mpi, group, error)
    LOGICAL, INTENT(IN)                      :: init_mpi
    INTEGER, INTENT(IN)                      :: group
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_init_lib', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: core, istat, mynode, node, &
                                                num_threads, numnode, unit_num
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(ma_mp_type)                         :: ma_env

    logger => cp_error_get_logger(error)

    IF ( isconfigured ) RETURN
    IF ( has_ma ) THEN
        ma_interface = ma_int_none
        IF (ma_hwloc) THEN
                ma_interface = ma_int_hwloc
        ELSE
          IF (ma_libnuma) THEN
             ma_interface = ma_int_libnuma
          ENDIF
        ENDIF
        num_threads=1

!$omp parallel
!$ num_threads = omp_get_num_threads()
!$omp end parallel

      CALL mp_environ(numnode, mynode, group)
      IF(init_mpi) THEN
        ma_env%numproc = numnode
        ma_env%mp_group = group
        ma_env%myid = r_pid
        ma_env%myproc = mynode
        ma_env%all_proc = .TRUE.
      ELSE ! the communicator is not usefull
        ma_env%numproc = numnode
        ma_env%mp_group = group
        ma_env%myid = r_pid
        ma_env%myproc = mynode
        ma_env%all_proc = .FALSE.
      ENDIF
       istat = ma_init_machine(core,node)
       IF (istat /= 0) THEN
           CALL stop_program(routineN,moduleN,__LINE__,&
               "Initialization of HWLOC/LIBNUMA failed!")
       ENDIF

       istat = ma_init_affinity(ma_env, num_threads,core,node)
       IF (istat /= 0) THEN
           CALL stop_program(routineN,moduleN,__LINE__,&
               "Initialization of AFFINITY failed!")
       ENDIF

       unit_num = cp_logger_get_default_io_unit(logger)

       !check the number of MPI and GPUs
        IF (ma_MPI_nGpus()) THEN
          IF (unit_num > 0 .AND. cp_logger_would_log(logger, cp_note_level)) THEN
             WRITE(unit_num,'(A)') "WARNING: Number of MPI greater than GPUs"
             WRITE(unit_num,'()')
           ENDIF
        ENDIF
        isconfigured = .TRUE.
    ENDIF


  END SUBROUTINE cp_ma_init_lib

! *****************************************************************************
!> \brief Performs the default binding of process
! *****************************************************************************
  SUBROUTINE cp_ma_default_affinity(isallproc)
    LOGICAL                                  :: isallproc

   IF(isallproc) THEN
      CALL ma_set_default_affinity()
      isdefault = isallproc
   ENDIF
  END SUBROUTINE cp_ma_default_affinity

! *****************************************************************************
!> \brief Finalizes Machine Architecture
! *****************************************************************************
  SUBROUTINE cp_ma_finalize_lib (error)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_finalize_lib', &
      routineP = moduleN//':'//routineN

    IF(has_ma) THEN
      CALL ma_finalize_machine()
      CALL ma_finalize_affinity()
    ELSE IF (has_ma_topology) THEN
      CALL ma_destroy_topology ()
    ENDIF
  END SUBROUTINE cp_ma_finalize_lib


END MODULE cp_ma_interface
