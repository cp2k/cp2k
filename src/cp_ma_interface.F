!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
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
  USE cp2k_info,                       ONLY: r_host_name,&
                                             r_pid
  USE cp_error_handling
  USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit,&
                                             cp_logger_type,&
                                             cp_logger_would_log,&
                                             cp_note_level,&
                                             cp_warning_level
  USE cp_output_handling
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             int_size
  USE ma_config
  USE machine_architecture
  USE machine_architecture_types
  USE message_passing,                 ONLY: mp_environ,&
                                             mp_sum
  USE string_utilities,                ONLY: ascii_to_string,&
                                             string_to_ascii
  USE termination,                     ONLY: stop_memory

!$ USE OMP_LIB

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_ma_interface'

  ! Interface to libma
  PUBLIC :: cp_ma_config, cp_ma_init_lib, cp_ma_finalize_lib
  PUBLIC :: cp_ma_run_on, cp_ma_thread_run_on 
  PUBLIC :: cp_ma_current_thread_run, cp_ma_set_neighbors
  PUBLIC :: cp_ma_thread_sched, cp_ma_verify_place 
  PUBLIC :: cp_ma_mpi_sched
  PUBLIC :: cp_ma_print_machine
  PUBLIC :: cp_ma_get_cuDev 

CONTAINS

! *****************************************************************************
!> \brief Configures options for Machine Architecture  
! *****************************************************************************
  SUBROUTINE cp_ma_config(root_section, error)
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_config', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=1)                         :: mpi_sched, use_sched
    LOGICAL                                  :: print_branch, print_full, &
                                                print_proc, print_resume, &
                                                print_thread, print_thread_cur
    TYPE(section_vals_type), POINTER         :: ma_section

!TODO 
!CALL cp_error_set (routineN, error_handle, error)

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
   
    CALL ma_set_conf_print_branch (print_branch)
    CALL ma_set_conf_print_resume (print_resume)
    CALL ma_set_conf_print_full (print_full)
    CALL ma_set_conf_print_proc (print_proc)
    CALL ma_set_conf_print_thread (print_thread)
    CALL ma_set_conf_print_thread_cur (print_thread_cur)
    CALL ma_set_conf_mempol (use_mempol)
    CALL ma_set_conf_sched (use_sched)
    CALL ma_set_conf_mpiSched (mpi_sched)
  END SUBROUTINE cp_ma_config

! *****************************************************************************
!> \brief Get the cuda devices for the MPIs
! *****************************************************************************
  SUBROUTINE cp_ma_get_cuDev(cuDev, myRank)
    INTEGER                                  :: cuDev, myRank

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_get_cuDev', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_num
    TYPE(cp_logger_type), POINTER            :: logger

    IF ( has_ma ) THEN
            unit_num = cp_logger_get_default_io_unit(logger)
            cuDev = ma_getcuDev(myRank)

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
          WRITE(unit_num,'(T2,A)') "SCHED | No thread scheduling"
      CASE (linear)
          WRITE(unit_num,'()')
          WRITE(unit_num,'(T2,A)') "SCHED | Applying linear thread scheduling"
      CASE DEFAULT
          WRITE(unit_num,'()')
          WRITE(unit_num,'(T2,A)') "SCHED | No valid thread scheduling"
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
          CALL cp_assert (.FALSE.,&
          cp_warning_level, cp_unimplemented_error_nr, routineP,&
          "Need hwloc/libNuma library to query machine architecture.")
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
    INTEGER                                  :: unit_nr

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_run_on', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=10)                        :: mempolicy
    INTEGER                                  :: ipe, istat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: all_core, all_memory, &
                                                all_node, all_node_mem, &
                                                all_pid
    LOGICAL                                  :: print_proc

    print_proc = ma_get_conf_print_proc() 
    IF ( has_ma ) THEN 
      ALLOCATE (all_core(para_env%num_pe),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"all_core",&
            para_env%num_pe*int_size)
      END IF
      all_core(:) = 0

      ALLOCATE (all_node(para_env%num_pe),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"all_node",&
            para_env%num_pe*int_size)
      END IF
      all_node(:) = 0

      ALLOCATE (all_memory(para_env%num_pe),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"all_memory",&
            para_env%num_pe*int_size)
      END IF
      all_memory(:) = 0

      ALLOCATE (all_node_mem(para_env%num_pe),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"all_node_mem",&
            para_env%num_pe*int_size)
      END IF
      all_node_mem(:) = 0   

      !Where process is running
      CALL ma_running_on(all_core(para_env%mepos+1), all_node(para_env%mepos+1))
 
   IF (print_proc) THEN
! ***** Print where process are running!

      CALL mp_sum(all_core,para_env%group)
      CALL mp_sum(all_node,para_env%group)

      ALLOCATE (all_pid(para_env%num_pe),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"all_pid",&
            para_env%num_pe*int_size)
      END IF
      all_pid(:) = 0
      all_pid(para_env%mepos+1) = r_pid
      CALL mp_sum(all_pid,para_env%group)

     IF ( ma_get_conf_mempol() .NE. 'N' ) THEN
      CALL ma_get_mempol(all_node_mem(para_env%mepos+1),all_memory(para_env%mepos+1))
      CALL mp_sum(all_memory,para_env%group)
      CALL mp_sum(all_node_mem,para_env%group)
     ENDIF      
 
      IF (unit_nr>0) THEN

        WRITE (unit_nr,"(T2,A)") "SCHED| Processes mapping"
        DO ipe=1,para_env%num_pe
          IF (all_node(ipe) .EQ. -1) THEN
            WRITE(unit_nr,"(T2,A,I8,A,I8,5X,30X,A,I8)")&
              "Process ", ipe-1, " ", all_pid(ipe), " runing on core ", all_core(ipe)
          ELSE
            WRITE (unit_nr,"(T2,A,I8,A,I8,5X,A,I6,7X,A,I8)")&
              "Process ", ipe-1, " ", all_pid(ipe)," running on NUMA node ",&
               all_node(ipe), " core ", all_core(ipe)
          ENDIF
        END DO
        
        WRITE (unit_nr,"(T2,A)") ""
        IF (ma_get_nnodes() .GT. 0 .AND. ma_get_conf_mempol() .NE. 'N') THEN 
           WRITE (unit_nr,"(T2,A)") "MEMORY| Processes memory mapping"
           DO ipe=1,para_env%num_pe
              CALL ma_string_mempol(all_memory(ipe), mempolicy) 
              IF (all_node_mem(ipe) .EQ. -1) THEN
                 WRITE (unit_nr,"(T2,A,I8,A,I8,5X,A,A,A)")&
                "Process ", ipe-1, " ", all_pid(ipe)," memory policy ",&
                 TRIM(mempolicy), " all nodes"
              ELSE 
                WRITE (unit_nr,"(T2,A,I8,A,I8,5X,A,A,7X,A,I8)")&
                "Process ", ipe-1, " ", all_pid(ipe)," memory policy ",&
                 TRIM(mempolicy), " node ", all_node_mem(ipe)
              ENDIF
           END DO
           WRITE (unit_nr,"(T2,A)") "" 
        ENDIF
      END IF

      DEALLOCATE (all_memory,STAT=istat)
      IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,"all_memory")
      END IF
      DEALLOCATE (all_node_mem,STAT=istat)
      IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,"all_node_memory")
      END IF
      DEALLOCATE (all_pid,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"all_pid")
      END IF

  ENDIF

      DEALLOCATE (all_core,STAT=istat)
      IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,"all_core")
      END IF
      DEALLOCATE (all_node,STAT=istat)
      IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,"all_node")
      END IF  
  ELSE
     CALL cp_assert (.NOT. print_proc,&
        cp_warning_level, cp_unimplemented_error_nr, routineP,&
        "Need hwloc/libnuma library to query machine architecture.")
  ENDIF
  END SUBROUTINE cp_ma_run_on

! *****************************************************************************
!> \brief Set the number of neighbors for a process. Neighbors is the process
!> \brief running on the same machine
! *****************************************************************************
  SUBROUTINE cp_ma_set_neighbors(ma_env, error)
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    TYPE(ma_mp_type)                         :: ma_env

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_set_neighbors', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: string
    INTEGER                                  :: istat, jpe, nr_occu
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: all_host

    IF ( has_ma ) THEN
        ALLOCATE (all_host(30,ma_env%numproc),STAT=istat)
        IF (istat /= 0) THEN
                CALL stop_memory(routineN,moduleN,__LINE__,"all_host",&
                30*ma_env%numproc*int_size)
        END IF
        all_host(:,:) = 0

        CALL string_to_ascii(r_host_name,all_host(:,ma_env%myproc+1))
        CALL mp_sum(all_host,ma_env%mp_group)

        nr_occu = 0
        DO jpe=1,ma_env%numproc
                CALL ascii_to_string(all_host(:,jpe),string)
                IF  (string .EQ. r_host_name) THEN
                nr_occu = nr_occu +1
                ENDIF
        END DO

        CALL ma_set_neighbors(nr_occu)
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
   IF (sched_thread .EQ. 'L' .OR. sched_thread .EQ. 'l' .OR. &
       sched_thread .EQ. 'N' .OR. sched_thread .EQ. 'n') THEN
     IF (print_thread) THEN
      CALL ma_thread_running_on(unit_num)
     ENDIF
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
      ELSE 
           mem_pol = -1        
      ENDIF
      CALL ma_set_mempol(mem_pol)
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

    sched_mpi = ma_get_conf_mpiSched()

    IF (ma_interface .NE. ma_int_none) THEN    
      IF (sched_mpi .EQ. 'L' .OR. sched_mpi .EQ. 'l') THEN
           pol_sched = linear
           CALL ma_sched_mpi(linear)
      ELSE
           pol_sched = nosched
           CALL ma_sched_mpi(nosched)
      ENDIF
    ELSE
      CALL cp_assert (pol_sched .LT. 0,&
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
           CALL ma_sched_threads(id, linear)
      ELSE
           pol_sched = nosched
           CALL ma_sched_threads(id, nosched)
      ENDIF
    ELSE
      CALL cp_assert (pol_sched .LT. 0,&
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


! *****************************************************************************
!> \brief Initializes Machine Architecture
! *****************************************************************************
  SUBROUTINE cp_ma_init_lib (group, error)
    INTEGER, INTENT(IN)                      :: group
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_init_lib', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat, mynode, num_threads, &
                                                numnode
    TYPE(ma_mp_type)                         :: ma_env

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
        CALL mp_environ (numnode, mynode, group)

        ma_env%myproc = mynode
        ma_env%numproc = numnode
        ma_env%mp_group = group

        istat = ma_init (ma_env, num_threads)
        CALL cp_ma_set_neighbors(ma_env, error)
    ELSE
        CALL cp_assert (.FALSE.,&
       cp_warning_level, cp_unimplemented_error_nr, routineP,&
       "Need hwloc/libnuma library to query machine architecture.")    
    ENDIF
    
  END SUBROUTINE cp_ma_init_lib

 
! *****************************************************************************
!> \brief Finalizes Machine Architecture
! *****************************************************************************
  SUBROUTINE cp_ma_finalize_lib (error)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ma_finalize_lib', &
      routineP = moduleN//':'//routineN

    CALL ma_finalize ()
  END SUBROUTINE cp_ma_finalize_lib


END MODULE cp_ma_interface
