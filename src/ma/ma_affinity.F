!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Affinity support for cp2k
!>          Interface to the HWLOC library
!> \author  Christiane Pousa Ribeiro
!> \date    2011-11-08
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2011-11-08
! *****************************************************************************
MODULE ma_affinity
  USE ISO_C_BINDING
  USE ma_config,                       ONLY: &
       hasnet, isdefault, ma_get_conf_mempol, ma_get_conf_mpisched, memory, &
       print_proc, print_thread, process, thread
  USE ma_errors,                       ONLY: ma_error_allocation,&
                                             ma_error_stop
  USE ma_kinds,                        ONLY: default_string_length,&
                                             int_4
  USE machine,                         ONLY: default_output_unit
  USE machine_architecture,            ONLY: &
       ma_get_core_node, ma_get_id, ma_get_mycore, ma_get_mynode, &
       ma_get_ncores, ma_get_netDev, ma_get_nnetDev, ma_get_nnodes, &
       ma_get_node_netDev, ma_get_proc_core, ma_get_proc_node, &
       ma_get_thread_id, ma_hw_get_mempol, ma_hw_set_mempol, &
       ma_my_first_core, ma_set_core, ma_set_first_core, ma_set_proc_core, &
       ma_set_proc_cores, ma_set_proc_node, ma_set_thread_allnodes, &
       ma_set_thread_cores, ma_set_thread_node, topology
  USE machine_architecture_types,      ONLY: &
       def, group, has_mpi, interleave, linear, local, ma_mp_type, &
       ma_process, manual, mpi, nosched, os, scatter, thread_inf, threads
  USE machine_architecture_utils,      ONLY: ascii_to_string,&
                                             integer_to_string,&
                                             ma_create_process,&
                                             ma_pack_threads,&
                                             ma_unpack_threads,&
                                             string_to_ascii
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_proc_name,&
                                             mp_sum

!$ USE OMP_LIB
 IMPLICIT NONE

 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ma_affinity'

  PUBLIC :: ma_init_affinity, ma_finalize_affinity

  PUBLIC :: ma_sched_threads, ma_sched_mpi

  PUBLIC :: ma_verify_place
  PUBLIC :: ma_thread_running_on, ma_current_thread_run

  PUBLIC :: ma_set_mempol

  PUBLIC :: ma_set_gpu_affinity
  PUBLIC :: ma_set_net_affinity

  PUBLIC :: ma_set_neighbors, ma_get_neighbors, ma_MPI_nGpus

  PUBLIC :: ma_print_proc_affinity, ma_set_default_affinity

! These are for Affinity module internal use.
!

  TYPE(ma_process)                               :: my_info
! The global view of the thread mapping on the machine
  TYPE(thread_inf), DIMENSION(:,:), ALLOCATABLE  :: thread_mapping


CONTAINS

! *****************************************************************************
!> \brief Initialize the affinity strategy information
!> \param ma_env            The parallel environment
!> \param num_threads       The number of threads of each process 
!> \param core ...
!> \param node ...
!> \retval istat ...
! *****************************************************************************
INTEGER FUNCTION ma_init_affinity(ma_env, num_threads, core, node) RESULT(istat)
  INTEGER                          :: num_threads, core, node
  TYPE(ma_mp_type)                 :: ma_env

  istat = -1
  ALLOCATE(thread_mapping(ma_env%numproc,num_threads),STAT=istat)
  IF (istat /= 0) CALL ma_error_stop(ma_error_allocation)
  CALL ma_create_process(my_info,ma_env,num_threads,core,node)
  CALL ma_set_neighbors(ma_env)
END FUNCTION ma_init_affinity

! *****************************************************************************
!> \brief Finalizes all the affinity strategy information
! *****************************************************************************
SUBROUTINE ma_finalize_affinity()
#if defined (__LIBNUMA) || defined (__HWLOC)
   INTEGER  :: istat
   DEALLOCATE(thread_mapping,STAT=istat)
   IF (istat /= 0) CALL ma_error_stop(ma_error_allocation)
#endif
END SUBROUTINE ma_finalize_affinity

! *****************************************************************************
!> \brief Default affinity strategy for CP2K. Place MPI process close to the
!> devices to avoid overheads. It maps process, threads, memory  and gpu.
! *****************************************************************************
SUBROUTINE ma_set_default_affinity()
#if defined (__LIBNUMA) || defined (__HWLOC)
! Place MPI closer to GPU
!$ INTEGER :: id
#if defined (__DBCSR_ACC) || defined (__PW_CUDA)
   INTEGER :: node
   IF (my_info%nr_neighbors .GT. topology%ngpus) THEN
     CALL ma_group_place(my_info%mp_info%myproc, mpi)
   ELSE
     node = ma_get_gpu_node(MOD(my_info%mp_info%myproc,topology%ngpus))
     CALL ma_numa_node_place(node,my_info%mp_info%myproc,mpi)
   ENDIF
   IF (my_info%nr_neighbors .GT. 1) THEN
      CALL ma_set_mempol(local)
   ELSE
      CALL ma_set_mempol(interleave)
   END IF
#else
   CALL ma_group_place(my_info%mp_info%myproc, mpi)
   IF (my_info%nr_neighbors .GT. 1) THEN
      CALL ma_set_mempol(local)
   ELSE
     IF ( my_info%nr_threads .GT. topology%ncores/topology%nnodes) THEN
       CALL ma_set_mempol(interleave)
     END IF
   END IF
#endif
!$omp parallel private(id)
!$ id = omp_get_thread_num()
!$ CALL ma_group_place(id,threads)
!$omp end parallel
#endif
isdefault = .TRUE.
END SUBROUTINE ma_set_default_affinity

! *****************************************************************************
!> \brief To convert a memory policy descriptor to a string
!> \param mempol          Memory policy descriptor 
!> \param mempolicy       String to hold the memory policy name 
! *****************************************************************************
SUBROUTINE ma_string_mempol(mempol, mempolicy)
    INTEGER, INTENT(IN)                      :: mempol
    CHARACTER(LEN=10)                        :: mempolicy

   SELECT CASE (mempol)
   CASE (os)
     mempolicy = 'DEFAULT'
   CASE (local)
     mempolicy = 'LOCAL'
   CASE (interleave)
     mempolicy = 'INTERLEAVE'
   CASE (manual)
     mempolicy = 'MANUAL'
   CASE DEFAULT
     mempolicy = 'NONE'
   END SELECT
END SUBROUTINE ma_string_mempol

! *****************************************************************************
!> \brief Set a memory policy for a process
!> \param mempol Memory policy descriptor
! *****************************************************************************
SUBROUTINE  ma_set_mempol (mempol)
    INTEGER                                  :: mempol

    INTEGER                                  :: istat, realnode, unit_nr

    istat = 0
    unit_nr = default_output_unit

  IF (mempol .EQ. manual) THEN
    realnode = memory(MOD(my_info%mp_info%myproc,my_info%nr_neighbors)+1)
  END IF

#if defined (__HWLOC) && !defined (__LIBNUMA)
 IF (.NOT. isdefault) THEN
  istat = ma_hw_set_mempol(mempol,realnode)
  IF (istat .LT. 0 .AND. my_info%mp_info%myproc .EQ. 0) THEN
   WRITE(unit_nr,'(T2,A,A)') "WARNING: Memory binding not supported", &
            " or can not be enforced - Using OS memory policy"
  END IF
 ENDIF
#endif
#if defined (__LIBNUMA) && !defined (__HWLOC)
 IF (.NOT. isdefault) THEN
  istat = ma_linux_set_mempol(mempol,realnode)
  IF (istat .LT. 0 .AND. my_info%mp_info%myproc .EQ. 0) THEN
   WRITE(unit_nr,'(T2,A,A)') "WARNING: Memory binding not supported", &
            " or can not be enforced - Using OS memory policy"
  END IF
 END IF
#endif

END SUBROUTINE ma_set_mempol

! *****************************************************************************
!> \brief Set a memory policy for a process
!> \param node            Node numa used by the policy 
!> \param mempol          Memory policy descriptor
! *****************************************************************************
SUBROUTINE  ma_get_mempol (node, mempol)
  INTEGER                :: node, mempol

#if defined (__HWLOC) && !defined (__LIBNUMA)
  CALL ma_hw_get_mempol(node,mempol)
#endif
#if defined (__LIBNUMA) && !defined (__HWLOC)
  CALL ma_linux_get_mempol(node,mempol)
#endif

 IF (mempol .EQ. manual) THEN
     node = memory(MOD(my_info%mp_info%myproc,my_info%nr_neighbors)+1)
 ENDIF

END SUBROUTINE ma_get_mempol

! *****************************************************************************
!> \brief Sets the number of neighbors for a process
!> \param ma_env parallel environment
! *****************************************************************************
  SUBROUTINE ma_set_neighbors(ma_env)
    TYPE(ma_mp_type)                         :: ma_env

    CHARACTER(LEN=default_string_length)     :: host_name, string
    INTEGER                                  :: istat, jpe, nr_occu
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: all_host

    IF (ma_env%all_proc) THEN

      ALLOCATE (all_host(default_string_length,ma_env%numproc),STAT=istat)
      IF (istat /= 0) CALL ma_error_stop(ma_error_allocation)

      all_host(:,:) = 0
      CALL mp_proc_name(host_name)

      CALL string_to_ascii(host_name,all_host(:,ma_env%myproc+1))
      CALL mp_sum(all_host,ma_env%mp_group)

      nr_occu = 0
      DO jpe=1,ma_env%numproc
         CALL ascii_to_string(all_host(:,jpe),string)
         IF  (string .EQ. host_name) THEN
              nr_occu = nr_occu +1
         ENDIF
      END DO

      my_info%nr_neighbors = nr_occu
   ELSE
      my_info%nr_neighbors = ma_get_ncores()
   ENDIF
  END SUBROUTINE ma_set_neighbors

! *****************************************************************************
!> \brief Gets the number of neighbors for a process
!> \retval nr_occu ...
! *****************************************************************************
  FUNCTION ma_get_neighbors () RESULT (nr_occu)
    INTEGER                                  :: nr_occu

    nr_occu = my_info%nr_neighbors

  END FUNCTION ma_get_neighbors

! *****************************************************************************
!> \brief Verify if the number of MPI processes is greater than GPUS
!> \retval is_greater ...
! *****************************************************************************
  FUNCTION ma_MPI_nGpus() RESULT(is_greater)
    LOGICAL                                  :: is_greater

    is_greater = .FALSE.
#if defined (__DBCSR_ACC) || defined (__PW_CUDA)
    IF (my_info%nr_neighbors .GT. topology%ngpus ) is_greater = .TRUE.
#endif
 END FUNCTION ma_MPI_nGpus

! *****************************************************************************
!> \brief Get the MPI/threads mapping on the machine cores
!> \param id thread OMP id - MPI rank 
!> \param sched_unit ...
! *****************************************************************************
  SUBROUTINE ma_def_place(id, sched_unit)
    INTEGER                                  :: id, sched_unit

  CALL ma_set_default_affinity()
  isdefault = .TRUE.

  END SUBROUTINE ma_def_place

! *****************************************************************************
!> \brief Get the MPI/threads mapping on the machine cores
!> \param id thread OMP id - MPI rank 
!> \param sched_unit ...
! *****************************************************************************
  SUBROUTINE ma_no_place(id, sched_unit)
    INTEGER                                  :: id, sched_unit

    INTEGER                                  :: core

    core = 0
    isdefault = .FALSE.
#if defined (__HWLOC)  && !defined (__LIBNUMA)
 IF ( sched_unit .EQ. threads ) THEN
  CALL ma_set_thread_cores(topology%ncores,0)
  my_info%threads_info(id+1)%core = ma_get_mycore()
  my_info%threads_info(id+1)%node = ma_get_mynode()
  my_info%threads_info(id+1)%id_omp = id
  my_info%threads_info(id+1)%id_real = ma_get_thread_id()
 ELSE
  CALL ma_set_proc_cores(topology%ncores,0)
  core = ma_get_proc_core(1)
  my_info%threads_info(1)%core = core
  my_info%threads_info(1)%node = ma_get_core_node(core)
  my_info%threads_info(1)%id_omp = id
  my_info%threads_info(1)%id_real = ma_get_id()
  my_info%mp_info%myid = my_info%threads_info(1)%id_real
  my_info%core = core
  my_info%node = my_info%threads_info(1)%node
 END IF
#endif
#if defined (__LIBNUMA) && !defined (__HWLOC)
 IF ( sched_unit .EQ. threads ) THEN
  CALL ma_set_proc_cores(topology%ncores,0)
  my_info%threads_info(id+1)%core = ma_linux_get_core()
  my_info%threads_info(id+1)%node = ma_linux_get_mynode()
  my_info%threads_info(id+1)%id_omp = id
  my_info%threads_info(id+1)%id_real = ma_linux_get_thread_id()
 ELSE
  CALL ma_set_proc_cores(topology%ncores,0)
  core = ma_linux_get_core()
  my_info%threads_info(1)%core = core
  my_info%threads_info(1)%node = ma_linux_get_core_node(core)
  my_info%threads_info(1)%id_omp = id
  my_info%threads_info(1)%id_real = ma_get_id()
  my_info%mp_info%myid = my_info%threads_info(1)%id_real
  my_info%core = core
  my_info%node = my_info%threads_info(1)%node
 END IF
#endif
  ! NOT USING NET AFFINITY
  hasnet = .FALSE.
  END SUBROUTINE ma_no_place

! *****************************************************************************
!> \brief Set the MPI/threads mapping on the machine cores - compact affinty
!> \param id thread OMP id - MPI rank
!> \param sched_unit ...
!> \note The strategy used here is a simple compact distrbution for threads of a
!> \note process. In this case, threads of a process take into account the process
!> \note affinity.
! *****************************************************************************

  SUBROUTINE ma_linear_place(id, sched_unit)
    INTEGER                                  :: id, sched_unit

    INTEGER                                  :: core

    core = 0
    isdefault = .FALSE.

#if defined (__HWLOC) && !defined (__LIBNUMA)
   IF (sched_unit .EQ. threads .AND. my_info%nr_threads .GT. 1) THEN
    core = MOD(my_info%mp_info%myproc*my_info%nr_threads+id,ma_get_ncores())
    CALL ma_set_core(core)
    my_info%threads_info(id+1)%core = core
    my_info%threads_info(id+1)%node = ma_get_mynode()
    my_info%threads_info(id+1)%id_omp = id
    my_info%threads_info(id+1)%id_real = ma_get_thread_id()
   ELSE
    IF(my_info%nr_threads .GT. 1) THEN
      core = MOD(my_info%mp_info%myproc*my_info%nr_threads,ma_get_ncores())
    ELSE
      core = (MOD(my_info%mp_info%myproc,ma_get_ncores()))
    ENDIF
    CALL ma_set_proc_core(core)
    my_info%threads_info(1)%core = core
    my_info%threads_info(1)%node = ma_get_proc_node()
    my_info%threads_info(1)%id_omp = id
    my_info%threads_info(1)%id_real = ma_get_id()
    my_info%mp_info%myid = my_info%threads_info(1)%id_real
    my_info%core = core
    my_info%node = my_info%threads_info(1)%node
   END IF
   !only if it is the Process/Thread - rank 0
   IF (( id .EQ. 0 ).AND. ( my_info%mp_info%myproc .EQ. 0 )) THEN
      CALL ma_set_net_affinity(sched_unit,id)
   ENDIF
#endif
#if defined (__LIBNUMA) && !defined (__HWLOC)
    IF (sched_unit .EQ. threads) THEN
      core = MOD(my_info%mp_info%myproc*my_info%nr_threads+id,ma_linux_ncores())
      CALL ma_linux_set_core(core)
      my_info%threads_info(id+1)%core = core
      my_info%threads_info(id+1)%node = ma_linux_get_mynode()
      my_info%threads_info(id+1)%id_omp = id
      my_info%threads_info(id+1)%id_real = ma_linux_get_thread_id()
    ELSE
     IF(my_info%nr_threads .GT. 1) THEN
       core = MOD(my_info%mp_info%myproc*my_info%nr_threads,ma_get_ncores())
     ELSE
       core = (MOD(my_info%mp_info%myproc,ma_get_ncores()))
    ENDIF
    CALL ma_linux_set_core(core)
    my_info%threads_info(1)%core = core
    my_info%threads_info(1)%node = ma_linux_get_mynode()
    my_info%threads_info(1)%id_omp = id
    my_info%threads_info(1)%id_real = ma_get_id()
    my_info%mp_info%myid = my_info%threads_info(1)%id_real
    my_info%core = core
    my_info%node = my_info%threads_info(1)%node
   ENDIF
   !only if it is the Process/Thread - rank 0
   IF (( id .EQ. 0 ).AND. ( my_info%mp_info%myproc .EQ. 0 )) THEN
        CALL ma_set_net_affinity(sched_unit,id)
   ENDIF
#endif
     END SUBROUTINE ma_linear_place

! *****************************************************************************
!> \brief Set the MPI/threads mapping on the machine cores - scatter
!> \param id thread OMP id - MPI rank 
!> \param sched_unit ...
!> \note The strategy used hier is a simple compact distrbution for threads of a
!> \note process. In this case, threads of a process take into account the process
!> \note affinity.
!> \note For MPI only process, the strategy is a scatter. In this case, the scheduling
!> \note provides as much as possible cache memory for MPI processes.
! *****************************************************************************

  SUBROUTINE ma_scatter_place(id, sched_unit)
    INTEGER                                  :: id, sched_unit

    INTEGER                                  :: core, distance

    core = 0
    distance = 0
    isdefault = .FALSE.

#if defined (__HWLOC) && !defined (__LIBNUMA)
 IF (sched_unit .EQ. threads) THEN
    core = my_info%core+id
    CALL ma_set_core(core)
    my_info%threads_info(id+1)%core = core
    my_info%threads_info(id+1)%node = ma_get_mynode()
    my_info%threads_info(id+1)%id_omp = id
    my_info%threads_info(id+1)%id_real = ma_get_thread_id()
   ELSE
    distance = topology%ncores/my_info%nr_neighbors
    core = (MOD(my_info%mp_info%myproc,ma_get_ncores()))*distance
    CALL ma_set_proc_core(core)
    my_info%threads_info(1)%core = core
    my_info%threads_info(1)%node = ma_get_proc_node()
    my_info%threads_info(1)%id_omp = id
    my_info%threads_info(1)%id_real = ma_get_id()
    my_info%mp_info%myid = my_info%threads_info(1)%id_real
    my_info%core = core
    my_info%node = my_info%threads_info(1)%node
   END IF
   !only if it is the Process/Thread - rank 0
   IF (( id .EQ. 0 ).AND. ( my_info%mp_info%myproc .EQ. 0 )) THEN
         CALL ma_set_net_affinity(sched_unit,id)
   ENDIF
#endif
#if defined (__LIBNUMA) && !defined (__HWLOC)
    IF (sched_unit .EQ. threads) THEN
      core = my_info%core+id
      CALL ma_linux_set_core(core)
      my_info%threads_info(id+1)%core = core
      my_info%threads_info(id+1)%node = ma_linux_get_mynode()
      my_info%threads_info(id+1)%id_omp = id
      my_info%threads_info(id+1)%id_real = ma_linux_get_thread_id()
    ELSE
       distance = topology%ncores/my_info%nr_neighbors
       core = (MOD(my_info%mp_info%myproc,ma_get_ncores()))*distance
    CALL ma_linux_set_core(core)
    my_info%threads_info(1)%core = core
    my_info%threads_info(1)%node = ma_linux_get_mynode()
    my_info%threads_info(1)%id_omp = id
    my_info%threads_info(1)%id_real = ma_get_id()
    my_info%mp_info%myid = my_info%threads_info(1)%id_real
    my_info%core = core
    my_info%node = my_info%threads_info(1)%node
   ENDIF
   !only if it is the Process/Thread - rank 0
   IF (( id .EQ. 0 ).AND. ( my_info%mp_info%myproc .EQ. 0 )) THEN
         CALL ma_set_net_affinity(sched_unit,id)
   ENDIF
#endif
     END SUBROUTINE ma_scatter_place

! *****************************************************************************
!> \brief Set the MPI/threads mapping on the machine cores - group
!> \param id thread OMP id - MPI rank
!> \param sched_unit ...
!> \note Group the MPI processes in a Hierarchical way following the machine topology
! *****************************************************************************

  SUBROUTINE ma_group_place(id, sched_unit)
    INTEGER                                  :: id, sched_unit

    INTEGER                                  :: core, distance, new_id, node

    core = 0
    distance = 0
    new_id = 0
    node = 0
    isdefault = .FALSE.

#if defined (__HWLOC) && !defined (__LIBNUMA)
 IF (sched_unit .EQ. threads) THEN
   IF (my_info%nr_threads .GT. (topology%ncores/topology%nnodes)) THEN
    CALL ma_set_thread_allnodes()
   ELSE
   ! threads inherit their parent CPU mask
   ! just get where they are running
    CALL ma_set_thread_node(my_info%node)
   ENDIF
    core = ma_get_mycore()
    my_info%threads_info(id+1)%core = core
    my_info%threads_info(id+1)%node = ma_get_core_node(core)
    my_info%threads_info(id+1)%id_omp = id
    my_info%threads_info(id+1)%id_real = ma_get_thread_id()
 ELSE
    !NUMA machine - integrated memory controllers
    IF ( topology%nnodes .GT. 1 ) THEN
     IF (my_info%nr_neighbors .LE. 1) THEN
       node = 0
     ELSE
       new_id = MOD(my_info%mp_info%myproc,my_info%nr_neighbors)
       IF ( my_info%nr_neighbors .LE. topology%nnodes) THEN
         node = new_id
       ELSE
         node = MOD(new_id/(my_info%nr_neighbors/topology%nnodes),topology%nnodes)
       ENDIF
     ENDIF
     CALL ma_set_proc_node(node)
    ELSE ! UMA machine
     distance = topology%ncores/my_info%nr_neighbors
     core = (MOD(my_info%mp_info%myproc,ma_get_ncores()))*distance
     node = 0
     CALL ma_set_thread_cores(distance, core)
    END IF
    my_info%threads_info(1)%core = ma_get_proc_core(1)
    my_info%threads_info(1)%node = node
    my_info%threads_info(1)%id_omp = 0 ! MPI process is always the master thread
    my_info%threads_info(1)%id_real = ma_get_id()
    my_info%mp_info%myid = my_info%threads_info(1)%id_real
    my_info%core = my_info%threads_info(1)%core
    my_info%node = node
 END IF
  !only if it is the Process/Thread - rank 0
  IF (( id .EQ. 0 ).AND. ( my_info%mp_info%myproc .EQ. 0 )) THEN
         CALL ma_set_net_affinity(sched_unit,id)
  ENDIF
#endif
#if defined (__LIBNUMA) && !defined (__HWLOC)
 IF (sched_unit .EQ. threads) THEN
   ! threads inherit their parent CPU mask
   ! just get where they are running
   CALL ma_linux_set_procnode(my_info%node)
    my_info%threads_info(id+1)%core = ma_linux_get_core()
    my_info%threads_info(id+1)%node = ma_linux_get_mynode()
    my_info%threads_info(id+1)%id_omp = id
    my_info%threads_info(id+1)%id_real = ma_linux_get_thread_id()
 ELSE
   !NUMA machine - integrated memory controllers
    IF ( topology%nnodes .GT. 1 ) THEN
     IF (my_info%nr_neighbors .LE. 1) THEN
       node = 0
     ELSE
       new_id = MOD(my_info%mp_info%myproc,my_info%nr_neighbors)
       IF ( my_info%nr_neighbors .LE. topology%nnodes) THEN
        node = new_id
       ELSE
         node = MOD(new_id/(my_info%nr_neighbors/topology%nnodes),topology%nnodes)
       ENDIF
     ENDIF
     CALL ma_set_proc_node(node)
    ELSE ! UMA machine
     distance = topology%ncores/my_info%nr_neighbors
     core = (MOD(my_info%mp_info%myproc,ma_get_ncores()))*distance
     node = 0
     CALL ma_set_proc_cores(distance,core)
    END IF
    my_info%threads_info(1)%core = ma_linux_get_core()
    my_info%threads_info(1)%node = node
    my_info%threads_info(1)%id_omp = 0
    my_info%threads_info(1)%id_real = ma_get_id()
    my_info%mp_info%myid = my_info%threads_info(1)%id_real
    my_info%core = my_info%threads_info(1)%core
    my_info%node = node
 END IF
 !only if it is the Process/Thread - rank 0
   IF (( id .EQ. 0 ).AND. ( my_info%mp_info%myproc .EQ. 0 )) THEN
     CALL ma_set_net_affinity(sched_unit,id)
   ENDIF
#endif
   END SUBROUTINE ma_group_place

! *****************************************************************************
!> \brief Set the MPI mapping on the machine cores - numa node
!> \param node ...
!> \param id thread OMP id - MPI rank
!> \param sched_unit ...
!> \note Group the MPI processes in a NUMA way following the machine topology
! *****************************************************************************
  SUBROUTINE ma_numa_node_place(node,id,sched_unit)
    INTEGER                                  :: node, id, sched_unit

    INTEGER                                  :: core

    core = 0
#if defined (__HWLOC) && !defined (__LIBNUMA)
 IF (sched_unit .EQ. threads) THEN
    core = MOD(my_info%mp_info%myproc*my_info%nr_threads+id,ma_get_ncores())
    CALL ma_set_core(core)
    my_info%threads_info(id+1)%core = ma_get_mycore()
    my_info%threads_info(id+1)%node = ma_get_mynode()
    my_info%threads_info(id+1)%id_omp = id
    my_info%threads_info(id+1)%id_real = ma_get_thread_id()
   ELSE
    CALL ma_set_proc_node(node)
    my_info%threads_info(1)%core = ma_get_proc_core(0)
    my_info%threads_info(1)%node = node
    my_info%threads_info(1)%id_omp = id
    my_info%threads_info(1)%id_real = ma_get_id()
    my_info%mp_info%myid = my_info%threads_info(1)%id_real
    my_info%core = my_info%threads_info(1)%core
    my_info%node = node
   END IF
   !only if it is the Process/Thread - rank 0
   IF (( id .EQ. 0 ).AND. ( my_info%mp_info%myproc .EQ. 0 )) THEN
         CALL ma_set_net_affinity(sched_unit,id)
   ENDIF
#endif
#if defined (__LIBNUMA) && !defined (__HWLOC)
 IF (sched_unit .EQ. threads) THEN
    core = MOD(my_info%mp_info%myproc*my_info%nr_threads+id,ma_linux_ncores())
    CALL ma_linux_set_core(core)
    my_info%threads_info(id+1)%core = core
    my_info%threads_info(id+1)%node = ma_linux_get_mynode()
    my_info%threads_info(id+1)%id_omp = id
    my_info%threads_info(id+1)%id_real = ma_linux_get_thread_id()
   ELSE
    CALL ma_set_proc_node(node)
    my_info%threads_info(1)%core = ma_linux_get_core()
    my_info%threads_info(1)%node = node
    my_info%threads_info(1)%id_omp = id
    my_info%threads_info(1)%id_real = ma_get_id()
    my_info%mp_info%myid = my_info%threads_info(1)%id_real
    my_info%core = my_info%threads_info(1)%core
    my_info%node = node
   END IF
   !only if it is the Process/Thread - rank 0
   IF (( id .EQ. 0 ).AND. ( my_info%mp_info%myproc .EQ. 0 )) THEN
         CALL ma_set_net_affinity(sched_unit,id)
   ENDIF
#endif

  END SUBROUTINE ma_numa_node_place

! *****************************************************************************
!> \brief Set the MPI/threads mapping on the machine cores - manual
!> \param id thread OMP id - MPI rank 
!> \param sched_unit ...
!> \note Group the MPI processes following the user inputs
! *****************************************************************************
  SUBROUTINE ma_manual_place(id,sched_unit)
    INTEGER                                  :: id, sched_unit

    INTEGER                                  :: core, node

    core = 0
    node = 0
    isdefault = .FALSE.

#if defined (__HWLOC) && !defined (__LIBNUMA)
 IF (sched_unit .EQ. threads) THEN
 ! offset to the parent MPI process
    core = my_info%core+thread(id+1)
    CALL ma_set_core(core)
    my_info%threads_info(id+1)%core = core
    my_info%threads_info(id+1)%node = ma_get_mynode()
    my_info%threads_info(id+1)%id_omp = id
    my_info%threads_info(id+1)%id_real = ma_get_thread_id()
   ELSE
    core = process(MOD(my_info%mp_info%myproc,my_info%nr_neighbors)+1)
    CALL ma_set_proc_core( core)
    node = ma_get_proc_node ()
    my_info%threads_info(1)%core = core
    my_info%threads_info(1)%node = node
    my_info%threads_info(1)%id_omp = id
    my_info%threads_info(1)%id_real = ma_get_id()
    my_info%mp_info%myid = my_info%threads_info(1)%id_real
    my_info%core = core
    my_info%node = node
   END IF
#endif
#if defined (__LIBNUMA) && !defined (__HWLOC)
 IF (sched_unit .EQ. threads) THEN
 !offset to the parent MPI process
    core = my_info%core+thread(id+1)
    CALL ma_linux_set_core(core)
    my_info%threads_info(id+1)%core = core
    my_info%threads_info(id+1)%node = ma_linux_get_mynode()
    my_info%threads_info(id+1)%id_omp = id
    my_info%threads_info(id+1)%id_real = ma_linux_get_thread_id()
   ELSE
    core = process(MOD(my_info%mp_info%myproc,my_info%nr_neighbors)+1)
    CALL ma_linux_set_core(core)
    node = ma_get_proc_node ()
    my_info%threads_info(1)%core = core
    my_info%threads_info(1)%node = node
    my_info%threads_info(1)%id_omp = id
    my_info%threads_info(1)%id_real = ma_get_id()
    my_info%mp_info%myid = my_info%threads_info(1)%id_real
    my_info%core = core
    my_info%node = node
   END IF
#endif
   END SUBROUTINE ma_manual_place

! *****************************************************************************
!> \brief Get the core and NUMA node of the process
!> \param my_core the core where the process is running
!> \param my_node the NUMA node where the process is running 
! *****************************************************************************
  SUBROUTINE ma_running_on (my_core, my_node)
      INTEGER                                  :: my_core, my_node
#if defined (__HWLOC) || defined (__LIBNUMA)
   IF ( ma_get_conf_mpiSched() .EQ. 'N' .OR. &
        ma_get_conf_mpiSched() .EQ. 'n' .OR. &
        ma_get_conf_mpiSched() .EQ. 'D' ) THEN
      my_core = ma_get_proc_core(1)
#if defined (__HWLOC)
    my_node = ma_get_core_node(my_core)
#endif
#if defined (__LIBNUMA)
      my_node = ma_linux_get_core_node(my_core)
#endif
      my_info%core = my_core
      my_info%node = my_node
    ELSE
     my_core = my_info%core
     my_node = my_info%node
    ENDIF
#endif
  END SUBROUTINE ma_running_on


! *****************************************************************************
!> \brief Update for all processes where the threads are running, also print
!> \param unit_num the output console 
! *****************************************************************************
  SUBROUTINE ma_thread_running_on(unit_num)
    INTEGER                                  :: unit_num

#if defined (__HWLOC) || defined (__LIBNUMA)
    INTEGER                                  :: id, istat, p, core
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: thread_send
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: thread_recv
    CHARACTER(len=4)                         :: first_core, second_core
    CHARACTER(LEN=default_string_length)     :: host_name
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: all_host
    CHARACTER(len=10)                        :: all_cores

  ! If not all MPI are present
  IF (.NOT. my_info%mp_info%all_proc) THEN
    CALL ma_show_threads(unit_num)
  ELSE
   IF (has_MPI) THEN

    ALLOCATE(thread_send(4,my_info%nr_threads))
    ALLOCATE(thread_recv(4,my_info%nr_threads,my_info%mp_info%numproc))

    thread_send(:,:) = 0
    thread_recv(:,:,:) = 0

    CALL ma_pack_threads(my_info,thread_send)
    CALL mp_allgather(thread_send, thread_recv, my_info%mp_info%mp_group)
    CALL ma_unpack_threads(my_info,thread_recv,thread_mapping)

    ALLOCATE (all_host(default_string_length,my_info%mp_info%numproc),STAT=istat)
    IF (istat /= 0) CALL ma_error_stop(ma_error_allocation)

    all_host(:,:) = 0
    CALL mp_proc_name(host_name)
    CALL string_to_ascii(host_name,all_host(:,my_info%mp_info%myproc+1))
    CALL mp_sum(all_host,my_info%mp_info%mp_group)
   ELSE
     host_name = "localhost"
   ENDIF

   IF(unit_num .GT. 0) THEN
   IF (my_info%nr_threads .GT. 1) THEN
     WRITE(unit_num,"(T2,A)") "OMP | Thread placement"
     DO p=1, my_info%mp_info%numproc
      IF (has_MPI) THEN
        CALL ascii_to_string(all_host(:,p),host_name)
      END IF
      WRITE(unit_num, "(T2,A,I4,A,A)") "Processs ", p-1, " running on ", TRIM(host_name)
      DO id=1, my_info%nr_threads
           core = thread_mapping(p,id)%core
           CALL integer_to_string(core,first_core)
           all_cores = first_core
           IF (thread_mapping(p,id)%node .EQ. -1) THEN
             WRITE(unit_num, "(T2,A,I4,2X,A,A)")   "OMP | Thread ",&
                   thread_mapping(p,id)%id_omp,&
                   " running on core ", all_cores
           ELSE
              WRITE(unit_num, "(T2,A,I4,2X,A,I4,2X,A,A)")   "OMP | Thread ",&
              thread_mapping(p,id)%id_omp," running on NUMA node ",&
              thread_mapping(p,id)%node, " core ", all_cores
           ENDIF
      END DO
     END DO
     WRITE(unit_num,'()')
   ENDIF
   ENDIF
  ENDIF
#endif

  END SUBROUTINE ma_thread_running_on

! *****************************************************************************
!> \brief Shows threads of each MPI process
!> \param unit_num ...
! *****************************************************************************
  SUBROUTINE ma_show_threads (unit_num)
    INTEGER                                  :: unit_num

    CHARACTER(len=10)                        :: all_cores
    CHARACTER(len=4)                         :: first_core
    INTEGER                                  :: core, id, p

  IF ( print_thread ) THEN
   IF(unit_num .GT. 0) THEN
   IF (my_info%nr_threads .GT. 1) THEN
     WRITE(unit_num,"()")
     WRITE(unit_num,"(T2,A)") "OMP | Thread placement"
      p=my_info%mp_info%numproc+1
      WRITE(unit_num, "(T2,A,I4)") "Processs ", p-1
      DO id=1, my_info%nr_threads
           core = thread_mapping(p,id)%core
           CALL integer_to_string(core,first_core)
           all_cores = first_core
           IF (thread_mapping(p,id)%node .EQ. -1) THEN
             WRITE(unit_num, "(T2,A,I4,2X,A,A)")   "OMP | Thread ",&
                   thread_mapping(p,id)%id_omp,&
                   " running on core ", all_cores
           ELSE
              WRITE(unit_num, "(T2,A,I4,2X,A,I4,2X,A,A)")   "OMP | Thread ",&
              thread_mapping(p,id)%id_omp," running on NUMA node ",&
              thread_mapping(p,id)%node, " core ", all_cores
       ENDIF
     END DO
     WRITE(unit_num,'()')
   ENDIF
   ENDIF
  ENDIF
  END SUBROUTINE ma_show_threads

! *****************************************************************************
!> \brief Print the current core and NUMA node where a MPI process is running
!> \param id the MPI process id (use cp2k_info, r_pid) 
!> \param unit_num the output console 
! *****************************************************************************
  SUBROUTINE ma_current_mpi_run(id,unit_num)
    INTEGER                                  :: id, unit_num

    INTEGER                                  :: core, node

     core = ma_get_proc_core(1)
     node = ma_get_proc_node()
     unit_num = default_output_unit
    IF(unit_num .GT. 0) THEN
      IF (node .EQ. -1 .AND. core .EQ. -1) THEN
        WRITE(unit_num, "(T2,A)") "WARNING: Compile with HWLOC"
      ELSE IF (node .EQ. -1) THEN
         WRITE(unit_num, "(T2,A,I6,A,I4,2X,A,I3)") "MPI | Process ", id,&
         " saved id ", my_info%mp_info%myproc, " running on core ", core
      ELSE
        WRITE(unit_num, "(T2,A,I6,A,I6,2X,A,I4,2X,A,I4)")   "MPI | Process ", &
             id, " saved id ", my_info%mp_info%myproc, &
             " running on NUMA node ", node, " core ", core
     ENDIF
   ENDIF
  END SUBROUTINE ma_current_mpi_run


! *****************************************************************************
!> \brief Print the current core and NUMA node where a thread is running
!> \param id Thread OMP id 
!> \param unit_num the output console 
! *****************************************************************************
  SUBROUTINE ma_current_thread_run(id, unit_num)
    INTEGER                                  :: id, unit_num

    INTEGER                                  :: core, node

     core = ma_get_mycore()
     node = ma_get_mynode()

   IF(unit_num .GT. 0) THEN
     IF (node .EQ. -1) THEN
         WRITE(unit_num, "(T2,A,I4,A,I4,2X,A,I3)") "OMP | Thread ", id,&
         " from process ", my_info%mp_info%myproc, " running on core ", core
     ELSE
        WRITE(unit_num, "(T2,A,I4,A,I4,2X,A,I4,2X,A,I4)")   "OMP | Thread ", id,&
        " from process ", my_info%mp_info%myproc, " running on NUMA node ", node, " core ", core
     ENDIF
   ENDIF
  END SUBROUTINE ma_current_thread_run

! *****************************************************************************
! Get the core for a thread
! *****************************************************************************
!> \brief ...
!> \param id ...
! *****************************************************************************
  SUBROUTINE ma_get_thread_run(id)
    INTEGER  :: id
#if defined (__HWLOC) || defined (__LIBNUMA)
  my_info%threads_info(id+1)%core = ma_get_mycore()
  my_info%threads_info(id+1)%node = ma_get_mynode()
  my_info%threads_info(id+1)%id_omp = id
#endif
#if defined (__HWLOC) && !defined (__LIBNUMA)
  my_info%threads_info(id+1)%id_real = ma_get_thread_id()
#endif
#if !defined (__HWLOC) && defined (__LIBNUMA)
  my_info%threads_info(id+1)%id_real = ma_linux_get_thread_id()
#endif
  END SUBROUTINE ma_get_thread_run


! *****************************************************************************
! Get/set the core for a thread
! *****************************************************************************
!> \brief ...
!> \param id ...
!> \param policy ...
! *****************************************************************************
  SUBROUTINE ma_sched_threads(id, policy)
    INTEGER                                  :: id, policy

     SELECT CASE (policy)
     CASE (nosched)
       CALL ma_no_place(id, threads)
     CASE (linear)
       CALL ma_linear_place(id, threads)
     CASE (scatter)
       CALL ma_scatter_place(id, threads)
     CASE (group)
       CALL ma_group_place(id, threads)
     CASE (manual)
       CALL ma_manual_place(id, threads)
     CASE (def)
       CALL ma_def_place(id, threads)
     END SELECT
  END SUBROUTINE ma_sched_threads

! *****************************************************************************
! Get/set the core for a mpi
! *****************************************************************************
!> \brief ...
!> \param policy ...
! *****************************************************************************
  SUBROUTINE ma_sched_mpi(policy)
    INTEGER                                  :: policy

     SELECT CASE (policy)
     CASE (nosched)
       CALL ma_no_place(my_info%mp_info%myproc, mpi)
     CASE (linear)
       CALL ma_linear_place(my_info%mp_info%myproc, mpi)
     CASE (scatter)
       CALL ma_scatter_place(my_info%mp_info%myproc, mpi)
     CASE (group)
       CALL ma_group_place(my_info%mp_info%myproc, mpi)
     CASE (manual)
       CALL ma_manual_place(my_info%mp_info%myproc, mpi)
     CASE (def)
       CALL ma_def_place(my_info%mp_info%myproc, mpi)
     END SELECT
  END SUBROUTINE ma_sched_mpi

! *****************************************************************************
!> \brief Verify the last and current core for a thread
!> \param id Thread OMP id 
!> \param unit_num the output console 
! *****************************************************************************
  SUBROUTINE ma_verify_place(id, unit_num)
    INTEGER                                  :: id, unit_num

    INTEGER                                  :: core, node, pe

     core = ma_get_mycore()
     node = ma_get_mynode()
     pe = my_info%mp_info%myproc

     IF (thread_mapping(pe+1,id+1)%core .NE. core) THEN
      WRITE(unit_num,'(T2,A)') "SCHED | Thread placement changed"
      IF (node .EQ. -1) THEN
         WRITE(unit_num, "(T2,A,I4,2X,A,I3)")   "OMP | Thread ", thread_mapping(pe+1,id+1)%id_omp,&
         " WAS running on core ", thread_mapping(pe+1,id+1)%core
         WRITE(unit_num, "(T2,A,I4,2X,A,I3)")   "OMP | Thread ", id,&
         " running on core ", core
      ELSE
        WRITE(unit_num, "(T2,A,I4,2X,A,I4,2X,A,I4)")   "OMP | Thread ", &
        thread_mapping(pe+1,id+1)%id_omp," WAS running on NUMA node ",&
        thread_mapping(pe+1,id+1)%node, " core ",thread_mapping(pe+1,id+1)%core
        WRITE(unit_num, "(T2,A,I4,2X,A,I4,2X,A,I4)")   "OMP | Thread ", id,&
        " running on NUMA node ", node, " core ", core
      ENDIF
     ENDIF
  END SUBROUTINE ma_verify_place

! *****************************************************************************
!> \brief Map the MPI process 0 to the NUMA node that holds the network card(s)
!> \param proc ...
!> \param id ...
! *****************************************************************************
  SUBROUTINE ma_set_net_affinity (proc, id)
    INTEGER(KIND=C_INT)                      :: proc, id

    INTEGER(KIND=C_INT)                      :: node

   node = my_info%node

#if defined (__HWLOC)
 IF(id .EQ. 0) THEN
   IF ( ma_get_nnetDev(node) .GT. 0 ) THEN
     CALL ma_set_first_core(node, proc)
   ELSE
     node = ma_get_node_netDev()
     CALL ma_set_first_core(node, proc)
     my_info%node = node
   ENDIF
 ENDIF
#elif defined (__LINUMA)
 IF(id .EQ. 0) THEN
   !usually node 0 and core 0 are the ones connect to I/O
   my_info%node = 0
   my_info%core = 0
   CALL ma_linux_set_core(my_info%core)
 ENDIF
#endif

 END SUBROUTINE ma_set_net_affinity

! *****************************************************************************
!> \brief Get the Network card affinity for the rank 0 process
!> \param cardId ...
! *****************************************************************************
  SUBROUTINE ma_get_net_affinity (cardId)
    INTEGER                                  :: cardId

!no net card affinity

    cardId = -1

#if defined (__HWLOC)
     cardId = ma_get_netDev(my_info%core, my_info%mp_info%myproc)
#endif
   END SUBROUTINE ma_get_net_affinity

! *****************************************************************************
!> \brief Set the GPU for a MPI process
!> \param myRank the MPI rank of the process 
!> \retval gpuId ...
! *****************************************************************************
 INTEGER FUNCTION ma_set_gpu_affinity (myRank) RESULT (gpuId)
   INTEGER(KIND=C_INT) ::  myRank, rank,core

   rank = my_info%mp_info%myproc
   core = my_info%core

   gpuId = 0

  IF (myRank .EQ. rank) THEN
#if defined (__DBCSR_ACC) || defined (__PW_CUDA)
#if defined (__HWLOC) || defined (__LIBNUMA)
   gpuId = ma_get_cuDev(core,myRank,ma_get_neighbors())
   my_info%gpu = gpuId
#endif
#endif
  ENDIF
 END FUNCTION ma_set_gpu_affinity

! *****************************************************************************
!> \brief Get the gpu affinity of a process
!> \param gpuId ...
! *****************************************************************************
  SUBROUTINE ma_get_gpu_affinity (gpuId)
    INTEGER                                  :: gpuId

    gpuId = -1

#if defined (__DBCSR_ACC) || defined (__PW_CUDA)
#if defined (__HWLOC) || defined (__LIBNUMA)
     gpuId = my_info%gpu
#endif
#endif
   END SUBROUTINE ma_get_gpu_affinity

! *****************************************************************************
!> \brief Print the affinity for all MPI processes
!> \param print_proc ...
!> \param print_mem ...
!> \param unit_nr ...
! *****************************************************************************
   SUBROUTINE ma_print_proc_affinity (print_proc,print_mem,unit_nr)
    LOGICAL                                  :: print_proc, print_mem
    INTEGER                                  :: unit_nr

#if defined (__HWLOC) || defined (__LIBNUMA)
    CHARACTER(len=1)                         :: mem_pol
    CHARACTER(len=10)                        :: all_cores, mempolicy
    CHARACTER(len=4)                         :: first_core, second_core
    INTEGER                                  :: core, ipe, istat
    INTEGER(KIND=int_4), ALLOCATABLE, &
      DIMENSION(:)                           :: all_core, &
                                                all_memory, all_node, &
                                                all_node_mem, all_pid
#if defined (__DBCSR_ACC) || defined (__PW_CUDA)
    INTEGER(KIND=int_4), ALLOCATABLE, &
      DIMENSION(:)                           :: all_gpu
#endif
#endif

#if defined (__HWLOC) || defined (__LIBNUMA)

    IF ((.NOT. my_info%mp_info%all_proc) .AND. print_proc) THEN
      CALL ma_current_mpi_run(my_info%mp_info%myid,unit_nr)
    ELSE
       ALLOCATE (all_core(my_info%mp_info%numproc),STAT=istat)
       IF (istat /= 0) THEN
         CALL ma_error_stop(ma_error_allocation)
       END IF
       all_core(:) = 0

       ALLOCATE (all_node(my_info%mp_info%numproc),STAT=istat)
       IF (istat /= 0) THEN
         CALL ma_error_stop(ma_error_allocation)
       END IF
       all_node(:) = 0

#if defined (__DBCSR_ACC) || defined (__PW_CUDA)

       ALLOCATE (all_gpu(my_info%mp_info%numproc),STAT=istat)
       IF (istat /= 0) THEN
         CALL ma_error_stop(ma_error_allocation)
       END IF
       all_gpu(:) = 0
       CALL ma_get_gpu_affinity(all_gpu(my_info%mp_info%myproc+1))
#endif

       ALLOCATE (all_memory(my_info%mp_info%numproc),STAT=istat)
       IF (istat /= 0) THEN
         CALL ma_error_stop(ma_error_allocation)
       END IF
       all_memory(:) = 0

       ALLOCATE (all_node_mem(my_info%mp_info%numproc),STAT=istat)
       IF (istat /= 0) THEN
         CALL ma_error_stop(ma_error_allocation)
       END IF
       all_node_mem(:) = 0

       !Where process is running

       CALL ma_running_on(all_core(my_info%mp_info%myproc+1), all_node(my_info%mp_info%myproc+1))

       IF (ma_get_conf_mpiSched() .EQ. 'G' .OR. &
           ma_get_conf_mpiSched() .EQ. 'g' .OR. &
           isdefault .AND. my_info%mp_info%numproc .NE. my_info%nr_neighbors) THEN
               all_core(my_info%mp_info%myproc+1)=ma_my_first_core()
       END IF

   IF (print_proc) THEN
 ! ***** Print where process are running!
#if defined (__DBCSR_ACC) || defined (__PW_CUDA)
       CALL mp_sum(all_gpu,my_info%mp_info%mp_group)
#endif

       CALL mp_sum(all_core,my_info%mp_info%mp_group)
       CALL mp_sum(all_node,my_info%mp_info%mp_group)

       ALLOCATE (all_pid(my_info%mp_info%numproc),STAT=istat)
       IF (istat /= 0) THEN
         CALL ma_error_stop(ma_error_allocation)
       END IF
       all_pid(:) = 0
       all_pid(my_info%mp_info%myproc+1) = my_info%mp_info%myid
       CALL mp_sum(all_pid,my_info%mp_info%mp_group)

      IF ( ma_get_conf_mempol() .NE. 'N' ) THEN
       CALL ma_get_mempol(all_node_mem(my_info%mp_info%myproc+1),all_memory(my_info%mp_info%myproc+1))
       CALL mp_sum(all_memory,my_info%mp_info%mp_group)
       CALL mp_sum(all_node_mem,my_info%mp_info%mp_group)
      ENDIF

       IF (unit_nr>0) THEN

       ! Inform the user about the network affinity
        IF (hasnet) THEN
           WRITE (unit_nr,"(T2,A)") "NETWORK| Affinity is on"
           WRITE (unit_nr,"()")
        ENDIF

        IF ( ma_get_conf_mpiSched() .EQ. 'N' .OR. &
             ma_get_conf_mpiSched() .EQ. 'n' .OR. &
             ma_get_conf_mpiSched() .EQ. 'D' ) THEN
           WRITE (unit_nr,"(T2,A)") "SCHED| Processes are now runing"
        ELSE
           WRITE (unit_nr,"(T2,A)") "SCHED| Processes mapping"
        ENDIF

         DO ipe=1,my_info%mp_info%numproc
           IF (ma_get_conf_mpiSched() .EQ. 'G' .OR. &
               ma_get_conf_mpiSched() .EQ. 'g' .OR. &
               isdefault .AND. my_info%mp_info%numproc .NE. my_info%nr_neighbors) THEN
             core = all_core(ipe)
             CALL integer_to_string(core,first_core)
             core = core + (topology%ncores/my_info%nr_neighbors)-1
             CALL integer_to_string(core,second_core)
             all_cores = TRIM(first_core)//"-"//TRIM(second_core)
           ELSE
             core = all_core(ipe)
             CALL integer_to_string(core,first_core)
             all_cores = first_core
           ENDIF
#if defined (__DBCSR_ACC) || defined (__PW_CUDA)
           IF (all_node(ipe) .EQ. -1) THEN
             WRITE(unit_nr,"(T2,A,I8,A,I8,5X,20X,A,A,A,I3)")&
                           "Process ", ipe-1, " ", all_pid(ipe),&
                           " on core ", TRIM(all_cores), &
                           "  attched to gpu ", all_gpu(ipe)
             ELSE
              WRITE (unit_nr,"(T2,A,I8,A,I8,5X,A,I4,2X,A,A,A,I3)")&
                             "Process ", ipe-1, " ", all_pid(ipe), &
                             " on NUMA node ",all_node(ipe),&
                             " core ",TRIM(all_cores), &
                             "  attached to gpu ",all_gpu(ipe)
          ENDIF

#else
          IF (all_node(ipe) .EQ. -1) THEN
             WRITE(unit_nr,"(T2,A,I8,A,I8,5X,30X,A,A)")&
               "Process ", ipe-1, " ", all_pid(ipe), " runing on core ",&
                TRIM(all_cores)
           ELSE
             WRITE (unit_nr,"(T2,A,I8,A,I8,5X,A,I6,7X,A,A)")&
               "Process ", ipe-1, " ", all_pid(ipe)," running on NUMA node ",&
                all_node(ipe), " core ", TRIM(all_cores)
           ENDIF
#endif
         END DO

         WRITE (unit_nr,"(T2,A)") ""
         IF (ma_get_nnodes() .GT. 0 .AND. ma_get_conf_mempol() .NE. 'N') THEN
            WRITE (unit_nr,"(T2,A)") "MEMORY| Processes memory mapping"
            mem_pol = ma_get_conf_mempol()

            DO ipe=1,my_info%mp_info%numproc
             IF (mem_pol .EQ. 'D') THEN
               WRITE (unit_nr,"(T2,A,I8,A,I8,5X,A)") &
                     "Process ", ipe-1, " ", all_pid(ipe)," OS memory policy"
             ELSE
               CALL ma_string_mempol(all_memory(ipe), mempolicy)
              IF (all_memory(ipe) .EQ. -1) THEN
                WRITE (unit_nr,"(T2,A,I8,A,I8,5X,A)") &
                "Process ", ipe-1, " ", all_pid(ipe)," memory policy not applied"
              ELSE IF (all_node_mem(ipe) .EQ. -1) THEN
                  WRITE (unit_nr,"(T2,A,I8,A,I8,5X,A,A,A)")&
                 "Process ", ipe-1, " ", all_pid(ipe)," memory policy ",&
                  TRIM(mempolicy), " all nodes"
               ELSE
                 WRITE (unit_nr,"(T2,A,I8,A,I8,5X,A,A,7X,A,I8)")&
                 "Process ", ipe-1, " ", all_pid(ipe)," memory policy ",&
                  TRIM(mempolicy), " node ", all_node_mem(ipe)
               ENDIF
             ENDIF
            END DO
            WRITE (unit_nr,"(T2,A)") ""
         ENDIF
       END IF

       DEALLOCATE (all_memory,STAT=istat)
       IF (istat /= 0) THEN
         CALL ma_error_stop(ma_error_allocation)
       END IF
       DEALLOCATE (all_node_mem,STAT=istat)
       IF (istat /= 0) THEN
         CALL ma_error_stop(ma_error_allocation)
       END IF
       DEALLOCATE (all_pid,STAT=istat)
       IF (istat /= 0) THEN
         CALL ma_error_stop(ma_error_allocation)
       END IF

   ENDIF

       DEALLOCATE (all_core,STAT=istat)
       IF (istat /= 0) THEN
         CALL ma_error_stop(ma_error_allocation)
       END IF
       DEALLOCATE (all_node,STAT=istat)
       IF (istat /= 0) THEN
         CALL ma_error_stop(ma_error_allocation)
       END IF
#if defined (__DBCSR_ACC) || defined (__PW_CUDA)
       DEALLOCATE (all_gpu,STAT=istat)
       IF (istat /= 0) THEN
         CALL ma_error_stop(ma_error_allocation)
       END IF
#endif
   ENDIF

#endif
   END SUBROUTINE ma_print_proc_affinity


END MODULE ma_affinity
