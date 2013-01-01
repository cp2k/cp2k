!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group   !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Parameter stack and queue methods.
!> \author  Urban Borstnik
!> \date    2011-06-17
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2011-06-17
! *****************************************************************************
MODULE dbcsr_pq_methods
  USE dbcsr_config,                    ONLY: has_omp_2_5,&
                                             mm_thread_workshare
  USE dbcsr_cuda_memory,               ONLY: dbcsr_cuda_dev_mem_hold,&
                                             dbcsr_cuda_dev_mem_release
  USE dbcsr_cuda_types,                ONLY: dbcsr_cuda_mem_type
  USE dbcsr_data_methods,              ONLY: dbcsr_data_get_size,&
                                             dbcsr_data_hold,&
                                             dbcsr_data_init,&
                                             dbcsr_data_release
  USE dbcsr_error_handling
  USE dbcsr_pq_types,                  ONLY: &
       dbcsr_pq_all_type, dbcsr_pq_type, dbcsr_ps_group_type, dbcsr_ps_obj, &
       dbcsr_ps_set_type, dbcsr_ps_state_empty, dbcsr_ps_state_filled, &
       dbcsr_ps_state_filling, dbcsr_ps_state_queued, dbcsr_ps_state_working, &
       dbcsr_ps_target_obj, dbcsr_ps_target_type, dbcsr_ps_type, &
       dbcsr_ps_width, max_regions, p_c_first, ps_host_global, &
       ps_host_parameter
  USE dbcsr_ptr_util,                  ONLY: memory_allocate,&
                                             memory_deallocate
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_memory_CUDA_host_pinned

  !$ USE OMP_LIB


  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_pq_methods'

  LOGICAL, PARAMETER :: debug_mod = .FALSE.
  LOGICAL, PARAMETER :: careful_mod = .FALSE.
  LOGICAL, PARAMETER :: full_target_lock = .TRUE.
  LOGICAL, PARAMETER :: print_flush_level = .FALSE.

  INTEGER, SAVE, PUBLIC :: s_sid = 0
  INTEGER, SAVE, PRIVATE :: sid = 0

  PUBLIC :: dbcsr_ps_create, dbcsr_ps_create_p, dbcsr_ps_destroy
  PUBLIC :: dbcsr_ps_init, dbcsr_ps_new, dbcsr_ps_new_p,&
            dbcsr_ps_hold, dbcsr_ps_release

  PUBLIC :: dbcsr_ps_add_data_ab, dbcsr_ps_rm_data_ab
  PUBLIC :: dbcsr_ps_add_target, dbcsr_ps_rm_target

  PUBLIC :: dbcsr_ps_add_data_cuda_ab, dbcsr_ps_rm_data_cuda_ab

  PUBLIC :: dbcsr_ps_empty, dbcsr_ps_size
  PUBLIC :: dbcsr_ps_view_open, dbcsr_ps_view_close

  PUBLIC :: dbcsr_psg_add_data_ab, dbcsr_psg_rm_data_ab

  PUBLIC :: dbcsr_psg_add_data_cuda_ab
  PUBLIC :: dbcsr_psg_rm_data_cuda_ab

  PUBLIC :: dbcsr_psg_view_open, dbcsr_psg_view_close
  PUBLIC :: dbcsr_psg_set_state, dbcsr_psg_get_state

  PUBLIC :: dbcsr_ps_target_init, dbcsr_ps_target_new
  PUBLIC :: dbcsr_ps_target_create, dbcsr_ps_target_destroy
  PUBLIC :: dbcsr_ps_target_hold, dbcsr_ps_target_release
  PUBLIC :: dbcsr_ps_target_rm_data, dbcsr_ps_target_add_data
  PUBLIC :: dbcsr_ps_target_rm_data_cuda, dbcsr_ps_target_add_data_cuda

  PUBLIC :: dbcsr_ps_target_lock_main, dbcsr_ps_target_unlock_main
  PUBLIC :: dbcsr_ps_target_lock_regions, dbcsr_ps_target_unlock_regions

  PUBLIC :: dbcsr_ps_target_new_regions

  PUBLIC :: dbcsr_ps_set_init, dbcsr_ps_set_finalize
  PUBLIC :: dbcsr_ps_set_create, dbcsr_ps_set_destroy
  PUBLIC :: dbcsr_ps_set_set_b_master
  PUBLIC :: dbcsr_ps_set_advance,&
            dbcsr_ps_set_get_group_p,&
            dbcsr_ps_set_get_n_working

  PUBLIC :: dbcsr_pq_create, dbcsr_pq_destroy
  PUBLIC :: dbcsr_pq_share, dbcsr_pq_unshare, dbcsr_pq_refresh
  PUBLIC :: dbcsr_pq_add_stack
  PUBLIC :: dbcsr_pq_get_own_stack, dbcsr_pq_get_any_stack
  PUBLIC :: dbcsr_pq_get_far_stack, dbcsr_pq_get_defined_stack
  PUBLIC :: dbcsr_pq_mark_done
  PUBLIC :: dbcsr_pq_flush_level_chg
  PUBLIC :: dbcsr_pq_handoff_level_chg

CONTAINS


! *****************************************************************************
!> \brief Creates a queue of parameter stacks
!> \param[out] queue              Queue of parameter stacks to create.
!> \param[in] nbins               Number of different bins to create.
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_pq_create (queue, nbins, error)
    TYPE(dbcsr_pq_type), INTENT(OUT)         :: queue
    INTEGER, INTENT(IN)                      :: nbins
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    ALLOCATE (queue%stacks(20, nbins))
    queue%nstacks_total = 0
    ALLOCATE (queue%nstacks(nbins))
    queue%nstacks(:) = 0
    queue%nbins = nbins
    !$ ALLOCATE (queue%queue_lock)
    !$ queue%queue_lock = 0
    !$ CALL OMP_INIT_LOCK (queue%queue_lock)
    NULLIFY (queue%all_queues)
    CALL dbcsr_pq_refresh (queue, flush_level=0, error=error)
  END SUBROUTINE dbcsr_pq_create


! *****************************************************************************
!> \brief Destroys a queue of parameter stacks
!> \param[out] queue              Queue of parameter stacks to create.
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_pq_destroy (queue, error)
    TYPE(dbcsr_pq_type), INTENT(INOUT)       :: queue
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    DEALLOCATE (queue%stacks)
    DEALLOCATE (queue%nstacks)
    queue%nstacks_total = 0
    !$ CALL OMP_DESTROY_LOCK (queue%queue_lock)
    !$ DEALLOCATE (queue%queue_lock)
  END SUBROUTINE dbcsr_pq_destroy

! *****************************************************************************
!> \brief Creates a queue of parameter stacks
!> \param[out] queue              Queue of parameter stacks to create.
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_pq_refresh (queue, flush_level, error)
    TYPE(dbcsr_pq_type), INTENT(INOUT)       :: queue
    INTEGER, INTENT(IN)                      :: flush_level
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    IF (careful_mod) THEN
       CALL dbcsr_assert (queue%nstacks_total, "EQ", 0,&
            dbcsr_fatal_level, dbcsr_internal_error, "dbcsr_pq_refresh",&
            "Local queue not empty", __LINE__, error=error)
    ENDIF
    queue%nstacks(:) = 0
    queue%nstacks_total = 0
    queue%done = .FALSE.
    queue%flush_level = flush_level
    queue%handoff = 0
    queue%nworking = 0
    !$OMP BARRIER
    !$OMP MASTER
    !$ IF (ASSOCIATED (queue%all_queues)) THEN
    !$    queue%all_queues%n_done = 0
    !$    queue%all_queues%n_working = OMP_GET_NUM_THREADS()
    !$    queue%all_queues%working(:) = .TRUE.
    !$ ENDIF
    !$OMP END MASTER
    !$OMP BARRIER
  END SUBROUTINE dbcsr_pq_refresh


! *****************************************************************************
!> \brief Shares queues among threads
!> \param[in,out] queue           Queue of parameter stacks local to thread
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_pq_share (queue, error)
    TYPE(dbcsr_pq_type), POINTER             :: queue
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: t
    TYPE(dbcsr_pq_all_type), POINTER, SAVE   :: all_queues

! This save here is just to get shared access...

    !$OMP MASTER
    ALLOCATE (all_queues)
    all_queues%n_done = 0
    all_queues%n_working = 1
    !$ all_queues%n_working = OMP_GET_NUM_THREADS()
    ALLOCATE (all_queues%queues(all_queues%n_working))
    ALLOCATE (all_queues%working(all_queues%n_working))
    all_queues%working(:) = .TRUE.
    !$OMP END MASTER
    !$OMP BARRIER
    queue%all_queues => all_queues
    t = 0
    !$ t = OMP_GET_THREAD_NUM ()
    all_queues%queues(t+1)%pq => queue
    !$OMP BARRIER
  END SUBROUTINE dbcsr_pq_share

! *****************************************************************************
!> \brief Shares queues among threads
!> \param[in,out] queue           Queue of parameter stacks local to thread
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_pq_unshare (queue, error)
    TYPE(dbcsr_pq_type), POINTER             :: queue
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    !$OMP BARRIER
    !$OMP MASTER
    DEALLOCATE (queue%all_queues%queues)
    DEALLOCATE (queue%all_queues%working)
    DEALLOCATE (queue%all_queues)
    !$OMP END MASTER
  END SUBROUTINE dbcsr_pq_unshare


! *****************************************************************************
!> \brief Marks a queue as done
!> \param[in,out] queue           Queue of parameter stacks local to thread
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_pq_mark_done (queue, error)
    TYPE(dbcsr_pq_type), INTENT(INOUT)       :: queue
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: t

    !$OMP FLUSH
    queue%done = .TRUE.
    t = 0
    !$ t = OMP_GET_THREAD_NUM ()
    queue%all_queues%working(t+1) = .FALSE.
    !$OMP ATOMIC
    queue%all_queues%n_done = queue%all_queues%n_done + 1
    !$OMP ATOMIC
    queue%all_queues%n_working = queue%all_queues%n_working - 1
    !$OMP FLUSH
  END SUBROUTINE dbcsr_pq_mark_done

! *****************************************************************************
!> \brief Changes the flush level of a queue
!> \see process_queue_preempt
!> \param[in,out] queue           Queue of parameter stacks
!> \param[in] chg                 Amount by which to change the flush level
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_pq_flush_level_chg (queue, chg, param_sets)
    TYPE(dbcsr_pq_type), INTENT(INOUT)       :: queue
    INTEGER, INTENT(IN)                      :: chg
    TYPE(dbcsr_ps_set_type), INTENT(IN)      :: param_sets

    queue%flush_level = MAX(0,&
                        MIN(param_sets%nbuffers*param_sets%group_size,&
                        queue%flush_level + chg))
    IF (print_flush_level) THEN
       !$ IF (.FALSE.) THEN
       WRITE(*,*)"GPU flush level", chg, queue%flush_level
       !$ ELSE
       !$ WRITE(*,*)"GPU flush level", OMP_GET_THREAD_NUM(), chg, queue%flush_level
       !$ ENDIF
    ENDIF
  END SUBROUTINE dbcsr_pq_flush_level_chg

! *****************************************************************************
!> \brief Changes the handoff level of a queue
!> \see process_queue_preempt
!> \param[in,out] queue           Queue of parameter stacks
!> \param[in] chg                 Amount by which to change the flush level
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_pq_handoff_level_chg (queue, chg, param_sets)
    TYPE(dbcsr_pq_type), INTENT(INOUT)       :: queue
    INTEGER, INTENT(IN)                      :: chg
    TYPE(dbcsr_ps_set_type), INTENT(IN)      :: param_sets

    queue%handoff = MAX(0,&
                        MIN(param_sets%nbuffers*param_sets%group_size,&
                        queue%handoff + chg))
    IF (print_flush_level) THEN
       !$ IF (.FALSE.) THEN
       WRITE(*,*)"GPU handoff lvl", chg, queue%handoff
       !$ ELSE
       !$ WRITE(*,*)"GPU handoff lvl", OMP_GET_THREAD_NUM(), chg, queue%handoff
       !$ ENDIF
    ENDIF
  END SUBROUTINE dbcsr_pq_handoff_level_chg


! *****************************************************************************
!> \brief Adds a stack to a queue
!>
!> Stacks with common m, n, and k values are added to the end while mixed
!> stacks are added to the beginning of the queue
!> \param[in,out] queue           Add stacks to this queue of parameter stacks
!> \param[in] stack               Stack to add to queue
!> \param[in,out] error           error
!> \note The queue parameter could be INTENT(INOUT) when OpenMP is not used.
! *****************************************************************************
  SUBROUTINE dbcsr_pq_add_stack (queue, stack, error)
    TYPE(dbcsr_pq_type), POINTER             :: queue
    TYPE(dbcsr_ps_obj), INTENT(IN)           :: stack
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_pq_add_stack', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bin, error_handle, i, new_n, &
                                                si
    TYPE(dbcsr_ps_obj), DIMENSION(:, :), &
      POINTER                                :: new_stacks

    IF (careful_mod) &
         CALL dbcsr_error_set(routineN, error_handle, error)
    !$ CALL OMP_SET_LOCK (queue%queue_lock)
    IF (mm_thread_workshare .AND. .NOT. has_omp_2_5) THEN
       !$OMP FLUSH
    ENDIF
    !$OMP CRITICAL
    s_sid = s_sid + 1
    si = s_sid
    !$OMP END CRITICAL
    IF (stack%s%stack_p .EQ. 0 .AND. debug_mod) THEN
       WRITE(*,*)"Adding 0 size"
    ENDIF
    IF (.FALSE.) THEN
       WRITE(*,*)routineN//" Adding stack size",&
            stack%s%stack_p, "ids", stack%s%id, si
    ENDIF
    !
    IF (.TRUE.) THEN
       CALL dbcsr_assert (stack%s%size_bin, "LE", queue%nbins,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Invalid size bin", __LINE__, error=error)
       CALL dbcsr_assert (stack%s%size_bin, "GE", 1,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Invalid size bin", __LINE__, error=error)
    ENDIF
    !
    bin = stack%s%size_bin
    new_n = queue%nstacks(bin) + 1
    !
    IF (new_n .GT. SIZE(queue%stacks, 1)) THEN
       ! Resize the stacks
       ALLOCATE (new_stacks(SIZE(queue%stacks, 1), queue%nbins))
       DO i = 1, queue%nbins
          new_stacks(1:queue%nstacks(bin), i) =&
               queue%stacks(1:queue%nstacks(bin), i)
       ENDDO
       DEALLOCATE (queue%stacks)
       ALLOCATE (queue%stacks(MAX(new_n, 2*SIZE(queue%stacks, 1)), queue%nbins))
       DO i = 1, queue%nbins
          queue%stacks(1:queue%nstacks(bin), i) = &
               new_stacks(1:queue%nstacks(bin), i)
       ENDDO
       DEALLOCATE (new_stacks)
    ENDIF
    !
    queue%stacks(new_n,bin) = stack
    queue%nstacks(bin) = queue%nstacks(bin) + 1
    queue%nstacks_total = queue%nstacks_total + 1
    IF (mm_thread_workshare .AND. .NOT. has_omp_2_5) THEN
       !$OMP FLUSH
    ENDIF
    !$ CALL OMP_UNSET_LOCK (queue%queue_lock)
    IF (careful_mod) &
         CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_pq_add_stack


! *****************************************************************************
!> \brief Gets a stack from my queue.
!> \note Currently makes this data structure act as a stack and not a queue.
!> \param[in,out] queue              Queue of parameter stacks to create.
!> \param[out] stack                 Stack to processes according to criteria
!> \param[out] found                 Whether a stack matching the criteria
!>                                   was found
!> \param[in] easy                   Stack should not be friendly to
!>                                   accelerators
!> \param[in,out] error           error
!> \note The queue parameter could be INTENT(INOUT) when OpenMP is not used.
! *****************************************************************************
  SUBROUTINE dbcsr_pq_get_own_stack (queue, stack, found, easy,&
       error)
    TYPE(dbcsr_pq_type), POINTER             :: queue
    TYPE(dbcsr_ps_obj), INTENT(OUT)          :: stack
    LOGICAL, INTENT(OUT)                     :: found
    LOGICAL, INTENT(IN)                      :: easy
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    IF (easy) THEN
       CALL dbcsr_pq_get_defined_stack (queue, stack, found, error=error)
    ELSE
       CALL dbcsr_pq_get_far_stack (queue, stack, found, error=error)
    ENDIF
  END SUBROUTINE dbcsr_pq_get_own_stack

! *****************************************************************************
!> \brief Gets a stack from my queue or others' queues.
!>
!> Other threads' stacks are taken if the other thread is still indexing or
!> the anyways argument is set.
!> \param[in,out] queue              Queue of parameter stacks to create.
!> \param[out] stack                 Stack to processes according to criteria
!> \param[out] found                 Whether a stack matching the criteria
!>                                   was found
!> \param[in] easy                   Stack should be friendly to accelerators
!> \param[in] anyways                Get stacks even if other thread is done
!>                                   indexing
!> \param[in,out] thread_ptr         thread counter for stack skimming
!> \param[in,out] error           error
!> \note The queue parameter could be INTENT(INOUT) when OpenMP is not used.
! *****************************************************************************
  SUBROUTINE dbcsr_pq_get_any_stack (queue, stack, found, easy, &
       anyways, thread_ptr, error)
    TYPE(dbcsr_pq_type), POINTER             :: queue
    TYPE(dbcsr_ps_obj), INTENT(OUT)          :: stack
    LOGICAL, INTENT(OUT)                     :: found
    LOGICAL, INTENT(IN)                      :: easy, anyways
    INTEGER, INTENT(INOUT)                   :: thread_ptr
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: i, myt, nt, t

    found = .FALSE.
    nt = 1
    myt = 0
    !$ myt = OMP_GET_THREAD_NUM ()
    !$ nt = OMP_GET_NUM_THREADS ()
    !$OMP FLUSH
    DO i = thread_ptr, thread_ptr + nt-1
       t = MOD (myt + i, nt)+1
       IF (anyways .OR. queue%all_queues%working(t)) THEN
          IF (easy) THEN
             CALL dbcsr_pq_get_defined_stack (queue%all_queues%queues(t)%pq,&
                  stack, found, error=error)
          ELSE
             CALL dbcsr_pq_get_far_stack (queue%all_queues%queues(t)%pq,&
                  stack, found, error=error)
          ENDIF
          IF (found) &
               EXIT
       ENDIF
    ENDDO
    thread_ptr = MOD(i, nt)+1
  END SUBROUTINE dbcsr_pq_get_any_stack


! *****************************************************************************
!> \brief Gets a stack from the queue from the end.
!> \param[in,out] queue              Queue of parameter stacks to create.
!> \param[in,out] error           error
!> \note The queue parameter could be INTENT(INOUT) when OpenMP is not used.
! *****************************************************************************
  SUBROUTINE dbcsr_pq_get_far_stack (queue, stack, found, error)
    TYPE(dbcsr_pq_type), POINTER             :: queue
    TYPE(dbcsr_ps_obj), INTENT(OUT)          :: stack
    LOGICAL, INTENT(OUT)                     :: found
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: bin

    !$ CALL OMP_SET_LOCK (queue%queue_lock)
    IF (mm_thread_workshare .AND. .NOT. has_omp_2_5) THEN
       !$OMP FLUSH
    ENDIF
    found = .FALSE.
    IF (queue%nstacks_total .GE. 1) THEN
       DO bin = queue%nbins, 1, -1
          IF (queue%nstacks(bin) .GT. 0) THEN
             stack = queue%stacks(queue%nstacks(bin),bin)
             queue%nstacks_total = queue%nstacks_total - 1
             queue%nstacks(bin) = queue%nstacks(bin) - 1
             found = .TRUE.
             EXIT
          ENDIF
       ENDDO
    ELSE
       found = .FALSE.
    ENDIF
    IF (mm_thread_workshare .AND. .NOT. has_omp_2_5) THEN
       !$OMP FLUSH
    ENDIF
    !$ CALL OMP_UNSET_LOCK (queue%queue_lock)
  END SUBROUTINE dbcsr_pq_get_far_stack



! *****************************************************************************
!> \brief Gets a stack with defined m, n, and k from the queue.
!> \param[in,out] queue              Queue of parameter stacks to create.
!> \param[in,out] error           error
!> \note The queue parameter could be INTENT(INOUT) when OpenMP is not used.
! *****************************************************************************
  SUBROUTINE dbcsr_pq_get_defined_stack (queue, stack, found, error)
    TYPE(dbcsr_pq_type), POINTER             :: queue
    TYPE(dbcsr_ps_obj), INTENT(OUT)          :: stack
    LOGICAL, INTENT(OUT)                     :: found
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: bin
    LOGICAL                                  :: dbg = .FALSE.

    !$ CALL OMP_SET_LOCK (queue%queue_lock)
    IF (mm_thread_workshare .AND. .NOT. has_omp_2_5) THEN
       !$OMP FLUSH
    ENDIF
    !$ IF (dbg) WRITE(*,*)'checking stacks', omp_get_thread_num(),&
    !$      queue%nstacks_total, queue%nstacks
    found = .FALSE.
    IF (queue%nstacks_total .GE. 1) THEN
       DO bin = 1, queue%nbins
          IF (queue%nstacks(bin) .GT. 0) THEN
             IF (queue%stacks(queue%nstacks(bin),bin)%s%defined_mnk) THEN
                stack = queue%stacks(queue%nstacks(bin),bin)
                queue%nstacks(bin) = queue%nstacks(bin) - 1
                queue%nstacks_total = queue%nstacks_total - 1
                found = .TRUE.
                EXIT
             ELSE
                EXIT
             ENDIF
          ENDIF
       ENDDO
    ELSE
       found = .FALSE.
    ENDIF
    IF (mm_thread_workshare .AND. .NOT. has_omp_2_5) THEN
       !$OMP FLUSH
    ENDIF
    !$ CALL OMP_UNSET_LOCK (queue%queue_lock)
    !$ IF (dbg .AND. found) &
    !$     write(*,*)"got stack", omp_get_thread_num(), "id", stack%s%sid
  END SUBROUTINE dbcsr_pq_get_defined_stack


  SUBROUTINE dbcsr_ps_target_init (t, error)
    TYPE(dbcsr_ps_target_obj), INTENT(OUT)   :: t
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    NULLIFY (t%t)
  END SUBROUTINE dbcsr_ps_target_init

  SUBROUTINE dbcsr_ps_target_new (t, error)
    TYPE(dbcsr_ps_target_obj), INTENT(OUT)   :: t
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    ALLOCATE (t%t)
    CALL dbcsr_ps_target_create (t%t, error)
    t%t%refcount = 1
  END SUBROUTINE dbcsr_ps_target_new

  SUBROUTINE dbcsr_ps_target_hold (t, error)
    TYPE(dbcsr_ps_target_obj), INTENT(INOUT) :: t
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    t%t%refcount = t%t%refcount + 1
  END SUBROUTINE dbcsr_ps_target_hold

  SUBROUTINE dbcsr_ps_target_release (t, error)
    TYPE(dbcsr_ps_target_obj), INTENT(INOUT) :: t
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    IF (ASSOCIATED (t%t)) THEN
       t%t%refcount = t%t%refcount - 1
       IF (t%t%refcount .EQ. 0) THEN
          CALL dbcsr_ps_target_destroy (t%t, error)
          DEALLOCATE (t%t)
          NULLIFY (t%t)
       ENDIF
    ENDIF
  END SUBROUTINE dbcsr_ps_target_release

  SUBROUTINE dbcsr_ps_target_create (t, error)
    TYPE(dbcsr_ps_target_type), INTENT(OUT)  :: t
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    !$ INTEGER :: i

    t%zero_first = 1
    t%zero_last = 0
    t%has_c_data = .FALSE.
    t%has_cuda_c_data = .FALSE.
    t%last_c_blk = 0
    NULLIFY (t%c_locks_dev)
    NULLIFY (t%stack_state_dev)
    NULLIFY (t%params_dev)
    t%owner = 0
    !$ t%owner = OMP_GET_THREAD_NUM ()
    !$ t%owner_wants_lock = .FALSE.
    !$ t%n_use = 0
    !$ CALL OMP_INIT_LOCK (t%target_lock)
    !$ ALLOCATE (t%region_locks(1024))
    !$ DO i = 1, SIZE(t%region_locks)
    !$    CALL OMP_INIT_LOCK (t%region_locks(i))
    !$ ENDDO
  END SUBROUTINE dbcsr_ps_target_create

  SUBROUTINE dbcsr_ps_target_destroy (t, error)
    TYPE(dbcsr_ps_target_type), &
      INTENT(INOUT)                          :: t
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    !$ INTEGER :: i

    IF (t%has_c_data) THEN
       CALL dbcsr_data_release (t%product_data_area)
    ENDIF
    IF (t%has_cuda_c_data) THEN
       CALL dbcsr_cuda_dev_mem_release (t%product_data_cuda, error=error)
    ENDIF
    IF (ASSOCIATED (t%c_locks_dev)) THEN
       CALL dbcsr_cuda_dev_mem_release (t%c_locks_dev, error=error)
    ENDIF
    ! the target of this pointer was allocated by the master, and shared between all threads
    !$OMP MASTER
    IF (ASSOCIATED (t%stack_state_dev)) THEN
       CALL dbcsr_cuda_dev_mem_release (t%stack_state_dev, error=error)
    ENDIF
    !$OMP END MASTER
    NULLIFY(t%stack_state_dev)
    IF (ASSOCIATED (t%params_dev)) THEN
       CALL dbcsr_cuda_dev_mem_release (t%params_dev, error=error)
    ENDIF
    !$ CALL OMP_DESTROY_LOCK (t%target_lock)
    !$ DO i = 1, SIZE(t%region_locks)
    !$    CALL OMP_DESTROY_LOCK (t%region_locks(i))
    !$ ENDDO
    !$ DEALLOCATE (t%region_locks)
  END SUBROUTINE dbcsr_ps_target_destroy

  SUBROUTINE dbcsr_ps_target_rm_data (t, error)
    TYPE(dbcsr_ps_target_type), &
      INTENT(INOUT)                          :: t
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    IF (t%has_c_data) THEN
       CALL dbcsr_data_release (t%product_data_area)
    ENDIF
  END SUBROUTINE dbcsr_ps_target_rm_data

  SUBROUTINE dbcsr_ps_target_rm_data_cuda (t, error)
    TYPE(dbcsr_ps_target_type), &
      INTENT(INOUT)                          :: t
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    IF (t%has_cuda_c_data) THEN
       CALL dbcsr_cuda_dev_mem_release (t%product_data_cuda, error=error)
    ENDIF
  END SUBROUTINE dbcsr_ps_target_rm_data_cuda

  SUBROUTINE dbcsr_ps_target_add_data (t, product_data, last_nze, error)
    TYPE(dbcsr_ps_target_type), &
      INTENT(INOUT)                          :: t
    TYPE(dbcsr_data_obj), INTENT(IN)         :: product_data
    INTEGER, INTENT(IN)                      :: last_nze
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    IF (t%has_c_data) THEN
       CALL dbcsr_ps_target_rm_data (t, error)
    ENDIF
    t%product_data_area = product_data
    CALL dbcsr_data_hold (t%product_data_area)
    t%zero_first = last_nze+1
    t%zero_last = last_nze
    t%has_c_data = .TRUE.
    !$ CALL dbcsr_ps_target_new_regions (t,&
    !$      dbcsr_data_get_size (t%product_data_area),&
    !$      error=error)
  END SUBROUTINE dbcsr_ps_target_add_data

  SUBROUTINE dbcsr_ps_target_add_data_cuda (t, product_cuda_data, error)
    TYPE(dbcsr_ps_target_type), &
      INTENT(INOUT)                          :: t
    TYPE(dbcsr_cuda_mem_type), POINTER       :: product_cuda_data
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    IF (t%has_cuda_c_data) THEN
       CALL dbcsr_ps_target_rm_data_cuda (t, error)
    ENDIF
    t%product_data_cuda => product_cuda_data
    CALL dbcsr_cuda_dev_mem_hold (t%product_data_cuda, error)
    t%has_cuda_c_data = .TRUE.
  END SUBROUTINE dbcsr_ps_target_add_data_cuda




  SUBROUTINE dbcsr_ps_init (stack, error)
    TYPE(dbcsr_ps_obj), INTENT(INOUT)        :: stack
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    NULLIFY (stack%s)
  END SUBROUTINE dbcsr_ps_init

  SUBROUTINE dbcsr_ps_hold (stack, error)
    TYPE(dbcsr_ps_obj), INTENT(INOUT)        :: stack
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    stack%s%refcount = stack%s%refcount + 1
  END SUBROUTINE dbcsr_ps_hold

  SUBROUTINE dbcsr_ps_release (stack, error)
    TYPE(dbcsr_ps_obj), INTENT(INOUT)        :: stack
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    IF (ASSOCIATED (stack%s)) THEN
       stack%s%refcount = stack%s%refcount - 1
       IF (stack%s%refcount .EQ. 0) THEN
          CALL dbcsr_ps_destroy (stack%s, error)
          DEALLOCATE (stack%s)
          NULLIFY (stack%s)
       ENDIF

    ENDIF
  END SUBROUTINE dbcsr_ps_release

  SUBROUTINE dbcsr_ps_new (stack, stack_size, mem_type, state, error)
    TYPE(dbcsr_ps_obj), INTENT(OUT)          :: stack
    INTEGER, INTENT(IN)                      :: stack_size, mem_type
    INTEGER, OPTIONAL, POINTER               :: state
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    ALLOCATE (stack%s)
    CALL dbcsr_ps_create (stack%s, stack_size, mem_type, state, error)
    stack%s%refcount = 1
  END SUBROUTINE dbcsr_ps_new

  SUBROUTINE dbcsr_ps_new_p (stack, master_stack, stack_size, stack_start, state, error)
    TYPE(dbcsr_ps_obj), INTENT(OUT)          :: stack
    TYPE(dbcsr_ps_obj), INTENT(IN)           :: master_stack
    INTEGER, INTENT(IN)                      :: stack_size, stack_start
    INTEGER, OPTIONAL, POINTER               :: state
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    ALLOCATE (stack%s)
    CALL dbcsr_ps_create_p (stack%s, master_stack%s, stack_size, stack_start, state, error)
    stack%s%refcount = 1
  END SUBROUTINE dbcsr_ps_new_p


! *****************************************************************************
!> \brief Creates an independent parameter stack
!>
!>        This just creates the structure; nothing of value is
!>        allocated or instantiated.
!> \param[out] stack          Parameter stack
!> \param[in] stack_size      Stack size
!> \param[in] mem_type        Memory type of stack data
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_ps_create(stack, stack_size, mem_type, state, error)
    TYPE(dbcsr_ps_type), INTENT(OUT)         :: stack
    INTEGER, INTENT(IN)                      :: stack_size, mem_type
    INTEGER, OPTIONAL, POINTER               :: state
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_ps_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set(routineN, error_handle, error)
    CALL ps_clear (stack, state, error)
    stack%mem_type = mem_type
      CALL memory_allocate (stack%parameters, stack_size*dbcsr_ps_width,&
         mem_type=mem_type, error=error)
    stack%own_data = .TRUE.
    stack%start_p = 1
    !$omp critical
    sid = sid+1
    stack%id = -1
    stack%id = sid
    !$omp end critical
    stack%id = -1
    stack%sid = -7
    IF (careful_mod) CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_ps_create


! *****************************************************************************
!> \brief Creates a parameter stack that points into an existing stack
!> \param[out] stack          Parameter stack
!> \param[in] master_stack    Existing parameter stack to point to
!> \param[in] stack_size      Stack size
!> \param[in] stack_start     Starting position in existing parameter stack
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_ps_create_p(stack, master_stack, stack_size, stack_start,&
       state, error)
    TYPE(dbcsr_ps_type), INTENT(OUT)         :: stack
    TYPE(dbcsr_ps_type), INTENT(IN)          :: master_stack
    INTEGER, INTENT(IN)                      :: stack_size, stack_start
    INTEGER, OPTIONAL, POINTER               :: state
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_ps_create_p', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set(routineN, error_handle, error)
    CALL ps_clear(stack, state, error)
    IF (careful_mod) THEN
       CALL dbcsr_assert ((stack_start+stack_size-1)*dbcsr_ps_width, "LE", &
            SIZE(master_stack%parameters),&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Can not point outside of master stack bounds.",&
            __LINE__, error=error)
    ENDIF
    stack%parameters => master_stack%parameters(&
         (stack_start-1)           *dbcsr_ps_width+1 : &
         (stack_start+stack_size-1)*dbcsr_ps_width)
    stack%own_data = .FALSE.
    stack%start_p = stack_start
    !$omp critical
    sid = sid+1
    stack%id = -1
    stack%id = sid
    !$omp end critical
    stack%sid = -7
    IF (careful_mod) CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_ps_create_p


! *****************************************************************************
!> \brief Initializes data in a stack
! *****************************************************************************
  SUBROUTINE ps_clear(stack, state, error)
    TYPE(dbcsr_ps_type), INTENT(OUT)         :: stack
    INTEGER, OPTIONAL, POINTER               :: state
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    stack%mem_type = 0
    stack%start_p = 1
    stack%stack_p = 0
    IF (PRESENT (state)) THEN
       stack%state => state
       stack%own_state = .FALSE.
    ELSE
       ALLOCATE (stack%state)
       stack%own_state = .TRUE.
    ENDIF
    stack%state = dbcsr_ps_state_empty
    NULLIFY (stack%parameters)
    stack%own_data = .FALSE.
    stack%driver = 0
    stack%driver_desc = 0
    stack%has_ab_data = .FALSE.
    stack%has_cuda_ab_data = .FALSE.
    CALL dbcsr_data_init (stack%left_data_area)
    CALL dbcsr_data_init (stack%right_data_area)
    stack%has_target = .FALSE.
    CALL dbcsr_ps_target_init (stack%t, error)
    stack%m = 0
    stack%n = 0
    stack%k = 0
    stack%defined_mnk = .FALSE.
    stack%size_bin = -1
  END SUBROUTINE ps_clear


! *****************************************************************************
!> \brief Destroys a parameter stack
!> \param[in,out] stack          Parameter stack
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_ps_destroy(stack, error)
    TYPE(dbcsr_ps_type), INTENT(INOUT)       :: stack
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_ps_destroy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set(routineN, error_handle, error)
    IF (careful_mod) CALL dbcsr_assert(stack%stack_p .EQ. 0,&
         dbcsr_warning_level, dbcsr_caller_error, routineN,&
         "Destroying a non-empty stack.",&
         __LINE__, error=error)
    IF (stack%own_data) THEN
     IF (stack%mem_type .NE. dbcsr_memory_CUDA_host_pinned) THEN
       CALL memory_deallocate (stack%parameters, stack%mem_type, error=error)
     ENDIF
    ENDIF
    CALL dbcsr_ps_rm_data_ab (stack, error)
    CALL dbcsr_ps_rm_data_cuda_ab (stack, error)
    CALL dbcsr_ps_rm_target (stack, error)
    IF (stack%own_state) THEN
       DEALLOCATE (stack%state)
    ENDIF
    stack%stack_p = -1
    IF (careful_mod) CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_ps_destroy


! *****************************************************************************
!> \brief Removes left and right matrix data areas from a stack
!> \param[in,out] stack           Parameter stack
!> \param[in] left_data           Left matrix data area
!> \param[in] right_data          Right matrix data area
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_ps_rm_data_ab(stack, error)
    TYPE(dbcsr_ps_type), INTENT(INOUT)       :: stack
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    IF (stack%has_ab_data) THEN
       CALL dbcsr_data_release (stack%left_data_area)
       CALL dbcsr_data_release (stack%right_data_area)
    ENDIF
    stack%has_ab_data = .FALSE.
  END SUBROUTINE dbcsr_ps_rm_data_ab

! *****************************************************************************
!> \brief Removes product matrix data area from a stack
!> \param[in,out] stack           Parameter stack
!> \param[in] left_data           Left matrix data area
!> \param[in] right_data          Right matrix data area
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_ps_rm_target(stack, error)
    TYPE(dbcsr_ps_type), INTENT(INOUT)       :: stack
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    IF (stack%has_target) THEN
       CALL dbcsr_ps_target_release(stack%t, error)
       stack%has_target = .FALSE.
    ENDIF
  END SUBROUTINE dbcsr_ps_rm_target


! *****************************************************************************
!> \brief Adds left and right matrix data areas to a stack
!> \param[in,out] stack           Parameter stack
!> \param[in] left_data           Left matrix data area
!> \param[in] right_data          Right matrix data area
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_ps_add_data_ab(stack, left_data, right_data, error)
    TYPE(dbcsr_ps_type), INTENT(INOUT)       :: stack
    TYPE(dbcsr_data_obj), INTENT(IN)         :: left_data, right_data
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    IF (stack%has_ab_data) THEN
       CALL dbcsr_ps_rm_data_ab (stack, error)
    ENDIF
    stack%left_data_area = left_data
    CALL dbcsr_data_hold (stack%left_data_area)
    stack%right_data_area = right_data
    CALL dbcsr_data_hold (stack%right_data_area)
    stack%has_ab_data = .TRUE.
  END SUBROUTINE dbcsr_ps_add_data_ab

! *****************************************************************************
!> \brief Adds left and right matrix data areas to a stack
!> \param[in,out] stack           Parameter stack
!> \param[in] product_data        Product matrix data area
!> \param[in] zero_first          Size of product matrix data area
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_ps_add_target(stack, t, error)
    TYPE(dbcsr_ps_type), INTENT(INOUT)       :: stack
    TYPE(dbcsr_ps_target_obj), INTENT(IN)    :: t
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    stack%t = t
    CALL dbcsr_ps_target_hold (stack%t, error)
    stack%has_target = .TRUE.
  END SUBROUTINE dbcsr_ps_add_target


! *****************************************************************************
!> \brief Removes left and right matrix CUDA data areas from a stack
!> \param[in,out] stack           Parameter stack
!> \param[in] left_data           Left matrix data area
!> \param[in] right_data          Right matrix data area
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_ps_rm_data_cuda_ab(stack, error)
    TYPE(dbcsr_ps_type), INTENT(INOUT)       :: stack
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    IF (stack%has_cuda_ab_data) THEN
       CALL dbcsr_cuda_dev_mem_release (stack%left_data_cuda, error=error)
       CALL dbcsr_cuda_dev_mem_release (stack%right_data_cuda, error=error)
    ENDIF
    stack%has_cuda_ab_data = .FALSE.
  END SUBROUTINE dbcsr_ps_rm_data_cuda_ab

! *****************************************************************************
!> \brief Adds left and right matrix data areas to a stack
!> \param[in,out] stack           Parameter stack
!> \param[in] left_data           Left matrix data area
!> \param[in] right_data          Right matrix data area
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_ps_add_data_cuda_ab(stack, left_cuda_data, right_cuda_data, error)
    TYPE(dbcsr_ps_type), INTENT(INOUT)       :: stack
    TYPE(dbcsr_cuda_mem_type), POINTER       :: left_cuda_data, &
                                                right_cuda_data
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    IF (stack%has_cuda_ab_data) THEN
       CALL dbcsr_ps_rm_data_cuda_ab (stack, error)
    ENDIF
    stack%left_data_cuda => left_cuda_data
    CALL dbcsr_cuda_dev_mem_hold (stack%left_data_cuda, error=error)
    stack%right_data_cuda => right_cuda_data
    CALL dbcsr_cuda_dev_mem_hold (stack%right_data_cuda, error=error)
    stack%has_cuda_ab_data = .TRUE.
  END SUBROUTINE dbcsr_ps_add_data_cuda_ab

! *****************************************************************************
!> \brief Removes left and right matrix data areas from a stack group
!> \param[in,out] group           Parameter stack group
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_psg_rm_data_ab(group, error)
    TYPE(dbcsr_ps_group_type), INTENT(INOUT) :: group
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: i

    CALL dbcsr_ps_rm_data_ab (group%master%s, error)
    DO i = 1, SIZE (group%stacks)
       CALL dbcsr_ps_rm_data_ab (group%stacks(i)%s, error)
    ENDDO
  END SUBROUTINE dbcsr_psg_rm_data_ab

! *****************************************************************************
!> \brief Adds left and right matrix data areas to a stack group
!> \param[in,out] group           Parameter stack group
!> \param[in] left_data           Left matrix data area
!> \param[in] right_data          Right matrix data area
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_psg_add_data_ab(group, left_data, right_data, error)
    TYPE(dbcsr_ps_group_type), INTENT(INOUT) :: group
    TYPE(dbcsr_data_obj), INTENT(IN)         :: left_data, right_data
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: i

    CALL dbcsr_ps_add_data_ab (group%master%s, left_data, right_data, error)
    DO i = 1, SIZE (group%stacks)
       CALL dbcsr_ps_add_data_ab (group%stacks(i)%s, left_data, right_data, error)
    ENDDO
  END SUBROUTINE dbcsr_psg_add_data_ab

! *****************************************************************************
!> \brief Removes left and right matrix CUDA data areas from a stack group
!> \param[in,out] group           Parameter stack group
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_psg_rm_data_cuda_ab(group, error)
    TYPE(dbcsr_ps_group_type), INTENT(INOUT) :: group
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: i

    CALL dbcsr_ps_rm_data_cuda_ab (group%master%s, error)
    DO i = 1, SIZE (group%stacks)
       CALL dbcsr_ps_rm_data_cuda_ab (group%stacks(i)%s, error)
    ENDDO
  END SUBROUTINE dbcsr_psg_rm_data_cuda_ab

! *****************************************************************************
!> \brief Adds left and right matrix data areas to a stack group
!> \param[in,out] group           Parameter stack group
!> \param[in] left_data           Left matrix data area
!> \param[in] right_data          Right matrix data area
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_psg_add_data_cuda_ab(group, left_cuda_data, right_cuda_data, error)
    TYPE(dbcsr_ps_group_type), INTENT(INOUT) :: group
    TYPE(dbcsr_cuda_mem_type), POINTER       :: left_cuda_data, &
                                                right_cuda_data
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: i

    CALL dbcsr_ps_add_data_cuda_ab (group%master%s, left_cuda_data, right_cuda_data, error)
    DO i = 1, SIZE (group%stacks)
       CALL dbcsr_ps_add_data_cuda_ab (group%stacks(i)%s, left_cuda_data, right_cuda_data, error)
    ENDDO
  END SUBROUTINE dbcsr_psg_add_data_cuda_ab


! *****************************************************************************
!> \brief Direct access to paramater stack
!> \par View details
!>      A pointer to the parameters array is returned.  Other data is
!>      returned as copies.  Therefore any changes must later be committed with
!>      dbcsr_ps_view_close.
!> \param[in] stack               Parameter stack
!> \param[out] parameters         Pointer to parameter stack
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_ps_view_open(stack, parameters,&
       stack_p, zero_first, zero_last, last_c_blk, error)
    TYPE(dbcsr_ps_type), INTENT(INOUT)       :: stack
    INTEGER, DIMENSION(:), POINTER           :: parameters
    INTEGER, INTENT(OUT)                     :: stack_p, zero_first, &
                                                zero_last, last_c_blk
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    parameters => stack%parameters
    stack_p = stack%stack_p
    zero_first = stack%t%t%zero_first
    zero_last = stack%t%t%zero_last
    last_c_blk = stack%t%t%last_c_blk
    stack%state = dbcsr_ps_state_filling
  END SUBROUTINE dbcsr_ps_view_open

! *****************************************************************************
!> \brief Direct access into stack
!> \param[in] stack               Parameter stack
!> \param[in] parameters          Pointer to parameter stack
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_ps_view_close(stack, parameters,&
       stack_p, zero_first, zero_last, last_c_blk, error)
    TYPE(dbcsr_ps_type), INTENT(INOUT)       :: stack
    INTEGER, DIMENSION(:), POINTER           :: parameters
    INTEGER, INTENT(IN)                      :: stack_p, zero_first, &
                                                zero_last, last_c_blk
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    stack%stack_p = stack_p
    stack%t%t%zero_first = zero_first
    stack%t%t%zero_last = zero_last
    stack%t%t%last_c_blk = last_c_blk
    stack%state = dbcsr_ps_state_filled
  END SUBROUTINE dbcsr_ps_view_close

! *****************************************************************************
!> \brief Returns whether a stack is empty
!> \result True if the stack is empty, otherwise False.
! *****************************************************************************
  LOGICAL FUNCTION dbcsr_ps_empty(stack)
    TYPE(dbcsr_ps_type), INTENT(INOUT) :: stack
    dbcsr_ps_empty = stack%stack_p .LT. 0
  END FUNCTION dbcsr_ps_empty

! *****************************************************************************
!> \brief Returns the size of the stack
!> \result The size of the stack
! *****************************************************************************
  INTEGER FUNCTION dbcsr_ps_size(stack)
    TYPE(dbcsr_ps_type), INTENT(INOUT) :: stack
    dbcsr_ps_size = SIZE (stack%parameters) / dbcsr_ps_width
  END FUNCTION dbcsr_ps_size


! *****************************************************************************
!> \brief Direct access to paramater stack group
!> \par View details
!>      A pointer to the parameters array is returned.  Other data is
!>      returned as copies.  Therefore any changes must later be committed with
!>      dbcsr_ps_view_close.
!> \param[in] group               Parameter stack group
!> \param[out] parameters         Pointers to the master parameter stack
!> \param[out] starts             Individual stack starting positions in the
!>                                master stack
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_psg_view_open(group, parameters,&
       stack_p, zero_first, zero_last, last_c_blk, starts, error)
    TYPE(dbcsr_ps_group_type), INTENT(INOUT) :: group
    INTEGER, DIMENSION(:), POINTER           :: parameters
    INTEGER, DIMENSION(:), INTENT(OUT)       :: stack_p
    INTEGER, INTENT(OUT)                     :: zero_first, zero_last, &
                                                last_c_blk
    INTEGER, DIMENSION(:), INTENT(OUT)       :: starts
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: i
    LOGICAL                                  :: dbg = .FALSE.

    parameters => group%master%s%parameters
    CALL dbcsr_psg_set_state (group, dbcsr_ps_state_filling, error)
    DO i = 1, SIZE (group%stacks)
       stack_p(i) = group%stacks(i)%s%stack_p
       starts(i) = group%stacks(i)%s%start_p
    END DO
    zero_first = group%master%s%t%t%zero_first
    zero_last  = group%master%s%t%t%zero_last
    last_c_blk = group%master%s%t%t%last_c_blk
    IF (dbg) &
       WRITE(*,*)"view open", stack_p, ";", zero_first, zero_last, last_c_blk, ";", starts
  END SUBROUTINE dbcsr_psg_view_open

! *****************************************************************************
!> \brief Close Direct access into stack group
!> \param[in] group               Parameter stack group
!> \param[in] parameters          Pointer to the master parameter stack
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_psg_view_close(group, parameters,&
       stack_p, zero_first, zero_last, last_c_blk, error)
    TYPE(dbcsr_ps_group_type), INTENT(INOUT) :: group
    INTEGER, DIMENSION(:), POINTER           :: parameters
    INTEGER, DIMENSION(:), INTENT(IN)        :: stack_p
    INTEGER, INTENT(IN)                      :: zero_first, zero_last, &
                                                last_c_blk
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: i
    LOGICAL                                  :: dbg = .FALSE.

    IF (dbg) &
         WRITE(*,*)"view clos", stack_p, ";", zero_first, zero_last, last_c_blk
    CALL dbcsr_psg_set_state (group, dbcsr_ps_state_filled, error)
    group%master%s%t%t%zero_first = zero_first
    group%master%s%t%t%zero_last  = zero_last
    group%master%s%t%t%last_c_blk = last_c_blk
    DO i = 1, SIZE (group%stacks)
       group%stacks(i)%s%stack_p = stack_p(i)
    ENDDO
  END SUBROUTINE dbcsr_psg_view_close

! *****************************************************************************
!> \brief Initializes a set of stacks
!> \param[in] group_size              The size of stack groups
!> \param[in] n_stack_mem_regions     Number of memory region in the stack
!> \param[in] stack_size              Each stack size
!> \param[in,out] error               error
! *****************************************************************************
  SUBROUTINE dbcsr_ps_set_init (group_size,n_stack_mem_regions,stack_size, error)
    INTEGER                                  :: group_size,n_stack_mem_regions, &
                                                stack_size
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER                                  :: ithread, nthreads, mem_type, k
#if defined (__DBCSR_CUDA)

    ithread = 0 ; nthreads = 1
    mem_type = dbcsr_memory_CUDA_host_pinned
    !$ ithread = OMP_GET_THREAD_NUM ()
    !$ nthreads = OMP_GET_NUM_THREADS ()

    !$OMP MASTER
     NULLIFY(ps_host_global)
     NULLIFY(ps_host_parameter)
     ALLOCATE (ps_host_global(nthreads))
     ALLOCATE (ps_host_parameter(nthreads))
    !$OMP END MASTER
    !$OMP BARRIER

    ALLOCATE (ps_host_global(ithread+1)%ps_local(n_stack_mem_regions))
    ALLOCATE (ps_host_parameter(ithread+1)%mem_stack(n_stack_mem_regions))

    DO k = 1, n_stack_mem_regions

          CALL memory_allocate (ps_host_global(ithread+1)%ps_local(k)%data,group_size,&
               mem_type=mem_type, error=error)
          CALL memory_allocate (ps_host_parameter(ithread+1)%mem_stack(k)%data,&
               stack_size*dbcsr_ps_width, mem_type=mem_type, error=error)
    ENDDO
#endif
  END SUBROUTINE dbcsr_ps_set_init

! *****************************************************************************
!> \brief Finalize a set of stacks
!> \param[in,out] error               error
! *****************************************************************************
  SUBROUTINE dbcsr_ps_set_finalize (error)
  TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

#if defined (__DBCSR_CUDA)
  INTEGER                                  :: k, i

     IF (ASSOCIATED(ps_host_global)) THEN
      DO k = 1, SIZE(ps_host_global)
       DO i = 1, SIZE(ps_host_global(k)%ps_local)
         CALL memory_deallocate(ps_host_global(k)%ps_local(i)%data,&
             dbcsr_memory_CUDA_host_pinned, error=error)
      ENDDO
         DEALLOCATE(ps_host_global(k)%ps_local)
      ENDDO
      DEALLOCATE(ps_host_global)
     ENDIF
     IF (ASSOCIATED(ps_host_parameter)) THEN
        DO k = 1, SIZE(ps_host_parameter)
         DO i = 1, SIZE(ps_host_parameter(k)%mem_stack)
           CALL memory_deallocate(ps_host_parameter(k)%mem_stack(i)%data,&
                            dbcsr_memory_CUDA_host_pinned, error=error)
        ENDDO
         DEALLOCATE(ps_host_parameter(k)%mem_stack)
        ENDDO
        DEALLOCATE(ps_host_parameter)
     ENDIF
#endif

  END SUBROUTINE dbcsr_ps_set_finalize

! *****************************************************************************
!> \brief Creates a set of stacks
!> \param[out] stack_set   Set of stacks
!> \param[in] group_size   The size of stack groups
!> \param[in] nbuffers     Number of buffers in memory region
!> \param[in] nmemregs     Number of memory regions.
!> \param[in,out] error    error
! *****************************************************************************
  SUBROUTINE dbcsr_ps_set_create (stack_set, stack_size, product_target,&
       group_size, nbuffers, nmemregs,&
       mem_type, error)
    TYPE(dbcsr_ps_set_type), INTENT(OUT)     :: stack_set
    INTEGER, INTENT(IN)                      :: stack_size
    TYPE(dbcsr_ps_target_obj), INTENT(IN)    :: product_target
    INTEGER, INTENT(IN)                      :: group_size, nbuffers, &
                                                nmemregs, mem_type
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_ps_set_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, i, ithread, j, k
    TYPE(dbcsr_ps_obj)                       :: memreg_stack

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set(routineN, error_handle, error)
    ithread=0
    !$ ithread = OMP_GET_THREAD_NUM ()
    ALLOCATE (stack_set%groups (nbuffers, nmemregs))
    stack_set%nbuffers = nbuffers
    stack_set%nmemregs = nmemregs
    stack_set%next_buffer = 1
    stack_set%next_memreg = 1
    stack_set%group_size = group_size
    ALLOCATE (stack_set%all_states (nmemregs))
    DO k = 1, stack_set%nmemregs

       ! Allocates the state array for all (non-master) stacks in this
       ! memory region.
      IF( mem_type .EQ. dbcsr_memory_CUDA_host_pinned) THEN
        stack_set%all_states(k)%data=>ps_host_global(ithread+1)%ps_local(k)%data
      ELSE
       CALL memory_allocate (stack_set%all_states(k)%data,&
            group_size * nbuffers,&
            mem_type=mem_type, error=error)
      ENDIF
       stack_set%all_states(k)%mem_type = mem_type

       DO j = 1, stack_set%nbuffers
          CALL dbcsr_ps_init (stack_set%groups(j,k)%master, error)
          ALLOCATE (stack_set%groups(j,k)%stacks(group_size))
          DO i = 1, SIZE (stack_set%groups(j,k)%stacks)
             CALL dbcsr_ps_init (stack_set%groups(j,k)%stacks(i), error)
          ENDDO
       ENDDO
       ! A new master stack for this memory region is created.  All
       ! (non-master) parameter arrays will be set to point into the
       ! parameter array of this master stack.

      CALL dbcsr_ps_init (memreg_stack, error)
      IF( mem_type .EQ. dbcsr_memory_CUDA_host_pinned) THEN
       ALLOCATE (memreg_stack%s)
       CALL ps_clear(memreg_stack%s,error=error)
       memreg_stack%s%parameters=>ps_host_parameter(ithread+1)%mem_stack(k)%data
       memreg_stack%s%mem_type = mem_type
       memreg_stack%s%own_data = .TRUE.
       memreg_stack%s%start_p = 1
       !$omp critical
       sid = sid+1
       memreg_stack%s%id = -1
       memreg_stack%s%sid = sid
    !$omp end critical
      memreg_stack%s%id = -1
      memreg_stack%s%sid = -7
      memreg_stack%s%refcount = 1
     ELSE
       CALL dbcsr_ps_new (memreg_stack,&
            stack_size*group_size*nbuffers,&
            mem_type=mem_type, error=error)
     ENDIF
       CALL dbcsr_ps_add_target (memreg_stack%s, product_target, error)
       CALL dbcsr_ps_set_set_b_master (stack_set, k,&
            memreg_stack, error=error)
       CALL dbcsr_ps_release (memreg_stack, error=error)
    ENDDO
    IF (careful_mod) CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_ps_set_create

! *****************************************************************************
!> \brief Releases all resources held by a stack set.
! *****************************************************************************
  SUBROUTINE dbcsr_ps_set_destroy (stack_set, error)
    TYPE(dbcsr_ps_set_type), INTENT(INOUT)   :: stack_set
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_ps_set_destroy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, i, j, k

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handle, error)
    DO k = 1, stack_set%nmemregs
       DO j = 1, stack_set%nbuffers
          DO i = 1, SIZE (stack_set%groups(j,k)%stacks)
             CALL dbcsr_ps_release (stack_set%groups(j,k)%stacks(i), error)
          ENDDO
          DEALLOCATE (stack_set%groups(j,k)%stacks)
          CALL dbcsr_ps_release (stack_set%groups(j,k)%master, error)
       ENDDO
       IF ( stack_set%all_states(k)%mem_type .NE. dbcsr_memory_CUDA_host_pinned ) THEN
       CALL memory_deallocate (stack_set%all_states(k)%data,&
                               stack_set%all_states(k)%mem_type, error=error)
       ENDIF
    ENDDO
    DEALLOCATE (stack_set%all_states)
    DEALLOCATE (stack_set%groups)
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_ps_set_destroy


! *****************************************************************************
!> \brief Sets the stack groups for all buffers in a given memory region.
!>
!> All actual group stacks share the same memory of the master_stack.
!> \param[in,out] stack_set   Set of stacks
!> \param[in] memreg          Memory region
!> \param[in] master_stack    Master stack for this memory region
!> \param[in,out] error    error
! *****************************************************************************
  SUBROUTINE dbcsr_ps_set_set_b_master (stack_set, memreg,&
       master_stack, error)
    TYPE(dbcsr_ps_set_type), INTENT(INOUT)   :: stack_set
    INTEGER, INTENT(IN)                      :: memreg
    TYPE(dbcsr_ps_obj), INTENT(IN)           :: master_stack
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: buff_array_len, buffer, &
                                                error_handle, group_size, i, &
                                                len, nbuffers, start
    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_ps_set_set_b_master', &
      routineP = moduleN//':'//routineN

    INTEGER, POINTER                         :: state_p

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set(routineN, error_handle, error)
    !
    nbuffers = stack_set%nbuffers
    buff_array_len = SIZE(master_stack%s%parameters)/dbcsr_ps_width / nbuffers
    !
    group_size = SIZE (stack_set%groups(1, memreg)%stacks)
    len = buff_array_len / group_size
    DO buffer = 1, nbuffers
       stack_set%groups(buffer, memreg)%master = master_stack
       CALL dbcsr_ps_hold (stack_set%groups(buffer, memreg)%master, error)
       DO i = 1, group_size
          start = (i-1)*len+1 + (buffer-1)*buff_array_len
          state_p => stack_set%all_states(memreg)%data(&
                        (buffer-1) * group_size + &
                        i-1 + 1 )
          CALL dbcsr_ps_new_p (stack_set%groups(buffer, memreg)%stacks(i),&
               master_stack, len, start,&
               state=state_p, &
               error=error)
          CALL dbcsr_ps_add_target(stack_set%groups(buffer, memreg)%stacks(i)%s,&
               master_stack%s%t, error)
       ENDDO
    ENDDO
    IF (careful_mod) CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_ps_set_set_b_master


! *****************************************************************************
!> \brief Returns a pointer to a stack group.
!>
!> If buffer or memreg are unspecified, their values are automatically
!> determined.
!> \param[in,out] stack_set   Set of stacks
!> \param[in] buffer          (optional) Return stack group of this buffer
!> \param[in] memreg          (optional) Return stack group from this memory
!>                            region
!> \param[in] master_array    returns the actual stack array
!> \param[in] wait            wait for stack group to be free; default is TRUE
!> \param[in,out] error    error
!> \par master_array
!>      The master array is the pointer to the master array associated
!>      with the stack group.  The starting elemnent of this group or
!>      stack into this master array is obtained by opening a
!>      view. \see dbcsr_psg_view_open
! *****************************************************************************
  FUNCTION dbcsr_ps_set_get_group_p (stack_set, buffer, memreg,&
       master_array, wait, error) RESULT (ps_group)
    TYPE(dbcsr_ps_set_type), INTENT(INOUT)   :: stack_set
    INTEGER, INTENT(IN), OPTIONAL            :: buffer, memreg
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: master_array
    LOGICAL, INTENT(IN), OPTIONAL            :: wait
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    TYPE(dbcsr_ps_group_type), POINTER       :: ps_group

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_ps_set_get_group_p', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a, b, error_handle, state
    LOGICAL                                  :: my_wait

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set(routineN, error_handle, error)
    IF (PRESENT (buffer)) THEN
       b = buffer
    ELSE
       b = stack_set%next_buffer
       !stack_set%next_buffer = MOD(stack_set%next_buffer-1, stack_set%nbuffers)+1
    ENDIF
    IF (PRESENT (memreg)) THEN
       a = memreg
    ELSE
       a = stack_set%next_memreg
       !stack_set%next_memreg = MOD(stack_set%next_memreg-1, stack_set%nmemregs)+1
    ENDIF
    IF (PRESENT (wait)) THEN
       my_wait = wait
    ELSE
       my_wait = .TRUE.
    ENDIF
    IF (my_wait) THEN
       ! Busy loop until processing of the requested group has finished.
       !$OMP FLUSH
       state = dbcsr_psg_get_state (stack_set%groups(b, a), error)
       DO WHILE (state .GE. dbcsr_ps_state_queued)
          !write(*,*)routineN//" Waiting", OMP_GET_THREAD_NUM()
          !$OMP FLUSH
          state = dbcsr_psg_get_state (stack_set%groups(b, a), error)
       ENDDO
    ENDIF
    ps_group => stack_set%groups(b, a)
    IF (PRESENT (master_array)) THEN
       master_array => ps_group%master%s%parameters
    ENDIF
    IF (careful_mod) CALL dbcsr_error_stop(error_handle, error)
  END FUNCTION dbcsr_ps_set_get_group_p

! *****************************************************************************
!> \brief Gets the number of stacks in the entire set of stacks being
!>        processed.
!>
!> This is used to balance work between an accelerator and the host CPU.
! *****************************************************************************
  FUNCTION dbcsr_ps_set_get_n_working (stack_set, error) RESULT (n_working)
    TYPE(dbcsr_ps_set_type), INTENT(INOUT)   :: stack_set
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER                                  :: n_working

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_ps_set_get_n_working', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: mem_region
    INTEGER, DIMENSION(:), POINTER           :: region_states

    n_working = 0
    !$OMP FLUSH
    DO mem_region = 1, SIZE(stack_set%all_states)
       region_states => stack_set%all_states(mem_region)%data
       n_working = n_working + &
            COUNT (region_states .EQ. dbcsr_ps_state_working)
    ENDDO
  END FUNCTION dbcsr_ps_set_get_n_working


  !> \brief Sets state for a whole process stack group
  SUBROUTINE dbcsr_psg_set_state (ps_group, state, error)
    TYPE(dbcsr_ps_group_type), INTENT(INOUT) :: ps_group
    INTEGER, INTENT(IN)                      :: state
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: i

    ps_group%master%s%state = state
    DO i = 1, SIZE(ps_group%stacks)
       ps_group%stacks(i)%s%state = state
    ENDDO
  END SUBROUTINE dbcsr_psg_set_state

  !> \brief Updates the state the process stack group
  SUBROUTINE dbcsr_psg_update_state (ps_group, error)
    TYPE(dbcsr_ps_group_type), INTENT(INOUT) :: ps_group
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: i, state

    state = -1
    DO i = 1, SIZE(ps_group%stacks)
       state = MAX(state, ps_group%stacks(i)%s%state)
    ENDDO
    ps_group%master%s%state = state
  END SUBROUTINE dbcsr_psg_update_state

  !> \brief Returns the state of the process stack group
  FUNCTION dbcsr_psg_get_state (ps_group, error) RESULT (state)
    TYPE(dbcsr_ps_group_type), INTENT(INOUT) :: ps_group
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER                                  :: state

    CALL dbcsr_psg_update_state (ps_group, error)
    state = ps_group%master%s%state
  END FUNCTION dbcsr_psg_get_state



! *****************************************************************************
!> \brief Advances the counters for automatic group fetching.
!> \see dbcsr-ps_set_get_group_p
!> \param[in,out] stack_set   Set of stacks
!> \param[in] buffer          (optional) Return stack group of this buffer
!> \param[in] memreg          (optional) Return stack group from this memory
!>                            region
!> \param[in] master_array    returns the actual stack array
!> \param[in,out] error    error
!> \par master_array
!>      The master array is the pointer to the master array associated
!>      with the stack group.  The starting elemnent of this group or
!>      stack into this master array is obtained by opening a
!>      view. \see dbcsr_psg_view_open
! *****************************************************************************
  SUBROUTINE dbcsr_ps_set_advance (stack_set, advance_memreg, error)
    TYPE(dbcsr_ps_set_type), INTENT(INOUT)   :: stack_set
    LOGICAL, INTENT(IN)                      :: advance_memreg
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_ps_set_advance', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: last_mr

! First advance group; if it wraps then advance the buffer

    last_mr = stack_set%next_memreg
    IF (advance_memreg) THEN
         stack_set%next_memreg = MOD(stack_set%next_memreg, stack_set%nmemregs)+1
      ELSE
         stack_set%next_buffer = MOD(stack_set%next_buffer, stack_set%nbuffers)+1
      ENDIF
  END SUBROUTINE dbcsr_ps_set_advance

! *****************************************************************************
!> \brief Locks the target
! *****************************************************************************
  SUBROUTINE dbcsr_ps_target_lock_main (t, flush, error)
    TYPE(dbcsr_ps_target_type), POINTER      :: t
    LOGICAL, INTENT(IN), OPTIONAL            :: flush
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    LOGICAL, PARAMETER                       :: dbg = .FALSE., &
                                                owner_priority = .TRUE.

    IF (owner_priority) THEN
       !$ IF (t%owner .EQ. OMP_GET_THREAD_NUM()) THEN
       !$    t%owner_wants_lock = .TRUE.
       !$OMP FLUSH
       !$ ELSE
       !$    DO WHILE (t%owner_wants_lock)
               !$OMP FLUSH
       !$    ENDDO
       !$ ENDIF
    ENDIF
    !$ CALL OMP_SET_LOCK (t%target_lock)
    !$ IF (mm_thread_workshare .AND. .NOT. has_omp_2_5) THEN
       !$OMP FLUSH
    !$ ENDIF
    !$ if (dbg) write(*,*)"  locked main", OMP_GET_THREAD_NUM(), t%owner
    IF (owner_priority) THEN
       !$ IF (t%owner .EQ. OMP_GET_THREAD_NUM()) THEN
       !$    t%owner_wants_lock = .FALSE.
             !$OMP FLUSH
       !$ ENDIF
       !$ IF (PRESENT (FLUSH)) THEN
       !$    IF (FLUSH) THEN
                !$OMP FLUSH
       !$       DO WHILE (t%n_use .GT. 0)
       !$          if (dbg) WRITE(*,*)" waiting", OMP_GET_THREAD_NUM(), t%n_use
                   !$OMP FLUSH
       !$       ENDDO
       !$    ENDIF
       !$ ENDIF
    ENDIF
  END SUBROUTINE dbcsr_ps_target_lock_main


! *****************************************************************************
!> \brief Unlocks the target
! *****************************************************************************
  SUBROUTINE dbcsr_ps_target_unlock_main (t, error)
    TYPE(dbcsr_ps_target_type), POINTER      :: t
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    !$ IF (mm_thread_workshare .AND. .NOT. has_omp_2_5) THEN
       !$OMP FLUSH
    !$ ENDIF
    !$ CALL OMP_UNSET_LOCK (t%target_lock)
    !$ if (dbg) write(*,*)"unlocked main", OMP_GET_THREAD_NUM(), t%owner
  END SUBROUTINE dbcsr_ps_target_unlock_main


! *****************************************************************************
!> \brief Locks target areas used for this stack processing
!> \param[in,out] target_desc     Unlocks this target
!> \param[in] stack_size          Size of stack
!> \param[in] parameters          Stack of paramaters
!> \param[in] bit_shift           Bit shift used to map target matrix data
!>                                offsets to regions
!> \param[in,out] which_locks     Locks taken
!> \param[in] nregions            Number of regions
!> \param]in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_ps_target_lock_regions(target_desc, stack_size, parameters,&
       bit_shift, which_locks, nregions, error)
    TYPE(dbcsr_ps_target_type), POINTER      :: target_desc
    INTEGER, INTENT(IN)                      :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width, &
      stack_size), INTENT(IN)                :: parameters
    INTEGER, INTENT(IN)                      :: bit_shift
    LOGICAL, DIMENSION(:), INTENT(INOUT)     :: which_locks
    INTEGER, INTENT(OUT)                     :: nregions
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_ps_target_lock_regions', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: error_handle

    !$ INTEGER                               :: sp, sft

    IF (careful_mod) &
         CALL dbcsr_error_set(routineN, error_handle, error)
    !$ IF (dbg .AND. OMP_GET_THREAD_NUM () .NE. target_desc%owner) THEN
    !$    write(*,*)"  lckng rgn  ", OMP_GET_THREAD_NUM(), target_desc%owner
    !$ ENDIF

    !$ CALL dbcsr_ps_target_lock_main (target_desc, error=error)

    IF (.NOT. full_target_lock) THEN
      !$ nregions = SIZE(target_desc%region_locks)
         !$OMP ATOMIC
      !$ target_desc%n_use = target_desc%n_use + 1
         !$OMP FLUSH
      !$ IF (dbg) WRITE(*,*)routineN//"    added one", OMP_GET_THREAD_NUM(),&
      !$      target_desc%n_use, target_desc%owner
      !$ which_locks(:) = .FALSE.
      !$ sft = -bit_shift
      !$ DO sp = 1, stack_size
      !$    which_locks(1+ISHFT(parameters(p_c_first,sp), sft)) = .TRUE.
      !$ ENDDO
      !$ DO sp = 1, nregions
      !$    IF (which_locks(sp)) THEN
      !$       CALL OMP_SET_LOCK (target_desc%region_locks(sp))
      !$    ENDIF
      !$ ENDDO
      !$ IF (dbg) WRITE(*,*)OMP_GET_THREAD_NUM(), which_locks
      !$ CALL dbcsr_ps_target_unlock_main (target_desc, error=error)
    ENDIF
    IF (careful_mod) &
         CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_ps_target_lock_regions


! *****************************************************************************
!> \brief Unlocks target areas used for this stack processing
!> \param[in,out] target_desc     Unlocks this target
!> \param[in,out] which_locks     Locks to unlock--those that I have held
!> \param[in] nregions            Number of regions
!> \param]in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_ps_target_unlock_regions(target_desc, &
       which_locks, nregions, error)
    TYPE(dbcsr_ps_target_type), POINTER      :: target_desc
    LOGICAL, DIMENSION(:), INTENT(INOUT)     :: which_locks
    INTEGER, INTENT(IN)                      :: nregions
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_ps_target_unlock_regions', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: error_handle

    !$ INTEGER                               :: sp

    IF (careful_mod) &
         CALL dbcsr_error_set(routineN, error_handle, error)

    IF (full_target_lock) THEN
       !$   CALL dbcsr_ps_target_unlock_main (target_desc, error=error)
    ELSE
      !$ DO sp = 1, nregions
      !$    IF (which_locks(sp)) THEN
      !$       CALL OMP_UNSET_LOCK (target_desc%region_locks(sp))
      !$    ENDIF
      !$ ENDDO
         !$OMP FLUSH
         !$OMP ATOMIC
      !$ target_desc%n_use = target_desc%n_use - 1
      !$ IF (dbg) WRITE(*,*)routineN//" subbed one", OMP_GET_THREAD_NUM(),&
      !$      target_desc%n_use, target_desc%owner
      !$ IF (dbg) THEN
      !$   write(*,*)OMP_GET_THREAD_NUM(), which_locks
      !$   DO sp = 1, nregions
      !$      which_locks(sp) =  OMP_TEST_LOCK (target_desc%region_locks(sp))
      !$   ENDDO
      !$   WRITE(*,*)OMP_GET_THREAD_NUM(), which_locks
      !$   DO sp = 1, nregions
      !$     IF (which_locks(sp)) &
      !$          CALL OMP_UNSET_LOCK (target_desc%region_locks(sp))
      !$   ENDDO
      !$ ENDIF
    ENDIF
    IF (careful_mod) &
         CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_ps_target_unlock_regions


! *****************************************************************************
!> \brief Resizes the set of region locks used for the target region
!> \param[in,out] t        Resize locks of this target
!> \param[in] max_size     Size of the data area
!> \param]in,out] error    error
! *****************************************************************************
  SUBROUTINE dbcsr_ps_target_new_regions (t, max_size, error)
    TYPE(dbcsr_ps_target_type), &
      INTENT(INOUT)                          :: t
    INTEGER, INTENT(IN)                      :: max_size
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_ps_target_new_regions', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = debug_mod

    INTEGER                                  :: error_handle, i, new_n_regions

    IF (careful_mod .OR. dbg) &
         CALL dbcsr_error_set(routineN, error_handle, error)
    !!$ new_n_regions = MAX(64, NINT (LOG(REAL(max_size))/LOG(2.0)))
    !$ new_n_regions = MIN(max_regions,&
    !$                     MAX(1024,&
    !$                         NINT (LOG(REAL(max_size))/LOG(2.0))))
    !$ t%bit_shift = CEILING (LOG(REAL(max_size)/REAL(new_n_regions))/LOG(2.0))
    !$ IF (new_n_regions .GT. SIZE(t%region_locks)) THEN
    !$    IF (dbg) WRITE(*,*)"new_regions", max_size, new_n_regions, t%bit_shift
    !$    DO i = 1, SIZE(t%region_locks)
    !$       IF (dbg) THEN
    !$          IF (.NOT. OMP_TEST_LOCK (t%region_locks(i)))&
    !$               WRITE(*,*)'oops! lock is set'
    !$          CALL OMP_UNSET_LOCK (t%region_locks(i))
    !$       ENDIF
    !$       CALL OMP_DESTROY_LOCK (t%region_locks(i))
    !$    ENDDO
    !$    DEALLOCATE (t%region_locks)
    !$    ALLOCATE (t%region_locks(new_n_regions))
    !$    DO i = 1, SIZE(t%region_locks)
    !  !$       t%region_locks(i) = 0
    !$       CALL OMP_INIT_LOCK (t%region_locks(i))
    !$    ENDDO
    !$ ENDIF
    IF (careful_mod .OR. dbg) &
         CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_ps_target_new_regions


END MODULE dbcsr_pq_methods
