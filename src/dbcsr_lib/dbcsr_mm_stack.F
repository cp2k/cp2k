!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Stacks of small matrix multiplications
!> \author  Urban Borstnik
!> \date    2011-09-26
!> \version 0.9
!>
!> <b>Modification history:</b>
!  - 2011-09-26 Split dbcsr_internal_operations
! *****************************************************************************
MODULE dbcsr_mm_stack
  USE dbcsr_block_operations,          ONLY: dbcsr_data_clear
  USE dbcsr_config,                    ONLY: &
       detailed_timing, driver_is_async, kernel_timing, mm_async, mm_driver, &
       mm_driver_blas, mm_driver_cuda, mm_driver_matmul, mm_driver_plasma, &
       mm_driver_smm, mm_workshare
  USE dbcsr_cuda_device,               ONLY: dbcsr_cuda_device_sync,&
                                             dbcsr_cuda_stream_sync
  USE dbcsr_cuda_memory,               ONLY: dbcsr_cuda_dev_mem_alloc,&
                                             dbcsr_cuda_dev_mem_dealloc,&
                                             dbcsr_cuda_dev_mem_realloc,&
                                             dbcsr_cuda_dev_mem_zero
  USE dbcsr_cuda_methods,              ONLY: dbcsr_cuda_dev_mem_get_alloc
  USE dbcsr_cuda_operations,           ONLY: dbcsr_cuda_cp_dev_to_host,&
                                             dbcsr_cuda_cp_host_to_dev,&
                                             dbcsr_cuda_do_mm_stack
  USE dbcsr_cuda_types,                ONLY: dbcsr_cuda_mem_type
  USE dbcsr_data_methods,              ONLY: dbcsr_data_ensure_size,&
                                             dbcsr_data_get_size,&
                                             dbcsr_data_set_size_referenced
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: dp,&
                                             real_4,&
                                             real_8,&
                                             sp
  USE dbcsr_machine,                   ONLY: m_flush,&
                                             m_walltime
  USE dbcsr_pq_methods,                ONLY: &
       dbcsr_pq_add_stack, dbcsr_pq_flush_level_chg, dbcsr_pq_get_any_stack, &
       dbcsr_pq_get_far_stack, dbcsr_pq_get_own_stack, &
       dbcsr_pq_handoff_level_chg, dbcsr_ps_set_get_group_p, &
       dbcsr_ps_set_get_n_working, dbcsr_ps_target_lock_main, &
       dbcsr_ps_target_lock_regions, dbcsr_ps_target_new_regions, &
       dbcsr_ps_target_unlock_main, dbcsr_ps_target_unlock_regions, &
       dbcsr_psg_get_state, dbcsr_psg_set_state
  USE dbcsr_pq_types,                  ONLY: &
       dbcsr_pq_type, dbcsr_ps_group_type, dbcsr_ps_obj, dbcsr_ps_set_type, &
       dbcsr_ps_state_empty, dbcsr_ps_state_queued, dbcsr_ps_state_working, &
       dbcsr_ps_target_obj, dbcsr_ps_target_type, dbcsr_ps_type, &
       dbcsr_ps_width, max_regions, mult_timers, p_a_first, p_b_first, &
       p_c_first, p_k, p_m, p_n, tmr_get_a_stack, tmr_get_other_stack, &
       tmr_kernel_time, tmr_process_me, tmr_process_other_stack, &
       tmr_process_stack, tmr_process_target, tmr_wait_for_empty_stack, &
       tmr_wait_for_my_target, tmr_wait_for_other_target
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8
  USE dbcsr_util,                      ONLY: xtime_start,&
                                             xtime_stop

  !$ USE OMP_LIB


  IMPLICIT NONE

  PRIVATE


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_stack'

  CHARACTER(len=*), PARAMETER, PRIVATE :: int_print = "(10(1X,I7))"

  REAL, PARAMETER                      :: default_resize_factor = 1.618034

  INTEGER :: multrec_calls = 0

  LOGICAL, PARAMETER, PUBLIC :: show_levels = .FALSE.


  PUBLIC :: enqueue_ps_group
  PUBLIC :: process_queue_mine,&
            process_queue_others,&
            process_queue_preempt,&
            get_stack_or_process_queue


  LOGICAL, PARAMETER :: debug_mod  = .FALSE.
  LOGICAL, PARAMETER :: careful_mod = .FALSE.

  !> \var max_stack_block_size  The maximal block size to be specially
  !>                            treated.
  INTEGER, PARAMETER :: max_stack_block_size = HUGE (INT (0))

  LOGICAL, PARAMETER, PRIVATE :: verbose_acc = .FALSE.

  REAL(kind=dp), PRIVATE :: index_time
  LOGICAL, PRIVATE :: do_index_time = .FALSE.
  LOGICAL, PRIVATE :: print_index_time = .FALSE.
  LOGICAL, PRIVATE :: measure_idle = .FALSE.

  REAL(KIND=dp), PUBLIC :: t_calc_step, t_dev_idle, t_process_stack,&
       t_dev_sync

  !$ INTEGER(KIND=omp_lock_kind), SAVE :: accel_lock

CONTAINS

! *****************************************************************************
!> \brief Puts stacks ready to be processed into a queue
!>
!> All non-empty stacks in the param_group are added to the queue.
!> \par Non-worksharing
!>      When mm_workshare is FALSE the stacks are instead processessed
!>      immediately, bypassing the queue.
!> \param[in,out] queue        My queue, to which parameter stacks are added.
!> \param[in] param_group      Parameter stacks to add to queue
!> \param[in,out] error        error
!> \note The queue and param_groups parameters could be INTENT(INOUT)
!>       when OpenMP is not used.
! *****************************************************************************
  SUBROUTINE enqueue_ps_group (queue, param_group, error)
    TYPE(dbcsr_pq_type), POINTER             :: queue
    TYPE(dbcsr_ps_group_type), POINTER       :: param_group
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'enqueue_ps_group', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: careful = careful_mod, &
                                                dbg = debug_mod

    INTEGER                                  :: error_handle, i

    IF (careful_mod) &
         CALL dbcsr_error_set(routineN, error_handle, error)
    CALL dbcsr_psg_set_state (param_group, dbcsr_ps_state_queued, error)
    param_group%master%s%driver = mm_driver
    !$ if (dbg) write(*,*)"  lckng targt", OMP_GET_THREAD_NUM(),&
    !$      param_group%master%s%t%t%owner
    CALL xtime_start (mult_timers(tmr_process_target))
    CALL process_ps_target_low (param_group%master%s%t,&
         param_group%master%s%t%t%product_data_area,&
         param_group%master%s%t%t%zero_first,&
         param_group%master%s%t%t%zero_last,&
         param_group%master%s%t%t%last_c_blk,&
         param_group%master%s%driver,&
         param_group%master%s%t%t%product_data_cuda,&
         param_group%master%s%t%t%c_locks_dev,&
         error)
    CALL xtime_stop (mult_timers(tmr_process_target))
    DO i = SIZE(param_group%stacks), 1, -1
       param_group%stacks(i)%s%driver = mm_driver
       IF (param_group%stacks(i)%s%stack_p .EQ. 0) THEN
           param_group%stacks(i)%s%state = dbcsr_ps_state_empty
          IF (dbg) WRITE(*,*)routineN//" not adding empty stack"
       ELSE
          IF (dbg) WRITE(*,*)routineN//" adding stack len", param_group%stacks(i)%s%stack_p
          IF (mm_workshare .OR. mm_async) THEN
             CALL dbcsr_pq_add_stack (queue, param_group%stacks(i), error=error)
          ELSE
             ! If no worksharing is active, just process the stack immediately.
             param_group%stacks(i)%s%driver = mm_driver
             CALL xtime_start (mult_timers(tmr_process_stack))
             CALL process_ps_stack (param_group%stacks(i)%s, error=error)
             CALL xtime_stop (mult_timers(tmr_process_stack))
          ENDIF
       ENDIF
    ENDDO
    IF (careful_mod) &
         CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE enqueue_ps_group


! *****************************************************************************
!> \brief Processes stacks from my queue
!>
!> \param[in,out] queue         My queue
!> \param[in] driver            Driver to use for processing queue
!> \param[in,out] error         error
!> \note The queue parameter could be INTENT(INOUT) when OpenMP is not used.
! *****************************************************************************
  SUBROUTINE process_queue_mine (queue, driver, error)
    TYPE(dbcsr_pq_type), POINTER             :: queue
    INTEGER, INTENT(IN)                      :: driver
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'process_queue_mine', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: careful = careful_mod, &
                                                dbg = debug_mod

    INTEGER                                  :: cnt, error_handle, max_handoff
    LOGICAL                                  :: easy, found
    TYPE(dbcsr_ps_obj)                       :: stack

    IF (careful_mod) &
         CALL dbcsr_error_set(routineN, error_handle, error)
    easy = driver_is_async (driver)
    cnt = 0
    IF (easy) THEN
          max_handoff = queue%handoff - queue%nworking
    ENDIF
    CALL xtime_start (mult_timers(tmr_get_a_stack))
    DO
       CALL dbcsr_pq_get_own_stack (queue, stack, found,&
            easy = easy, error=error)
       IF (dbg) THEN
          !$omp master
          WRITE(*,*)"Found my stack"
          !$omp end master
       ENDIF
       IF (found) THEN
          stack%s%driver = driver
          CALL xtime_stop (mult_timers(tmr_get_a_stack))
          CALL xtime_start (mult_timers(tmr_process_stack))
          CALL process_ps_stack(stack%s, error=error)
          CALL xtime_stop (mult_timers(tmr_process_stack))
          CALL xtime_start (mult_timers(tmr_get_a_stack))
          cnt = cnt + 1
          IF (easy .AND. cnt .GT. max_handoff) EXIT
       ELSE
          EXIT
       ENDIF
    ENDDO
    CALL xtime_stop (mult_timers(tmr_get_a_stack))
    !
    IF (careful_mod) &
         CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE process_queue_mine


! *****************************************************************************
!> \brief Processes stacks from other threads
!>
!> \param[in,out] queue         My queue
!> \param[in] forever           Processes until all other threads declare
!>                              they are done.
!> \param[in] driver            Driver to use for processing queue
!> \param[in,out] error         error
!> \note The queue parameter could be INTENT(INOUT) when OpenMP is not used.
! *****************************************************************************
  SUBROUTINE process_queue_others (queue, forever, driver, error)
    TYPE(dbcsr_pq_type), POINTER             :: queue
    LOGICAL, INTENT(IN)                      :: forever
    INTEGER, INTENT(IN)                      :: driver
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'process_queue_others', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: careful = careful_mod, &
                                                dbg = debug_mod

    INTEGER                                  :: cnt, error_handle, nt, &
                                                old_thread_ptr, t, thread_ptr
    LOGICAL                                  :: easy, found, keep_cycling
    TYPE(dbcsr_ps_obj)                       :: stack

    IF (careful_mod) &
         CALL dbcsr_error_set(routineN, error_handle, error)
    !
    !
    nt = 1
    !$nt = omp_get_num_threads()
    easy = driver_is_async (driver)
    cnt = 1
    old_thread_ptr = 0
    thread_ptr = 1
    found = .TRUE.
    keep_cycling = found .OR.&
         (forever .AND. queue%all_queues%n_working .GT. 0)
    DO WHILE (keep_cycling)
       CALL xtime_start (mult_timers(tmr_get_other_stack))
       CALL dbcsr_pq_get_any_stack (queue, stack, found,&
            easy = easy,&
            anyways = easy,&
            thread_ptr = thread_ptr, error=error)
       CALL xtime_stop (mult_timers(tmr_get_other_stack))
       IF (thread_ptr .NE. old_thread_ptr) cnt = 1
       old_thread_ptr = thread_ptr
       IF (found) THEN
          IF (dbg) THEN
             !$omp master
             WRITE(*,*)"Found other stack owned by", stack%s%t%t%owner
             !$omp end master
          ENDIF
          stack%s%driver = driver
          CALL xtime_start (mult_timers(tmr_process_other_stack))
          CALL process_ps_stack(stack%s, error=error)
          CALL xtime_stop (mult_timers(tmr_process_other_stack))
          cnt = cnt + 1
          t = 0
          !$ t = omp_get_thread_num ()
          t = MOD(thread_ptr + t, nt)+1
          IF (easy .AND.&
               cnt .GT. queue%all_queues%queues(t)%pq%handoff&
                     - queue%all_queues%queues(t)%pq%nworking) THEN
             thread_ptr = thread_ptr + 1
             WRITE(*,*)'new thread', thread_ptr
          ENDIF
       ENDIF
       !$OMP FLUSH
       keep_cycling = found .OR.&
            (forever .AND. queue%all_queues%n_working .GT. 0)
    ENDDO
    !
    IF (careful_mod) &
         CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE process_queue_others


! *****************************************************************************
!> \brief Processes stacks while the accelerator processor is busy.
!>
!> Processing is not done if number of stacks being processed is less
!> than the queue's flush_level.
!> \see get_stack_or_process_queue
!> \see csr_multiply
!> \param[in,out] queue         My queue
!> \param[in,out] param_sets    Parameter group sets
!> \param[in] driver            Driver to use for processing queue
!> \param[in,out] error         error
!> \note The queue parameter could be INTENT(INOUT) when OpenMP is not used.
! *****************************************************************************
  SUBROUTINE process_queue_preempt (queue, param_set, driver, error)
    TYPE(dbcsr_pq_type), POINTER             :: queue
    TYPE(dbcsr_ps_set_type), INTENT(INOUT)   :: param_set
    INTEGER, INTENT(IN)                      :: driver
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'process_queue_preempt', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: careful = careful_mod, &
                                                dbg = debug_mod

    INTEGER                                  :: error_handle, iam, n
    LOGICAL                                  :: found
    TYPE(dbcsr_ps_obj)                       :: stack

    IF (careful_mod) &
         CALL dbcsr_error_set(routineN, error_handle, error)
    !
    ! Fetch & process only if no one else is done (i.e., has spare
    ! time to process my stacks).
    found = .TRUE.
    n = dbcsr_ps_set_get_n_working (param_set, error)
    !$ if (.FALSE. .AND. dbg) &
    !$     write(*,*)routineN, omp_get_thread_num(), n, queue%flush_level
    CALL xtime_start (mult_timers(tmr_get_a_stack))
    DO WHILE (n .GT. queue%flush_level .AND. found)
       CALL dbcsr_pq_get_far_stack (queue, stack, found, error=error)
       IF (found) THEN
          IF (dbg) THEN
             iam = 0
             !$ iam = omp_get_thread_num ()
             WRITE(*,*)"Preempt: Found own stack", iam
          ENDIF
          CALL xtime_stop (mult_timers(tmr_get_a_stack))
          CALL xtime_start (mult_timers(tmr_process_stack))
          stack%s%driver = driver
          CALL process_ps_stack(stack%s, error=error)
          CALL xtime_stop (mult_timers(tmr_process_stack))
          CALL xtime_start (mult_timers(tmr_get_a_stack))
       ENDIF
       IF (dbg) THEN
          n = dbcsr_ps_set_get_n_working (param_set, error)
          IF (n .EQ. 0) WRITE(*,*)"oops"
       ENDIF
    ENDDO
    CALL xtime_stop (mult_timers(tmr_get_a_stack))
    !
    IF (careful_mod) &
         CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE process_queue_preempt


! *****************************************************************************
!> \brief Gets a new stack group prepared for writing.
!> \par Dealing with busy stack groups.
!> The queue will be processed while the stack group is busy (it is
!> queued or in processing).
!>
!> If the stack_group is not immediately available the queue's
!> flush_level is lowered.  It thus becomes more favorable to have the
!> host CPU process stacks.
!> \see process_queue_preempt
!> \see csr_multiply
!> \param[in,out] param_sets    The parameter group sets
!> \param[in,out] queue         My queue
!> \param[out] stack_group      Retrieved stack group
!> \param[in] driver            Driver to use for processing queue
!> \param[in,out] error         error
!> \note The queue parameter could be INTENT(INOUT) when OpenMP is not used.
! *****************************************************************************
  SUBROUTINE get_stack_or_process_queue (param_sets, queue, stack_group,&
       driver, error)
    TYPE(dbcsr_ps_set_type), INTENT(INOUT)   :: param_sets
    TYPE(dbcsr_pq_type), POINTER             :: queue
    TYPE(dbcsr_ps_group_type), POINTER       :: stack_group
    INTEGER, INTENT(IN)                      :: driver
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    LOGICAL, PARAMETER                       :: dbg = debug_mod

    INTEGER                                  :: cnt, n, stack_state
    LOGICAL                                  :: found
    TYPE(dbcsr_ps_obj)                       :: stack

    CALL xtime_start (mult_timers(tmr_wait_for_empty_stack))
    stack_group => dbcsr_ps_set_get_group_p (param_sets, wait=.FALSE.,&
         error=error)
    IF (mm_workshare) THEN
       !$OMP FLUSH
    ENDIF
    ! If the parameter stack we get is still busy then
    ! just do more crunching.
    cnt = 1
    stack_state = dbcsr_psg_get_state (stack_group, error)
    IF (stack_state .GE. dbcsr_ps_state_queued .AND. mm_async) THEN
       ! The parameter stack is still busy.  The others processors are
       ! too slow so we lower the threshold for local computation.
       !
       !$ if (dbg) write(*,*)"still in queue:", queue%nstacks, omp_get_thread_num()
       n = dbcsr_ps_set_get_n_working (param_sets, error)
       CALL dbcsr_pq_flush_level_chg (queue, -(n/4),&
            param_sets)
       CALL dbcsr_pq_handoff_level_chg (queue, -MAX((param_sets%group_size*3)/2,1),&
            param_sets)
    ENDIF
    CALL xtime_start (mult_timers(tmr_get_a_stack))
    DO WHILE (stack_state .GE. dbcsr_ps_state_queued)
       ! Just process the queue until the parameter stack is clear.
       CALL dbcsr_pq_get_far_stack (queue, stack, found, error=error)
       IF (found) THEN
          CALL xtime_stop (mult_timers(tmr_get_a_stack))
          CALL xtime_stop (mult_timers(tmr_wait_for_empty_stack))
          CALL xtime_start (mult_timers(tmr_process_me))
          stack%s%driver = driver
          CALL process_ps_stack(stack%s, error=error)
          CALL xtime_stop (mult_timers(tmr_process_me))
          CALL xtime_start (mult_timers(tmr_wait_for_empty_stack))
          CALL xtime_start (mult_timers(tmr_get_a_stack))
       ENDIF
       IF (mm_workshare) THEN
          !$OMP FLUSH
       ENDIF
       cnt = cnt + 1
       stack_state = dbcsr_psg_get_state (stack_group, error)
    ENDDO
    CALL xtime_stop (mult_timers(tmr_get_a_stack))
    CALL xtime_stop (mult_timers(tmr_wait_for_empty_stack))
  END SUBROUTINE get_stack_or_process_queue


! *****************************************************************************
!> \brief Prepares the target for stack processing.
!> \note Resizes data area, zeros new data, resets device locks, ...
! *****************************************************************************
  SUBROUTINE process_ps_target_low (target_desc,&
       product_data_area,&
       zero_first, zero_last, nblks, driver, card_data, c_locks_dev, error)
    TYPE(dbcsr_ps_target_obj), INTENT(INOUT) :: target_desc
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: product_data_area
    INTEGER, INTENT(INOUT)                   :: zero_first, zero_last
    INTEGER, INTENT(IN)                      :: nblks, driver
    TYPE(dbcsr_cuda_mem_type), POINTER       :: card_data, c_locks_dev
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'process_ps_target_low', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: careful = careful_mod, &
                                                dbg = debug_mod, &
                                                info = .FALSE.
    REAL, PARAMETER                          :: resize_factor = 1.618034

    INTEGER                                  :: c_size, istat, maxs, tmp_i

!   ---------------------------------------------------------------------------
!
! Resize target data area if necessary.
! Here we want the actual allocation size.

    maxs = dbcsr_data_get_size(product_data_area)
    IF (zero_last .GT. maxs) THEN
       ! Resize the target data area if the stack references a target
       ! matrix data element outside of its current range.  When
       ! appropriate a lock on the data area is taken.
       IF (dbg) WRITE(*,*)"Target too small"
       maxs = zero_last
       IF (mm_workshare) THEN
          CALL xtime_start (mult_timers(tmr_wait_for_my_target))
          ! Lock the target and wait for others to unlock their
          ! regions.
          CALL dbcsr_ps_target_lock_main (target_desc%t,&
               flush=.TRUE.,&
               error=error)
          CALL xtime_stop (mult_timers(tmr_wait_for_my_target))
       ENDIF
       IF (dbg) &
            WRITE(*,*)routineN//" Resizing to", LOG(REAL(maxs))/LOG(10.0)
       CALL dbcsr_data_ensure_size (product_data_area,&
            maxs, factor=resize_factor, error=error)
       IF (mm_workshare) THEN
          ! Potentially resize locks.  The number of region locks is
          ! adjusted to the allocated size of the new data area (to
          ! prevent changing them whenever the last addressed block is
          ! changed).
          maxs = dbcsr_data_get_size(product_data_area)
          CALL dbcsr_ps_target_new_regions (target_desc%t,&
               maxs, error=error)
          CALL dbcsr_ps_target_unlock_main (target_desc%t, error)
       ENDIF
    ENDIF
    IF (dbg) WRITE(*,*)"Setting referenced size"
    CALL dbcsr_data_set_size_referenced (product_data_area, zero_last)
    IF (dbg) WRITE(*,*)"Done"
    !
    ! Zero new blocks
    IF (zero_last .GE. zero_first) THEN
       IF (dbg) WRITE(*,*)"Zeroing"
       CALL dbcsr_data_clear (product_data_area, lb=zero_first, ub=zero_last)
       IF (dbg) WRITE(*,*)"Done"
    ENDIF
    !
    ! Cuda on-device resizing
    cuda_process_target: IF (driver .EQ. mm_driver_cuda) THEN
       CALL dbcsr_cuda_stream_sync (target_desc%t%owner+1, error=error)
       IF (dbg) WRITE(*,*)"Getting allocated size..."
       !$ IF (dbg) WRITE(*,*)"thread, target", omp_get_thread_num(),&
       !$      target_desc%t%owner
       c_size = dbcsr_cuda_dev_mem_get_alloc(card_data)
       IF (dbg) WRITE(*,*)"done", c_size, "vs", zero_last
       ! Resize & zero product data if too big.
       IF (zero_last .GT. c_size) THEN
          IF (dbg .OR. info) WRITE(*,*)routineN//" reallocating c_dev",&
               c_size, zero_last
          IF (detailed_timing) THEN
             t_dev_sync = t_dev_sync - m_walltime()
             CALL dbcsr_cuda_device_sync(error=error)
             t_dev_sync = t_dev_sync + m_walltime()
          ENDIF
          tmp_i = MAX(zero_last, INT(REAL(zero_last)*resize_factor) )
          IF (verbose_acc) WRITE(*,*)routineN//" reallocating c_dev",&
               tmp_i, c_size
          !--- !$OMP CRITICAL (crit_cuda)
          !$ IF (mm_workshare) CALL OMP_SET_LOCK (accel_lock)
          CALL dbcsr_cuda_dev_mem_realloc (card_data, tmp_i, stat=istat)
          IF (istat /= 0) THEN
             ! If the resize, including pre-oversizing failed, try to
             ! resize it to the actual number of locks used.
             IF (verbose_acc) WRITE(*,*)routineN//" Running out of memory"
             CALL dbcsr_cuda_dev_mem_realloc (card_data, zero_last, stat=istat)
          ENDIF
          c_size = dbcsr_cuda_dev_mem_get_alloc(card_data)
          IF (dbg) WRITE(*,*)routineN//" zeroing c_dev"
          CALL dbcsr_cuda_dev_mem_zero(card_data,&
               first=zero_first, last=c_size, error=error)
          !$ IF (mm_workshare) CALL OMP_UNSET_LOCK (accel_lock)
          !--- !$OMP END CRITICAL (crit_cuda)
       ELSE IF (dbg) THEN
          WRITE(*,*)"not reallocating c_dev"
       ENDIF
       !
       ! Resize locks, which are equal to the new block count.
       IF (dbg) WRITE(*,*)routineN//" lock size..."
       IF (dbg) WRITE(*,*)routineN, dbcsr_cuda_dev_mem_get_alloc(c_locks_dev)
       IF (nblks .GT. dbcsr_cuda_dev_mem_get_alloc(c_locks_dev)) THEN
          IF (dbg) WRITE(*,*)routineN//" reallocing locks"
          maxs = dbcsr_cuda_dev_mem_get_alloc(c_locks_dev)
          IF (dbg) WRITE(*,*)routineN//" lock size", maxs
          IF (detailed_timing) THEN
             t_dev_sync = t_dev_sync - m_walltime()
             CALL dbcsr_cuda_device_sync(error=error)
             t_dev_sync = t_dev_sync + m_walltime()
          ENDIF
          IF (dbg) WRITE(*,*)routineN//" reallocing locks"
          IF (verbose_acc .OR. info) WRITE(*,*)routineN//" reallocating locks", maxs
          IF (dbg) WRITE(*,*)routineN//" deallocating..."
          !--- !$OMP CRITICAL (crit_cuda)
          !$ IF (mm_workshare) CALL OMP_SET_LOCK (accel_lock)
          CALL dbcsr_cuda_dev_mem_dealloc(c_locks_dev, error=error)
          tmp_i = INT(REAL(nblks*4,kind=dp)*default_resize_factor)
          !WRITE(*,*)routineN//" reallocating locks", tmp_i, maxs
          IF (dbg) WRITE(*,*)routineN//" allocating..."
          CALL dbcsr_cuda_dev_mem_alloc(c_locks_dev, tmp_i, stat=istat)
          IF (istat /= 0) THEN
             IF (verbose_acc) WRITE(*,*)routineN//" Trying smaller allocation"
             tmp_i = nblks
             IF (dbg) WRITE(*,*)routineN//" failed, trying smaller..."
             CALL dbcsr_cuda_dev_mem_alloc(c_locks_dev, tmp_i, error=error)
          ENDIF
          IF (dbg) WRITE(*,*)routineN//" zeroing..."
          CALL dbcsr_cuda_dev_mem_zero(c_locks_dev,&
               first=1, last=dbcsr_cuda_dev_mem_get_alloc(c_locks_dev),&
               error=error)
          !$ IF (mm_workshare) CALL OMP_UNSET_LOCK (accel_lock)
          !--- !$OMP END CRITICAL (crit_cuda)
          IF (.FALSE.) THEN
             t_dev_sync = t_dev_sync - m_walltime()
             CALL dbcsr_cuda_device_sync(error=error)
             t_dev_sync = t_dev_sync + m_walltime()
          ENDIF
          !write(*,*)routineN//" done zeroing"
       ENDIF
    ENDIF cuda_process_target
    IF (dbg) WRITE(*,*)routineN//" done"
    zero_first = zero_last + 1
  END SUBROUTINE process_ps_target_low


! *****************************************************************************
!> \brief Unpacks stack variables to call lower-level stack processing
!>        drivers.
!*****************************************************************************
  SUBROUTINE process_ps_stack(param_stack, error)
    TYPE(dbcsr_ps_type), INTENT(INOUT)       :: param_stack
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'process_ps_stack', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: careful = careful_mod

    INTEGER                                  :: error_handle, stack_size

    IF (careful_mod) &
         CALL dbcsr_error_set(routineN, error_handle, error)
    IF (careful_mod) THEN
       IF (param_stack%driver .EQ. mm_driver_cuda) THEN
          CALL dbcsr_assert (param_stack%has_cuda_ab_data,&
               dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
               "A or B matrix data not present in stack.",&
               __LINE__, error=error)
       ELSE
          CALL dbcsr_assert (param_stack%has_ab_data,&
               dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
               "A or B matrix data not present in stack.",&
               __LINE__, error=error)
       ENDIF
    ENDIF
    !if (dbg) &
    !     write(*,*)routineN, OMP_GET_THREAD_NUM(), param_stack%stack_p,&
    !     param_stack%id, param_stack%sid, param_stack%t%t%owner
    !$ IF (.FALSE.) &
    !$      WRITE(*,*)"process_ps_stack      stack size is", param_stack%stack_p,&
    !$      omp_get_thread_num(), param_stack%t%t%owner, param_stack%driver,&
    !$      param_stack%sid, param_stack%id
    ! This has to be reset before the stack may be marked as empty.
    stack_size = param_stack%stack_p
    param_stack%stack_p = 0
    CALL process_mm_stack(param_stack%parameters,&
         param_stack%driver,&
         stack_size,&
         param_stack%left_data_area, param_stack%right_data_area,&
         param_stack%t%t%product_data_area,&
         param_stack%t%t%product_data_cuda, param_stack%t%t%has_cuda_c_data,&
         param_stack%left_data_cuda, param_stack%right_data_cuda,&
         param_stack%state, param_stack%t%t%stack_state_dev,&
         param_stack%m, param_stack%n, param_stack%k,&
         param_stack%max_m, param_stack%max_n, param_stack%max_k,&
         param_stack%defined_mnk,&
         param_stack%t%t%c_locks_dev, param_stack%t%t%params_dev,&
         param_stack%t%t, &
         error=error)
    IF (careful_mod) &
         CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE process_ps_stack


! *****************************************************************************
!> \brief Calls the various drivers that process the stack.
!>
!> \param[in] params           Stack of GEMM parameters
!> \param[in] n                Number of parameters
!> \param[in] left_data_area   Left-matrix data
!> \param[in] right_data_area  Right-matrix data
!> \param[in,out] zero_first   Zero product data area starting from this
!>                             element
!> \param[in,out] zero_last    Zero product data area up to this element
!> \param[in] lastblk          Number of blocks in product
!> \param[in,out] product_data_area  Data for results
! *****************************************************************************
  SUBROUTINE process_mm_stack(params, driver,&
       stack_size, &
       left_data_area, right_data_area, product_data_area,&
       product_data_card, has_product_data_card,&
       a_dev, b_dev,&
       state, stack_state_dev,&
       m, n, k, max_m, max_n, max_k, defined_mnk,&
       c_locks_dev, params_dev, &
       target_desc, &
       error)
    INTEGER, INTENT(IN)                      :: stack_size, driver
    INTEGER, DIMENSION(1:dbcsr_ps_width, &
      stack_size), INTENT(IN)                :: params
    TYPE(dbcsr_data_obj), INTENT(IN)         :: left_data_area, &
                                                right_data_area
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: product_data_area
    TYPE(dbcsr_cuda_mem_type), POINTER       :: product_data_card
    LOGICAL, INTENT(IN)                      :: has_product_data_card
    TYPE(dbcsr_cuda_mem_type), POINTER       :: a_dev, b_dev
    INTEGER, POINTER                         :: state
    TYPE(dbcsr_cuda_mem_type), POINTER       :: stack_state_dev
    INTEGER, INTENT(IN)                      :: m, n, k, max_m, max_n, max_k
    LOGICAL, INTENT(IN)                      :: defined_mnk
    TYPE(dbcsr_cuda_mem_type), POINTER       :: c_locks_dev, params_dev
    TYPE(dbcsr_ps_target_type), POINTER      :: target_desc
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'process_mm_stack', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: careful = careful_mod, &
                                                dbg = .FALSE.
    REAL, PARAMETER                          :: resize_factor = 1.618034

    INTEGER                                  :: nregions, sp
    LOGICAL, DIMENSION(max_regions)          :: which_locks
    REAL(KIND=dp)                            :: rnd, t_tmp

!   ---------------------------------------------------------------------------

    IF (detailed_timing)  THEN
     t_tmp=m_walltime()
     !$OMP ATOMIC
     t_process_stack = t_process_stack - t_tmp
    ENDIF  

    state = dbcsr_ps_state_working
    IF (dbg) THEN
       CALL RANDOM_NUMBER (rnd)
       IF (rnd < 0.01_dp) THEN
          WRITE(*,*)routineN//" Stack size", stack_size, dbcsr_ps_width
          CALL print_gemm_parameters(params(:,1:stack_size))
       ENDIF
    ENDIF
    !
    ! Verify stack consistency.  Only the upper bound is verified.
    IF (careful) THEN
       DO sp = 1, stack_size
          CALL dbcsr_assert (params(p_a_first,sp)&
               + params(p_m,sp) * params(p_k,sp) - 1,&
               "LE", dbcsr_data_get_size (left_data_area),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "A data out of bounds.", __LINE__, error=error)
          CALL dbcsr_assert (params(p_b_first,sp)&
               + params(p_k,sp) * params(p_n,sp) - 1,&
               "LE", dbcsr_data_get_size (right_data_area),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "B data out of bounds.", __LINE__, error=error)
          CALL dbcsr_assert (params(p_c_first,sp)&
               + params(p_m,sp) * params(p_n,sp) - 1,&
               "LE", dbcsr_data_get_size (product_data_area),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "C data out of bounds.", __LINE__, error=error)
       ENDDO
    ENDIF
    !
    ! Locks the target (data) area if there could be an access
    ! conflict there.
    !
    !$ IF (mm_workshare .AND. .NOT. driver_is_async (driver)) THEN
    !$    CALL xtime_start (mult_timers(tmr_wait_for_other_target))
    !$    CALL dbcsr_ps_target_lock_regions (target_desc,&
    !$    stack_size, params,&
    !$    target_desc%bit_shift,&
    !$    which_locks, nregions,&
    !$    error=error)
    !$    CALL xtime_stop (mult_timers(tmr_wait_for_other_target))
    !$ ENDIF
    IF (mm_workshare) THEN
       !$OMP FLUSH
    ENDIF
    IF (.NOT. mm_async) &
       CALL xtime_start(mult_timers(tmr_kernel_time))
    SELECT CASE (driver)
    CASE (mm_driver_matmul)
       SELECT CASE (product_data_area%d%data_type)
       CASE (dbcsr_type_real_4)
          CALL internal_process_mm_stack_s (params, &
               stack_size, &
               left_data_area%d%r_sp, right_data_area%d%r_sp, product_data_area%d%r_sp,&
               error=error)
       CASE (dbcsr_type_real_8)
          CALL internal_process_mm_stack_d (params,&
               stack_size,&
               left_data_area%d%r_dp, right_data_area%d%r_dp, product_data_area%d%r_dp,&
               error=error)
       CASE (dbcsr_type_complex_4)
          CALL internal_process_mm_stack_c (params,&
               stack_size,&
               left_data_area%d%c_sp, right_data_area%d%c_sp, product_data_area%d%c_sp,&
               error=error)
       CASE (dbcsr_type_complex_8)
          CALL internal_process_mm_stack_z (params,&
               stack_size,&
               left_data_area%d%c_dp, right_data_area%d%c_dp, product_data_area%d%c_dp,&
               error=error)
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type",__LINE__,error)
       END SELECT
    CASE (mm_driver_smm)
       SELECT CASE (product_data_area%d%data_type)
       CASE (dbcsr_type_real_4)
          CALL smm_process_mm_stack_s (params, &
               stack_size, &
               left_data_area%d%r_sp, right_data_area%d%r_sp, product_data_area%d%r_sp,&
               error=error)
       CASE (dbcsr_type_real_8)
          CALL smm_process_mm_stack_d (params,&
               stack_size,&
               left_data_area%d%r_dp, right_data_area%d%r_dp, product_data_area%d%r_dp,&
               error=error)
       CASE (dbcsr_type_complex_4)
          CALL smm_process_mm_stack_c (params,&
               stack_size,&
               left_data_area%d%c_sp, right_data_area%d%c_sp, product_data_area%d%c_sp,&
               error=error)
       CASE (dbcsr_type_complex_8)
          CALL smm_process_mm_stack_z (params,&
               stack_size,&
               left_data_area%d%c_dp, right_data_area%d%c_dp, product_data_area%d%c_dp,&
               error=error)
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type",__LINE__,error)
       END SELECT
    CASE (mm_driver_plasma)
       SELECT CASE (product_data_area%d%data_type)
       CASE (dbcsr_type_real_4)
          CALL plasma_process_mm_stack_s (params,&
               stack_size,&
               left_data_area%d%r_sp, right_data_area%d%r_sp, product_data_area%d%r_sp,&
               error=error)
       CASE (dbcsr_type_real_8)
          CALL plasma_process_mm_stack_d (params,&
               stack_size,&
               left_data_area%d%r_dp, right_data_area%d%r_dp, product_data_area%d%r_dp,&
               error=error)
       CASE (dbcsr_type_complex_4)
          CALL plasma_process_mm_stack_c (params,&
               stack_size,&
               left_data_area%d%c_sp, right_data_area%d%c_sp, product_data_area%d%c_sp,&
               error=error)
       CASE (dbcsr_type_complex_8)
          CALL plasma_process_mm_stack_z (params,&
               stack_size,&
               left_data_area%d%c_dp, right_data_area%d%c_dp, product_data_area%d%c_dp,&
               error=error)
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type",__LINE__,error)
       END SELECT
    CASE (mm_driver_blas)
       SELECT CASE (product_data_area%d%data_type)
       CASE (dbcsr_type_real_4)
          CALL blas_process_mm_stack_s (params,&
               stack_size,&
               left_data_area%d%r_sp, right_data_area%d%r_sp, product_data_area%d%r_sp,&
               error=error)
       CASE (dbcsr_type_real_8)
          CALL blas_process_mm_stack_d (params,&
               stack_size,&
               left_data_area%d%r_dp, right_data_area%d%r_dp, product_data_area%d%r_dp,&
               error=error)
       CASE (dbcsr_type_complex_4)
          CALL blas_process_mm_stack_c (params,&
               stack_size,&
               left_data_area%d%c_sp, right_data_area%d%c_sp, product_data_area%d%c_sp,&
               error=error)
       CASE (dbcsr_type_complex_8)
          CALL blas_process_mm_stack_z (params,&
               stack_size,&
               left_data_area%d%c_dp, right_data_area%d%c_dp, product_data_area%d%c_dp,&
               error=error)
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type",__LINE__,error)
       END SELECT
    CASE (mm_driver_cuda)
       IF (.NOT. has_product_data_card) &
            CALL dbcsr_assert (.FALSE.,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "No C data on card is specified.",&
            __LINE__, error=error)
       CALL cuda_process_mm_stack (params,&
            stack_size,&
            a_dev, b_dev, product_data_card,&
            c_locks_dev,&
            params_dev,&
            m, n, k, max_m, max_n, max_k, defined_mnk,&
            state, stack_state_dev,&
            target_desc%owner,&
            error=error)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Invalid multiplication driver",__LINE__,error)
    END SELECT
    IF (.NOT. mm_async) &
       CALL xtime_stop(mult_timers(tmr_kernel_time))
    !
    ! These writes have to be ordered.  Specifically the state should not
    ! be set before the stack has been used and the results written.
    IF (mm_workshare) THEN
       !$OMP FLUSH
    ENDIF
    IF (.NOT. driver_is_async (driver)) THEN
       state = dbcsr_ps_state_empty
    ENDIF
    IF (mm_workshare) THEN
       !$OMP FLUSH
    ENDIF
    !$ IF (mm_workshare .AND. .NOT. driver_is_async (driver)) THEN
    !$      CALL dbcsr_ps_target_unlock_regions (target_desc,&
    !$      which_locks, nregions,&
    !$      error=error)
    !$ ENDIF
    IF (detailed_timing) THEN
     t_tmp=m_walltime()
    !$OMP ATOMIC
     t_process_stack = t_process_stack + t_tmp
    ENDIF 
  END SUBROUTINE process_mm_stack


! *****************************************************************************
!> \brief Processes MM stack using CUDA.
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
!> \param[in,out] error        error
! *****************************************************************************
  SUBROUTINE cuda_process_mm_stack(params,&
       stack_size,&
       data_a_dev, data_b_dev, data_c_dev,&
       c_locks,&
       params_dev,&
       m, n, k, max_m, max_n, max_k, defined_mnk,&
       state, stack_state_dev, owner, &
       error)
    INTEGER, INTENT(IN)                      :: stack_size
    INTEGER, &
      DIMENSION(dbcsr_ps_width*stack_size), &
      INTENT(IN), TARGET                     :: params
    TYPE(dbcsr_cuda_mem_type), INTENT(IN)    :: data_a_dev, data_b_dev
    TYPE(dbcsr_cuda_mem_type), INTENT(INOUT) :: data_c_dev, c_locks, &
                                                params_dev
    INTEGER, INTENT(IN)                      :: m, n, k, max_m, max_n, max_k
    LOGICAL, INTENT(IN)                      :: defined_mnk
    INTEGER, POINTER                         :: state
    TYPE(dbcsr_cuda_mem_type), INTENT(IN)    :: stack_state_dev
    INTEGER, INTENT(IN)                      :: owner
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cuda_process_mm_stack', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, k_max, m_max, &
                                                n_max, stream
    INTEGER, DIMENSION(:), POINTER           :: params_p
    REAL(kind=dp)                            :: kt
    REAL(kind=dp), SAVE                      :: index_time = 0

!   ---------------------------------------------------------------------------

    IF (stack_size .EQ. 0) THEN
       !$OMP FLUSH
       state = dbcsr_ps_state_empty
       !$OMP FLUSH
       WRITE(*,*)routineN//" Stack size is 0"
       !$ write(*,*)"from thread", omp_get_thread_num()
       RETURN
    ENDIF
    stream = owner+1
    !$ IF (.FALSE.) &
    !$      WRITE(*,*)routineN//" Stack size is", stack_size, omp_get_thread_num()
    IF (careful_mod) &
         CALL dbcsr_error_set (routineN, error_handle, error)
    IF (.NOT. kernel_timing) t_calc_step = t_calc_step - m_walltime()
    IF (kernel_timing) THEN
       index_time = index_time + m_walltime()
       WRITE(*,'(1X,A20,1(1X,EN12.3,1X,A))')"Index time", index_time, "s"
    ENDIF
    params_p => params
    IF (kernel_timing) kt = -m_walltime()
    IF (m .GT. 0) THEN ; m_max = m ; ELSE ; m_max = -max_m ; ENDIF
    IF (n .GT. 0) THEN ; n_max = n ; ELSE ; n_max = -max_n ; ENDIF
    IF (k .GT. 0) THEN ; k_max = k ; ELSE ; k_max = -max_k ; ENDIF
    IF (kernel_timing) THEN
       kt = kt + m_walltime()
       WRITE(*,'(1X,A20,1(1X,EN12.3,1X,A))')"Max size time", kt, "s"
    ENDIF
    IF (kernel_timing) kt = -m_walltime()
    IF (kernel_timing) &
         CALL dbcsr_cuda_device_sync (error=error)
    IF (measure_idle) THEN
       t_dev_idle = t_dev_idle - m_walltime()
       CALL dbcsr_cuda_device_sync (error=error)
       t_dev_idle = t_dev_idle + m_walltime()
    ENDIF
    IF (kernel_timing) THEN
       kt = kt + m_walltime()
       WRITE(*,'(1X,A20,1(1X,EN12.3,1X,A))')"Thread sync time", kt, "s"
    ENDIF
    !
    IF (kernel_timing) kt = -m_walltime()
    CALL dbcsr_cuda_cp_host_to_dev (params_p,&
         params_dev%d_i, dbcsr_ps_width*stack_size,&
         async=.TRUE., stream=stream,&
         error=error)
    IF (kernel_timing) THEN
       kt = kt + m_walltime()
       WRITE(*,'(1X,A20,1(1X,EN12.3,1X,A))')"Parameter copy time", kt, "s"
    ENDIF
    !
    ! Schedule the resetting of the state of this stack as soon as the
    ! parameter upload finishes.  N.B. This must be the last change to
    ! the stack: its stack_p must be set to 0 and its data must be
    ! used <em>before</em> the state is set to dbcsr_ps_state_empty.
    IF (kernel_timing) kt = -m_walltime()
    CALL dbcsr_cuda_cp_dev_to_host (stack_state_dev%d_i, state,&
         async=.TRUE., stream=stream, error=error)
    IF (kernel_timing) THEN
       kt = kt + m_walltime()
       WRITE(*,'(1X,A20,1(1X,EN12.3,1X,A))')"State update time", kt, "s"
    ENDIF
    !
    IF (kernel_timing) CALL m_flush(6)
    IF (kernel_timing) kt = -m_walltime()
    !
    !--- !$OMP CRITICAL (crit_cuda)
    !$ IF (mm_workshare) CALL OMP_SET_LOCK (accel_lock)
    CALL dbcsr_cuda_do_mm_stack (params_dev%d_i, stack_size, dbcsr_ps_width,&
         data_a_dev, data_b_dev, data_c_dev,&
         c_locks%d_i, ABS(m_max), ABS(n_max), ABS(k_max), defined_mnk,&
         stream=stream, error=error)
    !$ IF (mm_workshare) CALL OMP_UNSET_LOCK (accel_lock)
    !--- !$OMP END CRITICAL (crit_cuda)
    IF (.FALSE.) CALL dbcsr_cuda_device_sync (error=error)
    IF (kernel_timing) THEN
       CALL dbcsr_cuda_device_sync (error=error)
       kt = kt + m_walltime()
       WRITE(*,'(1X,A20,1(1X,EN12.3,1X,A))')"Kernel time", kt, "s"
       t_calc_step = t_calc_step + kt
    ENDIF
    !
    IF (kernel_timing) index_time = -m_walltime()
    IF (.NOT. kernel_timing) t_calc_step = t_calc_step + m_walltime()

    IF (careful_mod) &
         CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE cuda_process_mm_stack



  SUBROUTINE print_gemm_parameters(params)
    INTEGER, DIMENSION(:, :), INTENT(in)     :: params

    INTEGER                                  :: sp

    DO sp = 1, SIZE(params,2)
       WRITE(*,'(1X,A,1X,I7,":",3(1X,I4,","),".",3(1X,I12,","))')&
            "GEMM PARAMETERS",&
               sp,&
               params(p_m,sp),&
               params(p_k,sp),&
               params(p_n,sp),&
               params(p_a_first,sp),&
               params(p_b_first,sp),&
               params(p_c_first,sp)
    ENDDO
  END SUBROUTINE print_gemm_parameters


#include "dbcsr_mm_stack_d.F"
#include "dbcsr_mm_stack_z.F"
#include "dbcsr_mm_stack_s.F"
#include "dbcsr_mm_stack_c.F"


END MODULE dbcsr_mm_stack
