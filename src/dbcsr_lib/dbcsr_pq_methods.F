!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  Urban Borstnik and the CP2K developers group   !
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
  USE array_types,                     ONLY: array_hold,&
                                             array_i1d_obj
  USE dbcsr_cuda_memory,               ONLY: dbcsr_cuda_dev_mem_hold,&
                                             dbcsr_cuda_dev_mem_release
  USE dbcsr_cuda_types,                ONLY: dbcsr_cuda_mem_type
  USE dbcsr_data_methods,              ONLY: dbcsr_data_hold,&
                                             dbcsr_data_init,&
                                             dbcsr_data_release
  USE dbcsr_error_handling
  USE dbcsr_pq_types,                  ONLY: &
       dbcsr_pq_type, dbcsr_ps_fin_q_type, dbcsr_ps_fin_type, &
       dbcsr_ps_group_type, dbcsr_ps_obj, dbcsr_ps_set_type, &
       dbcsr_ps_state_empty, dbcsr_ps_state_filled, dbcsr_ps_state_filling, &
       dbcsr_ps_state_working, dbcsr_ps_target_obj, dbcsr_ps_target_type, &
       dbcsr_ps_type, dbcsr_ps_width
  USE dbcsr_ptr_util,                  ONLY: memory_allocate,&
                                             memory_deallocate
  USE dbcsr_types,                     ONLY: dbcsr_data_obj

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_pq_methods'

  LOGICAL, PARAMETER :: careful_mod = .FALSE.

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
  PUBLIC :: dbcsr_psg_set_state

  PUBLIC :: dbcsr_ps_target_init, dbcsr_ps_target_new
  PUBLIC :: dbcsr_ps_target_create, dbcsr_ps_target_destroy
  PUBLIC :: dbcsr_ps_target_hold, dbcsr_ps_target_release
  PUBLIC :: dbcsr_ps_target_rm_data, dbcsr_ps_target_add_data
  PUBLIC :: dbcsr_ps_target_rm_data_cuda, dbcsr_ps_target_add_data_cuda

  PUBLIC :: dbcsr_ps_set_create, dbcsr_ps_set_destroy
  PUBLIC :: dbcsr_ps_set_set_b_master
  PUBLIC :: dbcsr_ps_set_advance,&
            dbcsr_ps_set_get_group_p

  PUBLIC :: dbcsr_pq_create, dbcsr_pq_destroy
  PUBLIC :: dbcsr_pq_add_stack
  PUBLIC :: dbcsr_pq_get_a_stack, dbcsr_pq_get_defined_stack

  PUBLIC :: dbcsr_ps_fin_q_add, dbcsr_ps_fin_q_pop

CONTAINS


! *****************************************************************************
!> \brief Creates a queue of parameter stacks
!> \param[out] queue              Queue of parameter stacks to create.
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_pq_create (queue, error)
    TYPE(dbcsr_pq_type), INTENT(OUT)         :: queue
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    ALLOCATE (queue%stacks(1))
    queue%nstacks = 0
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
    queue%nstacks = 0
  END SUBROUTINE dbcsr_pq_destroy


! *****************************************************************************
!> \brief Adds a stack to a queue
!> \param[in,out] queue              Queue of parameter stacks to create.
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_pq_add_stack (queue, stack, error)
    TYPE(dbcsr_pq_type), INTENT(INOUT)       :: queue
    TYPE(dbcsr_ps_obj), INTENT(IN)           :: stack
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_pq_add_stack', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, new_n, old_s
    TYPE(dbcsr_ps_obj), DIMENSION(:), &
      POINTER                                :: new_queue

    IF (careful_mod) &
         CALL dbcsr_error_set(routineN, error_handle, error)
    new_n = queue%nstacks + 1
    old_s = SIZE (queue%stacks)
    IF (new_n .GT. old_s) THEN
       ALLOCATE (new_queue (MAX(new_n, old_s*2)))
       new_queue(1:queue%nstacks) = queue%stacks(1:queue%nstacks)
       DEALLOCATE (queue%stacks)
       queue%stacks => new_queue
    ENDIF
    queue%stacks(new_n) = stack
    queue%nstacks = new_n
    IF (careful_mod) &
         CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_pq_add_stack


! *****************************************************************************
!> \brief Gets a stack from the queue
!> \param[in,out] queue              Queue of parameter stacks to create.
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_pq_get_a_stack (queue, stack, found, error)
    TYPE(dbcsr_pq_type), INTENT(INOUT)       :: queue
    TYPE(dbcsr_ps_obj), INTENT(OUT)          :: stack
    LOGICAL, INTENT(OUT)                     :: found
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    IF (queue%nstacks .GE. 1) THEN
       stack = queue%stacks(queue%nstacks)
       queue%nstacks = queue%nstacks - 1
       found = .TRUE.
    ELSE
       found = .FALSE.
    ENDIF
  END SUBROUTINE dbcsr_pq_get_a_stack


! *****************************************************************************
!> \brief Gets a stack from the queue
!> \param[in,out] queue              Queue of parameter stacks to create.
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_pq_get_defined_stack (queue, stack, found, error)
    TYPE(dbcsr_pq_type), INTENT(INOUT)       :: queue
    TYPE(dbcsr_ps_obj), INTENT(OUT)          :: stack
    LOGICAL, INTENT(OUT)                     :: found
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: i
    TYPE(dbcsr_ps_type)                      :: s

    IF (queue%nstacks .GE. 1) THEN
       found = .FALSE.
       DO i = queue%nstacks, 1, -1
          s = queue%stacks(i)%s
          !if (s%m .gt. 0 .and. &
          !     s%n .gt. 0 .and. &
          !     s%k .gt. 0) then
          !if (s%defined_mnk) then
          IF (.TRUE.) THEN
             found = .TRUE.
             stack = queue%stacks(i)
             queue%stacks(i:queue%nstacks-1) = &
                  queue%stacks(i+1:queue%nstacks)
             queue%nstacks = queue%nstacks - 1
             RETURN
          ENDIF
       ENDDO
    ELSE
       found = .FALSE.
    ENDIF
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

    t%zero_first = 1
    t%zero_last = 0
    t%has_c_data = .FALSE.
    t%has_cuda_c_data = .FALSE.
    t%last_c_blk = 0
    NULLIFY (t%c_locks_dev)
    NULLIFY (t%stack_state_dev)
    NULLIFY (t%params_dev)
    CALL dbcsr_ps_fin_q_create (t%fin_queue, error)
  END SUBROUTINE dbcsr_ps_target_create

  SUBROUTINE dbcsr_ps_target_destroy (t, error)
    TYPE(dbcsr_ps_target_type), &
      INTENT(INOUT)                          :: t
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    IF (t%has_c_data) THEN
       CALL dbcsr_data_release (t%product_data_area)
    ENDIF
    IF (t%has_cuda_c_data) THEN
       CALL dbcsr_cuda_dev_mem_release (t%product_data_cuda, error=error)
    ENDIF
    IF (ASSOCIATED (t%c_locks_dev)) THEN
       CALL dbcsr_cuda_dev_mem_release (t%c_locks_dev, error=error)
    ENDIF
    IF (ASSOCIATED (t%stack_state_dev)) THEN
       CALL dbcsr_cuda_dev_mem_release (t%stack_state_dev, error=error)
    ENDIF
    IF (ASSOCIATED (t%params_dev)) THEN
       CALL dbcsr_cuda_dev_mem_release (t%params_dev, error=error)
    ENDIF
    CALL dbcsr_ps_fin_q_destroy (t%fin_queue, error)
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
    stack%start_p = 0
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
       CALL memory_deallocate (stack%parameters, stack%mem_type, error=error)
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

    INTEGER                                  :: error_handle, i, j, k
    TYPE(dbcsr_ps_obj)                       :: memreg_stack

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set(routineN, error_handle, error)
    ALLOCATE (stack_set%groups (nbuffers, nmemregs))
    stack_set%nbuffers = nbuffers
    stack_set%nmemregs = nmemregs
    stack_set%next_buffer = 1
    stack_set%next_memreg = 1
    ALLOCATE (stack_set%all_states (nmemregs))
    DO k = 1, stack_set%nmemregs
       ! Allocates the state array for all (non-master) stacks in this
       ! memory region.
       CALL memory_allocate (stack_set%all_states(k)%data,&
            group_size * nbuffers,&
            mem_type=mem_type, error=error)
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
       CALL dbcsr_ps_new (memreg_stack,&
            stack_size*group_size*nbuffers,&
            mem_type=mem_type, error=error)
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
       CALL memory_deallocate (stack_set%all_states(k)%data,&
                               stack_set%all_states(k)%mem_type, error=error)
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
!> \param[in,out] error    error
!> \par master_array
!>      The master array is the pointer to the master array associated
!>      with the stack group.  The starting elemnent of this group or
!>      stack into this master array is obtained by opening a
!>      view. \see dbcsr_psg_view_open
! *****************************************************************************
  FUNCTION dbcsr_ps_set_get_group_p (stack_set, buffer, memreg,&
       master_array, error) RESULT (ps_group)
    TYPE(dbcsr_ps_set_type), INTENT(INOUT)   :: stack_set
    INTEGER, INTENT(IN), OPTIONAL            :: buffer, memreg
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: master_array
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    TYPE(dbcsr_ps_group_type), POINTER       :: ps_group

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_ps_set_get_group_p', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a, b, error_handle

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
    ! Busy loop until processing of the requested group has finished.
    DO WHILE (dbcsr_psg_get_state (stack_set%groups(b, a), error) &
         .EQ. dbcsr_ps_state_working)
    ENDDO
    ps_group => stack_set%groups(b, a)
    IF (PRESENT (master_array)) THEN
       master_array => ps_group%master%s%parameters
    ENDIF
    IF (careful_mod) CALL dbcsr_error_stop(error_handle, error)
  END FUNCTION dbcsr_ps_set_get_group_p

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

    INTEGER                                  :: last_a

! First advance group; if it wraps then advance the buffer

    last_a = stack_set%next_memreg
    IF (advance_memreg) &
         stack_set%next_memreg = MOD(stack_set%next_memreg, stack_set%nmemregs)+1
    IF (stack_set%next_memreg .LE. last_a .OR. .NOT. advance_memreg) THEN
       stack_set%next_buffer = MOD(stack_set%next_buffer, stack_set%nbuffers)+1
    ENDIF
  END SUBROUTINE dbcsr_ps_set_advance


  SUBROUTINE dbcsr_ps_fin_q_create (fin_queue, error)
    TYPE(dbcsr_ps_fin_q_type), INTENT(OUT)   :: fin_queue
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    ALLOCATE (fin_queue%f(1))
    fin_queue%n = 0
  END SUBROUTINE dbcsr_ps_fin_q_create

  SUBROUTINE dbcsr_ps_fin_q_destroy (fin_queue, error)
    TYPE(dbcsr_ps_fin_q_type), INTENT(INOUT) :: fin_queue
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    fin_queue%n = 0
    DEALLOCATE (fin_queue%f)
  END SUBROUTINE dbcsr_ps_fin_q_destroy

  SUBROUTINE dbcsr_ps_fin_q_add (fin_queue, DATA, xlate, nxlate, error)
    TYPE(dbcsr_ps_fin_q_type), INTENT(INOUT) :: fin_queue
    TYPE(dbcsr_data_obj), INTENT(IN)         :: DATA
    TYPE(array_i1d_obj), INTENT(IN)          :: xlate
    INTEGER, INTENT(IN)                      :: nxlate
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_ps_fin_q_add', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: new_n, new_s, old_n, old_s
    TYPE(dbcsr_ps_fin_type), DIMENSION(:), &
      POINTER                                :: tmp_f

    old_n = fin_queue%n
    new_n = old_n + 1
    old_s = SIZE (fin_queue%f)
    IF (new_n .GT. old_s) THEN
       new_s = MAX(new_n, old_s*2)
       ALLOCATE (tmp_f (new_s))
       tmp_f(1:old_n) = fin_queue%f(1:old_n)
       DEALLOCATE (fin_queue%f)
       fin_queue%f => tmp_f
    ENDIF
    fin_queue%n = new_n
    fin_queue%f(new_n)%src_data = DATA
    fin_queue%f(new_n)%xlate = xlate
    fin_queue%f(new_n)%n = nxlate
    CALL dbcsr_data_hold (fin_queue%f(new_n)%src_data)
    CALL array_hold (fin_queue%f(new_n)%xlate)
  END SUBROUTINE dbcsr_ps_fin_q_add

  SUBROUTINE dbcsr_ps_fin_q_pop (fin_queue, DATA, xlate, nxlate, found, error)
    TYPE(dbcsr_ps_fin_q_type), INTENT(INOUT) :: fin_queue
    TYPE(dbcsr_data_obj), INTENT(OUT)        :: DATA
    TYPE(array_i1d_obj), INTENT(OUT)         :: xlate
    INTEGER, INTENT(OUT)                     :: nxlate
    LOGICAL, INTENT(OUT)                     :: found
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_ps_fin_q_pop', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: old_n

    old_n = fin_queue%n
    IF (old_n .GE. 1) THEN
       DATA  = fin_queue%f(old_n)%src_data
       xlate = fin_queue%f(old_n)%xlate
       nxlate = fin_queue%f(old_n)%n
       fin_queue%n = old_n - 1
       found = .TRUE.
    ELSE
       found = .FALSE.
    ENDIF
  END SUBROUTINE dbcsr_ps_fin_q_pop

END MODULE dbcsr_pq_methods
