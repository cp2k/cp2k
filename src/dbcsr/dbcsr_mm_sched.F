!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Fourth layer of the dbcsr matrix-matrix multiplication.
!>         It hides the differences between performing calculations on the
!>         accelerator device or on the CPU.
!>
!> \author  Urban Borstnik
!>
!> <b>Modification history:</b>
!>  - 2010-02-23 Moved from dbcsr_operations
!>  - 2011-11    Moved parameter-stack processing routines to
!>               dbcsr_mm_methods.
!>  - 2013-01    extensive refactoring (Ole Schuett)
! *****************************************************************************


MODULE dbcsr_mm_sched
  USE dbcsr_block_operations,          ONLY: dbcsr_data_clear
  USE dbcsr_config,                    ONLY: accdrv_avoid_after_busy,&
                                             accdrv_do_inhomogenous,&
                                             accdrv_min_flop_process,&
                                             default_resize_factor,&
                                             mm_driver,&
                                             mm_driver_acc
  USE dbcsr_data_methods,              ONLY: dbcsr_data_ensure_size,&
                                             dbcsr_data_get_size
  USE dbcsr_error_handling
  USE dbcsr_mm_accdrv,                 ONLY: dbcsr_mm_accdrv_barrier,&
                                             dbcsr_mm_accdrv_finalize,&
                                             dbcsr_mm_accdrv_init,&
                                             dbcsr_mm_accdrv_lib_finalize,&
                                             dbcsr_mm_accdrv_lib_init,&
                                             dbcsr_mm_accdrv_phaseout,&
                                             dbcsr_mm_accdrv_process,&
                                             dbcsr_mm_accdrv_type
  USE dbcsr_mm_hostdrv,                ONLY: dbcsr_mm_hostdrv_barrier,&
                                             dbcsr_mm_hostdrv_finalize,&
                                             dbcsr_mm_hostdrv_init,&
                                             dbcsr_mm_hostdrv_phaseout,&
                                             dbcsr_mm_hostdrv_process,&
                                             dbcsr_mm_hostdrv_type
  USE dbcsr_mm_types,                  ONLY: p_a_first,&
                                             p_b_first,&
                                             p_c_first,&
                                             p_k,&
                                             p_m,&
                                             p_n,&
                                             stack_descriptor_type
  USE dbcsr_toollib,                   ONLY: sort
  USE dbcsr_types,                     ONLY: dbcsr_type,&
                                             dbcsr_work_type
  USE kinds,                           ONLY: int_4,&
                                             int_8
  USE message_passing,                 ONLY: mp_abort,&
                                             mp_bcast,&
                                             mp_environ,&
                                             mp_max,&
                                             mp_sum

 !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_sched'

  PUBLIC :: dbcsr_mm_sched_type
  PUBLIC :: dbcsr_mm_sched_lib_init, dbcsr_mm_sched_lib_finalize
  PUBLIC :: dbcsr_mm_sched_phaseout
  PUBLIC :: dbcsr_mm_sched_init, dbcsr_mm_sched_finalize
  PUBLIC :: dbcsr_mm_sched_process
  PUBLIC :: dbcsr_mm_sched_begin_burst, dbcsr_mm_sched_end_burst
  PUBLIC :: dbcsr_mm_sched_barrier

  ! *****************************************************************************
  TYPE dbcsr_mm_sched_type
    PRIVATE
    TYPE(dbcsr_work_type), POINTER  :: product_wm
    TYPE(dbcsr_mm_accdrv_type)      :: accdrv
    TYPE(dbcsr_mm_hostdrv_type)     :: hostdrv
    LOGICAL                         :: avoid_accdrv = .FALSE.
    LOGICAL                         :: product_wm_cleared = .FALSE.
    INTEGER                         :: product_wm_orig_datasize = -1
  END TYPE dbcsr_mm_sched_type

  ! *****************************************************************************
  TYPE stats_type
    INTEGER(kind=int_8)   :: cpu_num_stacks = 0
    INTEGER(kind=int_8)   :: acc_num_stacks = 0
    INTEGER(kind=int_8)   :: acc_flop = 0
    INTEGER(kind=int_8)   :: cpu_flop = 0
  END TYPE stats_type

! *****************************************************************************
!>  Statistics Counter
! *****************************************************************************
  ! Threadprivate counters to collect statistics
  TYPE(stats_type), SAVE :: stats
  !$OMP THREADPRIVATE(stats)

  INTEGER(kind=int_8), DIMENSION(:,:), POINTER, SAVE :: stats_num_mnk_stacks
  !$OMP THREADPRIVATE(stats_num_mnk_stacks)

  CONTAINS


! *****************************************************************************
! \brief Initialize the library
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_sched_lib_init(error)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: ithread, nthreads

    nthreads = 1; ithread = 0
    !$ nthreads = OMP_GET_NUM_THREADS () ; ithread = OMP_GET_THREAD_NUM ()

    ! reset statistics counter
    stats = stats_type()
    ALLOCATE(stats_num_mnk_stacks(1,5))
    stats_num_mnk_stacks(1,:) = 0 ! entry for the default stack

    CALL dbcsr_mm_accdrv_lib_init(error)


  END SUBROUTINE dbcsr_mm_sched_lib_init


! *****************************************************************************
! \brief Finalize the library and prints DBCSR statistics
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_sched_lib_finalize(group, output_unit, error)
    INTEGER, INTENT(IN)                      :: group, output_unit
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    TYPE(stats_type), SAVE                   :: report

   CALL dbcsr_mm_accdrv_lib_finalize(output_unit, error)

   ! Collect and output statistics ---------------------------------------------

   CALL stats_collect_from_threads(report, error)
   !$OMP BARRIER

   !$OMP MASTER
      CALL stats_collect_from_ranks(report, group, error)
      CALL stats_print_report(report, output_unit, error)
   !$OMP END MASTER
   DEALLOCATE(stats_num_mnk_stacks)
  END SUBROUTINE dbcsr_mm_sched_lib_finalize



! *****************************************************************************
!> \brief Makes sure that the product_wm is cleared.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE ensure_product_wm_cleared(this)
    TYPE(dbcsr_mm_sched_type), INTENT(INOUT) :: this

    INTEGER                                  :: allocated_datasize, &
                                                used_datasize

    IF(this%product_wm_cleared) RETURN

    ! The product's data_area could already contain some data.
    ! ( see: keep_product_data in dbcsr_operations.F )
    ! But this data might not occupy all the allocated memory in the data_area.
    ! Since, we don't want to keep track of unitialized memory we just zero it now.

    used_datasize = this%product_wm_orig_datasize
    allocated_datasize = dbcsr_data_get_size(this%product_wm%data_area)
    CALL dbcsr_data_clear(this%product_wm%data_area, lb=used_datasize+1, ub=allocated_datasize)
    this%product_wm_cleared = .TRUE.
  END SUBROUTINE ensure_product_wm_cleared


! *****************************************************************************
!> \brief Initializes a multiplication cycle for new set of C-blocks.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_sched_init(this, left, right, product_wm, error)
    TYPE(dbcsr_mm_sched_type), INTENT(INOUT) :: this
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_work_type), POINTER           :: product_wm
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_sched_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)

    this%product_wm => product_wm

    ! Clearing the product_wm takes too long, we gonna do it later and
    ! return now to allow for MPI to progress.
    ! We just have to remeber its datasize, in case it already contains data.
    this%product_wm_orig_datasize = this%product_wm%datasize

    CALL dbcsr_mm_hostdrv_init(this%hostdrv, left, right, product_wm, error)

    IF (mm_driver == mm_driver_acc) &
       CALL dbcsr_mm_accdrv_init(this%accdrv, left, right, product_wm, error)

    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_mm_sched_init

! *****************************************************************************
! \brief Signal approaching end of multiplication
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_sched_phaseout(this, error)
    TYPE(dbcsr_mm_sched_type), INTENT(INOUT) :: this
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CALL dbcsr_mm_hostdrv_phaseout(this%hostdrv, error)
    IF (mm_driver == mm_driver_acc) &
       CALL dbcsr_mm_accdrv_phaseout(this%accdrv, error)
  END SUBROUTINE dbcsr_mm_sched_phaseout


! *****************************************************************************
!> \brief Finalizes a multiplication cycle for a set of C-blocks.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_sched_finalize(this, error)
    TYPE(dbcsr_mm_sched_type), INTENT(INOUT) :: this
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_sched_finalize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)

    ! Just in case dbcsr_mm_sched_process was never called (really needed?)
    CALL ensure_product_wm_cleared(this)

    CALL dbcsr_mm_hostdrv_finalize(this%hostdrv, error)
    IF (mm_driver == mm_driver_acc) &
       CALL dbcsr_mm_accdrv_finalize(this%accdrv, error)

    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_mm_sched_finalize

! *****************************************************************************
! \brief Signal begin of a burst of calls to dbcsr_mm_sched_process.
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_sched_begin_burst(this, error)
    TYPE(dbcsr_mm_sched_type), INTENT(INOUT) :: this
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    this%avoid_accdrv = .FALSE.
  END SUBROUTINE dbcsr_mm_sched_begin_burst

! *****************************************************************************
! \brief Signal end of a burst of calls to dbcsr_mm_sched_process.
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_sched_end_burst(this, error)
    TYPE(dbcsr_mm_sched_type), INTENT(INOUT) :: this
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

!nothing to do here

  END SUBROUTINE dbcsr_mm_sched_end_burst

! *****************************************************************************
! \brief Signal that previous stacks should be processed first
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_sched_barrier(this, error)
    TYPE(dbcsr_mm_sched_type), INTENT(INOUT) :: this
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CALL dbcsr_mm_hostdrv_barrier(this%hostdrv, error)
    IF (mm_driver==mm_driver_acc) &
       CALL dbcsr_mm_accdrv_barrier(this%accdrv, error)

  END SUBROUTINE dbcsr_mm_sched_barrier

! *****************************************************************************
! \brief Processes a given stack.
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_sched_process(this, left, right, stack_data,&
                    stack_fillcount, stack_descr, error)
    TYPE(dbcsr_mm_sched_type), INTENT(INOUT) :: this
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    INTEGER, DIMENSION(:, :), POINTER        :: stack_data
    INTEGER, POINTER                         :: stack_fillcount
    TYPE(stack_descriptor_type), INTENT(IN)  :: stack_descr
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: sp, stacked_datasize
    INTEGER(kind=int_8)                      :: flop_per_entry, total_flop
    LOGICAL                                  :: success

    IF(stack_fillcount <= 0) CALL mp_abort("dbcsr_mm_sched_process: got empty stack")

    CALL ensure_product_wm_cleared(this)

    stacked_datasize   = this%product_wm%datasize
    CALL dbcsr_data_ensure_size(this%product_wm%data_area, stacked_datasize, &
                    factor=default_resize_factor, zero_pad=.TRUE., error=error)

     !!From here on there is no boundary checking due to assumed-SIZE-arguments.
     !!This is usefull to check stack parameters, BUT it works only for kind=dp
     IF(.FALSE.) THEN
        DO sp = 1, stack_fillcount
           IF(stack_data(p_a_first,sp) > SIZE(left%data_area%d%r_dp)) &
              STOP "left data out of range"
           IF(stack_data(p_b_first,sp) > SIZE(right%data_area%d%r_dp)) &
              STOP "right data out of range"
           IF(stack_data(p_c_first,sp) > SIZE(this%product_wm%data_area%d%r_dp)) THEN
              WRITE (*,*) "blub: ",stack_data(p_c_first,sp) , SIZE(this%product_wm%data_area%d%r_dp), &
                dbcsr_data_get_size(this%product_wm%data_area), stacked_datasize
              STOP "product data out of range"
           END IF
        END DO
     ENDIF


     IF(.FALSE.) THEN
       ! Check if homogenous stacks are indeed homogenous
       IF(stack_descr%defined_mnk) THEN
         DO sp = 1, stack_fillcount
           IF(stack_data(p_m,sp) /= stack_descr%m) CALL mp_abort("homogenous stacks check failed")
           IF(stack_data(p_n,sp) /= stack_descr%n) CALL mp_abort("homogenous stacks check failed")
           IF(stack_data(p_k,sp) /= stack_descr%k) CALL mp_abort("homogenous stacks check failed")
         END DO
       ENDIF
     ENDIF

     ! Submitting the stack for processing -------------------------------------
     flop_per_entry = 2*stack_descr%max_m*stack_descr%max_n*stack_descr%max_k
     total_flop = stack_fillcount * flop_per_entry

     IF (mm_driver==mm_driver_acc .AND. &
         flop_per_entry > accdrv_min_flop_process .AND. &
         (.NOT. this%avoid_accdrv) .AND. &
         (stack_descr%defined_mnk .OR. accdrv_do_inhomogenous) ) THEN
            CALL dbcsr_mm_accdrv_process(&
                    this%accdrv,&
                    left, right,&
                    params=stack_data,&
                    stack_size=stack_fillcount, &
                    stack_descr=stack_descr,&
                    success=success,&
                    error=error)

         IF(success) THEN
            ! update statistics
            stats%acc_num_stacks = stats%acc_num_stacks + 1
            stats%acc_flop = stats%acc_flop + total_flop
            CALL stats_add(m=stack_descr%m, n=stack_descr%n, k=stack_descr%k, &
                     stacksize_acc=INT(stack_fillcount,kind=int_8))
            RETURN
         ELSE
            this%avoid_accdrv = accdrv_avoid_after_busy
         ENDIF
     ENDIF


     !WRITE (*,*) "dbcsr_mm_sched_process: running hostdrv_process, stack_fillcount:", stack_fillcount

     CALL dbcsr_mm_hostdrv_process(&
                this%hostdrv,&
                left, right,&
                params=stack_data,&
                stack_size=stack_fillcount, &
                stack_descr=stack_descr,&
                success=success,&
                error=error)

     IF(success) THEN
        ! update statistics
        stats%cpu_num_stacks = stats%cpu_num_stacks + 1
        stats%cpu_flop = stats%cpu_flop + total_flop
        CALL stats_add(m=stack_descr%m, n=stack_descr%n, k=stack_descr%k, &
                 stacksize_cpu=INT(stack_fillcount,kind=int_8))
        RETURN
     ENDIF

     STOP "dbcsr_mm_sched_process_stack failed"


  END SUBROUTINE dbcsr_mm_sched_process


! ******************************************************************************
! \brief Helper-routine used by dbcsr_mm_sched_process to supply statistics.
! \author Ole Schuett
! ******************************************************************************
 SUBROUTINE stats_add(m,n,k, stacksize_cpu, stacksize_acc)
    INTEGER, INTENT(IN)                      :: m, n, k
    INTEGER(kind=int_8), OPTIONAL            :: stacksize_cpu, stacksize_acc

    INTEGER                                  :: i, s
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tmp

   DO i=1, SIZE(stats_num_mnk_stacks, 1)
       IF(stats_num_mnk_stacks(i,1)==m .AND. &
          stats_num_mnk_stacks(i,2)==n .AND. &
          stats_num_mnk_stacks(i,3)==k) THEN
            IF(PRESENT(stacksize_cpu))&
              stats_num_mnk_stacks(i,4) = stats_num_mnk_stacks(i,4) + stacksize_cpu
            IF(PRESENT(stacksize_acc))&
              stats_num_mnk_stacks(i,5) = stats_num_mnk_stacks(i,5) + stacksize_acc
            RETURN
       ENDIF
   END DO

   !not found, ok lets grow the list
   s = SIZE(stats_num_mnk_stacks, 1)
   tmp => stats_num_mnk_stacks
   ALLOCATE(stats_num_mnk_stacks(s+1, 5))
   stats_num_mnk_stacks(1:s,:) = tmp(:,:)
   stats_num_mnk_stacks(s+1,1) = m
   stats_num_mnk_stacks(s+1,2) = n
   stats_num_mnk_stacks(s+1,3) = k
   stats_num_mnk_stacks(s+1,4:5) = 0
   IF(PRESENT(stacksize_cpu))&
     stats_num_mnk_stacks(s+1,4) = stacksize_cpu
   IF(PRESENT(stacksize_acc))&
     stats_num_mnk_stacks(s+1,5) = stacksize_acc
   DEALLOCATE(tmp)
 END SUBROUTINE stats_add


! *****************************************************************************
! \brief Collects statistics from all OpenMP-threads
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE stats_collect_from_threads(report, error)
    TYPE(stats_type), INTENT(INOUT)          :: report
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: i, ithread, j, nthreads
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER, SAVE                          :: handover

   nthreads = 1; ithread = 0
   !$ nthreads = OMP_GET_NUM_THREADS () ; ithread = OMP_GET_THREAD_NUM ()

   !$OMP CRITICAL
   report%cpu_num_stacks = report%cpu_num_stacks + stats%cpu_num_stacks
   report%acc_num_stacks = report%acc_num_stacks + stats%acc_num_stacks
   report%acc_flop = report%acc_flop + stats%acc_flop
   report%cpu_flop = report%cpu_flop + stats%cpu_flop
   !$OMP END CRITICAL

   DO i=1, nthreads-1
     !$OMP BARRIER
     IF(ithread==i) handover => stats_num_mnk_stacks
     !$OMP BARRIER
     IF(ithread==0) THEN
       DO j=1, SIZE(handover, 1)
         CALL stats_add(m=INT(handover(j,1), kind=int_4),&
                        n=INT(handover(j,2), kind=int_4),&
                        k=INT(handover(j,3), kind=int_4),&
                        stacksize_cpu=handover(j,4),&
                        stacksize_acc=handover(j,5))
       END DO
     END IF
   END DO
  END SUBROUTINE stats_collect_from_threads


! *****************************************************************************
! \brief Collects statistics from all MPI-ranks
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE stats_collect_from_ranks(report, group, error)
    TYPE(stats_type), INTENT(INOUT)          :: report
    INTEGER, INTENT(IN)                      :: group
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: i, myrank, nranks, &
                                                sending_rank
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: mnk_collected
    INTEGER, DIMENSION(3)                    :: mnk

    CALL mp_environ(nranks, myrank, group)

    CALL mp_sum(report%acc_flop, group)
    CALL mp_sum(report%cpu_flop, group)
    CALL mp_sum(report%cpu_num_stacks, group)
    CALL mp_sum(report%acc_num_stacks, group)

    ! array mnk_collected is used as a logical-array, allows to use minloc
    ALLOCATE(mnk_collected(SIZE(stats_num_mnk_stacks,1)))
    mnk_collected = 0 ! init all to false

    ! broadcast stats of all mnk-combinations, which occured on any mpi rank
    DO
      ! each rank with uncollected stats tries to become the sending_rank
      sending_rank = -1
      IF(.NOT. ALL(mnk_collected==1)) sending_rank = myrank
      CALL mp_max(sending_rank, group)
      IF(sending_rank < 0) EXIT ! every rank got all mnk collected

      IF(sending_rank == myrank) THEN
        i = MINLOC(mnk_collected, dim=1)
        mnk = INT(stats_num_mnk_stacks(i,1:3), kind=int_4)
      ENDIF
      CALL mp_bcast(msg=mnk, source=sending_rank, gid=group)

      CALL stats_add(m=mnk(1), n=mnk(2), k=mnk(3), stacksize_cpu=0_int_8, stacksize_acc=0_int_8)
      DO i=1, SIZE(stats_num_mnk_stacks, 1)
         IF(ALL(stats_num_mnk_stacks(i,1:3)==mnk)) THEN
            IF(i<=SIZE(mnk_collected)) mnk_collected(i) = 1
            CALL mp_sum(stats_num_mnk_stacks(i,4:5), group)
         END IF
      END DO
    END DO

  END SUBROUTINE stats_collect_from_ranks


! *****************************************************************************
! \brief Prints collected statistics
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE stats_print_report(report, output_unit, error)
    TYPE(stats_type), INTENT(INOUT)          :: report
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: i, j
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: sort_idx, sort_key
    REAL                                     :: percent

    IF(output_unit<=0) RETURN

    WRITE (output_unit,"(1X,A,T47,A,T68,A,T77,A)") "COUNTER", "CPU", "ACC", "ACC%"

    percent = 100*REAL(report%acc_num_stacks)/REAL(MAX(INT(1,KIND=int_8),report%cpu_num_stacks+report%acc_num_stacks))
    WRITE (output_unit,'(A,T30,I20,1X,I20,T76,F5.1)') " number of processed stacks",&
       report%cpu_num_stacks, report%acc_num_stacks, percent

    percent = 100*REAL(stats_num_mnk_stacks(1,5))/REAL(MAX(INT(1,KIND=int_8),SUM(stats_num_mnk_stacks(1,4:5))))
    WRITE (output_unit,'(A,T30,I20,1X,I20,T76,F5.1)') " matmuls inhomo. stacks",&
       stats_num_mnk_stacks(1,4:5), percent

    percent = 100*REAL(SUM(stats_num_mnk_stacks(:,5)))/REAL(MAX(INT(1,KIND=int_8),SUM(stats_num_mnk_stacks(:,4:5))))
    WRITE (output_unit,'(A,T30,I20,1X,I20,T76,F5.1)') " matmuls total",&
       SUM(stats_num_mnk_stacks(:,4)), SUM(stats_num_mnk_stacks(:,5)), percent

    !sorting stat entries by flops per multiplication
    ALLOCATE(sort_key(SIZE(stats_num_mnk_stacks,1)-1))
    sort_key = 2*PRODUCT(stats_num_mnk_stacks(2:,1:3), DIM=2)
    ALLOCATE(sort_idx(SIZE(sort_key)))
    CALL sort(sort_key, SIZE(sort_key), sort_idx)

    DO i=1, SIZE(sort_idx)
     j = sort_idx(i)+1
     percent = 100*REAL(stats_num_mnk_stacks(j,5))/REAL(MAX(INT(1,KIND=int_8),SUM(stats_num_mnk_stacks(j,4:5))))
     WRITE (output_unit,"(A,I4,' x ',I4,' x ',I4,T30,I20,1X,I20,T76,F5.1)") &
       " flops",stats_num_mnk_stacks(j,1:3),stats_num_mnk_stacks(j,4:5)*2*PRODUCT(stats_num_mnk_stacks(j,1:3)), percent
    END DO

    percent = 100*REAL(report%acc_flop)/REAL(MAX(INT(1,KIND=int_8),report%cpu_flop+report%acc_flop))
    WRITE (output_unit,'(A,T30,I20,1X,I20,T76,F5.1)') " flops total",&
       report%cpu_flop, report%acc_flop, percent

  END SUBROUTINE stats_print_report

END MODULE dbcsr_mm_sched


