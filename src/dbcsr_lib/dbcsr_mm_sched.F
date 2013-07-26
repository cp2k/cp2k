!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Fourth layer of the dbcsr matrix-matrix multiplication.
!>         It hides the differences between performing calculations on the
!>         GPU or the CPU.
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
  USE dbcsr_config,                    ONLY: mm_driver,&
                                             mm_driver_cuda
  USE dbcsr_data_methods,              ONLY: dbcsr_data_ensure_size,&
                                             dbcsr_data_get_size,&
                                             dbcsr_data_get_type
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: int_4,&
                                             int_8
  USE dbcsr_message_passing,           ONLY: mp_abort
  USE dbcsr_mm_cudadrv,                ONLY: dbcsr_mm_cudadrv_finalize,&
                                             dbcsr_mm_cudadrv_init,&
                                             dbcsr_mm_cudadrv_lib_finalize,&
                                             dbcsr_mm_cudadrv_lib_init,&
                                             dbcsr_mm_cudadrv_process,&
                                             dbcsr_mm_cudadrv_type
  USE dbcsr_mm_hostdrv,                ONLY: dbcsr_mm_hostdrv_finalize,&
                                             dbcsr_mm_hostdrv_init,&
                                             dbcsr_mm_hostdrv_process,&
                                             dbcsr_mm_hostdrv_type
  USE dbcsr_mm_types,                  ONLY: p_a_first,&
                                             p_b_first,&
                                             p_c_first,&
                                             p_k,&
                                             p_m,&
                                             p_n,&
                                             stack_descriptor_type
  USE dbcsr_types,                     ONLY: dbcsr_type,&
                                             dbcsr_work_type

 !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_sched'

  PUBLIC :: dbcsr_mm_sched_type
  PUBLIC :: dbcsr_mm_sched_lib_init, dbcsr_mm_sched_lib_finalize
  PUBLIC :: dbcsr_mm_sched_init, dbcsr_mm_sched_finalize
  PUBLIC :: dbcsr_mm_sched_process

  ! *****************************************************************************
  TYPE dbcsr_mm_sched_type
    PRIVATE
    TYPE(dbcsr_work_type), POINTER  :: product_wm
    TYPE(dbcsr_mm_cudadrv_type) :: cudadrv
    TYPE(dbcsr_mm_hostdrv_type) :: hostdrv
  END TYPE dbcsr_mm_sched_type


! *****************************************************************************
!>  Statistics Counter
! *****************************************************************************
  ! Threadprivate counters to collect statistics
  INTEGER(kind=int_8), SAVE :: stats_num_cpu_stacks, stats_num_gpu_stacks
  !$OMP THREADPRIVATE(stats_num_cpu_stacks, stats_num_gpu_stacks)
  INTEGER(kind=int_8), SAVE :: stats_gpu_flop, stats_cpu_flop
  !$OMP THREADPRIVATE(stats_gpu_flop, stats_cpu_flop)

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
    stats_num_cpu_stacks = 0
    stats_num_gpu_stacks = 0
    stats_gpu_flop = 0
    stats_cpu_flop = 0
    ALLOCATE(stats_num_mnk_stacks(1,5))
    stats_num_mnk_stacks(1,:) = 0 ! entry for the default stack

     IF (mm_driver == mm_driver_cuda) &
        CALL dbcsr_mm_cudadrv_lib_init(error)


  END SUBROUTINE dbcsr_mm_sched_lib_init


! *****************************************************************************
! \brief Finalize the library and prints DBCSR statistics
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_sched_lib_finalize(output_unit, error)
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: i, ithread, j, nthreads
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER, SAVE                          :: handover
    INTEGER(kind=int_8), SAVE :: stats_cpu_flop_sum, stats_gpu_flop_sum, &
      stats_num_cpu_stacks_sum, stats_num_gpu_stacks_sum
    REAL                                     :: percent

   nthreads = 1; ithread = 0
   !$ nthreads = OMP_GET_NUM_THREADS () ; ithread = OMP_GET_THREAD_NUM ()

   IF (mm_driver == mm_driver_cuda) &
       CALL dbcsr_mm_cudadrv_lib_finalize(output_unit, error)

   ! Collect and output statistics ---------------------------------------------
   !$OMP MASTER
   stats_num_cpu_stacks_sum = 0
   stats_num_gpu_stacks_sum = 0
   stats_gpu_flop_sum = 0
   stats_cpu_flop_sum = 0
   !$OMP END MASTER
   !$OMP BARRIER

   !$OMP CRITICAL
   stats_num_cpu_stacks_sum = stats_num_cpu_stacks_sum + stats_num_cpu_stacks
   stats_num_gpu_stacks_sum = stats_num_gpu_stacks_sum + stats_num_gpu_stacks
   stats_gpu_flop_sum = stats_gpu_flop_sum + stats_gpu_flop
   stats_cpu_flop_sum = stats_cpu_flop_sum + stats_cpu_flop
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
                        stacksize_gpu=handover(j,5))
       END DO
     END IF
   END DO

   !$OMP MASTER
   IF(output_unit>0) THEN
      WRITE (output_unit,"(X,A,T47,A,T68,A,T77,A)") "COUNTER", "CPU", "GPU", "GPU%"

      percent = 100*REAL(stats_num_gpu_stacks_sum)/REAL(MAX(1,stats_num_cpu_stacks_sum+stats_num_gpu_stacks_sum))
      WRITE (output_unit,'(A,T30,I20,X,I20,T76,F5.1)') " number of processed stacks",&
         stats_num_cpu_stacks_sum, stats_num_gpu_stacks_sum, percent

      percent = 100*REAL(stats_num_mnk_stacks(1,5))/REAL(MAX(1,SUM(stats_num_mnk_stacks(1,4:5))))
      WRITE (output_unit,'(A,T30,I20,X,I20,T76,F5.1)') " mutmuls inhomo. stacks",&
         stats_num_mnk_stacks(1,4:5), percent

      DO i=2, SIZE(stats_num_mnk_stacks, 1)
       percent = 100*REAL(stats_num_mnk_stacks(i,5))/REAL(MAX(1,SUM(stats_num_mnk_stacks(i,4:5))))
       WRITE (output_unit,"(A,I4,' x ',I4,' x ',I4,T30,I20,X,I20,T76,F5.1)") &
         " matmuls",stats_num_mnk_stacks(i,:), percent
      END DO

      percent = 100*REAL(SUM(stats_num_mnk_stacks(:,5)))/REAL(MAX(1,SUM(stats_num_mnk_stacks(:,4:5))))
      WRITE (output_unit,'(A,T30,I20,X,I20,T76,F5.1)') " matmuls total",&
         SUM(stats_num_mnk_stacks(:,4)), SUM(stats_num_mnk_stacks(:,5)), percent

      percent = 100*REAL(stats_gpu_flop_sum)/REAL(MAX(1,stats_cpu_flop_sum+stats_gpu_flop_sum))
      WRITE (output_unit,'(A,T30,I20,X,I20,T76,F5.1)') " flops total",&
         stats_cpu_flop_sum, stats_gpu_flop_sum, percent

   END IF
   DEALLOCATE(stats_num_mnk_stacks)
   !$OMP END MASTER
  END SUBROUTINE dbcsr_mm_sched_lib_finalize




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

    INTEGER                                  :: allocated_datasize, &
                                                data_type, error_handler, &
                                                ithread, nthreads, &
                                                used_datasize

    CALL dbcsr_error_set(routineN, error_handler, error)


    nthreads = 1; ithread = 0
    !$ nthreads = OMP_GET_NUM_THREADS () ; ithread = OMP_GET_THREAD_NUM ()

    data_type = dbcsr_data_get_type(left%data_area)

    this%product_wm => product_wm
    used_datasize = this%product_wm%datasize
    allocated_datasize = dbcsr_data_get_size(this%product_wm%data_area)

    ! The product's data_area could already contain some data.
    ! ( see: keep_product_data in dbcsr_operations.F )
    ! But this data might not occupy all the allocated memory in the data_area.
    ! Since, we don't want to keep track of unitialized memory we just zero it now.
    CALL dbcsr_data_clear(this%product_wm%data_area, lb=used_datasize+1, ub=allocated_datasize)

    CALL dbcsr_mm_hostdrv_init(this%hostdrv, left, right, product_wm, error)

    IF (mm_driver == mm_driver_cuda) &
       CALL dbcsr_mm_cudadrv_init(this%cudadrv, left, right, product_wm, error)

    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_mm_sched_init


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

    CALL dbcsr_mm_hostdrv_finalize(this%hostdrv, error)
    IF (mm_driver == mm_driver_cuda) &
       CALL dbcsr_mm_cudadrv_finalize(this%cudadrv, error)

    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_mm_sched_finalize


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

    REAL, PARAMETER                          :: resize_factor = 1.618034

    INTEGER                                  :: sp, stacked_datasize
    INTEGER(kind=int_8)                      :: flop
    LOGICAL                                  :: success

    IF(stack_fillcount <= 0) CALL mp_abort("dbcsr_mm_sched_process: got empty stack")

    stacked_datasize   = this%product_wm%datasize
    CALL dbcsr_data_ensure_size(this%product_wm%data_area, stacked_datasize, &
                    factor=resize_factor, zero_pad=.TRUE., error=error)

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
           IF(stack_data(p_m,sp) /= stack_descr%m) CALL abort()
           IF(stack_data(p_n,sp) /= stack_descr%n) CALL abort()
           IF(stack_data(p_k,sp) /= stack_descr%k) CALL abort()
         END DO
       ENDIF
     ENDIF

     ! Submitting the stack for processing -------------------------------------

     flop = 2*stack_fillcount*stack_descr%max_m*stack_descr%max_n*stack_descr%max_k

     IF (mm_driver == mm_driver_cuda) THEN
         CALL dbcsr_mm_cudadrv_process(&
                    this%cudadrv,&
                    left, right,&
                    params=stack_data,&
                    stack_size=stack_fillcount, &
                    stack_descr=stack_descr,&
                    success=success,&
                    error=error)

         IF(success) THEN
            ! update statistics
            stats_num_gpu_stacks = stats_num_gpu_stacks + 1
            stats_gpu_flop = stats_gpu_flop + flop
            CALL stats_add(m=stack_descr%m, n=stack_descr%n, k=stack_descr%k, &
                     stacksize_gpu=INT(stack_fillcount,kind=int_8))
            RETURN
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
        stats_num_cpu_stacks = stats_num_cpu_stacks + 1
        stats_cpu_flop = stats_cpu_flop + flop
        CALL stats_add(m=stack_descr%m, n=stack_descr%n, k=stack_descr%k, &
                 stacksize_cpu=INT(stack_fillcount,kind=int_8))
        RETURN
     ENDIF

     STOP "dbcsr_mm_sched_process_stack failed"


  END SUBROUTINE dbcsr_mm_sched_process

 ! *****************************************************************************
 ! \brief Helper-routine used by dbcsr_mm_sched_process to supply statistics.
 ! *****************************************************************************
 SUBROUTINE stats_add(m,n,k, stacksize_cpu, stacksize_gpu)
    INTEGER, INTENT(IN)                      :: m, n, k
    INTEGER(kind=int_8), OPTIONAL            :: stacksize_cpu, stacksize_gpu

    INTEGER                                  :: i, s
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tmp

   DO i=1, SIZE(stats_num_mnk_stacks, 1)
       IF(stats_num_mnk_stacks(i,1)==m .AND. &
          stats_num_mnk_stacks(i,2)==n .AND. &
          stats_num_mnk_stacks(i,3)==k) THEN
            IF(PRESENT(stacksize_cpu))&
              stats_num_mnk_stacks(i,4) = stats_num_mnk_stacks(i,4) + stacksize_cpu
            IF(PRESENT(stacksize_gpu))&
              stats_num_mnk_stacks(i,5) = stats_num_mnk_stacks(i,5) + stacksize_gpu
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
   IF(PRESENT(stacksize_gpu))&
     stats_num_mnk_stacks(s+1,5) = stacksize_gpu
   DEALLOCATE(tmp)
 END SUBROUTINE stats_add
END MODULE dbcsr_mm_sched


