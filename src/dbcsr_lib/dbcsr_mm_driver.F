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


MODULE dbcsr_mm_driver
  USE dbcsr_block_operations,          ONLY: block_add,&
                                             dbcsr_data_clear
  USE dbcsr_config,                    ONLY: dbcsr_get_conf_nstacks,&
                                             mm_driver,&
                                             mm_driver_cuda,&
                                             mm_stack_size
  USE dbcsr_cuda_device,               ONLY: dbcsr_cuda_create_streams,&
                                             dbcsr_cuda_destroy_streams,&
                                             dbcsr_cuda_device_sync
  USE dbcsr_cuda_memory,               ONLY: dbcsr_cuda_dev_mem_alloc,&
                                             dbcsr_cuda_dev_mem_dealloc,&
                                             dbcsr_cuda_dev_mem_new,&
                                             dbcsr_cuda_dev_mem_realloc,&
                                             dbcsr_cuda_dev_mem_release,&
                                             dbcsr_cuda_dev_mem_zero,&
                                             dbcsr_cuda_host_mem_alloc,&
                                             dbcsr_cuda_host_mem_dealloc
  USE dbcsr_cuda_methods,              ONLY: dbcsr_cuda_dev_mem_get_alloc
  USE dbcsr_cuda_operations,           ONLY: dbcsr_cuda_cp_dev_to_host,&
                                             dbcsr_cuda_cp_host_to_dev
  USE dbcsr_cuda_types,                ONLY: dbcsr_cuda_mem_p_type,&
                                             dbcsr_cuda_mem_type
  USE dbcsr_data_methods,              ONLY: dbcsr_data_ensure_size,&
                                             dbcsr_data_get_size,&
                                             dbcsr_data_get_type,&
                                             dbcsr_data_init,&
                                             dbcsr_data_new,&
                                             dbcsr_data_release,&
                                             dbcsr_data_set_size_referenced
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: dp,&
                                             int_8,&
                                             real_8
  USE dbcsr_mm_stack,                  ONLY: process_mm_stack_cpu,&
                                             process_mm_stack_cuda
  USE dbcsr_mm_types,                  ONLY: dbcsr_ps_state_empty,&
                                             dbcsr_ps_state_working,&
                                             dbcsr_ps_width,&
                                             driver_carrier_type,&
                                             stack_descriptor_type
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_type,&
                                             dbcsr_type_int_4

 !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_driver'

  PUBLIC :: dbcsr_mm_driver_lib_init,   dbcsr_mm_driver_lib_finalize
  PUBLIC :: dbcsr_mm_driver_outer_init, dbcsr_mm_driver_outer_finalize
  PUBLIC :: dbcsr_mm_driver_inner_init, dbcsr_mm_driver_inner_finalize
  PUBLIC :: dbcsr_mm_driver_process_stack

  ! ===== Global Cuda Memory =====
  ! Allocating memory for cuda (on the card and host-pinned) is slow.
  ! Therefore, the memory is allocated once and stored in global variables.
  ! Variable cuda_mem_in_use should prevent concurred use of this global memory.
  LOGICAL cuda_mem_in_use
  !$OMP THREADPRIVATE(cuda_mem_in_use)
  TYPE(dbcsr_cuda_mem_type), POINTER :: cuda_dev_c_data, cuda_dev_c_locks
  !$OMP THREADPRIVATE(cuda_dev_c_data, cuda_dev_c_locks)
  TYPE(dbcsr_cuda_mem_p_type), DIMENSION(:), POINTER :: cuda_dev_stacklist_data
  !$OMP THREADPRIVATE(cuda_dev_stacklist_data)
  INTEGER, DIMENSION(:,:,:,:), POINTER     :: cuda_host_stackgrid_data
  !$OMP THREADPRIVATE(cuda_host_stackgrid_data)
  INTEGER, DIMENSION(:,:), POINTER         :: cuda_host_stackgrid_state
  !$OMP THREADPRIVATE(cuda_host_stackgrid_state)

  !These are shared between threads
  TYPE(dbcsr_cuda_mem_type), POINTER, SAVE :: cuda_dev_a_data => Null()
  TYPE(dbcsr_cuda_mem_type), POINTER, SAVE :: cuda_dev_b_data => Null()
  TYPE(dbcsr_cuda_mem_type), POINTER, SAVE :: cuda_dev_state_empty => Null()


  !INTEGER , SAVE :: num_cpu_stacks = 0
  !INTEGER , SAVE :: num_gpu_stacks = 0
  !INTEGER(kind=int_8), SAVE   :: gpu_flop = 0
  !INTEGER(kind=int_8), SAVE   :: cpu_flop = 0

  CONTAINS


! *****************************************************************************
! \brief Initialize the library
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_driver_lib_init(error)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: i, ithread, n_stack_buffers, &
                                                nstacks, nthreads
    INTEGER, DIMENSION(1), TARGET            :: tmp_arr
    INTEGER, DIMENSION(3)                    :: nxstacks
    INTEGER, DIMENSION(:), POINTER           :: tmp_ptr

    nthreads = 1; ithread = 0
    !$ nthreads = OMP_GET_NUM_THREADS () ; ithread = OMP_GET_THREAD_NUM ()

     IF (mm_driver == mm_driver_cuda) THEN
         CALL dbcsr_get_conf_nstacks(nxstacks, n_stack_buffers, error)
         nstacks = nxstacks(1) * nxstacks(2) * nxstacks(3) + 1

         !Ensure that lock_offset in dbcsr_cuda_calc.cu won't overflow (it's a int4).
         ! TODO: Check this code, before you start to used locking again.
         IF(4*mm_stack_size > 2**17)& ! The "4" due to stack_mm_mnk_sq5_d
            STOP "dbcsr_mm_driver_lib_init: mm_stack_size is too large"
         IF(nthreads*nstacks*n_stack_buffers > 2**15)&
            STOP "dbcsr_mm_driver_lib_init: nthreads*nstacks*n_stack_buffers is too large"

        ! 1. Things that are shared between threads.
        !-----------------------------------------------------------------------
        !$OMP BARRIER
        !$OMP MASTER
        CALL dbcsr_cuda_create_streams(nthreads*nstacks*n_stack_buffers, error)

        ! Constant dbcsr_ps_state_empty needed on card to signal completed transfers.
        tmp_arr = dbcsr_ps_state_empty
        tmp_ptr =>  tmp_arr
        CALL dbcsr_cuda_dev_mem_new(cuda_dev_state_empty, dbcsr_type_int_4, error=error)
        CALL dbcsr_cuda_dev_mem_alloc(cuda_dev_state_empty, 1, error=error)
        CALL dbcsr_cuda_cp_host_to_dev(host_mem=tmp_ptr, dev_mem=cuda_dev_state_empty%d_i, n=1,&
           async=.FALSE.,error=error)
        !$OMP END MASTER
        !$OMP BARRIER

        ! 2. Things that are local to each thread
        !------------------------------------------------------------------------
        ! Allocate the parameter stack on the card
        ALLOCATE(cuda_dev_stacklist_data(nstacks*n_stack_buffers))
        DO i=1, nstacks*n_stack_buffers
           NULLIFY(cuda_dev_stacklist_data(i)%p)
           CALL dbcsr_cuda_dev_mem_new(cuda_dev_stacklist_data(i)%p, dbcsr_type_int_4, error)
           CALL dbcsr_cuda_dev_mem_alloc(cuda_dev_stacklist_data(i)%p,&
                           dbcsr_ps_width*mm_stack_size, error=error)
        END DO

        ! Allocate the parameter stack on the host with pinned memory
        CALL dbcsr_cuda_host_mem_alloc(cuda_host_stackgrid_data, &
           dbcsr_ps_width, mm_stack_size, nstacks, n_stack_buffers, error=error)
        CALL dbcsr_cuda_host_mem_alloc(cuda_host_stackgrid_state, &
           nstacks, n_stack_buffers, error=error)

        ! Size of the C-Blocks unkown at this point, allocation is done late.
        NULLIFY(cuda_dev_c_data, cuda_dev_c_locks)

        !Setup safety guard against overlapping multiplications.
        cuda_mem_in_use = .FALSE.
    ENDIF

  END SUBROUTINE dbcsr_mm_driver_lib_init


! *****************************************************************************
! \brief Finalize the library
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_driver_lib_finalize(error)
     USE ISO_C_BINDING
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: i

   IF (mm_driver == mm_driver_cuda) THEN
       ! Releasing host pinned memory of stackgrid on the host
       CALL dbcsr_cuda_host_mem_dealloc(cuda_host_stackgrid_data, error)
       CALL dbcsr_cuda_host_mem_dealloc(cuda_host_stackgrid_state, error)

       ! releasing stuff on the card if necessary
       IF(ASSOCIATED(cuda_dev_c_data))&
          CALL dbcsr_cuda_dev_mem_release(cuda_dev_c_data, error=error)
       IF(ASSOCIATED(cuda_dev_c_locks))&
          CALL dbcsr_cuda_dev_mem_release(cuda_dev_c_locks, error=error)

       DO i=1, SIZE(cuda_dev_stacklist_data)
          CALL dbcsr_cuda_dev_mem_release(cuda_dev_stacklist_data(i)%p, error=error)
       END DO
       DEALLOCATE(cuda_dev_stacklist_data)

       !$OMP BARRIER
       !$OMP MASTER
       IF(ASSOCIATED(cuda_dev_a_data))&
          CALL dbcsr_cuda_dev_mem_release(cuda_dev_a_data, error=error)
       IF(ASSOCIATED(cuda_dev_b_data))&
          CALL dbcsr_cuda_dev_mem_release(cuda_dev_b_data, error=error)
       CALL dbcsr_cuda_dev_mem_release(cuda_dev_state_empty, error=error)
       CALL dbcsr_cuda_destroy_streams(error=error)

       !How much memory is still allocated on the card?
       !istat = dbcsr_cuda_dev_mem_info(mem_free, mem_avail)
       !WRITE (*,*) "after outer_finalize cuda mem: ",mem_free, mem_avail, istat
       !$OMP END MASTER
       !$OMP BARRIER

   END IF


   ! !$OMP MASTER
   !  WRITE (*,*) "num_stacks gpu:", num_gpu_stacks, " cpu: ",num_cpu_stacks, " ratio: ",REAL(num_gpu_stacks)/REAL(num_cpu_stacks)
   !  WRITE (*,*) "flop gpu:", gpu_flop, " cpu: ",cpu_flop, " ratio: ",REAL(gpu_flop)/REAL(cpu_flop)
   ! !$OMP END MASTER
  END SUBROUTINE dbcsr_mm_driver_lib_finalize




! *****************************************************************************
!> \brief Initializes a multiplication cycle for new set of C-blocks.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_driver_outer_init(left, right, product, drv_carrier,&
               fill_guess, nstacks, n_stack_buffers, stackgrid_data, stackgrid_state, error)
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_type), INTENT(INOUT)          :: product
    TYPE(driver_carrier_type), INTENT(INOUT) :: drv_carrier
    REAL(KIND=real_8)                        :: fill_guess
    INTEGER, INTENT(IN)                      :: nstacks, n_stack_buffers
    INTEGER, DIMENSION(:, :, :, :), POINTER  :: stackgrid_data
    INTEGER, DIMENSION(:, :), POINTER        :: stackgrid_state
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_driver_outer_init'

    INTEGER :: a_data_size, allocated_datasize, b_data_size, c_data_size, &
      c_locks_size, data_type, error_handler, ithread, nthreads, used_datasize

    CALL dbcsr_error_set(routineN, error_handler, error)

    nthreads = 1; ithread = 0
    !$ nthreads = OMP_GET_NUM_THREADS () ; ithread = OMP_GET_THREAD_NUM ()

    data_type = dbcsr_data_get_type(left%data_area)

    drv_carrier%product_wm => product%wms(ithread+1) !TODO: why do we need this ref?
    used_datasize = drv_carrier%product_wm%datasize
    allocated_datasize = dbcsr_data_get_size(drv_carrier%product_wm%data_area)

    ! The product's data_area could already contain some data.
    ! ( see: keep_product_data in dbcsr_operations.F )
    ! But this data might not occupy all the allocated memory in the data_area.
    ! Since, we don't want to keep track of unitialized memory we just zero it now.
    CALL dbcsr_data_clear(drv_carrier%product_wm%data_area, lb=used_datasize+1, ub=allocated_datasize)

    IF (mm_driver /= mm_driver_cuda) THEN
       ALLOCATE(stackgrid_data(dbcsr_ps_width, mm_stack_size, nstacks, n_stack_buffers))
       ALLOCATE(stackgrid_state(nstacks, n_stack_buffers))
    ELSE
       ! 0. Making shure there is only one multiplication going on at a time.
       !------------------------------------------------------------------------
       IF(.NOT. ASSOCIATED(cuda_dev_state_empty)) &
            STOP "dbcsr_mm_driver_outer_init: dbcsr not initialized"
       IF(cuda_mem_in_use) STOP "dbcsr_mm_driver: found overlapping multiplications"
       cuda_mem_in_use = .TRUE.

       ! 1. Things that are shared between threads.
       !------------------------------------------------------------------------
       !$OMP BARRIER
       !$OMP MASTER
       ! (re)Allocate Memory for A and B blocks on the Card
       a_data_size = dbcsr_data_get_size(left%data_area)
       b_data_size = dbcsr_data_get_size(right%data_area)
       CALL cuda_prepare_mem(cuda_dev_a_data, data_type, a_data_size, error=error)
       CALL cuda_prepare_mem(cuda_dev_b_data, data_type, b_data_size, error=error)
       !$OMP END MASTER
       !$OMP BARRIER


       ! 2. Things that are local to each thread
       !------------------------------------------------------------------------
       IF(SIZE(cuda_host_stackgrid_data, 1) /= dbcsr_ps_width  .OR. &
          SIZE(cuda_host_stackgrid_data, 2) /= mm_stack_size   .OR. &
          SIZE(cuda_host_stackgrid_data, 3) /= nstacks         .OR. &
          SIZE(cuda_host_stackgrid_data, 4) /= n_stack_buffers ) &
             STOP "dbcsr_mm_driver_outer_init: cuda_host_stackgrid_data wrong dims"

       IF(SIZE(cuda_host_stackgrid_state, 1) /= nstacks  .OR. &
          SIZE(cuda_host_stackgrid_state, 2) /= n_stack_buffers ) &
             STOP "dbcsr_mm_driver_outer_init: cuda_host_stackgrid_state wrong dims"

       stackgrid_data  => cuda_host_stackgrid_data
       stackgrid_state => cuda_host_stackgrid_state

       ! (re)Allocate Memory for C blocks on the Card and zero it.
       c_data_size = dbcsr_data_get_size(product%wms(ithread+1)%data_area)
       CALL cuda_prepare_mem(cuda_dev_c_data, data_type, c_data_size, zero=.TRUE., error=error)

       ! (re)Allocate and zero C locks on GPU.
       c_locks_size = INT(product%nblkcols_local,KIND=int_8)* &
                  INT(product%nblkrows_local,KIND=int_8) * fill_guess + 1
       CALL cuda_prepare_mem(cuda_dev_c_locks, dbcsr_type_int_4, c_locks_size, zero=.TRUE., error=error)

    ENDIF

    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_mm_driver_outer_init


! *****************************************************************************
! \brief Helper routine used by dbcsr_mm_driver_outer_init().
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE cuda_prepare_mem(dev_mem, data_type, data_size, zero, error)
    TYPE(dbcsr_cuda_mem_type), POINTER       :: dev_mem
    INTEGER, INTENT(IN)                      :: data_type, data_size
    LOGICAL, INTENT(IN), OPTIONAL            :: zero
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

     IF (.NOT. ASSOCIATED(dev_mem)) THEN
        !WRITE (*,*) "cuda_prepare_mem: allocting for the first time"
        CALL dbcsr_cuda_dev_mem_new(dev_mem, data_type, error=error)
        CALL dbcsr_cuda_dev_mem_alloc(dev_mem, data_size, error=error)
     END IF

     IF (dev_mem%data_type /= data_type) THEN
        !WRITE (*,*) "cuda_prepare_mem: changing datatype"
        CALL dbcsr_cuda_dev_mem_release(dev_mem, error=error)
        CALL dbcsr_cuda_dev_mem_new(dev_mem, data_type, error=error)
        CALL dbcsr_cuda_dev_mem_alloc(dev_mem, data_size, error=error)
     END IF

     IF (dbcsr_cuda_dev_mem_get_alloc(dev_mem) < data_size) THEN
        !WRITE (*,*) "cuda_prepare_mem: growing dev_mem"
        CALL dbcsr_cuda_dev_mem_dealloc(dev_mem, error=error)
        CALL dbcsr_cuda_dev_mem_alloc(dev_mem, data_size, error=error)
     END IF

     IF(PRESENT(zero)) THEN
        IF(zero) CALL dbcsr_cuda_dev_mem_zero(dev_mem, first=1, &
                    last=dbcsr_cuda_dev_mem_get_alloc(dev_mem), error=error)
     END IF
  END SUBROUTINE


! *****************************************************************************
!> \brief Initializes a multiplication cycle for a new set of A- and B-blocks.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_driver_inner_init(left, right, drv_carrier, error)
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(driver_carrier_type), INTENT(INOUT) :: drv_carrier
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_driver_inner_init'

    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)

    drv_carrier%left_data_area    = left%data_area
    drv_carrier%right_data_area   = right%data_area

    IF (mm_driver == mm_driver_cuda) THEN
       !$OMP BARRIER
       !$OMP MASTER
       !TODO: could be async if we use events to make the other streams wait
       !      for a and b to be uploaded.
       CALL dbcsr_cuda_cp_host_to_dev(host_mem=left%data_area,&
            dev_mem=cuda_dev_a_data, async = .FALSE., error=error)
       CALL dbcsr_cuda_cp_host_to_dev(host_mem=right%data_area,&
            dev_mem=cuda_dev_b_data, async = .FALSE., error=error)
       !$OMP END MASTER
       !$OMP BARRIER
    ENDIF

    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_mm_driver_inner_init


! *****************************************************************************
!> \brief Finalizes a multiplication cycle for a set of A- and B-blocks.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_driver_inner_finalize(drv_carrier, error)
    TYPE(driver_carrier_type), INTENT(INOUT) :: drv_carrier
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    IF (mm_driver == mm_driver_cuda) THEN
       !TODO: waiting for ALL streams - should only wait for those of this thread
       !$OMP BARRIER
       CALL dbcsr_cuda_device_sync(error)
       !$OMP BARRIER
    END IF
  END SUBROUTINE dbcsr_mm_driver_inner_finalize


! *****************************************************************************
!> \brief Finalizes a multiplication cycle for a set of C-blocks.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_driver_outer_finalize(drv_carrier, product,stackgrid_data, stackgrid_state, error)
    USE ISO_C_BINDING
    TYPE(driver_carrier_type), INTENT(INOUT) :: drv_carrier
    TYPE(dbcsr_type), INTENT(inout)          :: product
    INTEGER, DIMENSION(:, :, :, :), POINTER  :: stackgrid_data
    INTEGER, DIMENSION(:, :), POINTER        :: stackgrid_state
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: data_type, datasize, ithread
    TYPE(dbcsr_data_obj)                     :: tmp_data

   ithread = 0
   !$ ithread = omp_get_thread_num()

   IF (mm_driver /= mm_driver_cuda) THEN
      DEALLOCATE(stackgrid_data)
      DEALLOCATE(stackgrid_state)
   ELSE
       ! Host-pinned memory not dealloced, will be reused by next multiplication
       NULLIFY(stackgrid_data, stackgrid_state)

       ! Transfer C-data from GPU to host and adding it to host's result
       datasize = product%wms(ithread+1)%datasize
       data_type = dbcsr_data_get_type(product%wms(ithread+1)%data_area)

       CALL dbcsr_data_init(tmp_data)
       CALL dbcsr_data_new(tmp_data, data_type=data_type, data_size=datasize)
       CALL dbcsr_cuda_cp_dev_to_host(dev_mem=cuda_dev_c_data,&
            host_mem=tmp_data, error=error)
       CALL block_add(product%wms(ithread+1)%data_area, tmp_data,&
            len=datasize, error=error)
       CALL dbcsr_data_release(tmp_data)

       cuda_mem_in_use = .FALSE.
    ENDIF
  END SUBROUTINE dbcsr_mm_driver_outer_finalize


! *****************************************************************************
! \brief Processes a given stack.
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_driver_process_stack(drv_carrier, driver, stack_data,&
                    stack_fillcount, stack_state, stack_descr, stack_id, error)
    TYPE(driver_carrier_type), INTENT(INOUT) :: drv_carrier
    INTEGER, INTENT(IN)                      :: driver
    INTEGER, DIMENSION(:, :), POINTER        :: stack_data
    INTEGER, POINTER                         :: stack_fillcount, stack_state
    TYPE(stack_descriptor_type), POINTER     :: stack_descr
    INTEGER, INTENT(IN)                      :: stack_id
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    REAL, PARAMETER                          :: resize_factor = 1.618034

    INTEGER                                  :: allocated_datasize, nblks, &
                                                new_c_size, new_locks_size, &
                                                old_c_size, stacked_datasize

!INTEGER(kind=int_8)                      :: flop
!TODO: make difference clear between driver, mm_driver and mm_host_driver.
! collect statistics

    stack_descr%total_entries = stack_descr%total_entries + stack_fillcount

     stacked_datasize   = drv_carrier%product_wm%datasize
     allocated_datasize = dbcsr_data_get_size(drv_carrier%product_wm%data_area)

     ! tell the data_area how much of its memory is actually beeing used.
     CALL dbcsr_data_set_size_referenced(drv_carrier%product_wm%data_area, stacked_datasize)

     ! Resizing the product data (on host and card) ---------------------------
     IF (stacked_datasize > allocated_datasize) THEN

       ! Resize the target data area if the stack references a target
       ! matrix data element outside of its current range.  When
       ! appropriate a lock on the data area is taken.
       !WRITE (*,*) "dbcsr_mm_driver_process_stack: We need to resize c data"

       !IF (dbg) &
       !     WRITE(*,*)routineN//" Resizing to", LOG(REAL(maxs))/LOG(10.0)
       CALL dbcsr_data_ensure_size(drv_carrier%product_wm%data_area,&
           stacked_datasize, factor=resize_factor, zero_pad=.TRUE., error=error)

       IF(stacked_datasize > dbcsr_data_get_size(drv_carrier%product_wm%data_area))&
          STOP "dbcsr_mm_driver_process_stack: resizing failed"

       ! Cuda on-device resizing
       IF (mm_driver == mm_driver_cuda) THEN
          ! All kernels in the stream must be finished before
          ! reallocation is performed.
          !TODO: waiting for ALL streams - should only wait for those of this thread
          CALL dbcsr_cuda_device_sync(error)

          ! dbcsr_data_ensure_size oversizes the array - we need to know how much.
          old_c_size = dbcsr_cuda_dev_mem_get_alloc(cuda_dev_c_data)
          new_c_size = dbcsr_data_get_size(drv_carrier%product_wm%data_area)
          CALL dbcsr_cuda_dev_mem_realloc(cuda_dev_c_data, new_c_size, error=error)
          CALL dbcsr_cuda_dev_mem_zero(cuda_dev_c_data,&
               first=old_c_size, last=new_c_size, error=error)
        END IF
     END IF

     ! Resizing the product locks (only on the card) --------------------------
     IF (mm_driver == mm_driver_cuda) THEN
       ! Resize locks, which are equal to the new block count.
       nblks = drv_carrier%product_wm%lastblk 
       IF (nblks > dbcsr_cuda_dev_mem_get_alloc(cuda_dev_c_locks)) THEN
          !WRITE (*,*) "dbcsr_mm_driver_process_stack: We need to resize c_locks_cuda"
          ! All kernels in the stream must be finished before
          ! reallocation is performed.
          !TODO: waiting for ALL streams - should only wait for those of this thread
          CALL dbcsr_cuda_device_sync(error)
          new_locks_size = INT(REAL(nblks*4,kind=dp)*resize_factor)
          CALL dbcsr_cuda_dev_mem_realloc(cuda_dev_c_locks, new_locks_size, error=error)
          CALL dbcsr_cuda_dev_mem_zero(cuda_dev_c_locks,&
                                     first=1, last=new_locks_size, error=error)
        END IF
     END IF

     !!From here on there is no boundary checking due to assumed-SIZE-arguments.
     !!This is usefull to check stack parameters, BUT it works only for kind=dp 
     !DO sp = 1, stack_fillcount
     !   IF(stack_data(p_a_first,sp) > SIZE(drv_carrier%left_data_area%d%r_dp)) &
     !      STOP "left data out of range"
     !   IF(stack_data(p_b_first,sp) > SIZE(drv_carrier%right_data_area%d%r_dp)) &
     !      STOP "right data out of range"
     !   IF(stack_data(p_c_first,sp) > SIZE(drv_carrier%product_wm%data_area%d%r_dp)) THEN
     !      WRITE (*,*) "blub: ",stack_data(p_c_first,sp) , SIZE(drv_carrier%product_wm%data_area%d%r_dp), &
     !        dbcsr_data_get_size(drv_carrier%product_wm%data_area), stacked_datasize
     !      STOP "product data out of range"
     !   END IF
     !END DO

     ! Submitting the stack for processing -------------------------------------
     stack_state = dbcsr_ps_state_working

     ! flop = 2*stack_fillcount*stack_descr%max_m*stack_descr%max_n*stack_descr%max_k
     IF (driver == mm_driver_cuda) THEN
        ! !$OMP critical
        !    num_gpu_stacks = num_gpu_stacks + 1
        !    gpu_flop = gpu_flop + flop
        ! !$OMP end critical
        CALL process_mm_stack_cuda(params=stack_data,&
                                   stack_size=stack_fillcount,&
                                   data_a_dev=cuda_dev_a_data,&
                                   data_b_dev=cuda_dev_b_data,&
                                   data_c_dev=cuda_dev_c_data,&
                                   c_locks=cuda_dev_c_locks,&
                                   params_dev=cuda_dev_stacklist_data(stack_id)%p,&
                                   m=stack_descr%m,&
                                   n=stack_descr%n,&
                                   k=stack_descr%k,&
                                   max_m=stack_descr%max_m,&
                                   max_n=stack_descr%max_n,&
                                   max_k=stack_descr%max_k,&
                                   defined_mnk=stack_descr%defined_mnk,&
                                   state=stack_state,&
                                   stack_state_dev=cuda_dev_state_empty,&
                                   stack_id=stack_id,&
                                   error=error)
      ELSE
        ! !$OMP critical
        !    num_cpu_stacks = num_cpu_stacks + 1
        !    cpu_flop = cpu_flop + flop
        ! !$OMP end critical
        CALL process_mm_stack_cpu(params=stack_data,&
                                  driver=driver,&
                                  stack_size=stack_fillcount, &
                                  left_data_area=drv_carrier%left_data_area,&
                                  right_data_area=drv_carrier%right_data_area,&
                                  product_data_area=drv_carrier%product_wm%data_area,&
                                  state=stack_state,&
                                  m=stack_descr%m,&
                                  n=stack_descr%n,&
                                  k=stack_descr%k,&
                                  max_m=stack_descr%max_m,&
                                  max_n=stack_descr%max_n,&
                                  max_k=stack_descr%max_k,&
                                  defined_mnk=stack_descr%defined_mnk,&
                                  stack_id=stack_id,&
                                  error=error)
     END IF

     stack_fillcount = 0 !otherwise the stackcolumn gets processed again.


  END SUBROUTINE dbcsr_mm_driver_process_stack


END MODULE dbcsr_mm_driver

