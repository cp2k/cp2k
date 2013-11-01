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


MODULE dbcsr_mm_cudadrv
  USE ISO_C_BINDING
  USE dbcsr_block_operations,          ONLY: block_add
  USE dbcsr_config,                    ONLY: cudadrv_binning_binsize,&
                                             cudadrv_binning_nbins,&
                                             cudadrv_min_flop_sort,&
                                             cudadrv_nbuffers,&
                                             default_resize_factor,&
                                             mm_stack_size
  USE dbcsr_cuda_device,               ONLY: cudaSharedMemBankSizeEightByte,&
                                             dbcsr_cuda_set_shared_mem_config
  USE dbcsr_cuda_devmem,               ONLY: dbcsr_cuda_devmem_allocate,&
                                             dbcsr_cuda_devmem_allocated,&
                                             dbcsr_cuda_devmem_deallocate,&
                                             dbcsr_cuda_devmem_host2dev,&
                                             dbcsr_cuda_devmem_setzero,&
                                             dbcsr_cuda_devmem_type
  USE dbcsr_cuda_event,                ONLY: dbcsr_cuda_event_create,&
                                             dbcsr_cuda_event_destroy,&
                                             dbcsr_cuda_event_query,&
                                             dbcsr_cuda_event_record,&
                                             dbcsr_cuda_event_type,&
                                             dbcsr_cuda_stream_wait_event
  USE dbcsr_cuda_hostmem,              ONLY: dbcsr_cuda_hostmem_allocate,&
                                             dbcsr_cuda_hostmem_deallocate
  USE dbcsr_cuda_operations,           ONLY: dbcsr_cuda_do_mm_stack
  USE dbcsr_cuda_stream,               ONLY: dbcsr_cuda_stream_create,&
                                             dbcsr_cuda_stream_destroy,&
                                             dbcsr_cuda_stream_synchronize,&
                                             dbcsr_cuda_stream_type
  USE dbcsr_data_methods,              ONLY: dbcsr_data_dev2host,&
                                             dbcsr_data_ensure_size,&
                                             dbcsr_data_get_size,&
                                             dbcsr_data_get_type,&
                                             dbcsr_data_new,&
                                             dbcsr_data_release
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: default_string_length,&
                                             int_8
  USE dbcsr_mem_methods,               ONLY: dbcsr_mempool_destruct,&
                                             dbcsr_memtype_setup
  USE dbcsr_message_passing,           ONLY: mp_abort
  USE dbcsr_mm_types,                  ONLY: dbcsr_ps_width,&
                                             stack_descriptor_type
  USE dbcsr_toollib,                   ONLY: sort
  USE dbcsr_types,                     ONLY: dbcsr_data_area_type,&
                                             dbcsr_data_obj,&
                                             dbcsr_memtype_type,&
                                             dbcsr_type,&
                                             dbcsr_type_int_4,&
                                             dbcsr_work_type

 !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_cudadrv'

  LOGICAL, PARAMETER, PRIVATE :: careful_mod = .FALSE.


! *****************************************************************************
  TYPE dbcsr_mm_cudadrv_type
     PRIVATE
     TYPE(dbcsr_work_type), POINTER           :: product_wm => Null()
     TYPE(dbcsr_data_obj)                     :: c_buffer
  END TYPE dbcsr_mm_cudadrv_type

! *****************************************************************************

  PUBLIC :: dbcsr_mm_cudadrv_type
  PUBLIC :: dbcsr_mm_cudadrv_lib_init,   dbcsr_mm_cudadrv_lib_finalize
  PUBLIC :: dbcsr_mm_cudadrv_init, dbcsr_mm_cudadrv_finalize
  PUBLIC :: dbcsr_mm_cudadrv_process

  ! ===== Global Cuda Memory =====
  ! Allocating memory for cuda (on the card and host-pinned) is slow.
  ! Therefore, the memory is allocated once and stored in global variables.

  TYPE stack_buffer_type
     TYPE(dbcsr_cuda_devmem_type)           :: devmem
     INTEGER, DIMENSION(:,:), POINTER       :: hostmem
     TYPE(dbcsr_cuda_event_type)            :: ready
     TYPE(dbcsr_cuda_event_type)            :: calculated
     TYPE(dbcsr_cuda_stream_type)           :: stream
  END TYPE stack_buffer_type

  TYPE(stack_buffer_type), DIMENSION(:), ALLOCATABLE, TARGET  :: stack_buffers
  INTEGER :: curr_stack_buffer_idx = 0
  !$OMP THREADPRIVATE(stack_buffers, curr_stack_buffer_idx)

  TYPE(dbcsr_memtype_type), SAVE, PRIVATE :: memtype_cbuffer
  !$OMP THREADPRIVATE(memtype_cbuffer)

  CONTAINS


! *****************************************************************************
! \brief Initialize the library
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_cudadrv_lib_init(error)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

! nothing to do here

  END SUBROUTINE dbcsr_mm_cudadrv_lib_init


! *****************************************************************************
! \brief Finalize the library
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_cudadrv_lib_finalize(output_unit, error)
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

   IF(ALLOCATED(stack_buffers)) &
     CALL deallocate_stackbuffers(error)

   IF(ASSOCIATED(memtype_cbuffer%pool))&
     CALL dbcsr_mempool_destruct(memtype_cbuffer%pool, error)

   !How much memory is still allocated on the card?
   !istat = dbcsr_cuda_dev_mem_info(mem_free, mem_avail)
   !WRITE (*,*) "after finalize cuda mem: ",mem_free, mem_avail, istat

  END SUBROUTINE dbcsr_mm_cudadrv_lib_finalize


! *****************************************************************************
!> \brief Initializes a multiplication cycle for new set of C-blocks.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_cudadrv_init(this, left, right, product_wm, error)
    TYPE(dbcsr_mm_cudadrv_type), &
      INTENT(INOUT)                          :: this
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_work_type), POINTER           :: product_wm
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_cudadrv_init', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: poolname
    INTEGER                                  :: error_handler, i, ithread, &
                                                nthreads
    TYPE(dbcsr_data_obj)                     :: c_area

    CALL dbcsr_error_set(routineN, error_handler, error)

    nthreads = 1; ithread = 0
    !$ nthreads = OMP_GET_NUM_THREADS () ; ithread = OMP_GET_THREAD_NUM ()

    ! Setup global data which is reused in between multiplications
    !------------------------------------------------------------------------
    !Each thread has its own memtype with its own mempool
     WRITE (poolname, "(A,I4,A)") "cudadrv_cbuffer (thread=", ithread,")"
     CALL dbcsr_memtype_setup(memtype_cbuffer, cuda_hostalloc=.TRUE.,&
          cuda_devalloc=.TRUE., pool_capacity=1, pool_name=poolname, error=error)

    ! Allocate the buffers for transfering the stack_data
    IF(ALLOCATED(stack_buffers)) THEN
       IF(SIZE(stack_buffers)/=cudadrv_nbuffers) &
          CALL deallocate_stackbuffers(error)
    ENDIF

    IF(.NOT.ALLOCATED(stack_buffers)) THEN
       ALLOCATE(stack_buffers(cudadrv_nbuffers))
       DO i=1, cudadrv_nbuffers
          CALL dbcsr_cuda_devmem_allocate(stack_buffers(i)%devmem,&
                   dbcsr_type_int_4, dbcsr_ps_width*mm_stack_size, error)
          CALL dbcsr_cuda_hostmem_allocate(stack_buffers(i)%hostmem, &
                   dbcsr_ps_width, mm_stack_size, error=error)
          CALL dbcsr_cuda_event_create(stack_buffers(i)%ready, error)
          CALL dbcsr_cuda_event_create(stack_buffers(i)%calculated, error)
          CALL dbcsr_cuda_stream_create(stack_buffers(i)%stream, error)
       END DO
    END IF

    CALL dbcsr_cuda_set_shared_mem_config(cudaSharedMemBankSizeEightByte)

    ! Setup things for this particular multipcation
    !------------------------------------------------------------------------
    this%product_wm => product_wm
    c_area = this%product_wm%data_area

    CALL dbcsr_data_new(this%c_buffer, data_type=dbcsr_data_get_type(c_area), &
          data_size=dbcsr_data_get_size(c_area), memory_type=memtype_cbuffer)

    CALL dbcsr_cuda_devmem_setzero(this%c_buffer%d%cuda_devmem,&
                             stream=this%c_buffer%d%cuda_stream, error=error)

    CALL dbcsr_cuda_event_record(this%c_buffer%d%cuda_ready,&
                             stream=this%c_buffer%d%cuda_stream, error=error)

    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_mm_cudadrv_init


! *****************************************************************************
! \brief Helper routine used by dbcsr_mm_cudadrv_init() and dbcsr_mm_cudadrv_lib_finalize()
! \author Ole Schuett
! *****************************************************************************
 SUBROUTINE deallocate_stackbuffers(error)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: i

   DO i=1, SIZE(stack_buffers)
      CALL dbcsr_cuda_devmem_deallocate(stack_buffers(i)%devmem, error)
      CALL dbcsr_cuda_hostmem_deallocate(stack_buffers(i)%hostmem, error)
      CALL dbcsr_cuda_event_destroy(stack_buffers(i)%ready, error)
      CALL dbcsr_cuda_event_destroy(stack_buffers(i)%calculated, error)
      CALL dbcsr_cuda_stream_destroy(stack_buffers(i)%stream, error)
   END DO
   DEALLOCATE(stack_buffers)
 END SUBROUTINE deallocate_stackbuffers


! *****************************************************************************
!> \brief Finalizes a multiplication cycle for a set of C-blocks.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_cudadrv_finalize(this, error)
    USE ISO_C_BINDING
    TYPE(dbcsr_mm_cudadrv_type), &
      INTENT(INOUT)                          :: this
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    TYPE(dbcsr_data_obj)                     :: c_area

! Transfer C-data from GPU to host and adding it to host's result

    CALL dbcsr_data_dev2host(this%c_buffer, error)
    !TODO: not async, yet
    CALL dbcsr_cuda_stream_synchronize(this%c_buffer%d%cuda_stream, error)

    c_area = this%product_wm%data_area
    CALL block_add(c_area, this%c_buffer, error=error)
    CALL dbcsr_data_release(this%c_buffer)

  END SUBROUTINE dbcsr_mm_cudadrv_finalize


! *****************************************************************************
! \brief Sort stack entries with respect to the c_id.
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE stack_sort(params_in, params_out, stack_size)
    INTEGER, INTENT(IN)                      :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width, &
      stack_size), INTENT(OUT)               :: params_out
    INTEGER, DIMENSION(dbcsr_ps_width, &
      stack_size), INTENT(IN)                :: params_in

    INTEGER                                  :: i
    INTEGER, DIMENSION(stack_size)           :: c_sort, c_sort_ind

! sort by the C-blocks

     c_sort = params_in(6,:stack_size)
     CALL sort(c_sort, stack_size, c_sort_ind)
     DO i=1,stack_size
        params_out(:, i) = params_in(:, c_sort_ind(i))
     END DO

  END SUBROUTINE stack_sort


! *****************************************************************************
! \brief Roughly order stacks with a cheaper Binning-scheme by Peter Messmer
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE stack_binning(params_in, params_out, stack_size)
    INTEGER, INTENT(IN)                      :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width, &
      stack_size), INTENT(OUT)               :: params_out
    INTEGER, DIMENSION(dbcsr_ps_width, &
      stack_size), INTENT(IN)                :: params_in

    INTEGER                                  :: bin_id, i, top
    INTEGER, &
      DIMENSION(cudadrv_binning_nbins)       :: bin_top
    INTEGER, DIMENSION(dbcsr_ps_width)       :: val
    INTEGER, DIMENSION(dbcsr_ps_width, &
      cudadrv_binning_binsize, &
      cudadrv_binning_nbins)                 :: bin_arr

     bin_top = 1
     top = 1
     DO i=1,stack_size
        val(:) = params_in(:,i)
        bin_id = 1 + MODULO(INT(val(6)*(val(6)+3),KIND=int_8), &
                            INT(cudadrv_binning_nbins,KIND=int_8))
        IF(bin_top(bin_id) > cudadrv_binning_binsize) THEN
           params_out(:, top:top+bin_top(bin_id)-2) = bin_arr(:, 1:bin_top(bin_id)-1, bin_id)
           top = top + bin_top(bin_id) - 1
           bin_top(bin_id) = 1
        ENDIF
        bin_arr(:, bin_top(bin_id), bin_id) =  val(:)
        bin_top(bin_id) = bin_top(bin_id) + 1
     END DO
     DO i=1,cudadrv_binning_nbins
       IF(bin_top(i) > 1) THEN
          params_out(:, top:top+bin_top(i)-2) = bin_arr(:, 1:bin_top(i)-1, i)
          top = top + bin_top(i)-1
       ENDIF
     END DO

  END SUBROUTINE  stack_binning

! *****************************************************************************
! \brief Processes a given stack using CUDA
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_cudadrv_process(this, left, right, params, stack_size, &
       stack_descr, success, error)
    TYPE(dbcsr_mm_cudadrv_type), &
      INTENT(INOUT)                          :: this
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    INTEGER, INTENT(IN)                      :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width, &
      stack_size), INTENT(INOUT)             :: params
    TYPE(stack_descriptor_type), INTENT(IN)  :: stack_descr
    LOGICAL, INTENT(OUT)                     :: success
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_cudadrv_process', &
      routineP = moduleN//':'//routineN

    INTEGER :: error_handle, error_handle2, flop_per_entry, &
      next_stack_buffer_idx, stacked_datasize
    TYPE(dbcsr_data_area_type), POINTER      :: a_area, b_area, c_area
    TYPE(stack_buffer_type), POINTER         :: stackbuf

    next_stack_buffer_idx = MOD(curr_stack_buffer_idx, SIZE(stack_buffers))+1
    stackbuf => stack_buffers(next_stack_buffer_idx)

    IF(.NOT. dbcsr_cuda_event_query(stackbuf%ready, error)) THEN
       !WRITE (*,*) "dbcsr_mm_cudadrv_process: next buffer not ready"
       success = .FALSE.
       RETURN
    END IF

    curr_stack_buffer_idx = next_stack_buffer_idx
    !WRITE (*,*) "dbcsr_mm_cudadrv_process: using buffer: ", curr_stack_buffer_idx

    CALL dbcsr_error_set (routineN, error_handle, error)
    stacked_datasize   = this%product_wm%datasize
    CALL dbcsr_data_ensure_size(this%c_buffer, stacked_datasize, &
                     factor=default_resize_factor, zero_pad=.TRUE., error=error)

    !===========================================================================
    ! sort the stack. Since this costs CPU time, only a good idea if the CPUs
    ! are not too busy, or GPU gain is very large
    IF (careful_mod) CALL dbcsr_error_set (routineN//"_sort", error_handle2, error)
    flop_per_entry = 2*stack_descr%max_m*stack_descr%max_n*stack_descr%max_k

    IF(flop_per_entry > cudadrv_min_flop_sort) THEN
       CALL stack_sort(params, stackbuf%hostmem, stack_size)
    ELSE
       CALL stack_binning(params, stackbuf%hostmem, stack_size)
    ENDIF

    IF (careful_mod) CALL dbcsr_error_stop (error_handle2, error)


    CALL dbcsr_cuda_devmem_host2dev(stackbuf%devmem, hostmem=stackbuf%hostmem, stream=stackbuf%stream, error=error)
    CALL dbcsr_cuda_event_record(stackbuf%ready, stream=stackbuf%stream, error=error)

    a_area => left%data_area%d
    b_area => right%data_area%d
    c_area => this%c_buffer%d

    !WRITE (*,*) "dbcsr_mm_cudadrv_process: a_area%memory_type ", a_area%memory_type
    !WRITE (*,*) "dbcsr_mm_cudadrv_process: b_area%memory_type ", b_area%memory_type
    !WRITE (*,*) "dbcsr_mm_cudadrv_process: c_area%memory_type ", c_area%memory_type

    IF(.NOT. dbcsr_cuda_devmem_allocated(a_area%cuda_devmem)) &
       CALL mp_abort("dbcsr_mm_cudadrv_process: a_area%cuda_devmem not allocated")
    IF(.NOT. dbcsr_cuda_devmem_allocated(b_area%cuda_devmem)) &
       CALL mp_abort("dbcsr_mm_cudadrv_process: b_area%cuda_devmem not allocated")
    IF(.NOT. dbcsr_cuda_devmem_allocated(c_area%cuda_devmem)) &
       CALL mp_abort("dbcsr_mm_cudadrv_process: c_area%cuda_devmem not allocated")

    CALL dbcsr_cuda_stream_wait_event(stackbuf%stream, a_area%cuda_ready, error)
    CALL dbcsr_cuda_stream_wait_event(stackbuf%stream, b_area%cuda_ready, error)
    CALL dbcsr_cuda_stream_wait_event(stackbuf%stream, c_area%cuda_ready, error)

    CALL dbcsr_cuda_do_mm_stack(stackbuf%devmem, stack_size,&
         a_data=a_area%cuda_devmem,&
         b_data=b_area%cuda_devmem,&
         c_data=c_area%cuda_devmem,&
         m_max=stack_descr%max_m,&
         n_max=stack_descr%max_n,&
         k_max=stack_descr%max_k,&
         def_mnk=stack_descr%defined_mnk,&
         stream=stackbuf%stream, &
         success=success, error=error)

    IF(success) THEN
      CALL dbcsr_cuda_event_record(stackbuf%calculated, stream=stackbuf%stream, error=error)
      CALL dbcsr_cuda_stream_wait_event(a_area%cuda_stream, stackbuf%calculated, error)
      CALL dbcsr_cuda_stream_wait_event(b_area%cuda_stream, stackbuf%calculated, error)
      CALL dbcsr_cuda_stream_wait_event(c_area%cuda_stream, stackbuf%calculated, error)
    ENDIF

    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_mm_cudadrv_process



END MODULE dbcsr_mm_cudadrv


