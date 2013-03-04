!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
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
  USE dbcsr_config,                    ONLY: driver_is_async,&
                                             mm_driver_blas,&
                                             mm_driver_matmul,&
                                             mm_driver_plasma,&
                                             mm_driver_smm
  USE dbcsr_cuda_operations,           ONLY: dbcsr_cuda_cp_dev_to_host,&
                                             dbcsr_cuda_cp_host_to_dev,&
                                             dbcsr_cuda_do_mm_stack
  USE dbcsr_cuda_types,                ONLY: dbcsr_cuda_mem_type
  USE dbcsr_data_methods,              ONLY: dbcsr_data_get_size
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: dp,&
                                             real_4,&
                                             real_8,&
                                             sp
  USE dbcsr_mm_types,                  ONLY: &
       dbcsr_ps_state_empty, dbcsr_ps_state_working, dbcsr_ps_width, &
       p_a_first, p_b_first, p_c_first, p_k, p_m, p_n
  USE dbcsr_toollib,                   ONLY: sort
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8

  !$ USE OMP_LIB


  IMPLICIT NONE

  PRIVATE


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_stack'

  CHARACTER(len=*), PARAMETER, PRIVATE :: int_print = "(10(1X,I7))"

  REAL, PARAMETER                      :: default_resize_factor = 1.618034

  INTEGER :: multrec_calls = 0

  LOGICAL, PARAMETER, PUBLIC :: show_levels = .FALSE.
  LOGICAL, PARAMETER :: verbose_acc = .FALSE.

  PUBLIC :: process_mm_stack_cpu, process_mm_stack_cuda

  LOGICAL, PARAMETER :: debug_mod  = .FALSE.
  LOGICAL, PARAMETER :: careful_mod = .FALSE.

CONTAINS


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
  SUBROUTINE process_mm_stack_cpu(params, driver,&
       stack_size, &
       left_data_area, right_data_area, product_data_area,&
       state, m, n, k, max_m, max_n, max_k, defined_mnk,&
       stack_id,&
       error)
    INTEGER, INTENT(IN)                      :: stack_size, driver
    INTEGER, DIMENSION(1:dbcsr_ps_width, &
      stack_size), INTENT(INOUT)             :: params
    TYPE(dbcsr_data_obj), INTENT(IN)         :: left_data_area, &
                                                right_data_area
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: product_data_area
    INTEGER, POINTER                         :: state
    INTEGER, INTENT(IN)                      :: m, n, k, max_m, max_n, max_k
    LOGICAL, INTENT(IN)                      :: defined_mnk
    INTEGER, INTENT(IN)                      :: stack_id
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'process_mm_stack_cpu', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: careful = careful_mod, &
                                                dbg = .FALSE.
    REAL, PARAMETER                          :: resize_factor = 1.618034

    INTEGER                                  :: sp
    REAL(KIND=dp)                            :: rnd

    state = dbcsr_ps_state_working
    IF (dbg) THEN
       CALL RANDOM_NUMBER (rnd)
       IF (rnd < 0.01_dp) THEN
          WRITE(*,*)routineN//" Stack size", stack_size, dbcsr_ps_width
          CALL print_gemm_parameters(params(:,1:stack_size))
       ENDIF
    ENDIF

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
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Invalid multiplication driver",__LINE__,error)
    END SELECT

    IF (.NOT. driver_is_async (driver)) THEN
       state = dbcsr_ps_state_empty
    ENDIF
  END SUBROUTINE process_mm_stack_cpu


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
  SUBROUTINE process_mm_stack_cuda(params,&
       stack_size,&
       data_a_dev, data_b_dev, data_c_dev,&
       c_locks,&
       params_dev,&
       m, n, k, max_m, max_n, max_k, defined_mnk,&
       state, stack_state_dev,&
       stack_id,&
       error)
    INTEGER, INTENT(IN)                      :: stack_size
    INTEGER, &
      DIMENSION(dbcsr_ps_width*stack_size), &
      INTENT(INOUT), TARGET                  :: params
    TYPE(dbcsr_cuda_mem_type), INTENT(IN)    :: data_a_dev, data_b_dev
    TYPE(dbcsr_cuda_mem_type), INTENT(INOUT) :: data_c_dev, c_locks, &
                                                params_dev
    INTEGER, INTENT(IN)                      :: m, n, k, max_m, max_n, max_k
    LOGICAL, INTENT(IN)                      :: defined_mnk
    INTEGER, POINTER                         :: state
    TYPE(dbcsr_cuda_mem_type), INTENT(IN)    :: stack_state_dev
    INTEGER, INTENT(IN)                      :: stack_id
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'process_mm_stack_cuda', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: flops_to_sort = 2*13*13*13

    INTEGER                                  :: error_handle, error_handle2, &
                                                i, ithread, k_max, m_max, &
                                                n_max, nthreads, stream
    INTEGER, DIMENSION(:), POINTER           :: params_p
    INTEGER, DIMENSION(stack_size)           :: c_sort, c_sort_ind
    INTEGER, &
      DIMENSION(dbcsr_ps_width, stack_size)  :: params_sort, params_tmp

    nthreads = 1; ithread = 0
    !$ nthreads = OMP_GET_NUM_THREADS () ; ithread = OMP_GET_THREAD_NUM ()

    IF (stack_size .EQ. 0) THEN
       state = dbcsr_ps_state_empty
       WRITE(*,*)routineN//" Stack size is 0"
       RETURN
    ENDIF

    ! ithread is 0-based, stack_id is 1-based, stream should be 1-based
    stream = ithread + nthreads*(stack_id-1) + 1
    CALL dbcsr_error_set (routineN, error_handle, error)

    params_p => params
    IF (m .GT. 0) THEN ; m_max = m ; ELSE ; m_max = -max_m ; ENDIF
    IF (n .GT. 0) THEN ; n_max = n ; ELSE ; n_max = -max_n ; ENDIF
    IF (k .GT. 0) THEN ; k_max = k ; ELSE ; k_max = -max_k ; ENDIF

    !===========================================================================
    ! sort the stack. Since this costs CPU time, only a good idea if the CPUs
    ! are not too busy, or GPU gain is very large
    CALL dbcsr_error_set (routineN//"_sort", error_handle2, error)
    !TODO: now we are allways sorting
    !IF (2*m*n*k>=flops_to_sort .AND. stack_size .GT. 1) THEN
      params_tmp =  RESHAPE(params, (/dbcsr_ps_width, stack_size/))
      ! sort by the C-blocks
      c_sort = params_tmp(6,:)
      CALL sort(c_sort, stack_size, c_sort_ind)
      DO i=1,stack_size
         params_sort(:, i)=params_tmp(:, c_sort_ind(i))
      END DO
      params = RESHAPE(params_sort, (/dbcsr_ps_width*stack_size/))
    !ENDIF
    CALL dbcsr_error_stop (error_handle2, error)
    !===========================================================================


    CALL dbcsr_cuda_cp_host_to_dev (params_p,&
         params_dev%d_i, dbcsr_ps_width*stack_size,&
         async=.TRUE., stream=stream,&
         error=error)

    ! Schedule the resetting of the state of this stack as soon as the
    ! parameter upload finishes.  N.B. This must be the last change to
    ! the stack: its stack_p must be set to 0 and its data must be
    ! used <em>before</em> the state is set to dbcsr_ps_state_empty.
    CALL dbcsr_cuda_cp_dev_to_host (stack_state_dev%d_i, state,&
         async=.TRUE., stream=stream, error=error)

    !TODO: dbcsr_ps_width is a constant, shouldn't need to pass it as argument 
    CALL dbcsr_cuda_do_mm_stack (params_dev%d_i, stack_size, dbcsr_ps_width,&
         data_a_dev, data_b_dev, data_c_dev,&
         c_locks%d_i, ABS(m_max), ABS(n_max), ABS(k_max), defined_mnk,&
         stream=stream, error=error)

    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE process_mm_stack_cuda



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
