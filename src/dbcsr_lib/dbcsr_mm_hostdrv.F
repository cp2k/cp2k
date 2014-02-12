!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
MODULE dbcsr_mm_hostdrv
  USE dbcsr_config,                    ONLY: mm_driver,&
                                             mm_driver_blas,&
                                             mm_driver_cuda,&
                                             mm_driver_matmul,&
                                             mm_driver_plasma,&
                                             mm_driver_smm,&
                                             mm_host_driver
  USE dbcsr_data_methods,              ONLY: dbcsr_data_get_size
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: dp,&
                                             real_4,&
                                             real_8,&
                                             sp
  USE dbcsr_mm_types,                  ONLY: dbcsr_ps_width,&
                                             p_a_first,&
                                             p_b_first,&
                                             p_c_first,&
                                             p_k,&
                                             p_m,&
                                             p_n,&
                                             stack_descriptor_type
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_type,&
                                             dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8,&
                                             dbcsr_work_type

  !$ USE OMP_LIB


  IMPLICIT NONE

  PRIVATE


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_hostdrv'

  CHARACTER(len=*), PARAMETER, PRIVATE :: int_print = "(10(1X,I7))"

  INTEGER :: multrec_calls = 0

  LOGICAL, PARAMETER, PUBLIC :: show_levels = .FALSE.
  LOGICAL, PARAMETER :: verbose_acc = .FALSE.

  PUBLIC :: dbcsr_mm_hostdrv_process
  PUBLIC :: dbcsr_mm_hostdrv_type
  PUBLIC :: dbcsr_mm_hostdrv_init, dbcsr_mm_hostdrv_finalize
  PUBLIC :: dbcsr_mm_hostdrv_phaseout, dbcsr_mm_hostdrv_barrier

  LOGICAL, PARAMETER :: debug_mod  = .FALSE.
  LOGICAL, PARAMETER :: careful_mod = .FALSE.

  TYPE dbcsr_mm_hostdrv_type
      TYPE(dbcsr_data_obj)          :: data_area
  END TYPE dbcsr_mm_hostdrv_type

CONTAINS

! *****************************************************************************
! \brief Initialize the library
! *****************************************************************************
 SUBROUTINE dbcsr_mm_hostdrv_init(this, left, right, product_wm, error)
    TYPE(dbcsr_mm_hostdrv_type), &
      INTENT(INOUT)                          :: this
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_work_type), POINTER           :: product_wm
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_hostdrv_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)

    this%data_area = product_wm%data_area

    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_mm_hostdrv_init


! *****************************************************************************
! \brief Signal approaching end of multiplication
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_hostdrv_phaseout(this, error)
    TYPE(dbcsr_mm_hostdrv_type), &
      INTENT(INOUT)                          :: this
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

! nothing to do

  END SUBROUTINE dbcsr_mm_hostdrv_phaseout


! *****************************************************************************
! \brief Signal that previous stacks should be processed first
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_hostdrv_barrier(this, error)
    TYPE(dbcsr_mm_hostdrv_type), &
      INTENT(INOUT)                          :: this
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

! nothing to do here

  END SUBROUTINE dbcsr_mm_hostdrv_barrier


! *****************************************************************************
! \brief Finalize the library
! *****************************************************************************
  SUBROUTINE dbcsr_mm_hostdrv_finalize(this, error)
    TYPE(dbcsr_mm_hostdrv_type), &
      INTENT(INOUT)                          :: this
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

! nothing to do here

  END SUBROUTINE dbcsr_mm_hostdrv_finalize

! *****************************************************************************
!> \brief Calls the various drivers that process the stack.
!>
!> \param[in] params           Stack of GEMM parameters
!> \param[in] n                Number of parameters
!> \param[in] ab%left   Left-matrix data
!> \param[in] ab%right  Right-matrix data
!> \param[in,out] zero_first   Zero product data area starting from this
!>                             element
!> \param[in,out] zero_last    Zero product data area up to this element
!> \param[in] lastblk          Number of blocks in product
!> \param[in,out] this%data_area  Data for results
! *****************************************************************************
  SUBROUTINE dbcsr_mm_hostdrv_process(this, left, right, params, stack_size, &
       stack_descr, success, error)
    TYPE(dbcsr_mm_hostdrv_type), &
      INTENT(INOUT)                          :: this
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    INTEGER, INTENT(IN)                      :: stack_size
    INTEGER, DIMENSION(1:dbcsr_ps_width, &
      stack_size), INTENT(INOUT)             :: params
    TYPE(stack_descriptor_type), INTENT(IN)  :: stack_descr
    LOGICAL, INTENT(OUT)                     :: success
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_hostdrv_process', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: careful = careful_mod, &
                                                dbg = .FALSE.

    INTEGER                                  :: error_handle, sp
    REAL(KIND=dp)                            :: rnd

    IF(mm_driver==mm_driver_cuda) & !for cpu-only runs this is called too often
       CALL dbcsr_error_set (routineN, error_handle, error)

    success = .TRUE. !host driver never failes...hopefully

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
               "LE", dbcsr_data_get_size (left%data_area),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "A data out of bounds.", __LINE__, error=error)
          CALL dbcsr_assert (params(p_b_first,sp)&
               + params(p_k,sp) * params(p_n,sp) - 1,&
               "LE", dbcsr_data_get_size (right%data_area),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "B data out of bounds.", __LINE__, error=error)
          CALL dbcsr_assert (params(p_c_first,sp)&
               + params(p_m,sp) * params(p_n,sp) - 1,&
               "LE", dbcsr_data_get_size (this%data_area),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "C data out of bounds.", __LINE__, error=error)
       ENDDO
    ENDIF

    SELECT CASE (mm_host_driver)
    CASE (mm_driver_matmul)
       SELECT CASE (this%data_area%d%data_type)
       CASE (dbcsr_type_real_4)
          CALL internal_process_mm_stack_s (params, &
               stack_size, &
               left%data_area%d%r_sp, right%data_area%d%r_sp, this%data_area%d%r_sp,&
               error=error)
       CASE (dbcsr_type_real_8)
          CALL internal_process_mm_stack_d (params,&
               stack_size,&
               left%data_area%d%r_dp, right%data_area%d%r_dp, this%data_area%d%r_dp,&
               error=error)
       CASE (dbcsr_type_complex_4)
          CALL internal_process_mm_stack_c (params,&
               stack_size,&
               left%data_area%d%c_sp, right%data_area%d%c_sp, this%data_area%d%c_sp,&
               error=error)
       CASE (dbcsr_type_complex_8)
          CALL internal_process_mm_stack_z (params,&
               stack_size,&
               left%data_area%d%c_dp, right%data_area%d%c_dp, this%data_area%d%c_dp,&
               error=error)
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type",__LINE__,error)
       END SELECT
    CASE (mm_driver_smm)
       SELECT CASE (this%data_area%d%data_type)
       CASE (dbcsr_type_real_4)
          CALL smm_process_mm_stack_s (stack_descr, params, &
               stack_size, &
               left%data_area%d%r_sp, right%data_area%d%r_sp, this%data_area%d%r_sp,&
               error=error)
       CASE (dbcsr_type_real_8)
          CALL smm_process_mm_stack_d (stack_descr, params,&
               stack_size,&
               left%data_area%d%r_dp, right%data_area%d%r_dp, this%data_area%d%r_dp,&
               error=error)
       CASE (dbcsr_type_complex_4)
          CALL smm_process_mm_stack_c (stack_descr, params,&
               stack_size,&
               left%data_area%d%c_sp, right%data_area%d%c_sp, this%data_area%d%c_sp,&
               error=error)
       CASE (dbcsr_type_complex_8)
          CALL smm_process_mm_stack_z (stack_descr, params,&
               stack_size,&
               left%data_area%d%c_dp, right%data_area%d%c_dp, this%data_area%d%c_dp,&
               error=error)
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type",__LINE__,error)
       END SELECT
    CASE (mm_driver_plasma)
       SELECT CASE (this%data_area%d%data_type)
       CASE (dbcsr_type_real_4)
          CALL plasma_process_mm_stack_s (params,&
               stack_size,&
               left%data_area%d%r_sp, right%data_area%d%r_sp, this%data_area%d%r_sp,&
               error=error)
       CASE (dbcsr_type_real_8)
          CALL plasma_process_mm_stack_d (params,&
               stack_size,&
               left%data_area%d%r_dp, right%data_area%d%r_dp, this%data_area%d%r_dp,&
               error=error)
       CASE (dbcsr_type_complex_4)
          CALL plasma_process_mm_stack_c (params,&
               stack_size,&
               left%data_area%d%c_sp, right%data_area%d%c_sp, this%data_area%d%c_sp,&
               error=error)
       CASE (dbcsr_type_complex_8)
          CALL plasma_process_mm_stack_z (params,&
               stack_size,&
               left%data_area%d%c_dp, right%data_area%d%c_dp, this%data_area%d%c_dp,&
               error=error)
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type",__LINE__,error)
       END SELECT
    CASE (mm_driver_blas)
       SELECT CASE (this%data_area%d%data_type)
       CASE (dbcsr_type_real_4)
          CALL blas_process_mm_stack_s (params,&
               stack_size,&
               left%data_area%d%r_sp, right%data_area%d%r_sp, this%data_area%d%r_sp,&
               error=error)
       CASE (dbcsr_type_real_8)
          CALL blas_process_mm_stack_d (params,&
               stack_size,&
               left%data_area%d%r_dp, right%data_area%d%r_dp, this%data_area%d%r_dp,&
               error=error)
       CASE (dbcsr_type_complex_4)
          CALL blas_process_mm_stack_c (params,&
               stack_size,&
               left%data_area%d%c_sp, right%data_area%d%c_sp, this%data_area%d%c_sp,&
               error=error)
       CASE (dbcsr_type_complex_8)
          CALL blas_process_mm_stack_z (params,&
               stack_size,&
               left%data_area%d%c_dp, right%data_area%d%c_dp, this%data_area%d%c_dp,&
               error=error)
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type",__LINE__,error)
       END SELECT
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Invalid multiplication driver",__LINE__,error)
    END SELECT


    IF(mm_driver==mm_driver_cuda) & !for cpu-only runs this is called too often
       CALL dbcsr_error_stop(error_handle, error)


  END SUBROUTINE dbcsr_mm_hostdrv_process


! *****************************************************************************
!> \brief Helper-routine used by dbcsr_mm_hostdrv_process to print debug info.
! *****************************************************************************
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


#include "dbcsr_mm_hostdrv_d.f90"
#include "dbcsr_mm_hostdrv_z.f90"
#include "dbcsr_mm_hostdrv_s.f90"
#include "dbcsr_mm_hostdrv_c.f90"


END MODULE dbcsr_mm_hostdrv
