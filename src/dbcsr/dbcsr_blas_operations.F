!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Wrappers to BLAS calls.
!> \author  Urban Borstnik
!> \date    2010-02-18
!> \version 0.9
!>
!> <b>Modification history:</b>
!> - 2010-02-18 Moved from dbcsr_util
! *****************************************************************************
MODULE dbcsr_blas_operations

  USE dbcsr_data_methods,              ONLY: dbcsr_data_clear_2d_pointer,&
                                             dbcsr_data_set_2d_pointer
  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
                                             dbcsr_error_set,&
                                             dbcsr_error_type,&
                                             dbcsr_fatal_level,&
                                             dbcsr_wrong_args_error
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_scalar_type,&
                                             dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8

  !$ USE OMP_LIB
  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_blas_operations'

  ! BLAS routines
  PUBLIC :: dbcsr_blas_copy, dbcsr_blas_gemm

  ! LAPACK routine
  PUBLIC :: dbcsr_lapack_larnv


CONTAINS

!\> Interface to LAPACK larnv
  SUBROUTINE dbcsr_lapack_larnv( idist, iseed, n, x, error )
    INTEGER, INTENT(IN)                      :: idist
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: iseed
    INTEGER, INTENT(IN)                      :: n
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: x
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_lapack_larnv', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------
!CALL dbcsr_error_set (routineN, error_handler, error)

    SELECT CASE (x%d%data_type)
    CASE (dbcsr_type_real_4)
       CALL slarnv( idist, iseed, n, x%d%r_sp)
    CASE (dbcsr_type_real_8)
       CALL dlarnv( idist, iseed, n, x%d%r_dp)
    CASE (dbcsr_type_complex_4)
       CALL clarnv( idist, iseed, n, x%d%c_sp)
    CASE (dbcsr_type_complex_8)
       CALL zlarnv( idist, iseed, n, x%d%c_dp)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            routineN, "Invalid data type",__LINE__,error)
    END SELECT

    !CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_lapack_larnv

!\> Interface to BLAS SCAL
  SUBROUTINE dbcsr_blas_scal (N, alpha, x, INCX)
    INTEGER, INTENT(IN)                      :: N
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: alpha
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: x
    INTEGER, INTENT(IN)                      :: INCX

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_blas_scal', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

    CALL dbcsr_assert (alpha%data_type, "EQ", x%d%data_type,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Incompatible data types",__LINE__,error)
    SELECT CASE (alpha%data_type)
    CASE (dbcsr_type_real_4)
       !CALL BL_SCAL (n, alpha%r_sp, x%d%r_sp, incx)
       CALL SSCAL (n, alpha%r_sp, x%d%r_sp, incx)
    CASE (dbcsr_type_real_8)
       !CALL BL_SCAL (n, alpha%r_dp, x%d%r_dp, incx)
       CALL DSCAL (n, alpha%r_dp, x%d%r_dp, incx)
    CASE (dbcsr_type_complex_4)
       !CALL BL_SCAL (n, alpha%c_sp, x%d%c_sp, incx)
       CALL CSCAL (n, alpha%c_sp, x%d%c_sp, incx)
    CASE (dbcsr_type_complex_8)
       !CALL BL_SCAL (n, alpha%c_dp, x%d%c_dp, incx)
       CALL ZSCAL (n, alpha%c_dp, x%d%c_dp, incx)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            routineN, "Invalid data type",__LINE__,error)
    END SELECT
  END SUBROUTINE dbcsr_blas_scal

!\> Interface to BLAS DGEMM
  SUBROUTINE dbcsr_blas_gemm(TRANSA,TRANSB,&
       M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    CHARACTER(LEN=1), INTENT(IN)             :: TRANSA, TRANSB
    INTEGER, INTENT(IN)                      :: M, N, K
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: alpha
    TYPE(dbcsr_data_obj), INTENT(IN)         :: a
    INTEGER, INTENT(IN)                      :: LDA
    TYPE(dbcsr_data_obj), INTENT(IN)         :: b
    INTEGER, INTENT(IN)                      :: LDB
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: beta
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: c
    INTEGER, INTENT(IN)                      :: LDC

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_blas_gemm', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: assertions = .FALSE.

    TYPE(dbcsr_data_obj)                     :: a_tmp, b_tmp
    TYPE(dbcsr_error_type)                   :: error

    IF (assertions) THEN
       CALL dbcsr_assert (alpha%data_type, "EQ", beta%data_type,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Incompatible data types, alpha/beta",__LINE__,error)
       CALL dbcsr_assert (alpha%data_type, "EQ", a%d%data_type,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Incompatible data types, alpha, a",__LINE__,error)
       CALL dbcsr_assert (alpha%data_type, "EQ", b%d%data_type,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Incompatible data types, alpha, b",__LINE__,error)
       CALL dbcsr_assert (alpha%data_type, "EQ", c%d%data_type,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Incompatible data types, alpha, c",__LINE__,error)
    ENDIF
    a_tmp = a
    CALL dbcsr_data_set_2d_pointer (a_tmp, lda)
    b_tmp = b
    CALL dbcsr_data_set_2d_pointer (b_tmp, ldb)
    CALL dbcsr_data_set_2d_pointer (c, ldc)
    SELECT CASE (alpha%data_type)
    CASE (dbcsr_type_real_4)
       !CALL BL_GEMM (TRANSA, TRANSB, M, N, K, alpha%r_sp,&
       !     a%d%r2_sp, LDA, b%d%r2_sp, LDB, beta%r_sp, c%d%r2_sp, LDC)
       CALL SGEMM (TRANSA, TRANSB, M, N, K, alpha%r_sp,&
            a%d%r2_sp, LDA, b%d%r2_sp, LDB, beta%r_sp, c%d%r2_sp, LDC)
    CASE (dbcsr_type_real_8)
       !CALL BL_GEMM (TRANSA, TRANSB, M, N, K, alpha%r_dp,&
       !     a%d%r2_dp, LDA, b%d%r2_dp, LDB, beta%r_dp, c%d%r2_dp, LDC)
       CALL DGEMM (TRANSA, TRANSB, M, N, K, alpha%r_dp,&
            a%d%r2_dp, LDA, b%d%r2_dp, LDB, beta%r_dp, c%d%r2_dp, LDC)
    CASE (dbcsr_type_complex_4)
       !CALL BL_GEMM (TRANSA, TRANSB, M, N, K, alpha%c_sp,&
       !     a%d%c2_sp, LDA, b%d%c2_sp, LDB, beta%c_sp, c%d%c2_sp, LDC)
       CALL CGEMM (TRANSA, TRANSB, M, N, K, alpha%c_sp,&
            a%d%c2_sp, LDA, b%d%c2_sp, LDB, beta%c_sp, c%d%c2_sp, LDC)
    CASE (dbcsr_type_complex_8)
       !CALL BL_GEMM (TRANSA, TRANSB, M, N, K, alpha%c_dp,&
       !     a%d%c2_dp, LDA, b%d%c2_dp, LDB, beta%c_dp, c%d%c2_dp, LDC)
       CALL ZGEMM (TRANSA, TRANSB, M, N, K, alpha%c_dp,&
            a%d%c2_dp, LDA, b%d%c2_dp, LDB, beta%c_dp, c%d%c2_dp, LDC)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            routineN, "Invalid data type",__LINE__,error)
    END SELECT
    CALL dbcsr_data_clear_2d_pointer (a_tmp)
    CALL dbcsr_data_clear_2d_pointer (b_tmp)
    CALL dbcsr_data_clear_2d_pointer (c)
  END SUBROUTINE dbcsr_blas_gemm

  SUBROUTINE dbcsr_blas_copy(N,X,INCX,Y,INCY)
    INTEGER, INTENT(IN)                      :: N
    TYPE(dbcsr_data_obj), INTENT(IN)         :: x
    INTEGER, INTENT(IN)                      :: INCX
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: y
    INTEGER, INTENT(IN)                      :: INCY

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_blas_copy', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

    CALL dbcsr_assert (x%d%data_type, "EQ", y%d%data_type,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Incompatible data types",__LINE__,error)
    SELECT CASE (x%d%data_type)
    CASE (dbcsr_type_real_4)
       CALL SCOPY (n, x%d%r_sp, incx, y%d%r_sp, incy)
    CASE (dbcsr_type_real_8)
       CALL DCOPY (n, x%d%r_dp, incx, y%d%r_dp, incy)
    CASE (dbcsr_type_complex_4)
       CALL CCOPY (n, x%d%c_sp, incx, y%d%c_sp, incy)
    CASE (dbcsr_type_complex_8)
       CALL ZCOPY (n, x%d%c_dp, incx, y%d%c_dp, incy)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            routineN, "Invalid data type",__LINE__,error)
    END SELECT
  END SUBROUTINE dbcsr_blas_copy

END MODULE dbcsr_blas_operations
