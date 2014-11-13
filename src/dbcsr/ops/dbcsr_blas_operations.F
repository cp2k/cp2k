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

  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
                                             dbcsr_error_set,&
                                             dbcsr_error_type,&
                                             dbcsr_fatal_level,&
                                             dbcsr_wrong_args_error
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8

  !$ USE OMP_LIB
  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_blas_operations'

  ! BLAS routines

  ! LAPACK routine
  PUBLIC :: dbcsr_lapack_larnv


CONTAINS

!\> Interface to LAPACK larnv
! *****************************************************************************
!> \brief ...
!> \param idist ...
!> \param iseed ...
!> \param n ...
!> \param x ...
!> \param error ...
! *****************************************************************************
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

END MODULE dbcsr_blas_operations
