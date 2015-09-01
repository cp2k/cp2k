!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
                                             dbcsr_error_type,&
                                             dbcsr_fatal_level,&
                                             dbcsr_wrong_args_error
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8
  USE kinds,                           ONLY: int_8

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_blas_operations'

  PUBLIC :: dbcsr_lapack_larnv, set_larnv_seed

CONTAINS

! *****************************************************************************
!> \brief generate a seed respecting the lapack constraints, 
!>        - values between 0..4095 (2**12-1)
!>        - iseed(4) odd
!>        also try to avoid iseed that are zero.
!>        Have but with a unique 2D mapping (irow,icol), and a 'mini-seed' ival
!> \param irow 1..nrow
!> \param nrow nrow
!> \param icol 1..ncol
!> \param ncol ncol
!> \param ival mini-seed
!> \param iseed a lapack compatible seed
! *****************************************************************************
  SUBROUTINE set_larnv_seed(irow,nrow,icol,ncol,ival,iseed)
    INTEGER, INTENT(IN)                      :: irow, nrow, icol, ncol, ival
    INTEGER, INTENT(OUT)                     :: iseed(4)

    CHARACTER(len=*), PARAMETER :: routineN = 'set_larnv_seed', &
      routineP = moduleN//':'//routineN

    INTEGER(KIND=int_8)                      :: map

    map=((irow-1+icol*INT(nrow,int_8))*(1+MODULO(ival,2**16)))*2+1+0*ncol ! ncol used
    iseed(4)=INT(MODULO(     map      ,2_int_8**12)) ; map=map/2_int_8**12 ; ! keep odd
    iseed(3)=INT(MODULO(IEOR(map,3541_int_8),2_int_8**12)) ; map=map/2_int_8**12 
    iseed(2)=INT(MODULO(IEOR(map,1153_int_8),2_int_8**12)) ; map=map/2_int_8**12 
    iseed(1)=INT(MODULO(IEOR(map,2029_int_8),2_int_8**12)) ; map=map/2_int_8**12 
  END SUBROUTINE set_larnv_seed
! *****************************************************************************
!> \brief fils an array with random numbers
!> \param idist type of distribution (1..3, see lapack docs)
!> \param iseed requires each int to be in the range 0..2**12, and the iseed(4) odd
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

  END SUBROUTINE dbcsr_lapack_larnv

END MODULE dbcsr_blas_operations
