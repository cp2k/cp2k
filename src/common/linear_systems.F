!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Provides interfaces to LAPACK routines for factorisation and
!>      linear system solving
!> \note
!>      We are using LAPACK interfaces, so please make sure in IBM/AIX you have
!>      the lapack library before essl: "xlf90 ... -llapack -lessl" !!!
!> \par History
!>      none
!> \author JGH (30-5-2001)
! *****************************************************************************
MODULE linear_systems

  
  USE kinds,                           ONLY: dp
  USE lapack,                          ONLY: lapack_sgesv
  USE termination,                     ONLY: stop_program
#include "cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'linear_systems'

  PUBLIC :: solve_system

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param mysize ...
!> \param eigenvectors ...
! *****************************************************************************
SUBROUTINE solve_system ( matrix, mysize, eigenvectors )

    REAL(KIND=dp), INTENT(INOUT)             :: matrix( :, : )
    INTEGER, INTENT(IN)                      :: mysize
    REAL(KIND=dp), INTENT(INOUT)             :: eigenvectors( :, : )

    CHARACTER(len=*), PARAMETER :: routineN = 'solve_system', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: info, lda, ldb, nrhs, &
                                                ipiv( mysize )

  lda = SIZE ( matrix, 1 )
  ldb = SIZE ( eigenvectors, 1 )
  nrhs = SIZE ( eigenvectors, 2 )

  CALL lapack_sgesv ( mysize, nrhs, matrix, lda, ipiv, &
                      eigenvectors, ldb, info )
  IF ( info /= 0 ) THEN
     CALL stop_program(routineN,moduleN,__LINE__,"Error in inversion")
  END IF

END SUBROUTINE solve_system

END MODULE linear_systems

