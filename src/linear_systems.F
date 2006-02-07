!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/linear_systems [1.0] *
!!
!!   NAME
!!     linear_systems
!!
!!   FUNCTION
!!     Provides interfaces to LAPACK routines for factorisation and
!!     linear system solving
!!
!!   AUTHOR
!!     JGH (30-5-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!     We are using LAPACK interfaces, so please make sure in IBM/AIX you have
!!     the lapack library before essl: "xlf90 ... -llapack -lessl" !!!
!!
!!   SOURCE
!******************************************************************************

MODULE linear_systems

  USE kinds,                           ONLY: dp
  USE lapack,                          ONLY: lapack_cpotrf,&
                                             lapack_ctrtri,&
                                             lapack_sgesv,&
                                             lapack_spotrf,&
                                             lapack_strtri
  USE termination,                     ONLY: stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: cholesky_inverse 
  PUBLIC :: solve_system 

  INTERFACE cholesky_inverse
     MODULE PROCEDURE ch_inv_pdsym
     MODULE PROCEDURE ch_inv_pdher
  END INTERFACE

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE ch_inv_pdsym ( matrix, mysize, storageform )


    REAL(KIND=dp), INTENT(INOUT)             :: matrix( :, : )
    INTEGER, INTENT(IN)                      :: mysize
    CHARACTER(LEN=*), INTENT(IN)             :: storageform

    CHARACTER                                :: diag = "N", uplo
    INTEGER                                  :: info, lda

!------------------------------------------------------------------------------

  IF ( storageform ( 1:5 ) == "Lower" .OR. &
       storageform ( 1:5 ) == "LOWER" .OR. &
       storageform ( 1:5 ) == "lower" ) THEN
     uplo = "L"
  ELSE IF ( storageform ( 1:5 ) == "Upper" .OR. &
       storageform ( 1:5 ) == "upper" .OR. &
       storageform ( 1:5 ) == "UPPER" ) THEN
     uplo = "U"
  ELSE
     CALL stop_program ( "ch_inv_pdsym", "unknown form of storage" )
  END IF

  lda = SIZE ( matrix, 1 )

  CALL lapack_spotrf ( uplo, mysize, matrix, lda, info )
  IF ( info /= 0 ) THEN
     CALL stop_program ( "ch_inv_pdsym", "error in Cholesky factorisation" )
  END IF
  CALL lapack_strtri ( uplo, diag, mysize, matrix, lda, info )
  IF ( info /= 0 ) THEN
     CALL stop_program ( "ch_inv_pdsym", "error in inversion" )
  END IF

END SUBROUTINE ch_inv_pdsym

!******************************************************************************

SUBROUTINE ch_inv_pdher ( matrix, mysize, storageform )


    COMPLEX(KIND=dp), INTENT(INOUT)          :: matrix( :, : )
    INTEGER, INTENT(IN)                      :: mysize
    CHARACTER(LEN=*), INTENT(IN)             :: storageform

    CHARACTER                                :: diag = "N", uplo
    INTEGER                                  :: info, lda

!------------------------------------------------------------------------------

  IF ( storageform ( 1:5 ) == "Lower" .OR. &
       storageform ( 1:5 ) == "LOWER" .OR. &
       storageform ( 1:5 ) == "lower" ) THEN
     uplo = "L"
  ELSE IF ( storageform ( 1:5 ) == "Upper" .OR. &
       storageform ( 1:5 ) == "upper" .OR. &
       storageform ( 1:5 ) == "UPPER" ) THEN
     uplo = "U"
  ELSE
     CALL stop_program ( "ch_inv_pdher", "unknown form of storage" )
  END IF

  lda = SIZE ( matrix, 1 )

  CALL lapack_cpotrf ( uplo, mysize, matrix, lda, info )
  IF ( info /= 0 ) THEN
     CALL stop_program ( "ch_inv_pdher", "error in Cholesky factorisation" )
  END IF
  CALL lapack_ctrtri ( uplo, diag, mysize, matrix, lda, info )
  IF ( info /= 0 ) THEN
     CALL stop_program ( "ch_inv_pdher", "error in inversion" )
  END IF

END SUBROUTINE ch_inv_pdher

!******************************************************************************
SUBROUTINE solve_system ( matrix, mysize, eigenvectors )


    REAL(KIND=dp), INTENT(INOUT)             :: matrix( :, : )
    INTEGER, INTENT(IN)                      :: mysize
    REAL(KIND=dp), INTENT(INOUT)             :: eigenvectors( :, : )

    INTEGER                                  :: info, lda, ldb, nrhs
    INTEGER                                  :: ipiv( mysize )

  lda = SIZE ( matrix, 1 )
  ldb = SIZE ( eigenvectors, 1 )
  nrhs = SIZE ( eigenvectors, 2 )

  CALL lapack_sgesv ( mysize, nrhs, matrix, lda, ipiv, &
                      eigenvectors, ldb, info )
  IF ( info /= 0 ) THEN
     CALL stop_program ( "solve_system", "error in inversion" )
  END IF

END SUBROUTINE solve_system

!******************************************************************************

END MODULE linear_systems

!******************************************************************************
