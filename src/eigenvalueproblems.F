!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/eigenvalueproblems [1.0] *
!!
!!   NAME
!!     eigenvalueproblems
!!
!!   FUNCTION
!!     Provides interfaces to LAPACK eigenvalue/SVD routines
!!
!!   AUTHOR
!!     APSI
!!
!!   MODIFICATION HISTORY
!!     JGH (26-5-2001): delay D/S C/Z problem to the lapack library call
!!
!!   NOTES
!!     We are using LAPACK interfaces, so please make sure in IBM/AIX you have
!!     the lapack library before essl: "xlf90 ... -llapack -lessl" !!!
!!
!!   SOURCE
!******************************************************************************

MODULE eigenvalueproblems

  USE kinds,                           ONLY: dp
  USE lapack,                          ONLY: lapack_cgesvd,&
                                             lapack_chpev,&
                                             lapack_chpgv,&
                                             lapack_sgesvd,&
                                             lapack_sspgv,&
                                             lapack_ssyev
  USE termination,                     ONLY: stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: diagonalise, diagonalise_general_packed, singular_values

  INTERFACE diagonalise
     MODULE PROCEDURE diagonalise_ssyev
     MODULE PROCEDURE diagonalise_chpev
  END INTERFACE

  INTERFACE diagonalise_general_packed
     MODULE PROCEDURE diagonalise_sspgv
     MODULE PROCEDURE diagonalise_chpgv
  END INTERFACE

  INTERFACE singular_values
     MODULE PROCEDURE cp2k_sgesvd
     MODULE PROCEDURE cp2k_cgesvd
  END INTERFACE

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE diagonalise_ssyev ( matrix, mysize, storageform, eigenvalues, &
     eigenvectors )


    REAL(KIND=dp), INTENT(IN)                :: matrix( :, : )
    INTEGER, INTENT(IN)                      :: mysize
    CHARACTER(LEN=*), INTENT(IN)             :: storageform
    REAL(KIND=dp), INTENT(OUT)               :: eigenvalues( : ), &
                                                eigenvectors( :, : )

    CHARACTER                                :: jobz = "V", uplo
    INTEGER                                  :: info, lda, lwork
    REAL(KIND=dp)                            :: work( 3 * mysize - 1 )

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
     CALL stop_program ( "diagonalise_ssyev", "unknown form of storage" )
  END IF

  lda = SIZE ( matrix, 1 )
  lwork = 3 * mysize - 1

  eigenvectors = matrix

  CALL lapack_ssyev ( jobz, uplo, mysize, eigenvectors, lda, eigenvalues, &
       work, lwork, info )
  IF ( info /= 0 ) THEN
     CALL stop_program ( "diagonalise_ssyev", "error in diagonalisation" )
  END IF

END SUBROUTINE diagonalise_ssyev

!******************************************************************************

SUBROUTINE diagonalise_chpev ( matrix, mysize, storageform, eigenvalues, &
     eigenvectors )


    COMPLEX(KIND=dp), INTENT(INOUT)          :: matrix( : )
    INTEGER, INTENT(IN)                      :: mysize
    CHARACTER(LEN=*), INTENT(IN)             :: storageform
    REAL(KIND=dp), INTENT(OUT)               :: eigenvalues( : )
    COMPLEX(KIND=dp), INTENT(OUT)            :: eigenvectors( :, : )

    CHARACTER                                :: jobz = "V", uplo
    INTEGER                                  :: info
    COMPLEX(KIND=dp)                         :: work( 2 * mysize - 1 )
    REAL(KIND=dp)                            :: rwork( 3 * mysize - 2 )

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
     CALL stop_program ( "diagonalise_chpev", "unknown form of storage" )
  END IF

  CALL lapack_chpev ( jobz, uplo, mysize, matrix, eigenvalues, &
       eigenvectors, mysize, work, rwork, info )
  IF ( info /= 0 ) THEN
     CALL stop_program ( "diagonalise_chpev", "error in diagonalisation" )
  END IF

END SUBROUTINE diagonalise_chpev

!******************************************************************************

SUBROUTINE diagonalise_sspgv ( matrix, overlap, mysize, TYPE, storageform, &
     eigenvalues, eigenvectors )


    REAL(KIND=dp), INTENT(INOUT)             :: matrix( : ), overlap( : )
    INTEGER, INTENT(IN)                      :: mysize, TYPE
    CHARACTER(LEN=*), INTENT(IN)             :: storageform
    REAL(KIND=dp), INTENT(OUT)               :: eigenvalues( : ), &
                                                eigenvectors( :, : )

    CHARACTER                                :: jobz = "V", uplo
    INTEGER                                  :: info
    REAL(KIND=dp)                            :: work( 2 * mysize - 1 )

!----------------------------------------------------------------------

  IF ( storageform ( 1:5 ) == "Lower" .OR. &
       storageform ( 1:5 ) == "LOWER" .OR. &
       storageform ( 1:5 ) == "lower" ) THEN
     uplo = "L"
  ELSE IF ( storageform ( 1:5 ) == "Upper" .OR. &
       storageform ( 1:5 ) == "upper" .OR. &
       storageform ( 1:5 ) == "UPPER" ) THEN
     uplo = "U"
  ELSE
     CALL stop_program ( "diagonalise_sspgv", "unknown form of storage" )
  END IF

  CALL lapack_sspgv ( TYPE, jobz, uplo, mysize, matrix, overlap, eigenvalues, &
       eigenvectors, mysize, work, info )
  IF ( info /= 0 ) THEN
     CALL stop_program ( "diagonalise_sspgv", "error in diagonalisation" )
  END IF

END SUBROUTINE diagonalise_sspgv

!******************************************************************************

SUBROUTINE diagonalise_chpgv ( matrix, overlap, mysize, TYPE, storageform, &
     eigenvalues, eigenvectors )


    COMPLEX(KIND=dp), INTENT(INOUT)          :: matrix( : ), overlap( : )
    INTEGER, INTENT(IN)                      :: mysize, TYPE
    CHARACTER(LEN=*), INTENT(IN)             :: storageform
    REAL(KIND=dp), INTENT(OUT)               :: eigenvalues( : )
    COMPLEX(KIND=dp), INTENT(OUT)            :: eigenvectors( :, : )

    CHARACTER                                :: jobz = "V", uplo
    INTEGER                                  :: info
    COMPLEX(KIND=dp)                         :: work( 2 * mysize - 1 )
    REAL(KIND=dp)                            :: rwork( 3 * mysize - 2 )

!----------------------------------------------------------------------

  IF ( storageform ( 1:5 ) == "Lower" .OR. &
       storageform ( 1:5 ) == "LOWER" .OR. &
       storageform ( 1:5 ) == "lower" ) THEN
     uplo = "L"
  ELSE IF ( storageform ( 1:5 ) == "Upper" .OR. &
       storageform ( 1:5 ) == "upper" .OR. &
       storageform ( 1:5 ) == "UPPER" ) THEN
     uplo = "U"
  ELSE
     CALL stop_program ( "diagonalise_chpgv", "unknown form of storage" )
  END IF

  CALL lapack_chpgv ( TYPE, jobz, uplo, mysize, matrix, overlap, eigenvalues, &
       eigenvectors, mysize, work, rwork, info )
  IF ( info /= 0 ) THEN
     CALL stop_program ( "diagonalise_chpgv", "error in diagonalisation" )
  END IF

END SUBROUTINE diagonalise_chpgv

!******************************************************************************

SUBROUTINE cp2k_sgesvd ( matrix, svalues, mrow, ncol, uvec, vtvec )


    REAL(KIND=dp), INTENT(IN)                :: matrix( :, : )
    REAL(KIND=dp), INTENT(OUT)               :: svalues( : )
    INTEGER, INTENT(IN)                      :: mrow, ncol
    REAL(KIND=dp), INTENT(OUT)               :: uvec( :, : ), vtvec( :, : )

    CHARACTER, PARAMETER                     :: jobu = "A", jobvt = "A"

    INTEGER                                  :: info, lda, ldu, ldvt, lwork
    REAL(KIND=dp)                            :: work( 25 * ( mrow + ncol ) )

!----------------------------------------------------------------------

  lwork = 25 * ( mrow + ncol )
  lda = SIZE ( matrix, 1 )
  ldu = SIZE ( uvec, 1 )
  ldvt = SIZE ( vtvec, 1 )

  CALL lapack_sgesvd( jobu, jobvt, mrow, ncol, matrix, lda, svalues, &
        uvec, ldu, vtvec, ldvt, work, lwork, info )
  IF ( info /= 0 ) THEN
     CALL stop_program ( "cp2k_sgesvd", &
                         "error in singular value decomposition" )
  END IF

END SUBROUTINE cp2k_sgesvd

!******************************************************************************

SUBROUTINE cp2k_cgesvd ( matrix, svalues, mrow, ncol, uvec, vtvec )


    COMPLEX(KIND=dp), INTENT(IN)             :: matrix( :, : )
    REAL(KIND=dp), INTENT(OUT)               :: svalues( : )
    INTEGER, INTENT(IN)                      :: mrow, ncol
    COMPLEX(KIND=dp), INTENT(OUT)            :: uvec( :, : ), vtvec( :, : )

    CHARACTER, PARAMETER                     :: jobu = "A", jobvt = "A"

    INTEGER                                  :: info, lda, ldu, ldvt, lwork
    COMPLEX(KIND=dp)                         :: work( 25 * ( mrow + ncol ) )
    REAL(KIND=dp)                            :: rwork( 25 * ( mrow + ncol ) )

!----------------------------------------------------------------------

  lwork = 25 * ( mrow + ncol )
  lda = SIZE ( matrix, 1 )
  ldu = SIZE ( uvec, 1 )
  ldvt = SIZE ( vtvec, 1 )

  CALL lapack_cgesvd( jobu, jobvt, mrow, ncol, matrix, lda, svalues, &
        uvec, ldu, vtvec, ldvt, work, lwork, rwork, info )
  IF ( info /= 0 ) THEN
     CALL stop_program ( "cp2k_cgesvd", &
                         "error in singular value decomposition" )
  END IF

END SUBROUTINE cp2k_cgesvd

!******************************************************************************

END MODULE eigenvalueproblems

!******************************************************************************
