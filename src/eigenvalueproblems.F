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
!!
!!   AUTHOR
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

MODULE eigenvalueproblems
  
  USE kinds, ONLY : dbl
  USE lapack, ONLY : lapack_ssyev, lapack_chpev, lapack_chpgv
  USE stop_program, ONLY : stop_prg
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: diagonalise, diagonalise_overlap
  
  INTERFACE diagonalise
     MODULE PROCEDURE diagonalise_ssyev
     MODULE PROCEDURE diagonalise_dsyev
     MODULE PROCEDURE diagonalise_chpev
     MODULE PROCEDURE diagonalise_zhpev
  END INTERFACE
  
  INTERFACE diagonalise_overlap
     MODULE PROCEDURE diagonalise_chpgv
     MODULE PROCEDURE diagonalise_zhpgv
  END INTERFACE
  
!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE diagonalise_ssyev ( matrix, mysize, storageform, eigenvalues, &
     eigenvectors )
  
  IMPLICIT NONE
  
! Arguments
  CHARACTER ( LEN = * ), INTENT ( IN ) :: storageform
  INTEGER, INTENT ( IN ) :: mysize
  REAL ( KIND = 4 ), INTENT ( IN ) :: matrix ( :, : )
  REAL ( KIND = 4 ), INTENT ( OUT ) :: eigenvalues ( : )
  REAL ( KIND = 4 ), INTENT ( OUT ) :: eigenvectors ( :, : )
  
! Local
  CHARACTER :: jobz = "V", uplo
  INTEGER :: info, itype
  
  INTEGER :: lwork, lda
  REAL ( KIND = 4 ) :: work ( 3 * mysize - 1 )
  
!------------------------------------------------------------------------------
  
  IF ( storageform ( 1:5 ) == "Lower" .OR. &
       storageform ( 1:5 ) == "LOWER" .OR. &
       storageform ( 1:5 ) == "lower" ) THEN
     uplo = "L"
  ELSE IF ( storageform ( 1:5 ) == "Upper" .OR. &
       storageform ( 1:5 ) == "upper" .OR. &
       storageform ( 1:5 ) == "UPPER" ) then
     uplo = "U"
  ELSE
     CALL stop_prg ( "diagonalise_ssyev", "unknown form of storage" )
  END IF
  
  lda = SIZE ( matrix, 1 )
  lwork = 3 * mysize - 1
  
!*apsi* Because 'ssyev' is in-place...
  eigenvectors = matrix

  CALL lapack_ssyev ( jobz, uplo, mysize, eigenvectors, lda, eigenvalues, &
       work, lwork, info )
  IF ( info /= 0 ) THEN
     CALL stop_prg ( "diagonalise_ssyev", "error in diagonalisation" )
  END IF
  
END SUBROUTINE diagonalise_ssyev

!******************************************************************************

SUBROUTINE diagonalise_dsyev ( matrix, mysize, storageform, eigenvalues, &
     eigenvectors )
  
  IMPLICIT NONE
  
! Arguments
  CHARACTER ( LEN = * ), INTENT ( IN ) :: storageform
  INTEGER, INTENT ( IN ) :: mysize
  REAL ( KIND = 8 ), INTENT ( IN ) :: matrix ( :, : )
  REAL ( KIND = 8 ), INTENT ( OUT ) :: eigenvalues ( : )
  REAL ( KIND = 8 ), INTENT ( OUT ) :: eigenvectors ( :, : )
  
! Local
  CHARACTER :: jobz = "V", uplo
  INTEGER :: info, itype
  
  INTEGER :: lwork, lda
  REAL ( KIND = 8 ) :: work ( 3 * mysize - 1 )
  
!------------------------------------------------------------------------------
  
  IF ( storageform ( 1:5 ) == "Lower" .OR. &
       storageform ( 1:5 ) == "LOWER" .OR. &
       storageform ( 1:5 ) == "lower" ) THEN
     uplo = "L"
  ELSE IF ( storageform ( 1:5 ) == "Upper" .OR. &
       storageform ( 1:5 ) == "upper" .OR. &
       storageform ( 1:5 ) == "UPPER" ) then
     uplo = "U"
  ELSE
     CALL stop_prg ( "diagonalise_dsyev", "unknown form of storage" )
  END IF
  
  lda = SIZE ( matrix, 1 )
  lwork = 3 * mysize - 1
  
!*apsi* Because 'dsyev' is in-place...
  eigenvectors = matrix

  CALL lapack_ssyev ( jobz, uplo, mysize, eigenvectors, lda, eigenvalues, &
       work, lwork, info )
  IF ( info /= 0 ) THEN
     CALL stop_prg ( "diagonalise_dsyev", "error in diagonalisation" )
  END IF
  
END SUBROUTINE diagonalise_dsyev

!******************************************************************************

SUBROUTINE diagonalise_chpev ( matrix, mysize, storageform, eigenvalues, &
     eigenvectors )
  
  IMPLICIT NONE
  
! Arguments
  CHARACTER ( LEN = * ), INTENT ( IN ) :: storageform
  INTEGER, INTENT ( IN ) :: mysize
  REAL ( KIND = 4 ), INTENT ( OUT ) :: eigenvalues ( : )
  COMPLEX ( KIND = 4 ), INTENT ( OUT ) :: eigenvectors ( :, : )
  COMPLEX ( KIND = 4 ), INTENT ( INOUT ) :: matrix ( : )
  
! Local
  CHARACTER :: jobz = "V", uplo
  INTEGER :: info, itype
  
  REAL ( KIND = 4 ) :: rwork ( 3 * mysize - 2 )
  COMPLEX ( KIND = 4 ) :: work ( 2 * mysize - 1 )
  
!------------------------------------------------------------------------------
  
  IF ( storageform ( 1:5 ) == "Lower" .OR. &
       storageform ( 1:5 ) == "LOWER" .OR. &
       storageform ( 1:5 ) == "lower" ) THEN
     uplo = "L"
  ELSE IF ( storageform ( 1:5 ) == "Upper" .OR. &
       storageform ( 1:5 ) == "upper" .OR. &
       storageform ( 1:5 ) == "UPPER" ) then
     uplo = "U"
  ELSE
     CALL stop_prg ( "diagonalise_chpev", "unknown form of storage" )
  END IF
  
  CALL lapack_chpev ( jobz, uplo, mysize, matrix, eigenvalues, &
       eigenvectors, mysize, work, rwork, info )
  IF ( info /= 0 ) THEN
     CALL stop_prg ( "diagonalise_chpev", "error in diagonalisation" )
  END IF
  
END SUBROUTINE diagonalise_chpev

!******************************************************************************

SUBROUTINE diagonalise_zhpev ( matrix, mysize, storageform, eigenvalues, &
     eigenvectors )
  
  IMPLICIT NONE
  
! Arguments
  CHARACTER ( LEN = * ), INTENT ( IN ) :: storageform
  INTEGER, INTENT ( IN ) :: mysize
  REAL ( KIND = 8 ), INTENT ( OUT ) :: eigenvalues ( : )
  COMPLEX ( KIND = 8 ), INTENT ( OUT ) :: eigenvectors ( :, : )
  COMPLEX ( KIND = 8 ), INTENT ( INOUT ) :: matrix ( : )
  
! Local
  CHARACTER :: jobz = "V", uplo
  INTEGER :: info, itype
  
  REAL ( KIND = 8 ) :: rwork ( 3 * mysize - 2 )
  COMPLEX ( KIND = 8 ) :: work ( 2 * mysize - 1 )
  
!------------------------------------------------------------------------------
  
  IF ( storageform ( 1:5 ) == "Lower" .OR. &
       storageform ( 1:5 ) == "LOWER" .OR. &
       storageform ( 1:5 ) == "lower" ) THEN
     uplo = "L"
  ELSE IF ( storageform ( 1:5 ) == "Upper" .OR. &
       storageform ( 1:5 ) == "upper" .OR. &
       storageform ( 1:5 ) == "UPPER" ) then
     uplo = "U"
  ELSE
     CALL stop_prg ( "diagonalise_chpev", "unknown form of storage" )
  END IF
  
  CALL lapack_chpev ( jobz, uplo, mysize, matrix, eigenvalues, &
       eigenvectors, mysize, work, rwork, info )
  IF ( info /= 0 ) THEN
     CALL stop_prg ( "diagonalise_zhpev", "error in diagonalisation" )
  END IF
  
END SUBROUTINE diagonalise_zhpev

!**********************************************************************

SUBROUTINE diagonalise_chpgv ( matrix, overlap, mysize, type, storageform, &
     eigenvalues, eigenvectors )
  
  implicit none
  
! Arguments
  CHARACTER ( LEN = * ), INTENT ( IN ) :: storageform
  INTEGER, INTENT ( IN ) :: type, mysize
  REAL ( KIND = 4 ), INTENT ( OUT ) :: eigenvalues ( : )
  COMPLEX ( KIND = 4 ), INTENT ( OUT ) :: eigenvectors ( :, : )
  COMPLEX ( KIND = 4 ), INTENT ( INOUT ) :: matrix ( : ), overlap ( : )
  
! Local
  CHARACTER :: jobz = "V", uplo
  integer :: info
  real ( kind = 4 ) :: rwork ( 3 * mysize - 2 )
  COMPLEX ( kind = 4 ) :: work ( 2 * mysize - 1 )
  
!----------------------------------------------------------------------
  
  IF ( storageform ( 1:5 ) == "Lower" .OR. &
       storageform ( 1:5 ) == "LOWER" .OR. &
       storageform ( 1:5 ) == "lower" ) THEN
     uplo = "L"
  ELSE IF ( storageform ( 1:5 ) == "Upper" .OR. &
       storageform ( 1:5 ) == "upper" .OR. &
       storageform ( 1:5 ) == "UPPER" ) then
     uplo = "U"
  ELSE
     stop "diagonalise_chpev, unknown form of storage"
  END IF
  
  CALL lapack_chpgv ( type, jobz, uplo, mysize, matrix, overlap, eigenvalues, &
       eigenvectors, mysize, work, rwork, info )
  IF ( info /= 0 ) THEN
     CALL stop_prg ( "diagonalise_chpgv", "error in diagonalisation" )
  END IF

END SUBROUTINE diagonalise_chpgv

!******************************************************************************

SUBROUTINE diagonalise_zhpgv ( matrix, overlap, mysize, type, storageform, &
     eigenvalues, eigenvectors )

  IMPLICIT NONE

! Arguments
  CHARACTER ( LEN = * ), INTENT ( IN ) :: storageform
  INTEGER, INTENT ( IN ) :: type, mysize
  REAL ( KIND = 8 ), INTENT ( OUT ) :: eigenvalues ( : )
  COMPLEX ( KIND = 8 ), INTENT ( OUT ) :: eigenvectors ( :, : )
  COMPLEX ( KIND = 8 ), INTENT ( INOUT ) :: matrix ( : )
  COMPLEX ( KIND = 8 ), INTENT ( INOUT ) :: overlap ( : )
  
! Local
  CHARACTER :: jobz = "V", uplo
  INTEGER :: info
  REAL ( KIND = 8 ) :: rwork ( 3 * mysize - 2 )
  COMPLEX ( KIND = 8 ) :: work ( 2 * mysize - 1 )

!------------------------------------------------------------------------------

  IF ( storageform ( 1:5 ) == "Lower" .OR. &
       storageform ( 1:5 ) == "LOWER" .OR. &
       storageform ( 1:5 ) == "lower" ) THEN
     uplo = "L"
  ELSE IF ( storageform ( 1:5 ) == "Upper" .OR. &
       storageform ( 1:5 ) == "upper" .OR. &
       storageform ( 1:5 ) == "UPPER" ) then
     uplo = "U"
  ELSE
     stop "diagonalise_chpev, unknown form of storage"
  END IF

  CALL lapack_chpgv ( type, jobz, uplo, mysize, matrix, overlap, eigenvalues, &
       eigenvectors, mysize, work, rwork, info )
  IF ( info /= 0 ) THEN
     CALL stop_prg ( "diagonalise_zhpgv", "error in diagonalisation" )
  END IF

END SUBROUTINE diagonalise_zhpgv

END MODULE eigenvalueproblems
