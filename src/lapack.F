!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Interface to the LAPACK F77 library
!> \par History
!>      JGH (26-5-2001): New flag __SGL
!>                     : Added SVD routines
!>      JGH (30-5-2001): Added Cholesky factorisation and triangular inverse
!>      JGH (23-1-2009): Added linear system solution by SVD
!> \author APSI
! *****************************************************************************
MODULE lapack
  USE f77_blas

  IMPLICIT NONE
  PUBLIC ! make all routines public

  INTERFACE lapack_ssyev
! *****************************************************************************
     SUBROUTINE dsyev ( jobz, uplo, n, a, lda, w, work, lwork, info )
       USE kinds, ONLY: dp
    CHARACTER                                :: jobz, uplo
    INTEGER                                  :: N, LDA
    REAL(KIND=dp)                            :: A( LDA, * ), W( * ), WORK( * )
    INTEGER                                  :: LWORK, INFO

     END SUBROUTINE dsyev
  END INTERFACE

  INTERFACE lapack_ssygv
! *****************************************************************************
     SUBROUTINE dsygv ( itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info )
       USE kinds, ONLY: dp
    INTEGER                                  :: itype
    CHARACTER                                :: jobz, uplo
    INTEGER                                  :: N, LDA
    REAL(KIND=dp)                            :: A( LDA, * )
    INTEGER                                  :: LDB
    REAL(KIND=dp)                            :: B( LDB, * ), W( * ), WORK( * )
    INTEGER                                  :: LWORK, INFO

     END SUBROUTINE dsygv
  END INTERFACE

  INTERFACE lapack_chpev
! *****************************************************************************
     SUBROUTINE zhpev ( jobz, uplo, n, ap, w, z, ldz, work, rwork, info )
       USE kinds, ONLY: dp
    CHARACTER                                :: jobz, uplo
    INTEGER                                  :: n
    COMPLEX(KIND=dp)                         :: ap( * )
    REAL(KIND=dp)                            :: w( * )
    INTEGER                                  :: ldz
    COMPLEX(KIND=dp)                         :: z( ldz, * ), work( * )
    REAL(KIND=dp)                            :: rwork( * )
    INTEGER                                  :: info

     END SUBROUTINE zhpev
  END INTERFACE

  INTERFACE lapack_sspgv
! *****************************************************************************
     SUBROUTINE dspgv ( itype, jobz, uplo, n, ap, bp, w, z, ldz, work, &
          info )
       USE kinds, ONLY: dp
    INTEGER                                  :: itype
    CHARACTER                                :: jobz, uplo
    INTEGER                                  :: n
    REAL(KIND=dp)                            :: ap( * ), bp( * ), w( * )
    INTEGER                                  :: ldz
    REAL(KIND=dp)                            :: z( ldz, * ), work( * )
    INTEGER                                  :: info

     END SUBROUTINE dspgv
  END INTERFACE

  INTERFACE lapack_chpgv
! *****************************************************************************
     SUBROUTINE zhpgv ( itype, jobz, uplo, n, ap, bp, w, z, ldz, work, &
          rwork, info )
       USE kinds, ONLY: dp
    INTEGER                                  :: itype
    CHARACTER                                :: jobz, uplo
    INTEGER                                  :: n
    COMPLEX(KIND=dp)                         :: ap( * ), bp( * )
    REAL(KIND=dp)                            :: w( * )
    INTEGER                                  :: ldz
    COMPLEX(KIND=dp)                         :: z( ldz, * ), work( * )
    REAL(KIND=dp)                            :: rwork( * )
    INTEGER                                  :: info

     END SUBROUTINE zhpgv
  END INTERFACE

  INTERFACE lapack_sgesvd
! *****************************************************************************
     SUBROUTINE dgesvd ( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, &
          work, lwork, info )
       USE kinds, ONLY: dp
    CHARACTER                                :: jobu, jobvt
    INTEGER                                  :: m, n, lda
    REAL(KIND=dp)                            :: a( lda, * ), s( * )
    INTEGER                                  :: ldu
    REAL(KIND=dp)                            :: u( ldu, * )
    INTEGER                                  :: ldvt
    REAL(KIND=dp)                            :: vt( ldvt, * ), work( * )
    INTEGER                                  :: lwork, info

     END SUBROUTINE dgesvd
  END INTERFACE

  INTERFACE lapack_cgesvd
! *****************************************************************************
     SUBROUTINE zgesvd ( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, &
          work, lwork, rwork, info )
       USE kinds, ONLY: dp
    CHARACTER                                :: jobu, jobvt
    INTEGER                                  :: m, n, lda
    COMPLEX(KIND=dp)                         :: a( lda, * )
    REAL(KIND=dp)                            :: s( * )
    INTEGER                                  :: ldu
    COMPLEX(KIND=dp)                         :: u( ldu, * )
    INTEGER                                  :: ldvt
    COMPLEX(KIND=dp)                         :: vt( ldvt, * ), work( * )
    INTEGER                                  :: lwork
    REAL(KIND=dp)                            :: rwork( * )
    INTEGER                                  :: info

     END SUBROUTINE zgesvd
  END INTERFACE

  INTERFACE lapack_spotrf
! *****************************************************************************
     SUBROUTINE dpotrf ( uplo, n, a, lda, info )
       USE kinds, ONLY: dp
    CHARACTER                                :: uplo
    INTEGER                                  :: n, lda
    REAL(KIND=dp)                            :: a( lda, * )
    INTEGER                                  :: info

     END SUBROUTINE dpotrf
  END INTERFACE

  INTERFACE lapack_cpotrf
! *****************************************************************************
     SUBROUTINE zpotrf ( uplo, n, a, lda, info )
       USE kinds, ONLY: dp
    CHARACTER                                :: uplo
    INTEGER                                  :: n, lda
    COMPLEX(KIND=dp)                         :: a( lda, * )
    INTEGER                                  :: info

     END SUBROUTINE zpotrf
  END INTERFACE

  INTERFACE lapack_strtri
! *****************************************************************************
     SUBROUTINE dtrtri ( uplo, diag, n, a, lda, info )
       USE kinds, ONLY: dp
    CHARACTER                                :: uplo, diag
    INTEGER                                  :: n, lda
    REAL(KIND=dp)                            :: a( lda, * )
    INTEGER                                  :: info

     END SUBROUTINE dtrtri
  END INTERFACE

  INTERFACE lapack_ctrtri
! *****************************************************************************
     SUBROUTINE ztrtri ( uplo, diag, n, a, lda, info )
       USE kinds, ONLY: dp
    CHARACTER                                :: uplo, diag
    INTEGER                                  :: n, lda
    COMPLEX(KIND=dp)                         :: a( lda, * )
    INTEGER                                  :: info

     END SUBROUTINE ztrtri
  END INTERFACE

  INTERFACE lapack_sgesv
! *****************************************************************************
     SUBROUTINE dgesv ( n, nrhs, a, lda, ipiv, b, ldb, info )
       USE kinds, ONLY: dp
    INTEGER                                  :: N, NRHS, LDA
    REAL(KIND=dp)                            :: A( LDA, * )
    INTEGER                                  :: IPIV(*), LDB
    REAL(KIND=dp)                            :: B( LDB, * )
    INTEGER                                  :: INFO

     END SUBROUTINE dgesv
  END INTERFACE

  INTERFACE lapack_sgelss
! *****************************************************************************
     SUBROUTINE dgelss ( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info )
       USE kinds, ONLY: dp
    INTEGER                                  :: m, n, nrhs, lda
    REAL(KIND=dp)                            :: a( lda, * )
    INTEGER                                  :: ldb
    REAL(KIND=dp)                            :: b( ldb, *), s(*), rcond
    INTEGER                                  :: rank
    REAL(KIND=dp)                            :: work( * )
    INTEGER                                  :: lwork, info

     END SUBROUTINE dgelss
  END INTERFACE

END MODULE lapack

