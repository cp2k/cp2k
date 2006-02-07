!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/lapack [1.0] *
!!
!!   NAME
!!     lapack
!!
!!   FUNCTION
!!     Interface to the LAPACK F77 library
!!
!!   AUTHOR
!!     APSI
!!
!!   MODIFICATION HISTORY
!!     JGH (26-5-2001): New flag __sp_lib
!!                    : Added SVD routines
!!     JGH (30-5-2001): Added Cholesky facotrisation and triangular inverse
!!
!!   SOURCE
!******************************************************************************

MODULE lapack
  IMPLICIT NONE
  PUBLIC ! make all routines public 
#if defined (__sp_lib)
  INTERFACE lapack_ssyev
     SUBROUTINE ssyev ( jobz, uplo, n, a, lda, w, work, lwork, info )
       USE kinds, ONLY: dp
    CHARACTER                                :: jobz, uplo
    INTEGER                                  :: N, LDA
    REAL(KIND=dp)                            :: A( LDA, * ), W( * ), WORK( * )
    INTEGER                                  :: LWORK, INFO

     END SUBROUTINE ssyev
  END INTERFACE

  INTERFACE lapack_chpev
     SUBROUTINE chpev ( jobz, uplo, n, ap, w, z, ldz, work, rwork, info )
       USE kinds, ONLY: dp
    CHARACTER                                :: jobz, uplo
    INTEGER                                  :: n
    COMPLEX(KIND=dp)                         :: ap( * )
    REAL(KIND=dp)                            :: w( * )
    INTEGER                                  :: ldz
    COMPLEX(KIND=dp)                         :: z( ldz, * ), work( * )
    REAL(KIND=dp)                            :: rwork( * )
    INTEGER                                  :: info

     END SUBROUTINE chpev
  END INTERFACE

  INTERFACE lapack_sspgv
     SUBROUTINE sspgv ( itype, jobz, uplo, n, ap, bp, w, z, ldz, work, &
          info )
       USE kinds, ONLY: dp
    INTEGER                                  :: itype
    CHARACTER                                :: jobz, uplo
    INTEGER                                  :: n
    COMPLEX(KIND=dp)                         :: ap( * ), bp( * )
    REAL(KIND=dp)                            :: w( * )
    INTEGER                                  :: ldz
    COMPLEX(KIND=dp)                         :: z( ldz, * ), work( * )
    INTEGER                                  :: info

     END SUBROUTINE sspgv
  END INTERFACE

  INTERFACE lapack_chpgv
     SUBROUTINE chpgv ( itype, jobz, uplo, n, ap, bp, w, z, ldz, work, &
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

     END SUBROUTINE chpgv
  END INTERFACE

  INTERFACE lapack_sgesvd
     SUBROUTINE sgesvd ( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, &
          work, lwork, info )
       USE kinds, ONLY: dp
       USE kinds, ONLY: dp
    CHARACTER                                :: jobu, jobvt
    INTEGER                                  :: m, n, lda
    REAL(KIND=dp)                            :: a( lda, * ), s( * )
    INTEGER                                  :: ldu
    REAL(KIND=dp)                            :: u( ldu, * )
    INTEGER                                  :: ldvt
    REAL(KIND=dp)                            :: vt( ldvt, * ), work( * )
    INTEGER                                  :: lwork, info

     END SUBROUTINE sgesvd
  END INTERFACE

  INTERFACE lapack_cgesvd
     SUBROUTINE cgesvd ( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, &
          work, lwork, rwork, info )
       USE kinds, ONLY: dp
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

     END SUBROUTINE cgesvd
  END INTERFACE

  INTERFACE lapack_spotrf
     SUBROUTINE spotrf ( uplo, n, a, lda, info )
       USE kinds, ONLY: dp
    CHARACTER                                :: uplo
    INTEGER                                  :: n, lda
    REAL(KIND=dp)                            :: a( lda, * )
    INTEGER                                  :: info

     END SUBROUTINE spotrf
  END INTERFACE

  INTERFACE lapack_cpotrf
     SUBROUTINE cpotrf ( uplo, n, a, lda, info )
       USE kinds, ONLY: dp
    CHARACTER                                :: uplo
    INTEGER                                  :: n, lda
    COMPLEX(KIND=dp)                         :: a( lda, * )
    INTEGER                                  :: info

     END SUBROUTINE cpotrf
  END INTERFACE

  INTERFACE lapack_strtri
     SUBROUTINE strtri ( uplo, diag, n, a, lda, info )
       USE kinds, ONLY: dp
    CHARACTER                                :: uplo, diag
    INTEGER                                  :: n, lda
    COMPLEX(KIND=dp)                         :: a( lda, * )
    INTEGER                                  :: info

     END SUBROUTINE strtri
  END INTERFACE

  INTERFACE lapack_ctrtri
     SUBROUTINE ctrtri ( uplo, diag, n, a, lda, info )
       USE kinds, ONLY: dp
    CHARACTER                                :: uplo, diag
    INTEGER                                  :: n, lda
    COMPLEX(KIND=dp)                         :: a( lda, * )
    INTEGER                                  :: info

     END SUBROUTINE ctrtri
  END INTERFACE

  INTERFACE lapack_sgesv
     SUBROUTINE sgesv ( n, nrhs, a, lda, ipiv, b, ldb, info )
       USE kinds, ONLY: dp
    INTEGER                                  :: N, NRHS, LDA
    REAL(KIND=dp)                            :: A( LDA, * )
    INTEGER                                  :: IPIV(*), LDB
    REAL(KIND=dp)                            :: B( LDB, * )
    INTEGER                                  :: INFO

     END SUBROUTINE sgesv
  END INTERFACE

#else

  INTERFACE lapack_ssyev
     SUBROUTINE dsyev ( jobz, uplo, n, a, lda, w, work, lwork, info )
       USE kinds, ONLY: dp
    CHARACTER                                :: jobz, uplo
    INTEGER                                  :: N, LDA
    REAL(KIND=dp)                            :: A( LDA, * ), W( * ), WORK( * )
    INTEGER                                  :: LWORK, INFO

     END SUBROUTINE dsyev
  END INTERFACE

  INTERFACE lapack_chpev
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
     SUBROUTINE dpotrf ( uplo, n, a, lda, info )
       USE kinds, ONLY: dp
    CHARACTER                                :: uplo
    INTEGER                                  :: n, lda
    REAL(KIND=dp)                            :: a( lda, * )
    INTEGER                                  :: info

     END SUBROUTINE dpotrf
  END INTERFACE

  INTERFACE lapack_cpotrf
     SUBROUTINE zpotrf ( uplo, n, a, lda, info )
       USE kinds, ONLY: dp
    CHARACTER                                :: uplo
    INTEGER                                  :: n, lda
    COMPLEX(KIND=dp)                         :: a( lda, * )
    INTEGER                                  :: info

     END SUBROUTINE zpotrf
  END INTERFACE

  INTERFACE lapack_strtri
     SUBROUTINE dtrtri ( uplo, diag, n, a, lda, info )
       USE kinds, ONLY: dp
    CHARACTER                                :: uplo, diag
    INTEGER                                  :: n, lda
    REAL(KIND=dp)                            :: a( lda, * )
    INTEGER                                  :: info

     END SUBROUTINE dtrtri
  END INTERFACE

  INTERFACE lapack_ctrtri
     SUBROUTINE ztrtri ( uplo, diag, n, a, lda, info )
       USE kinds, ONLY: dp
    CHARACTER                                :: uplo, diag
    INTEGER                                  :: n, lda
    COMPLEX(KIND=dp)                         :: a( lda, * )
    INTEGER                                  :: info

     END SUBROUTINE ztrtri
  END INTERFACE

  INTERFACE lapack_sgesv
     SUBROUTINE dgesv ( n, nrhs, a, lda, ipiv, b, ldb, info )
       USE kinds, ONLY: dp
    INTEGER                                  :: N, NRHS, LDA
    REAL(KIND=dp)                            :: A( LDA, * )
    INTEGER                                  :: IPIV(*), LDB
    REAL(KIND=dp)                            :: B( LDB, * )
    INTEGER                                  :: INFO

     END SUBROUTINE dgesv
  END INTERFACE
#endif

END MODULE lapack

!!*****
!******************************************************************************
