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

#if defined (__sp_lib)
  INTERFACE lapack_ssyev
     SUBROUTINE ssyev ( jobz, uplo, n, a, lda, w, work, lwork, info )
       USE kinds, ONLY : dbl
       CHARACTER :: jobz, uplo
       INTEGER :: INFO, LDA, LWORK, N
       REAL ( KIND = dbl ) :: A( LDA, * ), W( * ), WORK( * )
     END SUBROUTINE ssyev
  END INTERFACE

  INTERFACE lapack_chpev
     SUBROUTINE chpev ( jobz, uplo, n, ap, w, z, ldz, work, rwork, info )
       USE kinds, ONLY : dbl
       CHARACTER :: jobz, uplo
       INTEGER :: n, ldz, info
       REAL ( KIND = dbl ) :: w ( * ), rwork ( * )
       COMPLEX ( KIND = dbl ) :: ap ( * ), work ( * ), z ( ldz, * )
     END SUBROUTINE chpev
  END INTERFACE

  INTERFACE lapack_sspgv
     SUBROUTINE sspgv ( itype, jobz, uplo, n, ap, bp, w, z, ldz, work, &
          info )
       USE kinds, ONLY : dbl
       CHARACTER :: jobz, uplo
       INTEGER :: itype, n, ldz, info
       REAL ( KIND = dbl ) :: w ( * )
       COMPLEX ( KIND = dbl ) :: ap ( * ), bp ( * ), work ( * ), z ( ldz, * )
     END SUBROUTINE sspgv
  END INTERFACE

  INTERFACE lapack_chpgv
     SUBROUTINE chpgv ( itype, jobz, uplo, n, ap, bp, w, z, ldz, work, &
          rwork, info )
       USE kinds, ONLY : dbl
       CHARACTER :: jobz, uplo
       INTEGER :: itype, n, ldz, info
       REAL ( KIND = dbl ) :: w ( * ), rwork ( * )
       COMPLEX ( KIND = dbl ) :: ap ( * ), bp ( * ), work ( * ), z ( ldz, * )
     END SUBROUTINE chpgv
  END INTERFACE

  INTERFACE lapack_sgesvd
     SUBROUTINE sgesvd ( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, &
          work, lwork, info )
       USE kinds, ONLY : dbl
       USE kinds, ONLY : dbl
       CHARACTER :: jobz, jobvt
       INTEGER :: m, n, lda, ldu, ldvt, lwork, info
       REAL ( KIND = dbl ) :: a ( lda, * ), s ( * ), u ( ldu, * ), &
                              vt ( ldvt, * ), work ( * )
     END SUBROUTINE sgesvd
  END INTERFACE

  INTERFACE lapack_cgesvd
     SUBROUTINE cgesvd ( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, &
          work, lwork, rwork, info )
       USE kinds, ONLY : dbl
       USE kinds, ONLY : dbl
       CHARACTER :: jobz, jobvt
       INTEGER :: m, n, lda, ldu, ldvt, lwork, info
       REAL ( KIND = dbl ) :: s ( * ), rwork ( * )
       COMPLEX ( KIND = dbl ) :: a ( lda, * ), u ( ldu, * ), &
                                 vt ( ldvt, * ), work ( * )
     END SUBROUTINE cgesvd
  END INTERFACE

  INTERFACE lapack_spotrf
     SUBROUTINE spotrf ( uplo, n, a, lda, info )
       USE kinds, ONLY : dbl
       CHARACTER :: uplo
       INTEGER :: info, n, lda
       REAL ( KIND = dbl ) :: a( lda, * )
     END SUBROUTINE spotrf
  END INTERFACE

  INTERFACE lapack_cpotrf
     SUBROUTINE cpotrf ( uplo, n, a, lda, info )
       USE kinds, ONLY : dbl
       CHARACTER :: uplo
       INTEGER :: info, n, lda
       COMPLEX ( KIND = dbl ) :: a( lda, * )
     END SUBROUTINE cpotrf
  END INTERFACE

  INTERFACE lapack_strtri
     SUBROUTINE strtri ( uplo, diag, n, a, lda, info )
       USE kinds, ONLY : dbl
       CHARACTER :: uplo, diag
       INTEGER :: info, n, lda
       COMPLEX ( KIND = dbl ) :: a( lda, * )
     END SUBROUTINE strtri
  END INTERFACE

  INTERFACE lapack_ctrtri
     SUBROUTINE ctrtri ( uplo, diag, n, a, lda, info )
       USE kinds, ONLY : dbl
       CHARACTER :: uplo, diag
       INTEGER :: info, n, lda
       COMPLEX ( KIND = dbl ) :: a( lda, * )
     END SUBROUTINE ctrtri
  END INTERFACE

  INTERFACE lapack_sgesv
     SUBROUTINE sgesv ( n, nrhs, a, lda, ipiv, b, ldb, info )
       USE kinds, ONLY : dbl
       INTEGER :: INFO, LDA, LDB, N, NRHS
       INTEGER :: IPIV (*)
       REAL ( KIND = dbl ) :: A( LDA, * ), B( LDB, * )
     END SUBROUTINE sgesv
  END INTERFACE

#else

  INTERFACE lapack_ssyev
     SUBROUTINE dsyev ( jobz, uplo, n, a, lda, w, work, lwork, info )
       USE kinds, ONLY : dbl
       CHARACTER :: jobz, uplo
       INTEGER :: INFO, LDA, LWORK, N
       REAL ( KIND = dbl ) :: A( LDA, * ), W( * ), WORK( * )
     END SUBROUTINE dsyev
  END INTERFACE

  INTERFACE lapack_chpev
     SUBROUTINE zhpev ( jobz, uplo, n, ap, w, z, ldz, work, rwork, info )
       USE kinds, ONLY : dbl
       CHARACTER :: jobz, uplo
       INTEGER :: n, ldz, info
       REAL ( KIND = dbl ) :: w ( * ), rwork ( * )
       COMPLEX ( KIND = dbl ) :: ap ( * ), work ( * ), z ( ldz, * )
     END SUBROUTINE zhpev
  END INTERFACE

  INTERFACE lapack_sspgv
     SUBROUTINE dspgv ( itype, jobz, uplo, n, ap, bp, w, z, ldz, work, &
          info )
       USE kinds, ONLY : dbl
       CHARACTER :: jobz, uplo
       INTEGER :: itype, n, ldz, info
       REAL ( KIND = dbl ) :: w ( * )
       REAL ( KIND = dbl ) :: ap ( * ), bp ( * ), work ( * ), z ( ldz, * )
     END SUBROUTINE dspgv
  END INTERFACE

  INTERFACE lapack_chpgv
     SUBROUTINE zhpgv ( itype, jobz, uplo, n, ap, bp, w, z, ldz, work, &
          rwork, info )
       USE kinds, ONLY : dbl
       CHARACTER :: jobz, uplo
       INTEGER :: itype, n, ldz, info
       REAL ( KIND = dbl ) :: w ( * ), rwork ( * )
       COMPLEX ( KIND = dbl ) :: ap ( * ), bp ( * ), work ( * ), z ( ldz, * )
     END SUBROUTINE zhpgv
  END INTERFACE

  INTERFACE lapack_sgesvd
     SUBROUTINE dgesvd ( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, &
          work, lwork, info )
       USE kinds, ONLY : dbl
       CHARACTER :: jobu, jobvt
       INTEGER :: m, n, lda, ldu, ldvt, lwork, info
       REAL ( KIND = dbl ) :: a ( lda, * ), s ( * ), u ( ldu, * ), &
                              vt ( ldvt, * ), work ( * )
     END SUBROUTINE dgesvd
  END INTERFACE

  INTERFACE lapack_cgesvd
     SUBROUTINE zgesvd ( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, &
          work, lwork, rwork, info )
       USE kinds, ONLY : dbl
       CHARACTER :: jobu, jobvt
       INTEGER :: m, n, lda, ldu, ldvt, lwork, info
       REAL ( KIND = dbl ) :: s ( * ), rwork ( * )
       COMPLEX ( KIND = dbl ) :: a ( lda, * ), u ( ldu, * ), &
                                 vt ( ldvt, * ), work ( * )
     END SUBROUTINE zgesvd
  END INTERFACE

  INTERFACE lapack_spotrf
     SUBROUTINE dpotrf ( uplo, n, a, lda, info )
       USE kinds, ONLY : dbl
       CHARACTER :: uplo
       INTEGER :: info, n, lda
       REAL ( KIND = dbl ) :: a( lda, * )
     END SUBROUTINE dpotrf
  END INTERFACE

  INTERFACE lapack_cpotrf
     SUBROUTINE zpotrf ( uplo, n, a, lda, info )
       USE kinds, ONLY : dbl
       CHARACTER :: uplo
       INTEGER :: info, n, lda
       COMPLEX ( KIND = dbl ) :: a( lda, * )
     END SUBROUTINE zpotrf
  END INTERFACE

  INTERFACE lapack_strtri
     SUBROUTINE dtrtri ( uplo, diag, n, a, lda, info )
       USE kinds, ONLY : dbl
       CHARACTER :: uplo, diag
       INTEGER :: info, n, lda
       REAL ( KIND = dbl ) :: a( lda, * )
     END SUBROUTINE dtrtri
  END INTERFACE

  INTERFACE lapack_ctrtri
     SUBROUTINE ztrtri ( uplo, diag, n, a, lda, info )
       USE kinds, ONLY : dbl
       CHARACTER :: uplo, diag
       INTEGER :: info, n, lda
       COMPLEX ( KIND = dbl ) :: a( lda, * )
     END SUBROUTINE ztrtri
  END INTERFACE

  INTERFACE lapack_sgesv
     SUBROUTINE dgesv ( n, nrhs, a, lda, ipiv, b, ldb, info )
       USE kinds, ONLY : dbl
       INTEGER :: INFO, LDA, LDB, N, NRHS
       INTEGER :: IPIV (*)
       REAL ( KIND = dbl ) :: A( LDA, * ), B( LDB, * )
     END SUBROUTINE dgesv
  END INTERFACE
#endif

END MODULE lapack

!!*****
!******************************************************************************
