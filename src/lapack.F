!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#if ! defined ( __YMP ) && ! defined ( __T3E )

!! We are using LAPACK interfaces, so please make sure in IBM/AIX you have
!!    the lapack library before essl: "xlf90 ... -llapack -lessl" !!!

MODULE lapack
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: lapack_ssyev, lapack_chpev, lapack_chpgv
  
  INTERFACE lapack_ssyev
     SUBROUTINE ssyev ( jobz, uplo, n, a, lda, w, work, lwork, info )
       CHARACTER :: jobz, uplo
       INTEGER :: INFO, LDA, LWORK, N
       REAL ( KIND = 4 ) :: A( LDA, * ), W( * ), WORK( * )
     END SUBROUTINE ssyev
     
     SUBROUTINE dsyev ( jobz, uplo, n, a, lda, w, work, lwork, info )
       CHARACTER :: jobz, uplo
       INTEGER :: INFO, LDA, LWORK, N
       REAL ( KIND = 8 ) :: A( LDA, * ), W( * ), WORK( * )
     END SUBROUTINE dsyev
  END INTERFACE
  
  INTERFACE lapack_chpev
     SUBROUTINE chpev ( jobz, uplo, n, ap, w, z, ldz, work, rwork, info )
       CHARACTER :: jobz, uplo
       INTEGER :: n, ldz, info
       REAL ( KIND = 4 ) :: w ( * ), rwork ( * )
       COMPLEX ( KIND = 4 ) :: ap ( * ), work ( * ), z ( ldz, * )
     END SUBROUTINE chpev
     
     SUBROUTINE zhpev ( jobz, uplo, n, ap, w, z, ldz, work, rwork, info )
       CHARACTER :: jobz, uplo
       INTEGER :: n, ldz, info
       REAL ( KIND = 8 ) :: w ( * ), rwork ( * )
       COMPLEX ( KIND = 8 ) :: ap ( * ), work ( * ), z ( ldz, * )
     END SUBROUTINE zhpev
  END INTERFACE
  
  INTERFACE lapack_chpgv
     SUBROUTINE chpgv ( itype, jobz, uplo, n, ap, bp, w, z, ldz, work, &
          rwork, info )
       CHARACTER :: jobz, uplo
       INTEGER :: itype, n, ldz, info
       REAL ( KIND = 4 ) :: w ( * ), rwork ( * )
       COMPLEX ( KIND = 4 ) :: ap ( * ), bp ( * ), work ( * ), z ( ldz, * )
     END SUBROUTINE chpgv
     
     SUBROUTINE zhpgv ( itype, jobz, uplo, n, ap, bp, w, z, ldz, work, &
          rwork, info )
       CHARACTER :: jobz, uplo
       INTEGER :: itype, n, ldz, info
       REAL ( KIND = 8 ) :: w ( * ), rwork ( * )
       COMPLEX ( KIND = 8 ) :: ap ( * ), bp ( * ), work ( * ), z ( ldz, * )
     END SUBROUTINE zhpgv
  END INTERFACE
  
END MODULE lapack

#else

MODULE lapack
  
  INTERFACE lapack_ssyev
     SUBROUTINE gsyev ( jobz, uplo, n, a, lda, w, work, lwork, info )
       CHARACTER :: jobz, uplo
       INTEGER :: INFO, LDA, LWORK, N
       REAL ( KIND = 4 ) :: A( LDA, * ), W( * ), WORK( * )
     END SUBROUTINE gsyev
     
     SUBROUTINE ssyev ( jobz, uplo, n, a, lda, w, work, lwork, info )
       CHARACTER :: jobz, uplo
       INTEGER :: INFO, LDA, LWORK, N
       REAL ( KIND = 8 ) :: A( LDA, * ), W( * ), WORK( * )
     END SUBROUTINE ssyev
  END INTERFACE
  
  INTERFACE lapack_chpev
     SUBROUTINE ghpev ( jobz, uplo, n, ap, w, z, ldz, work, rwork, info )
       CHARACTER :: jobz, uplo
       INTEGER :: n, ldz, info
       REAL ( KIND = 4 ) :: w ( * ), rwork ( * )
       COMPLEX ( KIND = 4 ) :: ap ( * ), work ( * ), z ( ldz, * )
     END SUBROUTINE ghpev
     
     SUBROUTINE chpev ( jobz, uplo, n, ap, w, z, ldz, work, rwork, info )
       CHARACTER :: jobz, uplo
       INTEGER :: n, ldz, info
       REAL ( KIND = 8 ) :: w ( * ), rwork ( * )
       COMPLEX ( KIND = 8 ) :: ap ( * ), work ( * ), z ( ldz, * )
     END SUBROUTINE chpev
  END INTERFACE
  
  INTERFACE lapack_chpgv
     SUBROUTINE ghpgv ( itype, jobz, uplo, n, ap, bp, w, z, ldz, work, &
          rwork, info )
       CHARACTER :: jobz, uplo
       INTEGER :: itype, n, ldz, info
       REAL ( KIND = 4 ) :: w ( * ), rwork ( * )
       COMPLEX ( KIND = 4 ) :: ap ( * ), bp ( * ), work ( * ), z ( ldz, * )
     END SUBROUTINE ghpgv
     
     SUBROUTINE chpgv ( itype, jobz, uplo, n, ap, bp, w, z, ldz, work, &
          rwork, info )
       CHARACTER :: jobz, uplo
       INTEGER :: itype, n, ldz, info
       REAL ( KIND = 8 ) :: w ( * ), rwork ( * )
       COMPLEX ( KIND = 8 ) :: ap ( * ), bp ( * ), work ( * ), z ( ldz, * )
     END SUBROUTINE chpgv
  END INTERFACE

END MODULE lapack

#endif
