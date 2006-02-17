MODULE f77_blas_netlib
!
! Provides interfaces for BLAS as distributed
! (http://www.netlib.org/blas/) Jan 2006.
!
! follows the same naming convention as the quick ref
! http://www.netlib.org/lapack/lug/node145.html
! for the arguments, and similar to the f77_blas_generic module
!
! Note: All procedures have intents specified for their arguments.
!       These are the de facto intents of these arguments, as found
!       in the netlib source/comments. Nevertheless, since the netlib
!       blas implementation is in Fortran77, INTENTs are
!       not specified. Some compilers might produce an error/warning (if strict
!       standard conformance checking is enabled). In that case, one option
!       is to edit the files and remove all INTENT declarations. E.g. :
!       sed 's/\(, INTENT(INOUT)\)\|\(, INTENT(OUT)\)\|\(, INTENT(IN)\)//' f77_blas_netlib.f90 > new.f90
!
! this file is the public domain. Joost VandeVondele.
!

PRIVATE

PUBLIC :: SROTG, DROTG, CROTG, ZROTG, SROTMG, DROTMG, SROT, DROT, CSROT, ZDROT, SROTM, DROTM, &
          SSWAP, DSWAP, CSWAP, ZSWAP, SSCAL, DSCAL, CSCAL, ZSCAL, CSSCAL, ZDSCAL, &
          SCOPY, DCOPY, CCOPY, ZCOPY, SAXPY, DAXPY, CAXPY, ZAXPY, SDOT, DDOT,  &
          CDOTU, ZDOTU, CDOTC, ZDOTC, SNRM2, DNRM2, SCNRM2, DZNRM2, &
          SASUM, SCASUM, DASUM, DZASUM, ISAMAX, IDAMAX, ICAMAX, IZAMAX

PUBLIC :: SGEMV, DGEMV, CGEMV, ZGEMV, SGBMV, DGBMV, CGBMV, ZGBMV, &
          CHEMV, ZHEMV, CHBMV, ZHBMV, CHPMV, ZHPMV, SSYMV, DSYMV, &
          SSBMV, DSBMV, SSPMV, DSPMV, STRMV, DTRMV, CTRMV, ZTRMV, &
          STBMV, DTBMV, CTBMV, ZTBMV, STPMV, DTPMV, CTPMV, ZTPMV, &
          STRSV, DTRSV, CTRSV, ZTRSV, STBSV, DTBSV, CTBSV, ZTBSV, &  
          STPSV, DTPSV, CTPSV, ZTPSV, SGER, DGER, CGERU, ZGERU,   &
          CGERC, ZGERC, CHER, ZHER, CHPR, ZHPR, CHER2, ZHER2,     &
          CHPR2, ZHPR2, SSYR, DSYR, SSPR, DSPR, SSYR2, DSYR2,     &
          SSPR2, DSPR2

PUBLIC :: SGEMM, DGEMM, CGEMM, ZGEMM, SSYMM, DSYMM, CSYMM, ZSYMM, CHEMM, ZHEMM, &
          SSYRK, DSYRK, CSYRK, ZSYRK, CHERK, ZHERK, SSYR2K, DSYR2K, CSYR2K, ZSYR2K, CHER2K, ZHER2K, &
          STRMM, DTRMM, CTRMM, ZTRMM, STRSM, DTRSM, CTRSM, ZTRSM

PUBLIC :: SDSDOT, DSDOT, DCABS1


!
! Level 1 BLAS
!
INTERFACE
   SUBROUTINE SROTG(A,B,C,S)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    REAL(SP), INTENT(INOUT)                  :: A, B
    REAL(SP), INTENT(OUT)                    :: C, S
   END SUBROUTINE SROTG
END INTERFACE
INTERFACE
   SUBROUTINE DROTG(A,B,C,S)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    REAL(DP), INTENT(INOUT)                  :: A, B
    REAL(DP), INTENT(OUT)                    :: C, S
   END SUBROUTINE DROTG
END INTERFACE
INTERFACE
   SUBROUTINE CROTG(A,B,C,S)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    COMPLEX(SP), INTENT(INOUT)               :: A, B
    REAL(SP), INTENT(OUT)                    :: C
    COMPLEX(SP), INTENT(OUT)                 :: S
   END SUBROUTINE CROTG
END INTERFACE
INTERFACE
   SUBROUTINE ZROTG(A,B,C,S)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    COMPLEX(DP), INTENT(INOUT)               :: A, B
    REAL(DP), INTENT(OUT)                    :: C
    COMPLEX(DP), INTENT(OUT)                 :: S
   END SUBROUTINE ZROTG
END INTERFACE
INTERFACE
   SUBROUTINE SROTMG(D1,D2,A,B,PARAM)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    REAL(SP), INTENT(INOUT)                  :: D1, D2, A
    REAL(SP), INTENT(IN)                     :: B
    REAL(SP), INTENT(OUT)                    :: PARAM(5)
   END SUBROUTINE SROTMG
END INTERFACE
INTERFACE
   SUBROUTINE DROTMG(D1,D2,A,B,PARAM)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    REAL(DP), INTENT(INOUT)                  :: D1, D2, A
    REAL(DP), INTENT(IN)                     :: B
    REAL(DP), INTENT(OUT)                    :: PARAM(5)
   END SUBROUTINE DROTMG
END INTERFACE
INTERFACE
   SUBROUTINE SROT(N,X,INCX,Y,INCY,C,S)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(SP), INTENT(IN)                     :: C, S
   END SUBROUTINE SROT
END INTERFACE
INTERFACE
   SUBROUTINE DROT(N,X,INCX,Y,INCY,C,S)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(DP), INTENT(IN)                     :: C, S
   END SUBROUTINE DROT
END INTERFACE
INTERFACE
   SUBROUTINE CSROT(N,X,INCX,Y,INCY,C,S)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(SP), INTENT(IN)                     :: C, S
   END SUBROUTINE CSROT
END INTERFACE
INTERFACE
   SUBROUTINE ZDROT(N,X,INCX,Y,INCY,C,S)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(DP), INTENT(IN)                     :: C, S
   END SUBROUTINE ZDROT
END INTERFACE
INTERFACE
   SUBROUTINE SROTM(N,X,INCX,Y,INCY,PARAM)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(SP), INTENT(IN)                     :: PARAM(5)
   END SUBROUTINE SROTM
END INTERFACE
INTERFACE
   SUBROUTINE DROTM(N,X,INCX,Y,INCY,PARAM)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(DP), INTENT(IN)                     :: PARAM(5)
   END SUBROUTINE DROTM
END INTERFACE
INTERFACE
   SUBROUTINE SSWAP(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE SSWAP
END INTERFACE
INTERFACE
   SUBROUTINE DSWAP(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE DSWAP
END INTERFACE
INTERFACE
   SUBROUTINE CSWAP(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE CSWAP
END INTERFACE
INTERFACE
   SUBROUTINE ZSWAP(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE ZSWAP
END INTERFACE
INTERFACE
   SUBROUTINE SSCAL(N,ALPHA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: ALPHA
    REAL(SP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE SSCAL
END INTERFACE
INTERFACE
   SUBROUTINE DSCAL(N,ALPHA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: ALPHA
    REAL(DP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE DSCAL
END INTERFACE
INTERFACE
   SUBROUTINE CSCAL(N,ALPHA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: ALPHA
    COMPLEX(SP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE CSCAL
END INTERFACE
INTERFACE
   SUBROUTINE ZSCAL(N,ALPHA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: ALPHA
    COMPLEX(DP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE ZSCAL
END INTERFACE
INTERFACE
   SUBROUTINE CSSCAL(N,ALPHA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: ALPHA
    COMPLEX(SP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE CSSCAL
END INTERFACE
INTERFACE
   SUBROUTINE ZDSCAL(N,ALPHA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: ALPHA
    COMPLEX(DP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE ZDSCAL
END INTERFACE
INTERFACE
   SUBROUTINE SCOPY(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE SCOPY
END INTERFACE
INTERFACE
   SUBROUTINE DCOPY(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE DCOPY
END INTERFACE
INTERFACE
   SUBROUTINE CCOPY(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE CCOPY
END INTERFACE
INTERFACE
   SUBROUTINE ZCOPY(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE ZCOPY
END INTERFACE
INTERFACE
   SUBROUTINE SAXPY(N,ALPHA,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE SAXPY
END INTERFACE
INTERFACE
   SUBROUTINE DAXPY(N,ALPHA,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE DAXPY
END INTERFACE
INTERFACE
   SUBROUTINE CAXPY(N,ALPHA,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE CAXPY
END INTERFACE
INTERFACE
   SUBROUTINE ZAXPY(N,ALPHA,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE ZAXPY
END INTERFACE
INTERFACE
   FUNCTION SDOT(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(IN)                     :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(SP)                                 :: SDOT
   END FUNCTION SDOT
END INTERFACE
INTERFACE
   FUNCTION DDOT(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(IN)                     :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(DP)                                 :: DDOT
   END FUNCTION DDOT
END INTERFACE
INTERFACE
   FUNCTION CDOTU(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(IN)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    COMPLEX(SP)                              :: CDOTU
   END FUNCTION CDOTU
END INTERFACE
INTERFACE
   FUNCTION ZDOTU(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(IN)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    COMPLEX(DP)                              :: ZDOTU
   END FUNCTION ZDOTU
END INTERFACE
INTERFACE
   FUNCTION CDOTC(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(IN)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    COMPLEX(SP)                              :: CDOTC
   END FUNCTION CDOTC
END INTERFACE
INTERFACE
   FUNCTION ZDOTC(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(IN)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    COMPLEX(DP)                              :: ZDOTC
   END FUNCTION ZDOTC
END INTERFACE
INTERFACE
   FUNCTION SNRM2(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP)                                 :: SNRM2
   END FUNCTION SNRM2
END INTERFACE
INTERFACE
   FUNCTION DNRM2(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP)                                 :: DNRM2
   END FUNCTION DNRM2
END INTERFACE
INTERFACE
   FUNCTION SCNRM2(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP)                                 :: SCNRM2
   END FUNCTION SCNRM2
END INTERFACE
INTERFACE
   FUNCTION DZNRM2(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP)                                 :: DZNRM2
   END FUNCTION DZNRM2
END INTERFACE
INTERFACE
   FUNCTION SASUM(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP)                                 :: SASUM
   END FUNCTION SASUM
END INTERFACE
INTERFACE
   FUNCTION SCASUM(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP)                                 :: SCASUM
   END FUNCTION SCASUM
END INTERFACE
INTERFACE
   FUNCTION DASUM(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP)                                 :: DASUM
   END FUNCTION DASUM
END INTERFACE
INTERFACE
   FUNCTION DZASUM(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP)                                 :: DZASUM
   END FUNCTION DZASUM
END INTERFACE
INTERFACE
   FUNCTION ISAMAX(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    INTEGER                                  :: ISAMAX
   END FUNCTION ISAMAX
END INTERFACE
INTERFACE
   FUNCTION IDAMAX(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    INTEGER                                  :: IDAMAX
   END FUNCTION IDAMAX
END INTERFACE
INTERFACE
   FUNCTION ICAMAX(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    INTEGER                                  :: ICAMAX
   END FUNCTION ICAMAX
END INTERFACE
INTERFACE
   FUNCTION IZAMAX(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    INTEGER                                  :: IZAMAX
   END FUNCTION IZAMAX
END INTERFACE
!
! Level 1 BLAS
!
INTERFACE
   SUBROUTINE SGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: TRANS
    INTEGER, INTENT(IN)                      :: M, N
    REAL(SP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    REAL(SP), INTENT(IN)                     :: A(LDA,*), X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(IN)                     :: BETA
    REAL(SP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE SGEMV
END INTERFACE
INTERFACE
   SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: TRANS
    INTEGER, INTENT(IN)                      :: M, N
    REAL(DP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    REAL(DP), INTENT(IN)                     :: A(LDA,*), X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(IN)                     :: BETA
    REAL(DP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE DGEMV
END INTERFACE
INTERFACE
   SUBROUTINE CGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: TRANS
    INTEGER, INTENT(IN)                      :: M, N
    COMPLEX(SP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(SP), INTENT(IN)                  :: A(LDA,*), X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(IN)                  :: BETA
    COMPLEX(SP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE CGEMV
END INTERFACE
INTERFACE
   SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: TRANS
    INTEGER, INTENT(IN)                      :: M, N
    COMPLEX(DP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(DP), INTENT(IN)                  :: A(LDA,*), X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(IN)                  :: BETA
    COMPLEX(DP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE ZGEMV
END INTERFACE
INTERFACE
   SUBROUTINE SGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: TRANS
    INTEGER, INTENT(IN)                      :: M, N, KL, KU
    REAL(SP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    REAL(SP), INTENT(IN)                     :: A(LDA,*), X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(IN)                     :: BETA
    REAL(SP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE SGBMV
END INTERFACE
INTERFACE
   SUBROUTINE DGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: TRANS
    INTEGER, INTENT(IN)                      :: M, N, KL, KU
    REAL(DP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    REAL(DP), INTENT(IN)                     :: A(LDA,*), X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(IN)                     :: BETA
    REAL(DP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE DGBMV
END INTERFACE
INTERFACE
   SUBROUTINE CGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: TRANS
    INTEGER, INTENT(IN)                      :: M, N, KL, KU
    COMPLEX(SP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(SP), INTENT(IN)                  :: A(LDA,*), X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(IN)                  :: BETA
    COMPLEX(SP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE CGBMV
END INTERFACE
INTERFACE
   SUBROUTINE ZGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: TRANS
    INTEGER, INTENT(IN)                      :: M, N, KL, KU
    COMPLEX(DP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(DP), INTENT(IN)                  :: A(LDA,*), X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(IN)                  :: BETA
    COMPLEX(DP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE ZGBMV
END INTERFACE
INTERFACE
   SUBROUTINE CHEMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(SP), INTENT(IN)                  :: A(LDA,*), X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(IN)                  :: BETA
    COMPLEX(SP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE CHEMV
END INTERFACE
INTERFACE
   SUBROUTINE ZHEMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(DP), INTENT(IN)                  :: A(LDA,*), X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(IN)                  :: BETA
    COMPLEX(DP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE ZHEMV
END INTERFACE
INTERFACE
   SUBROUTINE CHBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N, K
    COMPLEX(SP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(SP), INTENT(IN)                  :: A(LDA,*), X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(IN)                  :: BETA
    COMPLEX(SP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE CHBMV
END INTERFACE
INTERFACE
   SUBROUTINE ZHBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N, K
    COMPLEX(DP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(DP), INTENT(IN)                  :: A(LDA,*), X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(IN)                  :: BETA
    COMPLEX(DP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE ZHBMV
END INTERFACE
INTERFACE
   SUBROUTINE CHPMV(UPLO,N,ALPHA,A,X,INCX,BETA,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: ALPHA, A(*), X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(IN)                  :: BETA
    COMPLEX(SP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE CHPMV
END INTERFACE
INTERFACE
   SUBROUTINE ZHPMV(UPLO,N,ALPHA,A,X,INCX,BETA,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: ALPHA, A(*), X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(IN)                  :: BETA
    COMPLEX(DP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE ZHPMV
END INTERFACE
INTERFACE
   SUBROUTINE SSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    REAL(SP), INTENT(IN)                     :: A(LDA,*), X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(IN)                     :: BETA
    REAL(SP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE SSYMV
END INTERFACE
INTERFACE
   SUBROUTINE DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    REAL(DP), INTENT(IN)                     :: A(LDA,*), X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(IN)                     :: BETA
    REAL(DP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE DSYMV
END INTERFACE
INTERFACE
   SUBROUTINE SSBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N, K
    REAL(SP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    REAL(SP), INTENT(IN)                     :: A(LDA,*), X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(IN)                     :: BETA
    REAL(SP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE SSBMV
END INTERFACE
INTERFACE
   SUBROUTINE DSBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N, K
    REAL(DP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    REAL(DP), INTENT(IN)                     :: A(LDA,*), X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(IN)                     :: BETA
    REAL(DP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE DSBMV
END INTERFACE
INTERFACE
   SUBROUTINE SSPMV(UPLO,N,ALPHA,A,X,INCX,BETA,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: ALPHA, A(*), X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(IN)                     :: BETA
    REAL(SP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE SSPMV
END INTERFACE
INTERFACE
   SUBROUTINE DSPMV(UPLO,N,ALPHA,A,X,INCX,BETA,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: ALPHA, A(*), X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(IN)                     :: BETA
    REAL(DP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE DSPMV
END INTERFACE
INTERFACE
   SUBROUTINE STRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, LDA
    REAL(SP), INTENT(IN)                     :: A(LDA,*)
    REAL(SP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE STRMV
END INTERFACE
INTERFACE
   SUBROUTINE DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, LDA
    REAL(DP), INTENT(IN)                     :: A(LDA,*)
    REAL(DP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE DTRMV
END INTERFACE
INTERFACE
   SUBROUTINE CTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, LDA
    COMPLEX(SP), INTENT(IN)                  :: A(LDA,*)
    COMPLEX(SP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE CTRMV
END INTERFACE
INTERFACE
   SUBROUTINE ZTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, LDA
    COMPLEX(DP), INTENT(IN)                  :: A(LDA,*)
    COMPLEX(DP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE ZTRMV
END INTERFACE
INTERFACE
   SUBROUTINE STBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, K, LDA
    REAL(SP), INTENT(IN)                     :: A(LDA,*)
    REAL(SP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE STBMV
END INTERFACE
INTERFACE
   SUBROUTINE DTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, K, LDA
    REAL(DP), INTENT(IN)                     :: A(LDA,*)
    REAL(DP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE DTBMV
END INTERFACE
INTERFACE
   SUBROUTINE CTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, K, LDA
    COMPLEX(SP), INTENT(IN)                  :: A(LDA,*)
    COMPLEX(SP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE CTBMV
END INTERFACE
INTERFACE
   SUBROUTINE ZTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, K, LDA
    COMPLEX(DP), INTENT(IN)                  :: A(LDA,*)
    COMPLEX(DP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE ZTBMV
END INTERFACE
INTERFACE
   SUBROUTINE STPMV(UPLO,TRANS,DIAG,N,A,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: A(*)
    REAL(SP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE STPMV
END INTERFACE
INTERFACE
   SUBROUTINE DTPMV(UPLO,TRANS,DIAG,N,A,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: A(*)
    REAL(DP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE DTPMV
END INTERFACE
INTERFACE
   SUBROUTINE CTPMV(UPLO,TRANS,DIAG,N,A,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: A(*)
    COMPLEX(SP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE CTPMV
END INTERFACE
INTERFACE
   SUBROUTINE ZTPMV(UPLO,TRANS,DIAG,N,A,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: A(*)
    COMPLEX(DP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE ZTPMV
END INTERFACE
INTERFACE
   SUBROUTINE STRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, LDA
    REAL(SP), INTENT(IN)                     :: A(LDA,*)
    REAL(SP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE STRSV
END INTERFACE
INTERFACE
   SUBROUTINE DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, LDA
    REAL(DP), INTENT(IN)                     :: A(LDA,*)
    REAL(DP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE DTRSV
END INTERFACE
INTERFACE
   SUBROUTINE CTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, LDA
    COMPLEX(SP), INTENT(IN)                  :: A(LDA,*)
    COMPLEX(SP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE CTRSV
END INTERFACE
INTERFACE
   SUBROUTINE ZTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, LDA
    COMPLEX(DP), INTENT(IN)                  :: A(LDA,*)
    COMPLEX(DP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE ZTRSV
END INTERFACE
INTERFACE
   SUBROUTINE STBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, K, LDA
    REAL(SP), INTENT(IN)                     :: A(LDA,*)
    REAL(SP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE STBSV
END INTERFACE
INTERFACE
   SUBROUTINE DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, K, LDA
    REAL(DP), INTENT(IN)                     :: A(LDA,*)
    REAL(DP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE DTBSV
END INTERFACE
INTERFACE
   SUBROUTINE CTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, K, LDA
    COMPLEX(SP), INTENT(IN)                  :: A(LDA,*)
    COMPLEX(SP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE CTBSV
END INTERFACE
INTERFACE
   SUBROUTINE ZTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, K, LDA
    COMPLEX(DP), INTENT(IN)                  :: A(LDA,*)
    COMPLEX(DP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE ZTBSV
END INTERFACE
INTERFACE
   SUBROUTINE STPSV(UPLO,TRANS,DIAG,N,A,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: A(*)
    REAL(SP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE STPSV
END INTERFACE
INTERFACE
   SUBROUTINE DTPSV(UPLO,TRANS,DIAG,N,A,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: A(*)
    REAL(DP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE DTPSV
END INTERFACE
INTERFACE
   SUBROUTINE CTPSV(UPLO,TRANS,DIAG,N,A,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: A(*)
    COMPLEX(SP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE CTPSV
END INTERFACE
INTERFACE
   SUBROUTINE ZTPSV(UPLO,TRANS,DIAG,N,A,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: A(*)
    COMPLEX(DP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE ZTPSV
END INTERFACE
INTERFACE
   SUBROUTINE SGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: M, N
    REAL(SP), INTENT(IN)                     :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(IN)                     :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY, LDA
    REAL(SP), INTENT(INOUT)                  :: A(LDA,*)
   END SUBROUTINE SGER
END INTERFACE
INTERFACE
   SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: M, N
    REAL(DP), INTENT(IN)                     :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(IN)                     :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY, LDA
    REAL(DP), INTENT(INOUT)                  :: A(LDA,*)
   END SUBROUTINE DGER
END INTERFACE
INTERFACE
   SUBROUTINE CGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: M, N
    COMPLEX(SP), INTENT(IN)                  :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(IN)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY, LDA
    COMPLEX(SP), INTENT(INOUT)               :: A(LDA,*)
   END SUBROUTINE CGERU
END INTERFACE
INTERFACE
   SUBROUTINE ZGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: M, N
    COMPLEX(DP), INTENT(IN)                  :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(IN)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY, LDA
    COMPLEX(DP), INTENT(INOUT)               :: A(LDA,*)
   END SUBROUTINE ZGERU
END INTERFACE
INTERFACE
   SUBROUTINE CGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: M, N
    COMPLEX(SP), INTENT(IN)                  :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(IN)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY, LDA
    COMPLEX(SP), INTENT(INOUT)               :: A(LDA,*)
   END SUBROUTINE CGERC
END INTERFACE
INTERFACE
   SUBROUTINE ZGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: M, N
    COMPLEX(DP), INTENT(IN)                  :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(IN)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY, LDA
    COMPLEX(DP), INTENT(INOUT)               :: A(LDA,*)
   END SUBROUTINE ZGERC
END INTERFACE
INTERFACE
   SUBROUTINE CHER(UPLO,N,ALPHA,X,INCX,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: ALPHA
    COMPLEX(SP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX, LDA
    COMPLEX(SP), INTENT(INOUT)               :: A(LDA,*)
   END SUBROUTINE CHER
END INTERFACE
INTERFACE
   SUBROUTINE ZHER(UPLO,N,ALPHA,X,INCX,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: ALPHA
    COMPLEX(DP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX, LDA
    COMPLEX(DP), INTENT(INOUT)               :: A(LDA,*)
   END SUBROUTINE ZHER
END INTERFACE
INTERFACE
   SUBROUTINE CHPR(UPLO,N,ALPHA,X,INCX,AP)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: ALPHA
    COMPLEX(SP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(INOUT)               :: AP(*)
   END SUBROUTINE CHPR
END INTERFACE
INTERFACE
   SUBROUTINE ZHPR(UPLO,N,ALPHA,X,INCX,AP)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: ALPHA
    COMPLEX(DP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(INOUT)               :: AP(*)
   END SUBROUTINE ZHPR
END INTERFACE
INTERFACE
   SUBROUTINE CHER2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(IN)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY, LDA
    COMPLEX(SP), INTENT(INOUT)               :: A(LDA,*)
   END SUBROUTINE CHER2
END INTERFACE
INTERFACE
   SUBROUTINE ZHER2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(IN)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY, LDA
    COMPLEX(DP), INTENT(INOUT)               :: A(LDA,*)
   END SUBROUTINE ZHER2
END INTERFACE
INTERFACE
   SUBROUTINE CHPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(IN)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    COMPLEX(SP), INTENT(INOUT)               :: AP(*)
   END SUBROUTINE CHPR2
END INTERFACE
INTERFACE
   SUBROUTINE ZHPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(IN)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    COMPLEX(DP), INTENT(INOUT)               :: AP(*)
   END SUBROUTINE ZHPR2
END INTERFACE
INTERFACE
   SUBROUTINE SSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX, LDA
    REAL(SP), INTENT(INOUT)                  :: A(LDA,*)
   END SUBROUTINE SSYR
END INTERFACE
INTERFACE
   SUBROUTINE DSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX, LDA
    REAL(DP), INTENT(INOUT)                  :: A(LDA,*)
   END SUBROUTINE DSYR
END INTERFACE
INTERFACE
   SUBROUTINE SSPR(UPLO,N,ALPHA,X,INCX,AP)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(INOUT)                  :: AP(*)
   END SUBROUTINE SSPR
END INTERFACE
INTERFACE
   SUBROUTINE DSPR(UPLO,N,ALPHA,X,INCX,AP)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(INOUT)                  :: AP(*)
   END SUBROUTINE DSPR
END INTERFACE
INTERFACE
   SUBROUTINE SSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(IN)                     :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY, LDA
    REAL(SP), INTENT(INOUT)                  :: A(LDA,*)
   END SUBROUTINE SSYR2
END INTERFACE
INTERFACE
   SUBROUTINE DSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(IN)                     :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY, LDA
    REAL(DP), INTENT(INOUT)                  :: A(LDA,*)
   END SUBROUTINE DSYR2
END INTERFACE
INTERFACE
   SUBROUTINE SSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(IN)                     :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(SP), INTENT(INOUT)                  :: AP(*)
   END SUBROUTINE SSPR2
END INTERFACE
INTERFACE
   SUBROUTINE DSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(IN)                     :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(DP), INTENT(INOUT)                  :: AP(*)
   END SUBROUTINE DSPR2
END INTERFACE
!
! Level 3 BLAS
!
INTERFACE
   SUBROUTINE SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: TRANSA, TRANSB
    INTEGER, INTENT(IN)                      :: M, N, K
    REAL(SP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    REAL(SP), INTENT(IN)                     :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    REAL(SP), INTENT(IN)                     :: B(LDB,*), BETA
    INTEGER, INTENT(IN)                      :: LDC
    REAL(SP), INTENT(INOUT)                  :: C(LDC,*)
   END SUBROUTINE SGEMM
END INTERFACE
INTERFACE
   SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: TRANSA, TRANSB
    INTEGER, INTENT(IN)                      :: M, N, K
    REAL(DP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    REAL(DP), INTENT(IN)                     :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    REAL(DP), INTENT(IN)                     :: B(LDB,*), BETA
    INTEGER, INTENT(IN)                      :: LDC
    REAL(DP), INTENT(INOUT)                  :: C(LDC,*)
   END SUBROUTINE DGEMM
END INTERFACE
INTERFACE
   SUBROUTINE CGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: TRANSA, TRANSB
    INTEGER, INTENT(IN)                      :: M, N, K
    COMPLEX(SP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(SP), INTENT(IN)                  :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    COMPLEX(SP), INTENT(IN)                  :: B(LDB,*), BETA
    INTEGER, INTENT(IN)                      :: LDC
    COMPLEX(SP), INTENT(INOUT)               :: C(LDC,*)
   END SUBROUTINE CGEMM
END INTERFACE
INTERFACE
   SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: TRANSA, TRANSB
    INTEGER, INTENT(IN)                      :: M, N, K
    COMPLEX(DP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(DP), INTENT(IN)                  :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    COMPLEX(DP), INTENT(IN)                  :: B(LDB,*), BETA
    INTEGER, INTENT(IN)                      :: LDC
    COMPLEX(DP), INTENT(INOUT)               :: C(LDC,*)
   END SUBROUTINE ZGEMM
END INTERFACE
INTERFACE
   SUBROUTINE SSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: SIDE, UPLO
    INTEGER, INTENT(IN)                      :: M, N
    REAL(SP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    REAL(SP), INTENT(IN)                     :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    REAL(SP), INTENT(IN)                     :: B(LDB,*), BETA
    INTEGER, INTENT(IN)                      :: LDC
    REAL(SP), INTENT(INOUT)                  :: C(LDC,*)
   END SUBROUTINE SSYMM
END INTERFACE
INTERFACE
   SUBROUTINE DSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: SIDE, UPLO
    INTEGER, INTENT(IN)                      :: M, N
    REAL(DP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    REAL(DP), INTENT(IN)                     :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    REAL(DP), INTENT(IN)                     :: B(LDB,*), BETA
    INTEGER, INTENT(IN)                      :: LDC
    REAL(DP), INTENT(INOUT)                  :: C(LDC,*)
   END SUBROUTINE DSYMM
END INTERFACE
INTERFACE
   SUBROUTINE CSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: SIDE, UPLO
    INTEGER, INTENT(IN)                      :: M, N
    COMPLEX(SP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(SP), INTENT(IN)                  :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    COMPLEX(SP), INTENT(IN)                  :: B(LDB,*), BETA
    INTEGER, INTENT(IN)                      :: LDC
    COMPLEX(SP), INTENT(INOUT)               :: C(LDC,*)
   END SUBROUTINE CSYMM
END INTERFACE
INTERFACE
   SUBROUTINE ZSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: SIDE, UPLO
    INTEGER, INTENT(IN)                      :: M, N
    COMPLEX(DP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(DP), INTENT(IN)                  :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    COMPLEX(DP), INTENT(IN)                  :: B(LDB,*), BETA
    INTEGER, INTENT(IN)                      :: LDC
    COMPLEX(DP), INTENT(INOUT)               :: C(LDC,*)
   END SUBROUTINE ZSYMM
END INTERFACE
INTERFACE
   SUBROUTINE CHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: SIDE, UPLO
    INTEGER, INTENT(IN)                      :: M, N
    COMPLEX(SP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(SP), INTENT(IN)                  :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    COMPLEX(SP), INTENT(IN)                  :: B(LDB,*), BETA
    INTEGER, INTENT(IN)                      :: LDC
    COMPLEX(SP), INTENT(INOUT)               :: C(LDC,*)
   END SUBROUTINE CHEMM
END INTERFACE
INTERFACE
   SUBROUTINE ZHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: SIDE, UPLO
    INTEGER, INTENT(IN)                      :: M, N
    COMPLEX(DP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(DP), INTENT(IN)                  :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    COMPLEX(DP), INTENT(IN)                  :: B(LDB,*), BETA
    INTEGER, INTENT(IN)                      :: LDC
    COMPLEX(DP), INTENT(INOUT)               :: C(LDC,*)
   END SUBROUTINE ZHEMM
END INTERFACE
INTERFACE
   SUBROUTINE SSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS
    INTEGER, INTENT(IN)                      :: N, K
    REAL(SP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    REAL(SP), INTENT(IN)                     :: A(LDA,*), BETA
    INTEGER, INTENT(IN)                      :: LDC
    REAL(SP), INTENT(INOUT)                  :: C(LDC,*)
   END SUBROUTINE SSYRK
END INTERFACE
INTERFACE
   SUBROUTINE DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS
    INTEGER, INTENT(IN)                      :: N, K
    REAL(DP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    REAL(DP), INTENT(IN)                     :: A(LDA,*), BETA
    INTEGER, INTENT(IN)                      :: LDC
    REAL(DP), INTENT(INOUT)                  :: C(LDC,*)
   END SUBROUTINE DSYRK
END INTERFACE
INTERFACE
   SUBROUTINE CSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS
    INTEGER, INTENT(IN)                      :: N, K
    COMPLEX(SP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(SP), INTENT(IN)                  :: A(LDA,*), BETA
    INTEGER, INTENT(IN)                      :: LDC
    COMPLEX(SP), INTENT(INOUT)               :: C(LDC,*)
   END SUBROUTINE CSYRK
END INTERFACE
INTERFACE
   SUBROUTINE ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS
    INTEGER, INTENT(IN)                      :: N, K
    COMPLEX(DP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(DP), INTENT(IN)                  :: A(LDA,*), BETA
    INTEGER, INTENT(IN)                      :: LDC
    COMPLEX(DP), INTENT(INOUT)               :: C(LDC,*)
   END SUBROUTINE ZSYRK
END INTERFACE
INTERFACE
   SUBROUTINE CHERK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS
    INTEGER, INTENT(IN)                      :: N, K
    REAL(SP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(SP), INTENT(IN)                  :: A(LDA,*)
    REAL(SP), INTENT(IN)                     :: BETA
    INTEGER, INTENT(IN)                      :: LDC
    COMPLEX(SP), INTENT(INOUT)               :: C(LDC,*)
   END SUBROUTINE CHERK
END INTERFACE
INTERFACE
   SUBROUTINE ZHERK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS
    INTEGER, INTENT(IN)                      :: N, K
    REAL(DP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(DP), INTENT(IN)                  :: A(LDA,*)
    REAL(DP), INTENT(IN)                     :: BETA
    INTEGER, INTENT(IN)                      :: LDC
    COMPLEX(DP), INTENT(INOUT)               :: C(LDC,*)
   END SUBROUTINE ZHERK
END INTERFACE
INTERFACE
   SUBROUTINE SSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS
    INTEGER, INTENT(IN)                      :: N, K
    REAL(SP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    REAL(SP), INTENT(IN)                     :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    REAL(SP), INTENT(IN)                     :: B(LDB,*), BETA
    INTEGER, INTENT(IN)                      :: LDC
    REAL(SP), INTENT(INOUT)                  :: C(LDC,*)
   END SUBROUTINE SSYR2K
END INTERFACE
INTERFACE
   SUBROUTINE DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS
    INTEGER, INTENT(IN)                      :: N, K
    REAL(DP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    REAL(DP), INTENT(IN)                     :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    REAL(DP), INTENT(IN)                     :: B(LDB,*), BETA
    INTEGER, INTENT(IN)                      :: LDC
    REAL(DP), INTENT(INOUT)                  :: C(LDC,*)
   END SUBROUTINE DSYR2K
END INTERFACE
INTERFACE
   SUBROUTINE CSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS
    INTEGER, INTENT(IN)                      :: N, K
    COMPLEX(SP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(SP), INTENT(IN)                  :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    COMPLEX(SP), INTENT(IN)                  :: B(LDB,*), BETA
    INTEGER, INTENT(IN)                      :: LDC
    COMPLEX(SP), INTENT(INOUT)               :: C(LDC,*)
   END SUBROUTINE CSYR2K
END INTERFACE
INTERFACE
   SUBROUTINE ZSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS
    INTEGER, INTENT(IN)                      :: N, K
    COMPLEX(DP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(DP), INTENT(IN)                  :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    COMPLEX(DP), INTENT(IN)                  :: B(LDB,*), BETA
    INTEGER, INTENT(IN)                      :: LDC
    COMPLEX(DP), INTENT(INOUT)               :: C(LDC,*)
   END SUBROUTINE ZSYR2K
END INTERFACE
INTERFACE
   SUBROUTINE CHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS
    INTEGER, INTENT(IN)                      :: N, K
    COMPLEX(SP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(SP), INTENT(IN)                  :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    COMPLEX(SP), INTENT(IN)                  :: B(LDB,*)
    REAL(SP), INTENT(IN)                     :: BETA
    INTEGER, INTENT(IN)                      :: LDC
    COMPLEX(SP), INTENT(INOUT)               :: C(LDC,*)
   END SUBROUTINE CHER2K
END INTERFACE
INTERFACE
   SUBROUTINE ZHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS
    INTEGER, INTENT(IN)                      :: N, K
    COMPLEX(DP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(DP), INTENT(IN)                  :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    COMPLEX(DP), INTENT(IN)                  :: B(LDB,*)
    REAL(DP), INTENT(IN)                     :: BETA
    INTEGER, INTENT(IN)                      :: LDC
    COMPLEX(DP), INTENT(INOUT)               :: C(LDC,*)
   END SUBROUTINE ZHER2K
END INTERFACE
INTERFACE
   SUBROUTINE STRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: SIDE, UPLO, TRANSA, DIAG
    INTEGER, INTENT(IN)                      :: M, N
    REAL(SP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    REAL(SP), INTENT(IN)                     :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    REAL(SP), INTENT(INOUT)                  :: B(LDB,*)
   END SUBROUTINE STRMM
END INTERFACE
INTERFACE
   SUBROUTINE DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: SIDE, UPLO, TRANSA, DIAG
    INTEGER, INTENT(IN)                      :: M, N
    REAL(DP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    REAL(DP), INTENT(IN)                     :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    REAL(DP), INTENT(INOUT)                  :: B(LDB,*)
   END SUBROUTINE DTRMM
END INTERFACE
INTERFACE
   SUBROUTINE CTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: SIDE, UPLO, TRANSA, DIAG
    INTEGER, INTENT(IN)                      :: M, N
    COMPLEX(SP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(SP), INTENT(IN)                  :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    COMPLEX(SP), INTENT(INOUT)               :: B(LDB,*)
   END SUBROUTINE CTRMM
END INTERFACE
INTERFACE
   SUBROUTINE ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: SIDE, UPLO, TRANSA, DIAG
    INTEGER, INTENT(IN)                      :: M, N
    COMPLEX(DP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(DP), INTENT(IN)                  :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    COMPLEX(DP), INTENT(INOUT)               :: B(LDB,*)
   END SUBROUTINE ZTRMM
END INTERFACE
INTERFACE
   SUBROUTINE STRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: SIDE, UPLO, TRANSA, DIAG
    INTEGER, INTENT(IN)                      :: M, N
    REAL(SP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    REAL(SP), INTENT(IN)                     :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    REAL(SP), INTENT(INOUT)                  :: B(LDB,*)
   END SUBROUTINE STRSM
END INTERFACE
INTERFACE
   SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: SIDE, UPLO, TRANSA, DIAG
    INTEGER, INTENT(IN)                      :: M, N
    REAL(DP), INTENT(IN)                     :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    REAL(DP), INTENT(IN)                     :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    REAL(DP), INTENT(INOUT)                  :: B(LDB,*)
   END SUBROUTINE DTRSM
END INTERFACE
INTERFACE
   SUBROUTINE CTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: SIDE, UPLO, TRANSA, DIAG
    INTEGER, INTENT(IN)                      :: M, N
    COMPLEX(SP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(SP), INTENT(IN)                  :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    COMPLEX(SP), INTENT(INOUT)               :: B(LDB,*)
   END SUBROUTINE CTRSM
END INTERFACE
INTERFACE
   SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: SIDE, UPLO, TRANSA, DIAG
    INTEGER, INTENT(IN)                      :: M, N
    COMPLEX(DP), INTENT(IN)                  :: ALPHA
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(DP), INTENT(IN)                  :: A(LDA,*)
    INTEGER, INTENT(IN)                      :: LDB
    COMPLEX(DP), INTENT(INOUT)               :: B(LDB,*)
   END SUBROUTINE ZTRSM
END INTERFACE
INTERFACE
   FUNCTION SDSDOT(N,ALPHA,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(IN)                     :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(SP)                                 :: SDSDOT
   END FUNCTION SDSDOT
END INTERFACE
INTERFACE
   FUNCTION DSDOT(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(IN)                     :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(DP)                                 :: DSDOT
   END FUNCTION DSDOT
END INTERFACE
INTERFACE
   FUNCTION DCABS1(Z)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    COMPLEX(DP), INTENT(IN)                  :: Z
    REAL(DP)                                 :: DCABS1
   END FUNCTION DCABS1
END INTERFACE
END MODULE f77_blas_netlib
