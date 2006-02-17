MODULE f77_blas_generic
!
! Provides generic interfaces for BLAS as distributed
! (http://www.netlib.org/blas/) Jan 2006.
! Similar to the module F77_LAPACK as provided by LAPACK95.
! All BLAS calls in a code can be replaced by these generic
! ones, without modification of the arguments.
!
! the calls change the first letter (which indicate the kind in blas)
! into bl_
!
! i.e.
! CALL DDOT(N,X,INCX,Y,INCY)
! becomes
! CALL BL_DOT(N,X,INCX,Y,INCY)
!
! dsdot, sdsdot and dcabs1 have not been included as they are
! not generic (the first two functions are also not part of the blas
! provided as part of LAPACK V3.0), but interfaces can be found in
! the module f77_blas_netlib
!
! follows the same naming convention as the quick ref
! http://www.netlib.org/lapack/lug/node145.html
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

PUBLIC :: BL_ROTG, BL_ROTMG, BL_ROT, BL_ROTM, BL_SWAP, BL_SCAL, BL_COPY, &
          BL_AXPY, BL_DOT, BL_DOTU, BL_DOTC, BL_NRM2, BL_ASUM, BL_IAMAX

PUBLIC :: BL_GEMV, BL_GBMV, BL_HEMV, BL_HBMV, BL_HPMV, BL_SYMV, BL_SBMV, &
          BL_SPMV, BL_TRMV, BL_TBMV, BL_TPMV, BL_TRSV, BL_TBSV, BL_TPSV, &
          BL_GER, BL_GERU, BL_GERC, BL_HER, BL_HPR, BL_HER2, BL_HPR2,    &
          BL_SYR, BL_SPR, BL_SYR2, BL_SPR2

PUBLIC :: BL_GEMM, BL_SYMM, BL_HEMM, BL_SYRK, BL_HERK, BL_SYR2K,         &
          BL_HER2K, BL_TRMM, BL_TRSM

!
! Level 1 BLAS
!
INTERFACE BL_ROTG
   SUBROUTINE SROTG(A,B,C,S)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    REAL(SP), INTENT(INOUT)                  :: A, B
    REAL(SP), INTENT(OUT)                    :: C, S
   END SUBROUTINE SROTG
   SUBROUTINE DROTG(A,B,C,S)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    REAL(DP), INTENT(INOUT)                  :: A, B
    REAL(DP), INTENT(OUT)                    :: C, S
   END SUBROUTINE DROTG
   SUBROUTINE CROTG(A,B,C,S)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    COMPLEX(SP), INTENT(INOUT)               :: A, B
    REAL(SP), INTENT(OUT)                    :: C
    COMPLEX(SP), INTENT(OUT)                 :: S
   END SUBROUTINE CROTG
   SUBROUTINE ZROTG(A,B,C,S)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    COMPLEX(DP), INTENT(INOUT)               :: A, B
    REAL(DP), INTENT(OUT)                    :: C
    COMPLEX(DP), INTENT(OUT)                 :: S
   END SUBROUTINE ZROTG
END INTERFACE BL_ROTG
INTERFACE BL_ROTMG
   SUBROUTINE SROTMG(D1,D2,A,B,PARAM)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    REAL(SP), INTENT(INOUT)                  :: D1, D2, A
    REAL(SP), INTENT(IN)                     :: B
    REAL(SP), INTENT(OUT)                    :: PARAM(5)
   END SUBROUTINE SROTMG
   SUBROUTINE DROTMG(D1,D2,A,B,PARAM)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    REAL(DP), INTENT(INOUT)                  :: D1, D2, A
    REAL(DP), INTENT(IN)                     :: B
    REAL(DP), INTENT(OUT)                    :: PARAM(5)
   END SUBROUTINE DROTMG
END INTERFACE BL_ROTMG
INTERFACE BL_ROT
   SUBROUTINE SROT(N,X,INCX,Y,INCY,C,S)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(SP), INTENT(IN)                     :: C, S
   END SUBROUTINE SROT
   SUBROUTINE DROT(N,X,INCX,Y,INCY,C,S)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(DP), INTENT(IN)                     :: C, S
   END SUBROUTINE DROT
   SUBROUTINE CSROT(N,X,INCX,Y,INCY,C,S)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(SP), INTENT(IN)                     :: C, S
   END SUBROUTINE CSROT
   SUBROUTINE ZDROT(N,X,INCX,Y,INCY,C,S)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(DP), INTENT(IN)                     :: C, S
   END SUBROUTINE ZDROT
END INTERFACE BL_ROT
INTERFACE BL_ROTM
   SUBROUTINE SROTM(N,X,INCX,Y,INCY,PARAM)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(SP), INTENT(IN)                     :: PARAM(5)
   END SUBROUTINE SROTM
   SUBROUTINE DROTM(N,X,INCX,Y,INCY,PARAM)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(DP), INTENT(IN)                     :: PARAM(5)
   END SUBROUTINE DROTM
END INTERFACE BL_ROTM
INTERFACE BL_SWAP
   SUBROUTINE SSWAP(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE SSWAP
   SUBROUTINE DSWAP(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE DSWAP
   SUBROUTINE CSWAP(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE CSWAP
   SUBROUTINE ZSWAP(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE ZSWAP
END INTERFACE BL_SWAP
INTERFACE BL_SCAL
   SUBROUTINE SSCAL(N,ALPHA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: ALPHA
    REAL(SP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE SSCAL
   SUBROUTINE DSCAL(N,ALPHA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: ALPHA
    REAL(DP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE DSCAL
   SUBROUTINE CSCAL(N,ALPHA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: ALPHA
    COMPLEX(SP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE CSCAL
   SUBROUTINE ZSCAL(N,ALPHA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: ALPHA
    COMPLEX(DP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE ZSCAL
   SUBROUTINE CSSCAL(N,ALPHA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: ALPHA
    COMPLEX(SP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE CSSCAL
   SUBROUTINE ZDSCAL(N,ALPHA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: ALPHA
    COMPLEX(DP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE ZDSCAL
END INTERFACE BL_SCAL
INTERFACE BL_COPY
   SUBROUTINE SCOPY(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE SCOPY
   SUBROUTINE DCOPY(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE DCOPY
   SUBROUTINE CCOPY(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE CCOPY
   SUBROUTINE ZCOPY(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE ZCOPY
END INTERFACE BL_COPY
INTERFACE BL_AXPY
   SUBROUTINE SAXPY(N,ALPHA,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE SAXPY
   SUBROUTINE DAXPY(N,ALPHA,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE DAXPY
   SUBROUTINE CAXPY(N,ALPHA,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE CAXPY
   SUBROUTINE ZAXPY(N,ALPHA,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
   END SUBROUTINE ZAXPY
END INTERFACE BL_AXPY
INTERFACE BL_DOT
   FUNCTION SDOT(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(IN)                     :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(SP)                                 :: SDOT
   END FUNCTION SDOT
   FUNCTION DDOT(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(IN)                     :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(DP)                                 :: DDOT
   END FUNCTION DDOT
END INTERFACE BL_DOT
INTERFACE BL_DOTU
   FUNCTION CDOTU(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(IN)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    COMPLEX(SP)                              :: CDOTU
   END FUNCTION CDOTU
   FUNCTION ZDOTU(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(IN)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    COMPLEX(DP)                              :: ZDOTU
   END FUNCTION ZDOTU
END INTERFACE BL_DOTU
INTERFACE BL_DOTC
   FUNCTION CDOTC(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(IN)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    COMPLEX(SP)                              :: CDOTC
   END FUNCTION CDOTC
   FUNCTION ZDOTC(N,X,INCX,Y,INCY)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(IN)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY
    COMPLEX(DP)                              :: ZDOTC
   END FUNCTION ZDOTC
END INTERFACE BL_DOTC
INTERFACE BL_NRM2
   FUNCTION SNRM2(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP)                                 :: SNRM2
   END FUNCTION SNRM2
   FUNCTION DNRM2(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP)                                 :: DNRM2
   END FUNCTION DNRM2
   FUNCTION SCNRM2(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP)                                 :: SCNRM2
   END FUNCTION SCNRM2
   FUNCTION DZNRM2(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP)                                 :: DZNRM2
   END FUNCTION DZNRM2
END INTERFACE BL_NRM2
INTERFACE BL_ASUM
   FUNCTION SASUM(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP)                                 :: SASUM
   END FUNCTION SASUM
   FUNCTION SCASUM(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP)                                 :: SCASUM
   END FUNCTION SCASUM
   FUNCTION DASUM(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP)                                 :: DASUM
   END FUNCTION DASUM
   FUNCTION DZASUM(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP)                                 :: DZASUM
   END FUNCTION DZASUM
END INTERFACE BL_ASUM
INTERFACE BL_IAMAX
   FUNCTION ISAMAX(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    INTEGER                                  :: ISAMAX
   END FUNCTION ISAMAX
   FUNCTION IDAMAX(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    INTEGER                                  :: IDAMAX
   END FUNCTION IDAMAX
   FUNCTION ICAMAX(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    INTEGER                                  :: ICAMAX
   END FUNCTION ICAMAX
   FUNCTION IZAMAX(N,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    INTEGER                                  :: IZAMAX
   END FUNCTION IZAMAX
END INTERFACE BL_IAMAX
!
! Level 1 BLAS
!
INTERFACE BL_GEMV
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
END INTERFACE BL_GEMV
INTERFACE BL_GBMV
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
END INTERFACE BL_GBMV
INTERFACE BL_HEMV
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
END INTERFACE BL_HEMV
INTERFACE BL_HBMV
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
END INTERFACE BL_HBMV
INTERFACE BL_HPMV
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
END INTERFACE BL_HPMV
INTERFACE BL_SYMV
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
END INTERFACE BL_SYMV
INTERFACE BL_SBMV
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
END INTERFACE BL_SBMV
INTERFACE BL_SPMV
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
END INTERFACE BL_SPMV
INTERFACE BL_TRMV
   SUBROUTINE STRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, LDA
    REAL(SP), INTENT(IN)                     :: A(LDA,*)
    REAL(SP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE STRMV
   SUBROUTINE DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, LDA
    REAL(DP), INTENT(IN)                     :: A(LDA,*)
    REAL(DP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE DTRMV
   SUBROUTINE CTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, LDA
    COMPLEX(SP), INTENT(IN)                  :: A(LDA,*)
    COMPLEX(SP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE CTRMV
   SUBROUTINE ZTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, LDA
    COMPLEX(DP), INTENT(IN)                  :: A(LDA,*)
    COMPLEX(DP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE ZTRMV
END INTERFACE BL_TRMV
INTERFACE BL_TBMV
   SUBROUTINE STBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, K, LDA
    REAL(SP), INTENT(IN)                     :: A(LDA,*)
    REAL(SP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE STBMV
   SUBROUTINE DTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, K, LDA
    REAL(DP), INTENT(IN)                     :: A(LDA,*)
    REAL(DP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE DTBMV
   SUBROUTINE CTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, K, LDA
    COMPLEX(SP), INTENT(IN)                  :: A(LDA,*)
    COMPLEX(SP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE CTBMV
   SUBROUTINE ZTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, K, LDA
    COMPLEX(DP), INTENT(IN)                  :: A(LDA,*)
    COMPLEX(DP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE ZTBMV
END INTERFACE BL_TBMV
INTERFACE BL_TPMV
   SUBROUTINE STPMV(UPLO,TRANS,DIAG,N,A,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: A(*)
    REAL(SP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE STPMV
   SUBROUTINE DTPMV(UPLO,TRANS,DIAG,N,A,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: A(*)
    REAL(DP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE DTPMV
   SUBROUTINE CTPMV(UPLO,TRANS,DIAG,N,A,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: A(*)
    COMPLEX(SP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE CTPMV
   SUBROUTINE ZTPMV(UPLO,TRANS,DIAG,N,A,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: A(*)
    COMPLEX(DP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE ZTPMV
END INTERFACE BL_TPMV
INTERFACE BL_TRSV
   SUBROUTINE STRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, LDA
    REAL(SP), INTENT(IN)                     :: A(LDA,*)
    REAL(SP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE STRSV
   SUBROUTINE DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, LDA
    REAL(DP), INTENT(IN)                     :: A(LDA,*)
    REAL(DP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE DTRSV
   SUBROUTINE CTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, LDA
    COMPLEX(SP), INTENT(IN)                  :: A(LDA,*)
    COMPLEX(SP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE CTRSV
   SUBROUTINE ZTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, LDA
    COMPLEX(DP), INTENT(IN)                  :: A(LDA,*)
    COMPLEX(DP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE ZTRSV
END INTERFACE BL_TRSV
INTERFACE BL_TBSV
   SUBROUTINE STBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, K, LDA
    REAL(SP), INTENT(IN)                     :: A(LDA,*)
    REAL(SP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE STBSV
   SUBROUTINE DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, K, LDA
    REAL(DP), INTENT(IN)                     :: A(LDA,*)
    REAL(DP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE DTBSV
   SUBROUTINE CTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, K, LDA
    COMPLEX(SP), INTENT(IN)                  :: A(LDA,*)
    COMPLEX(SP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE CTBSV
   SUBROUTINE ZTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N, K, LDA
    COMPLEX(DP), INTENT(IN)                  :: A(LDA,*)
    COMPLEX(DP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE ZTBSV
END INTERFACE BL_TBSV
INTERFACE BL_TPSV
   SUBROUTINE STPSV(UPLO,TRANS,DIAG,N,A,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: A(*)
    REAL(SP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE STPSV
   SUBROUTINE DTPSV(UPLO,TRANS,DIAG,N,A,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: A(*)
    REAL(DP), INTENT(INOUT)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE DTPSV
   SUBROUTINE CTPSV(UPLO,TRANS,DIAG,N,A,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(SP), INTENT(IN)                  :: A(*)
    COMPLEX(SP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE CTPSV
   SUBROUTINE ZTPSV(UPLO,TRANS,DIAG,N,A,X,INCX)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO, TRANS, DIAG
    INTEGER, INTENT(IN)                      :: N
    COMPLEX(DP), INTENT(IN)                  :: A(*)
    COMPLEX(DP), INTENT(INOUT)               :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
   END SUBROUTINE ZTPSV
END INTERFACE BL_TPSV
INTERFACE BL_GER
   SUBROUTINE SGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: M, N
    REAL(SP), INTENT(IN)                     :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(IN)                     :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY, LDA
    REAL(SP), INTENT(INOUT)                  :: A(LDA,*)
   END SUBROUTINE SGER
   SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: M, N
    REAL(DP), INTENT(IN)                     :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(IN)                     :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY, LDA
    REAL(DP), INTENT(INOUT)                  :: A(LDA,*)
   END SUBROUTINE DGER
END INTERFACE BL_GER
INTERFACE BL_GERU
   SUBROUTINE CGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: M, N
    COMPLEX(SP), INTENT(IN)                  :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(IN)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY, LDA
    COMPLEX(SP), INTENT(INOUT)               :: A(LDA,*)
   END SUBROUTINE CGERU
   SUBROUTINE ZGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: M, N
    COMPLEX(DP), INTENT(IN)                  :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(IN)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY, LDA
    COMPLEX(DP), INTENT(INOUT)               :: A(LDA,*)
   END SUBROUTINE ZGERU
END INTERFACE BL_GERU
INTERFACE BL_GERC
   SUBROUTINE CGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: M, N
    COMPLEX(SP), INTENT(IN)                  :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(IN)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY, LDA
    COMPLEX(SP), INTENT(INOUT)               :: A(LDA,*)
   END SUBROUTINE CGERC
   SUBROUTINE ZGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    INTEGER, INTENT(IN)                      :: M, N
    COMPLEX(DP), INTENT(IN)                  :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(IN)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCY, LDA
    COMPLEX(DP), INTENT(INOUT)               :: A(LDA,*)
   END SUBROUTINE ZGERC
END INTERFACE BL_GERC
INTERFACE BL_HER
   SUBROUTINE CHER(UPLO,N,ALPHA,X,INCX,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: ALPHA
    COMPLEX(SP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX, LDA
    COMPLEX(SP), INTENT(INOUT)               :: A(LDA,*)
   END SUBROUTINE CHER
   SUBROUTINE ZHER(UPLO,N,ALPHA,X,INCX,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: ALPHA
    COMPLEX(DP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX, LDA
    COMPLEX(DP), INTENT(INOUT)               :: A(LDA,*)
   END SUBROUTINE ZHER
END INTERFACE BL_HER
INTERFACE BL_HPR
   SUBROUTINE CHPR(UPLO,N,ALPHA,X,INCX,AP)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: ALPHA
    COMPLEX(SP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(INOUT)               :: AP(*)
   END SUBROUTINE CHPR
   SUBROUTINE ZHPR(UPLO,N,ALPHA,X,INCX,AP)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: ALPHA
    COMPLEX(DP), INTENT(IN)                  :: X(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(INOUT)               :: AP(*)
   END SUBROUTINE ZHPR
END INTERFACE BL_HPR
INTERFACE BL_HER2
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
END INTERFACE BL_HER2
INTERFACE BL_HPR2
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
END INTERFACE BL_HPR2
INTERFACE BL_SYR
   SUBROUTINE SSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX, LDA
    REAL(SP), INTENT(INOUT)                  :: A(LDA,*)
   END SUBROUTINE SSYR
   SUBROUTINE DSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX, LDA
    REAL(DP), INTENT(INOUT)                  :: A(LDA,*)
   END SUBROUTINE DSYR
END INTERFACE BL_SYR
INTERFACE BL_SPR
   SUBROUTINE SSPR(UPLO,N,ALPHA,X,INCX,AP)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(SP), INTENT(IN)                     :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(INOUT)                  :: AP(*)
   END SUBROUTINE SSPR
   SUBROUTINE DSPR(UPLO,N,ALPHA,X,INCX,AP)
    INTEGER, PARAMETER :: SP=KIND(0.0),DP=KIND(0.0D0)
    CHARACTER(LEN=1), INTENT(IN)             :: UPLO
    INTEGER, INTENT(IN)                      :: N
    REAL(DP), INTENT(IN)                     :: ALPHA, X(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(INOUT)                  :: AP(*)
   END SUBROUTINE DSPR
END INTERFACE BL_SPR
INTERFACE BL_SYR2
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
END INTERFACE BL_SYR2
INTERFACE BL_SPR2
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
END INTERFACE BL_SPR2
!
! Level 3 BLAS
!
INTERFACE BL_GEMM
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
END INTERFACE BL_GEMM
INTERFACE BL_SYMM
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
END INTERFACE BL_SYMM
INTERFACE BL_HEMM
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
END INTERFACE BL_HEMM
INTERFACE BL_SYRK
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
END INTERFACE BL_SYRK
INTERFACE BL_HERK
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
END INTERFACE BL_HERK
INTERFACE BL_SYR2K
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
END INTERFACE BL_SYR2K
INTERFACE BL_HER2K
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
END INTERFACE BL_HER2K
INTERFACE BL_TRMM
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
END INTERFACE BL_TRMM
INTERFACE BL_TRSM
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
END INTERFACE BL_TRSM
END MODULE f77_blas_generic
