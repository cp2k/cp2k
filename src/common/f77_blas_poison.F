! *****************************************************************************
MODULE f77_blas_poison
!
! this routine 'poisons' the BLAS routines, i.e.
! if this module is USEd, direct calls to BLAS routines are not
! possible (unless the USE statement is appropriate )
! since there will be an argument mismatch. This provides an easy way
! to keep a code clean of direct BLAS calls.
!
! This file can be useful together with the f77_blas_generic module
! which provides generic interfaces to the BLAS routines. In that case
! only generic calls to BLAS routines will be allowed.
!
! this file is the public domain. Joost VandeVondele.
!

PRIVATE

PUBLIC :: SROTG,  DROTG, CROTG,  ZROTG, SROTMG, DROTMG, SROT,   DROT,  &
          ZROT,   CSROT, ZDROT,  SROTM, DROTM,  SSWAP,  DSWAP,  CSWAP, &
          ZSWAP,  SSCAL, DSCAL,  CSCAL, ZSCAL,  CSSCAL, ZDSCAL, SCOPY, &
          DCOPY,  CCOPY, ZCOPY,  SAXPY, DAXPY,  CAXPY,  ZAXPY,  SDOT,  &
          DDOT,   CDOTU, ZDOTU,  CDOTC, ZDOTC,  SNRM2,  DNRM2,  SCNRM2,&
          DZNRM2, SASUM, SCASUM, DASUM, DZASUM, ISAMAX, IDAMAX, ICAMAX,&
          IZAMAX

PUBLIC :: SGEMV, DGEMV, CGEMV, ZGEMV, SGBMV, DGBMV, CGBMV, ZGBMV, &
          CHEMV, ZHEMV, CHBMV, ZHBMV, CHPMV, ZHPMV, SSYMV, DSYMV, &
          SSBMV, DSBMV, SSPMV, DSPMV, STRMV, DTRMV, CTRMV, ZTRMV, &
          STBMV, DTBMV, CTBMV, ZTBMV, STPMV, DTPMV, CTPMV, ZTPMV, &
          STRSV, DTRSV, CTRSV, ZTRSV, STBSV, DTBSV, CTBSV, ZTBSV, &
          STPSV, DTPSV, CTPSV, ZTPSV, SGER,  DGER,  CGERU, ZGERU, &
          CGERC, ZGERC, CHER,  ZHER,  CHPR,  ZHPR,  CHER2, ZHER2, &
          CHPR2, ZHPR2, SSYR,  DSYR,  SSPR,  DSPR,  SSYR2, DSYR2, &
          SSPR2, DSPR2

PUBLIC :: SGEMM,  DGEMM,  CGEMM,  ZGEMM,  SSYMM,  DSYMM,  CSYMM,  ZSYMM, &
          CHEMM,  ZHEMM,  SSYRK,  DSYRK,  CSYRK,  ZSYRK,  CHERK,  ZHERK, &
          SSYR2K, DSYR2K, CSYR2K, ZSYR2K, CHER2K, ZHER2K, STRMM,  DTRMM, &
          CTRMM,  ZTRMM,  STRSM,  DTRSM,  CTRSM,  ZTRSM,  SDSDOT, DSDOT, &
          DCABS1, LSAME,  SCABS1

CONTAINS

! *****************************************************************************
 SUBROUTINE SROTG
 END SUBROUTINE SROTG
! *****************************************************************************
 SUBROUTINE DROTG
 END SUBROUTINE DROTG
! *****************************************************************************
 SUBROUTINE CROTG
 END SUBROUTINE CROTG
! *****************************************************************************
 SUBROUTINE ZROTG
 END SUBROUTINE ZROTG
! *****************************************************************************
 SUBROUTINE SROTMG
 END SUBROUTINE SROTMG
! *****************************************************************************
 SUBROUTINE DROTMG
 END SUBROUTINE DROTMG
! *****************************************************************************
 SUBROUTINE SROT
 END SUBROUTINE SROT
! *****************************************************************************
 SUBROUTINE DROT
 END SUBROUTINE DROT
! *****************************************************************************
 SUBROUTINE ZROT
 END SUBROUTINE ZROT
! *****************************************************************************
 SUBROUTINE CSROT
 END SUBROUTINE CSROT
! *****************************************************************************
 SUBROUTINE ZDROT
 END SUBROUTINE ZDROT
! *****************************************************************************
 SUBROUTINE SROTM
 END SUBROUTINE SROTM
! *****************************************************************************
 SUBROUTINE DROTM
 END SUBROUTINE DROTM
! *****************************************************************************
 SUBROUTINE SSWAP
 END SUBROUTINE SSWAP
! *****************************************************************************
 SUBROUTINE DSWAP
 END SUBROUTINE DSWAP
! *****************************************************************************
 SUBROUTINE CSWAP
 END SUBROUTINE CSWAP
! *****************************************************************************
 SUBROUTINE ZSWAP
 END SUBROUTINE ZSWAP
! *****************************************************************************
 SUBROUTINE SSCAL
 END SUBROUTINE SSCAL
! *****************************************************************************
 SUBROUTINE DSCAL
 END SUBROUTINE DSCAL
! *****************************************************************************
 SUBROUTINE CSCAL
 END SUBROUTINE CSCAL
! *****************************************************************************
 SUBROUTINE ZSCAL
 END SUBROUTINE ZSCAL
! *****************************************************************************
 SUBROUTINE CSSCAL
 END SUBROUTINE CSSCAL
! *****************************************************************************
 SUBROUTINE ZDSCAL
 END SUBROUTINE ZDSCAL
! *****************************************************************************
 SUBROUTINE SCOPY
 END SUBROUTINE SCOPY
! *****************************************************************************
 SUBROUTINE DCOPY
 END SUBROUTINE DCOPY
! *****************************************************************************
 SUBROUTINE CCOPY
 END SUBROUTINE CCOPY
! *****************************************************************************
 SUBROUTINE ZCOPY
 END SUBROUTINE ZCOPY
! *****************************************************************************
 SUBROUTINE SAXPY
 END SUBROUTINE SAXPY
! *****************************************************************************
 SUBROUTINE DAXPY
 END SUBROUTINE DAXPY
! *****************************************************************************
 SUBROUTINE CAXPY
 END SUBROUTINE CAXPY
! *****************************************************************************
 SUBROUTINE ZAXPY
 END SUBROUTINE ZAXPY
! *****************************************************************************
 FUNCTION SDOT()
   SDOT=0
 END FUNCTION SDOT
! *****************************************************************************
 FUNCTION DDOT()
   DDOT=0
 END FUNCTION DDOT
! *****************************************************************************
 FUNCTION CDOTU()
   CDOTU=0
 END FUNCTION CDOTU
! *****************************************************************************
 FUNCTION ZDOTU()
   ZDOTU=0
 END FUNCTION ZDOTU
! *****************************************************************************
 FUNCTION CDOTC()
   CDOTC=0
 END FUNCTION CDOTC
! *****************************************************************************
 FUNCTION ZDOTC()
   ZDOTC=0
 END FUNCTION ZDOTC
! *****************************************************************************
 FUNCTION SNRM2()
   SNRM2=0
 END FUNCTION SNRM2
! *****************************************************************************
 FUNCTION DNRM2()
   DNRM2=0
 END FUNCTION DNRM2
! *****************************************************************************
 FUNCTION SCNRM2()
   SCNRM2=0
 END FUNCTION SCNRM2
! *****************************************************************************
 FUNCTION DZNRM2()
   DZNRM2=0
 END FUNCTION DZNRM2
! *****************************************************************************
 FUNCTION SASUM()
   SASUM=0
 END FUNCTION SASUM
! *****************************************************************************
 FUNCTION SCASUM()
   SCASUM=0
 END FUNCTION SCASUM
! *****************************************************************************
 FUNCTION DASUM()
   DASUM=0
 END FUNCTION DASUM
! *****************************************************************************
 FUNCTION DZASUM()
   DZASUM=0
 END FUNCTION DZASUM
! *****************************************************************************
 FUNCTION ISAMAX()
   ISAMAX=0
 END FUNCTION ISAMAX
! *****************************************************************************
 FUNCTION IDAMAX()
   IDAMAX=0
 END FUNCTION IDAMAX
! *****************************************************************************
 FUNCTION ICAMAX()
   ICAMAX=0
 END FUNCTION ICAMAX
! *****************************************************************************
 FUNCTION IZAMAX()
   IZAMAX=0
 END FUNCTION IZAMAX
! *****************************************************************************
 SUBROUTINE SGEMV
 END SUBROUTINE SGEMV
! *****************************************************************************
 SUBROUTINE DGEMV
 END SUBROUTINE DGEMV
! *****************************************************************************
 SUBROUTINE CGEMV
 END SUBROUTINE CGEMV
! *****************************************************************************
 SUBROUTINE ZGEMV
 END SUBROUTINE ZGEMV
! *****************************************************************************
 SUBROUTINE SGBMV
 END SUBROUTINE SGBMV
! *****************************************************************************
 SUBROUTINE DGBMV
 END SUBROUTINE DGBMV
! *****************************************************************************
 SUBROUTINE CGBMV
 END SUBROUTINE CGBMV
! *****************************************************************************
 SUBROUTINE ZGBMV
 END SUBROUTINE ZGBMV
! *****************************************************************************
 SUBROUTINE CHEMV
 END SUBROUTINE CHEMV
! *****************************************************************************
 SUBROUTINE ZHEMV
 END SUBROUTINE ZHEMV
! *****************************************************************************
 SUBROUTINE CHBMV
 END SUBROUTINE CHBMV
! *****************************************************************************
 SUBROUTINE ZHBMV
 END SUBROUTINE ZHBMV
! *****************************************************************************
 SUBROUTINE CHPMV
 END SUBROUTINE CHPMV
! *****************************************************************************
 SUBROUTINE ZHPMV
 END SUBROUTINE ZHPMV
! *****************************************************************************
 SUBROUTINE SSYMV
 END SUBROUTINE SSYMV
! *****************************************************************************
 SUBROUTINE DSYMV
 END SUBROUTINE DSYMV
! *****************************************************************************
 SUBROUTINE SSBMV
 END SUBROUTINE SSBMV
! *****************************************************************************
 SUBROUTINE DSBMV
 END SUBROUTINE DSBMV
! *****************************************************************************
 SUBROUTINE SSPMV
 END SUBROUTINE SSPMV
! *****************************************************************************
 SUBROUTINE DSPMV
 END SUBROUTINE DSPMV
! *****************************************************************************
 SUBROUTINE STRMV
 END SUBROUTINE STRMV
! *****************************************************************************
 SUBROUTINE DTRMV
 END SUBROUTINE DTRMV
! *****************************************************************************
 SUBROUTINE CTRMV
 END SUBROUTINE CTRMV
! *****************************************************************************
 SUBROUTINE ZTRMV
 END SUBROUTINE ZTRMV
! *****************************************************************************
 SUBROUTINE STBMV
 END SUBROUTINE STBMV
! *****************************************************************************
 SUBROUTINE DTBMV
 END SUBROUTINE DTBMV
! *****************************************************************************
 SUBROUTINE CTBMV
 END SUBROUTINE CTBMV
! *****************************************************************************
 SUBROUTINE ZTBMV
 END SUBROUTINE ZTBMV
! *****************************************************************************
 SUBROUTINE STPMV
 END SUBROUTINE STPMV
! *****************************************************************************
 SUBROUTINE DTPMV
 END SUBROUTINE DTPMV
! *****************************************************************************
 SUBROUTINE CTPMV
 END SUBROUTINE CTPMV
! *****************************************************************************
 SUBROUTINE ZTPMV
 END SUBROUTINE ZTPMV
! *****************************************************************************
 SUBROUTINE STRSV
 END SUBROUTINE STRSV
! *****************************************************************************
 SUBROUTINE DTRSV
 END SUBROUTINE DTRSV
! *****************************************************************************
 SUBROUTINE CTRSV
 END SUBROUTINE CTRSV
! *****************************************************************************
 SUBROUTINE ZTRSV
 END SUBROUTINE ZTRSV
! *****************************************************************************
 SUBROUTINE STBSV
 END SUBROUTINE STBSV
! *****************************************************************************
 SUBROUTINE DTBSV
 END SUBROUTINE DTBSV
! *****************************************************************************
 SUBROUTINE CTBSV
 END SUBROUTINE CTBSV
! *****************************************************************************
 SUBROUTINE ZTBSV
 END SUBROUTINE ZTBSV
! *****************************************************************************
 SUBROUTINE STPSV
 END SUBROUTINE STPSV
! *****************************************************************************
 SUBROUTINE DTPSV
 END SUBROUTINE DTPSV
! *****************************************************************************
 SUBROUTINE CTPSV
 END SUBROUTINE CTPSV
! *****************************************************************************
 SUBROUTINE ZTPSV
 END SUBROUTINE ZTPSV
! *****************************************************************************
 SUBROUTINE SGER
 END SUBROUTINE SGER
! *****************************************************************************
 SUBROUTINE DGER
 END SUBROUTINE DGER
! *****************************************************************************
 SUBROUTINE CGERU
 END SUBROUTINE CGERU
! *****************************************************************************
 SUBROUTINE ZGERU
 END SUBROUTINE ZGERU
! *****************************************************************************
 SUBROUTINE CGERC
 END SUBROUTINE CGERC
! *****************************************************************************
 SUBROUTINE ZGERC
 END SUBROUTINE ZGERC
! *****************************************************************************
 SUBROUTINE CHER
 END SUBROUTINE CHER
! *****************************************************************************
 SUBROUTINE ZHER
 END SUBROUTINE ZHER
! *****************************************************************************
 SUBROUTINE CHPR
 END SUBROUTINE CHPR
! *****************************************************************************
 SUBROUTINE ZHPR
 END SUBROUTINE ZHPR
! *****************************************************************************
 SUBROUTINE CHER2
 END SUBROUTINE CHER2
! *****************************************************************************
 SUBROUTINE ZHER2
 END SUBROUTINE ZHER2
! *****************************************************************************
 SUBROUTINE CHPR2
 END SUBROUTINE CHPR2
! *****************************************************************************
 SUBROUTINE ZHPR2
 END SUBROUTINE ZHPR2
! *****************************************************************************
 SUBROUTINE SSYR
 END SUBROUTINE SSYR
! *****************************************************************************
 SUBROUTINE DSYR
 END SUBROUTINE DSYR
! *****************************************************************************
 SUBROUTINE SSPR
 END SUBROUTINE SSPR
! *****************************************************************************
 SUBROUTINE DSPR
 END SUBROUTINE DSPR
! *****************************************************************************
 SUBROUTINE SSYR2
 END SUBROUTINE SSYR2
! *****************************************************************************
 SUBROUTINE DSYR2
 END SUBROUTINE DSYR2
! *****************************************************************************
 SUBROUTINE SSPR2
 END SUBROUTINE SSPR2
! *****************************************************************************
 SUBROUTINE DSPR2
 END SUBROUTINE DSPR2
! *****************************************************************************
 SUBROUTINE SGEMM
 END SUBROUTINE SGEMM
! *****************************************************************************
 SUBROUTINE DGEMM
 END SUBROUTINE DGEMM
! *****************************************************************************
 SUBROUTINE CGEMM
 END SUBROUTINE CGEMM
! *****************************************************************************
 SUBROUTINE ZGEMM
 END SUBROUTINE ZGEMM
! *****************************************************************************
 SUBROUTINE SSYMM
 END SUBROUTINE SSYMM
! *****************************************************************************
 SUBROUTINE DSYMM
 END SUBROUTINE DSYMM
! *****************************************************************************
 SUBROUTINE CSYMM
 END SUBROUTINE CSYMM
! *****************************************************************************
 SUBROUTINE ZSYMM
 END SUBROUTINE ZSYMM
! *****************************************************************************
 SUBROUTINE CHEMM
 END SUBROUTINE CHEMM
! *****************************************************************************
 SUBROUTINE ZHEMM
 END SUBROUTINE ZHEMM
! *****************************************************************************
 SUBROUTINE SSYRK
 END SUBROUTINE SSYRK
! *****************************************************************************
 SUBROUTINE DSYRK
 END SUBROUTINE DSYRK
! *****************************************************************************
 SUBROUTINE CSYRK
 END SUBROUTINE CSYRK
! *****************************************************************************
 SUBROUTINE ZSYRK
 END SUBROUTINE ZSYRK
! *****************************************************************************
 SUBROUTINE CHERK
 END SUBROUTINE CHERK
! *****************************************************************************
 SUBROUTINE ZHERK
 END SUBROUTINE ZHERK
! *****************************************************************************
 SUBROUTINE SSYR2K
 END SUBROUTINE SSYR2K
! *****************************************************************************
 SUBROUTINE DSYR2K
 END SUBROUTINE DSYR2K
! *****************************************************************************
 SUBROUTINE CSYR2K
 END SUBROUTINE CSYR2K
! *****************************************************************************
 SUBROUTINE ZSYR2K
 END SUBROUTINE ZSYR2K
! *****************************************************************************
 SUBROUTINE CHER2K
 END SUBROUTINE CHER2K
! *****************************************************************************
 SUBROUTINE ZHER2K
 END SUBROUTINE ZHER2K
! *****************************************************************************
 SUBROUTINE STRMM
 END SUBROUTINE STRMM
! *****************************************************************************
 SUBROUTINE DTRMM
 END SUBROUTINE DTRMM
! *****************************************************************************
 SUBROUTINE CTRMM
 END SUBROUTINE CTRMM
! *****************************************************************************
 SUBROUTINE ZTRMM
 END SUBROUTINE ZTRMM
! *****************************************************************************
 SUBROUTINE STRSM
 END SUBROUTINE STRSM
! *****************************************************************************
 SUBROUTINE DTRSM
 END SUBROUTINE DTRSM
! *****************************************************************************
 SUBROUTINE CTRSM
 END SUBROUTINE CTRSM
! *****************************************************************************
 SUBROUTINE ZTRSM
 END SUBROUTINE ZTRSM
! *****************************************************************************
 FUNCTION SDSDOT()
   SDSDOT=0
 END FUNCTION SDSDOT
! *****************************************************************************
 FUNCTION DSDOT()
   DSDOT=0
 END FUNCTION DSDOT
! *****************************************************************************
 FUNCTION DCABS1()
   DCABS1=0
 END FUNCTION DCABS1
! *****************************************************************************
 FUNCTION LSAME()
   LSAME=0
 END FUNCTION LSAME
! *****************************************************************************
 FUNCTION SCABS1()
   SCABS1=0
 END FUNCTION SCABS1
END MODULE f77_blas_poison
