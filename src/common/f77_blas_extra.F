! *****************************************************************************
MODULE f77_blas_extra
!
! Provides extra generic interfaces for BLAS-like routines not distributed
! with BLAS (http://www.netlib.org/blas/) Jan 2006.
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
IMPLICIT NONE
PRIVATE

PUBLIC :: BL_COPY

INTERFACE BL_COPY
  MODULE PROCEDURE SDCOPY,DSCOPY,ZCCOPY,CZCOPY
END INTERFACE BL_COPY

CONTAINS

!
! Level 1 BLAS, based on netlib blas
!
! *****************************************************************************
!> \brief ...
!> \param N ...
!> \param X ...
!> \param INCX ...
!> \param Y ...
!> \param INCY ...
! *****************************************************************************
   SUBROUTINE SDCOPY(N,X,INCX,Y,INCY)
    INTEGER, INTENT(IN)                      :: N
    INTEGER, PARAMETER                       :: DP = KIND(0.0D0), &
                                                SP = KIND(0.0)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(DP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(SP), INTENT(IN)                     :: X(*)

    INTEGER                                  :: i, ix, iy, m, mp1

!

      IF(n.le.0)RETURN
      IF(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      IF(incx.lt.0)ix = (-n+1)*incx + 1
      IF(incy.lt.0)iy = (-n+1)*incy + 1
      DO 10 i = 1,n
        y(iy) = x(ix)
        ix = ix + incx
        iy = iy + incy
   10 CONTINUE
      RETURN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = MOD(n,7)
      IF( m .EQ. 0 ) go to 40
      DO 30 i = 1,m
        y(i) = x(i)
   30 CONTINUE
      IF( n .LT. 7 ) RETURN
   40 mp1 = m + 1
      DO 50 i = mp1,n,7
        y(i) = x(i)
        y(i + 1) = x(i + 1)
        y(i + 2) = x(i + 2)
        y(i + 3) = x(i + 3)
        y(i + 4) = x(i + 4)
        y(i + 5) = x(i + 5)
        y(i + 6) = x(i + 6)
   50 CONTINUE
      RETURN

   END SUBROUTINE SDCOPY

! *****************************************************************************
!> \brief ...
!> \param N ...
!> \param X ...
!> \param INCX ...
!> \param Y ...
!> \param INCY ...
! *****************************************************************************
   SUBROUTINE DSCOPY(N,X,INCX,Y,INCY)
    INTEGER, INTENT(IN)                      :: N
    INTEGER, PARAMETER                       :: DP = KIND(0.0D0), &
                                                SP = KIND(0.0)
    INTEGER, INTENT(IN)                      :: INCY
    REAL(SP), INTENT(INOUT)                  :: Y(*)
    INTEGER, INTENT(IN)                      :: INCX
    REAL(DP), INTENT(IN)                     :: X(*)

    INTEGER                                  :: i, ix, iy, m, mp1

!

      IF(n.le.0)RETURN
      IF(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      IF(incx.lt.0)ix = (-n+1)*incx + 1
      IF(incy.lt.0)iy = (-n+1)*incy + 1
      DO 10 i = 1,n
        y(iy) = x(ix)
        ix = ix + incx
        iy = iy + incy
   10 CONTINUE
      RETURN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = MOD(n,7)
      IF( m .EQ. 0 ) go to 40
      DO 30 i = 1,m
        y(i) = x(i)
   30 CONTINUE
      IF( n .LT. 7 ) RETURN
   40 mp1 = m + 1
      DO 50 i = mp1,n,7
        y(i) = x(i)
        y(i + 1) = x(i + 1)
        y(i + 2) = x(i + 2)
        y(i + 3) = x(i + 3)
        y(i + 4) = x(i + 4)
        y(i + 5) = x(i + 5)
        y(i + 6) = x(i + 6)
   50 CONTINUE
      RETURN

   END SUBROUTINE DSCOPY

! *****************************************************************************
!> \brief ...
!> \param N ...
!> \param X ...
!> \param INCX ...
!> \param Y ...
!> \param INCY ...
! *****************************************************************************
   SUBROUTINE CZCOPY(N,X,INCX,Y,INCY)
    INTEGER, INTENT(IN)                      :: N
    INTEGER, PARAMETER                       :: DP = KIND(0.0D0), &
                                                SP = KIND(0.0)
    INTEGER, INTENT(IN)                      :: INCY
    COMPLEX(DP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(SP), INTENT(IN)                  :: X(*)

    INTEGER                                  :: i, ix, iy, m, mp1

!

      IF(n.le.0)RETURN
      IF(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      IF(incx.lt.0)ix = (-n+1)*incx + 1
      IF(incy.lt.0)iy = (-n+1)*incy + 1
      DO 10 i = 1,n
        y(iy) = x(ix)
        ix = ix + incx
        iy = iy + incy
   10 CONTINUE
      RETURN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = MOD(n,7)
      IF( m .EQ. 0 ) go to 40
      DO 30 i = 1,m
        y(i) = x(i)
   30 CONTINUE
      IF( n .LT. 7 ) RETURN
   40 mp1 = m + 1
      DO 50 i = mp1,n,7
        y(i) = x(i)
        y(i + 1) = x(i + 1)
        y(i + 2) = x(i + 2)
        y(i + 3) = x(i + 3)
        y(i + 4) = x(i + 4)
        y(i + 5) = x(i + 5)
        y(i + 6) = x(i + 6)
   50 CONTINUE
      RETURN

   END SUBROUTINE CZCOPY

! *****************************************************************************
!> \brief ...
!> \param N ...
!> \param X ...
!> \param INCX ...
!> \param Y ...
!> \param INCY ...
! *****************************************************************************
   SUBROUTINE ZCCOPY(N,X,INCX,Y,INCY)
    INTEGER, INTENT(IN)                      :: N
    INTEGER, PARAMETER                       :: DP = KIND(0.0D0), &
                                                SP = KIND(0.0)
    INTEGER, INTENT(IN)                      :: INCY
    COMPLEX(SP), INTENT(INOUT)               :: Y(*)
    INTEGER, INTENT(IN)                      :: INCX
    COMPLEX(DP), INTENT(IN)                  :: X(*)

    INTEGER                                  :: i, ix, iy, m, mp1

!

      IF(n.le.0)RETURN
      IF(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      IF(incx.lt.0)ix = (-n+1)*incx + 1
      IF(incy.lt.0)iy = (-n+1)*incy + 1
      DO 10 i = 1,n
        y(iy) = x(ix)
        ix = ix + incx
        iy = iy + incy
   10 CONTINUE
      RETURN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = MOD(n,7)
      IF( m .EQ. 0 ) go to 40
      DO 30 i = 1,m
        y(i) = x(i)
   30 CONTINUE
      IF( n .LT. 7 ) RETURN
   40 mp1 = m + 1
      DO 50 i = mp1,n,7
        y(i) = x(i)
        y(i + 1) = x(i + 1)
        y(i + 2) = x(i + 2)
        y(i + 3) = x(i + 3)
        y(i + 4) = x(i + 4)
        y(i + 5) = x(i + 5)
        y(i + 6) = x(i + 6)
   50 CONTINUE
      RETURN

   END SUBROUTINE ZCCOPY

END MODULE f77_blas_extra
