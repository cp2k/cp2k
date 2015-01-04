!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculate spherical harmonics
!> \note
!>      Spherical Harmonics
!>          Numerical Stability up to L=15
!>          Accuracy > 1.E-12 up to L=15 tested
!>          Definition is consistent with orbital_transformation_matrices
!>      Clebsch-Gordon Coefficients
!>          Tested up to l=7 (i.e. L=14)
!> \todo
!>      1) Check if this definition is consistent with the
!>         Slater-Koster module
!> \par History
!>      JGH 28-Feb-2002 : Change of sign convention (-1^m)
!>      JGH  1-Mar-2002 : Clebsch-Gordon Coefficients
!> \author JGH 6-Oct-2000
! *****************************************************************************
MODULE spherical_harmonics

  
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE mathconstants,                   ONLY: dfac,&
                                             fac,&
                                             maxfac,&
                                             pi
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: y_lm
  PUBLIC :: clebsch_gordon_deallocate, clebsch_gordon_init,&
            clebsch_gordon
  PUBLIC :: legendre, dlegendre

  INTERFACE y_lm
     MODULE PROCEDURE rvy_lm, rry_lm, cvy_lm, ccy_lm
  END INTERFACE

  INTERFACE clebsch_gordon
     MODULE PROCEDURE clebsch_gordon_real, clebsch_gordon_complex
  END INTERFACE

  INTERFACE clebsch_gordon_getm
     MODULE PROCEDURE getm
  END INTERFACE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'spherical_harmonics'

  REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: cg_table
  INTEGER :: lmax=-1
  REAL(KIND=dp) :: osq2, sq2

CONTAINS

! Clebsch-Gordon Coefficients

! *****************************************************************************
!> \brief ...
!> \param l1 ...
!> \param m1 ...
!> \param l2 ...
!> \param m2 ...
!> \param clm ...
! *****************************************************************************
  SUBROUTINE clebsch_gordon_complex (l1,m1,l2,m2,clm)
    INTEGER, INTENT(IN)                      :: l1, m1, l2, m2
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: clm

    CHARACTER(len=*), PARAMETER :: routineN = 'clebsch_gordon_complex', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icase, ind, l, lm, lp, n

     l = l1 + l2
     IF (l > lmax) CALL clebsch_gordon_init ( l )
     n = l/2 + 1
     IF (n > SIZE(clm)) CALL stop_program(routineN,moduleN,__LINE__,&
                                          "Array too small")

     IF ( (m1 >= 0 .AND. m2 >= 0) .OR. (m1 < 0 .AND. m2 < 0) ) THEN
       icase = 1
     ELSE
       icase = 2
     END IF
     ind = order(l1,m1,l2,m2)

     DO lp = MOD(l,2),l,2
       lm = lp/2+1
       clm(lm) = cg_table(ind,lm,icase)
     END DO

  END SUBROUTINE clebsch_gordon_complex

! *****************************************************************************
!> \brief ...
!> \param l1 ...
!> \param m1 ...
!> \param l2 ...
!> \param m2 ...
!> \param rlm ...
! *****************************************************************************
  SUBROUTINE clebsch_gordon_real (l1,m1,l2,m2,rlm)
    INTEGER, INTENT(IN)                      :: l1, m1, l2, m2
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: rlm

    CHARACTER(len=*), PARAMETER :: routineN = 'clebsch_gordon_real', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icase1, icase2, ind, l, lm, &
                                                lp, mm(2), n
    REAL(KIND=dp)                            :: xsi

     l = l1 + l2
     IF (l > lmax) CALL clebsch_gordon_init ( l )
     n = l/2 + 1
     IF (n > SIZE(rlm,1)) CALL stop_program(routineN,moduleN,__LINE__,&
                                            "Array too small")

     ind = order(l1,m1,l2,m2)
     mm = getm(m1,m2)
     IF ( (m1 >= 0 .AND. m2 >= 0) .OR. (m1 < 0 .AND. m2 < 0) ) THEN
       icase1 = 1
       icase2 = 2
     ELSE
       icase1 = 2
       icase2 = 1
     END IF

     DO lp = MOD(l,2),l,2
       lm = lp/2+1
       xsi = get_factor(m1,m2,mm(1))
       rlm(lm,1) = xsi*cg_table(ind,lm,icase1)
       xsi = get_factor(m1,m2,mm(2))
       rlm(lm,2) = xsi*cg_table(ind,lm,icase2)
     END DO

  END SUBROUTINE clebsch_gordon_real

! *****************************************************************************
!> \brief ...
!> \param m1 ...
!> \param m2 ...
!> \retval m ...
! *****************************************************************************
  FUNCTION getm(m1,m2) RESULT(m)
    INTEGER, INTENT(IN)                      :: m1, m2
    INTEGER, DIMENSION(2)                    :: m

    INTEGER                                  :: mm, mp

     mp = m1 + m2
     mm = m1 - m2
     IF ( m1*m2 < 0  .OR. (m1*m2==0 .AND. (m1<0 .OR. m2<0))) THEN
       mp = -ABS(mp)
       mm = -ABS(mm)
     ELSE
       mp = ABS(mp)
       mm = ABS(mm)
     END IF
     m(1) = mp
     m(2) = mm
  END FUNCTION getm

! *****************************************************************************
!> \brief ...
!> \param m1 ...
!> \param m2 ...
!> \param m ...
!> \retval f ...
! *****************************************************************************
  FUNCTION get_factor(m1,m2,m) RESULT(f)
    INTEGER, INTENT(IN)                      :: m1, m2, m
    REAL(KIND=dp)                            :: f

    INTEGER                                  :: mx, my

     f = 1.0_dp
     IF ( ABS(m1)>=ABS(m2) ) THEN
       mx = m1
       my = m2
     ELSE
       mx = m2
       my = m1
     ENDIF
     IF ( mx*my == 0 ) THEN
       f = 1.0_dp
     ELSE IF ( m == 0 ) THEN
       IF(ABS(mx) /= ABS(my)) WRITE(*,'(A,3I6)') " 1) Illegal Case ",m1,m2,m
       IF(mx*my > 0) THEN
         f = 1.0_dp
       ELSE
         f = 0.0_dp
       END IF
     ELSE IF ( ABS(mx)+ABS(my) == m ) THEN
       f = osq2
       IF(mx<0) f = -osq2
     ELSE IF ( ABS(mx)+ABS(my) == -m ) THEN
       f = osq2
     ELSE IF ( ABS(mx)-ABS(my) == -m ) THEN
       IF(mx*my > 0) WRITE(*,'(A,3I6)') " 2) Illegal Case ",m1,m2,m
       IF( mx > 0 ) f = -osq2
       IF( mx < 0 ) f = osq2
     ELSE IF ( ABS(mx)-ABS(my) == m ) THEN
       IF(mx*my < 0) WRITE(*,'(A,3I6)') " 3) Illegal Case ",m1,m2,m
       f = osq2
     ELSE
       WRITE(*,'(A,3I6)') " 4) Illegal Case ",m1,m2,m
     END IF
  END FUNCTION get_factor

! *****************************************************************************
!> \brief ...
! *****************************************************************************
  SUBROUTINE clebsch_gordon_deallocate()
    CHARACTER(len=*), PARAMETER :: routineN = 'clebsch_gordon_deallocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr

    CALL timeset(routineN,handle)


     IF ( ALLOCATED ( cg_table ) ) THEN
        DEALLOCATE ( cg_table, STAT=ierr )
        IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"cg_table")
     END IF
    CALL timestop(handle)

  END SUBROUTINE clebsch_gordon_deallocate

! *****************************************************************************
!> \brief ...
!> \param l ...
! *****************************************************************************
  SUBROUTINE clebsch_gordon_init ( l )
    INTEGER, INTENT(IN)                      :: l

    CHARACTER(len=*), PARAMETER :: routineN = 'clebsch_gordon_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i1, i2, ierr, ix, iy, &
                                                l1, l2, lp, m, m1, m2, ml, &
                                                mp, n

     CALL timeset(routineN,handle)

     sq2 = SQRT(2.0_dp)
     osq2 = 1.0_dp/sq2

     IF ( l < 0 ) CALL stop_program(routineN,moduleN,__LINE__,"l < 0")
     IF ( ALLOCATED ( cg_table ) ) THEN
        DEALLOCATE ( cg_table, STAT=ierr )
        IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"cg_table")
     END IF
     ! maximum size of table
     n = (l**4 + 6*l**3 + 15*l**2 + 18*l + 8)/8
     m = l+1
     ALLOCATE ( cg_table(n,m,2), STAT=ierr )
     IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "cg_table",dp_size*2*n*m)
     lmax = l

     DO l1=0,lmax
       DO m1=0,l1
         iy = (l1*(l1+1))/2 + m1 + 1
         DO l2=l1,lmax
           ml=0
           IF ( l1==l2 ) ml = m1
           DO m2=ml,l2
             ix = (l2*(l2+1))/2 + m2 + 1
             i1 = (ix*(ix-1))/2 + iy
             DO lp=MOD(l1+l2,2),l1+l2,2
               i2 = lp/2 + 1
               mp = m2 + m1
               cg_table(i1,i2,1) = cgc(l1,m1,l2,m2,lp,mp)
               mp = ABS ( m2 - m1 )
               IF ( m2 >= m1 ) THEN
                  cg_table(i1,i2,2) = cgc(l1,m1,lp,mp,l2,m2)
               ELSE
                  cg_table(i1,i2,2) = cgc(l2,m2,lp,mp,l1,m1)
               END IF
             END DO
           END DO
         END DO
       END DO
     END DO

     CALL timestop(handle)

  END SUBROUTINE clebsch_gordon_init

! *****************************************************************************
!> \brief ...
!> \param l1 ...
!> \param m1 ...
!> \param l2 ...
!> \param m2 ...
!> \param lp ...
!> \param mp ...
!> \retval cgc ...
! *****************************************************************************
  FUNCTION cgc(l1,m1,l2,m2,lp,mp)
    INTEGER, INTENT(IN)                      :: l1, m1, l2, m2, lp, mp
    REAL(KIND=dp)                            :: cgc

    CHARACTER(len=*), PARAMETER :: routineN = 'cgc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: la, lb, ll, ma, mb, s, t, &
                                                tmax, tmin, z1, z2, z3, z4, z5
    REAL(KIND=dp)                            :: f1, f2, pref

! m1 >= 0; m2 >= 0; mp >= 0

    IF ( m1<0 .OR. m2<0 .OR. mp<0 ) THEN
       WRITE(*,*) l1,l2,lp
       WRITE(*,*) m1,m2,mp
       CALL stop_program(routineN,moduleN,__LINE__,"Illegal input values")
    END IF
    IF ( l2 < l1 ) THEN
      la = l2
      ma = m2
      lb = l1
      mb = m1
    ELSE
      la = l1
      ma = m1
      lb = l2
      mb = m2
    END IF

    IF ( MOD(la+lb+lp,2) == 0 .AND. la+lb >= lp .AND. lp >= lb-la &
         .AND. lb-mb >= 0 ) THEN
       ll = (2*lp+1) * (2*la+1) * (2*lb+1)
       pref = 1.0_dp/SQRT(4.0_dp*pi) * 0.5_dp * SQRT ( REAL(ll,dp) * &
              (sfac(lp-mp)/sfac(lp+mp)) * &
              (sfac(la-ma)/sfac(la+ma)) * (sfac(lb-mb)/sfac(lb+mb)) )
       s = (la+lb+lp)/2
       tmin = MAX ( 0, -lb+la-mp )
       tmax = MIN ( lb+la-mp, lp-mp, la-ma )
       f1 = REAL( 2 * (-1)**(s-lb-ma),KIND=dp) * (sfac(lb+mb)/sfac(lb-mb)) * &
            sfac(la+ma)/(sfac(s-lp)*sfac(s-lb)) * sfac(2*s-2*la)/sfac(s-la) * &
            (sfac(s)/sfac(2*s+1))
       f2 = 0.0_dp
       DO t=tmin,tmax
         z1 = lp+mp+t
         z2 = la+lb-mp-t
         z3 = lp-mp-t
         z4 = lb-la+mp+t
         z5 = la-ma-t
         f2 = f2 + (-1)**t * (sfac(z1)/(sfac(t)*sfac(z3)))* (sfac(z2)/(sfac(z4)*sfac(z5)))
       END DO
       cgc = pref * f1 * f2
    ELSE
       cgc = 0.0_dp
    END IF

  END FUNCTION cgc

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \retval fval ...
! *****************************************************************************
  FUNCTION sfac ( n ) RESULT (fval)
    INTEGER                                  :: n
    REAL(KIND=dp)                            :: fval

    INTEGER                                  :: i

    IF ( n > maxfac ) THEN
       fval = fac(maxfac)
       DO i = maxfac+1, n
          fval = REAL(i,dp)*fval
       END DO
    ELSE IF ( n >= 0 ) THEN
       fval = fac(n)
    ELSE
       fval = 1.0_dp
    END IF
  END FUNCTION sfac

! *****************************************************************************
!> \brief ...
!> \param l1 ...
!> \param m1 ...
!> \param l2 ...
!> \param m2 ...
!> \retval ind ...
! *****************************************************************************
  FUNCTION order(l1,m1,l2,m2) RESULT(ind)
    INTEGER, INTENT(IN)                      :: l1, m1, l2, m2
    INTEGER                                  :: ind

    INTEGER                                  :: i1, i2, ix, iy

    i1 = (l1*(l1+1))/2 + ABS(m1) + 1
    i2 = (l2*(l2+1))/2 + ABS(m2) + 1
    ix = MAX(i1,i2)
    iy = MIN(i1,i2)
    ind = (ix*(ix-1))/2+iy
  END FUNCTION order

! Calculation of Spherical Harmonics

! *****************************************************************************
!> \brief ...
!> \param r ...
!> \param y ...
!> \param l ...
!> \param m ...
! *****************************************************************************
  SUBROUTINE rvy_lm ( r, y, l, m )
!
! Real Spherical Harmonics
!                   _                   _
!                  |  [(2l+1)(l-|m|)!]   |1/2 m         cos(m p)   m>=0
!  Y_lm ( t, p ) = |---------------------|   P_l(cos(t))
!                  |[2Pi(1+d_m0)(l+|m|)!]|              sin(|m| p) m<0
!
! Input: r == (x,y,z) : normalised    x^2 + y^2 + z^2 = 1
!
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: y
    INTEGER, INTENT(IN)                      :: l, m

    CHARACTER(len=*), PARAMETER :: routineN = 'rvy_lm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: cp, lmm, lpm, pf, plm, rxy, &
                                                sp, t, z

  SELECT CASE ( l )
  CASE ( : -1 )
    CALL stop_program(routineN,moduleN,__LINE__,"Negative l value")
  CASE ( 0 )
    pf = SQRT ( 1.0_dp / ( 4.0_dp * pi ) )
    IF ( m /= 0 )  CALL stop_program(routineN,moduleN,__LINE__,&
                                     "l = 0 and m value out of bounds")
    y(:) = pf
  CASE ( 1 )
    SELECT CASE ( m )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "l = 1 and m value out of bounds")
    CASE ( 1 )
      pf = SQRT ( 3.0_dp / ( 4.0_dp * pi ) )
      y ( : ) = pf * r ( 1, : )
    CASE ( 0 )
      pf = SQRT ( 3.0_dp / ( 4.0_dp * pi ) )
      y ( : ) = pf * r ( 3, : )
    CASE ( -1 )
      pf = SQRT ( 3.0_dp / ( 4.0_dp * pi ) )
      y ( : ) = pf * r ( 2, : )
    END SELECT
  CASE ( 2 )
    SELECT CASE ( m )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "l = 2 and m value out of bounds")
    CASE ( 2 )
      pf = SQRT ( 15.0_dp / ( 16.0_dp * pi ) )
      y ( : ) = pf * ( r ( 1, : ) * r ( 1, : ) - r ( 2, : ) * r ( 2, : ) )
    CASE ( 1 )
      pf = SQRT ( 15.0_dp / ( 4.0_dp * pi ) )
      y ( : ) = pf * r ( 3, : ) * r ( 1, : )
    CASE ( 0 )
      pf = SQRT ( 5.0_dp / ( 16.0_dp * pi ) )
      y ( : ) = pf * ( 3.0_dp * r ( 3, : ) * r ( 3, : ) - 1.0_dp )
    CASE ( -1 )
      pf = SQRT ( 15.0_dp / ( 4.0_dp * pi ) )
      y ( : ) = pf * r ( 3, : ) * r ( 2, : )
    CASE ( -2 )
      pf = SQRT ( 15.0_dp / ( 16.0_dp * pi ) )
      y ( : ) = pf * 2.0_dp * r ( 1, : ) * r ( 2, : )
    END SELECT
  CASE ( 3 )
    SELECT CASE ( m )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "l = 3 and m value out of bounds")
    CASE ( 3 )
      pf = SQRT ( 35.0_dp / ( 32.0_dp * pi ) )
      y ( : ) = pf * r ( 1, : ) * ( r ( 1, : )**2 - 3.0_dp * r ( 2, : )**2 )
    CASE ( 2 )
      pf = SQRT ( 105.0_dp / ( 16.0_dp * pi ) )
      y ( : ) = pf * r ( 3, : ) * ( r ( 1, : )**2 - r ( 2, : )**2 )
    CASE ( 1 )
      pf = SQRT ( 21.0_dp / ( 32.0_dp * pi ) )
      y ( : ) = pf * r ( 1, : ) * ( 5.0_dp * r ( 3, : ) * r ( 3, : ) - 1.0_dp )
    CASE ( 0 )
      pf = SQRT ( 7.0_dp / ( 16.0_dp * pi ) )
      y ( : ) = pf * r ( 3, : ) * ( 5.0_dp * r ( 3, : ) * r ( 3, : ) - 3.0_dp )
    CASE ( -1 )
      pf = SQRT ( 21.0_dp / ( 32.0_dp * pi ) )
      y ( : ) = pf * r ( 2, : ) * ( 5.0_dp * r ( 3, : ) * r ( 3, : ) - 1.0_dp )
    CASE ( -2 )
      pf = SQRT ( 105.0_dp / ( 16.0_dp * pi ) )
      y ( : ) = pf * 2.0_dp * r ( 1, : ) * r ( 2, : ) * r ( 3, : )
    CASE ( -3 )
      pf = SQRT ( 35.0_dp / ( 32.0_dp * pi ) )
      y ( : ) = pf * r ( 2, : ) * ( 3.0_dp * r ( 1, : )**2 - r ( 2, : )**2 )
    END SELECT
  CASE DEFAULT
    IF ( m < -l .OR. m > l ) CALL stop_program(routineN,moduleN,__LINE__,&
                                               "m value out of bounds")
    lpm = fac ( l + ABS ( m ) )
    lmm = fac ( l - ABS ( m ) )
    IF ( m == 0 ) THEN
      t = 4.0_dp * pi
    ELSE
      t = 2.0_dp * pi
    END IF
    IF ( ABS ( lpm ) < EPSILON ( 1.0_dp ) ) THEN
      pf = REAL ( 2*l+1,KIND=dp)  / t
    ELSE
      pf = ( REAL ( 2*l+1,KIND=dp) * lmm ) / ( t * lpm )
    ENDIF
    pf = SQRT ( pf )
    DO i = 1, SIZE ( r ,2 )
      z = r ( 3, i )
!      plm = legendre ( z, l, m )
      plm = legendre ( z, l, ABS ( m ) )
      IF ( m == 0 ) THEN
        y ( i ) = pf * plm
      ELSE
        rxy = SQRT ( r ( 1, i ) ** 2  +  r ( 2, i ) ** 2 )
        IF ( rxy < EPSILON ( 1.0_dp ) ) THEN
          y ( i ) = 0.0_dp
        ELSE
          cp = r ( 1, i ) / rxy
          sp = r ( 2, i ) / rxy
          IF ( m > 0 ) THEN
            y ( i ) = pf * plm * cosn ( m, cp, sp )
          ELSE
            y ( i ) = pf * plm * sinn ( -m , cp, sp )
          END IF
        END IF
      END IF
    END DO
  END SELECT

  END SUBROUTINE rvy_lm

! *****************************************************************************
!> \brief ...
!> \param r ...
!> \param y ...
!> \param l ...
!> \param m ...
! *****************************************************************************
  SUBROUTINE rry_lm ( r, y, l, m )
!
! Real Spherical Harmonics
!                   _                   _
!                  |  [(2l+1)(l-|m|)!]   |1/2 m         cos(m p)   m>=0
!  Y_lm ( t, p ) = |---------------------|   P_l(cos(t))
!                  |[2Pi(1+d_m0)(l+|m|)!]|              sin(|m| p) m<0
!
! Input: r == (x,y,z) : normalised    x^2 + y^2 + z^2 = 1
!
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: r
    REAL(KIND=dp), INTENT(OUT)               :: y
    INTEGER, INTENT(IN)                      :: l, m

    CHARACTER(len=*), PARAMETER :: routineN = 'rry_lm', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: cp, lmm, lpm, pf, plm, rxy, &
                                                sp, t, z

  SELECT CASE ( l )
  CASE ( : -1 )
    CALL stop_program(routineN,moduleN,__LINE__,"Negative l value")
  CASE ( 0 )
    pf = SQRT ( 1.0_dp / ( 4.0_dp * pi ) )
    IF ( m /= 0 )  CALL stop_program(routineN,moduleN,__LINE__,&
                                     "l = 0 and m value out of bounds")
    y = pf
  CASE ( 1 )
    SELECT CASE ( m )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "l = 1 and m value out of bounds")
    CASE ( 1 )
      pf = SQRT ( 3.0_dp / ( 4.0_dp * pi ) )
      y = pf * r ( 1 )
    CASE ( 0 )
      pf = SQRT ( 3.0_dp / ( 4.0_dp * pi ) )
      y = pf * r ( 3 )
    CASE ( -1 )
      pf = SQRT ( 3.0_dp / ( 4.0_dp * pi ) )
      y = pf * r ( 2 )
    END SELECT
  CASE ( 2 )
    SELECT CASE ( m )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "l = 2 and m value out of bounds")
    CASE ( 2 )
      pf = SQRT ( 15.0_dp / ( 16.0_dp * pi ) )
      y = pf * ( r ( 1 ) * r ( 1 ) - r ( 2 ) * r ( 2 ) )
    CASE ( 1 )
      pf = SQRT ( 15.0_dp / ( 4.0_dp * pi ) )
      y = pf * r ( 3 ) * r ( 1 )
    CASE ( 0 )
      pf = SQRT ( 5.0_dp / ( 16.0_dp * pi ) )
      y = pf * ( 3.0_dp * r ( 3 ) * r ( 3 ) - 1.0_dp )
    CASE ( -1 )
      pf = SQRT ( 15.0_dp / ( 4.0_dp * pi ) )
      y = pf * r ( 3 ) * r ( 2 )
    CASE ( -2 )
      pf = SQRT ( 15.0_dp / ( 16.0_dp * pi ) )
      y = pf * 2.0_dp * r ( 1 ) * r ( 2 )
    END SELECT
  CASE ( 3 )
    SELECT CASE ( m )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 3 and m value out of bounds")
    CASE ( 3 )
      pf = SQRT ( 35.0_dp / ( 32.0_dp * pi ) )
      y = pf * r ( 1 ) * ( r ( 1 )**2 - 3.0_dp * r ( 2 )**2 )
    CASE ( 2 )
      pf = SQRT ( 105.0_dp / ( 16.0_dp * pi ) )
      y = pf * r ( 3 ) * ( r ( 1 )**2 - r ( 2 )**2 )
    CASE ( 1 )
      pf = SQRT ( 21.0_dp / ( 32.0_dp * pi ) )
      y = pf * r ( 1 ) * ( 5.0_dp * r ( 3 ) * r ( 3 ) - 1.0_dp )
    CASE ( 0 )
      pf = SQRT ( 7.0_dp / ( 16.0_dp * pi ) )
      y = pf * r ( 3 ) * ( 5.0_dp * r ( 3 ) * r ( 3 ) - 3.0_dp )
    CASE ( -1 )
      pf = SQRT ( 21.0_dp / ( 32.0_dp * pi ) )
      y = pf * r ( 2 ) * ( 5.0_dp * r ( 3 ) * r ( 3 ) - 1.0_dp )
    CASE ( -2 )
      pf = SQRT ( 105.0_dp / ( 16.0_dp * pi ) )
      y = pf * 2.0_dp * r ( 1 ) * r ( 2 ) * r ( 3 )
    CASE ( -3 )
      pf = SQRT ( 35.0_dp / ( 32.0_dp * pi ) )
      y = pf * r ( 2 ) * ( 3.0_dp * r ( 1 )**2 - r ( 2 )**2 )
    END SELECT
  CASE DEFAULT
    IF ( m < -l .OR. m > l ) CALL stop_program(routineN,moduleN,__LINE__,&
                                               "m value out of bounds")
    lpm = fac ( l + ABS ( m ) )
    lmm = fac ( l - ABS ( m ) )
    IF ( m == 0 ) THEN
      t = 4.0_dp * pi
    ELSE
      t = 2.0_dp * pi
    END IF
    IF ( ABS ( lpm ) < EPSILON ( 1.0_dp ) ) THEN
      pf = REAL ( 2*l+1,KIND=dp)  / t
    ELSE
      pf = ( REAL ( 2*l+1,KIND=dp) * lmm ) / ( t * lpm )
    ENDIF
    pf = SQRT ( pf )
    z = r ( 3 )
    plm = legendre ( z, l, m )
    IF ( m == 0 ) THEN
      y = pf * plm
    ELSE
      rxy = SQRT ( r ( 1 ) ** 2  +  r ( 2 ) ** 2 )
      IF ( rxy < EPSILON ( 1.0_dp ) ) THEN
        y = 0.0_dp
      ELSE
        cp = r ( 1 ) / rxy
        sp = r ( 2 ) / rxy
        IF ( m > 0 ) THEN
          y = pf * plm * cosn ( m, cp, sp )
        ELSE
          y = pf * plm * sinn ( -m , cp, sp )
        END IF
      END IF
    END IF
  END SELECT

  END SUBROUTINE rry_lm

! *****************************************************************************
!> \brief ...
!> \param r ...
!> \param y ...
!> \param l ...
!> \param m ...
! *****************************************************************************
  SUBROUTINE cvy_lm ( r, y, l, m )
!
! Complex Spherical Harmonics
!                   _                _
!                  | [(2l+1)(l-|m|)!] |1/2 m
!  Y_lm ( t, p ) = |------------------|   P_l(cos(t)) Exp[ i m p ]
!                  |  [4Pi(l+|m|)!]|  |
!
! Input: r == (x,y,z) : normalised    x^2 + y^2 + z^2 = 1
!
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r
    COMPLEX(KIND=dp), DIMENSION(:), &
      INTENT(OUT)                            :: y
    INTEGER, INTENT(IN)                      :: l, m

    CHARACTER(len=*), PARAMETER :: routineN = 'cvy_lm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, n
    REAL(KIND=dp)                            :: cp, lmm, lpm, pf, plm, rxy, &
                                                sp, t, ym, yp, z

  n = SIZE ( r ,2 )
  SELECT CASE ( l )
  CASE ( : -1 )
    CALL stop_program(routineN,moduleN,__LINE__,"Negative l value")
  CASE ( 0 )
    pf = SQRT ( 1.0_dp / ( 4.0_dp * pi ) )
    IF ( m /= 0 )  CALL stop_program(routineN,moduleN,__LINE__,&
                                     "l = 0 and m value out of bounds")
    y(:) = pf
  CASE ( 1 )
    SELECT CASE ( m )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 1 and m value out of bounds")
    CASE ( 1 )
      pf = SQRT ( 3.0_dp / ( 8.0_dp * pi ) )
      DO i = 1, n
        yp = r ( 1, i )
        ym = r ( 2, i )
        y ( i ) = pf * CMPLX ( yp, ym,KIND=dp)
      END DO
    CASE ( 0 )
      pf = SQRT ( 3.0_dp / ( 4.0_dp * pi ) )
      y ( : ) = pf * r ( 3, : )
    CASE ( -1 )
      pf = SQRT ( 3.0_dp / ( 8.0_dp * pi ) )
      DO i = 1, n
        yp = r ( 1, i )
        ym = r ( 2, i )
        y ( i ) = pf * CMPLX ( yp, -ym,KIND=dp)
      END DO
    END SELECT
  CASE ( 2 )
    SELECT CASE ( m )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 2 and m value out of bounds")
    CASE ( 2 )
      pf = SQRT ( 15.0_dp / ( 32.0_dp * pi ) )
      DO i = 1, n
        yp = ( r ( 1, i ) * r ( 1, i ) - r ( 2, i ) * r ( 2, i ) )
        ym = 2.0_dp * r ( 1, i ) * r ( 2, i )
        y ( i ) = pf * CMPLX ( yp, ym,KIND=dp)
      END DO
    CASE ( 1 )
      pf = SQRT ( 15.0_dp / ( 8.0_dp * pi ) )
      DO i = 1, n
        yp = r ( 3, i ) * r ( 1, i )
        ym = r ( 3, i ) * r ( 2, i )
        y ( i ) = pf * CMPLX ( yp, ym,KIND=dp)
      END DO
    CASE ( 0 )
      pf = SQRT ( 5.0_dp / ( 16.0_dp * pi ) )
      y ( : ) = pf * ( 3.0_dp * r ( 3, : ) * r ( 3, : ) - 1.0_dp )
    CASE ( -1 )
      pf = SQRT ( 15.0_dp / ( 8.0_dp * pi ) )
      DO i = 1, n
        yp = r ( 3, i ) * r ( 1, i )
        ym = r ( 3, i ) * r ( 2, i )
        y ( i ) = pf * CMPLX ( yp, -ym,KIND=dp)
      END DO
    CASE ( -2 )
      pf = SQRT ( 15.0_dp / ( 32.0_dp * pi ) )
      DO i = 1, n
        yp = ( r ( 1, i ) * r ( 1, i ) - r ( 2, i ) * r ( 2, i ) )
        ym = 2.0_dp * r ( 1, i ) * r ( 2, i )
        y ( i ) = pf * CMPLX ( yp, -ym,KIND=dp)
      END DO
    END SELECT
  CASE ( 3 )
    SELECT CASE ( m )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 3 and m value out of bounds")
    CASE ( 3 )
      pf = SQRT ( 35.0_dp / ( 64.0_dp * pi ) )
      DO i = 1, n
        yp = r ( 1, i ) * ( r ( 1, i )**2 - 3.0_dp * r ( 2, i )**2 )
        ym = r ( 2, i ) * ( 3.0_dp * r ( 1, i )**2 - r ( 2, i )**2 )
        y ( i ) = pf * CMPLX ( yp, ym,KIND=dp)
      END DO
    CASE ( 2 )
      pf = SQRT ( 105.0_dp / ( 32.0_dp * pi ) )
      DO i = 1, n
        yp = r ( 3, i ) * ( r ( 1, i )**2 - r ( 2, i )**2 )
        ym = 2.0_dp * r ( 1, i ) * r ( 2, i ) * r ( 3, i )
        y ( i ) = pf * CMPLX ( yp, ym,KIND=dp)
      END DO
    CASE ( 1 )
      pf = SQRT ( 21.0_dp / ( 64.0_dp * pi ) )
      DO i = 1, n
        yp = r ( 1, i ) * ( 5.0_dp * r ( 3, i ) * r ( 3, i ) - 1.0_dp )
        ym = r ( 2, i ) * ( 5.0_dp * r ( 3, i ) * r ( 3, i ) - 1.0_dp )
        y ( i ) = pf * CMPLX ( yp, ym,KIND=dp)
      END DO
    CASE ( 0 )
      pf = SQRT ( 7.0_dp / ( 16.0_dp * pi ) )
      y ( : ) = pf * r ( 3, : ) * ( 5.0_dp * r ( 3, : ) * r ( 3, : ) - 3.0_dp )
    CASE ( -1 )
      pf = SQRT ( 21.0_dp / ( 64.0_dp * pi ) )
      DO i = 1, n
        yp = r ( 1, i ) * ( 5.0_dp * r ( 3, i ) * r ( 3, i ) - 1.0_dp )
        ym = r ( 2, i ) * ( 5.0_dp * r ( 3, i ) * r ( 3, i ) - 1.0_dp )
        y ( i ) = pf * CMPLX ( yp, -ym,KIND=dp)
      END DO
    CASE ( -2 )
      pf = SQRT ( 105.0_dp / ( 32.0_dp * pi ) )
      DO i = 1, n
        yp = r ( 3, i ) * ( r ( 1, i )**2 - r ( 2, i )**2 )
        ym = 2.0_dp * r ( 1, i ) * r ( 2, i ) * r ( 3, i )
        y ( i ) = pf * CMPLX ( yp, -ym,KIND=dp)
      END DO
    CASE ( -3 )
      pf = SQRT ( 35.0_dp / ( 64.0_dp * pi ) )
      DO i = 1, n
        yp = r ( 1, i ) * ( r ( 1, i )**2 - 3.0_dp * r ( 2, i )**2 )
        ym = r ( 2, i ) * ( 3.0_dp * r ( 1, i )**2 - r ( 2, i )**2 )
        y ( i ) = pf * CMPLX ( yp, -ym,KIND=dp)
      END DO
    END SELECT
  CASE DEFAULT
    IF ( m < -l .OR. m > l ) CALL stop_program(routineN,moduleN,__LINE__,&
                                               "m value out of bounds")
    lpm = fac ( l + ABS ( m ) )
    lmm = fac ( l - ABS ( m ) )
    t = 4.0_dp * pi
    IF ( ABS ( lpm ) < EPSILON ( 1.0_dp ) ) THEN
      pf = REAL ( 2*l+1,KIND=dp)  / t
    ELSE
      pf = ( REAL ( 2*l+1,KIND=dp) * lmm ) / ( t * lpm )
    ENDIF
    pf = SQRT ( pf )
    DO i = 1, n
      z = r ( 3, i )
      plm = legendre ( z, l, m )
      IF ( m == 0 ) THEN
        y ( i ) = pf * plm
      ELSE
        rxy = SQRT ( r ( 1, i ) ** 2  +  r ( 2, i ) ** 2 )
        IF ( rxy < EPSILON ( 1.0_dp ) ) THEN
          y ( i ) = 0.0_dp
        ELSE
          cp = r ( 1, i ) / rxy
          sp = r ( 2, i ) / rxy
          IF ( m > 0 ) THEN
            yp = cosn ( m, cp, sp )
            ym = sinn ( m, cp, sp )
            y ( i ) = pf * plm * CMPLX( yp, ym,KIND=dp)
          ELSE
            yp = cosn ( -m, cp, sp )
            ym = sinn ( -m, cp, sp )
            y ( i ) = pf * plm * CMPLX( yp, -ym,KIND=dp)
          END IF
        END IF
      END IF
    END DO
  END SELECT

  END SUBROUTINE cvy_lm

! *****************************************************************************
!> \brief ...
!> \param r ...
!> \param y ...
!> \param l ...
!> \param m ...
! *****************************************************************************
  SUBROUTINE ccy_lm ( r, y, l, m )
!
! Complex Spherical Harmonics
!                   _                _
!                  | [(2l+1)(l-|m|)!] |1/2 m
!  Y_lm ( t, p ) = |------------------|   P_l(cos(t)) Exp[ i m p ]
!                  |  [4Pi(l+|m|)!]|  |
!
! Input: r == (x,y,z) : normalised    x^2 + y^2 + z^2 = 1
!
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: r
    COMPLEX(KIND=dp), INTENT(OUT)            :: y
    INTEGER, INTENT(IN)                      :: l, m

    CHARACTER(len=*), PARAMETER :: routineN = 'ccy_lm', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: cp, lmm, lpm, pf, plm, rxy, &
                                                sp, t, ym, yp, z

  SELECT CASE ( l )
  CASE ( : -1 )
    CALL stop_program(routineN,moduleN,__LINE__,"Negative l value")
  CASE ( 0 )
    pf = SQRT ( 1.0_dp / ( 4.0_dp * pi ) )
    IF ( m /= 0 )  CALL stop_program(routineN,moduleN,__LINE__,&
                                     "l = 0 and m value out of bounds")
    y = pf
  CASE ( 1 )
    SELECT CASE ( m )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 1 and m value out of bounds")
    CASE ( 1 )
      pf = SQRT ( 3.0_dp / ( 8.0_dp * pi ) )
      yp = r ( 1 )
      ym = r ( 2 )
      y = pf * CMPLX ( yp, ym,KIND=dp)
    CASE ( 0 )
      pf = SQRT ( 3.0_dp / ( 4.0_dp * pi ) )
      y = pf * r ( 3 )
    CASE ( -1 )
      pf = SQRT ( 3.0_dp / ( 8.0_dp * pi ) )
      yp = r ( 1 )
      ym = r ( 2 )
      y = pf * CMPLX ( yp, -ym,KIND=dp)
    END SELECT
  CASE ( 2 )
    SELECT CASE ( m )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 2 and m value out of bounds")
    CASE ( 2 )
      pf = SQRT ( 15.0_dp / ( 32.0_dp * pi ) )
      yp = ( r ( 1 ) * r ( 1 ) - r ( 2 ) * r ( 2 ) )
      ym = 2.0_dp * r ( 1 ) * r ( 2 )
      y = pf * CMPLX ( yp, ym,KIND=dp)
    CASE ( 1 )
      pf = SQRT ( 15.0_dp / ( 8.0_dp * pi ) )
      yp = r ( 3 ) * r ( 1 )
      ym = r ( 3 ) * r ( 2 )
      y = pf * CMPLX ( yp, ym,KIND=dp)
    CASE ( 0 )
      pf = SQRT ( 5.0_dp / ( 16.0_dp * pi ) )
      y = pf * ( 3.0_dp * r ( 3 ) * r ( 3 ) - 1.0_dp )
    CASE ( -1 )
      pf = SQRT ( 15.0_dp / ( 8.0_dp * pi ) )
      yp = r ( 3 ) * r ( 1 )
      ym = r ( 3 ) * r ( 2 )
      y = pf * CMPLX ( yp, -ym,KIND=dp)
    CASE ( -2 )
      pf = SQRT ( 15.0_dp / ( 32.0_dp * pi ) )
      yp = ( r ( 1 ) * r ( 1 ) - r ( 2 ) * r ( 2 ) )
      ym = 2.0_dp * r ( 1 ) * r ( 2 )
      y = pf * CMPLX ( yp, -ym,KIND=dp)
    END SELECT
  CASE ( 3 )
    SELECT CASE ( m )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 3 and m value out of bounds")
    CASE ( 3 )
      pf = SQRT ( 35.0_dp / ( 64.0_dp * pi ) )
      yp = r ( 1 ) * ( r ( 1 )**2 - 3.0_dp * r ( 2 )**2 )
      ym = r ( 2 ) * ( 3.0_dp * r ( 1 )**2 - r ( 2 )**2 )
      y = pf * CMPLX ( yp, ym,KIND=dp)
    CASE ( 2 )
      pf = SQRT ( 105.0_dp / ( 32.0_dp * pi ) )
      yp = r ( 3 ) * ( r ( 1 )**2 - r ( 2 )**2 )
      ym = 2.0_dp * r ( 1 ) * r ( 2 ) * r ( 3 )
      y = pf * CMPLX ( yp, ym,KIND=dp)
    CASE ( 1 )
      pf = SQRT ( 21.0_dp / ( 64.0_dp * pi ) )
      yp = r ( 1 ) * ( 5.0_dp * r ( 3 ) * r ( 3 ) - 1.0_dp )
      ym = r ( 2 ) * ( 5.0_dp * r ( 3 ) * r ( 3 ) - 1.0_dp )
      y = pf * CMPLX ( yp, ym,KIND=dp)
    CASE ( 0 )
      pf = SQRT ( 7.0_dp / ( 16.0_dp * pi ) )
      y = pf * r ( 3 ) * ( 5.0_dp * r ( 3 ) * r ( 3 ) - 3.0_dp )
    CASE ( -1 )
      pf = SQRT ( 21.0_dp / ( 64.0_dp * pi ) )
      yp = r ( 1 ) * ( 5.0_dp * r ( 3 ) * r ( 3 ) - 1.0_dp )
      ym = r ( 2 ) * ( 5.0_dp * r ( 3 ) * r ( 3 ) - 1.0_dp )
      y = pf * CMPLX ( yp, -ym,KIND=dp)
    CASE ( -2 )
      pf = SQRT ( 105.0_dp / ( 32.0_dp * pi ) )
      yp = r ( 3 ) * ( r ( 1 )**2 - r ( 2 )**2 )
      ym = 2.0_dp * r ( 1 ) * r ( 2 ) * r ( 3 )
      y = pf * CMPLX ( yp, -ym,KIND=dp)
    CASE ( -3 )
      pf = SQRT ( 35.0_dp / ( 64.0_dp * pi ) )
      yp = r ( 1 ) * ( r ( 1 )**2 - 3.0_dp * r ( 2 )**2 )
      ym = r ( 2 ) * ( 3.0_dp * r ( 1 )**2 - r ( 2 )**2 )
      y = pf * CMPLX ( yp, -ym,KIND=dp)
    END SELECT
  CASE DEFAULT
    IF ( m < -l .OR. m > l ) CALL stop_program(routineN,moduleN,__LINE__,&
                                               "m value out of bounds")
    lpm = fac ( l + ABS ( m ) )
    lmm = fac ( l - ABS ( m ) )
    t = 4.0_dp * pi
    IF ( ABS ( lpm ) < EPSILON ( 1.0_dp ) ) THEN
      pf = REAL ( 2*l+1,KIND=dp)  / t
    ELSE
      pf = ( REAL ( 2*l+1,KIND=dp) * lmm ) / ( t * lpm )
    ENDIF
    pf = SQRT ( pf )
    z = r ( 3 )
    plm = legendre ( z, l, m )
    IF ( m == 0 ) THEN
      y = pf * plm
    ELSE
      rxy = SQRT ( r ( 1 ) ** 2  +  r ( 2 ) ** 2 )
      IF ( rxy < EPSILON ( 1.0_dp ) ) THEN
        y = 0.0_dp
      ELSE
        cp = r ( 1 ) / rxy
        sp = r ( 2 ) / rxy
        IF ( m > 0 ) THEN
          yp = cosn ( m, cp, sp )
          ym = sinn ( m , cp, sp )
          y = pf * plm * CMPLX( yp, ym,KIND=dp)
        ELSE
          yp = cosn ( -m, cp, sp )
          ym = sinn ( -m , cp, sp )
          y = pf * plm * CMPLX( yp, -ym,KIND=dp)
        END IF
      END IF
    END IF
  END SELECT

  END SUBROUTINE ccy_lm

! *****************************************************************************
!> \brief ...
!> \param x ...
!> \param l ...
!> \param m ...
!> \retval plm ...
! *****************************************************************************
  FUNCTION legendre ( x, l, m ) RESULT ( plm )

    REAL(KIND=dp), INTENT(IN)                :: x
    INTEGER, INTENT(IN)                      :: l, m
    REAL(KIND=dp)                            :: plm

    CHARACTER(len=*), PARAMETER :: routineN = 'legendre', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: il, im, mm
    REAL(KIND=dp)                            :: fact, pll, pmm, pmmp1, somx2

  IF ( ABS ( x ) > 1.0_dp ) CALL stop_program(routineN,moduleN,__LINE__,"x value > 1")
  SELECT CASE ( l )
  CASE ( : -1 )
    CALL stop_program(routineN,moduleN,__LINE__,"Negative l value")
  CASE ( 0 )
    plm = 1.0_dp
  CASE ( 1 )
    SELECT CASE ( ABS ( m ) )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 1 and m value out of bounds")
    CASE ( 1 )
      plm = SQRT ( 1.0_dp - x * x )
    CASE ( 0 )
      plm = x
    END SELECT
  CASE ( 2 )
    SELECT CASE ( ABS ( m ) )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 2 and m value out of bounds")
    CASE ( 2 )
      plm = 3.0_dp * ( 1.0_dp - x * x )
    CASE ( 1 )
      plm = 3.0_dp * x * SQRT ( 1.0_dp - x * x )
    CASE ( 0 )
      plm = 1.5_dp * x * x - 0.5_dp
    END SELECT
  CASE ( 3 )
    SELECT CASE ( ABS ( m ) )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 3 and m value out of bounds")
    CASE ( 3 )
      plm = 15.0_dp * ( 1.0_dp - x * x ) ** 1.5_dp
    CASE ( 2 )
      plm = 15.0_dp * x * ( 1.0_dp - x * x )
    CASE ( 1 )
      plm = ( 7.5_dp * x * x - 1.5_dp ) * SQRT ( 1.0_dp - x * x )
    CASE ( 0 )
      plm = 2.5_dp * x**3 - 1.5_dp * x
    END SELECT
  CASE ( 4 )
    SELECT CASE ( ABS ( m ) )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 4 and m value out of bounds")
    CASE ( 4 )
      plm = 105.0_dp * ( 1.0_dp - x * x ) ** 2
    CASE ( 3 )
      plm = 105.0_dp * x * ( 1.0_dp - x * x ) ** 1.5_dp
    CASE ( 2 )
      plm = ( 52.5_dp * x * x - 7.5_dp ) * ( 1.0_dp - x * x )
    CASE ( 1 )
      plm = ( 17.5_dp * x ** 3 - 7.5_dp * x ) * SQRT ( 1.0_dp - x * x )
    CASE ( 0 )
      plm = 4.375_dp * x ** 4 - 3.75_dp * x ** 2 + 0.375_dp
    END SELECT
  CASE ( 5 )
    SELECT CASE ( ABS ( m ) )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 5 and m value out of bounds")
    CASE ( 5 )
      plm = 945.0_dp * ( 1.0_dp - x * x ) ** 2.5_dp
    CASE ( 4 )
      plm = 945.0_dp * x * ( 1.0_dp - x * x ) ** 2
    CASE ( 3 )
      plm = -( -472.5_dp * x * x + 52.5_dp ) * ( 1.0_dp - x * x ) ** 1.5_dp
    CASE ( 2 )
      plm = ( 157.5_dp * x ** 3 - 52.5_dp * x ) * ( 1.0_dp - x * x )
    CASE ( 1 )
      plm = -( -39.375_dp * x ** 4 + 26.25_dp * x * x - &
              1.875_dp ) * SQRT ( 1.0_dp - x * x )
    CASE ( 0 )
      plm = 7.875_dp * x ** 5 - 8.75_dp * x ** 3 + 1.875_dp * x
    END SELECT
  CASE ( 6 )
    SELECT CASE ( ABS ( m ) )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 6 and m value out of bounds")
    CASE ( 6 )
      plm = 10395.0_dp * ( 1.0_dp - x * x ) ** 3
    CASE ( 5 )
      plm = 10395.0_dp * x * ( 1.0_dp - x * x ) ** 2.5_dp
    CASE ( 4 )
      plm = ( 5197.5_dp * x * x - 472.5_dp ) * ( 1.0_dp - x * x ) ** 2
    CASE ( 3 )
      plm = -( - 1732.5_dp * x ** 3 + 472.5_dp * x ) * &
            ( 1.0_dp - x * x ) ** 1.5_dp
    CASE ( 2 )
      plm = ( 433.125_dp * x ** 4 - 236.25_dp * x ** 2 + &
              13.125_dp ) * ( 1.0_dp - x * x )
    CASE ( 1 )
      plm = -( -86.625_dp * x ** 5 + 78.75_dp * x ** 3 - &
              13.125_dp * x ) * SQRT ( 1.0_dp - x * x )
    CASE ( 0 )
      plm = 14.4375_dp * x ** 6 - 19.6875_dp * x ** 4 + &
            6.5625_dp * x ** 2 - 0.3125_dp
    END SELECT
  CASE DEFAULT
    mm = ABS ( m )
    IF ( mm > l ) CALL stop_program(routineN,moduleN,__LINE__,"m out of bounds")
! use recurence from numerical recipies
    pmm = 1.0_dp
    IF ( mm > 0 ) THEN
      somx2 = SQRT((1.0_dp-x)*(1.0_dp+x))
      fact = 1.0_dp
      DO im = 1, mm
        pmm = pmm * fact * somx2
        fact = fact + 2.0_dp
      END DO
    END IF
    IF ( l == mm ) THEN
      plm = pmm
    ELSE
      pmmp1 = x * REAL ( 2 * mm + 1,KIND=dp) * pmm
      IF ( l == mm + 1 ) THEN
        plm = pmmp1
      ELSE
        DO il = mm + 2, l
          pll = ( x * REAL ( 2 * il - 1,KIND=dp) * pmmp1 - &
                      REAL ( il + mm - 1,KIND=dp) * pmm ) / REAL ( il - mm,KIND=dp)
          pmm = pmmp1
          pmmp1 = pll
        END DO
        plm = pll
      END IF
    END IF
  END SELECT

  END FUNCTION legendre

! *****************************************************************************
!> \brief ...
!> \param x ...
!> \param l ...
!> \param m ...
!> \retval dplm ...
! *****************************************************************************
  FUNCTION dlegendre (x, l, m) RESULT (dplm)
    REAL(KIND=dp), INTENT(IN)                :: x
    INTEGER, INTENT(IN)                      :: l, m
    REAL(KIND=dp)                            :: dplm

    CHARACTER(len=*), PARAMETER :: routineN = 'dlegendre', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: il, im, mm
    REAL(KIND=dp)                            :: dpll, dpmm, dpmmp1, dsomx2, &
                                                fact, pll, plm, pmm, pmmp1, &
                                                somx2

  IF ( ABS ( x ) > 1.0_dp ) CALL stop_program(routineN,moduleN,__LINE__,"x value > 1")
  SELECT CASE ( l )
  CASE ( 0 )
    dplm = 0.0_dp
  CASE ( 1 )
    SELECT CASE ( ABS ( m ) )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 1 and m value out of bounds")
    CASE ( 1 )
      dplm = -x / SQRT ( 1.0_dp - x * x )
    CASE ( 0 )
      dplm = 1.0_dp
    END SELECT
  CASE ( 2 )
    SELECT CASE ( ABS ( m ) )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 2 and m value out of bounds")
    CASE ( 2 )
      dplm = -6.0_dp * x
    CASE ( 1 )
      dplm = 3.0_dp * SQRT ( 1.0_dp - x * x ) - 3.0_dp * x * x / SQRT ( 1.0_dp - x * x )
    CASE ( 0 )
      dplm = 3.0_dp * x
    END SELECT
  CASE ( 3 )
    SELECT CASE ( ABS ( m ) )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 3 and m value out of bounds")
    CASE ( 3 )
      dplm = -45.0_dp * SQRT( 1.0_dp - x * x ) * x
    CASE ( 2 )
      dplm = 15.0_dp * ( 1.0_dp - x * x ) - 30.0_dp * x * x
    CASE ( 1 )
      dplm = 15.0_dp * x * SQRT ( 1.0_dp - x * x ) - (x *  ( 7.5_dp * x * x - 1.5_dp ))/ SQRT ( 1.0_dp - x * x )
    CASE ( 0 )
      dplm = 7.5_dp * x * x - 1.5_dp
    END SELECT
  CASE ( 4 )
    SELECT CASE ( ABS ( m ) )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 4 and m value out of bounds")
    CASE ( 4 )
      dplm = - 420 * x * (1 - x*x)
    CASE ( 3 )
      dplm =  105.0_dp * ( ( 1.0_dp - x * x ) ** 1.5_dp - 3.0_dp * x * x * ( 1.0_dp - x * x ) ** 0.5_dp)
    CASE ( 2 )
      dplm = 105.0_dp * x * ( 1.0_dp - x * x ) -2.0_dp * x * ( 52.5_dp * x * x - 7.5_dp )
    CASE ( 1 )
      IF( ABS(x) - 1.0_dp < EPSILON(1.0_dp)) THEN
        dplm = 0.0_dp
      ELSE
        dplm = ( 17.5_dp * 3.0_dp *  x ** 2 - 7.5_dp ) * SQRT ( 1.0_dp - x * x ) - &
               x * ( 17.5_dp * x ** 3 - 7.5_dp * x ) / SQRT ( 1.0_dp - x * x )
      END IF
    CASE ( 0 )
      dplm = 4.375_dp * 4.0_dp *  x ** 3 - 3.75_dp * 2.0_dp * x
    END SELECT
  CASE ( 5 )
    SELECT CASE ( ABS ( m ) )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 5 and m value out of bounds")
    CASE ( 5 )
      dplm = -945.0_dp * 5.0_dp * x * ( 1.0_dp - x * x ) ** 1.5_dp
    CASE ( 4 )
      dplm = 945.0_dp * (( 1.0_dp - x * x ) ** 2 - 2.0_dp * x * x * ( 1.0_dp - x * x ))
    CASE ( 3 )
      dplm = 945.0_dp * x * ( 1.0_dp - x * x ) ** 1.5_dp - &
               3.0_dp * x * ( 472.5_dp * x * x - 52.5_dp ) * ( 1.0_dp - x * x ) ** 0.5_dp
    CASE ( 2 )
      dplm = (3.0_dp *157.5_dp * x ** 2 - 52.5_dp ) * ( 1.0_dp - x * x ) - &
                 ( 157.5_dp * x ** 3 - 52.5_dp * x ) * ( - 2.0_dp * x )
    CASE ( 1 )
      IF( ABS(x) - 1.0_dp < EPSILON(1.0_dp)) THEN
        dplm = 0.0_dp
      ELSE
        dplm = - (-39.375_dp * 4.0_dp * x * x * x + 2.0_dp * 26.25_dp*x ) * SQRT ( 1.0_dp - x * x )  +&
            x * ( -39.375_dp * x **4 + 26.25_dp * x * x - 1.875_dp ) / SQRT ( 1.0_dp - x * x )
      END IF
    CASE ( 0 )
      dplm = 5.0_dp * 7.875_dp * x ** 4 - 3.0_dp * 8.75_dp * x ** 2 + 1.875_dp
    END SELECT
  CASE ( 6 )
    SELECT CASE ( ABS ( m ) )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 6 and m value out of bounds")
    CASE ( 6 )
      dplm = -10395.0_dp * 6.0_dp * x * ( 1.0_dp - x * x ) ** 2
    CASE ( 5 )
      dplm = 10395.0_dp * ( ( 1.0_dp - x * x ) ** 2.5_dp - 5.0_dp * x * x * ( 1.0_dp - x * x ) ** 1.5_dp )
    CASE ( 4 )
      dplm = 2.0_dp * 5197.5_dp * x * ( 1.0_dp - x * x ) ** 2 - &
          4.0_dp * x * ( 5197.5_dp * x * x - 472.5_dp ) * ( 1.0_dp - x * x )
    CASE ( 3 )
      dplm = - (-3.0_dp * 1732.5_dp * x * x + 472.5_dp ) * ( 1.0_dp - x * x ) ** 1.5_dp + &
               ( - 1732.5_dp * x ** 3 + 472.5_dp * x ) * 3.0_dp * x * ( 1.0_dp - x * x ) ** 0.5_dp
    CASE ( 2 )
      dplm = (433.125_dp * 4.0_dp * x ** 3 -2.0_dp * 236.25_dp *x) * ( 1.0_dp - x * x ) - &
              2.0_dp * x * ( 433.125_dp * x ** 4 - 236.25_dp * x ** 2 + 13.125_dp )
    CASE ( 1 )
      IF( ABS(x) - 1.0_dp < EPSILON(1.0_dp)) THEN
        dplm = 0.0_dp
      ELSE
        dplm = -(-5.0_dp * 86.625_dp * x ** 4 + 3.0_dp *78.75_dp ** 2 - 13.125_dp ) * SQRT(1.0_dp-x*x) + &
            x * ( -86.625_dp * x ** 5 + 78.75_dp * x ** 3 - 13.125_dp * x ) / SQRT ( 1.0_dp - x * x )
      END IF
    CASE ( 0 )
      dplm = 14.4375_dp * 6.0_dp * x ** 5 - 19.6875_dp * 4.0_dp * x ** 3 + &
            6.5625_dp * 2.0_dp * x
    END SELECT
  CASE DEFAULT
!    CALL stop_program(routineN,moduleN,__LINE__,"l > 6 dplm not implemented")
    mm = ABS ( m )
    IF ( mm > l ) CALL stop_program(routineN,moduleN,__LINE__,"m out of bounds")
! use recurence from numerical recipies
    pmm = 1.0_dp
    dpmm = 0.0_dp
    IF ( mm > 0 ) THEN
      somx2 = SQRT((1.0_dp-x)*(1.0_dp+x))
      dsomx2 = ((1e0_dp - x)*(1e0_dp + x))**(0.5e0_dp*m-1)
      dpmm = -m*x*dsomx2*dfac(2*m-1)
      fact = 1.0_dp
      DO im = 1, mm
        pmm = -pmm * fact * somx2
        dpmm = -dpmm
        fact = fact + 2.0_dp
      END DO
    END IF
    IF ( l == mm ) THEN
      plm = pmm
      dplm = dpmm
    ELSE
      pmmp1 = x * REAL ( 2 * mm + 1,KIND=dp) * pmm
      dpmmp1 = (2*mm+1)*pmm + x*(2*mm+1)*dpmm
      IF ( l == mm + 1 ) THEN
        plm = pmmp1
        dplm = dpmmp1
      ELSE
        DO il = mm + 2, l
          pll = ( x * REAL ( 2 * il - 1,KIND=dp) * pmmp1 - &
                      REAL ( il + mm - 1,KIND=dp) * pmm ) / REAL ( il - mm,KIND=dp)
          dpll = ((2*il-1)*pmmp1 + x*(2*il-1)*dpmmp1 - (il+mm-1)*dpmm)/(il-mm)
          pmm = pmmp1
          pmmp1 = pll
          dpmm = dpmmp1
          dpmmp1 = dpll
        END DO
        plm = pll
        dplm = dpll
      END IF
    END IF

  END SELECT

  END FUNCTION dlegendre

! *****************************************************************************
!> \brief ...
!> \param x ...
!> \param l ...
!> \retval dPof1 ...
! *****************************************************************************
  FUNCTION dPof1(x, l)

    REAL(KIND=dp), INTENT(IN)                :: x
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp)                            :: dPof1

    CHARACTER(len=*), PARAMETER :: routineN = 'dPof1', &
      routineP = moduleN//':'//routineN

    IF ( ABS ( x ) -1.0_dp > EPSILON(1.0_dp)) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"|x| is not 1")
    END IF
    IF (x > 0.0_dp) THEN
      SELECT CASE ( l )
      CASE ( : -1 )
        CALL stop_program(routineN,moduleN,__LINE__,"Negative l value")
      CASE ( 0 )
        dPof1 = 0.0_dp
      CASE ( 1 )
        dPof1 = 1.0_dp
      CASE ( 2 )
        dPof1 = 3.0_dp
      CASE ( 3 )
        dPof1 = 6.0_dp
      CASE ( 4 )
        dPof1 = 10.0_dp
      CASE ( 5 )
        dPof1 = 15.0_dp
      CASE ( 6 )
        dPof1 = 21.0_dp
      CASE ( 7 :  )
        CALL stop_program(routineN,moduleN,__LINE__,"Not implemented")
      END SELECT
    ELSE
      SELECT CASE ( l )
      CASE ( : -1 )
        CALL stop_program(routineN,moduleN,__LINE__,"Negative l value")
      CASE ( 0 )
        dPof1 = 0.0_dp
      CASE ( 1 )
        dPof1 = 1.0_dp
      CASE ( 2 )
        dPof1 = -3.0_dp
      CASE ( 3 )
        dPof1 = 6.0_dp
      CASE ( 4 )
        dPof1 = -10.0_dp
      CASE ( 5 )
        dPof1 = 15.0_dp
      CASE ( 6 )
        dPof1 = -21.0_dp
      CASE ( 7 :  )
        CALL stop_program(routineN,moduleN,__LINE__,"Not implemented")
      END SELECT
    END IF

  END FUNCTION dPof1

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param k ...
!> \retval choose ...
! *****************************************************************************
  FUNCTION choose ( n, k )

    INTEGER, INTENT(IN)                      :: n, k
    REAL(KIND=dp)                            :: choose

  IF ( n >= k ) THEN
    choose = REAL ( NINT ( fac ( n ) / ( fac ( k ) * fac ( n - k ) ) ),KIND=dp)
  ELSE
    choose = 0.0_dp
  ENDIF

  END FUNCTION choose

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param c ...
!> \param s ...
!> \retval cosn ...
! *****************************************************************************
  FUNCTION cosn ( n, c, s )

    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), INTENT(IN)                :: c, s
    REAL(KIND=dp)                            :: cosn

    INTEGER                                  :: i, j

  cosn = 0.0_dp
  IF( ABS(c) < EPSILON ( 1.0_dp ) .OR. n==0 ) THEN
    IF ( MOD(n,2) == 0 ) THEN
      cosn = (-1.0_dp)**INT(n/2)
    ELSE
      cosn = 0.0_dp
    END IF
  ELSE
    DO i = n, 0, -2
      IF (i==0) THEN
        j = n
        IF (j==0) THEN
          cosn = cosn + choose ( n, j )
        ELSE IF ( MOD(j/2,2) == 0 ) THEN
          cosn = cosn + choose ( n, j ) *  s**j
        ELSE
          cosn = cosn - choose ( n, j ) *  s**j
        END IF
      ELSE
        j = n - i
        IF (j==0) THEN
          cosn = cosn + choose ( n, j ) * c**i
        ELSE IF ( MOD(j/2,2) == 0 ) THEN
          cosn = cosn + choose ( n, j ) * c**i * s**j
        ELSE
          cosn = cosn - choose ( n, j ) * c**i * s**j
        END IF
      END IF
    END DO
  END IF

  END FUNCTION cosn

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param c ...
!> \param s ...
!> \retval dcosn_dcp ...
! *****************************************************************************
  FUNCTION dcosn_dcp ( n, c, s )

    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), INTENT(IN)                :: c, s
    REAL(KIND=dp)                            :: dcosn_dcp

    INTEGER                                  :: i, j

  dcosn_dcp = 0.0_dp

  IF(  s < EPSILON ( 1.0_dp ) ) THEN
    dcosn_dcp = 0.0_dp
  ELSE
    DO i = n, 0, -2
      IF (i == 0) THEN
        dcosn_dcp = dcosn_dcp
      ELSE
        j = n - i
        IF ( MOD(j/2,2) == 0 ) THEN
          dcosn_dcp = dcosn_dcp + choose ( n, j ) * REAL(i,dp) * c**(i-1) * s**j
        ELSE
          dcosn_dcp = dcosn_dcp - choose ( n, j ) * REAL(i,dp) * c**(i-1) * s**j
        END IF
      END IF
    END DO
  END IF

  END FUNCTION dcosn_dcp

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param c ...
!> \param s ...
!> \retval dcosn_dsp ...
! *****************************************************************************
  FUNCTION dcosn_dsp ( n, c, s )

    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), INTENT(IN)                :: c, s
    REAL(KIND=dp)                            :: dcosn_dsp

    INTEGER                                  :: i, j

  dcosn_dsp = 0.0_dp
  IF( c < EPSILON ( 1.0_dp )  .OR. s < EPSILON ( 1.0_dp ) ) THEN
    dcosn_dsp = 0.0_dp
  ELSE
    DO i = n, 0, -2
       j = n - i
      IF (j == 0) THEN
        dcosn_dsp = dcosn_dsp
      ELSE
        IF ( MOD(j/2,2) == 0 ) THEN
          dcosn_dsp = dcosn_dsp + choose ( n, j ) * REAL(j,dp) * s**(j-1) * c**i
        ELSE
          dcosn_dsp = dcosn_dsp - choose ( n, j ) * REAL(j,dp) * s**(j-1) * c**i
        END IF
      END IF
    END DO
  END IF

  END FUNCTION dcosn_dsp

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param c ...
!> \param s ...
!> \retval sinn ...
! *****************************************************************************
  FUNCTION sinn ( n, c, s )

    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), INTENT(IN)                :: c, s
    REAL(KIND=dp)                            :: sinn

    INTEGER                                  :: i, j

  sinn = 0.0_dp

  IF( ABS(s) < EPSILON ( 1.0_dp ) .OR. n==0) THEN
    sinn = 0.0_dp
  ELSE IF( ABS(c) < EPSILON ( 1.0_dp ) ) THEN
    IF ( MOD(n,2) == 0 ) THEN
      sinn = 0.0_dp
    ELSE
      sinn = s*(-1.0_dp)**INT((n-1)/2)
    END IF
  ELSE
    DO i = n-1, 0, -2
      IF (i==0) THEN
        j = n
        IF ( MOD(j/2,2) == 0 ) THEN
          sinn = sinn + choose ( n, j ) *  s**j
        ELSE
          sinn = sinn - choose ( n, j ) *  s**j
        END IF
      ELSE
        j = n - i
        IF ( MOD(j/2,2) == 0 ) THEN
          sinn = sinn + choose ( n, j ) * c**i * s**j
        ELSE
          sinn = sinn - choose ( n, j ) * c**i * s**j
        END IF
      END IF
    END DO
  END IF

  END FUNCTION sinn

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param c ...
!> \param s ...
!> \retval dsinn_dcp ...
! *****************************************************************************
  FUNCTION dsinn_dcp ( n, c, s )

    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), INTENT(IN)                :: c, s
    REAL(KIND=dp)                            :: dsinn_dcp

    INTEGER                                  :: i, j

  dsinn_dcp = 0.0_dp

  IF( c < EPSILON ( 1.0_dp ) .OR. s < EPSILON ( 1.0_dp ) ) THEN
    dsinn_dcp = 0.0_dp
  ELSE
    DO i = n-1, 0, -2
      IF(i == 0) THEN
        dsinn_dcp = dsinn_dcp
      ELSE
        j = n - i
        IF ( MOD(j/2,2) == 0 ) THEN
          dsinn_dcp = dsinn_dcp + choose ( n, j ) * REAL(i,dp) * c**(i-1) * s**j
        ELSE
          dsinn_dcp = dsinn_dcp - choose ( n, j ) * REAL(i,dp) * c**(i-1) * s**j
        END IF
      END IF
    END DO
  END IF

  END FUNCTION dsinn_dcp

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param c ...
!> \param s ...
!> \retval dsinn_dsp ...
! *****************************************************************************
  FUNCTION dsinn_dsp ( n, c, s )

    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), INTENT(IN)                :: c, s
    REAL(KIND=dp)                            :: dsinn_dsp

    INTEGER                                  :: i, j

  dsinn_dsp = 0.0_dp

  IF( c < EPSILON ( 1.0_dp ) .OR. s < EPSILON ( 1.0_dp ) ) THEN
    dsinn_dsp = 0.0_dp
  ELSE
    DO i = n-1, 0, -2
      j = n - i
      IF (j == 0) THEN
        dsinn_dsp = dsinn_dsp
      ELSE
        IF ( MOD(j/2,2) == 0 ) THEN
          dsinn_dsp = dsinn_dsp + choose ( n, j ) * c**i * REAL(j,dp) * s**(j-1)
        ELSE
          dsinn_dsp = dsinn_dsp - choose ( n, j ) * c**i * REAL(j,dp) * s**(j-1)
        END IF
      END IF
    END DO
  END IF

  END FUNCTION dsinn_dsp

! *****************************************************************************
!> \brief ...
!> \param r ...
!> \param dy ...
!> \param l ...
!> \param m ...
! *****************************************************************************
  SUBROUTINE drvy_lm ( r, dy, l, m )
!
! Real Spherical Harmonics
!                   _                   _
!                  |  [(2l+1)(l-|m|)!]   |1/2 m         cos(m p)   m>=0
!  Y_lm ( t, p ) = |---------------------|   P_l(cos(t))
!                  |[2Pi(1+d_m0)(l+|m|)!]|              sin(|m| p) m<0
!
! Input: r == (x,y,z) : normalised    x^2 + y^2 + z^2 = 1
!
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: dy
    INTEGER, INTENT(IN)                      :: l, m

    CHARACTER(len=*), PARAMETER :: routineN = 'drvy_lm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: cp, dcp_dx, dcp_dy, dplm, &
                                                dsp_dx, dsp_dy, lmm, lpm, pf, &
                                                plm, rxy, sp, t, z

  SELECT CASE ( l )
  CASE ( : -1 )
    CALL stop_program(routineN,moduleN,__LINE__,"Negative l value")
  CASE ( 0 )
    pf = SQRT ( 1.0_dp / ( 4.0_dp * pi ) )
    IF ( m /= 0 )  CALL stop_program(routineN,moduleN,__LINE__,&
                                     "l = 0 and m value out of bounds")
    dy(1:3,:) = 0.0_dp
  CASE ( 1 )
    SELECT CASE ( m )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 1 and m value out of bounds")
    CASE ( 1 )
      pf = SQRT ( 3.0_dp / ( 4.0_dp * pi ) )
      dy (1, : ) = pf
      dy (2, : ) = 0.0_dp
      dy (3, : ) = 0.0_dp
    CASE ( 0 )
      pf = SQRT ( 3.0_dp / ( 4.0_dp * pi ) )
      dy (1, : ) = 0.0_dp
      dy (2, : ) = 0.0_dp
      dy (3, : ) = pf
    CASE ( -1 )
      pf = SQRT ( 3.0_dp / ( 4.0_dp * pi ) )
      dy (1, : ) = 0.0_dp
      dy (2, : ) = pf
      dy (3, : ) = 0.0_dp
    END SELECT
  CASE ( 2 )
    SELECT CASE ( m )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 2 and m value out of bounds")
    CASE ( 2 )
      pf = SQRT ( 15.0_dp / ( 16.0_dp * pi ) )
      dy (1, : ) =  pf * 2.0_dp * r ( 1, : )
      dy (2, : ) = -pf * 2.0_dp * r ( 2, : )
      dy (3, : ) =  0.0_dp
    CASE ( 1 )
      pf = SQRT ( 15.0_dp / ( 4.0_dp * pi ) )
      dy (1, : ) = pf * r ( 3, : )
      dy (2, : ) = 0.0_dp
      dy (3, : ) = pf * r ( 1, : )
    CASE ( 0 )
      pf = SQRT ( 5.0_dp / ( 16.0_dp * pi ) )
      dy (1, : ) = 0.0_dp
      dy (2, : ) = 0.0_dp
      dy (3, : ) = pf * 6.0_dp * r ( 3, : )
    CASE ( -1 )
      pf = SQRT ( 15.0_dp / ( 4.0_dp * pi ) )
      dy (1, : ) = 0.0_dp
      dy (2, : ) = pf * r ( 3, : )
      dy (3, : ) = pf * r ( 2, : )
    CASE ( -2 )
      pf = SQRT ( 15.0_dp / ( 16.0_dp * pi ) )
      dy (1, : ) = pf * 2.0_dp * r ( 2, : )
      dy (2, : ) = pf * 2.0_dp * r ( 1, : )
      dy (3, : ) = 0.0_dp
    END SELECT
  CASE ( 3 )
    SELECT CASE ( m )
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,"l = 3 and m value out of bounds")
    CASE ( 3 )
      pf = SQRT ( 35.0_dp / ( 32.0_dp * pi ) )
      dy (1, : ) =  pf * 3.0_dp * ( r ( 1, : )*r ( 1, : ) - r ( 2, : )* r ( 2, : ))
      dy (2, : ) = -pf * 6.0_dp * r ( 1, : ) * r ( 2, : )
      dy (3, : ) =  0.0_dp
    CASE ( 2 )
      pf = SQRT ( 105.0_dp / ( 16.0_dp * pi ) )
      dy (1, : ) =  pf * 2.0_dp * r ( 3, : ) * r ( 1, : )
      dy (2, : ) = -pf * 2.0_dp * r ( 3, : ) * r ( 2, : )
      dy (3, : ) =  pf * ( r ( 1, : )**2 - r ( 2, : )**2 )
    CASE ( 1 )
      pf = SQRT ( 21.0_dp / ( 32.0_dp * pi ) )
      dy (1, : ) =  pf * ( 5.0_dp * r ( 3, : ) * r ( 3, : ) - 1.0_dp )
      dy (2, : ) =  0.0_dp
      dy (3, : ) =  pf * 10.0_dp * r ( 1, : ) * r ( 3, : )
     CASE ( 0 )
      pf = SQRT ( 7.0_dp / ( 16.0_dp * pi ) )
      dy (1, : ) =  0.0_dp
      dy (2, : ) =  0.0_dp
      dy (3, : ) =  pf * ( 15.0_dp * r ( 3, : ) * r ( 3, : ) - 3.0_dp)
    CASE ( -1 )
      pf = SQRT ( 21.0_dp / ( 32.0_dp * pi ) )
      dy (1, : ) =  0.0_dp
      dy (2, : ) =  pf * ( 5.0_dp * r ( 3, : ) * r ( 3, : ) - 1.0_dp )
      dy (3, : ) =  pf * 10.0_dp * r ( 2, : ) * r ( 3, : )
    CASE ( -2 )
      pf = SQRT ( 105.0_dp / ( 16.0_dp * pi ) )
      dy (1, : ) =  pf * 2.0_dp * r ( 2, : ) * r ( 3, : )
      dy (2, : ) =  pf * 2.0_dp * r ( 1, : ) * r ( 3, : )
      dy (3, : ) =  pf * 2.0_dp * r ( 1, : ) * r ( 2, : )
    CASE ( -3 )
      pf = SQRT ( 35.0_dp / ( 32.0_dp * pi ) )
      dy (1, : ) =  pf * 6.0_dp * r ( 1, : ) * r ( 2, : )
      dy (2, : ) =  pf * 3.0_dp * ( r ( 1, : ) * r ( 1, : ) - r ( 2, : ) * r ( 2, : ))
      dy (3, : ) =  0.0_dp
    END SELECT
  CASE DEFAULT
    IF ( m < -l .OR. m > l ) CALL stop_program(routineN,moduleN,__LINE__,&
                                               "m value out of bounds")
    lpm = fac ( l + ABS ( m ) )
    lmm = fac ( l - ABS ( m ) )
    IF ( m == 0 ) THEN
      t = 4.0_dp * pi
    ELSE
      t = 2.0_dp * pi
    END IF
    IF ( ABS ( lpm ) < EPSILON ( 1.0_dp ) ) THEN
      pf = REAL ( 2*l+1,KIND=dp)  / t
    ELSE
      pf = ( REAL ( 2*l+1,KIND=dp) * lmm ) / ( t * lpm )
    ENDIF
    pf = SQRT ( pf )
    DO i = 1, SIZE ( r ,2 )
      z = r ( 3, i )
      dplm = dlegendre ( z, l, ABS ( m ) )
      plm = legendre ( z, l, m )

      IF ( m == 0 ) THEN
        dy (1, i ) = 0.0_dp
        dy (2, i ) = 0.0_dp
        dy (3, i ) = pf * dplm
      ELSE
        rxy = SQRT ( r ( 1, i ) ** 2  +  r ( 2, i ) ** 2 )
        IF ( rxy < EPSILON ( 1.0_dp ) ) THEN
          dy (1:3, i ) = 0.0_dp
        ELSE
          cp = r ( 1, i ) / rxy
          sp = r ( 2, i ) / rxy
          dcp_dx = (rxy-2.0_dp*r ( 1, i ) ** 2) / (rxy*rxy)
          dcp_dy = -2.0_dp*r ( 1, i ) * r ( 2, i ) / (rxy*rxy)
          dsp_dy = (rxy-2.0_dp*r ( 2, i ) ** 2) / (rxy*rxy)
          dsp_dx = -2.0_dp*r ( 1, i ) * r ( 2, i ) / (rxy*rxy)
          IF ( m > 0 ) THEN
            dy (1, i ) = pf * plm * (dcosn_dcp ( m, cp, sp ) * dcp_dx +&
                                     dcosn_dsp( m, cp, sp )  * dsp_dx)
            dy (2, i ) = pf * plm * (dcosn_dcp ( m, cp, sp ) * dcp_dy +&
                                     dcosn_dsp( m, cp, sp )  * dsp_dy)
            dy (3, i ) = pf * dplm * cosn ( m, cp, sp )
          ELSE
            dy (1, i ) = pf * plm * (dsinn_dcp ( -m, cp, sp ) * dcp_dx +&
                                     dsinn_dsp( -m, cp, sp )  * dsp_dx)
            dy (2, i ) = pf * plm * (dsinn_dcp ( -m, cp, sp ) * dcp_dy +&
                                     dsinn_dsp( -m, cp, sp )  * dsp_dy)
            dy (3, i ) = pf * dplm * sinn ( -m , cp, sp )
          END IF
        END IF
      END IF
    END DO
  END SELECT

  END SUBROUTINE drvy_lm

END MODULE spherical_harmonics
