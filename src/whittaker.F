!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates special integrals
!> \author JGH 10-08-2004
! *****************************************************************************
    MODULE whittaker

  USE f77_blas
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: dfac,&
                                             fac,&
                                             rootpi
  USE termination,                     ONLY: stop_program
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'whittaker'
  REAL(KIND=dp), PARAMETER :: epsilon=1.e-2_dp

  PRIVATE

  PUBLIC :: whittaker_c0a, whittaker_ci

CONTAINS

! *****************************************************************************
!> \brief int(y^(2+l) * exp(-alpha*y*y),y=0..x);
!>        wc(:)    :: output
!>        r(:)     :: coordinate
!>        expa(:)  :: exp(-alpha*r(:)**2)
!>        erfa(:)  :: erf(sqrt(alpha)*r(:))
!>        alpha    :: exponent
!>        l1, l2   :: L-quantum number
!>        n        :: number of points
!>
!> \param wc ...
!> \param r ...
!> \param expa ...
!> \param erfa ...
!> \param alpha ...
!> \param l1 ...
!> \param l2 ...
!> \param n ...
!> \author JGH 10-08-2004
! *****************************************************************************
  SUBROUTINE whittaker_c0a ( wc, r, expa, erfa, alpha, l1, l2, n )
    INTEGER, INTENT(IN)                      :: n, l2, l1
    REAL(KIND=dp), INTENT(IN)                :: alpha
    REAL(KIND=dp), DIMENSION(n)              :: erfa, expa, r, wc

    CHARACTER(len=*), PARAMETER :: routineN = 'whittaker_c0a', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, k, l
    REAL(dp)                                 :: t1, x, y

    l = l1 + l2

    IF ( MOD(l,2) /= 0 ) THEN
        CALL stop_program(moduleN,routineN,__LINE__,&
             "Total angular momentum has to be even")
    END IF
    IF ( l1 < l2 ) THEN
        CALL stop_program(moduleN,routineN,__LINE__,"l1 >= l2")
    END IF

    t1 = SQRT(alpha)
    y = REAL(l,dp)
    DO i = 1, n
      x  = r(i)
      IF ( t1*x < epsilon ) THEN
        wc(i) = x**l1 * ( x**2/(3._dp+y) - alpha*x**4/(5._dp+y) + &
                   alpha**2*x**6/(14._dp+2._dp*y) - &
                   alpha**3*x**8/(54._dp+6._dp*y) + &
                   alpha**4*x**10/(256._dp+24._dp*y)- &
                   alpha**5*x**12/120._dp/(13._dp+y))
      ELSE
        wc(i) = -rootpi*erfa(i)*alpha*dfac(l+1)
        DO k = 0, l/2
          wc(i) = wc(i) + expa(i)*x**(2*k+1)*t1**(2*k+3)*&
                  dfac(l+1)/dfac(2*k+1)*2**(k+1)
        END DO
        wc(i) = -wc(i)/2._dp**(l/2+2)/t1**(l+5)/x**(l2+1)
      END IF
    END DO

  END SUBROUTINE whittaker_c0a

! *****************************************************************************
!> \brief int(y^(2+l) * exp(-alpha*y*y),y=0..x);
!>        wc(:)    :: output
!>        r(:)     :: coordinate
!>        expa(:)  :: exp(-alpha*r(:)**2)
!>        erfa(:)  :: erf(sqrt(alpha)*r(:))
!>        alpha    :: exponent
!>        l        :: L-quantum number
!>        n        :: number of points
!>
!> \param wc ...
!> \param r ...
!> \param expa ...
!> \param erfa ...
!> \param alpha ...
!> \param l ...
!> \param n ...
!> \author JGH 10-08-2004
! *****************************************************************************
  SUBROUTINE whittaker_c0 ( wc, r, expa, erfa, alpha, l, n )
    INTEGER, INTENT(IN)                      :: n, l
    REAL(KIND=dp), INTENT(IN)                :: alpha
    REAL(KIND=dp), DIMENSION(n)              :: erfa, expa, r, wc

    CHARACTER(len=*), PARAMETER :: routineN = 'whittaker_c0', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, k
    REAL(dp) :: t1, t10, t11, t12, t13, t14, t16, t17, t18, t19, t2, t21, &
      t22, t23, t25, t28, t3, t30, t31, t34, t36, t39, t4, t41, t44, t45, &
      t46, t5, t50, t51, t52, t56, t6, t61, t7, t8, t9, x

    IF ( MOD(l,2) /= 0 ) THEN
        CALL stop_program(moduleN,routineN,__LINE__,&
             "Angular momentum has to be even")
    END IF

    SELECT CASE (l)

      CASE DEFAULT

        t1 = SQRT(alpha)
        DO i = 1, n
          x  = r(i)
          wc(i) = -rootpi*erfa(i)*alpha*dfac(l+1)
          DO k = 0, l/2
            wc(i) = wc(i) + expa(i)*x**(2*k+1)*t1**(2*k+3)*&
                    dfac(l+1)/dfac(2*k+1)*2**(k+1)
          END DO
          wc(i) = -wc(i)/2._dp**(l/2+2)/t1**(l+5)
        END DO

      CASE (0)

        t1 = SQRT(alpha)
        t2 = t1**2
        t11 = rootpi
        DO i = 1, n
          x  = r(i)
          t5 = x**2
          t7 = expa(i)
          t13 = erfa(i)
          t18 = -1._dp/t2/t1*(2._dp*x*t7*t1-t11*t13)/4._dp
          wc(i) = t18
        END DO

      CASE (2)

        t1 = SQRT(alpha)
        t2 = t1**2
        t3 = t2**2
        t17 = rootpi
        DO i = 1, n
          x  = r(i)
          t6 = x**2
          t9 = expa(i)
          t19 = erfa(i)
          t25 = -1._dp/t3/t1*(4._dp*t6*x*t9*t2*t1+6._dp*x*t9*t1-3*t17*t19)/8._dp
          wc(i) = t25
        END DO

      CASE (4)

        t1 = SQRT(alpha)
        t2 = t1**2
        t3 = t2*t1
        t4 = t2**2
        t23 = rootpi
        DO i = 1, n
          x  = r(i)
          t7 = x**2
          t8 = t7**2
          t11 = expa(i)
          t25 = erfa(i)
          t31 = -1._dp/t4/t3*(8._dp*t8*x*t11*t4*t1+20._dp*t7*x*t11*t3+30._dp*x*t11*t1-&
                15._dp*t23*t25)/16._dp
          wc(i) = t31
        END DO

      CASE (6)

        t8 = SQRT(alpha)
        t9 = t8**2
        t10 = t9**2
        t11 = t10**2
        t17 = t9*t8
        t28 = rootpi
        DO i = 1, n
          x  = r(i)
          t1 = x**2
          t2 = t1*x
          t3 = t1**2
          t6 = expa(i)
          t30 = erfa(i)
          t39 = -(16._dp*t3*t2*t6*t11*t8+56._dp*t3*x*t6*t10*t17+140._dp*t2*t6*t10*t8+&
                210._dp*x*t6*t17-105._dp*t28*t30*alpha)/t11/t17/32._dp
          wc(i) = t39
        END DO

      CASE (8)

        t8 = SQRT(alpha)
        t9 = t8**2
        t10 = t9*t8
        t11 = t9**2
        t12 = t11**2
        t34 = rootpi
        DO i = 1, n
          x  = r(i)
          t1 = x**2
          t2 = t1**2
          t3 = t2**2
          t6 = expa(i)
          t16 = t1*x
          t28 = t11*t8
          t36 = erfa(i)
          t45 = -(32._dp*t3*x*t6*t12*t10+144._dp*t2*t16*t6*t12*t8+504._dp*t2*x*t6*t11*t10+&
                1260._dp*t16*t6*t28+1890._dp*x*t6*t10-945._dp*t34*t36*alpha)/t12/t28/64._dp
          wc(i) = t45
        END DO

      CASE (10)

        t9 = SQRT(alpha)
        t10 = t9**2
        t11 = t10**2
        t12 = t11*t9
        t13 = t11**2
        t19 = t10*t9
        t30 = t11*t19
        t39 = rootpi
        DO i = 1, n
          x  = r(i)
          t1 = x**2
          t2 = t1*x
          t3 = t1**2
          t4 = t3**2
          t7 = expa(i)
          t41 = erfa(i)
          t50 = -(64._dp*t4*t2*t7*t13*t12+352._dp*t4*x*t7*t13*t19+&
                  1584._dp*t3*t2*t7*t13*t9+5544._dp*t3*x*t7*t30+&
                  13860._dp*t2*t7*t12+20790._dp*x*t7*t19-10395._dp*t39*t41*alpha)/&
                  t13/t30/128._dp
          wc(i) = t50
        END DO

      CASE (12)

        t9 = SQRT(alpha)
        t10 = t9**2
        t11 = t10*t9
        t12 = t10**2
        t13 = t12*t11
        t14 = t12**2
        t44 = rootpi
        DO i = 1, n
          x  = r(i)
          t1 = x**2
          t2 = t1**2
          t3 = t2*x
          t4 = t2**2
          t7 = expa(i)
         t18 = t1*x
         t21 = t12*t9
         t46 = erfa(i)
         t51 = t14**2
         t56 = -(128._dp*t4*t3*t7*t13*t14+832._dp*t4*t18*t7*t14*t21+&
                 4576._dp*t4*x*t7*t14*t11+20592._dp*t2*t18*t7*t14*t9+72072._dp*t3*t7*t13+&
                 180180._dp*t18*t7*t21+270270._dp*x*t7*t11-135135._dp*t44*t46*alpha)/&
                 t51/t9/256._dp
          wc(i) = t56
        END DO

      CASE (14)

        t10 = SQRT(alpha)
        t11 = t10**2
        t12 = t11**2
        t13 = t12**2
        t14 = t13**2
        t21 = t11*t10
        t22 = t12*t21
        t28 = t12*t10
        t50 = rootpi
        DO i = 1, n
          x  = r(i)
          t1 = x**2
          t2 = t1*x
          t3 = t1**2
          t4 = t3*t2
          t5 = t3**2
          t8 = expa(i)
          t18 = t3*x
          t52 = erfa(i)
          t61 = -(256._dp*t5*t4*t8*t14*t10+1920._dp*t5*t18*t8*t13*t22+&
                  12480._dp*t5*t2*t8*t13*t28+68640._dp*t5*x*t8*t13*t21+&
                  308880._dp*t4*t8*t13*t10+1081080._dp*t18*t8*t22+&
                  2702700._dp*t2*t8*t28+4054050._dp*x*t8*t21-&
                  2027025._dp*t50*t52*alpha)/t14/t21/512._dp
          wc(i) = t61
        END DO

    END SELECT

  END SUBROUTINE whittaker_c0

! *****************************************************************************
!> \brief int(y^(l+1) * exp(-alpha*y*y),y=x..infinity);
!>
!>                  (-1 - 1/2 l~)                          2
!>          1/2 alpha              GAMMA(1/2 l + 1, alpha x )
!>
!>
!>        wc(:)    :: output
!>        r(:)     :: coordinate
!>        expa(:)  :: exp(-alpha*r(:)**2)
!>        alpha    :: exponent
!>        l        :: L-quantum number
!>        n        :: number of points
!>
!> \param wc ...
!> \param r ...
!> \param expa ...
!> \param alpha ...
!> \param l ...
!> \param n ...
!> \author JGH 10-08-2004
! *****************************************************************************
  SUBROUTINE whittaker_ci ( wc, r, expa, alpha, l, n )
    INTEGER, INTENT(IN)                      :: n, l
    REAL(KIND=dp), INTENT(IN)                :: alpha
    REAL(KIND=dp), DIMENSION(n)              :: expa, r, wc

    CHARACTER(len=*), PARAMETER :: routineN = 'whittaker_ci', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, k
    REAL(dp)                                 :: t1, t10, t13, t14, t17, t18, &
                                                t2, t21, t25, t29, t3, t30, &
                                                t33, t4, t5, t6, t7, t8, t9, x

    IF ( MOD(l,2) /= 0 ) THEN
        CALL stop_program(moduleN,routineN,__LINE__,&
             "Angular momentum has to be even")
    END IF

    SELECT CASE (l)

      CASE DEFAULT

        DO i = 1, n
          x  = r(i)
          wc(i) = 0._dp
          DO k = 0, l/2
            wc(i) = wc(i) + alpha**k*x**(2*k)*fac(l/2)/fac(k)
          END DO
          wc(i) = 0.5_dp*wc(i)/alpha**(l/2+1)*expa(i)
        END DO

      CASE (0)

        DO i = 1, n
          x  = r(i)
          t2 = x**2
          t4 = expa(i)
          t6 = 1._dp/alpha*t4/2._dp
          wc(i) = t6
        END DO

      CASE (2)

        t6 = alpha**2
        DO i = 1, n
          x  = r(i)
          t1 = x**2
          t2 = alpha*t1
          t3 = expa(i)
          t9 = t3*(t2+1)/t6/2._dp
          wc(i) = t9
        END DO

      CASE (4)

        t5 = alpha**2
        DO i = 1, n
          x  = r(i)
          t1 = x**2
          t2 = alpha*t1
          t3 = expa(i)
          t4 = t1**2
          t13 = t3*(t4*t5+2._dp*t2+2._dp)/t5/alpha/2._dp
          wc(i) = t13
        END DO

      CASE (6)

        t6 = alpha**2
        t14 = t6**2
        DO i = 1, n
          x  = r(i)
          t1 = x**2
          t2 = alpha*t1
          t3 = expa(i)
          t4 = t1**2
          t17 = t3*(t4*t1*t6*alpha+3._dp*t4*t6+6._dp*t2+6._dp)/t14/2._dp
          wc(i) = t17
        END DO

      CASE (8)

        t6 = alpha**2
        t7 = t6**2
        DO i = 1, n
          x  = r(i)
          t1 = x**2
          t2 = alpha*t1
          t3 = expa(i)
          t4 = t1**2
          t5 = t4**2
          t21 = t3*(t5*t7+4._dp*t4*t1*t6*alpha+12._dp*t4*t6+24._dp*t2+24._dp)/t7/alpha/2._dp
          wc(i) = t21
        END DO

      CASE (10)

        t7 = alpha**2
        t8 = t7**2
        DO i = 1, n
          x  = r(i)
          t1 = x**2
          t2 = alpha*t1
          t3 = expa(i)
          t4 = t1**2
          t5 = t4**2
          t25 = t3*(t5*t1*t8*alpha+5._dp*t5*t8+20._dp*t4*t1*t7*alpha+60._dp*t4*t7+&
                120._dp*t2+120._dp)/t8/t7/2._dp
          wc(i) = t25
        END DO

      CASE (12)

        t7 = alpha**2
        t8 = t7**2
        t18 = t7*alpha
        DO i = 1, n
          x  = r(i)
          t1 = x**2
          t2 = alpha*t1
          t3 = expa(i)
          t4 = t1**2
          t5 = t4**2
          t29 = t3*(t5*t4*t8*t7+6._dp*t5*t1*t8*alpha+30._dp*t5*t8+120._dp*t4*t1*t18+&
                360._dp*t4*t7+720._dp*t2+720._dp)/t8/t18/2._dp
          wc(i) = t29
        END DO

      CASE (14)

        t8 = alpha**2
        t9 = t8*alpha
        t10 = t8**2
        t30 = t10**2
        DO i = 1, n
          x  = r(i)
          t1 = x**2
          t2 = alpha*t1
          t3 = expa(i)
          t4 = t1**2
          t5 = t4*t1
          t6 = t4**2
          t33 = t3*(t6*t5*t10*t9+7*t6*t4*t10*t8+42._dp*t6*t1*t10*alpha+&
                210._dp*t6*t10+840._dp*t5*t9+2520._dp*t4*t8+5040._dp*t2+5040._dp)/t30/2._dp
          wc(i) = t33
        END DO

    END SELECT

  END SUBROUTINE whittaker_ci

END MODULE whittaker

