!******************************************************************************
!!****f* cp2k/whittaker [1.0] *
!!
!!   NAME
!!     whittaker
!!
!!   FUNCTION
!!     Calculates integral
!!     Int(
!!
!!   AUTHOR
!!     JGH 10-08-2004
!!
!!****
!******************************************************************************
!
    MODULE whittaker
!
!------------------------------------------------------------------------------!

  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: rootpi
  USE termination,                     ONLY: stop_program

!
  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: module_name = "whittaker"

  PRIVATE

  PUBLIC :: whittaker_c0, whittaker_ci

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE whittaker_c0 ( wc, r, expa, erfa, alpha, l, n )

!
!   int(y^(2+l) * exp(-alpha*y*y),y=0..x);
!
!   wc(:)    :: output
!   r(:)     :: coordinate
!   expa(:)  :: exp(-alpha*r(:)**2)
!   erfa(:)  :: erf(sqrt(alpha)*r(:))
!   alpha    :: exponent
!   l        :: L-quantum number
!   n        :: number of points
!
    REAL(KIND=dp), DIMENSION(n)              :: wc, r, expa, erfa
    REAL(KIND=dp), INTENT (IN)               :: alpha
    INTEGER, INTENT (IN)                     :: l, n

    CHARACTER(LEN=*), PARAMETER :: routine_name = "whittaker_c0"

    REAL ( dp ) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,&
                   t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,&
                   t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,&
                   t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,&
                   t60,t61
    REAL ( dp ) :: x
    INTEGER     :: i

! -----------------------------------------------------------------------------

    IF ( MOD(l,2) /= 0 ) THEN
        CALL stop_program(module_name,routine_name,__LINE__,&
             "Angular momentum has to be even")
    END IF

    SELECT CASE (l)

      CASE DEFAULT

        CALL stop_program(module_name,routine_name,__LINE__,&
             "Angular momentum out of range 0-14")

      CASE (0)

        t1 = sqrt(alpha)
        t2 = t1**2
        t11 = rootpi                      
        DO i = 1, n
          x  = r(i)
          t5 = x**2
          t7 = expa(i)
          t13 = erfa(i)    
          t18 = -1/t2/t1*(2*x*t7*t1-t11*t13)/4
          wc(i) = t18
        END DO

      CASE (2)

        t1 = sqrt(alpha)
        t2 = t1**2
        t3 = t2**2
        t17 = rootpi                      
        DO i = 1, n
          x  = r(i)
          t6 = x**2
          t9 = expa(i)        
          t19 = erfa(i)    
          t25 = -1/t3/t1*(4*t6*x*t9*t2*t1+6*x*t9*t1-3*t17*t19)/8
          wc(i) = t25
        END DO

      CASE (4)

        t1 = sqrt(alpha)
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
          t31 = -1/t4/t3*(8*t8*x*t11*t4*t1+20*t7*x*t11*t3+30*x*t11*t1-&
                15*t23*t25)/16
          wc(i) = t31
        END DO

      CASE (6)

        t8 = sqrt(alpha)
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
          t39 = -(16*t3*t2*t6*t11*t8+56*t3*x*t6*t10*t17+140*t2*t6*t10*t8+&
                210*x*t6*t17-105*t28*t30*alpha)/t11/t17/32
          wc(i) = t39
        END DO

      CASE (8)

        t8 = sqrt(alpha)
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
          t45 = -(32*t3*x*t6*t12*t10+144*t2*t16*t6*t12*t8+504*t2*x*t6*t11*t10+&
                1260*t16*t6*t28+1890*x*t6*t10-945*t34*t36*alpha)/t12/t28/64
          wc(i) = t45
        END DO

      CASE (10)

        t9 = sqrt(alpha)
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
          t50 = -(64*t4*t2*t7*t13*t12+352*t4*x*t7*t13*t19+&
                  1584*t3*t2*t7*t13*t9+5544*t3*x*t7*t30+&
                  13860*t2*t7*t12+20790*x*t7*t19-10395*t39*t41*alpha)/&
                  t13/t30/128
          wc(i) = t50
        END DO

      CASE (12)

        t9 = sqrt(alpha)
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
         t56 = -(128*t4*t3*t7*t13*t14+832*t4*t18*t7*t14*t21+&
                 4576*t4*x*t7*t14*t11+20592*t2*t18*t7*t14*t9+72072*t3*t7*t13+&
                 180180*t18*t7*t21+270270*x*t7*t11-135135*t44*t46*alpha)/&
                 t51/t9/256
          wc(i) = t56
        END DO

      CASE (14)

        t10 = sqrt(alpha)
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
          t61 = -(256*t5*t4*t8*t14*t10+1920*t5*t18*t8*t13*t22+&
                  12480*t5*t2*t8*t13*t28+68640*t5*x*t8*t13*t21+&
                  308880*t4*t8*t13*t10+1081080*t18*t8*t22+&
                  2702700*t2*t8*t28+4054050*x*t8*t21-&
                  2027025*t50*t52*alpha)/t14/t21/512
          wc(i) = t61
        END DO

    END SELECT

  END SUBROUTINE whittaker_c0

! *****************************************************************************

  SUBROUTINE whittaker_ci ( wc, r, expa, alpha, l, n )

!
!   int(y^(l+1) * exp(-alpha*y*y),y=x..infinity);
!
!             (-1 - 1/2 l~)                          2
!     1/2 alpha              GAMMA(1/2 l + 1, alpha x )
!
!
!   wc(:)    :: output
!   r(:)     :: coordinate
!   expa(:)  :: exp(-alpha*r(:)**2)
!   alpha    :: exponent
!   l        :: L-quantum number
!   n        :: number of points
!
    REAL(KIND=dp), DIMENSION(n)              :: wc, r, expa
    REAL(KIND=dp), INTENT (IN)               :: alpha
    INTEGER, INTENT (IN)                     :: l, n

    CHARACTER(LEN=*), PARAMETER :: routine_name = "whittaker_ci"

    REAL ( dp ) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,&
                   t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,&
                   t32,t33
    REAL ( dp ) :: x
    INTEGER     :: i

! -----------------------------------------------------------------------------

    IF ( MOD(l,2) /= 0 ) THEN
        CALL stop_program(module_name,routine_name,__LINE__,&
             "Angular momentum has to be even")
    END IF

    SELECT CASE (l)

      CASE DEFAULT

        CALL stop_program(module_name,routine_name,__LINE__,&
             "Angular momentum out of range 0-14")

      CASE (0)

        DO i = 1, n
          x  = r(i)
          t2 = x**2
          t4 = expa(i)           
          t6 = 1/alpha*t4/2
          wc(i) = t6
        END DO

      CASE (2)

        t6 = alpha**2
        DO i = 1, n
          x  = r(i)
          t1 = x**2
          t2 = alpha*t1
          t3 = expa(i)           
          t9 = t3*(t2+1)/t6/2
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
          t13 = t3*(t4*t5+2*t2+2)/t5/alpha/2
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
          t17 = t3*(t4*t1*t6*alpha+3*t4*t6+6*t2+6)/t14/2
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
          t21 = t3*(t5*t7+4*t4*t1*t6*alpha+12*t4*t6+24*t2+24)/t7/alpha/2
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
          t25 = t3*(t5*t1*t8*alpha+5*t5*t8+20*t4*t1*t7*alpha+60*t4*t7+&
                120*t2+120)/t8/t7/2
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
          t29 = t3*(t5*t4*t8*t7+6*t5*t1*t8*alpha+30*t5*t8+120*t4*t1*t18+&
                360*t4*t7+720*t2+720)/t8/t18/2
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
          t33 = t3*(t6*t5*t10*t9+7*t6*t4*t10*t8+42*t6*t1*t10*alpha+&
                210*t6*t10+840*t5*t9+2520*t4*t8+5040*t2+5040)/t30/2
          wc(i) = t33
        END DO

    END SELECT

  END SUBROUTINE whittaker_ci

! *****************************************************************************

END MODULE whittaker

! *****************************************************************************
