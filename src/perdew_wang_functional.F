!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****m* cp2k/perdew_wang_functional *
!!
!!   NAME
!!     perdew_wang_functional
!!
!!   FUNCTION
!!     Calculates the correlation energy according to the
!!     functional due to Perdew and Wang
!!
!!   SYNOPSIS
!!     
!!
!!   MODIFICATION HISTORY
!!     , TCH, working version
!!
!!   SEE ALSO
!!     
!!
!!****

MODULE perdew_wang_functional

  USE functionals_utilities, ONLY : calc_fx, calc_rs, calc_rs_d, set_util, calc_z
  USE kinds, ONLY : dbl
  USE termination, ONLY : stop_program
  
  IMPLICIT NONE

  LOGICAL :: initialized = .FALSE.
  REAL(DBL), DIMENSION(-1:1) :: A, a1, b1, b2, b3, b4
  REAL(DBL), PARAMETER :: epsilon = 5.E-13

  PRIVATE

  PUBLIC :: perdew_wang

  INTERFACE perdew_wang
     MODULE PROCEDURE pw_info, pw_init, pw_lda, pw_lsd
  END INTERFACE

CONTAINS

! *****************************************************************************

  SUBROUTINE pw_init(method, cutoff)
    
    IMPLICIT NONE

    CHARACTER(3), INTENT(IN) :: method
    REAL(DBL), INTENT(IN) :: cutoff
    
    CALL set_util(cutoff)
    
    initialized = .FALSE.

    ! values for -ac are the same for all methods
    A(-1)  = 0.016887
    a1(-1) = 0.11125 
    b1(-1) = 10.357  
    b2(-1) = 3.6231
    b3(-1) = 0.88026
    b4(-1) = 0.49671

    SELECT CASE (method)

    CASE DEFAULT
       CALL stop_program("perdew_zunger_functional/pz_init", "unknown method")

    CASE ('PWO')
       A(0)  = 0.031091_dbl ; A(1)  = 0.015545_dbl
       a1(0) = 0.21370_dbl  ; a1(1) = 0.20548_dbl
       b1(0) = 7.5957_dbl   ; b1(1) = 14.1189_dbl
       b2(0) = 3.5876_dbl   ; b2(1) = 6.1977_dbl
       b3(0) = 1.6382_dbl   ; b3(1) = 3.3662_dbl
       b4(0) = 0.49294_dbl  ; b4(1) = 0.62517_dbl

    CASE ('DMC')
       A(0)  = 0.031091_dbl ;  A(1)  = 0.015545_dbl
       a1(0) = 0.026481_dbl ;  a1(1) = 0.022465_dbl
       b1(0) = 7.5957_dbl   ;  b1(1) = 14.1189_dbl
       b2(0) = 3.5876_dbl   ;  b2(1) = 6.1977_dbl
       b3(0) = -0.46647_dbl ;  b3(1) = -0.56043_dbl
       b4(0) = 0.13354_dbl  ;  b4(1) = 0.11313_dbl
       
    CASE ('VMC') 
       A(0)  = 0.031091_dbl ; A(1)  = 0.015545_dbl
       a1(0) = -0.002257_dbl; a1(1) = -0.009797_dbl
       b1(0) = 7.5957_dbl   ; b1(1) = 14.1189_dbl
       b2(0) = 3.5876_dbl   ; b2(1) = 6.1977_dbl
       b3(0) = -0.52669_dbl ; b3(1) = -0.91381_dbl
       b4(0) = 0.03755_dbl  ; b4(1) = 0.01538_dbl
       
    END SELECT
          
    initialized = .TRUE.

  END SUBROUTINE pw_init
  
! *****************************************************************************

  SUBROUTINE pw_lda(rho, ec, order)

    IMPLICIT NONE

    REAL(DBL), DIMENSION(1:), INTENT(IN) :: rho
    REAL(DBL), DIMENSION(1:,0:), INTENT(OUT) :: ec
    INTEGER, INTENT(IN), OPTIONAL :: order

    INTEGER :: k, n, m, order_ = 0
    LOGICAL, DIMENSION(0:3) :: calc = .FALSE.

    REAL(DBL), DIMENSION(:,:), ALLOCATABLE :: rs, ec0, ac
    
    IF (PRESENT(order)) order_ = order

    !! check order_ with an assertion here instead
    IF (order_ >= 0) THEN
       calc(0:order_) = .TRUE.
       m = order_ + 1
    ELSE
       order_ = -1 * order_
       calc(order_) = .TRUE.
       m = 1
    END IF

    n = SIZE(rho)
    IF (SIZE(ec,1) < n) &
         CALL stop_program("perdew_wang_functional/pw_lda", "SIZE(ec,1) too small")
    IF (SIZE(ec,2) < m) &
         CALL stop_program("perdew_wang_functional/pw_lda", "SIZE(ec,2) too small")

    IF (.NOT.initialized) CALL pw_init('PWO', epsilon)

    ALLOCATE(rs(1:n,0:order_), ec0(1:n,0:order_))

    CALL calc_rs_d(rho, rs(:,0:order_), order_)

    DO k=1, n
       CALL calc_g(rs(k,0), 0, ec0(k,0:0), order_)
    END DO

    m = 0
    IF (calc(0)) THEN
       ec(:,0) = ec0(:,0)
       m = m + 1
    END IF
    IF (calc(1)) THEN
       ec(:,m) = ec0(:,1)*rs(:,1)
       m =  m + 1
    END IF
    IF (calc(2)) THEN
       ec(:,m) = ec0(:,2)*rs(:,1)**2 + ec0(:,1)*rs(:,2)
       m = m + 1
    END IF
    IF (calc(3)) THEN
       ec(:,m) = ec0(:,3)*rs(:,1)**3 + ec0(:,2)*3._dbl*rs(:,1)*rs(:,2) + ec0(:,1)*rs(:,3)
    END IF

    DEALLOCATE(rs, ec0)

  END SUBROUTINE pw_lda
  
! *****************************************************************************

  SUBROUTINE pw_lsd(a, b, ec, order)

    IMPLICIT NONE

    REAL(DBL), DIMENSION(1:), INTENT(IN) :: a, b
    REAL(DBL), DIMENSION(1:,0:), INTENT(OUT) :: ec
    INTEGER, INTENT(IN), OPTIONAL :: order

    INTEGER :: k, l, n, m, order_ = 0
    LOGICAL, DIMENSION(0:3) :: calc = .FALSE.
    REAL(DBL), DIMENSION(:,:), ALLOCATABLE :: f
    REAL(DBL), DIMENSION(1:SIZE(a)) :: rs
    REAL(DBL) :: fpp, rho
    REAL(DBL) :: tr, tz, trr, trz, tzz, trrr, trrz, trzz, tzzz
    REAL(DBL), DIMENSION(:), ALLOCATABLE :: ac, e0, e1, r
    REAL(DBL), DIMENSION(:,:), ALLOCATABLE :: z

    !! d^2f(0)/dz^2 (second derivative of f at z=0)
    fpp = 1._dbl/(4._dbl/9._dbl*2._dbl/(2._dbl**(4._dbl/3._dbl) - 2))

    IF (PRESENT(order)) order_ = order
    IF (ABS(order_) > 3) CALL stop_program("perdew_wang_functional/pw_lsd", &
         "order must be between -3 and 3")

    IF (.NOT.initialized) CALL pw_init('PWO', epsilon)

    IF (order_ > 0) THEN
       calc(0:order_) = .TRUE.
       m = SUM((/ (k+1,k=0,order_) /))
    ELSE
       order_ = -1 * order_
       calc(order_) = .TRUE.
       m = 1
    END IF   

    n = SIZE(a)
    IF (SIZE(ec,1) < n) &
         CALL stop_program("perdew_wang_functional/pw_lda", "SIZE(ec,1) too small")
    IF (SIZE(ec,2) < m) &
         CALL stop_program("perdew_wang_functional/pw_lda", "SIZE(ec,2) too small")

    ALLOCATE (ac(0:order_), &
         e0(0:order_), &
         e1(0:order_), &
         r(0:order_), &
         z(0:order_,0:order_))
    
    ALLOCATE (f(1:n,0:order_))
 
    CALL calc_rs(a+b, rs)
    CALL calc_fx(a, b, f, order_)

    DO k=1, n

       r(0) = rs(k); rho = a(k) + b(k)

       CALL calc_g(r(0), -1, ac, order_)
       CALL calc_g(r(0), 0, e0, order_)
       CALL calc_g(r(0), 1, e1, order_)
       CALL calc_z(a(k), b(k), z, order_)
       
       !! calculate first partial derivatives
       IF (order_ >= 1) THEN

          r(1) = (-1._dbl/3._dbl)*r(0)/rho

          tr = e0(1) &
               + fpp*ac(1)*f(k,0) &
               - fpp*ac(1)*f(k,0)*z(0,0)**4 &
               + (e1(1)-e0(1))*f(k,0)*z(0,0)**4

          tz = fpp*ac(0)*f(k,1) &
               - fpp*ac(0)*f(k,1)*z(0,0)**4 &
               - fpp*ac(0)*f(k,0)*4._dbl*z(0,0)**3 &
               + (e1(0)-e0(0))*f(k,1)*z(0,0)**4 &
               + (e1(0)-e0(0))*f(k,0)*4._dbl*z(0,0)**3
          
       END IF
       
       !! calculate second partial derivatives
       IF (order_ >= 2) THEN
          
          r(2) = (-4._dbl/3._dbl)*r(1)/rho
          
          trr = e0(2) &
               + fpp*ac(2)*f(k,0) &
               - fpp*ac(2)*f(k,0)*z(0,0)**4 &
               + (e1(2)-e0(2))*f(k,0)*z(0,0)**4
          
          trz = fpp*ac(1)*f(k,1) &
               - fpp*ac(1)*f(k,1)*z(0,0)**4 &
               - fpp*ac(1)*f(k,0)*4._dbl*z(0,0)**3 &
               + (e1(1)-e0(1))*f(k,1)*z(0,0)**4 &
               + (e1(1)-e0(1))*f(k,0)*4._dbl*z(0,0)**3
          
          tzz = fpp*ac(0)*f(k,2) &
               - fpp*ac(0)*f(k,2)*z(0,0)**4 &
               - fpp*ac(0)*f(k,1)*8._dbl*z(0,0)**3 &
               - fpp*ac(0)*f(k,0)*12._dbl*z(0,0)**2 &
               + (e1(0)-e0(0))*f(k,2)*z(0,0)**4 &
               + (e1(0)-e0(0))*f(k,1)*8._dbl*z(0,0)**3 &
               + (e1(0)-e0(0))*f(k,0)*12._dbl*z(0,0)**2
       END IF

       !! calculate third derivatives
       IF (order_ >= 3) THEN

          r(3) = (-7._dbl/3._dbl)*r(2)/rho

          trrr = e0(3) &
               + fpp*ac(3)*f(k,0) &
               - fpp*ac(3)*f(k,0)*z(0,0)**4 &
               + (e1(3)-e0(3))*f(k,0)*z(0,0)**4

          trrz = fpp*ac(2)*f(k,1) &
               - fpp*ac(2)*f(k,1)*z(0,0)**4 &
               - fpp*ac(2)*f(k,0)*4._dbl*z(0,0)**3 &
               + (e1(2)-e0(2))*f(k,1)*z(0,0)**4 &
               + (e1(2)-e0(2))*f(k,0)*4._dbl*z(0,0)**3

          trzz = fpp*ac(1)*f(k,2) &
               - fpp*ac(1)*f(k,2)*z(0,0)**4 &
               - fpp*ac(1)*f(k,1)*8._dbl*z(0,0)**3 &
               - fpp*ac(1)*f(k,0)*12._dbl*z(0,0)**2 &
               + (e1(1)-e0(1))*f(k,2)*z(0,0)**4 &
               + (e1(1)-e0(1))*f(k,1)*8._dbl*z(0,0)**3 &
               + (e1(1)-e0(1))*f(k,0)*12._dbl*z(0,0)**2

          tzzz = fpp*ac(0)*f(k,3) &
               - fpp*ac(0)*f(k,3)*z(0,0)**4 &
               - fpp*ac(0)*f(k,2)*12._dbl*z(0,0)**3 &
               - fpp*ac(0)*f(k,1)*36._dbl*z(0,0)**2 &
               - fpp*ac(0)*f(k,0)*24._dbl*z(0,0) &
               + (e1(0)-e0(0))*f(k,3)*z(0,0)**4 &
               + (e1(0)-e0(0))*f(k,2)*12._dbl*z(0,0)**3 &
               + (e1(0)-e0(0))*f(k,1)*36._dbl*z(0,0)**2 &
               + (e1(0)-e0(0))*f(k,0)*24._dbl*z(0,0)
          
       END IF


       m=0
       IF (calc(0)) THEN
          
          ec(k,0) = e0(0) &
               + fpp*ac(0)*f(k,0)*(1._dbl - z(0,0)**4) &
               + (e1(0)-e0(0))*f(k,0)*z(0,0)**4
          m = m + 1

       END IF

       IF (calc(1)) THEN

          ec(k,m) = tr*r(1) + tz*z(1,0)
          ec(k,m+1) = tr*r(1) + tz*z(0,1)

          m = m + 2
       END IF
       
       IF (calc(2)) THEN
          
          ec(k,m) = &
               trr*r(1)**2 + 2._dbl*trz*r(1)*z(1,0) + tr*r(2) &
               + tzz*z(1,0)**2 + tz*z(2,0)
          ec(k,m+1) = &
               trr*r(1)**2 + trz*r(1)*(z(0,1)+z(1,0)) + tr*r(2) &
               + tzz*z(1,0)*z(0,1) + tz*z(1,1)
          ec(k,m+2) = &
               trr*r(1)**2 + 2._dbl*trz*r(1)*z(0,1) + tr*r(2) &
               + tzz*z(0,1)**2 + tz*z(0,2)

          m = m + 3
       END IF
       
       IF (calc(3)) THEN
          
          ec(k,m) = &
               trrr*r(1)**3 + 3._dbl*trrz*r(1)**2*z(1,0) &
               + 3._dbl*trr*r(1)*r(2) + 3._dbl*trz*r(2)*z(1,0) + tr*r(3) &
               + 3._dbl*trzz*r(1)*z(1,0)**2 + tzzz*z(1,0)**3 &
               + 3._dbl*trz*r(1)*z(2,0) &
               + 3._dbl*tzz*z(1,0)*z(2,0) + tz*z(3,0)
          ec(k,m+1) = &
               trrr*r(1)**3 + trrz*r(1)**2*(2._dbl*z(1,0)+z(0,1)) &
               + 2._dbl*trzz*r(1)*z(1,0)*z(0,1) &
               + 2._dbl*trz*(r(2)*z(1,0)+r(1)*z(1,1)) &
               + 3._dbl*trr*r(2)*r(1) + trz*r(2)*z(0,1) + tr*r(3) &
               + trzz*r(1)*z(1,0)**2 + tzzz*z(1,0)**2*z(0,1) &
               + 2._dbl*tzz*z(1,0)*z(1,1) &
               + trz*r(1)*z(2,0) + tzz*z(2,0)*z(0,1) + tz*z(2,1)
          ec(k,m+2) = &
               trrr*r(1)**3 + trrz*r(1)**2*(2._dbl*z(0,1)+z(1,0)) &
               + 2._dbl*trzz*r(1)*z(0,1)*z(1,0) &
               + 2._dbl*trz*(r(2)*z(0,1)+r(1)*z(1,1)) &
               + 3._dbl*trr*r(2)*r(1) + trz*r(2)*z(1,0) + tr*r(3) &
               + trzz*r(1)*z(0,1)**2 + tzzz*z(0,1)**2*z(1,0) &
               + 2._dbl*tzz*z(0,1)*z(1,1) &
               + trz*r(1)*z(0,2) + tzz*z(0,2)*z(1,0) + tz*z(1,2)
          ec(k,m+3) = &
               trrr*r(1)**3 + 3._dbl*trrz*r(1)**2*z(0,1) &
               + 3._dbl*trr*r(1)*r(2) + 3._dbl*trz*r(2)*z(0,1) + tr*r(3) &
               + 3._dbl*trzz*r(1)*z(0,1)**2 + tzzz*z(0,1)**3 &
               + 3._dbl*trz*r(1)*z(0,2) &
               + 3._dbl*tzz*z(0,1)*z(0,2) + tz*z(0,3)

       END IF

    END DO
    
    DEALLOCATE (ac, e0, e1, f, z)

  END SUBROUTINE pw_lsd

! *****************************************************************************

  SUBROUTINE calc_g(r, z, g, order)

!   ---------------------------------------------------------------------------
!   Calculates g and its derivatives wrt r up to 3rd order, where:
!
!   g = 2A(1+a1*r)ln(1+1/(2A(b1*r^1/2 + b2*r + b3*r^(3/2) + b4*r^2))).
!   ---------------------------------------------------------------------------

    IMPLICIT NONE

    REAL(DBL), INTENT(IN) :: r
    INTEGER, INTENT(IN) :: z
    REAL(DBL), DIMENSION(0:), INTENT(OUT) :: g
    INTEGER, INTENT(IN) :: order

    INTEGER :: k, n, order_
    REAL(DBL) :: A_, a1_, b1_, b2_, b3_, b4_, sr, rsr, rr
    REAL(DBL) :: t3, t11, t12, t14, t15, t16, t20, t22
    REAL(DBL) :: t40, t44, t45, t47, t48, t55, t56
        
    A_ = A(z); a1_ = a1(z)
    b1_ = b1(z); b2_ = b2(z); b3_ = b3(z); b4_ = b4(z)
    
    sr = SQRT(r)
    rsr = r*sr
    rr = r*r

    t3 = 1._dbl+a1_*r
    t11 = b1_*sr + b2_*r + b3_*rsr + b4_*rr
    t12 = t11**2
    t15 = 1._dbl + 0.5_dbl/A_/t11
    t16 = LOG(t15)
    t20 = 0.5_dbl*b1_/sr + b2_ + 1.5_dbl*b3_*sr + 2._dbl*b4_*r

    ! order 0 must always be calculated
    g(0) = -2._dbl*A_*t3*t16

    IF (order >= 1) THEN
       
       g(1) = -2._dbl*A_*a1_*t16 + t3*t20/(t12*t15)
       
    END IF
    
    IF (order >= 2) THEN

       t40 = -0.25_dbl*b1_/rsr + 0.75_dbl*b3_/sr + 2._dbl*b4_

       g(2) = 2._dbl*a1_*t20/(t12*t15) &
            - 2._dbl*(t20**2)*t3/(t12*t11*t15) &
            + t3*t40/(t12*t15) &
            + 0.5_dbl*t3*(t20**2)/(A_*(t12**2)*(t15**2))

    END IF

    IF (order >= 3) THEN

       t14 = 1._dbl/t12/t11
       t22 = t20**2
       t56 = t22*t20
       t47 = t15**2
       t48 = 1._dbl/t47
          
       t44 = t12**2
       t45 = 1._dbl/t44
       t55 = t3*t45

       g(3) = &
            - 6._dbl*a1_*t14*t22/t15 &
            + 3._dbl*a1_*t40/(t15*t12) &
            + 1.5_dbl*a1_*t45*t22*t48/A_ &
            + 6._dbl*t55*t56/t15 &
            - 6._dbl*t3*t14*t20*t40/t15 &
            - 3._dbl*t3*t56*t48/(A_*t44*t11) &
            + 0.375_dbl*t3*(b1_/(rr*sr)-b3_/rsr)/(t12*t15) &
            + 1.5_dbl*t55*t40*t48*t20/A_ &
            + 0.5_dbl*t3*t56/((A_**2)*t44*t12*t47*t15)
    END IF

  END SUBROUTINE calc_g

! *****************************************************************************

  SUBROUTINE pw_info ( reference, shortform )

    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: reference
    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: shortform

    IF ( PRESENT ( reference ) ) THEN
       reference = "J. P. Perdew and Yue Wang," &
                   //" Phys. Rev. B 45, 13244 (1992)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "J. P. Perdew et al., PRB 45, 13244 (1992)"
    END IF

  END SUBROUTINE pw_info

! *****************************************************************************

END MODULE perdew_wang_functional
