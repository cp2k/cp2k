!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****m* cp2k/perdew_zunger_funcional.F *
!!
!!   NAME
!!     perdew_zunger_functional
!!
!!   FUNCTION
!!     
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

MODULE perdew_zunger_functional
  
  USE functionals_utilities, ONLY : calc_fx, calc_rs, set_util
  USE kinds, ONLY : dbl
  USE termination, ONLY : stop_program
  
  IMPLICIT NONE

  LOGICAL :: initialized = .FALSE.
  REAL(DBL), DIMENSION(0:1) :: A = (0,0), B = (0,0), C = (0,0), D = (0,0), &
                               b1 = (0,0), b2 = (0,0), g = (0,0)

  REAL(DBL), PARAMETER :: epsilon = 5.E-13

  PRIVATE
  
  PUBLIC :: perdew_zunger
  
  INTERFACE perdew_zunger
     MODULE PROCEDURE pz, pz_init
  END INTERFACE

CONTAINS

! *****************************************************************************

  SUBROUTINE pz ( rhoa, rhob, ec, o_order )
    
    IMPLICIT NONE
    
    REAL(DBL), DIMENSION(:), INTENT(IN) :: rhoa, rhob
    REAL(DBL), DIMENSION(:,0:), INTENT(OUT) :: ec
    INTEGER, INTENT(IN), OPTIONAL :: o_order

    INTEGER :: k, m, order
    LOGICAL, DIMENSION(0:3) :: calc
    REAL(DBL), DIMENSION(1:SIZE(rhoa),0:3) :: rs
    REAL(DBL), DIMENSION(1:SIZE(rhoa),0:3) :: f, ec0, ec1

    calc = .FALSE.

    IF (PRESENT(o_order)) THEN
       order = o_order
    ELSE
       order = 0
    END IF

    IF (ABS(order) <= 3) THEN
       IF (order >= 0) THEN
          calc(0:order) = .TRUE.
          m = SUM((/ (m+1,m=0,order) /))
       ELSE
          order = -1 * order
          calc(order) = .TRUE.
          m = order + 1
       END IF
    ELSE
       CALL stop_program("perdew_zunger_functional/pz", &
            "parameter order must be between -3 and 3")
    END IF

    IF (SIZE(ec,2) /= m) & 
         CALL stop_program("perdew_zunger_functional/pz", &
         "array ec too small")

    IF (.NOT.initialized) THEN
       CALL pz_init(epsilon, 'DMC')
    END IF

    CALL calc_fx(rhoa, rhob, f, order)
    CALL calc_rs(rhoa + rhob, rs(:,0))
    CALL calc_g(rs(:,0), ec0, ec1, order)

    IF (order >= 1) rs(:,1) = (-1._dbl/3._dbl)*rs(:,0)/(rhoa(:) + rhob(:))
    IF (order >= 2) rs(:,2) = (-4._dbl/3._dbl)*rs(:,1)/(rhoa(:) + rhob(:))
    IF (order >= 3) rs(:,3) = (-7._dbl/3._dbl)*rs(:,2)/(rhoa(:) + rhob(:))

    m = 0
    IF (calc(0)) THEN
       ec(:,m) = ec0(:,0) + (ec1(:,0) - ec0(:,0))*f(:,0)
       m = m + 1
    END IF
    IF (calc(1)) THEN
       ec(:,m) = ec0(:,1)*rs(:,1) &
            + (ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,0) &
            + (ec1(:,0) - ec0(:,0))*f(:,1)*(2._dbl*rhob/(rhoa + rhob)**2)
       ec(:,m+1) = ec0(:,1)*rs(:,1) &
            + (ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,0) &
            + (ec1(:,0) - ec0(:,0))*f(:,1)*(-2._dbl*rhoa/(rhoa + rhob)**2)
       m = m + 2
    END IF
    IF (calc(2)) THEN
       ec(:,m) = ec0(:,1)*rs(:,2) &
            + ec0(:,2)*rs(:,1)**2 &
            + (ec1(:,2) - ec0(:,2))*rs(:,1)**2*f(:,0) &
            + (ec1(:,1) - ec0(:,1))*rs(:,2)*f(:,0) &
            + 2._dbl*(ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,1)*(2._dbl*rhob/(rhoa + rhob)**2) &
            + (ec1(:,0) - ec0(:,0))*f(:,2)*(2._dbl*rhob/(rhoa + rhob)**2)**2 &
            + (ec1(:,0) - ec0(:,0))*f(:,1)*(-4._dbl*rhob/(rhoa + rhob)**3)
       ec(:,m+1) = ec0(:,1)*rs(:,2) &
            + ec0(:,2)*rs(:,1)**2 &
            + (ec1(:,2) - ec0(:,2))*rs(:,1)**2*f(:,0) &
            + (ec1(:,1) - ec0(:,1))*rs(:,2)*f(:,0) &
            + (ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,1)*(2._dbl*(rhob-rhoa)/(rhoa + rhob)**2) &
            + (ec1(:,0) - ec0(:,0))*f(:,2)*(-4._dbl*(rhob*rhoa)/(rhoa + rhob)**4) &
            + (ec1(:,0) - ec0(:,0))*f(:,1)*(2._dbl*(rhoa-rhob)/(rhoa + rhob)**3)
       ec(:,m) = ec0(:,1)*rs(:,2) &
            + ec0(:,2)*rs(:,1)**2 &
            + (ec1(:,2) - ec0(:,2))*rs(:,1)**2*f(:,0) &
            + (ec1(:,1) - ec0(:,1))*rs(:,2)*f(:,0) &
            + 2._dbl*(ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,1)*(-2._dbl*rhoa/(rhoa + rhob)**2) &
            + (ec1(:,0) - ec0(:,0))*f(:,2)*(-2._dbl*rhoa/(rhoa + rhob)**2)**2 &
            + (ec1(:,0) - ec0(:,0))*f(:,1)*(4._dbl*rhoa/(rhoa + rhob)**3)
       m = m + 3
    END IF
    IF (calc(3)) THEN
       ec(:,m) = 0._dbl &
            + ec0(:,2)*rs(:,1)*rs(:,2) &
            + ec0(:,1)*rs(:,3)
       ec(:,m) = ec(:,m) &
            + ec0(:,3)*rs(:,1)**3 &
            + ec0(:,2)*2._dbl*rs(:,1)*rs(:,2) 
       ec(:,m) = ec(:,m) &
            + (ec1(:,3) - ec0(:,3))*rs(:,1)**3*f(:,0) &
            + (ec1(:,2) - ec0(:,2))*2._dbl*rs(:,1)*rs(:,2)*f(:,0) &
            + (ec1(:,2) - ec0(:,2))*rs(:,1)**2*f(:,1)*(2._dbl*rhob/(rhoa + rhob)**2)
       ec(:,m) = ec(:,m) & 
            + (ec1(:,2) - ec0(:,2))*rs(:,1)*rs(:,2)*f(:,0) &
            + (ec1(:,1) - ec0(:,1))*rs(:,3)*f(:,0) &
            + (ec1(:,1) - ec0(:,1))*rs(:,2)*f(:,1)*(2._dbl*rhob/(rhoa + rhob)**2)
       ec(:,m) = ec(:,m) &
            + 2._dbl*(ec1(:,2) - ec0(:,2))*rs(:,1)**2*f(:,1)*(2._dbl*rhob/(rhoa + rhob)**2) &
            + 2._dbl*(ec1(:,1) - ec0(:,1))*rs(:,2)*f(:,1)*(2._dbl*rhob/(rhoa + rhob)**2) &
            + 2._dbl*(ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,2)*(2._dbl*rhob/(rhoa + rhob)**2)**2 &
            + 2._dbl*(ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,1)*(-4._dbl*rhob/(rhoa + rhob)**3)
       ec(:,m) = ec(:,m) &
            + (ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,2)*(2._dbl*rhob/(rhoa + rhob)**2)**2 &
            + (ec1(:,0) - ec0(:,0))*f(:,3)*(2._dbl*rhob/(rhoa + rhob)**2)**3 &
            + (ec1(:,0) - ec0(:,0))*f(:,2)*(-16._dbl*rhob**2/(rhoa + rhob)**5)
       ec(:,m) = ec(:,m) &
            + (ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,1)*(-4._dbl*rhob/(rhoa + rhob)**3) &
            + (ec1(:,0) - ec0(:,0))*f(:,2)*(-8._dbl*rhob**2/(rhoa + rhob)**5) &
            + (ec1(:,0) - ec0(:,0))*f(:,1)*(12._dbl*rhob/(rhoa + rhob)**4)

       ec(:,m+1) = 0._dbl &
            + ec0(:,2)*rs(:,1)*rs(:,2) &
            + ec0(:,1)*rs(:,3)
       ec(:,m+1) = ec(:,m+1) &
            + ec0(:,3)*rs(:,1)**3 &
            + ec0(:,2)*2._dbl*rs(:,1)*rs(:,2) 
       ec(:,m+1) = ec(:,m+1) &
            + (ec1(:,3) - ec0(:,3))*rs(:,1)**3*f(:,0) &
            + (ec1(:,2) - ec0(:,2))*2._dbl*rs(:,1)*rs(:,2)*f(:,0) &
            + (ec1(:,2) - ec0(:,2))*rs(:,1)**2*f(:,1)*(-2._dbl*rhoa/(rhoa + rhob)**2)
       ec(:,m+1) = ec(:,m+1) & 
            + (ec1(:,2) - ec0(:,2))*rs(:,1)*rs(:,2)*f(:,0) &
            + (ec1(:,1) - ec0(:,1))*rs(:,3)*f(:,0) &
            + (ec1(:,1) - ec0(:,1))*rs(:,2)*f(:,1)*(-2._dbl*rhoa/(rhoa + rhob)**2)
       ec(:,m+1) = ec(:,m+1) &
            + (ec1(:,2) - ec0(:,2))*rs(:,1)**2*f(:,1)*(2._dbl*(rhob-rhoa)/(rhoa + rhob)**2) &
            + (ec1(:,1) - ec0(:,1))*rs(:,2)*f(:,1)*(2._dbl*(rhob-rhoa)/(rhoa + rhob)**2) &
            + (ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,2)*(-4._dbl*(rhob-rhoa)*rhoa/(rhoa + rhob)**4) &
            + (ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,1)*(-2._dbl*(3._dbl*rhob-rhoa)/(rhoa + rhob)**3) &
       ec(:,m+1) = ec(:,m+1) &
            + (ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,2)*(2._dbl*rhob/(rhoa + rhob)**2)**2 &
            + (ec1(:,0) - ec0(:,0))*f(:,3)*(-8._dbl*rhoa*rhob**2/(rhoa + rhob)**6) &
            + (ec1(:,0) - ec0(:,0))*f(:,2)*(8._dbl*rhob*(rhoa-rhob)/(rhoa + rhob)**5)
       ec(:,m+1) = ec(:,m+1) &
            + (ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,1)*(-4._dbl*rhob/(rhoa + rhob)**3) &
            + (ec1(:,0) - ec0(:,0))*f(:,2)*(8._dbl*rhoa*rhob/(rhoa + rhob)**5) &
            + (ec1(:,0) - ec0(:,0))*f(:,1)*(-4._dbl*(rhoa-2._dbl*rhob)/(rhoa + rhob)**4)

       ec(:,m+2) = 0._dbl &
            + ec0(:,2)*rs(:,1)*rs(:,2) &
            + ec0(:,1)*rs(:,3)
       ec(:,m+2) = ec(:,m+2) &
            + ec0(:,3)*rs(:,1)**3 &
            + ec0(:,2)*2._dbl*rs(:,1)*rs(:,2) 
       ec(:,m+2) = ec(:,m+2) &
            + (ec1(:,3) - ec0(:,3))*rs(:,1)**3*f(:,0) &
            + (ec1(:,2) - ec0(:,2))*2._dbl*rs(:,1)*rs(:,2)*f(:,0) &
            + (ec1(:,2) - ec0(:,2))*rs(:,1)**2*f(:,1)*(2._dbl*rhob/(rhoa + rhob)**2)
       ec(:,m+2) = ec(:,m+2) & 
            + (ec1(:,2) - ec0(:,2))*rs(:,1)*rs(:,2)*f(:,0) &
            + (ec1(:,1) - ec0(:,1))*rs(:,3)*f(:,0) &
            + (ec1(:,1) - ec0(:,1))*rs(:,2)*f(:,1)*(2._dbl*rhob/(rhoa + rhob)**2)
       ec(:,m+2) = ec(:,m+2) &
            + (ec1(:,2) - ec0(:,2))*rs(:,1)**2*f(:,1)*(2._dbl*(rhob-rhoa)/(rhoa + rhob)**2) &
            + (ec1(:,1) - ec0(:,1))*rs(:,2)*f(:,1)*(2._dbl*(rhob-rhoa)/(rhoa + rhob)**2) &
            + (ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,2)*(4._dbl*(rhob-rhoa)*rhob/(rhoa + rhob)**4) &
            + (ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,1)*(2._dbl*(3._dbl*rhoa-rhob)/(rhoa + rhob)**3) &
       ec(:,m+2) = ec(:,m+2) &
            + (ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,2)*(-2._dbl*rhoa/(rhoa + rhob)**2)**2 &
            + (ec1(:,0) - ec0(:,0))*f(:,3)*(-8._dbl*rhoa*rhob**2/(rhoa + rhob)**6) &
            + (ec1(:,0) - ec0(:,0))*f(:,2)*(-8._dbl*rhoa*(rhoa-rhob)/(rhoa + rhob)**5)
       ec(:,m+2) = ec(:,m+2) &
            + (ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,1)*(4._dbl*rhoa/(rhoa + rhob)**3) &
            + (ec1(:,0) - ec0(:,0))*f(:,2)*(8._dbl*rhob*rhoa/(rhoa + rhob)**5) &
            + (ec1(:,0) - ec0(:,0))*f(:,1)*(4._dbl*(rhob-2._dbl*rhoa)/(rhoa + rhob)**4)

       ec(:,m+3) = 0._dbl &
            + ec0(:,2)*rs(:,1)*rs(:,2) &
            + ec0(:,1)*rs(:,3)
       ec(:,m+3) = ec(:,m+3) &
            + ec0(:,3)*rs(:,1)**3 &
            + ec0(:,2)*2._dbl*rs(:,1)*rs(:,2) 
       ec(:,m+3) = ec(:,m+3) &
            + (ec1(:,3) - ec0(:,3))*rs(:,1)**3*f(:,0) &
            + (ec1(:,2) - ec0(:,2))*2._dbl*rs(:,1)*rs(:,2)*f(:,0) &
            + (ec1(:,2) - ec0(:,2))*rs(:,1)**2*f(:,1)*(-2._dbl*rhoa/(rhoa + rhob)**2)
       ec(:,m+3) = ec(:,m+3) & 
            + (ec1(:,2) - ec0(:,2))*rs(:,1)*rs(:,2)*f(:,0) &
            + (ec1(:,1) - ec0(:,1))*rs(:,3)*f(:,0) &
            + (ec1(:,1) - ec0(:,1))*rs(:,2)*f(:,1)*(-2._dbl*rhoa/(rhoa + rhob)**2)
       ec(:,m+3) = ec(:,m+3) &
            + 2._dbl(ec1(:,2) - ec0(:,2))*rs(:,1)**2*f(:,1)*(-2._dbl*rhoa/(rhoa + rhob)**2) &
            + 2._dbl*(ec1(:,1) - ec0(:,1))*rs(:,2)*f(:,1)*(-2._dbl*rhoa/(rhoa + rhob)**2) &
            + 2._dbl*(ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,2)*(-2._dbl*rhoa/(rhoa + rhob)**2)**2 &
            + 2._dbl*(ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,1)*(4._dbl*rhoa/(rhoa + rhob)**3)
       ec(:,m+3) = ec(:,m+3) &
            + (ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,2)*(-2._dbl*rhoa/(rhoa + rhob)**2)**2 &
            + (ec1(:,0) - ec0(:,0))*f(:,3)*(-2._dbl*rhoa/(rhoa + rhob)**2)**3 &
            + (ec1(:,0) - ec0(:,0))*f(:,2)*(-16._dbl*rhoa**2/(rhoa + rhob)**5)
       ec(:,m+3) = ec(:,m+3) &
            + (ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,1)*(4._dbl*rhoa/(rhoa + rhob)**3) &
            + (ec1(:,0) - ec0(:,0))*f(:,2)*(-8._dbl*rhoa**2/(rhoa + rhob)**5) &
            + (ec1(:,0) - ec0(:,0))*f(:,1)*(-12._dbl*rhoa/(rhoa + rhob)**4)

    END IF

  END SUBROUTINE pz

! *****************************************************************************

  SUBROUTINE pz_init ( cutoff, o_method )
    
    IMPLICIT NONE

    REAL(DBL), INTENT(IN) :: cutoff
    CHARACTER(3), INTENT(IN), OPTIONAL :: o_method
    
    CHARACTER(3) :: method

    CALL set_util(cutoff)

    IF (PRESENT(o_method)) THEN
       method = o_method
    ELSE
       method = 'DMC'
    END IF
    
    initialized = .FALSE.

    SELECT CASE (method)

    CASE DEFAULT
       CALL stop_program("perdew_zunger_functional/pz_init", "unknown method")

       !! insert original pz values

    CASE ('DMC')
       g(0)  = -0.103756_dbl; g(1)  = -0.065951_dbl
       b1(0) = 0.56371_dbl;   b1(1) = 1.11846_dbl
       b2(0) = 0.27358_dbl;   b2(1) = 0.18797_dbl
       A(0)  = 0.031091_dbl;  A(1)  = 0.015545_dbl
       B(0)  = -0.046644_dbl; B(1)  = -0.025599_dbl
       C(0)  = -0.00419_dbl;  C(1)  = -0.00329_dbl
       D(0)  = -0.00983_dbl;  D(1)  = -0.00300_dbl
       
    CASE ('VMC') 
       g(0)  = -0.093662_dbl; g(1)  = -0.055331_dbl
       b1(0) = 0.49453_dbl;   b1(1) = 0.93766_dbl
       b2(0) = 0.25534_dbl;   b2(1) = 0.14829_dbl
       A(0)  = 0.031091_dbl;  A(1)  = 0.015545_dbl
       B(0)  = -0.046644_dbl; B(1)  = -0.025599_dbl
       C(0)  = -0.00884_dbl;  C(1)  = -0.00677_dbl
       D(0)  = -0.00688_dbl;  D(1)  = -0.00093_dbl
       
    END SELECT
          
    initialized = .TRUE.
    
  END SUBROUTINE pz_init

! *****************************************************************************

  SUBROUTINE pz_0 ( rhoa, rhob, rs, f, g0, g1, ec )

    IMPLICIT NONE

    REAL(DBL), DIMENSION(:), INTENT(IN) :: rhoa, rhob, rs
    REAL(DBL), DIMENSION(:,0:), INTENT(IN) :: f, g0, g1
    REAL(DBL), DIMENSION(:,0:), INTENT(OUT) :: ec
    
    INTEGER :: k, np, sp

    REAL(DBL), DIMENSION(2) :: uec ! uninterpolated ec

    np = SIZE(rhoa)
    
    ec(:,0) = g0(:,0) + (g1(:,0) - g0(:,0))*f(:,0)

  END SUBROUTINE pz_0
  
! *****************************************************************************

  SUBROUTINE calc_g(rs, g0, g1, o_order)

    IMPLICIT NONE

    REAL(DBL), DIMENSION(:), INTENT(IN) :: rs
    REAL(DBL), DIMENSION(:,0:), INTENT(OUT) :: g0, g1
    INTEGER, INTENT(IN), OPTIONAL :: o_order

    INTEGER :: k, np, order
    REAL(DBL) :: r, sr

    IF (SIZE(g0,1) /= SIZE(rs) .OR. SIZE(g0,2) /= ABS(order)) &
         CALL stop_program("perdew_zunger_functional/calc_g", "size of g0() to small")
    IF (SIZE(g1,1) /= SIZE(rs) .OR. SIZE(g1,2) /= ABS(order)) &
         CALL stop_program("perdew_zunger_functional/calc_g", "size of g1() to small")

    order = 0
    IF (PRESENT(o_order)) order = o_order
    IF (order > 3 .OR. order < 0) CALL stop_program(&
         "perdew_zunger_functional/calc_g", "order must be between 0 and 3")
    
    np = SIZE(rs)

    DO k=1, np
       r = rs(k); sr = SQRT(r)
       IF (r >= 1.0_dbl) THEN
          g0(k,0) = g(0) / (1.0_dbl + b1(0)*sr + b2(0)*r)
          g1(k,0) = g(1) / (1.0_dbl + b1(1)*sr + b2(1)*r)
       ELSE
          g0(k,0) = A(0)*LOG(r) + B(0) + C(0)*r*LOG(r) + D(0)*r
          g1(k,0) = A(1)*LOG(r) + B(1) + C(1)*r*LOG(r) + D(1)*r
       END IF
    END DO

    IF (order >= 1) THEN
       DO k=1, np
          r = rs(k); sr = SQRT(r)
          IF (r >= 1.0_dbl) THEN
             g0(k,1) = -g(0) * (b1(0)/(2.0_dbl*sr) + b2(0)) / &
                  (1.0_dbl + b1(0)*sr + b2(0)*r)**2
             g1(k,1) = -g(1) * (b1(1)/(2.0_dbl*sr) + b2(1)) / &
                  (1.0_dbl + b1(1)*sr + b2(1)*r)**2
          ELSE
             g0(k,1) = A(0)/r + C(0)*LOG(r) + C(0) + D(0)
             g1(k,1) = A(1)/r + C(1)*LOG(r) + C(1) + D(1)
          END IF
       END DO
    END IF

    IF (order >= 2) THEN
       DO k=1, np
          r = rs(k); sr = SQRT(r)
          IF (rs(k) >= 1.0_dbl) THEN
             g0(k,2) = &
                  2._dbl * g(0) * (b1(0)/(2.0_dbl*sr) + b2(0))**2 / &
                  (1._dbl + b1(0)*sr + b2(0)*r)**3 &
                  + (1._dbl/4._dbl)*g(0)*b1(0) / &
                  ((1._dbl + b1(0)*sr + b2(0)*r)**2 * r**(3./2.))
             g1(k,2) = &
                  2._dbl * g(1) * (b1(1)/(2.0_dbl*sr) + b2(1))**2 / &
                  (1._dbl + b1(1)*sr + b2(1)*r)**3 &
                  + (1._dbl/4._dbl)*g(1)*b1(1) / &
                  ((1._dbl + b1(1)*sr + b2(1)*r)**2 * r**(3./2.))
          ELSE
             g0(k,2) = -A(0)/r**2 + C(0)/r
             g1(k,2) = -A(1)/r**2 + C(1)/r
          END IF
       END DO
    END IF

    IF (order >= 3) THEN
       DO k=1, np
          r = rs(k); sr = SQRT(r)
          IF (rs(k) >= 1.0_dbl) THEN
             g0(k,3) = &
                  -6.0_dbl * g(0) * (b1(0)/(2.0_dbl*sr) + b2(0))**3 / &
                  (1.0_dbl + b1(0)*sr + b2(0)*r)**4 &
                  -(3._dbl/2._dbl) * g(0) * (b1(0)/(2.0_dbl*sr) + b2(0)) * b1(0) / &
                  ((1._dbl + b1(0)*sr + b2(0)*r)**3 * r**(3./2.)) &
                  -(3._dbl/8._dbl) * g(0) * b1(0) / &
                  ((1._dbl + b1(0)*sr + b2(0)*r)**2 * r**(5./2.))
             g1(k,3) = &
                  -6.0_dbl * g(1) * (b1(1)/(2.0_dbl*sr) + b2(1))**3 / &
                  (1.0_dbl + b1(1)*sr + b2(1)*r)**4 &
                  -(3._dbl/2._dbl) * g(1) * (b1(1)/(2.0_dbl*sr) + b2(1)) * b1(1) / &
                  ((1._dbl + b1(1)*sr + b2(1)*r)**3 * r**(3./2.)) &
                  -(3._dbl/8._dbl) * g(1) * b1(1) / &
                  ((1._dbl + b1(1)*sr + b2(1)*r)**2 * r**(5./2.))
          ELSE
             g0(k,3) = A(0)/(2.0_dbl*r**3) - C(0)/r**2
             g0(k,3) = A(1)/(2.0_dbl*r**3) - C(1)/r**2
          END IF
       END DO
    END IF

  END SUBROUTINE calc_g

! *****************************************************************************



! *****************************************************************************

END MODULE PERDEW_ZUNGER_FUNCTIONAL
