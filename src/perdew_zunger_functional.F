!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****m* cp2k/perdew_zunger_funcional *
!!
!!   NAME
!!     perdew_zunger_functional
!!
!!   FUNCTION
!!     Calculate the Perdew-Zunger correlation potential and its
!!     derivatives with respect to the spin-up and spin-down
!!     densities up to 3rd order.
!!
!!   SYNOPSIS
!!     perdew_zunger ( method, cutoff )
!!     perdew_zunger ( rhoa, rhob, ec, order )
!!
!!   MODIFICATION HISTORY
!!     22-FEB-2002, TCH, working version
!!
!!   SEE ALSO
!!     functional_utilities
!!
!!****

MODULE perdew_zunger_functional
  
  USE functionals_utilities, ONLY : calc_fx, calc_rs, set_util
  USE kinds, ONLY : dbl
  USE termination, ONLY : stop_program
  
  IMPLICIT NONE

  LOGICAL :: initialized = .FALSE.
  REAL(DBL), DIMENSION(0:1) :: A = (0,0), B = (0,0), C = (0,0), D = (0,0), &
                               b1 = (0,0), b2 = (0,0), ga = (0,0)

  REAL(DBL), PARAMETER :: epsilon = 5.E-13

  PRIVATE
  
  PUBLIC :: perdew_zunger
  
  INTERFACE perdew_zunger
     MODULE PROCEDURE pz_info, pz_init, pz_lsd, pz_lda
  END INTERFACE

CONTAINS

! *****************************************************************************

!!****f* cp2k/perdew_zunger_functional/pz *
!!
!!   NAME
!!     perdew_zunger
!!
!!   FUNCTION
!!     Does the computation.
!!
!!   ARGUMENTS
!!     
!!
!!   MODIFICATION HISTORY
!!     , TCH, working version
!!
!!   SEE ALSO
!!     
!!
!!****
  SUBROUTINE pz_lsd ( rhoa, rhob, ec, order )
    
    IMPLICIT NONE
    
    REAL(DBL), DIMENSION(:), INTENT(IN) :: rhoa, rhob
    REAL(DBL), DIMENSION(:,0:), INTENT(OUT) :: ec
    INTEGER, INTENT(IN), OPTIONAL :: order

    INTEGER :: m, r_order
    LOGICAL, DIMENSION(0:3) :: calc
    REAL(DBL), DIMENSION(1:SIZE(rhoa),0:3) :: rs
    REAL(DBL), DIMENSION(1:SIZE(rhoa),0:3) :: f, ec0, ec1
    REAL(DBL), DIMENSION(1:SIZE(rhoa),0:3,0:3) :: z
    
#ifdef DEBUG
    integer :: k,n,errno
    n=size(rhoa)
    
    open (UNIT=10, FILE='log', STATUS='REPLACE', ACTION='WRITE', IOSTAT=errno)
#endif

    calc = .FALSE.

    r_order = 0
    IF (PRESENT(order)) r_order = order

    IF (ABS(r_order) <= 3) THEN
       IF (r_order >= 0) THEN
          calc(0:r_order) = .TRUE.
          m = SUM((/ (m+1,m=0,r_order) /))
       ELSE
          r_order = -1 * r_order
          calc(r_order) = .TRUE.
          m = r_order + 1
       END IF
    ELSE
       CALL stop_program("perdew_zunger_functional/pz", &
            "parameter order must be between -3 and 3")
    END IF

    IF (SIZE(ec,2) < m) & 
         CALL stop_program("perdew_zunger_functional/pz", &
         "array ec too small")

    IF (.NOT.initialized) THEN
       CALL pz_init('DMC', epsilon)
    END IF

    CALL calc_fx(rhoa, rhob, f(:,0:r_order), r_order)
    CALL calc_rs(rhoa + rhob, rs(:,0))
    CALL calc_g(rs(:,0), 0, ec0(:,0:r_order), r_order)
    CALL calc_g(rs(:,0), 1, ec1(:,0:r_order), r_order)
    CALL calc_z(rhoa, rhob, z, r_order)

    IF (r_order >= 1) rs(:,1) = (-1._dbl/3._dbl)*rs(:,0)/(rhoa(:) + rhob(:))
    IF (r_order >= 2) rs(:,2) = (-4._dbl/3._dbl)*rs(:,1)/(rhoa(:) + rhob(:))
    IF (r_order >= 3) rs(:,3) = (-7._dbl/3._dbl)*rs(:,2)/(rhoa(:) + rhob(:))

    m = 0
    IF (calc(0)) THEN
       ec(:,m) = ec0(:,0) + (ec1(:,0) - ec0(:,0))*f(:,0)
       m = m + 1
    END IF
    IF (calc(1)) THEN
       ec(:,m) = 0._dbl &
            + ec0(:,1)*rs(:,1) &
            + (ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,0) &
            + (ec1(:,0) - ec0(:,0))*f(:,1)*z(:,1,0)

       ec(:,m+1) = 0._dbl &
            + ec0(:,1)*rs(:,1) &
            + (ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,0) &
            + (ec1(:,0) - ec0(:,0))*f(:,1)*z(:,0,1)

       m = m + 2
    END IF
    IF (calc(2)) THEN
       ec(:,m) = 0._dbl &
            + ec0(:,2)*rs(:,1)**2 &
            + ec0(:,1)*rs(:,2) &
            + (ec1(:,2) - ec0(:,2))*rs(:,1)**2*f(:,0) &
            + (ec1(:,1) - ec0(:,1))*rs(:,2)*f(:,0) &
            + 2._dbl*(ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,1)*z(:,1,0) &
            + (ec1(:,0) - ec0(:,0))*f(:,2)*z(:,1,0)**2 &
            + (ec1(:,0) - ec0(:,0))*f(:,1)*z(:,2,0)

       ec(:,m+1) = 0._dbl &
            + ec0(:,2)*rs(:,1)**2 &
            + ec0(:,1)*rs(:,2) &
            + (ec1(:,2) - ec0(:,2))*rs(:,1)**2*f(:,0) &
            + (ec1(:,1) - ec0(:,1))*rs(:,2)*f(:,0) &
            + (ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,1)*(z(:,1,0)+z(:,0,1)) &
            + (ec1(:,0) - ec0(:,0))*f(:,2)*z(:,1,0)*z(:,0,1) &
            + (ec1(:,0) - ec0(:,0))*f(:,1)*z(:,1,1)

       ec(:,m) = 0._dbl &
            + ec0(:,1)*rs(:,2) &
            + ec0(:,2)*rs(:,1)**2 &
            + (ec1(:,2) - ec0(:,2))*rs(:,1)**2*f(:,0) &
            + (ec1(:,1) - ec0(:,1))*rs(:,2)*f(:,0) &
            + 2._dbl*(ec1(:,1) - ec0(:,1))*rs(:,1)*f(:,1)*z(:,0,1) &
            + (ec1(:,0) - ec0(:,0))*f(:,2)*z(:,0,1)**2 &
            + (ec1(:,0) - ec0(:,0))*f(:,1)*z(:,0,2)

       m = m + 3
    END IF
    IF (calc(3)) THEN
       ec(:,m) = 0._dbl &
            + ec0(:,2)*rs(:,1)*rs(:,2) &
            + ec0(:,1)*rs(:,3) &
            + ec0(:,3)*rs(:,1)**3 &
            + ec0(:,2)*2._dbl*rs(:,1)*rs(:,2)  &
            + (ec1(:,3) - ec0(:,3))*rs(:,1)**3*f(:,0) &
            + (ec1(:,2) - ec0(:,2))*2._dbl*rs(:,1)*rs(:,2)*f(:,0) &
            + (ec1(:,2) - ec0(:,2))*rs(:,1)**2*f(:,1)*z(:,1,0) &
            + (ec1(:,2)-ec0(:,2)) * rs(:,1) * rs(:,2) * f(:,0) &
            + (ec1(:,1)-ec0(:,1)) * rs(:,3) * f(:,0) &
            + (ec1(:,1)-ec0(:,1)) * rs(:,2) * f(:,1) * z(:,1,0) &
            + 2._dbl * (ec1(:,2)-ec0(:,2)) * rs(:,1)**2 * f(:,1) * z(:,1,0) &
            + 2._dbl * (ec1(:,1)-ec0(:,1)) * rs(:,2) * f(:,1) * z(:,1,0) &
            + 2._dbl * (ec1(:,1)-ec0(:,1)) * rs(:,1) * f(:,2) * z(:,1,0)**2 &
            + 2._dbl * (ec1(:,1)-ec0(:,1)) * rs(:,1) * f(:,1) * z(:,2,0) &
            + (ec1(:,1)-ec0(:,1)) * rs(:,1) * f(:,2) * z(:,1,0)**2 &
            + (ec1(:,0)-ec0(:,0)) * f(:,3) * z(:,1,0)**3 &
            + (ec1(:,0)-ec0(:,0)) * f(:,2) * 2._dbl * z(:,1,0)*z(:,2,0) &
            + (ec1(:,1)-ec0(:,1)) * rs(:,1) * f(:,1) * z(:,2,0) &
            + (ec1(:,0)-ec0(:,0)) * f(:,2) * z(:,1,0) * z(:,2,0) &
            + (ec1(:,0)-ec0(:,0)) * f(:,1) * z(:,3,0)

       ec(:,m+1) = 0._dbl &
            + ec0(:,2) * rs(:,1) * rs(:,2) &
            + ec0(:,1) * rs(:,3) &
            + ec0(:,3) * rs(:,1)**3 &
            + ec0(:,2) * 2._dbl * rs(:,1) * rs(:,2) &
            + (ec1(:,3)-ec0(:,3)) * rs(:,1)**3 * f(:,0) &
            + (ec1(:,2)-ec0(:,2)) * 2._dbl * rs(:,1) * rs(:,2) * f(:,0) &
            + (ec1(:,2)-ec0(:,2)) * rs(:,1)**2 * f(:,1) * z(:,0,1) &
            + (ec1(:,2)-ec0(:,2)) * rs(:,1) * rs(:,2) * f(:,0) &
            + (ec1(:,1)-ec0(:,1)) * rs(:,3) * f(:,0) &
            + (ec1(:,1)-ec0(:,1)) * rs(:,2) * f(:,1) * z(:,0,1) &
            + 2._dbl * (ec1(:,2)-ec0(:,2)) * rs(:,1)**2 * f(:,1) * z(:,1,0) &
            + 2._dbl * (ec1(:,1)-ec0(:,1)) * rs(:,2) * f(:,1) * z(:,1,0) &
            + 2._dbl * (ec1(:,1)-ec0(:,1)) * rs(:,1) * f(:,2) * z(:,0,1) * z(:,1,0) &
            + 2._dbl * (ec1(:,1)-ec0(:,1)) * rs(:,1) * f(:,1) * z(:,1,1) &
            + (ec1(:,1)-ec0(:,1)) * rs(:,1) * f(:,2) * z(:,1,0)**2 &
            + (ec1(:,0)-ec0(:,0)) * f(:,3) * z(:,0,1) * z(:,1,0)**2 &
            + (ec1(:,0)-ec0(:,0)) * f(:,2) * 2._dbl * z(:,1,0) * z(:,1,1) &
            + (ec1(:,1)-ec0(:,1)) * rs(:,1) * f(:,1) * z(:,2,0) &
            + (ec1(:,0)-ec0(:,0)) * f(:,2) * z(:,0,1) * z(:,2,0) &
            + (ec1(:,0)-ec0(:,0)) * f(:,1) * z(:,2,1)

       ec(:,m+2) = 0._dbl &
            + ec0(:,2) * rs(:,1) * rs(:,2) &
            + ec0(:,1) * rs(:,3) &
            + ec0(:,3) * rs(:,1)**3 &
            + ec0(:,2) * 2._dbl * rs(:,1) * rs(:,2) &
            + (ec1(:,3)-ec0(:,3)) * rs(:,1)**3 * f(:,0) &
            + (ec1(:,2)-ec0(:,2)) * 2._dbl * rs(:,1) * rs(:,2) * f(:,0) &
            + (ec1(:,2)-ec0(:,2)) * rs(:,1)**2 * f(:,1) * z(:,1,0) &
            + (ec1(:,2)-ec0(:,2)) * rs(:,1) * rs(:,2) * f(:,0) &
            + (ec1(:,1)-ec0(:,1)) * rs(:,3) * f(:,0) &
            + (ec1(:,1)-ec0(:,1)) * rs(:,2) * f(:,1) * z(:,1,0) &
            + 2._dbl * (ec1(:,2)-ec0(:,2)) * rs(:,1)**2 * f(:,1) * z(:,0,1) &
            + 2._dbl * (ec1(:,1)-ec0(:,1)) * rs(:,2) * f(:,1) * z(:,0,1) &
            + 2._dbl * (ec1(:,1)-ec0(:,1)) * rs(:,1) * f(:,2) * z(:,1,0) * z(:,0,1) &
            + 2._dbl * (ec1(:,1)-ec0(:,1)) * rs(:,1) * f(:,1) * z(:,1,1) &
            + (ec1(:,1)-ec0(:,1)) * rs(:,1) * f(:,2) * z(:,0,1)**2 &
            + (ec1(:,0)-ec0(:,0)) * f(:,3) * z(:,1,0) * z(:,0,1)**2 &
            + (ec1(:,0)-ec0(:,0)) * f(:,2) * 2._dbl * z(:,0,1) * z(:,1,1) &
            + (ec1(:,1)-ec0(:,1)) * rs(:,1) * f(:,1) * z(:,0,2) &
            + (ec1(:,0)-ec0(:,0)) * f(:,2) * z(:,1,0) * z(:,0,2) &
            + (ec1(:,0)-ec0(:,0)) * f(:,1) * z(:,1,2)

       ec(:,m+3) = 0._dbl &
            + ec0(:,2) * rs(:,1) * rs(:,2) &
            + ec0(:,1) * rs(:,3) &
            + ec0(:,3) * rs(:,1)**3 &
            + ec0(:,2) * 2._dbl * rs(:,1) * rs(:,2)  &
            + (ec1(:,3)-ec0(:,3)) * rs(:,1)**3 * f(:,0) &
            + (ec1(:,2)-ec0(:,2)) * 2._dbl * rs(:,1) * rs(:,2) * f(:,0) &
            + (ec1(:,2)-ec0(:,2)) * rs(:,1)**2 * f(:,1) * z(:,0,1) &
            + (ec1(:,2)-ec0(:,2)) * rs(:,1) * rs(:,2) * f(:,0) &
            + (ec1(:,1)-ec0(:,1)) * rs(:,3) * f(:,0) &
            + (ec1(:,1)-ec0(:,1)) * rs(:,2) * f(:,1) * z(:,0,1) &
            + 2._dbl * (ec1(:,2)-ec0(:,2)) * rs(:,1)**2 * f(:,1) * z(:,0,1) &
            + 2._dbl * (ec1(:,1)-ec0(:,1)) * rs(:,2) * f(:,1) * z(:,0,1) &
            + 2._dbl * (ec1(:,1)-ec0(:,1)) * rs(:,1) * f(:,2) * z(:,0,1)**2 &
            + 2._dbl * (ec1(:,1)-ec0(:,1)) * rs(:,1) * f(:,1) * z(:,0,2) &
            + (ec1(:,1)-ec0(:,1)) * rs(:,1) * f(:,2) * z(:,0,1)**2 &
            + (ec1(:,0)-ec0(:,0)) * f(:,3) * z(:,0,1)**3 &
            + (ec1(:,0)-ec0(:,0)) * f(:,2) * 2._dbl * z(:,0,1) * z(:,0,2) &
            + (ec1(:,1)-ec0(:,1)) * rs(:,1)* f(:,1) * z(:,0,2) &
            + (ec1(:,0)-ec0(:,0)) * f(:,2) * z(:,0,1) * z(:,0,2) &
            + (ec1(:,0)-ec0(:,0)) * f(:,1) * z(:,0,3)
    END IF

#ifdef DEBUG_PZ_LSD
!!!    DO k=1, n
!!!       WRITE (*, '(1X,ES12.5,4(4X,ES12.5))') rs(k,0), ec0(k,0), ec0(k,1), ec0(k,2), ec0(k,3)
!!!    END DO
#endif    

#ifdef DEBUG
    close(10)
#endif

  END SUBROUTINE pz_lsd

! *****************************************************************************

  SUBROUTINE pz_lda ( rho, ec, order )
    
    IMPLICIT NONE
    
    REAL(DBL), DIMENSION(:), INTENT(IN) :: rho
    REAL(DBL), DIMENSION(:,0:), INTENT(OUT) :: ec
    INTEGER, INTENT(IN), OPTIONAL :: order

    INTEGER :: m, r_order
    LOGICAL, DIMENSION(0:3) :: calc
    REAL(DBL), DIMENSION(1:SIZE(rho),0:3) :: rs
    REAL(DBL), DIMENSION(1:SIZE(rho),0:3) :: ec0

#ifdef DEBUG
    integer :: k,n,errno
    n=size(rho)
    
    open (UNIT=10, FILE='log', STATUS='REPLACE', ACTION='WRITE', IOSTAT=errno)
#endif

    calc = .FALSE.

    r_order = 0
    IF (PRESENT(order)) r_order = order

    IF (ABS(r_order) <= 3) THEN
       IF (r_order >= 0) THEN
          calc(0:r_order) = .TRUE.
          m = r_order + 1
       ELSE
          r_order = -1 * r_order
          calc(r_order) = .TRUE.
          m = 1
       END IF
    ELSE
       CALL stop_program("perdew_zunger_functional/pz", &
            "parameter order must be between -3 and 3")
    END IF

    IF (SIZE(ec,2) < m) & 
         CALL stop_program("perdew_zunger_functional/pz", &
         "array ec too small")

    IF (.NOT.initialized) THEN
       CALL pz_init('DMC', epsilon)
    END IF

    CALL calc_rs(rho, rs(:,0))
    CALL calc_g(rs(:,0), 0, ec0, r_order)

    IF (r_order >= 1) rs(:,1) = (-1._dbl/3._dbl) * rs(:,0) / rho(:)
    IF (r_order >= 2) rs(:,2) = (-4._dbl/3._dbl) * rs(:,1) / rho(:)
    IF (r_order >= 3) rs(:,3) = (-7._dbl/3._dbl) * rs(:,2) / rho(:) 

    m = 0
    IF (calc(0)) THEN
       ec(:,m) = ec0(:,0)
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

#ifdef DEBUG_PZ_LDA
    do k=1, n
       write (10,'(1X,ES12.5,5(4X,ES12.5))') rho(k), rs(k,0), ec(k,0), ec(k,1), ec0(k,1), rs(k,1)
    end do
#endif

#ifdef DEBUG
    close(10)
#endif

  END SUBROUTINE pz_lda

! *****************************************************************************

  SUBROUTINE pz_init ( method, cutoff )
    
    IMPLICIT NONE

    CHARACTER(3), INTENT(IN) :: method
    REAL(DBL), INTENT(IN) :: cutoff
    
    CALL set_util(cutoff)
    
    initialized = .FALSE.

    SELECT CASE (method)

    CASE DEFAULT
       CALL stop_program("perdew_zunger_functional/pz_init", "unknown method")
       
    CASE ('ORIG')
       ga(0)  = -0.1423_dbl;  ga(1)  = -0.0843_dbl
       b1(0) = 1.0529_dbl;    b1(1) = 1.3981_dbl
       b2(0) = 0.3334_dbl;    b2(1) = 0.2611_dbl
       A(0)  = 0.0311_dbl;    A(1)  = 0.01555_dbl
       B(0)  = -0.048_dbl;    B(1)  = -0.0269_dbl
       C(0)  = -0.0020_dbl;   C(1)  = -0.0007_dbl
       D(0)  = -0.0116_dbl;   D(1)  = -0.0048_dbl

    CASE ('DMC')
       ga(0)  = -0.103756_dbl; ga(1)  = -0.065951_dbl
       b1(0) = 0.56371_dbl;   b1(1) = 1.11846_dbl
       b2(0) = 0.27358_dbl;   b2(1) = 0.18797_dbl
       A(0)  = 0.031091_dbl;  A(1)  = 0.015545_dbl
       B(0)  = -0.046644_dbl; B(1)  = -0.025599_dbl
       C(0)  = -0.00419_dbl;  C(1)  = -0.00329_dbl
       D(0)  = -0.00983_dbl;  D(1)  = -0.00300_dbl
       
    CASE ('VMC') 
       ga(0)  = -0.093662_dbl; ga(1)  = -0.055331_dbl
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

  SUBROUTINE calc_g(rs, z, g, order)

    IMPLICIT NONE

    REAL(DBL), DIMENSION(:,0:), INTENT(OUT) :: g
    REAL(DBL), DIMENSION(:), INTENT(IN) :: rs
    INTEGER, INTENT(IN) :: z, order

    INTEGER :: k, np
    REAL(DBL) :: r, sr

    IF (SIZE(g,1) /= SIZE(rs) .OR. SIZE(g,2) < ABS(order)) &
         CALL stop_program("perdew_zunger_functional/calc_g", "size of g0() to small")

    IF (order > 3 .OR. order < 0) CALL stop_program(&
         "perdew_zunger_functional/calc_g", "order must be between 0 and 3")
    
    np = SIZE(rs)

    DO k=1, np
       r = rs(k)
       IF (r >= 1.0_dbl) THEN
          sr = SQRT(r)
          g(k,0) = ga(z) / (1.0_dbl + b1(z)*sr + b2(z)*r)
       ELSE
          g(k,0) = A(z)*LOG(r) + B(z) + C(z)*r*LOG(r) + D(z)*r
       END IF
    END DO

    IF (order >= 1) THEN
       DO k=1, np
          r = rs(k)
          IF (r >= 1.0_dbl) THEN
             sr = SQRT(r)
             g(k,1) = -ga(z) * (b1(z)/(2.0_dbl*sr) + b2(z)) / &
                  (1.0_dbl + b1(z)*sr + b2(z)*r)**2
          ELSE
             g(k,1) = A(z)/r + C(z)*LOG(r) + C(z) + D(z)
          END IF
       END DO
    END IF

    IF (order >= 2) THEN
       DO k=1, np
          r = rs(k)
          IF (rs(k) >= 1.0_dbl) THEN
             sr = SQRT(r)
             g(k,2) = &
                  2._dbl * ga(z) * (b1(z)/(2.0_dbl*sr) + b2(z))**2 &
                  / (1._dbl + b1(z)*sr + b2(z)*r)**3 &
                  + ga(z) * b1(z) &
                  / (4.0_dbl * (1._dbl + b1(z)*sr + b2(z)*r)**2 * sr**3)
          ELSE
             g(k,2) = -A(z)/r**2 + C(z)/r
          END IF
       END DO
    END IF

    IF (order >= 3) THEN
       DO k=1, np
          r = rs(k)
          IF (rs(k) >= 1.0_dbl) THEN
             sr = SQRT(r)
             g(k,3) = &
                  -6.0_dbl * ga(z) * (b1(z)/(2.0_dbl*sr) + b2(z))**3 / &
                  (1.0_dbl + b1(z)*sr + b2(z)*r)**4 &
                  -(3._dbl/2._dbl) * ga(z) * (b1(z)/(2.0_dbl*sr) + b2(z)) * b1(z) / &
                  ((1._dbl + b1(z)*sr + b2(z)*r)**3 * sr**3) &
                  -(3._dbl/8._dbl) * ga(z) * b1(z) / &
                  ((1._dbl + b1(z)*sr + b2(z)*r)**2 * sr**5)
          ELSE
             g(k,3) = 2.0_dbl*A(z)/r**3 - C(z)/r**2
          END IF
       END DO
    END IF

  END SUBROUTINE calc_g

! *****************************************************************************

  SUBROUTINE calc_z ( a, b, z, order )

    IMPLICIT NONE

    REAL(DBL), DIMENSION(:), INTENT(IN) :: a, b
    REAL(DBL), DIMENSION(:,0:,0:), INTENT(OUT) :: z
    INTEGER, INTENT(IN) :: order

    REAL(DBL), DIMENSION(1:SIZE(a)) :: c,d

    c = a+b

    z(:,0,0) = (a-b)/c
    IF (order >= 1) THEN
       d = c**2
       z(:,1,0) = 2._dbl*b/d
       z(:,0,1) = -2._dbl*a/d
    END IF
    IF (order >= 2) THEN
       d = c**3
       z(:,2,0) = -4._dbl*b/d
       z(:,1,1) = 2._dbl*(a-b)/d
       z(:,0,2) = 4._dbl*a/d
    END IF
    IF (order >= 3) THEN
       d = c**4
       z(:,3,0) = 12._dbl*b/d
       z(:,2,1) = -4._dbl*(a-2._dbl*b)/d
       z(:,1,2) = -4._dbl*(2._dbl*a-b)/d
       z(:,0,3) = -12._dbl*a/d
    END IF

  END SUBROUTINE calc_z
    
! *****************************************************************************

  SUBROUTINE pz_info ( reference, shortform )

    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: reference
    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: shortform

    IF ( PRESENT ( reference ) ) THEN
       reference = "J. P. Perdew and Alex Zunger," &
                   //" Phys. Rev. B 23, 5048 (1981)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "J. P. Perdew et al., PRB 23, 5048 (1981)"
    END IF

  END SUBROUTINE pz_info

! *****************************************************************************

END MODULE PERDEW_ZUNGER_FUNCTIONAL
