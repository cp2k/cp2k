!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****m* cp2k/xc_perdew_zunger *
!!
!!   NAME
!!     xc_perdew_zunger
!!
!!   FUNCTION
!!     Calculate the Perdew-Zunger correlation potential and 
!!     energy density and ist derivatives with respect to
!!     the spin-up and spin-down densities up to 3rd order.
!!
!!   SYNOPSIS
!!     perdew_zunger ( method, cutoff )
!!     perdew_zunger ( reference, shortform )
!!     perdew_zunger ( rho, ec, order )
!!     perdew_zunger ( rhoa, rhob, ec, order )
!!     perdew_zunger_energy_density ( rho, ec, order )
!!     perdew_zunger_energy_density ( rhoa, rhob, ec, order )
!!
!!   MODIFICATION HISTORY
!!     18-MAR-2002, TCH, working version
!!
!!   SEE ALSO
!!     functionals_utilities
!!
!!****

MODULE xc_perdew_zunger
  USE kinds,                           ONLY: dbl
  USE mathconstants,                   ONLY: zero
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type
  USE xc_functionals_utilities,        ONLY: calc_fx,&
                                             calc_rs,&
                                             calc_z,&
                                             set_util

  IMPLICIT NONE

  LOGICAL :: initialized = .FALSE.
  REAL(DBL), DIMENSION(0:1) :: A = (0,0), B = (0,0), C = (0,0), D = (0,0), &
                               b1 = (0,0), b2 = (0,0), ga = (0,0)

  REAL(DBL), PARAMETER :: epsilon = 5.E-13
  REAL(DBL) :: eps_rho

  PRIVATE
  
  PUBLIC :: perdew_zunger, perdew_zunger_energy_density
  
  INTERFACE perdew_zunger
     MODULE PROCEDURE pz_info, pz_init, pz_lsd, pz_lda
  END INTERFACE

  INTERFACE perdew_zunger_energy_density
     MODULE PROCEDURE pz_info, pz_init, pz_lsd_ed, pz_lda_ed
  END INTERFACE

CONTAINS

! *****************************************************************************

!!****f* cp2k/xc_perdew_zunger/pz_info *
!!
!!   NAME
!!     perdew_zunger
!!
!!   FUNCTION
!!     Return some info on the functionals.
!!
!!   ARGUMENTS
!!     reference - CHARACTER(*), INTENT(INOUT), OPTIONAL - full reference
!!     shortform - CHARACTER(*), INTENT(INOUT), OPTIONAL - short reference
!!
!!   MODIFICATION HISTORY
!!     18-MAR-2002, TCH, working version
!!
!!   SEE ALSO
!!     
!!
!!****
  SUBROUTINE pz_info ( reference, shortform )

    CHARACTER(LEN=*), INTENT(INOUT), OPTIONAL  :: reference, shortform

    IF ( PRESENT ( reference ) ) THEN
       reference = "J. P. Perdew and Alex Zunger," &
                   //" Phys. Rev. B 23, 5048 (1981)"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "J. P. Perdew et al., PRB 23, 5048 (1981)"
    END IF

  END SUBROUTINE pz_info

! *****************************************************************************

!!****f* cp2k/xc_perdew_zunger/pz_lda *
!!
!!   NAME
!!     perdew_zunger
!!
!!   FUNCTION
!!     Calculate the correlation energy and its derivatives
!!     wrt to rho (the electron density) up to 3rd order. This
!!     is the LDA version of the Perdew-Zunger correlation energy
!!     If no order argument is given, then the routine calculates
!!     just the energy.
!!
!!   ARGUMENTS
!!     rho - REAL(DBL), DIMENSION(1:n) - the electron density
!!     ec - REAL(DBL), DIMENSION(1:n,0:) - the potential and derivatives
!!     order - INTEGER, OPTIONAL - order of derivatives to calculate
!!
!!     order must lie between -3 and 3. If it is negative then only
!!     that order will be calculated, otherwise all derivatives up to
!!     that order will be calculated.
!!     ec must be large enough to accomodate the derivatives you ask
!!     for. If order is 2 then SIZE(ec,2) must be at least3 (0,1 and 2).
!!     If order is -3 then SIZE(ec,2) must be at least 1.
!!
!!   MODIFICATION HISTORY
!!     18-MAR-2002, TCH, working version
!!
!!   SEE ALSO
!!     functionals_utilities
!!
!!****
  SUBROUTINE pz_lda ( derivative_set, rho, ec, order )
    
    
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(DBL), DIMENSION(:), INTENT(IN)      :: rho
    REAL(DBL), DIMENSION(:, 0:), INTENT(INOUT) :: ec
    INTEGER, INTENT(IN), OPTIONAL            :: order

    INTEGER                                  :: k, m, n, order_, timer_handle
    LOGICAL, DIMENSION(0:3)                  :: calc
    REAL(DBL), DIMENSION(0:3)                :: ed

    CALL timeset('perdew_zunger lda', 'I', ' ', timer_handle)

    n = SIZE(rho)
    
    calc = .FALSE.

    order_ = 0
    IF (PRESENT(order)) order_ = order

    IF (ABS(order_) <= 3) THEN
       IF (order_ >= 0) THEN
          calc(0:order_) = .TRUE.
          m = order_ + 1
       ELSE
          order_ = -1 * order_
          calc(order_) = .TRUE.
          m = 1
       END IF
    ELSE
       CALL stop_program("xc_perdew_zunger/pz_lda", &
            "parameter order must be between -3 and 3")
    END IF

    IF (SIZE(ec,1) < n .OR. SIZE(ec,2) < m) & 
         CALL stop_program("xc_perdew_zunger/pz_lda", &
         "array ec too small")

    IF (.NOT.initialized) THEN
       CALL pz_init('ORIG', epsilon)
    END IF

!$omp parallel do private ( k, ed, m ) shared(order_)
    DO k=1, n
       
       IF ( rho(k) > eps_rho ) THEN

         CALL pz_lda_ed_loc(rho(k), ed, order_)
       
         m = 0
         IF (calc(0)) THEN
            ec(k,m) = ec(k,m) + rho(k)*ed(0)
            m = m + 1
         END IF
         IF (calc(1)) THEN
            ec(k,m) = ec(k,m) + ed(0) + rho(k)*ed(1)
            m =  m + 1
         END IF
         IF (calc(2)) THEN
            ec(k,m) = ec(k,m) + 2._dbl*ed(1) + rho(k)*ed(2)
            m = m + 1
         END IF
         IF (calc(3)) THEN
            ec(k,m) = ec(k,m) + 3._dbl*ed(2) + rho(k)*ed(3)
         END IF

       END IF

    END DO
!$omp end parallel do

    CALL timestop(zero, timer_handle)

  END SUBROUTINE pz_lda

! *****************************************************************************

!!****f* cp2k/xc_perdew_zunger/pw_lda_ed *
!!
!!   NAME
!!     perdew_zunger_energy_density
!!
!!   FUNCTION
!!     Calculate the correlation energy density and its derivatives
!!     wrt to rho (the electron density) up to 3rd order. This
!!     is the LDA version of the Perdew-Zunger correlation energy density.
!!     If no order argument is given, then the routine calculates just
!!     the 0th order
!!
!!   ARGUMENTS
!!     rho - REAL(DBL), DIMENSION(1:n) - the electron density
!!     ed - REAL(DBL), DIMENSION(1:n,0:) - the energy density and derivatives
!!     order - INTEGER, OPTIONAL - order of derivatives to calculate
!!
!!     order must lie between -3 and 3. If it is negative then only
!!     that order will be calculated, otherwise all derivatives up to
!!     that order will be calculated.
!!     ec must be large enough to accomodate the derivatives you ask
!!     for. If order is 2 then SIZE(ec,2) must be at least3 (0,1 and 2).
!!     If order is -3 then SIZE(ec,2) must be at least 1.
!!
!!   MODIFICATION HISTORY
!!     04-APR-2002, TCH, working version
!!
!!   SEE ALSO
!!     functionals_utilities
!!
!!****
  SUBROUTINE pz_lda_ed(rho, ed, order)


    REAL(DBL), DIMENSION(1:), INTENT(IN)     :: rho
    REAL(DBL), DIMENSION(1:, 0:), &
      INTENT(INOUT)                            :: ed
    INTEGER, INTENT(IN), OPTIONAL            :: order

    INTEGER                                  :: k, m, n, order_ = 0

    n = SIZE(rho)
    IF (PRESENT(order)) order_ = order

    IF (ABS(order_) <= 3) THEN
       IF (order_ >= 0) THEN
          m = order_ + 1
       ELSE
          m = 1
       END IF
    ELSE
       CALL stop_program("xc_perdew_zunger/pz_lda_ed", &
            "parameter order must be between -3 and 3")
    END IF

    IF (SIZE(ed,1) < n) &
         CALL stop_program("xc_perdew_zunger/pz_lda_ed", "SIZE(ec,1) too small")
    IF (SIZE(ed,2) < m) &
         CALL stop_program("xc_perdew_zunger/pz_lda_ed", "SIZE(ec,2) too small")
    
    IF (.NOT.initialized) CALL pz_init('ORIG', epsilon)

!$omp parallel do private ( k ) shared(order_)
    DO k=1, n
       IF ( rho(k) > eps_rho ) THEN
          CALL pz_lda_ed_loc(rho(k), ed(k,:), order_)
       END IF
    END DO

  END SUBROUTINE pz_lda_ed
  
! *****************************************************************************

!!****f* cp2k/xc_perdew_zunger/pz_lsd *
!!
!!   NAME
!!     perdew_zunger
!!
!!   FUNCTION
!!     Calculate the correlation energy and its derivatives
!!     wrt to rho (the electron density) up to 3rd order. This
!!     is the LSD version of the Perdew-Zunger correlation energy
!!     If no order argument is given, then the routine calculates
!!     just the energy.
!!
!!   ARGUMENTS
!!     rhoa - REAL(DBL), DIMENSION(1:n) - the density of spin up electrons
!!     rhob - REAL(DBL), DIMENSION(1:n) - the density of spin up electrons
!!     ec - REAL(DBL), DIMENSION(1:n,0:) - the potential and derivatives
!!     order - INTEGER, OPTIONAL - order of derivatives to calculate
!!
!!     order must lie between -3 and 3. If it is negative then only
!!     that order will be calculated, otherwise all derivatives up to
!!     that order will be calculated.
!!     ec must be large enough to accomodate the derivatives you ask
!!     for. If order is 2 then SIZE(ec,2) must be at least 6 (1 for 0th
!!     order, 2 for 1st order and 3 for 2nd order)
!!     If order is -3 then SIZE(ec,2) must be at least 4 (derivatives wrt
!!     to aaa, aab, abb and bbb)
!!
!!   MODIFICATION HISTORY
!!     18-MAR-2002, TCH, working version
!!
!!   SEE ALSO
!!     functionals_utilities
!!
!!****
  SUBROUTINE pz_lsd ( derivative_set, a, b, ec, order )
    
    
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    REAL(DBL), DIMENSION(:), INTENT(IN)      :: a, b
    REAL(DBL), DIMENSION(:, 0:), INTENT(INOUT) :: ec
    INTEGER, INTENT(IN), OPTIONAL            :: order

    INTEGER                                  :: k, m, n, order_, timer_handle
    LOGICAL, DIMENSION(0:3)                  :: calc
    REAL(DBL)                                :: rho
    REAL(DBL), DIMENSION(0:9)                :: ed

    CALL timeset('perdew_zunger lsd', 'I', ' ', timer_handle)

    n = SIZE(a)
    
    calc = .FALSE.

    order_ = 0
    IF (PRESENT(order)) order_ = order

    IF (ABS(order_) <= 3) THEN
       IF (order_ >= 0) THEN
          calc(0:order_) = .TRUE.
          m = SUM((/ (m+1,m=0,order_) /))
       ELSE
          order_ = -1 * order_
          calc(order_) = .TRUE.
          m = order_ + 1
       END IF
    ELSE
       CALL stop_program("xc_perdew_zunger/pz_lsd", &
            "parameter order must be between -3 and 3")
    END IF

    IF (SIZE(ec,1) < n .OR. SIZE(ec,2) < m) & 
         CALL stop_program("xc_perdew_zunger/pz_lsd", &
         "array ec too small")

    IF (.NOT.initialized) THEN
       CALL pz_init('ORIG', epsilon)
    END IF

!$omp parallel do private ( k, rho, ed, m ) shared(order_)
    DO k=1, n

       rho = a(k) + b(k)

       IF ( rho > eps_rho ) THEN

          CALL pz_lsd_ed_loc(a(k), b(k), ed, order_)

          m=0
          IF (calc(0)) THEN
             ec(k,m) = ec(k,m) + rho*ed(0)
             m = m + 1
          END IF

          IF (calc(1)) THEN
             ec(k,m)   = ec(k,m) + ed(0) + rho*ed(1)
             ec(k,m+1) = ec(k,m+1) + ed(0) + rho*ed(2)
             m = m + 2
          END IF
       
          IF (calc(2)) THEN
             ec(k,m)   = ec(k,m) + 2._dbl*ed(1) + rho*ed(3)
             ec(k,m+1) = ec(k,m+1) + ed(1) + ed(2) + rho*ed(4)
             ec(k,m+2) = ec(k,m+2) + 2._dbl*ed(2) + rho*ed(5)
             m = m + 3
          END IF
       
          IF (calc(3)) THEN
             ec(k,m)   = ec(k,m) + 3._dbl*ed(3) + rho*ed(6)
             ec(k,m+1) = ec(k,m+1) + 2._dbl*ed(4) + ed(3) + rho*ed(7)
             ec(k,m+2) = ec(k,m+2) + 2._dbl*ed(4) + ed(5) + rho*ed(8)
             ec(k,m+3) = ec(k,m+3) + 3._dbl*ed(5) + rho*ed(9)
          END IF

       END IF

    END DO
!$omp end parallel do
       
    CALL timestop(zero, timer_handle)

  END SUBROUTINE pz_lsd

! *****************************************************************************

!!****f* cp2k/xc_perdew_zunger/pw_lsd_ed *
!!
!!   NAME
!!     perdew_zunger_energy_density
!!
!!   FUNCTION
!!     Calculate the correlation energy density and its derivatives
!!     wrt to rho (the electron density) up to 3rd order. This
!!     is the LSD version of the Perdew-Zunger correlation energy density.
!!     If no order argument is given, then the routine calculates
!!     just the 0th order
!!
!!   ARGUMENTS
!!     rhoa - REAL(DBL), DIMENSION(1:n) - the density of spin up electrons
!!     rhob - REAL(DBL), DIMENSION(1:n) - the density of spin up electrons
!!     ed - REAL(DBL), DIMENSION(1:n,0:) - the energy density and derivatives
!!     order - INTEGER, OPTIONAL - order of derivatives to calculate
!!
!!     order must lie between -3 and 3. If it is negative then only
!!     that order will be calculated, otherwise all derivatives up to
!!     that order will be calculated.
!!     ec must be large enough to accomodate the derivatives you ask
!!     for. If order is 2 then SIZE(ec,2) must be at least 6 (1 for 0th
!!     order, 2 for 1st order and 3 for 2nd order)
!!     If order is -3 then SIZE(ec,2) must be at least 4 (derivatives wrt
!!     to aaa, aab, abb and bbb)
!!
!!   MODIFICATION HISTORY
!!     04-APR-2002, TCH, working version
!!
!!   SEE ALSO
!!     functionals_utilities
!!
!!****
  SUBROUTINE pz_lsd_ed(rhoa, rhob, ed, order)


    REAL(DBL), DIMENSION(1:), INTENT(IN)     :: rhoa, rhob
    REAL(DBL), DIMENSION(1:, 0:), &
      INTENT(INOUT)                            :: ed
    INTEGER, INTENT(IN), OPTIONAL            :: order

    INTEGER                                  :: k, m, n, order_ = 0

    IF (PRESENT(order)) order_ = order

    IF (ABS(order_) <= 3) THEN
       IF (order_ > 0) THEN
          m = SUM((/ (k+1,k=0,order_) /))
       ELSE
          m = 1
       END IF
    ELSE
       CALL stop_program("xc_perdew_zunger/pz", &
            "parameter order must be between -3 and 3")
    END IF

    n = SIZE(rhoa)
    IF (SIZE(ed,1) < n) &
         CALL stop_program("xc_perdew_zunger/pz_lsd_ed", "SIZE(ed,1) too small")
    IF (SIZE(ed,2) < m) &
         CALL stop_program("xc_perdew_zunger/pz_lsd_ed", "SIZE(ed,2) too small")

    IF (.NOT.initialized) CALL pz_init('ORIG', epsilon)

!$omp parallel do private ( k ) shared(order_)
    DO k=1, n
       IF ( rhoa(k) + rhob(k) > eps_rho ) THEN
          CALL pz_lsd_ed_loc(rhoa(k), rhob(k), ed(k,:), order_)
       END IF
    END DO

  END SUBROUTINE pz_lsd_ed

! *****************************************************************************

!!****f* cp2k/xc_perdew_zunger/pz_init *
!!
!!   NAME
!!     perdew_zunger
!!
!!   FUNCTION
!!     Initializes the functionals
!!
!!   ARGUMENTS
!!     method - CHARACTER(3) - name of the method used for parameters
!!     cutoff - REAL(DBL) - the cutoff density
!!
!!   MODIFICATION HISTORY
!!     18-MAR-2002, TCH, working version
!!
!!   SEE ALSO
!!     functionals_utilities
!!
!!****
  SUBROUTINE pz_init ( method, cutoff )
    

    CHARACTER(*), INTENT(IN)                 :: method
    REAL(DBL), INTENT(IN)                    :: cutoff

    CALL set_util(cutoff)

    eps_rho = cutoff
    
    initialized = .FALSE.

    SELECT CASE (method)

    CASE DEFAULT
       CALL stop_program("xc_perdew_zunger/pz_init", "unknown method")
       
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

  SUBROUTINE calc_g(r, z, g, order)


    REAL(DBL), INTENT(IN)                    :: r
    INTEGER, INTENT(IN)                      :: z
    REAL(DBL), DIMENSION(0:), INTENT(INOUT)    :: g
    INTEGER, INTENT(IN)                      :: order

    REAL(DBL)                                :: rr, rsr, sr

    IF (r >= 1.0_dbl) THEN

       sr = SQRT(r)
       ! order 0 must always be calculated
       g(0) = ga(z) / (1.0_dbl + b1(z)*sr + b2(z)*r)
       IF (order >= 1) THEN
          g(1) = -ga(z) * (b1(z)/(2.0_dbl*sr) + b2(z)) / &
               (1.0_dbl + b1(z)*sr + b2(z)*r)**2
       END IF
       IF (order >= 2) THEN
          rsr = r*sr
          g(2) = &
               2._dbl * ga(z) * (b1(z)/(2.0_dbl*sr) + b2(z))**2 &
               / (1._dbl + b1(z)*sr + b2(z)*r)**3 &
               + ga(z) * b1(z) &
               / (4.0_dbl * (1._dbl + b1(z)*sr + b2(z)*r)**2 * rsr)
       END IF
       IF (order >= 3) THEN
          g(3) = &
               -6.0_dbl * ga(z) * (b1(z)/(2.0_dbl*sr) + b2(z))**3 / &
               (1.0_dbl + b1(z)*sr + b2(z)*r)**4 &
               -(3._dbl/2._dbl) * ga(z) * (b1(z)/(2.0_dbl*sr) + b2(z)) * b1(z) / &
               ((1._dbl + b1(z)*sr + b2(z)*r)**3 * rsr) &
               -(3._dbl/8._dbl) * ga(z) * b1(z) / &
               ((1._dbl + b1(z)*sr + b2(z)*r)**2 * r*rsr)
       END IF

    ELSE

       ! order 0 must always be calculated
       g(0) = A(z)*LOG(r) + B(z) + C(z)*r*LOG(r) + D(z)*r
       IF (order >= 1) THEN
          g(1) = A(z)/r + C(z)*LOG(r) + C(z) + D(z)
       END IF
       IF (order >= 2) THEN
          rr = r*r
          g(2) = -A(z)/rr + C(z)/r
       END IF
       IF (order >= 3) THEN
          g(3) = 2.0_dbl*A(z)/(rr*r) - C(z)/rr
       END IF
       
    END IF

  END SUBROUTINE calc_g

! *****************************************************************************

  SUBROUTINE pz_lda_ed_loc ( rho, ed, order )
    
    
    REAL(DBL), INTENT(IN)                    :: rho
    REAL(DBL), DIMENSION(0:), INTENT(INOUT)    :: ed
    INTEGER, INTENT(IN)                      :: order

    INTEGER                                  :: m, order_
    LOGICAL, DIMENSION(0:3)                  :: calc
    REAL(DBL), DIMENSION(0:3)                :: e0, r

    calc = .FALSE.

    order_ = order
    IF (order_ >= 0) THEN
       calc(0:order_) = .TRUE.
    ELSE
       order_ = -1 * order_
       calc(order_) = .TRUE.
    END IF

    CALL calc_rs(rho, r(0))
    CALL calc_g(r(0), 0, e0, order_)

    IF (order_ >= 1) r(1) = (-1._dbl/3._dbl) * r(0) / rho
    IF (order_ >= 2) r(2) = (-4._dbl/3._dbl) * r(1) / rho
    IF (order_ >= 3) r(3) = (-7._dbl/3._dbl) * r(2) / rho

    m = 0
    IF (calc(0)) THEN
       ed(m) = e0(0)
       m = m + 1
    END IF
    IF (calc(1)) THEN
       ed(m) = e0(1)*r(1)
       m =  m + 1
    END IF
    IF (calc(2)) THEN
       ed(m) = e0(2)*r(1)**2 + e0(1)*r(2)
       m = m + 1
    END IF
    IF (calc(3)) THEN
       ed(m) = e0(3)*r(1)**3 + e0(2)*3._dbl*r(1)*r(2) + e0(1)*r(3)
    END IF

  END SUBROUTINE pz_lda_ed_loc

! *****************************************************************************

  SUBROUTINE pz_lsd_ed_loc ( a, b, ed, order )
    
    
    REAL(DBL), INTENT(IN)                    :: a, b
    REAL(DBL), DIMENSION(0:), INTENT(INOUT)    :: ed
    INTEGER, INTENT(IN), OPTIONAL            :: order

    INTEGER                                  :: m, order_
    LOGICAL, DIMENSION(0:3)                  :: calc
    REAL(DBL)                                :: rho, tr, trr, trrr, trrz, &
                                                trz, trzz, tz, tzz, tzzz
    REAL(DBL), DIMENSION(0:3)                :: e0, e1, f, r
    REAL(DBL), DIMENSION(0:3, 0:3)           :: z

    calc = .FALSE.

    order_ = 0
    IF (PRESENT(order)) order_ = order

    IF (order_ >= 0) THEN
       calc(0:order_) = .TRUE.
    ELSE
       order_ = -1 * order_
       calc(order_) = .TRUE.
    END IF

    rho = a + b

    CALL calc_fx(a, b, f(0:order_), order_)
    CALL calc_rs(rho, r(0))
    CALL calc_g(r(0), 0, e0(0:order_), order_)
    CALL calc_g(r(0), 1, e1(0:order_), order_)
    CALL calc_z(a, b, z(:,:), order_)

!! calculate first partial derivatives
    IF (order_ >= 1) THEN
       r(1) = (-1._dbl/3._dbl)*r(0)/rho
       tr = e0(1) + (e1(1)-e0(1))*f(0)
       tz = (e1(0)-e0(0))*f(1)
    END IF
    
!! calculate second partial derivatives
    IF (order_ >= 2) THEN
       r(2) = (-4._dbl/3._dbl)*r(1)/rho
       trr = e0(2) + (e1(2)-e0(2))*f(0)
       trz = (e1(1)-e0(1))*f(1)
       tzz = (e1(0)-e0(0))*f(2)
    END IF
    
!! calculate third derivatives
    IF (order_ >= 3) THEN
       r(3) = (-7._dbl/3._dbl)*r(2)/rho
       trrr = e0(3) + (e1(3)-e0(3))*f(0)
       trrz = (e1(2)-e0(2))*f(1)
       trzz = (e1(1)-e0(1))*f(2)
       tzzz = (e1(0)-e0(0))*f(3)
    END IF
    
    m=0
    IF (calc(0)) THEN
       ed(m) = e0(0) + (e1(0)-e0(0))*f(0)
       m = m + 1
    END IF
    
    IF (calc(1)) THEN
       ed(m)   = tr*r(1) + tz*z(1,0)
       ed(m+1) = tr*r(1) + tz*z(0,1)
       m = m + 2
    END IF
    
    IF (calc(2)) THEN
       ed(m)   = trr*r(1)**2 + 2._dbl*trz*r(1)*z(1,0) &
            + tr*r(2) + tzz*z(1,0)**2 + tz*z(2,0)
       ed(m+1) = trr*r(1)**2 + trz*r(1)*(z(0,1)+z(1,0)) &
            + tr*r(2) + tzz*z(1,0)*z(0,1) + tz*z(1,1)
       ed(m+2) = trr*r(1)**2 + 2._dbl*trz*r(1)*z(0,1) &
            + tr*r(2) + tzz*z(0,1)**2 + tz*z(0,2)
       m = m + 3
    END IF
       
    IF (calc(3)) THEN
       ed(m)   = trrr*r(1)**3 + 3._dbl*trrz*r(1)**2*z(1,0) &
            + 3._dbl*trr*r(1)*r(2) + 3._dbl*trz*r(2)*z(1,0) &
            + tr*r(3) + 3._dbl*trzz*r(1)*z(1,0)**2 &
            + tzzz*z(1,0)**3 + 3._dbl*trz*r(1)*z(2,0) &
            + 3._dbl*tzz*z(1,0)*z(2,0) + tz*z(3,0)
       ed(m+1) = trrr*r(1)**3 + trrz*r(1)**2*(2._dbl*z(1,0)+z(0,1)) &
            + 2._dbl*trzz*r(1)*z(1,0)*z(0,1) &
            + 2._dbl*trz*(r(2)*z(1,0)+r(1)*z(1,1)) &
            + 3._dbl*trr*r(2)*r(1) + trz*r(2)*z(0,1) + tr*r(3) &
            + trzz*r(1)*z(1,0)**2 + tzzz*z(1,0)**2*z(0,1) &
            + 2._dbl*tzz*z(1,0)*z(1,1) &
            + trz*r(1)*z(2,0) + tzz*z(2,0)*z(0,1) + tz*z(2,1)
       ed(m+2) = trrr*r(1)**3 + trrz*r(1)**2*(2._dbl*z(0,1)+z(1,0)) &
            + 2._dbl*trzz*r(1)*z(0,1)*z(1,0) &
            + 2._dbl*trz*(r(2)*z(0,1)+r(1)*z(1,1)) &
            + 3._dbl*trr*r(2)*r(1) + trz*r(2)*z(1,0) + tr*r(3) &
            + trzz*r(1)*z(0,1)**2 + tzzz*z(0,1)**2*z(1,0) &
            + 2._dbl*tzz*z(0,1)*z(1,1) &
            + trz*r(1)*z(0,2) + tzz*z(0,2)*z(1,0) + tz*z(1,2)
       ed(m+3) = trrr*r(1)**3 + 3._dbl*trrz*r(1)**2*z(0,1) &
            + 3._dbl*trr*r(1)*r(2) + 3._dbl*trz*r(2)*z(0,1) + tr*r(3) &
            + 3._dbl*trzz*r(1)*z(0,1)**2 + tzzz*z(0,1)**3 &
            + 3._dbl*trz*r(1)*z(0,2) &
            + 3._dbl*tzz*z(0,1)*z(0,2) + tz*z(0,3)
    END IF
    
  END SUBROUTINE pz_lsd_ed_loc

! *****************************************************************************

END MODULE XC_PERDEW_ZUNGER
