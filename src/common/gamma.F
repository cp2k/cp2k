!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of the incomplete Gamma function F_n(t) for multi-center
!>      integrals over Cartesian Gaussian functions.
!> \par History
!>      - restructured and cleaned (24.05.2004,MK)
!> \author Matthias Krack (07.01.1999)
! *****************************************************************************
MODULE gamma

  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE mathconstants,                   ONLY: ifac,&
                                             pi
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
#include "cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gamma'
  REAL(KIND = dp), PARAMETER  :: teps = 1.0E-13_dp

! *** Maximum n value of the tabulated F_n(t) function values ***

  INTEGER, SAVE :: current_nmax = -1

! *** F_n(t) table ***

  REAL(KIND = dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: ftable

! *** Public subroutines ***

  PUBLIC :: deallocate_md_ftable,&
            fgamma_0, fgamma_1,&
            fgamma_ref,&
            init_md_ftable

CONTAINS

! *****************************************************************************
!> \brief   Build a table of F_n(t) values in the range tmin <= t <= tmax
!>          with a stepsize of tdelta up to n equal to nmax.
!> \param nmax ...
!> \param tmin ...
!> \param tmax ...
!> \param tdelta ...
!> \date    11.01.1999
!> \par Parameters
!>       - nmax  : Maximum n value of F_n(t).
!>       - tdelta: Difference between two consecutive t abcissas (step size).
!>       - tmax  : Maximum t value.
!>       - tmin  : Minimum t value.
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE create_md_ftable(nmax,tmin,tmax,tdelta)

    INTEGER, INTENT(IN)                      :: nmax
    REAL(KIND=dp), INTENT(IN)                :: tmin, tmax, tdelta

    CHARACTER(len=*), PARAMETER :: routineN = 'create_md_ftable', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat, itab, itabmax, &
                                                itabmin, n
    REAL(KIND=dp)                            :: t

    IF (current_nmax > -1) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "An incomplete Gamma function table is already "//&
                        " allocated. Use the init routine for an update")
    END IF

    IF (nmax < 0) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "A negative n value for the initialization of the "//&
                        "incomplete Gamma function is invalid")
    END IF

!   *** Check arguments ***

    IF ((tmax < 0.0_dp).OR.&
        (tmin < 0.0_dp).OR.&
        (tdelta <= 0.0_dp).OR.&
        (tmin > tmax)) THEN
      CALL stop_program(routineN,moduleN,__LINE__,"Invalid arguments")
    END IF

    n = nmax + 6

    itabmin = FLOOR(tmin/tdelta)
    itabmax = CEILING((tmax - tmin)/tdelta)

    ALLOCATE (ftable(0:n,itabmin:itabmax),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineN,moduleN,__LINE__,"ftable",&
                       (n + 1)*(itabmax - itabmin + 1)*dp_size)
    END IF
    ftable = 0.0_dp

!   *** Fill table ***

    DO itab=itabmin,itabmax
      t = REAL(itab,dp)*tdelta
      ftable(0:n,itab) = fgamma_ref(n,t)
    END DO

!   *** Save initialization status ***

    current_nmax = nmax

  END SUBROUTINE create_md_ftable

! *****************************************************************************
!> \brief   Deallocate the table of F_n(t) values.
!> \date    24.05.2004
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_md_ftable()

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_md_ftable', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    IF (current_nmax > -1) THEN

      DEALLOCATE (ftable,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"ftable")
      END IF

      current_nmax = -1

    END IF

  END SUBROUTINE deallocate_md_ftable

! *****************************************************************************
!> \brief   Calculation of the incomplete Gamma function F(t) for multicenter
!>          integrals over Gaussian functions. f returns a vector with all
!>          F_n(t) values for 0 <= n <= nmax.
!> \param nmax ...
!> \param t ...
!> \param f ...
!> \date    08.01.1999,
!> \par History
!>          09.06.1999, MK : Changed from a FUNCTION to a SUBROUTINE
!> \par Literature
!>       L. E. McMurchie, E. R. Davidson, J. Comp. Phys. 26, 218 (1978)
!> \par Parameters
!>       - f   : The incomplete Gamma function F_n(t).
!>       - nmax: Maximum n value of F_n(t).
!>       - t   : Argument of the incomplete Gamma function.
!>       - kmax: Maximum number of iterations.
!>       - expt: Exponential term in the upward recursion of F_n(t).
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE fgamma_0(nmax,t,f)

    INTEGER, INTENT(IN)                      :: nmax
    REAL(KIND=dp), INTENT(IN)                :: t
    REAL(KIND=dp), DIMENSION(0:nmax), &
      INTENT(OUT)                            :: f

    INTEGER                                  :: itab, k, n
    REAL(KIND=dp)                            :: expt, g, tdelta, tmp, ttab

!   *** Calculate F(t) ***

    IF (t < teps) THEN

!     *** Special cases: t = 0 ***

      DO n=0,nmax
        f(n) = 1.0_dp/REAL(2*n + 1,dp)
      END DO

    ELSE IF (t <= 12.0_dp) THEN

!     *** 0 < t < 12 -> Taylor expansion ***

      tdelta = 0.1_dp

!     *** Pretabulation of the F_n(t) function ***
!     *** for the Taylor series expansion      ***

      IF (nmax > current_nmax) THEN
        CALL init_md_ftable(nmax)
      END IF


      itab = NINT(t/tdelta)
      ttab = REAL(itab,dp)*tdelta

      f(nmax) = ftable(nmax,itab)

      tmp = 1.0_dp
      DO k=1,6
        tmp = tmp * (ttab - t)
        f(nmax) = f(nmax) + ftable(nmax+k,itab)*tmp*ifac(k)
      END DO

      expt = EXP(-t)

!     *** Use the downward recursion relation to ***
!     *** generate the remaining F_n(t) values   ***

      DO n=nmax-1,0,-1
        f(n) = (2.0_dp*t*f(n+1) + expt)/REAL(2*n + 1,dp)
      END DO

    ELSE

!     *** t > 12 ***

      IF (t <= 15.0_dp) THEN

!       *** 12 < t <= 15 -> Four term polynom expansion ***

        g = 0.4999489092_dp - 0.2473631686_dp/t +&
            0.321180909_dp/t**2 - 0.3811559346_dp/t**3
        f(0) = 0.5_dp*SQRT(pi/t) - g*EXP(-t)/t

      ELSE IF (t <= 18.0_dp) THEN

!       *** 15 < t <= 18 -> Three term polynom expansion ***

        g = 0.4998436875_dp - 0.24249438_dp/t + 0.24642845_dp/t**2
        f(0) = 0.5_dp*SQRT(pi/t) - g*EXP(-t)/t

      ELSE IF (t <= 24.0_dp) THEN

!       *** 18 < t <= 24 -> Two term polynom expansion ***

        g = 0.499093162_dp - 0.2152832_dp/t
        f(0) = 0.5_dp*SQRT(pi/t) - g*EXP(-t)/t

      ELSE IF (t <= 30.0_dp) THEN

!       *** 24 < t <= 30 -> One term polynom expansion ***

        g = 0.49_dp
        f(0) = 0.5_dp*SQRT(pi/t) - g*EXP(-t)/t

      ELSE

!       *** t > 30 -> Asymptotic formula ***

        f(0) = 0.5_dp*SQRT(pi/t)

      END IF

      IF (t > REAL(2*nmax + 36,dp)) THEN
        expt = 0.0_dp
      ELSE
        expt = EXP(-t)
      END IF

!     *** Use the upward recursion relation to ***
!     *** generate the remaining F_n(t) values ***

      DO n=1,nmax
        f(n) = 0.5_dp*(REAL(2*n - 1,dp)*f(n-1) - expt)/t
      END DO

    END IF

  END SUBROUTINE fgamma_0

! *****************************************************************************
!> \brief   Calculation of the incomplete Gamma function F(t) for multicenter
!>          integrals over Gaussian functions. f returns a vector with all
!>          F_n(t) values for 0 <= n <= nmax.
!> \param nmax ...
!> \param t ...
!> \param f ...
!> \date    08.01.1999
!> \par Literature
!>       L. E. McMurchie, E. R. Davidson, J. Comp. Phys. 26, 218 (1978)
!> \par Parameters
!>       - f   : The incomplete Gamma function F_n(t).
!>       - nmax: Maximum n value of F_n(t).
!>       - t   : Argument of the incomplete Gamma function.
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE fgamma_1(nmax,t,f)

    INTEGER, INTENT(IN)                      :: nmax
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: t
    REAL(KIND=dp), &
      DIMENSION(SIZE(t, 1), 0:nmax), &
      INTENT(OUT)                            :: f

    INTEGER                                  :: i, itab, k, n
    REAL(KIND=dp)                            :: expt, g, tdelta, tmp, ttab

    DO i=1,SIZE(t,1)

!     *** Calculate F(t) ***

      IF (t(i) < teps) THEN

!       *** Special cases: t = 0 ***

        DO n=0,nmax
          f(i,n) = 1.0_dp/REAL(2*n + 1,dp)
        END DO

      ELSE IF (t(i) <= 12.0_dp) THEN

!       *** 0 < t < 12 -> Taylor expansion ***

        tdelta = 0.1_dp

!       *** Pretabulation of the F_n(t) function ***
!       *** for the Taylor series expansion      ***

        IF (nmax > current_nmax) THEN
         CALL init_md_ftable(nmax)
        END IF

        itab = NINT(t(i)/tdelta)
        ttab = REAL(itab,dp)*tdelta

        f(i,nmax) = ftable(nmax,itab)

        tmp = 1.0_dp
        DO k=1,6
          tmp = tmp * (ttab - t(i))
          f(i,nmax) = f(i,nmax) + ftable(nmax+k,itab)*tmp*ifac(k)
        END DO

        expt = EXP(-t(i))

!       *** Use the downward recursion relation to ***
!       *** generate the remaining F_n(t) values   ***

        DO n=nmax-1,0,-1
          f(i,n) = (2.0_dp*t(i)*f(i,n+1) + expt)/REAL(2*n + 1,dp)
        END DO

      ELSE

!       *** t > 12 ***

        IF (t(i) <= 15.0_dp) THEN

!         *** 12 < t <= 15 -> Four term polynom expansion ***

          g = 0.4999489092_dp - 0.2473631686_dp/t(i) +&
              0.321180909_dp/t(i)**2 - 0.3811559346_dp/t(i)**3
          f(i,0) = 0.5_dp*SQRT(pi/t(i)) - g*EXP(-t(i))/t(i)

        ELSE IF (t(i) <= 18.0_dp) THEN

!         *** 15 < t <= 18 -> Three term polynom expansion ***

          g = 0.4998436875_dp - 0.24249438_dp/t(i) + 0.24642845_dp/t(i)**2
          f(i,0) = 0.5_dp*SQRT(pi/t(i)) - g*EXP(-t(i))/t(i)

        ELSE IF (t(i) <= 24.0_dp) THEN

!         *** 18 < t <= 24 -> Two term polynom expansion ***

          g = 0.499093162_dp - 0.2152832_dp/t(i)
          f(i,0) = 0.5_dp*SQRT(pi/t(i)) - g*EXP(-t(i))/t(i)

        ELSE IF (t(i) <= 30.0_dp) THEN

!         *** 24 < t <= 30 -> One term polynom expansion ***

          g = 0.49_dp
          f(i,0) = 0.5_dp*SQRT(pi/t(i)) - g*EXP(-t(i))/t(i)

        ELSE

!         *** t > 30 -> Asymptotic formula ***

          f(i,0) = 0.5_dp*SQRT(pi/t(i))

        END IF

        IF (t(i) > REAL(2*nmax + 36,dp)) THEN
          expt = 0.0_dp
        ELSE
          expt = EXP(-t(i))
        END IF

!       *** Use the upward recursion relation to ***
!       *** generate the remaining F_n(t) values ***

        DO n=1,nmax
          f(i,n) = 0.5_dp*(REAL(2*n - 1,dp)*f(i,n-1) - expt)/t(i)
        END DO

      END IF

    END DO

  END SUBROUTINE fgamma_1

! *****************************************************************************
!> \brief   Calculation of the incomplete Gamma function F_n(t) using a
!>          spherical Bessel function expansion. fgamma_ref returns a
!>          vector with all F_n(t) values for 0 <= n <= nmax.
!>          For t values greater than 50 an asymptotic formula is used.
!>          This function is expected to return accurate F_n(t) values
!>          for any combination of n and t, but the calculation is slow
!>          and therefore the function may only be used for a pretabulation
!>          of F_n(t) values or for reference calculations.
!> \param nmax ...
!> \param t ...
!> \retval f ...
!> \date    07.01.1999
!> \par Literature
!>        F. E. Harris, Int. J. Quant. Chem. 23, 1469 (1983)
!> \par Parameters
!>       - expt   : Exponential term in the downward recursion of F_n(t).
!>       - factor : Prefactor of the Bessel function expansion.
!>       - nmax   : Maximum n value of F_n(t).
!>       - p      : Product of the Bessel function quotients.
!>       - r      : Quotients of the Bessel functions.
!>       - sumterm: One term in the sum over products of Bessel functions.
!>       - t      : Argument of the incomplete Gamma function.
!> \author  MK
!> \version 1.0
! *****************************************************************************
  FUNCTION fgamma_ref(nmax,t) RESULT(f)

    INTEGER, INTENT(IN)                      :: nmax
    REAL(KIND=dp), INTENT(IN)                :: t
    REAL(KIND=dp), DIMENSION(0:nmax)         :: f

    CHARACTER(len=*), PARAMETER :: routineN = 'fgamma_ref', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: kmax = 50
    REAL(KIND=dp), PARAMETER                 :: eps = EPSILON(0.0_dp)

    INTEGER                                  :: j, k, n
    REAL(KIND=dp)                            :: expt, factor, p, sumterm, &
                                                sumtot, term
    REAL(KIND=dp), DIMENSION(kmax+10)        :: r

!   ------------------------------------------------------------------
!   *** Initialization ***

    f(:) = 0.0_dp

    IF (t < teps) THEN

!     *** Special case: t = 0 => analytic expression ***

      DO n=0,nmax
        f(n) = 1.0_dp/REAL(2*n + 1,dp)
      END DO

    ELSE IF (t <= 50.0_dp) THEN

!     *** Initialize ratios of Bessel functions ***

      r(kmax+10) = 0.0_dp

      DO j=kmax+9,1,-1
        r(j) = -t/(REAL(4*j + 2,dp) - t*r(j+1))
      END DO

      factor = 2.0_dp*SINH(0.5_dp*t)*EXP(-0.5_dp*t)/t

      DO n=0,nmax

!       *** Initialize iteration ***

        sumtot = factor/REAL(2*n + 1,dp)
        term = 1.0_dp

!       *** Begin the summation and recursion ***

        DO k=1,kmax

          term = term*REAL(2*n - 2*k + 1,dp)/REAL(2*n + 2*k + 1,dp)

!         *** Product of Bessel function quotients ***

          p = 1.0_dp

          DO j=1,k
            p = p*r(j)
          END DO

          sumterm = factor*term*p*REAL(2*k + 1,dp)/REAL(2*n + 1,dp)

          IF (ABS(sumterm) < eps) THEN

!           *** Iteration converged ***

            EXIT

          ELSE IF (k == kmax) THEN

!           *** No convergence with kmax iterations ***

            CALL stop_program(routineN,moduleN,__LINE__,&
                              "Maximum number of iterations reached")

          ELSE

!           *** Add the current term to the sum and continue the iteration ***

            sumtot = sumtot + sumterm

          END IF

        END DO

        f(n) = sumtot

      END DO

    ELSE

!     *** Use asymptotic formula for t > 50 ***

      f(0) = 0.5_dp*SQRT(pi/t)

!     *** Use the upward recursion relation to ***
!     *** generate the remaining F_n(t) values ***

      expt = EXP(-t)

      DO n=1,nmax
        f(n) = 0.5_dp*(REAL(2*n - 1,dp)*f(n-1) - expt)/t
      END DO

    END IF

  END FUNCTION fgamma_ref

! *****************************************************************************
!> \brief   Initalize a table of F_n(t) values in the range 0 <= t <= 12 with
!>            a stepsize of 0.1 up to n equal to nmax for the Taylor series
!>            expansion used by McMurchie-Davidson (MD).
!> \param nmax ...
!> \date    10.06.1999
!> \par Parameters
!>       - nmax   : Maximum n value of F_n(t).
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE init_md_ftable(nmax)
    INTEGER, INTENT(IN)                      :: nmax

    CHARACTER(len=*), PARAMETER :: routineN = 'init_md_ftable', &
      routineP = moduleN//':'//routineN

    IF (nmax < 0) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "A negative n value for the initialization of the "//&
                        "incomplete Gamma function is invalid")
    END IF

!   *** Check, if the current initialization is sufficient ***

    IF (nmax > current_nmax) THEN

      CALL deallocate_md_ftable()

!     *** Pretabulation of the F_n(t) function ***
!     *** for the Taylor series expansion      ***

      CALL create_md_ftable(nmax,0.0_dp,12.0_dp,0.1_dp)

    END IF

  END SUBROUTINE init_md_ftable

END MODULE gamma
