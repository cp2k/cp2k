!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/gamma [1.0] *
!!
!!   NAME
!!     gamma
!!
!!   FUNCTION
!!     Calculation of the incomplete Gamma function F_n(t) for multicenter
!!     integrals over Cartesian Gaussian functions.
!!
!!   AUTHOR
!!     Matthias Krack (07.01.1999)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE gamma

! *****************************************************************************

! Index:

! SUBROUTINE build_ftable(nmax,tmin,tmax,tdelta)
! SUBROUTINE fgamma_0(nmax,t,f)
! SUBROUTINE fgamma_1(nmax,t,f)
! SUBROUTINE init_md_ftable(nmax)

! FUNCTION fgamma_ref(nmax,t) RESULT(f)

! *****************************************************************************

! Global variables:

! init_md_ftable_done: If .TRUE., the pretabulation for the F_n(t) calculation
!                      after McMurchie-Davidson is done.
! nmax_ftable        : Maximum n value of the tabulated F_n(t) function values.
! teps               : If t is lower than teps, then t is set to zero.

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  USE mathconstants, ONLY: fac,pi
  USE termination,   ONLY: stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameter ***

  REAL(wp), PARAMETER :: teps = 1.0E-13_wp

! *** Global variables ***

  INTEGER :: nmax_ftable = 0
  LOGICAL :: init_md_ftable_done = .FALSE.

  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ftable

! *** Public subroutines ***

  PUBLIC :: fgamma,&
            init_md_ftable

! *****************************************************************************

  INTERFACE fgamma
    MODULE PROCEDURE fgamma_0,fgamma_1
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE build_ftable(nmax,tmin,tmax,tdelta)

!   Purpose: Build a table of F_n(t) values in the range tmin <= t <= tmax
!            with a stepsize of tdelta up to n equal to nmax.

!   History: - Creation (11.01.1999, Matthias Krack)

!   ***************************************************************************

!   List of local variables:

!   nmax  : Maximum n value of F_n(t).
!   tdelta: Difference between two consecutive t abcissas (step size).
!   tmax  : Maximum t value.
!   tmin  : Minimum t value.

!   ***************************************************************************

    REAL(wp), INTENT(IN) :: tdelta,tmax,tmin
    INTEGER, INTENT(IN)  :: nmax

!   *** Local variables ***

    REAL(wp) :: t
    INTEGER  :: itab,itabmax,itabmin

!   ---------------------------------------------------------------------------

!   *** Check arguments ***

    IF ((tmax < 0.0_wp).OR.&
        (tmin < 0.0_wp).OR.&
        (tdelta <= 0.0_wp).OR.&
        (tmin > tmax).OR.&
        (nmax < 0)) THEN
      CALL stop_program("SUBROUTINE build_ftable (MODULE gamma)",&
                        "Invalid arguments")
    END IF

!   *** Memory allocation ***

    IF (ALLOCATED(ftable)) DEALLOCATE (ftable)

    itabmin = FLOOR(tmin/tdelta)
    itabmax = CEILING((tmax - tmin)/tdelta)

    ALLOCATE (ftable(0:nmax,itabmin:itabmax))

!   *** Initialization ***

    ftable(:,:) = 0.0_wp

!   *** Fill table ***

    DO itab=itabmin,itabmax
      t = REAL(itab,wp)*tdelta
      ftable(0:nmax,itab) = fgamma_ref(nmax,t)
    END DO

    nmax_ftable = nmax

  END SUBROUTINE build_ftable

! *****************************************************************************

  SUBROUTINE fgamma_0(nmax,t,f)

!   Purpose: Calculation of the incomplete Gamma function F(t) for multicenter
!            integrals over Gaussian functions. f returns a vector with all
!            F_n(t) values for 0 <= n <= nmax.

!   Literature: L. E. McMurchie, E. R. Davidson, J. Comp. Phys. 26, 218 (1978)

!   History: - Creation (08.01.1999, Matthias Krack)
!            - Changed from a FUNCTION to a SUBROUTINE
!              (09.06.1999, Matthias Krack)

!   ***************************************************************************

!   f   : The incomplete Gamma function F_n(t).
!   nmax: Maximum n value of F_n(t).
!   t   : Argument of the incomplete Gamma function.

!   List of local parameters:

!   kmax: Maximum number of iterations.

!   List of local variables:

!   expt: Exponential term in the upward recursion of F_n(t).

!   ***************************************************************************

    REAL(wp), INTENT(IN)                     :: t
    INTEGER, INTENT(IN)                      :: nmax
    REAL(wp), DIMENSION(0:nmax), INTENT(OUT) :: f

!   *** Local variables ***

    REAL(wp) :: expt,g,tdelta,ttab
    INTEGER  :: itab,k,n

!   ---------------------------------------------------------------------------

!   *** Calculate F(t) ***

    IF (t < teps) THEN

!     *** Special cases: t = 0 ***

      DO n=0,nmax
        f(n) = 1.0_wp/REAL(2*n + 1,wp)
      END DO

    ELSE IF (t <= 12.0_wp) THEN

!     *** 0 < t < 12 -> Taylor expansion ***

      tdelta = 0.1_wp

!     *** Pretabulation of the F_n(t) function ***
!     *** for the Taylor series expansion      ***

      IF ((.NOT.init_md_ftable_done).OR.(nmax + 6 > nmax_ftable)) THEN
        CALL init_md_ftable(nmax)
      END IF

      itab = NINT(t/tdelta)
      ttab = REAL(itab,wp)*tdelta

      f(nmax) = ftable(nmax,itab)

      DO k=1,6
        f(nmax) = f(nmax) + ftable(nmax+k,itab)*(ttab - t)**k/fac(k)
      END DO

      expt = EXP(-t)

!     *** Use the downward recursion relation to ***
!     *** generate the remaining F_n(t) values   ***

      DO n=nmax-1,0,-1
        f(n) = (2.0_wp*t*f(n+1) + expt)/REAL(2*n + 1,wp)
      END DO

    ELSE

!     *** t > 12 ***

      IF (t <= 15.0_wp) THEN

!       *** 12 < t <= 15 -> Four term polynom expansion ***

        g = 0.4999489092_wp - 0.2473631686_wp/t +&
            0.321180909_wp/t**2 - 0.3811559346_wp/t**3
        f(0) = 0.5_wp*SQRT(pi/t) - g*EXP(-t)/t

      ELSE IF (t <= 18.0_wp) THEN

!       *** 15 < t <= 18 -> Three term polynom expansion ***

        g = 0.4998436875_wp - 0.24249438_wp/t + 0.24642845_wp/t**2
        f(0) = 0.5_wp*SQRT(pi/t) - g*EXP(-t)/t

      ELSE IF (t <= 24.0_wp) THEN

!       *** 18 < t <= 24 -> Two term polynom expansion ***

        g = 0.499093162_wp - 0.2152832_wp/t
        f(0) = 0.5_wp*SQRT(pi/t) - g*EXP(-t)/t

      ELSE IF (t <= 30.0_wp) THEN

!       *** 24 < t <= 30 -> One term polynom expansion ***

        g = 0.49_wp
        f(0) = 0.5_wp*SQRT(pi/t) - g*EXP(-t)/t

      ELSE

!       *** t > 30 -> Asymptotic formula ***

        f(0) = 0.5_wp*SQRT(pi/t)

      END IF

      IF (t > REAL(2*nmax + 36,wp)) THEN
        expt = 0.0_wp
      ELSE
        expt = EXP(-t)
      END IF

!     *** Use the upward recursion relation to ***
!     *** generate the remaining F_n(t) values ***

      DO n=1,nmax
        f(n) = 0.5_wp*(REAL(2*n - 1,wp)*f(n-1) - expt)/t
      END DO

    END IF

  END SUBROUTINE fgamma_0

! *****************************************************************************

  SUBROUTINE fgamma_1(nmax,t,f)

!   Purpose: Calculation of the incomplete Gamma function F(t) for multicenter
!            integrals over Gaussian functions. f returns a vector with all
!            F_n(t) values for 0 <= n <= nmax.

!   Literature: L. E. McMurchie, E. R. Davidson, J. Comp. Phys. 26, 218 (1978)

!   History: - Creation (08.01.1999, Matthias Krack)

!   ***************************************************************************

!   f   : The incomplete Gamma function F_n(t).
!   nmax: Maximum n value of F_n(t).
!   t   : Argument of the incomplete Gamma function.

!   ***************************************************************************

    INTEGER, INTENT(IN)                                :: nmax
    REAL(wp), DIMENSION(:), INTENT(IN)                 :: t
    REAL(wp), DIMENSION(SIZE(t,1),0:nmax), INTENT(OUT) :: f

!   *** Local variables ***

    REAL(wp) :: expt,g,tdelta,ttab
    INTEGER  :: i,itab,k,n

!   ---------------------------------------------------------------------------

    DO i=1,SIZE(t,1)

!     *** Calculate F(t) ***

      IF (t(i) < teps) THEN

!       *** Special cases: t = 0 ***

        DO n=0,nmax
          f(i,n) = 1.0_wp/REAL(2*n + 1,wp)
        END DO

      ELSE IF (t(i) <= 12.0_wp) THEN

!       *** 0 < t < 12 -> Taylor expansion ***

        tdelta = 0.1_wp

!       *** Pretabulation of the F_n(t) function ***
!       *** for the Taylor series expansion      ***

        IF ((.NOT.init_md_ftable_done).OR.(nmax + 6 > nmax_ftable)) THEN
          CALL init_md_ftable(nmax)
        END IF

        itab = NINT(t(i)/tdelta)
        ttab = REAL(itab,wp)*tdelta

        f(i,nmax) = ftable(nmax,itab)

        DO k=1,6
          f(i,nmax) = f(i,nmax) + ftable(nmax+k,itab)*(ttab - t(i))**k/fac(k)
        END DO

        expt = EXP(-t(i))

!       *** Use the downward recursion relation to ***
!       *** generate the remaining F_n(t) values   ***

        DO n=nmax-1,0,-1
          f(i,n) = (2.0_wp*t(i)*f(i,n+1) + expt)/REAL(2*n + 1,wp)
        END DO

      ELSE

!       *** t > 12 ***

        IF (t(i) <= 15.0_wp) THEN

!         *** 12 < t <= 15 -> Four term polynom expansion ***

          g = 0.4999489092_wp - 0.2473631686_wp/t(i) +&
              0.321180909_wp/t(i)**2 - 0.3811559346_wp/t(i)**3
          f(i,0) = 0.5_wp*SQRT(pi/t(i)) - g*EXP(-t(i))/t(i)

        ELSE IF (t(i) <= 18.0_wp) THEN

!         *** 15 < t <= 18 -> Three term polynom expansion ***

          g = 0.4998436875_wp - 0.24249438_wp/t(i) + 0.24642845_wp/t(i)**2
          f(i,0) = 0.5_wp*SQRT(pi/t(i)) - g*EXP(-t(i))/t(i)

        ELSE IF (t(i) <= 24.0_wp) THEN

!         *** 18 < t <= 24 -> Two term polynom expansion ***

          g = 0.499093162_wp - 0.2152832_wp/t(i)
          f(i,0) = 0.5_wp*SQRT(pi/t(i)) - g*EXP(-t(i))/t(i)

        ELSE IF (t(i) <= 30.0_wp) THEN

!         *** 24 < t <= 30 -> One term polynom expansion ***

          g = 0.49_wp
          f(i,0) = 0.5_wp*SQRT(pi/t(i)) - g*EXP(-t(i))/t(i)

        ELSE

!         *** t > 30 -> Asymptotic formula ***

          f(i,0) = 0.5_wp*SQRT(pi/t(i))

        END IF

        IF (t(i) > REAL(2*nmax + 36,wp)) THEN
          expt = 0.0_wp
        ELSE
          expt = EXP(-t(i))
        END IF

!       *** Use the upward recursion relation to ***
!       *** generate the remaining F_n(t) values ***

        DO n=1,nmax
          f(i,n) = 0.5_wp*(REAL(2*n - 1,wp)*f(i,n-1) - expt)/t(i)
        END DO

      END IF

    END DO

  END SUBROUTINE fgamma_1

! *****************************************************************************

  FUNCTION fgamma_ref(nmax,t) RESULT(f)

!   Purpose: Calculation of the incomplete Gamma function F_n(t) using a
!            spherical Bessel function expansion. fgamma_ref returns a
!            vector with all F_n(t) values for 0 <= n <= nmax.
!            For t values greater than 50 an asymptotic formula is used.
!            This function is expected to return accurate F_n(t) values
!            for any combination of n and t, but the calculation is slow
!            and therefore the function may only be used for a pretabulation
!            of F_n(t) values or for reference calculations.

!   Literature: F. E. Harris, Int. J. Quant. Chem. 23, 1469 (1983)

!   History: - Creation (07.01.1999, Matthias Krack)

!   ***************************************************************************

!   List of local parameters:

!   eps : Numerical convergence criterion.
!   kmax: Maximum number of iterations.

!   List of local variables:

!   expt   : Exponential term in the downward recursion of F_n(t).
!   factor : Prefactor of the Bessel function expansion.
!   nmax   : Maximum n value of F_n(t).
!   p      : Product of the Bessel function quotients.
!   r      : Quotients of the Bessel functions.
!   sumterm: One term in the sum over products of Bessel functions.
!   t      : Argument of the incomplete Gamma function.

!   ***************************************************************************

    REAL(wp), INTENT(IN)        :: t
    INTEGER, INTENT(IN)         :: nmax
    REAL(wp), DIMENSION(0:nmax) :: f

!   *** Local parameters ***

    REAL(wp), PARAMETER :: eps = EPSILON(0.0_wp)
    INTEGER, PARAMETER  :: kmax = 50

!   *** Local variables ***

    REAL(wp) :: expt,factor,p,sumterm,sumtot,term
    INTEGER  :: j,k,n

    REAL(wp), DIMENSION(kmax+10) :: r

!   ------------------------------------------------------------------

!   *** Initialization ***

    f(:) = 0.0_wp

    IF (t < teps) THEN

!     *** Special case: t = 0 => analytic expression ***

      DO n=0,nmax
        f(n) = 1.0_wp/REAL(2*n + 1,wp)
      END DO

    ELSE IF (t <= 50.0_wp) THEN

!     *** Initialize ratios of Bessel functions ***

      r(kmax+10) = 0.0_wp

      DO j=kmax+9,1,-1
        r(j) = -t/(REAL(4*j + 2,wp) - t*r(j+1))
      END DO

      factor = 2.0_wp*SINH(0.5_wp*t)*EXP(-0.5_wp*t)/t

      DO n=0,nmax

!     *** Initialize iteration ***

        sumtot = factor/REAL(2*n + 1,wp)
        term = 1.0_wp

!       *** Begin the summation and recursion ***

        DO k=1,kmax

          term = term*REAL(2*n - 2*k + 1,wp)/REAL(2*n + 2*k + 1,wp)

!         *** Product of Bessel function quotients ***

          p = 1.0_wp

          DO j=1,k
            p = p*r(j)
          END DO

          sumterm = factor*term*p*REAL(2*k + 1,wp)/REAL(2*n + 1,wp)

          IF (ABS(sumterm) < eps) THEN

!           *** Iteration converged ***

            EXIT

          ELSE IF (k == kmax) THEN

!           *** No convergence with kmax iterations ***

            STOP "*** ERROR in FUNCTION fgamma_ref in MODULE gamma ***"

          ELSE

!           *** Add the current term to the sum and continue the iteration ***

            sumtot = sumtot + sumterm

          END IF

        END DO

        f(n) = sumtot

      END DO

    ELSE

!     *** Use asymptotic formula for t > 50 ***

      f(0) = 0.5_wp*SQRT(pi/t)

!     *** Use the upward recursion relation to ***
!     *** generate the remaining F_n(t) values ***

      expt = EXP(-t)

      DO n=1,nmax
        f(n) = 0.5_wp*(REAL(2*n - 1,wp)*f(n-1) - expt)/t
      END DO

    END IF

  END FUNCTION fgamma_ref

! *****************************************************************************

  SUBROUTINE init_md_ftable(nmax)

!   Purpose: Build a table of F_n(t) values in the range 0 <= t <= 12 with a
!            stepsize of 0.1 up to n equal to nmax for the Taylor series
!            expansion used by McMurchie-Davidson (MD).

!   History: - Creation (10.06.1999, Matthias Krack)

!   ***************************************************************************

!   nmax : Maximum n value of F_n(t).

!   ***************************************************************************

    INTEGER, INTENT(IN) :: nmax

!   ---------------------------------------------------------------------------

!   *** Pretabulation of the F_n(t) function ***
!   *** for the Taylor series expansion      ***

    CALL build_ftable(nmax+6,0.0_wp,12.0_wp,0.1_wp)

!   *** Table initialization done ***

    init_md_ftable_done = .TRUE.

  END SUBROUTINE init_md_ftable

! *****************************************************************************

END MODULE gamma
