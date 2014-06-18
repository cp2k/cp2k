!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Methods to perform on the fly statistical analysis of data
!>      -) Schiferl and Wallace, J. Chem. Phys. 83 (10) 1985
!> \author Teodoro Laino (01.2007) [tlaino]
!> \par History
!>         - Teodoro Laino (10.2008) [tlaino] - University of Zurich
!>           module made publicly available
! *****************************************************************************
MODULE statistical_methods
  USE erf_fn,                          ONLY: erf
  USE f77_blas
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE parallel_rng_types,              ONLY: next_random_number
  USE util,                            ONLY: sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'statistical_methods'
  LOGICAL, PARAMETER                   :: debug_this_module = .FALSE.
  INTEGER, PARAMETER, PUBLIC           :: min_sample_size = 20
  PUBLIC :: sw_test,&
            k_test,&
            vn_test

CONTAINS

! *****************************************************************************
!> \brief Shapiro - Wilk's test or W-statistic to test normality of a distribution
!>      R94  APPL. STATIST. (1995) VOL.44, NO.4
!>      Calculates the Shapiro-Wilk W test and its significance level
!> \par History
!>      Teodoro Laino (02.2007) [tlaino]
! *****************************************************************************
  SUBROUTINE sw_test ( ix, n, w, pw, error)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ix
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), INTENT(OUT)               :: w, pw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sw_test', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER :: bf1 = 0.8378_dp, c1(6) = (/ 0.000000_dp, &
      0.221157_dp,-0.147981_dp,-2.071190_dp,4.434685_dp, -2.706056_dp /), c2(6&
      ) = (/ 0.000000_dp,0.042981_dp,-0.293762_dp,-1.752461_dp,5.682633_dp, &
      -3.582633_dp /), &
      c3(4) = (/ 0.544000_dp,-0.399780_dp,0.025054_dp,-0.6714E-3_dp /), &
      c4(4) = (/ 1.3822_dp, -0.77857_dp, 0.062767_dp, -0.0020322_dp /), &
      c5(4) = (/ -1.5861_dp, -0.31082_dp, -0.083751_dp, 0.0038915_dp /), &
      c6(3) = (/ -0.4803_dp, -0.082676_dp, 0.0030302_dp /), &
      c7(2) = (/ 0.164_dp, 0.533_dp /), c8(2) = (/ 0.1736_dp, 0.315_dp /), &
      c9(2) = (/ 0.256_dp, -0.00635_dp /)
    REAL(KIND=dp), PARAMETER :: g(2) = (/ -2.273_dp, 0.459_dp /), &
      one = 1.0_dp, pi6 = 1.909859_dp, qtr = 0.25_dp, &
      small = EPSILON(0.0_dp), sqrth = 0.70711_dp, stqr = 1.047198_dp, &
      th = 0.375_dp, three = 3.0_dp, two = 2.0_dp, xx90 = 0.556_dp, &
      xx95 = 0.622_dp, z90 = 1.2816_dp, z95 = 1.6449_dp, z99 = 2.3263_dp, &
      zero = 0.0_dp, zm = 1.7509_dp, zss = 0.56268_dp

    INTEGER                                  :: i, i1, j, n2, output_unit, &
                                                stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: itmp
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a1, a2, an, an25, asa, fac, &
                                                gamma, m, range, rsn, s, sa, &
                                                sax, ssa, ssassx, ssumm2, &
                                                ssx, summ2, sx, w1, xi, xsx, &
                                                xx, y
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: a, x
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)
    ! Check for  N < 3
    IF (n  < 3 .OR. n > 5000)   THEN
       failure = .TRUE.
       IF (output_unit>0) WRITE(output_unit,'(A)')&
            "Shapiro Wilk test: number of points less than 3 or greated than 5000."
       IF (output_unit>0) WRITE(output_unit,'(A)')&
            "Not able to perform the test!"
    END IF
    ! Sort input array of data in ascending order
    IF (MOD(n,2)==0) THEN
       n2 = n/2
    ELSE
       n2 = (n-1)/2
    END IF
    ALLOCATE(x(n),stat=stat)
    CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    ALLOCATE(itmp(n),stat=stat)
    CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    ALLOCATE(a(n2),stat=stat)
    CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    x(:) = ix
    CALL sort(x, n, itmp)
    ! Check for zero range
    IF (.NOT.failure ) THEN
       range = x(n) - x(1)
       IF (range < small) failure = .TRUE.
       IF (failure.AND.(output_unit>0)) THEN
         WRITE(output_unit,'(A)')"Shapiro Wilk test: two data points are numerically identical."
         WRITE(output_unit,'(A)')"Not able to perform the test!"
       END IF
    END IF
    IF (.NOT.failure) THEN
       pw = one
       w  = one
       an = n
       ! Calculates coefficients for the test
       IF (n == 3) THEN
          a(1) = sqrth
       ELSE
          an25 = an + qtr
          summ2 = zero
          DO i = 1, n2
             CALL ppnd7((i - th)/an25, a(i))
             summ2 = summ2 + a(i) ** 2
          END DO
          summ2 = summ2 * two
          ssumm2 = SQRT(summ2)
          rsn = one / SQRT(an)
          a1 = poly(c1, 6, rsn) - a(1) / ssumm2
          ! Normalize coefficients
          IF (n > 5) THEN
             i1 = 3
             a2 = -a(2)/ssumm2 + poly(c2,6,rsn)
             fac = SQRT((summ2 - two * a(1)**2 - two * a(2)**2)/(one - two * a1**2 - two * a2**2))
             a(1) = a1
             a(2) = a2
          ELSE
             i1 = 2
             fac = SQRT((summ2 - two * a(1)**2)/ (one - two * a1**2))
             a(1) = a1
          END IF
          DO i = i1, n2
             a(i) = -a(i)/fac
          END DO
       END IF
       ! scaled X
       xx = x(1)/range
       sx = xx
       sa = -a(1)
       j = n - 1
       DO i = 2, n
          xi = x(i)/range
          sx = sx + xi
          IF (i /= j) sa = sa + SIGN(1, i - j) * a(MIN(i, j))
          xx = xi
          j = j - 1
       END DO
       ! Calculate W statistic as squared correlation
       ! between data and coefficients
       sa = sa/n
       sx = sx/n
       ssa = zero
       ssx = zero
       sax = zero
       j = n
       DO i = 1, n
          IF (i /= j) THEN
             asa = SIGN(1, i - j) * a(MIN(i, j)) - sa
          ELSE
             asa = -sa
          END IF
          xsx = x(i)/range - sx
          ssa = ssa + asa * asa
          ssx = ssx + xsx * xsx
          sax = sax + asa * xsx
          j = j - 1
       END DO
       ! W1 equals (1-W) calculated to avoid excessive rounding error
       ! for W very near 1 (a potential problem in very large samples)
       ssassx = SQRT(ssa * ssx)
       w1 = (ssassx - sax) * (ssassx + sax)/(ssa * ssx)
       w = one - w1
       ! Calculate significance level for W (exact for N=3)
       IF (n == 3) THEN
          pw = pi6 * (ASIN(SQRT(w)) - stqr)
       ELSE
          y  = LOG(w1)
          xx = LOG(an)
          m = zero
          s = one
          IF (n <= 11) THEN
             gamma = poly(g, 2, an)
             IF (y >= gamma) THEN
                pw = small
             ELSE
                y = -LOG(gamma - y)
                m = poly(c3, 4, an)
                s = EXP(poly(c4, 4, an))
                pw = alnorm((y - m)/s, .TRUE.)
             END IF
          ELSE
             m = poly(c5, 4, xx)
             s = EXP(poly(c6, 3, xx))
             pw = alnorm((y - m)/s, .TRUE.)
          END IF
       END IF
    END IF
    DEALLOCATE(x,stat=stat)
    CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    DEALLOCATE(itmp,stat=stat)
    CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    DEALLOCATE(a,stat=stat)
    CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)

  END SUBROUTINE sw_test

! *****************************************************************************
!> \brief Produces the normal deviate Z corresponding to a given lower tail area of P
!>      Z is accurate to about 1 part in 10**7.
!>      AS241  APPL. STATIST. (1988) VOL. 37, NO. 3, 477- 484.
!> \par History
!>      Original version by Alain J. Miller - 1996
!>      Teodoro Laino (02.2007) [tlaino]
! *****************************************************************************
  SUBROUTINE ppnd7 (p, normal_dev)
    REAL(KIND=dp), INTENT(IN)                :: p
    REAL(KIND=dp), INTENT(OUT)               :: normal_dev

    REAL(KIND=dp), PARAMETER :: a0 = 3.3871327179E+00_dp, &
      a1 = 5.0434271938E+01_dp, a2 = 1.5929113202E+02_dp, &
      a3 = 5.9109374720E+01_dp, b1 = 1.7895169469E+01_dp, &
      b2 = 7.8757757664E+01_dp, b3 = 6.7187563600E+01_dp, &
      c0 = 1.4234372777E+00_dp, c1 = 2.7568153900E+00_dp, &
      c2 = 1.3067284816E+00_dp, c3 = 1.7023821103E-01_dp, &
      const1 = 0.180625_dp, const2 = 1.6_dp, d1 = 7.3700164250E-01_dp, &
      d2 = 1.2021132975E-01_dp, e0 = 6.6579051150E+00_dp, &
      e1 = 3.0812263860E+00_dp, e2 = 4.2868294337E-01_dp, &
      e3 = 1.7337203997E-02_dp, f1 = 2.4197894225E-01_dp, &
      f2 = 1.2258202635E-02_dp, half = 0.5_dp, one = 1.0_dp
    REAL(KIND=dp), PARAMETER :: split1 = 0.425_dp, split2 = 5.0_dp, &
      zero = 0.0_dp

    REAL(KIND=dp)                            :: q, r

    q = p - half
    IF (ABS(q) <= split1) THEN
       r = const1 - q * q
       normal_dev = q * (((a3 * r + a2) * r + a1) * r + a0) / &
            (((b3 * r + b2) * r + b1) * r + one)
       RETURN
    ELSE
       IF (q < zero) THEN
          r = p
       ELSE
          r = one - p
       END IF
       IF (r <= zero) THEN
          normal_dev = zero
          RETURN
       END IF
       r = SQRT(-LOG(r))
       IF (r <= split2) THEN
          r = r - const2
          normal_dev = (((c3 * r + c2) * r + c1) * r + c0) / ((d2 * r + d1) * r + one)
       ELSE
          r = r - split2
          normal_dev = (((e3 * r + e2) * r + e1) * r + e0) / ((f2 * r + f1) * r + one)
       END IF
       IF (q < zero) normal_dev = - normal_dev
       RETURN
    END IF
  END SUBROUTINE ppnd7

! *****************************************************************************
!> \brief Evaluates the tail area of the standardised normal curve
!>      from x to infinity if upper is .true. or
!>      from minus infinity to x if upper is .false.
!>      AS66 Applied Statistics (1973) vol.22, no.3
!> \par History
!>      Original version by Alain J. Miller - 1996
!>      Teodoro Laino (02.2007) [tlaino]
! *****************************************************************************
  FUNCTION alnorm(x, upper) RESULT(fn_val)
    REAL(KIND=dp), INTENT(IN)                :: x
    LOGICAL, INTENT(IN)                      :: upper
    REAL(KIND=dp)                            :: fn_val

    REAL(KIND=dp), PARAMETER :: a1 = 5.75885480458_dp, a2 = 2.62433121679_dp, &
      a3 = 5.92885724438_dp, b1 = -29.8213557807_dp, b2 = 48.6959930692_dp, &
      c1 = -3.8052E-8_dp, c2 = 3.98064794E-4_dp, c3 = -0.151679116635_dp, &
      c4 = 4.8385912808_dp, c5 = 0.742380924027_dp, c6 = 3.99019417011_dp, &
      con = 1.28_dp, d1 = 1.00000615302_dp, d2 = 1.98615381364_dp, &
      d3 = 5.29330324926_dp, d4 = -15.1508972451_dp, d5 = 30.789933034_dp, &
      half = 0.5_dp, ltone = 7.0_dp, one = 1.0_dp, p = 0.398942280444_dp, &
      q = 0.39990348504_dp, r = 0.398942280385_dp, utzero = 18.66_dp, &
      zero = 0.0_dp

    LOGICAL                                  :: up
    REAL(KIND=dp)                            :: y, z

    up = upper
    z = x
    IF(z >=  zero) GO TO 10
    up = .NOT. up
    z = -z
10  IF(z <= ltone .OR. up .AND. z <= utzero) GO TO 20
    fn_val = zero
    GO TO 40
20  y = half*z*z
    IF(z > con) GO TO 30

    fn_val = half - z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))))
    GO TO 40
30  fn_val = r*EXP(-y)/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))))
40  IF(.NOT. up) fn_val = one - fn_val

  END FUNCTION alnorm

! *****************************************************************************
!> \brief Calculates the algebraic polynomial of order nored-1 with
!>      array of coefficients c.  Zero order coefficient is c(1)
!>      AS 181.2   Appl. Statist.  (1982) Vol. 31, No. 2
!> \par History
!>      Original version by Alain J. Miller - 1996
!>      Teodoro Laino (02.2007) [tlaino]
! *****************************************************************************
  FUNCTION poly(c, nord, x) RESULT(fn_val)

    REAL(KIND=dp), INTENT(IN)                :: c(:)
    INTEGER, INTENT(IN)                      :: nord
    REAL(KIND=dp), INTENT(IN)                :: x
    REAL(KIND=dp)                            :: fn_val

    INTEGER                                  :: i, j, n2
    REAL(KIND=dp)                            :: p

    fn_val = c(1)
    IF (nord == 1) RETURN
    p = x*c(nord)
    IF (nord == 2) GO TO 20
    n2 = nord - 2
    j = n2 + 1
    DO i = 1,n2
       p = (p + c(j))*x
       j = j - 1
    END DO
20  fn_val = fn_val + p
  END FUNCTION poly

! *****************************************************************************
!> \brief Kandall's test for correlation
!> \note
!>      tau:  Kendall's Tau
!>      z:    number of std devs from 0 of tau
!>      prob: tau's probability
!> \par History
!>      Teodoro Laino (02.2007) [tlaino]
! *****************************************************************************
  SUBROUTINE k_test(xdata,istart,n,tau,z,prob,error)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xdata
    INTEGER, INTENT(IN)                      :: istart, n
    REAL(KIND=dp)                            :: tau, z, prob
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'k_test', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: is, j, k, nt
    REAL(KIND=dp)                            :: a1, var

    nt  = n-istart+1
    IF (nt.GE.min_sample_size) THEN
       is=0
       DO j=istart,n-1
          DO k=j+1,n
             a1=xdata(j)-xdata(k)
             IF      (a1.GT.0.0_dp)THEN
                is=is+1
             ELSE IF (a1.LT.0.0_dp) THEN
                is=is-1
             END IF
          END DO
       END DO
       tau = REAL(is,KIND=dp)
       var = REAL(nt,KIND=dp)*REAL(nt-1,KIND=dp)*REAL(2*nt+5,KIND=dp)/18.0_dp
       z   = tau/SQRT(var)
       prob= erf(ABS(z)/SQRT(2.0_dp))
    ELSE
       tau  = 0.0_dp
       z    = 0.0_dp
       prob = 1.0_dp
    END IF
  END SUBROUTINE k_test

! *****************************************************************************
!> \brief Von Neumann test for serial correlation
!> \par History
!>      Teodoro Laino (02.2007) [tlaino]
! *****************************************************************************
  SUBROUTINE vn_test(xdata,n,r,u,prob,error)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xdata
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp)                            :: r, u, prob
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'vn_test', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: q, s, var, x

    IF (n.GE.min_sample_size) THEN
       x = 0.0_dp
       q = 0.0_dp
       s = 0.0_dp
       DO i = 1, n-1
          x = x + xdata(i)
          q = q + (xdata(i+1)-xdata(i))**2
       END DO
       x = x + xdata(n)
       x = x/REAL(n,KIND=dp)
       DO i = 1, n
          s = s + (xdata(i)-x)**2
       END DO
       s = s / REAL(n-1,KIND=dp)
       q = q / REAL(2*(n-1),KIND=dp)
       r = q/s
       var = SQRT(1.0_dp/REAL(n+1,KIND=dp)*(1.0_dp+1.0_dp/REAL(n-1,KIND=dp)))
       u = (r-1.0_dp)/var
       prob= erf(ABS(u)/SQRT(2.0_dp))
    ELSE
       r = 0.0_dp
       u = 0.0_dp
       prob = 1.0_dp
    END IF

  END SUBROUTINE vn_test

! *****************************************************************************
!> \brief Performs tests on statistical methods
!>      Debug use only
!> \par History
!>      Teodoro Laino (02.2007) [tlaino]
! *****************************************************************************
  SUBROUTINE tests(xdata, globenv, error)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xdata
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tests', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, n, stat
    LOGICAL                                  :: failure
    REAL(KINd=dp)                            :: prob, pw, r, tau, u, w, z
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ydata

    failure = .FALSE.
    IF (debug_this_module) THEN
       n = 50 ! original sample size
       NULLIFY(xdata)
       ALLOCATE(xdata(n), stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       DO i = 1, 10
          xdata(i) = 5.0_dp - REAL(i,KIND=dp)/2.0_dp + 0.1*&
                     next_random_number(globenv%gaussian_rng_stream,error=error)
          WRITE(3,*)xdata(i)
       END DO
       DO i = 11, n
          xdata(i) = 0.1*next_random_number(globenv%gaussian_rng_stream,error=error)
       END DO

       ! Test for trend
       DO i = 1, n
          CALL k_test(xdata,i,n,tau,z,prob,error)
          IF (prob<=0.2_dp) EXIT
       END DO
       WRITE(*,*)"Mann-Kendall test",i

       ! Test for normality distribution and for serial correlation
       DO i = 1,n
         ALLOCATE(ydata(n-i+1), stat=stat)
         CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
         ydata = xdata(i:n)
         CALL sw_test ( ydata, n-i+1, w, pw, error)
         CALL vn_test(ydata,n-i+1,r,u,prob,error)
         WRITE(*,*)"Shapiro Wilks test",i,w,pw
         WRITE(*,*)"Von Neu",i,r,u,prob
         DEALLOCATE(ydata, stat=stat)
         CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END DO

       DEALLOCATE(xdata, stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF
  END SUBROUTINE tests

END MODULE statistical_methods

