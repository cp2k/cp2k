!-----------------------------------------------------------------------
SUBROUTINE WRITE_DOCUMENTATION()
  write(6,'(A)') &
" scriptmini : a fortran program to minimise a script                  ", &
"               or external program.                                   ", &
"                                                                      ", &
"     Written by Joost VandeVondele                                    ", &
"                                                                      ", &
"     the script is treated as a black box, that given an              ", &
"     input vector x returns the function value f.                     ", &
"                                                                      ", &
" usage :                                                              ", &
"                                                                      ", &
"     ./scriptmini                                                     ", &
"                                                                      ", &
" inputs :                                                             ", &
"                                                                      ", &
"     -) the script should be called  scriptmini_eval                  ", &
"     -) the script should read the n real variables                   ", &
"        from  scriptmini_eval.in                                      ", &
"     -) the script should write the function value (1 real number)    ", &
"        to 'scriptmini_eval.out'                                      ", &
"     -) input of scriptmini is 'scriptmini.in' with the format:       ", &
"        N                                                             ", &
"        rhobeg rhoend                                                 ", &
"        maxfun                                                        ", &
"        iprint                                                        ", &
"        x[1] x[2] x[3] ... x[N]                                       ", &
"                                                                      ", &
"        where:                                                        ", &
"        N      : integer : is the number of variables                 ", &
"        rhobeg : real    : initial trust region radius, +- 10% of the ", &
"                           largest expected change in the variables   ", &
"        rhoend : real    : final trust region radius, +- the final    ", &
"                           uncertainty in the variables               ", &
"        maxfun : integer : the maximum number of calls to             ", &
"                           scriptmini_eval [O(10*N**2)]               ", &
"        iprint : integer : output level (0-3)                         ", &
"                           0 : no output at all (!)                   ", &
"                           ...                                        ", &
"                           3 : info at every step (recommended)       ", &
"        x[...] : real    : the initial values of the variables        "
END SUBROUTINE WRITE_DOCUMENTATION
!-------------------------------------------------------------------------------

MODULE Powell_Optimize

! Code converted using TO_F90 by Alan Miller
! Date: 2002-11-09  Time: 16:58:08

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

PRIVATE
PUBLIC  :: uobyqa


CONTAINS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uobyqa.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
SUBROUTINE uobyqa(n, x, rhobeg, rhoend, iprint, maxfun)

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(:)
REAL (dp), INTENT(IN)      :: rhobeg
REAL (dp), INTENT(IN)      :: rhoend
INTEGER, INTENT(IN)        :: iprint
INTEGER, INTENT(IN)        :: maxfun

!     This subroutine seeks the least value of a function of many variables,
!     by a trust region method that forms quadratic models by interpolation.
!     The algorithm is described in "UOBYQA: unconstrained optimization by
!     quadratic approximation" by M.J.D. Powell, Report DAMTP 2000/NA14,
!     University of Cambridge. The arguments of the subroutine are as follows.

!     N must be set to the number of variables and must be at least two.
!     Initial values of the variables must be set in X(1),X(2),...,X(N). They
!       will be changed to the values that give the least calculated F.
!     RHOBEG and RHOEND must be set to the initial and final values of a trust
!       region radius, so both must be positive with RHOEND<=RHOBEG. Typically
!       RHOBEG should be about one tenth of the greatest expected change to a
!       variable, and RHOEND should indicate the accuracy that is required in
!       the final values of the variables.
!     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
!       amount of printing.  Specifically, there is no output if IPRINT=0 and
!       there is output only at the return if IPRINT=1.  Otherwise, each new
!       value of RHO is printed, with the best vector of variables so far and
!       the corresponding value of the objective function. Further, each new
!       value of F with its variables are output if IPRINT=3.
!     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
!     The array W will be used for working space. Its length must be at least
!       ( N**4 + 8*N**3 + 23*N**2 + 42*N + max [ 2*N**2 + 4, 18*N ] ) / 4.

!     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
!     the value of the objective function for the variables X(1),X(2),...,X(N).

INTEGER  :: npt

!     Partition the working space array, so that different parts of it can be
!     treated separately by the subroutine that performs the main calculation.

npt = (n*n + 3*n + 2) / 2
CALL uobyqb(n, x, rhobeg, rhoend, iprint, maxfun, npt)
RETURN
END SUBROUTINE uobyqa

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uobyqb.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE uobyqb(n, x, rhobeg, rhoend, iprint, maxfun, npt)

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(:)
REAL (dp), INTENT(IN)      :: rhobeg
REAL (dp), INTENT(IN)      :: rhoend
INTEGER, INTENT(IN)        :: iprint
INTEGER, INTENT(IN)        :: maxfun
INTEGER, INTENT(IN)        :: npt

INTERFACE
  SUBROUTINE calfun(n, x, f)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)    :: n
    REAL (dp), INTENT(IN)  :: x(:)
    REAL (dp), INTENT(OUT) :: f
  END SUBROUTINE calfun
END INTERFACE

! The following arrays were previously passed as arguments:

REAL (dp)  :: xbase(n), xopt(n), xnew(n), xpt(npt,n), pq(npt-1)
REAL (dp)  :: pl(npt,npt-1), h(n,n), g(n), d(n), vlag(npt), w(npt)

!     The arguments N, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical to
!       the corresponding arguments in SUBROUTINE UOBYQA.

!     NPT is set by UOBYQA to (N*N+3*N+2)/2 for the above dimension statement.
!     XBASE will contain a shift of origin that reduces the contributions from
!       rounding errors to values of the model and Lagrange functions.
!     XOPT will be set to the displacement from XBASE of the vector of
!       variables that provides the least calculated F so far.
!     XNEW will be set to the displacement from XBASE of the vector of
!       variables for the current calculation of F.
!     XPT will contain the interpolation point coordinates relative to XBASE.
!     PQ will contain the parameters of the quadratic model.
!     PL will contain the parameters of the Lagrange functions.
!     H will provide the second derivatives that TRSTEP and LAGMAX require.
!     G will provide the first derivatives that TRSTEP and LAGMAX require.
!     D is reserved for trial steps from XOPT, except that it will contain
!       diagonal second derivatives during the initialization procedure.
!     VLAG will contain the values of the Lagrange functions at a new point X.
!     The array W will be used for working space.

REAL (dp)  :: half = 0.5_dp, one = 1.0_dp, tol = 0.01_dp, two = 2.0_dp
REAL (dp)  :: zero = 0.0_dp
REAL (dp)  :: ddknew, delta, detrat, diff, distest, dnorm, errtol, estim
REAL (dp)  :: evalue, f, fbase, fopt, fsave, ratio, rho, rhosq, sixthm
REAL (dp)  :: sum, sumg, sumh, temp, tempa, tworsq, vmax, vquad, wmult
INTEGER    :: i, ih, ip, iq, iw, j, jswitch, k, knew, kopt, ksave, ktemp
INTEGER    :: nf, nftest, nnp, nptm

!     Set some constants.

nnp = n + n + 1
nptm = npt - 1
nftest = MAX(maxfun,1)

!     Initialization. NF is the number of function calculations so far.

rho = rhobeg
rhosq = rho * rho
nf = 0
DO  i = 1, n
  xbase(i) = x(i)
  xpt(1:npt,i) = zero
END DO
pl(1:npt,1:nptm) = zero

!     The branch to label 120 obtains a new value of the objective function
!     and then there is a branch back to label 50, because the new function
!     value is needed to form the initial quadratic model. The least function
!     value so far and its index are noted below.

50 x(1:n) = xbase(1:n) + xpt(nf+1,1:n)
GO TO 150

70 IF (nf == 1) THEN
  fopt = f
  kopt = nf
  fbase = f
  j = 0
  jswitch = -1
  ih = n
ELSE
  IF (f < fopt) THEN
    fopt = f
    kopt = nf
  END IF
END IF

!     Form the gradient and diagonal second derivatives of the initial
!     quadratic model and Lagrange functions.

IF (nf <= nnp) THEN
  jswitch = -jswitch
  IF (jswitch > 0) THEN
    IF (j >= 1) THEN
      ih = ih + j
      IF (w(j) < zero) THEN
        d(j) = (fsave+f-two*fbase) / rhosq
        pq(j) = (fsave-f) / (two*rho)
        pl(1,ih) = -two / rhosq
        pl(nf-1,j) = half / rho
        pl(nf-1,ih) = one / rhosq
      ELSE
        pq(j) = (4.0D0*fsave-3.0D0*fbase-f) / (two*rho)
        d(j) = (fbase+f-two*fsave) / rhosq
        pl(1,j) = -1.5D0 / rho
        pl(1,ih) = one / rhosq
        pl(nf-1,j) = two / rho
        pl(nf-1,ih) = -two / rhosq
      END IF
      pq(ih) = d(j)
      pl(nf,j) = -half / rho
      pl(nf,ih) = one / rhosq
    END IF
    
!     Pick the shift from XBASE to the next initial interpolation point
!     that provides diagonal second derivatives.
    
    IF (j < n) THEN
      j = j + 1
      xpt(nf+1,j) = rho
    END IF
  ELSE
    fsave = f
    IF (f < fbase) THEN
      w(j) = rho
      xpt(nf+1,j) = two * rho
    ELSE
      w(j) = -rho
      xpt(nf+1,j) = -rho
    END IF
  END IF
  IF (nf < nnp) GO TO 50
  
!     Form the off-diagonal second derivatives of the initial quadratic model.
  
  ih = n
  ip = 1
  iq = 2
END IF
ih = ih + 1
IF (nf > nnp) THEN
  temp = one / (w(ip)*w(iq))
  tempa = f - fbase - w(ip) * pq(ip) - w(iq) * pq(iq)
  pq(ih) = (tempa - half*rhosq*(d(ip)+d(iq))) * temp
  pl(1,ih) = temp
  iw = ip + ip
  IF (w(ip) < zero) iw = iw + 1
  pl(iw,ih) = -temp
  iw = iq + iq
  IF (w(iq) < zero) iw = iw + 1
  pl(iw,ih) = -temp
  pl(nf,ih) = temp
  
!     Pick the shift from XBASE to the next initial interpolation point
!     that provides off-diagonal second derivatives.
  
  ip = ip + 1
END IF
IF (ip == iq) THEN
  ih = ih + 1
  ip = 1
  iq = iq + 1
END IF
IF (nf < npt) THEN
  xpt(nf+1,ip) = w(ip)
  xpt(nf+1,iq) = w(iq)
  GO TO 50
END IF

!     Set parameters to begin the iterations for the current RHO.

sixthm = zero
delta = rho
80 tworsq = (two*rho) ** 2
rhosq = rho * rho

!     Form the gradient of the quadratic model at the trust region centre.

90 knew = 0
ih = n
DO  j = 1, n
  xopt(j) = xpt(kopt,j)
  g(j) = pq(j)
  DO  i = 1, j
    ih = ih + 1
    g(i) = g(i) + pq(ih) * xopt(j)
    IF (i < j) g(j) = g(j) + pq(ih) * xopt(i)
    h(i,j) = pq(ih)
  END DO
END DO

!     Generate the next trust region step and test its length.  Set KNEW
!     to -1 if the purpose of the next F will be to improve conditioning,
!     and also calculate a lower bound on the Hessian term of the model Q.

CALL trstep(n, g, h, delta, tol, d, evalue)
temp = zero
DO i = 1, n
  temp = temp + d(i)**2
END DO
dnorm = MIN(delta,SQRT(temp))
errtol = -one
IF (dnorm < half*rho) THEN
  knew = -1
  errtol = half * evalue * rho * rho
  IF (nf <= npt+9) errtol = zero
  GO TO 290
END IF

!     Calculate the next value of the objective function.

130 DO  i = 1, n
  xnew(i) = xopt(i) + d(i)
  x(i) = xbase(i) + xnew(i)
END DO
150 IF (nf >= nftest) THEN
  IF (iprint > 0) WRITE(*, 5000)
  GO TO 420
END IF
nf = nf + 1
CALL calfun(n, x, f)
IF (iprint == 3) THEN
  WRITE(*, 5100) nf, f, x(1:n)
END IF
IF (nf <= npt) GO TO 70
IF (knew == -1) GO TO 420

!     Use the quadratic model to predict the change in F due to the step D,
!     and find the values of the Lagrange functions at the new point.

vquad = zero
ih = n
DO  j = 1, n
  w(j) = d(j)
  vquad = vquad + w(j) * pq(j)
  DO  i = 1, j
    ih = ih + 1
    w(ih) = d(i) * xnew(j) + d(j) * xopt(i)
    IF (i == j) w(ih) = half * w(ih)
    vquad = vquad + w(ih) * pq(ih)
  END DO
END DO
DO  k = 1, npt
  temp = zero
  DO  j = 1, nptm
    temp = temp + w(j) * pl(k,j)
  END DO
  vlag(k) = temp
END DO
vlag(kopt) = vlag(kopt) + one

!     Update SIXTHM, which is a lower bound on one sixth of the greatest
!     third derivative of F.

diff = f - fopt - vquad
sum = zero
DO  k = 1, npt
  temp = zero
  DO  i = 1, n
    temp = temp + (xpt(k,i)-xnew(i)) ** 2
  END DO
  temp = SQRT(temp)
  sum = sum + ABS(temp*temp*temp*vlag(k))
END DO
sixthm = MAX(sixthm, ABS(diff)/sum)

!     Update FOPT and XOPT if the new F is the least value of the objective
!     function so far.  Then branch if D is not a trust region step.

fsave = fopt
IF (f < fopt) THEN
  fopt = f
  xopt(1:n) = xnew(1:n)
END IF
ksave = knew
IF (knew <= 0) THEN
  
!     Pick the next value of DELTA after a trust region step.
  
  IF (vquad >= zero) THEN
    IF (iprint > 0) WRITE(*, 5200)
    GO TO 420
  END IF
  ratio = (f-fsave) / vquad
  IF (ratio <= 0.1D0) THEN
    delta = half * dnorm
  ELSE IF (ratio <= 0.7D0) THEN
    delta = MAX(half*delta,dnorm)
  ELSE
    delta = MAX(delta, 1.25D0*dnorm, dnorm+rho)
  END IF
  IF (delta <= 1.5D0*rho) delta = rho
  
!     Set KNEW to the index of the next interpolation point to be deleted.
  
  ktemp = 0
  detrat = zero
  IF (f >= fsave) THEN
    ktemp = kopt
    detrat = one
  END IF
  DO  k = 1, npt
    sum = zero
    DO  i = 1, n
      sum = sum + (xpt(k,i)-xopt(i)) ** 2
    END DO
    temp = ABS(vlag(k))
    IF (sum > rhosq) temp = temp * (sum/rhosq) ** 1.5D0
    IF (temp > detrat .AND. k /= ktemp) THEN
      detrat = temp
      ddknew = sum
      knew = k
    END IF
  END DO
  IF (knew == 0) GO TO 290
END IF

!     Replace the interpolation point that has index KNEW by the point XNEW,
!     and also update the Lagrange functions and the quadratic model.

DO  i = 1, n
  xpt(knew,i) = xnew(i)
END DO
temp = one / vlag(knew)
DO  j = 1, nptm
  pl(knew,j) = temp * pl(knew,j)
  pq(j) = pq(j) + diff * pl(knew,j)
END DO
DO  k = 1, npt
  IF (k /= knew) THEN
    temp = vlag(k)
    DO  j = 1, nptm
      pl(k,j) = pl(k,j) - temp * pl(knew,j)
    END DO
  END IF
END DO

!     Update KOPT if F is the least calculated value of the objective function.
!     Then branch for another trust region calculation.  The case KSAVE > 0
!     indicates that a model step has just been taken.

IF (f < fsave) THEN
  kopt = knew
  GO TO 90
END IF
IF (ksave > 0) GO TO 90
IF (dnorm > two*rho) GO TO 90
IF (ddknew > tworsq) GO TO 90

!     Alternatively, find out if the interpolation points are close
!     enough to the best point so far.

290 DO  k = 1, npt
  w(k) = zero
  DO  i = 1, n
    w(k) = w(k) + (xpt(k,i)-xopt(i)) ** 2
  END DO
END DO
320 knew = -1
distest = tworsq
DO  k = 1, npt
  IF (w(k) > distest) THEN
    knew = k
    distest = w(k)
  END IF
END DO

!     If a point is sufficiently far away, then set the gradient and Hessian
!     of its Lagrange function at the centre of the trust region, and find
!     half the sum of squares of components of the Hessian.

IF (knew > 0) THEN
  ih = n
  sumh = zero
  DO  j = 1, n
    g(j) = pl(knew,j)
    DO  i = 1, j
      ih = ih + 1
      temp = pl(knew,ih)
      g(j) = g(j) + temp * xopt(i)
      IF (i < j) THEN
        g(i) = g(i) + temp * xopt(j)
        sumh = sumh + temp * temp
      END IF
      h(i,j) = temp
    END DO
    sumh = sumh + half * temp * temp
  END DO
  
!     If ERRTOL is positive, test whether to replace the interpolation point
!     with index KNEW, using a bound on the maximum modulus of its Lagrange
!     function in the trust region.
  
  IF (errtol > zero) THEN
    w(knew) = zero
    sumg = zero
    DO  i = 1, n
      sumg = sumg + g(i) ** 2
    END DO
    estim = rho * (SQRT(sumg)+rho*SQRT(half*sumh))
    wmult = sixthm * distest ** 1.5D0
    IF (wmult*estim <= errtol) GO TO 320
  END IF
  
!     If the KNEW-th point may be replaced, then pick a D that gives a large
!     value of the modulus of its Lagrange function within the trust region.
!     Here the vector XNEW is used as temporary working space.
  
  CALL lagmax(n, g, h, rho, d, xnew, vmax)
  IF (errtol > zero) THEN
    IF (wmult*vmax <= errtol) GO TO 320
  END IF
  GO TO 130
END IF
IF (dnorm > rho) GO TO 90

!     Prepare to reduce RHO by shifting XBASE to the best point so far,
!     and make the corresponding changes to the gradients of the Lagrange
!     functions and the quadratic model.

IF (rho > rhoend) THEN
  ih = n
  DO  j = 1, n
    xbase(j) = xbase(j) + xopt(j)
    DO  k = 1, npt
      xpt(k,j) = xpt(k,j) - xopt(j)
    END DO
    DO  i = 1, j
      ih = ih + 1
      pq(i) = pq(i) + pq(ih) * xopt(j)
      IF (i < j) THEN
        pq(j) = pq(j) + pq(ih) * xopt(i)
        DO  k = 1, npt
          pl(k,j) = pl(k,j) + pl(k,ih) * xopt(i)
        END DO
      END IF
      DO  k = 1, npt
        pl(k,i) = pl(k,i) + pl(k,ih) * xopt(j)
      END DO
    END DO
  END DO
  
!     Pick the next values of RHO and DELTA.
  
  delta = half * rho
  ratio = rho / rhoend
  IF (ratio <= 16.0D0) THEN
    rho = rhoend
  ELSE IF (ratio <= 250.0D0) THEN
    rho = SQRT(ratio) * rhoend
  ELSE
    rho = 0.1D0 * rho
  END IF
  delta = MAX(delta,rho)
  IF (iprint >= 2) THEN
    IF (iprint >= 3) WRITE(*, 5300)
    WRITE(*, 5400) rho, nf
    WRITE(*, 5500) fopt, xbase(1:n)
  END IF
  GO TO 80
END IF

!     Return from the calculation, after another Newton-Raphson step, if
!     it is too short to have been tried before.

IF (errtol >= zero) GO TO 130
420 IF (fopt <= f) THEN
  DO  i = 1, n
    x(i) = xbase(i) + xopt(i)
  END DO
  f = fopt
END IF
IF (iprint >= 1) THEN
  WRITE(*, 5600) nf
  WRITE(*, 5500) f, x(1:n)
END IF
RETURN

5000 FORMAT (/T5, 'Return from UOBYQA because CALFUN has been',  &
    ' called MAXFUN times')
5100 FORMAT (/T5, 'Function number',i6,'    F =', g18.10,  &
    '    The corresponding X is:'/ (t3, 5g15.6))
5200 FORMAT (/T5, 'Return from UOBYQA because a trust',  &
    ' region step has failed to reduce Q')
5300 FORMAT (' ')
5400 FORMAT (/T5, 'New RHO =', g11.4, '     Number of function values =',i6)
5500 FORMAT (T5, 'Least value of F =', g23.15,  &
    '         The corresponding X is:'/ (t3, 5g15.6))
5600 FORMAT (/T5, 'At the return from UOBYQA',  &
    '     Number of function values =', i6)
END SUBROUTINE uobyqb

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% trstep.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE trstep(n, g, h, delta, tol, d, evalue)

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: g(:)
REAL (dp), INTENT(IN OUT)  :: h(:,:)
REAL (dp), INTENT(IN)      :: delta
REAL (dp), INTENT(IN)      :: tol
REAL (dp), INTENT(OUT)     :: d(:)
REAL (dp), INTENT(OUT)     :: evalue

!     N is the number of variables of a quadratic objective function, Q say.
!     G is the gradient of Q at the origin.
!     H is the Hessian matrix of Q.  Only the upper triangular and diagonal
!       parts need be set.  The lower triangular part is used to store the
!       elements of a Householder similarity transformation.
!     DELTA is the trust region radius, and has to be positive.
!     TOL is the value of a tolerance from the open interval (0,1).
!     D will be set to the calculated vector of variables.

!     EVALUE will be set to the least eigenvalue of H if and only if D is a
!     Newton-Raphson step.  Then EVALUE will be positive, but otherwise it
!     will be set to zero.

!     Let MAXRED be the maximum of Q(0)-Q(D) subject to ||D|| <= DELTA,
!     and let ACTRED be the value of Q(0)-Q(D) that is actually calculated.
!     We take the view that any D is acceptable if it has the properties

!             ||D|| <= DELTA  and  ACTRED <= (1-TOL)*MAXRED.

!     The calculation of D is done by the method of Section 2 of the paper
!     by MJDP in the 1997 Dundee Numerical Analysis Conference Proceedings,
!     after transforming H to tridiagonal form.

!     The arrays GG, TD, TN, W, PIV and Z will be used for working space.
REAL (dp)  :: gg(n), td(n), tn(n), w(n), piv(n), z(n)

REAL (dp)  :: delsq, dhd, dnorm, dsq, dtg, dtz, gam, gnorm, gsq, hnorm
REAL (dp)  :: par, parl, parlest, paru, paruest, phi, phil, phiu, pivksv
REAL (dp)  :: pivot, posdef, scale, shfmax, shfmin, shift, slope, sum
REAL (dp)  :: tdmin, temp, tempa, tempb, wsq, wwsq, wz, zsq
INTEGER    :: i, iterc, j, jp, k, kp, kpp, ksav, ksave, nm
REAL (dp)  :: one = 1.0_dp, two = 2.0_dp, zero = 0.0_dp

!     Initialization.

delsq = delta * delta
evalue = zero
nm = n - 1
DO  i = 1, n
  d(i) = zero
  td(i) = h(i,i)
  DO  j = 1, i
    h(i,j) = h(j,i)
  END DO
END DO

!     Apply Householder transformations to obtain a tridiagonal matrix that
!     is similar to H, and put the elements of the Householder vectors in
!     the lower triangular part of H.  Further, TD and TN will contain the
!     diagonal and other nonzero elements of the tridiagonal matrix.

DO  k = 1, nm
  kp = k + 1
  sum = zero
  IF (kp < n) THEN
    kpp = kp + 1
    DO  i = kpp, n
      sum = sum + h(i,k) ** 2
    END DO
  END IF
  IF (sum == zero) THEN
    tn(k) = h(kp,k)
    h(kp,k) = zero
  ELSE
    temp = h(kp,k)
    tn(k) = SIGN(SQRT(sum+temp*temp),temp)
    h(kp,k) = -sum / (temp+tn(k))
    temp = SQRT(two/(sum+h(kp,k)**2))
    DO  i = kp, n
      w(i) = temp * h(i,k)
      h(i,k) = w(i)
      z(i) = td(i) * w(i)
    END DO
    wz = zero
    DO  j = kp, nm
      jp = j + 1
      DO  i = jp, n
        z(i) = z(i) + h(i,j) * w(j)
        z(j) = z(j) + h(i,j) * w(i)
      END DO
      wz = wz + w(j) * z(j)
    END DO
    wz = wz + w(n) * z(n)
    DO  j = kp, n
      td(j) = td(j) + w(j) * (wz*w(j)-two*z(j))
      IF (j < n) THEN
        jp = j + 1
        DO  i = jp, n
          h(i,j) = h(i,j) - w(i) * z(j) - w(j) * (z(i)-wz*w(i))
        END DO
      END IF
    END DO
  END IF
END DO

!     Form GG by applying the similarity transformation to G.

gsq = zero
DO  i = 1, n
  gg(i) = g(i)
  gsq = gsq + g(i) ** 2
END DO
gnorm = SQRT(gsq)
DO  k = 1, nm
  kp = k + 1
  sum = zero
  DO  i = kp, n
    sum = sum + gg(i) * h(i,k)
  END DO
  DO  i = kp, n
    gg(i) = gg(i) - sum * h(i,k)
  END DO
END DO

!     Begin the trust region calculation with a tridiagonal matrix by
!     calculating the norm of H. Then treat the case when H is zero.

hnorm = ABS(td(1)) + ABS(tn(1))
tdmin = td(1)
tn(n) = zero
DO  i = 2, n
  temp = ABS(tn(i-1)) + ABS(td(i)) + ABS(tn(i))
  hnorm = MAX(hnorm,temp)
  tdmin = MIN(tdmin,td(i))
END DO
IF (hnorm == zero) THEN
  IF (gnorm == zero) GO TO 420
  scale = delta / gnorm
  DO  i = 1, n
    d(i) = -scale * gg(i)
  END DO
  GO TO 380
END IF

!     Set the initial values of PAR and its bounds.

parl = MAX(zero, -tdmin, gnorm/delta-hnorm)
parlest = parl
par = parl
paru = zero
paruest = zero
posdef = zero
iterc = 0

!     Calculate the pivots of the Cholesky factorization of (H+PAR*I).

160 iterc = iterc + 1
ksav = 0
piv(1) = td(1) + par
k = 1
170 IF (piv(k) > zero) THEN
  piv(k+1) = td(k+1) + par - tn(k) ** 2 / piv(k)
ELSE
  IF (piv(k) < zero .OR. tn(k) /= zero) GO TO 180
  ksav = k
  piv(k+1) = td(k+1) + par
END IF
k = k + 1
IF (k < n) GO TO 170
IF (piv(k) >= zero) THEN
  IF (piv(k) == zero) ksav = k
  
!     Branch if all the pivots are positive, allowing for the case when
!     G is zero.
  
  IF (ksav == 0 .AND. gsq > zero) GO TO 250
  IF (gsq == zero) THEN
    IF (par == zero) GO TO 380
    paru = par
    paruest = par
    IF (ksav == 0) GO TO 210
  END IF
  k = ksav
END IF

!     Set D to a direction of nonpositive curvature of the given tridiagonal
!     matrix, and thus revise PARLEST.

180 d(k) = one
IF (ABS(tn(k)) <= ABS(piv(k))) THEN
  dsq = one
  dhd = piv(k)
ELSE
  temp = td(k+1) + par
  IF (temp <= ABS(piv(k))) THEN
    d(k+1) = SIGN(one,-tn(k))
    dhd = piv(k) + temp - two * ABS(tn(k))
  ELSE
    d(k+1) = -tn(k) / temp
    dhd = piv(k) + tn(k) * d(k+1)
  END IF
  dsq = one + d(k+1) ** 2
END IF
190 IF (k > 1) THEN
  k = k - 1
  IF (tn(k) /= zero) THEN
    d(k) = -tn(k) * d(k+1) / piv(k)
    dsq = dsq + d(k) ** 2
    GO TO 190
  END IF
  d(1:k) = zero
END IF
parl = par
parlest = par - dhd / dsq

!     Terminate with D set to a multiple of the current D if the following
!     test suggests that it suitable to do so.

210 temp = paruest
IF (gsq == zero) temp = temp * (one-tol)
IF (paruest > zero .AND. parlest >= temp) THEN
  dtg = DOT_PRODUCT( d(1:n), gg(1:n) )
  scale = -SIGN(delta/SQRT(dsq),dtg)
  d(1:n) = scale * d(1:n)
  GO TO 380
END IF

!     Pick the value of PAR for the next iteration.

240 IF (paru == zero) THEN
  par = two * parlest + gnorm / delta
ELSE
  par = 0.5D0 * (parl+paru)
  par = MAX(par,parlest)
END IF
IF (paruest > zero) par = MIN(par,paruest)
GO TO 160

!     Calculate D for the current PAR in the positive definite case.

250 w(1) = -gg(1) / piv(1)
DO  i = 2, n
  w(i) = (-gg(i)-tn(i-1)*w(i-1)) / piv(i)
END DO
d(n) = w(n)
DO  i = nm, 1, -1
  d(i) = w(i) - tn(i) * d(i+1) / piv(i)
END DO

!     Branch if a Newton-Raphson step is acceptable.

dsq = zero
wsq = zero
DO  i = 1, n
  dsq = dsq + d(i) ** 2
  wsq = wsq + piv(i) * w(i) ** 2
END DO
IF (par /= zero .OR. dsq > delsq) THEN
  
!     Make the usual test for acceptability of a full trust region step.
  
  dnorm = SQRT(dsq)
  phi = one / dnorm - one / delta
  temp = tol * (one+par*dsq/wsq) - dsq * phi * phi
  IF (temp >= zero) THEN
    scale = delta / dnorm
    DO  i = 1, n
      d(i) = scale * d(i)
    END DO
    GO TO 380
  END IF
  IF (iterc >= 2 .AND. par <= parl) GO TO 380
  IF (paru > zero .AND. par >= paru) GO TO 380
  
!     Complete the iteration when PHI is negative.
  
  IF (phi < zero) THEN
    parlest = par
    IF (posdef == one) THEN
      IF (phi <= phil) GO TO 380
      slope = (phi-phil) / (par-parl)
      parlest = par - phi / slope
    END IF
    slope = one / gnorm
    IF (paru > zero) slope = (phiu-phi) / (paru-par)
    temp = par - phi / slope
    IF (paruest > zero) temp = MIN(temp,paruest)
    paruest = temp
    posdef = one
    parl = par
    phil = phi
    GO TO 240
  END IF
  
!     If required, calculate Z for the alternative test for convergence.
  
  IF (posdef == zero) THEN
    w(1) = one / piv(1)
    DO  i = 2, n
      temp = -tn(i-1) * w(i-1)
      w(i) = (SIGN(one,temp)+temp) / piv(i)
    END DO
    z(n) = w(n)
    DO  i = nm, 1, -1
      z(i) = w(i) - tn(i) * z(i+1) / piv(i)
    END DO
    wwsq = zero
    zsq = zero
    dtz = zero
    DO  i = 1, n
      wwsq = wwsq + piv(i) * w(i) ** 2
      zsq = zsq + z(i) ** 2
      dtz = dtz + d(i) * z(i)
    END DO
    
!     Apply the alternative test for convergence.
    
    tempa = ABS(delsq-dsq)
    tempb = SQRT(dtz*dtz+tempa*zsq)
    gam = tempa / (SIGN(tempb,dtz)+dtz)
    temp = tol * (wsq+par*delsq) - gam * gam * wwsq
    IF (temp >= zero) THEN
      DO  i = 1, n
        d(i) = d(i) + gam * z(i)
      END DO
      GO TO 380
    END IF
    parlest = MAX(parlest,par-wwsq/zsq)
  END IF
  
!     Complete the iteration when PHI is positive.
  
  slope = one / gnorm
  IF (paru > zero) THEN
    IF (phi >= phiu) GO TO 380
    slope = (phiu-phi) / (paru-par)
  END IF
  parlest = MAX(parlest,par-phi/slope)
  paruest = par
  IF (posdef == one) THEN
    slope = (phi-phil) / (par-parl)
    paruest = par - phi / slope
  END IF
  paru = par
  phiu = phi
  GO TO 240
END IF

!     Set EVALUE to the least eigenvalue of the second derivative matrix if
!     D is a Newton-Raphson step. SHFMAX will be an upper bound on EVALUE.

shfmin = zero
pivot = td(1)
shfmax = pivot
DO  k = 2, n
  pivot = td(k) - tn(k-1) ** 2 / pivot
  shfmax = MIN(shfmax,pivot)
END DO

!     Find EVALUE by a bisection method, but occasionally SHFMAX may be
!     adjusted by the rule of false position.

ksave = 0
350 shift = 0.5D0 * (shfmin+shfmax)
k = 1
temp = td(1) - shift

360 IF (temp > zero) THEN
  piv(k) = temp
  IF (k < n) THEN
    temp = td(k+1) - shift - tn(k) ** 2 / temp
    k = k + 1
    GO TO 360
  END IF
  shfmin = shift
ELSE
  IF (k < ksave) GO TO 370
  IF (k == ksave) THEN
    IF (pivksv == zero) GO TO 370
    IF (piv(k)-temp < temp-pivksv) THEN
      pivksv = temp
      shfmax = shift
    ELSE
      pivksv = zero
      shfmax = (shift*piv(k) - shfmin*temp) / (piv(k)-temp)
    END IF
  ELSE
    ksave = k
    pivksv = temp
    shfmax = shift
  END IF
END IF
IF (shfmin <= 0.99D0*shfmax) GO TO 350
370 evalue = shfmin

!     Apply the inverse Householder transformations to D.

380 nm = n - 1
DO  k = nm, 1, -1
  kp = k + 1
  sum = zero
  DO  i = kp, n
    sum = sum + d(i) * h(i,k)
  END DO
  DO  i = kp, n
    d(i) = d(i) - sum * h(i,k)
  END DO
END DO

!     Return from the subroutine.

420 RETURN
END SUBROUTINE trstep

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% lagmax.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE lagmax(n, g, h, rho, d, v, vmax)

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: g(:)
REAL (dp), INTENT(OUT)  :: h(:,:)
REAL (dp), INTENT(IN)   :: rho
REAL (dp), INTENT(OUT)  :: d(:)
REAL (dp), INTENT(OUT)  :: v(:)
REAL (dp), INTENT(OUT)  :: vmax

!     N is the number of variables of a quadratic objective function, Q say.
!     G is the gradient of Q at the origin.
!     H is the symmetric Hessian matrix of Q. Only the upper triangular and
!       diagonal parts need be set.
!     RHO is the trust region radius, and has to be positive.
!     D will be set to the calculated vector of variables.
!     The array V will be used for working space.
!     VMAX will be set to |Q(0)-Q(D)|.

!     Calculating the D that maximizes |Q(0)-Q(D)| subject to ||D|| <= RHO
!     requires of order N**3 operations, but sometimes it is adequate if
!     |Q(0)-Q(D)| is within about 0.9 of its greatest possible value. This
!     subroutine provides such a solution in only of order N**2 operations,
!     where the claim of accuracy has been tested by numerical experiments.

REAL (dp)  :: half = 0.5_dp, one = 1.0_dp, zero = 0.0_dp
REAL (dp)  :: dd, dhd, dlin, dsq, gd, gg, ghg, gnorm, halfrt, hmax, ratio
REAL (dp)  :: scale, sum, sumv, temp, tempa, tempb, tempc, tempd, tempv
REAL (dp)  :: vhg, vhv, vhw, vlin, vmu, vnorm, vsq, vv, wcos, whw, wsin, wsq
INTEGER    :: i, j, k

!     Preliminary calculations.

halfrt = SQRT(half)

!     Pick V such that ||HV|| / ||V|| is large.

hmax = zero
DO  i = 1, n
  sum = zero
  DO  j = 1, n
    h(j,i) = h(i,j)
    sum = sum + h(i,j) ** 2
  END DO
  IF (sum > hmax) THEN
    hmax = sum
    k = i
  END IF
END DO
DO  j = 1, n
  v(j) = h(k,j)
END DO

!     Set D to a vector in the subspace spanned by V and HV that maximizes
!     |(D,HD)|/(D,D), except that we set D=HV if V and HV are nearly parallel.
!     The vector that has the name D at label 60 used to be the vector W.

vsq = zero
vhv = zero
dsq = zero
DO  i = 1, n
  vsq = vsq + v(i) ** 2
  d(i) = DOT_PRODUCT( h(i,1:n), v(1:n) )
  vhv = vhv + v(i) * d(i)
  dsq = dsq + d(i) ** 2
END DO
IF (vhv*vhv <= 0.9999D0*dsq*vsq) THEN
  temp = vhv / vsq
  wsq = zero
  DO  i = 1, n
    d(i) = d(i) - temp * v(i)
    wsq = wsq + d(i) ** 2
  END DO
  whw = zero
  ratio = SQRT(wsq/vsq)
  DO  i = 1, n
    temp = DOT_PRODUCT( h(i,1:n), d(1:n) )
    whw = whw + temp * d(i)
    v(i) = ratio * v(i)
  END DO
  vhv = ratio * ratio * vhv
  vhw = ratio * wsq
  temp = half * (whw-vhv)
  temp = temp + SIGN(SQRT(temp**2+vhw**2),whw+vhv)
  DO  i = 1, n
    d(i) = vhw * v(i) + temp * d(i)
  END DO
END IF

!     We now turn our attention to the subspace spanned by G and D. A multiple
!     of the current D is returned if that choice seems to be adequate.

gg = zero
gd = zero
dd = zero
dhd = zero
DO  i = 1, n
  gg = gg + g(i) ** 2
  gd = gd + g(i) * d(i)
  dd = dd + d(i) ** 2
  sum = DOT_PRODUCT( h(i,1:n), d(1:n) )
  dhd = dhd + sum * d(i)
END DO
temp = gd / gg
vv = zero
scale = SIGN(rho/SQRT(dd),gd*dhd)
DO  i = 1, n
  v(i) = d(i) - temp * g(i)
  vv = vv + v(i) ** 2
  d(i) = scale * d(i)
END DO
gnorm = SQRT(gg)
IF (gnorm*dd <= 0.5D-2*rho*ABS(dhd) .OR. vv/dd <= 1.0D-4) THEN
  vmax = ABS(scale*(gd + half*scale*dhd))
  GO TO 170
END IF

!     G and V are now orthogonal in the subspace spanned by G and D.  Hence
!     we generate an orthonormal basis of this subspace such that (D,HV) is
!     negligible or zero, where D and V will be the basis vectors.

ghg = zero
vhg = zero
vhv = zero
DO  i = 1, n
  sum = DOT_PRODUCT( h(i,1:n), g(1:n) )
  sumv = DOT_PRODUCT( h(i,1:n), v(1:n) )
  ghg = ghg + sum * g(i)
  vhg = vhg + sumv * g(i)
  vhv = vhv + sumv * v(i)
END DO
vnorm = SQRT(vv)
ghg = ghg / gg
vhg = vhg / (vnorm*gnorm)
vhv = vhv / vv
IF (ABS(vhg) <= 0.01D0*MAX(ABS(ghg),ABS(vhv))) THEN
  vmu = ghg - vhv
  wcos = one
  wsin = zero
ELSE
  temp = half * (ghg-vhv)
  vmu = temp + SIGN(SQRT(temp**2+vhg**2),temp)
  temp = SQRT(vmu**2+vhg**2)
  wcos = vmu / temp
  wsin = vhg / temp
END IF
tempa = wcos / gnorm
tempb = wsin / vnorm
tempc = wcos / vnorm
tempd = wsin / gnorm
DO  i = 1, n
  d(i) = tempa * g(i) + tempb * v(i)
  v(i) = tempc * v(i) - tempd * g(i)
END DO

!     The final D is a multiple of the current D, V, D+V or D-V. We make the
!     choice from these possibilities that is optimal.

dlin = wcos * gnorm / rho
vlin = -wsin * gnorm / rho
tempa = ABS(dlin) + half * ABS(vmu+vhv)
tempb = ABS(vlin) + half * ABS(ghg-vmu)
tempc = halfrt * (ABS(dlin)+ABS(vlin)) + 0.25D0 * ABS(ghg+vhv)
IF (tempa >= tempb .AND. tempa >= tempc) THEN
  tempd = SIGN(rho,dlin*(vmu+vhv))
  tempv = zero
ELSE IF (tempb >= tempc) THEN
  tempd = zero
  tempv = SIGN(rho,vlin*(ghg-vmu))
ELSE
  tempd = SIGN(halfrt*rho,dlin*(ghg+vhv))
  tempv = SIGN(halfrt*rho,vlin*(ghg+vhv))
END IF
DO  i = 1, n
  d(i) = tempd * d(i) + tempv * v(i)
END DO
vmax = rho * rho * MAX(tempa,tempb,tempc)
170 RETURN
END SUBROUTINE lagmax

END MODULE Powell_Optimize

!-------------------------------------------------------------------------------
!
! Main program scriptmini
!
! reads input and starts optimisation
!
!-------------------------------------------------------------------------------
PROGRAM scriptmini

USE Powell_Optimize
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp)  :: rhobeg, rhoend
REAL(dp), DIMENSION(:), ALLOCATABLE :: x
INTEGER    :: iprint, maxfun, n,istat

OPEN(UNIT=17,FILE="scriptmini.in",FORM="FORMATTED",STATUS="OLD",IOSTAT=ISTAT)
IF (ISTAT.NE.0) THEN
   CALL write_documentation()
   STOP " Unable to open scriptmini.in "
END IF
READ(17,*,IOSTAT=ISTAT) N
IF (ISTAT.NE.0) THEN
   CALL write_documentation()
   STOP " Unable to read N in scriptmini.in "
END IF
ALLOCATE(x(N))
READ(17,*,IOSTAT=ISTAT) rhobeg,rhoend
IF (ISTAT.NE.0) THEN
   CALL write_documentation()
   STOP " Unable to read rhobeg,rhoend in scriptmini.in "
END IF
READ(17,*,IOSTAT=ISTAT) maxfun
IF (ISTAT.NE.0) THEN
   CALL write_documentation()
   STOP " Unable to read maxfun in scriptmini.in "
END IF
READ(17,*,IOSTAT=ISTAT) iprint
IF (ISTAT.NE.0) THEN
   CALL write_documentation()
   STOP " Unable to read iprint in scriptmini.in "
END IF
READ(17,*,IOSTAT=ISTAT) x
IF (ISTAT.NE.0) THEN
   CALL write_documentation()
   STOP " Unable to read x in scriptmini.in "
END IF

CALL uobyqa (n, x, rhobeg, rhoend, iprint, maxfun)

DEALLOCATE(x)

END PROGRAM scriptmini

!-------------------------------------------------------------------------------
!
! calfun: the actual evaluation of the script
!
!
!-------------------------------------------------------------------------------
SUBROUTINE calfun(n, x, f)

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: f

character(LEN=40) :: format

! write variables on a single line in the file
WRITE(format,'(A1,I4.4,A12)') '(',n,'(1X,F20.10))'
OPEN(UNIT=17,FILE="scriptmini_eval.in")
WRITE(UNIT=17,FMT=format) x(1:n)
CLOSE(UNIT=17)

! execute scriptmini_eval
CALL system("./scriptmini_eval")

! read value of the energy back
OPEN(UNIT=17,FILE="scriptmini_eval.out")
READ(UNIT=17,FMT=*) f
CLOSE(UNIT=17)

END SUBROUTINE calfun
