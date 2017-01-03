!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Simple splines
!> Splines are fully specified by the interpolation points, except that
!> at the ends, we have the freedom to prescribe the second derivatives.
!> If we know a derivative at an end (exactly), then best is to impose that.
!> Otherwise, it is better to use the "consistent" end conditions: the second
!> derivative is determined such that it is smooth.
!>
!> High level API: spline3, spline3ders.
!> Low level API: the rest of public soubroutines.
!>
!> Use the high level API to obtain cubic spline fit with consistent boundary
!> conditions and optionally the derivatives. Use the low level API if more fine
!> grained control is needed.
!>
!> This module is based on a code written by John E. Pask, LLNL.
!> \par History
!>      Adapted for use in CP2K  (30.12.2016,JGH)
!> \author JGH
! **************************************************************************************************
MODULE splines

   USE kinds,                           ONLY: dp
   USE lapack,                          ONLY: lapack_sgbsv
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'splines'

   PUBLIC :: spline3, spline3ders

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \param y ...
!> \param xnew ...
!> \retval ynew ...
! **************************************************************************************************
   FUNCTION spline3(x, y, xnew) RESULT(ynew)
      ! Takes the function values 'y' on the grid 'x' and returns new values 'ynew'
      ! at the given grid 'xnew' using cubic splines interpolation with such
      ! boundary conditions so that the 2nd derivative is consistent with the
      ! interpolating cubic.
      REAL(dp), INTENT(in)                               :: x(:), y(:), xnew(:)
      REAL(dp)                                           :: ynew(SIZE(xnew))

      INTEGER                                            :: i, ip
      REAL(dp)                                           :: c(0:4, SIZE(x)-1)

      ! get spline parameters: 2nd derivs at ends determined by cubic interpolation
      CALL spline3pars(x, y, [2, 2], [0._dp, 0._dp], c)

      ip = 0
      DO i = 1, SIZE(xnew)
         ip = iixmin(xnew(i), x, ip)
         ynew(i) = poly3(xnew(i), c(:, ip))
      END DO
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \param y ...
!> \param xnew ...
!> \param ynew ...
!> \param dynew ...
!> \param d2ynew ...
! **************************************************************************************************
   SUBROUTINE spline3ders(x, y, xnew, ynew, dynew, d2ynew)
      ! Just like 'spline', but also calculate 1st and 2nd derivatives
      REAL(dp), INTENT(in)                               :: x(:), y(:), xnew(:)
      REAL(dp), INTENT(out), OPTIONAL                    :: ynew(:), dynew(:), d2ynew(:)

      INTEGER                                            :: i, ip
      REAL(dp)                                           :: c(0:4, SIZE(x)-1)

      CALL spline3pars(x, y, [2, 2], [0._dp, 0._dp], c)

      ip = 0
      DO i = 1, SIZE(xnew)
         ip = iixmin(xnew(i), x, ip)
         IF (PRESENT(ynew)) ynew(i) = poly3(xnew(i), c(:, ip))
         IF (PRESENT(dynew)) dynew(i) = dpoly3(xnew(i), c(:, ip))
         IF (PRESENT(d2ynew)) d2ynew(i) = d2poly3(xnew(i), c(:, ip))
      END DO
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param xi ...
!> \param yi ...
!> \param bctype ...
!> \param bcval ...
!> \param c ...
! **************************************************************************************************
   SUBROUTINE spline3pars(xi, yi, bctype, bcval, c)
      ! Returns parameters c defining cubic spline interpolating x-y data xi, yi, with
      ! boundary conditions specified by bcytpe, bcvals
      REAL(dp), INTENT(in)                               :: xi(:), yi(:)
      INTEGER, INTENT(in)                                :: bctype(2)
      REAL(dp), INTENT(in)                               :: bcval(2)
      REAL(dp), INTENT(out)                              :: c(0:, :)

      INTEGER                                            :: i, i2, info, ipiv(4), j, n
      REAL(dp) :: Ae(4, 4), be(4), bemat(4, 1), c1, c2, c3, c4, ce(4), d2p1, d2pn, x0, xe(4), &
         ye(4), hi(SIZE(c, 2)), cs(2*SIZE(c, 2)), bs(2*SIZE(c, 2)), bmat(2*SIZE(c, 2), 1), &
         As(5, 2*SIZE(c, 2))
      INTEGER                                            :: ipiv2(2*SIZE(c, 2))

      ! x values of data
      ! y values of data
      ! type of boundary condition at each end:
      ! bctype(1) = type at left end, bctype(2) = type at right end.
      ! 1 = specified 2nd derivative, 2 = 2nd derivative consistent with interpolating cubic.
      ! boundary condition values at each end:
      ! bcval(1) = value at left end, bcval(2) = value at right end
      ! parameters defining spline: c(i,j) = ith parameter of jth
      ! spline polynomial, p_j = sum_{i=1}^4 c(i,j) (x-c(0,j))^(i-1), j = 1..n-1, n = # of data pts.
      ! dimensions: c(0:4,1:n-1)
      ! spline eq. matrix -- LAPACK band form
      ! spline eq. rhs vector
      ! spline eq. solution vector
      ! spline intervals
      ! end-cubic eq. matrix
      ! end-cubic eq. rhs vector
      ! end-cubic eq. solution vector
      ! x,y values at ends
      ! 2nd derivatives at ends
      ! expansion center
      ! expansion coefficients
      ! number of data points
      ! lapack variables

      ! check input parameters
      IF (bctype(1) < 1 .OR. bctype(1) > 2) CALL stop_error("spline3pars error: bctype /= 1 or 2.")
      IF (bctype(2) < 1 .OR. bctype(2) > 2) CALL stop_error("spline3pars error: bctype /= 1 or 2.")
      IF (SIZE(c, 1) /= 5) CALL stop_error("spline3pars error: size(c,1) /= 5.")
      IF (SIZE(c, 2) /= SIZE(xi)-1) CALL stop_error("spline3pars error: size(c,2) /= size(xi)-1.")
      IF (SIZE(xi) /= SIZE(yi)) CALL stop_error("spline3pars error: size(xi) /= size(yi)")

      ! To get rid of compiler warnings:
      d2p1 = 0
      d2pn = 0

      ! initializations
      n = SIZE(xi)
      DO i = 1, n-1
         hi(i) = xi(i+1)-xi(i)
      END DO

      ! compute interpolating-cubic 2nd derivs at ends, if required
      ! left end
      IF (bctype(1) == 2) THEN
         IF (n < 4) CALL stop_error("spline3pars error: n < 4")
         xe = xi(1:4)
         ye = yi(1:4)
         x0 = xe(1) ! center at end
         DO i = 1, 4
            DO j = 1, 4
               Ae(i, j) = (xe(i)-x0)**(j-1)
            END DO
         END DO
         Ae(:, 1) = 1 ! set 0^0 = 1
         be = ye; bemat(:, 1) = be
         CALL dgesv(4, 1, Ae, 4, ipiv, bemat, 4, info)
         IF (info /= 0) CALL stop_error("spline3pars error: dgesv error.")
         ce = bemat(:, 1)
         d2p1 = 2*ce(3)
      END IF
      ! right end
      IF (bctype(2) == 2) THEN
         IF (n < 4) CALL stop_error("spline3pars error: n < 4")
         xe = xi(n-3:n)
         ye = yi(n-3:n)
         x0 = xe(4) ! center at end
         DO i = 1, 4
            DO j = 1, 4
               Ae(i, j) = (xe(i)-x0)**(j-1)
            END DO
         END DO
         Ae(:, 1) = 1 ! set 0^0 = 1
         be = ye; bemat(:, 1) = be
         CALL dgesv(4, 1, Ae, 4, ipiv, bemat, 4, info)
         IF (info /= 0) CALL stop_error("spline3pars error: dgesv error.")
         ce = bemat(:, 1)
         d2pn = 2*ce(3)
      END IF

      ! set 2nd derivs at ends
      IF (bctype(1) == 1) d2p1 = bcval(1)
      IF (bctype(2) == 1) d2pn = bcval(2)

      ! construct spline equations -- LAPACK band form
      ! basis: phi1 = -(x-x_i)/h_i, phi2 = (x-x_{i+1})/h_i, phi3 = phi1^3-phi1, phi4 = phi2^3-phi2
      ! on interval [x_i,x_{i+1}] of length h_i = x_{i+1}-x_i
      !A=0  ! full matrix
      As = 0
      ! left end condition
      As(4, 1) = 6/hi(1)**2
      bs(1) = d2p1
      ! internal knot conditions
      DO i = 2, n-1
         i2 = 2*(i-1)
         As(5, i2-1) = 1/hi(i-1)
         As(4, i2) = 2/hi(i-1)
         As(3, i2+1) = 2/hi(i)
         As(2, i2+2) = 1/hi(i)
         As(5, i2) = 1/hi(i-1)**2
         As(4, i2+1) = -1/hi(i)**2
         bs(i2) = (yi(i+1)-yi(i))/hi(i)-(yi(i)-yi(i-1))/hi(i-1)
         bs(i2+1) = 0
      END DO
      ! right end condition
      As(4, 2*(n-1)) = 6/hi(n-1)**2
      bs(2*(n-1)) = d2pn

      ! solve spline equations -- LAPACK band form
      bmat(:, 1) = bs
      CALL lapack_sgbsv(2*(n-1), 1, 2, 1, As, 5, ipiv2, bmat, 2*(n-1), info)
      IF (info /= 0) CALL stop_error("spline3pars error: dgbsv error.")
      cs = bmat(:, 1)

      ! transform to (x-x0)^(i-1) basis and return
      DO i = 1, n-1
         ! coefficients in spline basis:
         c1 = yi(i)
         c2 = yi(i+1)
         c3 = cs(2*i-1)
         c4 = cs(2*i)
         ! coefficients in (x-x0)^(i-1) basis
         c(0, i) = xi(i)
         c(1, i) = c1
         c(2, i) = -(c1-c2+2*c3+c4)/hi(i)
         c(3, i) = 3*c3/hi(i)**2
         c(4, i) = (-c3+c4)/hi(i)**3
      END DO
   END SUBROUTINE

!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \param xi ...
!> \param c ...
!> \param val ...
!> \param der ...
! **************************************************************************************************
   SUBROUTINE spline3valder(x, xi, c, val, der)
      ! Returns value and 1st derivative of spline defined by knots xi and parameters c
      ! returned by spline3pars
      REAL(dp), INTENT(in)                               :: x, xi(:), c(0:, :)
      REAL(dp), INTENT(out)                              :: val, der

      INTEGER                                            :: i1, n

      ! point at which to evaluate spline
      ! spline knots (x values of data)
      ! spline parameters: c(i,j) = ith parameter of jth
      ! spline polynomial, p_j = sum_{i=1}^4 c(i,j) (x-c(0,j))^(i-1), j = 1..n-1, n = # of data pts.
      ! dimensions: c(0:4,1:n-1)
      ! value of spline at x
      ! 1st derivative of spline at x
      ! number of knots

      ! initialize, check input parameters
      n = SIZE(xi)
      IF (SIZE(c, 1) /= 5) CALL stop_error("spline3 error: size(c,1) /= 5.")
      IF (SIZE(c, 2) /= SIZE(xi)-1) CALL stop_error("spline3 error: size(c,2) /= size(xi)-1.")
      ! find interval containing x
      i1 = iix(x, xi)
      ! return value and derivative
      val = poly3(x, c(:, i1))
      der = dpoly3(x, c(:, i1))
   END SUBROUTINE

!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \param xi ...
!> \retval i1 ...
! **************************************************************************************************
   INTEGER FUNCTION iix(x, xi) RESULT(i1)
      ! Returns index i of interval [xi(i),xi(i+1)] containing x in mesh xi,
      ! with intervals indexed by left-most points.
      ! N.B.: x outside [x1,xn] are indexed to nearest end.
      ! Uses bisection, except if "x" lies in the first or second elements (which is
      ! often the case)
      REAL(dp), INTENT(in) :: x ! target value
      REAL(dp), INTENT(in) :: xi(:) ! mesh, xi(i) < xi(i+1)
      INTEGER n ! number of mesh points
      INTEGER i2, ic

      n = SIZE(xi)
      i1 = 1
      IF (n < 2) THEN
         CALL stop_error("error in iix: n < 2")
      ELSEIF (n == 2) THEN
         i1 = 1
      ELSEIF (n == 3) THEN
         IF (x <= xi(2)) THEN ! first element
            i1 = 1
         ELSE
            i1 = 2
         END IF
      ELSEIF (x <= xi(1)) THEN ! left end
         i1 = 1
      ELSEIF (x <= xi(2)) THEN ! first element
         i1 = 1
      ELSEIF (x <= xi(3)) THEN ! second element
         i1 = 2
      ELSEIF (x >= xi(n)) THEN ! right end
         i1 = n-1
      ELSE
         ! bisection: xi(i1) <= x < xi(i2)
         i1 = 3; i2 = n
         DO
            IF (i2-i1 == 1) EXIT
            ic = i1+(i2-i1)/2
            IF (x >= xi(ic)) THEN
               i1 = ic
            ELSE
               i2 = ic
            ENDIF
         END DO
      END IF
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \param xi ...
!> \param i_min ...
!> \retval ip ...
! **************************************************************************************************
   INTEGER FUNCTION iixmin(x, xi, i_min) RESULT(ip)
      ! Just like iix, but assumes that x >= xi(i_min)
      REAL(dp), INTENT(in) :: x, xi(:)
      INTEGER, INTENT(in) :: i_min
      IF (i_min >= 1 .AND. i_min <= SIZE(xi)-1) THEN
         ip = iix(x, xi(i_min:))+i_min-1
      ELSE
         ip = iix(x, xi)
      END IF
   END FUNCTION

!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \param n ...
!> \param x1 ...
!> \param xn ...
!> \retval iixun ...
! **************************************************************************************************
   FUNCTION iixun(x, n, x1, xn)
      ! Returns index i of interval [x(i),x(i+1)] containing x in uniform mesh defined by
      !   x(i) = x1 + (i-1)/(n-1)*(xn-x1), i = 1 .. n,
      ! with intervals indexed by left-most points.
      ! N.B.: x outside [x1,xn] are indexed to nearest end.
      REAL(dp), INTENT(in)                               :: x
      INTEGER, INTENT(in)                                :: n
      REAL(dp), INTENT(in)                               :: x1, xn
      INTEGER                                            :: iixun

      INTEGER                                            :: i

      ! index i of interval [x(i),x(i+1)] containing x
      ! target value
      ! number of mesh points
      ! initial point of mesh
      ! final point of mesh

      ! compute index
      i = INT((x-x1)/(xn-x1)*(n-1))+1
      ! reset if ouside 1..n
      IF (i < 1) i = 1
      IF (i > n-1) i = n-1
      iixun = i
   END FUNCTION

!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \param n ...
!> \param x1 ...
!> \param alpha ...
!> \param beta ...
!> \retval iixexp ...
! **************************************************************************************************
   FUNCTION iixexp(x, n, x1, alpha, beta)
      ! Returns index i of interval [x(i),x(i+1)] containing x in exponential mesh defined by
      !   x(i) = x1 + alpha [ exp(beta(i-1)) - 1 ], i = 1 .. n,
      ! where alpha = (x(n) - x(1))/[ exp(beta(n-1)) - 1 ],
      ! beta = log(r)/(n-2), r = (x(n)-x(n-1))/(x(2)-x(1)) = ratio of last to first interval,
      ! and intervals indexed by left-most points.
      ! N.B.: x outside [x1,xn] are indexed to nearest end.
      REAL(dp), INTENT(in)                               :: x
      INTEGER, INTENT(in)                                :: n
      REAL(dp), INTENT(in)                               :: x1, alpha, beta
      INTEGER                                            :: iixexp

      INTEGER                                            :: i

      ! index i of interval [x(i),x(i+1)] containing x
      ! target value
      ! number of mesh points
      ! initial point of mesh
      ! mesh parameter:
      !   x(i) = x1 + alpha [ exp(beta(i-1)) - 1 ], i = 1 .. n,
      ! where alpha = (x(n) - x(1))/[ exp(beta(n-1)) - 1 ],
      ! beta = log(r)/(n-2), r = (x(n)-x(n-1))/(x(2)-x(1)) = ratio of last to first interval,
      ! mesh parameter

      ! compute index
      i = INT(LOG((x-x1)/alpha+1)/beta)+1
      ! reset if outside 1..n
      IF (i < 1) i = 1
      IF (i > n-1) i = n-1
      iixexp = i
   END FUNCTION

!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \param c ...
!> \retval poly3 ...
! **************************************************************************************************
   FUNCTION poly3(x, c)
      ! returns value of polynomial \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
      REAL(dp), INTENT(in)                               :: x, c(0:)
      REAL(dp)                                           :: poly3

      REAL(dp)                                           :: dx

      ! point at which to evaluate polynomial
      ! coefficients: poly = \sum_{i=1}^4 c(i) (x-c(0))^(i-1)

      dx = x-c(0)
      poly3 = c(1)+c(2)*dx+c(3)*dx**2+c(4)*dx**3
   END FUNCTION

!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \param c ...
!> \retval dpoly3 ...
! **************************************************************************************************
   FUNCTION dpoly3(x, c)
      ! returns 1st derivative of polynomial \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
      REAL(dp), INTENT(in)                               :: x, c(0:)
      REAL(dp)                                           :: dpoly3

      REAL(dp)                                           :: dx

      ! point at which to evaluate polynomial
      ! coefficients: poly = \sum_{i=1}^4 c(i) (x-c(0))^(i-1)

      dx = x-c(0)
      dpoly3 = c(2)+2*c(3)*dx+3*c(4)*dx**2
   END FUNCTION

!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \param c ...
!> \retval d2poly3 ...
! **************************************************************************************************
   FUNCTION d2poly3(x, c)
      ! returns 2nd derivative of polynomial \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
      REAL(dp), INTENT(in)                               :: x, c(0:)
      REAL(dp)                                           :: d2poly3

      REAL(dp)                                           :: dx

      ! point at which to evaluate polynomial
      ! coefficients: poly = \sum_{i=1}^4 c(i) (x-c(0))^(i-1)

      dx = x-c(0)
      d2poly3 = 2*c(3)+6*c(4)*dx
   END FUNCTION

!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief ...
!> \param msg ...
! **************************************************************************************************
   SUBROUTINE stop_error(msg)
      ! Aborts the program
      CHARACTER(LEN=*)                                   :: msg

! Message to print on stdout
      CPABORT(msg)
   END SUBROUTINE

END MODULE splines
