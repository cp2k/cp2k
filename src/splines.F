!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE splines
  
! routines for handling splines
!   init_splinexy: allocates x and y vectors for splines
!   init_spline:   generate table for spline (allocate spl%y2)
!   spline:        return value of spline for given abscissa (optional:also y1)
!   spline_1:      return value of 1. derivative of spline for given abscissa
!   spline_int:    return value of integral on given interval of spline
!   kill_spline:   destructor ( spl%x,y und/oder spl%y2)
  
! NB: splines are always ``natural splines'', i.e. values of first
!     derivative at the end-points cannot be specified
!-----------------------------------------------------------------------------!
  
  USE kinds, ONLY : dbl
  USE stop_program, ONLY : stop_prg, stop_memory
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: spline_data, init_splinexy, init_spline, spline, spline_1, &
       spline_int, kill_spline, splineh
  
  TYPE spline_data
     REAL ( dbl ), POINTER :: x ( : ), y ( : )
     REAL ( dbl ), POINTER :: y2 ( : )
     INTEGER :: n, pos
     REAL ( dbl ) :: h, invh, h26, h16
  END TYPE spline_data
  
CONTAINS

!******************************************************************************

SUBROUTINE init_splinexy ( spl, nn )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( spline_data ), INTENT ( INOUT ) :: spl
  INTEGER, INTENT ( IN ) :: nn
  
! Locals
  INTEGER :: err
  
!------------------------------------------------------------------------------
  
  spl % n = nn
  IF ( ASSOCIATED ( spl % x ) ) THEN
     DEALLOCATE ( spl % x, STAT = err )
     IF ( err /= 0 ) CALL stop_memory ( 'init_splinexy', 'spl%x')
     NULLIFY ( spl % x )
  END IF
  IF ( ASSOCIATED ( spl % y ) ) THEN
     DEALLOCATE ( spl % y, STAT = err )
     IF ( err /= 0 ) CALL stop_memory ( 'init_splinexy', 'spl%y')
     NULLIFY ( spl % y )
  END IF
  
  ALLOCATE ( spl % x ( 1:nn ), STAT = err )
  IF ( err /= 0 ) CALL stop_memory ( 'init_splinexy', 'spl%x', nn )
  
  ALLOCATE ( spl%y(1:nn), STAT = err )
  IF ( err /= 0 ) CALL stop_memory ( 'init_splinexy', 'spl%y', nn)
  
END SUBROUTINE init_splinexy

!******************************************************************************

!   endpt: 's': regular spacing
!          'l': left; 'r': right; 'b': both = specify 1-deriv for
!                                                 which endpoints

SUBROUTINE init_spline ( spl, endpt, y1a, y1b )
  
  IMPLICIT NONE
  
! Arguments
  TYPE (spline_data), INTENT ( INOUT ) :: spl
  CHARACTER ( LEN = * ), INTENT ( IN ), OPTIONAL :: endpt
  REAL ( dbl ), INTENT ( IN ), OPTIONAL :: y1a, y1b
  
! Locals
  INTEGER :: err, i, k
  REAL ( dbl ) :: p, qn, sig, un, y1l, y1r
  REAL ( dbl ), POINTER :: ww ( : )
  CHARACTER ( LEN = 8 ) :: ep
  LOGICAL :: reg, lep, rep
  
!------------------------------------------------------------------------------
  
  IF ( .NOT. PRESENT ( endpt ) ) THEN
     ep = ''
  ELSE
     ep = endpt
  END IF
  
  reg = (scan(ep,'sS')>0)
  lep = (scan(ep,'lL')>0) .OR. (scan(ep,'bB')>0)
  rep = (scan(ep,'rR')>0) .OR. (scan(ep,'bB')>0)
  IF (lep) y1l = y1a
  IF (rep .AND. .NOT. lep) y1r = y1a
  IF (rep .AND. lep) y1r = y1b
  IF ((lep .OR. rep) .AND. .NOT. PRESENT ( y1a)) &
       CALL stop_prg ( 'init_spline','first deriv. at end-point missing')
  IF (lep .AND. rep .AND. .NOT. PRESENT ( y1b)) CALL stop_prg( &
       'init_spline','first deriv. at end-point missing')
  
  spl%pos = 1
  IF ( ASSOCIATED ( spl%y2)) THEN
     DEALLOCATE (spl%y2,STAT=err)
     IF (err /= 0 ) CALL stop_memory ( 'init_spline','spl%y2')
     NULLIFY (spl%y2)
  END IF
  ALLOCATE (spl%y2(1:spl%n),STAT=err)
  IF ( err /= 0 ) CALL stop_memory ( 'init_spline', 'spl%y2', spl % n )
  ALLOCATE (ww(1:spl%n),STAT=err)
  IF ( err /= 0 ) CALL stop_memory ( 'init_spline', 'ww', spl % n )
  
  IF (lep) THEN
     spl%y2(1) = -0.5_dbl
     ww(1) = ( 3.0_dbl &
          /(spl%x(2)-spl%x(1)))*((spl%y(2)-spl%y(1))/(spl%x(2)-spl% &
          x(1))-y1l)
  ELSE
     spl%y2(1) = 0
     ww(1) = 0.0_dbl
  END IF
  
  DO i = 2, spl%n - 1
     sig = (spl%x(i)-spl%x(i-1))/(spl%x(i+1)-spl%x(i-1))
     p = sig*spl%y2(i-1) + 2.0_dbl
     spl%y2(i) = (sig-1.)/p
     ww(i) = (6.0_dbl*((spl%y(i+1)-spl%y(i))/(spl%x(i+1)-spl%x(i))-(spl%y(i)- &
          spl%y(i-1))/(spl%x(i)-spl%x(i-1)))/(spl%x(i+ &
          1)-spl%x(i-1))-sig*ww(i-1))/p
  END DO
  
  IF ( rep ) THEN
     qn = 0.5_dbl
     un = ( 3.0_dbl &
          /(spl%x(spl%n)-spl%x(spl%n-1)))*(y1r-(spl%y(spl%n)-spl%y(spl &
          %n-1))/(spl%x(spl%n)-spl%x(spl%n-1)))
  ELSE
     qn = 0
     un = 0
  END IF
  
  spl%y2(spl%n) = (un-qn*ww(spl%n-1))/(qn*spl%y2(spl%n-1)+1.)
  DO k = spl%n - 1, 1, -1
     spl%y2(k) = spl%y2(k)*spl%y2(k+1) + ww(k)
  END DO
  
  DEALLOCATE (ww,STAT=err)
  IF (err /= 0 ) CALL stop_memory ( 'init_spline','ww')
  
  IF (reg) THEN
     spl%h = (spl%x(spl%n)-spl%x(1))/(spl%n-1.)
     spl%h16 = spl%h/6.
     spl%h26 = spl%h**2/6.
     spl%invh = 1/spl%h
  ELSE
     spl%h = 0
     spl%invh = 0
  END IF
  
END SUBROUTINE init_spline

!******************************************************************************

FUNCTION interv ( spl, xx )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( spline_data ), INTENT ( IN ) :: spl
  REAL ( dbl ), INTENT ( IN ) :: xx
  
! Locals
  INTEGER :: interv
  INTEGER :: khi, klo, i, p, n, k
  
!------------------------------------------------------------------------------
  
  IF (spl%h /= 0 ) THEN
     i = (xx-spl%x(1))*spl%invh + 1
     IF (i<1 .OR. i>spl%n) CALL stop_prg ( 'interv', &
          'illegal x-value passed to interv')
     interv = i
     RETURN
  END IF
  
  p = spl%pos
  IF (p>=spl%n) p = spl%n/2
  i = 0
  n = spl%n
  
! check if interval is close to previous interval
  IF ( xx < spl % x ( p + 1 ) ) THEN
     IF ( xx >= spl % x ( p ) ) THEN
        i = spl%pos
     ELSE IF ( p > 1 .AND. xx >= spl % x ( p - 1 ) ) THEN
        i = p - 1
     ELSE
        klo = 1
        khi = p + 1
     END IF
  ELSE IF ( p + 2 <= n .AND. xx < spl % x ( p + 2 ) ) THEN
     i = p + 1
  ELSE
     klo = p + 1
     khi = n
  END IF
  
  IF ( i == 0 ) THEN ! perform binary search
     IF ( xx < spl % x ( 1 ) .OR. xx > spl % x ( n ) ) &
          CALL stop_prg ( 'interv', 'xx value out of spline-range' )
     DO WHILE ( khi - klo > 1 )
        k = ( khi + klo ) / 2
        IF ( spl % x ( k ) > xx ) THEN
           khi = k
        ELSE
           klo = k
        END IF
     END DO
     i = klo
  END IF
  
  interv = i
  
END FUNCTION interv

!******************************************************************************

FUNCTION spline_slow ( spl, xx, y1 )
  
  IMPLICIT NONE
  
! Arguments
  TYPE (spline_data), INTENT ( INOUT ) :: spl
  REAL ( dbl ), INTENT ( IN ) :: xx
  REAL ( dbl ), INTENT ( OUT ), OPTIONAL :: y1
  
! Locals
  REAL ( dbl ) :: spline_slow
  INTEGER :: khi, klo
  REAL ( dbl ) :: a, b, h, invh
  
!------------------------------------------------------------------------------
  
  spl%pos = interv(spl,xx)
  klo = spl%pos
  khi = spl%pos + 1
  
  IF (spl%h /= 0 ) THEN
     h = spl%h
     invh = spl%invh
  ELSE
     h = spl%x(khi) - spl%x(klo)
     invh = spl%invh
     IF ( h == 0.0_dbl ) CALL stop_prg ( 'spline', 'bad spl%x input' )
  END IF
  
  a = (spl%x(khi)-xx)*invh
  b = 1 - a
  spline_slow = a*spl%y(klo) + b*spl%y(khi) + ((a**3-a)*spl%y2(klo)+(b** &
       3-b)*spl%y2(khi))*(h**2)/6.0_dbl
  IF ( PRESENT ( y1 ) ) y1 = (spl%y(khi)-spl%y(klo))*invh + &
       ((1.0_dbl-3*a**2)*spl%y2(klo)+(3*b**2-1.0_dbl)*spl%y2(khi))*h/6.0_dbl
  
END FUNCTION spline_slow

!******************************************************************************

FUNCTION spline ( spl, xx, y1 )
  
  IMPLICIT NONE
  
! Return value
  REAL ( dbl ) :: spline
  
! Arguments
  TYPE ( spline_data ), INTENT ( INOUT ) :: spl
  REAL ( dbl ), INTENT ( IN ) :: xx
  REAL ( dbl ), INTENT ( OUT ), OPTIONAL :: y1
  
! Locals
  INTEGER :: khi, klo, i
  REAL ( dbl ) :: a, b, h, invh, t, ylo, yhi, y2lo, y2hi, d, d0
  
!------------------------------------------------------------------------------
  
  h = spl % h
  invh = spl%invh
  IF (h /= 0 ) THEN
!     d=xx-spl%x(1); i=INT(d*spl%invh); d0=REAL(i)*h; i=i+1
     i = (xx-spl%x(1))*invh + 1
     IF (i<1 .OR. i>spl%n) CALL stop_prg ( 'spline','illegal x-value')
  ELSE
     i = interv(spl,xx)
     spl%pos = 1
  END IF
  
  a = (spl%x(i+1)-xx)*invh
  b = 1 - a
  t = -a*b
!  b=(d-d0)*invh; a=1-b; t=-a*b
  ylo = spl%y(i)
  yhi = spl%y(i+1)
  y2lo = spl%y2(i)
  y2hi = spl%y2(i+1)
  spline = a*ylo + b*yhi + ((a+1)*y2lo+(b+1)*y2hi)*t*spl%h26
  IF ( PRESENT ( y1)) y1 = (yhi-ylo)*invh + ((1.-3*a*a)*y2lo+(3*b*b-1.)* &
       y2hi)*spl%h16
  
END FUNCTION spline

!******************************************************************************

FUNCTION splineh ( spl, xx, y1 )
  
  IMPLICIT NONE
  
! Return value
  REAL ( dbl ) :: splineh

! Arguments
  TYPE (spline_data), INTENT ( IN ) :: spl
  REAL ( dbl ), INTENT ( IN ) :: xx
  REAL ( dbl ), INTENT ( OUT ) :: y1
  
! Locals
  INTEGER :: khi, klo, i
  REAL ( dbl ) :: a, b, h, invh, t, ylo, yhi, y2lo, y2hi, d, d0
  
!------------------------------------------------------------------------------
  
! fast spline for pair potentials without checks
  h = spl%h
  invh = spl%invh
  d=xx-spl%x(1); i=INT(d*spl%invh); d0=REAL(i)*h; i=i+1
  i = (xx-spl%x(1))*invh + 1
  
  a = (spl%x(i+1)-xx)*invh
  b = 1 - a
  t = -a*b
!    b=(d-d0)*invh; a=1-b; t=-a*b
  ylo = spl%y(i)
  yhi = spl%y(i+1)
  y2lo = spl%y2(i)
  y2hi = spl%y2(i+1)
  splineh = a*ylo + b*yhi + ((a+1)*y2lo+(b+1)*y2hi)*t*spl%h26
  y1 = (yhi-ylo)*invh + ((1.-3*a*a)*y2lo+(3*b*b-1.)*y2hi)*spl%h16
  
END FUNCTION splineh

!******************************************************************************

FUNCTION spline_1 ( spl, xx )
  
  IMPLICIT NONE
  
! Return value
  REAL ( dbl ) :: spline_1
  
! Arguments
  TYPE ( spline_data ), INTENT ( INOUT ) :: spl
  REAL ( dbl ), INTENT ( IN ) :: xx
  
! Locals
  INTEGER :: khi, klo
  REAL ( dbl ) :: a, b, h
  
!------------------------------------------------------------------------------
  
  spl % pos = interv ( spl, xx )
  klo = spl % pos
  khi = spl % pos + 1
  
  h = spl % x ( khi ) - spl % x ( klo )
  IF ( h == 0.0_dbl ) CALL stop_prg ( 'spline', 'bad spl%x input' )
  a = ( spl % x ( khi ) - xx ) / h
  b = 1.0_dbl - a
  spline_1 = ( spl%y(khi)-spl%y(klo) ) / h &
       + ( ( 1.0_dbl - 3.0_dbl * a ** 2 ) * spl % y2(klo) &
       + ( 3.0_dbl * b ** 2 - 1.0_dbl ) * spl % y2(khi) ) * h / 6.0_dbl
  
END FUNCTION spline_1

!******************************************************************************

FUNCTION stamm ( spl, p, x )
  
  IMPLICIT NONE
  
! Return value
  REAL ( dbl ) :: stamm
  
! Arguments
  TYPE (spline_data), INTENT ( IN ) :: spl
  INTEGER, INTENT ( IN ) :: p
  REAL ( dbl ), INTENT ( IN ) :: x
  
! Locals
  REAL ( dbl ) :: a, b, aa, bb, h
  
!------------------------------------------------------------------------------
  
  h = spl%x(p+1) - spl%x(p)
  b = (x-spl%x(p))/h
  a = 1 - b
  aa = a**2
  bb = b**2
  stamm = 0.5*h*(bb*spl%y(p+1)-aa*spl%y(p)) + h**3/12.*(aa*(1-0.5*aa)* &
       spl%y2(p)-bb*(1-0.5*bb)*spl%y2(p+1))
  
END FUNCTION stamm

!******************************************************************************

FUNCTION spline_int(spl,x0,x1)
  
  IMPLICIT NONE
  
! Return value
  REAL ( dbl ) :: spline_int
  
! Arguments
  TYPE (spline_data), INTENT ( INOUT ) :: spl
  REAL ( dbl ), INTENT ( IN ) :: x0, x1
  
! Local
  INTEGER :: j, pa, pb
  REAL ( dbl ) :: h, vorz, xa, xb, i1, i2
  
!------------------------------------------------------------------------------
  
  vorz = 1.0_dbl
  xa = MIN ( x0, x1 )
  xb = MAX ( x0, x1 )
  IF ( x0 > x1 ) vorz = -1.0_dbl
  IF ( xa < spl % x ( 1 ) .OR. xb > spl % x ( spl % n ) ) &
       CALL stop_prg ( 'spline_int', 'illegal integration range' )
  
  pa = interv ( spl, xa )
  pb = interv ( spl, xb )
  
  IF ( pa == pb ) THEN
     spline_int = vorz * ( stamm(spl,pa,xb) - stamm(spl,pa,xa) )
     RETURN
  END IF
  
  i1 = 0
  i2 = 0
  DO j = pa + 1, pb - 1
     h = spl%x(j+1) - spl%x(j)
     i1 = i1 + h*(spl%y(j)+spl%y(j+1))
     i2 = i2 + h**3*(spl%y2(j)+spl%y2(j+1))
  END DO
  h = spl%x(pa+1) - spl%x(pa)
  i1 = i1 + h*spl%y(pa+1)
  i2 = i2 + h**3*spl%y2(pa+1)
  h = spl%x(pb+1) - spl%x(pb)
  i1 = i1 + h*spl%y(pb)
  i2 = i2 + h**3*spl%y2(pb)
  
  spline_int = vorz*(i1/2.-i2/24.+stamm(spl,pb,xb)-stamm(spl,pa,xa))
  
END FUNCTION spline_int

!******************************************************************************

!   deallocate splines
!   what=='a' or not present: deallocate all (spl%x, spl%y, spl%y2)
!   what=='d': deallocate only data (spl%x, spl%y)
!   what=='2': deallocate only table of 2. derivatives (spl%y2)

SUBROUTINE kill_spline ( spl, what )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( spline_data ), INTENT ( INOUT ) :: spl
  CHARACTER, INTENT ( IN ), OPTIONAL :: what
  
! Locals
  CHARACTER :: w
  INTEGER :: err
  
!------------------------------------------------------------------------------
  
  IF ( PRESENT ( what ) ) THEN
     w = what
  ELSE
     w = 'a'
  END IF
  
  SELECT CASE ( w )
  CASE ( 'd', 'D' )
     IF ( ASSOCIATED ( spl%x)) THEN
        DEALLOCATE (spl%x,STAT=err)
        IF (err /= 0 ) CALL stop_memory ( 'kill_spline', 'spl%x' )
        NULLIFY (spl%x)
     END IF
     IF ( ASSOCIATED ( spl%y)) THEN
        DEALLOCATE (spl%y,STAT=err)
        IF (err /= 0 ) CALL stop_memory ( 'kill_spline', 'spl%y' )
        NULLIFY (spl%y)
     END IF
     
  CASE ( '2')
     IF ( ASSOCIATED ( spl%y2)) THEN
        DEALLOCATE (spl%y2,STAT=err)
        IF (err /= 0 ) CALL stop_memory ( 'kill_spline', 'spl%y2' )
        NULLIFY (spl%y2)
     END IF
     
  CASE ( 'a', 'A' )
     IF ( ASSOCIATED ( spl%x)) THEN
        DEALLOCATE (spl%x,STAT=err)
        IF (err /= 0 ) CALL stop_memory ( 'kill_spline', 'spl%x' )
        NULLIFY (spl%x)
     END IF
     IF ( ASSOCIATED ( spl%y)) THEN
        DEALLOCATE (spl%y,STAT=err)
        IF (err /= 0 ) CALL stop_memory ( 'kill_spline', 'spl%y' )
        NULLIFY (spl%y)
     END IF
     IF ( ASSOCIATED ( spl%y2)) THEN
        DEALLOCATE (spl%y2,STAT=err)
        IF (err /= 0 ) CALL stop_memory ( 'kill_spline', 'spl%y2' )
        NULLIFY (spl%y2)
     END IF
  END SELECT
  
END SUBROUTINE kill_spline

!******************************************************************************

END MODULE splines
