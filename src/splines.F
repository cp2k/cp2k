!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****h* cp2k/splines [1.0] *
!!
!!   NAME
!!     splines
!!
!!   FUNCTION
!!     routines for handling splines
!!
!!   AUTHOR
!!     various
!!
!!   MODIFICATION HISTORY
!!     2001-09-21-HAF added this doc entry and changed formatting
!!
!!   SOURCE
!******************************************************************************

MODULE splines

  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dbl
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE pair_potential_types,            ONLY: pair_potential_type

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='fist_environment_types'
  INTEGER, PRIVATE, SAVE :: last_spline_env_id_nr=0

  PUBLIC :: spline_env_release, spline_env_retain, spline_env_set, &
            spline_environment_type
  PUBLIC :: spline_env_create
  PUBLIC :: spline_data      ! the data structure for splines
  PUBLIC :: init_splinexy    ! allocates x and y vectors for splines
  PUBLIC :: init_spline      ! generate table for spline (allocates y2)
  PUBLIC :: spline           ! return value of spline for given abscissa
                             ! (optional:also y1)
  PUBLIC :: spline_1         ! return value of 1. derivative of spline
                             ! for given abscissa
  PUBLIC :: spline_int       ! return value of integral on given interval
                             ! of spline
  PUBLIC :: kill_spline      ! destructor ( spl%x,y und/oder spl%y2)
  PUBLIC :: splineh          ! return value of spline and 1. derivative
                             ! without checks (fast routine for pair_potential)

!******************************************************************************
!!****
!!****s* splines/spline_data [1.0] *
!!
!!   NAME
!!     spline_data
!!
!!   FUNCTION
!!     Data-structure that holds all needed information about
!!     a specific spline interpolation.
!!
!!   AUTHOR
!!     unknown
!!
!!   MODIFICATION HISTORY
!!     2001-09-19-HAF added this doc entry and changed formatting
!!
!!   SOURCE
!******************************************************************************

  TYPE spline_data
     REAL ( dbl ), POINTER :: x  ( : )    ! the x values
     REAL ( dbl ), POINTER :: y  ( : )    ! the function values y(x)
     REAL ( dbl ), POINTER :: y2 ( : )    ! the 2nd derivative via interpolation
     INTEGER               :: n           ! dimension of above arrays
     INTEGER               :: pos         ! index for arrays for last access
                                          ! not used if uniform increments
     REAL ( dbl )          :: h           ! uniform increment of x if applicable
     REAL ( dbl )          :: invh        ! inverse of h
     REAL ( dbl )          :: h26         ! 1/6 * h**2 if uniform increments
                                          ! 1/6 otherwise
     REAL ( dbl )          :: x1          ! starting x value if uniform incr.
     REAL ( dbl )          :: xn          ! end x value if uniform incr.
  END TYPE spline_data

  TYPE spline_environment_type
     INTEGER :: ref_count, id_nr
     TYPE ( spline_data ), POINTER, DIMENSION ( : ) :: spline_table
     TYPE ( pair_potential_type ), POINTER, DIMENSION ( :, : ) :: potparm
     INTEGER, POINTER, DIMENSION ( :, : ) :: spltab
     INTEGER  :: nspline_points
  END TYPE spline_environment_type
!******************************************************************************
!!*****

CONTAINS

!******************************************************************************
!!****
!!****s* splines/spline_env_set [1.0] *
!!
!!   NAME
!!     spline_env_set
!!
!!   FUNCTION
!!     set spline_env
!!   AUTHOR
!!     unknown
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE spline_env_set ( spline_env, potparm, spline_table, spltab,  &
                            nspline_points, error)
  TYPE(spline_environment_type), POINTER :: spline_env
  TYPE ( pair_potential_type ), POINTER, OPTIONAL :: potparm ( :, : )
  TYPE ( spline_data ), POINTER, OPTIONAL :: spline_table ( : )
  INTEGER, POINTER, OPTIONAL :: spltab ( :, : )
  INTEGER, OPTIONAL :: nspline_points
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  CHARACTER(len=*), PARAMETER :: routineN='spline_env_set',&
       routineP=moduleN//':'//routineN

  IF ( PRESENT ( potparm ) ) spline_env % potparm => potparm
  IF ( PRESENT ( spline_table ) ) spline_env % spline_table => spline_table
  IF ( PRESENT ( spltab ) ) spline_env % spltab => spltab
  IF ( PRESENT ( nspline_points ) ) spline_env % nspline_points = nspline_points
END SUBROUTINE spline_env_set
!******************************************************************************
!!****
!!****s* splines/spline_env_release [1.0] *
!!
!!   NAME
!!     spline_env_release
!!
!!   FUNCTION
!!     releases spline_env
!!   AUTHOR
!!     unknown
!!
!!   MODIFICATION HISTORY
!!     2001-09-19-HAF added this doc entry and changed formatting
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE spline_env_release(spline_env,error)
  TYPE(spline_environment_type), POINTER :: spline_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='spline_env_release',&
       routineP=moduleN//':'//routineN
  INTEGER :: i, stat
  failure=.FALSE.

  IF (ASSOCIATED(spline_env)) THEN
     CPPrecondition(spline_env%ref_count>0,cp_failure_level,routineP,error,failure)
     spline_env%ref_count=spline_env%ref_count-1
     IF (spline_env%ref_count<1) THEN
        DEALLOCATE (  spline_env % spltab, stat = stat )
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DO i = 1, SIZE (  spline_env % spline_table )
          CALL kill_spline ( spline_env % spline_table ( i )  )
        END DO
        DEALLOCATE (  spline_env % spline_table, stat = stat )
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE ( spline_env % potparm, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        NULLIFY ( spline_env % spltab )
        NULLIFY ( spline_env % spline_table )
        NULLIFY ( spline_env % potparm ) 
     END IF
  END IF
  NULLIFY ( spline_env )
END SUBROUTINE spline_env_release
!***************************************************************************

!******************************************************************************
!!****
!!****s* splines/spline_env_retain [1.0] *
!!
!!   NAME
!!     spline_env_retain
!!
!!   FUNCTION
!!     retains spline_env
!!    
!!
!!   AUTHOR
!!     unknown
!!
!!   MODIFICATION HISTORY
!!     2001-09-19-HAF added this doc entry and changed formatting
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE spline_env_retain(spline_env,error)
  TYPE(spline_environment_type), POINTER :: spline_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
 
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='spline_env_retain',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
!!  
  CPPrecondition(ASSOCIATED(spline_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(spline_env%ref_count>0,cp_failure_level,routineP,error,failure)
     spline_env%ref_count=spline_env%ref_count+1
  END IF
END SUBROUTINE spline_env_retain

!******************************************************************************
!!****
!!****s* splines/spline_env_create [1.0] *
!!
!!   NAME
!!     spline_env_create
!!
!!   FUNCTION
!!     Data-structure that holds all needed information about
!!     a specific spline interpolation.
!!
!!   AUTHOR
!!     unknown
!!
!!   MODIFICATION HISTORY
!!     2001-09-19-HAF added this doc entry and changed formatting
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE spline_env_create ( spline_env, error )
  TYPE(spline_environment_type), POINTER :: spline_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='spline_env_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat

  failure=.FALSE.
  
  ALLOCATE(spline_env, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
    spline_env%ref_count=1
    last_spline_env_id_nr=last_spline_env_id_nr+1
    spline_env%id_nr=last_spline_env_id_nr
    NULLIFY ( spline_env % spline_table )
    NULLIFY ( spline_env % spltab )
    NULLIFY ( spline_env % potparm )
    spline_env % nspline_points = 2500
  END IF
END SUBROUTINE spline_env_create

!******************************************************************************
!!****m* splines/init_splinexy [1.0] *
!!
!!   NAME
!!     init_splinexy
!!
!!   FUNCTION
!!     allocates storage for function table to be interpolated
!!     both x and y are allocated
!!
!!   AUTHOR
!!     unknown
!!
!!   INPUTS
!!     spl : spline_data structure to be initialized
!!     nn  : integer number of datapoints, that the function table will hold
!!
!!   MODIFICATION HISTORY
!!     2001-09-21-HAF added this doc entry and changed formatting
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE init_splinexy ( spl, nn )


    TYPE(spline_data), INTENT(INOUT)         :: spl
    INTEGER, INTENT(IN)                      :: nn

    INTEGER                                  :: err

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

  ALLOCATE ( spl % y ( 1:nn ), STAT = err )
  IF ( err /= 0 ) CALL stop_memory ( 'init_splinexy', 'spl%y', nn)

END SUBROUTINE init_splinexy

!******************************************************************************
!!****
!!****m* splines/init_spline [1.0] *
!!
!!   NAME
!!     init_spline
!!
!!   FUNCTION
!!     allocates storage for y2 table 
!!     calculates y2 table and other spline parameters
!!
!!   AUTHOR
!!     unknown
!!
!!   INPUTS
!!     spl : spline_data structure to be initialized
!!                      spl%y() must hold the function values
!!                      spl%x() must hold the absissa values in increasing
!!                      order OR if dx (below) is given, spl%x(1) must hold
!!                      the starting (left-most) point.
!!     dx  : (OPTIONAL) if present, only x(1) need to be valid, all other
!!                      x(i) are assumed to be x(1)+dx*(i-1)
!!                      (spline evaluations will also be faster)
!!     y1a : (OPTIONAL) if present, the 1-deriv of the left endpoint
!!                      if not present, natural spline condition at this end
!!                      (2-deriv == 0)
!!     y1b : (OPTIONAL) if present, the 1-deriv of the right endpoint
!!                      if not present, natural spline condition at this end
!!                      (2-deriv == 0)
!!
!!   EXAMPLES
!!     CALL init_spline(spline,dx=0.1_dbl)
!!     CALL init_spline(spline,y1b=0.0_dbl)
!!     CALL init_spline(spline,0.1_dbl,0.0_dbl,0.0_dbl)
!!
!!   NOTES
!!     if dx is given, the x array will be used as y2 array instead of
!!     allocating a new array. (y2 will become x, and x will be nullified)
!!
!!   MODIFICATION HISTORY
!!     2001-09-21-HAF added this doc entry and changed formatting
!!     2001-09-24-HAF changed interface and re-written
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE init_spline ( spl, dx, y1a, y1b )


    TYPE(spline_data), INTENT(INOUT)         :: spl
    REAL(dbl), INTENT(IN), OPTIONAL          :: dx, y1a, y1b

    INTEGER                                  :: err, i, n
    REAL(dbl)                                :: p, s
    REAL(dbl), POINTER                       :: ww( : )

!------------------------------------------------------------------------------

  IF ( ASSOCIATED ( spl%y2)) THEN
     DEALLOCATE (spl%y2,STAT=err)
     IF (err /= 0 ) CALL stop_memory ( 'init_spline','spl%y2')
     NULLIFY (spl%y2)
  END IF

  n = spl%n

  IF (PRESENT(dx)) THEN
     spl % x1   = spl%x(1)
     spl % xn   = spl%x(1)+ (n-1)*dx
     spl % pos  = 1
     spl % h    = dx
     spl % invh = 1.0_dbl/dx
     spl % h26  = dx**2 / 6.0_dbl
     spl % y2  => spl % x
     NULLIFY (spl%x)
     ALLOCATE (ww(1:n),STAT=err)
     IF ( err /= 0 ) CALL stop_memory ( 'init_spline', 'ww', spl % n )
     IF (PRESENT(y1a)) THEN
        spl % y2(1) = -0.5_dbl
        ww(1) = 3.0_dbl * ((spl%y(2)-spl%y(1))/dx - y1a)/dx
     ELSE
        spl % y2(1) = 0.0_dbl
        ww(1)       = 0.0_dbl
     END IF
     DO i = 2, n - 1
        s = 0.5
        p = 0.5_dbl * spl%y2(i-1) + 2.0_dbl
        spl%y2(i) = - 0.5_dbl / p
        ww(i) = ( 3.0_dbl*( spl%y(i+1)-2.0_dbl*spl%y(i)+spl%y(i-1) )/(dx*dx) &
                 - 0.5_dbl*ww(i-1) )/p
     END DO
     IF (PRESENT( y1b )) THEN
        spl%y2(n) = ( 3.0_dbl * (y1b - (spl%y(n)-spl%y(n-1)) / dx ) / dx - &
                      0.5_dbl * ww(n-1) ) / ( 0.5_dbl * spl%y2(n-1) + 1.0_dbl )
     ELSE
        spl%y2(n) = 0.0_dbl
     END IF
     DO i = n-1, 1, -1
        spl%y2(i) = spl%y2(i) * spl%y2(i+1) + ww(i)
     END DO
     DEALLOCATE (ww,STAT=err)
     IF (err /= 0 ) CALL stop_memory ( 'init_spline','ww')
  ELSE
     ALLOCATE (spl%y2(1:n),STAT=err)
     IF ( err /= 0 ) CALL stop_memory ( 'init_spline', 'spl%y2', spl % n )
     ALLOCATE (ww(1:n),STAT=err)
     IF ( err /= 0 ) CALL stop_memory ( 'init_spline', 'ww', spl % n )
     IF (PRESENT(y1a)) THEN
        spl % y2(1) = -0.5_dbl
        ww(1) = 3.0_dbl * ((spl%y(2)-spl%y(1))/(spl%x(2)-spl%x(1)) - y1a)/ &
                (spl%x(2)-spl%x(1))
     ELSE
        spl % y2(1) = 0.0_dbl
        ww(1)       = 0.0_dbl
     END IF
     DO i = 2, n - 1
        s = ( spl%x(i) - spl%x(i-1) ) / ( spl%x(i+1) - spl%x(i-1) )
        p = s * spl%y2(i-1) + 2.0_dbl
        spl%y2(i) = (s - 1.0_dbl) / p
        ww(i) = (6.0_dbl*( (spl%y(i+1)-spl%y(i)) / (spl%x(i+1)-spl%x(i)) - &
                           (spl%y(i)-spl%y(i-1)) / (spl%x(i)-spl%x(i-1)) ) / &
                         (spl%x(i+1) - spl%x(i-1)) - s*ww(i-1) )/p
     END DO
     IF (PRESENT( y1b )) THEN
        spl%y2(n) = ( 3.0_dbl * (y1b - (spl%y(n)-spl%y(n-1)) / &
                                       (spl%x(n)-spl%x(n-1)) ) / &
                      (spl%x(n)-spl%x(n-1)) - 0.5_dbl * ww(n-1) ) / &
                    ( 0.5_dbl * spl%y2(n-1) + 1.0_dbl )
     ELSE
        spl%y2(n) = 0.0_dbl
     END IF
     DO i = n-1, 1, -1
        spl%y2(i) = spl%y2(i) * spl%y2(i+1) + ww(i)
     END DO
     DEALLOCATE (ww,STAT=err)
     IF (err /= 0 ) CALL stop_memory ( 'init_spline','ww')
     spl % pos = 1
     spl % h    = 0.0_dbl
     spl % invh = 0.0_dbl
     spl % h26  = 1.0_dbl / 6.0_dbl
     spl % x1   = spl % x(1)
     spl % xn   = spl % x(n)
  END IF

END SUBROUTINE init_spline

!******************************************************************************
!!****
!!****i* splines/interv [1.2] *
!!
!!   NAME
!!     interv
!!
!!   FUNCTION
!!     calculates index to function table array for given absissa
!!
!!   AUTHOR
!!     unknown
!!
!!   INPUTS
!!     spl : spline_data structure
!!     xx  : absissa value to find nearest index for
!!
!!   OUTPUT
!!     nearest index for function tables
!!
!!   MODIFICATION HISTORY
!!     2001-09-24-HAF added this doc entry and changed formatting
!!     2001-09-24-HAF re-written to use mixed binary-interpolation search
!!                    if all else fails
!!     2002-11-18-HAF fixed return value in case of xx = x(n)
!!
!!   SOURCE
!******************************************************************************

FUNCTION interv ( spl, xx )


    TYPE(spline_data), INTENT(IN)            :: spl
    REAL(dbl), INTENT(IN)                    :: xx
    INTEGER                                  :: interv

    INTEGER                                  :: i, k, khi, klo, n, p
    REAL(dbl)                                :: xhi, xk, xlo

! Return value
!------------------------------------------------------------------------------

  IF (spl%h /= 0 ) THEN
     i = (xx-spl%x1)*spl%invh + 1
     IF (i<1 .OR. i>spl%n-1) THEN
        IF (xx > spl%xn) CALL stop_program ( 'interv', &
             'illegal x-value passed to interv')
        i = spl%n-1
     END IF
     interv = i
     RETURN
  END IF

  p = spl%pos
  IF ( (p >= spl%n) .OR. (p < 1) ) p = (spl%n+1)/2
  i = 0
  n = spl%n

! check if interval is close to previous interval
!
! How efficient is this? I would think in an MD simulation the xx
! are pretty random since it is not usually the same distance the next
! time step, but a different particle pair distance, etc... - HAF
!
  IF ( xx < spl % x ( p + 1 ) ) THEN
     IF ( xx >= spl % x ( p ) ) THEN
        i = p
     ELSE IF ( (p > 1) .AND. (xx >= spl % x ( p - 1 )) ) THEN
        i = p - 1
     ELSE
        klo = 1
        khi = p - 1
     END IF
  ELSE IF ( (p + 2 <= n) .AND. (xx < spl % x ( p + 2 )) ) THEN
     i = p + 1
  ELSE
     klo = p + 2
     khi = n
  END IF

  IF ( i == 0 ) THEN ! perform general search (binary + interpolation)
                     ! lots of other possibilities here though

     IF ( (xx < spl % x ( 1 )) .OR. (xx > spl % x ( n )) ) &
          CALL stop_program ( 'interv', 'xx value out of spline-range' )

     xlo = spl % x ( klo )
     xhi = spl % x ( khi )
     DO WHILE ( khi - klo > 1 )

        ! linear interpolation prediction:
        k = klo + (xx - xlo)*(khi-klo)/(xhi-xlo)
        IF (spl % x(k) >= xx) THEN
           khi = k
           xhi = spl % x(k)
        ELSE IF (spl % x(k+1) <= xx) THEN
           klo = k+1
           xlo = spl % x(k+1)
        ELSE
           klo = k
           EXIT
        END IF
        IF ( khi - klo <= 1 ) EXIT

        ! binary search:
        k = ( khi + klo ) / 2
        xk = spl % x ( k )
        IF ( xk > xx ) THEN
           khi = k
           xhi = xk
        ELSE
           klo = k
           xlo = xk
        END IF
     END DO
     i = MIN(klo,n-1)

  END IF

  interv = i

END FUNCTION interv

!******************************************************************************
!!****
!!****f* splines/spline [1.0] *
!!
!!   NAME
!!     spline
!!
!!   FUNCTION
!!     calculates the spline interpolated value at a given point
!!
!!   AUTHOR
!!     unknown
!!
!!   INPUTS
!!     spl : spline_data structure
!!     xx  : absissa value
!!
!!   OUTPUT
!!     spline interpolated value at xx
!!     OPTIONAL y1 : 1. derivative at xx
!!     OPTIONAL y2 : 2. derivative at xx
!!     OPTIONAL y3 : 3. derivative at xx
!!
!!   MODIFICATION HISTORY
!!     2001-09-25-HAF added this doc entry and changed formatting
!!
!!   SOURCE
!******************************************************************************

FUNCTION spline ( spl, xx, y1, y2, y3 )


! Return value
    TYPE(spline_data), INTENT(INOUT)         :: spl
    REAL(dbl), INTENT(IN)                    :: xx
    REAL(dbl), INTENT(OUT), OPTIONAL         :: y1, y2, y3
    REAL(dbl)                                :: spline

    INTEGER                                  :: i
    REAL(dbl)                                :: a, b, h, h26, invh, y2hi, &
                                                y2lo, yhi, ylo

!------------------------------------------------------------------------------

  if (spl % h /= 0.0_dbl) THEN
     ! index calculation from interv for speed reasons duplicated in this
     ! simple case:
     i = (xx-spl%x1)*spl%invh + 1
     IF (i<1 .OR. i>spl%n) THEN
        IF (xx > spl%xn) CALL stop_program ( 'spline', &
             'illegal x-value passed to spline')
        i = spl % n - 1
     END IF
     h    = spl % h
     invh = spl%invh
     h26  = spl%h26
     a = (spl%x1-xx)*invh+i
  ELSE
     ! general non-uniform case:
     i = interv(spl,xx)
     spl % pos = i
     h = spl%x(i+1)-spl%x(i)
     invh = 1.0_dbl/h
     h26  = h*h*spl%h26
     a = (spl%x(i+1)-xx)*invh
  END IF
  b = 1.0_dbl - a
  ylo = spl%y(i)
  yhi = spl%y(i+1)
  y2lo = spl%y2(i)
  y2hi = spl%y2(i+1)
  spline = a*ylo + b*yhi - ((a+1)*y2lo+(b+1)*y2hi)*a*b*h26
  IF ( PRESENT (y1) ) &
    y1 = invh*( yhi - ylo + ( (1-3*a*a)*y2lo + (3*b*b-1)*y2hi ) * h26 )
  IF ( PRESENT (y2) ) &
    y2 = a*y2lo + b*y2hi
  IF ( PRESENT (y3) ) &
    y3 = invh * ( y2hi - y2lo )
  RETURN

END FUNCTION spline

!******************************************************************************
!!****
!!****f* splines/splineh [1.0] *
!!
!!   NAME
!!     splineh
!!
!!   FUNCTION
!!     calculates the spline interpolated value at a given point
!!     and the first derivative
!!
!!   AUTHOR
!!     unknown
!!
!!   INPUTS
!!     spl : spline_data structure
!!     xx  : absissa value
!!
!!   OUTPUT
!!     spline interpolated value at xx
!!     y1 : 1. derivative at xx
!!
!!   NOTES
!!     the spline MUST have uniform x values and xx MUST be
!!     in the interpolation interval. No checks are done to ensure
!!     either condition.
!!
!!   MODIFICATION HISTORY
!!     2001-09-25-HAF added this doc entry and changed formatting
!!
!!   SOURCE
!******************************************************************************

FUNCTION splineh ( spl, xx, y1 )


! Return value
    TYPE(spline_data), INTENT(IN)            :: spl
    REAL(dbl), INTENT(IN)                    :: xx
    REAL(dbl), INTENT(OUT)                   :: y1
    REAL(dbl)                                :: splineh

    INTEGER                                  :: i
    REAL(dbl)                                :: a, b, h26, invh, y2hi, y2lo, &
                                                yhi, ylo

!------------------------------------------------------------------------------
! fast spline for pair potentials without checks

  h26 = spl%h26
  invh = spl%invh
  i = (xx-spl%x1)*invh + 1

  a = (spl%x1-xx)*invh + i
  b = 1 - a
  ylo = spl%y(i)
  yhi = spl%y(i+1)
  y2lo = spl%y2(i)
  y2hi = spl%y2(i+1)
  splineh = a*ylo + b*yhi - ((a+1)*y2lo+(b+1)*y2hi)*a*b*h26
  y1 = invh*((yhi-ylo) + ((1.-3*a*a)*y2lo+(3*b*b-1.)*y2hi)*h26)

END FUNCTION splineh

!******************************************************************************
!!****
!!****f* splines/spline_1 [1.0] *
!!
!!   NAME
!!     spline_1
!!
!!   FUNCTION
!!     calculates the 1. derivate of a spline at a given point
!!
!!   AUTHOR
!!     unknown
!!
!!   INPUTS
!!     spl : spline_data structure
!!     xx  : absissa value to find nearest index for
!!
!!   OUTPUT
!!     spline_1 first derivate of the spline interpolation at xx
!!
!!   MODIFICATION HISTORY
!!     2001-09-25-HAF added this doc entry and changed formatting
!!
!!   SOURCE
!******************************************************************************

FUNCTION spline_1 ( spl, xx )


! Return value
    TYPE(spline_data), INTENT(INOUT)         :: spl
    REAL(dbl), INTENT(IN)                    :: xx
    REAL(dbl)                                :: spline_1

    INTEGER                                  :: i
    REAL(dbl)                                :: a, b, h, h26, invh

!------------------------------------------------------------------------------

  if (spl % h /= 0.0_dbl) THEN
     ! index calculation from interv for speed reasons duplicated in this
     ! simple case:
     i = (xx-spl%x1)*spl%invh + 1
     IF (i<1 .OR. i>spl%n) THEN
        IF ( xx > spl%xn ) CALL stop_program ( 'spline', &
             'illegal x-value passed to spline')
        i = spl % n - 1
     END IF
     h    = spl % h
     invh = spl%invh
     h26  = spl%h26
     a = (spl%x1-xx)*invh+i
  ELSE
     ! general non-uniform case:
     i = interv(spl,xx)
     spl % pos = i
     h = spl%x(i+1)-spl%x(i)
     invh = 1.0_dbl/h
     h26  = h*h*spl%h26
     a = (spl%x(i+1)-xx)*invh
  END IF
  b = 1.0_dbl - a
  spline_1 = invh * ( spl%y(i+1)-spl%y(i) + &
                      ( (1-3*a*a)*spl%y2(i) + (3*b*b-1)*spl%y2(i+1) ) * h26 )
  RETURN

END FUNCTION spline_1

!******************************************************************************
!!****
!!****i* splines/stamm [1.0] *
!!
!!   NAME
!!     stamm
!!
!!   FUNCTION
!!     calculates the value of the undetermined integral of a spline at
!!     a given point for a given interval
!!
!!   AUTHOR
!!     unknown
!!
!!   INPUTS
!!     spl : spline_data structure
!!     p   : index of the spline interval to integrate
!!     xx  : absissa value to find nearest index for
!!
!!   OUTPUT
!!     value of the undetermined integral at xx for the interval p
!!
!!   MODIFICATION HISTORY
!!     2001-09-25-HAF added this doc entry and changed formatting
!!
!!   SOURCE
!******************************************************************************

FUNCTION stamm ( spl, p, xx )


! Return value
    TYPE(spline_data), INTENT(IN)            :: spl
    INTEGER, INTENT(IN)                      :: p
    REAL(dbl), INTENT(IN)                    :: xx
    REAL(dbl)                                :: stamm

    REAL(dbl)                                :: a, aa, b, bb, h, h26

!------------------------------------------------------------------------------

  h = spl%h
  IF (h == 0.0_dbl) THEN
     h = spl%x(p+1) - spl%x(p)
     h26 = h*h*spl%h26
     b = (xx-spl%x(p))/h
  ELSE
     h26 = spl%h26
     b = (xx-spl%x1)*spl%invh - p + 1
  END IF
  a = 1 - b
  aa = a**2
  bb = b**2
  stamm = 0.5*h*( bb*spl%y(p+1) - aa*spl%y(p) + &
                  h26 * ( aa*(1-0.5*aa)*spl%y2(p) - &
                          bb*(1-0.5*bb)*spl%y2(p+1) ) )

END FUNCTION stamm

!******************************************************************************
!!****
!!****f* splines/spline_int [1.0] *
!!
!!   NAME
!!     spline_int
!!
!!   FUNCTION
!!     calculates the integral of a spline for two given bounds
!!
!!   AUTHOR
!!     unknown
!!
!!   INPUTS
!!     spl : spline_data structure
!!     x0  : lower bound
!!     x1  : upper bound
!!
!!   OUTPUT
!!     the integrarl of the spline between x0 and x1
!!     /int_{x0}^{x1} spline(x) dx
!!
!!   MODIFICATION HISTORY
!!     2001-09-25-HAF added this doc entry and changed formatting
!!
!!   SOURCE
!******************************************************************************

FUNCTION spline_int( spl, x0, x1 )


! Return value
    TYPE(spline_data), INTENT(INOUT)         :: spl
    REAL(dbl), INTENT(IN)                    :: x0, x1
    REAL(dbl)                                :: spline_int

    INTEGER                                  :: j, pa, pb
    REAL(dbl)                                :: h, i1, i2, vorz, xa, xb

!------------------------------------------------------------------------------

  vorz = 1.0_dbl
  xa = MIN ( x0, x1 )
  xb = MAX ( x0, x1 )
  IF ( x0 > x1 ) vorz = -1.0_dbl

  h = spl % h
  if (h /= 0.0_dbl) THEN

     IF ( (xa < spl % x1) .OR. (xb > spl % xn) ) &
        CALL stop_program ( 'spline_int', 'illegal integration range' )

     ! index calculation from interv for speed reasons duplicated in this
     ! simple case:
     pa = MIN(INT((xa-spl%x1)*spl%invh + 1), spl%n - 1)
     pb = MIN(INT((xb-spl%x1)*spl%invh + 1), spl%n - 1)

     IF ( pa == pb ) THEN
        spline_int = vorz * ( stamm(spl,pa,xb) - stamm(spl,pa,xa) )
        RETURN
     END IF

     i1 = 0
     i2 = 0

     DO j = pa + 1, pb
        i1 = i1 + spl%y(j)
        i2 = i2 + spl%y2(j)
     END DO

     spline_int = vorz*(h*(i1-0.5*spl%h26*i2)+stamm(spl,pb,xb)-stamm(spl,pa,xa))

     RETURN

  ELSE

     ! general non-uniform case:

     IF ( xa < spl % x ( 1 ) .OR. xb > spl % x ( spl % n ) ) &
        CALL stop_program ( 'spline_int', 'illegal integration range' )

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
        i1 = i1 + h    * ( spl%y(j)  + spl%y(j+1)  )
        i2 = i2 + h**3 * ( spl%y2(j) + spl%y2(j+1) )
     END DO

     h = spl%x(pa+1) - spl%x(pa)
     i1 = i1 + h*spl%y(pa+1)
     i2 = i2 + h**3*spl%y2(pa+1)
     h = spl%x(pb+1) - spl%x(pb)
     i1 = i1 + h*spl%y(pb)
     i2 = i2 + h**3*spl%y2(pb)

     spline_int = vorz*(i1/2.-i2/24.+stamm(spl,pb,xb)-stamm(spl,pa,xa))

     RETURN

  END IF

END FUNCTION spline_int

!******************************************************************************
!!****
!!****m* splines/kill_spline [1.0] *
!!
!!   NAME
!!     kill_spline
!!
!!   FUNCTION
!!     deallocates (selected) storage of given spline data structure
!!
!!   AUTHOR
!!     unknown
!!
!!   INPUTS
!!     spl   : spline_data structure
!!     what  : OPTIONAL what should be deallocated:
!!                not present or 'a' : deallocate all (spl%x spl%y and spl%y2)
!!                'd'                : deallocate only data (spl%x and spl%y)
!!                '2'                : deallocate only 2. derivatives (spl%y2)
!!
!!   MODIFICATION HISTORY
!!     2001-09-26-HAF added this doc entry and changed formatting
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE kill_spline ( spl, what )


    TYPE(spline_data), INTENT(INOUT)         :: spl
    CHARACTER, INTENT(IN), OPTIONAL          :: what

    CHARACTER                                :: w
    INTEGER                                  :: err

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

  CASE DEFAULT
     CALL stop_program ( 'kill_spline', 'illegal deallocation option' )

  END SELECT

END SUBROUTINE kill_spline

!******************************************************************************
!!****

END MODULE splines
