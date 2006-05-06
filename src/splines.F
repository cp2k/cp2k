!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
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

  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_memory,&
       stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'splines'
  INTEGER, PRIVATE, SAVE :: last_spline_env_id_nr=0
  INTEGER, PRIVATE, SAVE :: last_spline_data_id_nr=0

  PUBLIC :: spline_env_release, spline_env_retain, spline_environment_type
  PUBLIC :: spline_env_create, spline_data_p_type
  PUBLIC :: spline_data_create, spline_data_copy
  PUBLIC :: spline_data_retain, spline_data_release
  PUBLIC :: spline_data_type ! the data structure for spline table
  PUBLIC :: init_splinexy    ! allocates x and y vectors for splines
  PUBLIC :: init_spline      ! generate table for spline (allocates y2)
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

  TYPE spline_data_type
     INTEGER :: ref_count, id_nr
     REAL (KIND=dp), POINTER :: x  ( : )    ! the x values
     REAL (KIND=dp), POINTER :: y  ( : )    ! the function values y(x)
     REAL (KIND=dp), POINTER :: y2 ( : )    ! the 2nd derivative via interpolation
     INTEGER               :: n             ! dimension of above arrays
     INTEGER               :: pos           ! index for arrays for last access
     ! not used if uniform increments
     REAL (KIND=dp)          :: h           ! uniform increment of x if applicable
     REAL (KIND=dp)          :: invh        ! inverse of h
     REAL (KIND=dp)          :: h26         ! 1/6 * h**2 if uniform increments
     ! 1/6 otherwise
     REAL (KIND=dp)          :: x1          ! starting x value if uniform incr.
     REAL (KIND=dp)          :: xn          ! end x value if uniform incr.
  END TYPE spline_data_type

  TYPE spline_data_p_type
     TYPE ( spline_data_type ), POINTER :: spline_data
  END TYPE spline_data_p_type

  TYPE spline_environment_type
     INTEGER :: ref_count, id_nr
     TYPE ( spline_data_p_type ), POINTER, DIMENSION ( : ) :: spline_data_p
     INTEGER, POINTER, DIMENSION ( :, : ) :: spltab
  END TYPE spline_environment_type
  !******************************************************************************
  !!*****

CONTAINS

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
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE spline_env_release(spline_env,error)
    TYPE(spline_environment_type), POINTER   :: spline_env
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_env_release', &
         routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure
    TYPE(spline_data_type), POINTER          :: spline_data

    failure=.FALSE.

    IF (ASSOCIATED(spline_env)) THEN
       CPPrecondition(spline_env%ref_count>0,cp_failure_level,routineP,error,failure)
       spline_env%ref_count=spline_env%ref_count-1
       IF (spline_env%ref_count<1) THEN
          DEALLOCATE (  spline_env % spltab, stat = stat )
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DO i = 1, SIZE (  spline_env % spline_data_p )
             spline_data => spline_env % spline_data_p ( i ) % spline_data
             CALL spline_data_release ( spline_data, error )
          END DO
          DEALLOCATE (  spline_env % spline_data_p, stat = stat )
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          NULLIFY ( spline_env % spltab )
          NULLIFY ( spline_env % spline_data_p )
          DEALLOCATE (  spline_env , stat = stat )
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)        
       END IF
    END IF
    NULLIFY ( spline_env )
  END SUBROUTINE spline_env_release
  !***************************************************************************

  !******************************************************************************
  !!****
  !!****s* splines/spline_data_release [1.0] *
  !!
  !!   NAME
  !!     spline_data_release
  !!
  !!   FUNCTION
  !!     releases spline_data
  !!    
  !!
  !!   AUTHOR
  !!     CJM
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE spline_data_release (spline_data,error)
    TYPE(spline_data_type), POINTER          :: spline_data
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_data_release', &
         routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(spline_data)) THEN
       CPPrecondition(spline_data%ref_count>0,cp_failure_level,routineP,error,failure)
       spline_data%ref_count=spline_data%ref_count-1
       IF (spline_data%ref_count<1) THEN
          IF ( ASSOCIATED ( spline_data%x)) THEN
             DEALLOCATE (spline_data%x,STAT=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             NULLIFY (spline_data%x)
          END IF
          IF ( ASSOCIATED ( spline_data%y)) THEN
             DEALLOCATE (spline_data%y,STAT=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             NULLIFY (spline_data%y)
          END IF
          IF ( ASSOCIATED ( spline_data%y2)) THEN
             DEALLOCATE (spline_data%y2,STAT=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             NULLIFY (spline_data%y2)
          END IF
          DEALLOCATE (  spline_data, stat = stat )
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY ( spline_data )
  END SUBROUTINE spline_data_release
  !******************************************************************************
  !!****
  !!****s* splines/spline_data_retain [1.0] *
  !!
  !!   NAME
  !!     spline_env_retain
  !!
  !!   FUNCTION
  !!     retains spline_env
  !!    
  !!
  !!   AUTHOR
  !!     CJM
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE spline_data_retain(spline_data,error)
    TYPE(spline_data_type), POINTER          :: spline_data
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_data_retain', &
         routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    !!  
    CPPrecondition(ASSOCIATED(spline_data),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(spline_data%ref_count>0,cp_failure_level,routineP,error,failure)
       spline_data%ref_count=spline_data%ref_count+1
    END IF
  END SUBROUTINE spline_data_retain

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
    TYPE(spline_environment_type), POINTER   :: spline_env
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_env_retain', &
         routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
  SUBROUTINE spline_env_create ( spline_env, ntype, error )
    TYPE(spline_environment_type), POINTER   :: spline_env
    INTEGER, INTENT(IN)                      :: ntype
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_env_create', &
         routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ntab, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(spline_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY ( spline_env % spline_data_p )
       NULLIFY ( spline_env % spltab )
       spline_env%ref_count=1
       last_spline_env_id_nr=last_spline_env_id_nr+1
       spline_env%id_nr=last_spline_env_id_nr
       ! Allocate the number of spline data tables (upper triangular)
       ntab = ( ntype * ntype + ntype ) / 2
       ALLOCATE ( spline_env % spline_data_p ( ntab ),  stat=stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       ALLOCATE ( spline_env % spltab ( ntype, ntype ), stat=stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       DO i = 1, ntab
          NULLIFY ( spline_env % spline_data_p ( i ) % spline_data )
          CALL spline_data_create ( spline_env % spline_data_p ( i ) % spline_data )
       END DO
    END IF
  END SUBROUTINE spline_env_create
  !******************************************************************************
  !!****
  !!****s* splines/spline_data_copy [1.0] *
  !!
  !!   NAME
  !!     spline_data_copy
  !!
  !!   FUNCTION
  !!     Copy Data-structure that constains spline table
  !!
  !!   AUTHOR
  !!     teo 11.2005
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE spline_data_copy ( spline_data_source, spline_data_dest, error )
    TYPE(spline_data_type), POINTER          :: spline_data_source, &
         spline_data_dest
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_data_copy', &
         routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPostcondition(ASSOCIATED(spline_data_source),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (.NOT.ASSOCIATED(spline_data_dest)) CALL spline_data_create(spline_data_dest,error)

       spline_data_dest%ref_count = spline_data_source%ref_count
       spline_data_dest%id_nr     = spline_data_source%id_nr
       spline_data_dest%n         = spline_data_source%n
       spline_data_dest%pos       = spline_data_source%pos
       spline_data_dest%h         = spline_data_source%h   
       spline_data_dest%invh      = spline_data_source%invh
       spline_data_dest%h26       = spline_data_source%h26 
       spline_data_dest%x1        = spline_data_source%x1  
       spline_data_dest%xn        = spline_data_source%xn  

       IF (ASSOCIATED(spline_data_source%x)) THEN
          ALLOCATE(spline_data_dest%x(SIZE(spline_data_source%x)), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          spline_data_dest%x = spline_data_source%x
       END IF
       IF (ASSOCIATED(spline_data_source%y)) THEN
          ALLOCATE(spline_data_dest%y(SIZE(spline_data_source%y)), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          spline_data_dest%y = spline_data_source%y
       END IF
       IF (ASSOCIATED(spline_data_source%y2)) THEN
          ALLOCATE(spline_data_dest%y2(SIZE(spline_data_source%y2)), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          spline_data_dest%y2 = spline_data_source%y2
       END IF
    END IF
  END SUBROUTINE spline_data_copy
  !******************************************************************************
  !!****
  !!****s* splines/spline_data_create [1.0] *
  !!
  !!   NAME
  !!     spline_data_create
  !!
  !!   FUNCTION
  !!     Data-structure that constains spline table
  !!
  !!   AUTHOR
  !!     unknown
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE spline_data_create ( spline_data, error )
    TYPE(spline_data_type), POINTER          :: spline_data
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_data_create', &
         routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(spline_data, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       spline_data%ref_count=1
       last_spline_data_id_nr=last_spline_data_id_nr+1
       spline_data%id_nr=last_spline_data_id_nr
       NULLIFY ( spline_data % x )
       NULLIFY ( spline_data % y )
       NULLIFY ( spline_data % y2 )
    END IF
  END SUBROUTINE spline_data_create

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


    TYPE(spline_data_type), POINTER          :: spl
    INTEGER, INTENT(IN)                      :: nn

    INTEGER                                  :: err

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
  !!     dx  :            x(i) are assumed to be x(1)+dx*(i-1)
  !!                      (spline evaluations will also be faster)
  !!     y1a : (OPTIONAL) if present, the 1-deriv of the left endpoint
  !!                      if not present, natural spline condition at this end
  !!                      (2-deriv == 0)
  !!     y1b : (OPTIONAL) if present, the 1-deriv of the right endpoint
  !!                      if not present, natural spline condition at this end
  !!                      (2-deriv == 0)
  !!
  !!   EXAMPLES
  !!     CALL init_spline(spline,dx=0.1_dp)
  !!     CALL init_spline(spline,y1b=0.0_dp)
  !!     CALL init_spline(spline,0.1_dp,0.0_dp,0.0_dp)
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

    TYPE(spline_data_type), POINTER          :: spl
    REAL(KIND=dp), INTENT(IN)                :: dx 
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: y1a, y1b

    INTEGER                                  :: err, i, n
    REAL(KIND=dp)                            :: p, s
    REAL(KIND=dp), POINTER                   :: ww( : )

    IF ( ASSOCIATED ( spl%y2)) THEN
       DEALLOCATE (spl%y2,STAT=err)
       IF (err /= 0 ) CALL stop_memory ( 'init_spline','spl%y2')
       NULLIFY (spl%y2)
    END IF

    n          = spl%n
    spl % x1   = spl%x(1)
    spl % xn   = spl%x(1)+ (n-1)*dx
    spl % pos  = 1
    spl % h    = dx
    spl % invh = 1.0_dp/dx
    spl % h26  = dx**2 / 6.0_dp
    spl % y2   => spl % x
    NULLIFY (spl%x)
    ALLOCATE (ww(1:n),STAT=err)
    IF ( err /= 0 ) CALL stop_memory ( 'init_spline', 'ww', spl % n )
    IF (PRESENT(y1a)) THEN
       spl % y2(1) = -0.5_dp
       ww(1) = 3.0_dp * ((spl%y(2)-spl%y(1))/dx - y1a)/dx
    ELSE
       spl % y2(1) = 0.0_dp
       ww(1)       = 0.0_dp
    END IF
    DO i = 2, n - 1
       s = 0.5_dp
       p = 0.5_dp * spl%y2(i-1) + 2.0_dp
       spl%y2(i) = - 0.5_dp / p
       ww(i) = ( 3.0_dp*( spl%y(i+1)-2.0_dp*spl%y(i)+spl%y(i-1) )/(dx*dx) &
            - 0.5_dp*ww(i-1) )/p
    END DO
    IF (PRESENT( y1b )) THEN
       spl%y2(n) = ( 3.0_dp * (y1b - (spl%y(n)-spl%y(n-1)) / dx ) / dx - &
            0.5_dp * ww(n-1) ) / ( 0.5_dp * spl%y2(n-1) + 1.0_dp )
    ELSE
       spl%y2(n) = 0.0_dp
    END IF
    DO i = n-1, 1, -1
       spl%y2(i) = spl%y2(i) * spl%y2(i+1) + ww(i)
    END DO
    DEALLOCATE (ww,STAT=err)
    IF (err /= 0 ) CALL stop_memory ( 'init_spline','ww')

  END SUBROUTINE init_spline

  !******************************************************************************
  !!****
  !!****f* splines/splineh [1.0] *
  !!
  !!   NAME
  !!     splineh
  !!
  !!   FUNCTION
  !!     calculates the spline interpolated value at a given point
  !!     and the first derivative WITHOUT checks
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

  FUNCTION splineh ( spl, xxi, y1 )
    ! Return value
    TYPE(spline_data_type), POINTER          :: spl
    REAL(KIND=dp), INTENT(IN)                :: xxi
    REAL(KIND=dp), INTENT(OUT)               :: y1
    REAL(KIND=dp)                            :: splineh

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: a, b, h26, invh, y2hi, y2lo, &
         yhi, ylo, xx, x4

    xx = 1.0_dp/xxi
    h26 = spl%h26
    invh = spl%invh
    i = (xx-spl%x1)*invh + 1

    a = (spl%x1-xx)*invh + REAL(i,kind=dp)
    b = 1.0_dp - a
    ylo = spl%y(i)
    yhi = spl%y(i+1)
    y2lo = spl%y2(i)
    y2hi = spl%y2(i+1)
    splineh = a*ylo + b*yhi - ((a+1)*y2lo+(b+1)*y2hi)*a*b*h26
    y1 = invh*((yhi-ylo) + ((1.-3*a*a)*y2lo+(3*b*b-1.)*y2hi)*h26)
    ! take into account the mapping w.r.t. 1/r^2
    x4 = xx * xx
    y1 = - y1 * x4
  END FUNCTION splineh

END MODULE splines
