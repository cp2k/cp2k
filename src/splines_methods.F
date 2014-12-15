!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines for handling splines
!> \par History
!>      2001-09-21-HAF added this doc entry and changed formatting
!> \author various
! *****************************************************************************
MODULE splines_methods

  
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE splines_types,                   ONLY: spline_data_p_type,&
                                             spline_data_type,&
                                             spline_factor_type
  USE termination,                     ONLY: stop_memory
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'splines_methods'

  PUBLIC :: init_splinexy    ! allocates x and y vectors for splines
  PUBLIC :: init_spline      ! generate table for spline (allocates y2)
  PUBLIC :: potential_s      ! return value of spline and 1. derivative
                             ! without checks (fast routine for pair_potential)
  PUBLIC :: spline_value     ! return value of spline and 1. derivative
                             ! without check (without assumption of 1/x^2 grid)

CONTAINS

! *****************************************************************************
!> \brief allocates storage for function table to be interpolated
!>      both x and y are allocated
!> \param spl spline_data structure to be initialized
!> \param nn integer number of datapoints, that the function table will hold
!> \par History
!>      2001-09-21-HAF added this doc entry and changed formatting
!> \author unknown
! *****************************************************************************
  SUBROUTINE init_splinexy ( spl, nn )

    TYPE(spline_data_type), POINTER          :: spl
    INTEGER, INTENT(IN)                      :: nn

    CHARACTER(len=*), PARAMETER :: routineN = 'init_splinexy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: err

    spl % n = nn

    IF ( ASSOCIATED ( spl % y ) ) THEN
       DEALLOCATE ( spl % y, STAT = err )
       IF (err /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"spl%y")
       NULLIFY ( spl % y )
    END IF

    IF ( ASSOCIATED ( spl % y2 ) ) THEN
       DEALLOCATE ( spl % y2, STAT = err )
       IF (err /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"spl%y2")
       NULLIFY ( spl % y2 )
    END IF

    ALLOCATE ( spl % y ( 1:nn ), STAT = err )
    IF (err /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                   "spl%y",dp_size*nn)

    ALLOCATE ( spl % y2 ( 1:nn ), STAT = err )
    IF (err /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                   "spl%y2",dp_size*nn)

  END SUBROUTINE init_splinexy

! *****************************************************************************
!> \brief allocates storage for y2 table
!>      calculates y2 table and other spline parameters
!> \param spl spline_data structure to be initialized
!>                       spl%y() must hold the function values
!>                       spl%x() must hold the absissa values in increasing
!>                       order OR if dx (below) is given, spl%x(1) must hold
!>                       the starting (left-most) point.
!> \param dx x(i) are assumed to be x(1)+dx*(i-1)
!>                       (spline evaluations will also be faster)
!>      y1a : (OPTIONAL) if present, the 1-deriv of the left endpoint
!>                       if not present, natural spline condition at this end
!>                       (2-deriv == 0)
!>      y1b : (OPTIONAL) if present, the 1-deriv of the right endpoint
!>                       if not present, natural spline condition at this end
!>                       (2-deriv == 0)
!> \param y1a ...
!> \param y1b ...
!> \par Examples
!>      CALL init_spline(spline,dx=0.1_dp)
!>      CALL init_spline(spline,y1b=0.0_dp)
!>      CALL init_spline(spline,0.1_dp,0.0_dp,0.0_dp)
!> \par History
!>      2001-09-21-HAF added this doc entry and changed formatting
!>      2001-09-24-HAF changed interface and re-written
!> \author unknown
!> \note
!>      if dx is given, the x array will be used as y2 array instead of
!>      allocating a new array. (y2 will become x, and x will be nullified)
! *****************************************************************************
  SUBROUTINE init_spline(spl, dx, y1a, y1b)

    TYPE(spline_data_type), POINTER          :: spl
    REAL(KIND=dp), INTENT(IN)                :: dx
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: y1a, y1b

    CHARACTER(len=*), PARAMETER :: routineN = 'init_spline', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: err, i, n
    REAL(KIND=dp)                            :: p, s
    REAL(KIND=dp), POINTER                   :: ww( : )

    n        = spl%n
    spl%xn   = spl%x1+ (n-1)*dx
    spl%h    = dx
    spl%invh = 1.0_dp/dx
    spl%h26  = dx**2 / 6.0_dp
    ALLOCATE (ww(1:n),STAT=err)
    IF (err /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                   "ww",dp_size*n)
    IF (PRESENT(y1a)) THEN
       spl%y2(1) = -0.5_dp
       ww(1) = 3.0_dp * ((spl%y(2)-spl%y(1))/dx - y1a)/dx
    ELSE
       spl%y2(1) = 0.0_dp
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
    IF (err /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"ww")

  END SUBROUTINE init_spline

! *****************************************************************************
!> \brief calculates the potential interpolated with splines value at a given point
!>      and the first derivative. Checks included to avoid just segfaulting!!
!> \param spl_p spline_data structure
!> \param xxi absissa value
!> \param y1 1. derivative at xx
!> \param spl_f ...
!> \param logger ...
!> \retval potential_s ...
!> \par Output
!>      spline interpolated value at xx
!> \par History
!>      2001-09-25-HAF added this doc entry and changed formatting
!> \author unknown
!> \note
!>      the spline MUST have uniform x values and xx MUST be
!>      in the interpolation interval. No checks are done to ensure
!>      either condition.
! *****************************************************************************
  FUNCTION potential_s (spl_p, xxi, y1, spl_f, logger )
    TYPE(spline_data_p_type), DIMENSION(:), &
      POINTER                                :: spl_p
    REAL(KIND=dp), INTENT(IN)                :: xxi
    REAL(KIND=dp), INTENT(OUT)               :: y1
    TYPE(spline_factor_type), POINTER        :: spl_f
    TYPE(cp_logger_type), POINTER            :: logger
    REAL(KIND=dp)                            :: potential_s

    REAL(KIND=dp), PARAMETER                 :: f13 = 1.0_dp/3.0_dp

    INTEGER                                  :: i, output_unit
    REAL(KIND=dp)                            :: a, b, h26, invh, x4, xx, xx0, &
                                                y2hi, y2lo, yhi, ylo, yy

    xx0  = 1.0_dp/xxi
    xx   = spl_f%rscale(1)*xx0
    x4   = xx * xx
    h26  = spl_p(1)%spline_data%h26
    invh = spl_p(1)%spline_data%invh
    IF (xx >= spl_p(1)%spline_data%xn) THEN
       ! In case the value is not on the spline let's print a warning and give the value
       ! for the smaller point available in the spline..
       ! This should happen in very few cases though..
       output_unit= cp_logger_get_default_unit_nr(logger)
       yy = spl_p(1)%spline_data%xn-spl_p(1)%spline_data%h
       WRITE(output_unit,FMT='(/,80("*"),/,"*",1X,"Value of r in Input =",F11.6,'//&
            '" not in the spline range. Using =",F11.6,T80,"*",/,80("*"))')SQRT(1.0_dp/xx), SQRT(1.0_dp/yy)
       xx = yy
    END IF
    i = INT((xx-spl_p(1)%spline_data%x1)*invh + 1)
    a = (spl_p(1)%spline_data%x1-xx)*invh + REAL(i,kind=dp)
    b = 1.0_dp - a

    ylo  = spl_p(1)%spline_data%y(i)
    yhi  = spl_p(1)%spline_data%y(i+1)
    y2lo = spl_p(1)%spline_data%y2(i)
    y2hi = spl_p(1)%spline_data%y2(i+1)
    potential_s = (a*ylo + b*yhi - ((a+1.0_dp)*y2lo+(b+1.0_dp)*y2hi)*a*b*h26)*spl_f%fscale(1)
    y1 = invh*((yhi-ylo) + ((f13-a*a)*y2lo-(f13-b*b)*y2hi)*3.0_dp*h26)
    y1 = 2.0_dp * y1 * x4 * spl_f%dscale(1)

    potential_s = potential_s + spl_f%cutoff
  END FUNCTION potential_s

! *****************************************************************************
!> \brief calculates the spline value at a given point
!>        (and possibly the first derivative) WITHOUT checks
!>        and without any funny scaling business, or weird
!>        1/x^2 grid assumptions
!>
!> \param spl ...
!> \param xx ...
!> \param y1 ...
!> \retval spline_value ...
!> \author HAF
! *****************************************************************************
  FUNCTION spline_value ( spl, xx, y1 )
    ! Return value
    TYPE(spline_data_type), POINTER          :: spl
    REAL(KIND=dp), INTENT(IN)                :: xx
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: y1
    REAL(KIND=dp)                            :: spline_value

    REAL(KIND=dp), PARAMETER                 :: f13 = 1.0_dp/3.0_dp

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: a, b, h26, invh, y2hi, y2lo, &
                                                yhi, ylo

    h26 = spl%h26
    invh = spl%invh
    i = INT((xx-spl%x1)*invh + 1)

    a = (spl%x1-xx)*invh + REAL(i,kind=dp)
    b = 1.0_dp - a
    ylo = spl%y(i)
    yhi = spl%y(i+1)
    y2lo = spl%y2(i)
    y2hi = spl%y2(i+1)
    spline_value = a*ylo + b*yhi - ((a+1.0_dp)*y2lo+(b+1.0_dp)*y2hi)*a*b*h26
    IF (PRESENT(y1)) y1 = invh*((yhi-ylo) + &
                           ((f13-a*a)*y2lo-(f13-b*b)*y2hi)*3.0_dp*h26)
  END FUNCTION spline_value

END MODULE splines_methods
