!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/xc_xalpha [1.0] *
!!
!!   NAME
!!     xc_xalpha
!!
!!   FUNCTION
!!     Calculate the local exchange functional
!!
!!   NOTES
!!     Order of derivatives is: LDA 0; 1; 2; 3; 
!!                              LSD 0; a  b; aa bb; aaa bbb;
!!
!!   AUTHOR
!!     JGH (17.02.2002)
!!
!!   MODIFICATION HISTORY
!!     JGH (26.02.2003) : OpenMP enabled
!!
!!*** *************************************************************************

MODULE xc_xalpha

! *****************************************************************************

  USE kinds, ONLY: dbl
  USE xc_derivative_types, ONLY: xc_derivative_set_type
  USE xc_functionals_utilities, ONLY: set_util, calc_rho13, setup_calculation
  USE termination, ONLY: stop_memory, stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  REAL(dbl), PARAMETER :: pi = 3.14159265358979323846264338_dbl
  REAL(dbl), PARAMETER :: f13 = 1._dbl/3._dbl, &
                          f23 = 2._dbl*f13, &
                          f43 = 4._dbl*f13

  PUBLIC :: xalpha

  REAL(dbl) :: xparam, flda, flsd
  REAL(dbl) :: eps_rho

! *****************************************************************************

  INTERFACE xalpha
    MODULE PROCEDURE xalpha_lda, xalpha_lsd, xalpha_init, xalpha_info
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE xalpha_init ( cutoff, parameter )

    REAL(dbl), INTENT(IN) :: cutoff
    REAL(dbl), INTENT(IN), OPTIONAL :: parameter

    eps_rho = cutoff
    CALL set_util ( cutoff )
    IF ( PRESENT ( parameter ) ) THEN
      xparam = parameter
    ELSE
      xparam = 2._dbl / 3._dbl
    END IF

    flda = -9._dbl/8._dbl * xparam * (3._dbl/pi)**f13
    flsd = flda * 2._dbl**f13

  END SUBROUTINE xalpha_init

! *****************************************************************************

  SUBROUTINE xalpha_info ( reference, shortform )

    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: reference
    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: shortform

    IF ( PRESENT ( reference ) ) THEN
       WRITE (reference,'(A,F8.4)') &
         "Dirac/Slater local exchange; parameter=",xparam
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       WRITE (shortform,'(A,F8.4)') "Dirac/Slater exchange",xparam
    END IF

  END SUBROUTINE xalpha_info

! *****************************************************************************

  SUBROUTINE xalpha_lda ( derivative_set, rho, pot, order )

!   ---------------------------------------------------------------------------

    TYPE(xc_derivative_set_type), POINTER :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

!   *** Local variables ***

    LOGICAL   :: calc(0:4)
    REAL(dbl), DIMENSION(:), ALLOCATABLE :: r13
    INTEGER :: n, error, m, mp, mt, mx(0:3,2)

!   ---------------------------------------------------------------------------

    CALL setup_calculation(order,mx,calc,100)
    mp = SIZE ( pot, 2 )
    m = ABS(order)
    mt = mx(m,2)
    IF ( mt > mp ) CALL stop_program ( "xalpha_lda", &
       "Size of POT array too small" )

    n = SIZE ( rho )
    ALLOCATE ( r13(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "xalpha_lda", "r13", n )
    
    CALL calc_rho13 ( rho, r13 )
    IF ( calc(0) ) THEN
       CALL xalpha_lda_0 ( rho, r13, pot(:,1) )
    END IF
    IF ( calc(1) ) THEN
       m = mx(1,1)
       CALL xalpha_lda_1 ( rho, r13, pot(:,m) )
    END IF
    IF ( calc(2) ) THEN
       m = mx(2,1)
       CALL xalpha_lda_2 ( rho, r13, pot(:,m) )
    END IF
    IF ( calc(3) ) THEN
       m = mx(3,1)
       CALL xalpha_lda_3 ( rho, r13, pot(:,m) )
    END IF

    DEALLOCATE ( r13, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "xalpha_lda", "rs" )

  END SUBROUTINE xalpha_lda

! *****************************************************************************

  SUBROUTINE xalpha_lsd ( derivative_set, rhoa, rhob, pot, order )

!   ---------------------------------------------------------------------------

    TYPE(xc_derivative_set_type), POINTER :: derivative_set
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

!   *** Local variables ***

    LOGICAL   :: calc(0:4)
    REAL(dbl), DIMENSION(:), ALLOCATABLE :: r13a, r13b
    INTEGER :: n, error, m, mp, mt, mx(0:3,2),m1,m2

!   ---------------------------------------------------------------------------

    CALL setup_calculation(order,mx,calc,200)
    mp = SIZE ( pot, 2 )
    m = ABS(order)
    mt = mx(m,2)
    IF ( mt > mp ) CALL stop_program ( "xalpha_lsd", &
       "Size of POT array too small" )

    n = SIZE ( rhoa )
    ALLOCATE ( r13a(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "xalpha_lsd", "r13a", n )
    ALLOCATE ( r13b(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "xalpha_lsd", "r13b", n )

    CALL calc_rho13 ( rhoa, r13a )
    CALL calc_rho13 ( rhob, r13b )
    IF ( calc(0) ) THEN
       CALL xalpha_lsd_0 ( rhoa, rhob, r13a, r13b, pot(:,1) )
    END IF
    IF ( calc(1) ) THEN
       m1 = mx(1,1)
       m2 = mx(1,2)
       CALL xalpha_lsd_1 ( rhoa, rhob, r13a, r13b, pot(:,m1:m2) )
    END IF
    IF ( calc(2) ) THEN
       m1 = mx(2,1)
       m2 = mx(2,2)
       CALL xalpha_lsd_2 ( rhoa, rhob, r13a, r13b, pot(:,m1:m2) )
    END IF
    IF ( calc(3) ) THEN
       m1 = mx(3,1)
       m2 = mx(3,2)
       CALL xalpha_lsd_3 ( rhoa, rhob, r13a, r13b, pot(:,m1:m2) )
    END IF

    DEALLOCATE ( r13a, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "xalpha_lsd", "r13a" )
    DEALLOCATE ( r13b, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "xalpha_lsd", "r13b" )

  END SUBROUTINE xalpha_lsd

! *****************************************************************************

  SUBROUTINE xalpha_lda_0 ( rho, r13, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, r13
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "xalpha_lda_0", "Inconsistent array sizes" )

!$omp parallel do private(ip)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         pot(ip) = pot(ip) + flda*r13(ip)*rho(ip)

      END IF

    END DO

  END SUBROUTINE xalpha_lda_0

! *****************************************************************************

  SUBROUTINE xalpha_lda_1 ( rho, r13, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, r13
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL ( dbl ) :: f
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "xalpha_lda_1", "Inconsistent array sizes" )

    f = f43 * flda

!$omp parallel do private(ip)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         pot(ip) = pot(ip) + f * r13(ip)

      END IF

    END DO

  END SUBROUTINE xalpha_lda_1

! *****************************************************************************

  SUBROUTINE xalpha_lda_2 ( rho, r13, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, r13
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL ( dbl ) :: f
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "xalpha_lda_2", "Inconsistent array sizes" )

    f = f13 * f43 * flda

!$omp parallel do private(ip)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         pot(ip) = pot(ip) + f * r13(ip) / rho(ip)

      END IF

    END DO

  END SUBROUTINE xalpha_lda_2

! *****************************************************************************

  SUBROUTINE xalpha_lda_3 ( rho, r13, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, r13
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL ( dbl ) :: f
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "xalpha_lda_3", "Inconsistent array sizes" )

    f = -f23 * f13 * f43 * flda

!$omp parallel do private(ip)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         pot(ip) = pot(ip) + f * r13(ip) / ( rho(ip) * rho(ip) )

      END IF

    END DO

  END SUBROUTINE xalpha_lda_3

! *****************************************************************************

  SUBROUTINE xalpha_lsd_0 ( rhoa, rhob, r13a, r13b, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, r13a, r13b
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: f
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "xalpha_lsd_0", "Inconsistent array sizes" )

!$omp parallel do private(ip)
    DO ip = 1, npoints

      IF ( rhoa(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + flsd * r13a(ip) * rhoa(ip)
      END IF
      IF ( rhob(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + flsd * r13b(ip) * rhob(ip)
      END IF

    END DO

  END SUBROUTINE xalpha_lsd_0

! *****************************************************************************

  SUBROUTINE xalpha_lsd_1 ( rhoa, rhob, r13a, r13b, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, r13a, r13b
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: f
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "xalpha_lsd_1", "Inconsistent array sizes" )

    f = f43 * flsd

!$omp parallel do private(ip)
    DO ip = 1, npoints

      IF ( rhoa(ip) > eps_rho ) THEN
         pot(ip,1) = pot(ip,1) + f * r13a(ip)
      END IF
      IF ( rhob(ip) > eps_rho ) THEN
         pot(ip,2) = pot(ip,2) + f * r13b(ip)
      END IF

    END DO

  END SUBROUTINE xalpha_lsd_1

! *****************************************************************************

  SUBROUTINE xalpha_lsd_2 ( rhoa, rhob, r13a, r13b, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, r13a, r13b
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: f
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "xalpha_lsd_2", "Inconsistent array sizes" )

    f = f13 * f43 * flsd

!$omp parallel do private(ip)
    DO ip = 1, npoints

      IF ( rhoa(ip) > eps_rho ) THEN
         pot(ip,1) = pot(ip,1) + f * r13a(ip)/rhoa(ip)
      END IF
      IF ( rhob(ip) > eps_rho ) THEN
         pot(ip,2) = pot(ip,2) + f * r13b(ip)/rhob(ip)
      END IF

    END DO

  END SUBROUTINE xalpha_lsd_2

! *****************************************************************************

  SUBROUTINE xalpha_lsd_3 ( rhoa, rhob, r13a, r13b, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, r13a, r13b
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: f
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "xalpha_lsd_3", "Inconsistent array sizes" )

    f = -f23 * f13 * f43 * flsd

!$omp parallel do private(ip)
    DO ip = 1, npoints

      IF ( rhoa(ip) > eps_rho ) THEN
         pot(ip,1) = pot(ip,1) + f * r13a(ip)/(rhoa(ip)*rhoa(ip))
      END IF
      IF ( rhob(ip) > eps_rho ) THEN
         pot(ip,2) = pot(ip,2) + f * r13b(ip)/(rhob(ip)*rhob(ip))
      END IF

    END DO

  END SUBROUTINE xalpha_lsd_3

! *****************************************************************************

END MODULE xc_xalpha

! *****************************************************************************
