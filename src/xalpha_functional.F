!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/xalpha_functional [1.0] *
!!
!!   NAME
!!     xalpha_functional
!!
!!   FUNCTION
!!     Calculate the local exchange functional
!!
!!   AUTHOR
!!     JGH (17.02.2002)
!!
!!   MODIFICATION HISTORY
!!     - none
!!
!!   SOURCE
!******************************************************************************

MODULE xalpha_functional

! *****************************************************************************

  USE kinds, ONLY: dbl
  USE functionals_utilities, ONLY: set_util, calc_rho13
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

  SUBROUTINE xalpha_lda ( rho, pot, order )

!   ---------------------------------------------------------------------------

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

!   *** Local variables ***

    LOGICAL   :: calc(0:4), allcalc
    REAL(dbl), DIMENSION(:), ALLOCATABLE :: r13
    INTEGER :: n, error, m, mp, mt

!   ---------------------------------------------------------------------------

    calc = .FALSE.
    mp = SIZE ( pot, 2 )
    allcalc = .TRUE.
    IF ( PRESENT ( order ) ) THEN
      IF ( order >= 0 ) THEN
        calc(0:order) = .TRUE.
        mt = order+1
      ELSE
        calc(-order) = .TRUE.
        allcalc = .FALSE.
        mt = 1
      END IF
    ELSE
       calc(0) = .TRUE.
       calc(1) = .TRUE.
       mt = 2
    END IF
    IF ( mt > mp ) CALL stop_program ( "xalpha_lda", &
       "Size of POT array too small" )

    n = SIZE ( rho )
    ALLOCATE ( r13(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "xalpha_lda", "r13", n )
    
    m = 1
    CALL calc_rho13 ( rho, r13 )
    IF ( calc(0) ) THEN
       CALL xalpha_lda_0 ( rho, r13, pot(:,1) )
    END IF
    IF ( calc(1) ) THEN
       IF ( allcalc ) m = 2
       CALL xalpha_lda_1 ( rho, r13, pot(:,m) )
    END IF
    IF ( calc(2) ) THEN
       IF ( allcalc ) m = 3
       CALL xalpha_lda_2 ( rho, r13, pot(:,m) )
    END IF
    IF ( calc(3) ) THEN
       IF ( allcalc ) m = 4
       CALL xalpha_lda_3 ( rho, r13, pot(:,m) )
    END IF

    DEALLOCATE ( r13, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "xalpha_lda", "rs" )

  END SUBROUTINE xalpha_lda

! *****************************************************************************

  SUBROUTINE xalpha_lsd ( rhoa, rhob, pot, order )

!   ---------------------------------------------------------------------------

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: pot
    INTEGER, INTENT(IN), OPTIONAL :: order

!   *** Local variables ***

    LOGICAL   :: calc(0:4), allcalc
    REAL(dbl), DIMENSION(:), ALLOCATABLE :: r13a, r13b
    INTEGER :: n, error, m, mp, mt

!   ---------------------------------------------------------------------------

    calc = .FALSE.
    mp = SIZE ( pot, 2 )
    allcalc = .TRUE.
    IF ( PRESENT ( order ) ) THEN
      IF ( order >= 0 ) THEN
        calc(0:order) = .TRUE.
        SELECT CASE (order)
        CASE (0)
          mt = 1
        CASE (1)
          mt = 3
        CASE (2)
          mt = 6
        CASE (3)
          mt = 10
        CASE DEFAULT
          mt = mp+100
        END SELECT
      ELSE
        calc(-order) = .TRUE.
        allcalc = .FALSE.
        mt = -order+1
      END IF
    ELSE
       calc(0) = .TRUE.
       calc(1) = .TRUE.
       mt = 3
    END IF
    IF ( mt > mp ) CALL stop_program ( "xalpha_lsd", &
       "Size of POT array too small" )

    n = SIZE ( rhoa )
    ALLOCATE ( r13a(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "xalpha_lsd", "r13a", n )
    ALLOCATE ( r13b(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "xalpha_lsd", "r13b", n )

    m = 1
    CALL calc_rho13 ( rhoa, r13a )
    CALL calc_rho13 ( rhob, r13b )
    IF ( calc(0) ) THEN
       CALL xalpha_lsd_0 ( rhoa, rhob, r13a, r13b, pot(:,1) )
    END IF
    IF ( calc(1) ) THEN
       IF ( allcalc ) m = 2
       CALL xalpha_lsd_1 ( rhoa, rhob, r13a, r13b, pot(:,m:m+1) )
    END IF
    IF ( calc(2) ) THEN
       IF ( allcalc ) m = 4
       CALL xalpha_lsd_2 ( rhoa, rhob, r13a, r13b, pot(:,m:m+2) )
    END IF
    IF ( calc(3) ) THEN
       IF ( allcalc ) m = 7
       CALL xalpha_lsd_3 ( rhoa, rhob, r13a, r13b, pot(:,m:m+3) )
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

    DO ip = 1, npoints

      IF ( rho(ip) < eps_rho ) THEN

         pot(ip) = 0.0_dbl

      ELSE

         pot(ip) = flda*r13(ip)*rho(ip)

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

    DO ip = 1, npoints

      IF ( rho(ip) < eps_rho ) THEN

         pot(ip) = 0.0_dbl

      ELSE

         pot(ip) = f * r13(ip)

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

    DO ip = 1, npoints

      IF ( rho(ip) < eps_rho ) THEN

         pot(ip) = 0.0_dbl

      ELSE

         pot(ip) = f * r13(ip) / rho(ip)

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
       CALL stop_program ( "xalpha_lda_2", "Inconsistent array sizes" )

    f = -f23 * f13 * f43 * flda

    DO ip = 1, npoints

      IF ( rho(ip) < eps_rho ) THEN

         pot(ip) = 0.0_dbl

      ELSE

         pot(ip) = f * r13(ip) / ( rho(ip) * rho(ip) )

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

    DO ip = 1, npoints

      IF ( rhoa(ip) < eps_rho ) THEN
         pot(ip) = 0.0_dbl
      ELSE
         pot(ip) = flsd * r13a(ip) * rhoa(ip)
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

    DO ip = 1, npoints

      IF ( rhoa(ip) < eps_rho ) THEN
         pot(ip,1) = 0.0_dbl
      ELSE
         pot(ip,1) = f * r13a(ip)
      END IF
      IF ( rhob(ip) < eps_rho ) THEN
         pot(ip,2) = 0.0_dbl
      ELSE
         pot(ip,2) = f * r13b(ip)
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

    DO ip = 1, npoints

      pot(ip,2) = 0._dbl
      IF ( rhoa(ip) < eps_rho ) THEN
         pot(ip,1) = 0.0_dbl
      ELSE
         pot(ip,1) = f * r13a(ip)/rhoa(ip)
      END IF
      IF ( rhob(ip) < eps_rho ) THEN
         pot(ip,3) = 0.0_dbl
      ELSE
         pot(ip,3) = f * r13b(ip)/rhob(ip)
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

    DO ip = 1, npoints

      pot(ip,2) = 0._dbl
      pot(ip,3) = 0._dbl
      IF ( rhoa(ip) < eps_rho ) THEN
         pot(ip,1) = 0.0_dbl
      ELSE
         pot(ip,1) = f * r13a(ip)/(rhoa(ip)*rhoa(ip))
      END IF
      IF ( rhob(ip) < eps_rho ) THEN
         pot(ip,4) = 0.0_dbl
      ELSE
         pot(ip,4) = f * r13b(ip)/(rhob(ip)*rhob(ip))
      END IF

    END DO

  END SUBROUTINE xalpha_lsd_3

! *****************************************************************************

END MODULE xalpha_functional

! *****************************************************************************
