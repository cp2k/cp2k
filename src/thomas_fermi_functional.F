!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/thomas_fermi_functional [1.0] *
!!
!!   NAME
!!     thomas_fermi_functional
!!
!!   FUNCTION
!!     Calculate the Thomas-Fermi kinetic energy functional
!!
!!   NOTES
!!     Order of derivatives is: LDA 0; 1; 2; 3;
!!                              LSD 0; a  b; aa bb; aaa bbb;
!!
!!   AUTHOR
!!     JGH (18.02.2002)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!******************************************************************************

MODULE thomas_fermi_functional

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
                          f43 = 4._dbl*f13, &
                          f53 = 5._dbl*f13

  PUBLIC :: thomas_fermi

  REAL(dbl) :: cf, flda, flsd
  REAL(dbl) :: eps_rho

! *****************************************************************************

  INTERFACE thomas_fermi
    MODULE PROCEDURE thomas_fermi_lda, thomas_fermi_lsd, &
                     thomas_fermi_init, thomas_fermi_info
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE thomas_fermi_init ( cutoff )

    REAL(dbl), INTENT(IN) :: cutoff

    eps_rho = cutoff
    CALL set_util ( cutoff )

    cf = 0.3_dbl*(3._dbl*pi*pi)**f23
    flda = cf
    flsd = flda * 2._dbl**f23

  END SUBROUTINE thomas_fermi_init

! *****************************************************************************

  SUBROUTINE thomas_fermi_info ( reference, shortform )

    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: reference
    CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: shortform

    IF ( PRESENT ( reference ) ) THEN
      reference = "Thomas-Fermi kinetic energy functional: see Parr and Yang"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Thomas-Fermi kinetic energy functional"
    END IF

  END SUBROUTINE thomas_fermi_info

! *****************************************************************************

  SUBROUTINE thomas_fermi_lda ( rho, pot, order )

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
    IF ( ABS ( order ) > 3 ) CALL stop_program ( "thomas_fermi_lda", &
       "order of derivative too high" )
    IF ( mt > mp ) CALL stop_program ( "thomas_fermi_lda", &
       "Size of POT array too small" )

    n = SIZE ( rho )
    ALLOCATE ( r13(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "thomas_fermi_lda", "r13", n )
    
    m = 1
    CALL calc_rho13 ( rho, r13 )
    IF ( calc(0) ) THEN
       CALL thomas_fermi_lda_0 ( rho, r13, pot(:,1) )
    END IF
    IF ( calc(1) ) THEN
       IF ( allcalc ) m = 2
       CALL thomas_fermi_lda_1 ( rho, r13, pot(:,m) )
    END IF
    IF ( calc(2) ) THEN
       IF ( allcalc ) m = 3
       CALL thomas_fermi_lda_2 ( rho, r13, pot(:,m) )
    END IF
    IF ( calc(3) ) THEN
       IF ( allcalc ) m = 4
       CALL thomas_fermi_lda_3 ( rho, r13, pot(:,m) )
    END IF

    DEALLOCATE ( r13, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "thomas_fermi_lda", "rs" )

  END SUBROUTINE thomas_fermi_lda

! *****************************************************************************

  SUBROUTINE thomas_fermi_lsd ( rhoa, rhob, pot, order )

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
          mt = 5
        CASE (3)
          mt = 7
        CASE DEFAULT
          mt = mp+100
        END SELECT
      ELSE
        calc(-order) = .TRUE.
        allcalc = .FALSE.
        mt = 2
      END IF
    ELSE
       calc(0) = .TRUE.
       calc(1) = .TRUE.
       mt = 3
    END IF
    IF ( ABS ( order ) > 3 ) CALL stop_program ( "thomas_fermi_lsd", &
       "order of derivative too high" )
    IF ( mt > mp ) CALL stop_program ( "thomas_fermi_lsd", &
       "Size of POT array too small" )

    n = SIZE ( rhoa )
    ALLOCATE ( r13a(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "thomas_fermi_lsd", "r13a", n )
    ALLOCATE ( r13b(n), STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "thomas_fermi_lsd", "r13b", n )

    m = 1
    CALL calc_rho13 ( rhoa, r13a )
    CALL calc_rho13 ( rhob, r13b )
    IF ( calc(0) ) THEN
       CALL thomas_fermi_lsd_0 ( rhoa, rhob, r13a, r13b, pot(:,1) )
    END IF
    IF ( calc(1) ) THEN
       IF ( allcalc ) m = 2
       CALL thomas_fermi_lsd_1 ( rhoa, rhob, r13a, r13b, pot(:,m:m+1) )
    END IF
    IF ( calc(2) ) THEN
       IF ( allcalc ) m = 4
       CALL thomas_fermi_lsd_2 ( rhoa, rhob, r13a, r13b, pot(:,m:m+1) )
    END IF
    IF ( calc(3) ) THEN
       IF ( allcalc ) m = 6
       CALL thomas_fermi_lsd_3 ( rhoa, rhob, r13a, r13b, pot(:,m:m+1) )
    END IF

    DEALLOCATE ( r13a, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "thomas_fermi_lsd", "r13a" )
    DEALLOCATE ( r13b, STAT=error )
    IF ( error /= 0 ) CALL stop_memory ( "thomas_fermi_lsd", "r13b" )

  END SUBROUTINE thomas_fermi_lsd

! *****************************************************************************

  SUBROUTINE thomas_fermi_lda_0 ( rho, r13, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, r13
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "thomas_fermi_lda_0", "Inconsistent array sizes" )

    DO ip = 1, npoints

      IF ( rho(ip) < eps_rho ) THEN

         pot(ip) = 0.0_dbl

      ELSE

         pot(ip) = flda*r13(ip)*r13(ip)*rho(ip)

      END IF

    END DO

  END SUBROUTINE thomas_fermi_lda_0

! *****************************************************************************

  SUBROUTINE thomas_fermi_lda_1 ( rho, r13, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, r13
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL ( dbl ) :: f
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "thomas_fermi_lda_1", "Inconsistent array sizes" )

    f = f53 * flda

    DO ip = 1, npoints

      IF ( rho(ip) < eps_rho ) THEN

         pot(ip) = 0.0_dbl

      ELSE

         pot(ip) = f * r13(ip)*r13(ip)

      END IF

    END DO

  END SUBROUTINE thomas_fermi_lda_1

! *****************************************************************************

  SUBROUTINE thomas_fermi_lda_2 ( rho, r13, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, r13
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL ( dbl ) :: f
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "thomas_fermi_lda_2", "Inconsistent array sizes" )

    f = f23 * f53 * flda

    DO ip = 1, npoints

      IF ( rho(ip) < eps_rho ) THEN

         pot(ip) = 0.0_dbl

      ELSE

         pot(ip) = f / r13(ip)

      END IF

    END DO

  END SUBROUTINE thomas_fermi_lda_2

! *****************************************************************************

  SUBROUTINE thomas_fermi_lda_3 ( rho, r13, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho, r13
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL ( dbl ) :: f
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rho )
    IF ( npoints > SIZE(pot) ) &
       CALL stop_program ( "thomas_fermi_lda_3", "Inconsistent array sizes" )

    f = -f13 * f23 * f53 * flda

    DO ip = 1, npoints

      IF ( rho(ip) < eps_rho ) THEN

         pot(ip) = 0.0_dbl

      ELSE

         pot(ip) = f / ( r13(ip) * rho(ip) )

      END IF

    END DO

  END SUBROUTINE thomas_fermi_lda_3

! *****************************************************************************

  SUBROUTINE thomas_fermi_lsd_0 ( rhoa, rhob, r13a, r13b, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, r13a, r13b
    REAL(dbl), DIMENSION(:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: f
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "thomas_fermi_lsd_0", "Inconsistent array sizes" )

    DO ip = 1, npoints

      IF ( rhoa(ip) < eps_rho ) THEN
         pot(ip) = 0.0_dbl
      ELSE
         pot(ip) = flsd * r13a(ip) * r13a(ip) * rhoa(ip)
      END IF
      IF ( rhob(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + flsd * r13b(ip) * r13b(ip) * rhob(ip)
      END IF

    END DO

  END SUBROUTINE thomas_fermi_lsd_0

! *****************************************************************************

  SUBROUTINE thomas_fermi_lsd_1 ( rhoa, rhob, r13a, r13b, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, r13a, r13b
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: f
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "thomas_fermi_lsd_1", "Inconsistent array sizes" )

    f = f53 * flsd

    DO ip = 1, npoints

      IF ( rhoa(ip) < eps_rho ) THEN
         pot(ip,1) = 0.0_dbl
      ELSE
         pot(ip,1) = f * r13a(ip) * r13a(ip)
      END IF
      IF ( rhob(ip) < eps_rho ) THEN
         pot(ip,2) = 0.0_dbl
      ELSE
         pot(ip,2) = f * r13b(ip) * r13b(ip)
      END IF

    END DO

  END SUBROUTINE thomas_fermi_lsd_1

! *****************************************************************************

  SUBROUTINE thomas_fermi_lsd_2 ( rhoa, rhob, r13a, r13b, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, r13a, r13b
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: f
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "thomas_fermi_lsd_2", "Inconsistent array sizes" )

    f = f23 * f53 * flsd

    DO ip = 1, npoints

      IF ( rhoa(ip) < eps_rho ) THEN
         pot(ip,1) = 0.0_dbl
      ELSE
         pot(ip,1) = f / r13a(ip)
      END IF
      IF ( rhob(ip) < eps_rho ) THEN
         pot(ip,2) = 0.0_dbl
      ELSE
         pot(ip,2) = f / r13b(ip)
      END IF

    END DO

  END SUBROUTINE thomas_fermi_lsd_2

! *****************************************************************************

  SUBROUTINE thomas_fermi_lsd_3 ( rhoa, rhob, r13a, r13b, pot )

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob, r13a, r13b
    REAL(dbl), DIMENSION(:,:), INTENT(OUT)  :: pot

    INTEGER :: ip, npoints
    REAL(dbl) :: f
!   ---------------------------------------------------------------------------

    npoints = SIZE ( rhoa )
    IF ( npoints > SIZE(pot,1) ) &
       CALL stop_program ( "thomas_fermi_lsd_3", "Inconsistent array sizes" )

    f = -f13 * f23 * f53 * flsd

    DO ip = 1, npoints

      IF ( rhoa(ip) < eps_rho ) THEN
         pot(ip,1) = 0.0_dbl
      ELSE
         pot(ip,1) = f / ( r13a(ip)*rhoa(ip))
      END IF
      IF ( rhob(ip) < eps_rho ) THEN
         pot(ip,2) = 0.0_dbl
      ELSE
         pot(ip,2) = f / ( r13b(ip)*rhob(ip))
      END IF

    END DO

  END SUBROUTINE thomas_fermi_lsd_3

! *****************************************************************************

END MODULE thomas_fermi_functional

! *****************************************************************************
